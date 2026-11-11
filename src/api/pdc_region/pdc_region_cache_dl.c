#include <time.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <mpi.h>

#include "pdc_utlist.h"
#include "pdc_config.h"
#include "pdc_id_pkg.h"
#include "pdc_obj.h"
#include "pdc_malloc.h"
#include "pdc_region.h"
#include "pdc_region_pkg.h"
#include "pdc_region_cache.h"
#include "pdc_region_cache_dl.h"
#include "pdc_region_prefetch.h"
#include "pdc_client_connect.h"

MPI_Comm client_cache_world_comm;
MPI_Comm client_cache_node_comm;

pdc_client_info client_info;

// New cache item, not in the list yet; its data is in slot_idx
pdc_object_data *
pdc_region_dl_new_item(pdcid_t obj_id, int ndim, uint64_t unit, uint64_t *offset, uint64_t *size,
                       uint64_t reg_buf_size, int slot_idx)
{
    pdc_object_data *item = (pdc_object_data *)PDC_calloc(1, sizeof(pdc_object_data));

    item->obj_id       = obj_id;
    item->unit         = unit;
    item->reg_ndim     = ndim;
    item->reg_buf_size = reg_buf_size;
    item->slot_idx     = slot_idx;
    item->target_rank  = -1;
    memcpy(item->reg_offset, offset, ndim * sizeof(uint64_t));
    memcpy(item->reg_size, size, ndim * sizeof(uint64_t));
    return item;
}

perr_t
pdc_region_dl_prepend(pdc_object_data *obj_cache_item)
{
    FUNC_ENTER(NULL);

    perr_t ret_value = SUCCEED;

    // Update metadata list information
    // Current item will be prepended to the head node
    obj_cache_item->next = client_info.local_cache_list_head;
    obj_cache_item->prev = NULL;

    // If list already exists, current head node's previous index should also point to itself
    if (client_info.local_cache_list_head != NULL) {
        obj_cache_item->next->prev = obj_cache_item;
    }

    // The new node becomes the head node of the object cache list
    client_info.local_cache_list_head = obj_cache_item;

    // If list did not exist, update the tail node metadata
    if (client_info.local_cache_list_tail == NULL) {
        client_info.local_cache_list_tail = obj_cache_item;
    }

    client_info.cached_item_num++;

done:
    fflush(stdout);
    FUNC_LEAVE(ret_value);
}

perr_t
pdc_region_dl_delete(pdc_object_data *obj_cache_item)
{
    FUNC_ENTER(NULL);

    perr_t ret_value = SUCCEED;

    // Update current item prev item's next item index
    if (obj_cache_item->prev != NULL) {
        obj_cache_item->prev->next = obj_cache_item->next;
    }
    else {
        client_info.local_cache_list_head = obj_cache_item->next; // if the item is the head
    }

    // Update current item next item's previous item index
    if (obj_cache_item->next != NULL) {
        obj_cache_item->next->prev = obj_cache_item->prev;
    }
    else {
        client_info.local_cache_list_tail = obj_cache_item->prev; // if the item is the tail
    }

    client_info.cached_item_num--;

done:
    fflush(stdout);
    FUNC_LEAVE(ret_value);
}

void
init_free_stack()
{
    for (int i = 0; i < MAX_SLOTS_PER_NODE - 1; i++) {
        client_info.header->free_list[i] = i + 1;
    }
    client_info.header->free_list[MAX_SLOTS_PER_NODE - 1] = SLOT_INVALID;
    client_info.header->free_stack_head                   = 0;
    client_info.header->pending_inserts                   = 0;
}

char *
get_data_ptr(int shared_index)
{
    size_t offset = (size_t)shared_index * (size_t)MAX_ITEM_SIZE;

    return client_info.node_shared_data_base + offset;
}

int
pop_free_slot()
{
    int current_head, next_head;

    MPI_Win_lock(MPI_LOCK_EXCLUSIVE, 0, 0, client_info.node_shared_data_win);

    MPI_Get(&current_head, 1, MPI_INT, 0, offsetof(SharedMemoryHeader, free_stack_head), 1, MPI_INT,
            client_info.node_shared_data_win);

    if (current_head == SLOT_INVALID) {
        MPI_Win_unlock(0, client_info.node_shared_data_win);
        return SLOT_INVALID;
    }

    MPI_Get(&next_head, 1, MPI_INT, 0, offsetof(SharedMemoryHeader, free_list) + current_head * sizeof(int),
            1, MPI_INT, client_info.node_shared_data_win);

    MPI_Put(&next_head, 1, MPI_INT, 0, offsetof(SharedMemoryHeader, free_stack_head), 1, MPI_INT,
            client_info.node_shared_data_win);

    MPI_Win_unlock(0, client_info.node_shared_data_win);

    // printf("[RANK %d] pop_free_slot index: %d\n", client_info.world_rank, current_head);
    // fflush(stdout);

    return current_head;
}

void
push_free_slot(int slot_idx)
{
    if (slot_idx < 0 || slot_idx >= MAX_SLOTS_PER_NODE)
        return;

    int current_head;

    MPI_Win_lock(MPI_LOCK_EXCLUSIVE, 0, 0, client_info.node_shared_data_win);

    MPI_Get(&current_head, 1, MPI_INT, 0, offsetof(SharedMemoryHeader, free_stack_head), 1, MPI_INT,
            client_info.node_shared_data_win);

    MPI_Put(&current_head, 1, MPI_INT, 0, offsetof(SharedMemoryHeader, free_list) + slot_idx * sizeof(int), 1,
            MPI_INT, client_info.node_shared_data_win);

    MPI_Put(&slot_idx, 1, MPI_INT, 0, offsetof(SharedMemoryHeader, free_stack_head), 1, MPI_INT,
            client_info.node_shared_data_win);

    MPI_Win_unlock(0, client_info.node_shared_data_win);
}

int
count_free_slots()
{
    int count = 0, current;

    MPI_Win_lock(MPI_LOCK_SHARED, 0, 0, client_info.node_shared_data_win);
    MPI_Fetch_and_op(NULL, &current, MPI_INT, 0, offsetof(SharedMemoryHeader, free_stack_head), MPI_NO_OP,
                     client_info.node_shared_data_win);
    MPI_Win_flush(0, client_info.node_shared_data_win);

    while (current != SLOT_INVALID && count < MAX_SLOTS_PER_NODE) {
        count++;
        current = client_info.header->free_list[current];
    }

    MPI_Win_unlock(0, client_info.node_shared_data_win);

    return count;
}

perr_t
pdc_region_dl_init()
{
    FUNC_ENTER(NULL);

    perr_t ret_value = SUCCEED;
    int *  node_world_ranks;
    int    mpi_alloc_error;

    client_info.world_rank = pdc_client_mpi_rank_g;
    client_info.world_size = pdc_client_mpi_size_g;

    MPI_Comm_dup(MPI_COMM_WORLD, &client_cache_world_comm); // Duplicate MPI_COMM_WORLD for client cache

    // 1. Node-level shared memory client_init setup
    MPI_Comm_split_type(client_cache_world_comm, MPI_COMM_TYPE_SHARED, client_info.world_rank, MPI_INFO_NULL,
                        &client_cache_node_comm);
    MPI_Comm_rank(client_cache_node_comm, &client_info.node_rank);
    MPI_Comm_size(client_cache_node_comm, &client_info.node_size);

    // 2. Broadcast node manager's rank to group clients by node
    if (client_info.node_rank == 0)
        client_info.node_manager_rank = client_info.world_rank;

    MPI_Bcast(&client_info.node_manager_rank, 1, MPI_INT, 0, client_cache_node_comm);

    // 3. Create a map for distinguishing inter-node and intra-node client group
    client_info.rank_to_node_id_map = (int *)PDC_malloc(client_info.world_size * sizeof(int));
    if (!client_info.rank_to_node_id_map)
        PGOTO_ERROR(FAIL, "pdc_region_dl_init - rank_to_node_id_map memory allocation failed");

    MPI_Allgather(&client_info.node_manager_rank, 1, MPI_INT, client_info.rank_to_node_id_map, 1, MPI_INT,
                  client_cache_world_comm);

    printf("[RANK %d] pdc_region_dl_init - step 3: rank_to_node_id_map\n", client_info.world_rank);
    fflush(stdout);

    // 4. Create local node map
    node_world_ranks = (int *)PDC_malloc(client_info.node_size * sizeof(int));
    if (!node_world_ranks)
        PGOTO_ERROR(FAIL, "pdc_region_dl_init - node_world_ranks memory allocation failed");

    client_info.world_to_node_rank_map = (int *)PDC_calloc(client_info.world_size, sizeof(int));
    if (!client_info.world_to_node_rank_map)
        PGOTO_ERROR(FAIL, "pdc_region_dl_init - world_to_node_rank_map memory allocation failed");

    MPI_Allgather(&client_info.world_rank, 1, MPI_INT, node_world_ranks, 1, MPI_INT, client_cache_node_comm);
    for (int i = 0; i < client_info.node_size; ++i) {
        client_info.world_to_node_rank_map[node_world_ranks[i]] = i;
    }

    PDC_free(node_world_ranks);

    printf("[RANK %d] pdc_region_dl_init - step 4: world_to_node_rank_map memory allocation\n",
           client_info.world_rank);
    fflush(stdout);

    // 5. Create shared memory window for node shared data
    MPI_Aint total_bytes =
        (MPI_Aint)sizeof(SharedMemoryHeader) + ((MPI_Aint)MAX_ITEM_SIZE * (MPI_Aint)MAX_SLOTS_PER_NODE);
    MPI_Aint local_alloc_size = (client_info.node_rank == 0) ? total_bytes : 0;

    mpi_alloc_error =
        MPI_Win_allocate_shared(local_alloc_size, 1, MPI_INFO_NULL, client_cache_node_comm,
                                &client_info.node_shared_base, &client_info.node_shared_data_win);

    printf("[RANK %d] pdc_region_dl_init - step 5: Create shared memory window\n", client_info.world_rank);
    fflush(stdout);

    if (mpi_alloc_error != MPI_SUCCESS)
        MPI_Abort(client_cache_world_comm, 1);

    if (client_info.node_rank != 0) {
        MPI_Aint size_out;
        int      disp_unit;

        MPI_Win_shared_query(client_info.node_shared_data_win, 0, &size_out, &disp_unit,
                             &client_info.node_shared_base);
    }

    printf("[RANK %d] pdc_region_dl_init - step 6: Query shared memory window\n", client_info.world_rank);
    fflush(stdout);

    client_info.header                = (SharedMemoryHeader *)client_info.node_shared_base;
    client_info.node_shared_data_base = (char *)client_info.node_shared_base + sizeof(SharedMemoryHeader);

    if (client_info.node_rank == 0) {
        init_free_stack();
        // printf("pdc_region_dl_init - step 7 -1 : Init the node shared memory\n", client_info.world_rank);
        // fflush(stdout);
        // memset(client_info.node_shared_data_base, 0, (size_t)MAX_SLOTS_PER_NODE * (size_t)MAX_ITEM_SIZE);
    }

    printf("[RANK %d] pdc_region_dl_init - step 7: Init the node shared memory\n", client_info.world_rank);
    fflush(stdout);

    MPI_Barrier(client_cache_node_comm);
    MPI_Barrier(client_cache_world_comm);

    client_info.local_cache_list_head = NULL;
    client_info.local_cache_list_tail = NULL;
    client_info.cached_item_num       = 0;
    client_info.client_cache_init     = 1;

done:
    FUNC_LEAVE(ret_value);
}

// TODO:
// Management of overalapping regions
perr_t
pdc_region_dl_insert(pdcid_t obj_id, int ndim, uint64_t unit, uint64_t *offset, uint64_t *size, void *buf,
                     uint64_t read_size)
{
    FUNC_ENTER(NULL);

    perr_t ret_value = SUCCEED;

    pdc_object_data *obj_cache_item = NULL;
    double           start          = MPI_Wtime();
    int              slot_idx;
    char *           data_ptr;

    // Check if the client cache has been initialized
    if (!client_info.client_cache_init)
        PGOTO_ERROR(FAIL, "pdc_region_dl_insert - object cache list not initialized");

    // Evict item if exceeding cache size
    // while (client_info.cached_item_num > MAX_ITEM_NUM)
    //     pdc_region_cache_evict();

    slot_idx = pop_free_slot();

    while (slot_idx == SLOT_INVALID) {
        // While an overlapped exchange runs, slots come back as sends complete;
        // evicting blindly could drop an item that is still being sent.
        if (pdc_region_dl_overlap_active())
            pdc_region_dl_overlap_make_room();
        else
            pdc_region_cache_evict();
        slot_idx = pop_free_slot();
    }
    if (pdc_region_dl_overlap_active())
        pdc_region_dl_overlap_record_insert();

    obj_cache_item = pdc_region_dl_new_item(
        obj_id, ndim, unit, offset, size, (read_size > MAX_ITEM_SIZE) ? MAX_ITEM_SIZE : read_size, slot_idx);
    if (!obj_cache_item)
        PGOTO_ERROR(FAIL, "PDC region cache - obj_cache_item memory allocation failed");
    obj_cache_item->consumed = 1; // inserted right after being read from the server

    data_ptr = get_data_ptr(slot_idx);
    memcpy(data_ptr, buf, obj_cache_item->reg_buf_size * sizeof(char));

    ret_value = pdc_region_dl_prepend(obj_cache_item);

    pdc_region_cache_timelog(start, "pdc_region_dl_insert - prepend new item to obj_cache list time");

done:
    FUNC_LEAVE(ret_value);
}

// TODO:
// Need to manage the object's offset cache information
// Need to manage when partial region is contained within
int
pdc_region_dl_local_search(pdcid_t obj_id, int ndim, uint64_t unit, uint64_t *offset, uint64_t *size,
                           void *buf, uint64_t read_size)
{
    FUNC_ENTER(NULL);

    perr_t ret_value = SUCCEED;

    uint64_t *       overlap_offset, *overlap_size;
    int              region_copy = 1;
    char *           data_ptr;
    int              is_cached = 0;
    pdc_object_data *obj_cache_iter;
    double           start, tmp_start, total_start = MPI_Wtime();

    if (!client_info.client_cache_init)
        PGOTO_ERROR(FAIL, "pdc_region_dl_search - object cache list not initialized");

    // Find if region is cached into local object list cache
    // If region contained, return the rank that contains the region
    obj_cache_iter = client_info.local_cache_list_head;

    while (obj_cache_iter != NULL) {
        if (obj_cache_iter->obj_id == obj_id) {
            is_cached = detect_region_contained(offset, size, obj_cache_iter->reg_offset,
                                                obj_cache_iter->reg_size, ndim);

            // If region contained, memcpy cached region data into transfer_request->buf
            if (is_cached) {
                start = MPI_Wtime();

                // Detect the offset range that is overlapped
                PDC_region_overlap_detect(ndim, offset, size, obj_cache_iter->reg_offset,
                                          obj_cache_iter->reg_size, &overlap_offset, &overlap_size);

                // memcpy the overlapped region
                data_ptr = get_data_ptr(obj_cache_iter->slot_idx);

                tmp_start = MPI_Wtime();

                for (int i = 0; i < ndim; ++i) {
                    if (offset[i] != obj_cache_iter->reg_offset[i]) {
                        region_copy = 0;
                        break;
                    }
                    if (size[i] != obj_cache_iter->reg_size[i]) {
                        region_copy = 0;
                        break;
                    }
                }

                if (region_copy) {
                    memcpy(buf, data_ptr, obj_cache_iter->reg_buf_size);
                    if (client_info.world_rank == 0)
                        printf("[RANK %d] Read entire region for obj_id: %lld\n", client_info.world_rank,
                               obj_id);
                }
                else {
                    memcpy_overlap_subregion(obj_cache_iter->reg_ndim, unit, data_ptr,
                                             obj_cache_iter->reg_offset, obj_cache_iter->reg_size, buf,
                                             offset, size, overlap_offset, overlap_size);
                }

                pdc_region_cache_timelog(tmp_start, "pdc_region_dl_local_search - memcpy data to buf");

                obj_cache_iter->consumed = 1;

                // Follow the LRU policy
                ret_value = pdc_region_dl_delete(obj_cache_iter);
                ret_value = pdc_region_dl_prepend(obj_cache_iter);

                free(overlap_offset);

                pdc_region_cache_timelog(start, "pdc_region_dl_local_search - local cache hit time");

                break;
            }
        }

        obj_cache_iter = obj_cache_iter->next;
    }

    if (!is_cached)
        printf("[RANK %d] cache miss for obj_id: %lld\n", client_info.world_rank, obj_id);
    else
        printf("[RANK %d] cache hit for obj_id: %lld\n", client_info.world_rank, obj_id);

    pdc_region_cache_timelog(total_start, "pdc_region_dl_local_search search time");

done:
    FUNC_LEAVE(is_cached);
}

perr_t
pdc_region_dl_update(pdcid_t obj_id, int ndim, uint64_t unit, uint64_t *offset, uint64_t *size, void *buf)
{
    FUNC_ENTER(NULL);

    perr_t ret_value = SUCCEED;

    pdc_object_data *obj_cache_iter, *obj_cache_item;
    int              is_overlapped = 0;

    if (client_info.cached_item_num == 0)
        goto done;

    // Find overlapping regions from the head of the list
    obj_cache_iter = client_info.local_cache_list_head;

    while (obj_cache_iter != NULL) {
        obj_cache_item = obj_cache_iter;
        obj_cache_iter = obj_cache_iter->next;

        if (obj_cache_item->obj_id == obj_id) {
            // Compare offset and offset + size and see if there is an overlap
            is_overlapped =
                check_overlap(ndim, offset, size, obj_cache_item->reg_offset, obj_cache_item->reg_size);

            // If there is overlapping region remove item from list
            if (is_overlapped) {
                // Delete the overlapped object item
                ret_value = pdc_region_dl_delete(obj_cache_item);

                push_free_slot(obj_cache_item->slot_idx);

                free(obj_cache_item);
            }
        }
    }

done:
    FUNC_LEAVE(ret_value);
}

perr_t
pdc_region_dl_evict()
{
    FUNC_ENTER(NULL);

    perr_t           ret_value = SUCCEED;
    pdc_object_data *obj_cache_item;

    if (client_info.cached_item_num == 0)
        goto done;

    // Evict the least recently used item that is not in the middle of being sent
    obj_cache_item = client_info.local_cache_list_tail;
    while (obj_cache_item != NULL && obj_cache_item->transfer_state == 1)
        obj_cache_item = obj_cache_item->prev;
    if (obj_cache_item == NULL)
        PGOTO_ERROR(FAIL, "pdc_region_dl_evict - every cached item is being sent");

    ret_value = pdc_region_dl_delete(obj_cache_item);

    push_free_slot(obj_cache_item->slot_idx);

    free(obj_cache_item);

done:
    FUNC_LEAVE(ret_value);
}

perr_t
pdc_region_dl_finalize()
{
    FUNC_ENTER(NULL);

    perr_t           ret_value = SUCCEED;
    pdc_object_data *obj_cache_iter, *obj_cache_item;
    double           start = MPI_Wtime();

    if (!client_info.client_cache_init)
        goto done;

    MPI_Barrier(client_cache_world_comm); // Make sure no shared memory is currently used

    free(client_info.world_to_node_rank_map);
    free(client_info.rank_to_node_id_map); // Free the rank_to_node_id_map

    obj_cache_iter = client_info.local_cache_list_tail;

    while (obj_cache_iter != NULL) {
        obj_cache_item = obj_cache_iter;
        obj_cache_iter = obj_cache_iter->prev;

        ret_value = pdc_region_dl_delete(obj_cache_item);
        free(obj_cache_item);
    }

    MPI_Win_free(&client_info.node_shared_data_win);
    MPI_Comm_free(&client_cache_node_comm);
    MPI_Comm_free(&client_cache_world_comm);

    pdc_region_cache_timelog(start, "pdc_region_dl_finalize - finalization time");

done:
    FUNC_LEAVE(ret_value);
}
