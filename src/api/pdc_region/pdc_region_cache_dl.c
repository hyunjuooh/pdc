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
        pdc_region_cache_evict();
        slot_idx = pop_free_slot();
    }

    obj_cache_item = (pdc_object_data *)PDC_malloc(sizeof(pdc_object_data));
    if (!obj_cache_item)
        PGOTO_ERROR(FAIL, "PDC region cache - obj_cache_item memory allocation failed");

    // Create object data
    obj_cache_item->obj_id             = obj_id;
    obj_cache_item->unit               = unit;
    obj_cache_item->reg_ndim           = ndim;
    obj_cache_item->reg_buf_size       = (read_size > MAX_ITEM_SIZE) ? MAX_ITEM_SIZE : read_size;
    obj_cache_item->target_rank        = -1;
    obj_cache_item->data_exchange_type = 0;
    obj_cache_item->slot_idx           = slot_idx;

    memcpy(obj_cache_item->reg_offset, offset, ndim * sizeof(uint64_t));
    memcpy(obj_cache_item->reg_size, size, ndim * sizeof(uint64_t));

    data_ptr = get_data_ptr(slot_idx);
    memcpy(data_ptr, buf, obj_cache_item->reg_buf_size * sizeof(char));

    obj_cache_item->prev = NULL;
    obj_cache_item->next = NULL;

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
pdc_region_dl_prepare_data_exchange(pdcid_t *global_prefetch_list, uint64_t *offset, uint64_t *size,
                                    int obj_prefetch_list_len)
{
    FUNC_ENTER(NULL);

    perr_t ret_value = SUCCEED;

    pdc_object_data *obj_cache_iter;
    int              i, j, prefetch_list_idx = 0, region_compare = 0, is_contained = 0;
    double           start = MPI_Wtime();

    if (offset != NULL) {
        region_compare = 1;
    }

    // Set the target rank for local cached object item
    obj_cache_iter = client_info.local_cache_list_head;

    while (obj_cache_iter != NULL) {
        for (i = 0; i < client_info.world_size; i++) {
            for (j = 0; j < obj_prefetch_list_len; j++) {
                if (obj_cache_iter->obj_id == global_prefetch_list[prefetch_list_idx]) {
                    if (!region_compare) {
                        obj_cache_iter->target_rank = i;
                        break;
                    }
                    else {
                        is_contained =
                            detect_region_contained(&offset[i], &size[i], obj_cache_iter->reg_offset,
                                                    obj_cache_iter->reg_size, obj_cache_iter->reg_ndim);

                        // If the prefetch list is contained within the local object cache list
                        // and if the target_rank is not the same as node_manager_rank
                        if (is_contained) {
                            obj_cache_iter->target_rank = i;
                            break;
                        }
                    }
                }

                prefetch_list_idx++;
            }
            if (obj_cache_iter->target_rank != -1) {
                break;
            }
        }

        is_contained      = 0;
        obj_cache_iter    = obj_cache_iter->next;
        prefetch_list_idx = 0;
    }

    // For debugging purpose
    // obj_cache_iter = client_info.local_cache_list_head;
    // while (obj_cache_iter != NULL) {
    //     printf("[RANK %d] prepare_data_exchange: object_id %lld , target_rank: %d\n",
    //     client_info.world_rank, obj_cache_iter->obj_id, obj_cache_iter->target_rank);

    //     fflush(stdout);

    //     obj_cache_iter = obj_cache_iter->next;
    // }

    pdc_region_cache_timelog(start, "pdc_region_dl_prepare_data_exchange - global prefetch list preparation");

done:
    fflush(stdout);
    FUNC_LEAVE(ret_value);
}

// Inter-node item header, exchanged before the data so every receiver knows
// what arrives from each sender (and in which order) before any data moves.
typedef struct {
    pdcid_t  obj_id;
    uint64_t unit;
    uint64_t reg_offset[3];
    uint64_t reg_size[3];
    uint64_t reg_buf_size;
    int      reg_ndim;
    int      padding;
} pdc_inter_item_hdr;

// Receive buffer for inter-node items that find no free slot. Allocated once,
// on first use, and only ever used by one blocking receive at a time.
static char *dl_drop_buf = NULL;

perr_t
pdc_region_dl_data_exchange(pdcid_t *global_prefetch_list, int obj_prefetch_list_len)
{
    FUNC_ENTER(NULL);

    perr_t ret_value = SUCCEED;

    pdc_object_data *exchange_head;
    double           start = MPI_Wtime(), tmp_timer;
    double           intra_time = 0.0, inter_meta_time = 0.0, inter_data_time = 0.0, delete_time = 0.0;
    int              node_size = client_info.node_size, world_size = client_info.world_size;
    int              old_cached_item_num = client_info.cached_item_num;
    int              chunk_size          = old_cached_item_num / NUM_CHUNKS;
    int              n_intra_sent = 0, n_intra_recv = 0, n_inter_sent = 0, n_inter_recv = 0, n_dropped = 0;

    // Items are exchanged in NUM_CHUNKS rounds so each node only needs free
    // slots for one round's incoming items before the round's sent items are freed.
    int *intra_send_counts = (int *)PDC_calloc(node_size, sizeof(int));
    int *intra_recv_counts = (int *)PDC_calloc(node_size, sizeof(int));
    int *intra_sdispls     = (int *)PDC_calloc(node_size, sizeof(int));
    int *intra_rdispls     = (int *)PDC_calloc(node_size, sizeof(int));
    int *inter_send_counts = (int *)PDC_calloc(world_size, sizeof(int));
    int *inter_recv_counts = (int *)PDC_calloc(world_size, sizeof(int));
    int *inter_sdispls     = (int *)PDC_calloc(world_size, sizeof(int));
    int *inter_rdispls     = (int *)PDC_calloc(world_size, sizeof(int));
    int *byte_counts_s     = (int *)PDC_calloc(world_size, sizeof(int));
    int *byte_counts_r     = (int *)PDC_calloc(world_size, sizeof(int));
    int *byte_displs_s     = (int *)PDC_calloc(world_size, sizeof(int));
    int *byte_displs_r     = (int *)PDC_calloc(world_size, sizeof(int));
    int *fill              = (int *)PDC_calloc(world_size > node_size ? world_size : node_size, sizeof(int));

    MPI_Datatype mpi_intra_transfer_unit;
    MPI_Type_contiguous(INTRA_TRANSFER_UNIT_SIZE, MPI_BYTE, &mpi_intra_transfer_unit);
    MPI_Type_commit(&mpi_intra_transfer_unit);

    exchange_head = client_info.local_cache_list_head;

    for (int c = 0; c < NUM_CHUNKS; c++) {
        int n_c = (c == NUM_CHUNKS - 1) ? old_cached_item_num - c * chunk_size : chunk_size;

        // Collect this round's items (the next n_c old items) and classify them.
        pdc_object_data **items     = NULL;
        int               n_intra   = 0;
        int               n_inter   = 0;
        if (n_c > 0)
            items = (pdc_object_data **)PDC_malloc(n_c * sizeof(pdc_object_data *));
        for (int i = 0; i < n_c && exchange_head != NULL; i++) {
            items[i]      = exchange_head;
            exchange_head = exchange_head->next;
        }

        memset(intra_send_counts, 0, node_size * sizeof(int));
        memset(inter_send_counts, 0, world_size * sizeof(int));
        for (int i = 0; i < n_c; i++) {
            int trg = items[i]->target_rank;
            if (trg < 0) {
                // Not in any rank's prefetch hint: dropped, its slot is freed below.
                items[i]->data_exchange_type = -1;
            }
            else if (client_info.rank_to_node_id_map[trg] == client_info.node_manager_rank) {
                items[i]->data_exchange_type = 0; // same node: hand over the slot
                intra_send_counts[client_info.world_to_node_rank_map[trg]]++;
                n_intra++;
            }
            else {
                items[i]->data_exchange_type = 1; // other node: send the data
                inter_send_counts[trg]++;
                n_inter++;
            }
        }

        // ---- Intra-node: send item metadata incl. slot index, no data copy ----
        tmp_timer = MPI_Wtime();
        MPI_Alltoall(intra_send_counts, 1, MPI_INT, intra_recv_counts, 1, MPI_INT, client_cache_node_comm);

        int intra_recv_total = 0, sdisp = 0;
        for (int i = 0; i < node_size; i++) {
            intra_sdispls[i] = sdisp;
            intra_rdispls[i] = intra_recv_total;
            sdisp += intra_send_counts[i];
            intra_recv_total += intra_recv_counts[i];
        }

        char *intra_send_buf = (char *)PDC_malloc((size_t)n_intra * INTRA_TRANSFER_UNIT_SIZE + 1);
        char *intra_recv_buf = (char *)PDC_malloc((size_t)intra_recv_total * INTRA_TRANSFER_UNIT_SIZE + 1);

        memset(fill, 0, node_size * sizeof(int));
        for (int i = 0; i < n_c; i++) {
            if (items[i]->data_exchange_type != 0)
                continue;
            int    trg_nrank = client_info.world_to_node_rank_map[items[i]->target_rank];
            char * p   = intra_send_buf + ((size_t)intra_sdispls[trg_nrank] + fill[trg_nrank]++) *
                                         INTRA_TRANSFER_UNIT_SIZE;
            size_t off = 0;
            memcpy(p + off, &items[i]->obj_id, sizeof(pdcid_t));
            off += sizeof(pdcid_t);
            memcpy(p + off, &items[i]->unit, sizeof(uint64_t));
            off += sizeof(uint64_t);
            memcpy(p + off, &items[i]->reg_ndim, sizeof(int));
            off += sizeof(int);
            memcpy(p + off, &items[i]->slot_idx, sizeof(int));
            off += sizeof(int);
            memcpy(p + off, items[i]->reg_offset, sizeof(uint64_t) * 3);
            off += sizeof(uint64_t) * 3;
            memcpy(p + off, items[i]->reg_size, sizeof(uint64_t) * 3);
            off += sizeof(uint64_t) * 3;
            memcpy(p + off, &items[i]->reg_buf_size, sizeof(uint64_t));
        }

        MPI_Alltoallv(intra_send_buf, intra_send_counts, intra_sdispls, mpi_intra_transfer_unit,
                      intra_recv_buf, intra_recv_counts, intra_rdispls, mpi_intra_transfer_unit,
                      client_cache_node_comm);

        for (int k = 0; k < intra_recv_total; k++) {
            pdc_object_data *item = (pdc_object_data *)PDC_malloc(sizeof(pdc_object_data));
            char *           p    = intra_recv_buf + (size_t)k * INTRA_TRANSFER_UNIT_SIZE;
            size_t           off  = 0;
            memcpy(&item->obj_id, p + off, sizeof(pdcid_t));
            off += sizeof(pdcid_t);
            memcpy(&item->unit, p + off, sizeof(uint64_t));
            off += sizeof(uint64_t);
            memcpy(&item->reg_ndim, p + off, sizeof(int));
            off += sizeof(int);
            memcpy(&item->slot_idx, p + off, sizeof(int));
            off += sizeof(int);
            memcpy(item->reg_offset, p + off, sizeof(uint64_t) * 3);
            off += sizeof(uint64_t) * 3;
            memcpy(item->reg_size, p + off, sizeof(uint64_t) * 3);
            off += sizeof(uint64_t) * 3;
            memcpy(&item->reg_buf_size, p + off, sizeof(uint64_t));
            item->target_rank        = -1;
            item->data_exchange_type = 0;
            pdc_region_dl_prepend(item);
        }
        free(intra_send_buf);
        free(intra_recv_buf);
        n_intra_sent += n_intra;
        n_intra_recv += intra_recv_total;
        intra_time += MPI_Wtime() - tmp_timer;

        // ---- Inter-node step 1: exchange item headers (small) ----
        tmp_timer = MPI_Wtime();
        MPI_Alltoall(inter_send_counts, 1, MPI_INT, inter_recv_counts, 1, MPI_INT, client_cache_world_comm);

        int inter_recv_total = 0;
        sdisp                = 0;
        for (int i = 0; i < world_size; i++) {
            inter_sdispls[i] = sdisp;
            inter_rdispls[i] = inter_recv_total;
            byte_counts_s[i] = inter_send_counts[i] * (int)sizeof(pdc_inter_item_hdr);
            byte_counts_r[i] = inter_recv_counts[i] * (int)sizeof(pdc_inter_item_hdr);
            byte_displs_s[i] = sdisp * (int)sizeof(pdc_inter_item_hdr);
            byte_displs_r[i] = inter_recv_total * (int)sizeof(pdc_inter_item_hdr);
            sdisp += inter_send_counts[i];
            inter_recv_total += inter_recv_counts[i];
        }

        pdc_inter_item_hdr *send_hdr  = (pdc_inter_item_hdr *)PDC_calloc(n_inter + 1, sizeof(pdc_inter_item_hdr));
        pdc_inter_item_hdr *recv_hdr  = (pdc_inter_item_hdr *)PDC_calloc(inter_recv_total + 1, sizeof(pdc_inter_item_hdr));
        pdc_object_data ** send_item = (pdc_object_data **)PDC_malloc((n_inter + 1) * sizeof(pdc_object_data *));

        // send_hdr/send_item are grouped by target rank, in a fixed order; the
        // item's position within its target's group is its MPI tag below.
        memset(fill, 0, world_size * sizeof(int));
        for (int i = 0; i < n_c; i++) {
            if (items[i]->data_exchange_type != 1)
                continue;
            int                 trg = items[i]->target_rank;
            int                 k   = inter_sdispls[trg] + fill[trg]++;
            pdc_inter_item_hdr *h   = &send_hdr[k];
            h->obj_id               = items[i]->obj_id;
            h->unit                 = items[i]->unit;
            h->reg_ndim             = items[i]->reg_ndim;
            h->reg_buf_size         = items[i]->reg_buf_size;
            memcpy(h->reg_offset, items[i]->reg_offset, sizeof(uint64_t) * 3);
            memcpy(h->reg_size, items[i]->reg_size, sizeof(uint64_t) * 3);
            send_item[k] = items[i];
        }

        MPI_Alltoallv(send_hdr, byte_counts_s, byte_displs_s, MPI_BYTE, recv_hdr, byte_counts_r, byte_displs_r,
                      MPI_BYTE, client_cache_world_comm);
        inter_meta_time += MPI_Wtime() - tmp_timer;

        // ---- Inter-node step 2: data straight from/into shared-window slots ----
        tmp_timer          = MPI_Wtime();
        MPI_Request *rreq  = (MPI_Request *)PDC_malloc((inter_recv_total + 1) * sizeof(MPI_Request));
        MPI_Request *sreq  = (MPI_Request *)PDC_malloc((n_inter + 1) * sizeof(MPI_Request));
        int *        rslot = (int *)PDC_malloc((inter_recv_total + 1) * sizeof(int));

        // Reserve a slot for every incoming item before any data moves, and post
        // the receive straight into it. Tag = position within the sender's group,
        // so a receive can only ever match the item its header describes.
        for (int s = 0; s < world_size; s++) {
            for (int j = 0; j < inter_recv_counts[s]; j++) {
                int k    = inter_rdispls[s] + j;
                rslot[k] = pop_free_slot();
                if (rslot[k] == SLOT_INVALID) {
                    rreq[k] = MPI_REQUEST_NULL; // received into dl_drop_buf after the others
                    continue;
                }
                MPI_Irecv(get_data_ptr(rslot[k]), (int)recv_hdr[k].reg_buf_size, MPI_BYTE, s, j,
                          client_cache_world_comm, &rreq[k]);
            }
        }

        for (int t = 0; t < world_size; t++) {
            for (int j = 0; j < inter_send_counts[t]; j++) {
                int k = inter_sdispls[t] + j;
                MPI_Isend(get_data_ptr(send_item[k]->slot_idx), (int)send_hdr[k].reg_buf_size, MPI_BYTE, t, j,
                          client_cache_world_comm, &sreq[k]);
            }
        }

        MPI_Waitall(inter_recv_total, rreq, MPI_STATUSES_IGNORE);

        // Items with no free slot on this node: drain them one at a time and drop
        // them, as before (the senders' Isends must still be matched).
        for (int s = 0; s < world_size; s++) {
            for (int j = 0; j < inter_recv_counts[s]; j++) {
                int k = inter_rdispls[s] + j;
                if (rslot[k] != SLOT_INVALID)
                    continue;
                if (dl_drop_buf == NULL)
                    dl_drop_buf = (char *)PDC_malloc(MAX_ITEM_SIZE);
                MPI_Recv(dl_drop_buf, (int)recv_hdr[k].reg_buf_size, MPI_BYTE, s, j, client_cache_world_comm,
                         MPI_STATUS_IGNORE);
                n_dropped++;
            }
        }

        // Only now do received items become visible in the cache list.
        for (int k = 0; k < inter_recv_total; k++) {
            if (rslot[k] == SLOT_INVALID)
                continue;
            pdc_object_data *item = (pdc_object_data *)PDC_malloc(sizeof(pdc_object_data));
            item->obj_id          = recv_hdr[k].obj_id;
            item->unit            = recv_hdr[k].unit;
            item->reg_ndim        = recv_hdr[k].reg_ndim;
            item->reg_buf_size    = recv_hdr[k].reg_buf_size;
            memcpy(item->reg_offset, recv_hdr[k].reg_offset, sizeof(uint64_t) * 3);
            memcpy(item->reg_size, recv_hdr[k].reg_size, sizeof(uint64_t) * 3);
            item->slot_idx           = rslot[k];
            item->target_rank        = -1;
            item->data_exchange_type = 0;
            pdc_region_dl_prepend(item);
        }

        // Sent slots may only be freed once their sends have completed.
        MPI_Waitall(n_inter, sreq, MPI_STATUSES_IGNORE);

        n_inter_sent += n_inter;
        n_inter_recv += inter_recv_total;
        free(rreq);
        free(sreq);
        free(rslot);
        free(send_hdr);
        free(recv_hdr);
        free(send_item);
        inter_data_time += MPI_Wtime() - tmp_timer;

        // ---- Remove this round's old items ----
        // Same-node items keep their slot (the receiver now owns it); slots of
        // items sent to another node, or dropped, return to the node's free list.
        tmp_timer = MPI_Wtime();
        for (int i = 0; i < n_c; i++) {
            if (items[i]->data_exchange_type != 0)
                push_free_slot(items[i]->slot_idx);
            pdc_region_dl_delete(items[i]);
            free(items[i]);
        }
        free(items);
        delete_time += MPI_Wtime() - tmp_timer;
    }

    MPI_Type_free(&mpi_intra_transfer_unit);

    free(intra_send_counts);
    free(intra_recv_counts);
    free(intra_sdispls);
    free(intra_rdispls);
    free(inter_send_counts);
    free(inter_recv_counts);
    free(inter_sdispls);
    free(inter_rdispls);
    free(byte_counts_s);
    free(byte_counts_r);
    free(byte_displs_s);
    free(byte_displs_r);
    free(fill);

    printf("[RANK %d] data exchange: rounds %d, intra-node sent %d recv %d, inter-node sent %d recv %d, "
           "dropped (no free slot) %d | intra %.3fs, inter headers %.3fs, inter data %.3fs, delete %.3fs\n",
           client_info.world_rank, NUM_CHUNKS, n_intra_sent, n_intra_recv, n_inter_sent, n_inter_recv,
           n_dropped, intra_time, inter_meta_time, inter_data_time, delete_time);

    pdc_region_cache_timelog(start, "pdc_region_dl_data_exchange - data exchange total time");

done:
    fflush(stdout);
    FUNC_LEAVE(ret_value);
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

    // Evict the last item following LRU policy
    obj_cache_item = client_info.local_cache_list_tail;

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
