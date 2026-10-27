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

// Overlapped data exchange (defined at the end of this file)
static int  overlap_exchange_active(void);
static void overlap_exchange_make_room(void);
static void overlap_exchange_record_insert(void);

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
        if (overlap_exchange_active())
            overlap_exchange_make_room();
        else
            pdc_region_cache_evict();
        slot_idx = pop_free_slot();
    }
    if (overlap_exchange_active())
        overlap_exchange_record_insert();

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
    obj_cache_item->consumed           = 1; // inserted right after being read from the server
    obj_cache_item->transfer_state     = 0;

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
            item->consumed           = 0;
            item->transfer_state     = 0;
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
            item->consumed           = 0;
            item->transfer_state     = 0;
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

/**************************************************************************/
/* Data exchange overlapped with computation                              */
/*                                                                        */
/* An item is forwarded to its next-epoch owner as soon as its current    */
/* owner has read it this epoch (server read in epoch 1, cache hit after  */
/* that), so the exchange runs while the application computes.           */
/**************************************************************************/

// Sent ahead of every item. A same-node item hands over its slot; for an item
// on another node the data follows as a second message.
typedef struct {
    pdcid_t  obj_id;
    uint64_t unit;
    uint64_t reg_offset[3];
    uint64_t reg_size[3];
    uint64_t reg_buf_size;
    int      reg_ndim;
    int      slot_idx;
    int      valid; // 0: the sender does not have the item, nothing follows
    int      padding;
} overlap_item_header;

// Where an object is read this epoch and who reads it next epoch
typedef struct {
    pdcid_t obj_id;
    int     current_owner;
    int     current_position; // position in the current owner's list, used for message tags
    int     next_owner;
} overlap_object_record;

typedef struct {
    pdcid_t             obj_id;
    int                 target_rank;
    int                 tag;
    int                 same_node;
    int                 state; // 0 waiting to be read, 1 being sent, 2 done
    pdc_object_data *   item;
    overlap_item_header header;
    MPI_Request         requests[2]; // header, data
} overlap_outgoing_item;

typedef struct {
    pdcid_t             obj_id;
    int                 source_rank;
    int                 tag;
    int                 same_node;
    int                 state; // 0 header posted, 1 waiting for a free slot, 2 data posted, 3 done
    int                 slot_idx;
    overlap_item_header header;
    MPI_Request         header_request;
    MPI_Request         data_request;
} overlap_incoming_item;

static struct {
    int                    active;
    int                    comm_created;
    MPI_Comm               comm;
    overlap_object_record *records;
    int                    num_records;
    pdcid_t *              next_ids_sorted; // this rank's next-epoch objects
    int                    num_next_ids;
    overlap_outgoing_item *outgoing;
    int                    num_outgoing;
    int                    num_outgoing_done;
    overlap_incoming_item *incoming;
    int                    num_incoming;
    int                    num_incoming_done;
    int                    pending_inserts; // this rank's share of the node's pending_inserts
    int                    num_same_node_sent;
    int                    num_other_node_sent;
    int                    num_same_node_received;
    int                    num_other_node_received;
    int                    num_missing;
} overlap_exchange;

static int
overlap_exchange_active(void)
{
    return overlap_exchange.active;
}

static int
compare_object_records(const void *a, const void *b)
{
    pdcid_t x = ((const overlap_object_record *)a)->obj_id;
    pdcid_t y = ((const overlap_object_record *)b)->obj_id;
    return (x > y) - (x < y);
}

static int
compare_object_ids(const void *a, const void *b)
{
    pdcid_t x = *(const pdcid_t *)a;
    pdcid_t y = *(const pdcid_t *)b;
    return (x > y) - (x < y);
}

static overlap_object_record *
find_object_record(pdcid_t obj_id)
{
    overlap_object_record key;
    key.obj_id = obj_id;
    return (overlap_object_record *)bsearch(&key, overlap_exchange.records, overlap_exchange.num_records,
                                            sizeof(overlap_object_record), compare_object_records);
}

static int
read_next_epoch(pdcid_t obj_id)
{
    return bsearch(&obj_id, overlap_exchange.next_ids_sorted, overlap_exchange.num_next_ids, sizeof(pdcid_t),
                   compare_object_ids) != NULL;
}

// This epoch's copy of an object (not one received for the next epoch)
static pdc_object_data *
find_current_item(pdcid_t obj_id)
{
    pdc_object_data *obj_cache_iter;

    for (obj_cache_iter = client_info.local_cache_list_head; obj_cache_iter != NULL;
         obj_cache_iter = obj_cache_iter->next) {
        if (obj_cache_iter->obj_id == obj_id && obj_cache_iter->transfer_state != 2)
            return obj_cache_iter;
    }
    return NULL;
}

static void
node_pending_inserts_add(int delta)
{
    int previous;

    MPI_Win_lock(MPI_LOCK_EXCLUSIVE, 0, 0, client_info.node_shared_data_win);
    MPI_Fetch_and_op(&delta, &previous, MPI_INT, 0, offsetof(SharedMemoryHeader, pending_inserts), MPI_SUM,
                     client_info.node_shared_data_win);
    MPI_Win_unlock(0, client_info.node_shared_data_win);
}

static int
node_pending_inserts_get(void)
{
    int value;

    MPI_Win_lock(MPI_LOCK_SHARED, 0, 0, client_info.node_shared_data_win);
    MPI_Fetch_and_op(NULL, &value, MPI_INT, 0, offsetof(SharedMemoryHeader, pending_inserts), MPI_NO_OP,
                     client_info.node_shared_data_win);
    MPI_Win_unlock(0, client_info.node_shared_data_win);
    return value;
}

// Called by pdc_region_dl_insert: one of this epoch's server reads has its slot now
static void
overlap_exchange_record_insert(void)
{
    if (overlap_exchange.pending_inserts > 0) {
        overlap_exchange.pending_inserts--;
        node_pending_inserts_add(-1);
    }
}

static void
overlap_start_send(overlap_outgoing_item *outgoing_item, pdc_object_data *item)
{
    overlap_item_header *header = &outgoing_item->header;

    memset(header, 0, sizeof(overlap_item_header));
    header->obj_id = outgoing_item->obj_id;
    if (item != NULL) {
        header->unit         = item->unit;
        header->reg_ndim     = item->reg_ndim;
        header->reg_buf_size = item->reg_buf_size;
        header->slot_idx     = item->slot_idx;
        header->valid        = 1;
        memcpy(header->reg_offset, item->reg_offset, sizeof(uint64_t) * 3);
        memcpy(header->reg_size, item->reg_size, sizeof(uint64_t) * 3);
        item->transfer_state = 1;
    }

    MPI_Isend(header, sizeof(overlap_item_header), MPI_BYTE, outgoing_item->target_rank,
              2 * outgoing_item->tag, overlap_exchange.comm, &outgoing_item->requests[0]);
    outgoing_item->requests[1] = MPI_REQUEST_NULL;
    if (item != NULL && !outgoing_item->same_node)
        MPI_Isend(get_data_ptr(item->slot_idx), (int)item->reg_buf_size, MPI_BYTE, outgoing_item->target_rank,
                  2 * outgoing_item->tag + 1, overlap_exchange.comm, &outgoing_item->requests[1]);

    outgoing_item->item  = item;
    outgoing_item->state = 1;
}

// Forward every item read since the last call
static void
overlap_start_sends_for_read_items(void)
{
    for (int i = 0; i < overlap_exchange.num_outgoing; i++) {
        overlap_outgoing_item *outgoing_item = &overlap_exchange.outgoing[i];
        if (outgoing_item->state != 0)
            continue;
        pdc_object_data *item = find_current_item(outgoing_item->obj_id);
        if (item != NULL && item->consumed && item->transfer_state == 0)
            overlap_start_send(outgoing_item, item);
    }
}

// Returns how many sends completed
static int
overlap_test_sends(void)
{
    int completed = 0, flag;

    for (int i = 0; i < overlap_exchange.num_outgoing; i++) {
        overlap_outgoing_item *outgoing_item = &overlap_exchange.outgoing[i];
        if (outgoing_item->state != 1)
            continue;
        MPI_Testall(2, outgoing_item->requests, &flag, MPI_STATUSES_IGNORE);
        if (!flag)
            continue;

        if (outgoing_item->item != NULL) {
            // A same-node receiver now owns the slot; data sent to another node
            // has its own copy there, so the slot returns to this node's free list.
            if (outgoing_item->same_node)
                overlap_exchange.num_same_node_sent++;
            else {
                push_free_slot(outgoing_item->item->slot_idx);
                overlap_exchange.num_other_node_sent++;
            }
            pdc_region_dl_delete(outgoing_item->item);
            free(outgoing_item->item);
            outgoing_item->item = NULL;
        }
        outgoing_item->state = 2;
        overlap_exchange.num_outgoing_done++;
        completed++;
    }
    return completed;
}

static void
overlap_add_received_item(overlap_item_header *header, int slot_idx)
{
    pdc_object_data *item = (pdc_object_data *)PDC_malloc(sizeof(pdc_object_data));

    item->obj_id       = header->obj_id;
    item->unit         = header->unit;
    item->reg_ndim     = header->reg_ndim;
    item->reg_buf_size = header->reg_buf_size;
    memcpy(item->reg_offset, header->reg_offset, sizeof(uint64_t) * 3);
    memcpy(item->reg_size, header->reg_size, sizeof(uint64_t) * 3);
    item->slot_idx           = slot_idx;
    item->target_rank        = -1;
    item->data_exchange_type = 0;
    item->consumed           = 0;
    item->transfer_state     = 2; // for the next epoch: never sent or read this epoch
    pdc_region_dl_prepend(item);
}

// Returns how many receives changed state
static int
overlap_test_receives(void)
{
    int changed = 0, flag, available = -1;

    for (int i = 0; i < overlap_exchange.num_incoming; i++) {
        overlap_incoming_item *incoming_item = &overlap_exchange.incoming[i];

        if (incoming_item->state == 0) {
            MPI_Test(&incoming_item->header_request, &flag, MPI_STATUS_IGNORE);
            if (!flag)
                continue;
            changed++;
            if (!incoming_item->header.valid) {
                incoming_item->state = 3; // read from the servers next epoch instead
                overlap_exchange.num_incoming_done++;
                overlap_exchange.num_missing++;
                continue;
            }
            if (incoming_item->same_node) {
                overlap_add_received_item(&incoming_item->header, incoming_item->header.slot_idx);
                incoming_item->state = 3;
                overlap_exchange.num_incoming_done++;
                overlap_exchange.num_same_node_received++;
                continue;
            }
            incoming_item->state = 1;
        }

        if (incoming_item->state == 1) {
            // Leave enough free slots for this epoch's remaining server reads on the node
            if (available < 0)
                available = count_free_slots() - node_pending_inserts_get();
            if (available <= 0)
                continue;
            incoming_item->slot_idx = pop_free_slot();
            if (incoming_item->slot_idx == SLOT_INVALID) {
                available = 0;
                continue;
            }
            available--;
            MPI_Irecv(get_data_ptr(incoming_item->slot_idx), (int)incoming_item->header.reg_buf_size, MPI_BYTE,
                      incoming_item->source_rank, 2 * incoming_item->tag + 1, overlap_exchange.comm,
                      &incoming_item->data_request);
            incoming_item->state = 2;
            changed++;
        }

        if (incoming_item->state == 2) {
            MPI_Test(&incoming_item->data_request, &flag, MPI_STATUS_IGNORE);
            if (!flag)
                continue;
            overlap_add_received_item(&incoming_item->header, incoming_item->slot_idx);
            incoming_item->state = 3;
            overlap_exchange.num_incoming_done++;
            overlap_exchange.num_other_node_received++;
            changed++;
        }
    }
    return changed;
}

// Called by pdc_region_dl_insert when no slot is free: completed sends free
// slots; otherwise evict (eviction never touches an item that is being sent).
static void
overlap_exchange_make_room(void)
{
    if (overlap_test_sends() > 0)
        return;
    pdc_region_cache_evict();
}

// Drop cached items this rank does not read next epoch and that are not being sent
static void
overlap_free_unneeded_items(void)
{
    pdc_object_data *obj_cache_iter = client_info.local_cache_list_head, *obj_cache_item;

    while (obj_cache_iter != NULL) {
        obj_cache_item = obj_cache_iter;
        obj_cache_iter = obj_cache_iter->next;
        if (obj_cache_item->transfer_state != 0 || read_next_epoch(obj_cache_item->obj_id))
            continue;
        pdc_region_dl_delete(obj_cache_item);
        push_free_slot(obj_cache_item->slot_idx);
        free(obj_cache_item);
    }
}

perr_t
pdc_region_dl_overlap_begin(pdcid_t *cur_ids, int cur_len, pdcid_t *next_ids, int next_len)
{
    FUNC_ENTER(NULL);

    perr_t           ret_value = SUCCEED;
    int              world_size = client_info.world_size, my_rank = client_info.world_rank;
    int              total_cur = 0, total_next = 0;
    pdc_object_data *obj_cache_iter;
    double           start = MPI_Wtime();

    if (!client_info.client_cache_init)
        PGOTO_ERROR(FAIL, "pdc_region_dl_overlap_begin - client cache not initialized");
    if (overlap_exchange.active)
        PGOTO_ERROR(FAIL, "pdc_region_dl_overlap_begin - previous exchange has not ended");

    if (!overlap_exchange.comm_created) {
        MPI_Comm_dup(client_cache_world_comm, &overlap_exchange.comm);
        overlap_exchange.comm_created = 1;
    }

    // Every rank's lists, so each rank knows the current and next owner of every object
    int *cur_counts  = (int *)PDC_calloc(world_size, sizeof(int));
    int *next_counts = (int *)PDC_calloc(world_size, sizeof(int));
    int *cur_displs  = (int *)PDC_calloc(world_size, sizeof(int));
    int *next_displs = (int *)PDC_calloc(world_size, sizeof(int));

    MPI_Allgather(&cur_len, 1, MPI_INT, cur_counts, 1, MPI_INT, overlap_exchange.comm);
    MPI_Allgather(&next_len, 1, MPI_INT, next_counts, 1, MPI_INT, overlap_exchange.comm);
    for (int r = 0; r < world_size; r++) {
        cur_displs[r]  = total_cur;
        next_displs[r] = total_next;
        total_cur += cur_counts[r];
        total_next += next_counts[r];
    }

    pdcid_t *all_cur  = (pdcid_t *)PDC_malloc((total_cur + 1) * sizeof(pdcid_t));
    pdcid_t *all_next = (pdcid_t *)PDC_malloc((total_next + 1) * sizeof(pdcid_t));
    MPI_Allgatherv(cur_ids, cur_len, MPI_UINT64_T, all_cur, cur_counts, cur_displs, MPI_UINT64_T,
                   overlap_exchange.comm);
    MPI_Allgatherv(next_ids, next_len, MPI_UINT64_T, all_next, next_counts, next_displs, MPI_UINT64_T,
                   overlap_exchange.comm);

    overlap_exchange.records = (overlap_object_record *)PDC_malloc((total_cur + 1) * sizeof(overlap_object_record));
    overlap_exchange.num_records = total_cur;
    for (int r = 0; r < world_size; r++) {
        for (int i = 0; i < cur_counts[r]; i++) {
            overlap_object_record *record = &overlap_exchange.records[cur_displs[r] + i];
            record->obj_id                = all_cur[cur_displs[r] + i];
            record->current_owner         = r;
            record->current_position      = i;
            record->next_owner            = -1;
        }
    }
    qsort(overlap_exchange.records, total_cur, sizeof(overlap_object_record), compare_object_records);
    for (int r = 0; r < world_size; r++) {
        for (int i = 0; i < next_counts[r]; i++) {
            overlap_object_record *record = find_object_record(all_next[next_displs[r] + i]);
            if (record != NULL)
                record->next_owner = r;
        }
    }

    overlap_exchange.num_next_ids    = next_len;
    overlap_exchange.next_ids_sorted = (pdcid_t *)PDC_malloc((next_len + 1) * sizeof(pdcid_t));
    memcpy(overlap_exchange.next_ids_sorted, next_ids, next_len * sizeof(pdcid_t));
    qsort(overlap_exchange.next_ids_sorted, next_len, sizeof(pdcid_t), compare_object_ids);

    // Items delivered by the previous exchange are this epoch's data now
    for (obj_cache_iter = client_info.local_cache_list_head; obj_cache_iter != NULL;
         obj_cache_iter = obj_cache_iter->next) {
        obj_cache_iter->consumed       = 0;
        obj_cache_iter->transfer_state = 0;
    }

    // Outgoing: this rank's objects that another rank reads next epoch
    overlap_exchange.outgoing =
        (overlap_outgoing_item *)PDC_calloc(cur_len + 1, sizeof(overlap_outgoing_item));
    overlap_exchange.num_outgoing = 0;
    for (int i = 0; i < cur_len; i++) {
        overlap_object_record *record = find_object_record(cur_ids[i]);
        if (record == NULL || record->current_owner != my_rank || record->current_position != i ||
            record->next_owner < 0 || record->next_owner == my_rank)
            continue;
        overlap_outgoing_item *outgoing_item = &overlap_exchange.outgoing[overlap_exchange.num_outgoing++];
        outgoing_item->obj_id                = record->obj_id;
        outgoing_item->target_rank           = record->next_owner;
        outgoing_item->tag                   = i;
        outgoing_item->same_node =
            client_info.rank_to_node_id_map[record->next_owner] == client_info.node_manager_rank;
        outgoing_item->requests[0] = MPI_REQUEST_NULL;
        outgoing_item->requests[1] = MPI_REQUEST_NULL;
    }

    // Incoming: objects this rank reads next epoch that another rank holds now.
    // Headers are small, so all their receives are posted right away.
    overlap_exchange.incoming =
        (overlap_incoming_item *)PDC_calloc(next_len + 1, sizeof(overlap_incoming_item));
    overlap_exchange.num_incoming = 0;
    for (int i = 0; i < next_len; i++) {
        overlap_object_record *record = find_object_record(next_ids[i]);
        if (record == NULL || record->next_owner != my_rank || record->current_owner == my_rank)
            continue;
        overlap_incoming_item *incoming_item = &overlap_exchange.incoming[overlap_exchange.num_incoming++];
        incoming_item->obj_id                = record->obj_id;
        incoming_item->source_rank           = record->current_owner;
        incoming_item->tag                   = record->current_position;
        incoming_item->same_node =
            client_info.rank_to_node_id_map[record->current_owner] == client_info.node_manager_rank;
        incoming_item->slot_idx     = SLOT_INVALID;
        incoming_item->data_request = MPI_REQUEST_NULL;
        MPI_Irecv(&incoming_item->header, sizeof(overlap_item_header), MPI_BYTE, incoming_item->source_rank,
                  2 * incoming_item->tag, overlap_exchange.comm, &incoming_item->header_request);
    }

    // This epoch's server reads still to come (epoch 1: all of them); receives
    // for next-epoch data must leave this many slots free on the node.
    overlap_exchange.pending_inserts = 0;
    for (int i = 0; i < cur_len; i++) {
        if (find_current_item(cur_ids[i]) == NULL)
            overlap_exchange.pending_inserts++;
    }
    if (overlap_exchange.pending_inserts > 0)
        node_pending_inserts_add(overlap_exchange.pending_inserts);

    overlap_exchange.num_outgoing_done       = 0;
    overlap_exchange.num_incoming_done       = 0;
    overlap_exchange.num_same_node_sent      = 0;
    overlap_exchange.num_other_node_sent     = 0;
    overlap_exchange.num_same_node_received  = 0;
    overlap_exchange.num_other_node_received = 0;
    overlap_exchange.num_missing             = 0;
    overlap_exchange.active                  = 1;

    free(cur_counts);
    free(next_counts);
    free(cur_displs);
    free(next_displs);
    free(all_cur);
    free(all_next);

    pdc_region_cache_timelog(start, "pdc_region_dl_overlap_begin - total time");

done:
    fflush(stdout);
    FUNC_LEAVE(ret_value);
}

perr_t
pdc_region_dl_overlap_progress(void)
{
    FUNC_ENTER(NULL);

    perr_t ret_value = SUCCEED;

    if (!overlap_exchange.active)
        goto done;

    overlap_start_sends_for_read_items();
    overlap_test_sends();
    overlap_test_receives();

done:
    FUNC_LEAVE(ret_value);
}

perr_t
pdc_region_dl_overlap_end(void)
{
    FUNC_ENTER(NULL);

    perr_t ret_value = SUCCEED;
    double start     = MPI_Wtime(), last_change = MPI_Wtime();
    int    warned    = 0, open_outgoing, open_incoming;

    if (!overlap_exchange.active)
        goto done;

    // No more server reads this epoch: release this rank's reserved slots
    if (overlap_exchange.pending_inserts > 0) {
        node_pending_inserts_add(-overlap_exchange.pending_inserts);
        overlap_exchange.pending_inserts = 0;
    }

    // Items not read this epoch still go to their next owner; items this rank
    // never had are reported, so the receiver reads them from the servers.
    for (int i = 0; i < overlap_exchange.num_outgoing; i++) {
        overlap_outgoing_item *outgoing_item = &overlap_exchange.outgoing[i];
        if (outgoing_item->state == 0) {
            pdc_object_data *item = find_current_item(outgoing_item->obj_id);
            overlap_start_send(outgoing_item, (item != NULL && item->transfer_state == 0) ? item : NULL);
        }
    }

    open_outgoing = overlap_exchange.num_outgoing - overlap_exchange.num_outgoing_done;
    open_incoming = overlap_exchange.num_incoming - overlap_exchange.num_incoming_done;

    while (overlap_exchange.num_outgoing_done < overlap_exchange.num_outgoing ||
           overlap_exchange.num_incoming_done < overlap_exchange.num_incoming) {
        if (overlap_test_sends() + overlap_test_receives() > 0) {
            last_change = MPI_Wtime();
            continue;
        }
        // Receives may be waiting for a slot: free what nobody needs any more
        overlap_free_unneeded_items();
        if (!warned && MPI_Wtime() - last_change > 60.0) {
            printf("[RANK %d] overlapped exchange: no progress for 60 s (%d sends, %d receives open); "
                   "MAX_SLOTS_PER_NODE may be too small for the next epoch's items\n",
                   client_info.world_rank, overlap_exchange.num_outgoing - overlap_exchange.num_outgoing_done,
                   overlap_exchange.num_incoming - overlap_exchange.num_incoming_done);
            fflush(stdout);
            warned = 1;
        }
        usleep(50);
    }

    // What is left in the cache is this rank's data for the next epoch
    overlap_free_unneeded_items();
    for (pdc_object_data *obj_cache_iter = client_info.local_cache_list_head; obj_cache_iter != NULL;
         obj_cache_iter = obj_cache_iter->next) {
        obj_cache_iter->consumed       = 0;
        obj_cache_iter->transfer_state = 0;
    }

    printf("[RANK %d] overlapped exchange: sent same-node %d other-node %d, received same-node %d "
           "other-node %d, missing %d | still open when the epoch ended: %d sends, %d receives, "
           "end wait %.3fs\n",
           client_info.world_rank, overlap_exchange.num_same_node_sent, overlap_exchange.num_other_node_sent,
           overlap_exchange.num_same_node_received, overlap_exchange.num_other_node_received,
           overlap_exchange.num_missing, open_outgoing, open_incoming, MPI_Wtime() - start);

    free(overlap_exchange.records);
    free(overlap_exchange.next_ids_sorted);
    free(overlap_exchange.outgoing);
    free(overlap_exchange.incoming);
    overlap_exchange.records         = NULL;
    overlap_exchange.next_ids_sorted = NULL;
    overlap_exchange.outgoing        = NULL;
    overlap_exchange.incoming        = NULL;
    overlap_exchange.active          = 0;

    pdc_region_cache_timelog(start, "pdc_region_dl_overlap_end - wait time");

done:
    fflush(stdout);
    FUNC_LEAVE(ret_value);
}
