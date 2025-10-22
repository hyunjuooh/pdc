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

#define NUM_CHUNKS         8
#define TRANSFER_UNIT_SIZE (sizeof(pdcid_t) + sizeof(int) + sizeof(uint64_t) * 8 + MAX_ITEM_SIZE)

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

perr_t
pdc_region_dl_init()
{
    perr_t ret_value = SUCCEED;
    int    mpi_alloc_error;
    int *  node_world_ranks;

    FUNC_ENTER(NULL);

    client_info.world_rank = pdc_client_mpi_rank_g;
    client_info.world_size = pdc_client_mpi_size_g;

    // Duplicate MPI_COMM_WORLD for client cache
    MPI_Comm_dup(MPI_COMM_WORLD, &client_cache_world_comm);

    // Node-level shared memory client_init setup
    MPI_Comm_split_type(client_cache_world_comm, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL,
                        &client_cache_node_comm);
    MPI_Comm_rank(client_cache_node_comm, &client_info.node_rank);
    MPI_Comm_size(client_cache_node_comm, &client_info.node_size);

    // Broadcast node manager's rank to group clients by node
    if (client_info.node_rank == 0)
        client_info.node_manager_rank = client_info.world_rank;

    MPI_Bcast(&client_info.node_manager_rank, 1, MPI_INT, 0, client_cache_node_comm);

    // Create a map for distinguishing inter-node and intra-node client group
    client_info.rank_to_node_id_map = (int *)PDC_malloc(client_info.world_size * sizeof(int));
    if (!client_info.rank_to_node_id_map)
        PGOTO_ERROR(FAIL, "pdc_region_dl_init - rank_to_node_id_map memory allocation failed");

    MPI_Allgather(&client_info.node_manager_rank, 1, MPI_INT, client_info.rank_to_node_id_map, 1, MPI_INT,
                  client_cache_world_comm);

    // Create local node map
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

    free(node_world_ranks);

    client_info.local_cache_list_head = NULL;
    client_info.local_cache_list_tail = NULL;
    client_info.cached_item_num       = 0;
    client_info.client_cache_init     = 1;

done:
    fflush(stdout);
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

    // Check if the client cache has been initialized
    if (!client_info.client_cache_init)
        PGOTO_ERROR(FAIL, "pdc_region_dl_insert - object cache list not initialized");

    // Evict item if exceeding cache size
    while (client_info.cached_item_num >= MAX_ITEM_NUM)
        pdc_region_cache_evict();


    obj_cache_item = (pdc_object_data *) PDC_malloc(sizeof(pdc_object_data));
    if (!obj_cache_item)
        PGOTO_ERROR(FAIL, "PDC region cache - obj_cache_item memory allocation failed");

    // Create object data
    obj_cache_item->obj_id       = obj_id;
    obj_cache_item->unit         = unit;
    obj_cache_item->reg_ndim     = ndim;
    obj_cache_item->reg_buf_size = read_size;
    obj_cache_item->target_rank  = -1;

    memcpy(obj_cache_item->reg_offset, offset, ndim * sizeof(uint64_t));
    memcpy(obj_cache_item->reg_size, size, ndim * sizeof(uint64_t));
    memcpy(obj_cache_item->reg_buf, buf, read_size * sizeof(char));

    obj_cache_item->prev = NULL;
    obj_cache_item->next = NULL;

    ret_value = pdc_region_dl_prepend(obj_cache_item);

    pdc_region_cache_timelog(start, "pdc_region_dl_insert - prepend new item to obj_cache list time");
done:
    fflush(stdout);
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
    int              is_cached = 0;
    pdc_object_data *obj_cache_iter;
    double           start, total_start = MPI_Wtime();

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
                memcpy_overlap_subregion(obj_cache_iter->reg_ndim, unit, obj_cache_iter->reg_buf,
                                         obj_cache_iter->reg_offset, obj_cache_iter->reg_size, buf, offset,
                                         size, overlap_offset, overlap_size);

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

    pdc_region_cache_timelog(total_start, "pdc_region_dl_local_search search time");

done:
    fflush(stdout);
    FUNC_LEAVE(is_cached);
}

/*int
pdc_region_dl_node_search(pdcid_t obj_id, int ndim, uint64_t unit, uint64_t *offset, uint64_t *size, void
*buf, uint64_t read_size)
{
    perr_t ret_value = SUCCEED;

    struct pdc_object_list  *obj_cache_iter;
    struct pdc_object_data  *obj_cache_data;
    uint64_t *               overlap_offset, *overlap_size;
    int                      i, is_cached = 0;
    double                   start, total_start = MPI_Wtime();

    FUNC_ENTER(NULL);

    for (int r = 0; r < client_info.node_size; r++) {
        if (r == client_info.node_rank)
            continue;

        for (int s = 0; s < MAX_ITEM_NUM; s++) {
            if (client_info.node_shm_base[r * MAX_ITEM_NUM + s].obj_id == obj_id) {
                obj_cache_data = &client_info.node_shm_base[r * MAX_ITEM_NUM + s];

                is_cached = detect_region_contained(offset, size, obj_cache_data->reg_offset,
                                                    obj_cache_data->reg_size, ndim);

                // If region contained, memcpy cached region data into transfer_request->buf
                if (is_cached) {
                    start = MPI_Wtime();

                    // Detect the offset range that is overlapped
                    PDC_region_overlap_detect(ndim, offset, size, obj_cache_data->reg_offset,
                                              obj_cache_data->reg_size, &overlap_offset, &overlap_size);

                    // memcpy the overlapped region
                    memcpy_overlap_subregion(obj_cache_data->reg_ndim, unit, obj_cache_data->reg_buf,
                                             obj_cache_data->reg_offset, obj_cache_data->reg_size, buf,
offset, size, overlap_offset, overlap_size);

                    free(overlap_offset);
                    pdc_region_cache_timelog(start, "pdc_region_dl_node_search - node cache hit time");

                    break;
                }
            }
        }

        if (is_cached) break;
    }

    pdc_region_cache_timelog(total_start, "pdc_region_dl_node_search search time");

done:
    fflush(stdout);
    FUNC_LEAVE(is_cached);
}*/

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

    pdc_region_cache_timelog(start, "pdc_region_dl_prepare_data_exchange - global prefetch list preparation");

done:
    fflush(stdout);
    FUNC_LEAVE(ret_value);
}

perr_t
pdc_region_dl_data_exchange(pdcid_t *global_prefetch_list, int obj_prefetch_list_len)
{
    FUNC_ENTER(NULL);

    perr_t ret_value = SUCCEED;

    pdc_object_data *obj_cache_iter, *exchange_head;
    char *           intra_node_send_buf, *inter_node_send_buf, *temp_intra_recv_buf, *temp_inter_recv_buf;
    double           start = MPI_Wtime();

    // Step 1: Calculating the max buffer sizes and allocate reusable data exchange buffers
    int    chunk_size           = client_info.cached_item_num / NUM_CHUNKS;
    size_t max_intra_send_chunk = 0, max_inter_send_chunk = 0, local_max_send_chunk = 0;;

    MPI_Datatype mpi_transfer_unit;
    MPI_Type_contiguous(TRANSFER_UNIT_SIZE, MPI_BYTE, &mpi_transfer_unit);
    MPI_Type_commit(&mpi_transfer_unit);

    for (int c = 0; c < NUM_CHUNKS; c++) {
        size_t current_intra_send = 0, current_inter_send = 0;

        int start_idx = c * chunk_size;
        int end_idx   = (c == NUM_CHUNKS - 1) ? client_info.cached_item_num : (c + 1) * chunk_size;

        int i = 0;
        for (obj_cache_iter = client_info.local_cache_list_head; obj_cache_iter != NULL;
             obj_cache_iter = obj_cache_iter->next) {
            if (i >= start_idx && i < end_idx) {
                if (client_info.rank_to_node_id_map[obj_cache_iter->target_rank] ==
                    client_info.node_manager_rank)
                    current_intra_send++;
                else
                    current_inter_send++;
            }
            i++;
        }
        if (current_intra_send > max_intra_send_chunk)
            max_intra_send_chunk = current_intra_send;
        if (current_inter_send > max_inter_send_chunk)
            max_inter_send_chunk = current_inter_send;
    }

    if (max_intra_send_chunk > max_inter_send_chunk)
        local_max_send_chunk = max_intra_send_chunk;
    else
        local_max_send_chunk = max_inter_send_chunk;

    size_t global_max_chunk_size = 0;
    MPI_Allreduce(&local_max_send_chunk, &global_max_chunk_size, 1, MPI_UNSIGNED_LONG, MPI_MAX, client_cache_world_comm);

    intra_node_send_buf = (char *)PDC_malloc(global_max_chunk_size * TRANSFER_UNIT_SIZE + 1);
    inter_node_send_buf = (char *)PDC_malloc(global_max_chunk_size * TRANSFER_UNIT_SIZE + 1);
    temp_intra_recv_buf =
        (char *)PDC_malloc(global_max_chunk_size * client_info.node_size * TRANSFER_UNIT_SIZE + 1);
    temp_inter_recv_buf =
        (char *)PDC_malloc(global_max_chunk_size * client_info.world_size * TRANSFER_UNIT_SIZE + 1);

    MPI_Barrier(client_cache_world_comm);

    // Step 2: Data exchange loop
    exchange_head = client_info.local_cache_list_head;
    for (int c = 0; c < NUM_CHUNKS; c++) {
        int start_idx = 0;
        int end_idx   = chunk_size;
        
        if (c == NUM_CHUNKS - 1) {
            if (client_info.cached_item_num % chunk_size != 0)
                end_idx = client_info.cached_item_num % chunk_size;
        }
        
        // Step 2-1: Intra-node shuffle for current chunk
        int  intra_node_item_count  = 0;
        int *intra_node_send_counts = (int *)PDC_calloc(client_info.node_size, sizeof(int));
        int  i                      = 0;
        for (obj_cache_iter = exchange_head; obj_cache_iter != NULL; obj_cache_iter = obj_cache_iter->next) {
            if (i >= start_idx && i < end_idx) {
                if (client_info.rank_to_node_id_map[obj_cache_iter->target_rank] ==
                    client_info.node_manager_rank) {
                    int target_nrank = client_info.world_to_node_rank_map[obj_cache_iter->target_rank];
                    intra_node_send_counts[target_nrank]++;
                    intra_node_item_count++;
                }
            }
            i++;
        }

        int *intra_node_recv_counts = (int *)PDC_calloc(client_info.node_size, sizeof(int));
        MPI_Alltoall(intra_node_send_counts, 1, MPI_INT, intra_node_recv_counts, 1, MPI_INT,
                     client_cache_node_comm);

        size_t total_intra_recv_items_chunk = 0;
        for (i = 0; i < client_info.node_size; i++) {
            // printf("[RANK %d] data_Exchange - intra_node_send_counts %d intra_node_recv_counts %d\n", pdc_client_mpi_rank_g, intra_node_send_counts[i], intra_node_recv_counts[i], i);
            total_intra_recv_items_chunk += (size_t)intra_node_recv_counts[i];
        }

        int *sdispls_intra     = (int *)PDC_malloc(client_info.node_size * sizeof(int));
        int *rdispls_intra     = (int *)PDC_malloc(client_info.node_size * sizeof(int));
        int current_sdisp_items = 0;
        int current_rdisp_items = 0;
        
        for (i = 0; i < client_info.node_size; i++) {
            sdispls_intra[i] = current_sdisp_items;
            rdispls_intra[i] = current_rdisp_items;

            current_sdisp_items += intra_node_send_counts[i];
            current_rdisp_items += intra_node_recv_counts[i];
        }

        // Packing send buffer
        if (intra_node_item_count > 0) {
            int *temp_offsets_intra = (int *)PDC_calloc(client_info.node_size, sizeof(int));
            int  item_idx           = 0;

            for (obj_cache_iter = exchange_head; obj_cache_iter != NULL;obj_cache_iter = obj_cache_iter->next) {
                if (item_idx >= start_idx && item_idx < end_idx) {
                    int trg_rank = obj_cache_iter->target_rank;

                    if (client_info.rank_to_node_id_map[trg_rank] == client_info.node_manager_rank) {
                        int  trg_nrank   = client_info.world_to_node_rank_map[trg_rank];
                        size_t offset    = ((size_t)sdispls_intra[trg_nrank] + (size_t)temp_offsets_intra[trg_nrank]) * TRANSFER_UNIT_SIZE;
                        char * send_ptr  = intra_node_send_buf + offset;
                        size_t current_offset = 0;

                        memcpy(send_ptr + current_offset, &obj_cache_iter->obj_id, sizeof(pdcid_t));
                        current_offset += sizeof(pdcid_t);
                        memcpy(send_ptr + current_offset, &obj_cache_iter->unit, sizeof(uint64_t));
                        current_offset += sizeof(uint64_t);
                        memcpy(send_ptr + current_offset, &obj_cache_iter->reg_ndim, sizeof(int));
                        current_offset += sizeof(int);
                        memcpy(send_ptr + current_offset, obj_cache_iter->reg_offset, sizeof(uint64_t) * 3);
                        current_offset += (sizeof(uint64_t) * 3);
                        memcpy(send_ptr + current_offset, obj_cache_iter->reg_size, sizeof(uint64_t) * 3);
                        current_offset += (sizeof(uint64_t) * 3);
                        memcpy(send_ptr + current_offset, &obj_cache_iter->reg_buf_size, sizeof(uint64_t));
                        current_offset += sizeof(uint64_t);
                        memcpy(send_ptr + current_offset, obj_cache_iter->reg_buf, MAX_ITEM_SIZE);

                        temp_offsets_intra[trg_nrank]++;
                    }
                }

                item_idx++;
            }
            free(temp_offsets_intra);
        }

        // Collective call for intra node shuffle
        MPI_Alltoallv(intra_node_send_buf, intra_node_send_counts, sdispls_intra, mpi_transfer_unit,
                      temp_intra_recv_buf, intra_node_recv_counts, rdispls_intra, mpi_transfer_unit,
                      client_cache_node_comm);

        // Unpacking recv buffer
        if (total_intra_recv_items_chunk > 0) {
            for (int sender_nrank = 0; sender_nrank < client_info.node_size; sender_nrank++) {
                int items_from_sender = intra_node_recv_counts[sender_nrank];
                char* base_recv_ptr = temp_intra_recv_buf + ((size_t)rdispls_intra[sender_nrank] * TRANSFER_UNIT_SIZE);

                for (int j = 0; j < items_from_sender; ++j) {
                    pdc_object_data *obj_cache_item = (pdc_object_data *)PDC_malloc(sizeof(pdc_object_data));
                    char* recv_ptr = base_recv_ptr + (j * (size_t)TRANSFER_UNIT_SIZE);
                    size_t           current_offset = 0;

                    memcpy(&obj_cache_item->obj_id, recv_ptr + current_offset, sizeof(pdcid_t));
                    current_offset += sizeof(pdcid_t);
                    memcpy(&obj_cache_item->unit, recv_ptr + current_offset, sizeof(uint64_t));
                    current_offset += sizeof(uint64_t);
                    memcpy(&obj_cache_item->reg_ndim, recv_ptr + current_offset, sizeof(int));
                    current_offset += sizeof(int);
                    memcpy(obj_cache_item->reg_offset, recv_ptr + current_offset, sizeof(uint64_t) * 3);
                    current_offset += (sizeof(uint64_t) * 3);
                    memcpy(obj_cache_item->reg_size, recv_ptr + current_offset, sizeof(uint64_t) * 3);
                    current_offset += (sizeof(uint64_t) * 3);
                    memcpy(&obj_cache_item->reg_buf_size, recv_ptr + current_offset, sizeof(uint64_t));
                    current_offset += sizeof(uint64_t);
                    memcpy(obj_cache_item->reg_buf, recv_ptr + current_offset, MAX_ITEM_SIZE);
    
                    obj_cache_item->target_rank = -1;
    
                    pdc_region_dl_prepend(obj_cache_item);
                }
                
            }
        }

        free(sdispls_intra);
        free(rdispls_intra);
        free(intra_node_send_counts);
        free(intra_node_recv_counts);

        // Stage 3: Inter-node data exchange for current chunk
        int  inter_node_item_count  = 0;
        int *inter_node_send_counts = (int *)PDC_calloc(client_info.world_size, sizeof(int));
        i                           = 0;
        
        for (obj_cache_iter = exchange_head; obj_cache_iter != NULL; obj_cache_iter = obj_cache_iter->next) {
            if (i >= start_idx && i < end_idx) {
                if (client_info.rank_to_node_id_map[obj_cache_iter->target_rank] !=
                    client_info.node_manager_rank) {
                    inter_node_send_counts[obj_cache_iter->target_rank]++;
                    inter_node_item_count++;
                }
            }
            i++;
        }

        int *inter_node_recv_counts = (int *)PDC_calloc(client_info.world_size, sizeof(int));
        MPI_Alltoall(inter_node_send_counts, 1, MPI_INT, inter_node_recv_counts, 1, MPI_INT,
                     client_cache_world_comm);

        size_t total_inter_recv_items_chunk = 0;
        for (int j = 0; j < client_info.world_size; j++) {
            total_inter_recv_items_chunk += inter_node_recv_counts[j];
        }

        int *sdispls_inter     = (int *)PDC_malloc(client_info.world_size * sizeof(int));
        int *rdispls_inter     = (int *)PDC_malloc(client_info.world_size * sizeof(int));
        
        current_sdisp_items = 0;
        current_rdisp_items = 0;

        for (int j = 0; j < client_info.world_size; ++j) {
            sdispls_inter[j] = current_sdisp_items;
            rdispls_inter[j] = current_rdisp_items;

            current_sdisp_items += inter_node_send_counts[j];
            current_rdisp_items += inter_node_recv_counts[j];
        }

        if (inter_node_item_count > 0) {
            int *temp_offsets_inter = (int *)PDC_calloc(client_info.world_size, sizeof(int));
            int  item_idx           = 0;
            for (obj_cache_iter = exchange_head; obj_cache_iter != NULL;
                 obj_cache_iter = obj_cache_iter->next) {
                if (item_idx >= start_idx && item_idx < end_idx) {
                    int trg_rank = obj_cache_iter->target_rank;
                    if (client_info.rank_to_node_id_map[trg_rank] != client_info.node_manager_rank) {
                        size_t offset = ((size_t)sdispls_inter[trg_rank] + (size_t)temp_offsets_inter[trg_rank]) * TRANSFER_UNIT_SIZE;
                        char * send_ptr       = inter_node_send_buf + offset;
                        size_t current_offset = 0;
                        
                        memcpy(send_ptr + current_offset, &obj_cache_iter->obj_id, sizeof(pdcid_t));
                        current_offset += sizeof(pdcid_t);
                        memcpy(send_ptr + current_offset, &obj_cache_iter->unit, sizeof(uint64_t));
                        current_offset += sizeof(uint64_t);
                        memcpy(send_ptr + current_offset, &obj_cache_iter->reg_ndim, sizeof(int));
                        current_offset += sizeof(int);
                        memcpy(send_ptr + current_offset, obj_cache_iter->reg_offset, sizeof(uint64_t) * 3);
                        current_offset += (sizeof(uint64_t) * 3);
                        memcpy(send_ptr + current_offset, obj_cache_iter->reg_size, sizeof(uint64_t) * 3);
                        current_offset += (sizeof(uint64_t) * 3);
                        memcpy(send_ptr + current_offset, &obj_cache_iter->reg_buf_size, sizeof(uint64_t));
                        current_offset += sizeof(uint64_t);
                        memcpy(send_ptr + current_offset, obj_cache_iter->reg_buf, MAX_ITEM_SIZE);

                        temp_offsets_inter[trg_rank]++;
                    }
                }
                item_idx++;
            }
            free(temp_offsets_inter);
        }

        MPI_Alltoallv(inter_node_send_buf, inter_node_send_counts, sdispls_inter, mpi_transfer_unit,
                      temp_inter_recv_buf, inter_node_recv_counts, rdispls_inter, mpi_transfer_unit,
                      client_cache_world_comm);

        // Unpack inter-node exchange
        if (total_inter_recv_items_chunk > 0) {
            for (int sender_wrank = 0; sender_wrank < client_info.world_size; sender_wrank++) {
                int items_from_sender = inter_node_recv_counts[sender_wrank];
                char* base_recv_ptr = temp_inter_recv_buf + ((size_t)rdispls_inter[sender_wrank] * TRANSFER_UNIT_SIZE);

                for (int j = 0; j < items_from_sender; ++j) {
                    pdc_object_data *obj_cache_item = (pdc_object_data *) PDC_malloc(sizeof(pdc_object_data));
                    char* recv_ptr = base_recv_ptr + (j * (size_t)TRANSFER_UNIT_SIZE);
                    size_t           current_offset = 0;

                    memcpy(&obj_cache_item->obj_id, recv_ptr + current_offset, sizeof(pdcid_t));
                    current_offset += sizeof(pdcid_t);
                    memcpy(&obj_cache_item->unit, recv_ptr + current_offset, sizeof(uint64_t));
                    current_offset += sizeof(uint64_t);
                    memcpy(&obj_cache_item->reg_ndim, recv_ptr + current_offset, sizeof(int));
                    current_offset += sizeof(int);
                    memcpy(obj_cache_item->reg_offset, recv_ptr + current_offset, sizeof(uint64_t) * 3);
                    current_offset += (sizeof(uint64_t) * 3);
                    memcpy(obj_cache_item->reg_size, recv_ptr + current_offset, sizeof(uint64_t) * 3);
                    current_offset += (sizeof(uint64_t) * 3);
                    memcpy(&obj_cache_item->reg_buf_size, recv_ptr + current_offset, sizeof(uint64_t));
                    current_offset += sizeof(uint64_t);
                    memcpy(obj_cache_item->reg_buf, recv_ptr + current_offset, MAX_ITEM_SIZE);
    
                    pdc_region_dl_prepend(obj_cache_item);
                }
            }
        }

        free(sdispls_inter);
        free(rdispls_inter);
        free(inter_node_send_counts);
        free(inter_node_recv_counts);

        i = 0;
        while (i >= start_idx && i < end_idx && exchange_head != NULL) {
            obj_cache_iter = exchange_head;
            exchange_head = obj_cache_iter->next;
            pdc_region_dl_delete(obj_cache_iter);
            free(obj_cache_iter);
            
            i++;
        }
    }

    // For debugging purpose
    obj_cache_iter = client_info.local_cache_list_head;
    while (obj_cache_iter != NULL) {
        if (client_info.world_rank == 0)
            printf("[RANK %d] object_id %lld offset %lld\n", client_info.world_rank, obj_cache_iter->obj_id, obj_cache_iter->reg_offset[0]);

        fflush(stdout);

        obj_cache_iter = obj_cache_iter->next;
    }

    free(intra_node_send_buf);
    free(inter_node_send_buf);
    free(temp_intra_recv_buf);
    free(temp_inter_recv_buf);

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

                free(obj_cache_item);
            }
        }
    }

done:
    fflush(stdout);
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

    free(obj_cache_item);

done:
    fflush(stdout);
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

    MPI_Comm_free(&client_cache_node_comm);
    MPI_Comm_free(&client_cache_world_comm);

    pdc_region_cache_timelog(start, "pdc_region_dl_finalize - finalization time");

done:
    fflush(stdout);
    FUNC_LEAVE(ret_value);
}
