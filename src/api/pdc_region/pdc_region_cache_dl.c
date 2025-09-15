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

// #define MAX_CACHE_SIZE 4294967296
// #define MAX_CACHE_SIZE 268435456
#define MAX_ITEM_SIZE  134214281
#define MAX_ITEM_NUM   256
#define INVALID_ID     -1
#define TAG_PLACEMENT_ORDER 100
#define TAG_LOCAL_PLAN 200

MPI_Comm client_cache_world_comm;
MPI_Comm client_cache_node_comm;

pdc_client_info client_info;

int pop_free_slot() {
    int free_index = client_info.free_slot_indices[--client_info.free_slot_count];
    
    if (client_info.free_slot_count <= 0) return -2;
    
    return free_index;
}

void push_free_slot(int slot_index) {
    if (client_info.free_slot_count < MAX_ITEM_NUM)
        client_info.free_slot_indices[client_info.free_slot_count++] = slot_index;
}

int list_find_by_id(pdcid_t obj_id) {
    pdc_object_list *obj_cache_iter;

    obj_cache_iter = client_info.list_head;
    while (obj_cache_iter!= NULL) {
        if (obj_cache_iter->obj_id == obj_id) {
            return 1;
        }
        obj_cache_iter = obj_cache_iter->next;
    }

    return 0;
}

perr_t
pdc_region_dl_prepend(pdc_object_list* obj_cache_item)
{
    perr_t ret_value = SUCCEED;

    FUNC_ENTER(NULL);

    // Update the obj_cache_item prev, next
    // Current item will be prepended to the head node
    obj_cache_item->next = client_info.list_head;
    obj_cache_item->prev = NULL;

    // If list already exists, current head node's previous index should also point to itself
    if (client_info.list_head != NULL) {
        obj_cache_item->next->prev = obj_cache_item;
    }

    // The new node becomes the head node of the object cache list
    client_info.list_head = obj_cache_item;

    // If list did not exist, update the tail node metadata
    if (client_info.list_tail == NULL) {
        client_info.list_tail = obj_cache_item;
    }

done:
    fflush(stdout);
    FUNC_LEAVE(ret_value);
}

perr_t
pdc_region_dl_delete(pdc_object_list* obj_cache_item)
{
    perr_t ret_value = SUCCEED;

    FUNC_ENTER(NULL);

    // Update current item prev item's next item index
    if (obj_cache_item->prev != NULL) {
        obj_cache_item->prev->next = obj_cache_item->next;
    } else {
        client_info.list_head = obj_cache_item->next; // if the item is the head
    }

    // Update current item next item's previous item index
    if (obj_cache_item->next != NULL) {
        obj_cache_item->next->prev = obj_cache_item->prev;
    } else {
        client_info.list_tail = obj_cache_item->prev; // if the item is the tail
    }

done:
    fflush(stdout);
    FUNC_LEAVE(ret_value);
}



perr_t
pdc_region_dl_init()
{
    perr_t ret_value = SUCCEED; 
    int mpi_alloc_error;

    FUNC_ENTER(NULL);

    client_info.world_rank = pdc_client_mpi_rank_g;
    client_info.world_size = pdc_client_mpi_size_g;
    
    // Duplicate MPI_COMM_WORLD for client cache
    MPI_Comm_dup(MPI_COMM_WORLD, &client_cache_world_comm);

    // Node-level shared memory client_init setup
    MPI_Comm_split_type(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL, &client_cache_node_comm);
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
    
    MPI_Allgather(&client_info.node_manager_rank, 1, MPI_INT, client_info.rank_to_node_id_map, 1, 
                  MPI_INT, client_cache_world_comm);
    
    // Node manager will allocate the shared memory
    client_info.node_shm_base = NULL;
    client_info.node_shm_size = (client_info.node_rank == 0) ? client_info.node_size * MAX_ITEM_NUM * sizeof(pdc_object_data) : 0;
    mpi_alloc_error = MPI_Win_allocate_shared(client_info.node_shm_size, 1, MPI_INFO_NULL, client_cache_node_comm, 
                                              &client_info.node_shm_base, &client_info.node_shm_win);
    if (mpi_alloc_error != MPI_SUCCESS)
        PGOTO_ERROR(FAIL, "pdc_region_dl_init - node shared memory allocation failed");

    // Query the shared memory base pointer
    if (client_info.node_rank != 0) {
        int displ_unit;
        MPI_Win_shared_query(client_info.node_shm_win, 0, &client_info.node_shm_size, &displ_unit, &client_info.node_shm_base);
    }

    // Get current client's shared memory base
    client_info.local_shm_base = &client_info.node_shm_base[client_info.node_rank * MAX_ITEM_NUM];
    for (int i=0; i < MAX_ITEM_NUM; i++) 
        client_info.local_shm_base[i].obj_id = INVALID_ID;

    // Initialize free slot information
    client_info.free_slot_count = MAX_ITEM_NUM;
    for (int i=0; i < MAX_ITEM_NUM; i++) {
        client_info.free_slot_indices[i] = MAX_ITEM_NUM - 1 - i;
    }

    client_info.list_head = NULL;
    client_info.list_tail = NULL;
    client_info.list_size = 0;
    client_info.data_exchange_item_num = 0;
    client_info.client_init = 1;
    
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
    perr_t                  ret_value = SUCCEED;
    
    int                     new_slot_idx;
    pdc_object_data        *obj_cache_data;
    pdc_object_list        *obj_cache_item;
    double                  start;

    FUNC_ENTER(NULL);

    start = MPI_Wtime();

    if (!client_info.client_init)
        PGOTO_ERROR(FAIL, "pdc_region_dl_insert - object cache list not initialized");
    
    do {
        new_slot_idx = pop_free_slot();
        if (new_slot_idx == -2) {
            pdc_region_cache_evict();
            break;
        }
    } while (new_slot_idx == INVALID_ID);

    if (client_info.free_slot_count <= 0) {
	//printf("do statement: new slot idx %d\n", new_slot_idx);
        pdc_region_cache_evict();
        new_slot_idx = pop_free_slot();
    }
    
    obj_cache_data = &client_info.local_shm_base[new_slot_idx];
    if (obj_cache_data->obj_id != INVALID_ID)
        PGOTO_ERROR(FAIL, "pdc_region_dl_insert - The slot is not valid");
    
    // Create object data
    obj_cache_data->obj_id   = obj_id;
    obj_cache_data->unit     = unit;
    obj_cache_data->reg_ndim = ndim;
    obj_cache_data->reg_buf_size = read_size;

    memcpy(obj_cache_data->reg_offset, offset, ndim * sizeof(uint64_t));
    memcpy(obj_cache_data->reg_size, size, ndim * sizeof(uint64_t));
    memcpy(obj_cache_data->reg_buf, buf, read_size * sizeof(char));

    // Add the obj_cache_data and prepend it to the cache list
    obj_cache_item = (struct pdc_object_list *)PDC_malloc(sizeof(struct pdc_object_list));
    if (!obj_cache_item)
        PGOTO_ERROR(FAIL, "pdc_region_dl_insert - obj_cache_item memory allocation failed");
    
    obj_cache_item->obj_id = obj_id;
    obj_cache_item->slot_index = new_slot_idx;
    
    ret_value = pdc_region_dl_prepend(obj_cache_item);

    client_info.list_size++;

    pdc_region_cache_timelog(start, "pdc_region_dl_insert - prepend new item to obj_cache list time");
done:
    fflush(stdout);
    FUNC_LEAVE(ret_value);
}

// TODO:
// Need to manage the object's offset cache information
// Need to manage when partial region is contained within
int
pdc_region_dl_local_search(pdcid_t obj_id, int ndim, uint64_t unit, uint64_t *offset, uint64_t *size, void *buf,
                     uint64_t read_size)
{
    perr_t ret_value = SUCCEED;

    struct pdc_object_list  *obj_cache_iter;
    struct pdc_object_data  *obj_cache_data;
    uint64_t *               overlap_offset, *overlap_size;
    int                      i, is_cached = 0;
    double                   start, total_start = MPI_Wtime();

    FUNC_ENTER(NULL);

    if (!client_info.client_init)
        PGOTO_ERROR(FAIL, "pdc_region_dl_search - object cache list not initialized");

    // Find if region is cached into local object list cache
    // If region contained, return the rank that contains the region
    obj_cache_iter = client_info.list_head;

    while (obj_cache_iter != NULL) {
        if (obj_cache_iter->obj_id == obj_id) {
            obj_cache_data = &client_info.local_shm_base[obj_cache_iter->slot_index];
            
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
                                         obj_cache_data->reg_offset, obj_cache_data->reg_size, buf, offset,
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

int
pdc_region_dl_node_search(pdcid_t obj_id, int ndim, uint64_t unit, uint64_t *offset, uint64_t *size, void *buf,
                     uint64_t read_size)
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
                                             obj_cache_data->reg_offset, obj_cache_data->reg_size, buf, offset,
                                             size, overlap_offset, overlap_size);
    
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
}

perr_t
pdc_region_dl_prepare_data_exchange(pdcid_t *global_prefetch_list, uint64_t *offset, uint64_t *size, 
                                    int obj_prefetch_list_len)
{
    FUNC_ENTER(NULL);

    perr_t ret_value = SUCCEED;

    pdc_object_list  *obj_cache_iter;
    pdc_object_data  *obj_cache_data;
    
    int    i, j, prefetch_list_idx = 0, region_compare = 0, is_contained = 0;
    double start = MPI_Wtime();

    client_info.local_data_exchange = 
                (DataExchangeLocation *)PDC_malloc(client_info.list_size * sizeof(DataExchangeLocation));

    if (offset != NULL) {
        region_compare = 1;
    }

    // Set the target rank for local cached object item
    obj_cache_iter = client_info.list_head; 
    
    while (obj_cache_iter != NULL) {
        for (i = 0; i < client_info.world_size; i++) {
            for (j = 0; j < obj_prefetch_list_len; j++) {
                if (obj_cache_iter->obj_id == global_prefetch_list[prefetch_list_idx]) {
        		    if (!region_compare) {
                        client_info.local_data_exchange[client_info.data_exchange_item_num].obj_id 
                            = obj_cache_iter->obj_id;
                        client_info.local_data_exchange[client_info.data_exchange_item_num].slot_index 
                            = obj_cache_iter->slot_index;
                        client_info.local_data_exchange[client_info.data_exchange_item_num].source_rank 
                            = client_info.world_rank;
                        client_info.local_data_exchange[client_info.data_exchange_item_num].target_rank 
                            = i;
                        
            			break;
        		    }
        		    else {
                        obj_cache_data = &client_info.local_shm_base[obj_cache_iter->slot_index];
                        
            			is_contained = detect_region_contained(&offset[i], &size[i], obj_cache_data->reg_offset,
                    					                       obj_cache_data->reg_size, obj_cache_data->reg_ndim);
			
                		// If the prefetch list is contained within the local object cache list
                        // and if the target_rank is not the same as node_manager_rank
            			if (is_contained) {
                            client_info.local_data_exchange[client_info.data_exchange_item_num].obj_id 
                                = obj_cache_iter->obj_id;
                            client_info.local_data_exchange[client_info.data_exchange_item_num].slot_index 
                                = obj_cache_iter->slot_index;
            			    client_info.local_data_exchange[client_info.data_exchange_item_num].source_rank 
                                = client_info.world_rank;
                            client_info.local_data_exchange[client_info.data_exchange_item_num].target_rank 
                                = i;

                            // printf("[RANK %d] pdc_region_prefetch_by_objid - object list needs to be send to rank %d\n",
                            //         pdc_client_mpi_rank_g, 
                            //         client_info.local_data_exchange[client_info.data_exchange_item_num].target_rank);
                            // fflush(stdout);
                            
            			    break;
            			}
        		    }
        		}
        		prefetch_list_idx++;
    	    }
            
            if (is_contained) {
                client_info.data_exchange_item_num++;

                break;
            }
    	}
        
    	is_contained = 0;
    	obj_cache_iter = obj_cache_iter->next;
    	prefetch_list_idx = 0;
    }

    pdc_region_cache_timelog(start,
                             "pdc_region_dl_prepare_data_exchange - global prefetch list preparation");

done:
    fflush(stdout);
    FUNC_LEAVE(ret_value);
}

// perr_t
// pdc_region_dl_rebuild_list() 
// {
//     FUNC_ENTER(NULL);
    
//     perr_t ret_value = SUCCEED;
//     double start = MPI_Wtime();

//     pdc_object_list  *obj_cache_iter, *obj_cache_item;

//     // Step 1: Initialize free slot information
//     client_info.free_slot_count = MAX_ITEM_NUM;
//     for (int i=0; i < MAX_ITEM_NUM; i++) {
//         client_info.free_slot_indices[i] = MAX_ITEM_NUM - 1 - i;
//     }

//     // Step 2: Free existing list after data exchange
//     obj_cache_iter = client_info.list_tail;

//     while (obj_cache_iter != NULL) {
//         obj_cache_item = obj_cache_iter;
//         obj_cache_iter = obj_cache_iter->prev;

//         ret_value = pdc_region_dl_delete(obj_cache_item);
//         free(obj_cache_item);

//         client_info.list_size--;
//     }
    
//     obj_cache_item = NULL;
    
//     // Step 3: Insert new items
//     for (int i=0; i < MAX_ITEM_NUM; i++) {
//         if (client_info.local_shm_base[i].obj_id != INVALID_ID) { 
//             int tmp_slot_index;
//             // Add the obj_cache_data and prepend it to the cache list
//             obj_cache_item = (struct pdc_object_list *)PDC_malloc(sizeof(struct pdc_object_list));
//             if (!obj_cache_item)
//                 PGOTO_ERROR(FAIL, "pdc_region_dl_rebuild_list - obj_cache_item memory allocation failed");
            
//             obj_cache_item->obj_id = client_info.local_shm_base[i].obj_id;
//             obj_cache_item->slot_index = i;
            
//             ret_value = pdc_region_dl_prepend(obj_cache_item);
//             tmp_slot_index = pop_free_slot();

//             printf("[RANK %d] pdc_region_dl_rebuild_list - tmp_slot_index %d, shared memory index: %d\n",
//                 pdc_client_mpi_rank_g, tmp_slot_index, i);
//             fflush(stdout);
        
//             client_info.list_size++;
//         }
//     }

//     pdc_region_cache_timelog(start, "pdc_region_dl_rebuild_list - total time");
    
// done:
//     fflush(stdout);
//     FUNC_LEAVE(ret_value);
    
// }

perr_t
pdc_region_dl_data_exchange(pdcid_t *global_prefetch_list, int obj_prefetch_list_len)
{
    FUNC_ENTER(NULL);
    
    perr_t ret_value = SUCCEED;
    
    double start = MPI_Wtime();

    MPI_Comm node_manager_comm;
    int node_rank = (client_info.node_rank == 0) ? 0 : MPI_UNDEFINED;
    MPI_Comm_split(client_cache_world_comm, node_rank, client_info.world_rank, &node_manager_comm);

    int num_managers = 0;
    if (node_manager_comm != MPI_COMM_NULL) {
        MPI_Comm_size(node_manager_comm, &num_managers);
        MPI_Win_create(client_info.node_shm_win, client_info.node_size * MAX_ITEM_NUM * sizeof(pdc_object_data), 
                       1, MPI_INFO_NULL, node_manager_comm, &client_info.node_manager_win);

        printf("[RANK %d] pdc_region_dl_data_exchange - node_manager_comm created\n", pdc_client_mpi_rank_g);
        fflush(stdout);
    }

    // Node exchange only done by node managers
    if (client_info.node_rank == 0) {
        printf("[RANK %d] pdc_region_dl_data_exchange - node manager entered\n",
                pdc_client_mpi_rank_g);
        fflush(stdout);

        DataExchangeLocation **node_plans = PDC_malloc(client_info.node_size * sizeof(DataExchangeLocation *));
        int *node_plan_sizes_cnt = (int *) PDC_malloc(client_info.node_size * sizeof(int));
        int total_node_plan_size = client_info.list_size;
        
        node_plans[0] = client_info.local_data_exchange;
        node_plan_sizes_cnt[0] = client_info.list_size;

        // Step 1-1: Aggregate the data exchange plan to the node manager
        for (int i=1; i < client_info.node_size; i++) {
            MPI_Status status;
            int        size_in_bytes;

            MPI_Probe(client_info.world_rank + i, TAG_LOCAL_PLAN, client_cache_world_comm, &status);
            
            MPI_Get_count(&status, MPI_BYTE, &size_in_bytes);
            node_plans[i] = PDC_malloc(size_in_bytes);
            MPI_Recv(node_plans[i], size_in_bytes, MPI_BYTE, client_info.world_rank + i, 
                     TAG_LOCAL_PLAN, client_cache_world_comm, &status);
            
            node_plan_sizes_cnt[i] = size_in_bytes / sizeof(DataExchangeLocation);
            total_node_plan_size += node_plan_sizes_cnt[i];
        }

        printf("[RANK %d] pdc_region_dl_data_exchange - here 1 total_node_plan_size: %d, num_managers: %d\n",
                pdc_client_mpi_rank_g, total_node_plan_size, num_managers);
        fflush(stdout);

        // Step 1-2: Flatten the global plan into a single array
        DataExchangeLocation *node_full_plan = 
                              (DataExchangeLocation *) PDC_malloc(total_node_plan_size * sizeof(DataExchangeLocation));
        
        int current_pos = 0;
        for (int i=0; i < client_info.node_size; i++) {
            memcpy(&node_full_plan[current_pos], node_plans[i], node_plan_sizes_cnt[i] * sizeof(DataExchangeLocation));
            current_pos += node_plan_sizes_cnt[i];
            if (i > 0) free(node_plans[i]);
        }

        free(node_plans); free(node_plan_sizes_cnt);

        printf("[RANK %d] pdc_region_dl_data_exchange - here 2: num_managers: %d\n", pdc_client_mpi_rank_g, num_managers);
        fflush(stdout);

        // Step 1-3: Prepare and share global_plan between node managers only
        // int global_plan_size = num_managers * total_node_plan_size;
        // DataExchangeLocation *global_plan = 
        //                       (DataExchangeLocation *) PDC_malloc(global_plan_size * sizeof(DataExchangeLocation));
        
        // MPI_Allgather(node_full_plan, total_node_plan_size * sizeof(DataExchangeLocation), MPI_BYTE,
        //               global_plan, total_node_plan_size * sizeof(DataExchangeLocation), MPI_BYTE,
        //               node_manager_comm);

        int *recvcounts = (int *) PDC_malloc(num_managers * sizeof(int));
        int *displs = (int *) PDC_malloc(num_managers * sizeof(int));
        MPI_Allgather(&total_node_plan_size, 1, MPI_INT, recvcounts, 1, MPI_INT, node_manager_comm);

        int global_plan_size = 0;
        for (int i=0; i < num_managers; i++) {
            displs[i] = global_plan_size;
            global_plan_size += recvcounts[i];
        }
        
        DataExchangeLocation *global_plan = 
                              (DataExchangeLocation *) PDC_malloc(global_plan_size * sizeof(DataExchangeLocation));

        for (int i=0; i < num_managers; i++) {
            recvcounts[i] *= sizeof(DataExchangeLocation);
            displs[i] *= sizeof(DataExchangeLocation);
        }

        MPI_Allgatherv(node_full_plan, total_node_plan_size * sizeof(DataExchangeLocation), MPI_BYTE,
                       global_plan, recvcounts, displs, MPI_BYTE, node_manager_comm);
        
        free(recvcounts); free(displs); free(node_full_plan);
        
        printf("[RANK %d] pdc_region_dl_data_exchange - here 3\n", pdc_client_mpi_rank_g);
        fflush(stdout);

        // Step 2: Prepare inter-node placement orders
        MPI_Aint vacated_slots_disp[client_info.node_size * MAX_ITEM_NUM];
        // MPI_Aint *vacated_slots_disp = (MPI_Aint *)PDC_malloc(client_info.node_size * MAX_ITEM_NUM * sizeof(MPI_Aint));
        int vacated_cnt = 0;

        // Step 2-1: Calculate the vacated slots wihtin the node
        for (int i=0; i < global_plan_size; i++) {
            if (client_info.rank_to_node_id_map[global_plan[i].source_rank] == client_info.world_rank &&
                client_info.rank_to_node_id_map[global_plan[i].target_rank] != client_info.world_rank) {
                int source_node_rank = global_plan[i].source_rank - client_info.world_rank; // node_rank within the node
                int slot_index = global_plan[i].slot_index;

                vacated_slots_disp[vacated_cnt++] = 
                    (source_node_rank * MAX_ITEM_NUM + slot_index) * sizeof(pdc_object_data);
            }
        }

        printf("[RANK %d] pdc_region_dl_data_exchange - here 4 vacated_cnt: %d\n",
                pdc_client_mpi_rank_g, vacated_cnt);
        fflush(stdout);

        // Step 2-2: Static mapping of the slot_space
        int num_items_to_send = 0;
        
        for (int i=0; i < global_plan_size; i++) {
            if (client_info.rank_to_node_id_map[global_plan[i].source_rank] == client_info.world_rank &&
                client_info.rank_to_node_id_map[global_plan[i].source_rank] != 
                client_info.rank_to_node_id_map[global_plan[i].target_rank]) {
                num_items_to_send++;
            }
        }

        printf("[RANK %d] pdc_region_dl_data_exchange - here 5 num_items_to_send: %d\n",
                pdc_client_mpi_rank_g, num_items_to_send);
        fflush(stdout);

        // Step 2-3: Build a temporary outgoing item buffer
        pdc_object_data *outgoing_buffer = (pdc_object_data *) PDC_malloc(num_items_to_send * sizeof(pdc_object_data));
        if (!outgoing_buffer)
            PGOTO_ERROR(FAIL, "pdc_region_dl_data_exchange - outgoing_buffer memory allocation failed");
        
        int outgoing_idx = 0;

        for (int i=0; i < global_plan_size; i++) {
            if (client_info.rank_to_node_id_map[global_plan[i].source_rank] == client_info.world_rank &&
                client_info.rank_to_node_id_map[global_plan[i].source_rank] != 
                client_info.rank_to_node_id_map[global_plan[i].target_rank]) {
                int source_node_rank = global_plan[i].source_rank - client_info.world_rank;
                pdc_object_data *source_ptr = 
                                &client_info.node_shm_base[source_node_rank * MAX_ITEM_NUM + global_plan[i].slot_index];
                
                memcpy(&outgoing_buffer[outgoing_idx], source_ptr, sizeof(pdc_object_data));
                outgoing_idx++;
            }
        }

        MPI_Barrier(node_manager_comm);

        // Step 2-4: For each item destined for this node, assign it one of the vacaated slot
        MPI_Request requests[global_plan_size];
        int req_cnt = 0;
        
        for (int i=0; i < global_plan_size; i++) {
            if (client_info.rank_to_node_id_map[global_plan[i].source_rank] != client_info.world_rank &&
                client_info.rank_to_node_id_map[global_plan[i].target_rank] == client_info.world_rank) {
                if (vacated_cnt > 0) {
                    MPI_Aint dest_disp = vacated_slots_disp[--vacated_cnt];
                    DataPlacementOrder order = {global_plan[i].obj_id, global_plan[i].source_rank, global_plan[i].slot_index, global_plan[i].target_rank, dest_disp};
                    int source_manager = client_info.rank_to_node_id_map[global_plan[i].source_rank];
                    
                    MPI_Isend(&order, sizeof(DataPlacementOrder), MPI_BYTE, source_manager, 
                              TAG_PLACEMENT_ORDER, client_cache_world_comm, &requests[req_cnt++]);
                } 
            }
        }

        MPI_Waitall(req_cnt, requests, MPI_STATUSES_IGNORE);
        MPI_Barrier(node_manager_comm);

        printf("[RANK %d] pdc_region_dl_data_exchange - here 6\n", pdc_client_mpi_rank_g);
        fflush(stdout);

        int *manager_world_ranks = (int *) PDC_malloc(num_managers * sizeof(int));
        MPI_Allgather(&client_info.world_rank, 1, MPI_INT, manager_world_ranks, 1, MPI_INT, node_manager_comm);

        // Step 3: Inter-node data movement
        DataPlacementOrder *orders = (DataPlacementOrder *) PDC_malloc(num_items_to_send * sizeof(DataPlacementOrder));
        
        for (int i=0; i < num_items_to_send; i++) {
            MPI_Status status;
            
            MPI_Recv(&orders[i], sizeof(DataPlacementOrder), MPI_BYTE, MPI_ANY_SOURCE, TAG_PLACEMENT_ORDER, 
                     client_cache_world_comm, &status);
        }

        // MPI_Win_fence(0, client_info.node_manager_win);
        
        for (int i=0; i < num_items_to_send; i++) {
            DataPlacementOrder *order = &orders[i];

            pdc_object_data *swap_buffer = NULL;
            for (int j=0; j < num_items_to_send; j++) {
                if (outgoing_buffer[j].obj_id == order->obj_id) {
                    swap_buffer = &outgoing_buffer[j];
                    break;
                }
            }

            printf("[RANK %d] pdc_region_dl_data_exchange - Inter-node data movement: object_id %d\n",
                pdc_client_mpi_rank_g, swap_buffer->obj_id);
            fflush(stdout);

            if (swap_buffer) {
                int target_manager_rank = -1;
                for (int j=0; j < num_managers; j++) {
                    if(manager_world_ranks[j] == client_info.rank_to_node_id_map[order->target_rank]) {
                        target_manager_rank = j;
                        break;
                    }
                }

                printf("[RANK %d] pdc_region_dl_data_exchange - Inter-node data movement: target_manager_rank %d\n",
                pdc_client_mpi_rank_g, target_manager_rank);
                fflush(stdout);

                if (target_manager_rank != -1) {
                    MPI_Win_lock(MPI_LOCK_SHARED, target_manager_rank, 0, client_info.node_manager_win);
                    MPI_Put(swap_buffer, sizeof(pdc_object_data), MPI_BYTE, target_manager_rank, order->target_disp, 
                            sizeof(pdc_object_data), MPI_BYTE, client_info.node_manager_win);
                    MPI_Win_unlock(target_manager_rank, client_info.node_manager_win);
                }
            }
        }

        // MPI_Win_fence(0, client_info.node_manager_win);

        printf("[RANK %d] pdc_region_dl_data_exchange - here 7\n",
                pdc_client_mpi_rank_g);
        fflush(stdout);

        if (node_manager_comm != MPI_COMM_NULL) {
            MPI_Barrier(node_manager_comm);
        // MPI_Comm_free(node_manager_comm);
        }
        
        free(global_plan);
        free(outgoing_buffer);
        free(manager_world_ranks);
        free(orders);
    } else {
        printf("[RANK %d] pdc_region_dl_data_exchange - workers entered\n",
                pdc_client_mpi_rank_g);
        fflush(stdout);
        
        MPI_Send(client_info.local_data_exchange, client_info.list_size * sizeof(DataExchangeLocation), MPI_BYTE,
                 client_info.node_manager_rank, TAG_LOCAL_PLAN, client_cache_world_comm);
    }

    MPI_Barrier(client_cache_world_comm);

    printf("[RANK %d] pdc_region_dl_data_exchange - here 7\n",
            pdc_client_mpi_rank_g);
    fflush(stdout);

    // if (node_manager_comm != MPI_COMM_NULL) {
    //     // MPI_Barrier(node_manager_comm);
    //     MPI_Comm_free(node_manager_comm);
    // }

    printf("[RANK %d] pdc_region_dl_data_exchange - here 8\n",
            pdc_client_mpi_rank_g);
    fflush(stdout);

    //ret_value = pdc_region_dl_rebuild_list();

    free(client_info.local_data_exchange);
    
    pdc_region_cache_timelog(start, "pdc_region_dl_data_exchange - total time");
    
done:
    fflush(stdout);
    FUNC_LEAVE(ret_value);
}

perr_t
pdc_region_dl_update(pdcid_t obj_id, int ndim, uint64_t unit, uint64_t *offset, uint64_t *size, void *buf)
{
    FUNC_ENTER(NULL);
    
    perr_t ret_value = SUCCEED;

    pdc_object_list  *obj_cache_iter, *obj_cache_item;
    pdc_object_data  *obj_cache_data;
    
    int                      is_overlapped = 0;

    if (client_info.list_size == 0)
        goto done;

    // Find overlapping regions from the head of the list
    obj_cache_iter = client_info.list_head;

    while (obj_cache_iter != NULL) {
        obj_cache_item = obj_cache_iter;
        obj_cache_iter = obj_cache_iter->next;

        if (obj_cache_item->obj_id == obj_id) {
            obj_cache_data = &client_info.local_shm_base[obj_cache_item->slot_index];
            
            // Compare offset and offset + size and see if there is an overlap
            is_overlapped =
                check_overlap(ndim, offset, size, obj_cache_data->reg_offset, obj_cache_data->reg_size);

            // If there is overlapping region remove item from list
            if (is_overlapped) {
                // Delete the overlapped object item
                ret_value = pdc_region_dl_delete(obj_cache_item);

                obj_cache_data->obj_id = INVALID_ID;
                push_free_slot(obj_cache_item->slot_index);
                free(obj_cache_item);

                client_info.list_size--;
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
    perr_t ret_value = SUCCEED;

    pdc_object_list  *obj_cache_item;
    pdc_object_data  *obj_cache_data;

    FUNC_ENTER(NULL);

    if (client_info.list_size == 0) {
        goto done;
    }

    // Evict the last item following LRU policy
    obj_cache_item = client_info.list_tail;
    obj_cache_data = &client_info.local_shm_base[obj_cache_item->slot_index];

    ret_value = pdc_region_dl_delete(obj_cache_item);

    obj_cache_data->obj_id = INVALID_ID;
    push_free_slot(obj_cache_item->slot_index);
    free(obj_cache_item);

    client_info.list_size--;
    
done:
    fflush(stdout);
    FUNC_LEAVE(ret_value);
}

perr_t
pdc_region_dl_finalize()
{
    perr_t             ret_value = SUCCEED;
    pdc_object_list *obj_cache_iter, *obj_cache_item;
    double             start     = MPI_Wtime();

    FUNC_ENTER(NULL);

    if (!client_info.client_init)
        return ret_value;

    MPI_Barrier(client_cache_world_comm);   // Make sure no shared memory is currently used

    obj_cache_iter = client_info.list_tail; // Free the obj_cache_list of the client_info

    while (obj_cache_iter != NULL) {
        obj_cache_item = obj_cache_iter;
        obj_cache_iter = obj_cache_iter->prev;

        ret_value = pdc_region_dl_delete(obj_cache_item);

        free(obj_cache_item);

        client_info.list_size--;
    }

    if (client_info.node_shm_win != MPI_WIN_NULL) {
        MPI_Win_free(&client_info.node_shm_win); // Free the shared memory of the client_info
    }

    free(client_info.rank_to_node_id_map);  // Free the rank_to_node_id_map

    MPI_Comm_free(&client_cache_node_comm);
    MPI_Comm_free(&client_cache_world_comm);

    pdc_region_cache_timelog(start, "pdc_region_dl_finalize - finalization time");

done:
    fflush(stdout);
    FUNC_LEAVE(ret_value);
}
