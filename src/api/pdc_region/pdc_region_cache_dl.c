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

MPI_Comm client_cache_world_comm;
MPI_Comm client_cache_node_comm;

pdc_client_info client_info;

int pop_free_slot() {
    int free_index = client_info.free_slot_indices[--client_info.free_slot_count];
    
    if (client_info.free_slot_count <= 0) return -1;
    
    return free_index;
}

void push_free_slot(int slot_index) {
    if (client_info.free_slot_count < MAX_ITEM_NUM)
        client_info.free_slot_indices[client_info.free_slot_count++] = slot_index;
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
    
    MPI_Allgather(&client_info.node_manager_rank, 1, MPI_INT, client_info.rank_to_node_id_map, 1, MPI_INT, client_cache_world_comm);
    
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
    client_info.swap_buffer.obj_id = INVALID_ID;
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

    new_slot_idx = pop_free_slot();
    if (new_slot_idx == -1)
        PGOTO_ERROR(FAIL, "pdc_region_dl_insert - No empty slot available for new item");

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
    obj_cache_item->target_rank = INVALID_ID;
    
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
                pdc_region_cache_timelog(start, "pdc_region_dl_search - local cache hit time");

                break;
            }
        }

        obj_cache_iter = obj_cache_iter->next;
    }

    pdc_region_cache_timelog(total_start, "pdc_region_dl_search search time");

done:
    fflush(stdout);
    FUNC_LEAVE(is_cached);
}

perr_t
pdc_region_dl_prepare_data_exchange(char **global_prefetch_list, uint64_t *offset, uint64_t *size,
                                    int *global_list_len)
{
    perr_t ret_value = SUCCEED;
    int    i, j, item_idx, prefetch_list_idx = 0, region_compare = 0, is_contained = 0;
    double start;

    FUNC_ENTER(NULL);

    if (offset != NULL) {
        region_compare = 1;
    }

    // // Set the target rank for local cached object item
    // item_idx = obj_cache_list_metadata->head_idx;

    // while (item_idx != -1) {
    //     for (i = 0; i < pdc_client_mpi_size_g; i++) {
    //         for (j = 0; j < global_list_len[i]; j++) {
    //             // Object name detected in the prefetch list
    //             if (strcmp(obj_cache_list[item_idx].obj_name, global_prefetch_list[prefetch_list_idx]) == 0) {
    //                 // If region_compare = 0, it requires the whole region
    //                 if (!region_compare) {
    //                     obj_cache_list[item_idx].target_rank = i;
    //                     break;
    //                 }
    //                 else {
    //                     is_contained = detect_region_contained(
    //                         offset, size, obj_cache_list[item_idx].reg_offset,
    //                         obj_cache_list[item_idx].reg_size, obj_cache_list[item_idx].reg_ndim);

    //                     // If the prefetch list is contained within the local object cache list
    //                     if (is_contained) {
    //                         obj_cache_list[item_idx].target_rank = i;
    //                         break;
    //                     }
    //                 }
    //             }
    //             prefetch_list_idx++;
    //         }
    //         if (obj_cache_list[item_idx].target_rank != -1) {
    //             break;
    //         }
    //     }
    //     item_idx          = obj_cache_list[item_idx].next;
    //     prefetch_list_idx = 0;
    // }

    pdc_region_cache_timelog(start,
                             "pdc_region_dl_collect_global_metadata - global metadata collection time");

done:
    fflush(stdout);
    FUNC_LEAVE(ret_value);
}

perr_t
pdc_region_dl_data_pack(char *recv_buf, int item_idx, int total_size)
{
    perr_t ret_value = SUCCEED;
    int    position  = 0;

    FUNC_ENTER(NULL);

    // MPI_Pack(obj_cache_list[item_idx].obj_name, MAX_NAME_LEN, MPI_CHAR, recv_buf, total_size, &position,
    //          client_cache_comm);
    // MPI_Pack(obj_cache_list[item_idx].unit, 1, MPI_UINT64_T, recv_buf, total_size, &position,
    //          client_cache_comm);

    // MPI_Pack(obj_cache_list[item_idx].reg_ndim, 1, MPI_INT, recv_buf, total_size, &position,
    //          client_cache_comm);
    // MPI_Pack(obj_cache_list[item_idx].reg_offset, 3, MPI_UINT64_T, recv_buf, total_size, &position,
    //          client_cache_comm);
    // MPI_Pack(obj_cache_list[item_idx].reg_size, 3, MPI_UINT64_T, recv_buf, total_size, &position,
    //          client_cache_comm);

    // // MPI_Pack(obj_cache_list[item_idx].buf_offset, 1, MPI_UINT64_T, recv_buf, total_size, &position,
    // // client_cache_comm);
    // MPI_Pack(obj_cache_list[item_idx].buf_size, 1, MPI_UINT64_T, recv_buf, total_size, &position,
    //          client_cache_comm);

    // if (!region_buf_backward) {
    //     MPI_Pack(region_buf + obj_cache_list[item_idx].buf_offset, obj_cache_list[item_idx].buf_size, MPI_INT,
    //              recv_buf, total_size, &position, client_cache_comm);
    // }
    // else {
    //     MPI_Pack(region_buf_end_ptr - obj_cache_list[item_idx].buf_offset, obj_cache_list[item_idx].buf_size,
    //              MPI_INT, recv_buf, total_size, &position, client_cache_comm);
    // }

done:
    fflush(stdout);
    FUNC_LEAVE(ret_value);
}

perr_t
pdc_region_dl_data_exchange(char **global_prefetch_list, int obj_prefetch_list_len)
{
    perr_t ret_value      = SUCCEED;
    char * receive_buffer = NULL, *send_buffer = NULL;
    int    sends_completed_count    = 0;
    int    receives_completed_count = 0;
    int    sends_initiated_idx      = 0, item_idx;
    double start;

    FUNC_ENTER(NULL);

    // // Set the receive buffer size
    // MPI_Allreduce(&local_max_obj_size, &global_max_obj_size, 1, MPI_UINT64_T, MPI_MAX, client_cache_comm);

    // int total_size = 0, temp_size;

    // MPI_Pack_size(MAX_NAME_LEN, MPI_CHAR, client_cache_comm, &temp_size);
    // total_size += temp_size;
    // MPI_Pack_size(1, MPI_UINT64_T, client_cache_comm, &temp_size);
    // total_size += temp_size;

    // MPI_Pack_size(1, MPI_INT, client_cache_comm, &temp_size);
    // total_size += temp_size;
    // MPI_Pack_size(3, MPI_UINT64_T, client_cache_comm, &temp_size);
    // total_size += temp_size;
    // MPI_Pack_size(3, MPI_UINT64_T, client_cache_comm, &temp_size);
    // total_size += temp_size;

    // // MPI_Pack_size(1, MPI_UINT64_T, client_cache_comm, &temp_size); total_size += temp_size;
    // MPI_Pack_size(1, MPI_UINT64_T, client_cache_comm, &temp_size);
    // total_size += temp_size;

    // MPI_Pack_size(global_max_obj_size, MPI_CHAR, MPI_COMM_WORLD, &temp_size);
    // total_size += temp_size;

    // if (global_max_obj_size > 0) {
    //     send_buffer    = (char *)malloc(total_size);
    //     receive_buffer = (char *)malloc(total_size);
    //     if (!receive_buffer || !send_buffer) {
    //         fprintf(stderr, "[Rank %d] Failed to allocate receive_buffer (size %d)\n", my_rank,
    //                 global_max_obj_size);
    //         MPI_Abort(client_cache_comm, 1);
    //     }
    // }
    // else {
    //     PGOTO_ERROR(FAIL, "pdc_region_dl_data_exchange - global_max_obj_size not specified");
    // }

    // // Post initial receive if data exchange expected
    // if (obj_prefetch_list_len > 0 && global_max_obj_size > 0) {
    //     MPI_Irecv(receive_buffer, global_max_obj_size, MPI_BYTE, MPI_ANY_SOURCE, MPI_ANY_TAG,
    //               client_cache_comm, &current_recv_request);
    //     // printf("[Rank %d] Posted initial MPI_Irecv.\n", my_rank);
    // }

    // int global_all_done = 0;
    // while (!global_all_done) {
    //     int local_sends_done    = (sends_completed_count == obj_cache_list_item_num);
    //     int local_receives_done = (receives_completed_count == obj_prefetch_list_len);

    //     // --- Try to initiate sends ---
    //     sends_initiated_idx = obj_cache_list_metadata->head_idx;
    //     if (!local_sends_done && (sends_initiated_idx != -1)) {
    //         // SegmentInfo *current_send_segment = &segments_to_send[sends_initiated_idx];

    //         // MPI_Pack the node item to send
    //         pdc_region_dl_data_pack(send_buffer, item_idx, total_size);

    //         if (!current_send_segment->is_initiated) {
    //             // printf("[Rank %d] Initiating send for segment ID %d to Rank %d (Size %d, Tag %d)\n",
    //             //       my_rank, current_send_segment->original_segment_id,
    //             //       current_send_segment->partner_rank, current_send_segment->size,
    //             //       current_send_segment->tag);
    //             MPI_Isend(current_send_segment->data, current_send_segment->size, MPI_BYTE,
    //                       current_send_segment->partner_rank, current_send_segment->tag, MPI_COMM_WORLD,
    //                       &current_send_segment->request);
    //             current_send_segment->is_initiated = 1;
    //             sends_initiated_idx++;
    //         }
    //     }

    //     // --- Check for completed sends ---
    //     if (!local_sends_done) {
    //         for (int i = 0; i < sends_initiated_idx; ++i) {
    //             SegmentInfo *seg = &segments_to_send[i];
    //             if (seg->is_initiated && !seg->is_completed) {
    //                 int        flag = 0;
    //                 MPI_Status send_status;
    //                 MPI_Test(&seg->request, &flag, &send_status);
    //                 if (flag) {
    //                     // printf("[Rank %d] Send for segment ID %d to Rank %d completed.\n", my_rank,
    //                     // seg->original_segment_id, seg->partner_rank);
    //                     free(seg->data); // Free the sent data buffer
    //                     seg->data         = NULL;
    //                     seg->is_completed = 1;
    //                     sends_completed_count++;
    //                 }
    //             }
    //         }
    //     }

    //     // --- Check for completed receives ---
    //     if (!local_receives_done && current_recv_request != MPI_REQUEST_NULL) {
    //         int        flag = 0;
    //         MPI_Status recv_status;
    //         MPI_Test(&current_recv_request, &flag, &recv_status);

    //         if (flag) { // A segment has been received
    //             int actual_recv_size;
    //             MPI_Get_count(&recv_status, MPI_BYTE, &actual_recv_size);

    //             // Find which expected segment this corresponds to (for verification)
    //             // This matching is crucial. For simplicity, we assume source/tag identify it.
    //             int matched_expected_idx = -1;
    //             for (int i = 0; i < num_segments_to_receive; ++i) {
    //                 if (!expected_segments_info[i].is_completed &&
    //                     expected_segments_info[i].partner_rank == recv_status.MPI_SOURCE &&
    //                     expected_segments_info[i].tag == recv_status.MPI_TAG) {
    //                     matched_expected_idx = i;
    //                     break;
    //                 }
    //             }

    //             if (matched_expected_idx != -1) {
    //                 expected_segments_info[matched_expected_idx].data = (char *)malloc(actual_recv_size);
    //                 if (!expected_segments_info[matched_expected_idx].data) {
    //                     fprintf(stderr, "[Rank %d] Failed to allocate for received data store\n", my_rank);
    //                     MPI_Abort(MPI_COMM_WORLD, 1);
    //                 }
    //                 memcpy(expected_segments_info[matched_expected_idx].data, receive_buffer,
    //                        actual_recv_size);
    //                 expected_segments_info[matched_expected_idx].size = actual_recv_size; // Store actual size
    //                 expected_segments_info[matched_expected_idx].is_completed =
    //                     1; // Mark as processed from buffer

    //                 printf("[Rank %d] Received segment from Rank %d, Tag %d (Size %d, Checksum %d). Expected "
    //                        "ID %d. Stored.\n",
    //                        my_rank, recv_status.MPI_SOURCE, recv_status.MPI_TAG, actual_recv_size,
    //                        calculate_checksum(expected_segments_info[matched_expected_idx].data,
    //                                           actual_recv_size),
    //                        expected_segments_info[matched_expected_idx].original_segment_id);

    //                 receives_completed_count++;
    //             }
    //             else {
    //                 // This is an unexpected message or a logic error in matching
    //                 fprintf(stderr,
    //                         "[Rank %d] Received UNEXPECTED segment from Rank %d, Tag %d, Size %d. Aborting "
    //                         "or ignoring.\n",
    //                         my_rank, recv_status.MPI_SOURCE, recv_status.MPI_TAG, actual_recv_size);
    //                 // Potentially handle this error, for now, we'll just note it.
    //                 // If not aborting, this message effectively "uses up" the Irecv.
    //             }

    //             current_recv_request = MPI_REQUEST_NULL; // Receive buffer is now free for another Irecv

    //             // Post next receive if more are expected and buffer is free
    //             if (receives_completed_count < num_segments_to_receive && global_max_segment_size > 0) {
    //                 MPI_Irecv(receive_buffer, global_max_segment_size, MPI_BYTE, MPI_ANY_SOURCE, MPI_ANY_TAG,
    //                           MPI_COMM_WORLD, &current_recv_request);
    //                 // printf("[Rank %d] Posted next MPI_Irecv.\n", my_rank);
    //             }
    //         }
    //     }
    //     else if (!local_receives_done && current_recv_request == MPI_REQUEST_NULL &&
    //              num_segments_to_receive > 0 && global_max_segment_size > 0) {
    //         // This case: current_recv_request was null (e.g. after processing one, or initially if first post
    //         // failed/not done) AND we still expect more.
    //         MPI_Irecv(receive_buffer, global_max_segment_size, MPI_BYTE, MPI_ANY_SOURCE, MPI_ANY_TAG,
    //                   MPI_COMM_WORLD, &current_recv_request);
    //         // printf("[Rank %d] Re-posted MPI_Irecv because current was NULL and expecting more.\n",
    //         // my_rank);
    //     }

    //     // --- Global Termination Check ---
    //     local_sends_done    = (sends_completed_count == num_segments_to_send);
    //     local_receives_done = (receives_completed_count == num_segments_to_receive);
    //     int local_all_done  = (local_sends_done && local_receives_done);

    //     MPI_Allreduce(&local_all_done, &global_all_done, 1, MPI_INT, MPI_LAND, MPI_COMM_WORLD);
    // }

    // printf("[Rank %d] Shuffling loop completed. Sent: %d/%d. Received: %d/%d.\n", my_rank,
    //        sends_completed_count, num_segments_to_send, receives_completed_count, num_segments_to_receive);

done:
    fflush(stdout);
    FUNC_LEAVE(ret_value);
}

item_delete_info
pdc_region_dl_update(pdcid_t obj_id, int ndim, uint64_t unit, uint64_t *offset, uint64_t *size, void *buf)
{
    perr_t ret_value = SUCCEED;

    pdc_object_list  *obj_cache_iter, *obj_cache_item;
    pdc_object_data  *obj_cache_data;
    
    size_t                   updated_size  = 0;
    int                      is_overlapped = 0, updated_item_num = 0, one_item_obj_list = 0;
    item_delete_info         result;

    FUNC_ENTER(NULL);

    if (client_info.list_size == 0) {
        result.deleted_size     = updated_size;
        result.deleted_item_num = updated_item_num;
        return result;
    }

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
                updated_size += obj_cache_data->reg_buf_size;

                // Delete the overlapped object item
                ret_value = pdc_region_dl_delete(obj_cache_item);

                obj_cache_data->obj_id = INVALID_ID;
                push_free_slot(obj_cache_item->slot_index);
                free(obj_cache_item);

                client_info.list_size--;
                updated_item_num++;
            }
        }
    }

    result.deleted_size     = updated_size;
    result.deleted_item_num = updated_item_num;

done:
    fflush(stdout);
    FUNC_LEAVE(result);
}

item_delete_info
pdc_region_dl_evict()
{
    perr_t ret_value = SUCCEED;

    pdc_object_list  *obj_cache_item;
    pdc_object_data  *obj_cache_data;
    
    size_t           deleted_size = 0;
    item_delete_info result;

    FUNC_ENTER(NULL);

    if (client_info.list_size == 0) {
        result.deleted_size     = deleted_size;
        result.deleted_item_num = 0;
        return result;
    }

    // Evict the last item following LRU policy
    obj_cache_item = client_info.list_tail;
    obj_cache_data = &client_info.local_shm_base[obj_cache_item->slot_index];

    ret_value = pdc_region_dl_delete(obj_cache_item);

    obj_cache_data->obj_id = INVALID_ID;
    push_free_slot(obj_cache_item->slot_index);
    free(obj_cache_item);

    client_info.list_size--;

    deleted_size += obj_cache_data->reg_buf_size;

    result.deleted_size     = deleted_size;
    result.deleted_item_num = 1;

done:
    fflush(stdout);
    FUNC_LEAVE(result);
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
