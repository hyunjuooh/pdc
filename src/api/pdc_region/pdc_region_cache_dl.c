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
#define MAX_CACHE_SIZE 268435456
#define MAX_ITEM_NUM   1000

MPI_Comm client_cache_comm;

int init_object_cache       = 0;
int obj_cache_list_item_num = 0; // Tracks the number of objects cached within the list

// Object cache list that will be used for object cache management
struct pdc_object_cache *        obj_cache_list;
struct pdc_object_list_metadata *obj_cache_list_metadata;

char * metadata_base_ptr = NULL;
size_t metadata_size;

// Data exchange related variables
char *global_metadata_list = NULL;
uint64_t global_max_obj_size = 0;
uint64_t local_max_obj_size = 0;

// For global cache creation
MPI_Win  shared_obj_cache_win = MPI_WIN_NULL;
MPI_Aint buffer_win_size;
// char *   shared_buf;

char *   region_buf;
char *   region_buf_end_ptr;
uint64_t region_buf_offset;
int      region_buf_backward = 0;
int      region_buf_order[MAX_ITEM_NUM];

// Return the next index for obj cache item insertion
int
obj_cache_new_item_index()
{
    perr_t ret_value = SUCCEED;
    int    new_item_idx;

    FUNC_ENTER(NULL);

    if (obj_cache_list_metadata->free_idx == -1)
        PGOTO_ERROR(FAIL, "obj_cache_new_item_index - PDC client cache obtaining new item index failed");

    new_item_idx                      = obj_cache_list_metadata->free_idx;
    obj_cache_list_metadata->free_idx = obj_cache_list[new_item_idx].next;

done:
    fflush(stdout);
    FUNC_LEAVE(new_item_idx);
}

perr_t
pdc_region_dl_list_init()
{
    perr_t ret_value = SUCCEED;

    int    i;
    double start = MPI_Wtime();

    FUNC_ENTER(NULL);

    if (metadata_base_ptr == NULL)
        PGOTO_ERROR(FAIL, "pdc_region_dl_list_init - metadata_base_ptr not initialized");

    // Initialize obj_cache_list_metadata index
    obj_cache_list_metadata->head_idx        = -1;
    obj_cache_list_metadata->tail_idx        = -1;
    obj_cache_list_metadata->free_idx        = 0;
    obj_cache_list_metadata->cached_item_num = 0;

    // Initialize prev, next index for obj_cache_list
    for (i = 0; i < MAX_ITEM_NUM - 1; i++) {
        obj_cache_list[i].prev = -1;
        obj_cache_list[i].next = i + 1;
    }

    // The last item does not have next index
    obj_cache_list[MAX_ITEM_NUM - 1].next = -1;

    // Initialize offset and cached item number
    region_buf_offset       = 0;
    obj_cache_list_item_num = 0;

    // From pdc_region_cache
    total_buf_size = 0;
    total_item_num = 0;

    pdc_region_cache_timelog(start, "pdc_region_dl_list_init - total time");

done:
    fflush(stdout);
    FUNC_LEAVE(ret_value);
}

perr_t
pdc_region_dl_prepend(int item_idx)
{
    perr_t ret_value = SUCCEED;

    FUNC_ENTER(NULL);

    // If it was deleted for prepend
    if (obj_cache_list_metadata->free_idx == item_idx) {
        obj_cache_list_metadata->free_idx = obj_cache_list[item_idx].next;
    }

    // Update metadata list information
    // Current item will be prepended to the head node
    obj_cache_list[item_idx].next = obj_cache_list_metadata->head_idx;
    obj_cache_list[item_idx].prev = -1;

    // If list already exists, current head node's previous index should also point to itself
    if (obj_cache_list_metadata->head_idx != -1) {
        obj_cache_list[obj_cache_list[item_idx].next].prev = item_idx;
    }

    // The new node becomes the head node of the object cache list
    obj_cache_list_metadata->head_idx = item_idx;

    // If list did not exist, update the tail node metadata
    if (obj_cache_list_metadata->tail_idx == -1) {
        obj_cache_list_metadata->tail_idx = item_idx;
    }

done:
    fflush(stdout);
    FUNC_LEAVE(ret_value);
}

perr_t
pdc_region_dl_delete(int item_idx)
{
    perr_t ret_value = SUCCEED;

    FUNC_ENTER(NULL);

    // Update current item prev item's next item index
    if (obj_cache_list[item_idx].prev != -1) {
        obj_cache_list[obj_cache_list[item_idx].prev].next = obj_cache_list[item_idx].next;
    }
    else {
        obj_cache_list_metadata->head_idx = obj_cache_list[item_idx].next; // if the item is the head
    }

    // Update current item next item's previous item index
    if (obj_cache_list[item_idx].next != -1) {
        obj_cache_list[obj_cache_list[item_idx].next].prev = obj_cache_list[item_idx].prev;
    }
    else {
        obj_cache_list_metadata->tail_idx = obj_cache_list[item_idx].prev; // if the item is the tail
    }

    obj_cache_list[item_idx].next     = obj_cache_list_metadata->free_idx;
    obj_cache_list_metadata->free_idx = item_idx;

done:
    fflush(stdout);
    FUNC_LEAVE(ret_value);
}

perr_t
pdc_region_dl_init()
{
    perr_t ret_value = SUCCEED;

    int mpi_alloc_error, i;
    double start;

    FUNC_ENTER(NULL);

    // Duplicate MPI_COMM_WORLD for client cache
    MPI_Comm_dup(MPI_COMM_WORLD, &client_cache_comm);

    // Region buffer memory allocation
    start           = MPI_Wtime();
    buffer_win_size = MAX_CACHE_SIZE * sizeof(char);

    region_buf = (char *)PDC_malloc(MAX_CACHE_SIZE * sizeof(char));
    region_buf_end_ptr = region_buf + (MAX_CACHE_SIZE * sizeof(char));

    pdc_region_cache_timelog(start, "pdc_region_dl_init - Init region memory buffer");

    // mpi_alloc_error = MPI_Win_create(region_buf, buffer_win_size, sizeof(char), MPI_INFO_NULL,
    //                                  client_cache_comm, &shared_obj_cache_win);

    // if (mpi_alloc_error != MPI_SUCCESS)
    //     PGOTO_ERROR(
    //         FAIL,
    //         "pdc_region_dl_collect_global_metadata - MPI shared allocation for shared_obj_cache_win failed");

    pdc_region_cache_timelog(start, "pdc_region_dl_init - create window");

    // Metadata memory allocation
    start = MPI_Wtime();

    metadata_size     = sizeof(pdc_object_list_metadata) + MAX_ITEM_NUM * sizeof(pdc_object_cache);
    metadata_base_ptr = (char *)PDC_malloc(metadata_size);

    obj_cache_list_metadata = (pdc_object_list_metadata *)metadata_base_ptr;
    obj_cache_list = (pdc_object_cache *)((char *)metadata_base_ptr + sizeof(pdc_object_list_metadata));

    pdc_region_cache_timelog(start, "pdc_region_dl_init - Init metadata memory buffer");

    // global_metadata_list = (char *)PDC_malloc(metadata_size * pdc_client_mpi_size_g);
    // if (!global_metadata_list)
    //     PGOTO_ERROR(FAIL, "global metadata memory allocation failed");

    ret_value = pdc_region_dl_list_init();

    init_object_cache = 1;

done:
    fflush(stdout);
    FUNC_LEAVE(ret_value);
}

// TODO:
// Management of overalapping regions
perr_t
pdc_region_dl_insert(char *obj_name, int ndim, uint64_t unit, uint64_t *offset, uint64_t *size, void *buf,
                     uint64_t read_size)
{
    perr_t ret_value = SUCCEED;

    int    new_item_idx;
    double start;

    FUNC_ENTER(NULL);

    if (!init_object_cache)
        PGOTO_ERROR(FAIL, "pdc_region_dl_insert - object cache list not initialized");

    new_item_idx = obj_cache_new_item_index();
    if (new_item_idx == -1)
        PGOTO_ERROR(FAIL, "pdc_region_dl_insert - new item index not found");

    // Create obj_cache_list item
    // obj_cache_list[new_item_idx].obj_id   = obj_id;
    snprintf(obj_cache_list[new_item_idx].obj_name, sizeof(obj_cache_list[new_item_idx].obj_name), "%s",
             obj_name);

    obj_cache_list[new_item_idx].unit     = unit;
    obj_cache_list[new_item_idx].reg_ndim = ndim;
    obj_cache_list[new_item_idx].buf_size = read_size;
    
    // Initialize data exchange metadata info
    obj_cache_list[item_idx].target_rank = -1;
    obj_cache_list[item_idx].is_initiated = 0;
    obj_cache_list[item_idx].is_completed = 0;
    obj_cache_list[item_idx].request = MPI_REQUEST_NULL;
    
    obj_cache_list[new_item_idx].reg_offset[0] = offset[0];
    obj_cache_list[new_item_idx].reg_size[0]   = size[0];

    if (ndim > 1) {
        obj_cache_list[new_item_idx].reg_offset[1] = offset[1];
        obj_cache_list[new_item_idx].reg_size[1]   = size[1];
    }
    if (ndim > 2) {
        obj_cache_list[new_item_idx].reg_offset[2] = offset[2];
        obj_cache_list[new_item_idx].reg_size[2]   = size[2];
    }

    start = MPI_Wtime();

    if (!region_buf_backward) {
        // memcpy region data to obj_cache_list item
        obj_cache_list[new_item_idx].buf_offset = region_buf_offset;
        memcpy(region_buf + region_buf_offset, buf, sizeof(char) * read_size);
        region_buf_offset += (sizeof(char) * read_size);
    } else {
        region_buf_offset += (sizeof(char) * read_size);
        memcpy(region_buf_end_ptr - region_buf_offset, buf, sizeof(char) * read_size);
        obj_cache_list[new_item_idx].buf_offset = region_buf_offset;
    }

    if (local_max_obj_size < read_size) {
        local_max_obj_size = read_size;
    }

    pdc_region_cache_timelog(start, "pdc_region_dl_insert - memcpy region data to obj_cache_list item time");

    start = MPI_Wtime();

    // Prepend the item to the list and update metadata
    ret_value = pdc_region_dl_prepend(new_item_idx);

    pdc_region_cache_timelog(start, "pdc_region_dl_insert - prepend new item to obj_cache list time");

    obj_cache_list_item_num++;
    obj_cache_list_metadata->cached_item_num = obj_cache_list_item_num;

done:
    fflush(stdout);
    FUNC_LEAVE(ret_value);
}

// TODO:
// Need to manage the object's offset cache information
// Need to manage when partial region is contained within
int
pdc_region_dl_search(char *obj_name, int ndim, uint64_t unit, uint64_t *offset, uint64_t *size, void *buf,
                     uint64_t read_size)
{
    perr_t ret_value = SUCCEED;

    uint64_t *overlap_offset, *overlap_size;
    int       i, item_idx, mpi_alloc_error;
    int       is_cached = 0;

    double start = MPI_Wtime();

    FUNC_ENTER(NULL);

    if (!init_object_cache)
        PGOTO_ERROR(FAIL, "pdc_region_dl_search - object cache list not initialized");

    // Find if region is cached into local object list cache
    // If region contained, return the rank that contains the region
    item_idx = obj_cache_list_metadata->head_idx;

    while (item_idx != -1) {
        // if (obj_cache_list[item_idx].obj_id == obj_id) {
        if (strcmp(obj_cache_list[item_idx].obj_name, obj_name) == 0) {
            is_cached = detect_region_contained(offset, size, obj_cache_list[item_idx].reg_offset,
                                                obj_cache_list[item_idx].reg_size, ndim);

            if (is_cached) {
                // If region contained, memcpy cached region data into transfer_request->buf
                // Detect the offset range that is overlapped
                start = MPI_Wtime();
                PDC_region_overlap_detect(ndim, offset, size, obj_cache_list[item_idx].reg_offset,
                                          obj_cache_list[item_idx].reg_size, &overlap_offset, &overlap_size);

                if (!region_buf_backward) {
                    // memcpy the overlapped region
                    memcpy_overlap_subregion(
                        obj_cache_list[item_idx].reg_ndim, unit, region_buf + obj_cache_list[item_idx].buf_offset,
                        obj_cache_list[item_idx].reg_offset, obj_cache_list[item_idx].reg_size, buf, offset, size,
                        overlap_offset, overlap_size);
                } else {
                    memcpy_overlap_subregion(
                        obj_cache_list[item_idx].reg_ndim, unit, region_buf_end_ptr + obj_cache_list[item_idx].buf_offset,
                        obj_cache_list[item_idx].reg_offset, obj_cache_list[item_idx].reg_size, buf, offset, size,
                        overlap_offset, overlap_size);
                }

                // Follow the LRU policy
                pdc_region_dl_delete(item_idx);
                pdc_region_dl_prepend(item_idx);

                free(overlap_offset);

                pdc_region_cache_timelog(start, "pdc_region_dl_search - local cache hit time");

                break;
            }
        }

        item_idx = obj_cache_list[item_idx].next;
    }

    pdc_region_cache_timelog(start, "pdc_region_dl_search search time");

done:
    fflush(stdout);
    FUNC_LEAVE(is_cached);
}

perr_t
pdc_region_dl_prepare_data_exchange(char **global_prefetch_list, uint64_t *offset, uint64_t *size, int *global_list_len)
{
    perr_t                   ret_value = SUCCEED;
    int                      i, j, item_idx, prefetch_list_idx=0, region_compare=0, is_contained=0;
    double                   start;

    FUNC_ENTER(NULL);

    if (offset != NULL) {
        region_compare = 1;
    }

    // Set the target rank for local cached object item
    item_idx = obj_cache_list_metadata->head_idx;
    
    while (item_idx != -1) {
        for (i = 0; i < pdc_client_mpi_size_g; i++) {
            for (j = 0; j < global_list_len[i]; j++) {
                // Object name detected in the prefetch list
                if (strcmp(obj_cache_list[item_idx].obj_name, global_prefetch_list[prefetch_list_idx]) == 0) {
                    // If region_compare = 0, it requires the whole region
                    if (!region_compare) {
                        obj_cache_list[item_idx].target_rank = i;
                        break;
                    } else {
                        is_contained = detect_region_contained(offset, size, obj_cache_list[item_idx].reg_offset,
                                                obj_cache_list[item_idx].reg_size, obj_cache_list[item_idx].reg_ndim);

                        // If the prefetch list is contained within the local object cache list
                        if (is_contained) {
                            obj_cache_list[item_idx].target_rank = i;
                            break;
                        }
                    }
                }
                prefetch_list_idx++;
            }
            if (obj_cache_list[item_idx].target_rank != -1) {
                break;
            }
        }
        item_idx = obj_cache_list[item_idx].next;
        prefetch_list_idx = 0;
    }

    pdc_region_cache_timelog(start,
                             "pdc_region_dl_collect_global_metadata - global metadata collection time");

done:
    fflush(stdout);
    FUNC_LEAVE(ret_value);
}

perr_t
pdc_region_dl_data_pack(char *recv_buf, int item_idx, int total_size) {
    perr_t                   ret_value = SUCCEED;
    int                      position = 0;

    FUNC_ENTER(NULL);

    MPI_Pack(obj_cache_list[item_idx].obj_name, MAX_NAME_LEN, MPI_CHAR, recv_buf, total_size, &position, client_cache_comm);
    MPI_Pack(obj_cache_list[item_idx].unit, 1, MPI_UINT64_T, recv_buf, total_size, &position, client_cache_comm);

    MPI_Pack(obj_cache_list[item_idx].reg_ndim, 1, MPI_INT, recv_buf, total_size, &position, client_cache_comm);
    MPI_Pack(obj_cache_list[item_idx].reg_offset, 3, MPI_UINT64_T, recv_buf, total_size, &position, client_cache_comm);
    MPI_Pack(obj_cache_list[item_idx].reg_size, 3, MPI_UINT64_T, recv_buf, total_size, &position, client_cache_comm);

    // MPI_Pack(obj_cache_list[item_idx].buf_offset, 1, MPI_UINT64_T, recv_buf, total_size, &position, client_cache_comm);
    MPI_Pack(obj_cache_list[item_idx].buf_size, 1, MPI_UINT64_T, recv_buf, total_size, &position, client_cache_comm);

    if (!region_buf_backward) {
        MPI_Pack(region_buf + obj_cache_list[item_idx].buf_offset, obj_cache_list[item_idx].buf_size, MPI_INT, recv_buf, total_size, &position, client_cache_comm);
    } else {
	MPI_Pack(region_buf_end_ptr - obj_cache_list[item_idx].buf_offset, obj_cache_list[item_idx].buf_size, MPI_INT, recv_buf, total_size, &position, client_cache_comm);
    }

done:
    fflush(stdout);
    FUNC_LEAVE(ret_value);
}

perr_t
pdc_region_dl_data_exchange(char **global_prefetch_list, int obj_prefetch_list_len)
{
    perr_t                   ret_value = SUCCEED;
    char                    *receive_buffer = NULL, *send_buffer = NULL;
    int                      sends_completed_count = 0;
    int                      receives_completed_count = 0;
    int                      sends_initiated_idx = 0, item_idx;
    double                   start;

    FUNC_ENTER(NULL);

    // Set the receive buffer size
    MPI_Allreduce(&local_max_obj_size, &global_max_obj_size, 1, MPI_UINT64_T, MPI_MAX, client_cache_comm);

    int total_size = 0, temp_size;

    MPI_Pack_size(MAX_NAME_LEN, MPI_CHAR, client_cache_comm, &temp_size); total_size += temp_size;
    MPI_Pack_size(1, MPI_UINT64_T, client_cache_comm, &temp_size); total_size += temp_size;
    
    MPI_Pack_size(1, MPI_INT, client_cache_comm, &temp_size); total_size += temp_size;
    MPI_Pack_size(3, MPI_UINT64_T, client_cache_comm, &temp_size); total_size += temp_size;
    MPI_Pack_size(3, MPI_UINT64_T, client_cache_comm, &temp_size); total_size += temp_size;

    // MPI_Pack_size(1, MPI_UINT64_T, client_cache_comm, &temp_size); total_size += temp_size;
    MPI_Pack_size(1, MPI_UINT64_T, client_cache_comm, &temp_size); total_size += temp_size;

    MPI_Pack_size(global_max_obj_size, MPI_CHAR, MPI_COMM_WORLD, &temp_size); total_size += temp_size;

    if (global_max_obj_size > 0) {
	send_buffer = (char *)malloc(total_size);
        receive_buffer = (char *)malloc(total_size);
        if (!receive_buffer || !send_buffer) {
            fprintf(stderr, "[Rank %d] Failed to allocate receive_buffer (size %d)\n", my_rank, global_max_obj_size);
            MPI_Abort(client_cache_comm, 1);
        }
    } else {
        PGOTO_ERROR(FAIL, "pdc_region_dl_data_exchange - global_max_obj_size not specified");
    }

    // Post initial receive if data exchange expected
    if (obj_prefetch_list_len > 0 && global_max_obj_size > 0) {
        MPI_Irecv(receive_buffer, global_max_obj_size, MPI_BYTE,
                  MPI_ANY_SOURCE, MPI_ANY_TAG, client_cache_comm, &current_recv_request);
        // printf("[Rank %d] Posted initial MPI_Irecv.\n", my_rank);
    }

    int global_all_done = 0;
    while (!global_all_done) {
        int local_sends_done = (sends_completed_count == obj_cache_list_item_num);
        int local_receives_done = (receives_completed_count == obj_prefetch_list_len);

        // --- Try to initiate sends ---
	sends_initiated_idx = obj_cache_list_metadata->head_idx;
        if (!local_sends_done && (sends_initiated_idx != -1)) {
            // SegmentInfo *current_send_segment = &segments_to_send[sends_initiated_idx];
	    
	    // MPI_Pack the node item to send
	    pdc_region_dl_data_pack(send_buffer, item_idx, total_size);

            if (!current_send_segment->is_initiated) {
                 // printf("[Rank %d] Initiating send for segment ID %d to Rank %d (Size %d, Tag %d)\n",
                 //       my_rank, current_send_segment->original_segment_id, current_send_segment->partner_rank,
                 //       current_send_segment->size, current_send_segment->tag);
                MPI_Isend(current_send_segment->data, current_send_segment->size, MPI_BYTE,
                          current_send_segment->partner_rank, current_send_segment->tag,
                          MPI_COMM_WORLD, &current_send_segment->request);
                current_send_segment->is_initiated = 1;
                sends_initiated_idx++;
            }
        }

        // --- Check for completed sends ---
        if (!local_sends_done) {
            for (int i = 0; i < sends_initiated_idx; ++i) {
                SegmentInfo *seg = &segments_to_send[i];
                if (seg->is_initiated && !seg->is_completed) {
                    int flag = 0;
                    MPI_Status send_status;
                    MPI_Test(&seg->request, &flag, &send_status);
                    if (flag) {
                        // printf("[Rank %d] Send for segment ID %d to Rank %d completed.\n", my_rank, seg->original_segment_id, seg->partner_rank);
                        free(seg->data); // Free the sent data buffer
                        seg->data = NULL;
                        seg->is_completed = 1;
                        sends_completed_count++;
                    }
                }
            }
        }

        // --- Check for completed receives ---
        if (!local_receives_done && current_recv_request != MPI_REQUEST_NULL) {
            int flag = 0;
            MPI_Status recv_status;
            MPI_Test(&current_recv_request, &flag, &recv_status);

            if (flag) { // A segment has been received
                int actual_recv_size;
                MPI_Get_count(&recv_status, MPI_BYTE, &actual_recv_size);

                // Find which expected segment this corresponds to (for verification)
                // This matching is crucial. For simplicity, we assume source/tag identify it.
                int matched_expected_idx = -1;
                for(int i=0; i < num_segments_to_receive; ++i) {
                    if (!expected_segments_info[i].is_completed &&
                        expected_segments_info[i].partner_rank == recv_status.MPI_SOURCE &&
                        expected_segments_info[i].tag == recv_status.MPI_TAG) {
                        matched_expected_idx = i;
                        break;
                    }
                }

                if (matched_expected_idx != -1) {
                    expected_segments_info[matched_expected_idx].data = (char *)malloc(actual_recv_size);
                    if (!expected_segments_info[matched_expected_idx].data) {
                        fprintf(stderr, "[Rank %d] Failed to allocate for received data store\n", my_rank);
                        MPI_Abort(MPI_COMM_WORLD, 1);
                    }
                    memcpy(expected_segments_info[matched_expected_idx].data, receive_buffer, actual_recv_size);
                    expected_segments_info[matched_expected_idx].size = actual_recv_size; // Store actual size
                    expected_segments_info[matched_expected_idx].is_completed = 1; // Mark as processed from buffer

                    printf("[Rank %d] Received segment from Rank %d, Tag %d (Size %d, Checksum %d). Expected ID %d. Stored.\n",
                           my_rank, recv_status.MPI_SOURCE, recv_status.MPI_TAG, actual_recv_size,
                           calculate_checksum(expected_segments_info[matched_expected_idx].data, actual_recv_size),
                           expected_segments_info[matched_expected_idx].original_segment_id);

                    receives_completed_count++;
                } else {
                    // This is an unexpected message or a logic error in matching
                    fprintf(stderr, "[Rank %d] Received UNEXPECTED segment from Rank %d, Tag %d, Size %d. Aborting or ignoring.\n",
                           my_rank, recv_status.MPI_SOURCE, recv_status.MPI_TAG, actual_recv_size);
                    // Potentially handle this error, for now, we'll just note it.
                    // If not aborting, this message effectively "uses up" the Irecv.
                }


                current_recv_request = MPI_REQUEST_NULL; // Receive buffer is now free for another Irecv

                // Post next receive if more are expected and buffer is free
                if (receives_completed_count < num_segments_to_receive && global_max_segment_size > 0) {
                    MPI_Irecv(receive_buffer, global_max_segment_size, MPI_BYTE,
                              MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &current_recv_request);
                    // printf("[Rank %d] Posted next MPI_Irecv.\n", my_rank);
                }
            }
        } else if (!local_receives_done && current_recv_request == MPI_REQUEST_NULL && num_segments_to_receive > 0 && global_max_segment_size > 0) {
             // This case: current_recv_request was null (e.g. after processing one, or initially if first post failed/not done)
             // AND we still expect more.
             MPI_Irecv(receive_buffer, global_max_segment_size, MPI_BYTE,
                       MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &current_recv_request);
             // printf("[Rank %d] Re-posted MPI_Irecv because current was NULL and expecting more.\n", my_rank);
        }


        // --- Global Termination Check ---
        local_sends_done = (sends_completed_count == num_segments_to_send);
        local_receives_done = (receives_completed_count == num_segments_to_receive);
        int local_all_done = (local_sends_done && local_receives_done);

        MPI_Allreduce(&local_all_done, &global_all_done, 1, MPI_INT, MPI_LAND, MPI_COMM_WORLD);
    }

    printf("[Rank %d] Shuffling loop completed. Sent: %d/%d. Received: %d/%d.\n",
           my_rank, sends_completed_count, num_segments_to_send,
           receives_completed_count, num_segments_to_receive);
    
done:
    fflush(stdout);
    FUNC_LEAVE(ret_value);
}


item_delete_info
pdc_region_dl_update(char *obj_name, int ndim, uint64_t unit, uint64_t *offset, uint64_t *size, void *buf)
{
    perr_t ret_value = SUCCEED;

    size_t           updated_size  = 0;
    int              is_overlapped = 0, item_idx, tmp_item_idx, updated_item_num = 0;
    item_delete_info result;

    FUNC_ENTER(NULL);

    if (obj_cache_list_item_num == 0) {
        result.deleted_size     = updated_size;
        result.deleted_item_num = updated_item_num;
        return result;
    }

    // Find overlapping regions from the head of the list
    item_idx = obj_cache_list_metadata->head_idx;

    while (item_idx != -1) {
        tmp_item_idx = obj_cache_list[item_idx].next;
        // if (obj_cache_list[item_idx].obj_id == obj_id) {
        if (strcmp(obj_cache_list[item_idx].obj_name, obj_name) == 0) {
            // Compare offset and offset + size and see if there is an overlap
            is_overlapped = check_overlap(ndim, offset, size, obj_cache_list[item_idx].reg_offset,
                                          obj_cache_list[item_idx].reg_size);

            // If there is overlapping region remove item from list
            if (is_overlapped) {
                updated_size += obj_cache_list[item_idx].buf_size;

                ret_value = pdc_region_dl_delete(item_idx);

                if (ret_value != SUCCEED)
                    PGOTO_ERROR(FAIL, "pdc_region_dl_update - Deleting item failed");

                obj_cache_list_item_num--;
                obj_cache_list_metadata->cached_item_num = obj_cache_list_item_num;
                updated_item_num++;
            }
        }

        item_idx = tmp_item_idx;
    }

    result.deleted_size     = updated_size;
    result.deleted_item_num = updated_item_num;

done:
    fflush(stdout);
    FUNC_LEAVE(result);
}

item_delete_info
pdc_region_dl_evict_by_size(size_t required_size)
{
    perr_t ret_value = SUCCEED;

    size_t           deleted_size = 0;
    int              item_idx, deleted_item_num = 0;
    item_delete_info result;

    FUNC_ENTER(NULL);

    if (obj_cache_list_item_num == 0) {
        result.deleted_size     = deleted_size;
        result.deleted_item_num = 0;
        return result;
    }

    item_idx = obj_cache_list_metadata->tail_idx;

    // Delete item from end of list (following LRU policy)
    while (item_idx != -1) {
        required_size -= obj_cache_list[item_idx].buf_size;
        deleted_size += obj_cache_list[item_idx].buf_size;

        // Update the region_buf_offset for global data buffer
        // Applies only when evicting the last item
        region_buf_offset -= obj_cache_list[item_idx].buf_size;

        ret_value = pdc_region_dl_delete(item_idx);

        if (ret_value != SUCCEED)
            PGOTO_ERROR(FAIL, "pdc_region_dl_evict - Deleting item failed");

        obj_cache_list_item_num--;
        obj_cache_list_metadata->cached_item_num = obj_cache_list_item_num;
        deleted_item_num++;

        if (required_size < MAX_CACHE_SIZE) {
            break;
        }

        item_idx = obj_cache_list_metadata->tail_idx;
    }

    result.deleted_size     = deleted_size;
    result.deleted_item_num = deleted_item_num;

done:
    fflush(stdout);
    FUNC_LEAVE(result);
}

item_delete_info
pdc_region_dl_evict_by_num()
{
    perr_t ret_value = SUCCEED;

    size_t           deleted_size = 0;
    int              item_idx;
    item_delete_info result;

    FUNC_ENTER(NULL);

    if (obj_cache_list_item_num == 0) {
        result.deleted_size     = deleted_size;
        result.deleted_item_num = 0;
        return result;
    }

    // Delete one item from end of list (following LRU policy)
    item_idx = obj_cache_list_metadata->tail_idx;

    // Update the region_buf_offset for global data buffer
    // Applies only when evicting the last item
    region_buf_offset -= obj_cache_list[item_idx].buf_size;

    ret_value = pdc_region_dl_delete(item_idx);

    if (ret_value != SUCCEED)
        PGOTO_ERROR(FAIL, "pdc_region_dl_evict - Deleting item failed");

    deleted_size += obj_cache_list[item_idx].buf_size;
    obj_cache_list_item_num--;

    obj_cache_list_metadata->cached_item_num = obj_cache_list_item_num;

    result.deleted_size     = deleted_size;
    result.deleted_item_num = 1;

done:
    fflush(stdout);
    FUNC_LEAVE(result);
}

perr_t
pdc_region_dl_clean_list()
{
    perr_t   ret_value = SUCCEED;
    int      item_idx;
    char *   tmp_region_buf;
    uint64_t tmp_region_buf_offset = 0;

    FUNC_ENTER(NULL);

    if (obj_cache_list_item_num == 0)
        return ret_value;

    item_idx = obj_cache_list_metadata->head_idx;

    if (init_object_cache) {
        tmp_region_buf = (char *)PDC_malloc(MAX_CACHE_SIZE);

        // Update the region buffer data offset info
        memcpy(tmp_region_buf, region_buf, MAX_CACHE_SIZE);

        region_buf_offset = 0;

        while (item_idx != -1) {
            tmp_region_buf_offset               = obj_cache_list[item_idx].buf_offset;
            obj_cache_list[item_idx].buf_offset = region_buf_offset;

            memcpy(region_buf + region_buf_offset, tmp_region_buf + tmp_region_buf_offset,
                   sizeof(char) * obj_cache_list[item_idx].buf_size);

            region_buf_offset += (sizeof(char) * obj_cache_list[item_idx].buf_size);

            item_idx = obj_cache_list[item_idx].next;
        }

        free(tmp_region_buf);
    }

done:
    fflush(stdout);
    FUNC_LEAVE(ret_value);
}

perr_t
pdc_region_dl_finalize()
{
    perr_t ret_value = SUCCEED;
    double start     = MPI_Wtime();

    FUNC_ENTER(NULL);

    if (!init_object_cache)
        return ret_value;

    // Make sure no shared memory is currently used
    MPI_Barrier(client_cache_comm);

    // ret_value = pdc_region_global_metadata_free();

    // if (shared_obj_cache_win != MPI_WIN_NULL) {
    //     MPI_Win_free(&shared_obj_cache_win);
    //     free(shared_buf);
    //     free(global_metadata_list);
    // }

    if (shared_obj_cache_win != MPI_WIN_NULL) {
        MPI_Win_free(&shared_obj_cache_win);
    }
    
    free(global_metadata_list);
    free(metadata_base_ptr);
    free(region_buf);

    MPI_Comm_free(&client_cache_comm);

    init_object_cache = 0;

    pdc_region_cache_timelog(start, "pdc_region_dl_finalize - finalization time");

done:
    fflush(stdout);
    FUNC_LEAVE(ret_value);
}
