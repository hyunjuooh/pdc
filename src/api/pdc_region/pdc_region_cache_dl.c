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

#define MAX_CACHE_SIZE 34359738368
// #define MAX_CACHE_SIZE 4294967296
// #define MAX_CACHE_SIZE 268435456

MPI_Comm client_cache_comm;

// Initialization flag
int init_object_cache;

// Tracks the number of objects cached within the list
int cached_item_num;
int data_exchange_item_num;
int exist_item_num;

// Data exchange
uint64_t global_max_obj_size;
uint64_t local_max_obj_size;

// Object cache list that will be used for object cache management
struct pdc_object_cache *obj_cache_list, *obj_cache_list_end;

perr_t
pdc_region_dl_prepend(struct pdc_object_cache *obj_cache_item)
{
    perr_t ret_value = SUCCEED;

    FUNC_ENTER(NULL);

    // Update metadata list information
    // Current item will be prepended to the head node
    obj_cache_item->next = obj_cache_list;
    obj_cache_item->prev = NULL;

    // If list already exists, current head node's previous index should also point to itself
    if (obj_cache_list != NULL) {
        obj_cache_item->next->prev = obj_cache_item;
    }

    // The new node becomes the head node of the object cache list
    obj_cache_list = obj_cache_item;

    // If list did not exist, update the tail node metadata
    if (obj_cache_list_end == NULL) {
        obj_cache_list_end = obj_cache_item;
    }

done:
    fflush(stdout);
    FUNC_LEAVE(ret_value);
}

perr_t
pdc_region_dl_delete(struct pdc_object_cache *obj_cache_item)
{
    perr_t ret_value = SUCCEED;

    FUNC_ENTER(NULL);

    // Update current item prev item's next item index
    if (obj_cache_item->prev != NULL) {
        obj_cache_item->prev->next = obj_cache_item->next;
    }
    else {
        obj_cache_list = obj_cache_item->next; // if the item is the head
    }

    // Update current item next item's previous item index
    if (obj_cache_item->next != NULL) {
        obj_cache_item->next->prev = obj_cache_item->prev;
    }
    else {
        obj_cache_list_end = obj_cache_item->prev; // if the item is the tail
    }

done:
    fflush(stdout);
    FUNC_LEAVE(ret_value);
}

perr_t
pdc_region_dl_init()
{
    perr_t ret_value = SUCCEED;

    int    mpi_alloc_error, i;
    double start = MPI_Wtime();

    FUNC_ENTER(NULL);

    // Duplicate MPI_COMM_WORLD for client cache
    MPI_Comm_dup(MPI_COMM_WORLD, &client_cache_comm);

    obj_cache_list     = NULL;
    obj_cache_list_end = NULL;

    init_object_cache   = 1;
    cached_item_num     = 0;
    global_max_obj_size = 0;
    local_max_obj_size  = 0;
    exist_item_num      = 0;

    pdc_region_cache_timelog(start, "pdc_region_dl_init - Total time");

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

    struct pdc_object_cache *obj_cache_item = NULL;
    double                   start, total_start = MPI_Wtime();

    FUNC_ENTER(NULL);

    if (!init_object_cache)
        PGOTO_ERROR(FAIL, "pdc_region_dl_insert - object cache list not initialized");

    obj_cache_item = (struct pdc_object_cache *)PDC_malloc(sizeof(struct pdc_object_cache));
    if (!obj_cache_item)
        PGOTO_ERROR(FAIL, "PDC region cache - obj_cache_item memory allocation failed");

    // Create obj_cache_list item
    // obj_cache_list[new_item_idx].obj_id   = obj_id;
    snprintf(obj_cache_item->obj_name, sizeof(obj_cache_item->obj_name), "%s", obj_name);

    obj_cache_item->unit     = unit;
    obj_cache_item->reg_ndim = ndim;
    obj_cache_item->buf_size = read_size;

    obj_cache_item->reg_offset[0] = offset[0];
    obj_cache_item->reg_size[0]   = size[0];

    if (ndim > 1) {
        obj_cache_item->reg_offset[1] = offset[1];
        obj_cache_item->reg_size[1]   = size[1];
    }
    if (ndim > 2) {
        obj_cache_item->reg_offset[2] = offset[2];
        obj_cache_item->reg_size[2]   = size[2];
    }

    obj_cache_item->target_rank  = -1;
    obj_cache_item->is_initiated = 0;
    obj_cache_item->is_completed = 0;
    obj_cache_item->request      = MPI_REQUEST_NULL;

    obj_cache_item->prev = NULL;
    obj_cache_item->next = NULL;

    start = MPI_Wtime();

    // memcpy region data to obj_cache_list item
    obj_cache_item->buf = (char *)PDC_malloc(sizeof(char) * read_size);
    memcpy(obj_cache_item->buf, buf, sizeof(char) * read_size);

    pdc_region_cache_timelog(start, "pdc_region_dl_insert - memcpy region data to obj_cache_list item time");

    ret_value = pdc_region_dl_prepend(obj_cache_item);

    cached_item_num++;

    if (read_size > local_max_obj_size)
        local_max_obj_size = read_size;

    pdc_region_cache_timelog(total_start, "pdc_region_dl_insert - Total time");

done:
    fflush(stdout);
    FUNC_LEAVE(ret_value);
}

int
pdc_region_dl_search(char *obj_name, int ndim, uint64_t unit, uint64_t *offset, uint64_t *size, void *buf,
                     uint64_t read_size)
{
    perr_t ret_value = SUCCEED;

    struct pdc_object_cache *obj_cache_iter;
    uint64_t *               overlap_offset, *overlap_size;
    int                      i, is_cached = 0, one_item_obj_list = 0;
    double                   start, total_start = MPI_Wtime();

    FUNC_ENTER(NULL);

    if (!init_object_cache)
        PGOTO_ERROR(FAIL, "pdc_region_dl_search - object cache list not initialized");

    // Find if region is cached into local object list cache
    // If region contained, return the rank that contains the region
    obj_cache_iter = obj_cache_list;

    while (obj_cache_iter != NULL) {
        if (strcmp(obj_cache_iter->obj_name, obj_name) == 0) {
            is_cached = detect_region_contained(offset, size, obj_cache_iter->reg_offset,
                                                obj_cache_iter->reg_size, ndim);

            // If region contained, memcpy cached region data into transfer_request->buf
            if (is_cached) {
                start = MPI_Wtime();

                // Detect the offset range that is overlapped
                PDC_region_overlap_detect(ndim, offset, size, obj_cache_iter->reg_offset,
                                          obj_cache_iter->reg_size, &overlap_offset, &overlap_size);

                // memcpy the overlapped region
                memcpy_overlap_subregion(obj_cache_iter->reg_ndim, unit, obj_cache_iter->buf,
                                         obj_cache_iter->reg_offset, obj_cache_iter->reg_size, buf, offset,
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

    pdc_region_cache_timelog(total_start, "pdc_region_dl_search - total time");

done:
    fflush(stdout);
    FUNC_LEAVE(is_cached);
}

perr_t
pdc_region_dl_prepare_data_exchange(char **global_prefetch_list, uint64_t *offset, uint64_t *size,
                                    int *global_list_len)
{
    perr_t ret_value = SUCCEED;

    struct pdc_object_cache *obj_cache_iter;
    int                      i, j, prefetch_list_idx = 0, region_compare = 0, is_contained = 0;
    double                   start = MPI_Wtime();

    FUNC_ENTER(NULL);

    MPI_Barrier(client_cache_comm);

    if (offset != NULL)
        region_compare = 1;

    data_exchange_item_num = 0;
    exist_item_num         = 0;

    // Set the target rank for local cached object item
    obj_cache_iter = obj_cache_list;

    while (obj_cache_iter != NULL) {
        for (i = 0; i < pdc_client_mpi_size_g; i++) {
            for (j = 0; j < global_list_len[i]; j++) {
                // Object name detected in the prefetch list
                if (strcmp(obj_cache_iter->obj_name, global_prefetch_list[prefetch_list_idx]) == 0) {
                    // If region_compare = 0, it requires the whole region
                    if (!region_compare) {
                        obj_cache_iter->target_rank = i;
                        break;
                    }
                    else {
                        is_contained =
                            detect_region_contained(&offset[i], &size[i], obj_cache_iter->reg_offset,
                                                    obj_cache_iter->reg_size, obj_cache_iter->reg_ndim);

                        // If the prefetch list is contained within the local object cache list
                        if (is_contained) {
                            obj_cache_iter->target_rank = i;
                            break;
                        }
                    }
                }
                prefetch_list_idx++;
            }

            if (obj_cache_iter->target_rank != -1) {
                data_exchange_item_num++;
                if (obj_cache_iter->target_rank == pdc_client_mpi_rank_g) {
                    exist_item_num++;
                }
                break;
            }
        }

        is_contained      = 0;
        obj_cache_iter    = obj_cache_iter->next;
        prefetch_list_idx = 0;
    }

    // for debugging purpose
    obj_cache_iter = obj_cache_list;
    while (obj_cache_iter != NULL) {
        printf("rank: %d, obj_name: %s, target_rank: %d\n", pdc_client_mpi_rank_g, obj_cache_iter->obj_name,
               obj_cache_iter->target_rank);
        obj_cache_iter = obj_cache_iter->next;
    }

    pdc_region_cache_timelog(start, "pdc_region_dl_prepare_data_exchange - total time");

done:
    fflush(stdout);
    FUNC_LEAVE(ret_value);
}

uint64_t
pdc_region_dl_data_exchange_size()
{
    uint64_t total_size = 0;
    int      temp_size;

    double start = MPI_Wtime();

    FUNC_ENTER(NULL);

    MPI_Pack_size(MAX_NAME_LEN, MPI_CHAR, client_cache_comm, &temp_size);
    total_size += temp_size;
    MPI_Pack_size(1, MPI_UINT64_T, client_cache_comm, &temp_size);
    total_size += temp_size;

    MPI_Pack_size(1, MPI_INT, client_cache_comm, &temp_size);
    total_size += temp_size;
    MPI_Pack_size(3, MPI_UINT64_T, client_cache_comm, &temp_size);
    total_size += temp_size;
    MPI_Pack_size(3, MPI_UINT64_T, client_cache_comm, &temp_size);
    total_size += temp_size;

    MPI_Pack_size(1, MPI_UINT64_T, client_cache_comm, &temp_size);
    total_size += temp_size;
    MPI_Pack_size(global_max_obj_size, MPI_CHAR, MPI_COMM_WORLD, &temp_size);
    total_size += temp_size;

    pdc_region_cache_timelog(start, "pdc_region_dl_data_exchange_size - total time");

done:
    fflush(stdout);
    FUNC_LEAVE(total_size);
}

int
pdc_region_dl_data_pack(char *send_buf, struct pdc_object_cache *obj_cache_item, int total_size)
{
    perr_t ret_value = SUCCEED;

    int    position = 0;
    double start    = MPI_Wtime();

    FUNC_ENTER(NULL);

    MPI_Pack(obj_cache_item->obj_name, MAX_NAME_LEN, MPI_CHAR, send_buf, total_size, &position,
             client_cache_comm);
    MPI_Pack(&(obj_cache_item->unit), 1, MPI_UINT64_T, send_buf, total_size, &position, client_cache_comm);

    MPI_Pack(&(obj_cache_item->reg_ndim), 1, MPI_INT, send_buf, total_size, &position, client_cache_comm);
    MPI_Pack(obj_cache_item->reg_offset, 3, MPI_UINT64_T, send_buf, total_size, &position, client_cache_comm);
    MPI_Pack(obj_cache_item->reg_size, 3, MPI_UINT64_T, send_buf, total_size, &position, client_cache_comm);

    MPI_Pack(&(obj_cache_item->buf_size), 1, MPI_UINT64_T, send_buf, total_size, &position,
             client_cache_comm);
    MPI_Pack(obj_cache_item->buf, obj_cache_item->buf_size, MPI_CHAR, send_buf, total_size, &position,
             client_cache_comm);

    pdc_region_cache_timelog(start, "pdc_region_dl_data_pack - total time");

done:
    fflush(stdout);
    FUNC_LEAVE(position);
}

perr_t
pdc_region_dl_data_unpack(char *recv_buf, struct pdc_object_cache *obj_cache_item, int total_size)
{
    perr_t ret_value = SUCCEED;

    int    position = 0;
    double start    = MPI_Wtime();

    FUNC_ENTER(NULL);

    MPI_Unpack(recv_buf, total_size, &position, obj_cache_item->obj_name, MAX_NAME_LEN, MPI_CHAR,
               client_cache_comm);
    MPI_Unpack(recv_buf, total_size, &position, &(obj_cache_item->unit), 1, MPI_UINT64_T, client_cache_comm);

    MPI_Unpack(recv_buf, total_size, &position, &(obj_cache_item->reg_ndim), 1, MPI_INT, client_cache_comm);
    MPI_Unpack(recv_buf, total_size, &position, obj_cache_item->reg_offset, 3, MPI_UINT64_T,
               client_cache_comm);
    MPI_Unpack(recv_buf, total_size, &position, obj_cache_item->reg_size, 3, MPI_UINT64_T, client_cache_comm);

    MPI_Unpack(recv_buf, total_size, &position, &(obj_cache_item->buf_size), 1, MPI_UINT64_T,
               client_cache_comm);

    obj_cache_item->buf = (char *)PDC_malloc(sizeof(char) * obj_cache_item->buf_size);
    MPI_Unpack(recv_buf, total_size, &position, obj_cache_item->buf, obj_cache_item->buf_size, MPI_CHAR,
               client_cache_comm);

    obj_cache_item->target_rank = -1;

    pdc_region_cache_timelog(start, "pdc_region_dl_data_unpack - total time");

done:
    fflush(stdout);
    FUNC_LEAVE(ret_value);
}

perr_t
pdc_region_dl_data_exchange(int obj_prefetch_list_len)
{
    perr_t ret_value = SUCCEED;

    char *                   receive_buffer = NULL, *send_buffer = NULL;
    struct pdc_object_cache *obj_cache_iter, *obj_cache_item, *recv_obj_item;

    MPI_Request current_recv_request = MPI_REQUEST_NULL;

    int      global_all_done          = 0;
    int      global_send_completed    = 0;
    int      global_receive_completed = 0;
    int      sends_completed_count    = 0;
    int      receives_completed_count = 0;
    
    // int    num_segments_to_receive =  obj_prefetch_list_len - exist_item_num;
    int      num_segments_to_receive = obj_prefetch_list_len;
    uint64_t total_size;
    int      position            = 0;
    int      send_buffer_flag    = 0;
    int      local_sends_done    = (sends_completed_count == data_exchange_item_num);
    int      local_receives_done = (receives_completed_count == num_segments_to_receive);
    int      local_all_done      = (local_sends_done && local_receives_done);

    double start = MPI_Wtime();

    FUNC_ENTER(NULL);

    // Set the receive buffer size
    MPI_Allreduce(&local_max_obj_size, &global_max_obj_size, 1, MPI_UINT64_T, MPI_MAX, client_cache_comm);

    pdc_region_cache_timelog(start, "pdc_region_dl_data_exchange - global_max_obj_size");

    total_size = pdc_region_dl_data_exchange_size();

    if (global_max_obj_size > 0) {
        send_buffer    = (char *)PDC_malloc(total_size);
        receive_buffer = (char *)PDC_malloc(total_size);
        recv_obj_item  = (struct pdc_object_cache *)PDC_malloc(sizeof(struct pdc_object_cache));

        if (!receive_buffer || !send_buffer) {
            fprintf(stderr, "[Rank %d] Failed to allocate receive_buffer (size %d)\n", pdc_client_mpi_rank_g,
                    global_max_obj_size);
            MPI_Abort(client_cache_comm, 1);
        }
    }
    else {
        PGOTO_ERROR(FAIL, "pdc_region_dl_data_exchange - global_max_obj_size not specified");
    }

    pdc_region_cache_timelog(start, "pdc_region_dl_data_exchange - malloc buffer");

    // Post initial receive if data exchange expected
    if (num_segments_to_receive > 0 && global_max_obj_size > 0) {
        MPI_Irecv(receive_buffer, total_size, MPI_PACKED, MPI_ANY_SOURCE, MPI_ANY_TAG, client_cache_comm,
                  &current_recv_request);
    }

    pdc_region_cache_timelog(start, "pdc_region_dl_data_exchange - post initial receive");

    obj_cache_iter = obj_cache_list;

    while (!global_all_done) {
        // while (!local_all_done) {
        // int local_sends_done    = (sends_completed_count == data_exchange_item_num);
        // int local_receives_done = (receives_completed_count == num_segments_to_receive);
        local_sends_done    = (sends_completed_count == data_exchange_item_num);
        local_receives_done = (receives_completed_count == num_segments_to_receive);

        printf("[RANK %d] sends_completed_count: %d, receives_completed_count: %d, data_exchange_item_num: "
               "%d, num_segments_to_receive: %d\n",
               pdc_client_mpi_rank_g, sends_completed_count, receives_completed_count, data_exchange_item_num,
               num_segments_to_receive);
        fflush(stdout);

        // Initiate send
        if (!local_sends_done && (obj_cache_iter != NULL) && (obj_cache_iter->target_rank != -1) &&
            !(obj_cache_iter->is_initiated) && !send_buffer_flag) {
            printf("[RANK %d] obj_name: %s sent, offset: %d\n", pdc_client_mpi_rank_g,
                   obj_cache_iter->obj_name, obj_cache_iter->reg_offset[0]);
            fflush(stdout);

            // MPI_Pack the node item to send
            position = pdc_region_dl_data_pack(send_buffer, obj_cache_iter, total_size);
            MPI_Isend(send_buffer, position, MPI_PACKED, obj_cache_iter->target_rank, 0, client_cache_comm,
                      &obj_cache_iter->request);
            obj_cache_iter->is_initiated = 1;
            send_buffer_flag             = 1;
        }
        pdc_region_cache_timelog(start, "pdc_region_dl_data_exchange - initiate send");

        // printf("[RANK %d] send_buffer_flag: %d, obj_name: %s, initiated: %d\n", pdc_client_mpi_rank_g,
        //         send_buffer_flag, obj_cache_iter->obj_name, obj_cache_iter->is_initiated);
        // fflush(stdout);

        // --- Check for completed send ---
        // If send_buffer is currently used make sure that it was sent completely
        if (send_buffer_flag && (obj_cache_iter != NULL) && obj_cache_iter->is_initiated) {
            int        flag = 0;
            MPI_Status send_status;
            MPI_Test(&obj_cache_iter->request, &flag, &send_status);

            printf("[RANK %d] Check completed sends flag: %d, obj_name: %s\n", pdc_client_mpi_rank_g, flag,
                   obj_cache_iter->obj_name);
            fflush(stdout);

            if (flag) {
                printf("[RANK %d] Send completed obj_name: %s deleted offset: %d\n", pdc_client_mpi_rank_g,
                       obj_cache_iter->obj_name, obj_cache_iter->reg_offset[0]);
                fflush(stdout);

                obj_cache_item = obj_cache_iter;
                obj_cache_iter = obj_cache_iter->next;

                // Delete the item
                ret_value = pdc_region_dl_delete(obj_cache_item);
                pdc_region_cache_delete(obj_cache_item->buf_size, 1);

                free(obj_cache_item->buf);
                free(obj_cache_item);

                cached_item_num--;
                // obj_cache_item->is_completed=1;
                send_buffer_flag = 0;
                sends_completed_count++;
            }
        }

        pdc_region_cache_timelog(start, "pdc_region_dl_data_exchange - check completed send");

        // --- Check for completed receives ---
        // if (!local_receives_done && (current_recv_request != MPI_REQUEST_NULL)) {
        if (!local_receives_done && (current_recv_request != MPI_REQUEST_NULL)) {
            int        flag = 0;
            MPI_Status recv_status;
            MPI_Test(&current_recv_request, &flag, &recv_status);

            if (flag) { // A segment has been received
                ret_value = pdc_region_dl_data_unpack(receive_buffer, recv_obj_item, total_size);

                ret_value = pdc_region_cache_insert(1, recv_obj_item->obj_name, recv_obj_item->reg_ndim,
                                                    recv_obj_item->unit, recv_obj_item->reg_offset,
                                                    recv_obj_item->reg_size, recv_obj_item->buf);

                printf(
                    "[RANK %d] Received obj_name: %s, offset: %d, head_obj_name: %s, head_obj_offset: %d\n",
                    pdc_client_mpi_rank_g, recv_obj_item->obj_name, recv_obj_item->reg_offset[0],
                    obj_cache_list->obj_name, obj_cache_list->reg_offset[0]);
                fflush(stdout);

                receives_completed_count++;
                current_recv_request = MPI_REQUEST_NULL; // Receive buffer is now free for another Irecv

                // Post next receive if more are expected and buffer is free
                if (receives_completed_count < num_segments_to_receive && global_max_obj_size > 0) {
                    MPI_Irecv(receive_buffer, total_size, MPI_PACKED, MPI_ANY_SOURCE, MPI_ANY_TAG,
                              client_cache_comm, &current_recv_request);
                }
            }
            pdc_region_cache_timelog(start, "pdc_region_dl_data_exchange - check completed receives");
        }
        else if (!local_receives_done && current_recv_request == MPI_REQUEST_NULL &&
                 obj_prefetch_list_len > 0 && global_max_obj_size > 0) {
            // This case: current_recv_request was null (e.g. after processing one, or initially if first post
            // failed/not done) AND we still expect more.
            MPI_Irecv(receive_buffer, total_size, MPI_PACKED, MPI_ANY_SOURCE, MPI_ANY_TAG, client_cache_comm,
                      &current_recv_request);
        }

        // --- Global Termination Check ---
        local_sends_done    = (sends_completed_count == data_exchange_item_num);
        local_receives_done = (receives_completed_count == num_segments_to_receive);
        // int local_all_done  = (local_sends_done && local_receives_done);
        local_all_done = (local_sends_done && local_receives_done);

        // MPI_Allreduce(&local_all_done, &global_all_done, 1, MPI_INT, MPI_LAND, client_cache_comm);
        // MPI_Allreduce(&local_sends_done, &global_all_done, 1, MPI_INT, MPI_LAND, client_cache_comm);

        MPI_Allreduce(&data_exchange_item_num, &global_send_completed, 1, MPI_INT, MPI_SUM, client_cache_comm);
        MPI_Allreduce(&receives_completed_count, &global_receive_completed, 1, MPI_INT, MPI_SUM, client_cache_comm);

        global_all_done = (global_send_completed == global_receive_completed);
    }

    printf("[Rank %d] Shuffling loop completed. Sent: %d/%d. Received: %d/%d.\n", pdc_client_mpi_rank_g,
           sends_completed_count, data_exchange_item_num, receives_completed_count, obj_prefetch_list_len);

    MPI_Barrier(client_cache_comm);

    free(receive_buffer);
    free(send_buffer);
    free(recv_obj_item->buf);
    free(recv_obj_item);

done:
    fflush(stdout);
    FUNC_LEAVE(ret_value);
}

item_delete_info
pdc_region_dl_update(char *obj_name, int ndim, uint64_t unit, uint64_t *offset, uint64_t *size, void *buf)
{
    perr_t ret_value = SUCCEED;

    struct pdc_object_cache *obj_cache_iter, *obj_cache_item;
    size_t                   updated_size  = 0;
    int                      is_overlapped = 0, updated_item_num = 0, one_item_obj_list = 0;
    item_delete_info         result;

    FUNC_ENTER(NULL);

    if (cached_item_num == 0) {
        result.deleted_size     = updated_size;
        result.deleted_item_num = updated_item_num;
        return result;
    }

    // Find overlapping regions from the head of the list
    obj_cache_iter = obj_cache_list;

    while (obj_cache_iter != NULL) {
        obj_cache_item = obj_cache_iter;
        obj_cache_iter = obj_cache_iter->next;

        if (strcmp(obj_cache_item->obj_name, obj_name) == 0) {
            // Compare offset and offset + size and see if there is an overlap
            is_overlapped =
                check_overlap(ndim, offset, size, obj_cache_item->reg_offset, obj_cache_item->reg_size);

            // If there is overlapping region remove item from list
            if (is_overlapped) {
                updated_size += obj_cache_item->buf_size;

                // Delete the overlapped object item
                ret_value = pdc_region_dl_delete(obj_cache_item);

                free(obj_cache_item->buf);
                free(obj_cache_item);

                cached_item_num--;
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
pdc_region_dl_evict(size_t required_size)
{
    perr_t ret_value = SUCCEED;

    struct pdc_object_cache *obj_cache_iter, *obj_cache_item;
    size_t                   deleted_size     = 0;
    int                      deleted_item_num = 0;
    item_delete_info         result;

    FUNC_ENTER(NULL);

    if (cached_item_num == 0) {
        result.deleted_size     = deleted_size;
        result.deleted_item_num = 0;
        return result;
    }

    obj_cache_iter = obj_cache_list_end;

    // Delete item from end of list (following LRU policy)
    while (obj_cache_iter != NULL) {
        obj_cache_item = obj_cache_iter;
        obj_cache_iter = obj_cache_iter->prev;

        required_size -= obj_cache_item->buf_size;
        deleted_size += obj_cache_item->buf_size;

        // Delete the last item of the list and free the buffer
        ret_value = pdc_region_dl_delete(obj_cache_item);

        free(obj_cache_item->buf);
        free(obj_cache_item);

        cached_item_num--;
        deleted_item_num++;

        if (required_size < MAX_CACHE_SIZE) {
            break;
        }
    }

    result.deleted_size     = deleted_size;
    result.deleted_item_num = deleted_item_num;

done:
    fflush(stdout);
    FUNC_LEAVE(result);
}

perr_t
pdc_region_dl_finalize()
{
    perr_t ret_value = SUCCEED;

    struct pdc_object_cache *obj_cache_iter, *obj_cache_item;
    int                      i;
    double                   start = MPI_Wtime();

    FUNC_ENTER(NULL);

    if (!init_object_cache)
        return ret_value;

    // Make sure no shared memory is currently used
    MPI_Barrier(client_cache_comm);

    obj_cache_iter = obj_cache_list_end;

    while (obj_cache_iter != NULL) {
        obj_cache_item = obj_cache_iter;
        obj_cache_iter = obj_cache_iter->prev;

        ret_value = pdc_region_dl_delete(obj_cache_item);

        free(obj_cache_item->buf);
        free(obj_cache_item);

        cached_item_num--;
    }

    MPI_Comm_free(&client_cache_comm);

    pdc_region_cache_timelog(start, "pdc_region_dl_finalize - finalization time");

done:
    fflush(stdout);
    FUNC_LEAVE(ret_value);
}
