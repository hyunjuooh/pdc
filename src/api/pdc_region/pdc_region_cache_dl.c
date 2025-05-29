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

#define MAX_CACHE_SIZE 4294967296
// #define MAX_CACHE_SIZE 268435456

MPI_Comm client_cache_comm;

// Initialization flag
int init_object_cache = 0;

// Tracks the number of objects cached within the list
int cached_obj_num = 0;
int max_obj_size = 0;
int global_max_obj_size = 0;

// Object cache list that will be used for object cache management
struct pdc_object_cache *        obj_cache_list, *obj_cache_list_end;

perr_t
pdc_region_dl_init()
{
    perr_t ret_value = SUCCEED;

    int mpi_alloc_error, i;
    double start = MPI_Wtime();

    FUNC_ENTER(NULL);

    // Duplicate MPI_COMM_WORLD for client cache
    MPI_Comm_dup(MPI_COMM_WORLD, &client_cache_comm);

    obj_cache_list = NULL;
    obj_cache_list_end = NULL;

    pdc_region_cache_timelog(start, "pdc_region_dl_init - Total time");

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

    struct pdc_object_cache *obj_cache_item = NULL;
    double start, total_start = MPI_Wtime();

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
    obj_cache_item->target_rank = -1;

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

    start = MPI_Wtime();

    // memcpy region data to obj_cache_list item
    obj_cache_item->buf = (char *)PDC_malloc(sizeof(char) * read_size);
    memcpy(obj_cache_item->buf, buf, sizeof(char) * read_size);

    pdc_region_cache_timelog(start, "pdc_region_dl_insert - memcpy region data to obj_cache_list item time");

    if (obj_cache_list == NULL) {
        DL_PREPEND(obj_cache_list, obj_cache_item);
        obj_cache_list_end = obj_cache_item;
    }
    else {
        DL_PREPEND(obj_cache_list, obj_cache_item);
    }

    cached_obj_num++;

    if (read_size > max_obj_length)
        max_obj_length = read_size;

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
    uint64_t                *overlap_offset, *overlap_size;
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
            is_cached = detect_region_contained(offset, size, obj_cache_iter->reg_offset, obj_cache_iter->reg_size, ndim);

            // If region contained, memcpy cached region data into transfer_request->buf
            if (is_cached) {
                start = MPI_Wtime();
                
                // Detect the offset range that is overlapped
                PDC_region_overlap_detect(ndim, offset, size, obj_cache_iter->reg_offset,
                                          obj_cache_iter->reg_size, &overlap_offset, &overlap_size);

                // memcpy the overlapped region
                memcpy_overlap_subregion(
                    obj_cache_iter->reg_ndim, unit, obj_cache_iter->buf,
                    obj_cache_iter->reg_offset, obj_cache_iter->reg_size, buf, offset, size,
                    overlap_offset, overlap_size);

                // Follow the LRU policy
                // Update the obj_cache_list_end information
                if (obj_cache_iter == obj_cache_list_end) {
                    if (obj_cache_list_end->prev) {
                        obj_cache_list_end = obj_cache_list_end->prev;
                    }
                    else {
                        one_item_obj_list = 1;
                    }
                }

                // Move the recently searched object to the front of the list
                if (!one_item_obj_list) {
                    DL_DELETE(obj_cache_list, obj_cache_iter);
                    DL_PREPEND(obj_cache_list, obj_cache_iter);
                }

                free(overlap_offset);

                pdc_region_cache_timelog(start, "pdc_region_dl_search - local cache hit time");

                break;
            }
        }

        obj_cache_iter = obj_cache_iter->next;
    }

    pdc_region_cache_timelog(total_start, "pdc_region_dl_search - Total time");

done:
    fflush(stdout);
    FUNC_LEAVE(is_cached);
}

perr_t
pdc_region_dl_prepare_data_exchange(char **global_prefetch_list, int *global_list_len, int *global_list_item_len)
{
    perr_t ret_value = SUCCEED;

    int    mpi_alloc_error, i, j, item_idx = 0;
    struct pdc_object_cache *obj_cache_iter;
    double start;

    FUNC_ENTER(NULL);

    // Set the target rank for local cached object item
    obj_cache_iter = obj_cache_list;
    while (obj_cache_iter != NULL) { {
        for (i = 0; i < pdc_client_mpi_size_g; i++) {
            for (j = 0; j < global_list_len[i]; j++) {
                if (strcmp(obj_cache_iter->obj_name, global_prefetch_list[item_idx]) == 0) {
                    obj_cache_iter->target_rank = i;
                    break;
                }
                item_idx++;
            }
            if (obj_cache_iter->target_rank != -1) {
                break;
            }
        }
        obj_cache_iter = obj_cache_iter->next;
        item_idx = 0;
    }

    pdc_region_cache_timelog(start,
                             "pdc_region_dl_collect_global_metadata - global metadata collection time");

done:
    fflush(stdout);
    FUNC_LEAVE(ret_value);
}

int
pdc_region_dl_global_data_exchange(int recv_num)
{
    perr_t                    ret_value = SUCCEED;
    int                       is_cached=0;
    char                     *data_exchange_buf;
    double                    start;

    FUNC_ENTER(NULL);

    // Find the larget object size for data exchange buffer
    MPI_Allreduce(&max_obj_size, &global_max_obj_size, 1, MPI_INT, MPI_MAX, client_cache_comm);

    // Allocate data exchange buffer
    data_exchange_buf = (char *)PDC_malloc(global_max_obj_size);
    if (!data_exchange_buf)
        PGOTO_ERROR(FAIL, "pdc_region_dl_global_data_exchange - data exchange memory allocation failed");

    
    
    

done:
    fflush(stdout);
    FUNC_LEAVE(is_cached);
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

    if (cached_obj_num == 0) {
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
            is_overlapped = check_overlap(ndim, offset, size, obj_cache_item->reg_offset,
                                          obj_cache_item->reg_size);

            // If there is overlapping region remove item from list
            if (is_overlapped) {
                updated_size += obj_cache_item->buf_size;

                if (obj_cache_item == obj_cache_list_end) {
                    if (obj_cache_list_end->prev) {
                        obj_cache_list_end = obj_cache_list_end->prev;
                    }
                }

                // Delete the overlapped object item
                DL_DELETE(obj_cache_list, obj_cache_item);

                free(obj_cache_item->buf);
                free(obj_cache_item);

                cached_obj_num--;
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
    size_t           deleted_size = 0;
    int              deleted_item_num = 0;
    item_delete_info result;

    FUNC_ENTER(NULL);

    if (cached_obj_num == 0) {
        result.deleted_size     = deleted_size;
        result.deleted_item_num = 0;
        return result;
    }

    obj_cache_iter = obj_cache_list_end;

    // Delete item from end of list (following LRU policy)
    while (obj_cache_iter != NULL) {
        obj_cache_item = obj_cache_iter;

        // Update the list end item
        obj_cache_iter     = obj_cache_iter->prev;
        obj_cache_list_end = obj_cache_iter;
        
        required_size -= obj_cache_item->buf_size;
        deleted_size += obj_cache_item->buf_size;

        // Delete the last item of the list and free the buffer
        DL_DELETE(obj_cache_list, obj_cache_item);

        free(obj_cache_item->buf);
        free(obj_cache_item);

        cached_obj_num--;
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
    double start     = MPI_Wtime();
    struct pdc_object_cache *obj_cache_iter, *obj_cache_item;

    FUNC_ENTER(NULL);

    if (!init_object_cache)
        return ret_value;

    // Make sure no shared memory is currently used
    MPI_Barrier(client_cache_comm);

    obj_cache_iter = obj_cache_list;

    while (obj_cache_iter != NULL) {
        obj_cache_item = obj_cache_iter;

        // Update the list end item
        obj_cache_iter     = obj_cache_iter->prev;
        obj_cache_list_end = obj_cache_iter;

        // Delete the last item of the list and free the buffer
        DL_DELETE(obj_cache_list, obj_cache_item);

        free(obj_cache_item->buf);
        free(obj_cache_item);
    }

    MPI_Comm_free(&client_cache_comm);

    init_object_cache = 0;
    cached_obj_num = 0

    pdc_region_cache_timelog(start, "pdc_region_dl_finalize - finalization time");

done:
    fflush(stdout);
    FUNC_LEAVE(ret_value);
}
