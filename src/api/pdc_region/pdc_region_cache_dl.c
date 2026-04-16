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
#define MAX_CACHE_SIZE 34359738368
#define MAX_ITEM_NUM   1000

MPI_Comm client_cache_comm;

int init_object_cache       = 0;
int obj_cache_list_item_num = 0; // Tracks the number of objects cached within the list

// Object cache list that will be used for object cache management
struct pdc_object_cache *        obj_cache_list;
struct pdc_object_list_metadata *obj_cache_list_metadata;

char * metadata_base_ptr = NULL;
size_t metadata_size;

char *global_metadata_list = NULL;

// For global cache creation
MPI_Win  shared_obj_cache_win = MPI_WIN_NULL;
MPI_Aint buffer_win_size;
char *   shared_buf;
char *   region_buf;
uint64_t region_buf_offset;

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
    // MPI_Aint buffer_win_size;

    double start;

    FUNC_ENTER(NULL);

    // Duplicate MPI_COMM_WORLD for client cache
    MPI_Comm_dup(MPI_COMM_WORLD, &client_cache_comm);

    // Region buffer memory allocation
    start           = MPI_Wtime();
    buffer_win_size = MAX_CACHE_SIZE * sizeof(char);

    region_buf = (char *)PDC_malloc(MAX_CACHE_SIZE * sizeof(char));

    pdc_region_cache_timelog(start, "pdc_region_dl_init - Init region memory buffer");

    // Metadata memory allocation
    start = MPI_Wtime();

    metadata_size     = sizeof(pdc_object_list_metadata) + MAX_ITEM_NUM * sizeof(pdc_object_cache);
    metadata_base_ptr = (char *)PDC_malloc(metadata_size);

    obj_cache_list_metadata = (pdc_object_list_metadata *)metadata_base_ptr;
    obj_cache_list = (pdc_object_cache *)((char *)metadata_base_ptr + sizeof(pdc_object_list_metadata));

    pdc_region_cache_timelog(start, "pdc_region_dl_init - Init metadata memory buffer");

    ret_value = pdc_region_dl_list_init();

    init_object_cache = 1;

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
    obj_cache_list[new_item_idx].obj_id = obj_id;
    // snprintf(obj_cache_list[new_item_idx].obj_name, sizeof(obj_cache_list[new_item_idx].obj_name), "%s",
    //          obj_name);

    obj_cache_list[new_item_idx].unit     = unit;
    obj_cache_list[new_item_idx].reg_ndim = ndim;
    obj_cache_list[new_item_idx].buf_size = read_size;

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

    // memcpy region data to obj_cache_list item
    obj_cache_list[new_item_idx].buf_offset = region_buf_offset;
    memcpy(region_buf + region_buf_offset, buf, sizeof(char) * read_size);
    region_buf_offset += (sizeof(char) * read_size);

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
pdc_region_dl_search(pdcid_t obj_id, int ndim, uint64_t unit, uint64_t *offset, uint64_t *size, void *buf,
                     uint64_t read_size)
{
    perr_t ret_value = SUCCEED;

    uint64_t *overlap_offset, *overlap_size;
    int       i, item_idx, mpi_alloc_error;
    int       is_cached = 0;

    int region_copy = 1;

    double start = MPI_Wtime();

    FUNC_ENTER(NULL);

    if (!init_object_cache)
        PGOTO_ERROR(FAIL, "pdc_region_dl_search - object cache list not initialized");

    // Find if region is cached into local object list cache
    // If region contained, return the rank that contains the region
    item_idx = obj_cache_list_metadata->head_idx;

    while (item_idx != -1) {
        if (obj_cache_list[item_idx].obj_id == obj_id) {
            // if (strcmp(obj_cache_list[item_idx].obj_name, obj_name) == 0) {
            is_cached = detect_region_contained(offset, size, obj_cache_list[item_idx].reg_offset,
                                                obj_cache_list[item_idx].reg_size, ndim);

            if (is_cached) {
                // If region contained, memcpy cached region data into transfer_request->buf
                // Detect the offset range that is overlapped
                start = MPI_Wtime();
                PDC_region_overlap_detect(ndim, offset, size, obj_cache_list[item_idx].reg_offset,
                                          obj_cache_list[item_idx].reg_size, &overlap_offset, &overlap_size);

                for (int i = 0; i < ndim; ++i) {
                    if (offset[i] != obj_cache_list[item_idx].reg_offset[i]) {
                        region_copy = 0;
                        break;
                    }
                    if (size[i] != obj_cache_list[item_idx].reg_size[i]) {
                        region_copy = 0;
                        break;
                    }
                }

                if (region_copy) {
                    memcpy(buf, region_buf + obj_cache_list[item_idx].buf_offset,
                           obj_cache_list[item_idx].buf_size);
                    // if (client_info.world_rank == 0)
                    //     printf("[RANK %d] Read entire region for obj_id: %lld\n", client_info.world_rank,
                    //            obj_id);
                }
                else {
                    // memcpy the overlapped region
                    memcpy_overlap_subregion(
                        obj_cache_list[item_idx].reg_ndim, unit,
                        region_buf + obj_cache_list[item_idx].buf_offset, obj_cache_list[item_idx].reg_offset,
                        obj_cache_list[item_idx].reg_size, buf, offset, size, overlap_offset, overlap_size);
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
pdc_region_dl_collect_global_metadata()
{
    perr_t ret_value = SUCCEED;

    int    mpi_alloc_error, item_idx;
    double start;

    FUNC_ENTER(NULL);

    // If the shared_buf was not freed prior to new allocation free that first
    if (shared_buf != NULL)
        ret_value = pdc_region_global_metadata_free();

    // Create window to expose local cache to other ranks
    start      = MPI_Wtime();
    shared_buf = (char *)PDC_malloc(MAX_CACHE_SIZE * sizeof(char));

    mpi_alloc_error = MPI_Win_create(shared_buf, buffer_win_size, sizeof(char), MPI_INFO_NULL,
                                     client_cache_comm, &shared_obj_cache_win);

    if (mpi_alloc_error != MPI_SUCCESS)
        PGOTO_ERROR(
            FAIL,
            "pdc_region_dl_collect_global_metadata - MPI shared allocation for shared_obj_cache_win failed");

    pdc_region_cache_timelog(start, "pdc_region_dl_collect_global_metadata - create window");

    // Copy local region into the shared window
    start = MPI_Wtime();
    memcpy(shared_buf, region_buf, MAX_CACHE_SIZE);
    pdc_region_cache_timelog(start, "pdc_region_dl_collect_global_metadata - memcpy time");

    // Shared global metadata memory allocation
    start                = MPI_Wtime();
    global_metadata_list = (char *)PDC_malloc(metadata_size * pdc_client_mpi_size_g);
    if (!global_metadata_list)
        PGOTO_ERROR(FAIL, "global metadata memory allocation failed");

    // Collect the metadata information from all ranks
    mpi_alloc_error = MPI_Allgather(metadata_base_ptr, metadata_size, MPI_CHAR, global_metadata_list,
                                    metadata_size, MPI_CHAR, client_cache_comm);
    if (mpi_alloc_error != MPI_SUCCESS)
        PGOTO_ERROR(FAIL, "pdc_region_dl_search - MPI Allgather failed");

    pdc_region_cache_timelog(start,
                             "pdc_region_dl_collect_global_metadata - global metadata collection time");

done:
    fflush(stdout);
    FUNC_LEAVE(ret_value);
}

perr_t
pdc_region_global_metadata_free()
{
    perr_t ret_value = SUCCEED;
    double start     = MPI_Wtime();

    FUNC_ENTER(NULL);

    MPI_Barrier(client_cache_comm);

    if (shared_obj_cache_win != MPI_WIN_NULL) {
        MPI_Win_free(&shared_obj_cache_win);
        free(shared_buf);
        free(global_metadata_list);
    }

    pdc_region_cache_timelog(start, "pdc_region_global_metadata_free - total time");

done:
    fflush(stdout);
    FUNC_LEAVE(ret_value);
}

int
pdc_region_dl_global_search(pdcid_t obj_id, uint64_t *offset, uint64_t *size)
{
    perr_t ret_value = SUCCEED;

    int                       i, item_idx, is_cached = 0, is_contained = 1;
    void *                    current_metadata_list;
    pdc_object_list_metadata *current_metadata;
    pdc_object_cache *        current_cache_list;
    char *                    buf;
    double                    start;

    FUNC_ENTER(NULL);

    if (global_metadata_list == NULL) {
        if (pdc_client_mpi_rank_g == 0)
            printf("[RANK %d] pdc_region_dl_global_search - global metadata list not created %p\n",
                   pdc_client_mpi_rank_g, (void *)global_metadata_list);
        goto done;
    }

    // MPI_Win_fence(0, shared_obj_cache_win);

    for (i = 0; i < pdc_client_mpi_size_g; i++) {
        current_metadata_list = (char *)global_metadata_list + (i * metadata_size);

        current_metadata = (pdc_object_list_metadata *)current_metadata_list;
        current_cache_list =
            (pdc_object_cache *)((char *)current_metadata_list + sizeof(pdc_object_list_metadata));

        item_idx = current_metadata->head_idx;

        while (item_idx != -1) {
            if (current_cache_list[item_idx].obj_id == obj_id) {
                // if (strcmp(current_cache_list[item_idx].obj_name, obj_name) == 0) {
                if (offset != NULL) {
                    is_contained = detect_region_contained(
                        offset, size, current_cache_list[item_idx].reg_offset,
                        current_cache_list[item_idx].reg_size, current_cache_list[item_idx].reg_ndim);
                }

                if (is_contained) {
                    start = MPI_Wtime();

                    buf = (char *)PDC_malloc(current_cache_list[item_idx].buf_size);

                    MPI_Win_lock(MPI_LOCK_SHARED, i, 0, shared_obj_cache_win);

                    MPI_Get(buf, current_cache_list[item_idx].buf_size, MPI_BYTE, i,
                            current_cache_list[item_idx].buf_offset, current_cache_list[item_idx].buf_size,
                            MPI_BYTE, shared_obj_cache_win);

                    MPI_Win_flush(i, shared_obj_cache_win);
                    MPI_Win_unlock(i, shared_obj_cache_win);

                    // If the region exist insert it to local cache list
                    ret_value = pdc_region_cache_insert(
                        current_cache_list[item_idx].obj_id, current_cache_list[item_idx].reg_ndim,
                        current_cache_list[item_idx].unit, current_cache_list[item_idx].reg_offset,
                        current_cache_list[item_idx].reg_size, buf);

                    free(buf);

                    is_cached = 1;

                    pdc_region_cache_timelog(start, "pdc_region_dl_search - global cache hit time");

                    break;
                }
            }

            item_idx = current_cache_list[item_idx].next;
        }

        if (is_cached)
            break;
    }

    // MPI_Win_fence(0, shared_obj_cache_win);

done:
    fflush(stdout);
    FUNC_LEAVE(is_cached);
}

item_delete_info
pdc_region_dl_update(pdcid_t obj_id, int ndim, uint64_t unit, uint64_t *offset, uint64_t *size, void *buf)
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
        if (obj_cache_list[item_idx].obj_id == obj_id) {
            // if (strcmp(obj_cache_list[item_idx].obj_name, obj_name) == 0) {
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

    ret_value = pdc_region_global_metadata_free();

    free(metadata_base_ptr);
    free(region_buf);

    MPI_Comm_free(&client_cache_comm);

    init_object_cache = 0;

    pdc_region_cache_timelog(start, "pdc_region_dl_finalize - finalization time");

done:
    fflush(stdout);
    FUNC_LEAVE(ret_value);
}