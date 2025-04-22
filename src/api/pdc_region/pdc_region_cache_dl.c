#include <time.h>
#include <stdlib.h>
#include <unistd.h>
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
#include "pdc_client_connect.h"

#define MAX_CACHE_SIZE 4294967296
#define MAX_ITEM_NUM   1000

MPI_Comm client_cache_comm;

int init_object_cache       = 0;
int obj_cache_list_item_num = 0; // Tracks the number of objects cached within the list

// Object cache list that will be used for object cache management
struct pdc_object_cache *        obj_cache_list;
struct pdc_object_list_metadata *obj_cache_list_metadata;

char * metadata_base_ptr;
size_t metadata_size;

void *global_metadata_list           = NULL;
int   global_metadata_list_collected = 0;

// For global cache creation
MPI_Win  shared_obj_cache_win;
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
pdc_region_dl_prepend(int item_idx)
{
    perr_t ret_value = SUCCEED;

    FUNC_ENTER(NULL);

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

    int      mpi_alloc_error, i;
    MPI_Aint buffer_win_size;

    double start;

    FUNC_ENTER(NULL);

    obj_cache_list_item_num = 0;

    // Duplicate MPI_COMM_WORLD for client cache
    MPI_Comm_dup(MPI_COMM_WORLD, &client_cache_comm);

    // Shared memory buffer allocation
    start             = MPI_Wtime();
    region_buf_offset = 0;
    buffer_win_size   = MAX_CACHE_SIZE * sizeof(char);

    shared_buf = (char *)PDC_malloc(MAX_CACHE_SIZE * sizeof(char));
    region_buf = (char *)PDC_malloc(MAX_CACHE_SIZE * sizeof(char));

    mpi_alloc_error = MPI_Win_create(shared_buf, buffer_win_size, sizeof(char), MPI_INFO_NULL,
                                     client_cache_comm, &shared_obj_cache_win);

    if (mpi_alloc_error != MPI_SUCCESS)
        PGOTO_ERROR(FAIL, "PDC region dl init - MPI shared allocation for shared_obj_cache_win failed");

    pdc_region_cache_timelog(start, "pdc_region_dl_init - Init shared region memory buffer");

    // Shared memory metadata allocation
    start = MPI_Wtime();

    metadata_size     = sizeof(pdc_object_list_metadata) + MAX_ITEM_NUM * sizeof(pdc_object_cache);
    metadata_base_ptr = (char *)PDC_malloc(metadata_size);

    pdc_region_cache_timelog(start, "pdc_region_dl_init - Init metadata memory buffer");

    // Continuous memory
    obj_cache_list_metadata = (pdc_object_list_metadata *)metadata_base_ptr;
    obj_cache_list = (pdc_object_cache *)((char *)metadata_base_ptr + sizeof(pdc_object_list_metadata));

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

    if (!init_object_cache) {
        // ret_value = pdc_region_dl_init();
        PGOTO_ERROR(FAIL, "pdc_region_dl_insert - object cache list not initialized");
    }

    new_item_idx = obj_cache_new_item_index();
    if (new_item_idx == -1)
        PGOTO_ERROR(FAIL, "pdc_region_dl_insert - new item index not found");

    // Create obj_cache_list item
    obj_cache_list[new_item_idx].obj_id   = obj_id;
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

perr_t
pdc_region_dl_collect_global_metadata()
{
    perr_t ret_value = SUCCEED;

    int    mpi_alloc_error;
    double start;

    FUNC_ENTER(NULL);

    // if (global_metadata_list != NULL) {
    //     free(global_metadata_list);
    //     global_metadata_list = NULL;
    // }

    start = MPI_Wtime();
    memcpy(shared_buf, region_buf, MAX_CACHE_SIZE);
    pdc_region_cache_timelog(start, "pdc_region_dl_collect_global_metadata - memcpy time");

    // Allgather the metadata information
    mpi_alloc_error = MPI_Allgather(metadata_base_ptr, metadata_size, MPI_BYTE, global_metadata_list,
                                    metadata_size, MPI_BYTE, client_cache_comm);
    if (mpi_alloc_error != MPI_SUCCESS)
        PGOTO_ERROR(FAIL, "pdc_region_dl_search - MPI Allgather failed");

    global_metadata_list_collected = 1;

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

    double start = MPI_Wtime();

    FUNC_ENTER(NULL);

    if (!init_object_cache) {
        return is_cached;
    }

    if (!obj_cache_list_item_num) {
        return is_cached;
    }

    // Find if region is cached into local object list cache
    // If region contained, return the rank that contains the region
    item_idx = obj_cache_list_metadata->head_idx;

    while (item_idx != -1) {
        if (obj_cache_list[item_idx].obj_id == obj_id) {
            is_cached = detect_region_contained(offset, size, obj_cache_list[item_idx].reg_offset,
                                                obj_cache_list[item_idx].reg_size, ndim);

            if (is_cached) {
                // If region contained, memcpy cached region data into transfer_request->buf
                // Detect the offset range that is overlapped
                start = MPI_Wtime();
                PDC_region_overlap_detect(ndim, offset, size, obj_cache_list[item_idx].reg_offset,
                                          obj_cache_list[item_idx].reg_size, &overlap_offset, &overlap_size);

                // memcpy the overlapped region
                memcpy_overlap_subregion(
                    obj_cache_list[item_idx].reg_ndim, unit, region_buf + obj_cache_list[item_idx].buf_offset,
                    obj_cache_list[item_idx].reg_offset, obj_cache_list[item_idx].reg_size, buf, offset, size,
                    overlap_offset, overlap_size);

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

int
pdc_region_dl_global_search(pdcid_t obj_id)
{
    perr_t ret_value = SUCCEED;

    int                       i, item_idx, is_cached = 0;
    void *                    current_metadata_list;
    pdc_object_list_metadata *current_metadata;
    pdc_object_cache *        current_cache_list;
    char *                    buf;

    FUNC_ENTER(NULL);

    if (!global_metadata_list_collected) {
        if (pdc_client_mpi_rank_g == 0)
            printf("[RANK %d] pdc_region_dl_search - global metadata list not created %p\n",
                   pdc_client_mpi_rank_g, (void *)global_metadata_list);
        goto done;
    }

    for (i = 0; i < pdc_client_mpi_size_g; i++) {
        if (i == pdc_client_mpi_rank_g)
            continue;

        current_metadata_list = (char *)global_metadata_list + (i * metadata_size);

        current_metadata = (pdc_object_list_metadata *)current_metadata_list;
        current_cache_list =
            (pdc_object_cache *)((char *)current_metadata_list + sizeof(pdc_object_list_metadata));

        while (item_idx != -1) {
            if (current_cache_list[item_idx].obj_id == obj_id) {
                buf = (char *)PDC_malloc(current_cache_list[item_idx].buf_size);

                MPI_Win_lock(MPI_LOCK_SHARED, i, 0, shared_obj_cache_win);

                MPI_Get(buf, current_cache_list[item_idx].buf_size, MPI_BYTE, i,
                        current_cache_list[item_idx].buf_offset, current_cache_list[item_idx].buf_size,
                        MPI_BYTE, shared_obj_cache_win);

                MPI_Win_flush(i, shared_obj_cache_win);
                MPI_Win_unlock(i, shared_obj_cache_win);

                ret_value = pdc_region_cache_insert(
                    current_cache_list[item_idx].obj_id, current_cache_list[item_idx].reg_ndim,
                    current_cache_list[item_idx].unit, current_cache_list[item_idx].reg_offset,
                    current_cache_list[item_idx].reg_size, buf);

                is_cached = 1;

                free(buf);
                break;
            }

            item_idx = current_cache_list[item_idx].next;
        }
    }
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

    // Delete item from end of list (following LRU policy)
    while (item_idx != -1) {
        item_idx = obj_cache_list_metadata->tail_idx;

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

    ret_value = pdc_region_dl_delete(item_idx);

    // Update the region_buf_offset for global data buffer
    // Applies only when evicting the last item
    region_buf_offset -= obj_cache_list[item_idx].buf_size;

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

    MPI_Win_free(&shared_obj_cache_win);

    free(metadata_base_ptr);
    free(region_buf);
    free(shared_buf);

    if (global_metadata_list != NULL)
        free(global_metadata_list);

    MPI_Comm_free(&client_cache_comm);

    init_object_cache = 0;

    pdc_region_cache_timelog(start, "pdc_region_dl_finalize - finalization time");

done:
    fflush(stdout);
    FUNC_LEAVE(ret_value);
}