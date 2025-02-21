#include <time.h>
#include <stdlib.h>
#include <unistd.h>
#include "pdc_utlist.h"
#include "pdc_config.h"
#include "pdc_id_pkg.h"
#include "pdc_obj.h"
#include "pdc_obj_pkg.h"
#include "pdc_malloc.h"
#include "pdc_prop_pkg.h"
#include "pdc_region.h"
#include "pdc_region_pkg.h"
#include "pdc_region_cache.h"
#include "pdc_region_cache_dl.h"
#include "pdc_obj_pkg.h"
#include "pdc_interface.h"
#include "pdc_transforms_pkg.h"
#include "pdc_client_connect.h"
#include "pdc_analysis_pkg.h"
#include <mpi.h>

// Temporary defined
#define MAX_CACHE_SIZE 34359738368

// Object cache list that will be used for object cache management
static struct pdc_object_cache *obj_cache_list, *obj_cache_list_end;

perr_t
pdc_region_dl_init()
{
    perr_t ret_value = SUCCEED;

    FUNC_ENTER(NULL);

    obj_cache_list     = NULL;
    obj_cache_list_end = NULL;

done:
    fflush(stdout);
    FUNC_LEAVE(ret_value);
}

// TODO:
// Implement when the region is overlapping - partially containing
// Need to manage the object's offset cache information
// Currently considering fully contained case
int
pdc_region_dl_search(pdcid_t obj_id, int ndim, uint64_t unit, uint64_t *offset, uint64_t *size, void *buf)
{
    perr_t ret_value = SUCCEED;

    struct pdc_object_cache *obj_cache_iter;
    struct pdc_region_cache *reg_cache_iter;
    uint64_t *               overlap_offset, *overlap_size;
    int                      region_contained = 0;

    obj_cache_iter = obj_cache_list;

    // Navigate through the object list
    while (obj_cache_iter != NULL) {
        if (obj_cache_iter->obj_id == obj_id) {
            reg_cache_iter = obj_cache_iter->reg_cache_list;

            // Compare offset and offset + size and see if the region is contained
            // If region cache exists, memcpy the data into the buf
            while (reg_cache_iter != NULL) {
                // Check if the region is fully contained within the region list
                region_contained =
                    detect_region_contained(offset, size, reg_cache_iter->reg_offset,
                                            reg_cache_iter->reg_size, reg_cache_iter->reg_ndim);

                // Get the offset and size information of overlapped region part
                PDC_region_overlap_detect(ndim, offset, size, reg_cache_iter->reg_offset,
                                          reg_cache_iter->reg_size, &overlap_offset, &overlap_size);

                // Copy the overlapped part into the provided transfer_request buffer
                memcpy_overlap_subregion(reg_cache_iter->reg_ndim, unit, reg_cache_iter->buf,
                                         reg_cache_iter->reg_offset, reg_cache_iter->reg_size, buf, offset,
                                         size, overlap_offset, overlap_size);

                // Update region cache list according to LRU policy
                if (region_contained) {
                    ret_value = pdc_region_dl_LRU(reg_cache_iter, obj_cache_iter->reg_cache_list,
                                                  obj_cache_iter->reg_cache_list_end, 2);
                    free(overlap_offset);
                    break;
                }
                reg_cache_iter = reg_cache_iter->next;
            }

            // Update object cache list according to LRU policy
            if (region_contained) {
                ret_value = pdc_region_dl_LRU(obj_cache_iter, obj_cache_list, obj_cache_list_end, 1);
                break;
            }
        }
        obj_cache_iter = obj_cache_iter->next;
    }

    return region_contained;
}

// Insert the region to the list
perr_t
pdc_region_dl_insert(pdcid_t obj_id, int ndim, uint64_t unit, uint64_t *offset, uint64_t *size, void *buf,
                     uint64_t read_size)
{
    perr_t ret_value = SUCCEED;

    struct pdc_object_cache *obj_cache_iter, *obj_cache_item = NULL;
    struct pdc_region_cache *reg_cache_item;

    int obj_list_end_init = 0;

    double reg_list_start, reg_memcpy;

    FUNC_ENTER(NULL);

    obj_cache_iter = obj_cache_list;
    if (obj_cache_list == NULL) {
        obj_list_end_init = 1;
    }

    // Searching for obj_id in the object cache list
    while (obj_cache_iter != NULL) {
        if (obj_cache_iter->obj_id == obj_id) {
            obj_cache_item = obj_cache_iter;
            break;
        }
        obj_cache_iter = obj_cache_iter->next;
    }

    // If the obj_id does not exist add object item to the object cache list
    if (obj_cache_item == NULL) {
        obj_cache_item = (struct pdc_object_cache *)PDC_malloc(sizeof(struct pdc_object_cache));
        if (!obj_cache_item)
            PGOTO_ERROR(FAIL, "PDC region cache - obj_cache_item memory allocation failed");

        obj_cache_item->obj_id         = obj_id;
        obj_cache_item->reg_cache_list = NULL;

        if (obj_cache_list == NULL) {
            obj_cache_list_end = obj_cache_list;
        }

        DL_PREPEND(obj_cache_list, obj_cache_item);
        if (obj_list_end_init) {
            obj_cache_list_end = obj_cache_list;
        }
    }

    // Check if there are overlapping parts of region lists

    // Insert the region to the list
    // Check if the region cache list exists for the obj_id
    // If it does not exists create the list and insert the region
    reg_list_start = MPI_Wtime();

    reg_cache_item = (struct pdc_region_cache *)PDC_malloc(sizeof(struct pdc_region_cache));
    if (!reg_cache_item)
        PGOTO_ERROR(FAIL, "PDC region cache - reg_cache_item memory allocation failed");

    reg_cache_item->reg_ndim = ndim;

    // memcpy offset and size continuously
    reg_cache_item->reg_offset = (uint64_t *)PDC_malloc(sizeof(uint64_t) * ndim * 2);
    if (!reg_cache_item->reg_offset)
        PGOTO_ERROR(FAIL, "PDC region cache - reg_cache_item->reg_offset memory allocation failed");

    reg_cache_item->reg_size = reg_cache_item->reg_offset + ndim;

    reg_cache_item->buf = (char *)PDC_malloc(sizeof(char) * read_size);
    if (!reg_cache_item->buf)
        PGOTO_ERROR(FAIL, "PDC region cache - reg_cache_item->buf memory allocation failed");

    reg_memcpy = MPI_Wtime();

    memcpy(reg_cache_item->reg_offset, offset, sizeof(uint64_t) * ndim);
    memcpy(reg_cache_item->reg_size, size, sizeof(uint64_t) * ndim);
    memcpy(reg_cache_item->buf, buf, sizeof(char) * read_size);

    pdc_region_cache_timelog(5, reg_memcpy, "pdc_region_cache_insert reg_memcpy time");

    if (obj_cache_item->reg_cache_list == NULL) {
        DL_PREPEND(obj_cache_item->reg_cache_list, reg_cache_item);
        obj_cache_item->reg_cache_list_end = obj_cache_item->reg_cache_list;
    }
    else {
        DL_PREPEND(obj_cache_item->reg_cache_list, reg_cache_item);
    }

    pdc_region_cache_timelog(5, reg_list_start, "pdc_region_cache_insert reg_list_insert time");

done:
    fflush(stdout);
    FUNC_LEAVE(ret_value);
}

int
pdc_region_dl_update(pdcid_t obj_id, int ndim, uint64_t unit, uint64_t *offset, uint64_t *size, void *buf)
{

    struct pdc_object_cache *obj_cache_iter, *obj_cache_item;
    struct pdc_region_cache *reg_cache_iter, *reg_cache_item;
    uint64_t *               overlap_offset, *overlap_size;
    int                      is_overlapped     = 0;
    int                      one_item_reg_list = 0, one_item_obj_list = 0;

    obj_cache_iter = obj_cache_list;

    // Navigate through the object list
    while (obj_cache_iter != NULL) {
        if (obj_cache_iter->obj_id == obj_id) {
            reg_cache_iter = obj_cache_iter->reg_cache_list;

            while (reg_cache_iter != NULL) {
                // Compare offset and offset + size and see if there is an overlap
                is_overlapped =
                    check_overlap(ndim, offset, size, reg_cache_iter->reg_offset, reg_cache_iter->reg_size);

                reg_cache_item = reg_cache_iter;
                reg_cache_iter = reg_cache_iter->next;

                if (is_overlapped) {
                    // If the target region is the end region, update the reg_cache_list_end
                    if (reg_cache_item == obj_cache_iter->reg_cache_list_end) {
                        if (reg_cache_item->prev) {
                            obj_cache_iter->reg_cache_list_end = reg_cache_item->prev;
                        }
                        else {
                            one_item_reg_list = 1;
                        }
                    }

                    // Delete the last item of the list and free the buffer
                    DL_DELETE(obj_cache_iter->reg_cache_list, reg_cache_item);

                    free(reg_cache_item->reg_offset);
                    free(reg_cache_item->buf);
                    free(reg_cache_item);
                }
            }
        }

        obj_cache_item = obj_cache_iter;
        obj_cache_iter = obj_cache_iter->next;

        // If the target object had one region but deleted, delete the object item as well
        if (one_item_reg_list) {
            // If the target object is the end object, update the obj_cache_list_end
            if (obj_cache_item == obj_cache_list_end) {
                if (obj_cache_list_end->prev) {
                    obj_cache_list_end = obj_cache_list_end->prev;
                }
            }
            // Delete the empty object item
            DL_DELETE(obj_cache_list, obj_cache_item);
            free(obj_cache_item);
        }
    }

    return is_overlapped;
}

size_t
pdc_region_dl_evict(size_t required_size)
{
    struct pdc_object_cache *obj_cache_iter, *obj_cache_item;
    struct pdc_region_cache *reg_cache_item;
    size_t                   deleted_size = 0;

    obj_cache_iter = obj_cache_list_end;

    // From the end of the object list, free the object item of the list until it matches the required size
    while (obj_cache_iter != NULL) {
        obj_cache_item = obj_cache_iter;

        // From the end of the region list, free the region item of the list until it matches the required
        // size
        do {
            reg_cache_item                     = obj_cache_iter->reg_cache_list_end;
            obj_cache_iter->reg_cache_list_end = obj_cache_iter->reg_cache_list_end->prev;

            required_size -= sizeof(reg_cache_item->buf);
            deleted_size += sizeof(reg_cache_item->buf);

            // Delete the last item of the list and free the buffer
            DL_DELETE(obj_cache_iter->reg_cache_list, reg_cache_item);

            free(reg_cache_item->reg_offset);
            free(reg_cache_item->buf);
            free(reg_cache_item);

            if (required_size < MAX_CACHE_SIZE) {
                break;
            }
        } while (reg_cache_item != NULL);

        if (required_size < MAX_CACHE_SIZE) {
            break;
        }

        obj_cache_iter     = obj_cache_iter->prev;
        obj_cache_list_end = obj_cache_iter;

        // Delete the empty object item
        DL_DELETE(obj_cache_list, obj_cache_item);
        free(obj_cache_item);
    }

    return deleted_size;
}

// Update object cache list according to LRU policy
perr_t
pdc_region_dl_LRU(void *target_item, void *target_list, void *target_list_end, int type)
{
    perr_t                   ret_value = SUCCEED;
    struct pdc_object_cache *target_obj;
    struct pdc_region_cache *target_reg, *target_reg_list, *target_reg_list_end;

    // Indication for one item list
    int one_item_list = 0;

    FUNC_ENTER(NULL);

    // Object list update according to LRU policy
    if (type == 1) {
        target_obj = (pdc_object_cache *)target_item;

        // Update the obj_cache_list_end information
        if (target_obj == obj_cache_list_end) {
            if (obj_cache_list_end->prev) {
                obj_cache_list_end = obj_cache_list_end->prev;
            }
            else {
                one_item_list = 1;
            }
        }

        // Move the recently searched object to the front of the list
        if (!one_item_list) {
            DL_DELETE(obj_cache_list, target_obj);
            DL_PREPEND(obj_cache_list, target_obj);
        }
    }

    // Region list update according to LRU policy
    if (type == 2) {
        target_reg          = (pdc_region_cache *)target_item;
        target_reg_list     = (pdc_region_cache *)target_list;
        target_reg_list_end = (pdc_region_cache *)target_list_end;

        // Move the recently searched region into the front of the list
        if (target_reg == target_reg_list_end) {
            if (target_reg_list_end->prev) {
                target_reg_list_end = target_reg_list_end->prev;
            }
            else {
                one_item_list = 1;
            }
        }
        if (!one_item_list) {
            DL_DELETE(target_reg_list, target_reg);
            DL_PREPEND(target_reg_list, target_reg);
        }
    }

done:
    fflush(stdout);
    FUNC_LEAVE(ret_value);
}