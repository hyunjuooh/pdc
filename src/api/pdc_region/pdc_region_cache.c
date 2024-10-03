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
#include "pdc_obj_pkg.h"
#include "pdc_interface.h"
#include "pdc_transforms_pkg.h"
#include "pdc_client_connect.h"
#include "pdc_analysis_pkg.h"
#include <mpi.h>

// Temporary defined variable
#define MAX_CACHE_SIZE 34359738368

static size_t                   total_buf_size;
static struct pdc_object_cache *obj_cache_list, *obj_cache_list_end;

// Initialization of global variables
perr_t
pdc_region_cache_init()
{
    perr_t ret_value = SUCCEED;

    FUNC_ENTER(NULL);

    total_buf_size = 0;
    obj_cache_list = NULL;

done:
    fflush(stdout);
    FUNC_LEAVE(ret_value);
}

// TODO:
// Implement when the region is overlapping - partially containing
// Need to manage the object's offset cache information
// Currently considering fully contained case
int
pdc_region_cache_search(pdcid_t obj_id, int ndim, size_t unit, uint64_t *offset, uint64_t *size, void *buf)
{
    // perr_t                   ret_value = SUCCEED;
    struct pdc_object_cache *obj_cache_iter;
    struct pdc_region_cache *reg_cache_iter;
    uint64_t *               overlap_offset, *overlap_size;
    int                      region_contained = 0;

    obj_cache_iter = obj_cache_list;

    // FUNC_ENTER(NULL);

    // Navigate through the object list
    while (obj_cache_iter != NULL) {
        if (obj_cache_iter->obj_id == obj_id) {
            reg_cache_iter = obj_cache_iter->reg_cache_list;

            // Compare offset and offset + size and see if the region is contained
            while (reg_cache_iter != NULL) {
                // Check if the region is fully contained within the region list
                region_contained =
                    detect_region_contained(offset, size, reg_cache_iter->reg_offset,
                                            reg_cache_iter->reg_size, reg_cache_iter->reg_ndim);

                // Currently considering fully contained case
                if (region_contained) {
                    // Get the offset and size information of overlapped region part
                    PDC_region_overlap_detect(ndim, offset, size, reg_cache_iter->reg_offset,
                                              reg_cache_iter->reg_size, &overlap_offset, &overlap_size);

                    // Copy the overlapped part into the provided transfer_request buffer
                    memcpy_overlap_subregion(reg_cache_iter->reg_ndim, unit, reg_cache_iter->buf,
                                             reg_cache_iter->reg_offset, reg_cache_iter->reg_size, buf,
                                             offset, size, overlap_offset, overlap_size);

                    // Move the recently searched region into the front of the list
                    DL_DELETE(obj_cache_iter->reg_cache_list, reg_cache_iter);
                    DL_PREPEND(obj_cache_iter->reg_cache_list, reg_cache_iter);

                    free(overlap_offset);
                    free(overlap_size);

                    break;
                }

                reg_cache_iter = reg_cache_iter->next;
            }

            // Update the obj_cache_list_end information
            if (obj_cache_iter == obj_cache_list_end) {
                obj_cache_list_end = obj_cache_list_end->prev;
            }

            // Move the recently searched object to the front of the list
            DL_DELETE(obj_cache_list, obj_cache_iter);
            DL_PREPEND(obj_cache_list, obj_cache_iter);

            break;
        }

        obj_cache_iter = obj_cache_iter->next;
    }

    return region_contained;

// done:
//     fflush(stdout);
//     FUNC_LEAVE(ret_value);
}

// Insert the region to the list
perr_t
pdc_region_cache_insert(pdcid_t obj_id, int ndim, uint64_t *offset, uint64_t *size, void *buf)
{
    perr_t ret_value = SUCCEED;

    struct pdc_object_cache *obj_cache_iter, *obj_cache_item = NULL;
    struct pdc_region_cache *reg_cache_item;

    FUNC_ENTER(NULL);

    // Check if there is remaining buffer size to insert region
    // If there is no remaining capacity, free the buffer according to LRU policy
    if (total_buf_size + sizeof(buf) > MAX_CACHE_SIZE) {
        pdc_region_cache_evict(total_buf_size + sizeof(buf));
    }

    // If there is no object list, generate the list and insert the item
    if (obj_cache_list == NULL) {
        obj_cache_item = (struct pdc_object_cache *)PDC_malloc(sizeof(struct pdc_object_cache));
        if (!obj_cache_item)
            PGOTO_ERROR(0, "PDC region cache - obj_cache_item memory allocation failed");

        obj_cache_item->obj_id         = obj_id;
        obj_cache_item->reg_cache_list = NULL;

        DL_PREPEND(obj_cache_list, obj_cache_item);
        obj_cache_list_end = obj_cache_list;
    }
    else {
        // Check if the specific obj_id exists in the object cache list
        obj_cache_iter = obj_cache_list;

        while (obj_cache_iter != NULL) {
            if (obj_cache_iter->obj_id == obj_id) {
                obj_cache_item = obj_cache_iter;
                break;
            }
            obj_cache_iter = obj_cache_iter->next;
        }

        // If it does not exists create the list prior to insertion
        if (obj_cache_item == NULL) {
            obj_cache_item = (struct pdc_object_cache *)PDC_malloc(sizeof(struct pdc_object_cache));
            if (!obj_cache_item)
                PGOTO_ERROR(0, "PDC region cache - obj_cache_item memory allocation failed");

            obj_cache_item->obj_id         = obj_id;
            obj_cache_item->reg_cache_list = NULL;

            DL_PREPEND(obj_cache_list, obj_cache_item);
        }
    }

    // Check if there are overlapping parts of region lists

    // Insert the region to the list
    // Check if the region cache list exists for the obj_id
    // If it does not exists create the list and insert the region

    // pdc_malloc
    // check if allocation is successful or not
    reg_cache_item = (struct pdc_region_cache *)PDC_malloc(sizeof(struct pdc_region_cache));
    if (!reg_cache_item)
        PGOTO_ERROR(0, "PDC region cache - reg_cache_item memory allocation failed");

    reg_cache_item->reg_ndim = ndim;

    // memcpy offset and size continuously
    reg_cache_item->reg_offset = (uint64_t *)PDC_malloc(sizeof(uint64_t) * ndim * 2);
    if (!reg_cache_item->reg_offset)
        PGOTO_ERROR(0, "PDC region cache - reg_cache_item->reg_offset memory allocation failed");

    reg_cache_item->reg_size = reg_cache_item->reg_offset + ndim;

    reg_cache_item->buf = (char *)PDC_malloc(sizeof(char) * sizeof(buf));
    if (!reg_cache_item->reg_offset)
        PGOTO_ERROR(0, "PDC region cache - reg_cache_item->reg_size memory allocation failed");

    memcpy(reg_cache_item->reg_offset, offset, sizeof(uint64_t) * ndim);
    memcpy(reg_cache_item->reg_size, size, sizeof(uint64_t) * ndim);
    memcpy(reg_cache_item->buf, buf, sizeof(buf));

    if (obj_cache_item->reg_cache_list == NULL) {
        DL_PREPEND(obj_cache_item->reg_cache_list, reg_cache_item);
        obj_cache_item->reg_cache_list_end = obj_cache_item->reg_cache_list;
    }
    else {
        DL_PREPEND(obj_cache_item->reg_cache_list, reg_cache_item);
    }

    total_buf_size += sizeof(buf);

done:
    fflush(stdout);
    FUNC_LEAVE(ret_value);
}

// // Check the overlapping part of the regions
// perr_t pdc_region_cache_overlap(pdcid_t obj_id, int ndim, uint64_t *offset, uint64_t *size, void *buf){

// }

// Evict the region cache and object cache according to LRU policy
perr_t
pdc_region_cache_evict(size_t required_size)
{
    perr_t                   ret_value = SUCCEED;
    struct pdc_object_cache *obj_cache_iter, *obj_cache_item;
    struct pdc_region_cache *reg_cache_item;

    FUNC_ENTER(NULL);

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

            // Delete the last item of the list and free the buffer
            DL_DELETE(obj_cache_iter->reg_cache_list, reg_cache_item);

            free(reg_cache_item->reg_offset);
            free(reg_cache_item->reg_size);
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

done:
    fflush(stdout);
    FUNC_LEAVE(ret_value);
}
