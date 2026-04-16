/*
 * Copyright Notice for
 * Proactive Data Containers (PDC) Software Library and Utilities
 * -----------------------------------------------------------------------------

 *** Copyright Notice ***

 * Proactive Data Containers (PDC) Copyright (c) 2017, The Regents of the
 * University of California, through Lawrence Berkeley National Laboratory,
 * UChicago Argonne, LLC, operator of Argonne National Laboratory, and The HDF
 * Group (subject to receipt of any required approvals from the U.S. Dept. of
 * Energy).  All rights reserved.

 * If you have questions about your rights to use or distribute this software,
 * please contact Berkeley Lab's Innovation & Partnerships Office at  IPO@lbl.gov.

 * NOTICE.  This Software was developed under funding from the U.S. Department of
 * Energy and the U.S. Government consequently retains certain rights. As such, the
 * U.S. Government has been granted for itself and others acting on its behalf a
 * paid-up, nonexclusive, irrevocable, worldwide license in the Software to
 * reproduce, distribute copies to the public, prepare derivative works, and
 * perform publicly and display publicly, and to permit other to do so.
 */

#ifndef PDC_REGION_CACHE_DL_H
#define PDC_REGION_CACHE_DL_H

#include "pdc_public.h"
#include "pdc_obj.h"

#define MAX_NAME_LEN 100

/**************************/
/* Library Private Struct */
/**************************/

typedef struct pdc_object_cache {
    // PDC Object information
    pdcid_t obj_id;
    // char     obj_name[MAX_NAME_LEN];
    uint64_t unit;

    // Remote region information
    int      reg_ndim;
    uint64_t reg_offset[3];
    uint64_t reg_size[3];

    // Region data information
    uint64_t buf_offset;
    uint64_t buf_size;

    // Double linked list for cached object list
    int next;
    int prev;
} pdc_object_cache;

typedef struct pdc_object_list_metadata {
    int head_idx;
    int tail_idx;
    int free_idx;
    int cached_item_num;
} pdc_object_list_metadata;

typedef struct item_delete_info {
    size_t deleted_size;
    int    deleted_item_num;
} item_delete_info;

/**************************************************************************/
/* Private Functions for Linked List Structure Client-side Region Caching */
/**************************************************************************/

perr_t pdc_region_dl_init();

perr_t pdc_region_dl_list_init();

int pdc_region_dl_search(pdcid_t obj_id, int ndim, uint64_t unit, uint64_t *offset, uint64_t *size, void *buf,
                         uint64_t read_size);

perr_t pdc_region_dl_collect_global_metadata();

perr_t pdc_region_global_metadata_free();

int pdc_region_dl_global_search(pdcid_t obj_id, uint64_t *offset, uint64_t *size);

perr_t pdc_region_dl_insert(pdcid_t obj_id, int ndim, uint64_t unit, uint64_t *offset, uint64_t *size,
                            void *buf, uint64_t read_size);

item_delete_info pdc_region_dl_update(pdcid_t obj_id, int ndim, uint64_t unit, uint64_t *offset,
                                      uint64_t *size, void *buf);

item_delete_info pdc_region_dl_evict_by_size(size_t required_size);

item_delete_info pdc_region_dl_evict_by_num();

perr_t pdc_region_dl_clean_list();

perr_t pdc_region_dl_finalize();

#endif /* PDC_REGION_CACHE_DL_H */