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

/*typedef struct pdc_object_cache {
    // PDC Object information
    // pdcid_t  obj_id;
    char     obj_name[MAX_NAME_LEN];
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
*/

// Linked list version structure
typedef struct pdc_object_cache {
    // PDC Object information
    char     obj_name[MAX_NAME_LEN];
    uint64_t unit;

    // Remote region information
    int      reg_ndim;
    uint64_t reg_offset[3];
    uint64_t reg_size[3];

    // Region data information
    uint64_t buf_size;
    char *   buf;

    // Data exchange
    int target_rank;

    // Double linked list for cached object list
    struct pdc_object_cache *prev;
    struct pdc_object_cache *next;
} pdc_object_cache;

typedef struct item_delete_info {
    size_t deleted_size;
    int    deleted_item_num;
} item_delete_info;

/**************************************************************************/
/* Private Functions for Linked List Structure Client-side Region Caching */
/**************************************************************************/

perr_t pdc_region_dl_init();

int pdc_region_dl_search(char *obj_name, int ndim, uint64_t unit, uint64_t *offset, uint64_t *size, void *buf,
                         uint64_t read_size);

perr_t pdc_region_dl_prepare_data_exchange(char **global_prefetch_list, int *global_list_len,
                                           int *global_list_item_len);

int pdc_region_dl_global_data_exchange(int recv_num);

perr_t pdc_region_dl_insert(char *obj_name, int ndim, uint64_t unit, uint64_t *offset, uint64_t *size,
                            void *buf, uint64_t read_size);

item_delete_info pdc_region_dl_update(char *obj_name, int ndim, uint64_t unit, uint64_t *offset,
                                      uint64_t *size, void *buf);

item_delete_info pdc_region_dl_evict(size_t required_size);

perr_t pdc_region_dl_finalize();

#endif /* PDC_REGION_CACHE_DL_H */
