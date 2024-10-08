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

#ifndef PDC_REGION_CACHE_H
#define PDC_REGION_CACHE_H

#include "pdc_public.h"
#include "pdc_obj.h"

/**************************/
/* Library Private Struct */
/**************************/

typedef struct pdc_object_cache {
    // PDC Object information
    pdcid_t obj_id;
    // int         obj_ndim;

    // Cached region list for this object
    struct pdc_region_cache *reg_cache_list, *reg_cache_list_end;

    // Double linked list for cached object list
    struct pdc_object_cache *prev;
    struct pdc_object_cache *next;
} pdc_object_cache;

typedef struct pdc_region_cache {
    // Region information(remote region)
    int       reg_ndim;
    uint64_t *reg_offset;
    uint64_t *reg_size;

    // Region Buffer
    char *buf;

    struct pdc_region_cache *prev;
    struct pdc_region_cache *next;
} pdc_region_cache;

/****************************************************/
/* Private Functions for Client-side Region Caching */
/****************************************************/

perr_t pdc_region_cache_init();

int pdc_region_cache_search(pdcid_t obj_id, int ndim, uint64_t unit, uint64_t *offset, uint64_t *size,
                            void *buf);

perr_t pdc_region_cache_insert(pdcid_t obj_id, int ndim, uint64_t unit, uint64_t *offset, uint64_t *size, void *buf);

perr_t pdc_region_cache_evict(size_t required_size);

#endif /* PDC_REGION_CACHE_H */