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

#include <mpi.h>
#include "pdc_public.h"
#include "pdc_obj.h"

// #define MAX_NAME_LEN 100
#define MAX_ITEM_SIZE  134214281
#define MAX_ITEM_NUM   256

/**************************/
/* Library Private Struct */
/**************************/

typedef struct pdc_object_data {
    // PDC Object information
    pdcid_t  obj_id;
    uint64_t unit;

    // Remote region information
    int      reg_ndim;

    uint64_t reg_offset[3];
    uint64_t reg_size[3];
    uint64_t reg_buf_size;

    // Region data
    char reg_buf[MAX_ITEM_SIZE];
} pdc_object_data;

typedef struct pdc_object_list {
    pdcid_t  obj_id;
    int      slot_index;
    int      target_rank;

    struct pdc_object_list* prev;
    struct pdc_object_list* next;
} pdc_object_list;

// Specifies client specific information
typedef struct pdc_client_info {
    // Global and local rank info
    int world_rank;
    int world_size;

    int node_rank;
    int node_size;
    int node_manager_rank;

    MPI_Win  node_shm_win;
    MPI_Aint node_shm_size;

    pdc_object_data* node_shm_base; // Shared memory base pointer
    pdc_object_data* local_shm_base; // Current client's shared memory base pointer
    pdc_object_data  swap_buffer; // Buffer for data exchange

    pdc_object_list* list_head; // Linked list head
    pdc_object_list* list_tail; // Linked list head
    int list_size;

    int free_slot_indices[MAX_ITEM_NUM];
    int free_slot_count;

    int client_init;
    
    int* rank_to_node_id_map;
} pdc_client_info;

typedef struct item_delete_info {
    size_t deleted_size;
    int    deleted_item_num;
} item_delete_info;

/**************************************************************************/
/* Private Functions for Linked List Structure Client-side Region Caching */
/**************************************************************************/

perr_t pdc_region_dl_init();

int pdc_region_dl_local_search(pdcid_t obj_id, int ndim, uint64_t unit, uint64_t *offset, uint64_t *size, void *buf,
                               uint64_t read_size);

perr_t pdc_region_dl_insert(pdcid_t obj_id, int ndim, uint64_t unit, uint64_t *offset, uint64_t *size, void *buf,
                            uint64_t read_size);

item_delete_info pdc_region_dl_update(pdcid_t obj_id, int ndim, uint64_t unit, uint64_t *offset, uint64_t *size, void *buf);

item_delete_info pdc_region_dl_evict();

perr_t pdc_region_dl_finalize();

#endif /* PDC_REGION_CACHE_DL_H */
