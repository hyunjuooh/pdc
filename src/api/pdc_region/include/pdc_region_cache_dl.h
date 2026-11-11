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

#define NUM_CHUNKS 21
// #define NUM_CHUNKS               4

// Overlapped exchange: most items a rank sends to one other rank in one batch
// exchange (one per training step, so at least the batch size; more wait for
// the next one)
#define BATCH_ITEMS_PER_DEST 8

// 1GB data generation for bdcats
// #define MAX_ITEM_SIZE 33554432

// #define MAX_ITEM_SIZE      67108864
// #define MAX_SLOTS_PER_NODE 234

// #define MAX_ITEM_SIZE      134217728

#define MAX_ITEM_SIZE      157286400
#define MAX_SLOTS_PER_NODE 950
// #define MAX_SLOTS_PER_NODE 712
// #define MAX_SLOTS_PER_NODE 475

#define SLOT_INVALID -1

/**************************/
/* Library Private Struct */
/**************************/
typedef struct {
    int free_stack_head;
    // Server reads still to come on this node in the current epoch while an
    // overlapped exchange runs; next-epoch receives must leave this many slots free.
    int  pending_inserts;
    int  free_list[MAX_SLOTS_PER_NODE];
    char padding[64 - sizeof(int)];
} SharedMemoryHeader;

typedef struct pdc_object_data {
    // PDC Object information
    pdcid_t  obj_id;
    uint64_t unit;

    // Remote region information
    int reg_ndim;

    uint64_t reg_offset[3];
    uint64_t reg_size[3];
    uint64_t reg_buf_size;

    // Region data
    int slot_idx;
    // char reg_buf[MAX_ITEM_SIZE];

    // Index info
    int target_rank;
    int data_exchange_type; // 0 is intra, 1 is inter

    // Overlapped data exchange state
    int consumed;       // read in the current epoch (server read or cache hit)
    int transfer_state; // 0 idle, 1 being sent, 2 received for the next epoch

    struct pdc_object_data *prev;
    struct pdc_object_data *next;
} pdc_object_data;

// Specifies client specific information
typedef struct pdc_client_info {
    // Global and local rank info
    int world_rank;
    int world_size;

    int node_rank;
    int node_size;
    int node_manager_rank;

    int cached_item_num;

    pdc_object_data *local_cache_list_head;
    pdc_object_data *local_cache_list_tail;

    MPI_Win             node_shared_data_win;
    SharedMemoryHeader *header;
    void *              node_shared_base;
    char *              node_shared_data_base;

    int *rank_to_node_id_map;
    int *world_to_node_rank_map;
    int *target_ranks;

    int client_cache_init;
} pdc_client_info;

/**************************************************************************/
/* Private Functions for Linked List Structure Client-side Region Caching */
/**************************************************************************/

perr_t pdc_region_dl_init();

int pdc_region_dl_local_search(pdcid_t obj_id, int ndim, uint64_t unit, uint64_t *offset, uint64_t *size,
                               void *buf, uint64_t read_size);

perr_t pdc_region_dl_insert(pdcid_t obj_id, int ndim, uint64_t unit, uint64_t *offset, uint64_t *size,
                            void *buf, uint64_t read_size);

perr_t pdc_region_dl_update(pdcid_t obj_id, int ndim, uint64_t unit, uint64_t *offset, uint64_t *size,
                            void *buf);

perr_t pdc_region_dl_evict();

perr_t pdc_region_dl_finalize();

// Cache list and the node's slot pool, shared with the data exchange
extern pdc_client_info client_info;
extern MPI_Comm        client_cache_world_comm;
extern MPI_Comm        client_cache_node_comm;

pdc_object_data *pdc_region_dl_new_item(pdcid_t obj_id, int ndim, uint64_t unit, uint64_t *offset,
                                        uint64_t *size, uint64_t reg_buf_size, int slot_idx);
perr_t           pdc_region_dl_prepend(pdc_object_data *obj_cache_item);
perr_t           pdc_region_dl_delete(pdc_object_data *obj_cache_item);
char *           get_data_ptr(int shared_index);
int              pop_free_slot();
void             push_free_slot(int slot_idx);
int              count_free_slots();

/**************************************************************************/
/* Data exchange between epochs (pdc_region_dl_exchange.c)                */
/**************************************************************************/

perr_t pdc_region_dl_prepare_data_exchange(pdcid_t *global_prefetch_list, uint64_t *offset, uint64_t *size,
                                           int obj_prefetch_list_len);

perr_t pdc_region_dl_data_exchange(pdcid_t *global_prefetch_list, int obj_prefetch_list_len);

// Exchange that overlaps with computation: begin at the start of an epoch with
// this epoch's and next epoch's object lists (metadata ids, in read order),
// progress after every step, end before the next epoch starts.
perr_t pdc_region_dl_overlap_begin(pdcid_t *cur_ids, int cur_len, pdcid_t *next_ids, int next_len);
perr_t pdc_region_dl_overlap_progress(void);
perr_t pdc_region_dl_overlap_end(void);

// Used by pdc_region_dl_insert while an overlapped exchange runs
int  pdc_region_dl_overlap_active(void);
void pdc_region_dl_overlap_make_room(void);
void pdc_region_dl_overlap_record_insert(void);

#endif /* PDC_REGION_CACHE_DL_H */
