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

// Data exchange between epochs: each cached item moves to the rank that reads
// it next epoch, either blocking between epochs (pdc_region_dl_data_exchange)
// or overlapped with computation (pdc_region_dl_overlap_begin/progress/end).

/**************************************************************************/
/* Shared by both exchanges                                               */
/**************************************************************************/

// A cached item's metadata (pdc_object_data without its data and local list
// state), sent ahead of every item. A same-node item hands over its slot; for
// an item on another node the data follows as a second message.
typedef struct {
    pdcid_t  obj_id;
    uint64_t unit;
    uint64_t reg_offset[3];
    uint64_t reg_size[3];
    uint64_t reg_buf_size;
    int      reg_ndim;
    int      slot_idx;
    int      valid; // overlapped exchange: 0 the sender does not have the item, nothing follows
    int      padding;
} pdc_object_metadata;

static void
metadata_from_item(pdc_object_metadata *metadata, pdc_object_data *item)
{
    memset(metadata, 0, sizeof(pdc_object_metadata));
    metadata->obj_id       = item->obj_id;
    metadata->unit         = item->unit;
    metadata->reg_ndim     = item->reg_ndim;
    metadata->reg_buf_size = item->reg_buf_size;
    metadata->slot_idx     = item->slot_idx;
    metadata->valid        = 1;
    memcpy(metadata->reg_offset, item->reg_offset, sizeof(uint64_t) * 3);
    memcpy(metadata->reg_size, item->reg_size, sizeof(uint64_t) * 3);
}

// A received item, its data in slot_idx, becomes visible in the cache list
static void
add_received_item(pdc_object_metadata *metadata, int slot_idx, int transfer_state)
{
    pdc_object_data *item = pdc_region_dl_new_item(metadata->obj_id, metadata->reg_ndim, metadata->unit,
                                                   metadata->reg_offset, metadata->reg_size,
                                                   metadata->reg_buf_size, slot_idx);
    item->transfer_state = transfer_state;
    pdc_region_dl_prepend(item);
}

// Fills displs from counts and returns the total
static int
counts_to_displs(int n, const int *counts, int *displs)
{
    int total = 0;

    for (int i = 0; i < n; i++) {
        displs[i] = total;
        total += counts[i];
    }
    return total;
}

/**************************************************************************/
/* Blocking exchange between epochs                                       */
/**************************************************************************/

perr_t
pdc_region_dl_prepare_data_exchange(pdcid_t *global_prefetch_list, uint64_t *offset, uint64_t *size,
                                    int obj_prefetch_list_len)
{
    FUNC_ENTER(NULL);

    perr_t ret_value = SUCCEED;

    pdc_object_data *obj_cache_iter;
    int              i, j, prefetch_list_idx = 0, region_compare = 0, is_contained = 0;
    double           start = MPI_Wtime();

    if (offset != NULL) {
        region_compare = 1;
    }

    // Set the target rank for local cached object item
    obj_cache_iter = client_info.local_cache_list_head;

    while (obj_cache_iter != NULL) {
        for (i = 0; i < client_info.world_size; i++) {
            for (j = 0; j < obj_prefetch_list_len; j++) {
                if (obj_cache_iter->obj_id == global_prefetch_list[prefetch_list_idx]) {
                    if (!region_compare) {
                        obj_cache_iter->target_rank = i;
                        break;
                    }
                    else {
                        is_contained =
                            detect_region_contained(&offset[i], &size[i], obj_cache_iter->reg_offset,
                                                    obj_cache_iter->reg_size, obj_cache_iter->reg_ndim);

                        // If the prefetch list is contained within the local object cache list
                        // and if the target_rank is not the same as node_manager_rank
                        if (is_contained) {
                            obj_cache_iter->target_rank = i;
                            break;
                        }
                    }
                }

                prefetch_list_idx++;
            }
            if (obj_cache_iter->target_rank != -1) {
                break;
            }
        }

        is_contained      = 0;
        obj_cache_iter    = obj_cache_iter->next;
        prefetch_list_idx = 0;
    }

    pdc_region_cache_timelog(start, "pdc_region_dl_prepare_data_exchange - global prefetch list preparation");

done:
    fflush(stdout);
    FUNC_LEAVE(ret_value);
}

// Metadata of this round's items of one exchange type (0 same node, 1 other
// node), grouped by destination in a fixed order: dest's group starts at displs[dest].
static void
pack_metadata(pdc_object_data **items, int n_items, int type, int n_dest, const int *displs, int *fill,
             pdc_object_metadata *metadata, pdc_object_data **packed_items)
{
    memset(fill, 0, n_dest * sizeof(int));
    for (int i = 0; i < n_items; i++) {
        if (items[i]->data_exchange_type != type)
            continue;
        int trg  = items[i]->target_rank;
        int dest = (type == 0) ? client_info.world_to_node_rank_map[trg] : trg;
        int k    = displs[dest] + fill[dest]++;
        metadata_from_item(&metadata[k], items[i]);
        if (packed_items != NULL)
            packed_items[k] = items[i];
    }
}

// Exchanges the per-destination counts, then the metadata. Returns the received
// metadata, grouped by sender (sender s's group starts at rdispls[s]).
static pdc_object_metadata *
exchange_metadata(MPI_Comm comm, int comm_size, MPI_Datatype metadata_type, pdc_object_metadata *send_meta,
                 int *send_counts, int *sdispls, int *recv_counts, int *rdispls, int *recv_total)
{
    pdc_object_metadata *recv_meta;

    MPI_Alltoall(send_counts, 1, MPI_INT, recv_counts, 1, MPI_INT, comm);
    *recv_total = counts_to_displs(comm_size, recv_counts, rdispls);

    recv_meta = (pdc_object_metadata *)PDC_calloc(*recv_total + 1, sizeof(pdc_object_metadata));
    MPI_Alltoallv(send_meta, send_counts, sdispls, metadata_type, recv_meta, recv_counts, rdispls,
                  metadata_type, comm);
    return recv_meta;
}

// Receive buffer for inter-node items that find no free slot. Allocated once,
// on first use, and only ever used by one blocking receive at a time.
static char *dl_drop_buf = NULL;

perr_t
pdc_region_dl_data_exchange(pdcid_t *global_prefetch_list, int obj_prefetch_list_len)
{
    FUNC_ENTER(NULL);

    perr_t ret_value = SUCCEED;

    pdc_object_data *exchange_head;
    double           start = MPI_Wtime(), tmp_timer;
    double           intra_time = 0.0, inter_meta_time = 0.0, inter_data_time = 0.0, delete_time = 0.0;
    int              node_size = client_info.node_size, world_size = client_info.world_size;
    int              old_cached_item_num = client_info.cached_item_num;
    int              chunk_size          = old_cached_item_num / NUM_CHUNKS;
    int              n_intra_sent = 0, n_intra_recv = 0, n_inter_sent = 0, n_inter_recv = 0, n_dropped = 0;

    // Items are exchanged in NUM_CHUNKS rounds so each node only needs free
    // slots for one round's incoming items before the round's sent items are freed.
    // Per-destination counts and displacements, all sized world_size (>= node_size).
    int *counts            = (int *)PDC_calloc(9 * (size_t)world_size, sizeof(int));
    int *intra_send_counts = counts, *intra_recv_counts = counts + world_size;
    int *intra_sdispls = counts + 2 * world_size, *intra_rdispls = counts + 3 * world_size;
    int *inter_send_counts = counts + 4 * world_size, *inter_recv_counts = counts + 5 * world_size;
    int *inter_sdispls = counts + 6 * world_size, *inter_rdispls = counts + 7 * world_size;
    int *fill          = counts + 8 * world_size;

    MPI_Datatype metadata_type;
    MPI_Type_contiguous(sizeof(pdc_object_metadata), MPI_BYTE, &metadata_type);
    MPI_Type_commit(&metadata_type);

    exchange_head = client_info.local_cache_list_head;

    for (int c = 0; c < NUM_CHUNKS; c++) {
        int n_c = (c == NUM_CHUNKS - 1) ? old_cached_item_num - c * chunk_size : chunk_size;

        // Collect this round's items (the next n_c old items) and classify them.
        pdc_object_data **items   = (pdc_object_data **)PDC_malloc((n_c + 1) * sizeof(pdc_object_data *));
        int               n_intra = 0;
        int               n_inter = 0;
        for (int i = 0; i < n_c && exchange_head != NULL; i++) {
            items[i]      = exchange_head;
            exchange_head = exchange_head->next;
        }

        memset(intra_send_counts, 0, node_size * sizeof(int));
        memset(inter_send_counts, 0, world_size * sizeof(int));
        for (int i = 0; i < n_c; i++) {
            int trg = items[i]->target_rank;
            if (trg < 0) {
                // Not in any rank's prefetch hint: dropped, its slot is freed below.
                items[i]->data_exchange_type = -1;
            }
            else if (client_info.rank_to_node_id_map[trg] == client_info.node_manager_rank) {
                items[i]->data_exchange_type = 0; // same node: hand over the slot
                intra_send_counts[client_info.world_to_node_rank_map[trg]]++;
                n_intra++;
            }
            else {
                items[i]->data_exchange_type = 1; // other node: send the data
                inter_send_counts[trg]++;
                n_inter++;
            }
        }

        // ---- Intra-node: send item metadata incl. slot index, no data copy ----
        tmp_timer = MPI_Wtime();
        int                  intra_recv_total;
        pdc_object_metadata *intra_send_meta =
            (pdc_object_metadata *)PDC_calloc(n_intra + 1, sizeof(pdc_object_metadata));
        counts_to_displs(node_size, intra_send_counts, intra_sdispls);
        pack_metadata(items, n_c, 0, node_size, intra_sdispls, fill, intra_send_meta, NULL);

        pdc_object_metadata *intra_recv_meta =
            exchange_metadata(client_cache_node_comm, node_size, metadata_type, intra_send_meta,
                              intra_send_counts, intra_sdispls, intra_recv_counts, intra_rdispls,
                              &intra_recv_total);
        for (int k = 0; k < intra_recv_total; k++)
            add_received_item(&intra_recv_meta[k], intra_recv_meta[k].slot_idx, 0);

        free(intra_send_meta);
        free(intra_recv_meta);
        n_intra_sent += n_intra;
        n_intra_recv += intra_recv_total;
        intra_time += MPI_Wtime() - tmp_timer;

        // ---- Inter-node step 1: exchange item metadata (small) ----
        tmp_timer = MPI_Wtime();
        int                  inter_recv_total;
        pdc_object_metadata *send_meta =
            (pdc_object_metadata *)PDC_calloc(n_inter + 1, sizeof(pdc_object_metadata));
        pdc_object_data **send_item =
            (pdc_object_data **)PDC_malloc((n_inter + 1) * sizeof(pdc_object_data *));

        // send_meta/send_item are grouped by target rank, in a fixed order; the
        // item's position within its target's group is its MPI tag below.
        counts_to_displs(world_size, inter_send_counts, inter_sdispls);
        pack_metadata(items, n_c, 1, world_size, inter_sdispls, fill, send_meta, send_item);

        pdc_object_metadata *recv_meta =
            exchange_metadata(client_cache_world_comm, world_size, metadata_type, send_meta,
                              inter_send_counts, inter_sdispls, inter_recv_counts, inter_rdispls,
                              &inter_recv_total);
        inter_meta_time += MPI_Wtime() - tmp_timer;

        // ---- Inter-node step 2: data straight from/into shared-window slots ----
        tmp_timer          = MPI_Wtime();
        MPI_Request *rreq  = (MPI_Request *)PDC_malloc((inter_recv_total + 1) * sizeof(MPI_Request));
        MPI_Request *sreq  = (MPI_Request *)PDC_malloc((n_inter + 1) * sizeof(MPI_Request));
        int *        rslot = (int *)PDC_malloc((inter_recv_total + 1) * sizeof(int));

        // Reserve a slot for every incoming item before any data moves, and post
        // the receive straight into it. Tag = position within the sender's group,
        // so a receive can only ever match the item its metadata describes.
        for (int s = 0; s < world_size; s++) {
            for (int j = 0; j < inter_recv_counts[s]; j++) {
                int k    = inter_rdispls[s] + j;
                rslot[k] = pop_free_slot();
                if (rslot[k] == SLOT_INVALID) {
                    rreq[k] = MPI_REQUEST_NULL; // received into dl_drop_buf after the others
                    continue;
                }
                MPI_Irecv(get_data_ptr(rslot[k]), (int)recv_meta[k].reg_buf_size, MPI_BYTE, s, j,
                          client_cache_world_comm, &rreq[k]);
            }
        }

        for (int t = 0; t < world_size; t++) {
            for (int j = 0; j < inter_send_counts[t]; j++) {
                int k = inter_sdispls[t] + j;
                MPI_Isend(get_data_ptr(send_item[k]->slot_idx), (int)send_meta[k].reg_buf_size, MPI_BYTE, t,
                          j, client_cache_world_comm, &sreq[k]);
            }
        }

        MPI_Waitall(inter_recv_total, rreq, MPI_STATUSES_IGNORE);

        // Items with no free slot on this node: drain them one at a time and drop
        // them, as before (the senders' Isends must still be matched).
        for (int s = 0; s < world_size; s++) {
            for (int j = 0; j < inter_recv_counts[s]; j++) {
                int k = inter_rdispls[s] + j;
                if (rslot[k] != SLOT_INVALID)
                    continue;
                if (dl_drop_buf == NULL)
                    dl_drop_buf = (char *)PDC_malloc(MAX_ITEM_SIZE);
                MPI_Recv(dl_drop_buf, (int)recv_meta[k].reg_buf_size, MPI_BYTE, s, j, client_cache_world_comm,
                         MPI_STATUS_IGNORE);
                n_dropped++;
            }
        }

        // Only now do received items become visible in the cache list.
        for (int k = 0; k < inter_recv_total; k++) {
            if (rslot[k] != SLOT_INVALID)
                add_received_item(&recv_meta[k], rslot[k], 0);
        }

        // Sent slots may only be freed once their sends have completed.
        MPI_Waitall(n_inter, sreq, MPI_STATUSES_IGNORE);

        n_inter_sent += n_inter;
        n_inter_recv += inter_recv_total;
        free(rreq);
        free(sreq);
        free(rslot);
        free(send_meta);
        free(recv_meta);
        free(send_item);
        inter_data_time += MPI_Wtime() - tmp_timer;

        // ---- Remove this round's old items ----
        // Same-node items keep their slot (the receiver now owns it); slots of
        // items sent to another node, or dropped, return to the node's free list.
        tmp_timer = MPI_Wtime();
        for (int i = 0; i < n_c; i++) {
            if (items[i]->data_exchange_type != 0)
                push_free_slot(items[i]->slot_idx);
            pdc_region_dl_delete(items[i]);
            free(items[i]);
        }
        free(items);
        delete_time += MPI_Wtime() - tmp_timer;
    }

    MPI_Type_free(&metadata_type);
    free(counts);

    printf("[RANK %d] data exchange: rounds %d, intra-node sent %d recv %d, inter-node sent %d recv %d, "
           "dropped (no free slot) %d | intra %.3fs, inter metadata %.3fs, inter data %.3fs, delete %.3fs\n",
           client_info.world_rank, NUM_CHUNKS, n_intra_sent, n_intra_recv, n_inter_sent, n_inter_recv,
           n_dropped, intra_time, inter_meta_time, inter_data_time, delete_time);

    pdc_region_cache_timelog(start, "pdc_region_dl_data_exchange - data exchange total time");

done:
    fflush(stdout);
    FUNC_LEAVE(ret_value);
}

/**************************************************************************/
/* Data exchange overlapped with computation: the exchange above in small */
/* pieces, one batch exchange per training step (per progress call) for   */
/* the items read in that step's batch. Each batch exchange:              */
/*   1. the items' metadata, in one fixed-size block per destination,     */
/*      goes out in one MPI_Ialltoall                                     */
/*      (MPI_Alltoall with PDC_OVERLAP_METADATA=blocking)                 */
/*   2. once it has arrived, same-node items hand over their slot and     */
/*      items for another node are sent slot to slot (Isend / Irecv)      */
/* and the data moves while the application computes. Batch exchanges are */
/* collective: every rank calls progress once per step, as often as every */
/* other rank.                                                            */
/**************************************************************************/

typedef struct {
    int same_node_sent, other_node_sent, same_node_received, other_node_received, missing;
} overlap_stats;

// This rank's share of the node's pending_inserts: this epoch's server reads still
// to come, for which receives of next-epoch data must leave free slots
static int overlap_pending_inserts;

static void
print_overlap_summary(overlap_stats *stats, int open_sends, int open_receives, double end_wait)
{
    printf("[RANK %d] overlapped exchange: sent same-node %d other-node %d, received same-node %d "
           "other-node %d, missing %d | still open when the epoch ended: %d sends, %d receives, "
           "end wait %.3fs\n",
           client_info.world_rank, stats->same_node_sent, stats->other_node_sent, stats->same_node_received,
           stats->other_node_received, stats->missing, open_sends, open_receives, end_wait);
}

// This epoch's copy of an object (not one received for the next epoch)
static pdc_object_data *
find_current_item(pdcid_t obj_id)
{
    pdc_object_data *obj_cache_iter;

    for (obj_cache_iter = client_info.local_cache_list_head; obj_cache_iter != NULL;
         obj_cache_iter = obj_cache_iter->next) {
        if (obj_cache_iter->obj_id == obj_id && obj_cache_iter->transfer_state != 2)
            return obj_cache_iter;
    }
    return NULL;
}

// Every cached item becomes this epoch's data: not read or sent yet
static void
reset_item_states(void)
{
    pdc_object_data *obj_cache_iter;

    for (obj_cache_iter = client_info.local_cache_list_head; obj_cache_iter != NULL;
         obj_cache_iter = obj_cache_iter->next) {
        obj_cache_iter->consumed       = 0;
        obj_cache_iter->transfer_state = 0;
    }
}

static void
node_pending_inserts_add(int delta)
{
    int previous;

    MPI_Win_lock(MPI_LOCK_EXCLUSIVE, 0, 0, client_info.node_shared_data_win);
    MPI_Fetch_and_op(&delta, &previous, MPI_INT, 0, offsetof(SharedMemoryHeader, pending_inserts), MPI_SUM,
                     client_info.node_shared_data_win);
    MPI_Win_unlock(0, client_info.node_shared_data_win);
}

static int
node_pending_inserts_get(void)
{
    int value;

    MPI_Win_lock(MPI_LOCK_SHARED, 0, 0, client_info.node_shared_data_win);
    MPI_Fetch_and_op(NULL, &value, MPI_INT, 0, offsetof(SharedMemoryHeader, pending_inserts), MPI_NO_OP,
                     client_info.node_shared_data_win);
    MPI_Win_unlock(0, client_info.node_shared_data_win);
    return value;
}

// Every rank's id list, concatenated in rank order: rank r's starts at displs[r]
static pdcid_t *
allgather_ids(MPI_Comm comm, pdcid_t *ids, int len, int *counts, int *displs, int *total)
{
    pdcid_t *all_ids;

    MPI_Allgather(&len, 1, MPI_INT, counts, 1, MPI_INT, comm);
    *total  = counts_to_displs(client_info.world_size, counts, displs);
    all_ids = (pdcid_t *)PDC_malloc((*total + 1) * sizeof(pdcid_t));
    MPI_Allgatherv(ids, len, MPI_UINT64_T, all_ids, counts, displs, MPI_UINT64_T, comm);
    return all_ids;
}

// A batch exchange's data messages are tagged (id % BATCH_TAG_RANGE) * BATCH_ITEMS_PER_DEST
// + position in the block, which stays below the MPI minimum tag limit of 32767
#define BATCH_TAG_RANGE 4000

// The exchange of one batch's items; open until all of its transfers have finished
typedef struct batch_exchange {
    int         id;
    int         metadata_arrived;
    MPI_Request metadata_request;
    int         num_open; // sends and receives not finished yet

    // world_size blocks of BATCH_ITEMS_PER_DEST entries, block r for rank r. An entry
    // with valid = 1 is an item; once the metadata has arrived, one still in transfer.
    pdc_object_metadata *send_metadata;
    pdc_object_metadata *recv_metadata;
    pdc_object_data **   send_items;
    MPI_Request *        send_requests;
    MPI_Request *        recv_requests;
    int *                recv_slots; // SLOT_INVALID: waiting for a free slot

    struct batch_exchange *next;
} batch_exchange;

// Who reads an object next epoch
typedef struct {
    pdcid_t obj_id;
    int     owner;
} next_epoch_owner;

static struct {
    int               active;
    int               comm_created;
    MPI_Comm          comm;
    int               blocking_metadata;
    next_epoch_owner *owners; // every rank's next-epoch objects, sorted by obj_id
    int               num_owners;
    int               num_batches; // batch exchanges started this epoch
    batch_exchange *  first_batch; // open batch exchanges, oldest first
    batch_exchange *  last_batch;
    overlap_stats     stats;
} overlap_exchange;

// Which items a batch exchange takes
enum { BATCH_READ_ITEMS, BATCH_ALL_ITEMS, BATCH_NO_ITEMS };

static int
compare_owners(const void *a, const void *b)
{
    pdcid_t x = ((const next_epoch_owner *)a)->obj_id;
    pdcid_t y = ((const next_epoch_owner *)b)->obj_id;
    return (x > y) - (x < y);
}

// Rank that reads obj_id next epoch, -1 if none
static int
next_owner(pdcid_t obj_id)
{
    next_epoch_owner key, *found;

    key.obj_id = obj_id;
    found      = (next_epoch_owner *)bsearch(&key, overlap_exchange.owners, overlap_exchange.num_owners,
                                        sizeof(next_epoch_owner), compare_owners);
    return found != NULL ? found->owner : -1;
}

static int
read_next_epoch(pdcid_t obj_id)
{
    return next_owner(obj_id) == client_info.world_rank;
}

// Drop cached items this rank does not read next epoch and that are not being sent
static void
free_unneeded_items(void)
{
    pdc_object_data *obj_cache_iter = client_info.local_cache_list_head, *obj_cache_item;

    while (obj_cache_iter != NULL) {
        obj_cache_item = obj_cache_iter;
        obj_cache_iter = obj_cache_iter->next;
        if (obj_cache_item->transfer_state != 0 || read_next_epoch(obj_cache_item->obj_id))
            continue;
        pdc_region_dl_delete(obj_cache_item);
        push_free_slot(obj_cache_item->slot_idx);
        free(obj_cache_item);
    }
}

static int
same_node(int rank)
{
    return client_info.rank_to_node_id_map[rank] == client_info.node_manager_rank;
}

static int
batch_exchange_tag(batch_exchange *batch, int k)
{
    return (batch->id % BATCH_TAG_RANGE) * BATCH_ITEMS_PER_DEST + k % BATCH_ITEMS_PER_DEST;
}

// Starts the next batch exchange; returns how many items did not fit (a later one takes them)
static int
batch_exchange_start(int which)
{
    int            n           = client_info.world_size * BATCH_ITEMS_PER_DEST, left = 0;
    int            block_bytes = BATCH_ITEMS_PER_DEST * (int)sizeof(pdc_object_metadata);
    int *          fill        = (int *)PDC_calloc(client_info.world_size, sizeof(int));
    batch_exchange *batch       = (batch_exchange *)PDC_calloc(1, sizeof(batch_exchange));

    batch->id            = overlap_exchange.num_batches++;
    batch->send_metadata = (pdc_object_metadata *)PDC_calloc(n, sizeof(pdc_object_metadata));
    batch->recv_metadata = (pdc_object_metadata *)PDC_calloc(n, sizeof(pdc_object_metadata));
    batch->send_items    = (pdc_object_data **)PDC_calloc(n, sizeof(pdc_object_data *));
    batch->send_requests = (MPI_Request *)PDC_malloc(n * sizeof(MPI_Request));
    batch->recv_requests = (MPI_Request *)PDC_malloc(n * sizeof(MPI_Request));
    batch->recv_slots    = (int *)PDC_malloc(n * sizeof(int));
    for (int k = 0; k < n; k++) {
        batch->send_requests[k] = MPI_REQUEST_NULL;
        batch->recv_requests[k] = MPI_REQUEST_NULL;
        batch->recv_slots[k]    = SLOT_INVALID;
    }

    for (pdc_object_data *item = client_info.local_cache_list_head; item != NULL && which != BATCH_NO_ITEMS;
         item                  = item->next) {
        if (item->transfer_state != 0 || (which == BATCH_READ_ITEMS && !item->consumed))
            continue;
        int target = next_owner(item->obj_id);
        if (target < 0 || target == client_info.world_rank)
            continue;
        if (fill[target] == BATCH_ITEMS_PER_DEST) {
            left++;
            continue;
        }
        int k = target * BATCH_ITEMS_PER_DEST + fill[target]++;
        metadata_from_item(&batch->send_metadata[k], item);
        batch->send_items[k] = item;
        item->transfer_state = 1; // being sent: never evicted
    }
    free(fill);

    if (overlap_exchange.blocking_metadata) {
        MPI_Alltoall(batch->send_metadata, block_bytes, MPI_BYTE, batch->recv_metadata, block_bytes, MPI_BYTE,
                     overlap_exchange.comm);
        batch->metadata_request = MPI_REQUEST_NULL;
    }
    else
        MPI_Ialltoall(batch->send_metadata, block_bytes, MPI_BYTE, batch->recv_metadata, block_bytes,
                      MPI_BYTE, overlap_exchange.comm, &batch->metadata_request);

    if (overlap_exchange.last_batch != NULL)
        overlap_exchange.last_batch->next = batch;
    else
        overlap_exchange.first_batch = batch;
    overlap_exchange.last_batch = batch;
    return left;
}

// A batch exchange's metadata has arrived: same-node items hand over their slot, data for
// other nodes starts moving (its receives are posted once a slot is free)
static void
batch_exchange_metadata_arrived(batch_exchange *batch)
{
    int n = client_info.world_size * BATCH_ITEMS_PER_DEST;

    batch->metadata_arrived = 1;
    for (int k = 0; k < n; k++) {
        int peer = k / BATCH_ITEMS_PER_DEST;

        if (batch->send_metadata[k].valid) {
            pdc_object_data *item = batch->send_items[k];
            if (same_node(peer)) {
                // The receiver now owns the slot
                pdc_region_dl_delete(item);
                free(item);
                batch->send_metadata[k].valid = 0;
                overlap_exchange.stats.same_node_sent++;
            }
            else {
                MPI_Isend(get_data_ptr(item->slot_idx), (int)item->reg_buf_size, MPI_BYTE, peer,
                          batch_exchange_tag(batch, k), overlap_exchange.comm, &batch->send_requests[k]);
                batch->num_open++;
            }
        }

        if (batch->recv_metadata[k].valid) {
            if (same_node(peer)) {
                add_received_item(&batch->recv_metadata[k], batch->recv_metadata[k].slot_idx, 2);
                batch->recv_metadata[k].valid = 0;
                overlap_exchange.stats.same_node_received++;
            }
            else
                batch->num_open++;
        }
    }
}

// Tests a batch exchange's open sends and receives, posts receives that have a slot now;
// returns how many finished or started
static int
batch_exchange_test(batch_exchange *batch, int *available)
{
    int n = client_info.world_size * BATCH_ITEMS_PER_DEST, changed = 0, flag;

    for (int k = 0; k < n; k++) {
        if (batch->send_metadata[k].valid) {
            MPI_Test(&batch->send_requests[k], &flag, MPI_STATUS_IGNORE);
            if (flag) {
                // The data has its own copy on the other node: the slot returns to this node's free list
                pdc_object_data *item = batch->send_items[k];
                push_free_slot(item->slot_idx);
                pdc_region_dl_delete(item);
                free(item);
                batch->send_metadata[k].valid = 0;
                batch->num_open--;
                overlap_exchange.stats.other_node_sent++;
                changed++;
            }
        }

        pdc_object_metadata *received = &batch->recv_metadata[k];
        if (!received->valid)
            continue;
        if (batch->recv_slots[k] == SLOT_INVALID) {
            // Leave enough free slots for this epoch's remaining server reads on the node
            if (*available < 0)
                *available = count_free_slots() - node_pending_inserts_get();
            if (*available <= 0)
                continue;
            batch->recv_slots[k] = pop_free_slot();
            if (batch->recv_slots[k] == SLOT_INVALID) {
                *available = 0;
                continue;
            }
            (*available)--;
            MPI_Irecv(get_data_ptr(batch->recv_slots[k]), (int)received->reg_buf_size, MPI_BYTE,
                      k / BATCH_ITEMS_PER_DEST, batch_exchange_tag(batch, k), overlap_exchange.comm,
                      &batch->recv_requests[k]);
            changed++;
        }
        MPI_Test(&batch->recv_requests[k], &flag, MPI_STATUS_IGNORE);
        if (flag) {
            add_received_item(received, batch->recv_slots[k], 2);
            received->valid = 0;
            batch->num_open--;
            overlap_exchange.stats.other_node_received++;
            changed++;
        }
    }
    return changed;
}

static void
batch_exchange_free(batch_exchange *batch)
{
    free(batch->send_metadata);
    free(batch->recv_metadata);
    free(batch->send_items);
    free(batch->send_requests);
    free(batch->recv_requests);
    free(batch->recv_slots);
    free(batch);
}

// Advances every open batch exchange, oldest first, and frees the finished ones;
// returns how many sends and receives finished or started
static int
batch_exchanges_advance(void)
{
    int            changed = 0, flag, available = -1;
    batch_exchange *batch = overlap_exchange.first_batch, *previous = NULL, *next;

    while (batch != NULL) {
        next = batch->next;
        if (!batch->metadata_arrived) {
            MPI_Test(&batch->metadata_request, &flag, MPI_STATUS_IGNORE);
            if (flag) {
                batch_exchange_metadata_arrived(batch);
                changed++;
            }
        }
        if (batch->metadata_arrived) {
            changed += batch_exchange_test(batch, &available);
            if (batch->num_open == 0) {
                if (previous != NULL)
                    previous->next = next;
                else
                    overlap_exchange.first_batch = next;
                if (overlap_exchange.last_batch == batch)
                    overlap_exchange.last_batch = previous;
                batch_exchange_free(batch);
                batch = next;
                continue;
            }
        }
        previous = batch;
        batch    = next;
    }
    return changed;
}

// Sends and receives not finished: items not yet in a batch exchange count as open sends;
// receives are known once a batch exchange's metadata has arrived
static void
batch_exchanges_count_open(int *open_sends, int *open_receives)
{
    int n = client_info.world_size * BATCH_ITEMS_PER_DEST;

    *open_sends    = 0;
    *open_receives = 0;
    for (pdc_object_data *item = client_info.local_cache_list_head; item != NULL; item = item->next) {
        int target = next_owner(item->obj_id);
        if (item->transfer_state == 0 && target >= 0 && target != client_info.world_rank)
            (*open_sends)++;
    }
    for (batch_exchange *batch = overlap_exchange.first_batch; batch != NULL; batch = batch->next) {
        for (int k = 0; k < n; k++) {
            *open_sends += batch->send_metadata[k].valid;
            if (batch->metadata_arrived)
                *open_receives += batch->recv_metadata[k].valid;
        }
    }
}

perr_t
pdc_region_dl_overlap_begin(pdcid_t *cur_ids, int cur_len, pdcid_t *next_ids, int next_len)
{
    FUNC_ENTER(NULL);

    perr_t ret_value = SUCCEED;
    int    world_size = client_info.world_size, total_next;
    double start = MPI_Wtime();

    if (!client_info.client_cache_init)
        PGOTO_ERROR(FAIL, "pdc_region_dl_overlap_begin - client cache not initialized");
    if (overlap_exchange.active)
        PGOTO_ERROR(FAIL, "pdc_region_dl_overlap_begin - previous exchange has not ended");

    if (!overlap_exchange.comm_created) {
        const char *metadata = getenv("PDC_OVERLAP_METADATA");

        MPI_Comm_dup(client_cache_world_comm, &overlap_exchange.comm);
        overlap_exchange.comm_created      = 1;
        overlap_exchange.blocking_metadata = metadata != NULL && strcmp(metadata, "blocking") == 0;
        if (client_info.world_rank == 0)
            printf("[RANK 0] overlapped exchange mode: batch exchanges, metadata %s\n",
                   overlap_exchange.blocking_metadata ? "MPI_Alltoall" : "MPI_Ialltoall");
    }

    // Who reads every object next epoch
    int *    counts   = (int *)PDC_calloc(2 * (size_t)world_size, sizeof(int));
    int *    displs   = counts + world_size;
    pdcid_t *all_next = allgather_ids(overlap_exchange.comm, next_ids, next_len, counts, displs, &total_next);

    overlap_exchange.owners     = (next_epoch_owner *)PDC_malloc((total_next + 1) * sizeof(next_epoch_owner));
    overlap_exchange.num_owners = total_next;
    for (int r = 0; r < world_size; r++) {
        for (int i = 0; i < counts[r]; i++) {
            overlap_exchange.owners[displs[r] + i].obj_id = all_next[displs[r] + i];
            overlap_exchange.owners[displs[r] + i].owner  = r;
        }
    }
    qsort(overlap_exchange.owners, total_next, sizeof(next_epoch_owner), compare_owners);
    free(counts);
    free(all_next);

    // Items delivered by the previous exchange are this epoch's data now
    reset_item_states();

    // This epoch's server reads still to come (epoch 1: all of them); receives
    // for next-epoch data must leave this many slots free on the node.
    overlap_pending_inserts = 0;
    for (int i = 0; i < cur_len; i++) {
        if (find_current_item(cur_ids[i]) == NULL)
            overlap_pending_inserts++;
    }
    if (overlap_pending_inserts > 0)
        node_pending_inserts_add(overlap_pending_inserts);

    overlap_exchange.num_batches = 0;
    memset(&overlap_exchange.stats, 0, sizeof(overlap_exchange.stats));
    overlap_exchange.active = 1;

    pdc_region_cache_timelog(start, "pdc_region_dl_overlap_begin - total time");

done:
    fflush(stdout);
    FUNC_LEAVE(ret_value);
}

perr_t
pdc_region_dl_overlap_progress(void)
{
    FUNC_ENTER(NULL);

    perr_t ret_value = SUCCEED;

    if (!overlap_exchange.active)
        goto done;

    batch_exchange_start(BATCH_READ_ITEMS);
    batch_exchanges_advance();

done:
    FUNC_LEAVE(ret_value);
}

perr_t
pdc_region_dl_overlap_end(void)
{
    FUNC_ENTER(NULL);

    perr_t ret_value = SUCCEED;
    double start = MPI_Wtime(), last_change = MPI_Wtime();
    int    warned = 0, open_sends, open_receives, max_batches, left, max_left;

    if (!overlap_exchange.active)
        goto done;

    // No more server reads this epoch: release this rank's reserved slots
    if (overlap_pending_inserts > 0) {
        node_pending_inserts_add(-overlap_pending_inserts);
        overlap_pending_inserts = 0;
    }

    batch_exchanges_advance();
    batch_exchanges_count_open(&open_sends, &open_receives);

    // Batch exchanges pair up across ranks in the order they are started: a rank that
    // started fewer adds empty ones. (On client_cache_world_comm, so these
    // collectives never interleave with the batch exchanges' own.)
    MPI_Allreduce(&overlap_exchange.num_batches, &max_batches, 1, MPI_INT, MPI_MAX, client_cache_world_comm);
    while (overlap_exchange.num_batches < max_batches)
        batch_exchange_start(BATCH_NO_ITEMS);

    // Items not read this epoch still go to their next owner
    do {
        left = batch_exchange_start(BATCH_ALL_ITEMS);
        MPI_Allreduce(&left, &max_left, 1, MPI_INT, MPI_MAX, client_cache_world_comm);
    } while (max_left > 0);

    while (overlap_exchange.first_batch != NULL) {
        if (batch_exchanges_advance() > 0) {
            last_change = MPI_Wtime();
            continue;
        }
        // Receives may be waiting for a slot: free what nobody needs any more
        free_unneeded_items();
        if (!warned && MPI_Wtime() - last_change > 60.0) {
            int s, r;
            batch_exchanges_count_open(&s, &r);
            printf("[RANK %d] overlapped exchange: no progress for 60 s (%d sends, %d receives open); "
                   "MAX_SLOTS_PER_NODE may be too small for the next epoch's items\n",
                   client_info.world_rank, s, r);
            fflush(stdout);
            warned = 1;
        }
        usleep(50);
    }

    // What is left in the cache is this rank's data for the next epoch
    free_unneeded_items();
    reset_item_states();

    // Objects this rank reads next epoch that did not arrive are read from the servers
    for (int i = 0; i < overlap_exchange.num_owners; i++) {
        if (overlap_exchange.owners[i].owner == client_info.world_rank &&
            find_current_item(overlap_exchange.owners[i].obj_id) == NULL)
            overlap_exchange.stats.missing++;
    }

    print_overlap_summary(&overlap_exchange.stats, open_sends, open_receives, MPI_Wtime() - start);

    free(overlap_exchange.owners);
    overlap_exchange.owners     = NULL;
    overlap_exchange.num_owners = 0;
    overlap_exchange.active     = 0;

    pdc_region_cache_timelog(start, "pdc_region_dl_overlap_end - wait time");

done:
    fflush(stdout);
    FUNC_LEAVE(ret_value);
}

/**************************************************************************/
/* Hooks pdc_region_dl_insert calls while an exchange runs                */
/**************************************************************************/

int
pdc_region_dl_overlap_active(void)
{
    return overlap_exchange.active;
}

// No slot is free: finished sends free slots; otherwise evict (eviction never
// touches an item that is being sent)
void
pdc_region_dl_overlap_make_room(void)
{
    if (batch_exchanges_advance() > 0)
        return;
    pdc_region_cache_evict();
}

// One of this epoch's server reads has its slot now
void
pdc_region_dl_overlap_record_insert(void)
{
    if (overlap_pending_inserts > 0) {
        overlap_pending_inserts--;
        node_pending_inserts_add(-1);
    }
}
