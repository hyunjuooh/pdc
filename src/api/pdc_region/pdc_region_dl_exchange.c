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
/* Data exchange overlapped with computation                              */
/*                                                                        */
/* An item is forwarded to its next-epoch owner as soon as its current    */
/* owner has read it this epoch (server read in epoch 1, cache hit after  */
/* that), so the exchange runs while the application computes.           */
/**************************************************************************/

// Where an object is read this epoch and who reads it next epoch
typedef struct {
    pdcid_t obj_id;
    int     current_owner;
    int     current_position; // position in the current owner's list, used for message tags
    int     next_owner;
} overlap_object_record;

typedef struct {
    pdcid_t             obj_id;
    int                 target_rank;
    int                 tag;
    int                 same_node;
    int                 state; // 0 waiting to be read, 1 being sent, 2 done
    pdc_object_data *   item;
    pdc_object_metadata metadata;
    MPI_Request         requests[2]; // metadata, data
} overlap_outgoing_item;

typedef struct {
    pdcid_t             obj_id;
    int                 source_rank;
    int                 tag;
    int                 same_node;
    int                 state; // 0 metadata posted, 1 waiting for a free slot, 2 data posted, 3 done
    int                 slot_idx;
    pdc_object_metadata metadata;
    MPI_Request         metadata_request;
    MPI_Request         data_request;
} overlap_incoming_item;

static struct {
    int                    active;
    int                    comm_created;
    MPI_Comm               comm;
    overlap_object_record *records;
    int                    num_records;
    pdcid_t *              next_ids_sorted; // this rank's next-epoch objects
    int                    num_next_ids;
    overlap_outgoing_item *outgoing;
    int                    num_outgoing;
    int                    num_outgoing_done;
    overlap_incoming_item *incoming;
    int                    num_incoming;
    int                    num_incoming_done;
    int                    pending_inserts; // this rank's share of the node's pending_inserts
    struct {
        int same_node_sent, other_node_sent, same_node_received, other_node_received, missing;
    } stats;
} overlap_exchange;

int
pdc_region_dl_overlap_active(void)
{
    return overlap_exchange.active;
}

static int
compare_object_records(const void *a, const void *b)
{
    pdcid_t x = ((const overlap_object_record *)a)->obj_id;
    pdcid_t y = ((const overlap_object_record *)b)->obj_id;
    return (x > y) - (x < y);
}

static int
compare_object_ids(const void *a, const void *b)
{
    pdcid_t x = *(const pdcid_t *)a;
    pdcid_t y = *(const pdcid_t *)b;
    return (x > y) - (x < y);
}

static overlap_object_record *
find_object_record(pdcid_t obj_id)
{
    overlap_object_record key;
    key.obj_id = obj_id;
    return (overlap_object_record *)bsearch(&key, overlap_exchange.records, overlap_exchange.num_records,
                                            sizeof(overlap_object_record), compare_object_records);
}

static int
read_next_epoch(pdcid_t obj_id)
{
    return bsearch(&obj_id, overlap_exchange.next_ids_sorted, overlap_exchange.num_next_ids, sizeof(pdcid_t),
                   compare_object_ids) != NULL;
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

// Called by pdc_region_dl_insert: one of this epoch's server reads has its slot now
void
pdc_region_dl_overlap_record_insert(void)
{
    if (overlap_exchange.pending_inserts > 0) {
        overlap_exchange.pending_inserts--;
        node_pending_inserts_add(-1);
    }
}

// item == NULL: this rank does not have the object, the metadata says so (valid = 0)
static void
overlap_start_send(overlap_outgoing_item *outgoing_item, pdc_object_data *item)
{
    pdc_object_metadata *metadata = &outgoing_item->metadata;

    if (item != NULL) {
        metadata_from_item(metadata, item);
        item->transfer_state = 1;
    }
    else {
        memset(metadata, 0, sizeof(pdc_object_metadata));
        metadata->obj_id = outgoing_item->obj_id;
    }

    MPI_Isend(metadata, sizeof(pdc_object_metadata), MPI_BYTE, outgoing_item->target_rank,
              2 * outgoing_item->tag, overlap_exchange.comm, &outgoing_item->requests[0]);
    outgoing_item->requests[1] = MPI_REQUEST_NULL;
    if (item != NULL && !outgoing_item->same_node)
        MPI_Isend(get_data_ptr(item->slot_idx), (int)item->reg_buf_size, MPI_BYTE, outgoing_item->target_rank,
                  2 * outgoing_item->tag + 1, overlap_exchange.comm, &outgoing_item->requests[1]);

    outgoing_item->item  = item;
    outgoing_item->state = 1;
}

// Forward every item read since the last call
static void
overlap_start_sends_for_read_items(void)
{
    for (int i = 0; i < overlap_exchange.num_outgoing; i++) {
        overlap_outgoing_item *outgoing_item = &overlap_exchange.outgoing[i];
        if (outgoing_item->state != 0)
            continue;
        pdc_object_data *item = find_current_item(outgoing_item->obj_id);
        if (item != NULL && item->consumed && item->transfer_state == 0)
            overlap_start_send(outgoing_item, item);
    }
}

// Returns how many sends completed
static int
overlap_test_sends(void)
{
    int completed = 0, flag;

    for (int i = 0; i < overlap_exchange.num_outgoing; i++) {
        overlap_outgoing_item *outgoing_item = &overlap_exchange.outgoing[i];
        if (outgoing_item->state != 1)
            continue;
        MPI_Testall(2, outgoing_item->requests, &flag, MPI_STATUSES_IGNORE);
        if (!flag)
            continue;

        if (outgoing_item->item != NULL) {
            // A same-node receiver now owns the slot; data sent to another node
            // has its own copy there, so the slot returns to this node's free list.
            if (outgoing_item->same_node)
                overlap_exchange.stats.same_node_sent++;
            else {
                push_free_slot(outgoing_item->item->slot_idx);
                overlap_exchange.stats.other_node_sent++;
            }
            pdc_region_dl_delete(outgoing_item->item);
            free(outgoing_item->item);
            outgoing_item->item = NULL;
        }
        outgoing_item->state = 2;
        overlap_exchange.num_outgoing_done++;
        completed++;
    }
    return completed;
}

// Returns how many receives changed state
static int
overlap_test_receives(void)
{
    int changed = 0, flag, available = -1;

    for (int i = 0; i < overlap_exchange.num_incoming; i++) {
        overlap_incoming_item *incoming_item = &overlap_exchange.incoming[i];

        if (incoming_item->state == 0) {
            MPI_Test(&incoming_item->metadata_request, &flag, MPI_STATUS_IGNORE);
            if (!flag)
                continue;
            changed++;
            if (!incoming_item->metadata.valid) {
                incoming_item->state = 3; // read from the servers next epoch instead
                overlap_exchange.num_incoming_done++;
                overlap_exchange.stats.missing++;
                continue;
            }
            if (incoming_item->same_node) {
                add_received_item(&incoming_item->metadata, incoming_item->metadata.slot_idx, 2);
                incoming_item->state = 3;
                overlap_exchange.num_incoming_done++;
                overlap_exchange.stats.same_node_received++;
                continue;
            }
            incoming_item->state = 1;
        }

        if (incoming_item->state == 1) {
            // Leave enough free slots for this epoch's remaining server reads on the node
            if (available < 0)
                available = count_free_slots() - node_pending_inserts_get();
            if (available <= 0)
                continue;
            incoming_item->slot_idx = pop_free_slot();
            if (incoming_item->slot_idx == SLOT_INVALID) {
                available = 0;
                continue;
            }
            available--;
            MPI_Irecv(get_data_ptr(incoming_item->slot_idx), (int)incoming_item->metadata.reg_buf_size,
                      MPI_BYTE, incoming_item->source_rank, 2 * incoming_item->tag + 1, overlap_exchange.comm,
                      &incoming_item->data_request);
            incoming_item->state = 2;
            changed++;
        }

        if (incoming_item->state == 2) {
            MPI_Test(&incoming_item->data_request, &flag, MPI_STATUS_IGNORE);
            if (!flag)
                continue;
            add_received_item(&incoming_item->metadata, incoming_item->slot_idx, 2);
            incoming_item->state = 3;
            overlap_exchange.num_incoming_done++;
            overlap_exchange.stats.other_node_received++;
            changed++;
        }
    }
    return changed;
}

// Called by pdc_region_dl_insert when no slot is free: completed sends free
// slots; otherwise evict (eviction never touches an item that is being sent).
void
pdc_region_dl_overlap_make_room(void)
{
    if (overlap_test_sends() > 0)
        return;
    pdc_region_cache_evict();
}

// Drop cached items this rank does not read next epoch and that are not being sent
static void
overlap_free_unneeded_items(void)
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

// Every rank's id list, concatenated in rank order: rank r's starts at displs[r]
static pdcid_t *
allgather_ids(pdcid_t *ids, int len, int *counts, int *displs, int *total)
{
    pdcid_t *all_ids;

    MPI_Allgather(&len, 1, MPI_INT, counts, 1, MPI_INT, overlap_exchange.comm);
    *total  = counts_to_displs(client_info.world_size, counts, displs);
    all_ids = (pdcid_t *)PDC_malloc((*total + 1) * sizeof(pdcid_t));
    MPI_Allgatherv(ids, len, MPI_UINT64_T, all_ids, counts, displs, MPI_UINT64_T, overlap_exchange.comm);
    return all_ids;
}

perr_t
pdc_region_dl_overlap_begin(pdcid_t *cur_ids, int cur_len, pdcid_t *next_ids, int next_len)
{
    FUNC_ENTER(NULL);

    perr_t ret_value = SUCCEED;
    int    world_size = client_info.world_size, my_rank = client_info.world_rank;
    int    total_cur, total_next;
    double start = MPI_Wtime();

    if (!client_info.client_cache_init)
        PGOTO_ERROR(FAIL, "pdc_region_dl_overlap_begin - client cache not initialized");
    if (overlap_exchange.active)
        PGOTO_ERROR(FAIL, "pdc_region_dl_overlap_begin - previous exchange has not ended");

    if (!overlap_exchange.comm_created) {
        MPI_Comm_dup(client_cache_world_comm, &overlap_exchange.comm);
        overlap_exchange.comm_created = 1;
    }

    // Every rank's lists, so each rank knows the current and next owner of every object
    int *cur_counts  = (int *)PDC_calloc(4 * (size_t)world_size, sizeof(int));
    int *next_counts = cur_counts + world_size;
    int *cur_displs  = cur_counts + 2 * world_size;
    int *next_displs = cur_counts + 3 * world_size;

    pdcid_t *all_cur  = allgather_ids(cur_ids, cur_len, cur_counts, cur_displs, &total_cur);
    pdcid_t *all_next = allgather_ids(next_ids, next_len, next_counts, next_displs, &total_next);

    overlap_exchange.records =
        (overlap_object_record *)PDC_malloc((total_cur + 1) * sizeof(overlap_object_record));
    overlap_exchange.num_records = total_cur;
    for (int r = 0; r < world_size; r++) {
        for (int i = 0; i < cur_counts[r]; i++) {
            overlap_object_record *record = &overlap_exchange.records[cur_displs[r] + i];
            record->obj_id                = all_cur[cur_displs[r] + i];
            record->current_owner         = r;
            record->current_position      = i;
            record->next_owner            = -1;
        }
    }
    qsort(overlap_exchange.records, total_cur, sizeof(overlap_object_record), compare_object_records);
    for (int r = 0; r < world_size; r++) {
        for (int i = 0; i < next_counts[r]; i++) {
            overlap_object_record *record = find_object_record(all_next[next_displs[r] + i]);
            if (record != NULL)
                record->next_owner = r;
        }
    }

    overlap_exchange.num_next_ids    = next_len;
    overlap_exchange.next_ids_sorted = (pdcid_t *)PDC_malloc((next_len + 1) * sizeof(pdcid_t));
    memcpy(overlap_exchange.next_ids_sorted, next_ids, next_len * sizeof(pdcid_t));
    qsort(overlap_exchange.next_ids_sorted, next_len, sizeof(pdcid_t), compare_object_ids);

    // Items delivered by the previous exchange are this epoch's data now
    reset_item_states();

    // Outgoing: this rank's objects that another rank reads next epoch
    overlap_exchange.outgoing =
        (overlap_outgoing_item *)PDC_calloc(cur_len + 1, sizeof(overlap_outgoing_item));
    overlap_exchange.num_outgoing = 0;
    for (int i = 0; i < cur_len; i++) {
        overlap_object_record *record = find_object_record(cur_ids[i]);
        if (record == NULL || record->current_owner != my_rank || record->current_position != i ||
            record->next_owner < 0 || record->next_owner == my_rank)
            continue;
        overlap_outgoing_item *outgoing_item = &overlap_exchange.outgoing[overlap_exchange.num_outgoing++];
        outgoing_item->obj_id                = record->obj_id;
        outgoing_item->target_rank           = record->next_owner;
        outgoing_item->tag                   = i;
        outgoing_item->same_node =
            client_info.rank_to_node_id_map[record->next_owner] == client_info.node_manager_rank;
        outgoing_item->requests[0] = MPI_REQUEST_NULL;
        outgoing_item->requests[1] = MPI_REQUEST_NULL;
    }

    // Incoming: objects this rank reads next epoch that another rank holds now.
    // Metadata messages are small, so all their receives are posted right away.
    overlap_exchange.incoming =
        (overlap_incoming_item *)PDC_calloc(next_len + 1, sizeof(overlap_incoming_item));
    overlap_exchange.num_incoming = 0;
    for (int i = 0; i < next_len; i++) {
        overlap_object_record *record = find_object_record(next_ids[i]);
        if (record == NULL || record->next_owner != my_rank || record->current_owner == my_rank)
            continue;
        overlap_incoming_item *incoming_item = &overlap_exchange.incoming[overlap_exchange.num_incoming++];
        incoming_item->obj_id                = record->obj_id;
        incoming_item->source_rank           = record->current_owner;
        incoming_item->tag                   = record->current_position;
        incoming_item->same_node =
            client_info.rank_to_node_id_map[record->current_owner] == client_info.node_manager_rank;
        incoming_item->slot_idx     = SLOT_INVALID;
        incoming_item->data_request = MPI_REQUEST_NULL;
        MPI_Irecv(&incoming_item->metadata, sizeof(pdc_object_metadata), MPI_BYTE, incoming_item->source_rank,
                  2 * incoming_item->tag, overlap_exchange.comm, &incoming_item->metadata_request);
    }

    // This epoch's server reads still to come (epoch 1: all of them); receives
    // for next-epoch data must leave this many slots free on the node.
    overlap_exchange.pending_inserts = 0;
    for (int i = 0; i < cur_len; i++) {
        if (find_current_item(cur_ids[i]) == NULL)
            overlap_exchange.pending_inserts++;
    }
    if (overlap_exchange.pending_inserts > 0)
        node_pending_inserts_add(overlap_exchange.pending_inserts);

    overlap_exchange.num_outgoing_done = 0;
    overlap_exchange.num_incoming_done = 0;
    memset(&overlap_exchange.stats, 0, sizeof(overlap_exchange.stats));
    overlap_exchange.active = 1;

    free(cur_counts);
    free(all_cur);
    free(all_next);

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

    overlap_start_sends_for_read_items();
    overlap_test_sends();
    overlap_test_receives();

done:
    FUNC_LEAVE(ret_value);
}

perr_t
pdc_region_dl_overlap_end(void)
{
    FUNC_ENTER(NULL);

    perr_t ret_value = SUCCEED;
    double start     = MPI_Wtime(), last_change = MPI_Wtime();
    int    warned    = 0, open_outgoing, open_incoming;

    if (!overlap_exchange.active)
        goto done;

    // No more server reads this epoch: release this rank's reserved slots
    if (overlap_exchange.pending_inserts > 0) {
        node_pending_inserts_add(-overlap_exchange.pending_inserts);
        overlap_exchange.pending_inserts = 0;
    }

    // Items not read this epoch still go to their next owner; items this rank
    // never had are reported, so the receiver reads them from the servers.
    for (int i = 0; i < overlap_exchange.num_outgoing; i++) {
        overlap_outgoing_item *outgoing_item = &overlap_exchange.outgoing[i];
        if (outgoing_item->state == 0) {
            pdc_object_data *item = find_current_item(outgoing_item->obj_id);
            overlap_start_send(outgoing_item, (item != NULL && item->transfer_state == 0) ? item : NULL);
        }
    }

    open_outgoing = overlap_exchange.num_outgoing - overlap_exchange.num_outgoing_done;
    open_incoming = overlap_exchange.num_incoming - overlap_exchange.num_incoming_done;

    while (overlap_exchange.num_outgoing_done < overlap_exchange.num_outgoing ||
           overlap_exchange.num_incoming_done < overlap_exchange.num_incoming) {
        if (overlap_test_sends() + overlap_test_receives() > 0) {
            last_change = MPI_Wtime();
            continue;
        }
        // Receives may be waiting for a slot: free what nobody needs any more
        overlap_free_unneeded_items();
        if (!warned && MPI_Wtime() - last_change > 60.0) {
            printf("[RANK %d] overlapped exchange: no progress for 60 s (%d sends, %d receives open); "
                   "MAX_SLOTS_PER_NODE may be too small for the next epoch's items\n",
                   client_info.world_rank, overlap_exchange.num_outgoing - overlap_exchange.num_outgoing_done,
                   overlap_exchange.num_incoming - overlap_exchange.num_incoming_done);
            fflush(stdout);
            warned = 1;
        }
        usleep(50);
    }

    // What is left in the cache is this rank's data for the next epoch
    overlap_free_unneeded_items();
    reset_item_states();

    printf("[RANK %d] overlapped exchange: sent same-node %d other-node %d, received same-node %d "
           "other-node %d, missing %d | still open when the epoch ended: %d sends, %d receives, "
           "end wait %.3fs\n",
           client_info.world_rank, overlap_exchange.stats.same_node_sent,
           overlap_exchange.stats.other_node_sent, overlap_exchange.stats.same_node_received,
           overlap_exchange.stats.other_node_received, overlap_exchange.stats.missing, open_outgoing,
           open_incoming, MPI_Wtime() - start);

    free(overlap_exchange.records);
    free(overlap_exchange.next_ids_sorted);
    free(overlap_exchange.outgoing);
    free(overlap_exchange.incoming);
    overlap_exchange.records         = NULL;
    overlap_exchange.next_ids_sorted = NULL;
    overlap_exchange.outgoing        = NULL;
    overlap_exchange.incoming        = NULL;
    overlap_exchange.active          = 0;

    pdc_region_cache_timelog(start, "pdc_region_dl_overlap_end - wait time");

done:
    fflush(stdout);
    FUNC_LEAVE(ret_value);
}
