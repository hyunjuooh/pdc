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
#include "pdc_obj_pkg.h"
#include "pdc_interface.h"
#include "pdc_transforms_pkg.h"
#include "pdc_client_connect.h"
#include "pdc_analysis_pkg.h"
#include <mpi.h>

typedef struct pdc_region_prefetch {
    struct pdc_region_info  *region_prefetch_info;
    struct pdc_region_prefetch *prev;
    struct pdc_region_prefetch *next;
} pdc_region_prefetch;

typedef struct transfer_request_pending{
    pdcid_t transfer_request_id;
    pdcid_t remote_reg_id;
    struct transfer_request_pending *prev;
    struct transfer_request_pending *next;
} transfer_request_pending;

static struct pdc_region_prefetch          *region_prefetch_list;
static struct transfer_request_pending     *transfer_request_pending_list;
static size_t total_buffer_size;

struct pdc_region_prefetch                 *global_region_prefetch_list = &region_prefetch_list;

int PDCregion_prefetch_init()
{
    total_buffer_size = 0;
    region_prefetch_list = NULL;
    transfer_request_pending_list = NULL;

    return 0;
}

int PDCregion_prefetch_status(pdcid_t reg_id){

    struct transfer_request_pending   *transfer_request_pending_iter = NULL;
    pdc_transfer_status_t      completed;
    perr_t                     ret;

    transfer_request_pending_iter = transfer_request_pending_list;

    while (transfer_request_pending_iter != NULL){
        if (transfer_request_pending_iter->remote_reg_id == reg_id){
            ret = PDCregion_transfer_wait(transfer_request_pending_iter->transfer_request_id);
            if (ret != SUCCEED) {
                printf("Failed region transfer start\n");
            }

            ret = PDCregion_transfer_close(transfer_request_pending_iter->transfer_request_id);
            if (ret != SUCCEED) {
                printf("Fail to region transfer close @ line %d\n", __LINE__);
            }

            DL_DELETE(transfer_request_pending_list, transfer_request_pending_iter);
            break;
        }

        transfer_request_pending_iter = transfer_request_pending_iter->next;
    }

    /*if (transfer_request_pending_iter == NULL){
        printf("No transfer request made for specific remote region id");
        return 0;
    }*/

    return 1;
}

void* PDCregion_prefetch_search(pdcid_t reg_id)
{
    struct pdc_region_prefetch        *reg_prefetch_iter;

    reg_prefetch_iter = region_prefetch_list;

    PDCregion_prefetch_status(reg_id);

    while (reg_prefetch_iter != NULL) {
        if (reg_prefetch_iter->region_prefetch_info->local_id == reg_id){
            if (reg_prefetch_iter->region_prefetch_info->local_prefetched == 1) {
                return reg_prefetch_iter->region_prefetch_info->buf;
            }
        }
        reg_prefetch_iter = reg_prefetch_iter->next;
    }

    return NULL;
}

// Prefetch specified region into the client local memory
// reg_id      target region to prefetch
// obj_id      object which contains the target region
// remote_reg  used for transfer

pdcid_t PDCregion_prefetch(size_t buf_size, pdcid_t local_reg_id, pdcid_t remote_reg_id, pdcid_t obj_id) {
    pdcid_t                    ret_value = 0;
    struct _pdc_id_info        *objinfo = NULL;
    struct _pdc_id_info        *local_reginfo = NULL;
    struct _pdc_id_info        *remote_reginfo = NULL;
    struct _pdc_obj_info       *obj;
    struct pdc_region_info     *local_reg, *remote_reg;
    //int                        *buf = (int *)malloc(buf_size);
    char                       *buf = (char *)malloc(sizeof(char) * buf_size);

    // used for prefetch management
    struct pdc_region_prefetch        *region_prefetch_item;
    struct pdc_region_info            *region_prefetch_info;

    struct transfer_request_pending    *transfer_request_item;

    int ndim, i;

    // indicate if the transfer was successful or not
    perr_t                     ret;

    // used for transfer
    pdcid_t                    transfer_request;

    FUNC_ENTER(NULL);

    local_reginfo = PDC_find_id(local_reg_id);
    local_reg = (struct pdc_region_info *)(local_reginfo->obj_ptr);

    if (local_reginfo == NULL)
        printf("No local_regioninfo found\n");

    remote_reginfo = PDC_find_id(remote_reg_id);
    remote_reg = (struct pdc_region_info *)(remote_reginfo->obj_ptr);

    if (remote_reginfo == NULL)
        printf("No remote_reginfo found\n");

    objinfo = PDC_find_id(obj_id);
    obj = (struct _pdc_obj_info *)(objinfo->obj_ptr);

    if (objinfo == NULL)
        printf("No objinfo found\n");

    //check if the region is already in the local memory
    if (PDCregion_prefetch_search(remote_reg_id) != NULL)
        goto done;

    // if region is not prefetched, transfer request from server and append to local region cache list
    memset(buf, 0, buf_size);

    transfer_request = PDCregion_transfer_create(buf, PDC_READ, obj_id, local_reg_id, remote_reg_id);

    ret = PDCregion_transfer_start(transfer_request);
    if (ret != SUCCEED) {
        printf("Failed region transfer start\n");
        goto done;
    }

    transfer_request_item = (struct transfer_request_pending *)malloc(sizeof(struct transfer_request_pending));
    transfer_request_item->transfer_request_id = transfer_request;
    transfer_request_item->remote_reg_id = remote_reg_id;
    DL_APPEND(transfer_request_pending_list, transfer_request_item);

    region_prefetch_item = (struct pdc_region_prefetch *)malloc(sizeof(struct pdc_region_prefetch));

    // memory allocation for region_prefetch_info
    region_prefetch_item->region_prefetch_info = (struct pdc_region_info *)malloc(sizeof(struct pdc_region_info));
    region_prefetch_info = region_prefetch_item->region_prefetch_info;

    region_prefetch_info->local_id = remote_reg_id;
    region_prefetch_info->ndim = remote_reg->ndim;
    ndim = remote_reg->ndim;
    region_prefetch_info->offset = (uint64_t *)malloc(sizeof(uint64_t) * ndim * 2);
    region_prefetch_info->size = region_prefetch_info->offset + ndim;
    region_prefetch_info->buf = buf;
    region_prefetch_info->unit = remote_reg->unit;
    
    // memory copy for offset, size
    memcpy(region_prefetch_info->offset, remote_reg->offset, sizeof(uint64_t) * ndim);
    memcpy(region_prefetch_info->size, remote_reg->size, sizeof(uint64_t) * ndim);

    // change metadata information to indicate region location
    region_prefetch_info->local_prefetched = remote_reg->local_prefetched = 1;

    // add the pdc region that was read from the transfer request to the pdc_region_cache_list
    DL_APPEND(region_prefetch_list, region_prefetch_item);

    total_buffer_size += buf_size;

done:
    fflush(stdout);
    FUNC_LEAVE(ret_value);
}

int
PDCregion_prefetch_free()
{
   struct pdc_region_prefetch  *region_prefetch_iter, *region_prefetch_tmp;

   region_prefetch_iter = region_prefetch_list;
   while (region_prefetch_iter != NULL) {
       free(region_prefetch_iter->region_prefetch_info);
       region_prefetch_tmp = region_prefetch_iter;
       region_prefetch_iter = region_prefetch_iter->next;
       free(region_prefetch_tmp);
   }

   return 0;
}
