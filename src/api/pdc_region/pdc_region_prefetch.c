#include <time.h>
#include <stdlib.h>
#include <unistd.h>
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
#include "pdc_client_connect.h"

pdcid_t *obj_prefetch_list;
int      obj_prefetch_list_len;

perr_t
PDCregion_receive_prefetch_hint(const pdcid_t *arr, int obj_array_len)
{
    perr_t ret_value = SUCCEED;

    int     i, item_idx = 0;
    int     current_rank = pdc_client_mpi_rank_g;
    pdcid_t obj_item;

    FUNC_ENTER(NULL);

    // for (i = 0; i < obj_array_len; i++) {
    //     printf("Received from C: %s\n", arr[i]);
    // }

    obj_prefetch_list_len = obj_array_len;
    obj_prefetch_list     = (pdcid_t *)PDC_malloc(obj_prefetch_list_len * sizeof(pdcid_t));
    memcpy(obj_prefetch_list, arr, obj_prefetch_list_len);

    if (pdc_client_mpi_rank_g == 0)
        printf("[RANK %d] PDCregion_receive_prefetch_hint - object list item length %d \n",
               pdc_client_mpi_rank_g, obj_prefetch_list_len);

    // // Get the list of object id for each rank that should be prefetched
    // for (i = 0; i < obj_prefetch_list_len; i++) {
    //     obj_item = PDCobj_open()
    //     obj_prefetch_list[item_idx] = obj_array[i];
    //     item_idx++;
    // }

done:
    fflush(stdout);
    FUNC_LEAVE(ret_value);
}

perr_t
PDCregion_prefetch_by_objid()
{
    perr_t ret_value = SUCCEED;
    int    i, is_cached;

    FUNC_ENTER(NULL);

    if (obj_prefetch_list == NULL) {
        if (pdc_client_mpi_rank_g == 0)
            printf("[RANK %d] pdc_region_prefetch_by_objid - object list not created\n",
                   pdc_client_mpi_rank_g);

        goto done;
    }

    ret_value = pdc_region_dl_collect_global_metadata();

    for (i = 0; i < obj_prefetch_list_len; i++) {
        is_cached = pdc_region_dl_global_search(obj_prefetch_list[i]);
    }

    free(obj_prefetch_list);

done:
    fflush(stdout);
    FUNC_LEAVE(ret_value);
}
