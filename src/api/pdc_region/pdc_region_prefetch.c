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

pdcid_t * obj_prefetch_list;
uint64_t *reg_offset_list;
int       reg_dim;
int       obj_prefetch_list_len;

perr_t
PDCregion_receive_prefetch_hint(const pdcid_t *arr, const pdcid_t *arr2, int obj_array_len)
{
    perr_t ret_value = SUCCEED;

    struct _pdc_id_info *   objinfo2, *reginfo2;
    struct _pdc_obj_info *  obj2;
    struct pdc_region_info *reg2;
    uint64_t *              ptr;
    int                     i;

    FUNC_ENTER(NULL);

    obj_prefetch_list_len = obj_array_len;
    obj_prefetch_list     = (pdcid_t *)PDC_malloc(obj_prefetch_list_len * sizeof(pdcid_t));

    // If arr2 is null, it will read the whole region of the object
    if (arr2 != NULL) {
        reginfo2 = PDC_find_id(arr2[0]);
        if (reginfo2 == NULL)
            PGOTO_ERROR(FAIL, "cannot locate remote region ID");

        reg2 = (struct pdc_region_info *)(reginfo2->obj_ptr);

        reg_dim = reg2->ndim;

        reg_offset_list = (uint64_t *)PDC_malloc(sizeof(uint64_t) * reg_dim * 2 * obj_prefetch_list_len);
        ptr             = reg_offset_list;
    }

    // Convert received object id
    for (i = 0; i < obj_prefetch_list_len; i++) {
        objinfo2 = PDC_find_id(arr[i]);
        if (objinfo2 == NULL)
            PGOTO_ERROR(FAIL, "cannot locate remote object ID");
        obj2 = (struct _pdc_obj_info *)(objinfo2->obj_ptr);

        obj_prefetch_list[i] = obj2->obj_info_pub->meta_id;

        if (reg_offset_list != NULL) {
            reginfo2 = PDC_find_id(arr2[i]);
            if (reginfo2 == NULL)
                PGOTO_ERROR(FAIL, "cannot locate remote region ID");

            reg2 = (struct pdc_region_info *)(reginfo2->obj_ptr);
            memcpy(ptr, reg2->offset, sizeof(uint64_t) * reg_dim);
            ptr += reg_dim;

            memcpy(ptr, reg2->size, sizeof(uint64_t) * reg_dim);
            ptr += reg_dim;
        }

        /*if (pdc_client_mpi_rank_g == 0)
            printf("[RANK %d] Object id item %d finished\n",
                   pdc_client_mpi_rank_g, i);
        fflush(stdout);*/
    }

done:
    fflush(stdout);
    FUNC_LEAVE(ret_value);
}

perr_t
PDCregion_prefetch_by_objid()
{
    perr_t    ret_value = SUCCEED;
    uint64_t *offset, *size, *ptr;
    int       i, is_cached;

    FUNC_ENTER(NULL);

    if (obj_prefetch_list == NULL) {
        if (pdc_client_mpi_rank_g == 0)
            printf("[RANK %d] pdc_region_prefetch_by_objid - object list not created\n",
                   pdc_client_mpi_rank_g);

        goto done;
    }

    ret_value = pdc_region_dl_collect_global_metadata();

    ret_value = pdc_region_dl_list_init();

    ptr    = reg_offset_list;
    offset = (uint64_t *)PDC_malloc(sizeof(uint64_t) * reg_dim);
    size   = (uint64_t *)PDC_malloc(sizeof(uint64_t) * reg_dim);

    for (i = 0; i < obj_prefetch_list_len; i++) {
        if (ptr != NULL) {
            memcpy(offset, ptr, sizeof(uint64_t) * reg_dim);
            ptr += reg_dim;

            memcpy(size, ptr, sizeof(uint64_t) * reg_dim);
            ptr += reg_dim;

            if (pdc_client_mpi_rank_g == 0)
                printf("[RANK %d] pdc_region_prefetch_by_objid with offset, object id %d\n",
                       pdc_client_mpi_rank_g, obj_prefetch_list[i]);
            fflush(stdout);

            is_cached = pdc_region_dl_global_search(obj_prefetch_list[i], offset, size);
        }
        else {
            if (pdc_client_mpi_rank_g == 0)
                printf("[RANK %d] pdc_region_prefetch_by_objid without offset\n", pdc_client_mpi_rank_g);
            is_cached = pdc_region_dl_global_search(obj_prefetch_list[i], NULL, NULL);
        }
    }

    ret_value = pdc_region_global_metadata_free();

    free(offset);
    free(size);
    free(obj_prefetch_list);

    if (reg_offset_list != NULL)
        free(reg_offset_list);

done:
    fflush(stdout);
    FUNC_LEAVE(ret_value);
}
