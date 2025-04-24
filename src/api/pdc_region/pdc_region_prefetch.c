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
int *    obj_idx_array;

perr_t
PDCregion_receive_prefetch_hint(int *obj_array, int obj_array_len)
{
    perr_t ret_value = SUCCEED;
    int    i, sample_per_rank, item_idx = 0;
    char   obj_filename[2000];

    FUNC_ENTER(NULL);

    if (pdc_client_mpi_rank_g == 0) {
        printf("PDCregion_receive_prefetch_hint: ");
        fflush(stdout);
        // for (i = 0; i < obj_array_len; i++) {
        //     printf("%s, \n", filename[obj_array[i]]);
        //     snprintf(obj_filename, sizeof(obj_filename), "%s%s", filename[obj_array[i]], "-records");
        //     printf("%s, \n", obj_filename);
        //     fflush(stdout);
        // }
    }

    obj_prefetch_list_len = obj_array_len;
    obj_prefetch_list     = (pdcid_t *)PDC_malloc(obj_prefetch_list_len * sizeof(pdcid_t));

    sample_per_rank = obj_prefetch_list_len / pdc_client_mpi_size_g;

    // Get the list of object id for each rank that should be prefetched
    for (i = pdc_client_mpi_rank_g * sample_per_rank; i < (pdc_client_mpi_rank_g + 1) * sample_per_rank;
         i++) {
        snprintf(obj_filename, sizeof(obj_filename), "%s%s", filename[obj_array[i]], "-records");
        obj_prefetch_list[item_idx] = PDCobj_open(obj_filename, pdc_id);
        PDCobj_close(obj_prefetch_list[item_idx]);

        item_idx++;
    }

done:
    fflush(stdout);
    FUNC_LEAVE(ret_value);
}

perr_t
PDCregion_receive_dataset(const char **dataset_array, int obj_array_len)
{
    perr_t ret_value = SUCCEED;
    int    i;

    FUNC_ENTER(NULL);

    if (pdc_client_mpi_rank_g == 0) {
        printf("PDCregion_receive_prefetch_hint: ");
        fflush(stdout);
        for (i = 0; i < obj_array_len; i++) {
            printf("%s, \n", dataset_array[i]);
            fflush(stdout);
        }
    }

done:
    fflush(stdout);
    FUNC_LEAVE(ret_value);
}

// perr_t
// PDCregion_receive_prefetch_hint(const char **obj_array, int obj_array_len)
// {
//     perr_t ret_value = SUCCEED;
//     int    i, sample_per_rank, item_idx = 0;

//     FUNC_ENTER(NULL);

//     if (pdc_client_mpi_rank_g == 0) {
//         for (i = 0; i < obj_array_len; i++) {
//             printf("PDCregion_receive_prefetch_hint %d: %s\n", i, obj_array[i]);
//             fflush(stdout);
//         }
//     }

//     // obj_prefetch_list_len = obj_array_len;
//     // obj_prefetch_list     = (pdcid_t *)PDC_malloc(obj_prefetch_list_len * sizeof(pdcid_t));

//     // sample_per_rank = obj_prefetch_list_len / pdc_client_mpi_size_g;

//     // // Get the list of object id for each rank that should be prefetched
//     // for (i = pdc_client_mpi_rank_g * sample_per_rank; i < (pdc_client_mpi_rank_g + 1) * sample_per_rank;
//     //      i++) {
//     //     obj_prefetch_list[item_idx] = PDCobj_open(obj_array[i], pdc_id);
//     //     PDCobj_close(obj_prefetch_list[item_idx]);

//     //     item_idx++;
//     // }

// done:
//     fflush(stdout);
//     FUNC_LEAVE(ret_value);
// }

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
