#include <time.h>
#include <stdlib.h>
#include <string.h>
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
#include "pdc_region_prefetch.h"
#include "pdc_region_cache_dl.h"
#include "pdc_client_connect.h"

// pdcid_t * obj_prefetch_list = NULL;
char **obj_prefetch_list;
char **global_prefetch_list;
int *  global_list_len;
int *  global_list_item_len;

uint64_t *reg_offset_list;

int reg_dim;
int obj_prefetch_list_len;

perr_t
pdc_region_prefetch_init()
{
    perr_t ret_value = SUCCEED;

    FUNC_ENTER(NULL);

    obj_prefetch_list     = NULL;
    reg_offset_list       = NULL;
    reg_dim               = 1;
    obj_prefetch_list_len = 0;

done:
    fflush(stdout);
    FUNC_LEAVE(ret_value);
}

perr_t
PDCregion_print_prefetch_list()
{
    perr_t ret_value = SUCCEED;
    int    i;

    FUNC_ENTER(NULL);

    printf("[RANK %d] Prefetch list item number %d\n ", pdc_client_mpi_rank_g, obj_prefetch_list_len);

done:
    fflush(stdout);
    FUNC_LEAVE(ret_value);
}

perr_t
PDCregion_receive_prefetch_hint(char *arr[], pdcid_t *arr2, int obj_array_len)
{
    perr_t ret_value = SUCCEED;

    struct _pdc_id_info *   reginfo2;
    struct pdc_region_info *reg2;

    uint64_t *ptr;
    int       i;
    double    start;

    FUNC_ENTER(NULL);

    start = MPI_Wtime();

    obj_prefetch_list_len = obj_array_len;
    obj_prefetch_list     = (char *)PDC_malloc(obj_prefetch_list_len * sizeof(char *));

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
        // obj_prefetch_list[i] = arr[i];
        obj_prefetch_list[i] = strdup(arr[i]);
        if (obj_prefetch_list[i] == NULL) {
            PGOTO_ERROR(FAIL, "Fail to copy object name");
        }

        if (arr2 != NULL) {
            reginfo2 = PDC_find_id(arr2[i]);
            if (reginfo2 == NULL)
                PGOTO_ERROR(FAIL, "cannot locate remote region ID");

            reg2 = (struct pdc_region_info *)(reginfo2->obj_ptr);
            memcpy(ptr, reg2->offset, sizeof(uint64_t) * reg_dim);
            ptr += reg_dim;

            memcpy(ptr, reg2->size, sizeof(uint64_t) * reg_dim);
            ptr += reg_dim;
        }
    }

    // PDCregion_print_prefetch_list();

    pdc_region_cache_timelog(start, "PDCregion_receive_prefetch_hint - Total time");

done:
    fflush(stdout);
    FUNC_LEAVE(ret_value);
}

perr_t
pdc_region_prepare_global_prefetch_list()
{
    perr_t ret_value = SUCCEED;
    int *  displs, *list_item_len, *global_bytes;
    int    local_bytes = 0, pos = 0, total_item_num = 0, total_bytes = 0;
    char * sendlist, *recvlist;
    int    i;
    double start;

    FUNC_ENTER(NULL);

    // Collect all list length globally
    global_list_len = (int *)PDC_malloc(pdc_client_mpi_size_g * sizeof(int));
    MPI_Allgather(&obj_prefetch_list_len, 1, MPI_INT, global_list_len, 1, MPI_INT, MPI_COMM_WORLD);

    displs    = (int *)PDC_malloc(pdc_client_mpi_size_g * sizeof(int));
    displs[0] = 0;
    for (i = 0; i < pdc_client_mpi_size_g; i++) {
        total_item_num += global_list_len[i];
        if (i > 0) {
            displs[i] = displs[i - 1] + global_list_len[i - 1];
        }
    }

    list_item_len = (int *)PDC_malloc(obj_prefetch_list_len * sizeof(int));
    for (int i = 0; i < obj_prefetch_list_len; i++) {
        // Include null terminator for string
        list_item_len[i] = strlen(obj_prefetch_list[i]) + 1;
    }

    global_list_item_len = (int *)PDC_malloc(total_item_num * sizeof(int));
    MPI_Allgatherv(list_item_len, total_item_num, MPI_INT, global_list_item_len, global_list_len, displs,
                   MPI_INT, MPI_COMM_WORLD);

    // Flatten prefetch obj name list for collection
    for (i = 0; i < obj_prefetch_list_len; i++) {
        local_bytes += list_item_len[i];
    }

    sendlist = (char *)PDC_malloc(local_bytes);
    for (i = 0; i < obj_prefetch_list_len; i++) {
        memcpy(sendlist + pos, obj_prefetch_list[i], list_item_len[i]);
        pos += list_item_len[i];
    }

    // Prepare data size and displacements for collection of flatten prefetch list
    global_bytes = (int *)PDC_malloc(pdc_client_mpi_size_g * sizeof(int));
    MPI_Allgather(&local_bytes, 1, MPI_INT, global_bytes, 1, MPI_INT, MPI_COMM_WORLD);

    displs[0] = 0;
    for (i = 0; i < pdc_client_mpi_size_g; i++) {
        total_bytes += global_bytes[i];
        if (i > 0) {
            displs[i] = displs[i - 1] + global_bytes[i - 1];
        }
    }

    // Collect flatten prefetch list
    recvlist = (char *)PDC_malloc(total_bytes);
    MPI_Allgatherv(sendlist, local_bytes, MPI_CHAR, recvlist, total_bytes, displs, MPI_CHAR, MPI_COMM_WORLD);

    // Create a global list for each rank
    global_prefetch_list = (char *)PDC_malloc(total_item_num * sizeof(char *));

    pos = 0;
    for (i = 0; i < total_item_num; i++) {
        global_prefetch_list[i] = (char *)PDC_malloc(global_list_item_len[i]);
        memcpy(global_prefetch_list[i], recvlist + pos, global_list_item_len[i]);
        pos += global_list_item_len[i];
    }

    free(displs);
    free(list_item_len);
    free(sendlist);
    free(global_bytes);
    free(recvlist);

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
    double    start;

    FUNC_ENTER(NULL);

    start = MPI_Wtime();

    if (obj_prefetch_list == NULL) {
        if (pdc_client_mpi_rank_g == 0)
            printf("[RANK %d] pdc_region_prefetch_by_objid - object list not created\n",
                   pdc_client_mpi_rank_g);

        goto done;
    }

    ret_value = pdc_region_prepare_global_prefetch_list();
    ret_value =
        pdc_region_dl_prepare_data_exchange(global_prefetch_list, global_list_len, global_list_item_len);
    ret_value = pdc_region_dl_global_data_exchange(obj_prefetch_list_len);

    // ret_value = pdc_region_dl_list_init();

    // ptr    = reg_offset_list;
    // offset = (uint64_t *)PDC_malloc(sizeof(uint64_t) * reg_dim);
    // size   = (uint64_t *)PDC_malloc(sizeof(uint64_t) * reg_dim);

    // for (i = 0; i < obj_prefetch_list_len; i++) {
    //     if (ptr != NULL) {
    //         memcpy(offset, ptr, sizeof(uint64_t) * reg_dim);
    //         ptr += reg_dim;

    //         memcpy(size, ptr, sizeof(uint64_t) * reg_dim);
    //         ptr += reg_dim;

    //         is_cached = pdc_region_dl_global_search(obj_prefetch_list[i], offset, size);
    //     }
    //     else {
    //         is_cached = pdc_region_dl_global_search(obj_prefetch_list[i], NULL, NULL);
    //     }
    // }

    // ret_value = pdc_region_global_metadata_free();

    // free(offset);
    // free(size);

    // for (i = 0; i < obj_prefetch_list_len; i++) {
    //     free(obj_prefetch_list[i]);
    // }
    // free(obj_prefetch_list);

    // obj_prefetch_list = NULL;

    // if (reg_offset_list != NULL) {
    //     free(reg_offset_list);
    //     reg_offset_list = NULL;
    // }

    pdc_region_cache_timelog(start, "PDCregion_prefetch_by_objid - Total time");

done:
    fflush(stdout);
    FUNC_LEAVE(ret_value);
}