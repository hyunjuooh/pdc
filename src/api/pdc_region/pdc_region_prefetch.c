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
char **   obj_prefetch_list;
uint64_t *reg_offset_list;
uint64_t *reg_size_list;

int reg_dim;
int obj_prefetch_list_len;

// Global prefetch list
char **   global_obj_prefetch_list;
int *     global_list_len;
uint64_t *global_offset_list;
uint64_t *global_size_list;
int       global_list_size;

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

    uint64_t *ptr, *ptr2;
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

        // reg_offset_list = (uint64_t *)PDC_malloc(sizeof(uint64_t) * reg_dim * obj_prefetch_list_len);
        // reg_size_list   = (uint64_t *)PDC_malloc(sizeof(uint64_t) * reg_dim * obj_prefetch_list_len);
        reg_offset_list = (uint64_t *)PDC_malloc(sizeof(uint64_t) * reg_dim);
        reg_size_list   = (uint64_t *)PDC_malloc(sizeof(uint64_t) * reg_dim);

        ptr  = reg_offset_list;
        ptr2 = reg_size_list;
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
            memcpy(reg_offset_list, reg2->offset, sizeof(uint64_t) * reg_dim);
            // memcpy(ptr, reg2->offset, sizeof(uint64_t) * reg_dim);
            // ptr += reg_dim;

            memcpy(reg_size_list, reg2->size, sizeof(uint64_t) * reg_dim);
            // memcpy(ptr2, reg2->size, sizeof(uint64_t) * reg_dim);
            // ptr2 += reg_dim;
        }
    }

    // PDCregion_print_prefetch_list();

    // printf("[RANK %d] PDCregion_receive_prefetch_hint: pid=%d, var_a=%d, &var_a=%p, "
    //        "PDCregion_print_prefetch_list=%p\n",
    //        pdc_client_mpi_rank_g, getpid(), obj_prefetch_list_len, (void *)&obj_prefetch_list_len,
    //        (void *)PDCregion_print_prefetch_list);

    pdc_region_cache_timelog(start, "PDCregion_receive_prefetch_hint - Total time");

done:
    fflush(stdout);
    FUNC_LEAVE(ret_value);
}

void
flatten_strings(char **src, int count, char **flat_buf, int **lengths, int *total_bytes)
{
    *lengths     = malloc(count * sizeof(int));
    *total_bytes = 0;
    for (int i = 0; i < count; i++) {
        (*lengths)[i] = strlen(src[i]) + 1; // include null terminator
        *total_bytes += (*lengths)[i];
    }

    *flat_buf  = malloc(*total_bytes);
    int offset = 0;
    for (int i = 0; i < count; i++) {
        memcpy(*flat_buf + offset, src[i], (*lengths)[i]);
        offset += (*lengths)[i];
    }
}

perr_t
pdc_region_prepare_global_prefetch_list()
{
    perr_t ret_value = SUCCEED;
    double start;

    FUNC_ENTER(NULL);

    start = MPI_Wtime();
    printf("Rank %d entered global_prefetch_list creation\n", pdc_client_mpi_rank_g);
    fflush(stdout);
    
    // Gather how many strings each rank has
    global_list_len = (int *)PDC_malloc(pdc_client_mpi_size_g * sizeof(int));
    MPI_Allgather(&obj_prefetch_list_len, 1, MPI_INT, global_list_len, 1, MPI_INT, MPI_COMM_WORLD);

    printf("Rank %d Gather how many strings each rank has\n", pdc_client_mpi_rank_g);
    fflush(stdout);

    // Flatten local strings and gather lengths
    char *flat_data;
    int * lengths, total_bytes;
    flatten_strings(obj_prefetch_list, obj_prefetch_list_len, &flat_data, &lengths, &total_bytes);

    printf("Rank %d Flatten local strings and gather lengths\n", pdc_client_mpi_rank_g);
    fflush(stdout);

    // Gather string lengths from all ranks
    int  total_strings       = 0;
    int *recv_lengths        = NULL;
    int *recv_lengths_displs = (int *)PDC_malloc(pdc_client_mpi_size_g * sizeof(int));
    int *recv_lengths_counts = (int *)PDC_malloc(pdc_client_mpi_size_g * sizeof(int));
    recv_lengths_displs[0]   = 0;
    for (int i = 0; i < pdc_client_mpi_size_g; i++) {
        recv_lengths_counts[i] = global_list_len[i];
        if (i > 0)
            recv_lengths_displs[i] = recv_lengths_displs[i - 1] + global_list_len[i - 1];
        total_strings += global_list_len[i];
    }
    recv_lengths = (int *)PDC_malloc(total_strings * sizeof(int));
    MPI_Allgatherv(lengths, obj_prefetch_list_len, MPI_INT, recv_lengths, recv_lengths_counts,
                   recv_lengths_displs, MPI_INT, MPI_COMM_WORLD);

    printf("Rank %d Gather how many strings each rank has\n", pdc_client_mpi_rank_g);
    fflush(stdout);

    // Gather flattened data
    int *local_bytes = (int *)PDC_malloc(sizeof(int));
    *local_bytes     = total_bytes;
    int *recv_bytes  = (int *)PDC_malloc(pdc_client_mpi_size_g * sizeof(int));
    MPI_Allgather(local_bytes, 1, MPI_INT, recv_bytes, 1, MPI_INT, MPI_COMM_WORLD);

    int *recv_bytes_displs = (int *)PDC_malloc(pdc_client_mpi_size_g * sizeof(int));
    recv_bytes_displs[0]   = 0;
    for (int i = 1; i < pdc_client_mpi_size_g; i++) {
        recv_bytes_displs[i] = recv_bytes_displs[i - 1] + recv_bytes[i - 1];
    }

    int total_recv_bytes =
        recv_bytes_displs[pdc_client_mpi_size_g - 1] + recv_bytes[pdc_client_mpi_size_g - 1];
    char *global_flat = (char *)PDC_malloc(total_recv_bytes);

    MPI_Allgatherv(flat_data, total_bytes, MPI_CHAR, global_flat, recv_bytes, recv_bytes_displs, MPI_CHAR,
                   MPI_COMM_WORLD);

    printf("Rank %d Gather flattened data\n", pdc_client_mpi_rank_g);
    fflush(stdout);

    // Reconstruct global char** list
    global_obj_prefetch_list = (char *)PDC_malloc(total_strings * sizeof(char *));
    int offset               = 0;
    for (int i = 0; i < total_strings; i++) {
        global_obj_prefetch_list[i] = strdup(global_flat + offset);
        offset += recv_lengths[i];
    }

    global_list_size = total_strings;

    printf("Rank %d Reconstruct global char** list\n", pdc_client_mpi_rank_g);
    fflush(stdout);

    if (reg_offset_list != NULL) {
        global_offset_list = (uint64_t *)PDC_malloc(reg_dim * pdc_client_mpi_size_g * sizeof(uint64_t));
        global_size_list   = (uint64_t *)PDC_malloc(reg_dim * pdc_client_mpi_size_g * sizeof(uint64_t));

        MPI_Allgather(reg_offset_list, reg_dim, MPI_UINT64_T, global_offset_list, reg_dim, MPI_UINT64_T,
                      MPI_COMM_WORLD);
        MPI_Allgather(reg_size_list, reg_dim, MPI_UINT64_T, global_size_list, reg_dim, MPI_UINT64_T,
                      MPI_COMM_WORLD);
    }
    else {
        global_offset_list = NULL;
        global_size_list   = NULL;
        printf("Rank %d No region offset list \n", pdc_client_mpi_rank_g);
        fflush(stdout);
    }

    // Print gathered result at each rank
    // if (pdc_client_mpi_rank_g == 0) {
    //     printf("Rank %d received:\n", pdc_client_mpi_rank_g);
    //     // for (int i = 0; i < total_strings; i++) {
    //     //     printf("  [%d] %s\n", i, global_obj_prefetch_list[i]);
    //     // }

    //     // for (int i = 0; i < pdc_client_mpi_size_g; i++) {
    //     //     printf("  [%d] %d\n", i, global_list_len[i]);
    //     // }

    //     printf("  From rank %d: list1 =", pdc_client_mpi_rank_g);
    //     for (int i = 0; i < pdc_client_mpi_size_g; i++) {
    //         printf(" %" PRIu64, global_offset_list[i]);
    //     }
    //     printf(" | list2 =");
    //     for (int i = 0; i < pdc_client_mpi_size_g; i++) {
    //         printf(" %" PRIu64, global_size_list[i]);
    //     }
    //     printf("\n");
    //     fflush(stdout);
    // }

    free(flat_data);
    free(lengths);
    free(local_bytes);
    free(recv_lengths);
    free(recv_lengths_displs);
    free(recv_lengths_counts);
    free(recv_bytes);
    free(recv_bytes_displs);
    free(global_flat);

done:
    fflush(stdout);
    FUNC_LEAVE(ret_value);
}

perr_t
PDCregion_prefetch_by_objid()
{
    perr_t    ret_value = SUCCEED;
    uint64_t *offset, *size, *ptr, *ptr2;
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
    ret_value = pdc_region_dl_prepare_data_exchange(global_obj_prefetch_list, global_offset_list,
                                                    global_size_list, global_list_len);

    ret_value = pdc_region_dl_data_exchange(obj_prefetch_list_len);

    for (i = 0; i < obj_prefetch_list_len; i++) {
        free(obj_prefetch_list[i]);
    }

    for (i = 0; i < global_list_size; i++) {
        free(global_obj_prefetch_list[i]);
    }

    free(obj_prefetch_list);
    free(global_obj_prefetch_list);
    free(global_offset_list);
    free(global_size_list);
    free(reg_offset_list);
    free(reg_size_list);

    obj_prefetch_list        = NULL;
    global_obj_prefetch_list = NULL;
    global_offset_list       = NULL;
    global_size_list         = NULL;
    reg_offset_list          = NULL;
    reg_size_list            = NULL;

    pdc_region_cache_timelog(start, "PDCregion_prefetch_by_objid - Total time");

done:
    fflush(stdout);
    FUNC_LEAVE(ret_value);
}