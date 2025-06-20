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
#include "pdc_region_prefetch.h"
#include "pdc_client_connect.h"

// #define MAX_CACHE_SIZE 268435456
// #define MAX_CACHE_SIZE 4294967296
#define MAX_CACHE_SIZE 34359738368
#define MAX_ITEM_NUM   1000

size_t total_buf_size = 0;
int    total_item_num = 0;

pdcid_t pdc_id;

// Initialization of global variables
perr_t
pdc_region_cache_init(pdcid_t pdcid)
{
    perr_t ret_value = SUCCEED;
    double start;

    FUNC_ENTER(NULL);

    start = MPI_Wtime();

    total_buf_size = 0;
    pdc_id         = pdcid;

    ret_value = pdc_region_dl_init();
    ret_value = pdc_region_prefetch_init();

    pdc_region_cache_timelog(start, "pdc_region_cache_init - total time");

done:
    fflush(stdout);
    FUNC_LEAVE(ret_value);
}

// TODO:
// Currently considering fully contained case
int
pdc_region_cache_search(char *obj_name, int ndim, uint64_t unit, uint64_t *offset, uint64_t *size, void *buf)
{
    int      region_contained = 0;
    uint64_t read_size        = 0;
    double   start;

    start = MPI_Wtime();

    // Calculation of read size
    if (ndim >= 1)
        read_size = unit * size[0];
    if (ndim >= 2)
        read_size *= size[1];
    if (ndim >= 3)
        read_size *= size[2];

    // Search on doubly linked list
    region_contained = pdc_region_dl_search(obj_name, ndim, unit, offset, size, buf, read_size);

    pdc_region_cache_timelog(start, "pdc_region_cache_search - total time");

    return region_contained;
}

// Insert the region to the list
perr_t
pdc_region_cache_insert(int data_exchange, char *obj_name, int ndim, uint64_t unit, uint64_t *offset,
                        uint64_t *size, void *buf)
{
    perr_t ret_value = SUCCEED;

    uint64_t read_size = 0;
    double   start;

    FUNC_ENTER(NULL);

    start = MPI_Wtime();

    // Calculation of read size
    if (ndim >= 1)
        read_size = unit * size[0];
    if (ndim >= 2)
        read_size *= size[1];
    if (ndim >= 3)
        read_size *= size[2];

    // Check if there is remaining buffer size to insert region
    // If there is no remaining capacity, free the buffer according to LRU policy
    // If it is during data_exchange skip the eviction since it deletes item internally
    if ((total_buf_size + read_size > MAX_CACHE_SIZE) && !data_exchange) {
        pdc_region_cache_evict(total_buf_size + read_size);
    }

    pdc_region_cache_timelog(start, "pdc_region_cache_insert - pdc_region_cache_evict time");

    ret_value = pdc_region_dl_insert(obj_name, ndim, unit, offset, size, buf, read_size);

    if (ret_value == SUCCEED) {
        total_buf_size += read_size;
        total_item_num += 1;
    }

    pdc_region_cache_timelog(start, "pdc_region_cache_insert - total time");

done:
    fflush(stdout);
    FUNC_LEAVE(ret_value);
}

// Check if there are overlapping parts when PDC_WRITE transfer_request is executed
// If there is evict that item since it is out of date
perr_t
pdc_region_cache_update(char *obj_name, int ndim, uint64_t unit, uint64_t *offset, uint64_t *size, void *buf)
{
    perr_t ret_value = SUCCEED;

    item_delete_info update_result;
    double           start = MPI_Wtime();

    FUNC_ENTER(NULL);

    update_result = pdc_region_dl_update(obj_name, ndim, unit, offset, size, buf);
    fflush(stdout);

    pdc_region_cache_timelog(start, "pdc_region_cache_update - total time");

    total_buf_size -= update_result.deleted_size;
    total_item_num -= update_result.deleted_item_num;

done:
    fflush(stdout);
    FUNC_LEAVE(ret_value);
}

// Evict the region cache and object cache according to LRU policy
perr_t
pdc_region_cache_evict(size_t required_size)
{
    perr_t ret_value = SUCCEED;

    item_delete_info eviction_result;
    double           start;

    start = MPI_Wtime();

    FUNC_ENTER(NULL);

    eviction_result = pdc_region_dl_evict(required_size);

    total_buf_size -= eviction_result.deleted_size;
    total_item_num -= eviction_result.deleted_item_num;

    pdc_region_cache_timelog(start, "pdc_region_cache_evict - pdc_region_dl_evict time");

done:
    fflush(stdout);
    FUNC_LEAVE(ret_value);
}

perr_t
pdc_region_cache_delete(size_t deleted_size, int deleted_item_num)
{
    perr_t ret_value = SUCCEED;

    double start;

    start = MPI_Wtime();

    FUNC_ENTER(NULL);

    total_buf_size -= deleted_size;
    total_item_num -= deleted_item_num;

    pdc_region_cache_timelog(start, "pdc_region_cache_evict - pdc_region_dl_evict time");

done:
    fflush(stdout);
    FUNC_LEAVE(ret_value);
}

perr_t
pdc_region_cache_finalize()
{
    perr_t ret_value = SUCCEED;

    FUNC_ENTER(NULL);

    ret_value = pdc_region_dl_finalize();

done:
    fflush(stdout);
    FUNC_LEAVE(ret_value);
}

void
pdc_region_cache_timelog(double start_time, const char *message)
{
    int        rank_limit = 5;
    double     end_time;
    time_t     cur_time = time(NULL);
    struct tm *log_time = localtime(&cur_time);

    end_time = MPI_Wtime();
    if (pdc_client_mpi_rank_g <= rank_limit) {
        cur_time = time(NULL);
        log_time = localtime(&cur_time);
        printf("[CACHE_LOG] [%02d:%02d:%02d] [RANK %d] [TOTAL_RANK %d] | %s : %f\n", log_time->tm_hour,
               log_time->tm_min, log_time->tm_sec, pdc_client_mpi_rank_g, pdc_client_mpi_size_g, message,
               end_time - start_time);
        fflush(stdout);
    }
}