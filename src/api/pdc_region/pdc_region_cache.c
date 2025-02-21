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
#include "pdc_region_cache.h"
#include "pdc_region_cache_dl.h"
#include "pdc_obj_pkg.h"
#include "pdc_interface.h"
#include "pdc_transforms_pkg.h"
#include "pdc_client_connect.h"
#include "pdc_analysis_pkg.h"
#include <mpi.h>

// Temporary defined
#define MAX_CACHE_SIZE 34359738368

static size_t total_buf_size;

// Initialization of global variables
perr_t
pdc_region_cache_init()
{
    perr_t ret_value = SUCCEED;

    FUNC_ENTER(NULL);

    total_buf_size = 0;

done:
    fflush(stdout);
    FUNC_LEAVE(ret_value);
}

// TODO:
// Implement when the region is overlapping - partially containing
// Need to manage the object's offset cache information
// Currently considering fully contained case
int
pdc_region_cache_search(pdcid_t obj_id, int ndim, uint64_t unit, uint64_t *offset, uint64_t *size, void *buf)
{
    int region_contained = 0;

    double     start, end;
    time_t     cur_time = time(NULL);
    struct tm *log_time = localtime(&cur_time);

    start = MPI_Wtime();

    // Search on doubly linked list
    region_contained = pdc_region_dl_search(obj_id, ndim, unit, offset, size, buf);

    pdc_region_cache_timelog(5, start, "pdc_region_cache_search time");

    return region_contained;
}

// Insert the region to the list
perr_t
pdc_region_cache_insert(pdcid_t obj_id, int ndim, uint64_t unit, uint64_t *offset, uint64_t *size, void *buf)
{
    perr_t   ret_value = SUCCEED;
    uint64_t read_size = 0;

    double start;

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
    if (total_buf_size + read_size > MAX_CACHE_SIZE) {
        pdc_region_cache_evict(total_buf_size + sizeof(buf));
    }

    ret_value = pdc_region_dl_insert(obj_id, ndim, unit, offset, size, buf, read_size);

    if (ret_value == SUCCEED) {
        total_buf_size += read_size;
    }

    pdc_region_cache_timelog(5, start, "pdc_region_cache_insert total time");

done:
    fflush(stdout);
    FUNC_LEAVE(ret_value);
}

// Check if there are overlapping parts when PDC_WRITE transfer_request is executed
// If there is evict that item since it is out of date
// TODO:
// Add checking if the ndim matches the region item ndim
perr_t
pdc_region_cache_update(pdcid_t obj_id, int ndim, uint64_t unit, uint64_t *offset, uint64_t *size, void *buf)
{
    perr_t ret_value = SUCCEED;
    double start;
    int    updated = 0;

    FUNC_ENTER(NULL);

    start = MPI_Wtime();

    updated = pdc_region_dl_update(obj_id, ndim, unit, offset, size, buf);

    if (updated) {
        total_buf_size -= updated;
    }

    pdc_region_cache_timelog(5, start, "pdc_region_cache_update time");

done:
    fflush(stdout);
    FUNC_LEAVE(ret_value);
}

// Evict the region cache and object cache according to LRU policy
perr_t
pdc_region_cache_evict(size_t required_size)
{
    perr_t ret_value = SUCCEED;
    size_t deleted_size;

    double start;

    time_t     cur_time = time(NULL);
    struct tm *log_time = localtime(&cur_time);

    start = MPI_Wtime();

    FUNC_ENTER(NULL);

    deleted_size = pdc_region_dl_evict(required_size);

    // total_buf_size update
    total_buf_size -= deleted_size;

    pdc_region_cache_timelog(5, start, "pdc_region_cache_evict time");

done:
    fflush(stdout);
    FUNC_LEAVE(ret_value);
}

void
pdc_region_cache_timelog(int rank_limit, double start_time, const char *message)
{
    double     end_time;
    time_t     cur_time = time(NULL);
    struct tm *log_time = localtime(&cur_time);

    end_time = MPI_Wtime();
    if (pdc_client_mpi_rank_g <= rank_limit) {
        cur_time = time(NULL);
        log_time = localtime(&cur_time);
        printf("[CACHE_LOG] [%02d:%02d:%02d] [RANK %d] | %s : %f\n", log_time->tm_hour, log_time->tm_min,
               log_time->tm_sec, pdc_client_mpi_rank_g, message, end_time - start_time);
    }
}