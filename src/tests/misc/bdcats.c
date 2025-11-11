/*
 * Copyright Notice for
 * Proactive Data Containers (PDC) Software Library and Utilities
 * -----------------------------------------------------------------------------

 *** Copyright Notice ***

 * Proactive Data Containers (PDC) Copyright (c) 2017, The Regents of the
 * University of California, through Lawrence Berkeley National Laboratory,
 * UChicago Argonne, LLC, operator of Argonne National Laboratory, and The HDF
 * Group (subject to receipt of any required approvals from the U.S. Dept. of
 * Energy).  All rights reserved.

 * If you have questions about your rights to use or distribute this software,
 * please contact Berkeley Lab's Innovation & Partnerships Office at  IPO@lbl.gov.

 * NOTICE.  This Software was developed under funding from the U.S. Department of
 * Energy and the U.S. Government consequently retains certain rights. As such, the
 * U.S. Government has been granted for itself and others acting on its behalf a
 * paid-up, nonexclusive, irrevocable, worldwide license in the Software to
 * reproduce, distribute copies to the public, prepare derivative works, and
 * perform publicly and display publicly, and to permit other to do so.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <getopt.h>
#include <time.h>
#include <sys/time.h>
#include <math.h>
#include <inttypes.h>
#include "pdc.h"
#include "test_helper.h"

#define NPARTICLES 8388608

double
uniform_random_number()
{
    return (((double)rand()) / ((double)(RAND_MAX)));
}

void
print_usage()
{
    LOG_JUST_PRINT("Usage: srun -n ./bdcats #particles\n");
}

int
main(int argc, char **argv)
{
    int     rank = 0, size = 1;
    pdcid_t pdc_id, cont_id;
    pdcid_t obj_xx, obj_yy, obj_zz, obj_pxx, obj_pyy, obj_pzz, obj_id11, obj_id22;
    pdcid_t region_x, region_y, region_z, region_px, region_py, region_pz, region_id1, region_id2;
    pdcid_t region_xx, region_yy, region_zz, region_pxx, region_pyy, region_pzz, region_id11, region_id22;
    int     ret_value = TSUCCEED;

    float *   x, *y, *z;
    float *   px, *py, *pz;
    int *     id1, *id2;
    uint64_t  numparticles;
    int       ndim = 1;
    uint64_t *offset;
    uint64_t *offset_remote;
    uint64_t *mysize;

    pdcid_t transfer_request_x, transfer_request_y, transfer_request_z, transfer_request_px,
        transfer_request_py, transfer_request_pz, transfer_request_id1, transfer_request_id2;

#ifdef ENABLE_MPI
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
#endif

    numparticles = NPARTICLES;
    if (argc == 2) {
        numparticles = atoll(argv[1]);
        if (rank == 0)
            LOG_INFO("Writing %" PRIu64 " number of particles with %d clients.\n", numparticles, size);
    }

    x = (float *)malloc(numparticles * sizeof(float));
    y = (float *)malloc(numparticles * sizeof(float));
    z = (float *)malloc(numparticles * sizeof(float));

    px = (float *)malloc(numparticles * sizeof(float));
    py = (float *)malloc(numparticles * sizeof(float));
    pz = (float *)malloc(numparticles * sizeof(float));

    id1 = (int *)malloc(numparticles * sizeof(int));
    id2 = (int *)malloc(numparticles * sizeof(int));

    // create a pdc
    TASSERT((pdc_id = PDCinit("pdc")) != 0, "Call to PDCinit succeeded", "Call to PDCinit failed");
    // open a container
    TASSERT((cont_id = PDCcont_open("c1", pdc_id)) > 0, "Call to PDCcont_open succeeded",
            "Call to PDCcont_open failed");
    // open objects
    TASSERT((obj_xx = PDCobj_open("obj-var-xx", pdc_id)) != 0, "Call to PDCobj_open succeeded",
            "Call to PDCobj_open failed");
    TASSERT((obj_yy = PDCobj_open("obj-var-yy", pdc_id)) != 0, "Call to PDCobj_open succeeded",
            "Call to PDCobj_open failed");
    TASSERT((obj_zz = PDCobj_open("obj-var-zz", pdc_id)) != 0, "Call to PDCobj_open succeeded",
            "Call to PDCobj_open failed");
    TASSERT((obj_pxx = PDCobj_open("obj-var-pxx", pdc_id)) != 0, "Call to PDCobj_open succeeded",
            "Call to PDCobj_open failed");
    TASSERT((obj_pyy = PDCobj_open("obj-var-pyy", pdc_id)) != 0, "Call to PDCobj_open succeeded",
            "Call to PDCobj_open failed");
    TASSERT((obj_pzz = PDCobj_open("obj-var-pyy", pdc_id)) != 0, "Call to PDCobj_open succeeded",
            "Call to PDCobj_open failed");
    TASSERT((obj_id11 = PDCobj_open("obj-var-pyy", pdc_id)) != 0, "Call to PDCobj_open succeeded",
            "Call to PDCobj_open failed");
    TASSERT((obj_id22 = PDCobj_open("id22", pdc_id)) != 0, "Call to PDCobj_open succeeded",
            "Call to PDCobj_open failed");

    offset           = (uint64_t *)malloc(sizeof(uint64_t) * ndim);
    offset_remote    = (uint64_t *)malloc(sizeof(uint64_t) * ndim);
    mysize           = (uint64_t *)malloc(sizeof(uint64_t) * ndim);
    offset[0]        = 0;
    offset_remote[0] = rank * numparticles;
    mysize[0]        = numparticles;

    // create regions
    TASSERT((region_x = PDCregion_create(ndim, offset, mysize)) != 0, "Call to PDCregion_create succeeded",
            "Call to PDCregion_create failed");
    TASSERT((region_y = PDCregion_create(ndim, offset, mysize)) != 0, "Call to PDCregion_create succeeded",
            "Call to PDCregion_create failed");
    TASSERT((region_z = PDCregion_create(ndim, offset, mysize)) != 0, "Call to PDCregion_create succeeded",
            "Call to PDCregion_create failed");
    TASSERT((region_px = PDCregion_create(ndim, offset, mysize)) != 0, "Call to PDCregion_create succeeded",
            "Call to PDCregion_create failed");
    TASSERT((region_py = PDCregion_create(ndim, offset, mysize)) != 0, "Call to PDCregion_create succeeded",
            "Call to PDCregion_create failed");
    TASSERT((region_pz = PDCregion_create(ndim, offset, mysize)) != 0, "Call to PDCregion_create succeeded",
            "Call to PDCregion_create failed");
    TASSERT((region_id1 = PDCregion_create(ndim, offset, mysize)) != 0, "Call to PDCregion_create succeeded",
            "Call to PDCregion_create failed");
    TASSERT((region_id2 = PDCregion_create(ndim, offset, mysize)) != 0, "Call to PDCregion_create succeeded",
            "Call to PDCregion_create failed");

    TASSERT((region_xx = PDCregion_create(ndim, offset_remote, mysize)) != 0, "Call to PDCregion_create succeeded",
            "Call to PDCregion_create failed");
    TASSERT((region_yy = PDCregion_create(ndim, offset_remote, mysize)) != 0, "Call to PDCregion_create succeeded",
            "Call to PDCregion_create failed");
    TASSERT((region_zz = PDCregion_create(ndim, offset_remote, mysize)) != 0, "Call to PDCregion_create succeeded",
            "Call to PDCregion_create failed");
    TASSERT((region_pxx = PDCregion_create(ndim, offset_remote, mysize)) != 0, "Call to PDCregion_create succeeded",
            "Call to PDCregion_create failed");
    TASSERT((region_pyy = PDCregion_create(ndim, offset_remote, mysize)) != 0, "Call to PDCregion_create succeeded",
            "Call to PDCregion_create failed");
    TASSERT((region_pzz = PDCregion_create(ndim, offset_remote, mysize)) != 0, "Call to PDCregion_create succeeded",
            "Call to PDCregion_create failed");
    TASSERT((region_id11 = PDCregion_create(ndim, offset_remote, mysize)) != 0, "Call to PDCregion_create succeeded",
            "Call to PDCregion_create failed");
    TASSERT((region_id22 = PDCregion_create(ndim, offset_remote, mysize)) != 0, "Call to PDCregion_create succeeded",
            "Call to PDCregion_create failed");

#ifdef ENABLE_MPI
    MPI_Barrier(MPI_COMM_WORLD);
#endif

    TASSERT((transfer_request_x = PDCregion_transfer_create(&x[0], PDC_READ, obj_xx, region_x, region_xx)) !=
                0,
            "Call to PDCregion_transfer_create succeeded", "Call to PDCregion_transfer_create failed");
    TASSERT((transfer_request_y = PDCregion_transfer_create(&y[0], PDC_READ, obj_yy, region_y, region_yy)) !=
                0,
            "Call to PDCregion_transfer_create succeeded", "Call to PDCregion_transfer_create failed");
    TASSERT((transfer_request_z = PDCregion_transfer_create(&z[0], PDC_READ, obj_zz, region_z, region_zz)) !=
                0,
            "Call to PDCregion_transfer_create succeeded", "Call to PDCregion_transfer_create failed");
    TASSERT((transfer_request_px =
                 PDCregion_transfer_create(&px[0], PDC_READ, obj_pxx, region_px, region_pxx)) != 0,
            "Call to PDCregion_transfer_create succeeded", "Call to PDCregion_transfer_create failed");
    TASSERT((transfer_request_py =
                 PDCregion_transfer_create(&py[0], PDC_READ, obj_pyy, region_py, region_pyy)) != 0,
            "Call to PDCregion_transfer_create succeeded", "Call to PDCregion_transfer_create failed");
    TASSERT((transfer_request_pz =
                 PDCregion_transfer_create(&pz[0], PDC_READ, obj_pzz, region_pz, region_pzz)) != 0,
            "Call to PDCregion_transfer_create succeeded", "Call to PDCregion_transfer_create failed");
    TASSERT((transfer_request_id1 =
                 PDCregion_transfer_create(&id1[0], PDC_READ, obj_id11, region_id1, region_id11)) != 0,
            "Call to PDCregion_transfer_create succeeded", "Call to PDCregion_transfer_create failed");
    TASSERT((transfer_request_id2 =
                 PDCregion_transfer_create(&id2[0], PDC_READ, obj_id22, region_id2, region_id22)) != 0,
            "Call to PDCregion_transfer_create succeeded", "Call to PDCregion_transfer_create failed");

    TASSERT(PDCregion_transfer_start(transfer_request_x) >= 0, "Call to PDCregion_transfer_start succeeded",
            "Call to PDCregion_transfer_start failed");
    TASSERT(PDCregion_transfer_start(transfer_request_y) >= 0, "Call to PDCregion_transfer_start succeeded",
            "Call to PDCregion_transfer_start failed");
    TASSERT(PDCregion_transfer_start(transfer_request_z) >= 0, "Call to PDCregion_transfer_start succeeded",
            "Call to PDCregion_transfer_start failed");
    TASSERT(PDCregion_transfer_start(transfer_request_px) >= 0, "Call to PDCregion_transfer_start succeeded",
            "Call to PDCregion_transfer_start failed");
    TASSERT(PDCregion_transfer_start(transfer_request_py) >= 0, "Call to PDCregion_transfer_start succeeded",
            "Call to PDCregion_transfer_start failed");
    TASSERT(PDCregion_transfer_start(transfer_request_pz) >= 0, "Call to PDCregion_transfer_start succeeded",
            "Call to PDCregion_transfer_start failed");
    TASSERT(PDCregion_transfer_start(transfer_request_id1) >= 0, "Call to PDCregion_transfer_start succeeded",
            "Call to PDCregion_transfer_start failed");
    TASSERT(PDCregion_transfer_start(transfer_request_id2) >= 0, "Call to PDCregion_transfer_start succeeded",
            "Call to PDCregion_transfer_start failed");

    TASSERT(PDCregion_transfer_wait(transfer_request_x) >= 0, "Call to PDCregion_transfer_wait succeeded",
            "Call to PDCregion_transfer_wait failed");
    TASSERT(PDCregion_transfer_wait(transfer_request_y) >= 0, "Call to PDCregion_transfer_wait succeeded",
            "Call to PDCregion_transfer_wait failed");
    TASSERT(PDCregion_transfer_wait(transfer_request_z) >= 0, "Call to PDCregion_transfer_wait succeeded",
            "Call to PDCregion_transfer_wait failed");
    TASSERT(PDCregion_transfer_wait(transfer_request_px) >= 0, "Call to PDCregion_transfer_wait succeeded",
            "Call to PDCregion_transfer_wait failed");
    TASSERT(PDCregion_transfer_wait(transfer_request_py) >= 0, "Call to PDCregion_transfer_wait succeeded",
            "Call to PDCregion_transfer_wait failed");
    TASSERT(PDCregion_transfer_wait(transfer_request_pz) >= 0, "Call to PDCregion_transfer_wait succeeded",
            "Call to PDCregion_transfer_wait failed");
    TASSERT(PDCregion_transfer_wait(transfer_request_id1) >= 0, "Call to PDCregion_transfer_wait succeeded",
            "Call to PDCregion_transfer_wait failed");
    TASSERT(PDCregion_transfer_wait(transfer_request_id2) >= 0, "Call to PDCregion_transfer_wait succeeded",
            "Call to PDCregion_transfer_wait failed");

    TASSERT(PDCregion_transfer_close(transfer_request_x) >= 0, "Call to PDCregion_transfer_close succeeded",
            "Call to PDCregion_transfer_close failed");
    TASSERT(PDCregion_transfer_close(transfer_request_y) >= 0, "Call to PDCregion_transfer_close succeeded",
            "Call to PDCregion_transfer_close failed");
    TASSERT(PDCregion_transfer_close(transfer_request_z) >= 0, "Call to PDCregion_transfer_close succeeded",
            "Call to PDCregion_transfer_close failed");
    TASSERT(PDCregion_transfer_close(transfer_request_px) >= 0, "Call to PDCregion_transfer_close succeeded",
            "Call to PDCregion_transfer_close failed");
    TASSERT(PDCregion_transfer_close(transfer_request_py) >= 0, "Call to PDCregion_transfer_close succeeded",
            "Call to PDCregion_transfer_close failed");
    TASSERT(PDCregion_transfer_close(transfer_request_pz) >= 0, "Call to PDCregion_transfer_close succeeded",
            "Call to PDCregion_transfer_close failed");
    TASSERT(PDCregion_transfer_close(transfer_request_id1) >= 0, "Call to PDCregion_transfer_close succeeded",
            "Call to PDCregion_transfer_close failed");
    TASSERT(PDCregion_transfer_close(transfer_request_id2) >= 0, "Call to PDCregion_transfer_close succeeded",
            "Call to PDCregion_transfer_close failed");

    PDC_timing_report("read");

    TASSERT(PDCobj_close(obj_xx) >= 0, "Call to PDCobj_close succeeded", "Call to PDCobj_close failed");
    TASSERT(PDCobj_close(obj_yy) >= 0, "Call to PDCobj_close succeeded", "Call to PDCobj_close failed");
    TASSERT(PDCobj_close(obj_zz) >= 0, "Call to PDCobj_close succeeded", "Call to PDCobj_close failed");
    TASSERT(PDCobj_close(obj_pxx) >= 0, "Call to PDCobj_close succeeded", "Call to PDCobj_close failed");
    TASSERT(PDCobj_close(obj_pyy) >= 0, "Call to PDCobj_close succeeded", "Call to PDCobj_close failed");
    TASSERT(PDCobj_close(obj_pzz) >= 0, "Call to PDCobj_close succeeded", "Call to PDCobj_close failed");
    TASSERT(PDCobj_close(obj_id11) >= 0, "Call to PDCobj_close succeeded", "Call to PDCobj_close failed");
    TASSERT(PDCobj_close(obj_id22) >= 0, "Call to PDCobj_close succeeded", "Call to PDCobj_close failed");
    TASSERT(PDCregion_close(region_x) >= 0, "Call to PDCregion_close succeeded",
            "Call to PDCregion_close failed");
    TASSERT(PDCregion_close(region_y) >= 0, "Call to PDCregion_close succeeded",
            "Call to PDCregion_close failed");
    TASSERT(PDCregion_close(region_z) >= 0, "Call to PDCregion_close succeeded",
            "Call to PDCregion_close failed");
    TASSERT(PDCregion_close(region_px) >= 0, "Call to PDCregion_close succeeded",
            "Call to PDCregion_close failed");
    TASSERT(PDCregion_close(region_py) >= 0, "Call to PDCregion_close succeeded",
            "Call to PDCregion_close failed");
    TASSERT(PDCregion_close(region_pz) >= 0, "Call to PDCregion_close succeeded",
            "Call to PDCregion_close failed");
    TASSERT(PDCregion_close(region_id1) >= 0, "Call to PDCregion_close succeeded",
            "Call to PDCregion_close failed");
    TASSERT(PDCregion_close(region_id2) >= 0, "Call to PDCregion_close succeeded",
            "Call to PDCregion_close failed");
    TASSERT(PDCregion_close(region_xx) >= 0, "Call to PDCregion_close succeeded",
            "Call to PDCregion_close failed");
    TASSERT(PDCregion_close(region_yy) >= 0, "Call to PDCregion_close succeeded",
            "Call to PDCregion_close failed");
    TASSERT(PDCregion_close(region_zz) >= 0, "Call to PDCregion_close succeeded",
            "Call to PDCregion_close failed");
    TASSERT(PDCregion_close(region_pxx) >= 0, "Call to PDCregion_close succeeded",
            "Call to PDCregion_close failed");
    TASSERT(PDCregion_close(region_pyy) >= 0, "Call to PDCregion_close succeeded",
            "Call to PDCregion_close failed");
    TASSERT(PDCregion_close(region_pzz) >= 0, "Call to PDCregion_close succeeded",
            "Call to PDCregion_close failed");
    TASSERT(PDCregion_close(region_id11) >= 0, "Call to PDCregion_close succeeded",
            "Call to PDCregion_close failed");
    TASSERT(PDCregion_close(region_id22) >= 0, "Call to PDCregion_close succeeded",
            "Call to PDCregion_close failed");
    // close a container
    TASSERT(PDCcont_close(cont_id) >= 0, "Call to PDCcont_close succeeded", "Call to PDCcont_close failed");
    TASSERT(PDCclose(pdc_id) >= 0, "Call to PDCclose succeeded", "Call to PDCclose failed");

    free(offset);
    free(offset_remote);
    free(mysize);
    free(x);
    free(y);
    free(z);
    free(px);
    free(py);
    free(pz);
    free(id1);
    free(id2);

done:
#ifdef ENABLE_MPI
    MPI_Finalize();
#endif

    return ret_value;
}
