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

#define NPARTICLES 8388608

double
uniform_random_number()
{
    return (((double)rand()) / ((double)(RAND_MAX)));
}

void
print_usage()
{
    printf("Usage: srun -n ./bdcats #particles #transfer_request\n");
}

int
main(int argc, char *argv[])
{
    int     rank = 0, size = 1, i = 0, j = 0;
    int     num_transfer_request = 0;
    double  t0, t1, start, end;
    pdcid_t pdc_id, cont_id;
    pdcid_t obj_xx, obj_yy, obj_zz, obj_pxx, obj_pyy, obj_pzz, obj_id11, obj_id22;
    pdcid_t region_x, region_y, region_z, region_px, region_py, region_pz, region_id1, region_id2;
    pdcid_t region_xx, region_yy, region_zz, region_pxx, region_pyy, region_pzz, region_id11, region_id22;
    perr_t  ret;

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

    if (argc >= 2) {
        numparticles = atoll(argv[1]);
        if (rank == 0)
            printf("Writing %" PRIu64 " number of particles with %d clients.\n", numparticles, size);
    }
    else {
        numparticles = NPARTICLES;
    }

    num_transfer_request = atoi(argv[2]);

    x = (float *)malloc(numparticles * sizeof(float));
    y = (float *)malloc(numparticles * sizeof(float));
    z = (float *)malloc(numparticles * sizeof(float));

    px = (float *)malloc(numparticles * sizeof(float));
    py = (float *)malloc(numparticles * sizeof(float));
    pz = (float *)malloc(numparticles * sizeof(float));

    id1 = (int *)malloc(numparticles * sizeof(int));
    id2 = (int *)malloc(numparticles * sizeof(int));

    // create a pdc
    pdc_id = PDCinit("pdc");

    // open a container
    cont_id = PDCcont_open("c1", pdc_id);
    if (cont_id == 0)
        printf("Fail to create container @ line  %d!\n", __LINE__);

    // open objects
    obj_xx = PDCobj_open("obj-var-xx", pdc_id);
    if (obj_xx == 0) {
        printf("Error when open object %s\n", "obj-var-xx");
        exit(-1);
    }
    obj_yy = PDCobj_open("obj-var-yy", pdc_id);
    if (obj_yy == 0) {
        printf("Error when open object %s\n", "obj-var-yy");
        exit(-1);
    }
    obj_zz = PDCobj_open("obj-var-zz", pdc_id);
    if (obj_zz == 0) {
        printf("Error when open object %s\n", "obj-var-zz");
        exit(-1);
    }
    obj_pxx = PDCobj_open("obj-var-pxx", pdc_id);
    if (obj_pxx == 0) {
        printf("Error when open object %s\n", "obj-var-pxx");
        exit(-1);
    }
    obj_pyy = PDCobj_open("obj-var-pyy", pdc_id);
    if (obj_pyy == 0) {
        printf("Error when open object %s\n", "obj-var-pyy");
        exit(-1);
    }
    obj_pzz = PDCobj_open("obj-var-pzz", pdc_id);
    if (obj_pzz == 0) {
        printf("Error when open object %s\n", "obj-var-pzz");
        exit(-1);
    }
    obj_id11 = PDCobj_open("id11", pdc_id);
    if (obj_id11 == 0) {
        printf("Error when open object %s\n", "id11");
        exit(-1);
    }
    obj_id22 = PDCobj_open("id22", pdc_id);
    if (obj_id22 == 0) {
        printf("Error when open object %s\n", "id22");
        exit(-1);
    }

    offset           = (uint64_t *)malloc(sizeof(uint64_t) * ndim);
    offset_remote    = (uint64_t *)malloc(sizeof(uint64_t) * ndim);
    mysize           = (uint64_t *)malloc(sizeof(uint64_t) * ndim);
    offset[0]        = 0;
    offset_remote[0] = rank * numparticles;
    mysize[0]        = numparticles;

    // create a region
    region_x   = PDCregion_create(ndim, offset, mysize);
    region_y   = PDCregion_create(ndim, offset, mysize);
    region_z   = PDCregion_create(ndim, offset, mysize);
    region_px  = PDCregion_create(ndim, offset, mysize);
    region_py  = PDCregion_create(ndim, offset, mysize);
    region_pz  = PDCregion_create(ndim, offset, mysize);
    region_id1 = PDCregion_create(ndim, offset, mysize);
    region_id2 = PDCregion_create(ndim, offset, mysize);

    region_xx   = PDCregion_create(ndim, offset_remote, mysize);
    region_yy   = PDCregion_create(ndim, offset_remote, mysize);
    region_zz   = PDCregion_create(ndim, offset_remote, mysize);
    region_pxx  = PDCregion_create(ndim, offset_remote, mysize);
    region_pyy  = PDCregion_create(ndim, offset_remote, mysize);
    region_pzz  = PDCregion_create(ndim, offset_remote, mysize);
    region_id11 = PDCregion_create(ndim, offset_remote, mysize);
    region_id22 = PDCregion_create(ndim, offset_remote, mysize);

#ifdef ENABLE_MPI
    MPI_Barrier(MPI_COMM_WORLD);
    t0    = MPI_Wtime();
    start = MPI_Wtime();
#endif

    printf("PDC read whole region and cache it into client-side cache\n");

    transfer_request_x   = PDCregion_transfer_create(&x[0], PDC_READ, obj_xx, region_x, region_xx);
    transfer_request_y   = PDCregion_transfer_create(&y[0], PDC_READ, obj_yy, region_y, region_yy);
    transfer_request_z   = PDCregion_transfer_create(&z[0], PDC_READ, obj_zz, region_z, region_zz);
    transfer_request_px  = PDCregion_transfer_create(&px[0], PDC_READ, obj_pxx, region_px, region_pxx);
    transfer_request_py  = PDCregion_transfer_create(&py[0], PDC_READ, obj_pyy, region_py, region_pyy);
    transfer_request_pz  = PDCregion_transfer_create(&pz[0], PDC_READ, obj_pzz, region_pz, region_pzz);
    transfer_request_id1 = PDCregion_transfer_create(&id1[0], PDC_READ, obj_id11, region_id1, region_id11);
    transfer_request_id2 = PDCregion_transfer_create(&id2[0], PDC_READ, obj_id22, region_id2, region_id22);

    ret = PDCregion_transfer_start(transfer_request_x);
    ret = PDCregion_transfer_start(transfer_request_y);
    ret = PDCregion_transfer_start(transfer_request_z);
    ret = PDCregion_transfer_start(transfer_request_px);
    ret = PDCregion_transfer_start(transfer_request_py);
    ret = PDCregion_transfer_start(transfer_request_pz);
    ret = PDCregion_transfer_start(transfer_request_id1);
    ret = PDCregion_transfer_start(transfer_request_id2);

#ifdef ENABLE_MPI
    MPI_Barrier(MPI_COMM_WORLD);
    t1 = MPI_Wtime();
    if (rank == 0) {
        printf("[CACHE_LOG] Transfer start time: %f\n", t1 - t0);
    }
#endif

    ret = PDCregion_transfer_wait(transfer_request_x);
    ret = PDCregion_transfer_wait(transfer_request_y);
    ret = PDCregion_transfer_wait(transfer_request_z);
    ret = PDCregion_transfer_wait(transfer_request_px);
    ret = PDCregion_transfer_wait(transfer_request_py);
    ret = PDCregion_transfer_wait(transfer_request_pz);
    ret = PDCregion_transfer_wait(transfer_request_id1);
    ret = PDCregion_transfer_wait(transfer_request_id2);

#ifdef ENABLE_MPI
    MPI_Barrier(MPI_COMM_WORLD);
    t0 = MPI_Wtime();
    if (rank == 0) {
        printf("[CACHE_LOG] Transfer wait time: %f\n", t0 - t1);
    }
#endif

    ret = PDCregion_transfer_close(transfer_request_x);
    ret = PDCregion_transfer_close(transfer_request_y);
    ret = PDCregion_transfer_close(transfer_request_z);
    ret = PDCregion_transfer_close(transfer_request_px);
    ret = PDCregion_transfer_close(transfer_request_py);
    ret = PDCregion_transfer_close(transfer_request_pz);
    ret = PDCregion_transfer_close(transfer_request_id1);
    ret = PDCregion_transfer_close(transfer_request_id2);

#ifdef ENABLE_MPI
    MPI_Barrier(MPI_COMM_WORLD);
    t1 = MPI_Wtime();
    if (rank == 0) {
        printf("[CACHE_LOG] Transfer close time: %f\n", t1 - t0);
    }
#endif

    if (PDCregion_close(region_x) < 0)
        printf("fail to close region region_x\n");

    if (PDCregion_close(region_y) < 0)
        printf("fail to close region region_y\n");

    if (PDCregion_close(region_z) < 0)
        printf("fail to close region region_z\n");

    if (PDCregion_close(region_px) < 0)
        printf("fail to close region region_px\n");

    if (PDCregion_close(region_py) < 0)
        printf("fail to close region region_py\n");

    if (PDCobj_close(region_pz) < 0)
        printf("fail to close region region_pz\n");

    if (PDCobj_close(region_id1) < 0)
        printf("fail to close region region_id1\n");

    if (PDCobj_close(region_id2) < 0)
        printf("fail to close region region_id2\n");

    if (PDCregion_close(region_xx) < 0)
        printf("fail to close region region_xx\n");

    if (PDCregion_close(region_yy) < 0)
        printf("fail to close region region_yy\n");

    if (PDCregion_close(region_zz) < 0)
        printf("fail to close region region_zz\n");

    if (PDCregion_close(region_pxx) < 0)
        printf("fail to close region region_pxx\n");

    if (PDCregion_close(region_pyy) < 0)
        printf("fail to close region region_pyy\n");

    if (PDCregion_close(region_pzz) < 0)
        printf("fail to close region region_pzz\n");

    if (PDCregion_close(region_id11) < 0)
        printf("fail to close region region_id11\n");

    if (PDCregion_close(region_id22) < 0)
        printf("fail to close region region_id22\n");

#ifdef ENABLE_MPI
    MPI_Barrier(MPI_COMM_WORLD);
    end = MPI_Wtime();
    if (rank == 0) {
        printf("[CACHE_LOG] Total execution time - Cache Insertion: %f\n", end - start);
    }
    start = MPI_Wtime();
    printf("PDC execute multiple transfer request regarding partial regions\n");
#endif

    mysize[0] = numparticles / num_transfer_request;

    // create a region
    region_x   = PDCregion_create(ndim, offset, mysize);
    region_y   = PDCregion_create(ndim, offset, mysize);
    region_z   = PDCregion_create(ndim, offset, mysize);
    region_px  = PDCregion_create(ndim, offset, mysize);
    region_py  = PDCregion_create(ndim, offset, mysize);
    region_pz  = PDCregion_create(ndim, offset, mysize);
    region_id1 = PDCregion_create(ndim, offset, mysize);
    region_id2 = PDCregion_create(ndim, offset, mysize);

    for (i = 0; i < num_transfer_request; i++) {
        offset_remote[0] = rank * numparticles + (mysize[0] * i);

        region_xx   = PDCregion_create(ndim, offset_remote, mysize);
        region_yy   = PDCregion_create(ndim, offset_remote, mysize);
        region_zz   = PDCregion_create(ndim, offset_remote, mysize);
        region_pxx  = PDCregion_create(ndim, offset_remote, mysize);
        region_pyy  = PDCregion_create(ndim, offset_remote, mysize);
        region_pzz  = PDCregion_create(ndim, offset_remote, mysize);
        region_id11 = PDCregion_create(ndim, offset_remote, mysize);
        region_id22 = PDCregion_create(ndim, offset_remote, mysize);

#ifdef ENABLE_MPI
        MPI_Barrier(MPI_COMM_WORLD);
        t0 = MPI_Wtime();
#endif

        transfer_request_x  = PDCregion_transfer_create(&x[0], PDC_READ, obj_xx, region_x, region_xx);
        transfer_request_y  = PDCregion_transfer_create(&y[0], PDC_READ, obj_yy, region_y, region_yy);
        transfer_request_z  = PDCregion_transfer_create(&z[0], PDC_READ, obj_zz, region_z, region_zz);
        transfer_request_px = PDCregion_transfer_create(&px[0], PDC_READ, obj_pxx, region_px, region_pxx);
        transfer_request_py = PDCregion_transfer_create(&py[0], PDC_READ, obj_pyy, region_py, region_pyy);
        transfer_request_pz = PDCregion_transfer_create(&pz[0], PDC_READ, obj_pzz, region_pz, region_pzz);
        transfer_request_id1 =
            PDCregion_transfer_create(&id1[0], PDC_READ, obj_id11, region_id1, region_id11);
        transfer_request_id2 =
            PDCregion_transfer_create(&id2[0], PDC_READ, obj_id22, region_id2, region_id22);

        ret = PDCregion_transfer_start(transfer_request_x);
        ret = PDCregion_transfer_start(transfer_request_y);
        ret = PDCregion_transfer_start(transfer_request_z);
        ret = PDCregion_transfer_start(transfer_request_px);
        ret = PDCregion_transfer_start(transfer_request_py);
        ret = PDCregion_transfer_start(transfer_request_pz);
        ret = PDCregion_transfer_start(transfer_request_id1);
        ret = PDCregion_transfer_start(transfer_request_id2);

#ifdef ENABLE_MPI
        MPI_Barrier(MPI_COMM_WORLD);
        t1 = MPI_Wtime();
        if (rank == 0) {
            printf("[CACHE_LOG] Transfer start time: %f\n", t1 - t0);
        }
#endif
        ret = PDCregion_transfer_wait(transfer_request_x);
        ret = PDCregion_transfer_wait(transfer_request_y);
        ret = PDCregion_transfer_wait(transfer_request_z);
        ret = PDCregion_transfer_wait(transfer_request_px);
        ret = PDCregion_transfer_wait(transfer_request_py);
        ret = PDCregion_transfer_wait(transfer_request_pz);
        ret = PDCregion_transfer_wait(transfer_request_id1);
        ret = PDCregion_transfer_wait(transfer_request_id2);

#ifdef ENABLE_MPI
        MPI_Barrier(MPI_COMM_WORLD);
        t0 = MPI_Wtime();
        if (rank == 0) {
            printf("[CACHE_LOG] Transfer wait time: %f\n", t0 - t1);
        }
#endif

        ret = PDCregion_transfer_close(transfer_request_x);
        ret = PDCregion_transfer_close(transfer_request_y);
        ret = PDCregion_transfer_close(transfer_request_z);
        ret = PDCregion_transfer_close(transfer_request_px);
        ret = PDCregion_transfer_close(transfer_request_py);
        ret = PDCregion_transfer_close(transfer_request_pz);
        ret = PDCregion_transfer_close(transfer_request_id1);
        ret = PDCregion_transfer_close(transfer_request_id2);

#ifdef ENABLE_MPI
        MPI_Barrier(MPI_COMM_WORLD);
        t1 = MPI_Wtime();
        if (rank == 0) {
            printf("[CACHE_LOG] Transfer close time: %f\n", t1 - t0);
        }
#endif

        // Check if data written previously has been correctly read.
        for (j = offset_remote[0]; j < offset_remote[0] + mysize[0]; ++j) {
            if (id1[j] != j) {
                printf("wrong value %d!=%d @ line %d\n", id1[j], j, __LINE__);
                break;
            }
            if (id2[j] != j * 2) {
                printf("wrong value %d!=%d @ line %d\n", id2[j], j, __LINE__);
                break;
            }
        }

        if (PDCregion_close(region_xx) < 0)
            printf("fail to close region region_xx\n");

        if (PDCregion_close(region_yy) < 0)
            printf("fail to close region region_yy\n");

        if (PDCregion_close(region_zz) < 0)
            printf("fail to close region region_zz\n");

        if (PDCregion_close(region_pxx) < 0)
            printf("fail to close region region_pxx\n");

        if (PDCregion_close(region_pyy) < 0)
            printf("fail to close region region_pyy\n");

        if (PDCregion_close(region_pzz) < 0)
            printf("fail to close region region_pzz\n");

        if (PDCregion_close(region_id11) < 0)
            printf("fail to close region region_id11\n");

        if (PDCregion_close(region_id22) < 0)
            printf("fail to close region region_id22\n");

#ifdef ENABLE_MPI
        MPI_Barrier(MPI_COMM_WORLD);
#endif
    }

    if (PDCregion_close(region_x) < 0)
        printf("fail to close region region_xx\n");

    if (PDCregion_close(region_y) < 0)
        printf("fail to close region region_yy\n");

    if (PDCregion_close(region_z) < 0)
        printf("fail to close region region_zz\n");

    if (PDCregion_close(region_px) < 0)
        printf("fail to close region region_pxx\n");

    if (PDCregion_close(region_py) < 0)
        printf("fail to close region region_pyy\n");

    if (PDCregion_close(region_pz) < 0)
        printf("fail to close region region_pzz\n");

    if (PDCregion_close(region_id1) < 0)
        printf("fail to close region region_id11\n");

    if (PDCregion_close(region_id2) < 0)
        printf("fail to close region region_id22\n");

    if (PDCobj_close(obj_xx) < 0)
        printf("fail to close obj_xx\n");

    if (PDCobj_close(obj_yy) < 0)
        printf("fail to close object obj_yy\n");

    if (PDCobj_close(obj_zz) < 0)
        printf("fail to close object obj_zz\n");

    if (PDCobj_close(obj_pxx) < 0)
        printf("fail to close object obj_pxx\n");

    if (PDCobj_close(obj_pyy) < 0)
        printf("fail to close object obj_pyy\n");

    if (PDCobj_close(obj_pzz) < 0)
        printf("fail to close object obj_pzz\n");

    if (PDCobj_close(obj_id11) < 0)
        printf("fail to close object obj_id11\n");

    if (PDCobj_close(obj_id22) < 0)
        printf("fail to close object obj_id22\n");

    // close a container
    if (PDCcont_close(cont_id) < 0)
        printf("fail to close container c1\n");

    if (PDCclose(pdc_id) < 0)
        printf("fail to close PDC\n");

    free(x);
    free(y);
    free(z);
    free(px);
    free(py);
    free(pz);
    free(id1);
    free(id2);
    free(offset);
    free(offset_remote);
    free(mysize);

#ifdef ENABLE_MPI
    end = MPI_Wtime();
    if (rank == 0) {
        printf("[CACHE_LOG] Total execution time: %f\n", end - start);
    }
    MPI_Finalize();
#endif

    return 0;
}
