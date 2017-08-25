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

#ifdef ENABLE_MPI
#include "mpi.h"
#endif

#include "pdc.h"


int main(int argc, char **argv) {
    pdcid_t pdc, create_prop1, create_prop2, create_prop;
    PDC_prop_type type;
    int rank = 0, size = 1;
    
#ifdef ENABLE_MPI
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
#endif
    // create a pdc
    pdc = PDC_init("pdc");

    // create an object property
    create_prop1 = PDCprop_create(PDC_OBJ_CREATE, pdc);
    if(create_prop1 <= 0) {
        printf("Fail to create @ line %d\n", __LINE__);
    }
    // create another object property
    create_prop2 = PDCprop_create(PDC_OBJ_CREATE, pdc);
    if(create_prop2 <= 0) {
        printf("Fail to create @ line %d\n", __LINE__);
    }

    if(PDCprop_close(create_prop1)<0)
        printf("Fail to close property @ line %d\n", __LINE__);
    else
        printf("successfully close property # %lld\n", create_prop1);
    if(PDCprop_close(create_prop2)<0)
        printf("Fail to close property @ line %d\n", __LINE__);
    else
        printf("successfully close property # %lld\n", create_prop2);

    // create a container property
    create_prop = PDCprop_create(PDC_CONT_CREATE, pdc);
    if(create_prop > 0) {
        if(type == PDC_CONT_CREATE)
            printf("Create a container property, id is %lld\n", create_prop);
        else if(type == PDC_OBJ_CREATE)
            printf("Create an object property, id is %lld\n", create_prop);
    }
    else
        printf("Fail to create @ line  %d!\n", __LINE__);

    // close property
   if(PDCprop_close(create_prop)<0)
       printf("Fail to close property @ line %d\n", __LINE__);
   else
       printf("successfully close property # %lld\n", create_prop);

    // close a pdc
    if(PDC_close(pdc) < 0)
       printf("fail to close PDC\n");
    else
       printf("PDC is closed\n");

#ifdef ENABLE_MPI
    MPI_Finalize();
#endif
    return 0;
}