/* Copyright (c) 2018, Mahendra Kumar Verma
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 * 1. Redistributions of source code must retain the above copyright
 *  notice, this list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright
 *  notice, this list of conditions and the following disclaimer in the
 *  documentation and/or other materials provided with the distribution.
 * 3. All advertising materials mentioning features or use of this software
 *  must display the following acknowledgement:
 *  This product includes software developed by the Simulation and Modeling Lab,
 *  IIT Kanpur.
 * 4. Neither the name of the Simulation and Modeling Lab, IIT Kanpur nor the
 *  names of its contributors may be used to endorse or promote products
 *  derived from this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY Mahendra Kumar Verma ''AS IS'' AND ANY
 * EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL <COPYRIGHT HOLDER> BE LIABLE FOR ANY
 * DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 *SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 */


 /*! \file  mpidata.cc
  *
  * \brief  This file contains contructor and definitions of member functions
  *         of mpidata class.
  * \author Rishabh Sahu
  * \date   Mar 2018
  */



#include "mpidata.h"
#include <mpi.h>

/**
 ********************************************************************************************************************************************
 * \brief   Constructor of the mpidata class
 *
 *          The short constructor of mpidata class merely resizes the array of shareRanks and MPI_Request datatypes.
 *
 * \param   inputArray is the blitz array whose sub-arrays have to be created and synchronised across processors.
 * \param   xStag is boolean representing whether the scalar field provided is staggered in X direction.
 * \param   yStag is boolean representing whether the scalar field provided is staggered in Y direction.
 * \param   parallelData is a const reference to the global data contained in the parallel class.
 * \param   gridData is a const reference to the grid information.
 ********************************************************************************************************************************************
 */
mpidata::mpidata(blitz::Array<double, 3> inputArray, bool xStag, bool yStag, const parallel &parallelData, const grid &gridData_): dataField(inputArray), rankData(parallelData), gridData(gridData_) {
    shareRanks.resize(4);
    shareRanks=rankData.nearRanks;
    if (xStag){
        shareRanks(0) = MPI_PROC_NULL;
        shareRanks(1) = MPI_PROC_NULL;
    }
    if (yStag || gridData.isPlanar){
        shareRanks(2) = MPI_PROC_NULL;
        shareRanks(3) = MPI_PROC_NULL;
    }

    recvStatus.resize(4);
    recvRequest.resize(4);
}

/**
 ********************************************************************************************************************************************
 * \brief   Function to create the subarray MPI_Datatypes
 *
 *          Function calculates the required subarrays that need to be shared across MPI nodes and stores them in
 *          MPI_Datatypes that can simply be used whenever syncing is required.
 *
 *
 * \param   xStag is boolean representing whether the scalar field provided is staggered in X direction.
 * \param   yStag is boolean representing whether the scalar field provided is staggered in Y direction.
 * \param   zStag is boolean representing whether the scalar field provided is staggered in Z direction.
 ********************************************************************************************************************************************
 */
void mpidata::createSubarrays(const bool xStag, const bool yStag, const bool zStag) {

    blitz::TinyVector<int, 3> subSize;
    blitz::TinyVector<int, 3> saStarts;
    blitz::TinyVector<int, 3> globSize;

    globSize(0) = gridData.local_colloq_x;
    globSize(1) = gridData.local_colloq_y;
    globSize(2) = gridData.local_colloq_z;

    if (xStag){
        globSize(0) = globSize(0) - 1;
    }
    if (yStag){
        globSize(1) = globSize(1) - 1;
    }
    if (zStag){
        globSize(2) = globSize(2) - 1;
    }
    saStarts = 0;
    subSize = globSize;


    // Subarray for left side
    subSize(0) = 1;
    saStarts(0) = 1;
    MPI_Type_create_subarray(3, globSize.data(), subSize.data(), saStarts.data(), MPI_ORDER_C, MPI_DOUBLE_PRECISION, &sendSubarrayX0);
    MPI_Type_commit(&sendSubarrayX0);

    saStarts(0) = 0;
    MPI_Type_create_subarray(3, globSize.data(), subSize.data(), saStarts.data(), MPI_ORDER_C, MPI_DOUBLE_PRECISION, &recvSubarrayX0);
    MPI_Type_commit(&recvSubarrayX0);

    // Subarray for right side
    subSize(0) = 1;
    saStarts(0) = globSize(0) - 2;
    MPI_Type_create_subarray(3, globSize.data(), subSize.data(), saStarts.data(), MPI_ORDER_C, MPI_DOUBLE_PRECISION, &sendSubarrayX1);
    MPI_Type_commit(&sendSubarrayX1);

    saStarts(0) = globSize(0) - 1;
    MPI_Type_create_subarray(3, globSize.data(), subSize.data(), saStarts.data(), MPI_ORDER_C, MPI_DOUBLE_PRECISION, &recvSubarrayX1);
    MPI_Type_commit(&recvSubarrayX1);



    subSize = globSize;
    saStarts = 0;


    // Front and back is only required if grid is three dimensional.
    if (!gridData.isPlanar){
        // Subarray for front side
        subSize(1) = 1;
        saStarts(1) = 1;
        MPI_Type_create_subarray(3, globSize.data(), subSize.data(), saStarts.data(), MPI_ORDER_C, MPI_DOUBLE_PRECISION, &sendSubarrayY0);
        MPI_Type_commit(&sendSubarrayY0);


        saStarts(1) = 0;
        MPI_Type_create_subarray(3, globSize.data(), subSize.data(), saStarts.data(), MPI_ORDER_C, MPI_DOUBLE_PRECISION, &recvSubarrayY0);
        MPI_Type_commit(&recvSubarrayY0);

        // Subarray for back side
        subSize(1) = 1;
        saStarts(1) = globSize(1) - 2;
        MPI_Type_create_subarray(3, globSize.data(), subSize.data(), saStarts.data(), MPI_ORDER_C, MPI_DOUBLE_PRECISION, &sendSubarrayY1);
        MPI_Type_commit(&sendSubarrayY1);

        saStarts(1) = globSize(1) - 1;
        MPI_Type_create_subarray(3, globSize.data(), subSize.data(), saStarts.data(), MPI_ORDER_C, MPI_DOUBLE_PRECISION, &recvSubarrayY1);
        MPI_Type_commit(&recvSubarrayY1);
    }

}

/**
 ********************************************************************************************************************************************
 * \brief   Function to sync data across all MPI nodes.
 *
 *          This is the core function of the mpidata class.
 *          The end slices of each sub-domain recieves data from their corresponding neighbouring sub-domains,
 *          while the interior slices of each sub-domain sends data to their corresponding neighbouring sub-domains.
 *
 *          All the data slices are send as subarray MPI derived data-types created in the \ref createSubarrays function.
 *          As a result, \ref syncData must be called only after the subarrays have been created.
 *
 *          The data transfer is implemented here with a mixture of blocking and non-blocking communication calls.
 *          The recieves are non-blocking, while the sends are blocking. This combination prevents inter-processor deadlock.
 ********************************************************************************************************************************************
 */
void mpidata::syncData() {
    recvRequest = MPI_REQUEST_NULL;

    MPI_Irecv(dataField.dataFirst(), 1, recvSubarrayX0, shareRanks(0), 1, MPI_COMM_WORLD, &recvRequest(0));
    MPI_Irecv(dataField.dataFirst(), 1, recvSubarrayX1, shareRanks(1), 1, MPI_COMM_WORLD, &recvRequest(1));
    MPI_Send(dataField.dataFirst(), 1, sendSubarrayX0, shareRanks(0), 1, MPI_COMM_WORLD);
    MPI_Send(dataField.dataFirst(), 1, sendSubarrayX1, shareRanks(1), 1, MPI_COMM_WORLD);

    //Front and back sides only need to be synced if grid is three-dimensional
    if (!gridData.isPlanar){
        MPI_Irecv(dataField.dataFirst(), 1, recvSubarrayY0, shareRanks(2), 1, MPI_COMM_WORLD, &recvRequest(2));
        MPI_Irecv(dataField.dataFirst(), 1, recvSubarrayY1, shareRanks(3), 1, MPI_COMM_WORLD, &recvRequest(3));
        MPI_Send(dataField.dataFirst(), 1, sendSubarrayY0, shareRanks(2), 1, MPI_COMM_WORLD);
        MPI_Send(dataField.dataFirst(), 1, sendSubarrayY1, shareRanks(3), 1, MPI_COMM_WORLD);
    }


    MPI_Waitall(4, recvRequest.dataFirst(), recvStatus.dataFirst());
}
