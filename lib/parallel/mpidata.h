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
 * 4. Neitfileio.ccher the name of the Simulation and Modeling Lab, IIT Kanpur nor the
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


 /*! \file  mpidata.h
  *
  * \brief  This file contains header of mpidata class.
  * \author Rishabh Sahu
  * \date   Mar 2018
  */



#ifndef MPIDATA_H
#define MPIDATA_H

#include <blitz/array.h>
#include <mpi.h>

#include "parallel.h"
#include "grid.h"

class mpidata {
    private:
        /** MPI subarray datatype for the slice of data to be sent to the neighbouring sub-domain to left along x-direction */
        MPI_Datatype sendSubarrayX0;
        /** MPI subarray datatype for the slice of data to be sent to the neighbouring sub-domain to right along x-direction */
        MPI_Datatype sendSubarrayX1;
        /** MPI subarray datatype for the slice of data to be sent to the neighbouring sub-domain to left along y-direction (front side) */
        MPI_Datatype sendSubarrayY0;
        /** MPI subarray datatype for the slice of data to be sent to the neighbouring sub-domain to right along y-direction (rear side) */
        MPI_Datatype sendSubarrayY1;

        /** MPI subarray datatype for the slice of data to be received from the neighbouring sub-domain to left along x-direction */
        MPI_Datatype recvSubarrayX0;
        /** MPI subarray datatype for the slice of data to be received from the neighbouring sub-domain to right along x-direction */
        MPI_Datatype recvSubarrayX1;
        /** MPI subarray datatype for the slice of data to be received from the neighbouring sub-domain to left along y-direction (front side) */
        MPI_Datatype recvSubarrayY0;
        /** MPI subarray datatype for the slice of data to be received from the neighbouring sub-domain to right along y-direction (rear side) */
        MPI_Datatype recvSubarrayY1;

        /** An array of MPI_Request data-types necessary for obtaining output from the non-blocking receive MPI_Irecv in the syncData function. */
        blitz::Array<MPI_Request, 1> recvRequest;

        /** An array of MPI_Status data-types necessary for obtaining output from the non-blocking receive MPI_Irecv in the syncData function. */
        blitz::Array<MPI_Status, 1> recvStatus;

        /** Array of values of the data field which needs to be synchronised across processors */
        blitz::Array<double, 3> dataField;

        /** A const reference to the global variables stored in the parallel class to access rank data */
        const parallel &rankData;

        /** A const reference to the global variables stored in the parallel class to access rank data */
        const grid &gridData;

        /** Blitz array containing the IDs of MPI nodes to which data is to be shared */
        blitz::Array<int,1> shareRanks;

    public:
        mpidata(blitz::Array<double, 3> inputArray, bool xStag_, bool yStag_, const parallel &parallelData, const grid &gridData_);

        void createSubarrays(const bool xStag, const bool yStag, const bool zStag);

        void syncData();
};

/**
 ********************************************************************************************************************************************
 *  \class mpidata mpidata.h "lib/mpidata.h"
 *  \brief Class to store MPI derived datatypes for individual arrays.
 *
 *  Since the solver uses staggered grids, the data of scalar fields are stored in arrays of different sizes in each direction
 *  depending on whether the field is staggered or not in that direction.
 *  As a result, the limits of the sub-arrays to be sent across inter-processor boundaries is different for different arrays.
 *  Hence the <B>mpidata</B> class contains MPI_SUBARRAY derived datatypes to be initialized along with different fields in order
 *  to store their sub-arrays for inter-processor communication.
 ********************************************************************************************************************************************
 */

#endif
