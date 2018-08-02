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


 /*! \file  field.h
  *
  * \brief  This file contains header of field class.
  * \author Rishabh Sahu
  * \date   Mar 2018
  */




#ifndef FIELD_H
#define FIELD_H

#include <blitz/array.h>
#include "mpidata.h"
#include "grid.h"

class field {

    public:
        /** The blitz array that contains the data of scalar field. */
        blitz::Array<double, 3> F;

        /** Boolean representing whether field is scattered in X direction. */
        bool xStag;

        /** Boolean representing whether field is scattered in Y direction. */
        bool yStag;

        /** Boolean representing whether field is scattered in Z direction. */
        bool zStag;

        /** Array size of scalar field in X direction. */
        int local_Nx;

        /** Array size of scalar field in X direction. */
        int local_Ny;

        /** Array size of scalar field in Y direction. */
        int local_Nz;

        /** Array size of scalar field in Z direction. */
        const grid &gridData;

        /** MPI handle used to sync the field data between different MPI nodes. */
        mpidata *mpiHandle;


        field(const grid &gridData, bool xStag_, bool yStag_, bool zStag_);

};

/**
 ********************************************************************************************************************************************
 *  \class field field.h "lib/field.h"
 *  \brief  A class which abstracts concept of a scalar field for Yee's algorithm.
 *
 *          Yee's algorithm requires different lattice positions for different components of
 *          E field and H field. These positions are represented by three boolean variables that
 *          keep track of whether the field component is staggered in X, Y and Z direction.
 *          Additionally, the class contains information about the array size of scalar field and
 *          constant reference to grid information.
 ********************************************************************************************************************************************
 */


#endif
