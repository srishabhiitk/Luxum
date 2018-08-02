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


 /*! \file  field.cc
  *
  * \brief  This file contains contructor of field class.
  * \author Rishabh Sahu
  * \date   Mar 2018
  */



#include "field.h"

/**
 ********************************************************************************************************************************************
 * \brief   Constructor of the field class
 *
 *          Initializes the data members using the input parameters and grid information.
 *          Resizes the Blitz array to required size and then initializes each element of
 *          array to zero. Also, initializes the MPI handle by calling the mpidata class.
 *
 * \param   gridData is a constant reference to grid information.
 * \param   xStag_ is a boolean representing whether the field is staggered in X direction.
 * \param   yStag_ is a boolean representing whether the field is staggered in Y direction.
 * \param   zStag_ is a boolean representing whether the field is staggered in Z direction.
 ********************************************************************************************************************************************
 */
field::field(const grid &gridData, bool xStag_, bool yStag_, bool zStag_): gridData(gridData){
  local_Nx=gridData.local_colloq_x;
  local_Ny=gridData.local_colloq_y;
  local_Nz=gridData.local_colloq_z;
  xStag=xStag_;
  yStag=yStag_;
  zStag=zStag_;


  if (xStag_){
    local_Nx=local_Nx-1;
  }
  if (yStag_){
    local_Ny=local_Ny-1;
  }
  if (zStag_){
    local_Nz=local_Nz-1;
  }

  F.resize(local_Nx,local_Ny,local_Nz);
  #pragma omp parallel for num_threads(gridData.inputData.n_threads) schedule(dynamic)
  for (int i=0;i<local_Nx;i++){
      for (int j=0;j<local_Ny;j++){
          for (int k=0;k<local_Nz;k++){
            F(i,j,k) = 0;
          }
      }
  }

  mpiHandle = new mpidata(F, xStag, yStag, gridData.rankData, gridData);
  // The call to following function marks the array elements that would need to shared
  // to different MPI nodes if syncing of this scalar field is required.
  mpiHandle->createSubarrays(xStag, yStag, zStag);

}
