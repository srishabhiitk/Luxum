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


 /*! \file  grid.cc
  *
  * \brief  This file contains contructor and definitions of member functions
  *         of grid class.
  * \author Rishabh Sahu
  * \date   Mar 2018
  */



#include "grid.h"

/**
 ********************************************************************************************************************************************
 * \brief   Constructor of the grid class
 *
 *          Initializes the datamembers of grid class such the number of colloquated points,
 *          stating and ending index of colloquated points for specific MPI nodes. This is
 *          calculated using the information from object of parallel class which contains
 *          relevant information like number of available processors in different dimensions.
 *
 * \param   solParam is a const reference to the global data contained in the reader class
 * \param   parallelData is a reference to the object of parallel class
 ********************************************************************************************************************************************
 */
grid::grid(const reader &solParam, parallel &parallelData): inputData(solParam),
                                                            rankData(parallelData) {


    // Storing some data from reader class.
    int Nx = inputData.Nx;
    int Ny = inputData.Ny;
    int npX = inputData.npX;
    int npY = inputData.npY;

    int xRank = rankData.xRank;
    int yRank = rankData.yRank;

    isPlanar = inputData.isPlanar;


    // Calculating staring and ending index accoring to number of processors available.
    colloq_start_index_x = xRank*(Nx-2)/npX;
    colloq_start_index_y = yRank*(Ny-2)/npY;
    colloq_start_index_z = 0;
    colloq_end_index_x = colloq_start_index_x+(Nx-2)/npX+1;
    colloq_end_index_y = colloq_start_index_y+(Ny-2)/npY+1;
    colloq_end_index_z = inputData.Nz - 1;

    // Calculating number of colloquated points per MPI node accoring to number of processors available.
    local_colloq_x = (Nx-2)/npX+2;
    local_colloq_y = (Ny-2)/npY+2;
    local_colloq_z = inputData.Nz;


    // If grid is two dimensional, certain changes to above used formula is required.
    if (isPlanar){
        local_colloq_y = 1;
        colloq_start_index_y = 0;
        colloq_end_index_y = 0;
    }

}
