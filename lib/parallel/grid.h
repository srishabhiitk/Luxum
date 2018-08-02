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


 /*! \file  grid.h
  *
  * \brief  This file contains header of grid class.
  * \author Rishabh Sahu
  * \date   Mar 2018
  */



#ifndef GRID_H
#define GRID_H

#include "reader.h"
#include "parallel.h"

class grid {
    public:
        /** A const reference to the global variables stored in the reader class to access user set parameters */
        const reader &inputData;

        /** A reference to the global variables stored in the parallel class to access MPI related parameters */
        const parallel &rankData;

        /** Stores the number of colloquated points in X-direction for a specific MPI node.  */
        int local_colloq_x;

        /** Stores the number of colloquated points in Y-direction for a specific MPI node.  */
        int local_colloq_y;

        /** Stores the number of colloquated points in Z-direction for a specific MPI node.  */
        int local_colloq_z;

        /** Stores the starting global index in X-direction for colloquated points for a specific MPI node. */
        int colloq_start_index_x;

        /** Stores the starting global index in Y-direction for colloquated points for a specific MPI node. */
        int colloq_start_index_y;

        /** Stores the starting global index in Z-direction for colloquated points for a specific MPI node. */
        int colloq_start_index_z;

        /** Stores the ending global index in X-direction for colloquated points for a specific MPI node. */
        int colloq_end_index_x;

        /** Stores the ending global index in Y-direction for colloquated points for a specific MPI node. */
        int colloq_end_index_y;

        /** Stores the ending global index in Z-direction for colloquated points for a specific MPI node. */
        int colloq_end_index_z;

        /** Boolean to store whether the grid is two-dimensional or not. */
        bool isPlanar;




        grid(const reader &solParam, parallel &parallelData);
};

/**
 ********************************************************************************************************************************************
 *  \class grid grid.h "lib/grid.h"
 *  \brief  Contains information about local colloquated indices for specific MPI nodes.
 *
 *          The datamembers of this class contain information about how and which points
 *          are located in specific MPI node with respect to full global arrays.
 ********************************************************************************************************************************************
 */

#endif
