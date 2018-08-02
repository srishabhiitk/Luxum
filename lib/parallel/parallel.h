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


 /*! \file  parallel.h
  *
  * \brief  This file contains header of parallel class.
  * \author Rishabh Sahu
  * \date   Mar 2018
  */



#ifndef PARALLEL_H
#define PARALLEL_H

#include <blitz/array.h>
#include <mpi.h>
#include "reader.h"

class parallel {
    private:
        void assignRanks();
        void getNeighbours();

    public:
        // ALL THE INTEGERS USED BELOW ARE POSITIVE. STILL IT IS BETTER TO USE int INSTEAD OF unsigned int [1]
        /** The MPI rank of each sub-domain */
        int rank;

        /** The total number of cores available for computation */
        int nProc;

        /** npX and npY indicates the number of sub-domain divisions along the X and Y directions respectively */
        //@{
        const int npX, npY;
        //@}

        /** xRank and yRank indicates the rank in terms of sub-domain divisions along the X and Y directions respectively.
         *  Like the global rank variable, these values also start from 0 to npX -1 and npY - 1 respectively. */
        //@{
        int xRank, yRank;
        //@}

        /** Array of ranks of the 4 neighbouring sub-domains - Left, Right, Front, Back */
        blitz::Array<int, 1> nearRanks;

        parallel(const reader &iDat);

        inline int findRank(int xR, int yR);
};

/**
 ********************************************************************************************************************************************
 *  \class parallel parallel.h "lib/parallel.h"
 *          \brief      Class for all the global variables and functions related to parallelization.
 *
 *                      After MPI_Init, every process has its own rank. Moreover, after performing domain decomposition, each process has its own
 *                      xRank and yRank to identify its position within the global computational domain.
 *                      These data, along with the data to identify the neighbouring processes for inter-domain communication are stored in the
 *                      <B>parallel</B> class.
 *                      This class is initialized only once at the start of the solver.
 ********************************************************************************************************************************************
 */

#endif
