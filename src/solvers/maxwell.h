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


 /*! \file  maxwell.h
  *
  * \brief  Header file for class maxwell.
  * \author Rishabh Sahu
  * \date   Mar 2018
  */





#ifndef MAXWELL_H
#define MAXWELL_H

#include <iostream>
#include <math.h>
#include <blitz/array.h>
#include "parallel.h"
#include "reader.h"
#include "grid.h"
#include "vfield.h"
#include "fileio.h"
#include <sstream>
#include <sys/time.h>

class maxwell {
    public:
        /** Stores the number of openMP threads per MPI node. */
        int n_threads;

        /** Number of total timesteps to be executed. */
        int num_timesteps;

        /** Stores value of electric permeability of vaccum. */
        double epsilon0;

        /** Stores value of magnetic permeability of vaccum. */
        double mew0;

        /** Stores courant factor of Yee's algorithm. */
        double S;

        /** Stores value of speed of light calculated from 'epsilon0' and 'mew0'. */
        double c;

        /** Resolution in time. */
        double dt;

        /** Resolution in X-direction. */
        double dx;

        /** Resolution in Y-direction. */
        double dy;

        /** Resolution in Z-direction. */
        double dz;

        /** PML parameter: Maximum value of a. */
        double a_max;

        /** PML parameter: Maximum value of k. */
        double k_max;

        /** PML parameter: Polynomial exponent for sigma grading. */
        double m_pml;

        /** PML parameter: Polynomial exponent for a function grading. */
        double ma_pml;

        /** PML parameter: Optimal value of sigma max for 5 to 10 cell thick PML. */
        double sigma_optimal;

        /** PML parameter: Maximum value of sigma function. */
        double sigma_max;

        /** PML parameter: Length of PML in units similar to that of dx. */
        double d_pml;

        /** Plane-wave source parameter: Width of Gaussian pulse with flat wave-front. */
        double sigma;

        /** Plane-wave source parameter: Offset value of Gaussian pulse with flat wave-front.
        * It can be used to delay or advance the pulse arrival time in Total Field zone. */
        double init;

        /** Wavelength of light. */
        double lambda;

        /** Poiter to object of reader class. */
        reader &inputData;

        /** Poiter to object of parallel class. */
        parallel &mpi;

        /** Poiter to object of grid class. */
        grid &gridData;

        /** Boolean storing whether to use plane-wave source or not. */
        bool usePW;

        /** Boolean storing whether to use PML boundary condition or not. */
        bool usePML;

        /** Boolean storing whether to use point source or not. */
        bool usePS;

        /** Coefficient representing X-direction of H-field of plane-wave. */
        double H_x_dir;

        /** Coefficient representing Y-direction of H-field of plane-wave. */
        double H_y_dir;

        /** Coefficient representing Z-direction of H-field of plane-wave. */
        double H_z_dir;

        /** Coefficient representing X-direction of E-field of plane-wave. */
        double E_x_dir;

        /** Coefficient representing Y-direction of E-field of plane-wave. */
        double E_y_dir;

        /** Coefficient representing Z-direction of E-field of plane-wave. */
        double E_z_dir;

        /** Coefficient representing X-direction of wavevector of plane-wave. */
        double kx;

        /** Coefficient representing Y-direction of wavevector of plane-wave. */
        double ky;

        /** Coefficient representing Z-direction of wavevector of plane-wave. */
        double kz;

        /** Global X-index marking start of Total-field zone for plane-wave. */
        int p1;

        /** Global Y-index marking start of Total-field zone for plane-wave. */
        int q1;

        /** Global Z-index marking start of Total-field zone for plane-wave. */
        int r1;

        /** Global X-index marking end of Total-field zone for plane-wave. */
        int p2;

        /** Global Y-index marking end of Total-field zone for plane-wave. */
        int q2;

        /** Global Z-index marking end of Total-field zone for plane-wave. */
        int r2;

        /** Stores the local X-index from where plane-wave has to start for specific MPI-node. */
        int start_x;

        /** Stores the local Y-index from where plane-wave has to start for specific MPI-node. */
        int start_y;

        /** Stores the local Z-index from where plane-wave has to start for specific MPI-node. */
        int start_z;

        /** Stores the local X-index till where plane-wave has to go for specific MPI-node. */
        int end_x;

        /** Stores the local Y-index till where plane-wave has to go for specific MPI-node. */
        int end_y;

        /** Stores the local Z-index till where plane-wave has to go for specific MPI-node. */
        int end_z;

        /** Amplitude of point source. */
        int amp_source;

        /** Angular-frequency of point source. */
        int omega_source;

        /** Global X-index of point source location. */
        int Nx_source;

        /** Global Y-index of point source location. */
        int Ny_source;

        /** Global Z-index of point source location. */
        int Nz_source;


        maxwell(reader &_inputData, parallel &_mpi, grid &_gridData);
        void solve();
        void solve3d();
        void solve3d_pml();
        void solve_planar();
        void solve_planar_pml();
        void plane_wave_execute(vfield *curl, int timestep);
        void plane_wave_execute(vfield *curl1, vfield *curl2, int timestep);
        void plane_wave_initialise();
        int sigma_fn(int dim, bool is_E, blitz::Array<double,1> w);
        int k_fn(int dim, bool is_E, blitz::Array<double,1> w);
        int a_fn(int dim, bool is_E, blitz::Array<double,1> w);
        int b_fn(int dim, bool is_E, blitz::Array<double,1> w);
        int c_fn(int dim, bool is_E, blitz::Array<double,1> w);
        bool check_global_limits(blitz::TinyVector<int,3> v, bool xstag, bool ystag);
        blitz::TinyVector<int,3> global_to_local(blitz::TinyVector<int,3> glob);
};

/**
 ********************************************************************************************************************************************
 *  \class maxwell maxwell.h "lib/maxwell.h"
 *  \brief  Contains many functions that implement Yee's algorithm along with various sources and boundary conditions.
 *
 *  The class initializes various variables according to user's choice as read by reader class.
 *  The class has various functions that implement Yee's algorithm for two and three dimensions along with different boundary condition.
 *  The class also contains various functions that collectively implement plane-wave source using Total-field scattered-field formulation.
 ********************************************************************************************************************************************
 */

#endif
