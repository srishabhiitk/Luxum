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


 /*! \file  reader.h
  *
  * \brief  This file contains header of reader class.
  * \author Rishabh Sahu
  * \date   Mar 2018
  */



#ifndef READER_H
#define READER_H

#include <string>
#include <fstream>
#include <blitz/array.h>
#include <yaml-cpp/yaml.h>

class reader {
    public:
        /** Number of processors in X direction. */
        int npX;

        /** Number of processors in Y direction. */
        int npY;

        /** Resolution in X-direction. */
        double dx;

        /** Resolution in Y-direction. */
        double dy;

        /** Resolution in Z-direction. */
        double dz;

        /** Stores courant factor of Yee's algorithm. */
        double S;

        /** Number of total timesteps to be executed. */
        int num_timesteps;

        /** Number of grid points in X-direction. */
        int Nx;

        /** Number of grid points in Y-direction. */
        int Ny;

        /** Number of grid points in Z-direction. */
        int Nz;

        /** Stores the number of openMP threads per MPI node. */
        int n_threads;

        /** Boolean to store whether the solution is planar or two-dimensional. */
        bool isPlanar;

        /** Boolean storing whether to use plane-wave source or not. */
        bool usePW;

        /** Boolean storing whether to use PML boundary condition or not. */
        bool usePML;

        /** Boolean storing whether to use point source or not. */
        bool usePS;

        /** Coefficient representing X-direction of wavevector of plane-wave. */
        double kx_pw;

        /** Coefficient representing Y-direction of wavevector of plane-wave. */
        double ky_pw;

        /** Coefficient representing Z-direction of wavevector of plane-wave. */
        double kz_pw;

        /** Coefficient representing X-direction of E-field of plane-wave. */
        double Ex_pw;

        /** Coefficient representing Y-direction of E-field of plane-wave. */
        double Ey_pw;

        /** Coefficient representing Z-direction of E-field of plane-wave. */
        double Ez_pw;

        /** Offset for plane-wave pulse. Can be used to hasten or delay arrival of pulse. */
        double off_pw;

        /** Global X-index marking start of Total-field zone for plane-wave. */
        int x1_pw;

        /** Global Y-index marking start of Total-field zone for plane-wave. */
        int y1_pw;

        /** Global Z-index marking start of Total-field zone for plane-wave. */
        int z1_pw;

        /** Global X-index marking end of Total-field zone for plane-wave. */
        int x2_pw;

        /** Global Y-index marking end of Total-field zone for plane-wave. */
        int y2_pw;

        /** Global Z-index marking end of Total-field zone for plane-wave. */
        int z2_pw;

        /** Global X-index of point source location. */
        int Nx_source;

        /** Global Y-index of point source location. */
        int Ny_source;

        /** Global Z-index of point source location. */
        int Nz_source;

        /** Amplitude of point source. */
        double amp_source;

        /** Angular-frequency of point source. */
        double omega_source;

        /** PML parameter: Length of PML in units similar to that of dx. */
        double d_pml;

        /** Plane-wave source parameter: Width of Gaussian pulse with flat wave-front. */
        double sigma;

        /** PML parameter: Polynomial exponent for sigma grading. */
        double sigma_m;

        /** PML parameter: Maximum value of kappa. */
        double kappa;

        /** PML parameter: Maximum value of a. */
        double a_max;

        /** PML parameter: Polynomial exponent for a function grading. */
        double a_m;


        reader();
        void readYAML();
        void checkData();

};

/**
 ********************************************************************************************************************************************
 *  \class reader reader.h "lib/reader.h"
 *          \brief      Contains all the global variables set by the user through the yaml file
 *
 *                      The class reads the paramters.yaml file and stores all the simulation paramters in publicly accessible constants.
 *                      The class also has a function to check the consistency of the user set paramters and throw exceptions.
 *                      The class is best initialized as a constant to prevent inadvertent tampering of the global variables it contains.
 ********************************************************************************************************************************************
 */

#endif
