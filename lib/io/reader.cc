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


 /*! \file  reader.cc
  *
  * \brief  This file contains contructor and definitions of member functions
  *         of reader class.
  * \author Rishabh Sahu
  * \date   Mar 2018
  */



#include <iostream>
#include "reader.h"


/**
 ********************************************************************************************************************************************
 * \brief   Contructor to reader class
 *
 *          The brief constructor calls two functions - one to read the YAML file and other to
 *          validate the read arguments and abort the program if any dispute is found.
 ********************************************************************************************************************************************
 */
reader::reader() {
    readYAML();
    checkData();
}

/**
 ********************************************************************************************************************************************
 * \brief   Function to open the yaml file and read the parameters
 *
 *          The function opens the parameters.yaml file and reads the simulation parameters into its member variables that are publicly
 *          accessible.
 ********************************************************************************************************************************************
 */
void reader::readYAML() {
    std::ifstream inFile;
    inFile.open("parameters.yaml", std::ifstream::in);

    YAML::Node yamlNode;
    YAML::Parser parser(inFile);


    parser.GetNextDocument(yamlNode);
    yamlNode["isPlanar"] >> isPlanar;
    // isPlanar = true;


    parser.GetNextDocument(yamlNode);
    yamlNode["X Length"] >> Nx;
    yamlNode["Y Length"] >> Ny;
    yamlNode["Z Length"] >> Nz;
    yamlNode["X Resolution"] >> dx;
    yamlNode["Y Resolution"] >> dy;
    yamlNode["Z Resolution"] >> dz;
    yamlNode["Courant factor"] >> S;
    yamlNode["Number of timesteps"] >> num_timesteps;

    parser.GetNextDocument(yamlNode);
    yamlNode["X Number of Procs"] >> npX;
    yamlNode["Y Number of Procs"] >> npY;
    yamlNode["Number of Threads per node"] >> n_threads;

    parser.GetNextDocument(yamlNode);
    yamlNode["usePointSource"] >> usePS;
    yamlNode["Source X Pos"] >> Nx_source;
    yamlNode["Source Y Pos"] >> Ny_source;
    yamlNode["Source Z Pos"] >> Nz_source;
    yamlNode["Source Amplitude"] >> amp_source;
    yamlNode["Source Omega"] >> omega_source;
    yamlNode["usePW"] >> usePW;
    yamlNode["kx_pw"] >> kx_pw;
    yamlNode["ky_pw"] >> ky_pw;
    yamlNode["kz_pw"] >> kz_pw;
    yamlNode["Ex_pw"] >> Ex_pw;
    yamlNode["Ey_pw"] >> Ey_pw;
    yamlNode["Ez_pw"] >> Ez_pw;
    yamlNode["off_pw"] >> off_pw;
    yamlNode["x1_pw"] >> x1_pw;
    yamlNode["y1_pw"] >> y1_pw;
    yamlNode["z1_pw"] >> z1_pw;
    yamlNode["x2_pw"] >> x2_pw;
    yamlNode["y2_pw"] >> y2_pw;
    yamlNode["z2_pw"] >> z2_pw;


    parser.GetNextDocument(yamlNode);
    yamlNode["usePML"] >> usePML;
    yamlNode["Thickness of PML"] >> d_pml;
    yamlNode["sigma_max/sigma_optimal"] >> sigma;
    yamlNode["m value of sigma"] >> sigma_m;
    yamlNode["kappa max"] >> kappa;
    yamlNode["max value of a"] >> a_max;
    yamlNode["m value for a"] >> a_m;


    inFile.close();
}

/**
 ********************************************************************************************************************************************
 * \brief   Function to perform a check on the consistency of user-set parameters
 *
 *          In order to catch potential errors early on, a few basic checks are performed here to validate the paramters set
 *          by the user.
 *          Additional checks to be performed on the paramters can be added to this function if necessary.
 ********************************************************************************************************************************************
 */
 void reader::checkData() {

     if ((Nx-2)%npX!=0) {
         std::cout << "ERROR: Condition-- (Number of points - 2) mod Number of processors --is not satisfied in X-direction. Aborting" << std::endl;
         exit(0);
     }
     if ((Ny-2)%npY!=0) {
         std::cout << "ERROR: Condition-- (Number of points - 2) mod Number of processors --is not satisfied in Y-direction. Aborting" << std::endl;
         exit(0);
     }

     if (isPlanar==true && Ny!=1){
         std::cout << "ERROR: Please put number of grid points in Y-direction as 1 for planar solutions.";
         exit(0);
     }

     if (isPlanar==true && npY!=1){
         std::cout << "ERROR: Please put number of processors in Y-direction as 1 for planar solutions.";
         exit(0);
     }

 }
