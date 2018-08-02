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


 /*! \file  vfield.cc
  *
  * \brief  This file contains contructor and definitions of member functions
  *         of vfield class.
  * \author Rishabh Sahu
  * \date   Mar 2018
  */



#include <iostream>
#include "vfield.h"


/**
 ********************************************************************************************************************************************
 * \brief   Constructor of the vfield class
 *
 *          Initializes the data members, especially the pointers to scalar field objects are
 *          intialised according to need.
 *
 * \param   gridData is a constant reference to grid information.
 * \param   isFaceCentered_ is a boolean representing whether the vector field is face centered.
 ********************************************************************************************************************************************
 */
vfield::vfield(const grid &gridData, bool isFaceCentered_): gridData(gridData) {
    isFaceCentered = isFaceCentered_;
    isPlanar = gridData.isPlanar;

    if (isPlanar){
        if (isFaceCentered_){
            Vx = new field(gridData, false, false, true);
            Vy = NULL;
            Vz = new field(gridData, true, false, false);
        }
        else{
            Vx = NULL;
            Vy = new field(gridData, false, false, false);
            Vz = NULL;
        }
    }
    else{
        if (isFaceCentered_){
            Vx = new field(gridData, false, true, true);
            Vy = new field(gridData, true, false, true);
            Vz = new field(gridData, true, true, false);
        }
        else{
            Vx = new field(gridData, true, false, false);
            Vy = new field(gridData, false, true, false);
            Vz = new field(gridData, false, false, true);
        }
    }

}


/**
 ********************************************************************************************************************************************
 * \brief   Calculates curl of vector field of corresponding to object calling it.
 *
 *          The function takes one argument of vfield class where it stores the curl.
 *          The curl of a face centered vector field will be edge centered and vice-versa.
 *          The curl is only calculated if the vfield object passed and vfield object calling
 *          are of different types.
 *
 * \param   curl is a pointer to vfield object where the calculated curl will be stored.
 ********************************************************************************************************************************************
 */
void vfield::curl_3d(vfield *curl){
    int n_threads = gridData.inputData.n_threads;
    bool isFaceCentered_c = curl->isFaceCentered;
    if (isFaceCentered!=isFaceCentered_c){

        double dx = gridData.inputData.dx;
        double dy = gridData.inputData.dy;
        double dz = gridData.inputData.dz;

        int Nx = gridData.local_colloq_x;
        int Ny = gridData.local_colloq_y;
        int Nz = gridData.local_colloq_z;



        if (!isFaceCentered){
            #pragma omp parallel for num_threads(n_threads) schedule(dynamic)
            for (int i=0;i<Nx;i++){
                for (int j=0;j<(Ny-1);j++){
                    for (int k=0;k<(Nz-1);k++){
                      curl->Vx->F(i,j,k) = (Vz->F(i,j+1,k)-Vz->F(i,j,k))/dy - (Vy->F(i,j,k+1)-Vy->F(i,j,k))/dz;
                    }
                }
            }
            #pragma omp parallel for num_threads(n_threads) schedule(dynamic)
            for (int j=0;j<Ny;j++){
                for (int i=0;i<(Nx-1);i++){
                    for (int k=0;k<(Nz-1);k++){
                      curl->Vy->F(i,j,k) = (Vx->F(i,j,k+1)-Vx->F(i,j,k))/dz - (Vz->F(i+1,j,k)-Vz->F(i,j,k))/dx;
                    }
                }
            }
            #pragma omp parallel for num_threads(n_threads) schedule(dynamic)
            for (int k=0;k<Nz;k++){
                for (int j=0;j<(Ny-1);j++){
                    for (int i=0;i<(Nx-1);i++){
                      curl->Vz->F(i,j,k) = (Vy->F(i+1,j,k)-Vy->F(i,j,k))/dx - (Vx->F(i,j+1,k)-Vx->F(i,j,k))/dy;
                    }
                }
            }
        }
        else{
            #pragma omp parallel for num_threads(n_threads) schedule(dynamic)
            for (int i=0;i<(Nx-1);i++){
                for (int j=1;j<(Ny-1);j++){
                    for (int k=1;k<(Nz-1);k++){
                      curl->Vx->F(i,j,k) = (Vz->F(i,j,k)-Vz->F(i,j-1,k))/dy - (Vy->F(i,j,k)-Vy->F(i,j,k-1))/dz;
                    }
                }
            }
            #pragma omp parallel for num_threads(n_threads) schedule(dynamic)
            for (int j=0;j<(Ny-1);j++){
                for (int i=1;i<(Nx-1);i++){
                    for (int k=1;k<(Nz-1);k++){
                      curl->Vy->F(i,j,k) = (Vx->F(i,j,k)-Vx->F(i,j,k-1))/dz - (Vz->F(i,j,k)-Vz->F(i-1,j,k))/dx;
                    }
                }
            }
            #pragma omp parallel for num_threads(n_threads) schedule(dynamic)
            for (int k=0;k<(Nz-1);k++){
                for (int j=1;j<(Ny-1);j++){
                    for (int i=1;i<(Nx-1);i++){
                      curl->Vz->F(i,j,k) = (Vy->F(i,j,k)-Vy->F(i-1,j,k))/dx - (Vx->F(i,j,k)-Vx->F(i,j-1,k))/dy;
                    }
                }
            }
        }
    }
    else{
    std::cout<<"First two parameter should be staggered and unstaggered."<<std::endl;
    }

}


/**
 ********************************************************************************************************************************************
 * \brief   Calculates curl of vector field of corresponding to object calling it.
 *
 *          The function takes two argument of vfield class where it stores two parts of the curl.
 *          The two parts of the curl correspond to two final values corresponding to different dimensions
 *          which need to be subtracted to get the original curl in a given dimension.
 *          The first parameter contains the positive part of the curl.
 *          The second parameter contains the negative part of the curl.
 *          Actual curl = curl1 - curl2
 *
 *          The curl of a face centered vector field will be edge centered and vice-versa.
 *          The curl is only calculated if the vfield object passed and vfield object calling
 *          are of different types.
 *
 * \param   curl1 is a pointer to vfield object where positive part of the calculated curl will be stored.
 * \param   curl2 is a pointer to vfield object where negative part of the calculated curl will be stored.
 ********************************************************************************************************************************************
 */
void vfield::curl_3d_adv(vfield *curl1, vfield *curl2){
    int n_threads = gridData.inputData.n_threads;
    bool isFaceCentered_c1 = curl1->isFaceCentered;
    bool isFaceCentered_c2 = curl2->isFaceCentered;
    if (isFaceCentered!=isFaceCentered_c1 && isFaceCentered_c1==isFaceCentered_c2){

        double dx = gridData.inputData.dx;
        double dy = gridData.inputData.dy;
        double dz = gridData.inputData.dz;

        int Nx = gridData.local_colloq_x;
        int Ny = gridData.local_colloq_y;
        int Nz = gridData.local_colloq_z;



        if (!isFaceCentered){
            #pragma omp parallel for num_threads(n_threads) schedule(dynamic)
            for (int i=0;i<Nx;i++){
                for (int j=0;j<(Ny-1);j++){
                    for (int k=0;k<(Nz-1);k++){
                      curl1->Vx->F(i,j,k) = (Vz->F(i,j+1,k)-Vz->F(i,j,k))/dy;
                      curl2->Vx->F(i,j,k) = (Vy->F(i,j,k+1)-Vy->F(i,j,k))/dz;
                    }
                }
            }
            #pragma omp parallel for num_threads(n_threads) schedule(dynamic)
            for (int j=0;j<Ny;j++){
                for (int i=0;i<(Nx-1);i++){
                    for (int k=0;k<(Nz-1);k++){
                      curl1->Vy->F(i,j,k) = (Vx->F(i,j,k+1)-Vx->F(i,j,k))/dz;
                      curl2->Vy->F(i,j,k) = (Vz->F(i+1,j,k)-Vz->F(i,j,k))/dx;
                    }
                }
            }
            #pragma omp parallel for num_threads(n_threads) schedule(dynamic)
            for (int k=0;k<Nz;k++){
                for (int j=0;j<(Ny-1);j++){
                    for (int i=0;i<(Nx-1);i++){
                      curl1->Vz->F(i,j,k) = (Vy->F(i+1,j,k)-Vy->F(i,j,k))/dx;
                      curl2->Vz->F(i,j,k) = (Vx->F(i,j+1,k)-Vx->F(i,j,k))/dy;
                    }
                }
            }
        }
        else{
            #pragma omp parallel for num_threads(n_threads) schedule(dynamic)
            for (int i=0;i<(Nx-1);i++){
                for (int j=1;j<(Ny-1);j++){
                    for (int k=1;k<(Nz-1);k++){
                      curl1->Vx->F(i,j,k) = (Vz->F(i,j,k)-Vz->F(i,j-1,k))/dy;
                      curl2->Vx->F(i,j,k) = (Vy->F(i,j,k)-Vy->F(i,j,k-1))/dz;
                    }
                }
            }
            #pragma omp parallel for num_threads(n_threads) schedule(dynamic)
            for (int j=0;j<(Ny-1);j++){
                for (int i=1;i<(Nx-1);i++){
                    for (int k=1;k<(Nz-1);k++){
                      curl1->Vy->F(i,j,k) = (Vx->F(i,j,k)-Vx->F(i,j,k-1))/dz;
                      curl2->Vy->F(i,j,k) = (Vz->F(i,j,k)-Vz->F(i-1,j,k))/dx;
                    }
                }
            }
            #pragma omp parallel for num_threads(n_threads) schedule(dynamic)
            for (int k=0;k<(Nz-1);k++){
                for (int j=1;j<(Ny-1);j++){
                    for (int i=1;i<(Nx-1);i++){
                      curl1->Vz->F(i,j,k) = (Vy->F(i,j,k)-Vy->F(i-1,j,k))/dx;
                      curl2->Vz->F(i,j,k) = (Vx->F(i,j,k)-Vx->F(i,j-1,k))/dy;
                    }
                }
            }
        }
    }
    else{
    std::cout<<"Error: Passed parameters should both be either face-centered or not face-centered. The calling object should have opposite respective property!"<<std::endl;
    }

}





/**
 ********************************************************************************************************************************************
 * \brief   Calculates curl of vector field of corresponding to object calling it for planar grid.
 *
 *          The function takes one argument of vfield class where it stores the curl.
 *          The curl of a face centered vector field will be edge centered and vice-versa.
 *          The curl is only calculated if the vfield object passed and vfield object calling
 *          are of different types.
 *
 * \param   curl is a pointer to vfield object where the calculated curl will be stored.
 ********************************************************************************************************************************************
 */
void vfield::curl_planar(vfield *curl){
    int n_threads = gridData.inputData.n_threads;

    double dx = gridData.inputData.dx;
    double dz = gridData.inputData.dz;

    int Nx = gridData.local_colloq_x;
    int Nz = gridData.local_colloq_z;


    if (curl->Vy==NULL && Vx==NULL && Vz==NULL){
        #pragma omp parallel for num_threads(n_threads) schedule(dynamic)
        for (int i=0; i<(Nx); i++){
            for (int j=0; j<(Nz-1); j++){
                curl->Vx->F(i,0,j) = -(Vy->F(i,0,j+1)-Vy->F(i,0,j))/dz;
            }
        }
        #pragma omp parallel for num_threads(n_threads) schedule(dynamic)
        for (int i=0; i<(Nx-1); i++){
            for (int j=0; j<(Nz); j++){
                curl->Vz->F(i,0,j) = (Vy->F(i+1,0,j)-Vy->F(i,0,j))/dx;
            }
        }
    }
    else if (Vy==NULL && curl->Vx==NULL && curl->Vz==NULL){
        #pragma omp parallel for num_threads(n_threads) schedule(dynamic)
        for (int i=1; i<(Nx-1); i++){
            for (int j=1; j<(Nz-1); j++){
                curl->Vy->F(i,0,j) = (Vx->F(i,0,j)-Vx->F(i,0,j-1))/dz - (Vz->F(i,0,j)-Vz->F(i-1,0,j))/dx;
            }
        }
    }
    else{
        std::cout<<"The parameters passed to curl_planar function are not correctly defined for planar transverse fields!"<<std::endl;
    }
}


/**
 ********************************************************************************************************************************************
 * \brief   Calculates curl of vector field of corresponding to object calling it for planar grid.
 *
 *          The function takes two argument of vfield class where it stores two parts of the curl.
 *          The two parts of the curl correspond to two final values corresponding to different dimensions
 *          which need to be subtracted to get the original curl in a given dimension.
 *          The first parameter contains the positive part of the curl.
 *          The second parameter contains the negative part of the curl.
 *          Actual curl = curl1 - curl2
 *
 *          The curl of a face centered vector field will be edge centered and vice-versa.
 *          The curl is only calculated if the vfield object passed and vfield object calling
 *          are of different types.
 *
 * \param   curl1 is a pointer to vfield object where positive part of the calculated curl will be stored.
 * \param   curl2 is a pointer to vfield object where negative part of the calculated curl will be stored.
 ********************************************************************************************************************************************
 */
void vfield::curl_planar_adv(vfield *curl1, vfield *curl2){
    int n_threads = gridData.inputData.n_threads;

    double dx = gridData.inputData.dx;
    double dz = gridData.inputData.dz;

    int Nx = gridData.local_colloq_x;
    int Nz = gridData.local_colloq_z;


    if (curl1->Vy==NULL && curl2->Vy==NULL && Vx==NULL && Vz==NULL){
        #pragma omp parallel for num_threads(n_threads) schedule(dynamic)
        for (int i=0; i<(Nx); i++){
            for (int j=0; j<(Nz-1); j++){
                curl2->Vx->F(i,0,j) = (Vy->F(i,0,j+1)-Vy->F(i,0,j))/dz;
            }
        }
        #pragma omp parallel for num_threads(n_threads) schedule(dynamic)
        for (int i=0; i<(Nx-1); i++){
            for (int j=0; j<(Nz); j++){
                curl1->Vz->F(i,0,j) = (Vy->F(i+1,0,j)-Vy->F(i,0,j))/dx;
            }
        }
    }
    else if (Vy==NULL && curl1->Vx==NULL && curl1->Vz==NULL && curl2->Vx==NULL && curl2->Vz==NULL){
        #pragma omp parallel for num_threads(n_threads) schedule(dynamic)
        for (int i=1; i<(Nx-1); i++){
            for (int j=1; j<(Nz-1); j++){
                curl1->Vy->F(i,0,j) = (Vx->F(i,0,j)-Vx->F(i,0,j-1))/dz;
                curl2->Vy->F(i,0,j) = (Vz->F(i,0,j)-Vz->F(i-1,0,j))/dx;
            }
        }
    }
    else{
        std::cout<<"The parameters passed to curl_planar function are not correctly defined for planar transverse fields!"<<std::endl;
    }
}
