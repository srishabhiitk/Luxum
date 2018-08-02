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


 /*! \file  maxwell.cc
  *
  * \brief  File contains constructor for Maxwell class along with other useful
  *         functions to implement plane-waves and PML.
  * \author Rishabh Sahu
  * \date   Mar 2018
  */

#include "maxwell.h"


/**
 ********************************************************************************************************************************************
 * \brief   Constructor of the maxwell class
 *
 *          Initialises various variables and constants according to conditions ans values given in yaml file.
 *
 *          Also calculates some quantities such as direction of H field according to directions specified of
 *          wave vector 'k' and E field.
 *
 *          Also calls the function plane_wave_initialise(), in case plane-waves are required by the user.
 *
 * \param   _inputData is a reference to reader class object.
 * \param   _mpi is a reference to parallel class object.
 * \param   _gridData is a reference to grid class object.
 ********************************************************************************************************************************************
 */
maxwell::maxwell(reader &_inputData, parallel &_mpi, grid &_gridData): inputData(_inputData), mpi(_mpi), gridData(_gridData) {
    n_threads = inputData.n_threads;
    num_timesteps = inputData.num_timesteps;
    epsilon0 = 8.85*pow(10,-12);
    mew0 = M_PI*4*pow(10,-7);
    S = inputData.S;
    c = 1/pow(epsilon0*mew0,0.5);
    dt=S*inputData.dx/c;
    dx=inputData.dx;
    dy=inputData.dy;
    dz=inputData.dz;
    sigma=pow(10,-4)/0.5;
    init=inputData.off_pw;
    lambda=10*dx;
    usePW = inputData.usePW;
    usePML = inputData.usePML;
    usePS = inputData.usePS;

    if (usePML){
        a_max=inputData.a_max;
        k_max=inputData.kappa;
        m_pml=inputData.sigma_m;
        ma_pml=inputData.a_m;
        sigma_optimal=0.8*(m_pml+1)/dx/pow(mew0/epsilon0,0.5);
        sigma_max=sigma_optimal*inputData.sigma;
        d_pml=inputData.d_pml;
    }
    if (usePW){
        E_x_dir = inputData.Ex_pw;
        E_y_dir = inputData.Ey_pw;
        E_z_dir = inputData.Ez_pw;
        kx = inputData.kx_pw;
        ky = inputData.ky_pw;
        kz = inputData.kz_pw;
        //Calculate H using cross product between k and E
        H_x_dir = ky*E_z_dir - kz*E_y_dir;
        H_y_dir = kz*E_x_dir - kx*E_z_dir;
        H_z_dir = kx*E_y_dir - ky*E_x_dir;
        //Now normalise direction of H
        double H_abs = pow(pow(H_x_dir,2)+pow(H_y_dir,2)+pow(H_z_dir,2), 0.5);
        H_x_dir = H_x_dir/H_abs;
        H_y_dir = H_y_dir/H_abs;
        H_z_dir = H_z_dir/H_abs;


        // double theta = M_PI/2.0;
        // H_x_dir = sin(theta);
        // H_z_dir = -cos(theta);
        // E_y_dir = 1;
        // kx = cos(theta);
        // kz = sin(theta);
        // ky = 0;

        // H_x_dir = -1.0/sqrt(6);
        // H_y_dir = 2.0/sqrt(6);
        // H_z_dir = -1.0/sqrt(6);
        // E_x_dir = 1.0/sqrt(2);
        // E_y_dir = 0;
        // E_z_dir = -1.0/sqrt(2);
        // kx = 1.0/sqrt(3);
        // ky = 1.0/sqrt(3);
        // kz = 1.0/sqrt(3);

        // H_x_dir = 0.0;
        // H_y_dir = 0.0;
        // H_z_dir = -1.0;
        // E_x_dir = 1.0;
        // E_y_dir = 0.0;
        // E_z_dir = 0.0;
        // kx = 0.0;
        // ky = 1.0;
        // kz = 0.0;

        p1 = inputData.x1_pw;
        q1 = inputData.y1_pw;
        r1 = inputData.z1_pw;
        p2 = inputData.x2_pw;
        q2 = inputData.y2_pw;
        r2 = inputData.z2_pw;
        plane_wave_initialise();
    }
    if (usePS){
        amp_source=inputData.amp_source;
        omega_source=inputData.omega_source;
        Nx_source=inputData.Nx_source;
        Ny_source=inputData.Ny_source;
        Nz_source=inputData.Nz_source;
    }
}


/**
 ********************************************************************************************************************************************
 * \brief   A function which decides the routine to call depending on the choice provided by user
 *          in the yaml file.
 ********************************************************************************************************************************************
 */
void maxwell::solve(){
    if (gridData.isPlanar){
        if (usePML){
            solve_planar_pml();
        }
        else{
            solve_planar();
        }
    }
    else{
        if (usePML){
            solve3d_pml();
        }
        else{
            solve3d();
        }
    }
}






/**
 ********************************************************************************************************************************************
 * \brief   This function edits the curl array provided to it adding the plane-wave values at the boundary
 *          of total field and scattered field. The function also takes timestep number to calculate the
 *          value of plane-wave theoretically.
 *
 *          It is enough to add the values of plane wave at the boundary of an imaginary total field cube.
 *          This keeps the plane-wave sustained throughout the cube.
 *
 *\param    curl is the poniter to curl array which needs to edited so that it contains the plane-wave values.
 *\param    timestep is the number of time-step at which the solution is currently at.
 ********************************************************************************************************************************************
 */
void maxwell::plane_wave_execute(vfield *curl, int timestep){
    blitz::TinyVector<int,3> temp;


    if (gridData.isPlanar){

        if (curl->isFaceCentered){
            // Changes on YZ plane
            if ((p1-1)>=gridData.colloq_start_index_x && (p1-1)<=gridData.colloq_end_index_x && start_z!=-1){
                for (int k=start_z; k<=end_z; k++){
                    temp = (p1-1),0,k;
                    temp = global_to_local(temp);
                    curl->Vz->F(temp) = curl->Vz->F(temp) + exp(-pow(((p1*dx*kx + k*dz*kz) + init*dx -(timestep-0.5)*c*dt),2)/sigma)*E_y_dir/dx;
                }
            }
            if (p2>=gridData.colloq_start_index_x && p2<=gridData.colloq_end_index_x && start_z!=-1){
                for (int k=start_z; k<=end_z; k++){
                    temp = p2,0,k;
                    temp = global_to_local(temp);
                    curl->Vz->F(temp) = curl->Vz->F(temp) - exp(-pow(((p2*dx*kx + k*dz*kz) + init*dx -(timestep-0.5)*c*dt),2)/sigma)*E_y_dir/dx;
                }
            }



            // Changes on XY plane
            if ((r1-1)>=gridData.colloq_start_index_z && (r1-1)<=gridData.colloq_end_index_z && start_x!=-1){
                for (int i=start_x; i<=end_x; i++){
                    temp = i,0,(r1-1);
                    temp = global_to_local(temp);
                    curl->Vx->F(temp) = curl->Vx->F(temp) - exp(-pow(((i*dx*kx + r1*dz*kz) + init*dx -(timestep-0.5)*c*dt),2)/sigma)*E_y_dir/dz;
                }
            }
            if (r2>=gridData.colloq_start_index_z && r2<=gridData.colloq_end_index_z && start_x!=-1){
                for (int i=start_x; i<=end_x; i++){
                    temp = i,0,r2;
                    temp = global_to_local(temp);
                    curl->Vx->F(temp) = curl->Vx->F(temp) + exp(-pow(((i*dx*kx + r2*dz*kz) + init*dx -(timestep-0.5)*c*dt),2)/sigma)*E_y_dir/dz;
                }
            }
        }
        else{
            // *************  if curl is not face-centered  *****************

            // Changes on YZ plane
            if (p1>=gridData.colloq_start_index_x && p1<=gridData.colloq_end_index_x && start_z!=-1){
                for (int k=start_z; k<=end_z; k++){
                    temp = p1,0,k;
                    temp = global_to_local(temp);
                    curl->Vy->F(temp) = curl->Vy->F(temp) + exp(-pow((((p1-0.5)*dx*kx + k*dz*kz) + init*dx -(timestep)*c*dt),2)/sigma)*(-1)*H_z_dir/dx/c/mew0;
                }
            }
            if (p2>=gridData.colloq_start_index_x && p2<=gridData.colloq_end_index_x && start_z!=-1){
                for (int k=start_z; k<=end_z; k++){
                    temp = p2,0,k;
                    temp = global_to_local(temp);
                    curl->Vy->F(temp) = curl->Vy->F(temp) - exp(-pow((((p2+0.5)*dx*kx + k*dz*kz) + init*dx -(timestep)*c*dt),2)/sigma)*(-1)*H_z_dir/dx/c/mew0;
                }
            }



            // Changes on XY plane
            if (r1>=gridData.colloq_start_index_z && r1<=gridData.colloq_end_index_z && start_x!=-1){
                for (int i=start_x; i<=end_x; i++){
                    temp = i,0,r1;
                    temp = global_to_local(temp);
                    curl->Vy->F(temp) = curl->Vy->F(temp) - exp(-pow(((i*dx*kx + (r1-0.5)*dz*kz) + init*dx -(timestep)*c*dt),2)/sigma)*(-1)*H_x_dir/dz/c/mew0;
                }
            }
            if (r2>=gridData.colloq_start_index_z && r2<=gridData.colloq_end_index_z && start_x!=-1){
                for (int i=start_x; i<=end_x; i++){
                    temp = i,0,r2;
                    temp = global_to_local(temp);
                    curl->Vy->F(temp) = curl->Vy->F(temp) + exp(-pow(((i*dx*kx + (r2+0.5)*dz*kz) + init*dx -(timestep)*c*dt),2)/sigma)*(-1)*H_x_dir/dz/c/mew0;
                }
            }
        }



    }


    else{
        // **********************************************************************
        // ************************  For 3-D plane wave  ************************
        // **********************************************************************
        if (curl->isFaceCentered){
            // Changes on YZ plane
            if ((p1-1)>=gridData.colloq_start_index_x && (p1-1)<=gridData.colloq_end_index_x && start_y!=-1 && start_z!=-1){
                for (int j=start_y; j<=end_y-1; j++){
                    for (int k=start_z; k<=end_z; k++){
                        temp = (p1-1),j,k;
                        temp = global_to_local(temp);
                        curl->Vz->F(temp) = curl->Vz->F(temp) - exp(-pow(((p1*dx*kx + (j+0.5)*dy*ky + k*dz*kz) + init*dx -(timestep-0.5)*c*dt),2)/sigma)*E_y_dir/dx;
                    }
                }
                for (int j=start_y; j<=end_y; j++){
                    for (int k=start_z; k<=end_z-1; k++){
                        temp = (p1-1),j,k;
                        temp = global_to_local(temp);
                        curl->Vy->F(temp) = curl->Vy->F(temp) + exp(-pow(((p1*dx*kx + j*dy*ky + (k+0.5)*dz*kz) + init*dx -(timestep-0.5)*c*dt),2)/sigma)*E_z_dir/dx;
                    }
                }
            }
            if (p2>=gridData.colloq_start_index_x && p2<=gridData.colloq_end_index_x && start_y!=-1 && start_z!=-1){
                for (int j=start_y; j<=end_y-1; j++){
                    for (int k=start_z; k<=end_z; k++){
                        temp = p2,j,k;
                        temp = global_to_local(temp);
                        curl->Vz->F(temp) = curl->Vz->F(temp) + exp(-pow(((p2*dx*kx + (j+0.5)*dy*ky + k*dz*kz) + init*dx -(timestep-0.5)*c*dt),2)/sigma)*E_y_dir/dx;
                    }
                }
                for (int j=start_y; j<=end_y; j++){
                    for (int k=start_z; k<=end_z-1; k++){
                        temp = p2,j,k;
                        temp = global_to_local(temp);
                        curl->Vy->F(temp) = curl->Vy->F(temp) - exp(-pow(((p2*dx*kx + j*dy*ky + (k+0.5)*dz*kz) + init*dx -(timestep-0.5)*c*dt),2)/sigma)*E_z_dir/dx;
                    }
                }
            }


            // Changes on XZ plane
            if ((q1-1)>=gridData.colloq_start_index_y && (q1-1)<=gridData.colloq_end_index_y && start_x!=-1 && start_z!=-1){
                for (int i=start_x; i<=end_x-1; i++){
                    for (int k=start_z; k<=end_z; k++){
                        temp = i,(q1-1),k;
                        temp = global_to_local(temp);
                        curl->Vz->F(temp) = curl->Vz->F(temp) + exp(-pow((((i+0.5)*dx*kx + q1*dy*ky + k*dz*kz) + init*dx -(timestep-0.5)*c*dt),2)/sigma)*E_x_dir/dy;
                    }
                }
                for (int i=start_x; i<=end_x; i++){
                    for (int k=start_z; k<=end_z-1; k++){
                        temp = i,(q1-1),k;
                        temp = global_to_local(temp);
                        curl->Vx->F(temp) = curl->Vx->F(temp) - exp(-pow(((i*dx*kx + q1*dy*ky + (k+0.5)*dz*kz) + init*dx -(timestep-0.5)*c*dt),2)/sigma)*E_z_dir/dy;
                    }
                }
            }
            if (q2>=gridData.colloq_start_index_y && q2<=gridData.colloq_end_index_y && start_x!=-1 && start_z!=-1){
                for (int i=start_x; i<=end_x-1; i++){
                    for (int k=start_z; k<=end_z; k++){
                        temp = i,q2,k;
                        temp = global_to_local(temp);
                        curl->Vz->F(temp) = curl->Vz->F(temp) - exp(-pow((((i+0.5)*dx*kx + q2*dy*ky + k*dz*kz) + init*dx -(timestep-0.5)*c*dt),2)/sigma)*E_x_dir/dy;
                    }
                }
                for (int i=start_x; i<=end_x; i++){
                    for (int k=start_z; k<=end_z-1; k++){
                        temp = i,q2,k;
                        temp = global_to_local(temp);
                        curl->Vx->F(temp) = curl->Vx->F(temp) + exp(-pow(((i*dx*kx + q2*dy*ky + (k+0.5)*dz*kz) + init*dx -(timestep-0.5)*c*dt),2)/sigma)*E_z_dir/dy;
                    }
                }
            }



            // Changes on XY plane
            if ((r1-1)>=gridData.colloq_start_index_z && (r1-1)<=gridData.colloq_end_index_z && start_x!=-1 && start_y!=-1){
                for (int i=start_x; i<=end_x-1; i++){
                    for (int j=start_y; j<=end_y; j++){
                        temp = i,j,(r1-1);
                        temp = global_to_local(temp);
                        curl->Vy->F(temp) = curl->Vy->F(temp) - exp(-pow((((i+0.5)*dx*kx + j*dy*ky + r1*dz*kz) + init*dx -(timestep-0.5)*c*dt),2)/sigma)*E_x_dir/dz;
                    }
                }
                for (int i=start_x; i<=end_x; i++){
                    for (int j=start_y; j<=end_y-1; j++){
                        temp = i,j,(r1-1);
                        temp = global_to_local(temp);
                        curl->Vx->F(temp) = curl->Vx->F(temp) + exp(-pow(((i*dx*kx + (j+0.5)*dy*ky + r1*dz*kz) + init*dx -(timestep-0.5)*c*dt),2)/sigma)*E_y_dir/dz;
                    }
                }
            }
            if (r2>=gridData.colloq_start_index_z && r2<=gridData.colloq_end_index_z && start_x!=-1 && start_y!=-1){
                for (int i=start_x; i<=end_x-1; i++){
                    for (int j=start_y; j<=end_y; j++){
                        temp = i,j,r2;
                        temp = global_to_local(temp);
                        curl->Vy->F(temp) = curl->Vy->F(temp) + exp(-pow((((i+0.5)*dx*kx + j*dy*ky + r2*dz*kz) + init*dx -(timestep-0.5)*c*dt),2)/sigma)*E_x_dir/dz;
                    }
                }
                for (int i=start_x; i<=end_x; i++){
                    for (int j=start_y; j<=end_y-1; j++){
                        temp = i,j,r2;
                        temp = global_to_local(temp);
                        curl->Vx->F(temp) = curl->Vx->F(temp) - exp(-pow(((i*dx*kx + (j+0.5)*dy*ky + r2*dz*kz) + init*dx -(timestep-0.5)*c*dt),2)/sigma)*E_y_dir/dz;
                    }
                }
            }
        }
        else{
            // *************  if curl is not face-centered  *****************

            // Changes on YZ plane
            if (p1>=gridData.colloq_start_index_x && p1<=gridData.colloq_end_index_x && start_y!=-1 && start_z!=-1){
                for (int j=start_y; j<=end_y-1; j++){
                    for (int k=start_z; k<=end_z; k++){
                        temp = p1,j,k;
                        temp = global_to_local(temp);
                        curl->Vy->F(temp) = curl->Vy->F(temp) + exp(-pow((((p1-0.5)*dx*kx + (j+0.5)*dy*ky + k*dz*kz) + init*dx -(timestep)*c*dt),2)/sigma)*H_z_dir/dx/c/mew0;
                    }
                }
                for (int j=start_y; j<=end_y; j++){
                    for (int k=start_z; k<=end_z-1; k++){
                        temp = p1,j,k;
                        temp = global_to_local(temp);
                        curl->Vz->F(temp) = curl->Vz->F(temp) - exp(-pow((((p1-0.5)*dx*kx + j*dy*ky + (k+0.5)*dz*kz) + init*dx -(timestep)*c*dt),2)/sigma)*H_y_dir/dx/c/mew0;
                    }
                }
            }
            if (p2>=gridData.colloq_start_index_x && p2<=gridData.colloq_end_index_x && start_y!=-1 && start_z!=-1){
                for (int j=start_y; j<=end_y-1; j++){
                    for (int k=start_z; k<=end_z; k++){
                        temp = p2,j,k;
                        temp = global_to_local(temp);
                        curl->Vy->F(temp) = curl->Vy->F(temp) - exp(-pow((((p2+0.5)*dx*kx + (j+0.5)*dy*ky + k*dz*kz) + init*dx -(timestep)*c*dt),2)/sigma)*H_z_dir/dx/c/mew0;
                    }
                }
                for (int j=start_y; j<=end_y; j++){
                    for (int k=start_z; k<=end_z-1; k++){
                        temp = p2,j,k;
                        temp = global_to_local(temp);
                        curl->Vz->F(temp) = curl->Vz->F(temp) + exp(-pow((((p2+0.5)*dx*kx + j*dy*ky + (k+0.5)*dz*kz) + init*dx -(timestep)*c*dt),2)/sigma)*H_y_dir/dx/c/mew0;
                    }
                }
            }


            // Changes on XZ plane
            if (q1>=gridData.colloq_start_index_y && q1<=gridData.colloq_end_index_y && start_x!=-1 && start_z!=-1){
                for (int i=start_x; i<=end_x-1; i++){
                    for (int k=start_z; k<=end_z; k++){
                        temp = i,q1,k;
                        temp = global_to_local(temp);
                        curl->Vx->F(temp) = curl->Vx->F(temp) - exp(-pow((((i+0.5)*dx*kx + (q1-0.5)*dy*ky + k*dz*kz) + init*dx -(timestep)*c*dt),2)/sigma)*H_z_dir/dy/c/mew0;
                    }
                }
                for (int i=start_x; i<=end_x; i++){
                    for (int k=start_z; k<=end_z-1; k++){
                        temp = i,q1,k;
                        temp = global_to_local(temp);
                        curl->Vz->F(temp) = curl->Vz->F(temp) + exp(-pow(((i*dx*kx + (q1-0.5)*dy*ky + (k+0.5)*dz*kz) + init*dx -(timestep)*c*dt),2)/sigma)*H_x_dir/dy/c/mew0;
                    }
                }
            }
            if (q2>=gridData.colloq_start_index_y && q2<=gridData.colloq_end_index_y && start_x!=-1 && start_z!=-1){
                for (int i=start_x; i<=end_x-1; i++){
                    for (int k=start_z; k<=end_z; k++){
                        temp = i,q2,k;
                        temp = global_to_local(temp);
                        curl->Vx->F(temp) = curl->Vx->F(temp) + exp(-pow((((i+0.5)*dx*kx + (q2+0.5)*dy*ky + k*dz*kz) + init*dx -(timestep)*c*dt),2)/sigma)*H_z_dir/dy/c/mew0;
                    }
                }
                for (int i=start_x; i<=end_x; i++){
                    for (int k=start_z; k<=end_z-1; k++){
                        temp = i,q2,k;
                        temp = global_to_local(temp);
                        curl->Vz->F(temp) = curl->Vz->F(temp) - exp(-pow(((i*dx*kx + (q2+0.5)*dy*ky + (k+0.5)*dz*kz) + init*dx -(timestep)*c*dt),2)/sigma)*H_x_dir/dy/c/mew0;
                    }
                }
            }



            // Changes on XY plane
            if (r1>=gridData.colloq_start_index_z && r1<=gridData.colloq_end_index_z && start_x!=-1 && start_y!=-1){
                for (int i=start_x; i<=end_x-1; i++){
                    for (int j=start_y; j<=end_y; j++){
                        temp = i,j,r1;
                        temp = global_to_local(temp);
                        curl->Vx->F(temp) = curl->Vx->F(temp) + exp(-pow((((i+0.5)*dx*kx + j*dy*ky + (r1-0.5)*dz*kz) + init*dx -(timestep)*c*dt),2)/sigma)*H_y_dir/dz/c/mew0;
                    }
                }
                for (int i=start_x; i<=end_x; i++){
                    for (int j=start_y; j<=end_y-1; j++){
                        temp = i,j,r1;
                        temp = global_to_local(temp);
                        curl->Vy->F(temp) = curl->Vy->F(temp) - exp(-pow(((i*dx*kx + (j+0.5)*dy*ky + (r1-0.5)*dz*kz) + init*dx -(timestep)*c*dt),2)/sigma)*H_x_dir/dz/c/mew0;
                    }
                }
            }
            if (r2>=gridData.colloq_start_index_z && r2<=gridData.colloq_end_index_z && start_x!=-1 && start_y!=-1){
                for (int i=start_x; i<=end_x-1; i++){
                    for (int j=start_y; j<=end_y; j++){
                        temp = i,j,r2;
                        temp = global_to_local(temp);
                        curl->Vx->F(temp) = curl->Vx->F(temp) - exp(-pow((((i+0.5)*dx*kx + j*dy*ky + (r2+0.5)*dz*kz) + init*dx -(timestep)*c*dt),2)/sigma)*H_y_dir/dz/c/mew0;
                    }
                }
                for (int i=start_x; i<=end_x; i++){
                    for (int j=start_y; j<=end_y-1; j++){
                        temp = i,j,r2;
                        temp = global_to_local(temp);
                        curl->Vy->F(temp) = curl->Vy->F(temp) + exp(-pow(((i*dx*kx + (j+0.5)*dy*ky + (r2+0.5)*dz*kz) + init*dx -(timestep)*c*dt),2)/sigma)*H_x_dir/dz/c/mew0;
                    }
                }
            }
        }
    }  //End of else
}


/**
 ********************************************************************************************************************************************
 * \brief   This function is the overload to plane_wave_execute(vfield *curl, int timestep).
 *
 *          It does the exact same thing but because 3D solutions require two curl arrays instead of one, it
 *          takes both of them as input and adds the plane-wave values accordingly.
 *
 *\param    curl1 is the poniter to first curl array which needs to edited so that it contains the plane-wave values.
 *\param    curl2 is the poniter to second curl array which needs to edited so that it contains the plane-wave values.
 *\param    timestep is the number of time-step at which the solution is currently at.
 ********************************************************************************************************************************************
 */
void maxwell::plane_wave_execute(vfield *curl1, vfield *curl2, int timestep){
    blitz::TinyVector<int,3> temp;


    if (gridData.isPlanar){ //For 2-D plane wave
        if (curl1->isFaceCentered){
            // Changes on YZ plane
            if ((p1-1)>=gridData.colloq_start_index_x && (p1-1)<=gridData.colloq_end_index_x && start_z!=-1){
                for (int k=start_z; k<=end_z; k++){
                    temp = (p1-1),0,k;
                    temp = global_to_local(temp);
                    curl1->Vz->F(temp) = curl1->Vz->F(temp) + exp(-pow(((p1*dx*kx + k*dz*kz) + init*dx -(timestep-0.5)*c*dt),2)/sigma)*E_y_dir/dx;
                }

            }
            if (p2>=gridData.colloq_start_index_x && p2<=gridData.colloq_end_index_x && start_z!=-1){
                for (int k=start_z; k<=end_z; k++){
                    temp = p2,0,k;
                    temp = global_to_local(temp);
                    curl1->Vz->F(temp) = curl1->Vz->F(temp) - exp(-pow(((p2*dx*kx + k*dz*kz) + init*dx -(timestep-0.5)*c*dt),2)/sigma)*E_y_dir/dx;
                }

            }



            // Changes on XY plane
            if ((r1-1)>=gridData.colloq_start_index_z && (r1-1)<=gridData.colloq_end_index_z && start_x!=-1){
                for (int i=start_x; i<=end_x; i++){
                    temp = i,0,(r1-1);
                    temp = global_to_local(temp);
                    curl2->Vx->F(temp) = curl2->Vx->F(temp) + exp(-pow(((i*dx*kx + r1*dz*kz) + init*dx -(timestep-0.5)*c*dt),2)/sigma)*E_y_dir/dz;

                }
            }
            if (r2>=gridData.colloq_start_index_z && r2<=gridData.colloq_end_index_z && start_x!=-1){
                for (int i=start_x; i<=end_x; i++){
                    temp = i,0,r2;
                    temp = global_to_local(temp);
                    curl2->Vx->F(temp) = curl2->Vx->F(temp) - exp(-pow(((i*dx*kx + r2*dz*kz) + init*dx -(timestep-0.5)*c*dt),2)/sigma)*E_y_dir/dz;
                }
            }
        }
        else{
            // *************  if curl is not face-centered  *****************

            // Changes on YZ plane
            if (p1>=gridData.colloq_start_index_x && p1<=gridData.colloq_end_index_x && start_z!=-1){
                for (int k=start_z; k<=end_z; k++){
                    temp = p1,0,k;
                    temp = global_to_local(temp);
                    curl2->Vy->F(temp) = curl2->Vy->F(temp) - exp(-pow((((p1-0.5)*dx*kx + k*dz*kz) + init*dx -(timestep)*c*dt),2)/sigma)*(-1)*H_z_dir/dx/c/mew0;
                }
            }
            if (p2>=gridData.colloq_start_index_x && p2<=gridData.colloq_end_index_x && start_z!=-1){
                for (int k=start_z; k<=end_z; k++){
                    temp = p2,0,k;
                    temp = global_to_local(temp);
                    curl2->Vy->F(temp) = curl2->Vy->F(temp) + exp(-pow((((p2+0.5)*dx*kx + k*dz*kz) + init*dx -(timestep)*c*dt),2)/sigma)*(-1)*H_z_dir/dx/c/mew0;
                }
            }



            // Changes on XY plane
            if (r1>=gridData.colloq_start_index_z && r1<=gridData.colloq_end_index_z && start_x!=-1){
                for (int i=start_x; i<=end_x; i++){
                    temp = i,0,r1;
                    temp = global_to_local(temp);
                    curl1->Vy->F(temp) = curl1->Vy->F(temp) - exp(-pow(((i*dx*kx + (r1-0.5)*dz*kz) + init*dx -(timestep)*c*dt),2)/sigma)*(-1)*H_x_dir/dz/c/mew0;
                }
            }
            if (r2>=gridData.colloq_start_index_z && r2<=gridData.colloq_end_index_z && start_x!=-1){
                for (int i=start_x; i<=end_x; i++){
                    temp = i,0,r2;
                    temp = global_to_local(temp);
                    curl1->Vy->F(temp) = curl1->Vy->F(temp) + exp(-pow(((i*dx*kx + (r2+0.5)*dz*kz) + init*dx -(timestep)*c*dt),2)/sigma)*(-1)*H_x_dir/dz/c/mew0;
                }
            }
        }
    }
    else{
        // **********************************************************************
        // ************************  For 3-D plane wave  ************************
        // **********************************************************************
        if (curl1->isFaceCentered){
            // Changes on YZ plane
            if ((p1-1)>=gridData.colloq_start_index_x && (p1-1)<=gridData.colloq_end_index_x && start_y!=-1 && start_z!=-1){
                for (int j=start_y; j<=end_y-1; j++){
                    for (int k=start_z; k<=end_z; k++){
                        temp = (p1-1),j,k;
                        temp = global_to_local(temp);
                        curl1->Vz->F(temp) = curl1->Vz->F(temp) - exp(-pow(((p1*dx*kx + (j+0.5)*dy*ky + k*dz*kz) + init*dx -(timestep-0.5)*c*dt),2)/sigma)*E_y_dir/dx;
                    }
                }
                for (int j=start_y; j<=end_y; j++){
                    for (int k=start_z; k<=end_z-1; k++){
                        temp = (p1-1),j,k;
                        temp = global_to_local(temp);
                        curl2->Vy->F(temp) = curl2->Vy->F(temp) - exp(-pow(((p1*dx*kx + j*dy*ky + (k+0.5)*dz*kz) + init*dx -(timestep-0.5)*c*dt),2)/sigma)*E_z_dir/dx;
                    }
                }
            }
            if (p2>=gridData.colloq_start_index_x && p2<=gridData.colloq_end_index_x && start_y!=-1 && start_z!=-1){
                for (int j=start_y; j<=end_y-1; j++){
                    for (int k=start_z; k<=end_z; k++){
                        temp = p2,j,k;
                        temp = global_to_local(temp);
                        curl1->Vz->F(temp) = curl1->Vz->F(temp) + exp(-pow(((p2*dx*kx + (j+0.5)*dy*ky + k*dz*kz) + init*dx -(timestep-0.5)*c*dt),2)/sigma)*E_y_dir/dx;
                    }
                }
                for (int j=start_y; j<=end_y; j++){
                    for (int k=start_z; k<=end_z-1; k++){
                        temp = p2,j,k;
                        temp = global_to_local(temp);
                        curl2->Vy->F(temp) = curl2->Vy->F(temp) + exp(-pow(((p2*dx*kx + j*dy*ky + (k+0.5)*dz*kz) + init*dx -(timestep-0.5)*c*dt),2)/sigma)*E_z_dir/dx;
                    }
                }
            }


            // Changes on XZ plane
            if ((q1-1)>=gridData.colloq_start_index_y && (q1-1)<=gridData.colloq_end_index_y && start_x!=-1 && start_z!=-1){
                for (int i=start_x; i<=end_x-1; i++){
                    for (int k=start_z; k<=end_z; k++){
                        temp = i,(q1-1),k;
                        temp = global_to_local(temp);
                        curl2->Vz->F(temp) = curl2->Vz->F(temp) - exp(-pow((((i+0.5)*dx*kx + q1*dy*ky + k*dz*kz) + init*dx -(timestep-0.5)*c*dt),2)/sigma)*E_x_dir/dy;
                    }
                }
                for (int i=start_x; i<=end_x; i++){
                    for (int k=start_z; k<=end_z-1; k++){
                        temp = i,(q1-1),k;
                        temp = global_to_local(temp);
                        curl1->Vx->F(temp) = curl1->Vx->F(temp) - exp(-pow(((i*dx*kx + q1*dy*ky + (k+0.5)*dz*kz) + init*dx -(timestep-0.5)*c*dt),2)/sigma)*E_z_dir/dy;
                    }
                }
            }
            if (q2>=gridData.colloq_start_index_y && q2<=gridData.colloq_end_index_y && start_x!=-1 && start_z!=-1){
                for (int i=start_x; i<=end_x-1; i++){
                    for (int k=start_z; k<=end_z; k++){
                        temp = i,q2,k;
                        temp = global_to_local(temp);
                        curl2->Vz->F(temp) = curl2->Vz->F(temp) + exp(-pow((((i+0.5)*dx*kx + q2*dy*ky + k*dz*kz) + init*dx -(timestep-0.5)*c*dt),2)/sigma)*E_x_dir/dy;
                    }
                }
                for (int i=start_x; i<=end_x; i++){
                    for (int k=start_z; k<=end_z-1; k++){
                        temp = i,q2,k;
                        temp = global_to_local(temp);
                        curl1->Vx->F(temp) = curl1->Vx->F(temp) + exp(-pow(((i*dx*kx + q2*dy*ky + (k+0.5)*dz*kz) + init*dx -(timestep-0.5)*c*dt),2)/sigma)*E_z_dir/dy;
                    }
                }
            }



            // Changes on XY plane
            if ((r1-1)>=gridData.colloq_start_index_z && (r1-1)<=gridData.colloq_end_index_z && start_x!=-1 && start_y!=-1){
                for (int i=start_x; i<=end_x-1; i++){
                    for (int j=start_y; j<=end_y; j++){
                        temp = i,j,(r1-1);
                        temp = global_to_local(temp);
                        curl1->Vy->F(temp) = curl1->Vy->F(temp) - exp(-pow((((i+0.5)*dx*kx + j*dy*ky + r1*dz*kz) + init*dx -(timestep-0.5)*c*dt),2)/sigma)*E_x_dir/dz;
                    }
                }
                for (int i=start_x; i<=end_x; i++){
                    for (int j=start_y; j<=end_y-1; j++){
                        temp = i,j,(r1-1);
                        temp = global_to_local(temp);
                        curl2->Vx->F(temp) = curl2->Vx->F(temp) - exp(-pow(((i*dx*kx + (j+0.5)*dy*ky + r1*dz*kz) + init*dx -(timestep-0.5)*c*dt),2)/sigma)*E_y_dir/dz;
                    }
                }
            }
            if (r2>=gridData.colloq_start_index_z && r2<=gridData.colloq_end_index_z && start_x!=-1 && start_y!=-1){
                for (int i=start_x; i<=end_x-1; i++){
                    for (int j=start_y; j<=end_y; j++){
                        temp = i,j,r2;
                        temp = global_to_local(temp);
                        curl1->Vy->F(temp) = curl1->Vy->F(temp) + exp(-pow((((i+0.5)*dx*kx + j*dy*ky + r2*dz*kz) + init*dx -(timestep-0.5)*c*dt),2)/sigma)*E_x_dir/dz;
                    }
                }
                for (int i=start_x; i<=end_x; i++){
                    for (int j=start_y; j<=end_y-1; j++){
                        temp = i,j,r2;
                        temp = global_to_local(temp);
                        curl2->Vx->F(temp) = curl2->Vx->F(temp) + exp(-pow(((i*dx*kx + (j+0.5)*dy*ky + r2*dz*kz) + init*dx -(timestep-0.5)*c*dt),2)/sigma)*E_y_dir/dz;
                    }
                }
            }
        }
        else{
            // *************  if curl is not face-centered  *****************

            // Changes on YZ plane
            if (p1>=gridData.colloq_start_index_x && p1<=gridData.colloq_end_index_x && start_y!=-1 && start_z!=-1){
                for (int j=start_y; j<=end_y-1; j++){
                    for (int k=start_z; k<=end_z; k++){
                        temp = p1,j,k;
                        temp = global_to_local(temp);
                        curl2->Vy->F(temp) = curl2->Vy->F(temp) - exp(-pow((((p1-0.5)*dx*kx + (j+0.5)*dy*ky + k*dz*kz) + init*dx -(timestep)*c*dt),2)/sigma)*H_z_dir/dx/c/mew0;
                    }
                }
                for (int j=start_y; j<=end_y; j++){
                    for (int k=start_z; k<=end_z-1; k++){
                        temp = p1,j,k;
                        temp = global_to_local(temp);
                        curl1->Vz->F(temp) = curl1->Vz->F(temp) - exp(-pow((((p1-0.5)*dx*kx + j*dy*ky + (k+0.5)*dz*kz) + init*dx -(timestep)*c*dt),2)/sigma)*H_y_dir/dx/c/mew0;
                    }
                }
            }
            if (p2>=gridData.colloq_start_index_x && p2<=gridData.colloq_end_index_x && start_y!=-1 && start_z!=-1){
                for (int j=start_y; j<=end_y-1; j++){
                    for (int k=start_z; k<=end_z; k++){
                        temp = p2,j,k;
                        temp = global_to_local(temp);
                        curl2->Vy->F(temp) = curl2->Vy->F(temp) + exp(-pow((((p2+0.5)*dx*kx + (j+0.5)*dy*ky + k*dz*kz) + init*dx -(timestep)*c*dt),2)/sigma)*H_z_dir/dx/c/mew0;
                    }
                }
                for (int j=start_y; j<=end_y; j++){
                    for (int k=start_z; k<=end_z-1; k++){
                        temp = p2,j,k;
                        temp = global_to_local(temp);
                        curl1->Vz->F(temp) = curl1->Vz->F(temp) + exp(-pow((((p2+0.5)*dx*kx + j*dy*ky + (k+0.5)*dz*kz) + init*dx -(timestep)*c*dt),2)/sigma)*H_y_dir/dx/c/mew0;
                    }
                }
            }


            // Changes on XZ plane
            if (q1>=gridData.colloq_start_index_y && q1<=gridData.colloq_end_index_y && start_x!=-1 && start_z!=-1){
                for (int i=start_x; i<=end_x-1; i++){
                    for (int k=start_z; k<=end_z; k++){
                        temp = i,q1,k;
                        temp = global_to_local(temp);
                        curl1->Vx->F(temp) = curl1->Vx->F(temp) - exp(-pow((((i+0.5)*dx*kx + (q1-0.5)*dy*ky + k*dz*kz) + init*dx -(timestep)*c*dt),2)/sigma)*H_z_dir/dy/c/mew0;
                    }
                }
                for (int i=start_x; i<=end_x; i++){
                    for (int k=start_z; k<=end_z-1; k++){
                        temp = i,q1,k;
                        temp = global_to_local(temp);
                        curl2->Vz->F(temp) = curl2->Vz->F(temp) - exp(-pow(((i*dx*kx + (q1-0.5)*dy*ky + (k+0.5)*dz*kz) + init*dx -(timestep)*c*dt),2)/sigma)*H_x_dir/dy/c/mew0;
                    }
                }
            }
            if (q2>=gridData.colloq_start_index_y && q2<=gridData.colloq_end_index_y && start_x!=-1 && start_z!=-1){
                for (int i=start_x; i<=end_x-1; i++){
                    for (int k=start_z; k<=end_z; k++){
                        temp = i,q2,k;
                        temp = global_to_local(temp);
                        curl1->Vx->F(temp) = curl1->Vx->F(temp) + exp(-pow((((i+0.5)*dx*kx + (q2+0.5)*dy*ky + k*dz*kz) + init*dx -(timestep)*c*dt),2)/sigma)*H_z_dir/dy/c/mew0;
                    }
                }
                for (int i=start_x; i<=end_x; i++){
                    for (int k=start_z; k<=end_z-1; k++){
                        temp = i,q2,k;
                        temp = global_to_local(temp);
                        curl2->Vz->F(temp) = curl2->Vz->F(temp) + exp(-pow(((i*dx*kx + (q2+0.5)*dy*ky + (k+0.5)*dz*kz) + init*dx -(timestep)*c*dt),2)/sigma)*H_x_dir/dy/c/mew0;
                    }
                }
            }



            // Changes on XY plane
            if (r1>=gridData.colloq_start_index_z && r1<=gridData.colloq_end_index_z && start_x!=-1 && start_y!=-1){
                for (int i=start_x; i<=end_x-1; i++){
                    for (int j=start_y; j<=end_y; j++){
                        temp = i,j,r1;
                        temp = global_to_local(temp);
                        curl2->Vx->F(temp) = curl2->Vx->F(temp) - exp(-pow((((i+0.5)*dx*kx + j*dy*ky + (r1-0.5)*dz*kz) + init*dx -(timestep)*c*dt),2)/sigma)*H_y_dir/dz/c/mew0;
                    }
                }
                for (int i=start_x; i<=end_x; i++){
                    for (int j=start_y; j<=end_y-1; j++){
                        temp = i,j,r1;
                        temp = global_to_local(temp);
                        curl1->Vy->F(temp) = curl1->Vy->F(temp) - exp(-pow(((i*dx*kx + (j+0.5)*dy*ky + (r1-0.5)*dz*kz) + init*dx -(timestep)*c*dt),2)/sigma)*H_x_dir/dz/c/mew0;
                    }
                }
            }
            if (r2>=gridData.colloq_start_index_z && r2<=gridData.colloq_end_index_z && start_x!=-1 && start_y!=-1){
                for (int i=start_x; i<=end_x-1; i++){
                    for (int j=start_y; j<=end_y; j++){
                        temp = i,j,r2;
                        temp = global_to_local(temp);
                        curl2->Vx->F(temp) = curl2->Vx->F(temp) + exp(-pow((((i+0.5)*dx*kx + j*dy*ky + (r2+0.5)*dz*kz) + init*dx -(timestep)*c*dt),2)/sigma)*H_y_dir/dz/c/mew0;
                    }
                }
                for (int i=start_x; i<=end_x; i++){
                    for (int j=start_y; j<=end_y-1; j++){
                        temp = i,j,r2;
                        temp = global_to_local(temp);
                        curl1->Vy->F(temp) = curl1->Vy->F(temp) + exp(-pow(((i*dx*kx + (j+0.5)*dy*ky + (r2+0.5)*dz*kz) + init*dx -(timestep)*c*dt),2)/sigma)*H_x_dir/dz/c/mew0;
                    }
                }
            }
        }
    }
}


/**
 ********************************************************************************************************************************************
 * \brief   This function initialises certain variables according to the grid distribution over different
 *          different processors. These variables are crucial to functioning of plane_wave_execute().
 *
 *          The variables are initialised once and then used throughout all the timesteps without change.
 *          Hence, the function is called once by the constructor if plane-wave source is reuqired.
 *
 ********************************************************************************************************************************************
 */
void maxwell::plane_wave_initialise(){

    start_x=-1, end_x=-1, start_y=-1, end_y=-1, start_z=-1, end_z=-1;

    if (p1>gridData.colloq_end_index_x){
        start_x = -1;
        end_x = -1;
    }
    else if(p1>=gridData.colloq_start_index_x){
        start_x = p1;
        if (p2>gridData.colloq_end_index_x){
            end_x = gridData.colloq_end_index_x;
        }
        else{
            end_x = p2;
        }
    }
    else{
        if(p2>gridData.colloq_end_index_x){
            start_x = gridData.colloq_start_index_x;
            end_x = gridData.colloq_end_index_x;
        }
        else if (p2>=gridData.colloq_start_index_x){
            start_x = gridData.colloq_start_index_x;
            end_x = p2;
        }
        else{
            start_x = -1;
            end_x = -1;
        }
    }



    if (q1>gridData.colloq_end_index_y){
        start_y = -1;
        end_y = -1;
    }
    else if(q1>=gridData.colloq_start_index_y){
        start_y = q1;
        if (q2>gridData.colloq_end_index_y){
            end_y = gridData.colloq_end_index_y;
        }
        else{
            end_y = q2;
        }
    }
    else{
        if(q2>gridData.colloq_end_index_y){
            start_y = gridData.colloq_start_index_y;
            end_y = gridData.colloq_end_index_y;
        }
        else if (q2>=gridData.colloq_start_index_y){
            start_y = gridData.colloq_start_index_y;
            end_y = q2;
        }
        else{
            start_y = -1;
            end_y = -1;
        }
    }

    if (r1>gridData.colloq_end_index_z){
        start_z = -1;
        end_z = -1;
    }
    else if(r1>=gridData.colloq_start_index_z){
        start_z = r1;
        if (r2>gridData.colloq_end_index_z){
            end_z = gridData.colloq_end_index_z;
        }
        else{
            end_z = r2;
        }
    }
    else{
        if(r2>gridData.colloq_end_index_z){
            start_z = gridData.colloq_start_index_z;
            end_z = gridData.colloq_end_index_z;
        }
        else if (r2>=gridData.colloq_start_index_z){
            start_z = gridData.colloq_start_index_z;
            end_z = r2;
        }
        else{
            start_z = -1;
            end_z = -1;
        }
    }


    // std::cout<<"Rank: "<<mpi.rank<<", X-start: "<<start_x<<std::endl;
    // std::cout<<"Rank: "<<mpi.rank<<", X-end: "<<end_x<<std::endl;
    // std::cout<<"Rank: "<<mpi.rank<<", Y-start: "<<start_y<<std::endl;
    // std::cout<<"Rank: "<<mpi.rank<<", Y-end: "<<end_y<<std::endl;
    // std::cout<<"Rank: "<<mpi.rank<<", Z-start: "<<start_z<<std::endl;
    // std::cout<<"Rank: "<<mpi.rank<<", Z-end: "<<end_z<<std::endl;
    //
    // std::cout<<"Rank: "<<mpi.rank<<", Colloq X-start: "<<gridData.colloq_start_index_x<<std::endl;
    // std::cout<<"Rank: "<<mpi.rank<<", Colloq X-end: "<<gridData.colloq_end_index_x<<std::endl;
    // std::cout<<"Rank: "<<mpi.rank<<", Colloq Y-start: "<<gridData.colloq_start_index_y<<std::endl;
    // std::cout<<"Rank: "<<mpi.rank<<", Colloq Y-end: "<<gridData.colloq_end_index_y<<std::endl;
    // std::cout<<"Rank: "<<mpi.rank<<", Colloq Z-start: "<<gridData.colloq_start_index_z<<std::endl;
    // std::cout<<"Rank: "<<mpi.rank<<", Colloq Z-end: "<<gridData.colloq_end_index_z<<std::endl;
}


/**
 ********************************************************************************************************************************************
 * \brief   This small function simply checks whether a given global index specified by parameter v lies
 *          within the MPI node calling the function. This is checked using grid data and staggered nature
 *          of the scalar field whose index is being questioned.
 *
 *\param    v is a three element Blitz TinyVector which conains the global index in question.
 *\param    xstag is a boolean which contains whether the field is question is staggered in X.
 *\param    ystag is a boolean which contains whether the field is question is staggered in Y.
 *
 ********************************************************************************************************************************************
 */
bool maxwell::check_global_limits(blitz::TinyVector<int,3> v, bool xstag, bool ystag){
    int Nx_e = gridData.colloq_end_index_x;
    int Ny_e = gridData.colloq_end_index_y;
    if (xstag){
        Nx_e = Nx_e - 1;
    }
    if (ystag && (!gridData.isPlanar)){
        Ny_e = Ny_e - 1;
    }
    if ((v(0)-gridData.colloq_start_index_x)>=0 && (v(1)-gridData.colloq_start_index_y)>=0){
        if ((v(0)-Nx_e)<=0 && (v(1)-Ny_e)<=0){
            return true;
        }
    }
    return false;
}


/**
 ********************************************************************************************************************************************
 * \brief   This function returns a local index to the MPI node of global index, if the result from
 *          function check_global_limits() is true. This is used to position the point source in MPI
 *          node as user will provide the index in global index terms but it will be finally required
 *          in local index terms for a specific MPI node.
 *
 *\param    glob is a three element Blitz TinyVector which conains the global index required to be
 *          converted to local index.
 *
 ********************************************************************************************************************************************
 */
blitz::TinyVector<int,3> maxwell::global_to_local(blitz::TinyVector<int,3> glob){
    blitz::TinyVector<int,3> loc = glob;
    loc(0) = glob(0) - gridData.colloq_start_index_x;
    loc(1) = glob(1) - gridData.colloq_start_index_y;
    return loc;
}


/**
 ********************************************************************************************************************************************
 * \brief   This function calculates and stores the polynomial graded profile of loss term sigma in
 *          variable 'w'. The function can calculate the grading profile in all directions even for
 *          rectangular grid. The direction is specified by parameter 'dim' with its value being 0, 1,
 *          and 2 for X, Y, and Z direction respectively. The function is used to calculate profile for
 *          both E and H field. The parameter 'is_E' is set to true if profile is needed for E field
 *          otherwise false.
 *
 *\param    dim specifes the dimension in which grading profile is to be calculated.
 *\param    is_E is set to true if profile is required for E field, otherwise it is set to false.
 *\param    w is the variable where the calculated grading profile will be stored.
 *
 ********************************************************************************************************************************************
 */
int maxwell::sigma_fn(int dim, bool is_E, blitz::Array<double,1> w){
  double delta=0;
  int N=0;
  int N_local=0;
  int offset=0;
  if (dim==0){
    delta=gridData.inputData.dx;
    N_local=gridData.local_colloq_x;
    N=gridData.inputData.Nx;
    offset=gridData.colloq_start_index_x;
  }
  else if (dim==1){
    delta=gridData.inputData.dy;
    N_local=gridData.local_colloq_y;
    N=gridData.inputData.Ny;
    offset=gridData.colloq_start_index_y;
  }
  else if (dim==2){
    delta=gridData.inputData.dz;
    N_local=gridData.local_colloq_z;
    N=gridData.inputData.Nz;
    offset=0;
  }
  else{
    std::cout<<"Warning: Dim should be 1, 2 or 3. Returning Default!"<<std::endl;
    w=0;
    return 0;
  }

  blitz::firstIndex ii;
  blitz::Array<double,1> temp;
  if (is_E){
    temp.resize(N_local);
    temp=ii+offset;
  }
  else{
    temp.resize(N_local-1);
    temp=ii+0.5+offset;
  }

  double d=d_pml/delta;
  w=0.0;
  w=(temp<=(d-1))*pow((((d-1)-temp)/(d-1)),m_pml)*sigma_max+(temp>(d-1))*w;
  w=(temp>=(N-d))*pow((temp-N+d)/(d-1),m_pml)*sigma_max+(temp<(N-d))*w;
  return 0;
}


/**
 ********************************************************************************************************************************************
 * \brief   Similar to sigma_fn() this function also calculates another quantity useful for PML application,
 *          specifically for strectched PML. The profile for this function is also polynomial. The function
 *          calculates the profile for all directions in rectangular grid.
 *          The input parameters are similar to sigma_fn() and serve same function.
 *
 *\param    dim specifes the dimension in which grading profile is to be calculated.
 *\param    is_E is set to true if profile is required for E field, otherwise it is set to false.
 *\param    w is the variable where the calculated grading profile will be stored.
 *
 ********************************************************************************************************************************************
 */
int maxwell::k_fn(int dim, bool is_E, blitz::Array<double,1> w){
    double delta=0;
    int N=0;
    int N_local=0;
    int offset=0;
    if (dim==0){
      delta=gridData.inputData.dx;
      N_local=gridData.local_colloq_x;
      N=gridData.inputData.Nx;
      offset=gridData.colloq_start_index_x;
    }
    else if (dim==1){
      delta=gridData.inputData.dy;
      N_local=gridData.local_colloq_y;
      N=gridData.inputData.Ny;
      offset=gridData.colloq_start_index_y;
    }
    else if (dim==2){
      delta=gridData.inputData.dz;
      N_local=gridData.local_colloq_z;
      N=gridData.inputData.Nz;
      offset=0;
    }
    else{
      std::cout<<"Warning: Dim should be 1, 2 or 3. Returning Default!"<<std::endl;
      w=1.0;
      return 0;
    }

    blitz::firstIndex ii;
    blitz::Array<double,1> temp;
    if (is_E){
      temp.resize(N_local);
      temp=ii+offset;
    }
    else{
      temp.resize(N_local-1);
      temp=ii+0.5+offset;
    }


    double d=d_pml/delta;
    w=1;
    w=(temp<=(d-1))*(pow((((d-1)-temp)/(d-1)),m_pml)*(k_max-1)+1)+(temp>(d-1))*w;
    w=(temp>=(N-d))*(pow((temp-N+d)/(d-1),m_pml)*(k_max-1)+1)+(temp<(N-d))*w;
    return 0;

}


/**
 ********************************************************************************************************************************************
 * \brief   Similar to sigma_fn() this function also calculates another quantity useful for PML application,
 *          specifically for complex PML. The function calculates the profile for all directions in
 *          rectangular grid.
 *          The input parameters are similar to sigma_fn() and serve same function.
 *
 *\param    dim specifes the dimension in which grading profile is to be calculated.
 *\param    is_E is set to true if profile is required for E field, otherwise it is set to false.
 *\param    w is the variable where the calculated grading profile will be stored.
 *
 ********************************************************************************************************************************************
 */
int maxwell::a_fn(int dim, bool is_E, blitz::Array<double,1> w){
    double delta=0;
    int N=0;
    int N_local=0;
    int offset=0;
    if (dim==0){
      delta=gridData.inputData.dx;
      N_local=gridData.local_colloq_x;
      N=gridData.inputData.Nx;
      offset=gridData.colloq_start_index_x;
    }
    else if (dim==1){
      delta=gridData.inputData.dy;
      N_local=gridData.local_colloq_y;
      N=gridData.inputData.Ny;
      offset=gridData.colloq_start_index_y;
    }
    else if (dim==2){
      delta=gridData.inputData.dz;
      N_local=gridData.local_colloq_z;
      N=gridData.inputData.Nz;
      offset=0;
    }
    else{
      std::cout<<"Warning: Dim should be 1, 2 or 3. Returning Default!"<<std::endl;
      w=0;
      return 0;
    }

    blitz::firstIndex ii;
    blitz::Array<double,1> temp;
    if (is_E){
      temp.resize(N_local);
      temp=ii+offset;
    }
    else{
      temp.resize(N_local-1);
      temp=ii+0.5+offset;
    }


    double d=d_pml/delta;
    w=0.0;
    w=(temp<=(d-1))*pow((temp/(d-1)),ma_pml)*a_max+(temp>(d-1))*w;
    w=(temp>=(N-d))*pow((N-1-temp)/(d-1),ma_pml)*a_max+(temp<(N-d))*w;
    return 0;
}


/**
 ********************************************************************************************************************************************
 * \brief   Similar to sigma_fn() this function also calculates another quantity useful for PML application,
 *          specifically for complex PML. The function calculates the profile for all directions in
 *          rectangular grid.
 *          The input parameters are similar to sigma_fn() and serve same function.
 *
 *\param    dim specifes the dimension in which grading profile is to be calculated.
 *\param    is_E is set to true if profile is required for E field, otherwise it is set to false.
 *\param    w is the variable where the calculated grading profile will be stored.
 *
 ********************************************************************************************************************************************
 */
int maxwell::b_fn(int dim, bool is_E, blitz::Array<double,1> w){
  blitz::Array<double,1> temp_sigma, temp_a, temp_k;
  int N=0;
  if (dim==0){
    N=gridData.local_colloq_x;
  }
  else if (dim==1){
    N=gridData.local_colloq_y;
  }
  else if (dim==2){
    N=gridData.local_colloq_z;
  }
  else{
    std::cout<<"Warning: Dim should be 1, 2 or 3. Returning Default!"<<std::endl;
  }
  if (is_E){
    temp_sigma.resize(N);
    sigma_fn(dim,true,temp_sigma);
    temp_a.resize(N);
    a_fn(dim,true,temp_a);
    temp_k.resize(N);
    k_fn(dim,true,temp_k);
  }
  else{
    temp_sigma.resize(N-1);
    sigma_fn(dim,false,temp_sigma);
    temp_a.resize(N-1);
    a_fn(dim,false,temp_a);
    temp_k.resize(N-1);
    k_fn(dim,false,temp_k);
  }

  w=exp(-((temp_sigma/temp_k/epsilon0)-(temp_a/epsilon0))*dt);

  return 0;
}


/**
 ********************************************************************************************************************************************
 * \brief   Similar to sigma_fn() this function also calculates another quantity useful for PML application,
 *          specifically for complex PML. The function calculates the profile for all directions in
 *          rectangular grid.
 *          The input parameters are similar to sigma_fn() and serve same function.
 *
 *\param    dim specifes the dimension in which grading profile is to be calculated.
 *\param    is_E is set to true if profile is required for E field, otherwise it is set to false.
 *\param    w is the variable where the calculated grading profile will be stored.
 *
 ********************************************************************************************************************************************
 */
int maxwell::c_fn(int dim, bool is_E, blitz::Array<double,1> w){
  blitz::Array<double,1> temp_sigma, temp_a, temp_k;
  int N=0;
  if (dim==0){
    N=gridData.local_colloq_x;
  }
  else if (dim==1){
    N=gridData.local_colloq_y;
  }
  else if (dim==2){
    N=gridData.local_colloq_z;
  }
  else{
    std::cout<<"Warning: Dim should be 1, 2 or 3. Returning Default!"<<std::endl;
  }
  if (is_E){
    temp_sigma.resize(N);
    sigma_fn(dim,true,temp_sigma);
    temp_a.resize(N);
    a_fn(dim,true,temp_a);
    temp_k.resize(N);
    k_fn(dim,true,temp_k);
  }
  else{
    temp_sigma.resize(N-1);
    sigma_fn(dim,false,temp_sigma);
    temp_a.resize(N-1);
    a_fn(dim,false,temp_a);
    temp_k.resize(N-1);
    k_fn(dim,false,temp_k);
  }

  b_fn(dim,is_E,w);
  if (!is_E){
    N=N-1;
  }
  int i;
  for (i=0;i<N;i++){
    if (temp_sigma(i)==0){
      w(i)=0.0;
    }
    else{
      w(i)=temp_sigma(i)/(temp_sigma(i)*temp_k(i)+pow(temp_k(i),2)*temp_a(i))*(w(i)-1);
    }
  }
  return 0;
}
