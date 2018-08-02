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


 /*! \file  maxwell_two_dim_pml.cc
  *
  * \brief  This file contains funtion to run two dimensional simulation with
  *         Perfectly Matched Layer.
  * \author Rishabh Sahu
  * \date   Mar 2018
  */



#include "../maxwell.h"

/**
 ********************************************************************************************************************************************
 * \brief   Function that calculates the timestepping solution when grid specified is two-dimensional
 *          with absorbing boundary condition.
 *
 * \author  Rishabh Sahu
 * \date    4 March 2018
 ********************************************************************************************************************************************
 */
void maxwell::solve_planar_pml(){
    double Ca = 1.0, Da = 1.0;
    int t = 0;//timestep
    int Nx = gridData.local_colloq_x;
    int Nz = gridData.local_colloq_z;
    double local_max, max = 0;//Variables to keep track of maximum value of E_total.


    //Declaring required variables to store streched PML grading.
    blitz::Array<double,1> k_ex(Nx),k_ez(Nz);
    blitz::Array<double,1> k_hx(Nx-1),k_hz(Nz-1);

    k_fn(0,true,k_ex);
    k_fn(2,true,k_ez);
    k_fn(0,false,k_hx);
    k_fn(2,false,k_hz);


    //Declaring required variables to store complex PML grading.
    blitz::Array<double,1> b_ex(Nx),b_ez(Nz),c_ex(Nx),c_ez(Nz);
    blitz::Array<double,1> b_hx(Nx-1),b_hz(Nz-1),c_hx(Nx-1),c_hz(Nz-1);

    b_fn(0,true,b_ex);
    b_fn(2,true,b_ez);
    c_fn(0,true,c_ex);
    c_fn(2,true,c_ez);
    b_fn(0,false,b_hx);
    b_fn(2,false,b_hz);
    c_fn(0,false,c_hx);
    c_fn(2,false,c_hz);


    vfield *E = new vfield(gridData,false);
    vfield *H = new vfield(gridData,true);

    vfield *Cb = new vfield(gridData,false);
    vfield *Db = new vfield(gridData,true);

    vfield *Ecurl1 = new vfield(gridData,true);
    vfield *Ecurl2 = new vfield(gridData,true);
    vfield *Hcurl1 = new vfield(gridData,false);
    vfield *Hcurl2 = new vfield(gridData,false);
    //Resizing the nonessential variables to 0 size to prevent accidental usage.
    Ecurl1->Vx->F.resize(0,0,0);
    Ecurl2->Vz->F.resize(0,0,0);

    vfield *psi_E1 = new vfield(gridData,false);
    vfield *psi_E2 = new vfield(gridData,false);
    vfield *psi_H1 = new vfield(gridData,true);
    vfield *psi_H2 = new vfield(gridData,true);
    //Resizing the nonessential variables to 0 size to prevent accidental usage.
    psi_H1->Vx->F.resize(0,0,0);
    psi_H2->Vz->F.resize(0,0,0);


    /****** Following code reads a Cb file to initiale random structure with a refrcative index. *****/
    // fileio Cb_reader("wg2.h5", true, gridData, Cb->Vy);
    // Cb_reader.readHDF5("data");
    // Cb->Vy->F=(Cb->Vy->F*11 + 1)*epsilon0;
    // Cb->Vy->F=dt/Cb->Vy->F;


    Cb->Vy->F=dt/epsilon0;
    Db->Vx->F=dt/mew0;
    Db->Vz->F=dt/mew0;

    blitz::TinyVector<int,3> source;
    source = int(Nx_source), 0, int(Nz_source);

    fileio Ey_writer("Ey_data_2d.h5", false, gridData, E->Vy);

    //Main timestepping loop
    for (t=0; t<num_timesteps; t++){

        //Putting a point source if required
        if (usePS){
            if (check_global_limits(source, false, false)){
                E->Vz->F(global_to_local(source))=amp_source*sin(omega_source*t*dt);
            }
        }

        //Taking the maximum value of E_total from all MPI nodes.
        MPI_Reduce(&local_max, &max, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

        // Print the maximum value of E_total using one processor.
        // This keeps the solution in check by monitering for any exploding solutions.
        if (mpi.rank==0){
            std::cout<<"t: "<<t<<", max: "<<max<<std::endl;
        }

        E->curl_planar_adv(Ecurl1, Ecurl2);
        //Change the curl if plane-wave is required.
        if (usePW){
            plane_wave_execute(Ecurl1, Ecurl2, t);
        }
        #pragma omp parallel for num_threads(n_threads) schedule(dynamic)
        for (int i=0;i<Nx;i++){
            for (int j=0;j<(Nz-1);j++){
                psi_H2->Vx->F(i,0,j)=b_hz(j)*psi_H2->Vx->F(i,0,j)+c_hz(j)*Ecurl2->Vx->F(i,0,j);
                H->Vx->F(i,0,j)=Da*H->Vx->F(i,0,j)-Db->Vx->F(i,0,j)*(-Ecurl2->Vx->F(i,0,j)/k_hz(j)-psi_H2->Vx->F(i,0,j));
            }
        }
        #pragma omp parallel for num_threads(n_threads) schedule(dynamic)
        for (int i=0;i<(Nx-1);i++){
            for (int j=0;j<Nz;j++){
                psi_H1->Vz->F(i,0,j)=b_hx(i)*psi_H1->Vz->F(i,0,j)+c_hx(i)*Ecurl1->Vz->F(i,0,j);
                H->Vz->F(i,0,j)=Da*H->Vz->F(i,0,j)-Db->Vz->F(i,0,j)*(Ecurl1->Vz->F(i,0,j)/k_hx(i)+psi_H1->Vz->F(i,0,j));
            }
        }

        H->curl_planar_adv(Hcurl1, Hcurl2);
        //Change the curl if plane-wave is required.
        if (usePW){
            plane_wave_execute(Hcurl1, Hcurl2, t);
        }
        #pragma omp parallel for num_threads(n_threads) schedule(dynamic)
        for (int i=0;i<(Nx-1);i++){
            for (int j=0;j<(Nz-1);j++){
                psi_E1->Vy->F(i,0,j)=b_ez(j)*psi_E1->Vy->F(i,0,j)+c_ez(j)*Hcurl1->Vy->F(i,0,j);
                psi_E2->Vy->F(i,0,j)=b_ex(i)*psi_E2->Vy->F(i,0,j)+c_ex(i)*Hcurl2->Vy->F(i,0,j);
                E->Vy->F(i,0,j)=Ca*E->Vy->F(i,0,j)+Cb->Vy->F(i,0,j)*(Hcurl1->Vy->F(i,0,j)/k_ez(j)-Hcurl2->Vy->F(i,0,j)/k_ex(i)+psi_E1->Vy->F(i,0,j)-psi_E2->Vy->F(i,0,j));
            }
        }

        //In this case, syncing of only Ey is required!
        E->Vy->mpiHandle->syncData();

        local_max = blitz::max(E->Vy->F*E->Vy->F);

        //Writing the data after certain number of timesteps
        if (t%10==0){
          Ey_writer.writeHDF5(t);
        }

    }//end of main timestepping loop

    Ey_writer.closeFileio();
    // Cb_reader.closeFileio();
}
