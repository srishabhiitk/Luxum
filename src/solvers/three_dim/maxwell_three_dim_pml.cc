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


 /*! \file  maxwell_three_dim_pml.cc
  *
  * \brief  This file contains funtion to run three dimensional simulation with
  *         Perfectly Matched Layer.
  * \author Rishabh Sahu
  * \date   Mar 2018
  */



#include "../maxwell.h"


/**
 ********************************************************************************************************************************************
 * \brief   Function that calculates the timestepping solution when grid specified is three-dimensional
 *          with absorbing boundary condition.
 ********************************************************************************************************************************************
 */
void maxwell::solve3d_pml(){
    double Ca = 1.0, Da = 1.0;
    int t = 0; //timestep
    int Nx = gridData.local_colloq_x;
    int Ny = gridData.local_colloq_y;
    int Nz = gridData.local_colloq_z;
    double local_max, max = 0; //Variables to keep track of maximum value of E_total.


    //Declaring required variables to store streched PML grading.
    blitz::Array<double,1> k_ex(Nx),k_ey(Ny),k_ez(Nz);
    blitz::Array<double,1> k_hx(Nx-1),k_hy(Ny-1),k_hz(Nz-1);

    k_fn(0,true,k_ex);
    k_fn(1,true,k_ey);
    k_fn(2,true,k_ez);
    k_fn(0,false,k_hx);
    k_fn(1,false,k_hy);
    k_fn(2,false,k_hz);


    //Declaring required variables to store complex PML grading.
    blitz::Array<double,1> b_ex(Nx),b_ey(Ny),b_ez(Nz),c_ex(Nx),c_ey(Ny),c_ez(Nz);
    blitz::Array<double,1> b_hx(Nx-1),b_hy(Ny-1),b_hz(Nz-1),c_hx(Nx-1),c_hy(Ny-1),c_hz(Nz-1);

    b_fn(0,true,b_ex);
    b_fn(1,true,b_ey);
    b_fn(2,true,b_ez);
    c_fn(0,true,c_ex);
    c_fn(1,true,c_ey);
    c_fn(2,true,c_ez);
    b_fn(0,false,b_hx);
    b_fn(1,false,b_hy);
    b_fn(2,false,b_hz);
    c_fn(0,false,c_hx);
    c_fn(1,false,c_hy);
    c_fn(2,false,c_hz);



    vfield *E = new vfield(gridData,false);
    vfield *Etotal = new vfield(gridData,false);
    vfield *H = new vfield(gridData,true);

    vfield *Cb = new vfield(gridData,false);
    vfield *Db = new vfield(gridData,true);

    vfield *Ecurl1 = new vfield(gridData,true);
    vfield *Ecurl2 = new vfield(gridData,true);
    vfield *Hcurl1 = new vfield(gridData,false);
    vfield *Hcurl2 = new vfield(gridData,false);

    vfield *psi_E1 = new vfield(gridData,false);
    vfield *psi_E2 = new vfield(gridData,false);
    vfield *psi_H1 = new vfield(gridData,true);
    vfield *psi_H2 = new vfield(gridData,true);

    fileio Cb_reader_x("airplane_100_x.h5", true, gridData, Cb->Vx);
    fileio Cb_reader_y("airplane_100_y.h5", true, gridData, Cb->Vy);
    fileio Cb_reader_z("airplane_100_z.h5", true, gridData, Cb->Vz);

    fileio Ez_writer("Ez_data.h5", false, gridData, Etotal->Vz);


    Cb_reader_x.readHDF5("data");
    Cb_reader_y.readHDF5("data");
    Cb_reader_z.readHDF5("data");



    Cb->Vx->F=epsilon0*(1+Cb->Vx->F*1000000);
    Cb->Vy->F=epsilon0*(1+Cb->Vy->F*1000000);
    Cb->Vz->F=epsilon0*(1+Cb->Vz->F*1000000);

    Cb->Vx->F=dt/Cb->Vx->F;
    Cb->Vy->F=dt/Cb->Vy->F;
    Cb->Vz->F=dt/Cb->Vz->F;

    Db->Vx->F=dt/mew0;
    Db->Vy->F=dt/mew0;
    Db->Vz->F=dt/mew0;


    /****************  The following code can be used to put a solid block of some refractive index
                       in center of 100X100X100 grid of size 10X10X10.         *******************/
    // //Solid block at center
    // blitz::TinyVector<int,3> loBound;
    // loBound=45,45,45;
    // loBound = global_to_local(loBound);
    // blitz::TinyVector<int,3> upBound;
    // upBound=55,55,55;
    // upBound = global_to_local(upBound);
    // blitz::RectDomain<3> rd;
    // rd=blitz::RectDomain<3>(loBound,upBound);
    // Cb->Vx->F(rd) = dt/(epsilon0*6);
    // Cb->Vy->F(rd) = dt/(epsilon0*6);
    // Cb->Vz->F(rd) = dt/(epsilon0*6);




    blitz::TinyVector<int,3> source;
    source = int(Nx_source), int(Ny_source), int(Nz_source);


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

        E->curl_3d_adv(Ecurl1,Ecurl2);
        //Change the curl if plane-wave is required.
        if (usePW){
            plane_wave_execute(Ecurl1,Ecurl2,t);
        }

        #pragma omp parallel for num_threads(n_threads) schedule(dynamic)
        for (int i=0;i<Nx;i++){
          for (int j=0;j<(Ny-1);j++){
            for (int k=0;k<(Nz-1);k++){
              psi_H1->Vx->F(i,j,k)=b_hy(j)*psi_H1->Vx->F(i,j,k)+c_hy(j)*Ecurl1->Vx->F(i,j,k);
              psi_H2->Vx->F(i,j,k)=b_hz(k)*psi_H2->Vx->F(i,j,k)+c_hz(k)*Ecurl2->Vx->F(i,j,k);
              H->Vx->F(i,j,k)=Da*H->Vx->F(i,j,k)-Db->Vx->F(i,j,k)*(Ecurl1->Vx->F(i,j,k)/k_hy(j)-Ecurl2->Vx->F(i,j,k)/k_hz(k))-Db->Vx->F(i,j,k)*(psi_H1->Vx->F(i,j,k)-psi_H2->Vx->F(i,j,k));
            }
          }
        }

        #pragma omp parallel for num_threads(n_threads) schedule(dynamic)
        for (int i=0;i<(Nx-1);i++){
          for (int j=0;j<Ny;j++){
            for (int k=0;k<(Nz-1);k++){
              psi_H1->Vy->F(i,j,k)=b_hz(k)*psi_H1->Vy->F(i,j,k)+c_hz(k)*Ecurl1->Vy->F(i,j,k);
              psi_H2->Vy->F(i,j,k)=b_hx(i)*psi_H2->Vy->F(i,j,k)+c_hx(i)*Ecurl2->Vy->F(i,j,k);
              H->Vy->F(i,j,k)=Da*H->Vy->F(i,j,k)-Db->Vy->F(i,j,k)*(Ecurl1->Vy->F(i,j,k)/k_hz(k)-Ecurl2->Vy->F(i,j,k)/k_hx(i))-Db->Vy->F(i,j,k)*(psi_H1->Vy->F(i,j,k)-psi_H2->Vy->F(i,j,k));
            }
          }
        }

        #pragma omp parallel for num_threads(n_threads) schedule(dynamic)
        for (int i=0;i<(Nx-1);i++){
          for (int j=0;j<(Ny-1);j++){
            for (int k=0;k<Nz;k++){
              psi_H1->Vz->F(i,j,k)=b_hx(i)*psi_H1->Vz->F(i,j,k)+c_hx(i)*Ecurl1->Vz->F(i,j,k);
              psi_H2->Vz->F(i,j,k)=b_hy(j)*psi_H2->Vz->F(i,j,k)+c_hy(j)*Ecurl2->Vz->F(i,j,k);
              H->Vz->F(i,j,k)=Da*H->Vz->F(i,j,k)-Db->Vz->F(i,j,k)*(Ecurl1->Vz->F(i,j,k)/k_hx(i)-Ecurl2->Vz->F(i,j,k)/k_hy(j))-Db->Vz->F(i,j,k)*(psi_H1->Vz->F(i,j,k)-psi_H2->Vz->F(i,j,k));
            }
          }
        }




        H->curl_3d_adv(Hcurl1,Hcurl2);
        //Change the curl if plane-wave is required.
        if (usePW){
            plane_wave_execute(Hcurl1,Hcurl2,t);
        }

        #pragma omp parallel for num_threads(n_threads) schedule(dynamic)
        for (int i=0;i<(Nx-1);i++){
          for (int j=0;j<Ny;j++){
            for (int k=0;k<Nz;k++){
              psi_E1->Vx->F(i,j,k)=b_ey(j)*psi_E1->Vx->F(i,j,k)+c_ey(j)*Hcurl1->Vx->F(i,j,k);
              psi_E2->Vx->F(i,j,k)=b_ez(k)*psi_E2->Vx->F(i,j,k)+c_ez(k)*Hcurl2->Vx->F(i,j,k);
              E->Vx->F(i,j,k)=Ca*E->Vx->F(i,j,k)+Cb->Vx->F(i,j,k)*(Hcurl1->Vx->F(i,j,k)/k_ey(j)-Hcurl2->Vx->F(i,j,k)/k_ez(k))+Cb->Vx->F(i,j,k)*(psi_E1->Vx->F(i,j,k)-psi_E2->Vx->F(i,j,k));
            }
          }
        }

        #pragma omp parallel for num_threads(n_threads) schedule(dynamic)
        for (int i=0;i<Nx;i++){
          for (int j=0;j<(Ny-1);j++){
            for (int k=0;k<Nz;k++){
              psi_E1->Vy->F(i,j,k)=b_ez(k)*psi_E1->Vy->F(i,j,k)+c_ez(k)*Hcurl1->Vy->F(i,j,k);
              psi_E2->Vy->F(i,j,k)=b_ex(i)*psi_E2->Vy->F(i,j,k)+c_ex(i)*Hcurl2->Vy->F(i,j,k);
              E->Vy->F(i,j,k)=Ca*E->Vy->F(i,j,k)+Cb->Vy->F(i,j,k)*(Hcurl1->Vy->F(i,j,k)/k_ez(k)-Hcurl2->Vy->F(i,j,k)/k_ex(i))+Cb->Vy->F(i,j,k)*(psi_E1->Vy->F(i,j,k)-psi_E2->Vy->F(i,j,k));
            }
          }
        }

        #pragma omp parallel for num_threads(n_threads) schedule(dynamic)
        for (int i=0;i<Nx;i++){
          for (int j=0;j<Ny;j++){
            for (int k=0;k<(Nz-1);k++){
              psi_E1->Vz->F(i,j,k)=b_ex(i)*psi_E1->Vz->F(i,j,k)+c_ex(i)*Hcurl1->Vz->F(i,j,k);
              psi_E2->Vz->F(i,j,k)=b_ey(j)*psi_E2->Vz->F(i,j,k)+c_ey(j)*Hcurl2->Vz->F(i,j,k);
              E->Vz->F(i,j,k)=Ca*E->Vz->F(i,j,k)+Cb->Vz->F(i,j,k)*(Hcurl1->Vz->F(i,j,k)/k_ex(i)-Hcurl2->Vz->F(i,j,k)/k_ey(j))+Cb->Vz->F(i,j,k)*(psi_E1->Vz->F(i,j,k)-psi_E2->Vz->F(i,j,k));
            }
          }
        }




        E->Vz->mpiHandle->syncData();
        E->Vy->mpiHandle->syncData();
        E->Vx->mpiHandle->syncData();



        #pragma omp parallel for num_threads(n_threads) schedule(dynamic)
        for (int i=0;i<(Nx-1);i++){
          for (int j=0;j<(Ny-1);j++){
            for (int k=0;k<(Nz-1);k++){
              Etotal->Vz->F(i,j,k)=E->Vz->F(i,j,k)*E->Vz->F(i,j,k)+E->Vy->F(i,j,k)*E->Vy->F(i,j,k)+E->Vx->F(i,j,k)*E->Vx->F(i,j,k);
            }
          }
        }


        local_max = blitz::max(Etotal->Vz->F);

        //Writing the data after certain number of timesteps
        if ((t%3==0)){
            Ez_writer.writeHDF5(t);
        }

    }//end of main timestepping loop

    Ez_writer.closeFileio();

}
