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


 /*! \file  maxwell_three_dim.cc
  *
  * \brief  This file contains funtion to run three dimensional simulation without
  *         without Perfectly Matched Layer.
  * \author Rishabh Sahu
  * \date   Mar 2018
  */




#include "../maxwell.h"

/**
 ********************************************************************************************************************************************
 * \brief   Function that calculates the timestepping solution when grid specified is three-dimensional
 *          without absorbing boundary condition.
 ********************************************************************************************************************************************
 */
void maxwell::solve3d(){
    double Ca = 1.0, Da = 1.0;
    double Cb = dt/epsilon0;
    double Db = dt/mew0;

    int t = 0; //timestep
    int Nx = gridData.local_colloq_x;
    int Ny = gridData.local_colloq_y;
    int Nz = gridData.local_colloq_z;
    double local_max, max = 0; //Variables to keep track of maximum value of E_total.

    struct timeval begin,end;


    vfield *E = new vfield(gridData,false);
    vfield *H = new vfield(gridData,true);

    vfield *Ecurl = new vfield(gridData,true);
    vfield *Hcurl = new vfield(gridData,false);

    vfield *Etotal = new vfield(gridData,false);

    blitz::TinyVector<int,3> source;
    source = int(inputData.Nx_source), int(inputData.Ny_source), int(inputData.Nz_source);

    fileio Ez_writer("Ez_data2.h5", false, gridData, Etotal->Vz);

    gettimeofday(&begin,NULL); //Get start time before main loop begins

    // Main timestepping loop
    for (t=0; t<num_timesteps; t++){

        // Putting the point source
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


        E->curl_3d(Ecurl);
        //Change the curl if plane-wave is required.
        if (usePW){
            plane_wave_execute(Ecurl,t);
        }

        #pragma omp parallel for num_threads(n_threads) schedule(dynamic)
        for (int i=0;i<Nx;i++){
          for (int j=0;j<(Ny-1);j++){
            for (int k=0;k<(Nz-1);k++){
              H->Vx->F(i,j,k)=Da*H->Vx->F(i,j,k)-Db*Ecurl->Vx->F(i,j,k);
            }
          }
        }

        #pragma omp parallel for num_threads(n_threads) schedule(dynamic)
        for (int i=0;i<(Nx-1);i++){
          for (int j=0;j<Ny;j++){
            for (int k=0;k<(Nz-1);k++){
              H->Vy->F(i,j,k)=Da*H->Vy->F(i,j,k)-Db*Ecurl->Vy->F(i,j,k);
            }
          }
        }

        #pragma omp parallel for num_threads(n_threads) schedule(dynamic)
        for (int i=0;i<(Nx-1);i++){
          for (int j=0;j<(Ny-1);j++){
            for (int k=0;k<Nz;k++){
              H->Vz->F(i,j,k)=Da*H->Vz->F(i,j,k)-Db*Ecurl->Vz->F(i,j,k);
            }
          }
        }

        H->curl_3d(Hcurl);
        //Change the curl if plane-wave is required.
        if (usePW){
            plane_wave_execute(Hcurl,t);
        }

        #pragma omp parallel for num_threads(n_threads) schedule(dynamic)
        for (int i=0;i<(Nx-1);i++){
          for (int j=0;j<Ny;j++){
            for (int k=0;k<Nz;k++){
              E->Vx->F(i,j,k)=Ca*E->Vx->F(i,j,k)+Cb*Hcurl->Vx->F(i,j,k);
            }
          }
        }

        #pragma omp parallel for num_threads(n_threads) schedule(dynamic)
        for (int i=0;i<Nx;i++){
          for (int j=0;j<(Ny-1);j++){
            for (int k=0;k<Nz;k++){
              E->Vy->F(i,j,k)=Ca*E->Vy->F(i,j,k)+Cb*Hcurl->Vy->F(i,j,k);
            }
          }
        }

        #pragma omp parallel for num_threads(n_threads) schedule(dynamic)
        for (int i=0;i<Nx;i++){
          for (int j=0;j<Ny;j++){
            for (int k=0;k<(Nz-1);k++){
              E->Vz->F(i,j,k)=Ca*E->Vz->F(i,j,k)+Cb*Hcurl->Vz->F(i,j,k);
            }
          }
        }




        E->Vz->mpiHandle->syncData();
        E->Vy->mpiHandle->syncData();
        E->Vx->mpiHandle->syncData();



        //Has to be placed after syncing because sometimes in plane-wave simulation
        //curl gets updated by plane-wave where it is not updated (is at boundary)
        //it does not affect anything because things get reset after syncing but
        //curl specifically at that point keeps on getting edited! (NOT RESOLVED YET!)
        #pragma omp parallel for num_threads(n_threads) schedule(dynamic)
        for (int i=0;i<(Nx-1);i++){
          for (int j=0;j<(Ny-1);j++){
            for (int k=0;k<(Nz-1);k++){
              Etotal->Vz->F(i,j,k)=E->Vz->F(i,j,k)*E->Vz->F(i,j,k)+E->Vy->F(i,j,k)*E->Vy->F(i,j,k)+E->Vx->F(i,j,k)*E->Vx->F(i,j,k);
            }
          }
        }

        local_max = blitz::max(Etotal->Vz->F);


        if (t%10==0){
            Ez_writer.writeHDF5(t);
        }

    } // End of main timestepping loop

    gettimeofday(&end,NULL);
    std::cout<<"Time taken: "<<((end.tv_sec-begin.tv_sec)*1000000u + end.tv_usec - begin.tv_usec)/1.0e6<<" sec"<<std::endl;
    Ez_writer.closeFileio();
}
