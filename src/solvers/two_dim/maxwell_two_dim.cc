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


 /*! \file  maxwell_two_dim.cc
  *
  * \brief  This file contains funtion to run two dimensional simulation without
  *         Perfectly Matched Layer.
  * \author Rishabh Sahu
  * \date   Mar 2018
  */



#include "../maxwell.h"

/**
 ********************************************************************************************************************************************
 * \brief   Function that calculates the timestepping solution when grid specified is two-dimensional
 *          without absorbing boundary condition.
 ********************************************************************************************************************************************
 */
void maxwell::solve_planar(){
    double Ca = 1.0, Da = 1.0;
    double Db = dt/mew0;
    int t = 0; //timestep
    int Nx = gridData.local_colloq_x;
    int Nz = gridData.local_colloq_z;
    double local_max, max = 0;


    struct timeval begin,end;


    vfield *E = new vfield(gridData,false);
    vfield *Cb = new vfield(gridData,false);
    vfield *H = new vfield(gridData,true);

    vfield *Ecurl = new vfield(gridData,true);
    vfield *Hcurl = new vfield(gridData,false);

    blitz::TinyVector<int,3> source;
    source = int(Nx_source), 0, int(Nz_source);

    fileio Ey_writer("Ey_data_2d.h5", false, gridData, E->Vy);

    // fileio Cb_reader("lens_distort.h5", true, gridData, Cb->Vy);


    // Cb_reader.readHDF5("data");

    // Cb->Vy->F=(Cb->Vy->F*1.5+1)*epsilon0;
    Cb->Vy->F=dt/epsilon0;

    gettimeofday(&begin,NULL); //Get start time before main loop begins
    // Main timestepping loop
    for (t=0; t<num_timesteps; t++){

        // Putting the point source
        if (usePS){
            if (check_global_limits(source, false, false)){
                E->Vz->F(global_to_local(source))=amp_source*sin(omega_source*t*dt);
            }
        }


        E->curl_planar(Ecurl);
        //Change the curl if plane-wave is required.
        if (usePW){
            plane_wave_execute(Ecurl,t);
        }
        #pragma omp parallel for num_threads(n_threads) schedule(dynamic)
        for (int i=0;i<Nx;i++){
          for (int j=0;j<(Nz-1);j++){
             H->Vx->F(i,0,j)=Da*H->Vx->F(i,0,j)-Db*Ecurl->Vx->F(i,0,j);
          }
        }
        #pragma omp parallel for num_threads(n_threads) schedule(dynamic)
        for (int i=0;i<(Nx-1);i++){
          for (int j=0;j<Nz;j++){
             H->Vz->F(i,0,j)=Da*H->Vz->F(i,0,j)-Db*Ecurl->Vz->F(i,0,j);
          }
        }

        H->curl_planar(Hcurl);
        //Change the curl if plane-wave is required.
        if (usePW){
            plane_wave_execute(Hcurl,t);
        }
        #pragma omp parallel for num_threads(n_threads) schedule(dynamic)
        for (int i=0;i<(Nx-1);i++){
          for (int j=0;j<(Nz-1);j++){
             E->Vy->F(i,0,j)=Ca*E->Vy->F(i,0,j)+Cb->Vy->F(i,0,j)*Hcurl->Vy->F(i,0,j);
          }
        }

        //In this case, syncing of only Ey is required!
        E->Vy->mpiHandle->syncData();

        local_max = blitz::max(E->Vy->F*E->Vy->F);

        //Taking the maximum value of E_total from all MPI nodes.
        MPI_Reduce(&local_max, &max, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

        // Print the maximum value of E_total using one processor.
        // This keeps the solution in check by monitering for any exploding solutions.
        if (mpi.rank==0){
            std::cout<<"t: "<<t<<", max: "<<max<<std::endl;
        }

        // Writing data after certain time intervals
        if (t%20==0 && t>=100){
            Ey_writer.writeHDF5(t);
        }

    } // End of main timestepping loop



    gettimeofday(&end,NULL);
    if (mpi.rank==0){
        std::cout<<"Time taken: "<<((end.tv_sec-begin.tv_sec)*1000000u + end.tv_usec - begin.tv_usec)/1.0e6<<" sec"<<std::endl;
    }

    Ey_writer.closeFileio();
    // Cb_reader.closeFileio();
}
