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


 /*! \file  fileio.cc
  *
  * \brief  This file contains contructor and definitions of member functions
  *         of fileio class. 
  * \author Rishabh Sahu
  * \date   Mar 2018
  */



#include "fileio.h"

/**
 ********************************************************************************************************************************************
 * \brief   Constructor of the fileio class
 *
 *          The constructor initializes the datamembers according to the scalar field provided.
 *          The subarray that is contained in the specific MPI node is calculated according to the
 *          processor rank. The subarray is size is different for scalar field staggered differently
 *          and whether the object will be used to read or write the data.
 *
 * \param   fileName is the name of file whose data is to be read or in which data is to be written.
 * \param   read is a boolean which is set to true if object will be used to read HDF5 file.
 * \param   mesh is a constant reference to grid information.
 * \param   iField is a reference to the scalar field object whose data is to be read or written.
 ********************************************************************************************************************************************
 */
fileio::fileio(const char *fileName, bool read, const grid &mesh, field *iField): mesh(mesh), outField(iField) {
    // Create a property list for collectively opening a file by all processors
    plist_id = H5Pcreate(H5P_FILE_ACCESS);
    H5Pset_fapl_mpio(plist_id, MPI_COMM_WORLD, MPI_INFO_NULL);

    // First create a file handle with the path to the file
    if (read){
        fileHandle = H5Fopen(fileName, H5F_ACC_RDONLY, plist_id);
    }
    else{
        fileHandle = H5Fcreate(fileName, H5F_ACC_TRUNC, H5P_DEFAULT, plist_id);
    }


    // Close the property list for later reuse
    H5Pclose(plist_id);




    //Preparing source Dataspace that will be written to bigger Dataspace

    locSize[0] = outField->local_Nx;
    locSize[1] = outField->local_Ny;
    locSize[2] = outField->local_Nz;

    memSpace = H5Screate_simple(3, locSize, NULL);

    locOffset[0] = 0;
    locOffset[1] = 0;
    locOffset[2] = 0;

    if (!read){
        //in x-driection
        if (outField->xStag){
            locSize[0] = locSize[0]-1; //if field is stagged
            if (mesh.rankData.xRank==(mesh.rankData.npX-1)){
                locSize[0] = locSize[0]+1;
            }
        }
        else{
            locSize[0] = locSize[0]-2; //if field is colloquated
            locOffset[0] = 1;
            if (mesh.rankData.xRank==(mesh.rankData.npX-1)){
                locSize[0] = locSize[0]+1;
            }
            if (mesh.rankData.xRank==0){
                locSize[0] = locSize[0]+1;
                locOffset[0] = 0;
            }
        }

        //in y direction
        if (outField->yStag){
            locSize[1] = locSize[1]-1;
            if (mesh.rankData.yRank==(mesh.rankData.npY-1)){
                locSize[1] = locSize[1]+1;
            }
        }
        else{
            locSize[1] = locSize[1]-2;
            locOffset[1] = 1;
            if (mesh.rankData.yRank==(mesh.rankData.npY-1)){
                locSize[1] = locSize[1]+1;
            }
            if (mesh.rankData.yRank==0){
                locSize[1] = locSize[1]+1;
                locOffset[1] = 0;
            }
        }
    }

    status = H5Sselect_hyperslab(memSpace, H5S_SELECT_SET, locOffset, NULL, locSize, NULL);
    if (status) {
        if (mesh.rankData.rank == 0) {
            std::cout << "Error in creating hyperslab while writing data. Aborting" << std::endl;
        }
        exit(0);
    }



    //Now preparing the bigger dataspace in which the data from source will be written to
    gloSize[0] = mesh.inputData.Nx;
    gloSize[1] = mesh.inputData.Ny;
    gloSize[2] = mesh.inputData.Nz;

    if (outField->xStag){
        gloSize[0] = gloSize[0]-1;
    }
    if (outField->yStag){
        gloSize[1] = gloSize[1]-1;
    }
    if (outField->zStag){
        gloSize[2] = gloSize[2]-1;
    }

    fileSpace = H5Screate_simple(3, gloSize, NULL);

    gloOffset[0] = mesh.colloq_start_index_x + locOffset[0];
    gloOffset[1] = mesh.colloq_start_index_y + locOffset[1];
    gloOffset[2] = 0;

    status = H5Sselect_hyperslab(fileSpace, H5S_SELECT_SET, gloOffset, NULL, locSize, NULL);
    if (status) {
        if (mesh.rankData.rank == 0) {
            std::cout << "Error in creating hyperslab while writing data. Aborting" << std::endl;
        }
        exit(0);
    }
}



/**
 ********************************************************************************************************************************************
 * \brief   Function to write scalar field data to HDF5 file
 *
 *          The function write the data to scalar field. The necessary information is taken
 *          from datamember of the object calling the function.
 *
 * \param   t is dataset number which can be changed for subsequent writes at different time steps.
 ********************************************************************************************************************************************
 */
void fileio::writeHDF5(double t) {

    temp<<"t="<<t;
    datasetName = temp.str();

    // Create the dataset *for the file*, linking it to the file handle.
    // Correspondingly, it will use the *core* dataspace, as only the core has to be written excluding the pads
    dataSet = H5Dcreate2(fileHandle, datasetName.c_str(), H5T_NATIVE_DOUBLE, fileSpace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    // Create a property list to use collective data write
    plist_id = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);

    // Write the dataset. Most important thing to note is that the 3rd and 4th arguments represent the *source* and *destination* dataspaces.
    // The source here is the memSpace pointing to the memory buffer. Note that its view has been adjusted using hyperslab.
    // The destination is the fileSpace. Though the fileSpace is smaller than the memSpace,
    // only the appropriate hyperslab within the memSpace is transferred to the destination.
    status = H5Dwrite(dataSet, H5T_NATIVE_DOUBLE, memSpace, fileSpace, plist_id, outField->F.dataFirst());
    if (status) {
        if (mesh.rankData.rank == 0) {
            std::cout << "Error in writing output to HDF file. Aborting" << std::endl;
        }
        exit(0);
    }

    temp.str("");
    H5Pclose(plist_id);
    H5Dclose(dataSet);
}

/**
 ********************************************************************************************************************************************
 * \brief   Function to read data from HDF5 file
 *
 *          The function reads the dataset from HDF5 file and writes it to scalar field data.
 *          The necessary information is taken from datamember of the object calling the function.
 *
 * \param   datasetName is name of dataset which is to be read.
 ********************************************************************************************************************************************
 */
void fileio::readHDF5(std::string datasetName) {
    // Opening the dataset of the specified name
    dataSet = H5Dopen(fileHandle, datasetName.c_str(), H5P_DEFAULT);

    // Create a property list to use collective data write
    plist_id = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);

    // Write the dataset. Most important thing to note is that the 3rd and 4th arguments represent the *source* and *destination* dataspaces.
    // The source here is the memSpace pointing to the memory buffer. Note that its view has been adjusted using hyperslab.
    // The destination is the fileSpace. Though the fileSpace is smaller than the memSpace,
    // only the appropriate hyperslab within the memSpace is transferred to the destination.
    status = H5Dread(dataSet, H5T_NATIVE_DOUBLE, memSpace, fileSpace, plist_id, (void*)outField->F.data());
    if (status) {
        if (mesh.rankData.rank == 0) {
            std::cout << "Error in reading input from HDF file. Aborting" << std::endl;
        }
        exit(0);
    }

    temp.str("");
    H5Pclose(plist_id);
    H5Dclose(dataSet);
}

/**
 ********************************************************************************************************************************************
 * \brief   Function to close fileio object.
 *
 *          Call this before ending the program. This ensures a smoothly written HDF5 file.
 *          Not calling this may cause file to get correupted.
 ********************************************************************************************************************************************
 */
void fileio::closeFileio(){
    // CLOSE/RELEASE RESOURCES
    H5Sclose(memSpace);
    H5Sclose(fileSpace);
    H5Fclose(fileHandle);
}
