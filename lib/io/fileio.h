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


 /*! \file  fileio.h
  *
  * \brief  This file contains header of fileio class.
  * \author Rishabh Sahu
  * \date   Mar 2018
  */



#ifndef FILEIO_H
#define FILEIO_H

#include <blitz/array.h>
#include <iostream>
#include <string>
#include <sstream>
#include "field.h"
#include "grid.h"
#include "hdf5.h"
//using namespace H5;

class fileio {
    public:
        fileio(const char *fileName, bool read, const grid &mesh, field *iField);

        void writeHDF5(double t);
        void readHDF5(std::string datasetName);
        void closeFileio();

    private:
        /** A const reference to the grid information. */
        const grid &mesh;

        /** A reference to the scalar field object whose data is to be read or written. */
        field *outField;

        /** A reference to the scalar field object whose data is to be read or written. */
        std::stringstream temp;
        std::string datasetName;

        hsize_t gloSize[3];
        hsize_t locSize[3];

        hsize_t gloOffset[3];
        hsize_t locOffset[3];


        hid_t plist_id;

        hid_t fileHandle;

        hid_t dataSet;

        hid_t memSpace;
        hid_t fileSpace;

        herr_t status;



};

/**
 ********************************************************************************************************************************************
 *  \class fileio fileio.h "lib/fileio.h"
 *          \brief  Class to read and write data from and to HDF5 files.
 *
 *                  The objects of this class can be made using pointer to a scalar field object.
 *                  Once the rwHDF5 object is created, the member functions calls make it easy to
 *                  read and write data.
 *
 *                  writeHDF5() writes the current data stored in scalar field in the dataset name passed.
 *                  readHDF5() reads the dataset name passed and writes in scalar field object.
 ********************************************************************************************************************************************
 */

#endif
