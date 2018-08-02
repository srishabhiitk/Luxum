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


 /*! \file  vfield.h
  *
  * \brief  This file contains header of vfield class.
  * \author Rishabh Sahu
  * \date   Mar 2018
  */



#ifndef VFIELD_H
#define VFIELD_H

#include <blitz/array.h>
#include "field.h"
#include "grid.h"

class vfield {

    public:

        /** Pointer to scalar field object which contains X component. */
        field *Vx;

        /** Pointer to scalar field object which contains Y component. */
        field *Vy;

        /** Pointer to scalar field object which contains Z component. */
        field *Vz;

        /** A constant reference to the grid information. */
        const grid &gridData;

        /** Boolean representing whether the vector field is face centered. */
        bool isFaceCentered;

        /** Boolean representing whether the grid is planar or two dimensional. */
        bool isPlanar;



        vfield(const grid &gridData, bool isFaceCentered_);

        void curl_3d(vfield *curl);
        void curl_3d_adv(vfield *curl1,vfield *curl2);
        void curl_planar(vfield *curl);
        void curl_planar_adv(vfield *curl1, vfield *curl2);


};

/**
 ********************************************************************************************************************************************
 *  \class vfield vfield.h "lib/vfield.h"
 *  \brief  A class which abstracts concept of a vector field for Yee's algorithm.
 *
 *          Yee's algorithm requires E and H vector fields at specific lattice points - face center and
 *          edge center. This class abtracts these fields out. The type of field can be specified using
 *          input parameter whether it is face centered or not.
 *
 *          The class also defines actions that can be performed on these vector fields. For now, the
 *          actions are only limited curl functions, but more can be easily added.
 ********************************************************************************************************************************************
 */



#endif
