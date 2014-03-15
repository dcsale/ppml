      !-------------------------------------------------------------------------
      !  Module       :            ppm_module_hamjac_ext
      !-------------------------------------------------------------------------
      !
      !  Purpose      : procedure module for ppm_hamjac_ext
      !
      !  Remarks      : 
      !
      !  References   : 
      !
      !  Revisions    :
      !-------------------------------------------------------------------------
      !  $Log: ppm_module_hamjac_ext.f,v $
      !  Revision 1.1.1.1  2007/07/13 10:18:58  ivos
      !  CBL version of the PPM library
      !
      !  Revision 1.1  2005/07/25 00:34:07  ivos
      !  Initial check-in.
      !
      !-------------------------------------------------------------------------
      ! Copyright (c) 2010 CSE Lab (ETH Zurich), MOSAIC Group (ETH Zurich), 
      !                    Center for Fluid Dynamics (DTU)
      !
      !
      ! This file is part of the Parallel Particle Mesh Library (PPM).
      !
      ! PPM is free software: you can redistribute it and/or modify
      ! it under the terms of the GNU Lesser General Public License 
      ! as published by the Free Software Foundation, either 
      ! version 3 of the License, or (at your option) any later 
      ! version.
      !
      ! PPM is distributed in the hope that it will be useful,
      ! but WITHOUT ANY WARRANTY; without even the implied warranty of
      ! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
      ! GNU General Public License for more details.
      !
      ! You should have received a copy of the GNU General Public License
      ! and the GNU Lesser General Public License along with PPM. If not,
      ! see <http://www.gnu.org/licenses/>.
      !
      ! Parallel Particle Mesh Library (PPM)
      ! ETH Zurich
      ! CH-8092 Zurich, Switzerland

      !-------------------------------------------------------------------------

#define __SINGLE_PRECISION 1
#define __DOUBLE_PRECISION 2
#define __2D               3
#define __3D               4
#define __VEC              5
#define __SCA              6

      MODULE ppm_module_hamjac_ext

        !-----------------------------------------------------
        !  Interface
        !-----------------------------------------------------
        INTERFACE ppm_hamjac_ext_step

           MODULE PROCEDURE ppm_hamjac_ext_step_3ds
           MODULE PROCEDURE ppm_hamjac_ext_step_3dd
           MODULE PROCEDURE ppm_hamjac_ext_step_3dsv
           MODULE PROCEDURE ppm_hamjac_ext_step_3ddv

        END INTERFACE

        INTERFACE ppm_hamjac_ext

           MODULE PROCEDURE ppm_hamjac_ext_3ds
           MODULE PROCEDURE ppm_hamjac_ext_3dd
           MODULE PROCEDURE ppm_hamjac_ext_3dsv
           MODULE PROCEDURE ppm_hamjac_ext_3ddv

        END INTERFACE


      CONTAINS
#define __DIME  __3D
#define __MODE  __SCA
#define __KIND  __SINGLE_PRECISION
        ! 3D SCA SINGLE
#include "ppm_hamjac_ext_step_3d.f"
#undef __KIND
#define __KIND  __DOUBLE_PRECISION
        ! 3D SCA SINGLE
#include "ppm_hamjac_ext_step_3d.f"
#undef __KIND
#undef __MODE
#undef __DIME

#define __DIME  __3D
#define __MODE  __SCA
#define __KIND  __SINGLE_PRECISION
        ! 3D SCA SINGLE
#include "ppm_hamjac_ext_3d.f"
#undef __KIND
#define __KIND  __DOUBLE_PRECISION
        ! 3D SCA SINGLE
#include "ppm_hamjac_ext_3d.f"
#undef __KIND
#undef __MODE
#undef __DIME


#define __DIME  __3D
#define __MODE  __VEC
#define __KIND  __SINGLE_PRECISION
        ! 3D VEC SINGLE
#include "ppm_hamjac_ext_step_3d.f"
#undef __KIND
#define __KIND  __DOUBLE_PRECISION
        ! 3D VEC SINGLE
#include "ppm_hamjac_ext_step_3d.f"
#undef __KIND
#undef __MODE
#undef __DIME

#define __DIME  __3D
#define __MODE  __VEC
#define __KIND  __SINGLE_PRECISION
        ! 3D VEC SINGLE
#include "ppm_hamjac_ext_3d.f"
#undef __KIND
#define __KIND  __DOUBLE_PRECISION
        ! 3D VEC SINGLE
#include "ppm_hamjac_ext_3d.f"
#undef __KIND
#undef __MODE
#undef __DIME




      END MODULE ppm_module_hamjac_ext
        

        

