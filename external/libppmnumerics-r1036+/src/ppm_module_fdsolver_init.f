#ifdef __XLF
@PROCESS NOHOT
#endif
      !-------------------------------------------------------------------------
      !  Module       :                ppm_module_fdsolver_solve
      !-------------------------------------------------------------------------
      !
      !  Purpose      : This module includes the source code for field solver 
      !                        routines.
      !
      !  Remarks      :
      !
      !  References   :
      !
      !  Revisions    :
      !-------------------------------------------------------------------------
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
     
      !-------------------------------------------------------------------------
      !  Define types
      !-------------------------------------------------------------------------
#define __SINGLE_PRECISION  1
#define __DOUBLE_PRECISION  2
#define __COMPLEX           3
#define __DOUBLE_COMPLEX    4

#define __SFIELD            9
#define __VFIELD           10

#define __2D               22
#define __3D               33


   MODULE ppm_module_fdsolver_init

         !----------------------------------------------------------------------
         !  Define interface to ppm_fdsolver_init
         !----------------------------------------------------------------------
         INTERFACE ppm_fdsolver_init
            MODULE PROCEDURE ppm_fdsolver_init_2d_sca_s
            MODULE PROCEDURE ppm_fdsolver_init_2d_sca_d
            MODULE PROCEDURE ppm_fdsolver_init_2d_vec_s
            MODULE PROCEDURE ppm_fdsolver_init_2d_vec_d
            MODULE PROCEDURE ppm_fdsolver_init_3d_sca_s
            MODULE PROCEDURE ppm_fdsolver_init_3d_sca_d
            MODULE PROCEDURE ppm_fdsolver_init_3d_vec_s
            MODULE PROCEDURE ppm_fdsolver_init_3d_vec_d
         END INTERFACE
         !----------------------------------------------------------------------
         !  include the source 
         !----------------------------------------------------------------------
         CONTAINS

#define __MESH_DIM __2D
#define __DIM __SFIELD
#define __KIND __SINGLE_PRECISION
#include "ppm_fdsolver_init.f"
#undef  __KIND

#define __KIND __DOUBLE_PRECISION
#include "ppm_fdsolver_init.f"
#undef  __KIND
#undef __DIM


#define __DIM __VFIELD
#define __KIND __SINGLE_PRECISION
#include "ppm_fdsolver_init.f"
#undef  __KIND

#define __KIND __DOUBLE_PRECISION
#include "ppm_fdsolver_init.f"
#undef  __KIND
#undef __DIM
#undef __MESH_DIM




#define __MESH_DIM __3D
#define __DIM __SFIELD
#define __KIND __SINGLE_PRECISION
#include "ppm_fdsolver_init.f"
#undef  __KIND

#define __KIND __DOUBLE_PRECISION
#include "ppm_fdsolver_init.f"
#undef  __KIND
#undef __DIM


#define __DIM __VFIELD
#define __KIND __SINGLE_PRECISION
#include "ppm_fdsolver_init.f"
#undef  __KIND

#define __KIND __DOUBLE_PRECISION
#include "ppm_fdsolver_init.f"
#undef  __KIND
#undef __DIM
#undef __MESH_DIM


      END MODULE ppm_module_fdsolver_init
