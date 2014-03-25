      !-------------------------------------------------------------------------
      !  Module       :             ppm_module_gmm_march_bkwd
      !-------------------------------------------------------------------------
      !
      !  Purpose      : This module includes the source code for the 
      !                 backward marching routine of the marching method.
      !
      !  Remarks      :
      !
      !  References   :
      !
      !  Revisions    :
      !-------------------------------------------------------------------------
      !  $Log: ppm_module_gmm_march_bkwd.f,v $
      !  Revision 1.1.1.1  2007/07/13 10:18:58  ivos
      !  CBL version of the PPM library
      !
      !  Revision 1.1  2005/04/27 01:08:32  ivos
      !  Initial commit, but tested.
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
#define __SINGLE_PRECISION         1
#define __DOUBLE_PRECISION         2
#define __2D                       3
#define __3D                       4

      MODULE ppm_module_gmm_march_bkwd

         !----------------------------------------------------------------------
         !  Define interfaces to ppm_gmm_march_bkwd
         !----------------------------------------------------------------------
         INTERFACE ppm_gmm_march_bkwd
            MODULE PROCEDURE ppm_gmm_march_bkwd_2ds
            MODULE PROCEDURE ppm_gmm_march_bkwd_2dd
            MODULE PROCEDURE ppm_gmm_march_bkwd_3ds
            MODULE PROCEDURE ppm_gmm_march_bkwd_3dd
         END INTERFACE

         !----------------------------------------------------------------------
         !  include the source 
         !----------------------------------------------------------------------
         CONTAINS

#define __DIM __2D
#define __KIND __SINGLE_PRECISION
#include "ppm_gmm_march_bkwd.f"
#undef __KIND
#define __KIND __DOUBLE_PRECISION
#include "ppm_gmm_march_bkwd.f"
#undef __KIND
#undef __DIM

#define __DIM __3D
#define __KIND __SINGLE_PRECISION
#include "ppm_gmm_march_bkwd.f"
#undef __KIND
#define __KIND __DOUBLE_PRECISION
#include "ppm_gmm_march_bkwd.f"
#undef __KIND
#undef __DIM

      END MODULE ppm_module_gmm_march_bkwd
