      !-------------------------------------------------------------------------
      !  Module       :             ppm_module_gmm_reinitialize
      !-------------------------------------------------------------------------
      !
      !  Purpose      : This module includes the source code for the 
      !                 reinitialize routine of the marching method.
      !
      !  Remarks      :
      !
      !  References   :
      !
      !  Revisions    :
      !-------------------------------------------------------------------------
      !  $Log: ppm_module_gmm_reinitialize.f,v $
      !  Revision 1.1.1.1  2007/07/13 10:18:58  ivos
      !  CBL version of the PPM library
      !
      !  Revision 1.3  2005/04/27 01:06:14  ivos
      !  Convergence tests completed, cleaned up code, optmized code (Shark),
      !  and changed structure to allow faster compilation.
      !
      !  Revision 1.2  2005/04/21 04:48:23  ivos
      !  Cleaned interfaces and removed unnecessary overloaded versions.
      !
      !  Revision 1.1  2005/03/11 04:16:00  ivos
      !  Initial implementation.
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

      MODULE ppm_module_gmm_reinitialize

         !----------------------------------------------------------------------
         !  Define interfaces to ppm_gmm_reinitialize
         !----------------------------------------------------------------------
         INTERFACE ppm_gmm_reinitialize
            MODULE PROCEDURE ppm_gmm_reinitialize_2ds
            MODULE PROCEDURE ppm_gmm_reinitialize_2dd
            MODULE PROCEDURE ppm_gmm_reinitialize_3ds
            MODULE PROCEDURE ppm_gmm_reinitialize_3dd
         END INTERFACE

         !----------------------------------------------------------------------
         !  include the source 
         !----------------------------------------------------------------------
         CONTAINS

#define __DIM __2D
#define __KIND __SINGLE_PRECISION
#include "ppm_gmm_reinitialize.f"
#undef __KIND
#define __KIND __DOUBLE_PRECISION
#include "ppm_gmm_reinitialize.f"
#undef __KIND
#undef __DIM

#define __DIM __3D
#define __KIND __SINGLE_PRECISION
#include "ppm_gmm_reinitialize.f"
#undef __KIND
#define __KIND __DOUBLE_PRECISION
#include "ppm_gmm_reinitialize.f"
#undef __KIND
#undef __DIM

      END MODULE ppm_module_gmm_reinitialize
