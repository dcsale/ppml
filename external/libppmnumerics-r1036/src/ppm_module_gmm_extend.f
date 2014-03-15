      !-------------------------------------------------------------------------
      !  Module       :               ppm_module_gmm_extend
      !-------------------------------------------------------------------------
      !
      !  Purpose      : This module includes the source code for the 
      !                 extend routine of the marching method.
      !
      !  Remarks      :
      !
      !  References   :
      !
      !  Revisions    :
      !-------------------------------------------------------------------------
      !  $Log: ppm_module_gmm_extend.f,v $
      !  Revision 1.1.1.1  2007/07/13 10:18:58  ivos
      !  CBL version of the PPM library
      !
      !  Revision 1.4  2005/04/21 04:48:24  ivos
      !  Cleaned interfaces and removed unnecessary overloaded versions.
      !
      !  Revision 1.3  2005/03/12 04:08:36  ivos
      !  Misc bug fixes.
      !
      !  Revision 1.2  2005/03/11 21:10:08  ivos
      !  Added thresholded extensions and closest point transform.
      !
      !  Revision 1.1  2005/03/11 04:15:59  ivos
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
#define __VFIELD                   3
#define __SFIELD                   4
#define __2D                       5
#define __3D                       6
#define __YES                      7
#define __NO                       8

      MODULE ppm_module_gmm_extend

         !----------------------------------------------------------------------
         !  Modules
         !----------------------------------------------------------------------
         USE ppm_module_data, ONLY: ppm_kind_single, ppm_kind_double 
         PRIVATE :: ppm_kind_single, ppm_kind_double
         
         !----------------------------------------------------------------------
         !  Work memory
         !----------------------------------------------------------------------
         REAL(ppm_kind_single), DIMENSION(:,:,:  ), POINTER :: ext_wrk_2ds
         REAL(ppm_kind_single), DIMENSION(:,:,:,:), POINTER :: ext_wrk_3ds
         REAL(ppm_kind_double), DIMENSION(:,:,:  ), POINTER :: ext_wrk_2dd
         REAL(ppm_kind_double), DIMENSION(:,:,:,:), POINTER :: ext_wrk_3dd

         !----------------------------------------------------------------------
         !  Define interfaces to ppm_gmm_extend
         !----------------------------------------------------------------------
         INTERFACE ppm_gmm_extend
            ! 2d scalar and vector field versions with cutoff
            MODULE PROCEDURE ppm_gmm_extend_2d_tsca_s
            MODULE PROCEDURE ppm_gmm_extend_2d_tsca_d
            MODULE PROCEDURE ppm_gmm_extend_2d_tvec_s
            MODULE PROCEDURE ppm_gmm_extend_2d_tvec_d

            ! 3d scalar and vector field versions with cutoff
            MODULE PROCEDURE ppm_gmm_extend_3d_tsca_s
            MODULE PROCEDURE ppm_gmm_extend_3d_tsca_d
            MODULE PROCEDURE ppm_gmm_extend_3d_tvec_s
            MODULE PROCEDURE ppm_gmm_extend_3d_tvec_d

            ! 2d scalar and vector field versions with function pointer
            MODULE PROCEDURE ppm_gmm_extend_2d_ksca_s
            MODULE PROCEDURE ppm_gmm_extend_2d_ksca_d
            MODULE PROCEDURE ppm_gmm_extend_2d_kvec_s
            MODULE PROCEDURE ppm_gmm_extend_2d_kvec_d

            ! 3d scalar and vector field versions with function pointer
            MODULE PROCEDURE ppm_gmm_extend_3d_ksca_s
            MODULE PROCEDURE ppm_gmm_extend_3d_ksca_d
            MODULE PROCEDURE ppm_gmm_extend_3d_kvec_s
            MODULE PROCEDURE ppm_gmm_extend_3d_kvec_d
         END INTERFACE

         !----------------------------------------------------------------------
         !  include the source 
         !----------------------------------------------------------------------
         CONTAINS

#define __KICKOFF __NO
#define __DIM __2D
#define __TYPE __SFIELD
#define __KIND __SINGLE_PRECISION
#include "ppm_gmm_extend.f"
#undef __KIND
#define __KIND __DOUBLE_PRECISION
#include "ppm_gmm_extend.f"
#undef __KIND
#undef __TYPE

#define __TYPE __VFIELD
#define __KIND __SINGLE_PRECISION
#include "ppm_gmm_extend.f"
#undef __KIND
#define __KIND __DOUBLE_PRECISION
#include "ppm_gmm_extend.f"
#undef __KIND
#undef __TYPE
#undef __DIM

#define __DIM __3D
#define __TYPE __SFIELD
#define __KIND __SINGLE_PRECISION
#include "ppm_gmm_extend.f"
#undef __KIND
#define __KIND __DOUBLE_PRECISION
#include "ppm_gmm_extend.f"
#undef __KIND
#undef __TYPE

#define __TYPE __VFIELD
#define __KIND __SINGLE_PRECISION
#include "ppm_gmm_extend.f"
#undef __KIND
#define __KIND __DOUBLE_PRECISION
#include "ppm_gmm_extend.f"
#undef __KIND
#undef __TYPE
#undef __DIM
#undef __KICKOFF

#define __KICKOFF __YES
#define __DIM __2D
#define __TYPE __SFIELD
#define __KIND __SINGLE_PRECISION
#include "ppm_gmm_extend.f"
#undef __KIND
#define __KIND __DOUBLE_PRECISION
#include "ppm_gmm_extend.f"
#undef __KIND
#undef __TYPE

#define __TYPE __VFIELD
#define __KIND __SINGLE_PRECISION
#include "ppm_gmm_extend.f"
#undef __KIND
#define __KIND __DOUBLE_PRECISION
#include "ppm_gmm_extend.f"
#undef __KIND
#undef __TYPE
#undef __DIM

#define __DIM __3D
#define __TYPE __SFIELD
#define __KIND __SINGLE_PRECISION
#include "ppm_gmm_extend.f"
#undef __KIND
#define __KIND __DOUBLE_PRECISION
#include "ppm_gmm_extend.f"
#undef __KIND
#undef __TYPE

#define __TYPE __VFIELD
#define __KIND __SINGLE_PRECISION
#include "ppm_gmm_extend.f"
#undef __KIND
#define __KIND __DOUBLE_PRECISION
#include "ppm_gmm_extend.f"
#undef __KIND
#undef __TYPE
#undef __DIM
#undef __KICKOFF

      END MODULE ppm_module_gmm_extend
