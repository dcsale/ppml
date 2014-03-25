#ifdef __XLF
@PROCESS NOHOT
#endif
      !-------------------------------------------------------------------------
      !  Module       :            ppm_module_fdsolver_fft_bd
      !-------------------------------------------------------------------------
      !
      !  Purpose      : This module includes the source code for the fft
      !                 routines
      !
      !  Remarks      :
      !
      !  References   :
      !
      !  Revisions    :
      !-------------------------------------------------------------------------
      !  $Log: ppm_module_fdsolver_fft_bd.f,v $
      !  Revision 1.1.1.1  2007/07/13 10:18:58  ivos
      !  CBL version of the PPM library
      !
      !  Revision 1.1  2005/02/16 12:07:35  hiebers
      !  initial implementation
      !
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
#define __SINGLE_PRECISION           1
#define __DOUBLE_PRECISION           2
#define __SINGLE_PRECISION_COMPLEX   5
#define __DOUBLE_PRECISION_COMPLEX   6
#define __SINGLE_PRECISION_COMPLEX_Z 7
#define __DOUBLE_PRECISION_COMPLEX_Z 8

#define __SLAB                      10


      MODULE ppm_module_fdsolver_fft_bd

         !----------------------------------------------------------------------
         !  Define interface to ppm_fdsolver_fft_bd
         !----------------------------------------------------------------------
         INTERFACE ppm_fdsolver_fft_bd_slab
            MODULE PROCEDURE ppm_fdsolver_fft_bd_slab_3ds
            MODULE PROCEDURE ppm_fdsolver_fft_bd_slab_3dd
         END INTERFACE

         INTERFACE ppm_fdsolver_fft_bd
            MODULE PROCEDURE ppm_fdsolver_fft_bd_2ds
            MODULE PROCEDURE ppm_fdsolver_fft_bd_2dd
            MODULE PROCEDURE ppm_fdsolver_fft_bd_2dc
            MODULE PROCEDURE ppm_fdsolver_fft_bd_2dcc

            MODULE PROCEDURE ppm_fdsolver_fft_bd_3ds
            MODULE PROCEDURE ppm_fdsolver_fft_bd_3dd
            MODULE PROCEDURE ppm_fdsolver_fft_bd_3dc
            MODULE PROCEDURE ppm_fdsolver_fft_bd_3dcc
         END INTERFACE

         INTERFACE ppm_fdsolver_fft_bd_z
            MODULE PROCEDURE ppm_fdsolver_fft_bd_z_3dc
            MODULE PROCEDURE ppm_fdsolver_fft_bd_z_3dcc
         END INTERFACE

         !----------------------------------------------------------------------
         !  include the source 
         !----------------------------------------------------------------------
         CONTAINS

#define __CASE __SLAB
 
#define __KIND __SINGLE_PRECISION
#include "ppm_fdsolver_fft_bd_3d.f"
#undef __KIND

#define __KIND __DOUBLE_PRECISION
#include "ppm_fdsolver_fft_bd_3d.f"
#undef __KIND

#undef __CASE


#define __KIND __SINGLE_PRECISION
#include "ppm_fdsolver_fft_bd_2d.f"
#include "ppm_fdsolver_fft_bd_3d.f"
#undef __KIND

#define __KIND __DOUBLE_PRECISION
#include "ppm_fdsolver_fft_bd_2d.f"
#include "ppm_fdsolver_fft_bd_3d.f"
#undef __KIND

#define __KIND __SINGLE_PRECISION_COMPLEX
#include "ppm_fdsolver_fft_bd_2d.f"
#include "ppm_fdsolver_fft_bd_3d.f"
#undef __KIND

#define __KIND __DOUBLE_PRECISION_COMPLEX
#include "ppm_fdsolver_fft_bd_2d.f"
#include "ppm_fdsolver_fft_bd_3d.f"
#undef __KIND

#define __KIND __SINGLE_PRECISION_COMPLEX_Z
#include "ppm_fdsolver_fft_bd_3d.f"
#undef __KIND

#define __KIND __DOUBLE_PRECISION_COMPLEX_Z
#include "ppm_fdsolver_fft_bd_3d.f"
#undef __KIND

      END MODULE ppm_module_fdsolver_fft_bd
