#ifdef __XLF
@PROCESS NOHOT
#endif
      !-------------------------------------------------------------------------
      !  Module       :            ppm_module_util_fft_backward
      !-------------------------------------------------------------------------
      !
      !  Purpose      : This module includes the source code for the utility
      !                 routines.
      !
      !  Remarks      :
      !
      !  References   :
      !
      !  Revisions    :
      !-------------------------------------------------------------------------
      !  $Log: ppm_module_util_fft_backward.f,v $
      !  Revision 1.1.1.1  2007/07/13 10:19:00  ivos
      !  CBL version of the PPM library
      !
      !  Revision 1.2  2004/11/05 18:16:38  michaebe
      !  added xlf compiler directives
      !
      !  Revision 1.1  2004/07/26 07:30:12  ivos
      !  First commit after spitting the old modules into single-interface
      !  units.
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

      MODULE ppm_module_util_fft_backward

         !----------------------------------------------------------------------
         !  Define interface to ppm_util_fft_backward
         !----------------------------------------------------------------------
         INTERFACE ppm_util_fft_backward
            MODULE PROCEDURE ppm_util_fft_backward_2ds
            MODULE PROCEDURE ppm_util_fft_backward_2dd
            MODULE PROCEDURE ppm_util_fft_backward_2dc
            MODULE PROCEDURE ppm_util_fft_backward_2dcc

            MODULE PROCEDURE ppm_util_fft_backward_3ds
            MODULE PROCEDURE ppm_util_fft_backward_3dd
            MODULE PROCEDURE ppm_util_fft_backward_3dc
            MODULE PROCEDURE ppm_util_fft_backward_3dcc
         END INTERFACE

         !----------------------------------------------------------------------
         !  include the source 
         !----------------------------------------------------------------------
         CONTAINS
 
#define __KIND __SINGLE_PRECISION
#include "ppm_util_fft_backward_2d.f"
#include "ppm_util_fft_backward_3d.f"
#undef __KIND

#define __KIND __DOUBLE_PRECISION
#include "ppm_util_fft_backward_2d.f"
#include "ppm_util_fft_backward_3d.f"
#undef __KIND

#define __KIND __SINGLE_PRECISION_COMPLEX
#include "ppm_util_fft_backward_2d.f"
#include "ppm_util_fft_backward_3d.f"
#undef __KIND

#define __KIND __DOUBLE_PRECISION_COMPLEX
#include "ppm_util_fft_backward_2d.f"
#include "ppm_util_fft_backward_3d.f"
#undef __KIND

      END MODULE ppm_module_util_fft_backward
