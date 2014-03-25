      !-------------------------------------------------------------------------
      !  Module       :             ppm_module_comp_pp_mk_table
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
#define __SINGLE_PRECISION_COMPLEX 3
#define __DOUBLE_PRECISION_COMPLEX 4

#define __INTERNAL                 5
#define __USER_FUNCTION            6

      MODULE ppm_module_comp_pp_mk_table
      !!! [[ppm_comp_pp_mk_table]]
      !!! This module provides the routines concerned
      !!! with computing kernel interactions within a set of particles.
         !----------------------------------------------------------------------
         !  Define interfaces to the lookup table creation routine
         !----------------------------------------------------------------------
         INTERFACE ppm_comp_pp_mk_table
            MODULE PROCEDURE ppm_comp_pp_mk_table_si
            MODULE PROCEDURE ppm_comp_pp_mk_table_di
            MODULE PROCEDURE ppm_comp_pp_mk_table_sci
            MODULE PROCEDURE ppm_comp_pp_mk_table_dci
            MODULE PROCEDURE ppm_comp_pp_mk_table_su
            MODULE PROCEDURE ppm_comp_pp_mk_table_du
            MODULE PROCEDURE ppm_comp_pp_mk_table_scu
            MODULE PROCEDURE ppm_comp_pp_mk_table_dcu
         END INTERFACE

         !----------------------------------------------------------------------
         !  include the source 
         !----------------------------------------------------------------------
         CONTAINS

#define __KERNEL __INTERNAL
#define __KIND __SINGLE_PRECISION
#include "ppm_comp_pp_mk_table.f"
#undef __KIND

#define __KIND __DOUBLE_PRECISION
#include "ppm_comp_pp_mk_table.f"
#undef __KIND

#define __KIND __SINGLE_PRECISION_COMPLEX
#include "ppm_comp_pp_mk_table.f"
#undef __KIND

#define __KIND __DOUBLE_PRECISION_COMPLEX
#include "ppm_comp_pp_mk_table.f"
#undef __KIND
#undef __KERNEL

#define __KERNEL __USER_FUNCTION
#define __KIND __SINGLE_PRECISION
#include "ppm_comp_pp_mk_table.f"
#undef __KIND

#define __KIND __DOUBLE_PRECISION
#include "ppm_comp_pp_mk_table.f"
#undef __KIND

#define __KIND __SINGLE_PRECISION_COMPLEX
#include "ppm_comp_pp_mk_table.f"
#undef __KIND

#define __KIND __DOUBLE_PRECISION_COMPLEX
#include "ppm_comp_pp_mk_table.f"
#undef __KIND
#undef __KERNEL

      END MODULE ppm_module_comp_pp_mk_table
