      !-------------------------------------------------------------------------
      ! Module         :            ppm_module_mg_finalize
      !-------------------------------------------------------------------------
      !
      ! Purpose       : multigrid finalize module
      !               
      !
      ! Remarks       :
      !
      ! References    : 
      !
      ! Revisions     :
      !-------------------------------------------------------------------------
      !  $Log: ppm_module_mg_finalize.f,v $
      !  Revision 1.1.1.1  2007/07/13 10:19:00  ivos
      !  CBL version of the PPM library
      !
      !  Revision 1.1  2004/09/22 18:49:34  kotsalie
      !  MG new version
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


#define __SINGLE_PRECISION 1
#define __DOUBLE_PRECISION 2
#define __INTEGER          3
#define __LOGICAL          4
#define __2D               7
#define __3D               8
#define __SFIELD           9
#define __VFIELD          10

MODULE ppm_module_mg_finalize   
  !--------------------------------------------------------------------------
  !Modules
  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  ! INCLUDE THE SOURCES
  !-----------------------------------------------------------------------------

CONTAINS
#define __DIM __SFIELD
#define __MESH_DIM __2D
#define __KIND __SINGLE_PRECISION
#include "ppm_mg_finalize.f"
#undef __KIND

#define __KIND __DOUBLE_PRECISION
#include "ppm_mg_finalize.f"
#undef __KIND
#undef __MESH_DIM

#define __MESH_DIM __3D
#define __KIND __SINGLE_PRECISION
#include "ppm_mg_finalize.f"
#undef __KIND

#define __KIND __DOUBLE_PRECISION
#include "ppm_mg_finalize.f"
#undef __KIND
#undef __MESH_DIM
#undef __DIM

#define __DIM __VFIELD
#define __MESH_DIM __2D
#define __KIND __SINGLE_PRECISION
#include "ppm_mg_finalize.f"
#undef __KIND

#define __KIND __DOUBLE_PRECISION
#include "ppm_mg_finalize.f"
#undef __KIND
#undef __MESH_DIM

#define __MESH_DIM __3D
#define __KIND __SINGLE_PRECISION
#include "ppm_mg_finalize.f"
#undef __KIND

#define __KIND __DOUBLE_PRECISION
#include "ppm_mg_finalize.f"
#undef __KIND
#undef __MESH_DIM
#undef __DIM



END MODULE ppm_module_mg_finalize


