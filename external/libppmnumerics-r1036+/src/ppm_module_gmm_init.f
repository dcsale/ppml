      !-------------------------------------------------------------------------
      !  Module       :               ppm_module_gmm_init
      !-------------------------------------------------------------------------
      !
      !  Purpose      : This module includes the source code for the 
      !                 initialization routine of the marching method.
      !
      !  Remarks      :
      !
      !  References   :
      !
      !  Revisions    :
      !-------------------------------------------------------------------------
      !  $Log: ppm_module_gmm_init.f,v $
      !  Revision 1.1.1.1  2007/07/13 10:18:58  ivos
      !  CBL version of the PPM library
      !
      !  Revision 1.1  2005/03/10 01:37:15  ivos
      !  Initial check-in. BEWARE: Not tested in parallel yet!
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
     
      MODULE ppm_module_gmm_init

         !----------------------------------------------------------------------
         !  Define interfaces to ppm_gmm_init
         !----------------------------------------------------------------------
         INTERFACE ppm_gmm_init
            MODULE PROCEDURE ppm_gmm_init
         END INTERFACE

         !----------------------------------------------------------------------
         !  include the source 
         !----------------------------------------------------------------------
        CONTAINS

#include "ppm_gmm_init.f"

      END MODULE ppm_module_gmm_init
