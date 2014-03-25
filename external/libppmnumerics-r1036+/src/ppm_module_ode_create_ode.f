      !-------------------------------------------------------------------------
      !  Module       :            ppm_module_ode_create_ode
      !-------------------------------------------------------------------------
      !
      !  Purpose      : procedure module for ppm_ode_create_ode
      !
      !  Remarks      : 
      !
      !  References   : 
      !
      !  Revisions    :
      !-------------------------------------------------------------------------
      !  $Log: ppm_module_ode_create_ode.f,v $
      !  Revision 1.1.1.1  2007/07/13 10:19:00  ivos
      !  CBL version of the PPM library
      !
      !  Revision 1.2  2004/07/26 11:20:29  michaebe
      !  added dummy interface
      !
      !  Revision 1.1  2004/07/26 07:45:47  michaebe
      !  Procedure modules created in the course of atomization.
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


      MODULE ppm_module_ode_create_ode

        !-----------------------------------------------------
        !  Dummy interface
        !-----------------------------------------------------
        INTERFACE ppm_ode_create_ode
           MODULE PROCEDURE ppm_ode_create_ode
        END INTERFACE
        
      CONTAINS

#include "ppm_ode_create_ode.f"

      END MODULE ppm_module_ode_create_ode
