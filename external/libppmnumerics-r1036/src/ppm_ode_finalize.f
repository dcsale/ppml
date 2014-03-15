      !-------------------------------------------------------------------------
      !  Subroutine   :                 ppm_ode_finalize
      !-------------------------------------------------------------------------
      !
      !  Purpose      : garbage collection of the ode
      !
      !  Input        : 
      !
      !  Output       : info                   (I) return status
      !
      !  Remarks      : 
      !
      !  References   : 
      !
      !  Revisions    :
      !-------------------------------------------------------------------------
      !  $Log: ppm_ode_finalize.f,v $
      !  Revision 1.1.1.1  2007/07/13 10:19:00  ivos
      !  CBL version of the PPM library
      !
      !  Revision 1.7  2004/07/26 15:38:51  ivos
      !  Inserted missing USE statements to resolve undefined references
      !  at link stage.
      !
      !  Revision 1.6  2004/07/26 13:49:18  ivos
      !  Removed Routines sections from the header comment.
      !
      !  Revision 1.5  2004/07/26 11:33:03  michaebe
      !  inserted the use of the ppm ode data module.
      !
      !  Revision 1.4  2004/07/26 07:47:52  michaebe
      !  Atomized. Otherwise no changes.
      !
      !  Revision 1.3  2004/06/10 16:20:03  ivos
      !  Moved all cpp directtives to column 1. The NEC cpp did not recognize
      !  them otherwise!!!
      !
      !  Revision 1.2  2004/02/20 14:55:38  michaebe
      !  added the use of ppm_module_util.
      !
      !  Revision 1.1  2004/02/19 08:33:54  michaebe
      !  initial implementation.
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
      SUBROUTINE ppm_ode_finalize(info)
        !-----------------------------------------------------------------------
        !  Includes
        !-----------------------------------------------------------------------
#include "ppm_define.h"

        !-----------------------------------------------------------------------
        !  Modules
        !-----------------------------------------------------------------------
        USE ppm_module_data_ode
        USE ppm_module_data
        USE ppm_module_substart
        USE ppm_module_substop
        USE ppm_module_error
        USE ppm_module_alloc
        IMPLICIT NONE
        
        !-----------------------------------------------------------------------
        !  Arguments
        !-----------------------------------------------------------------------
        INTEGER,                    INTENT(  out) :: info
        INTEGER                                   :: iopt
        INTEGER, DIMENSION(3)                     :: tlda
        !-----------------------------------------------------------------------
        !  call substart
        !-----------------------------------------------------------------------
        CALL substart('ppm_ode_finalize',t0,info)
        !-----------------------------------------------------------------------
        !  check input arguments
        !-----------------------------------------------------------------------
        IF(ppm_debug.GT.0) THEN
           !--------------------------------------------------------------------
           ! check if ppm is initialized
           !--------------------------------------------------------------------
           IF(.NOT.ppm_initialized) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_ppm_noinit,'ppm_ode_finalize',&
                   & 'Please call ppm_init first!',__LINE__,info)
              GOTO 9999
           END IF
        END IF
        iopt = ppm_param_dealloc
        tlda(1) = ppm_max_mid_allocd
#include "ppm_ode_modalloc.h"        
9999    CONTINUE        
        !-----------------------------------------------------------------------
        ! substop
        !-----------------------------------------------------------------------
        CALL substop('ppm_ode_finalize',t0,info)
        RETURN
      END SUBROUTINE ppm_ode_finalize
