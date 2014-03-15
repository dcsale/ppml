      !-------------------------------------------------------------------------
      !  Function     :                 ppm_ode_alldone
      !-------------------------------------------------------------------------
      !
      !  Purpose      : will say .true. if all modes are integrated 
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
      !  $Log: ppm_ode_alldone.f,v $
      !  Revision 1.1.1.1  2007/07/13 10:19:00  ivos
      !  CBL version of the PPM library
      !
      !  Revision 1.6  2004/07/26 11:59:40  ivos
      !  Fixes to make it compile.
      !
      !  Revision 1.5  2004/07/26 11:33:03  michaebe
      !  inserted the use of the ppm ode data module.
      !
      !  Revision 1.4  2004/07/26 11:28:13  michaebe
      !  syntax error... removed.
      !
      !  Revision 1.3  2004/07/26 07:46:51  michaebe
      !  Atomized. Otherwise no changes.
      !
      !  Revision 1.2  2004/06/10 16:20:03  ivos
      !  Moved all cpp directtives to column 1. The NEC cpp did not recognize
      !  them otherwise!!!
      !
      !  Revision 1.1  2004/02/19 08:33:53  michaebe
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
      FUNCTION ppm_ode_alldone(info)
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
      IMPLICIT NONE
      !-----------------------------------------------------------------------
      !  Arguments
      !-----------------------------------------------------------------------
      INTEGER,                    INTENT(  OUT) :: info
      LOGICAL                                   :: ppm_ode_alldone
      !-----------------------------------------------------------------------
      !  call substart
      !-----------------------------------------------------------------------
      CALL substart('ppm_ode_alldone',t0,info)
        IF(ppm_debug.GT.0) THEN
           IF(.NOT.ppm_initialized) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_ppm_noinit,'ppm_ode_alldone',&
                   & 'Please call ppm_init first!',__LINE__,info)
              GOTO 9999
           END IF
        END IF
        IF(ANY(ppm_ode_state(1:ppm_max_mid) .NE. ppm_ode_state_finished)) THEN
           ppm_ode_alldone = .FALSE.
        ELSE
           ppm_ode_alldone = .TRUE.
        END IF
9999    CONTINUE        
      !-----------------------------------------------------------------------
      ! substop
      !-----------------------------------------------------------------------
      CALL substop('ppm_ode_alldone',t0,info)
      RETURN

      END FUNCTION ppm_ode_alldone
