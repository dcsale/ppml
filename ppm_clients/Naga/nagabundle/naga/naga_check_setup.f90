!------------------------------------------------------------------------------
! Subroutine :  naga_check_setup.f90
!------------------------------------------------------------------------------
!      _  __
!     / |/ /__ ____ ____ _
!    /    / _ `/ _ `/ _ `/
!   /_/|_/\_,_/\_, /\_,_/
!             /___/
!   A PPM Vortex-In-Cell client - 2011 (c)
!   Naga is distributed under the GNU Lesser General Public License version 3
!   Johannes Tophøj Rasmussen, johannes.tophoj@gmail.com
!   Technical University of Denmark
!   Department of Mechanical Engineering
!
!
! This routines checks the input parameters, e.g. if the initialised 
! domain,data etc. is valid
! NOTE: This routine needs to be improved extensively
!------------------------------------------------------------------------------

MODULE naga_mod_check_setup

IMPLICIT NONE

INTERFACE naga_check_setup
  MODULE PROCEDURE naga_check_setup
END INTERFACE

CONTAINS

SUBROUTINE naga_check_setup(info)

USE naga_mod_globals
USE naga_mod_say

IMPLICIT NONE
!-------------------------------------------------------------------------
! Arguments
!-------------------------------------------------------------------------
INTEGER, INTENT(OUT) :: info


!-------------------------------------------------------------------------
! Local variables
!-------------------------------------------------------------------------
REAL(MK)             :: t0


!----------------------------------------------------------------------------
! Initialise routine
!----------------------------------------------------------------------------
CALL substart('naga_check_setup',t0,info)


!----------------------------------------------------------------------------
! Check extent of patches (min vs max) and if they overlap bounds of parents
!----------------------------------------------------------------------------
IF (nlevels .LE. 0) THEN
  CALL naga_say(rank, 'naga_check_setup','No patches have been setup.')
  info = -1
  GOTO 9999
ENDIF
IF (maxpatches .LE. 0) THEN
  CALL naga_say(rank, 'naga_check_setup','No patches have been setup.')
  info = -1
  GOTO 9999
ENDIF


!----------------------------------------------------------------------------
! Check if the validation field is toggled for the validation test cases
!----------------------------------------------------------------------------
IF (validationfield .NE. 1) THEN
  IF (flowcase .EQ. 301 .OR. &
      flowcase .EQ. 302 .OR. &
      flowcase .EQ. 303 .OR. &
      flowcase .EQ. 304) THEN
      CALL naga_say(rank, 'naga_check_setup','The VALIDATIONFIELD must be set 1 for the validation flowcases.')
      info = -1
      GOTO 9999
  ENDIF
ENDIF



!----------------------------------------------------------------------------
! Return
!----------------------------------------------------------------------------
 9999 CONTINUE
CALL substop('naga_check_setup',t0,info)
RETURN


END SUBROUTINE naga_check_setup

END MODULE naga_mod_check_setup

