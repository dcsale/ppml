!------------------------------------------------------------------------------
! Subroutine :  naga_flowcase.f90
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
! This routines performs various actions during the simulation based on the
! flow case and ipos (determines position in the algorithm). So far only
! ipos = 1 has been implemented corresponding to the beginning of the RHS
! evaluation.
!------------------------------------------------------------------------------

MODULE naga_mod_flowcase

IMPLICIT NONE

INTERFACE naga_flowcase
  MODULE PROCEDURE naga_flowcase
END INTERFACE

CONTAINS

SUBROUTINE naga_flowcase(ipos,info)

USE naga_mod_globals
USE naga_mod_say

IMPLICIT NONE

!-------------------------------------------------------------------------
! Arguments
!-------------------------------------------------------------------------
INTEGER, INTENT(IN   )                                  :: ipos
INTEGER, INTENT(  OUT)                                  :: info


!-------------------------------------------------------------------------
! Local variables
!-------------------------------------------------------------------------
REAL(MK),PARAMETER              :: PI=ACOS(-1.0_mk)
REAL(MK)                        :: t0
INTEGER                         :: ilevel,ipatch


!----------------------------------------------------------------------------
! Initialise routine
!----------------------------------------------------------------------------
CALL substart('naga_flowcase',t0,info)
info = -1

!position before the RHS begins
IF (ipos .EQ. 1) THEN

  !Flowcase 32: disrupt y-component of uinfinity
  IF (flowcase .EQ. 32) THEN
    IF (time .GT. 3.0_MK .AND. time .LT. 4.0_MK) THEN
      uinfinity(2) = SIN(PI*(time-3.0_MK))
    ELSE
      uinfinity(2) = 0.0_MK
    ENDIF
  ENDIF


ELSE
  CALL naga_say(rank,'naga_flowcase','Undefined stage provided.')
  GOTO 9999
ENDIF



info = 0
!----------------------------------------------------------------------------
! Return
!----------------------------------------------------------------------------
 9999 CONTINUE
CALL substop('naga_flowcase',t0,info)
RETURN


END SUBROUTINE naga_flowcase

END MODULE naga_mod_flowcase

