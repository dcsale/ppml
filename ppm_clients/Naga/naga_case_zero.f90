!------------------------------------------------------------------------------
! Subroutine :  naga_case_zero.f90
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
! This routines sets up the initial ZERO vorticity field
!------------------------------------------------------------------------------

MODULE naga_mod_case_zero

IMPLICIT NONE

INTERFACE naga_case_zero
  MODULE PROCEDURE naga_case_zero
END INTERFACE

CONTAINS

SUBROUTINE naga_case_zero(info)

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
INTEGER              :: i,j,k
INTEGER              :: ilevel,ipatch,isubl,isub


!----------------------------------------------------------------------------
! Initialise routine
!----------------------------------------------------------------------------
CALL substart('naga_case_zero',t0,info)


!----------------------------------------------------------------------------
! Loop through the sub domain points
!----------------------------------------------------------------------------
DO ilevel=1,nlevels
  DO ipatch=1,npatches(ilevel)
    DO isub=1,topos(ilevel)%nsublist
      isubl=topos(ilevel)%isublist(isub)
      DO k=1-gstw(1),ptcset(ilevel,ipatch)%snx(3,isubl)+gstw(1)
        DO j=1-gstw(2),ptcset(ilevel,ipatch)%snx(2,isubl)+gstw(2)
          DO i=1-gstw(3),ptcset(ilevel,ipatch)%snx(1,isubl)+gstw(3)
            wf(ilevel,ipatch)%fld(1,i,j,k,isub) = 0.0_MK
            wf(ilevel,ipatch)%fld(2,i,j,k,isub) = 0.0_MK
            wf(ilevel,ipatch)%fld(3,i,j,k,isub) = 0.0_MK
          END DO !i
        END DO !j
      END DO !k
    END DO !isub
  END DO !ipatch
END DO !ilevel


!----------------------------------------------------------------------------
! Return
!----------------------------------------------------------------------------
 9999 CONTINUE
CALL substop('naga_case_zero',t0,info)
RETURN


END SUBROUTINE naga_case_zero

END MODULE naga_mod_case_zero

