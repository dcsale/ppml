!------------------------------------------------------------------------------
! Subroutine :  naga_structure_none.f90
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
! This routines sets up 0 chi/ubar fields
!------------------------------------------------------------------------------

MODULE naga_mod_structure_none

IMPLICIT NONE

INTERFACE naga_structure_none
  MODULE PROCEDURE naga_structure_none
END INTERFACE

CONTAINS

SUBROUTINE naga_structure_none(info)

USE naga_mod_globals

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
CALL substart('naga_structure_none',t0,info)


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
            chif(ilevel,ipatch)%fld(i,j,k,isub)    = 0.0_MK
            ubarf(ilevel,ipatch)%fld(1,i,j,k,isub) = 0.0_MK
            ubarf(ilevel,ipatch)%fld(2,i,j,k,isub) = 0.0_MK
            ubarf(ilevel,ipatch)%fld(3,i,j,k,isub) = 0.0_MK
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
CALL substop('naga_structure_none',t0,info)
RETURN


END SUBROUTINE naga_structure_none

END MODULE naga_mod_structure_none

