  !------------------------------------------------------------------------------
! Subroutine :  naga_add_uinfinity.f90
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
! This routines adds u infinity to the velocity field (imposes the free stream)
!------------------------------------------------------------------------------

MODULE naga_mod_add_uinfinity

IMPLICIT NONE

INTERFACE naga_add_uinfinity
  MODULE PROCEDURE naga_add_uinfinity
END INTERFACE

CONTAINS

SUBROUTINE naga_add_uinfinity(info)

USE naga_mod_globals
USE naga_mod_say

IMPLICIT NONE

!-------------------------------------------------------------------------
! Arguments
!-------------------------------------------------------------------------
INTEGER, INTENT(INOUT)                                  :: info


!-------------------------------------------------------------------------
! Local variables
!-------------------------------------------------------------------------
REAL(MK)                        :: t0
INTEGER                         :: ilevel,ipatch
INTEGER                         :: isub,isubl
INTEGER                         :: i,j,k


!----------------------------------------------------------------------------
! Initialise routine
!----------------------------------------------------------------------------
CALL substart('naga_add_uinfinity',t0,info)


!----------------------------------------------------------------------------
! Calculate velocity on mesh
!----------------------------------------------------------------------------
DO ilevel=1,nlevels
  DO ipatch=1,npatches(ilevel)
    DO isub=1,topos(ilevel)%nsublist
      isubl=topos(ilevel)%isublist(isub)
      DO k=1-gstw(3),ptcset(ilevel,ipatch)%snx(3,isubl)+gstw(3)
        DO j=1-gstw(2),ptcset(ilevel,ipatch)%snx(2,isubl)+gstw(2)
          DO i=1-gstw(1),ptcset(ilevel,ipatch)%snx(1,isubl)+gstw(1)
            uf(ilevel,ipatch)%fld(1,i,j,k,isub) = &
                         & uf(ilevel,ipatch)%fld(1,i,j,k,isub) + uinfinity(1)
            uf(ilevel,ipatch)%fld(2,i,j,k,isub) = &
                         & uf(ilevel,ipatch)%fld(2,i,j,k,isub) + uinfinity(2)
            uf(ilevel,ipatch)%fld(3,i,j,k,isub) = &
                         & uf(ilevel,ipatch)%fld(3,i,j,k,isub) + uinfinity(3)
          ENDDO
        ENDDO
      ENDDO
    ENDDO
  ENDDO
ENDDO


!----------------------------------------------------------------------------
! Return
!----------------------------------------------------------------------------
 9999 CONTINUE
CALL substop('naga_add_uinfinity',t0,info)
RETURN


END SUBROUTINE naga_add_uinfinity

END MODULE naga_mod_add_uinfinity

