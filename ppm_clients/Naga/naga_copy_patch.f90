!------------------------------------------------------------------------------
! Subroutine :  naga_copy_patch.f90
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
! This routine copies the geometry settings of one patch setup to another setup
! which has already been allocated
!------------------------------------------------------------------------------

MODULE naga_mod_copy_patch

IMPLICIT NONE

INTERFACE naga_copy_patch
  MODULE PROCEDURE naga_copy_patch
END INTERFACE

CONTAINS

SUBROUTINE naga_copy_patch(psetin,psetout,info)

USE naga_mod_globals
USE naga_mod_say

IMPLICIT NONE
!-------------------------------------------------------------------------
! Arguments
!-------------------------------------------------------------------------
TYPE(patch_setup),DIMENSION(:,:),POINTER               :: psetin
TYPE(patch_setup),DIMENSION(:,:),POINTER               :: psetout
INTEGER, INTENT(OUT)                        :: info


!-------------------------------------------------------------------------
! Local variables
!-------------------------------------------------------------------------
REAL(MK)                        :: t0
INTEGER                         :: ilevel,ipatch


!----------------------------------------------------------------------------
! Initialise routine
!----------------------------------------------------------------------------
CALL substart('naga_copy_patch',t0,info)


!-------------------------------------------------------------------------
! We now know the number of patches and levels. First declare the 
! patch_setting array as
! [number of levels] x [maximum number of patches on one level] but
! we only allocate the data of the patches that do exist
!-------------------------------------------------------------------------
IF(ASSOCIATED(psetout)) THEN
  DEALLOCATE(psetout,stat=info)
  IF (info .NE. 0) THEN
    CALL naga_say(rank, 'naga_copy_patch','Failed to deallocate patch.')
    GOTO 9999
  ENDIF
END IF

ALLOCATE(psetout(nlevels,maxpatches))

DO ilevel=1,nlevels
  DO ipatch=1,npatches(ilevel)
    psetout(ilevel,ipatch)%lvl = psetin(ilevel,ipatch)%lvl
    psetout(ilevel,ipatch)%par = psetin(ilevel,ipatch)%par
    psetout(ilevel,ipatch)%buf = psetin(ilevel,ipatch)%buf
    !min of patch
    psetout(ilevel,ipatch)%min(1) = psetin(ilevel,ipatch)%min(1)
    psetout(ilevel,ipatch)%min(2) = psetin(ilevel,ipatch)%min(2)
    psetout(ilevel,ipatch)%min(3) = psetin(ilevel,ipatch)%min(3)
    !max of patch
    psetout(ilevel,ipatch)%max(1) = psetin(ilevel,ipatch)%max(1)
    psetout(ilevel,ipatch)%max(2) = psetin(ilevel,ipatch)%max(2)
    psetout(ilevel,ipatch)%max(3) = psetin(ilevel,ipatch)%max(3)
    !resolution
    psetout(ilevel,ipatch)%gnx     = psetin(ilevel,ipatch)%gnx
    !grid spacing
    psetout(ilevel,ipatch)%dx(1)  = psetin(ilevel,ipatch)%dx(1)
    psetout(ilevel,ipatch)%dx(2)  = psetin(ilevel,ipatch)%dx(2)
    psetout(ilevel,ipatch)%dx(3)  = psetin(ilevel,ipatch)%dx(3)
  ENDDO
ENDDO


!----------------------------------------------------------------------------
! Return
!----------------------------------------------------------------------------
 9999 CONTINUE
CALL substop('naga_copy_patch',t0,info)
RETURN


END SUBROUTINE naga_copy_patch

END MODULE naga_mod_copy_patch

