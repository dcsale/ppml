!------------------------------------------------------------------------------
! Subroutine : naga_allocate_patch.f90
!------------------------------------------------------------------------------
! _ __
! / |/ /__ ____ ____ _
! / / _ `/ _ `/ _ `/
! /_/|_/\_,_/\_, /\_,_/
! /___/
! A PPM Vortex-In-Cell client - 2011 (c)
! Naga is distributed under the GNU Lesser General Public License version 3
! Johannes Tophøj Rasmussen, johannes.tophoj@gmail.com
! Technical University of Denmark
! Department of Mechanical Engineering
!
!
! This module contains routines for allocating vector and scalar fields
!
! p=particle
! f=field
! s=scalar
! v=vector
! c=complex
!------------------------------------------------------------------------------
MODULE naga_mod_allocate_patch
INTERFACE naga_allocate_patch
  MODULE PROCEDURE naga_allocate_patch_ps
  MODULE PROCEDURE naga_allocate_patch_pv
  MODULE PROCEDURE naga_allocate_patch_fv
  MODULE PROCEDURE naga_allocate_patch_fs
  MODULE PROCEDURE naga_allocate_patch_fvc
  MODULE PROCEDURE naga_allocate_patch_fsc
END INTERFACE
CONTAINS
SUBROUTINE naga_allocate_patch_ps(pset,patch,info)
USE naga_mod_globals
USE naga_mod_say
IMPLICIT NONE
!-----------------------------------------------------------------------------
! Arguments
!-----------------------------------------------------------------------------
TYPE(patch_setup),DIMENSION(:,:),POINTER :: pset
TYPE(patch_part_s),DIMENSION(:,:),POINTER :: patch
INTEGER, INTENT(OUT) :: info
!-----------------------------------------------------------------------------
! Local variables
!-----------------------------------------------------------------------------
REAL(MK) :: t0
INTEGER, DIMENSION(3) :: indl, indu
INTEGER :: ilevel, ipatch, isub, isubl
!-----------------------------------------------------------------------------
! Initialise routine
!-----------------------------------------------------------------------------
CALL substart('naga_allocate_patch',t0,info)
!-----------------------------------------------------------------------------
! Allocate the multi level, multi patch type (not the actual fields)
!-----------------------------------------------------------------------------
ALLOCATE(patch(nlevels,maxpatches))
!-----------------------------------------------------------------------------
! Return
!-----------------------------------------------------------------------------
 9999 CONTINUE
CALL substop('naga_allocate_patch',t0,info)
RETURN
END SUBROUTINE naga_allocate_patch_ps
SUBROUTINE naga_allocate_patch_pv(pset,patch,info)
USE naga_mod_globals
USE naga_mod_say
IMPLICIT NONE
!-----------------------------------------------------------------------------
! Arguments
!-----------------------------------------------------------------------------
TYPE(patch_setup),DIMENSION(:,:),POINTER :: pset
TYPE(patch_part_v),DIMENSION(:,:),POINTER :: patch
INTEGER, INTENT(OUT) :: info
!-----------------------------------------------------------------------------
! Local variables
!-----------------------------------------------------------------------------
REAL(MK) :: t0
INTEGER, DIMENSION(3) :: indl, indu
INTEGER :: ilevel, ipatch, isub, isubl
!-----------------------------------------------------------------------------
! Initialise routine
!-----------------------------------------------------------------------------
CALL substart('naga_allocate_patch',t0,info)
!-----------------------------------------------------------------------------
! Allocate the multi level, multi patch type (not the actual fields)
!-----------------------------------------------------------------------------
ALLOCATE(patch(nlevels,maxpatches))
!-----------------------------------------------------------------------------
! Return
!-----------------------------------------------------------------------------
 9999 CONTINUE
CALL substop('naga_allocate_patch',t0,info)
RETURN
END SUBROUTINE naga_allocate_patch_pv
SUBROUTINE naga_allocate_patch_fs(pset,patch,info)
USE naga_mod_globals
USE naga_mod_say
IMPLICIT NONE
!-----------------------------------------------------------------------------
! Arguments
!-----------------------------------------------------------------------------
TYPE(patch_setup),DIMENSION(:,:),POINTER :: pset
TYPE(patch_field_s),DIMENSION(:,:),POINTER :: patch
INTEGER, INTENT(OUT) :: info
!-----------------------------------------------------------------------------
! Local variables
!-----------------------------------------------------------------------------
REAL(MK) :: t0
INTEGER, DIMENSION(3) :: indl, indu
INTEGER :: ilevel, ipatch, isub, isubl
!-----------------------------------------------------------------------------
! Initialise routine
!-----------------------------------------------------------------------------
CALL substart('naga_allocate_patch',t0,info)
!-----------------------------------------------------------------------------
! Allocate the multi level, multi patch type (not the actual fields)
!-----------------------------------------------------------------------------
ALLOCATE(patch(nlevels,maxpatches))
!-----------------------------------------------------------------------------
! Now loop through all patches at all levels and swing the wand
!-----------------------------------------------------------------------------
DO ilevel=1,nlevels
  DO ipatch=1,npatches(ilevel)
    !-----------------------------------------------------------------------
    ! Deallocate the field if it is already allocated
    !-----------------------------------------------------------------------
    IF(ASSOCIATED(patch(ilevel,ipatch)%fld)) THEN
      DEALLOCATE(patch(ilevel,ipatch)%fld,stat=info)
      IF (info .NE. 0) THEN
        CALL naga_say(rank,'naga_allocate_patch','Failed to deallocate array.')
        GOTO 9999
      ENDIF
    END IF
    indl = 1 - pset(ilevel,ipatch)%gst(1)
    indl = 1 - pset(ilevel,ipatch)%gst(2)
    indl = 1 - pset(ilevel,ipatch)%gst(3)
    indu = 0
    DO isub=1,topos(ilevel)%nsublist
      isubl=topos(ilevel)%isublist(isub)
      !-----------------------------------------------------------------------
      ! Determine indicies for the field
      !-----------------------------------------------------------------------
      indu(1) = MAX(indu(1),&
              & pset(ilevel,ipatch)%snx(1,isubl) + pset(ilevel,ipatch)%gst(1))
      indu(2) = MAX(indu(2),&
              & pset(ilevel,ipatch)%snx(2,isubl) + pset(ilevel,ipatch)%gst(2))
      indu(3) = MAX(indu(3),&
              & pset(ilevel,ipatch)%snx(3,isubl) + pset(ilevel,ipatch)%gst(3))
    ENDDO !isub
    !-----------------------------------------------------------------------
    ! Allocate the field
    !-----------------------------------------------------------------------
    ALLOCATE(patch(ilevel,ipatch)%fld(&
    & indl(1):indu(1),indl(2):indu(2),indl(3):indu(3),&
    & 1:topos(ilevel)%nsublist),stat=info)
    IF (info .NE. 0) THEN
      CALL naga_say(rank,'naga_allocate_patch_fs','Failed to allocate array.')
      GOTO 9999
    ENDIF
  ENDDO !ipatch
ENDDO !ilevel
!-----------------------------------------------------------------------------
! Return
!-----------------------------------------------------------------------------
 9999 CONTINUE
CALL substop('naga_allocate_patch',t0,info)
RETURN
END SUBROUTINE naga_allocate_patch_fs
SUBROUTINE naga_allocate_patch_fv(pset,patch,info)
USE naga_mod_globals
USE naga_mod_say
IMPLICIT NONE
!-----------------------------------------------------------------------------
! Arguments
!-----------------------------------------------------------------------------
TYPE(patch_setup),DIMENSION(:,:),POINTER :: pset
TYPE(patch_field_v),DIMENSION(:,:),POINTER :: patch
INTEGER, INTENT(OUT) :: info
!-----------------------------------------------------------------------------
! Local variables
!-----------------------------------------------------------------------------
REAL(MK) :: t0
INTEGER, DIMENSION(3) :: indl, indu
INTEGER :: ilevel, ipatch, isub, isubl
!-----------------------------------------------------------------------------
! Initialise routine
!-----------------------------------------------------------------------------
CALL substart('naga_allocate_patch',t0,info)
!-----------------------------------------------------------------------------
! Allocate the multi level, multi patch type (not the actual fields)
!-----------------------------------------------------------------------------
ALLOCATE(patch(nlevels,maxpatches))
!-----------------------------------------------------------------------------
! Now loop through all patches at all levels and swing the wand
!-----------------------------------------------------------------------------
DO ilevel=1,nlevels
  DO ipatch=1,npatches(ilevel)
    !-----------------------------------------------------------------------
    ! Deallocate the field if it is already allocated
    !-----------------------------------------------------------------------
    IF(ASSOCIATED(patch(ilevel,ipatch)%fld)) THEN
      DEALLOCATE(patch(ilevel,ipatch)%fld,stat=info)
      IF (info .NE. 0) THEN
        CALL naga_say(rank,'naga_allocate_patch','Failed to deallocate array.')
        GOTO 9999
      ENDIF
    END IF
    indl = 1 - pset(ilevel,ipatch)%gst(1)
    indl = 1 - pset(ilevel,ipatch)%gst(2)
    indl = 1 - pset(ilevel,ipatch)%gst(3)
    indu = 0
    DO isub=1,topos(ilevel)%nsublist
      isubl=topos(ilevel)%isublist(isub)
      !-----------------------------------------------------------------------
      ! Determine indicies for the field
      !-----------------------------------------------------------------------
      indu(1) = MAX(indu(1),&
              & pset(ilevel,ipatch)%snx(1,isubl) + pset(ilevel,ipatch)%gst(1))
      indu(2) = MAX(indu(2),&
              & pset(ilevel,ipatch)%snx(2,isubl) + pset(ilevel,ipatch)%gst(2))
      indu(3) = MAX(indu(3),&
              & pset(ilevel,ipatch)%snx(3,isubl) + pset(ilevel,ipatch)%gst(3))
    ENDDO !isub
    !-----------------------------------------------------------------------
    ! Allocate the field
    !-----------------------------------------------------------------------
    ALLOCATE(patch(ilevel,ipatch)%fld(&
    & ncom,indl(1):indu(1),indl(2):indu(2),indl(3):indu(3),&
    & 1:topos(ilevel)%nsublist),stat=info)
    IF (info .NE. 0) THEN
      CALL naga_say(rank,'naga_allocate_patch_fs','Failed to allocate array.')
      GOTO 9999
    ENDIF
  ENDDO !ipatch
ENDDO !ilevel
!-----------------------------------------------------------------------------
! Return
!-----------------------------------------------------------------------------
 9999 CONTINUE
CALL substop('naga_allocate_patch',t0,info)
RETURN
END SUBROUTINE naga_allocate_patch_fv
SUBROUTINE naga_allocate_patch_fsc(pset,patch,info)
USE naga_mod_globals
USE naga_mod_say
IMPLICIT NONE
!-----------------------------------------------------------------------------
! Arguments
!-----------------------------------------------------------------------------
TYPE(patch_setup),DIMENSION(:,:),POINTER :: pset
TYPE(patch_field_sc),DIMENSION(:,:),POINTER :: patch
INTEGER, INTENT(OUT) :: info
!-----------------------------------------------------------------------------
! Local variables
!-----------------------------------------------------------------------------
REAL(MK) :: t0
INTEGER, DIMENSION(3) :: indl, indu
INTEGER :: ilevel, ipatch, isub, isubl
!-----------------------------------------------------------------------------
! Initialise routine
!-----------------------------------------------------------------------------
CALL substart('naga_allocate_patch',t0,info)
!-----------------------------------------------------------------------------
! Allocate the multi level, multi patch type (not the actual fields)
!-----------------------------------------------------------------------------
ALLOCATE(patch(nlevels,maxpatches))
!-----------------------------------------------------------------------------
! Now loop through all patches at all levels and swing the wand
!-----------------------------------------------------------------------------
DO ilevel=1,nlevels
  DO ipatch=1,npatches(ilevel)
    !-----------------------------------------------------------------------
    ! Deallocate the field if it is already allocated
    !-----------------------------------------------------------------------
    IF(ASSOCIATED(patch(ilevel,ipatch)%fld)) THEN
      DEALLOCATE(patch(ilevel,ipatch)%fld,stat=info)
      IF (info .NE. 0) THEN
        CALL naga_say(rank,'naga_allocate_patch','Failed to deallocate array.')
        GOTO 9999
      ENDIF
    END IF
    indl = 1 - pset(ilevel,ipatch)%gst(1)
    indl = 1 - pset(ilevel,ipatch)%gst(2)
    indl = 1 - pset(ilevel,ipatch)%gst(3)
    indu = 0
    DO isub=1,topos(ilevel)%nsublist
      isubl=topos(ilevel)%isublist(isub)
      !-----------------------------------------------------------------------
      ! Determine indicies for the field
      !-----------------------------------------------------------------------
      indu(1) = MAX(indu(1),&
              & pset(ilevel,ipatch)%snx(1,isubl) + pset(ilevel,ipatch)%gst(1))
      indu(2) = MAX(indu(2),&
              & pset(ilevel,ipatch)%snx(2,isubl) + pset(ilevel,ipatch)%gst(2))
      indu(3) = MAX(indu(3),&
              & pset(ilevel,ipatch)%snx(3,isubl) + pset(ilevel,ipatch)%gst(3))
    ENDDO !isub
    !-----------------------------------------------------------------------
    ! Allocate the field
    !-----------------------------------------------------------------------
    ALLOCATE(patch(ilevel,ipatch)%fld(&
    & indl(1):indu(1),indl(2):indu(2),indl(3):indu(3),&
    & 1:topos(ilevel)%nsublist),stat=info)
    IF (info .NE. 0) THEN
      CALL naga_say(rank,'naga_allocate_patch_fs','Failed to allocate array.')
      GOTO 9999
    ENDIF
  ENDDO !ipatch
ENDDO !ilevel
!-----------------------------------------------------------------------------
! Return
!-----------------------------------------------------------------------------
 9999 CONTINUE
CALL substop('naga_allocate_patch',t0,info)
RETURN
END SUBROUTINE naga_allocate_patch_fsc
SUBROUTINE naga_allocate_patch_fvc(pset,patch,info)
USE naga_mod_globals
USE naga_mod_say
IMPLICIT NONE
!-----------------------------------------------------------------------------
! Arguments
!-----------------------------------------------------------------------------
TYPE(patch_setup),DIMENSION(:,:),POINTER :: pset
TYPE(patch_field_vc),DIMENSION(:,:),POINTER :: patch
INTEGER, INTENT(OUT) :: info
!-----------------------------------------------------------------------------
! Local variables
!-----------------------------------------------------------------------------
REAL(MK) :: t0
INTEGER, DIMENSION(3) :: indl, indu
INTEGER :: ilevel, ipatch, isub, isubl
!-----------------------------------------------------------------------------
! Initialise routine
!-----------------------------------------------------------------------------
CALL substart('naga_allocate_patch',t0,info)
!-----------------------------------------------------------------------------
! Allocate the multi level, multi patch type (not the actual fields)
!-----------------------------------------------------------------------------
ALLOCATE(patch(nlevels,maxpatches))
!-----------------------------------------------------------------------------
! Now loop through all patches at all levels and swing the wand
!-----------------------------------------------------------------------------
DO ilevel=1,nlevels
  DO ipatch=1,npatches(ilevel)
    !-----------------------------------------------------------------------
    ! Deallocate the field if it is already allocated
    !-----------------------------------------------------------------------
    IF(ASSOCIATED(patch(ilevel,ipatch)%fld)) THEN
      DEALLOCATE(patch(ilevel,ipatch)%fld,stat=info)
      IF (info .NE. 0) THEN
        CALL naga_say(rank,'naga_allocate_patch','Failed to deallocate array.')
        GOTO 9999
      ENDIF
    END IF
    indl = 1 - pset(ilevel,ipatch)%gst(1)
    indl = 1 - pset(ilevel,ipatch)%gst(2)
    indl = 1 - pset(ilevel,ipatch)%gst(3)
    indu = 0
    DO isub=1,topos(ilevel)%nsublist
      isubl=topos(ilevel)%isublist(isub)
      !-----------------------------------------------------------------------
      ! Determine indicies for the field
      !-----------------------------------------------------------------------
      indu(1) = MAX(indu(1),&
              & pset(ilevel,ipatch)%snx(1,isubl) + pset(ilevel,ipatch)%gst(1))
      indu(2) = MAX(indu(2),&
              & pset(ilevel,ipatch)%snx(2,isubl) + pset(ilevel,ipatch)%gst(2))
      indu(3) = MAX(indu(3),&
              & pset(ilevel,ipatch)%snx(3,isubl) + pset(ilevel,ipatch)%gst(3))
    ENDDO !isub
    !-----------------------------------------------------------------------
    ! Allocate the field
    !-----------------------------------------------------------------------
    ALLOCATE(patch(ilevel,ipatch)%fld(&
    & ncom,indl(1):indu(1),indl(2):indu(2),indl(3):indu(3),&
    & 1:topos(ilevel)%nsublist),stat=info)
    IF (info .NE. 0) THEN
      CALL naga_say(rank,'naga_allocate_patch_fs','Failed to allocate array.')
      GOTO 9999
    ENDIF
  ENDDO !ipatch
ENDDO !ilevel
!-----------------------------------------------------------------------------
! Return
!-----------------------------------------------------------------------------
 9999 CONTINUE
CALL substop('naga_allocate_patch',t0,info)
RETURN
END SUBROUTINE naga_allocate_patch_fvc
END MODULE naga_mod_allocate_patch
