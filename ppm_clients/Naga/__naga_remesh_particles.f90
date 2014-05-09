!------------------------------------------------------------------------------
! Subroutine : naga_remesh_particles.f90
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
! This module contains routines for remeshing/creating particles
! Note: concentration/scalar transport seems not to be included here
!------------------------------------------------------------------------------
MODULE naga_mod_remesh_particles
INTERFACE naga_remesh_particles
  MODULE PROCEDURE naga_remesh_particles_v
  MODULE PROCEDURE naga_remesh_particles_vs
  MODULE PROCEDURE naga_remesh_particles_vv
END INTERFACE
CONTAINS
SUBROUTINE naga_remesh_particles_v(info)
USE naga_mod_globals
USE naga_mod_say
USE ppm_module_rmsh
IMPLICIT NONE
!-----------------------------------------------------------------------------
! Arguments
!-----------------------------------------------------------------------------
INTEGER, INTENT(OUT) :: info
!-----------------------------------------------------------------------------
! Local variables
!-----------------------------------------------------------------------------
REAL(MK) :: t0
INTEGER :: ilevel, ipatch
!-----------------------------------------------------------------------------
! Initialise routine
!-----------------------------------------------------------------------------
CALL substart('naga_remesh_particles',t0,info)
!-----------------------------------------------------------------------------
! Now loop through all patches at all levels and swing the wand
!-----------------------------------------------------------------------------
DO ilevel=1,nlevels
  DO ipatch=1,npatches(ilevel)
      !-----------------------------------------------------------------------
      ! Deallocate the particles if they are already allocated
      !-----------------------------------------------------------------------
      !-----------------------------------------------------------------------
      ! Arrays for first order time integration (and higher order)
      !-----------------------------------------------------------------------
      IF(ASSOCIATED(xp(ilevel,ipatch)%val)) THEN
        DEALLOCATE(xp(ilevel,ipatch)%val,stat=info)
        IF (info .NE. 0) THEN
          CALL naga_say(rank,'naga_remesh_particles','Failed to deallocate existing particles.')
          GOTO 9999
        ENDIF
      END IF
      IF(ASSOCIATED(wp(ilevel,ipatch)%val)) THEN
        DEALLOCATE(wp(ilevel,ipatch)%val,stat=info)
        IF (info .NE. 0) THEN
          CALL naga_say(rank,'naga_remesh_particles','Failed to deallocate existing particles.')
          GOTO 9999
        ENDIF
      END IF
      IF(ASSOCIATED(up(ilevel,ipatch)%val)) THEN
        DEALLOCATE(up(ilevel,ipatch)%val,stat=info)
        IF (info .NE. 0) THEN
          CALL naga_say(rank,'naga_remesh_particles','Failed to deallocate existing particles.')
          GOTO 9999
        ENDIF
      END IF
      IF(ASSOCIATED(dwp(ilevel,ipatch)%val)) THEN
        DEALLOCATE(dwp(ilevel,ipatch)%val,stat=info)
        IF (info .NE. 0) THEN
          CALL naga_say(rank,'naga_remesh_particles','Failed to deallocate existing particles.')
          GOTO 9999
        ENDIF
      END IF
      !-----------------------------------------------------------------------
      ! Arrays for second order time integration
      !-----------------------------------------------------------------------
      IF (timeintscheme .EQ. 2) THEN
        IF(ASSOCIATED(xp0(ilevel,ipatch)%val)) THEN
          DEALLOCATE(xp0(ilevel,ipatch)%val,stat=info)
          IF (info .NE. 0) THEN
            CALL naga_say(rank,'naga_remesh_particles','Failed to deallocate existing particles.')
            GOTO 9999
          ENDIF
        END IF
        IF(ASSOCIATED(wp0(ilevel,ipatch)%val)) THEN
          DEALLOCATE(wp0(ilevel,ipatch)%val,stat=info)
          IF (info .NE. 0) THEN
            CALL naga_say(rank,'naga_remesh_particles','Failed to deallocate existing particles.')
            GOTO 9999
          ENDIF
        END IF
        IF(ASSOCIATED(up0(ilevel,ipatch)%val)) THEN
          DEALLOCATE(up0(ilevel,ipatch)%val,stat=info)
          IF (info .NE. 0) THEN
            CALL naga_say(rank,'naga_remesh_particles','Failed to deallocate existing particles.')
            GOTO 9999
          ENDIF
        END IF
        IF(ASSOCIATED(dwp0(ilevel,ipatch)%val)) THEN
          DEALLOCATE(dwp0(ilevel,ipatch)%val,stat=info)
          IF (info .NE. 0) THEN
            CALL naga_say(rank,'naga_remesh_particles','Failed to deallocate existing particles.')
            GOTO 9999
          ENDIF
        END IF
      ENDIF
      !-----------------------------------------------------------------------
      ! Allocate the field
      !-----------------------------------------------------------------------
      CALL ppm_rmsh_create_part(ptcset(ilevel,ipatch)%topoid,&
                              & ptcset(ilevel,ipatch)%meshid,&
                              & xp(ilevel,ipatch)%val,&
                              & ptcset(ilevel,ipatch)%np,&
                              & wp(ilevel,ipatch)%val,&
                              & ncom,&
                              & wf(ilevel,ipatch)%fld,&
                              & (/-0.01_MK,HUGE(0.0_MK)/),& !to let concentration live
                              & info,&
                              & resetpos=.TRUE.&
                              &)
      IF (info .NE. 0) THEN
        CALL naga_say(rank,'naga_remesh_particles_s','Failed remesh particles.')
        GOTO 9999
      ENDIF
  ENDDO !ipatch
ENDDO !ilevel
!-----------------------------------------------------------------------------
! Allocate the arrays that were not allocated by ppm_rmsh_create_part
!-----------------------------------------------------------------------------
DO ilevel=1,nlevels
  DO ipatch=1,npatches(ilevel)
    ALLOCATE(up (ilevel,ipatch)%val(ncom,ptcset(ilevel,ipatch)%np), &
           & dwp (ilevel,ipatch)%val(ncom,ptcset(ilevel,ipatch)%np))
    IF (timeintscheme .EQ. 2) THEN
      ALLOCATE(xp0 (ilevel,ipatch)%val(ncom,ptcset(ilevel,ipatch)%np), &
             & wp0 (ilevel,ipatch)%val(ncom,ptcset(ilevel,ipatch)%np), &
             & up0 (ilevel,ipatch)%val(ncom,ptcset(ilevel,ipatch)%np), &
             & dwp0(ilevel,ipatch)%val(ncom,ptcset(ilevel,ipatch)%np))
    ENDIF
  ENDDO
ENDDO
!-----------------------------------------------------------------------------
! Return
!-----------------------------------------------------------------------------
 9999 CONTINUE
CALL substop('naga_remesh_particles',t0,info)
RETURN
END SUBROUTINE naga_remesh_particles_v
SUBROUTINE naga_remesh_particles_vs(field2,part2,info)
USE naga_mod_globals
USE naga_mod_say
USE ppm_module_rmsh
IMPLICIT NONE
!-----------------------------------------------------------------------------
! Arguments
!-----------------------------------------------------------------------------
TYPE(patch_field_s),DIMENSION(:,:),POINTER :: field2
TYPE(patch_part_s),DIMENSION(:,:),POINTER :: part2
INTEGER, INTENT(OUT) :: info
!-----------------------------------------------------------------------------
! Local variables
!-----------------------------------------------------------------------------
REAL(MK) :: t0
INTEGER :: ilevel, ipatch
!-----------------------------------------------------------------------------
! Initialise routine
!-----------------------------------------------------------------------------
CALL substart('naga_remesh_particles',t0,info)
!-----------------------------------------------------------------------------
! Now loop through all patches at all levels and swing the wand
!-----------------------------------------------------------------------------
DO ilevel=1,nlevels
  DO ipatch=1,npatches(ilevel)
      !-----------------------------------------------------------------------
      ! Deallocate the particles if they are already allocated
      !-----------------------------------------------------------------------
      !-----------------------------------------------------------------------
      ! Arrays for first order time integration (and higher order)
      !-----------------------------------------------------------------------
      IF(ASSOCIATED(xp(ilevel,ipatch)%val)) THEN
        DEALLOCATE(xp(ilevel,ipatch)%val,stat=info)
        IF (info .NE. 0) THEN
          CALL naga_say(rank,'naga_remesh_particles','Failed to deallocate existing particles.')
          GOTO 9999
        ENDIF
      END IF
      IF(ASSOCIATED(wp(ilevel,ipatch)%val)) THEN
        DEALLOCATE(wp(ilevel,ipatch)%val,stat=info)
        IF (info .NE. 0) THEN
          CALL naga_say(rank,'naga_remesh_particles','Failed to deallocate existing particles.')
          GOTO 9999
        ENDIF
      END IF
      IF(ASSOCIATED(up(ilevel,ipatch)%val)) THEN
        DEALLOCATE(up(ilevel,ipatch)%val,stat=info)
        IF (info .NE. 0) THEN
          CALL naga_say(rank,'naga_remesh_particles','Failed to deallocate existing particles.')
          GOTO 9999
        ENDIF
      END IF
      IF(ASSOCIATED(dwp(ilevel,ipatch)%val)) THEN
        DEALLOCATE(dwp(ilevel,ipatch)%val,stat=info)
        IF (info .NE. 0) THEN
          CALL naga_say(rank,'naga_remesh_particles','Failed to deallocate existing particles.')
          GOTO 9999
        ENDIF
      END IF
      IF(ASSOCIATED(part2(ilevel,ipatch)%val)) THEN
        DEALLOCATE(part2(ilevel,ipatch)%val,stat=info)
        IF (info .NE. 0) THEN
          CALL naga_say(rank,'naga_remesh_particles','Failed to deallocate existing particles.')
          GOTO 9999
        ENDIF
      END IF
      !-----------------------------------------------------------------------
      ! Arrays for second order time integration
      !-----------------------------------------------------------------------
      IF (timeintscheme .EQ. 2) THEN
        IF(ASSOCIATED(xp0(ilevel,ipatch)%val)) THEN
          DEALLOCATE(xp0(ilevel,ipatch)%val,stat=info)
          IF (info .NE. 0) THEN
            CALL naga_say(rank,'naga_remesh_particles','Failed to deallocate existing particles.')
            GOTO 9999
          ENDIF
        END IF
        IF(ASSOCIATED(wp0(ilevel,ipatch)%val)) THEN
          DEALLOCATE(wp0(ilevel,ipatch)%val,stat=info)
          IF (info .NE. 0) THEN
            CALL naga_say(rank,'naga_remesh_particles','Failed to deallocate existing particles.')
            GOTO 9999
          ENDIF
        END IF
        IF(ASSOCIATED(up0(ilevel,ipatch)%val)) THEN
          DEALLOCATE(up0(ilevel,ipatch)%val,stat=info)
          IF (info .NE. 0) THEN
            CALL naga_say(rank,'naga_remesh_particles','Failed to deallocate existing particles.')
            GOTO 9999
          ENDIF
        END IF
        IF(ASSOCIATED(dwp0(ilevel,ipatch)%val)) THEN
          DEALLOCATE(dwp0(ilevel,ipatch)%val,stat=info)
          IF (info .NE. 0) THEN
            CALL naga_say(rank,'naga_remesh_particles','Failed to deallocate existing particles.')
            GOTO 9999
          ENDIF
        END IF
      ENDIF
      !-----------------------------------------------------------------------
      ! Allocate the field
      !-----------------------------------------------------------------------
      CALL ppm_rmsh_create_part(ptcset(ilevel,ipatch)%topoid,&
                              & ptcset(ilevel,ipatch)%meshid,&
                              & xp(ilevel,ipatch)%val,&
                              & ptcset(ilevel,ipatch)%np,&
                              & wp(ilevel,ipatch)%val,&
                              & ncom,&
                              & wf(ilevel,ipatch)%fld,&
                              & (/-0.01_MK,HUGE(0.0_MK)/),& !to let concentration live
                              & info,&
                              & resetpos=.TRUE.&
                              & ,field_wp=field2(ilevel,ipatch)%fld,&
                              & wp =part2(ilevel,ipatch)%val &
                              &)
      IF (info .NE. 0) THEN
        CALL naga_say(rank,'naga_remesh_particles_s','Failed remesh particles.')
        GOTO 9999
      ENDIF
  ENDDO !ipatch
ENDDO !ilevel
!-----------------------------------------------------------------------------
! Allocate the arrays that were not allocated by ppm_rmsh_create_part
!-----------------------------------------------------------------------------
DO ilevel=1,nlevels
  DO ipatch=1,npatches(ilevel)
    ALLOCATE(up (ilevel,ipatch)%val(ncom,ptcset(ilevel,ipatch)%np), &
           & dwp (ilevel,ipatch)%val(ncom,ptcset(ilevel,ipatch)%np))
    IF (timeintscheme .EQ. 2) THEN
      ALLOCATE(xp0 (ilevel,ipatch)%val(ncom,ptcset(ilevel,ipatch)%np), &
             & wp0 (ilevel,ipatch)%val(ncom,ptcset(ilevel,ipatch)%np), &
             & up0 (ilevel,ipatch)%val(ncom,ptcset(ilevel,ipatch)%np), &
             & dwp0(ilevel,ipatch)%val(ncom,ptcset(ilevel,ipatch)%np))
    ENDIF
  ENDDO
ENDDO
!-----------------------------------------------------------------------------
! Return
!-----------------------------------------------------------------------------
 9999 CONTINUE
CALL substop('naga_remesh_particles',t0,info)
RETURN
END SUBROUTINE naga_remesh_particles_vs
SUBROUTINE naga_remesh_particles_vv(field2,part2,info)
USE naga_mod_globals
USE naga_mod_say
USE ppm_module_rmsh
IMPLICIT NONE
!-----------------------------------------------------------------------------
! Arguments
!-----------------------------------------------------------------------------
TYPE(patch_field_v),DIMENSION(:,:),POINTER :: field2
TYPE(patch_part_v),DIMENSION(:,:),POINTER :: part2
INTEGER, INTENT(OUT) :: info
!-----------------------------------------------------------------------------
! Local variables
!-----------------------------------------------------------------------------
REAL(MK) :: t0
INTEGER :: ilevel, ipatch
!-----------------------------------------------------------------------------
! Initialise routine
!-----------------------------------------------------------------------------
CALL substart('naga_remesh_particles',t0,info)
!-----------------------------------------------------------------------------
! Now loop through all patches at all levels and swing the wand
!-----------------------------------------------------------------------------
DO ilevel=1,nlevels
  DO ipatch=1,npatches(ilevel)
      !-----------------------------------------------------------------------
      ! Deallocate the particles if they are already allocated
      !-----------------------------------------------------------------------
      !-----------------------------------------------------------------------
      ! Arrays for first order time integration (and higher order)
      !-----------------------------------------------------------------------
      IF(ASSOCIATED(xp(ilevel,ipatch)%val)) THEN
        DEALLOCATE(xp(ilevel,ipatch)%val,stat=info)
        IF (info .NE. 0) THEN
          CALL naga_say(rank,'naga_remesh_particles','Failed to deallocate existing particles.')
          GOTO 9999
        ENDIF
      END IF
      IF(ASSOCIATED(wp(ilevel,ipatch)%val)) THEN
        DEALLOCATE(wp(ilevel,ipatch)%val,stat=info)
        IF (info .NE. 0) THEN
          CALL naga_say(rank,'naga_remesh_particles','Failed to deallocate existing particles.')
          GOTO 9999
        ENDIF
      END IF
      IF(ASSOCIATED(up(ilevel,ipatch)%val)) THEN
        DEALLOCATE(up(ilevel,ipatch)%val,stat=info)
        IF (info .NE. 0) THEN
          CALL naga_say(rank,'naga_remesh_particles','Failed to deallocate existing particles.')
          GOTO 9999
        ENDIF
      END IF
      IF(ASSOCIATED(dwp(ilevel,ipatch)%val)) THEN
        DEALLOCATE(dwp(ilevel,ipatch)%val,stat=info)
        IF (info .NE. 0) THEN
          CALL naga_say(rank,'naga_remesh_particles','Failed to deallocate existing particles.')
          GOTO 9999
        ENDIF
      END IF
      IF(ASSOCIATED(part2(ilevel,ipatch)%val)) THEN
        DEALLOCATE(part2(ilevel,ipatch)%val,stat=info)
        IF (info .NE. 0) THEN
          CALL naga_say(rank,'naga_remesh_particles','Failed to deallocate existing particles.')
          GOTO 9999
        ENDIF
      END IF
      !-----------------------------------------------------------------------
      ! Arrays for second order time integration
      !-----------------------------------------------------------------------
      IF (timeintscheme .EQ. 2) THEN
        IF(ASSOCIATED(xp0(ilevel,ipatch)%val)) THEN
          DEALLOCATE(xp0(ilevel,ipatch)%val,stat=info)
          IF (info .NE. 0) THEN
            CALL naga_say(rank,'naga_remesh_particles','Failed to deallocate existing particles.')
            GOTO 9999
          ENDIF
        END IF
        IF(ASSOCIATED(wp0(ilevel,ipatch)%val)) THEN
          DEALLOCATE(wp0(ilevel,ipatch)%val,stat=info)
          IF (info .NE. 0) THEN
            CALL naga_say(rank,'naga_remesh_particles','Failed to deallocate existing particles.')
            GOTO 9999
          ENDIF
        END IF
        IF(ASSOCIATED(up0(ilevel,ipatch)%val)) THEN
          DEALLOCATE(up0(ilevel,ipatch)%val,stat=info)
          IF (info .NE. 0) THEN
            CALL naga_say(rank,'naga_remesh_particles','Failed to deallocate existing particles.')
            GOTO 9999
          ENDIF
        END IF
        IF(ASSOCIATED(dwp0(ilevel,ipatch)%val)) THEN
          DEALLOCATE(dwp0(ilevel,ipatch)%val,stat=info)
          IF (info .NE. 0) THEN
            CALL naga_say(rank,'naga_remesh_particles','Failed to deallocate existing particles.')
            GOTO 9999
          ENDIF
        END IF
      ENDIF
      !-----------------------------------------------------------------------
      ! Allocate the field
      !-----------------------------------------------------------------------
      CALL ppm_rmsh_create_part(ptcset(ilevel,ipatch)%topoid,&
                              & ptcset(ilevel,ipatch)%meshid,&
                              & xp(ilevel,ipatch)%val,&
                              & ptcset(ilevel,ipatch)%np,&
                              & wp(ilevel,ipatch)%val,&
                              & ncom,&
                              & wf(ilevel,ipatch)%fld,&
                              & (/-0.01_MK,HUGE(0.0_MK)/),& !to let concentration live
                              & info,&
                              & resetpos=.TRUE.&
                              & ,field_wp=field2(ilevel,ipatch)%fld,&
                              & wp =part2(ilevel,ipatch)%val &
                              & ,lda2=ncom&
                              &)
      IF (info .NE. 0) THEN
        CALL naga_say(rank,'naga_remesh_particles_s','Failed remesh particles.')
        GOTO 9999
      ENDIF
  ENDDO !ipatch
ENDDO !ilevel
!-----------------------------------------------------------------------------
! Allocate the arrays that were not allocated by ppm_rmsh_create_part
!-----------------------------------------------------------------------------
DO ilevel=1,nlevels
  DO ipatch=1,npatches(ilevel)
    ALLOCATE(up (ilevel,ipatch)%val(ncom,ptcset(ilevel,ipatch)%np), &
           & dwp (ilevel,ipatch)%val(ncom,ptcset(ilevel,ipatch)%np))
    IF (timeintscheme .EQ. 2) THEN
      ALLOCATE(xp0 (ilevel,ipatch)%val(ncom,ptcset(ilevel,ipatch)%np), &
             & wp0 (ilevel,ipatch)%val(ncom,ptcset(ilevel,ipatch)%np), &
             & up0 (ilevel,ipatch)%val(ncom,ptcset(ilevel,ipatch)%np), &
             & dwp0(ilevel,ipatch)%val(ncom,ptcset(ilevel,ipatch)%np))
    ENDIF
  ENDDO
ENDDO
!-----------------------------------------------------------------------------
! Return
!-----------------------------------------------------------------------------
 9999 CONTINUE
CALL substop('naga_remesh_particles',t0,info)
RETURN
END SUBROUTINE naga_remesh_particles_vv
END MODULE naga_mod_remesh_particles
