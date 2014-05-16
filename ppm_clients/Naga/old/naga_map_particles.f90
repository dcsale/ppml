!------------------------------------------------------------------------------
! Subroutine :  naga_map_particles.f90
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
! This routine maps particles between sub domains
!------------------------------------------------------------------------------

MODULE naga_mod_map_particles

INTERFACE naga_map_particles
  MODULE PROCEDURE naga_map_particles
END INTERFACE

CONTAINS

SUBROUTINE naga_map_particles(istage,info)

USE naga_mod_globals
USE naga_mod_say
USE naga_mod_catch_particles
USE ppm_module_map_part_partial
USE ppm_module_map_part

IMPLICIT NONE

!-----------------------------------------------------------------------------
! Arguments
!-----------------------------------------------------------------------------
!istage is required such that e.g. xp0,wp0,etc are not mapped after final step
INTEGER, INTENT(IN)                               :: istage
INTEGER, INTENT(OUT)                              :: info


!-----------------------------------------------------------------------------
! Local variables
!-----------------------------------------------------------------------------
REAL(MK)              :: t0
INTEGER               :: ilevel, ipatch
INTEGER               :: newnp
INTEGER               :: i !@tmp
!!REAL(MK)              :: maxx,maxy,maxz
LOGICAL               :: ignoreunassigned


!-----------------------------------------------------------------------------
! Initialise routine
!-----------------------------------------------------------------------------
CALL substart('naga_map_particles',t0,info)


!-----------------------------------------------------------------------------
! Catch particles leaving the computational domain
!-----------------------------------------------------------------------------
IF (domainbc .EQ. 1) THEN
  CALL naga_catch_particles(xp,info) !@ THIS MIGHT NOT BE NECESSARY
  ignoreunassigned = .FALSE.
ELSE
  ignoreunassigned = .TRUE.
ENDIF


!-----------------------------------------------------------------------------
! Now loop through all patches at all levels and swing the wand
!-----------------------------------------------------------------------------
DO ilevel=1,nlevels
  DO ipatch=1,npatches(ilevel)
    !---------------------------------------------------------------------------
    ! push positions
    !---------------------------------------------------------------------------
    CALL ppm_map_part_partial(ptcset(ilevel,ipatch)%topoid,&
                            & xp(ilevel,ipatch)%val,&
                            & ptcset(ilevel,ipatch)%np,&
                            !& info)
                            & info,ignoreunassigned)

    !!DO i=1,ptcset(ilevel,ipatch)%np
      !!write(*,*) xp(ilevel,ipatch)%val(:,i)
    !!END DO
    IF (info .NE. 0) THEN
      CALL naga_say(rank,'naga_map_particles','Failed to initialise mapping.')
      GOTO 9999
    ENDIF
    !---------------------------------------------------------------------------
    ! push additional values
    !---------------------------------------------------------------------------
    ! vorticity
    CALL ppm_map_part_push   (wp(ilevel,ipatch)%val,&
                            & ncom,&
                            & ptcset(ilevel,ipatch)%np,&
                            & info)
    IF (info .NE. 0) THEN
      CALL naga_say(rank,'naga_map_particles','Failed to push vorticity particles.')
      GOTO 9999
    ENDIF
    ! When doing higher order integration schemes additional particle states exist
    ! which must be mapped as well
    IF ((timeintscheme .EQ. 2 .AND. istage .EQ. 1)) THEN
      ! position 0
      CALL ppm_map_part_push   (xp0(ilevel,ipatch)%val,&
                              & ncom,&
                              & ptcset(ilevel,ipatch)%np,&
                              & info)
      IF (info .NE. 0) THEN
        CALL naga_say(rank,'naga_map_particles','Failed to push position 0 particles.')
        GOTO 9999
      ENDIF
      ! vorticity 0
      CALL ppm_map_part_push   (wp0(ilevel,ipatch)%val,&
                              & ncom,&
                              & ptcset(ilevel,ipatch)%np,&
                              & info)
      IF (info .NE. 0) THEN
        CALL naga_say(rank,'naga_map_particles','Failed to push vorticity 0 particles.')
        GOTO 9999
      ENDIF
      ! velocity 0
      CALL ppm_map_part_push   (up0(ilevel,ipatch)%val,&
                              & ncom,&
                              & ptcset(ilevel,ipatch)%np,&
                              & info)
      IF (info .NE. 0) THEN
        CALL naga_say(rank,'naga_map_particles','Failed to push velocity 0 particles.')
        GOTO 9999
      ENDIF
      ! vorticity derivatives 0
      CALL ppm_map_part_push   (dwp0(ilevel,ipatch)%val,&
                              & ncom,&
                              & ptcset(ilevel,ipatch)%np,&
                              & info)
      IF (info .NE. 0) THEN
        CALL naga_say(rank,'naga_map_particles','Failed to push vorticity derivatives 0 particles.')
        GOTO 9999
      ENDIF
    ENDIF
    IF (concentration .EQ. 1) THEN
      ! concentration
      CALL ppm_map_part_push   (cp(ilevel,ipatch)%val,&
                              & ptcset(ilevel,ipatch)%np,&
                              & info)
      IF (info .NE. 0) THEN
        CALL naga_say(rank,'naga_map_particles','Failed to push concentration particles.')
        GOTO 9999
      ENDIF
    ENDIF
    !---------------------------------------------------------------------------
    ! send the data
    !---------------------------------------------------------------------------
    CALL ppm_map_part_send   (ptcset(ilevel,ipatch)%np,&
                            & newnp,&
                            & info)
    IF (info .NE. 0) THEN
      CALL naga_say(rank,'naga_map_particles','Failed to send particles.')
      GOTO 9999
    ENDIF
    !---------------------------------------------------------------------------
    ! receive
    !---------------------------------------------------------------------------
    IF (concentration .EQ. 1) THEN
      ! concentration
      CALL ppm_map_part_pop    (cp(ilevel,ipatch)%val,&
                              & ptcset(ilevel,ipatch)%np,&
                              & newnp,&
                              & info)
      IF (info .NE. 0) THEN
        GOTO 9999
        CALL naga_say(rank,'naga_map_particles','Failed to pop concentration particles.')
      ENDIF
    ENDIF
    ! When doing higher order integration schemes additional particle states exist
    ! which must be mapped as well
    IF ((timeintscheme .EQ. 2 .AND. istage .EQ. 1)) THEN
      ! vorticity derivatives 0
      CALL ppm_map_part_pop    (dwp0(ilevel,ipatch)%val,&
                              & ncom,&
                              & ptcset(ilevel,ipatch)%np,&
                              & newnp,&
                              & info)
      IF (info .NE. 0) THEN
        CALL naga_say(rank,'naga_map_particles','Failed to pop vorticity derivatives 0 particles.')
        GOTO 9999
      ENDIF
      ! velocity 0
      CALL ppm_map_part_pop    (up0(ilevel,ipatch)%val,&
                              & ncom,&
                              & ptcset(ilevel,ipatch)%np,&
                              & newnp,&
                              & info)
      IF (info .NE. 0) THEN
        CALL naga_say(rank,'naga_map_particles','Failed to pop velocity 0 particles.')
        GOTO 9999
      ENDIF
      ! vorticity 0
      CALL ppm_map_part_pop    (wp0(ilevel,ipatch)%val,&
                              & ncom,&
                              & ptcset(ilevel,ipatch)%np,&
                              & newnp,&
                              & info)
      IF (info .NE. 0) THEN
        CALL naga_say(rank,'naga_map_particles','Failed to pop vorticity 0 particles.')
        GOTO 9999
      ENDIF
      ! position 0
      CALL ppm_map_part_pop    (xp0(ilevel,ipatch)%val,&
                              & ncom,&
                              & ptcset(ilevel,ipatch)%np,&
                              & newnp,&
                              & info)
      IF (info .NE. 0) THEN
        CALL naga_say(rank,'naga_map_particles','Failed to pop position 0 particles.')
        GOTO 9999
      ENDIF
    ENDIF
    ! These guys are always there:
    ! vorticity
    CALL ppm_map_part_pop    (wp(ilevel,ipatch)%val,&
                            & ncom,&
                            & ptcset(ilevel,ipatch)%np,&
                            & newnp,&
                            & info)
    IF (info .NE. 0) THEN
      CALL naga_say(rank,'naga_map_particles','Failed to pop vorticity particles.')
      GOTO 9999
    ENDIF
    !---------------------------------------------------------------------------
    ! receive positions
    !---------------------------------------------------------------------------
    CALL ppm_map_part_pop    (xp(ilevel,ipatch)%val,&
                            & ndim,&
                            & ptcset(ilevel,ipatch)%np,&
                            & newnp,&
                            & info)
    ptcset(ilevel,ipatch)%np = newnp
    IF (info .NE. 0) THEN
      CALL naga_say(rank,'naga_map_particles','Failed to particle positions.')
      GOTO 9999
    ENDIF

  ENDDO !ipatch
ENDDO !ilevel


!-----------------------------------------------------------------------------
! Reallocate particles arrays (u and dw)
! This may not be the correct way to handle u and dw. Maybe they should be
! mapped later on.
! If particles are not remeshed all xp0,wp0,etc must also be resized to np
!-----------------------------------------------------------------------------
!@ This is a place to note - it might cause problems in the rk2+ integrations
DO ilevel=1,nlevels
  DO ipatch=1,npatches(ilevel)
    DEALLOCATE(up (ilevel,ipatch)%val)
    DEALLOCATE(dwp(ilevel,ipatch)%val)
    ALLOCATE(up (ilevel,ipatch)%val(ncom,ptcset(ilevel,ipatch)%np), &
             dwp(ilevel,ipatch)%val(ncom,ptcset(ilevel,ipatch)%np))
  ENDDO
ENDDO
IF (istage .NE. 1 .AND. iremesh .EQ. 0) THEN
  DO ilevel=1,nlevels
    DO ipatch=1,npatches(ilevel)
      DEALLOCATE(up0 (ilevel,ipatch)%val)
      DEALLOCATE(dwp0(ilevel,ipatch)%val)
      DEALLOCATE(xp0 (ilevel,ipatch)%val)
      DEALLOCATE(wp0 (ilevel,ipatch)%val)
      ALLOCATE(xp0 (ilevel,ipatch)%val(ncom,ptcset(ilevel,ipatch)%np), &
             & wp0 (ilevel,ipatch)%val(ncom,ptcset(ilevel,ipatch)%np), &
             & up0 (ilevel,ipatch)%val(ncom,ptcset(ilevel,ipatch)%np), &
             & dwp0(ilevel,ipatch)%val(ncom,ptcset(ilevel,ipatch)%np))
    ENDDO
  ENDDO
ENDIF

!-----------------------------------------------------------------------------
! Return
!-----------------------------------------------------------------------------
 9999 CONTINUE
CALL substop('naga_map_particles',t0,info)
RETURN


END SUBROUTINE naga_map_particles

END MODULE naga_mod_map_particles
