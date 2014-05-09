!------------------------------------------------------------------------------
! Subroutine : naga_catch_particles.f90
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
! This routines catches particles and reinserts them periodically into the
! computational domain
! Particles are being caught twice as truncation errors may reinsert them
! outside the domain at the 16th digit. Sounds unlike but with lots of
! particles it is unfortunately probable.
!------------------------------------------------------------------------------
MODULE naga_mod_catch_particles
IMPLICIT NONE
INTERFACE naga_catch_particles
  MODULE PROCEDURE naga_catch_particles
END INTERFACE
CONTAINS
SUBROUTINE naga_catch_particles(partx,info)
USE naga_mod_globals
USE naga_mod_say
IMPLICIT NONE
!-------------------------------------------------------------------------
! Arguments
!-------------------------------------------------------------------------
INTEGER, INTENT(INOUT) :: info
TYPE(patch_part_v),DIMENSION(:,:),POINTER :: partx
!-------------------------------------------------------------------------
! Local variables
!-------------------------------------------------------------------------
REAL(MK) :: t0
REAL(MK) :: gmin,gmax,glen
INTEGER :: j
INTEGER :: ilevel,ipatch
INTEGER :: repetitions
!----------------------------------------------------------------------------
! Initialise routine
!----------------------------------------------------------------------------
CALL substart('naga_catch_particles',t0,info)
info = 0
DO ilevel=1,nlevels
  DO ipatch=1,npatches(ilevel)
    IF (ilevel .NE. 1 .OR. &
      & ilevel .EQ. 1 .AND. domainbc .EQ. 1) THEN
      gmin = ptcset(ilevel,ipatch)%min(1)
      gmax = ptcset(ilevel,ipatch)%max(1)
      glen = gmax-gmin
      DO j=1,ptcset(ilevel,ipatch)%np
        DO repetitions=1,2
          IF (partx(ilevel,ipatch)%val(1,j) .GE. gmax) THEN
            partx(ilevel,ipatch)%val(1,j) = partx(ilevel,ipatch)%val(1,j) - glen
          ELSE IF (partx(ilevel,ipatch)%val(1,j) .LT. gmin) THEN
            partx(ilevel,ipatch)%val(1,j) = partx(ilevel,ipatch)%val(1,j) + glen
          ENDIF
        ENDDO !repetitions
      ENDDO
      gmin = ptcset(ilevel,ipatch)%min(2)
      gmax = ptcset(ilevel,ipatch)%max(2)
      glen = gmax-gmin
      DO j=1,ptcset(ilevel,ipatch)%np
        DO repetitions=1,2
          IF (partx(ilevel,ipatch)%val(2,j) .GE. gmax) THEN
            partx(ilevel,ipatch)%val(2,j) = partx(ilevel,ipatch)%val(2,j) - glen
          ELSE IF (partx(ilevel,ipatch)%val(2,j) .LT. gmin) THEN
            partx(ilevel,ipatch)%val(2,j) = partx(ilevel,ipatch)%val(2,j) + glen
          ENDIF
        ENDDO !repetitions
      ENDDO
      gmin = ptcset(ilevel,ipatch)%min(3)
      gmax = ptcset(ilevel,ipatch)%max(3)
      glen = gmax-gmin
      DO j=1,ptcset(ilevel,ipatch)%np
        DO repetitions=1,2
          IF (partx(ilevel,ipatch)%val(3,j) .GE. gmax) THEN
            partx(ilevel,ipatch)%val(3,j) = partx(ilevel,ipatch)%val(3,j) - glen
          ELSE IF (partx(ilevel,ipatch)%val(3,j) .LT. gmin) THEN
            partx(ilevel,ipatch)%val(3,j) = partx(ilevel,ipatch)%val(3,j) + glen
          ENDIF
        ENDDO !repetitions
      ENDDO
    ENDIF
  ENDDO !ipatch
ENDDO !ilevel
!----------------------------------------------------------------------------
! Return
!----------------------------------------------------------------------------
 9999 CONTINUE
CALL substop('naga_catch_particles',t0,info)
RETURN
END SUBROUTINE naga_catch_particles
END MODULE naga_mod_catch_particles
