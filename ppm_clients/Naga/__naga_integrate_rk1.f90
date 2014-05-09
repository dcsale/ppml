!------------------------------------------------------------------------------
! Subroutine : naga_integrate_rk1.f90
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
! This routines integrates particle positions and strengths in time using
! Euler / first order Runge-Kutta integration
!------------------------------------------------------------------------------
MODULE naga_mod_integrate_rk1
IMPLICIT NONE
INTERFACE naga_integrate_rk1
  MODULE PROCEDURE naga_integrate_rk1
END INTERFACE
CONTAINS
SUBROUTINE naga_integrate_rk1(info)
USE naga_mod_globals
USE naga_mod_say
USE naga_mod_rhs
USE naga_mod_map_particles
IMPLICIT NONE
!-------------------------------------------------------------------------
! Arguments
!-------------------------------------------------------------------------
INTEGER, INTENT(INOUT) :: info
!-------------------------------------------------------------------------
! Local variables
!-------------------------------------------------------------------------
REAL(MK) :: t0
INTEGER :: i,j
INTEGER :: ilevel,ipatch
INTEGER :: istage
!----------------------------------------------------------------------------
! Initialise routine
!----------------------------------------------------------------------------
CALL substart('naga_integrate_rk1',t0,info)
info = 0
!----------------------------------------------------------------------------
! Start time integration loop
!----------------------------------------------------------------------------
DO WHILE(time .LT. endtime)
  !----------------------------------------------------------------------------
  ! This could be the place for writing restart files
  !----------------------------------------------------------------------------
  IF (rank .EQ. 0) THEN
    write(*,*) 'itime',itime
  ENDIF
  !----------------------------------------------------------------------------
  ! First (and only) sub step
  !----------------------------------------------------------------------------
  istage = 1
  CALL naga_rhs(istage,info)
  IF (info .NE. 0) THEN
    CALL naga_say(rank, 'naga_integrate_rk1','Failed to compute RHS.')
    GOTO 9999
  END IF
  DO ilevel=1,nlevels
    DO ipatch=1,npatches(ilevel)
      DO j=1,ptcset(ilevel,ipatch)%np
        DO i=1,ncom
          wp(ilevel,ipatch)%val(i,j) = wp(ilevel,ipatch)%val(i,j) &
              + dtime*dwp(ilevel,ipatch)%val(i,j)
          xp(ilevel,ipatch)%val(i,j) = xp(ilevel,ipatch)%val(i,j) &
              + dtime*up(ilevel,ipatch)%val(i,j)
        ENDDO
      ENDDO
    ENDDO
  ENDDO
  !----------------------------------------------------------------------------
  ! Map and catch particles
  !----------------------------------------------------------------------------
  CALL naga_map_particles(istage,info)
  IF (info .NE. 0) THEN
    CALL naga_say(rank, 'naga_integrate_rk1','Failed to map particles.')
    GOTO 9999
  END IF
  !----------------------------------------------------------------------------
  ! Increment values
  !----------------------------------------------------------------------------
  itime = itime + 1
  time = time + dtime
  !----------------------------------------------------------------------------
  ! Stop criteria
  !----------------------------------------------------------------------------
  IF (maxitime .NE. 0) THEN
    IF (itime .GE. maxitime) THEN
      IF (rank .EQ. 0) THEN
        CALL naga_say(rank, 'naga_integrate_rk1','Maximum time step reached, exiting.')
      ENDIF
      GOTO 9999
    ENDIF
  ENDIF
  IF (forceabort) THEN
    IF (rank .EQ. 0) THEN
      CALL naga_say(rank, 'naga_integrate_rk1','ABORT file found, exiting.')
    ENDIF
    GOTO 9999
  ENDIF
ENDDO
IF (rank .EQ. 0) THEN
  CALL naga_say(rank, 'naga_integrate_rk1','Time integration complete, exiting.')
ENDIF
!----------------------------------------------------------------------------
! Return
!----------------------------------------------------------------------------
 9999 CONTINUE
CALL substop('naga_integrate_rk1',t0,info)
RETURN
END SUBROUTINE naga_integrate_rk1
END MODULE naga_mod_integrate_rk1
