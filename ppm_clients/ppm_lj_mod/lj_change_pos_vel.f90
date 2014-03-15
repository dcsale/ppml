SUBROUTINE lj_change_pos_vel(info)
!======================================================================!
! integrates in time yielding new positions and velocities
!======================================================================!
  USE ppm_module_core
  USE lj_module_global
  IMPLICIT NONE
  INCLUDE 'mpif.h'
  ! arguments
  INTEGER,INTENT(OUT)                   :: info
  ! local variables
  INTEGER                               :: ip

  !====================================================================!
  ! initialize
  info = 0 ! change if error occurs

  !====================================================================!
  ! time integration (forward Euler)

  IF (ndim .EQ. 2) THEN
     DO ip=1,Npart

        accp(1,ip) = SIGN(MIN(ABS(accp(1,ip)),amax),accp(1,ip))
        accp(2,ip) = SIGN(MIN(ABS(accp(2,ip)),amax),accp(2,ip))
        
        ! new positions:
        xp(1,ip) = xp(1,ip) + vp(1,ip)*dt + 0.5_prec*accp(1,ip)*dt**2
        xp(2,ip) = xp(2,ip) + vp(2,ip)*dt + 0.5_prec*accp(2,ip)*dt**2
        ! impose periodic boundaries:
        IF (xp(1,ip) .GT. maxdom(1)) xp(1,ip) = xp(1,ip) - lendom(1)
        IF (xp(2,ip) .GT. maxdom(2)) xp(2,ip) = xp(2,ip) - lendom(2)
        IF (xp(1,ip) .LT. mindom(1)) xp(1,ip) = xp(1,ip) + lendom(1)
        IF (xp(2,ip) .LT. mindom(2)) xp(2,ip) = xp(2,ip) + lendom(2)
        ! checker:
        IF (xp(1,ip) .GT. maxdom(1)) WRITE(*,*)  &
             'ERROR: Time step too big, ip,xp(1,ip)=',ip,xp(1,ip)
        IF (xp(2,ip) .GT. maxdom(2)) WRITE(*,*)  &
             'ERROR: Time step too big, ip,xp(2,ip)=',ip,xp(2,ip)
        IF (xp(1,ip) .LT. mindom(1)) WRITE(*,*)  &
             'ERROR: Time step too big, ip,xp(1,ip)=',ip,xp(1,ip)
        IF (xp(2,ip) .LT. mindom(2)) WRITE(*,*)  &
             'ERROR: Time step too big, ip,xp(2,ip)=',ip,xp(2,ip)
        ! new velocities:
        vp(1,ip) = vp(1,ip) + accp(1,ip)*dt
        vp(2,ip) = vp(2,ip) + accp(2,ip)*dt
        
        vp(1,ip) = SIGN(MIN(ABS(vp(1,ip)),vmax),vp(1,ip))
        vp(2,ip) = SIGN(MIN(ABS(vp(2,ip)),vmax),vp(2,ip))
        
     ENDDO
  ELSE
     DO ip=1,Npart

        accp(1,ip) = SIGN(MIN(ABS(accp(1,ip)),amax),accp(1,ip))
        accp(2,ip) = SIGN(MIN(ABS(accp(2,ip)),amax),accp(2,ip))
        accp(3,ip) = SIGN(MIN(ABS(accp(3,ip)),amax),accp(3,ip))
        
        ! new positions:
        xp(1,ip) = xp(1,ip) + vp(1,ip)*dt + 0.5_prec*accp(1,ip)*dt**2
        xp(2,ip) = xp(2,ip) + vp(2,ip)*dt + 0.5_prec*accp(2,ip)*dt**2
        xp(3,ip) = xp(3,ip) + vp(3,ip)*dt + 0.5_prec*accp(3,ip)*dt**2
        ! impose periodic boundaries:
        IF (xp(1,ip) .GT. maxdom(1)) xp(1,ip) = xp(1,ip) - lendom(1)
        IF (xp(2,ip) .GT. maxdom(2)) xp(2,ip) = xp(2,ip) - lendom(2)
        IF (xp(3,ip) .GT. maxdom(3)) xp(3,ip) = xp(3,ip) - lendom(3)
        IF (xp(1,ip) .LT. mindom(1)) xp(1,ip) = xp(1,ip) + lendom(1)
        IF (xp(2,ip) .LT. mindom(2)) xp(2,ip) = xp(2,ip) + lendom(2)
        IF (xp(3,ip) .LT. mindom(3)) xp(3,ip) = xp(3,ip) + lendom(3)
        ! checker:
        IF (xp(1,ip) .GT. maxdom(1)) WRITE(*,*)  &
             'ERROR: Time step too big, ip,xp(1,ip)=',ip,xp(1,ip)
        IF (xp(2,ip) .GT. maxdom(2)) WRITE(*,*)  &
             'ERROR: Time step too big, ip,xp(2,ip)=',ip,xp(2,ip)
        IF (xp(3,ip) .GT. maxdom(3)) WRITE(*,*)  &
             'ERROR: Time step too big, ip,xp(3,ip)=',ip,xp(3,ip)
        IF (xp(1,ip) .LT. mindom(1)) WRITE(*,*)  &
             'ERROR: Time step too big, ip,xp(1,ip)=',ip,xp(1,ip)
        IF (xp(2,ip) .LT. mindom(2)) WRITE(*,*)  &
             'ERROR: Time step too big, ip,xp(2,ip)=',ip,xp(2,ip)
        IF (xp(3,ip) .LT. mindom(3)) WRITE(*,*)  &
             'ERROR: Time step too big, ip,xp(3,ip)=',ip,xp(3,ip)
        ! new velocities:
        vp(1,ip) = vp(1,ip) + accp(1,ip)*dt
        vp(2,ip) = vp(2,ip) + accp(2,ip)*dt
        vp(3,ip) = vp(3,ip) + accp(3,ip)*dt
        
        vp(1,ip) = SIGN(MIN(ABS(vp(1,ip)),vmax),vp(1,ip))
        vp(2,ip) = SIGN(MIN(ABS(vp(2,ip)),vmax),vp(2,ip))
        vp(3,ip) = SIGN(MIN(ABS(vp(3,ip)),vmax),vp(3,ip))
        
     ENDDO
  ENDIF

9999  CONTINUE ! jump here upon error

END SUBROUTINE lj_change_pos_vel
