SUBROUTINE lj_set_particles(info)
!======================================================================!
! set initial particle positions and velocities
!======================================================================!
  USE lj_module_global
  IMPLICIT NONE
  ! arguments
  INTEGER,INTENT(OUT) :: info
  ! local variables
  INTEGER :: ip,seedsize,i
  INTEGER,DIMENSION(:),ALLOCATABLE :: seed
  REAL(prec) :: randvar,twopi,pi
  REAL(prec),DIMENSION(:),ALLOCATABLE :: x1,x2,x3,x4,x5,x6
  !====================================================================!
  ! initialization
  info = 0
  pi = ACOS(-1._prec)
  twopi = 2._prec*pi
  Npart = Npart_global/nproc
  IF (rank .EQ. 0) Npart = Npart + MOD(Npart_global,nproc)
  !====================================================================!
  ! allocation
  ALLOCATE(xp(ndim,Npart),STAT=info)
  IF (info .NE. 0) THEN
     WRITE(*,*)'ERROR in lj_set_particles: allocation failed.'
     GOTO 9999
  ENDIF
  ALLOCATE(vp(ndim,Npart),STAT=info)
  IF (info .NE. 0) THEN
     WRITE(*,*)'ERROR in lj_set_particles: allocation failed.'
     GOTO 9999
  ENDIF
  ALLOCATE(x1(Npart),STAT=info)
  IF (info .NE. 0) THEN
     WRITE(*,*)'ERROR in lj_set_particles: allocation failed.'
     GOTO 9999
  ENDIF
  ALLOCATE(x2(Npart),STAT=info)
  IF (info .NE. 0) THEN
     WRITE(*,*)'ERROR in lj_set_particles: allocation failed.'
     GOTO 9999
  ENDIF
  ALLOCATE(x3(Npart),STAT=info)
  IF (info .NE. 0) THEN
     WRITE(*,*)'ERROR in lj_set_particles: allocation failed.'
     GOTO 9999
  ENDIF
  ALLOCATE(x4(Npart),STAT=info)
  IF (info .NE. 0) THEN
     WRITE(*,*)'ERROR in lj_set_particles: allocation failed.'
     GOTO 9999
  ENDIF
  IF (ndim .EQ. 3) THEN
     ALLOCATE(x5(Npart),STAT=info)
     IF (info .NE. 0) THEN
        WRITE(*,*)'ERROR in lj_set_particles: allocation failed.'
        GOTO 9999
     ENDIF
     ALLOCATE(x6(Npart),STAT=info)
     IF (info .NE. 0) THEN
        WRITE(*,*)'ERROR in lj_set_particles: allocation failed.'
        GOTO 9999
     ENDIF
  ENDIF
  !====================================================================!
  ! set seed for random number generator
  CALL RANDOM_SEED(SIZE=seedsize)
  ALLOCATE(seed(seedsize))
  DO i=1,seedsize
     seed(i) = 12345*(11356*rank + 37)
  ENDDO
  CALL RANDOM_SEED(PUT=seed)
  DEALLOCATE(seed)
  !====================================================================!
  ! distribute particles randomly in computational domain
  CALL RANDOM_NUMBER(xp)
  IF (ndim .EQ. 2) THEN
     DO ip=1,Npart
        xp(1,ip) = xp(1,ip)*lendom(1) + mindom(1)
        xp(2,ip) = xp(2,ip)*lendom(2) + mindom(2)
     ENDDO
  ELSE
     DO ip=1,Npart
        xp(1,ip) = xp(1,ip)*lendom(1) + mindom(1)
        xp(2,ip) = xp(2,ip)*lendom(2) + mindom(2)
        xp(3,ip) = xp(3,ip)*lendom(3) + mindom(3)
     ENDDO
  ENDIF
  !====================================================================!
  ! choose random velocity components, normally distributed with
  ! variance sigma2
  CALL RANDOM_NUMBER(x1)
  CALL RANDOM_NUMBER(x2)
  CALL RANDOM_NUMBER(x3)
  CALL RANDOM_NUMBER(x4)
  IF (ndim .EQ. 2) THEN
     sigma2 = 3._prec*k_boltzmann*Temp/(2._prec*mass) ! in (Angstrom/s)^2
     DO ip=1,Npart
        vp(1,ip) = SQRT(-2._prec*sigma2*LOG(x1(ip)))*COS(twopi*x2(ip))
        vp(2,ip) = SQRT(-2._prec*sigma2*LOG(x3(ip)))*COS(twopi*x4(ip))
     ENDDO
  ELSE
     sigma2 = 3._prec*k_boltzmann*Temp/mass ! in (Angstrom/s)^2
     CALL RANDOM_NUMBER(x5)
     CALL RANDOM_NUMBER(x6)
     DO ip=1,Npart
        vp(1,ip) = SQRT(-2._prec*sigma2*LOG(x1(ip)))*COS(twopi*x2(ip))
        vp(2,ip) = SQRT(-2._prec*sigma2*LOG(x3(ip)))*COS(twopi*x4(ip))
        vp(3,ip) = SQRT(-2._prec*sigma2*LOG(x5(ip)))*COS(twopi*x6(ip))
     ENDDO
  ENDIF
  !====================================================================!
  ! deallocation
  DEALLOCATE(x1)
  DEALLOCATE(x2)
  DEALLOCATE(x3)
  DEALLOCATE(x4)
  IF (ndim .EQ. 3) THEN
     DEALLOCATE(x5)
     DEALLOCATE(x6)
  ENDIF
9999 CONTINUE ! jump here upon error
END SUBROUTINE lj_set_particles
