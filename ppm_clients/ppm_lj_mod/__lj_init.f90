SUBROUTINE lj_init(info)
!======================================================================!
! initialization of Lennard-Jones client
!======================================================================!
  USE ppm_module_core
  USE lj_module_global
  IMPLICIT NONE
  INCLUDE 'mpif.h'
  ! arguments
  INTEGER,INTENT(OUT) :: info
  ! local variables
  INTEGER :: ip
  !====================================================================!
  ! initialize some general variables
  info = 0 ! change if error occurs
  tolexp = LOG10(EPSILON(1._prec)) ! exponent of machine epsilon
  !====================================================================!
  ! initialize MPI
  comm = MPI_COMM_WORLD
  CALL MPI_Init(info)
  IF (info .NE. 0) THEN
     WRITE(*,*)'ERROR in lj_init: MPI_Init failed.'
     GOTO 9999
  ENDIF
  CALL MPI_Comm_Size(comm,nproc,info)
  IF (info .NE. 0) THEN
     WRITE(*,*)'ERROR in lj_init: MPI_Comm_Size failed.'
     GOTO 9999
  ENDIF
  CALL MPI_Comm_Rank(comm,rank,info)
  IF (info .NE. 0) THEN
     WRITE(*,*)'ERROR in lj_init: MPI_Comm_Rank failed.'
     GOTO 9999
  ENDIF
  !====================================================================!
  ! get input parameter from command line or from control file
  CALL define_args()
  CALL disable_ctrl() !comment out if a control file is provided
  CALL parse_args(info)
  if (info .NE. 0 ) THEN
      GOTO 9999
  endif
  !====================================================================!
  ! intialize the PPM library
  CALL ppm_init(ndim,prec,tolexp,comm,debug,info)
  IF (info .NE. 0) THEN
     WRITE(*,*)'ERROR in lj_init: ppm_init failed.'
     GOTO 9999
  ENDIF
  !====================================================================!
  ! define domain length in Angstrom = 10^-10 m
  lendom = maxdom - mindom
  !====================================================================!
  ! set particles (initial positions and velocities of the gas molecules)
  CALL lj_set_particles(info)
  IF (info .NE. 0) THEN
     WRITE(*,*)'ERROR in lj_init: lj_set_particles failed.'
     GOTO 9999
  ENDIF
  !======================================================================!
  ! create new topology, map particles, create/map ghosts, and build cell
  ! lists
  CALL lj_newtopo(info)
  IF (info .NE. 0) THEN
     WRITE(*,*)'ERROR in lj_init: lj_newtopo failed.'
     GOTO 9999
  ENDIF
  !====================================================================!
  ! compute a reasonable time stepsize (characteristic length /
  ! characteristic speed * safety factor) and bound velocities and
  ! accelerations for numeric stability, necessary mainly because of
  ! the (unrealistical) random initial positions of the molecules
  dt = 0.1_prec*SQRT(mass*sigma**2/epsil)
  vmax = 0.01*sigma/dt ! max. velocity (num. cutoff!)
  amax = vmax/dt ! max. acceleration (num. cutoff!)
  IF (ndim .EQ. 2) THEN
     DO ip=1,Mpart
        vp(1,ip) = SIGN(MIN(vmax,ABS(vp(1,ip))),vp(1,ip))
        vp(2,ip) = SIGN(MIN(vmax,ABS(vp(2,ip))),vp(2,ip))
     ENDDO
  ELSE
     DO ip=1,Mpart
        vp(1,ip) = SIGN(MIN(vmax,ABS(vp(1,ip))),vp(1,ip))
        vp(2,ip) = SIGN(MIN(vmax,ABS(vp(2,ip))),vp(2,ip))
        vp(3,ip) = SIGN(MIN(vmax,ABS(vp(3,ip))),vp(3,ip))
     ENDDO
  ENDIF
  !====================================================================!
  ! allocate accp
  IF (lsymm) THEN
     ALLOCATE(accp(ndim,Mpart))
     IF (info .NE. 0) THEN
        WRITE(*,*)'ERROR in lj_init: Allocation failed.'
        GOTO 9999
     ENDIF
  ELSE
     ALLOCATE(accp(ndim,Npart))
     IF (info .NE. 0) THEN
        WRITE(*,*)'ERROR in lj_init: Allocation failed.'
        GOTO 9999
     ENDIF
  ENDIF
9999 CONTINUE ! jump here upon error
END SUBROUTINE lj_init
