MODULE lj_module_global
!======================================================================!
! contains global variables (for all procedures using this module)
!======================================================================!
  USE ppm_module_core
  USE ppm_module_typedef
  USE ppm_module_ctrl
  IMPLICIT NONE
  !====================================================================!
  ! parameters for initializing
  INTEGER, PARAMETER :: prec = 8
  INTEGER, PARAMETER :: ndim = 3
  !====================================================================!
  ! general stuff
  INTEGER :: debug,tolexp,nproc,rank,comm
  !====================================================================!
  ! computational domain
  REAL(prec),DIMENSION(ndim) :: mindom,maxdom ! lower and upper
                              ! bound of the computational domain
  REAL(prec),DIMENSION(ndim) :: lendom ! side length of
                              ! computational domain
  !====================================================================!
  ! particles
  INTEGER :: Npart_global ! number of particles in
                              ! whole domain (sum of Npart of all
                              ! processors
  INTEGER :: Npart ! number of particles on each processor
                              ! without ghost particles
  INTEGER :: Mpart ! number of particles on each processor
                              ! including ghost particles
  REAL(prec),DIMENSION(:,:),POINTER :: xp ! particle positions
  REAL(prec),DIMENSION(:,:),POINTER :: vp ! particle velocities
  !====================================================================!
  ! topology
  INTEGER :: decomp,assig,topo_id
  INTEGER :: nsublist,nsubs
  INTEGER ,DIMENSION(:) ,POINTER :: bcdef,sub2proc,isublist
  REAL(prec),DIMENSION(:) ,POINTER :: cost,proc_speed
  REAL(prec),DIMENSION(:,:),POINTER :: minsub,maxsub
  !====================================================================!
  ! ghosts...
  INTEGER,PARAMETER :: isymm = 1
  LOGICAL :: lsymm
  REAL(prec) :: ghostsize
  REAL(prec) :: skin ! skin on top of cutoff
                                             ! radius when creating ghostlayers
  !====================================================================!
  ! cell lists
  INTEGER :: npartner
  INTEGER, DIMENSION(:) ,POINTER :: mycell
  INTEGER, DIMENSION(:,:),POINTER :: ipartner,jpartner
  TYPE (ppm_t_clist),DIMENSION(:) ,POINTER :: clist
  REAL(prec),DIMENSION(ndim) :: cellsize
  !====================================================================!
  ! particle interactions
  REAL(prec),DIMENSION(:,:),POINTER :: accp
  REAL(prec),DIMENSION(3) :: comp_accel_arg
  !====================================================================!
  ! thermodynamics
  REAL(prec),PARAMETER :: Rgas = 8.314E20_prec ! in kg*Angstrom^2/(K*mol*s^2)
  REAL(prec),PARAMETER :: k_boltzmann = 1.38E-3_prec ! in kg*Angstrom^2/(K*s^2)
  !====================================================================!
  ! gas properties (Argon dimer)
  REAL(prec),PARAMETER :: mass = 6.69E-26_prec ! kg
  REAL(prec),PARAMETER :: sigma = 3.40_prec ! LJ parameter in Angstrom
  REAL(prec),PARAMETER :: epsil = 1.654E-1_prec ! LJ parameter in kg*Angstrom^2/s^2
  !====================================================================!
  ! INITIAL state of gas
  REAL(prec) :: sigma2 ! variance of velocity components
  REAL(prec) :: Temp ! temperature in K
  !====================================================================!
  ! artificial constraints on velocity and acceleration (for stability)
  REAL(prec) :: vmax,amax
  !====================================================================!
  ! numeric parameters
  INTEGER :: nsteps,framerate
  INTEGER :: it
  REAL(prec) :: dt
  REAL(prec) :: cutoff ! in Angstrom
  CONTAINS
  SUBROUTINE define_args
      CALL arg_group('General parameters')
      CALL arg(mindom, 'mindom', &
          long_flag = '--mindom', &
          ctrl_name = 'mindom', &
          default = (/0._prec, 0._prec, 0._prec/))
      CALL arg(maxdom, 'maxdom', &
          long_flag = '--maxdom', &
          ctrl_name = 'maxdom', &
          default = (/200._prec, 200._prec, 200._prec/))
      CALL arg(Npart_global, 'Npart_global', &
          flag = '-N', &
          long_flag = '--Npart-global', &
          ctrl_name = 'Npart_global', &
          default = 195, &
          min = 10)
      CALL arg(debug, 'debug', &
          flag = '-d', &
          long_flag = '--debug', &
          ctrl_name = 'debug', &
          default = 0)
      CALL arg(Temp, 'Temperature', &
          flag = '-N', &
          long_flag = '--Temp-global', &
          ctrl_name = 'Temperature', &
          default = 50._prec, &
          min = 0._prec)
      CALL arg_group('Numerics')
      CALL arg(nsteps, 'nsteps', &
          flag = '-n', &
          long_flag = '--nsteps', &
          ctrl_name = 'nsteps', &
          default = 10000, &
          min = 1)
      CALL arg(framerate, 'framerate', &
          flag = '-f', &
          long_flag = '--framerate', &
          ctrl_name = 'framerate', &
          default = 10, &
          min = 0)
      CALL arg(cutoff, 'cutoff', &
          flag = '-c', &
          long_flag = '--cutoff', &
          ctrl_name = 'cutoff', &
          default = 10._prec, &
          min = 0._prec)
      CALL arg(skin, 'skin', &
          flag = '-s', &
          long_flag = '--skin', &
          ctrl_name = 'skin', &
          default = 0._prec, &
          min = 0._prec)
  END SUBROUTINE define_args
END MODULE lj_module_global
