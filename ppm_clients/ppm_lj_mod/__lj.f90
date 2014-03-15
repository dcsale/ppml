PROGRAM lj
!======================================================================!
! Lennard-Jones client: simulates the motion of gas molecules living in
! a force field governed by the Lennard-Jones potential
!
! Author: Birte Schrader and Philipp Boenhof
!-------------------------------------------------------------------------
! Parallel Particle Mesh Library (PPM)
! ETH Zurich
! CH-8092 Zurich, Switzerland
!-------------------------------------------------------------------------
  USE lj_module_global
  USE lj_module_writeout
  IMPLICIT NONE
  INCLUDE 'mpif.h'
  INTEGER :: info,ip
  CHARACTER(LEN=256) :: filename,myformat
  REAL(KIND(1.d0)) :: t1,t2
  !====================================================================!
  ! initialize
  info = 0
  CALL lj_init(info)
  IF (info .NE. 0) THEN
     WRITE(*,*)'ERROR in lj (main): lj_init failed.'
     GOTO 9999
  ENDIF
  WRITE(myformat,'(A,I1,A)'),'(A,I',CEILING(LOG10(REAL(nproc))),',A,I8,A)'
  !====================================================================!
  ! simulation
  DO it=1,nsteps
     t1 = MPI_WTIME(info)
     CALL lj_move_particles(info)
     IF (info .NE. 0) THEN
        WRITE(*,*)'ERROR in lj (main): lj_move_particles failed.'
        GOTO 9999
     ENDIF
     t2 = MPI_WTIME(info)
     WRITE(*,*)'lj_move_particles: ',t2-t1,' seconds'
     t1 = MPI_WTIME(info)
     CALL lj_remap(info)
     IF (info .NE. 0) THEN
        WRITE(*,*)'ERROR in lj (main): lj_remap failed.'
        GOTO 9999
     ENDIF
     t2 = MPI_WTIME(info)
     WRITE(*,*)'lj_remap: ',t2-t1,' seconds'
     t1 = MPI_WTIME(info)
     !=================================================================!
     ! writeout
     IF (MOD(it-1,framerate) .EQ. 0) THEN
        filename = 'processor.pdb'
        CALL lj_writeout(filename,info)
        IF (info .NE. 0) THEN
           WRITE(*,*)'ERROR in lj (main): lj_writeout failed.'
           GOTO 9999
        ENDIF
     ENDIF
     WRITE(*,myformat)'[',rank,'] (lj) : completed time step ',it
     t2 = MPI_WTIME(info)
     WRITE(*,*)'output: ',t2-t1,' seconds'
     t1 = MPI_WTIME(info)
  ENDDO
  !====================================================================!
  ! finalize
  CALL lj_finalize(info)
  IF (info .NE. 0) THEN
     WRITE(*,*)'ERROR in lj (main): lj_finalize failed.'
     GOTO 9999
  ENDIF
  t2 = MPI_WTIME(info)
  WRITE(*,*)'finalize: ',t2-t1,' seconds'
9999 CONTINUE ! jump here upon error
  IF (info .NE. 0) WRITE(*,*)'Controlled termination on error.'
END PROGRAM lj
