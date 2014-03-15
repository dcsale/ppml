SUBROUTINE lj_finalize(info)
!======================================================================!
! finalization of Lennard-Jones client
!======================================================================!
  USE lj_module_global
  USE ppm_module_core
  IMPLICIT NONE
  INCLUDE 'mpif.h'
  ! arguments
  INTEGER,INTENT(OUT) :: info
  ! local variables
  !====================================================================!
  ! initialize
  info = 0
  !====================================================================!
  ! deallocation
  DEALLOCATE(xp)
  DEALLOCATE(vp)
  DEALLOCATE(accp)
  !====================================================================!
  ! finanlize the PPM library
  CALL ppm_finalize(info)
  IF (info .NE. 0) THEN
     WRITE(*,*)'ERROR in lj_finalize: ppm_finalize failed.'
     GOTO 9999
  ENDIF
  !====================================================================!
  ! finalize MPI
  CALL MPI_Finalize(info)
  IF (info .NE. 0) THEN
     WRITE(*,*)'ERROR in lj_finalize: MPI_Finalize failed.'
     GOTO 9999
  ENDIF
9999 CONTINUE ! jump here upon error
END SUBROUTINE lj_finalize
