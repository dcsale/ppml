      !-------------------------------------------------------------------------
      !  Subroutine   :                     client_time
      !-------------------------------------------------------------------------
      !
      !  Purpose      : Retreives the current time (cpu). Uses either 
      !                 MPI_Wtime, etime or the f90 CPU_TIME intrinsic, 
      !                 based on the DEFINES that are set (see
      !                 ppm_define.h):
      !                     __MPI       uses MPI_Wtime
      !                     __ETIME     uses etime
      !                     none        uses CPU_TIME
      !
      !  Input        :
      !
      !  Input/output : 
      !
      !  Output       : timing    (F) current CPU clock time
      !
      !  Routines     : etime (C intrinsic)
      !
      !  Remarks      :
      !
      !  References   :
      !-------------------------------------------------------------------------
      !  Parallel Particle Mesh Library (PPM)
      !  ETH Zurich
      !  CH-8092 Zurich, Switzerland
      !-------------------------------------------------------------------------

      SUBROUTINE client_time(timing)

      !-------------------------------------------------------------------------
      !  Includes
      !-------------------------------------------------------------------------
      IMPLICIT NONE
#ifdef __MPI
      INCLUDE 'mpif.h'
#endif
      !-------------------------------------------------------------------------
      !  Arguments     
      !-------------------------------------------------------------------------
      REAL(8), INTENT(  OUT) :: timing
      !-------------------------------------------------------------------------
      !  Local variables 
      !-------------------------------------------------------------------------
#ifdef __ETIME
      REAL(8) :: array(2)
      !-------------------------------------------------------------------------
      !  Externals 
      !-------------------------------------------------------------------------
      REAL(8), EXTERNAL :: etime
#endif
      
#ifdef __MPI
      !-------------------------------------------------------------------------
      !  Call the MPI function MPI_Wtime
      !-------------------------------------------------------------------------
      timing = MPI_Wtime()
#else
      !-------------------------------------------------------------------------
      !  Call the C routine: etime to get the cpu time
      !-------------------------------------------------------------------------
#ifdef __ETIME
      timing = etime(array)
#else
      CALL CPU_TIME(timing)
#endif
#endif

      !-------------------------------------------------------------------------
      !  Return 
      !-------------------------------------------------------------------------
 9999 CONTINUE
      RETURN
      END SUBROUTINE client_time
