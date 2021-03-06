      !-------------------------------------------------------------------------
      !  Subroutine   :                    client_finalize
      !-------------------------------------------------------------------------
      !
      !  Purpose      : Clean up global memory and close open units. This
      !                 also finalizes ppm and mpi.
      !
      !  Input        : errstat       (I) error status. 0 if this is a
      !                                   normal termination. .NE. 0 if
      !                                   this is a program abort.
      !
      !  Input/output : 
      !
      !  Output       : info          (I) return status. 0 on success.
      !
      !  Routines     : substart
      !                 substop
      !                 ppm_finalize
      !                 MPI_Finalize
      !
      !  Remarks      :
      !
      !  References   :
      !-------------------------------------------------------------------------
      !  Parallel Particle Mesh Library (PPM)
      !  ETH Zurich
      !  CH-8092 Zurich, Switzerland
      !-------------------------------------------------------------------------

      SUBROUTINE client_finalize(errstat,info)

      !-------------------------------------------------------------------------
      !  Modules 
      !-------------------------------------------------------------------------
      USE client_global
      USE ppm_module_core
      USE ppm_module_neighlist
      IMPLICIT NONE

      !-------------------------------------------------------------------------
      !  Includes 
      !-------------------------------------------------------------------------
#ifdef __MPI
      INCLUDE 'mpif.h'
#endif
      !-------------------------------------------------------------------------
      !  Arguments     
      !-------------------------------------------------------------------------
      INTEGER, INTENT(IN   )   :: errstat
      INTEGER, INTENT(  OUT)   :: info
      !-------------------------------------------------------------------------
      !  Local variables 
      !-------------------------------------------------------------------------
      REAL(8)                  :: t0
      !-------------------------------------------------------------------------
      !  Externals 
      !-------------------------------------------------------------------------
      
      !-------------------------------------------------------------------------
      !  Initialise 
      !-------------------------------------------------------------------------
      CALL substart('client_finalize',t0,info)

      !-------------------------------------------------------------------------
      !  Deallocate global memory
      !-------------------------------------------------------------------------
      IF (ASSOCIATED(xp)) DEALLOCATE(xp,STAT=info)

      !-------------------------------------------------------------------------
      !  Deallocate neighbor lists
      !-------------------------------------------------------------------------
!       CALL ppm_clist_destroy(clist,info)
      IF (ASSOCIATED(Nm)) DEALLOCATE(Nm,STAT=info)

      !-------------------------------------------------------------------------
      !  Close log file
      !-------------------------------------------------------------------------
      CLOSE(client_log_unit,IOSTAT=info)

      !-------------------------------------------------------------------------
      !  Finalize PPM
      !-------------------------------------------------------------------------
      CALL ppm_finalize(info)

#ifdef __MPI
      IF (errstat .EQ. 0) THEN
          !---------------------------------------------------------------------
          !  Finalize MPI
          !---------------------------------------------------------------------
          CALL MPI_Finalize(info)
          IF (info .NE. 0) THEN
              WRITE(*,*) 'FAILED TO FINALIZE MPI. BAILING OUT!'
              GOTO 9999
          ENDIF
      ELSE
          !---------------------------------------------------------------------
          !  Abort MPI
          !---------------------------------------------------------------------
          CALL MPI_Abort(comm,errstat,info)
          IF (info .NE. 0) THEN
              WRITE(*,*) 'FAILED TO ABORT MPI!'
              GOTO 9999
          ENDIF
      ENDIF
#endif

      !-------------------------------------------------------------------------
      !  Return 
      !-------------------------------------------------------------------------
 9999 CONTINUE
      CALL substop('client_finalize',t0,info)
      RETURN
      END SUBROUTINE client_finalize
