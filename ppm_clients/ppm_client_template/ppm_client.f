      !-------------------------------------------------------------------------
      !  Program      :                      ppm_client
      !-------------------------------------------------------------------------
      !
      !  Purpose      : A simple client template
      !
      !  Note         : 
      !-------------------------------------------------------------------------
      !  Parallel Particle Mesh Library (PPM)
      !  ETH Zurich
      !  CH-8092 Zurich, Switzerland
      !-------------------------------------------------------------------------

      PROGRAM ppm_client

      !-------------------------------------------------------------------------
      !  Modules
      !-------------------------------------------------------------------------
      USE ppm_module_data
      USE client_global
      USE client_module_io
      IMPLICIT NONE
      !-------------------------------------------------------------------------
      !  Includes
      !-------------------------------------------------------------------------
#ifdef __MPI
      INCLUDE 'mpif.h'
#endif
      !-------------------------------------------------------------------------
      !  Local variables
      !-------------------------------------------------------------------------
      INTEGER               :: info
      INTEGER, DIMENSION(6) :: bcdef
      REAL(MK)              :: t0
      INTEGER               :: istep
      !-------------------------------------------------------------------------
      !  Externals
      !-------------------------------------------------------------------------

      !-------------------------------------------------------------------------
      !  Initialise
      !-------------------------------------------------------------------------
      bcdef(1:6) = ppm_param_bcdef_freespace
      info = 0
      CALL client_init(info)
      IF (info .EQ. exit_gracefully) THEN
         GOTO 7000
      ELSE IF (info .LT. 0) THEN
          CALL pwrite(rank,'ppm_client','Failed to initialize.',info)
          GOTO 8000
      ENDIF
      IF (info .GT. 0) GOTO 9999
      CALL substart('ppm_client',t0,info)

      !-------------------------------------------------------------------------
      !  Write initial diagnostics
      !-------------------------------------------------------------------------
      CALL client_write_diag(20,.TRUE.,info)
      IF (info .NE. 0) THEN
          CALL pwrite(rank,'ppm_client',   &
     &        'Writing intitial diagnostics failed.',info)
          GOTO 8000
      ENDIF

      !-------------------------------------------------------------------------
      !  Write initial condition to output file
      !-------------------------------------------------------------------------
      CALL client_write_output(0,20,info)
      IF (info .NE. 0) THEN
          CALL pwrite(rank,'ppm_client','Writing intitial condition failed.',&
     &                info)
          GOTO 8000
      ENDIF

      
      !-------------------------------------------------------------------------
      !  Do time stepping
      !-------------------------------------------------------------------------
      istep = 0
      DO WHILE (.FALSE.) ! fill in loop termination condition
         istep = istep + 1
         !----------------------------------------------------------------------
         !  Place here your main numerical code that is executed at every
         !  time step.
         !----------------------------------------------------------------------
 


         !----------------------------------------------------------------------
         !  Write diagnostics file if needed
         !----------------------------------------------------------------------
         IF (MOD(istep,freqdiag) .EQ. 0) THEN
             CALL client_write_diag(20,.FALSE.,info)
             IF (info .NE. 0) THEN
                 CALL pwrite(rank,'ppm_client','Writing diagnostics failed.',info)
                 GOTO 8000
             ENDIF
         ENDIF
      ENDDO               ! time loop




      !-------------------------------------------------------------------------
      !  Finalize
      !-------------------------------------------------------------------------
7000  IF (rank.EQ.0) THEN
         CALL pwrite(rank,'ppm_client','Done. Goodbye...',info)
      ENDIF
      CALL client_finalize(0,info)

      !-------------------------------------------------------------------------
      !  Skip error handler
      !-------------------------------------------------------------------------
      GOTO 9999

      !-------------------------------------------------------------------------
      !  Error handler
      !-------------------------------------------------------------------------
 8000 CONTINUE
      CALL client_finalize(1,info)
      IF (ASSOCIATED(xp)) DEALLOCATE(xp,STAT=info)

      !-------------------------------------------------------------------------
      !  End
      !-------------------------------------------------------------------------
 9999 CONTINUE
      CALL substop('ppm_client',t0,info)

      END PROGRAM ppm_client
