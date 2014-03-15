      !-------------------------------------------------------------------------
      !  Subroutine   :                  client_write_output
      !-------------------------------------------------------------------------
      !
      !  Purpose      : This routine writes an output file for
      !                 visualization.
      !
      !  Input        : istep      (I) time step number completed
      !                 iUnit      (I) IO unit to be used for file output
      !
      !                 The global variables xp,wp are read.
      !
      !  Input/output : 
      !
      !  Output       : info       (I) return status (0 if no error)
      !
      !  Routines     : ppm_io_open
      !                 ppm_io
      !                 ppm_io_close
      !
      !  Remarks      :
      !
      !  References   :
      !-------------------------------------------------------------------------
      !  Parallel Particle Mesh Library (PPM)
      !  ETH Zurich
      !  CH-8092 Zurich, Switzerland
      !-------------------------------------------------------------------------

      SUBROUTINE client_write_output(istep,iUnit,info)

      !-------------------------------------------------------------------------
      !  Modules
      !-------------------------------------------------------------------------
      USE client_global
      USE ppm_module_io
      IMPLICIT NONE
      !-------------------------------------------------------------------------
      !  Arguments     
      !-------------------------------------------------------------------------
      INTEGER                 , INTENT(IN   ) :: istep,iUnit
      INTEGER                 , INTENT(  OUT) :: info
      !-------------------------------------------------------------------------
      !  Local variables 
      !-------------------------------------------------------------------------
      REAL(8)                                 :: t0
      INTEGER , DIMENSION(4  )                :: ivec
      REAL(MK), DIMENSION(2  )                :: fvec
      INTEGER                                 :: ifmt,jUnit,i
      CHARACTER(LEN=MAXCHAR)                  :: filename,cbuf
      !-------------------------------------------------------------------------
      !  Externals 
      !-------------------------------------------------------------------------
      
      !-------------------------------------------------------------------------
      !  Initialise 
      !-------------------------------------------------------------------------
      CALL substart('client_write_output',t0,info)
      WRITE(filename,'(A,I8.8,A)') TRIM(outputfile),istep,'.out'













      !-------------------------------------------------------------------------
      !  Confirm
      !-------------------------------------------------------------------------
      IF (rank .EQ. 0) THEN
          WRITE(cbuf,'(2A)') 'Output written to file ',TRIM(filename)
          CALL pwrite(rank,'client_write_output',cbuf,info)
          WRITE(client_log_unit,'(A)') TRIM(cbuf)
      ENDIF

      !-------------------------------------------------------------------------
      !  Return 
      !-------------------------------------------------------------------------
 9999 CONTINUE
      CALL substop('client_write_output',t0,info)
      RETURN
      END SUBROUTINE client_write_output
