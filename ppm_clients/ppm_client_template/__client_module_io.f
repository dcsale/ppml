      !-------------------------------------------------------------------------
      ! Module : client_module_io
      !-------------------------------------------------------------------------
      !
      ! Purpose : This module contains all subroutines needed for
      ! file I/O.
      !
      ! Remarks :
      !
      ! References :
      !-------------------------------------------------------------------------
      ! Parallel Particle Mesh Library (PPM)
      ! ETH Zurich
      ! CH-8092 Zurich, Switzerland
      !-------------------------------------------------------------------------
      MODULE client_module_io
         !----------------------------------------------------------------------
         ! Define interface to client_write_output
         !----------------------------------------------------------------------
         INTERFACE client_write_output
             MODULE PROCEDURE client_write_output
         END INTERFACE
         !----------------------------------------------------------------------
         ! Define interface to client_write_diag
         !----------------------------------------------------------------------
         INTERFACE client_write_diag
             MODULE PROCEDURE client_write_diag
         END INTERFACE
         !----------------------------------------------------------------------
         ! Include the source
         !----------------------------------------------------------------------
         CONTAINS
      !-------------------------------------------------------------------------
      ! Subroutine : client_write_output
      !-------------------------------------------------------------------------
      !
      ! Purpose : This routine writes an output file for
      ! visualization.
      !
      ! Input : istep (I) time step number completed
      ! iUnit (I) IO unit to be used for file output
      !
      ! The global variables xp,wp are read.
      !
      ! Input/output :
      !
      ! Output : info (I) return status (0 if no error)
      !
      ! Routines : ppm_io_open
      ! ppm_io
      ! ppm_io_close
      !
      ! Remarks :
      !
      ! References :
      !-------------------------------------------------------------------------
      ! Parallel Particle Mesh Library (PPM)
      ! ETH Zurich
      ! CH-8092 Zurich, Switzerland
      !-------------------------------------------------------------------------
      SUBROUTINE client_write_output(istep,iUnit,info)
      !-------------------------------------------------------------------------
      ! Modules
      !-------------------------------------------------------------------------
      USE client_global
      USE ppm_module_io
      IMPLICIT NONE
      !-------------------------------------------------------------------------
      ! Arguments
      !-------------------------------------------------------------------------
      INTEGER , INTENT(IN ) :: istep,iUnit
      INTEGER , INTENT( OUT) :: info
      !-------------------------------------------------------------------------
      ! Local variables
      !-------------------------------------------------------------------------
      REAL(8) :: t0
      INTEGER , DIMENSION(4 ) :: ivec
      REAL(MK), DIMENSION(2 ) :: fvec
      INTEGER :: ifmt,jUnit,i
      CHARACTER(LEN=MAXCHAR) :: filename,cbuf
      !-------------------------------------------------------------------------
      ! Externals
      !-------------------------------------------------------------------------
      !-------------------------------------------------------------------------
      ! Initialise
      !-------------------------------------------------------------------------
      CALL substart('client_write_output',t0,info)
      WRITE(filename,'(A,I8.8,A)') TRIM(outputfile),istep,'.out'
      !-------------------------------------------------------------------------
      ! Confirm
      !-------------------------------------------------------------------------
      IF (rank .EQ. 0) THEN
          WRITE(cbuf,'(2A)') 'Output written to file ',TRIM(filename)
          CALL pwrite(rank,'client_write_output',cbuf,info)
          WRITE(client_log_unit,'(A)') TRIM(cbuf)
      ENDIF
      !-------------------------------------------------------------------------
      ! Return
      !-------------------------------------------------------------------------
 9999 CONTINUE
      CALL substop('client_write_output',t0,info)
      RETURN
      END SUBROUTINE client_write_output
      !-------------------------------------------------------------------------
      ! Subroutine : client_write_diag
      !-------------------------------------------------------------------------
      !
      ! Purpose : This routine writes diagnostics to a file.
      ! For the parallel version, all processors send
      ! their information to rank 0, who writes the file.
      !
      ! Input : Unit (I) IO unit to be used for file output
      ! lfirst (L) .TRUE. if this is the first
      ! CALL to this routine in the current
      ! run. If this is not a restart run,
      ! it will cause old diag files to be
      ! overwritten. Otherwise, data is
      ! appended.
      !
      ! The global variables xp,wp are read.
      !
      ! Input/output :
      !
      ! Output : info (I) return status (0 if no error)
      !
      ! Routines : ppm_io_open
      ! ppm_io
      ! ppm_io_close
      ! substart
      ! substop
      !
      ! Remarks :
      !
      ! References :
      !-------------------------------------------------------------------------
      ! Parallel Particle Mesh Library (PPM)
      ! ETH Zurich
      ! CH-8092 Zurich, Switzerland
      !-------------------------------------------------------------------------
      SUBROUTINE client_write_diag(iUnit,lfirst,info)
      !-------------------------------------------------------------------------
      ! Modules
      !-------------------------------------------------------------------------
      USE client_global
      USE ppm_module_io
      IMPLICIT NONE
      !-------------------------------------------------------------------------
      ! Arguments
      !-------------------------------------------------------------------------
      INTEGER , INTENT(IN ) :: iUnit
      LOGICAL , INTENT(IN ) :: lfirst
      INTEGER , INTENT( OUT) :: info
      !-------------------------------------------------------------------------
      ! Local variables
      !-------------------------------------------------------------------------
      REAL(8) :: t0
      CHARACTER(LEN=MAXCHAR) :: cbuf
      INTEGER :: ifmt
      !-------------------------------------------------------------------------
      ! Externals
      !-------------------------------------------------------------------------
      !-------------------------------------------------------------------------
      ! Initialise
      !-------------------------------------------------------------------------
      CALL substart('client_write_diag',t0,info)
      !-------------------------------------------------------------------------
      ! Check if diagnostics are enabled
      !-------------------------------------------------------------------------
      IF (LEN_TRIM(diagfile) .LT. 1) THEN
          CALL pwrite(rank,'client_write_diag', &
     & 'Diagnostics are disabled. Skipping.',info)
          GOTO 9999
      ENDIF
      !-------------------------------------------------------------------------
      ! Format
      !-------------------------------------------------------------------------
      IF (diagfmt(1:1) .EQ. 'f' .OR. diagfmt(1:1) .EQ. 'F') THEN
          ifmt = ppm_param_io_ascii
      ELSE
          ifmt = ppm_param_io_binary
      ENDIF
      !-------------------------------------------------------------------------
      ! Compute diagnostics and write to file using the ppm I/O routines
      !-------------------------------------------------------------------------
      !-------------------------------------------------------------------------
      ! Confirm
      !-------------------------------------------------------------------------
      IF (rank .EQ. 0) THEN
          IF (lfirst) THEN
              WRITE(cbuf,'(2A)') 'Diagnostics written to file ',TRIM(diagfile)
          ELSE
              WRITE(cbuf,'(2A)') 'Diagnostics appended to file ',TRIM(diagfile)
          ENDIF
          CALL pwrite(rank,'client_write_diag',cbuf,info)
          WRITE(client_log_unit,'(A)') TRIM(cbuf)
      ENDIF
      !-------------------------------------------------------------------------
      ! Return
      !-------------------------------------------------------------------------
 9999 CONTINUE
      CALL substop('client_write_diag',t0,info)
      RETURN
      END SUBROUTINE client_write_diag
      END MODULE client_module_io
