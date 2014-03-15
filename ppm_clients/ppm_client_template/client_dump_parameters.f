      !-------------------------------------------------------------------------
      !  Subroutine   :                 client_dump_parameters
      !-------------------------------------------------------------------------
      !
      !  Purpose      : Prints some of the parameters to stdout as well as
      !                 the log file.
      !
      !  Input        :
      !
      !  Input/output : 
      !
      !  Output       : info       (I) return status. 0 on success.
      !
      !  Routines     :
      !
      !  Remarks      :
      !
      !  References   :
      !-------------------------------------------------------------------------
      !  Parallel Particle Mesh Library (PPM)
      !  ETH Zurich
      !  CH-8092 Zurich, Switzerland
      !-------------------------------------------------------------------------

      SUBROUTINE client_dump_parameters(info)
      !-------------------------------------------------------------------------
      !  Modules 
      !-------------------------------------------------------------------------
      USE client_global
      IMPLICIT NONE
      !-------------------------------------------------------------------------
      !  Arguments     
      !-------------------------------------------------------------------------
      INTEGER, INTENT(   OUT)  :: info
      !-------------------------------------------------------------------------
      !  Local variables 
      !-------------------------------------------------------------------------
      INTEGER, DIMENSION(8)    :: tval
      CHARACTER(LEN = MAXCHAR) :: cbuf
      CHARACTER(LEN = 20)      :: cfmt
      REAL(8)                  :: t0
      INTEGER                  :: ilen
      !-------------------------------------------------------------------------
      !  Externals 
      !-------------------------------------------------------------------------
      
      !-------------------------------------------------------------------------
      !  Initialise 
      !-------------------------------------------------------------------------
      CALL substart('client_dump_parameters',t0,info)

      !-------------------------------------------------------------------------
      !  Root only
      !-------------------------------------------------------------------------
      IF (rank .EQ. 0) THEN
         !----------------------------------------------------------------------
         !  Get time and date 
         !----------------------------------------------------------------------
         CALL DATE_AND_TIME(values=tval)
         WRITE(cbuf,'(A,I4.4,A,5(I2.2,A))') '   Date and Time of Start:   ', &
     &      tval(1),'-',tval(2),'-',tval(3),' ',tval(5),':',tval(6),':',     &
     &      tval(7),' '
         ilen = LEN_TRIM(cbuf)
         WRITE(client_log_unit,'(A)') cbuf(1:ilen)
         WRITE(*,'(A)') cbuf(1:ilen)
         WRITE(cbuf,'(A,I4)') '   Number of processors:     ',nproc
         ilen = LEN_TRIM(cbuf)
         WRITE(client_log_unit,'(A)') cbuf(1:ilen)
         WRITE(*,'(A)') cbuf(1:ilen)
         !----------------------------------------------------------------------
         !  Control and debug (command line parameters)
         !----------------------------------------------------------------------
         WRITE(cbuf,'(2A)') '   Control file read:        ',TRIM(ctrlfile)
         ilen = LEN_TRIM(cbuf)
         WRITE(client_log_unit,'(A)') cbuf(1:ilen)
         WRITE(*,'(A)') cbuf(1:ilen)
         WRITE(cbuf,'(A,I1.1)') '   Debug level for this run: ',debug
         ilen = LEN_TRIM(cbuf)
         WRITE(client_log_unit,'(A)') cbuf(1:ilen)
         WRITE(*,'(A)') cbuf(1:ilen)
         WRITE(cbuf,'(2A)') '   Log file name:            ',TRIM(logfile)
         ilen = LEN_TRIM(cbuf)
         WRITE(client_log_unit,'(A)') cbuf(1:ilen)
         WRITE(*,'(A)') cbuf(1:ilen)
         !----------------------------------------------------------------------
         !  Run case name
         !----------------------------------------------------------------------
         WRITE(cbuf,'(A)') '---- CASE DETAILS ---------------------------------'
         ilen = LEN_TRIM(cbuf)
         WRITE(client_log_unit,'(A)') cbuf(1:ilen)
         WRITE(*,'(A)') cbuf(1:ilen)
         WRITE(cbuf,'(2A)') '   This is problem case      ',TRIM(casename)
         ilen = LEN_TRIM(cbuf)
         WRITE(client_log_unit,'(A)') cbuf(1:ilen)
         WRITE(*,'(A)') cbuf(1:ilen)
         !----------------------------------------------------------------------
         !  Domain decomposition
         !----------------------------------------------------------------------
         WRITE(cbuf,'(A)') '---- DOMAIN DECOMPOSITION -------------------------'
         ilen = LEN_TRIM(cbuf)
         WRITE(client_log_unit,'(A)') cbuf(1:ilen)
         WRITE(*,'(A)') cbuf(1:ilen)
         WRITE(cbuf,'(2A)') '   Domain decomposition used: octtree'
         ilen = LEN_TRIM(cbuf)
         WRITE(client_log_unit,'(A)') cbuf(1:ilen)
         WRITE(*,'(A)') cbuf(1:ilen)
         !----------------------------------------------------------------------
         !  Diagnostics
         !----------------------------------------------------------------------
         WRITE(cbuf,'(A)') '---- DIAGNOSTIC OUTPUT ----------------------------'
         ilen = LEN_TRIM(cbuf)
         WRITE(client_log_unit,'(A)') cbuf(1:ilen)
         WRITE(*,'(A)') cbuf(1:ilen)
         IF (LEN_TRIM(diagfile) .LT. 1) THEN
            WRITE(cbuf,'(2A)') '   Diagnostics file:         ','disabled'
            ilen = LEN_TRIM(cbuf)
            WRITE(client_log_unit,'(A)') cbuf(1:ilen)
            WRITE(*,'(A)') cbuf(1:ilen)
         ELSE
            WRITE(cbuf,'(A,I6)') '   Time steps between diag:  ',freqdiag
            ilen = LEN_TRIM(cbuf)
            WRITE(client_log_unit,'(A)') cbuf(1:ilen)
            WRITE(*,'(A)') cbuf(1:ilen)
            WRITE(cbuf,'(2A)') '   Diagnostics file:         ',TRIM(diagfile)
            ilen = LEN_TRIM(cbuf)
            WRITE(client_log_unit,'(A)') cbuf(1:ilen)
            WRITE(*,'(A)') cbuf(1:ilen)
         ENDIF
         !----------------------------------------------------------------------
         !  Output
         !----------------------------------------------------------------------
         WRITE(cbuf,'(A)') '---- RESULT DATA OUTPUT ---------------------------'
         ilen = LEN_TRIM(cbuf)
         WRITE(client_log_unit,'(A)') cbuf(1:ilen)
         WRITE(*,'(A)') cbuf(1:ilen)
         WRITE(cbuf,'(A,I6)') '   Time steps between output:',freqoutput
         ilen = LEN_TRIM(cbuf)
         WRITE(client_log_unit,'(A)') cbuf(1:ilen)
         WRITE(*,'(A)') cbuf(1:ilen)
         WRITE(cbuf,'(2A)') '   Output file name stub:    ',TRIM(outputfile)
         ilen = LEN_TRIM(cbuf)
         WRITE(client_log_unit,'(A)') cbuf(1:ilen)
         WRITE(*,'(A)') cbuf(1:ilen)
         !----------------------------------------------------------------------
         !  Flags
         !----------------------------------------------------------------------
         IF (checkmap) THEN
            WRITE(cbuf,'(2A)') '   Check mapping results:    ','yes'
         ELSE
            WRITE(cbuf,'(2A)') '   Check mapping results:    ','no'
         ENDIF
         ilen = LEN_TRIM(cbuf)
         WRITE(client_log_unit,'(A)') cbuf(1:ilen)
         WRITE(*,'(A)') cbuf(1:ilen)
         IF (probeproc) THEN
            WRITE(cbuf,'(2A)') '   Probe processor speeds:   ','yes'
         ELSE
            WRITE(cbuf,'(2A)') '   Probe processor speeds:   ','no'
         ENDIF
         ilen = LEN_TRIM(cbuf)
         WRITE(client_log_unit,'(A)') cbuf(1:ilen)
         WRITE(*,'(A)') cbuf(1:ilen)
         IF (dumpparticles) THEN
            WRITE(cbuf,'(2A)') '   Dump particle positions:  ','yes'
         ELSE
            WRITE(cbuf,'(2A)') '   Dump particle positions:  ','no'
         ENDIF
         ilen = LEN_TRIM(cbuf)
         WRITE(client_log_unit,'(A)') cbuf(1:ilen)
         WRITE(*,'(A)') cbuf(1:ilen)
      ENDIF

      !-------------------------------------------------------------------------
      !  Return 
      !-------------------------------------------------------------------------
 9999 CONTINUE
      CALL substop('client_dump_parameters',t0,info)
      RETURN
      END SUBROUTINE client_dump_parameters
