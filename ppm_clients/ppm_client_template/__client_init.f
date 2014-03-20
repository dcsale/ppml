      !-------------------------------------------------------------------------
      ! Subroutine : client_init
      !-------------------------------------------------------------------------
      !
      ! Purpose : Initialize the client and read all input files.
      !
      ! Input :
      !
      ! Input/output :
      !
      ! Output : info (I) return status. 0 on success.
      !
      ! Routines : client_read_input
      ! client_defaults
      ! client_read_ctrl
      !
      !
      ! Remarks :
      !
      ! References :
      !
      ! Revisions :
      !-------------------------------------------------------------------------
      !
      !-------------------------------------------------------------------------
      ! Parallel Particle Mesh Library (PPM)
      ! ETH Zurich
      ! CH-8092 Zurich, Switzerland
      !-------------------------------------------------------------------------
      SUBROUTINE client_init(info)
      !-------------------------------------------------------------------------
      ! Modules
      !-------------------------------------------------------------------------
      USE ppm_module_io
      USE ppm_module_core
      USE client_global
      USE client_module_io
      IMPLICIT NONE
      !-------------------------------------------------------------------------
      ! Includes
      !-------------------------------------------------------------------------
      !-------------------------------------------------------------------------
      ! Arguments
      !-------------------------------------------------------------------------
      INTEGER, INTENT( OUT) :: info
      !-------------------------------------------------------------------------
      ! Local variables
      !-------------------------------------------------------------------------
      INTEGER :: ilen,ios,tolexp
      REAL(MK), DIMENSION(3) :: bsize
      CHARACTER(LEN = MAXCHAR) :: cbuf
      REAL(8) :: t0
      !-------------------------------------------------------------------------
      ! Externals
      !-------------------------------------------------------------------------
      !-------------------------------------------------------------------------
      ! Initialise
      !-------------------------------------------------------------------------
      info = 0
      !-------------------------------------------------------------------------
      ! Initialise number of processors
      !-------------------------------------------------------------------------
      Nproc = 1
      rank = 0
      !-------------------------------------------------------------------------
      ! Parse args and ctrl file
      !-------------------------------------------------------------------------
      CALL define_args
      CALL parse_args(info)
      IF (info .NE. 0) GOTO 9999
      !-------------------------------------------------------------------------
      ! Open log file
      !-------------------------------------------------------------------------
      IF (Nproc .LT. 10) THEN
          WRITE(cbuf,'(2A,I1.1)') TRIM(logfile),'.',rank
      ELSEIF (Nproc .LT. 100) THEN
          WRITE(cbuf,'(2A,I2.2)') TRIM(logfile),'.',rank
      ELSEIF (Nproc .LT. 1000) THEN
          WRITE(cbuf,'(2A,I3.3)') TRIM(logfile),'.',rank
      ELSE
          WRITE(cbuf,'(2A,I5.5)') TRIM(logfile),'.',rank
      ENDIF
      OPEN(client_log_unit,FILE=cbuf,ACTION='WRITE',IOSTAT=ios)
      CLOSE(client_log_unit,STATUS='DELETE')
      OPEN(client_log_unit,FILE=cbuf,ACTION='WRITE',IOSTAT=ios)
      !-------------------------------------------------------------------------
      ! Initialise the ppm library
      !-------------------------------------------------------------------------
      tolexp = INT(LOG10(EPSILON(lmyeps)))+10
      CALL ppm_init(3,MK,tolexp,comm,debug,info,ppm_log_unit)
      IF (info .NE. 0) THEN
          CALL pwrite(rank,'client_init','Failed to initialize PPM library.',info)
          GOTO 9999
      ENDIF
      !-------------------------------------------------------------------------
      ! Set ppm log file unit
      !-------------------------------------------------------------------------
      CALL ppm_io_set_unit(6,0,ppm_log_unit,info)
      IF (info .NE. 0) GOTO 9999
      !-------------------------------------------------------------------------
      ! Write some of the parameters
      !-------------------------------------------------------------------------
      CALL client_dump_parameters(info)
      !-------------------------------------------------------------------------
      ! Nullify data pointers
      !-------------------------------------------------------------------------
      NULLIFY(xp)
      !-------------------------------------------------------------------------
      ! Here comes all initialization code
      ! - setting up topology
      ! - setup initial condition of particles/meshes:
      ! * reading from files
      ! * programmatically
      ! - mapping particles/meshes onto topology
      ! - doing any other pre-processing needed for initial condition
      !-------------------------------------------------------------------------
      !-------------------------------------------------------------------------
      ! End of init
      !-------------------------------------------------------------------------
      CALL pwrite(rank,'client_init','Initialization complete.',info)
      IF(rank .EQ. 0) WRITE(client_log_unit,'(A)') 'Initialization complete.'
      !-------------------------------------------------------------------------
      ! Return
      !-------------------------------------------------------------------------
 9999 CONTINUE
      CALL substop('client_init',t0,info)
      RETURN
      END SUBROUTINE client_init
