      !-------------------------------------------------------------------------
      !  Module       :                    client_global
      !-------------------------------------------------------------------------
      !
      !  Purpose      :  Declare global types and variables for ppm_rc.
      !                
      !  Remarks      :
      !
      !  References   :
      !
      !  Revisions    :
      !-------------------------------------------------------------------------
      !  Parallel Particle Mesh Library (PPM)
      !  ETH Zurich
      !  CH-8092 Zurich, Switzerland
      !-------------------------------------------------------------------------

      MODULE client_global
         !----------------------------------------------------------------------
         !  Modules
         !----------------------------------------------------------------------
         USE ppm_module_ctrl
         USE ppm_module_typedef
         USE ppm_module_substart
         USE ppm_module_substop
!          USE ppm_module_neighlist, ONLY: ppm_type_ptr_to_clist
!          PRIVATE :: ppm_type_ptr_to_clist

         !----------------------------------------------------------------------
         !  Global TYPEs
         !----------------------------------------------------------------------
         TYPE ptr_to_list
             INTEGER, DIMENSION(:), POINTER            :: list
         END TYPE

         !----------------------------------------------------------------------
         !  Header file
         !----------------------------------------------------------------------
         INCLUDE 'client_param.h'

         !----------------------------------------------------------------------
         !  Particle positions and properties
         !----------------------------------------------------------------------
         ! Particle positions. 1st index: x,y,z; 2nd: ipart
         REAL(MK), DIMENSION(:,:), POINTER             :: xp

         !----------------------------------------------------------------------
         !  Topology
         !----------------------------------------------------------------------
         REAL(MK), DIMENSION(3)                        :: min_phys,max_phys
         REAL(MK), DIMENSION(:  ), POINTER             :: proc_speed
         REAL(MK), DIMENSION(:  ), POINTER             :: sub_cost
         TYPE(ppm_t_topo), POINTER                     :: topo

         !----------------------------------------------------------------------
         !  Domain decomposition
         !----------------------------------------------------------------------
         INTEGER                                       :: ppm_decomp


         !----------------------------------------------------------------------
         !  Number of particles on local processor
         !----------------------------------------------------------------------
         ! without ghosts (i.e. real particles)
         INTEGER                                       :: Npart
         ! total particles including the ghosts
         INTEGER                                       :: Mpart

         !----------------------------------------------------------------------
         !  Program control flags
         !----------------------------------------------------------------------
         LOGICAL                                       :: checkmap
         LOGICAL                                       :: probeproc
         LOGICAL                                       :: dumpparticles

         !----------------------------------------------------------------------
         !  Float comparison tolerance
         !----------------------------------------------------------------------
         REAL(MK)                                      :: lmyeps

         !----------------------------------------------------------------------
         !  
         !----------------------------------------------------------------------
         INTEGER                                       :: stime

         !----------------------------------------------------------------------
         !  I/O and control
         !----------------------------------------------------------------------
         INTEGER                                       :: freqoutput,freqdiag
         CHARACTER(LEN=MAXCHAR)                        :: logfile,ctrlfile
         CHARACTER(LEN=MAXCHAR)                        :: inputfile,abortfile
         CHARACTER(LEN=MAXCHAR)                        :: outputfile
         CHARACTER(LEN=MAXCHAR)                        :: diagfile
         CHARACTER(LEN=MAXCHAR)                        :: prgname,diagfmt
         CHARACTER(LEN=MAXCHAR)                        :: casename
         INTEGER                                       :: debug
         INTEGER                                       :: ppm_log_unit
         INTEGER                                       :: client_log_unit

         !----------------------------------------------------------------------
         !  MPI stuff
         !----------------------------------------------------------------------
         INTEGER                                       :: rank,Nproc,comm
         INTEGER                                       :: MPPREC

         !----------------------------------------------------------------------
         !  Verlet and cell list memory
         !----------------------------------------------------------------------
         ! number of cells in all directions
         INTEGER, DIMENSION(:,:), POINTER              :: Nm
         ! cell list
!          TYPE(ppm_type_ptr_to_clist), DIMENSION(:), POINTER :: clist

      CONTAINS

         SUBROUTINE define_args
           PROCEDURE(string_func) diag_format
           PROCEDURE(string_func) trunc30

           CALL arg_group('Case')

           CALL arg(casename, 'casename',        &
                    ctrl_name    = 'casename',   &
                    default      = 'ppm_client', &
                    help         = 'Name of this case. "ppm_client" is&
                    & the default. Must be between 1 and 30 characters&
                    & long!',                    &
                    default_func = trunc30)

           CALL arg_group('Input')

           CALL arg(inputfile, 'inputfile',          &
                    ctrl_name    = 'inputfile',      &
                    default      = 'ppm_client.tif', &
                    help         = 'Name of data input file.', &
                    default_func = trunc30)

           CALL arg_group('Output')

           CALL arg(outputfile, 'outputfile',        &
                    ctrl_name    = 'outputfile',     &
                    default      = 'ppm_client.out', &
                    help         = 'Name stub of ouput files. <timeste&
                    &p>.out will be appended.',      &
                    default_func = trunc30)

           CALL arg(freqoutput, 'freqoutput', &
                    ctrl_name = 'freqoutput', &
                    default   = 10,           &
                    min       = 1,            &
                    help      = 'Write output file every # iterations.')

           CALL arg_group('Diagnostics')

           CALL arg(diagfile, 'diagfile',             &
                    ctrl_name    = 'diagnosticsfile', &
                    default      = 'ppm_client.diag', &
                    help         = 'Name of file to which diagnostics &
                    &will be APPENDED.',              &
                    default_func = trunc30)

           CALL arg(freqdiag, 'freqdiag',   &
                    ctrl_name = 'freqdiag', &
                    default   = 1,          &
                    min       = 1,          &
                    help      = 'Write diag file every # iterations.')

           CALL arg(diagfmt, 'diagfmt',                  &
                    ctrl_name    = 'diagnostics_format', &
                    default      = 'FORMATTED',          &
                    default_func = diag_format,          &
                    help         = 'Diagnostics file format. One of:\n * &
                    &FORMATTED    (i.e. ASCII text)\n * UNFORMATTED  (&
                    &i.e. binary)')

           CALL arg_group('Program control parameters')

           CALL arg(debug, 'debug',      &
                    flag      = '-d',    &
                    ctrl_name = 'debug', &
                    default   = 0,       &
                    help      = 'PPM Debug Level. 0,1 or 2. Set to 0 f&
                    &or normal run.')

           CALL arg(abortfile, 'abortfile',     &
                    ctrl_name    = 'abortfile', &
                    default      = 'ABORT',     &
                    help         = 'Name of the abort file. At every t&
                    &ime step, the program will check if a file with t&
                    &his name is present in the directory where it was&
                    & started. If such a file is found, the simulation&
                    & terminates gracefully, writing output files.', &
                    default_func = trunc30)

           CALL arg(logfile, 'logfile',              &
                    ctrl_name    = 'logfile',        &
                    default      = 'ppm_client.log', &
                    help         = 'Name of the log file used by ppm_c&
                    &lient',                         &
                    default_func = trunc30)

           CALL arg(checkmap, 'checkmap',   &
                    ctrl_name = 'checkmap', &
                    help      = 'Check the result of the global map? T&
                    &his only makes sense in a parallel run.', &
                    default   = .FALSE.)

           CALL arg(probeproc, 'probeproc',  &
                    ctrl_name = 'probeproc', &
                    help      = 'Probe relative processor speeds and u&
                    &se this info for load balancing? This only makes &
                    &sense in a parallel run.', &
                    default   = .FALSE.)

           CALL arg(dumpparticles, 'dumpparticles', &
                    ctrl_name = 'dumpparticles',    &
                    help      = 'Dump particle positions in ASCII file&
                    &s?\nWARNING: large files will be written!!', &
                    default   = .FALSE.)

           CALL arg_group('Domain')

           CALL arg(min_phys, 'min_phys', &
                    default = (/ 0.0_MK, 0.0_MK, 0.0_MK /))

           CALL arg(max_phys, 'max_phys', &
                    default = (/ 1.0_MK, 1.0_MK, 1.0_MK /))

           IF (MK .EQ. KIND(1.0D0)) THEN
              lmyeps = 1D-10
           ELSE
              lmyeps = 1E-5
           END IF

           ppm_decomp = -1
           ppm_log_unit = 99
           client_log_unit = 88

         END SUBROUTINE define_args

      END MODULE client_global

LOGICAL FUNCTION diag_format(diagfmt)
  CHARACTER(LEN=*), POINTER :: diagfmt
  IF (diagfmt(1:1) .EQ. 'F' .OR. diagfmt(1:1) .EQ. 'f') THEN
     diagfmt = 'FORMATTED'
  ELSE
     diagfmt = 'UNFORMATTED'
  END IF
  diag_format = .TRUE.
END FUNCTION diag_format

LOGICAL FUNCTION trunc30(string)
  CHARACTER(LEN=*), POINTER :: string
  INTEGER                   :: len
  trunc30 = .TRUE.
  len = LEN_TRIM(string)
  IF (len .LT. 1) THEN
     trunc30 = .FALSE.     ! if .false. module ctrl reverts to default
  ELSE IF (len .GT. 30) THEN
     string = string(1:30)
  END IF
END FUNCTION trunc30
