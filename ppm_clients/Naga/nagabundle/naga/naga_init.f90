!------------------------------------------------------------------------------
! Subroutine :  naga_init.f90
!------------------------------------------------------------------------------
!      _  __
!     / |/ /__ ____ ____ _
!    /    / _ `/ _ `/ _ `/
!   /_/|_/\_,_/\_, /\_,_/
!             /___/
!   A PPM Vortex-In-Cell client - 2011 (c)
!   Naga is distributed under the GNU Lesser General Public License version 3
!   Johannes Tophøj Rasmussen, johannes.tophoj@gmail.com
!   Technical University of Denmark
!   Department of Mechanical Engineering
!
!
! This routines initialises ppm, MPI, reads the inputs, checks, etc
!------------------------------------------------------------------------------

MODULE naga_mod_init

IMPLICIT NONE

INTERFACE naga_init
  MODULE PROCEDURE naga_init
END INTERFACE

CONTAINS

SUBROUTINE naga_init(info)

USE naga_mod_globals
USE naga_mod_say
USE naga_mod_defaults
USE naga_mod_check_setup
USE naga_mod_read_ctrl
USE ppm_module_init
USE ppm_module_io
#ifdef __DOTIMING
USE naga_mod_time
#endif

IMPLICIT NONE

!-------------------------------------------------------------------------
! Arguments
!-------------------------------------------------------------------------
INTEGER, INTENT(OUT) :: info


!-------------------------------------------------------------------------
! Local variables
!-------------------------------------------------------------------------
REAL(MK)             :: t0
INTEGER              :: info2
CHARACTER(LEN=256)   :: msg
CHARACTER(LEN=256)   :: nodename
INTEGER              :: name_len


!----------------------------------------------------------------------------
! Initialise routine
!----------------------------------------------------------------------------
CALL substart('naga_init',t0,info)


!-------------------------------------------------------------------------
! Set default values
!-------------------------------------------------------------------------
CALL naga_defaults(info)


!-------------------------------------------------------------------------
! Read values from Ctrl file
! First check if the control file has been specified on the command line
!-------------------------------------------------------------------------
IF (IARGC() .EQ. 1) THEN
  CALL GETARG(1,ctrlfile)
  WRITE(msg,*) 'Using control file :', ctrlfile(1:32)
  info = -1
  CALL naga_say(rank, 'naga_init',msg,info)
ENDIF
CALL naga_read_ctrl(ctrlfile,info)
IF (info .NE. 0) THEN
  CALL naga_say(-1, 'naga_init','Error reading control file, aborting.')
  GOTO 9999
ENDIF


!-------------------------------------------------------------------------
! Check settings
!-------------------------------------------------------------------------
CALL naga_check_setup(info)
IF (info .NE. 0) THEN
  CALL naga_say(-1, 'naga_init','Erroneous input, aborting.')
  GOTO 9999
ENDIF


!-------------------------------------------------------------------------
! Setup constants
!-------------------------------------------------------------------------
ppm_log_unit  = 99
naga_log_unit = 88


!-------------------------------------------------------------------------
! Start MPI
!-------------------------------------------------------------------------
comm = MPI_COMM_WORLD
CALL MPI_Init(info)
IF (info .NE. 0) THEN
  CALL naga_say(-1, 'naga_init','Failed to initialise MPI, aborting.')
  GOTO 9999
ENDIF
!-------------------------------------------------------------------------
! Get rank
!-------------------------------------------------------------------------
CALL MPI_Comm_Rank(comm,rank,info)
IF (info .NE. 0) THEN
  CALL naga_say(-1, 'naga_init','Failed to get rank, aborting.')
  GOTO 9999
ENDIF
!-------------------------------------------------------------------------
! Get total number of processors
!-------------------------------------------------------------------------
CALL MPI_Comm_Size(comm,Nproc,info)
IF (info .NE. 0) THEN
  CALL naga_say(rank, 'naga_init','Failed to get number of nodes, aborting.')
  GOTO 9999
ENDIF
!-------------------------------------------------------------------------
! Here we could get the name of the node. !But we wont
!-------------------------------------------------------------------------
CALL MPI_GET_PROCESSOR_NAME(nodename, name_len, info)
WRITE(*,*) 'rank: ',rank, ' starting on "',nodename(1:name_len),'"'


!-------------------------------------------------------------------------
! Initialise PPM. The arguments are:
!  - number of dimensions
!  - precision
!  - 10^threshold for what is considered zero
!  - debug level
!  - return state variable
!  - io unit for ppm log file
! Open log file
!-------------------------------------------------------------------------
CALL ppm_init(ndim,MK,-10,comm,debug,info,ppm_log_unit)
IF (info .NE. 0) THEN
  CALL naga_say(rank, 'naga_init','Failed to initialise PPM, aborting.')
  GOTO 9999
ENDIF
CALL ppm_io_set_unit(6,0,ppm_log_unit,info)
IF (info .NE. 0) THEN
  CALL naga_say(rank, 'naga_init','Failed to open log file, aborting.')
  GOTO 9999
ENDIF


!-------------------------------------------------------------------------
! Possibly print the parameters of the simulation to the screen
!-------------------------------------------------------------------------


#ifdef __DOTIMING
!-------------------------------------------------------------------------
!  Initialise timing routine
!-------------------------------------------------------------------------
CALL naga_time_initialise(info)
IF (info .NE. 0) THEN
  CALL naga_say(rank,'naga_init','Error initialising timings.')
  GOTO 9999
ENDIF
#endif


!----------------------------------------------------------------------------
! Return
!----------------------------------------------------------------------------
 9999 CONTINUE
CALL substop('naga_init',t0,info2)
RETURN


END SUBROUTINE naga_init

END MODULE naga_mod_init

