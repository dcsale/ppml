!------------------------------------------------------------------------------
! Subroutine :  naga_time.f90
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
! Times parts of the code. Is called as:
! CALL naga_time_initialise()
! CALL naga_time_str(id)
! CALL naga_time_stp(id)
! CALL naga_time_write()
! the id specifies a specific timer and all when called an id is started and
! stopped repeatedly the cumulative time is recorded.
!
! The calls can be contained within cpp ifdef statements such that timing is
! only performed when switched on at compile time.
!------------------------------------------------------------------------------

MODULE naga_mod_time

USE naga_mod_globals, ONLY: rank,runtag,iruntag,MK,substart,substop
USE naga_mod_say

IMPLICIT NONE

INTEGER,PARAMETER              :: slots=100
INTEGER, DIMENSION(slots)      :: starttime
REAL(MK),DIMENSION(slots)      :: totaltime
INTEGER                        :: count_max


CONTAINS


SUBROUTINE naga_time_initialise(info)

IMPLICIT NONE
!-------------------------------------------------------------------------
! Arguments
!-------------------------------------------------------------------------
INTEGER, INTENT(OUT)          :: info

!-------------------------------------------------------------------------
! Local variables
!-------------------------------------------------------------------------
REAL(MK)                      :: t0
INTEGER                       :: info2


!----------------------------------------------------------------------------
! Initialise routine
!----------------------------------------------------------------------------
CALL substart('naga_time_initialise',t0,info2)
info = 0

!----------------------------------------------------------------------------
! Set arrays to zero and get count_max
!----------------------------------------------------------------------------
starttime = 0
totaltime = 0
CALL SYSTEM_CLOCK(COUNT_MAX=count_max)


!----------------------------------------------------------------------------
! Return
!----------------------------------------------------------------------------
 9999 CONTINUE
CALL substop('naga_time_initialise',t0,info2)
RETURN


END SUBROUTINE naga_time_initialise


SUBROUTINE naga_time_str(id,info)

IMPLICIT NONE
!-------------------------------------------------------------------------
! Arguments
!-------------------------------------------------------------------------
INTEGER, INTENT(IN)           :: id
INTEGER, INTENT(OUT),OPTIONAL :: info

!-------------------------------------------------------------------------
! Local variables
!-------------------------------------------------------------------------
REAL(MK)                      :: t0
INTEGER                       :: info2

!----------------------------------------------------------------------------
! Initialise routine
!----------------------------------------------------------------------------
CALL substart('naga_time_str',t0,info2)

IF (PRESENT(info)) THEN
  info = 0
ENDIF


!----------------------------------------------------------------------------
! Get start time
!----------------------------------------------------------------------------
CALL SYSTEM_CLOCK(COUNT=starttime(id)) ! Start timing


!----------------------------------------------------------------------------
! Return
!----------------------------------------------------------------------------
 9999 CONTINUE
CALL substop('naga_time_str',t0,info2)
RETURN

END SUBROUTINE naga_time_str


SUBROUTINE naga_time_stp(id,info)

IMPLICIT NONE
!-------------------------------------------------------------------------
! Arguments
!-------------------------------------------------------------------------
INTEGER, INTENT(IN)           :: id
INTEGER, INTENT(OUT),OPTIONAL :: info

!-------------------------------------------------------------------------
! Local variables
!-------------------------------------------------------------------------
REAL(MK)                      :: t0
INTEGER                       :: info2
INTEGER                       :: count_end

!----------------------------------------------------------------------------
! Initialise routine
!----------------------------------------------------------------------------
CALL substart('naga_time_stp',t0,info2)

IF (PRESENT(info)) THEN
  info = 0
ENDIF


!----------------------------------------------------------------------------
! Get stop time and calculate and add difference
!----------------------------------------------------------------------------
CALL SYSTEM_CLOCK(COUNT=count_end)
IF (count_end .GE. starttime(id)) THEN
  totaltime(id) = totaltime(id)  + count_end -  starttime(id)
ELSE
  totaltime(id) = totaltime(id)  + count_end - (starttime(id) - count_max)
ENDIF


!----------------------------------------------------------------------------
! Return
!----------------------------------------------------------------------------
 9999 CONTINUE
CALL substop('naga_time_stp',t0,info2)
RETURN


END SUBROUTINE naga_time_stp


SUBROUTINE naga_time_write(info)

IMPLICIT NONE
!-------------------------------------------------------------------------
! Arguments
!-------------------------------------------------------------------------
INTEGER, INTENT(OUT)          :: info

!-------------------------------------------------------------------------
! Local variables
!-------------------------------------------------------------------------
REAL(MK)                      :: t0
INTEGER                       :: info2
INTEGER                       :: i
INTEGER                       :: count_rate
CHARACTER(256)                :: filename

!----------------------------------------------------------------------------
! Initialise routine
!----------------------------------------------------------------------------
CALL substart('naga_time_write',t0,info2)
info = 0

!----------------------------------------------------------------------------
! Open file
!----------------------------------------------------------------------------
WRITE(filename,'(A,A,I5.5,A)') runtag(1:iruntag), 'R', rank, '.timing'
OPEN(21,FILE=filename,IOSTAT=info,POSITION='append',STATUS='replace')
IF (info .NE. 0) THEN
  CALL naga_say(rank, 'naga_time_write','Failed to open output file.')
  GOTO 9999
ENDIF

!----------------------------------------------------------------------------
! Calculate physical time and write to file
!----------------------------------------------------------------------------
CALL SYSTEM_CLOCK(COUNT_RATE=count_rate)
DO i = 1,slots
  WRITE(21,'(I8,20E20.12)') i, REAL(totaltime(i),MK)/REAL(count_rate,MK)
END DO

CLOSE(21)


!----------------------------------------------------------------------------
! Return
!----------------------------------------------------------------------------
 9999 CONTINUE
CALL substop('naga_time_write',t0,info2)
RETURN

END SUBROUTINE naga_time_write


END MODULE naga_mod_time

