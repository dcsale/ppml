!------------------------------------------------------------------------------
! Subroutine : naga_output.f90
!------------------------------------------------------------------------------
! _ __
! / |/ /__ ____ ____ _
! / / _ `/ _ `/ _ `/
! /_/|_/\_,_/\_, /\_,_/
! /___/
! A PPM Vortex-In-Cell client - 2011 (c)
! Naga is distributed under the GNU Lesser General Public License version 3
! Johannes Tophøj Rasmussen, johannes.tophoj@gmail.com
! Technical University of Denmark
! Department of Mechanical Engineering
!
!
! Sets up, formats and writes diagnostics to file
! is called as
! naga_output(.TRUE./.FALSE.) to intialise (logical toggles header)
! naga_output('description', value) to add a value
! naga_output(fid) to write to file opened with file id fid
! NOTE: No comment character is added to header.
! in gnuplot use: set key autotitle columnhead
!------------------------------------------------------------------------------
MODULE naga_mod_output
USE naga_mod_globals
USE naga_mod_say
IMPLICIT NONE
CHARACTER(len=1024) :: descline
CHARACTER(len=1024) :: dataline
LOGICAL :: line0
INTEGER :: fields
INTEGER,PARAMETER :: flen=20
INTERFACE naga_output
  MODULE PROCEDURE naga_output_clear
  MODULE PROCEDURE naga_output_int
  MODULE PROCEDURE naga_output_real
  MODULE PROCEDURE naga_output_write
END INTERFACE
CONTAINS
!clear output buffers and toggle writing header
SUBROUTINE naga_output_clear(first,info)
IMPLICIT NONE
!-------------------------------------------------------------------------
! Arguments
!-------------------------------------------------------------------------
LOGICAL, INTENT(IN) :: first
INTEGER, INTENT(OUT),OPTIONAL :: info
!-------------------------------------------------------------------------
! Local variables
!-------------------------------------------------------------------------
REAL(MK) :: t0
INTEGER :: info2
!----------------------------------------------------------------------------
! Initialise routine
!----------------------------------------------------------------------------
CALL substart('naga_output',t0,info2)
!----------------------------------------------------------------------------
! Toggle header write and clear buffers
!----------------------------------------------------------------------------
IF (first .EQV. .TRUE.) THEN
  line0 = .TRUE.
  descline = ''
ELSE
  line0 = .FALSE.
END IF
dataline = ''
fields = 0
!----------------------------------------------------------------------------
! Return
!----------------------------------------------------------------------------
 9999 CONTINUE
CALL substop('naga_output',t0,info2)
RETURN
END SUBROUTINE naga_output_clear
!add integer value
SUBROUTINE naga_output_int(desc,value,info)
IMPLICIT NONE
!-------------------------------------------------------------------------
! Arguments
!-------------------------------------------------------------------------
CHARACTER(LEN=*), INTENT(IN) :: desc
INTEGER, INTENT(IN) :: value
INTEGER, INTENT(OUT),OPTIONAL :: info
!-------------------------------------------------------------------------
! Local variables
!-------------------------------------------------------------------------
REAL(MK) :: t0
INTEGER :: info2
!----------------------------------------------------------------------------
! Initialise routine
!----------------------------------------------------------------------------
CALL substart('naga_output',t0,info2)
!----------------------------------------------------------------------------
! Append value and possibly description
!----------------------------------------------------------------------------
IF (line0) THEN
  descline((fields*flen+1):((fields+1)*flen+1)) = desc
END IF
WRITE(dataline((fields*flen+1):((fields+1)*flen+1)),'(I20)') value
fields = fields + 1
!----------------------------------------------------------------------------
! Return
!----------------------------------------------------------------------------
 9999 CONTINUE
CALL substop('naga_output',t0,info2)
RETURN
END SUBROUTINE naga_output_int
!add real value
SUBROUTINE naga_output_real(desc,value,info)
IMPLICIT NONE
!-------------------------------------------------------------------------
! Arguments
!-------------------------------------------------------------------------
CHARACTER(LEN=*), INTENT(IN) :: desc
REAL(MK), INTENT(IN) :: value
INTEGER, INTENT(OUT),OPTIONAL :: info
!-------------------------------------------------------------------------
! Local variables
!-------------------------------------------------------------------------
REAL(MK) :: t0
INTEGER :: info2
!----------------------------------------------------------------------------
! Initialise routine
!----------------------------------------------------------------------------
CALL substart('naga_output',t0,info2)
!----------------------------------------------------------------------------
! Append value and possibly description
!----------------------------------------------------------------------------
IF (line0) THEN
  descline((fields*flen+1):((fields+1)*flen+1)) = desc
END IF
WRITE(dataline((fields*flen+1):((fields+1)*flen+1)),'(20E20.12)') value
fields = fields + 1
!----------------------------------------------------------------------------
! Return
!----------------------------------------------------------------------------
 9999 CONTINUE
CALL substop('naga_output',t0,info2)
RETURN
END SUBROUTINE naga_output_real
!write buffers to file
SUBROUTINE naga_output_write(fid,info)
IMPLICIT NONE
!-------------------------------------------------------------------------
! Arguments
!-------------------------------------------------------------------------
INTEGER, INTENT(IN) :: fid
INTEGER, INTENT(OUT),OPTIONAL :: info
!-------------------------------------------------------------------------
! Local variables
!-------------------------------------------------------------------------
REAL(MK) :: t0
INTEGER :: info2
!----------------------------------------------------------------------------
! Initialise routine
!----------------------------------------------------------------------------
CALL substart('naga_output',t0,info2)
!----------------------------------------------------------------------------
! Comment this back in to add comment character to description line
! (this will modify the description
! of the first variable - therefore always include a leading space in it)
! and write to file
!----------------------------------------------------------------------------
!descline(1:1) = '#'
IF (line0) THEN
  WRITE(fid,*) descline(1:(fields*flen))
ENDIF
WRITE(fid,*) dataline(1:(fields*flen))
!----------------------------------------------------------------------------
! Return
!----------------------------------------------------------------------------
 9999 CONTINUE
CALL substop('naga_output',t0,info2)
RETURN
END SUBROUTINE naga_output_write
END MODULE naga_mod_output
