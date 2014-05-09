!------------------------------------------------------------------------------
! Subroutine : naga_say.f90
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
! This routine outputs the arguments 'rank','caller' and 'mesg' in the format
! [rank](caller) : mesg
! rank is output with leading zeros. If rank is negative it, and, [] is omitted
! if info = -1 only rank 0 outputs
!@ make a similar routine with info being INTENT(IN) so it can be called without variable
!------------------------------------------------------------------------------
MODULE naga_mod_say
IMPLICIT NONE
INTERFACE naga_say
  MODULE PROCEDURE naga_say
END INTERFACE
CONTAINS
SUBROUTINE naga_say(rank,caller,mesg,info)
IMPLICIT NONE
!-------------------------------------------------------------------------
! Arguments
!-------------------------------------------------------------------------
INTEGER, INTENT(IN) :: rank
CHARACTER(LEN=*), INTENT(IN) :: caller,mesg
INTEGER, INTENT(OUT), OPTIONAL :: info
!-------------------------------------------------------------------------
! Local variables
!-------------------------------------------------------------------------
INTEGER :: ios
INTEGER :: icaller,imsg
IF (PRESENT(info)) THEN
  info = 0
ENDIF
!-------------------------------------------------------------------------
! Get length of messages
!-------------------------------------------------------------------------
icaller = LEN_TRIM(caller)
imsg = LEN_TRIM(mesg)
!-------------------------------------------------------------------------
! Define the print format
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
! Do the print
!-------------------------------------------------------------------------
IF (rank.LT.0) THEN !All print (without their id)
  WRITE(*,'(4A)',IOSTAT=ios) '(', caller(1:icaller), ') : ', &
  & mesg(1:imsg)
ELSE
  IF (rank .EQ. 0) THEN !Only root prints
    WRITE(*,'(A,I4.4,4A)',IOSTAT=ios) '[',rank,'](',caller(1:icaller),') : ', &
    & mesg(1:imsg)
    IF (PRESENT(info)) THEN
      info = ios
    ENDIF
  ELSE !All print (with id)
    WRITE(*,'(A,I4.4,4A)',IOSTAT=ios) '[',rank,'](',caller(1:icaller),') : ', &
    & mesg(1:imsg)
  ENDIF
ENDIF
!-------------------------------------------------------------------------
! Return
!-------------------------------------------------------------------------
 9999 CONTINUE
RETURN
END SUBROUTINE naga_say
END MODULE naga_mod_say
