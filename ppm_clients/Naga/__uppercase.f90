      !-------------------------------------------------------------------------
      ! Subroutine : UpperCase
      !-------------------------------------------------------------------------
      !
      ! Purpose : This routine converts a character string to upper case
      ! letters.
      !
      ! Input : ilen (I) length of input string
      !
      ! Input/output : string (C) input: string to be converted.
      ! output: upper case version of it
      !
      ! Output : into (I) return status
      !
      ! Routines :
      !
      ! Remarks :
      !
      ! References :
      !
      ! Revisions :
      !-------------------------------------------------------------------------
      ! $Log: UpperCase.f,v $
      ! Revision 1.1.1.1 2007/07/13 10:18:16 ivos
      ! CBL version of the PPM library
      !
      ! Revision 1.1 2004/06/02 15:28:40 ivos
      ! Initial implementation.
      !
      !-------------------------------------------------------------------------
      ! Jens Honore Walther
      ! Institute of Computational Science
      ! ETH Zentrum, Hirschengraben 84
      ! CH-8092 Zurich, Switzerland
      !-------------------------------------------------------------------------
      SUBROUTINE UpperCase(string,ilen,info)
      IMPLICIT NONE
      !-------------------------------------------------------------------------
      ! Arguments
      !-------------------------------------------------------------------------
      CHARACTER(LEN=*), INTENT(INOUT) :: string
      INTEGER , INTENT(IN ) :: ilen
      INTEGER , INTENT( OUT) :: info
      !-------------------------------------------------------------------------
      ! Local variables
      !-------------------------------------------------------------------------
      INTEGER :: i,j
      INTEGER :: i1,i2,i3,iadd
      !-------------------------------------------------------------------------
      ! Externals
      !-------------------------------------------------------------------------
      !-------------------------------------------------------------------------
      ! Initialise
      !-------------------------------------------------------------------------
      info = 0
      !-------------------------------------------------------------------------
      ! Convert to upper case
      !-------------------------------------------------------------------------
      i1 = IACHAR('a') - 1
      i2 = IACHAR('z') + 1
      i3 = IACHAR('A')
      iadd = i3 - i1 - 1
      DO i=1,ilen
         j = IACHAR(string(i:i))
         IF (j.GT.i1.AND.j.LT.i2) THEN
            string(i:i) = CHAR(j+iadd)
         ENDIF
      ENDDO
      !-------------------------------------------------------------------------
      ! Return
      !-------------------------------------------------------------------------
 9999 CONTINUE
      RETURN
      END SUBROUTINE UpperCase
