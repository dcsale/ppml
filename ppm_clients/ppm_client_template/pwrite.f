      !-------------------------------------------------------------------------
      !  Subroutine   :                     pwrite
      !-------------------------------------------------------------------------
      !
      !  Purpose      : Subroutine for parallel writing for debugging purposes.
      !
      !  Input        : rank      (I) MPI rank of writing processor
      !                 caller    (C) name of calling subroutine
      !                 cbuf      (C) message
      !
      !  Input/output : 
      !
      !  Output       : info      (I) error status
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

      SUBROUTINE pwrite(rank,caller,cbuf,info)

      IMPLICIT NONE
      !-------------------------------------------------------------------------
      !  Arguments     
      !-------------------------------------------------------------------------
      INTEGER          :: rank
      CHARACTER(LEN=*) :: caller,cbuf
      INTEGER          :: info
      !-------------------------------------------------------------------------
      !  Local variables 
      !-------------------------------------------------------------------------
      CHARACTER(LEN=256) :: cformat
      INTEGER            :: ios
      INTEGER            :: icaller,ibuf
      !-------------------------------------------------------------------------
      !  Externals 
      !-------------------------------------------------------------------------
      
      !-------------------------------------------------------------------------
      !  Initialise 
      !-------------------------------------------------------------------------

      !-------------------------------------------------------------------------
      !  Get length of messages
      !-------------------------------------------------------------------------
      icaller = LEN_TRIM(caller)
      ibuf    = LEN_TRIM(cbuf)

      !-------------------------------------------------------------------------
      !  Define the print format
      !-------------------------------------------------------------------------
      IF     (rank.LT.0) THEN
         cformat = '(4A)'
      ELSEIF (rank.LT.10) THEN
         cformat = '(A,I1,4A)' 
      ELSEIF (rank.LT.100) THEN
         cformat = '(A,I2,4A)' 
      ELSEIF (rank.LT.1000) THEN
         cformat = '(A,I3,4A)' 
      ELSE
         cformat = '(A,I4,4A)' 
      ENDIF 

      !-------------------------------------------------------------------------
      !  Do the print
      !-------------------------------------------------------------------------
      IF (rank.LT.0) THEN
         WRITE(*,cformat,IOSTAT=ios)                                    &
     &      '(',                                                        &
     &      caller(1:icaller),                                          &
     &      ') : ',                                                     &
     &      cbuf(1:ibuf)
      ELSE
         WRITE(*,cformat,IOSTAT=ios)                                    &
     &      '[',rank,'] (',                                             &
     &      caller(1:icaller),                                          &
     &      ') : ',                                                     &
     &      cbuf(1:ibuf)
      ENDIF 

      !-------------------------------------------------------------------------
      !  Return 
      !-------------------------------------------------------------------------
 9999 CONTINUE
      RETURN
      END SUBROUTINE pwrite
