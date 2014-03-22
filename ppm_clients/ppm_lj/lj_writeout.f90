#if __TYPE == __REAL
  SUBROUTINE lj_writeout_real(property,filename,info)
#elif __TYPE == __INTEGER
  SUBROUTINE lj_writeout_integer(property,filename,info)
#elif __TYPE == __RANK
  SUBROUTINE lj_writeout_rank(filename,info)
#endif
!======================================================================!
! write out particle positions and one of their properties in pdb-file
! data is collected from all processes and written out sequentially.
! if no property is given, the rank of the host processor is taken as a
! property
!======================================================================!
  USE lj_module_global, ONLY: xp,ndim,prec,nproc,rank,comm,Npart,  &
       nsteps,it,framerate
  IMPLICIT NONE
  INCLUDE 'mpif.h'
  ! arguments
#if __TYPE == __REAL
  REAL(prec),DIMENSION(:),INTENT(IN)    :: property
  REAL(prec),DIMENSION(:),ALLOCATABLE   :: propp_iproc
#elif __TYPE == __INTEGER
  INTEGER,DIMENSION(:),INTENT(IN)       :: property
  INTEGER,DIMENSION(:),ALLOCATABLE      :: propp_iproc
#endif
  CHARACTER(LEN = 256),INTENT(IN)       :: filename
  INTEGER,INTENT(OUT)                   :: info
  ! local variables
  INTEGER,DIMENSION(:),ALLOCATABLE      :: Npart_vec
  REAL(prec),DIMENSION(:,:),ALLOCATABLE :: xp_iproc
  INTEGER                               :: ip,count,iproc,i
  INTEGER,DIMENSION(MPI_STATUS_SIZE)    :: status
  CHARACTER                             :: altloc
  CHARACTER (LEN =  4 )                 :: atom_name
  CHARACTER                             :: chains
  CHARACTER (LEN =  2 )                 :: charge
  CHARACTER (LEN =  2 )                 :: element
  CHARACTER                             :: icode
  INTEGER                               :: resno
  REAL(8)                               :: occ
  CHARACTER (LEN =  3 )                 :: resname
  CHARACTER (LEN =  4 )                 :: segid

  !====================================================================!
  ! initialize
  info = 0
  atom_name = 'H   '
  altloc = ' '
  resname = '   '
  chains = ' '
  resno = 0
  icode = ' '
  occ = 1.0D+00
  segid = '    '
  element = '  '
  charge = '  '

  !====================================================================!
  ! open file
  IF (rank .EQ. 0) THEN
     IF (it .EQ. 1) THEN
        OPEN(23,FILE=filename,FORM='FORMATTED',STATUS='REPLACE',IOSTAT=info)
        IF (info .NE. 0) THEN
           WRITE(*,*)'ERROR in lj_writeout: opening file failed.'
           GOTO 9999
        ENDIF
     ENDIF
     WRITE (23,'(a6,i5)')'MODEL ',(it-1)/framerate + 1
  ENDIF

  !====================================================================!
  ! rank 0 collects data and writes out sequentially
  IF (rank .EQ. 0) THEN
     ALLOCATE(Npart_vec(nproc))
     IF (info .NE. 0) THEN
        WRITE(*,*)'ERROR in lj_writeout: Allocation failed.'
        GOTO 9999
     ENDIF
  ENDIF

  CALL MPI_Gather(Npart,1,MPI_INTEGER,Npart_vec,1,MPI_INTEGER,0,comm,info)
  IF (info .NE. 0) THEN
     WRITE(*,*)'ERROR in lj_writeout: MPI_Gather failed.'
     GOTO 9999
  ENDIF

  IF (rank .EQ. 0) THEN

     ip = 0
     DO iproc = 1,nproc

        IF (iproc .GT. 1) THEN

           ALLOCATE(xp_iproc(ndim,Npart_vec(iproc)))
           IF (info .NE. 0) THEN
              WRITE(*,*)'ERROR in lj_writeout: Allocation failed.'
              GOTO 9999
           ENDIF

#if __TYPE == __RANK

           count = ndim*Npart_vec(iproc)
           IF (prec .EQ. 4) THEN
              CALL MPI_Recv(xp_iproc,count,MPI_REAL,iproc-1,iproc-1,comm,status,info)
           ELSE
              CALL MPI_Recv(xp_iproc,count,MPI_DOUBLE_PRECISION,  &
                   iproc-1,iproc-1,comm,status,info)
           ENDIF

           IF (ndim .EQ. 2) THEN
           
              DO i=1,Npart_vec(iproc)
                 ip = ip + 1
                 WRITE (23, &
                      '(a6,i5,1x,a4,a1,a3,a1,1x,i4,a1,3x,3f8.3,2f6.2,6x,a4,a2,a2)' ) &
                      'ATOM  ', ip, atom_name, altloc, resname, chains, &
                      resno,icode, xp_iproc(1,i),xp_iproc(2,i),0., occ, REAL(iproc-1), segid,  &
                      element, charge        
              END DO

           ELSE
              
              DO i=1,Npart_vec(iproc)
                 ip = ip + 1
                 WRITE (23, &
                      '(a6,i5,1x,a4,a1,a3,a1,1x,i4,a1,3x,3f8.3,2f6.2,6x,a4,a2,a2)' ) &
                      'ATOM  ', ip, atom_name, altloc, resname, chains, &
                      resno,icode, xp_iproc(1,i),xp_iproc(2,i),xp_iproc(3,i), occ, REAL(iproc-1), segid,  &
                      element, charge        
              END DO

           ENDIF

           DEALLOCATE(xp_iproc)
#else
           ALLOCATE(propp_iproc(Npart_vec(iproc)))
           IF (info .NE. 0) THEN
              WRITE(*,*)'ERROR in lj_writeout: Allocation failed.'
              GOTO 9999
           ENDIF

           IF (prec .EQ. 4) THEN
              count = ndim*Npart_vec(iproc)
              CALL MPI_Recv(xp_iproc,count,MPI_REAL,iproc-1,iproc-1,comm,status,info)
#if __TYPE == __REAL
              count = Npart_vec(iproc)
              CALL MPI_Recv(propp_iproc,count,MPI_REAL,iproc-1,iproc-1,comm,status,info)
#endif
           ELSE
              count = ndim*Npart_vec(iproc)
              CALL MPI_Recv(xp_iproc,count,MPI_DOUBLE_PRECISION,  &
                   iproc-1,iproc-1,comm,status,info)
#if __TYPE == __REAL
              count = Npart_vec(iproc)
              CALL MPI_Recv(propp_iproc,count,MPI_DOUBLE_PRECISION,  &
                   iproc-1,iproc-1,comm,status,info)
#endif
           ENDIF

#if __TYPE == __INTEGER
              count = Npart_vec(iproc)
              CALL MPI_Recv(propp_iproc,count,MPI_INTEGER,  &
                   iproc-1,iproc-1,comm,status,info)
#endif

           IF (ndim .EQ. 2) THEN
           
#if __TYPE == __REAL
              DO i=1,Npart_vec(iproc)
                 ip = ip + 1
                 WRITE (23, &
                      '(a6,i5,1x,a4,a1,a3,a1,1x,i4,a1,3x,3f8.3,2f6.2,6x,a4,a2,a2)' ) &
                      'ATOM  ', ip, atom_name, altloc, resname, chains, &
                      resno,icode, xp_iproc(1,i),xp_iproc(2,i),0., occ, propp_iproc(i), segid,  &
                      element, charge        
              END DO
#else
              DO i=1,Npart_vec(iproc)
                 ip = ip + 1
                 WRITE (23, &
                      '(a6,i5,1x,a4,a1,a3,a1,1x,i4,a1,3x,3f8.3,2f6.2,6x,a4,a2,a2)' ) &
                      'ATOM  ', ip, atom_name, altloc, resname, chains, &
                      resno,icode, xp_iproc(1,i),xp_iproc(2,i),0., occ, REAL(propp_iproc(i)),  &
                      segid, element, charge        
              END DO
#endif
              
           ELSE
              
#if __TYPE == __REAL
              DO i=1,Npart_vec(iproc)
                 ip = ip + 1
                 WRITE (23, &
                      '(a6,i5,1x,a4,a1,a3,a1,1x,i4,a1,3x,3f8.3,2f6.2,6x,a4,a2,a2)' ) &
                      'ATOM  ', ip, atom_name, altloc, resname, chains, &
                      resno,icode, xp_iproc(1,i),xp_iproc(2,i),xp_iproc(3,i), occ, propp_iproc(i), segid,  &
                      element, charge        
              END DO
#else
              DO i=1,Npart_vec(iproc)
                 ip = ip + 1
                 WRITE (23, &
                      '(a6,i5,1x,a4,a1,a3,a1,1x,i4,a1,3x,3f8.3,2f6.2,6x,a4,a2,a2)' ) &
                      'ATOM  ', ip, atom_name, altloc, resname, chains, &
                      resno,icode, xp_iproc(1,i),xp_iproc(2,i),xp_iproc(3,i), occ, REAL(propp_iproc(i)),  &
                      segid, element, charge        
              END DO
#endif
              
           ENDIF

           DEALLOCATE(propp_iproc)
           DEALLOCATE(xp_iproc)
           
#endif

        ELSE

#if __TYPE == __RANK

           IF (ndim .EQ. 2) THEN
           
              DO i=1,Npart_vec(iproc)
                 ip = ip + 1
                 WRITE (23, &
                      '(a6,i5,1x,a4,a1,a3,a1,1x,i4,a1,3x,3f8.3,2f6.2,6x,a4,a2,a2)' ) &
                      'ATOM  ', ip, atom_name, altloc, resname, chains, &
                      resno,icode, xp(1,i),xp(2,i),0., occ, 0., segid,  &
                      element, charge        
              END DO

           ELSE
              
              DO i=1,Npart_vec(iproc)
                 ip = ip + 1
                 WRITE (23, &
                      '(a6,i5,1x,a4,a1,a3,a1,1x,i4,a1,3x,3f8.3,2f6.2,6x,a4,a2,a2)' ) &
                      'ATOM  ', ip, atom_name, altloc, resname, chains, &
                      resno,icode, xp(1,i),xp(2,i),xp(3,i), occ, 0., segid,  &
                      element, charge        
              END DO

           ENDIF

#else

           IF (ndim .EQ. 2) THEN
           
#if __TYPE == __REAL
              DO i=1,Npart_vec(iproc)
                 ip = ip + 1
                 WRITE (23, &
                      '(a6,i5,1x,a4,a1,a3,a1,1x,i4,a1,3x,3f8.3,2f6.2,6x,a4,a2,a2)' ) &
                      'ATOM  ', ip, atom_name, altloc, resname, chains, &
                      resno,icode, xp(1,i),xp(2,i),0., occ, property(i), segid,  &
                      element, charge        
              END DO
#else
              DO i=1,Npart_vec(iproc)
                 ip = ip + 1
                 WRITE (23, &
                      '(a6,i5,1x,a4,a1,a3,a1,1x,i4,a1,3x,3f8.3,2f6.2,6x,a4,a2,a2)' ) &
                      'ATOM  ', ip, atom_name, altloc, resname, chains, &
                      resno,icode, xp(1,i),xp(2,i),0., occ, REAL(property(i)),  &
                      segid, element, charge        
              END DO
#endif
              
           ELSE
              
#if __TYPE == __REAL
              DO i=1,Npart_vec(iproc)
                 ip = ip + 1
                 WRITE (23, &
                      '(a6,i5,1x,a4,a1,a3,a1,1x,i4,a1,3x,3f8.3,2f6.2,6x,a4,a2,a2)' ) &
                      'ATOM  ', ip, atom_name, altloc, resname, chains, &
                      resno,icode, xp(1,i),xp(2,i),xp(3,i), occ, property(i), segid,  &
                      element, charge        
              END DO
#else
              DO i=1,Npart_vec(iproc)
                 ip = ip + 1
                 WRITE (23, &
                      '(a6,i5,1x,a4,a1,a3,a1,1x,i4,a1,3x,3f8.3,2f6.2,6x,a4,a2,a2)' ) &
                      'ATOM  ', ip, atom_name, altloc, resname, chains, &
                      resno,icode, xp(1,i),xp(2,i),xp(3,i), occ, REAL(property(i)),  &
                      segid, element, charge        
              END DO
#endif
              
           ENDIF
           
#endif

        ENDIF

     ENDDO

  !====================================================================!
  ! terminate pdb-file and close file

     IF (rank .EQ. 0) THEN
        WRITE (23,'(a6)')'ENDMDL'
        IF (it .GE. nsteps-framerate+1) THEN
           CLOSE(23,IOSTAT=info)
           IF (info .NE. 0) THEN
              WRITE(*,*)'ERROR in lj_writeout: closing file failed.'
              GOTO 9999
           ENDIF
        ENDIF
     ENDIF

     DEALLOCATE(Npart_vec)
     
  ELSE

     count = ndim*Npart
     IF (prec .EQ. 4) THEN
        CALL MPI_Send(xp,count,MPI_REAL,0,rank,comm,info)
        IF (info .NE. 0) THEN
           WRITE(*,*)'ERROR in lj_writeout: MPI_Send failed.'
           GOTO 9999
        ENDIF
     ELSE
        CALL MPI_Send(xp,count,MPI_DOUBLE_PRECISION,0,rank,comm,info)
        IF (info .NE. 0) THEN
           WRITE(*,*)'ERROR in lj_writeout: MPI_Send failed.'
           GOTO 9999
        ENDIF
     ENDIF

     count = Npart
#if __TYPE == __REAL
     IF (prec .EQ. 4) THEN
        CALL MPI_Send(property,count,MPI_REAL,0,rank,comm,info)
        IF (info .NE. 0) THEN
           WRITE(*,*)'ERROR in lj_writeout: MPI_Send failed.'
           GOTO 9999
        ENDIF
     ELSE
        CALL MPI_Send(property,count,MPI_DOUBLE_PRECISION,0,rank,comm,info)
        IF (info .NE. 0) THEN
           WRITE(*,*)'ERROR in lj_writeout: MPI_Send failed.'
           GOTO 9999
        ENDIF
     ENDIF
#elif __TYPE == __INTEGER
        CALL MPI_Send(property,count,MPI_INTEGER,0,rank,comm,info)
        IF (info .NE. 0) THEN
           WRITE(*,*)'ERROR in lj_writeout: MPI_Send failed.'
           GOTO 9999
        ENDIF
#endif

  ENDIF

9999 CONTINUE ! jump here upon error

#if __TYPE == __REAL
  END SUBROUTINE lj_writeout_real
#elif __TYPE == __INTEGER
  END SUBROUTINE lj_writeout_integer
#else
  END SUBROUTINE lj_writeout_rank
#endif

