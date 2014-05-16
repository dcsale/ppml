!------------------------------------------------------------------------------
! Subroutine :  naga_stl_read.f90
!------------------------------------------------------------------------------
!      _  __
!     / |/ /__ ____ ____ _
!    /    / _ `/ _ `/ _ `/
!   /_/|_/\_,_/\_, /\_,_/
!             /___/
!   A PPM Vortex-In-Cell client - 2011 (c)
!   Johannes Tophøj Rasmussen, johannes.tophoj@gmail.com
!   Technical University of Denmark
!   Department of Mechanical Engineering
!
!
! This routines reads STL files and preprocesses them
!------------------------------------------------------------------------------

SUBROUTINE naga_stl_read(info)

  !!USE module_wvic
  !!USE ppm_module_write
  USE naga_mod_globals
  USE naga_mod_say
  !USE naga_mod_stl

  IMPLICIT NONE

  !-------------------------------------------------------------------------
  ! Arguments
  !-------------------------------------------------------------------------
  INTEGER, INTENT(inout)                :: info


  !-------------------------------------------------------------------------
  ! Local variables
  !-------------------------------------------------------------------------
  CHARACTER(LEN=156)       :: msg
  INTEGER                  :: tr
  INTEGER                  :: vertex,state
  !!INTEGER                  :: stlfilelen
  INTEGER                  :: iline
  INTEGER                  :: iunit
  INTEGER                  :: ilen,ios
  CHARACTER(LEN=256)       :: cbuf, tmpbuf
  !!CHARACTER(LEN=256)       :: cvalue,carg
  LOGICAL                  :: stlfileexists
  REAL(MK), DIMENSION(3)   :: vertex1,vertex2,vertex3,vecu,vecv,vecw
  REAL(MK)                 :: t0


  !----------------------------------------------------------------------------
  ! Initialise routine
  !----------------------------------------------------------------------------
  CALL substart('naga_stl_read',t0,info)


  !-----------------------------------------------------------------------------
  ! Definition of file unit
  !-----------------------------------------------------------------------------
  iunit = 30


  !-----------------------------------------------------------------------------
  ! open STL file
  !-----------------------------------------------------------------------------
  stlt(1)%lenfilename = LEN_TRIM(stlopt_filename)
  INQUIRE(FILE=stlopt_filename(1:stlt(1)%lenfilename),EXIST=stlfileexists)
  IF(.NOT. stlfileexists) THEN
    WRITE(msg,'(2A)')'No such STL file: ',stlopt_filename(1:stlt(1)%lenfilename)
    CALL naga_say(rank,'naga_stl_read',msg)
    info = 1
    GOTO 9999
  END IF

  OPEN(iunit,FILE=stlopt_filename(1:stlt(1)%lenfilename),IOSTAT=ios,ACTION='READ')
  IF(ios .NE. 0) THEN
    WRITE(msg,'(2A)')'Failed to open STL file: ',stlopt_filename(1:stlt(1)%lenfilename)
    CALL naga_say(rank,'naga_stl_read',msg)
    info = 1
    GOTO 9999
  END IF

  !-----------------------------------------------------------------------------
  ! Count facets
  !-----------------------------------------------------------------------------
  stlt(1)%tri_count = 0
  DO
    READ(iUnit,'(A)',END=50,ERR=200) cbuf
    cbuf=ADJUSTL(cbuf)
    ilen = LEN_TRIM(cbuf)
    CALL UpperCase(cbuf,5,info)
    IF('FACET' .EQ. cbuf(1:5)) THEN
      stlt(1)%tri_count = stlt(1)%tri_count + 1
    END IF
  END DO

50 REWIND(iUnit)
  !-----------------------------------------------------------------------------
  ! Allocate vector dotproduct arrays
  !-----------------------------------------------------------------------------
  ALLOCATE(stlt(1)%tri_norm(stlt(1)%tri_count,3))
  ALLOCATE(stlt(1)%tri_base(stlt(1)%tri_count,3))
  ALLOCATE(stlt(1)%tri_vecu(stlt(1)%tri_count,3))
  ALLOCATE(stlt(1)%tri_vecv(stlt(1)%tri_count,3))
  ALLOCATE(stlt(1)%tri_vecw(stlt(1)%tri_count,3))
  ALLOCATE(stlt(1)%tri_denom(stlt(1)%tri_count))
  ALLOCATE(stlt(1)%tri_udotv(stlt(1)%tri_count))
  ALLOCATE(stlt(1)%tri_udotu(stlt(1)%tri_count))
  ALLOCATE(stlt(1)%tri_vdotv(stlt(1)%tri_count))
  ALLOCATE(stlt(1)%tri_wdotw(stlt(1)%tri_count))

  ALLOCATE(stlt(1)%tri_udotu2d(stlt(1)%tri_count))
  ALLOCATE(stlt(1)%tri_vdotv2d(stlt(1)%tri_count))
  ALLOCATE(stlt(1)%tri_udotv2d(stlt(1)%tri_count))
  ALLOCATE(stlt(1)%tri_denom2d(stlt(1)%tri_count))
  !-----------------------------------------------------------------------------
  ! scan file
  !-----------------------------------------------------------------------------
  state = 0
  iline = 0
  tr=0
  vertex = 1
  DO
    iline = iline + 1        ! increment line
    READ(iUnit,'(A)',END=100,ERR=200) cbuf
    cbuf=ADJUSTL(cbuf)
    ilen = LEN_TRIM(cbuf)

    !--------------------------------------------------------------------------
    !  Skip comment or empty lines
    !--------------------------------------------------------------------------
    IF(ilen .GT. 0 .AND. cbuf(1:1) .NE. '#') THEN
      CALL UpperCase(cbuf,ilen,info)
      IF('SOLID' .EQ. cbuf(1:5)) THEN

        IF(state .NE. 0) THEN
          WRITE(msg,*)'STL input: misplaced SOLID, line ',iline
          CALL naga_say(rank,'naga_stl_read',msg)
          Info = 1
          GOTO 9999
        ELSE
          state=1
          CYCLE
        END IF
      END IF
      IF('FACET NORMAL' .EQ. cbuf(1:12)) THEN
        IF(state .NE. 1) THEN
          WRITE(msg,*)'STL input: misplaced FACET, line ',iline
          CALL naga_say(rank,'naga_stl_read',msg)
          Info = 1
          GOTO 9999
        ELSE
          tr=tr+1
          tmpbuf=cbuf(13:ilen)
          tmpbuf=ADJUSTL(tmpbuf)
          READ(tmpbuf,*,iostat=ios,err=200) stlt(1)%tri_norm(tr,:)
          state=2
          CYCLE
        END IF
      END IF
      IF('VERTEX' .EQ. cbuf(1:6)) THEN
        IF(state .NE. 2) THEN
          WRITE(msg,*)'STL input: misplaced VERTEX, line ',iline
          CALL naga_say(rank,'naga_stl_read',msg)
          Info = 1
          GOTO 9999
        ELSE IF (vertex .EQ. 1) THEN
          tmpbuf=cbuf(7:ilen)
          tmpbuf=ADJUSTL(tmpbuf)
          READ(tmpbuf,*,iostat=ios,err=200) vertex1
          vertex1 = vertex1 * stlopt_scale
          vertex1 = vertex1 + stlopt_translate
          IF (tr .EQ. 1) THEN
            stl(1)%bndminx = vertex1(1)
            stl(1)%bndmaxx = vertex1(1)
            stl(1)%bndminy = vertex1(2)
            stl(1)%bndmaxy = vertex1(2)
            stl(1)%bndminz = vertex1(3)
            stl(1)%bndmaxz = vertex1(3)
          END IF
          stl(1)%bndminx = MIN(stl(1)%bndminx,vertex1(1))
          stl(1)%bndmaxx = MAX(stl(1)%bndmaxx,vertex1(1))
          stl(1)%bndminy = MIN(stl(1)%bndminy,vertex1(2))
          stl(1)%bndmaxy = MAX(stl(1)%bndmaxy,vertex1(2))
          stl(1)%bndminz = MIN(stl(1)%bndminz,vertex1(3))
          stl(1)%bndmaxz = MAX(stl(1)%bndmaxz,vertex1(3))
          vertex=2
          CYCLE
        ELSE IF (vertex .EQ. 2) THEN
          tmpbuf=cbuf(7:ilen)
          tmpbuf=ADJUSTL(tmpbuf)
          READ(tmpbuf,*,iostat=ios,err=200) vertex2
          vertex2 = vertex2 * stlopt_scale
          vertex2 = vertex2 + stlopt_translate
          stl(1)%bndminx = MIN(stl(1)%bndminx,vertex2(1))
          stl(1)%bndmaxx = MAX(stl(1)%bndmaxx,vertex2(1))
          stl(1)%bndminy = MIN(stl(1)%bndminy,vertex2(2))
          stl(1)%bndmaxy = MAX(stl(1)%bndmaxy,vertex2(2))
          stl(1)%bndminz = MIN(stl(1)%bndminz,vertex2(3))
          stl(1)%bndmaxz = MAX(stl(1)%bndmaxz,vertex2(3))
          vertex=3
          CYCLE
        ELSE IF (vertex .EQ. 3) THEN
          tmpbuf=cbuf(7:ilen)
          tmpbuf=ADJUSTL(tmpbuf)
          READ(tmpbuf,*,iostat=ios,err=200) vertex3
          vertex3 = vertex3 * stlopt_scale
          vertex3 = vertex3 + stlopt_translate
          stl(1)%bndminx = MIN(stl(1)%bndminx,vertex3(1))
          stl(1)%bndmaxx = MAX(stl(1)%bndmaxx,vertex3(1))
          stl(1)%bndminy = MIN(stl(1)%bndminy,vertex3(2))
          stl(1)%bndmaxy = MAX(stl(1)%bndmaxy,vertex3(2))
          stl(1)%bndminz = MIN(stl(1)%bndminz,vertex3(3))
          stl(1)%bndmaxz = MAX(stl(1)%bndmaxz,vertex3(3))
          vertex=1
          stlt(1)%tri_base(tr,:) = vertex1
          vecu = vertex2 - vertex1
          vecv = vertex3 - vertex1
          vecw = vertex3 - vertex2
          stlt(1)%tri_vecu(tr,:) = vecu
          stlt(1)%tri_vecv(tr,:) = vecv
          stlt(1)%tri_vecw(tr,:) = vecw
          stlt(1)%tri_udotu(tr) = (vecu(1)*vecu(1)+vecu(2)*vecu(2)+vecu(3)*vecu(3))
          stlt(1)%tri_vdotv(tr) = (vecv(1)*vecv(1)+vecv(2)*vecv(2)+vecv(3)*vecv(3))
          stlt(1)%tri_wdotw(tr) = (vecw(1)*vecw(1)+vecw(2)*vecw(2)+vecw(3)*vecw(3))
          stlt(1)%tri_udotv(tr) = (vecu(1)*vecv(1)+vecu(2)*vecv(2)+vecu(3)*vecv(3))
          stlt(1)%tri_denom(tr)  = 1.0_MK/ &
                 & (stlt(1)%tri_udotu(tr)*stlt(1)%tri_vdotv(tr)-stlt(1)%tri_udotv(tr)**2)


          IF (stlopt_inout_direction .EQ. 1) THEN
            stlt(1)%tri_udotu2d(tr) = (vecu(2)*vecu(2)+vecu(3)*vecu(3))
            stlt(1)%tri_vdotv2d(tr) = (vecv(2)*vecv(2)+vecv(3)*vecv(3))
            stlt(1)%tri_udotv2d(tr) = (vecu(2)*vecv(2)+vecu(3)*vecv(3))
          ELSE IF (stlopt_inout_direction .EQ. 2) THEN
            stlt(1)%tri_udotu2d(tr) = (vecu(3)*vecu(3)+vecu(1)*vecu(1))
            stlt(1)%tri_vdotv2d(tr) = (vecv(3)*vecv(3)+vecv(1)*vecv(1))
            stlt(1)%tri_udotv2d(tr) = (vecu(3)*vecv(3)+vecu(1)*vecv(1))
          ELSE IF (stlopt_inout_direction .EQ. 3) THEN
            stlt(1)%tri_udotu2d(tr) = (vecu(1)*vecu(1)+vecu(2)*vecu(2))
            stlt(1)%tri_vdotv2d(tr) = (vecv(1)*vecv(1)+vecv(2)*vecv(2))
            stlt(1)%tri_udotv2d(tr) = (vecu(1)*vecv(1)+vecu(2)*vecv(2))
          ENDIF
          stlt(1)%tri_denom2d(tr) = 1.0_MK/ &
           & (stlt(1)%tri_udotu2d(tr)*stlt(1)%tri_vdotv2d(tr)-stlt(1)%tri_udotv2d(tr)**2)
          CYCLE
        END IF
      END IF
      IF('ENDFACET' .EQ. cbuf(1:8)) THEN
        IF(state .NE. 2) THEN
          WRITE(msg,*)'STL input: misplaced ENDFACET, line ',iline
          CALL naga_say(rank,'naga_stl_read',msg)
          Info = 1
          GOTO 9999
        ELSE
          state=1
          CYCLE
        END IF
      END IF
      IF('ENDSOLID' .EQ. cbuf(1:8)) THEN
        IF((state .NE. 1) .AND. (vertex .NE. 1)) THEN
          WRITE(msg,*)'STL input: misplaced ENDSOLID, line ',iline
          CALL naga_say(rank,'naga_stl_read',msg)
          Info = 1
          GOTO 9999
        ELSE
          state=0
          CYCLE
        END IF
      END IF
    END IF

  END DO   




  GOTO 9999



200 CONTINUE
  WRITE(msg,'(A,I5,2A)') 'Error reading line: ',iline, &
                       & ' of file: ',stlopt_filename(1:stlt(1)%lenfilename)
  CALL naga_say(rank,'naga_stl_read',msg)
  ilen = LEN_TRIM(cbuf)
  CALL naga_say(rank,'naga_stl_read',cbuf(1:ilen))
  Info = 1
  GOTO 9999

  !-----------------------------------------------------------------------------
  !  End of file
  !-----------------------------------------------------------------------------
100 Info = 0

  !-----------------------------------------------------------------------------
  !  Close file
  !-----------------------------------------------------------------------------
  CLOSE(iUnit)

  !-----------------------------------------------------------------------------
  !  Do last checks
  !-----------------------------------------------------------------------------
  IF (stlt(1)%tri_count .NE. tr) THEN
    WRITE(msg,*)' Error: Bad triangle count ', stlt(1)%tri_count, ' /= ',tr
    CALL naga_say(rank,'naga_stl_read',msg)
    Info = 1
    GOTO 9999
  END IF
  IF (stlt(1)%tri_count .LT. 4) THEN
    WRITE(msg,*)
    CALL naga_say(rank,'naga_stl_read','Error: At least 4 facets required')
    Info = 1
    GOTO 9999
  END IF

  IF ((stl(1)%bndminx .LT. ptcset(1,1)%min(1)) .OR. &
    & (stl(1)%bndmaxx .GT. ptcset(1,1)%max(1)) .OR. &
    & (stl(1)%bndminy .LT. ptcset(1,1)%min(2)) .OR. &
    & (stl(1)%bndmaxy .GT. ptcset(1,1)%max(2)) .OR. &
    & (stl(1)%bndminz .LT. ptcset(1,1)%min(3)) .OR. &
    & (stl(1)%bndmaxz .GT. ptcset(1,1)%max(3))) THEN
    IF (stlopt_check_bounding) THEN 
      CALL naga_say(rank,'naga_stl_read','Error: STL-object exceeds  computational domain:')
      WRITE(msg,*)' minx: ',stl(1)%bndminx,' miny: ',stl(1)%bndminy,' minz: ',stl(1)%bndminz, ''
      CALL naga_say(rank,'naga_stl_read',msg)
      WRITE(msg,*)' maxx: ',stl(1)%bndmaxx,' maxy: ',stl(1)%bndmaxy,' maxz: ',stl(1)%bndmaxz, ''
      CALL naga_say(rank,'naga_stl_read',msg)
      Info = 1
      GOTO 9999
    ELSE
      CALL naga_say(rank,'naga_stl_read','Warning: Object exceeds domain')
    END IF 
  END IF

  IF (rank .EQ. 0) THEN
    WRITE(msg,*) 'rank:',rank,'Succesfully read', stlt(1)%tri_count, ' facets from ', stlopt_filename(1:stlt(1)%lenfilename)
    CALL naga_say(rank,'naga_stl_read',msg)
    WRITE(msg,*)' minx:',stl(1)%bndminx,'miny:',stl(1)%bndminy, 'minz:',stl(1)%bndminz
    CALL naga_say(rank,'naga_stl_read',msg)
    WRITE(msg,*)' maxx:',stl(1)%bndmaxx,'maxy:',stl(1)%bndmaxy, 'maxz:',stl(1)%bndmaxz
    CALL naga_say(rank,'naga_stl_read',msg)
  END IF

  !----------------------------------------------------------------------------
  ! Return
  !----------------------------------------------------------------------------
 9999 CONTINUE
 CALL substop('naga_stl_read',t0,info)
 RETURN


END SUBROUTINE naga_stl_read


