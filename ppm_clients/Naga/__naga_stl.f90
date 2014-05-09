!------------------------------------------------------------------------------
! Subroutine : naga_stl.f90
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
! This module contains routines for initialising solids from STL input
!------------------------------------------------------------------------------
MODULE naga_mod_stl
USE naga_mod_globals
IMPLICIT NONE
INTERFACE naga_stl_read
  MODULE PROCEDURE naga_stl_read
END INTERFACE
INTERFACE naga_stl_init
  MODULE PROCEDURE naga_stl_init
END INTERFACE
CONTAINS
!------------------------------------------------------------------------------
! Subroutine : naga_stl_read.f90
!------------------------------------------------------------------------------
! _ __
! / |/ /__ ____ ____ _
! / / _ `/ _ `/ _ `/
! /_/|_/\_,_/\_, /\_,_/
! /___/
! A PPM Vortex-In-Cell client - 2011 (c)
! Johannes Tophøj Rasmussen, johannes.tophoj@gmail.com
! Technical University of Denmark
! Department of Mechanical Engineering
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
  INTEGER, INTENT(inout) :: info
  !-------------------------------------------------------------------------
  ! Local variables
  !-------------------------------------------------------------------------
  CHARACTER(LEN=156) :: msg
  INTEGER :: tr
  INTEGER :: vertex,state
  !!INTEGER :: stlfilelen
  INTEGER :: iline
  INTEGER :: iunit
  INTEGER :: ilen,ios
  CHARACTER(LEN=256) :: cbuf, tmpbuf
  !!CHARACTER(LEN=256) :: cvalue,carg
  LOGICAL :: stlfileexists
  REAL(MK), DIMENSION(3) :: vertex1,vertex2,vertex3,vecu,vecv,vecw
  REAL(MK) :: t0
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
    iline = iline + 1 ! increment line
    READ(iUnit,'(A)',END=100,ERR=200) cbuf
    cbuf=ADJUSTL(cbuf)
    ilen = LEN_TRIM(cbuf)
    !--------------------------------------------------------------------------
    ! Skip comment or empty lines
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
          stlt(1)%tri_denom(tr) = 1.0_MK/ &
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
  ! End of file
  !-----------------------------------------------------------------------------
100 Info = 0
  !-----------------------------------------------------------------------------
  ! Close file
  !-----------------------------------------------------------------------------
  CLOSE(iUnit)
  !-----------------------------------------------------------------------------
  ! Do last checks
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
!------------------------------------------------------------------------------
! Subroutine : naga_stl_init.f90
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
! This routines initialises the solid mask from triangulated surface STL files.
! Calculates intersections between (a line from the point to a
! reference point) and (the triangles). Calculates barycentric coordinates and
! determines whether the line intersects the triangle. If so the direction of
! interesection is stored and after checking all triangles it can be determined
! if the point is interior or exterior.
!------------------------------------------------------------------------------
SUBROUTINE naga_stl_init(info)
  USE naga_mod_globals
  USE naga_mod_stepfunction
  USE naga_mod_say
  IMPLICIT NONE
  INTEGER :: info
  INTEGER :: i,j,k,tr,isub,isubl
  INTEGER :: inout,intersections
  REAL(Mk), DIMENSION(3) :: vecT,vecP,vecS,vece
  REAL(Mk), DIMENSION(3) :: vecX
  REAL(Mk) :: Pdotu,Pdotv,a,b,c
  REAL(Mk) :: epsilon
  INTEGER :: ibmin,ibmax,jbmin,jbmax,kbmin,kbmax
  INTEGER :: imin,imax,jmin,jmax,kmin,kmax
  REAL(Mk) :: minimumdist,diag
  INTEGER :: npenextra
  REAL(Mk) :: udotx,vdotx,dist
  INTEGER :: ilevel,ipatch
  REAL(Mk) :: dx,dy,dz
  info = 0
  DO ilevel=1,nlevels
    DO ipatch=1,npatches(ilevel)
      dx = ptcset(ilevel,ipatch)%dx(1)
      dy = ptcset(ilevel,ipatch)%dx(2)
      dz = ptcset(ilevel,ipatch)%dx(3)
      epsilon = 1.0_MK / (sqrt(dx**2 + dy**2 + dz**2))
      chif(ilevel,ipatch)%fld = 0.0_MK
      diag = (ptcset(ilevel,ipatch)%max(1) - ptcset(ilevel,ipatch)%min(1))**2+ &
           & (ptcset(ilevel,ipatch)%max(2) - ptcset(ilevel,ipatch)%min(2))**2+ &
           & (ptcset(ilevel,ipatch)%max(3) - ptcset(ilevel,ipatch)%min(3))**2
      DO isub=1,topos(ilevel)%nsublist
        isubl=topos(ilevel)%isublist(isub)
        !---------------------------------------------------------------------------
        ! determine boundaries for the cells to be evaluated
        ! the boundaries are expanded by the number of cells equal to 50% of
        ! step1_interval cell diagonals
        ! ijkbmin will only exceed ijkbmax when a part of the CV is in the subdomain
        ! so that DO loops only initiate when CVs are present in the subdomain:
        !---------------------------------------------------------------------------
        npenextra = ceiling((step1_interval*0.5_Mk + step1_offset)/epsilon)
        ibmin = NINT((stl(1)%bndminx-ptcset(1,1)%min(1))/dx)+1 - &
              & (ptcset(ilevel,ipatch)%sistr(1,isubl)-1)-npenextra
        jbmin = NINT((stl(1)%bndminy-ptcset(1,1)%min(2))/dy)+1 - &
              & (ptcset(ilevel,ipatch)%sistr(2,isubl)-1)-npenextra
        kbmin = NINT((stl(1)%bndminz-ptcset(1,1)%min(3))/dz)+1 - &
              & (ptcset(ilevel,ipatch)%sistr(3,isubl)-1)-npenextra
        ibmax = NINT((stl(1)%bndmaxx-ptcset(1,1)%min(1))/dx)+1 - &
              & (ptcset(ilevel,ipatch)%sistr(1,isubl)-1)+npenextra
        jbmax = NINT((stl(1)%bndmaxy-ptcset(1,1)%min(2))/dy)+1 - &
              & (ptcset(ilevel,ipatch)%sistr(2,isubl)-1)+npenextra
        kbmax = NINT((stl(1)%bndmaxz-ptcset(1,1)%min(3))/dz)+1 - &
              & (ptcset(ilevel,ipatch)%sistr(3,isubl)-1)+npenextra
        imin=max(ibmin,1-gstw(1))
        jmin=max(jbmin,1-gstw(2))
        kmin=max(kbmin,1-gstw(3))
        imax=min(ibmax,ptcset(ilevel,ipatch)%snx(1,isubl)+gstw(1))
        jmax=min(jbmax,ptcset(ilevel,ipatch)%snx(2,isubl)+gstw(2))
        kmax=min(kbmax,ptcset(ilevel,ipatch)%snx(3,isubl)+gstw(3))
        IF (stlopt_inout_direction .EQ. 1) THEN
    DO k=kmin,kmax
      DO j=jmin,jmax
        DO i=imin,imax
          !---------------------------------------------------------------------
          ! First check if the point is inside or outside the solid
          !---------------------------------------------------------------------
          vecT(1) = REAL(i-1+(ptcset(ilevel,ipatch)%sistr(1,isubl)-1),mk)*dx+ptcset(1,1)%min(1)
          vecT(2) = REAL(j-1+(ptcset(ilevel,ipatch)%sistr(2,isubl)-1),mk)*dy+ptcset(1,1)%min(2)
          vecT(3) = REAL(k-1+(ptcset(ilevel,ipatch)%sistr(3,isubl)-1),mk)*dz+ptcset(1,1)%min(3)
          intersections = 0
          inout = 0
          DO tr=1,stlt(1)%tri_count
            !check if triangle is aligned with the point in the bounds_direction
            !perhaps add elseif verbose a more thorough test
            ! maybe change to some finite value
            !implement warning if stl points are too close to grid coordinates
            IF (stlt(1)%tri_norm(tr,1).NE.0.0_mk) THEN
              vecX(1) = vecT(2) - stlt(1)%tri_base(tr,2)
              vecX(2) = vecT(3) - stlt(1)%tri_base(tr,3)
              vecX(3) = vecT(1) - stlt(1)%tri_base(tr,1)
              !beware of vecX index
              udotx = stlt(1)%tri_vecu(tr,2)*vecX(1) + &
                    & stlt(1)%tri_vecu(tr,3)*vecX(2) !@+ &
                    !@& stlt(1)%tri_vecu(tr,1)*vecX(3)
              vdotx = stlt(1)%tri_vecv(tr,2)*vecX(1) + &
                    & stlt(1)%tri_vecv(tr,3)*vecX(2) !@+ &
                    !@& stlt(1)%tri_vecv(tr,1)*vecX(3)
              !get barycentric coordinates
              a = (udotx*stlt(1)%tri_vdotv2d(tr)-vdotx*stlt(1)%tri_udotv2d(tr))*stlt(1)%tri_denom2d(tr)
              b = (vdotx*stlt(1)%tri_udotu2d(tr)-udotx*stlt(1)%tri_udotv2d(tr))*stlt(1)%tri_denom2d(tr)
              c = (1.0_mk-a-b)
              IF (a.GE.0.0_mk .AND. b.GE.0.0_mk.AND.c.GE.0.0_mk) THEN
                dist = a*(stlt(1)%tri_base(tr,1)+stlt(1)%tri_vecu(tr,1)) + &
                     & b*(stlt(1)%tri_base(tr,1)+stlt(1)%tri_vecv(tr,1)) + &
                     & c* stlt(1)%tri_base(tr,1)
                IF (dist .LT. vecT(1)) THEN
                  intersections = intersections + 1
                  IF (stlopt_check_intersections) THEN
                    IF (stlt(1)%tri_norm(tr,1).GT.0.0_mk) THEN
                      inout = inout + 1
                    ELSE
                      inout = inout - 1
                    ENDIF
                  ENDIF
                ENDIF
              ENDIF
            ENDIF
          ENDDO !tr
          minimumdist = diag
          DO tr=1,stlt(1)%tri_count
            vecS=stlt(1)%tri_base(tr,:)-vecT
            !----------------------------------------------------------------
            ! x3 x0 is the projection of T onto the
            ! /\ triangle plane
            ! x0 / \ x1 is used as triangle base point
            ! / \ vecU = x2 - x1
            ! / \ vecV = x3 - x1
            ! /________\ vecW = x3 - x2
            ! x1 x2
            !
            ! x0 x1 T is the point to be evaluated
            ! - - - ======== - - - === is the triangle (in its plane)
            ! | / vecP = x0 - x1
            !vecP | / vecS vecS = x1 - T
            ! | /
            ! |/
            ! T
            !
            ! Meassure distance to triangle.
            ! This is not the vecP from figure 6.2 in thesis
            ! a and b are being reused. Initial definition is in thesis.
            ! tri_norm contains normal vector of triangle plane
            !----------------------------------------------------------------
            vecP = stlt(1)%tri_norm(tr,:) * &
            & ( stlt(1)%tri_norm(tr,1)*vecS(1) + stlt(1)%tri_norm(tr,2)*vecS(2) + &
            & stlt(1)%tri_norm(tr,3)*vecS(3)) - vecS
            Pdotu = vecP(1)*stlt(1)%tri_vecu(tr,1) + vecP(2)*stlt(1)%tri_vecu(tr,2) + &
            & vecP(3)*stlt(1)%tri_vecu(tr,3)
            Pdotv = vecP(1)*stlt(1)%tri_vecv(tr,1) + vecP(2)*stlt(1)%tri_vecv(tr,2) + &
            & vecP(3)*stlt(1)%tri_vecv(tr,3)
            a = (Pdotu*stlt(1)%tri_vdotv(tr) - Pdotv*stlt(1)%tri_udotv(tr)) * stlt(1)%tri_denom(tr)
            b = (Pdotv*stlt(1)%tri_udotu(tr) - Pdotu*stlt(1)%tri_udotv(tr)) * stlt(1)%tri_denom(tr)
            IF ((a .GE. 0.0_mk) .AND. (b .GE. 0.0_mk) .AND. &
            & ((a+b) .LE. 1.0_mk)) THEN
              !get distance normal to triangle
              a = (stlt(1)%tri_norm(tr,1)*vecS(1) + stlt(1)%tri_norm(tr,2)*vecS(2) + &
              & stlt(1)%tri_norm(tr,3)*vecS(3))**2
            ELSE
              !get distance to edge or corner
              a = vecP(1)*stlt(1)%tri_vecu(tr,1) + vecP(2)*stlt(1)%tri_vecu(tr,2) + &
                & vecP(3)*stlt(1)%tri_vecu(tr,3)
              b = vecP(1)*stlt(1)%tri_vecv(tr,1) + vecP(2)*stlt(1)%tri_vecv(tr,2) + &
                & vecP(3)*stlt(1)%tri_vecv(tr,3)
              c = (vecP(1)-stlt(1)%tri_vecu(tr,1))*stlt(1)%tri_vecw(tr,1) + &
                & (vecP(2)-stlt(1)%tri_vecu(tr,2))*stlt(1)%tri_vecw(tr,2) + &
                & (vecP(3)-stlt(1)%tri_vecu(tr,3))*stlt(1)%tri_vecw(tr,3)
              IF (a .GT. stlt(1)%tri_udotu(tr)) THEN !point is closest to x2
                !a = tri_udotu(tr)
                vece=-vecS-stlt(1)%tri_vecu(tr,:)
                a = SUM(vece**2)
              ELSEIF (a .LT. 0.0_mk) THEN !point is closest to x1
                !a = 0
                !vece=-vecS
                a = SUM(vecS**2)
              ELSE !point is closest to vecU
                vece=-vecS-stlt(1)%tri_vecu(tr,:)*a/stlt(1)%tri_udotu(tr)
                a = SUM(vece**2)
                !this may not be the cheapest way to do this
              ENDIF
              IF (b .GT. stlt(1)%tri_vdotv(tr)) THEN !point is closest to x3
                !b = tri_vdotv(tr)
                vece=-vecS-stlt(1)%tri_vecv(tr,:)
                b = SUM(vece**2)
              ELSEIF (b .LT. 0.0_mk) THEN !point is closest to x1
                !b = 0
                !vece=-vecS
                b = SUM(vecS**2)
              ELSE !point is closest to vecV
                vece=-vecS-stlt(1)%tri_vecv(tr,:)*b/stlt(1)%tri_vdotv(tr)
                b = SUM(vece**2)
              ENDIF
              IF (c .GT. stlt(1)%tri_wdotw(tr)) THEN !point is closest to x3
                !c = tri_wdotw(tr)
                vece=-vecS-stlt(1)%tri_vecv(tr,:)
                c = SUM(vece**2)
              ELSEIF (c .LT. 0.0_mk) THEN !point is closest to x2
                !c = 0
                vece=-vecS-stlt(1)%tri_vecu(tr,:)
                c = SUM(vecS**2)
              ELSE !point is closest to vecW
                vece=-vecS-stlt(1)%tri_vecu(tr,:)-stlt(1)%tri_vecw(tr,:)*c/stlt(1)%tri_wdotw(tr)
                c = SUM(vece**2)
              ENDIF
              a = min(a,b)
              a = min(a,c)
            ENDIF
            minimumdist = min(minimumdist,a)
          END DO !tr triangles
          !!chif(ilevel,ipatch)%fld(i,j,k,isub) = naga_stepfunc1(minimumdist*epsilon)
          minimumdist = sqrt(minimumdist)!JTR cleanup
          IF (stlopt_check_intersections) THEN
            IF (inout.LT.0) THEN
              minimumdist = sign(minimumdist,-1.0_mk)
            ENDIF
          ELSE IF(MOD(intersections,2) .NE. 0) THEN
            minimumdist = sign(minimumdist,-1.0_mk)
          ENDIF
          !The multiplication with epsilon can be done to the step function limits instead
          !Could be done either in this routine or when the limits are determined (better).
          !That would save a multiplication of epsilon, but would in the long run require
          !multiple epsilon. Leave it for now
          chif(ilevel,ipatch)%fld(i,j,k,isub) = naga_stepfunc1(minimumdist*epsilon)
        END DO !inner array index
      END DO !middle array index
    END DO !outer array index
        ELSE IF (stlopt_inout_direction .EQ. 2) THEN
    DO k=kmin,kmax
      DO j=jmin,jmax
        DO i=imin,imax
          !---------------------------------------------------------------------
          ! First check if the point is inside or outside the solid
          !---------------------------------------------------------------------
          vecT(1) = REAL(i-1+(ptcset(ilevel,ipatch)%sistr(1,isubl)-1),mk)*dx+ptcset(1,1)%min(1)
          vecT(2) = REAL(j-1+(ptcset(ilevel,ipatch)%sistr(2,isubl)-1),mk)*dy+ptcset(1,1)%min(2)
          vecT(3) = REAL(k-1+(ptcset(ilevel,ipatch)%sistr(3,isubl)-1),mk)*dz+ptcset(1,1)%min(3)
          intersections = 0
          inout = 0
          DO tr=1,stlt(1)%tri_count
            !check if triangle is aligned with the point in the bounds_direction
            !perhaps add elseif verbose a more thorough test
            ! maybe change to some finite value
            !implement warning if stl points are too close to grid coordinates
            IF (stlt(1)%tri_norm(tr,2).NE.0.0_mk) THEN
              vecX(1) = vecT(3) - stlt(1)%tri_base(tr,3)
              vecX(2) = vecT(1) - stlt(1)%tri_base(tr,1)
              vecX(3) = vecT(2) - stlt(1)%tri_base(tr,2)
              !beware of vecX index
              udotx = stlt(1)%tri_vecu(tr,3)*vecX(1) + &
                    & stlt(1)%tri_vecu(tr,1)*vecX(2) !@+ &
                    !@& stlt(1)%tri_vecu(tr,2)*vecX(3)
              vdotx = stlt(1)%tri_vecv(tr,3)*vecX(1) + &
                    & stlt(1)%tri_vecv(tr,1)*vecX(2) !@+ &
                    !@& stlt(1)%tri_vecv(tr,2)*vecX(3)
              !get barycentric coordinates
              a = (udotx*stlt(1)%tri_vdotv2d(tr)-vdotx*stlt(1)%tri_udotv2d(tr))*stlt(1)%tri_denom2d(tr)
              b = (vdotx*stlt(1)%tri_udotu2d(tr)-udotx*stlt(1)%tri_udotv2d(tr))*stlt(1)%tri_denom2d(tr)
              c = (1.0_mk-a-b)
              IF (a.GE.0.0_mk .AND. b.GE.0.0_mk.AND.c.GE.0.0_mk) THEN
                dist = a*(stlt(1)%tri_base(tr,2)+stlt(1)%tri_vecu(tr,2)) + &
                     & b*(stlt(1)%tri_base(tr,2)+stlt(1)%tri_vecv(tr,2)) + &
                     & c* stlt(1)%tri_base(tr,2)
                IF (dist .LT. vecT(2)) THEN
                  intersections = intersections + 1
                  IF (stlopt_check_intersections) THEN
                    IF (stlt(1)%tri_norm(tr,2).GT.0.0_mk) THEN
                      inout = inout + 1
                    ELSE
                      inout = inout - 1
                    ENDIF
                  ENDIF
                ENDIF
              ENDIF
            ENDIF
          ENDDO !tr
          minimumdist = diag
          DO tr=1,stlt(1)%tri_count
            vecS=stlt(1)%tri_base(tr,:)-vecT
            !----------------------------------------------------------------
            ! x3 x0 is the projection of T onto the
            ! /\ triangle plane
            ! x0 / \ x1 is used as triangle base point
            ! / \ vecU = x2 - x1
            ! / \ vecV = x3 - x1
            ! /________\ vecW = x3 - x2
            ! x1 x2
            !
            ! x0 x1 T is the point to be evaluated
            ! - - - ======== - - - === is the triangle (in its plane)
            ! | / vecP = x0 - x1
            !vecP | / vecS vecS = x1 - T
            ! | /
            ! |/
            ! T
            !
            ! Meassure distance to triangle.
            ! This is not the vecP from figure 6.2 in thesis
            ! a and b are being reused. Initial definition is in thesis.
            ! tri_norm contains normal vector of triangle plane
            !----------------------------------------------------------------
            vecP = stlt(1)%tri_norm(tr,:) * &
            & ( stlt(1)%tri_norm(tr,1)*vecS(1) + stlt(1)%tri_norm(tr,2)*vecS(2) + &
            & stlt(1)%tri_norm(tr,3)*vecS(3)) - vecS
            Pdotu = vecP(1)*stlt(1)%tri_vecu(tr,1) + vecP(2)*stlt(1)%tri_vecu(tr,2) + &
            & vecP(3)*stlt(1)%tri_vecu(tr,3)
            Pdotv = vecP(1)*stlt(1)%tri_vecv(tr,1) + vecP(2)*stlt(1)%tri_vecv(tr,2) + &
            & vecP(3)*stlt(1)%tri_vecv(tr,3)
            a = (Pdotu*stlt(1)%tri_vdotv(tr) - Pdotv*stlt(1)%tri_udotv(tr)) * stlt(1)%tri_denom(tr)
            b = (Pdotv*stlt(1)%tri_udotu(tr) - Pdotu*stlt(1)%tri_udotv(tr)) * stlt(1)%tri_denom(tr)
            IF ((a .GE. 0.0_mk) .AND. (b .GE. 0.0_mk) .AND. &
            & ((a+b) .LE. 1.0_mk)) THEN
              !get distance normal to triangle
              a = (stlt(1)%tri_norm(tr,1)*vecS(1) + stlt(1)%tri_norm(tr,2)*vecS(2) + &
              & stlt(1)%tri_norm(tr,3)*vecS(3))**2
            ELSE
              !get distance to edge or corner
              a = vecP(1)*stlt(1)%tri_vecu(tr,1) + vecP(2)*stlt(1)%tri_vecu(tr,2) + &
                & vecP(3)*stlt(1)%tri_vecu(tr,3)
              b = vecP(1)*stlt(1)%tri_vecv(tr,1) + vecP(2)*stlt(1)%tri_vecv(tr,2) + &
                & vecP(3)*stlt(1)%tri_vecv(tr,3)
              c = (vecP(1)-stlt(1)%tri_vecu(tr,1))*stlt(1)%tri_vecw(tr,1) + &
                & (vecP(2)-stlt(1)%tri_vecu(tr,2))*stlt(1)%tri_vecw(tr,2) + &
                & (vecP(3)-stlt(1)%tri_vecu(tr,3))*stlt(1)%tri_vecw(tr,3)
              IF (a .GT. stlt(1)%tri_udotu(tr)) THEN !point is closest to x2
                !a = tri_udotu(tr)
                vece=-vecS-stlt(1)%tri_vecu(tr,:)
                a = SUM(vece**2)
              ELSEIF (a .LT. 0.0_mk) THEN !point is closest to x1
                !a = 0
                !vece=-vecS
                a = SUM(vecS**2)
              ELSE !point is closest to vecU
                vece=-vecS-stlt(1)%tri_vecu(tr,:)*a/stlt(1)%tri_udotu(tr)
                a = SUM(vece**2)
                !this may not be the cheapest way to do this
              ENDIF
              IF (b .GT. stlt(1)%tri_vdotv(tr)) THEN !point is closest to x3
                !b = tri_vdotv(tr)
                vece=-vecS-stlt(1)%tri_vecv(tr,:)
                b = SUM(vece**2)
              ELSEIF (b .LT. 0.0_mk) THEN !point is closest to x1
                !b = 0
                !vece=-vecS
                b = SUM(vecS**2)
              ELSE !point is closest to vecV
                vece=-vecS-stlt(1)%tri_vecv(tr,:)*b/stlt(1)%tri_vdotv(tr)
                b = SUM(vece**2)
              ENDIF
              IF (c .GT. stlt(1)%tri_wdotw(tr)) THEN !point is closest to x3
                !c = tri_wdotw(tr)
                vece=-vecS-stlt(1)%tri_vecv(tr,:)
                c = SUM(vece**2)
              ELSEIF (c .LT. 0.0_mk) THEN !point is closest to x2
                !c = 0
                vece=-vecS-stlt(1)%tri_vecu(tr,:)
                c = SUM(vecS**2)
              ELSE !point is closest to vecW
                vece=-vecS-stlt(1)%tri_vecu(tr,:)-stlt(1)%tri_vecw(tr,:)*c/stlt(1)%tri_wdotw(tr)
                c = SUM(vece**2)
              ENDIF
              a = min(a,b)
              a = min(a,c)
            ENDIF
            minimumdist = min(minimumdist,a)
          END DO !tr triangles
          !!chif(ilevel,ipatch)%fld(i,j,k,isub) = naga_stepfunc1(minimumdist*epsilon)
          minimumdist = sqrt(minimumdist)!JTR cleanup
          IF (stlopt_check_intersections) THEN
            IF (inout.LT.0) THEN
              minimumdist = sign(minimumdist,-1.0_mk)
            ENDIF
          ELSE IF(MOD(intersections,2) .NE. 0) THEN
            minimumdist = sign(minimumdist,-1.0_mk)
          ENDIF
          !The multiplication with epsilon can be done to the step function limits instead
          !Could be done either in this routine or when the limits are determined (better).
          !That would save a multiplication of epsilon, but would in the long run require
          !multiple epsilon. Leave it for now
          chif(ilevel,ipatch)%fld(i,j,k,isub) = naga_stepfunc1(minimumdist*epsilon)
        END DO !inner array index
      END DO !middle array index
    END DO !outer array index
        ELSE IF (stlopt_inout_direction .EQ. 3) THEN
    DO k=kmin,kmax
      DO j=jmin,jmax
        DO i=imin,imax
          !---------------------------------------------------------------------
          ! First check if the point is inside or outside the solid
          !---------------------------------------------------------------------
          vecT(1) = REAL(i-1+(ptcset(ilevel,ipatch)%sistr(1,isubl)-1),mk)*dx+ptcset(1,1)%min(1)
          vecT(2) = REAL(j-1+(ptcset(ilevel,ipatch)%sistr(2,isubl)-1),mk)*dy+ptcset(1,1)%min(2)
          vecT(3) = REAL(k-1+(ptcset(ilevel,ipatch)%sistr(3,isubl)-1),mk)*dz+ptcset(1,1)%min(3)
          intersections = 0
          inout = 0
          DO tr=1,stlt(1)%tri_count
            !check if triangle is aligned with the point in the bounds_direction
            !perhaps add elseif verbose a more thorough test
            ! maybe change to some finite value
            !implement warning if stl points are too close to grid coordinates
            IF (stlt(1)%tri_norm(tr,3).NE.0.0_mk) THEN
              vecX(1) = vecT(1) - stlt(1)%tri_base(tr,1)
              vecX(2) = vecT(2) - stlt(1)%tri_base(tr,2)
              vecX(3) = vecT(3) - stlt(1)%tri_base(tr,3)
              !beware of vecX index
              udotx = stlt(1)%tri_vecu(tr,1)*vecX(1) + &
                    & stlt(1)%tri_vecu(tr,2)*vecX(2) !@+ &
                    !@& stlt(1)%tri_vecu(tr,3)*vecX(3)
              vdotx = stlt(1)%tri_vecv(tr,1)*vecX(1) + &
                    & stlt(1)%tri_vecv(tr,2)*vecX(2) !@+ &
                    !@& stlt(1)%tri_vecv(tr,3)*vecX(3)
              !get barycentric coordinates
              a = (udotx*stlt(1)%tri_vdotv2d(tr)-vdotx*stlt(1)%tri_udotv2d(tr))*stlt(1)%tri_denom2d(tr)
              b = (vdotx*stlt(1)%tri_udotu2d(tr)-udotx*stlt(1)%tri_udotv2d(tr))*stlt(1)%tri_denom2d(tr)
              c = (1.0_mk-a-b)
              IF (a.GE.0.0_mk .AND. b.GE.0.0_mk.AND.c.GE.0.0_mk) THEN
                dist = a*(stlt(1)%tri_base(tr,3)+stlt(1)%tri_vecu(tr,3)) + &
                     & b*(stlt(1)%tri_base(tr,3)+stlt(1)%tri_vecv(tr,3)) + &
                     & c* stlt(1)%tri_base(tr,3)
                IF (dist .LT. vecT(3)) THEN
                  intersections = intersections + 1
                  IF (stlopt_check_intersections) THEN
                    IF (stlt(1)%tri_norm(tr,3).GT.0.0_mk) THEN
                      inout = inout + 1
                    ELSE
                      inout = inout - 1
                    ENDIF
                  ENDIF
                ENDIF
              ENDIF
            ENDIF
          ENDDO !tr
          minimumdist = diag
          DO tr=1,stlt(1)%tri_count
            vecS=stlt(1)%tri_base(tr,:)-vecT
            !----------------------------------------------------------------
            ! x3 x0 is the projection of T onto the
            ! /\ triangle plane
            ! x0 / \ x1 is used as triangle base point
            ! / \ vecU = x2 - x1
            ! / \ vecV = x3 - x1
            ! /________\ vecW = x3 - x2
            ! x1 x2
            !
            ! x0 x1 T is the point to be evaluated
            ! - - - ======== - - - === is the triangle (in its plane)
            ! | / vecP = x0 - x1
            !vecP | / vecS vecS = x1 - T
            ! | /
            ! |/
            ! T
            !
            ! Meassure distance to triangle.
            ! This is not the vecP from figure 6.2 in thesis
            ! a and b are being reused. Initial definition is in thesis.
            ! tri_norm contains normal vector of triangle plane
            !----------------------------------------------------------------
            vecP = stlt(1)%tri_norm(tr,:) * &
            & ( stlt(1)%tri_norm(tr,1)*vecS(1) + stlt(1)%tri_norm(tr,2)*vecS(2) + &
            & stlt(1)%tri_norm(tr,3)*vecS(3)) - vecS
            Pdotu = vecP(1)*stlt(1)%tri_vecu(tr,1) + vecP(2)*stlt(1)%tri_vecu(tr,2) + &
            & vecP(3)*stlt(1)%tri_vecu(tr,3)
            Pdotv = vecP(1)*stlt(1)%tri_vecv(tr,1) + vecP(2)*stlt(1)%tri_vecv(tr,2) + &
            & vecP(3)*stlt(1)%tri_vecv(tr,3)
            a = (Pdotu*stlt(1)%tri_vdotv(tr) - Pdotv*stlt(1)%tri_udotv(tr)) * stlt(1)%tri_denom(tr)
            b = (Pdotv*stlt(1)%tri_udotu(tr) - Pdotu*stlt(1)%tri_udotv(tr)) * stlt(1)%tri_denom(tr)
            IF ((a .GE. 0.0_mk) .AND. (b .GE. 0.0_mk) .AND. &
            & ((a+b) .LE. 1.0_mk)) THEN
              !get distance normal to triangle
              a = (stlt(1)%tri_norm(tr,1)*vecS(1) + stlt(1)%tri_norm(tr,2)*vecS(2) + &
              & stlt(1)%tri_norm(tr,3)*vecS(3))**2
            ELSE
              !get distance to edge or corner
              a = vecP(1)*stlt(1)%tri_vecu(tr,1) + vecP(2)*stlt(1)%tri_vecu(tr,2) + &
                & vecP(3)*stlt(1)%tri_vecu(tr,3)
              b = vecP(1)*stlt(1)%tri_vecv(tr,1) + vecP(2)*stlt(1)%tri_vecv(tr,2) + &
                & vecP(3)*stlt(1)%tri_vecv(tr,3)
              c = (vecP(1)-stlt(1)%tri_vecu(tr,1))*stlt(1)%tri_vecw(tr,1) + &
                & (vecP(2)-stlt(1)%tri_vecu(tr,2))*stlt(1)%tri_vecw(tr,2) + &
                & (vecP(3)-stlt(1)%tri_vecu(tr,3))*stlt(1)%tri_vecw(tr,3)
              IF (a .GT. stlt(1)%tri_udotu(tr)) THEN !point is closest to x2
                !a = tri_udotu(tr)
                vece=-vecS-stlt(1)%tri_vecu(tr,:)
                a = SUM(vece**2)
              ELSEIF (a .LT. 0.0_mk) THEN !point is closest to x1
                !a = 0
                !vece=-vecS
                a = SUM(vecS**2)
              ELSE !point is closest to vecU
                vece=-vecS-stlt(1)%tri_vecu(tr,:)*a/stlt(1)%tri_udotu(tr)
                a = SUM(vece**2)
                !this may not be the cheapest way to do this
              ENDIF
              IF (b .GT. stlt(1)%tri_vdotv(tr)) THEN !point is closest to x3
                !b = tri_vdotv(tr)
                vece=-vecS-stlt(1)%tri_vecv(tr,:)
                b = SUM(vece**2)
              ELSEIF (b .LT. 0.0_mk) THEN !point is closest to x1
                !b = 0
                !vece=-vecS
                b = SUM(vecS**2)
              ELSE !point is closest to vecV
                vece=-vecS-stlt(1)%tri_vecv(tr,:)*b/stlt(1)%tri_vdotv(tr)
                b = SUM(vece**2)
              ENDIF
              IF (c .GT. stlt(1)%tri_wdotw(tr)) THEN !point is closest to x3
                !c = tri_wdotw(tr)
                vece=-vecS-stlt(1)%tri_vecv(tr,:)
                c = SUM(vece**2)
              ELSEIF (c .LT. 0.0_mk) THEN !point is closest to x2
                !c = 0
                vece=-vecS-stlt(1)%tri_vecu(tr,:)
                c = SUM(vecS**2)
              ELSE !point is closest to vecW
                vece=-vecS-stlt(1)%tri_vecu(tr,:)-stlt(1)%tri_vecw(tr,:)*c/stlt(1)%tri_wdotw(tr)
                c = SUM(vece**2)
              ENDIF
              a = min(a,b)
              a = min(a,c)
            ENDIF
            minimumdist = min(minimumdist,a)
          END DO !tr triangles
          !!chif(ilevel,ipatch)%fld(i,j,k,isub) = naga_stepfunc1(minimumdist*epsilon)
          minimumdist = sqrt(minimumdist)!JTR cleanup
          IF (stlopt_check_intersections) THEN
            IF (inout.LT.0) THEN
              minimumdist = sign(minimumdist,-1.0_mk)
            ENDIF
          ELSE IF(MOD(intersections,2) .NE. 0) THEN
            minimumdist = sign(minimumdist,-1.0_mk)
          ENDIF
          !The multiplication with epsilon can be done to the step function limits instead
          !Could be done either in this routine or when the limits are determined (better).
          !That would save a multiplication of epsilon, but would in the long run require
          !multiple epsilon. Leave it for now
          chif(ilevel,ipatch)%fld(i,j,k,isub) = naga_stepfunc1(minimumdist*epsilon)
        END DO !inner array index
      END DO !middle array index
    END DO !outer array index
        ENDIF
      END DO !isub
    END DO !ipatch
  END DO !ilevel
END SUBROUTINE naga_stl_init
END MODULE naga_mod_stl
