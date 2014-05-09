!------------------------------------------------------------------------------
! Subroutine : naga_read_ctrl.f90
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
! This routines reads the input/control file
! Please remember to define default values in naga_defaults
!------------------------------------------------------------------------------
MODULE naga_mod_read_ctrl
IMPLICIT NONE
INTERFACE naga_read_ctrl
  MODULE PROCEDURE naga_read_ctrl
END INTERFACE
CONTAINS
SUBROUTINE naga_read_ctrl(cfile,info)
USE naga_mod_globals
USE naga_mod_say
IMPLICIT NONE
!-------------------------------------------------------------------------
! Arguments
!-------------------------------------------------------------------------
CHARACTER(LEN=*), INTENT(IN) :: cfile
INTEGER, INTENT(OUT) :: info
!-------------------------------------------------------------------------
! Local variables
!-------------------------------------------------------------------------
REAL(MK) :: t0
INTEGER :: info2
INTEGER :: i,j,idx
INTEGER :: iunit, cfilelen, ios, ilen, iline
CHARACTER(LEN=MAXCHAR) :: entrybuf, entryarg, entryval
LOGICAL :: cfileexists
INTEGER :: layoutread
INTEGER :: ilevel,ipatch
INTEGER, DIMENSION(3) :: setbuf !buffer for patch settings
REAL(MK),DIMENSION(6) :: dimbuf !buffer for patch dims.
INTEGER, DIMENSION(3) :: resbuf !buffer for resolution
!----------------------------------------------------------------------------
! Initialise routine
!----------------------------------------------------------------------------
CALL substart('naga_read_ctrl',t0,info)
!----------------------------------------------------------------------------
! Define the unit number to read control file from
!----------------------------------------------------------------------------
iunit = 20
!----------------------------------------------------------------------------
! Define that the layout has not been read
!----------------------------------------------------------------------------
layoutread = 0
!----------------------------------------------------------------------------
! Check the control file
!----------------------------------------------------------------------------
cfilelen = LEN_TRIM(cfile)
IF (cfilelen.LT.1) THEN
  CALL naga_say(rank,'naga_read_ctrl','No control file given.')
  info = -1
  GOTO 9999
ENDIF
INQUIRE(FILE=cfile,EXIST=cfileexists)
IF (.NOT.cfileexists) THEN
  WRITE(entrybuf,'(2A)')'No such file: ',cfile(1:cfilelen)
  CALL naga_say(rank,'naga_read_ctrl',entrybuf)
  info = -1
  GOTO 9999
ENDIF
!-------------------------------------------------------------------------
! Open the file
!-------------------------------------------------------------------------
OPEN(iunit,FILE=cfile,IOSTAT=ios,ACTION='READ')
IF (ios.NE.0) THEN
  WRITE(entrybuf,'(2A)')'Failed to open control file: ',cfile(1:cfilelen)
  CALL naga_say(rank,'naga_read_ctrl',entrybuf)
  info = -1
  GOTO 9999
ENDIF
!-------------------------------------------------------------------------
! Scan file
!-------------------------------------------------------------------------
iline = 0
DO
  !----------------------------------------------------------------------
  ! Increment line counter
  !----------------------------------------------------------------------
  iline = iline + 1
  !----------------------------------------------------------------------
  ! Read line
  !----------------------------------------------------------------------
  READ(iUnit,'(A)',END=100,ERR=200) entrybuf
  ilen = LEN_TRIM(entrybuf)
  !----------------------------------------------------------------------
  ! Debug print of Ctrl file
  !----------------------------------------------------------------------
  IF (debug.GT.1) THEN
    CALL naga_say(rank,'naga_read_ctrl',entrybuf(1:ilen),info)
  ENDIF
  !----------------------------------------------------------------------
  ! Skip comment or empty lines·
  !----------------------------------------------------------------------
  IF (ilen.GT.0.AND.entrybuf(1:1).NE.'#') THEN
    !-------------------------------------------------------------------
    ! Remove space
    !-------------------------------------------------------------------
    j = 0
    DO i=1,ilen
      IF (entrybuf(i:i).NE.' ') THEN
        j = j + 1
        entrybuf(j:j) = entrybuf(i:i)
      ENDIF
    ENDDO
    !-------------------------------------------------------------------
    ! Update length of string
    !-------------------------------------------------------------------
    ilen = j
    !-------------------------------------------------------------------
    ! Find position of =
    !-------------------------------------------------------------------
    idx = INDEX(entrybuf,'=')
    !-------------------------------------------------------------------
    ! Exit if missing
    !-------------------------------------------------------------------
    IF (idx.LT.0) THEN
      WRITE(entrybuf,'(A,I5)')'Incorrect line: ',iline
      CALL naga_say(rank,'naga_read_ctrl',entrybuf,info)
      info = -1
      GOTO 9999
    ENDIF
    !-------------------------------------------------------------------
    ! Get argument and value
    !-------------------------------------------------------------------
    entryarg = ADJUSTL(entrybuf(1:idx-1))
    entryval = ADJUSTL(entrybuf(idx+1:ilen))
    !-------------------------------------------------------------------
    ! Convert to upper case
    !-------------------------------------------------------------------
    CALL UpperCase(entryarg,idx-1,info)
    !-------------------------------------------------------------------
    ! Update the appropriate input
    ! Please remember to define default values in naga_defaults
    !-------------------------------------------------------------------
    IF (entryarg.EQ.'DEBUG') THEN
      READ(entryval,*,IOSTAT=ios,ERR=200) debug
    ELSE IF (entryarg.EQ.'VALIDATIONFIELD') THEN
      READ(entryval,*,IOSTAT=ios,ERR=200) validationfield
    ELSE IF (entryarg.EQ.'FLOWCASE') THEN
      READ(entryval,*,IOSTAT=ios,ERR=200) flowcase
    ELSE IF (entryarg.EQ.'GHOSTWIDTH') THEN
      READ(entryval,*,IOSTAT=ios,ERR=200) gstw
    ELSE IF (entryarg.EQ.'TAG') THEN
      READ(entryval,*,IOSTAT=ios,ERR=200) runtag
      iruntag = LEN_TRIM(runtag)
    !!ELSE IF (entryarg.EQ.'LAYOUTFILE') THEN
      !!READ(entryval,*,IOSTAT=ios,ERR=200) layoutfile
      !!ilayoutfile = LEN_TRIM(layoutfile)
    ELSE IF (entryarg.EQ.'DT') THEN
      READ(entryval,*,IOSTAT=ios,ERR=200) dtime
    ELSE IF (entryarg.EQ.'ENDTIME') THEN
      READ(entryval,*,IOSTAT=ios,ERR=200) endtime
    ELSE IF (entryarg.EQ.'IREMESH') THEN
      READ(entryval,*,IOSTAT=ios,ERR=200) iremesh
    ELSE IF (entryarg.EQ.'IREPROJECT') THEN
      READ(entryval,*,IOSTAT=ios,ERR=200) ireproject
    ELSE IF (entryarg.EQ.'IDUMP_VORTICITY') THEN
      READ(entryval,*,IOSTAT=ios,ERR=200) idumpvrt
    ELSE IF (entryarg.EQ.'IDUMP_RHS') THEN
      READ(entryval,*,IOSTAT=ios,ERR=200) idumpdvrt
    ELSE IF (entryarg.EQ.'IDUMP_VELOCITY') THEN
      READ(entryval,*,IOSTAT=ios,ERR=200) idumpvel
    ELSE IF (entryarg.EQ.'IDUMP_MASK') THEN
      READ(entryval,*,IOSTAT=ios,ERR=200) idumpchi
    ELSE IF (entryarg.EQ.'IDUMP_CONCENTRATION') THEN
      READ(entryval,*,IOSTAT=ios,ERR=200) idumpconc
    ELSE IF (entryarg.EQ.'BOUNDARY_CONDITION') THEN
      READ(entryval,*,IOSTAT=ios,ERR=200) domainbc
    ELSE IF (entryarg.EQ.'UINFINITY') THEN
      READ(entryval,*,IOSTAT=ios,ERR=200) uinfinity
    ELSE IF (entryarg.EQ.'MAXITIME') THEN
      READ(entryval,*,IOSTAT=ios,ERR=200) maxitime
    ELSE IF (entryarg.EQ.'TIMEINTEGRATIONSCHEME') THEN
      READ(entryval,*,IOSTAT=ios,ERR=200) timeintscheme
    ELSE IF (entryarg.EQ.'PENALIZATIONSCHEME') THEN
      READ(entryval,*,IOSTAT=ios,ERR=200) penalizationscheme
    ELSE IF (entryarg.EQ.'CLEARINTERIORRHS') THEN
      READ(entryval,*,IOSTAT=ios,ERR=200) clearinteriorrhs
    ELSE IF (entryarg.EQ.'RHSSCHEME') THEN
      READ(entryval,*,IOSTAT=ios,ERR=200) rhsscheme
    ELSE IF (entryarg.EQ.'FREESPACEKERNEL') THEN
      READ(entryval,*,IOSTAT=ios,ERR=200) freespacekernel
    ELSE IF (entryarg.EQ.'VELOCITYSCHEME') THEN
      READ(entryval,*,IOSTAT=ios,ERR=200) velocityscheme
    ELSE IF (entryarg.EQ.'PENALIZATIONPARAMETER') THEN
      READ(entryval,*,IOSTAT=ios,ERR=200) penalizationparam
    ELSE IF (entryarg.EQ.'PENALIZATIONADAPT') THEN
      READ(entryval,*,IOSTAT=ios,ERR=200) penalizationadapt
    ELSE IF (entryarg.EQ.'PENALIZATION') THEN
      READ(entryval,*,IOSTAT=ios,ERR=200) penalization
    ELSE IF (entryarg.EQ.'CONCENTRATION') THEN
      READ(entryval,*,IOSTAT=ios,ERR=200) concentration
    ELSE IF (entryarg.EQ.'NU') THEN
      READ(entryval,*,IOSTAT=ios,ERR=200) nu
    ELSE IF (entryarg.EQ.'LESMODEL') THEN
      READ(entryval,*,IOSTAT=ios,ERR=200) lesmodel
    ELSE IF (entryarg.EQ.'TORUS_RADIUS1') THEN
      READ(entryval,*,IOSTAT=ios,ERR=200) torus_radius1
    ELSE IF (entryarg.EQ.'TORUS_RADIUS2') THEN
      READ(entryval,*,IOSTAT=ios,ERR=200) torus_radius2
    ELSE IF (entryarg.EQ.'SPHERE_RADIUS') THEN
      READ(entryval,*,IOSTAT=ios,ERR=200) sphere_radius
    ELSE IF (entryarg.EQ.'STEP1_INTERVAL') THEN
      READ(entryval,*,IOSTAT=ios,ERR=200) step1_interval
    ELSE IF (entryarg.EQ.'STEP1_LINEARFRACTION') THEN
      READ(entryval,*,IOSTAT=ios,ERR=200) step1_linearfraction
    ELSE IF (entryarg.EQ.'STEP1_OFFSET') THEN
      READ(entryval,*,IOSTAT=ios,ERR=200) step1_offset
    ELSE IF (entryarg.EQ.'STL_FILENAME') THEN
      READ(entryval,*,IOSTAT=ios,ERR=200) stlopt_filename
    ELSE IF (entryarg.EQ.'STL_SCALE') THEN
      READ(entryval,*,IOSTAT=ios,ERR=200) stlopt_scale
    ELSE IF (entryarg.EQ.'STL_TRANSLATE') THEN
      READ(entryval,*,IOSTAT=ios,ERR=200) stlopt_translate
    ELSE IF (entryarg.EQ.'STL_CHECK_BOUNDS') THEN
      READ(entryval,*,IOSTAT=ios,ERR=200) stlopt_check_bounding
    ELSE IF (entryarg.EQ.'STL_CHECK_INTERSECTIONS') THEN
      READ(entryval,*,IOSTAT=ios,ERR=200) stlopt_check_intersections
    ELSE IF (entryarg.EQ.'STL_INOUT_DIRECTION') THEN
      READ(entryval,*,IOSTAT=ios,ERR=200) stlopt_inout_direction
    !----------------------------------------------------------------------
    ! *** *** *** *** Some ASCII flashing *** *** *** ***
    ! READING THE LAYOUT (set-up of patches)
    ! Here we have few more lines as a lot of variables have to be set
    !----------------------------------------------------------------------
    ! First the dimensions of the layout needs to be set
    ELSE IF (entryarg.EQ.'LAYOUT_SIZE') THEN
      ! A check to see if we have already set the number of patches
      IF (layoutread .NE. 0) THEN
        WRITE(entrybuf,'(2A,I5)') 'Layout dimensions can only be set once, ',&
        &'and before the layout is being initialised. Error on line ',iline
        ilen = LEN_TRIM(entrybuf)
        CALL naga_say(rank,'naga_read_ctrl',entrybuf(1:ilen),info)
        info = -1
        GOTO 9999
      ! And the action:
      ELSE
        READ(entryval,*,IOSTAT=ios,ERR=200) nlevels,maxpatches
        ALLOCATE(npatches(nlevels))
        ALLOCATE(ptcset(nlevels,maxpatches))
        nlevels = 0
        npatches = 0
        layoutread = 1
      ENDIF
    ELSE IF (entryarg.EQ.'LAYOUT') THEN
      !first check if layout dimensions have been set
      IF (layoutread .NE. 1) THEN
        WRITE(entrybuf,'(2A)') 'Layout dimensions (LAYOUT_SIZE=nlevels,',&
        &'maxpatches) must be set before the patch layout is being initialised.'
        ilen = LEN_TRIM(entrybuf)
        CALL naga_say(rank,'naga_read_ctrl',entrybuf(1:ilen),info)
        info = -1
        GOTO 9999
      ENDIF
      !----------------------------------------------------------------------
      ! Use npatches to keep track of the next patch numbers on all levels,
      ! therby recounting npatches.
      ! Save the settings:
      !
      !----------------------------------------------------------------------
      READ(entryval,*,IOSTAT=ios,ERR=200) setbuf,dimbuf,resbuf
      !determine the number of levels (the highest level)
      nlevels = MAX(nlevels,setbuf(1))
      !update the present patch of the level
      npatches(setbuf(1)) = npatches(setbuf(1)) + 1
      !resolution level
      ptcset(setbuf(1),npatches(setbuf(1)))%lvl = setbuf(1)
      !patch parent
      ptcset(setbuf(1),npatches(setbuf(1)))%par = setbuf(2)
      !patch buffer layer width
      ptcset(setbuf(1),npatches(setbuf(1)))%buf = setbuf(3)
      !min of patch
      ptcset(setbuf(1),npatches(setbuf(1)))%min(1) = dimbuf(1)
      ptcset(setbuf(1),npatches(setbuf(1)))%min(2) = dimbuf(2)
      ptcset(setbuf(1),npatches(setbuf(1)))%min(3) = dimbuf(3)
      !max of patch
      ptcset(setbuf(1),npatches(setbuf(1)))%max(1) = dimbuf(4)
      ptcset(setbuf(1),npatches(setbuf(1)))%max(2) = dimbuf(5)
      ptcset(setbuf(1),npatches(setbuf(1)))%max(3) = dimbuf(6)
      !resolution
      ptcset(setbuf(1),npatches(setbuf(1)))%gnx = resbuf
    ELSE
      WRITE(entrybuf,'(A,I5,2A)')'Warning: Ignoring setting on line ',iline, &
          & ' : ', entryarg(1:idx-1)
      CALL naga_say(rank,'naga_read_ctrl',entrybuf,info)
    END IF
  END IF
END DO
!-------------------------------------------------------------------------
! Something went wrong
!-------------------------------------------------------------------------
  200 CONTINUE
WRITE(entrybuf,'(A,I5,2A)') 'Error reading line: ',iline, &
& ' of file: ',cfile(1:cfilelen)
ilen = LEN_TRIM(entrybuf)
CALL naga_say(rank,'naga_read_ctrl',entrybuf(1:ilen),info)
info = -1
GOTO 9999
!-------------------------------------------------------------------------
! End of file - post-process data
!-------------------------------------------------------------------------
  100 CONTINUE
DO ilevel=1,nlevels
  DO ipatch=1,npatches(ilevel)
    IF (ilevel .EQ. 1) THEN !check if level=1
      !compute grid spacing - subtract 1 grid point for periodic domains
      IF (domainbc .EQ. 0) THEN
        ptcset(ilevel,npatches(ilevel))%dx(1) = &
        & ( ptcset(ilevel,npatches(ilevel))%max(1) &
        & - ptcset(ilevel,npatches(ilevel))%min(1) ) &
        & / REAL(ptcset(ilevel,npatches(ilevel))%gnx (1)-1,MK)
        ptcset(ilevel,npatches(ilevel))%dx(2) = &
        & ( ptcset(ilevel,npatches(ilevel))%max(2) &
        & - ptcset(ilevel,npatches(ilevel))%min(2) ) &
        & / REAL(ptcset(ilevel,npatches(ilevel))%gnx (2)-1,MK)
        ptcset(ilevel,npatches(ilevel))%dx(3) = &
        & ( ptcset(ilevel,npatches(ilevel))%max(3) &
        & - ptcset(ilevel,npatches(ilevel))%min(3) ) &
        & / REAL(ptcset(ilevel,npatches(ilevel))%gnx (3)-1,MK)
      ELSEIF (domainbc .EQ. 1) THEN
        ptcset(ilevel,npatches(ilevel))%dx(1) = &
        & ( ptcset(ilevel,npatches(ilevel))%max(1) &
        & - ptcset(ilevel,npatches(ilevel))%min(1) ) &
        & / REAL(ptcset(ilevel,npatches(ilevel))%gnx (1)-1,MK)
        ptcset(ilevel,npatches(ilevel))%dx(2) = &
        & ( ptcset(ilevel,npatches(ilevel))%max(2) &
        & - ptcset(ilevel,npatches(ilevel))%min(2) ) &
        & / REAL(ptcset(ilevel,npatches(ilevel))%gnx (2)-1,MK)
        ptcset(ilevel,npatches(ilevel))%dx(3) = &
        & ( ptcset(ilevel,npatches(ilevel))%max(3) &
        & - ptcset(ilevel,npatches(ilevel))%min(3) ) &
        & / REAL(ptcset(ilevel,npatches(ilevel))%gnx (3)-1,MK)
      ELSE
        CALL naga_say(rank,'naga_read_ctrl','In control file: wrong BOUNDARY_CONDITION',info)
      ENDIF
    ELSE
      !for all child patches fit their bounds to parent meshes
      ptcset( ilevel , ipatch )%min(1) = &
      & ANINT( (ptcset( ilevel , ipatch )%min(1) &
      & - ptcset(ilevel-1,ptcset(ilevel,ipatch)%par)%min(1) ) &
      & / ptcset(ilevel-1,ptcset(ilevel,ipatch)%par)% dx(1),MK) &
      & * ptcset(ilevel-1,ptcset(ilevel,ipatch)%par)% dx(1)
      ptcset( ilevel , ipatch )%min(2) = &
      & ANINT( (ptcset( ilevel , ipatch )%min(2) &
      & - ptcset(ilevel-1,ptcset(ilevel,ipatch)%par)%min(2) ) &
      & / ptcset(ilevel-1,ptcset(ilevel,ipatch)%par)% dx(2),MK) &
      & * ptcset(ilevel-1,ptcset(ilevel,ipatch)%par)% dx(2)
      ptcset( ilevel , ipatch )%min(3) = &
      & ANINT( (ptcset( ilevel , ipatch )%min(3) &
      & - ptcset(ilevel-1,ptcset(ilevel,ipatch)%par)%min(3) ) &
      & / ptcset(ilevel-1,ptcset(ilevel,ipatch)%par)% dx(3),MK) &
      & * ptcset(ilevel-1,ptcset(ilevel,ipatch)%par)% dx(3)
      !max of patch
      ptcset( ilevel , ipatch )%max(1) = &
      & ANINT( (ptcset( ilevel , ipatch )%max(1) &
      & - ptcset(ilevel-1,ptcset(ilevel,ipatch)%par)%min(1) ) &
      & / ptcset(ilevel-1,ptcset(ilevel,ipatch)%par)% dx(1),MK) &
      & * ptcset(ilevel-1,ptcset(ilevel,ipatch)%par)% dx(1)
      ptcset( ilevel , ipatch )%max(2) = &
      & ANINT( (ptcset( ilevel , ipatch )%max(2) &
      & - ptcset(ilevel-1,ptcset(ilevel,ipatch)%par)%min(2) ) &
      & / ptcset(ilevel-1,ptcset(ilevel,ipatch)%par)% dx(2),MK) &
      & * ptcset(ilevel-1,ptcset(ilevel,ipatch)%par)% dx(2)
      ptcset( ilevel , ipatch )%max(3) = &
      & ANINT( (ptcset( ilevel , ipatch )%max(3) &
      & - ptcset(ilevel-1,ptcset(ilevel,ipatch)%par)%min(3) ) &
      & / ptcset(ilevel-1,ptcset(ilevel,ipatch)%par)% dx(3),MK) &
      & * ptcset(ilevel-1,ptcset(ilevel,ipatch)%par)% dx(3)
      !grid spacing
      ptcset(ilevel,ipatch)%dx(1) = &
      & ptcset(ilevel-1,ptcset(ilevel,ipatch)%par)%dx(1)/2.0_MK
      ptcset(ilevel,ipatch)%dx(2) = &
      & ptcset(ilevel-1,ptcset(ilevel,ipatch)%par)%dx(2)/2.0_MK
      ptcset(ilevel,ipatch)%dx(3) = &
      & ptcset(ilevel-1,ptcset(ilevel,ipatch)%par)%dx(3)/2.0_MK
      !resolution
      ptcset(ilevel,ipatch)%gnx(1) = NINT( &
      & ( ptcset(ilevel,ipatch)%max(1) &
      & - ptcset(ilevel,ipatch)%min(1)) &
      & / ptcset(ilevel,ipatch)% dx(1))
      ptcset(ilevel,ipatch)%gnx(2) = NINT( &
      & ( ptcset(ilevel,ipatch)%max(2) &
      & - ptcset(ilevel,ipatch)%min(2)) &
      & / ptcset(ilevel,ipatch)% dx(2))
      ptcset(ilevel,ipatch)%gnx(3) = NINT( &
      & ( ptcset(ilevel,ipatch)%max(3) &
      & - ptcset(ilevel,ipatch)%min(3)) &
      & / ptcset(ilevel,ipatch)% dx(3))
    ENDIF
  END DO
END DO
!-------------------------------------------------------------------------
! Close file
!-------------------------------------------------------------------------
CLOSE(iunit)
!----------------------------------------------------------------------------
! Return
!----------------------------------------------------------------------------
CALL substop('naga_read_ctrl',t0,info2)
 9999 CONTINUE
RETURN
END SUBROUTINE naga_read_ctrl
END MODULE naga_mod_read_ctrl
