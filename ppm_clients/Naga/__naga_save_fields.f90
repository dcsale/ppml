!------------------------------------------------------------------------------
! Subroutine : naga_save_fields.f90
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
! This routines calls the vtk routine for the fields depending of choice:
! choice is a factor of 2: vorticity
! choice is a factor of 3: velocity
! choice is a factor of 5: solid mask
! choice is a factor of 7: concentration
! choice is a factor of 11: vorticity RHS
!------------------------------------------------------------------------------
MODULE naga_mod_save_fields
IMPLICIT NONE
INTERFACE naga_save_fields
  MODULE PROCEDURE naga_save_fields
END INTERFACE
CONTAINS
SUBROUTINE naga_save_fields(info,choice)
USE naga_mod_globals
USE naga_mod_say
USE naga_mod_dump_vtk
IMPLICIT NONE
!-------------------------------------------------------------------------
! Arguments
!-------------------------------------------------------------------------
INTEGER, INTENT(OUT) :: info
INTEGER, INTENT(IN) :: choice
!-------------------------------------------------------------------------
! Local variables
!-------------------------------------------------------------------------
REAL(MK) :: t0
CHARACTER(len=256) :: msg
INTEGER :: ilevel,ipatch
LOGICAL :: dumpvrt,dumpvel,dumpchi,dumpcnc,dumpdwp
!----------------------------------------------------------------------------
! Initialise routine
!----------------------------------------------------------------------------
CALL substart('naga_save_fields',t0,info)
info = 0
!----------------------------------------------------------------------------
! Determine whether to dump or not
!----------------------------------------------------------------------------
dumpvrt = .FALSE.
dumpvel = .FALSE.
dumpchi = .FALSE.
dumpcnc = .FALSE.
dumpdwp = .FALSE.
IF (MOD(choice,2) .EQ. 0) THEN
  IF (idumpvrt .NE. 0) THEN
    IF (MOD(itime,idumpvrt) .EQ. 0) THEN
      dumpvrt = .TRUE.
    ENDIF
  ENDIF
  IF (forcedumpvrt) THEN
    dumpvrt = .TRUE.
  ENDIF
ENDIF
IF (MOD(choice,3) .EQ. 0) THEN
  IF (idumpvel .NE. 0) THEN
    IF (MOD(itime,idumpvel) .EQ. 0) THEN
      dumpvel = .TRUE.
    ENDIF
  ENDIF
  IF (forcedumpvel) THEN
    dumpvel = .TRUE.
  ENDIF
ENDIF
IF (MOD(choice,5) .EQ. 0) THEN
  IF (idumpchi .NE. 0) THEN
    IF (MOD(itime,idumpchi) .EQ. 0) THEN
      dumpchi = .TRUE.
    ENDIF
  ENDIF
  IF (forcedumpchi) THEN
    dumpchi = .TRUE.
  ENDIF
ENDIF
IF (MOD(choice,7) .EQ. 0 .AND. concentration .EQ. 1) THEN
  IF (idumpconc .NE. 0) THEN
    IF (MOD(itime,idumpconc) .EQ. 0) THEN
      dumpcnc = .TRUE.
    ENDIF
  ENDIF
  IF (forcedumpcnc) THEN
    dumpcnc = .TRUE.
  ENDIF
ENDIF
IF (MOD(choice,11) .EQ. 0) THEN
  IF (idumpdvrt .NE. 0) THEN
    IF (MOD(itime,idumpdvrt) .EQ. 0) THEN
      dumpdwp = .TRUE.
    ENDIF
  ENDIF
  IF (forcedumpdwp) THEN
    dumpdwp = .TRUE.
  ENDIF
ENDIF
!----------------------------------------------------------------------------
! Assuming all fields are ready we call the vtk routines (if it is time):
!----------------------------------------------------------------------------
DO ilevel=1,nlevels
  DO ipatch=1,npatches(ilevel)
    IF (dumpvrt) THEN
      msg = 'vrt'
      CALL naga_dump_vtk(ptcset(ilevel,ipatch)%topoid, &
      & ptcset(ilevel,ipatch)%meshid, wf(ilevel,ipatch)%fld,msg,itime,info)
      IF (info .NE. 0) THEN
        CALL naga_say(rank, 'naga_save_fields','Failed to save vorticity as vtk.')
        GOTO 9999
      ELSE
        IF (rank .EQ. 0) THEN
          CALL naga_say(rank, 'naga_save_fields','Saved vorticity as vtk.',info)
        ENDIF
      ENDIF
    ENDIF
    IF (dumpvel) THEN
      msg = 'vel'
      CALL naga_dump_vtk(ptcset(ilevel,ipatch)%topoid, &
      & ptcset(ilevel,ipatch)%meshid, uf(ilevel,ipatch)%fld,msg,itime,info)
      IF (info .NE. 0) THEN
        CALL naga_say(rank, 'naga_save_fields','Failed to save velocity as vtk.')
        GOTO 9999
      ELSE
        IF (rank .EQ. 0) THEN
          CALL naga_say(rank, 'naga_save_fields','Saved velocity as vtk.',info)
        ENDIF
      ENDIF
    ENDIF
    IF (dumpchi .AND. (penalization .NE. 0)) THEN
      msg = 'chi'
      CALL naga_dump_vtk(ptcset(ilevel,ipatch)%topoid, &
      & ptcset(ilevel,ipatch)%meshid, chif(ilevel,ipatch)%fld,msg,itime,info)
      IF (info .NE. 0) THEN
        CALL naga_say(rank, 'naga_save_fields','Failed to save mask as vtk.')
        GOTO 9999
      ELSE
        IF (rank .EQ. 0) THEN
          CALL naga_say(rank, 'naga_save_fields','Saved mask as vtk.',info)
        ENDIF
      ENDIF
    ENDIF
    IF (dumpcnc) THEN
      msg = 'conc'
      CALL naga_dump_vtk(ptcset(ilevel,ipatch)%topoid, &
      & ptcset(ilevel,ipatch)%meshid, cf(ilevel,ipatch)%fld,msg,itime,info)
      IF (info .NE. 0) THEN
        CALL naga_say(rank, 'naga_save_fields','Failed to save concentration as vtk.')
        GOTO 9999
      ELSE
        IF (rank .EQ. 0) THEN
          CALL naga_say(rank, 'naga_save_fields','Saved concentration as vtk.',info)
        ENDIF
      ENDIF
    ENDIF
    IF (dumpdwp) THEN
      msg = 'dvrt'
      CALL naga_dump_vtk(ptcset(ilevel,ipatch)%topoid, &
      & ptcset(ilevel,ipatch)%meshid, dwf(ilevel,ipatch)%fld,msg,itime,info)
      IF (info .NE. 0) THEN
        CALL naga_say(rank, 'naga_save_fields','Failed to save vorticity RHS as vtk.')
        GOTO 9999
      ELSE
        IF (rank .EQ. 0) THEN
          CALL naga_say(rank, 'naga_save_fields','Saved vorticity RHS as vtk.',info)
        ENDIF
      ENDIF
    ENDIF
  END DO !ipatch
END DO !ilevel
!----------------------------------------------------------------------------
! Return
!----------------------------------------------------------------------------
 9999 CONTINUE
CALL substop('naga_save_fields',t0,info)
RETURN
END SUBROUTINE naga_save_fields
END MODULE naga_mod_save_fields
