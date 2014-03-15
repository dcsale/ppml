!------------------------------------------------------------------------------
! Subroutine : naga_ghost_fields.f90
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
! This routine ghosts vorticity. Could be extended with arguments to toggle
! ghosting of velocity etc.
! The fields to be ghosted is determined by which prime integers the input
! argument 'choice' is divisible by:
! 2: Vorticity
! 3: Velocity
! 5: Solid mask
! 7: Concentration
! 11: Vorticity RHS
! 13: Strain rate
!------------------------------------------------------------------------------
MODULE naga_mod_ghost_fields
INTERFACE naga_ghost_fields
  MODULE PROCEDURE naga_ghost_fields
END INTERFACE
CONTAINS
SUBROUTINE naga_ghost_fields(info,choice)
USE naga_mod_globals
USE naga_mod_say
USE naga_mod_extrapolate
USE ppm_module_map_field_ghost
USE ppm_module_map_field
USE ppm_module_data
IMPLICIT NONE
!-----------------------------------------------------------------------------
! Arguments
!-----------------------------------------------------------------------------
INTEGER, INTENT(OUT) :: info
INTEGER, INTENT(IN),OPTIONAL :: choice
!-----------------------------------------------------------------------------
! Local variables
!-----------------------------------------------------------------------------
REAL(MK) :: t0
INTEGER :: ilevel,ipatch
INTEGER :: toggle
INTEGER :: nbase,nextra
!-----------------------------------------------------------------------------
! Initialise routine
!-----------------------------------------------------------------------------
CALL substart('naga_ghost_fields',t0,info)
IF (PRESENT(choice)) THEN
  toggle = choice
ELSE
  toggle = 2
ENDIF
DO ilevel=1,nlevels
  DO ipatch=1,npatches(ilevel)
    CALL ppm_map_field_ghost_get(ptcset(ilevel,ipatch)%topoid,&
       & ptcset(ilevel,ipatch)%meshid,gstw,info)
    !vorticity
    IF (MOD(toggle,2) .EQ. 0) THEN
      CALL ppm_map_field_push(ptcset(ilevel,ipatch)%topoid,&
         & ptcset(ilevel,ipatch)%meshid,wf(ilevel,ipatch)%fld,&
         & ncom,info)
      IF (info .NE. 0) THEN
        CALL naga_say(rank,'naga_ghost_fields','Failed to push vrt.')
        GOTO 9999
      ENDIF
    ENDIF
    !velocity
    IF (MOD(toggle,3) .EQ. 0) THEN
      CALL ppm_map_field_push(ptcset(ilevel,ipatch)%topoid,&
         & ptcset(ilevel,ipatch)%meshid,uf(ilevel,ipatch)%fld,&
         & ncom,info)
      IF (info .NE. 0) THEN
        CALL naga_say(rank,'naga_ghost_fields','Failed to push vel.')
        GOTO 9999
      ENDIF
    ENDIF
    !solid mask
    IF (MOD(toggle,5) .EQ. 0) THEN
      CALL ppm_map_field_push(ptcset(ilevel,ipatch)%topoid,&
         & ptcset(ilevel,ipatch)%meshid,chif(ilevel,ipatch)%fld,&
         & info)
      IF (info .NE. 0) THEN
        CALL naga_say(rank,'naga_ghost_fields','Failed to push chi.')
        GOTO 9999
      ENDIF
    ENDIF
    !concentration
    IF (MOD(toggle,7) .EQ. 0) THEN
      CALL ppm_map_field_push(ptcset(ilevel,ipatch)%topoid,&
         & ptcset(ilevel,ipatch)%meshid,cf(ilevel,ipatch)%fld,&
         & info)
      IF (info .NE. 0) THEN
        CALL naga_say(rank,'naga_ghost_fields','Failed to push conc.')
        GOTO 9999
      ENDIF
    ENDIF
    !vorticity time derivative
    IF (MOD(toggle,11) .EQ. 0) THEN
      CALL ppm_map_field_push(ptcset(ilevel,ipatch)%topoid,&
         & ptcset(ilevel,ipatch)%meshid,dwf(ilevel,ipatch)%fld,&
         & ncom,info)
      IF (info .NE. 0) THEN
        CALL naga_say(rank,'naga_ghost_fields','Failed to push dvrt.')
        GOTO 9999
      ENDIF
    ENDIF
    !strain rate
    IF (MOD(toggle,13) .EQ. 0) THEN
      CALL ppm_map_field_push(ptcset(ilevel,ipatch)%topoid,&
         & ptcset(ilevel,ipatch)%meshid,sf(ilevel,ipatch)%fld,&
         & info)
      IF (info .NE. 0) THEN
        CALL naga_say(rank,'naga_ghost_fields','Failed to push sf.')
        GOTO 9999
      ENDIF
    ENDIF
    !---------------------------------------------------------------------------
    !SYMMETRY AROUND THESE LINES
    !---------------------------------------------------------------------------
    CALL ppm_map_field_send(info)
    !strain rate
    IF (MOD(toggle,13) .EQ. 0) THEN
      CALL ppm_map_field_pop(ptcset(ilevel,ipatch)%topoid,&
         & ptcset(ilevel,ipatch)%meshid,sf(ilevel,ipatch)%fld,&
         & gstw,info)
      IF (info .NE. 0) THEN
        CALL naga_say(rank,'naga_ghost_fields','Failed to pop sf.')
        GOTO 9999
      ENDIF
    ENDIF
    !vorticity time derivative
    IF (MOD(toggle,11) .EQ. 0) THEN
      CALL ppm_map_field_pop(ptcset(ilevel,ipatch)%topoid,&
         & ptcset(ilevel,ipatch)%meshid,dwf(ilevel,ipatch)%fld,&
         & ncom,gstw,info)
      IF (info .NE. 0) THEN
        CALL naga_say(rank,'naga_ghost_fields','Failed to pop dvrt.')
        GOTO 9999
      ENDIF
    ENDIF
    !concentration
    IF (MOD(toggle,7) .EQ. 0) THEN
      CALL ppm_map_field_pop(ptcset(ilevel,ipatch)%topoid,&
         & ptcset(ilevel,ipatch)%meshid,cf(ilevel,ipatch)%fld,&
         & gstw,info)
      IF (info .NE. 0) THEN
        CALL naga_say(rank,'naga_ghost_fields','Failed to pop conc.')
        GOTO 9999
      ENDIF
    ENDIF
    !solid mask
    IF (MOD(toggle,5) .EQ. 0) THEN
      CALL ppm_map_field_pop(ptcset(ilevel,ipatch)%topoid,&
         & ptcset(ilevel,ipatch)%meshid,chif(ilevel,ipatch)%fld,&
         & gstw,info)
      IF (info .NE. 0) THEN
        CALL naga_say(rank,'naga_ghost_fields','Failed to pop chi.')
        GOTO 9999
      ENDIF
    ENDIF
    !velocity
    IF (MOD(toggle,3) .EQ. 0) THEN
      CALL ppm_map_field_pop(ptcset(ilevel,ipatch)%topoid,&
         & ptcset(ilevel,ipatch)%meshid,uf(ilevel,ipatch)%fld,&
         & ncom,gstw,info)
      IF (info .NE. 0) THEN
        CALL naga_say(rank,'naga_ghost_fields','Failed to pop vel.')
        GOTO 9999
      ENDIF
    ENDIF
    !vorticity
    IF (MOD(toggle,2) .EQ. 0) THEN
      CALL ppm_map_field_pop(ptcset(ilevel,ipatch)%topoid,&
         & ptcset(ilevel,ipatch)%meshid,wf(ilevel,ipatch)%fld,&
         & ncom,gstw,info)
      IF (info .NE. 0) THEN
        CALL naga_say(rank,'naga_ghost_fields','Failed to pop vrt.')
        GOTO 9999
      ENDIF
    ENDIF
  ENDDO
ENDDO
IF (info .NE. 0) THEN
  CALL naga_say(rank,'naga_ghost_fields','Failed to ghost fields.')
  GOTO 9999
ENDIF
!-----------------------------------------------------------------------------
! Now call extrapolation routine if freespace
!-----------------------------------------------------------------------------
IF (domainbc .EQ. 0) THEN
  nbase = 4
  nextra = 2
  DO ilevel=1,nlevels
    DO ipatch=1,npatches(ilevel)
      IF (MOD(toggle,2) .EQ. 0) THEN
        CALL naga_extrapolate(ptcset(ilevel,ipatch),&
        & wf(ilevel,ipatch)%fld,nextra,nbase,info)
      ENDIF
      !velocity
      IF (MOD(toggle,3) .EQ. 0) THEN
        CALL naga_extrapolate(ptcset(ilevel,ipatch),&
        & uf(ilevel,ipatch)%fld,nextra,nbase,info)
      ENDIF
      !solid mask
      IF (MOD(toggle,5) .EQ. 0) THEN
        CALL naga_extrapolate(ptcset(ilevel,ipatch),&
        & chif(ilevel,ipatch)%fld,nextra,nbase,info)
      ENDIF
      !concentration
      IF (MOD(toggle,7) .EQ. 0) THEN
        CALL naga_extrapolate(ptcset(ilevel,ipatch),&
        & cf(ilevel,ipatch)%fld,nextra,nbase,info)
      ENDIF
      !vorticity time derivative
      IF (MOD(toggle,11) .EQ. 0) THEN
        CALL naga_extrapolate(ptcset(ilevel,ipatch),&
        & dwf(ilevel,ipatch)%fld,nextra,nbase,info)
      ENDIF
      !strain rate
      IF (MOD(toggle,13) .EQ. 0) THEN
        CALL naga_extrapolate(ptcset(ilevel,ipatch),&
        & sf(ilevel,ipatch)%fld,nextra,nbase,info)
      ENDIF
    ENDDO
  ENDDO
ENDIF
!-----------------------------------------------------------------------------
! Return
!-----------------------------------------------------------------------------
 9999 CONTINUE
CALL substop('naga_ghost_fields',t0,info)
RETURN
END SUBROUTINE naga_ghost_fields
END MODULE naga_mod_ghost_fields
