!------------------------------------------------------------------------------
! Subroutine :  naga_rhs.f90
!------------------------------------------------------------------------------ 
!      _  __
!     / |/ /__ ____ ____ _
!    /    / _ `/ _ `/ _ `/
!   /_/|_/\_,_/\_, /\_,_/
!             /___/
!   A PPM Vortex-In-Cell client - 2011 (c)
!   Naga is distributed under the GNU Lesser General Public License version 3
!   Johannes Tophøj Rasmussen, johannes.tophoj@gmail.com
!   Technical University of Denmark
!   Department of Mechanical Engineering
!
!
! This routines computes the rhs for the ode
!------------------------------------------------------------------------------

MODULE naga_mod_rhs

IMPLICIT NONE

INTERFACE naga_rhs
  MODULE PROCEDURE naga_rhs
END INTERFACE

CONTAINS

SUBROUTINE naga_rhs(istage,info)

USE naga_mod_globals
USE naga_mod_say
USE naga_mod_dump_vtk
USE naga_mod_remesh_particles
USE naga_mod_add_uinfinity
USE naga_mod_diagnostics
USE naga_mod_strainrate
USE naga_mod_save_fields
USE naga_mod_penalize_explicit
USE naga_mod_penalize_interpolation
USE naga_mod_penalize_interpolation_update
USE naga_mod_ghost_fields
USE naga_mod_vorticity_rhs
USE naga_mod_flowcase
USE ppm_module_interp_p2m
USE ppm_module_interp_m2p
USE ppm_module_poisson
#ifdef __DOTIMING
USE naga_mod_time
#endif

IMPLICIT NONE

!-------------------------------------------------------------------------
! Arguments
!-------------------------------------------------------------------------
INTEGER, INTENT(IN   )                                  :: istage
INTEGER, INTENT(  OUT)                                  :: info


!-------------------------------------------------------------------------
! Local variables
!-------------------------------------------------------------------------
REAL(MK)                        :: t0
INTEGER                         :: ilevel,ipatch
INTEGER                         :: choice
LOGICAL                         :: doremesh,doreprojection
LOGICAL                         :: dostrainrate,dodiagnostics


!----------------------------------------------------------------------------
! Initialise routine
!----------------------------------------------------------------------------
CALL substart('naga_rhs',t0,info)


!----------------------------------------------------------------------------
! Determine which actions to take during this stage/time step:
!----------------------------------------------------------------------------
! If istage=1 and we wish to do so then: Perhaps do some reprojection and if
! we do fancy it do a bit o' the good ol' remeshing. Actually always remesh 
! after reprojecting
! !@ if possible it would be preferable to also remesh dw and u - and by this
! I do not mean remesh as such but just to create ordered particle values 
! directly from the mesh points. PPM does not support this...
!----------------------------------------------------------------------------
doremesh       = .FALSE.
doreprojection = .FALSE.
dostrainrate   = .FALSE.
dodiagnostics  = .FALSE.
IF (istage .EQ. 1) THEN
  !At first RK step do diagnostics and thus also strain rate
  dostrainrate   = .TRUE.
  dodiagnostics  = .TRUE.
  IF (ireproject .NE. 0) THEN
    IF (MOD(itime,ireproject) .EQ. 0) THEN
      doreprojection = .TRUE.
      doremesh       = .TRUE.
    ENDIF
  ENDIF

  IF (iremesh .NE. 0) THEN
    IF (MOD(itime,iremesh) .EQ. 0) THEN
      doremesh       = .TRUE.
    ENDIF
  ENDIF
  !Always remesh if we penalize (interpolation)
  IF (penalization .EQ. 1) THEN
    doremesh       = .TRUE.
  ENDIF
ENDIF

IF (lesmodel .EQ. 1) THEN
  dostrainrate   = .TRUE.
ENDIF

!----------------------------------------------------------------------------
! Flowcase specific actions in the beginning of the rhs schedule (pos #1)
!----------------------------------------------------------------------------
CALL naga_flowcase(1,info)
IF (info .NE. 0) THEN
  CALL naga_say(rank,'naga_rhs','Failed to do flowcase specific actions #1.')
  GOTO 9999
ENDIF

#ifdef __DOTIMING
!Example of timer start
CALL naga_time_str(1)
#endif
!----------------------------------------------------------------------------
! Interpolate particle vorticity to mesh
! ppm_map_fields_ghost_put is being called automatically by ppm_interp_p2m
!----------------------------------------------------------------------------
DO ilevel=1,nlevels
  DO ipatch=1,npatches(ilevel)
    wf(ilevel,ipatch)%fld = 0.0_MK !@
    CALL ppm_interp_p2m(ptcset(ilevel,ipatch)%topoid,&
                      & ptcset(ilevel,ipatch)%meshid,&
                      & xp(ilevel,ipatch)%val,&
                      & ptcset(ilevel,ipatch)%np,&
                      & wp(ilevel,ipatch)%val,&
                      & ncom,&
                      & ppm_param_rmsh_kernel_mp4,&
                      & gstw,&
                      & wf(ilevel,ipatch)%fld,&
                      & info)
    IF (info .NE. 0) THEN
      CALL naga_say(rank,'naga_rhs','Failed to interpolate vorticity to mesh.')
      GOTO 9999
    ENDIF
  ENDDO !ipatch
ENDDO !ilevel
#ifdef __DOTIMING
!Example of timer stop
CALL naga_time_stp(1)
#endif

!@ In time reprojection could be moved here

!----------------------------------------------------------------------------
! Interpolate particle concentration to mesh
!----------------------------------------------------------------------------
IF (concentration .EQ. 1) THEN
  DO ilevel=1,nlevels
    DO ipatch=1,npatches(ilevel)
      CALL ppm_interp_p2m(ptcset(ilevel,ipatch)%topoid,&
                        & ptcset(ilevel,ipatch)%meshid,&
                        & xp(ilevel,ipatch)%val,&
                        & ptcset(ilevel,ipatch)%np,&
                        & cp(ilevel,ipatch)%val,&
                        & ppm_param_rmsh_kernel_mp4,&
                        & gstw,&
                        & cf(ilevel,ipatch)%fld,&
                        & info)
      IF (info .NE. 0) THEN
        CALL naga_say(rank,'naga_rhs','Failed to interpolate vorticity to mesh.')
        GOTO 9999
      ENDIF
    ENDDO !ipatch
  ENDDO !ilevel
ENDIF



!----------------------------------------------------------------------------
! Calculate velocity on mesh
! Vorticity has not been ghosted but we dont need it for this
!----------------------------------------------------------------------------
DO ilevel=1,nlevels
  DO ipatch=1,npatches(ilevel)
    CALL ppm_poisson_solve(ptcset(ilevel,ipatch)%topoid,&
       & ptcset(ilevel,ipatch)%meshid,ppmpoisson(ilevel,ipatch),&
       & wf(ilevel,ipatch)%fld,uf(ilevel,ipatch)%fld,gstw,info)
    IF (info .NE. 0) THEN
      CALL naga_say(rank,'naga_rhs','Failed to compute velocity.')
      GOTO 9999
    ENDIF
  ENDDO
ENDDO
CALL naga_add_uinfinity(info)
IF (info .NE. 0) THEN
  CALL naga_say(rank,'naga_rhs','Failed to add free stream velocity.')
  GOTO 9999
ENDIF


!----------------------------------------------------------------------------
! Compute strain rate for use in diagnostics and/or LES model
!----------------------------------------------------------------------------
IF (dostrainrate) THEN
  CALL naga_strainrate(info)
  IF (info .NE. 0) THEN
    CALL naga_say(rank,'naga_rhs','Failed to compute strainrate.')
    GOTO 9999
  ENDIF
ENDIF

!----------------------------------------------------------------------------
! Do diagnostics - velocity and vorticity exist in their virgin format. 
! Velocity has been ghosted in ppm_poisson_solve. Vorticity has NOT been ghosted
! Also this needs to be done before explicit penalization
!----------------------------------------------------------------------------
IF (dodiagnostics) THEN
  CALL naga_diagnostics(info)
  IF (info .NE. 0) THEN
    CALL naga_say(rank,'naga_rhs','Failed to compute diagnostics.')
    GOTO 9999
  ENDIF
ENDIF


IF (istage .EQ. 1) THEN
  !----------------------------------------------------------------------------
  ! Dump field values to disk - vorticity, velocity, mask and concentration
  !----------------------------------------------------------------------------
  choice = 1
  choice=choice*2
  choice=choice*3
  choice=choice*5
  choice=choice*7
  CALL naga_save_fields(info,choice)
  IF (info .NE. 0) THEN
    CALL naga_say(rank, 'naga_integrate_rhs','Failed to save fields.')
    GOTO 9999
  END IF
ENDIF


!----------------------------------------------------------------------------
! If istage=1 do penalization (including additional velocity computation)
!----------------------------------------------------------------------------
IF (istage .EQ. 1 .AND. penalization .EQ. 1) THEN
  CALL naga_penalize_interpolation_update(info)

  !----------------------------------------------------------------------------
  ! Recalculate velocity
  ! The vorticity ghosts are not updated after the penalization but are not 
  ! required either.
  !----------------------------------------------------------------------------
  DO ilevel=1,nlevels
    DO ipatch=1,npatches(ilevel)
      CALL ppm_poisson_solve(ptcset(ilevel,ipatch)%topoid,&
      & ptcset(ilevel,ipatch)%meshid,ppmpoisson(ilevel,ipatch),&
      & wf(ilevel,ipatch)%fld,uf(ilevel,ipatch)%fld,gstw,info)
      IF (info .NE. 0) THEN
        CALL naga_say(rank,'naga_rhs','Failed to compute velocity.')
        GOTO 9999
      ENDIF
    ENDDO
  ENDDO
  CALL naga_add_uinfinity(info)
  IF (info .NE. 0) THEN
    CALL naga_say(rank,'naga_rhs','Failed to add free stream velocity.')
    GOTO 9999
  ENDIF
ENDIF


!----------------------------------------------------------------------------
! Perform vorticity reprojection to kill divergence of the field
! Vorticity ghosts have not been updated here (@) either - thats ok for now
! as the poisson routines handles the necessary extrapolation of vorticity.
!----------------------------------------------------------------------------
IF (doreprojection) THEN
!CALL naga_ghost_fields(info,choice=2) !@
  DO ilevel=1,nlevels
    DO ipatch=1,npatches(ilevel) 
      !CALL ppm_poisson_solve(ptcset(ilevel,ipatch)%topoid,&
      !& ptcset(ilevel,ipatch)%meshid,ppmpoisson(ilevel,ipatch),&
      !& wf(ilevel,ipatch)%fld,uf(ilevel,ipatch)%fld,gstw,info,&
      !& tmpoperation = ppm_poisson_opr_repr)

      ! to call the newer v1.2.2 Poisson solver
      !ppm_poisson_solve(topoid,
      !                  meshid,
      !                  ppmpoisson,
      !                  fieldin,
      !                  fieldout,
      !                  gstw,
      !                  info,
      !                  tmpcase)
      CALL ppm_poisson_solve(ptcset(ilevel,ipatch)%topoid,&
                            &ptcset(ilevel,ipatch)%meshid,&
                            &ppmpoisson(ilevel,ipatch),&
                            &wf(ilevel,ipatch)%fld,&
                            &uf(ilevel,ipatch)%fld,&
                            &gstw,&
                            &info,&
                            &ppm_poisson_grn_reprojec)
    ENDDO
  ENDDO
ENDIF


!----------------------------------------------------------------------------
! Remesh particles. Velocity and vorticity (and concentration) remains on mesh
!----------------------------------------------------------------------------
IF (doremesh) THEN
  DO ilevel=1,nlevels
    DO ipatch=1,npatches(ilevel)
      IF (concentration .NE. 1) THEN
        CALL naga_remesh_particles(info)
      ELSE
        ! Concentration field present
        CALL naga_remesh_particles(cf,cp,info)
      ENDIF
      IF (info .NE. 0) THEN
        CALL naga_say(rank, 'naga_setup','Failed to create vorticity particles.')
        GOTO 9999
      ENDIF
    ENDDO
  ENDDO
ENDIF



!-------------------------------------------------------------------------------
! No matter what, calculate vorticity rhs here
! - Ghost vorticity first
! - also strain rate if we do LES
! - and afterwards domega/dt as well
! naga_ghost_fields also extrapolates field values to ghost layer if freespace
!-------------------------------------------------------------------------------
choice = 2
IF (lesmodel .EQ. 1) THEN
  choice = choice * 13
ENDIF
CALL naga_ghost_fields(info,choice=choice)
CALL naga_vorticity_rhs(info)
!----------------------------------------------------------------------------
! If we do explicit penalization then this is the place
!----------------------------------------------------------------------------
IF (penalization .EQ. 2) THEN
  CALL naga_penalize_explicit(info)
END IF
CALL naga_ghost_fields(info,choice=11)



!----------------------------------------------------------------------------
! Dump field values to disk - vorticity derivative
!----------------------------------------------------------------------------
IF (istage .EQ. 1) THEN
  choice = 1
  choice=choice*11
  CALL naga_save_fields(info,choice)
  IF (info .NE. 0) THEN
    CALL naga_say(rank, 'naga_integrate_rhs','Failed to save fields.')
    GOTO 9999
  END IF
END IF


!----------------------------------------------------------------------------
! Interpolate RHS and velocities to particles
! !@ in time this should only be done when remeshing is not done.
!----------------------------------------------------------------------------
!IF ((istage .NE. 1) .OR. (iremesh .EQ. 0) .OR. (MOD(itime,iremesh) .NE. 0)) THEN
DO ilevel=1,nlevels
  DO ipatch=1,npatches(ilevel)
    CALL ppm_interp_m2p(ptcset(ilevel,ipatch)%topoid,&
                      & ptcset(ilevel,ipatch)%meshid,&
                      & xp(ilevel,ipatch)%val,&
                      & ptcset(ilevel,ipatch)%np,&
                      & dwp(ilevel,ipatch)%val,&
                      & ncom,&
                      & ppm_param_rmsh_kernel_mp4,&
                      & gstw,&
                      & dwf(ilevel,ipatch)%fld,&
                      & info)
    CALL ppm_interp_m2p(ptcset(ilevel,ipatch)%topoid,&
                      & ptcset(ilevel,ipatch)%meshid,&
                      & xp(ilevel,ipatch)%val,&
                      & ptcset(ilevel,ipatch)%np,&
                      & up(ilevel,ipatch)%val,&
                      & ncom,&
                      & ppm_param_rmsh_kernel_mp4,&
                      & gstw,&
                      & uf(ilevel,ipatch)%fld,&
                      & info)
    IF (info .NE. 0) THEN
      CALL naga_say(rank,'naga_rhs','Failed to interpolate time derivatives to particles.')
      GOTO 9999
    ENDIF
  ENDDO !ipatch
ENDDO !ilevel
!ENDIF


!----------------------------------------------------------------------------
! Return
!----------------------------------------------------------------------------
 9999 CONTINUE
CALL substop('naga_rhs',t0,info)
RETURN


END SUBROUTINE naga_rhs

END MODULE naga_mod_rhs

