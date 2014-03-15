!------------------------------------------------------------------------------
! Subroutine :  Subroutine :  naga_setup.f90
!------------------------------------------------------------------------------
!      _  __
!     / |/ /__ ____ ____ _
!    /    / _ `/ _ `/ _ `/
!   /_/|_/\_,_/\_, /\_,_/
!             /___/

!   A PPM Vortex-In-Cell client Vortex-In-Cell client - 2011 (c)
!   Naga is distributed under the GNU Lesser General Public License version 3
!   Johannes Tophøj Rasmussen, johannes.tophoj@gmail.com
!   Technical University of Denmark
!   Department of Mechanical Engineering
!
!
! This routines declares fields, particles etc and sets up the initial problem
!------------------------------------------------------------------------------

MODULE naga_mod_setup

IMPLICIT NONE

INTERFACE naga_setup
  MODULE PROCEDURE naga_setup
END INTERFACE

CONTAINS

SUBROUTINE naga_setup(info)

USE naga_mod_globals
USE naga_mod_say
USE naga_mod_allocate_patch
USE naga_mod_case_test
USE naga_mod_case_testreprojection
USE naga_mod_case_testperiodic
USE naga_mod_case_testring
USE naga_mod_case_inittaylorgreen
USE naga_mod_case_zero
USE naga_mod_case_torus
USE naga_mod_setup_topologies
USE naga_mod_get_topologies
USE naga_mod_copy_patch
USE naga_mod_remesh_particles
USE naga_mod_structure_sphere
USE naga_mod_structure_none
USE naga_mod_concentration_sin
USE naga_mod_stl
USE naga_mod_stepfunction
USE ppm_module_poisson

IMPLICIT NONE

!-------------------------------------------------------------------------
! Arguments
!-------------------------------------------------------------------------
INTEGER, INTENT(INOUT) :: info
!-------------------------------------------------------------------------
! Local variables
!-------------------------------------------------------------------------
REAL(MK)                        :: t0
INTEGER                         :: ilevel,ipatch
INTEGER                         :: poissonoperation
INTEGER                         :: poissonderivative
INTEGER                         :: poissonbc
INTEGER                         :: poissongreen
CHARACTER(len=256)              :: filename
LOGICAL                         :: file_exists


!----------------------------------------------------------------------------
! Initialise routine
!----------------------------------------------------------------------------
CALL substart('naga_setup',t0,info)


!-------------------------------------------------------------------------
! Delete diagnostics file and ABORT file if it already exists
! (diagnostics file should not be deleted when restarting (is implemented))
!-------------------------------------------------------------------------
IF (rank .EQ. 0) THEN
  WRITE(filename,'(A,A)') runtag(1:iruntag), '.diag'
  INQUIRE(FILE=filename,EXIST=file_exists)
  IF (file_exists) THEN
    OPEN(10,FILE=filename)
    CLOSE(10,STATUS='delete')
  ENDIF
  WRITE(filename,'(A)') 'ABORT'
  INQUIRE(FILE=filename,EXIST=file_exists)
  IF (file_exists) THEN
    OPEN(10,FILE=filename)
    CLOSE(10,STATUS='delete')
  ENDIF
ENDIF

!-------------------------------------------------------------------------
! Maybe this is the place to read restart files
!-------------------------------------------------------------------------


!-------------------------------------------------------------------------
! Setup a topology for each level and store the topologies in a vector
!-------------------------------------------------------------------------
CALL naga_setup_topologies(ptcset,ppm_param_decomp_cartesian,gstw,info)
CALL naga_get_topologies(ptcset,topos,info)


!----------------------------------------------------------------------------
! Now allocate the fields
!----------------------------------------------------------------------------
CALL naga_allocate_patch(ptcset,wf,info)
IF (info .NE. 0) THEN
  CALL naga_say(rank, 'naga_setup','Failed to allocate wf.')
  GOTO 9999
ENDIF
CALL naga_allocate_patch(ptcset,dwf,info)
IF (info .NE. 0) THEN
  CALL naga_say(rank, 'naga_setup','Failed to allocate dwf.')
  GOTO 9999
ENDIF
CALL naga_allocate_patch(ptcset,uf,info)
IF (info .NE. 0) THEN
  CALL naga_say(rank, 'naga_setup','Failed to allocate uf.')
  GOTO 9999
ENDIF
IF (penalization .GT. 0) THEN
  CALL naga_allocate_patch(ptcset,chif,info)
  IF (info .NE. 0) THEN
    CALL naga_say(rank, 'naga_setup','Failed to allocate chif.')
    GOTO 9999
  ENDIF
  CALL naga_allocate_patch(ptcset,ubarf,info)
  IF (info .NE. 0) THEN
    CALL naga_say(rank, 'naga_setup','Failed to allocate ubarf.')
    GOTO 9999
  ENDIF
  !----------------------------------------------------------------------------
  ! Set mask and ubar 0 - just in case.
  !----------------------------------------------------------------------------
  CALL naga_structure_none(info)
  IF (info .NE. 0) THEN
    CALL naga_say(rank, 'naga_setup','Failed to initialise 0 chi/ubar.')
    GOTO 9999
  ENDIF
  !----------------------------------------------------------------------------
  ! Initialise step function (stepfunc1)
  !----------------------------------------------------------------------------
  CALL naga_init_stepfunc1(info)
ENDIF
!If toggled a field for validation purposes is allocated
IF (validationfield .EQ. 1) THEN
  CALL naga_allocate_patch(ptcset,psif,info)
  IF (info .NE. 0) THEN
    CALL naga_say(rank, 'naga_setup','Failed to allocate psif.')
    GOTO 9999
  ENDIF
ENDIF

!If toggled, allocate concentration field patches
IF (concentration .EQ. 1) THEN
  CALL naga_allocate_patch(ptcset,cf,info)
  IF (info .NE. 0) THEN
    CALL naga_say(rank, 'naga_setup','Failed to allocate cf.')
    GOTO 9999
  ENDIF
ENDIF

!If toggled, allocate concentration field patches
IF (lesmodel .EQ. 1) THEN
  CALL naga_allocate_patch(ptcset,sf,info)
  IF (info .NE. 0) THEN
    CALL naga_say(rank, 'naga_setup','Failed to allocate sf.')
    GOTO 9999
  ENDIF
ENDIF


!----------------------------------------------------------------------------
! Allocate particle patch structs (the data arrays are allocated elsewhere)
!----------------------------------------------------------------------------
IF (timeintscheme .EQ. 1  .OR. &
  & timeintscheme .EQ. 2       ) THEN
  CALL naga_allocate_patch(ptcset,xp,info)
  IF (info .NE. 0) THEN
    CALL naga_say(rank, 'naga_setup','Failed to allocate xp.')
    GOTO 9999
  ENDIF
  CALL naga_allocate_patch(ptcset,wp,info)
  IF (info .NE. 0) THEN
    CALL naga_say(rank, 'naga_setup','Failed to allocate wp.')
    GOTO 9999
  ENDIF
  CALL naga_allocate_patch(ptcset,up,info)
  IF (info .NE. 0) THEN
    CALL naga_say(rank, 'naga_setup','Failed to allocate up.')
    GOTO 9999
  ENDIF
  CALL naga_allocate_patch(ptcset,dwp,info)
  IF (info .NE. 0) THEN
    CALL naga_say(rank, 'naga_setup','Failed to allocate dwp.')
    GOTO 9999
  ENDIF
ENDIF

IF (timeintscheme .EQ. 2) THEN
  CALL naga_allocate_patch(ptcset,xp0,info)
  IF (info .NE. 0) THEN
    CALL naga_say(rank, 'naga_setup','Failed to allocate xp.')
    GOTO 9999
  ENDIF
  CALL naga_allocate_patch(ptcset,wp0,info)
  IF (info .NE. 0) THEN
    CALL naga_say(rank, 'naga_setup','Failed to allocate wp.')
    GOTO 9999
  ENDIF
  CALL naga_allocate_patch(ptcset,up0,info)
  IF (info .NE. 0) THEN
    CALL naga_say(rank, 'naga_setup','Failed to allocate up.')
    GOTO 9999
  ENDIF
  CALL naga_allocate_patch(ptcset,dwp0,info)
  IF (info .NE. 0) THEN
    CALL naga_say(rank, 'naga_setup','Failed to allocate dwp.')
    GOTO 9999
  ENDIF
ENDIF

!If toggled allocate concentration particle patches
IF (concentration .EQ. 1) THEN
  CALL naga_allocate_patch(ptcset,cp,info)
  IF (info .NE. 0) THEN
    CALL naga_say(rank, 'naga_setup','Failed to allocate cp.')
    GOTO 9999
  ENDIF
ENDIF


!----------------------------------------------------------------------------
! Determine velocity derivatives and other settings for poisson solver
! (these settings may be changed flow-specifically)
!----------------------------------------------------------------------------
IF (velocityscheme .EQ. 0) THEN
  poissonoperation  = ppm_poisson_opr_none
  poissonderivative = ppm_poisson_drv_sp
ELSE IF (velocityscheme .EQ. 1) THEN
  poissonoperation  = ppm_poisson_opr_vel
  poissonderivative = ppm_poisson_drv_fd2
ELSE IF (velocityscheme .EQ. 2) THEN
  poissonoperation  = ppm_poisson_opr_vel
  poissonderivative = ppm_poisson_drv_fd4
ELSE IF (velocityscheme .EQ. 3) THEN
  poissonoperation  = ppm_poisson_opr_vel
  poissonderivative = ppm_poisson_drv_sp
ELSE
  CALL naga_say(rank,'Naga_setup','Undefined velocity scheme specified.')
  info = -1
  GOTO 9999
ENDIF

IF (domainbc .EQ. 0) THEN
  poissonbc    = ppm_poisson_bc_fre
  IF (freespacekernel .EQ. 0) THEN
    poissongreen = ppm_poisson_grn_pois_per
  ELSE IF (freespacekernel .EQ. 1) THEN
    poissongreen = ppm_poisson_grn_pois_fre
  ELSE IF (freespacekernel .EQ. 2) THEN
    poissongreen = ppm_poisson_grn_pois_blob2
  ELSE IF (freespacekernel .EQ. 3) THEN
    poissongreen = ppm_poisson_grn_pois_blob4
  ELSE IF (freespacekernel .EQ. 4) THEN
    poissongreen = ppm_poisson_grn_pois_blob6
  ELSE IF (freespacekernel .EQ. 5) THEN
    poissongreen = ppm_poisson_grn_pois_blob8
  ELSE IF (freespacekernel .EQ. 6) THEN
    poissongreen = ppm_poisson_grn_pois_blob10
  ELSE
    CALL naga_say(rank,'Naga_setup','Undefined Greens function specified.')
    info = -1
    GOTO 9999
  ENDIF
ELSE IF (domainbc .EQ. 1) THEN
  poissonbc    = ppm_poisson_bc_per
  poissongreen = ppm_poisson_grn_pois_per
ENDIF


!----------------------------------------------------------------------------
! Initialise the flow
!----------------------------------------------------------------------------
IF (.NOT. restarted) THEN
  SELECT CASE(flowcase)
  CASE(0)
    !-------------------------------------------------------------------------
    ! Dummy flow case
    !-------------------------------------------------------------------------
    IF (rank .EQ. 0) THEN
      CALL naga_say(rank,'naga_setup','No flowcase selected')
    ENDIF
    info = -1
    GOTO 9999
  CASE(1)
    !-------------------------------------------------------------------------
    ! Test flow case
    !-------------------------------------------------------------------------
    IF (rank .EQ. 0) THEN
      CALL naga_say(rank,'naga_setup','Flowcase 1 selected')
    ENDIF
    CALL naga_say(rank,'Naga_setup','Flowcase 1 not implemented.')
    STOP
  CASE(2)
    !-------------------------------------------------------------------------
    ! Vortex ring
    !-------------------------------------------------------------------------
    IF (rank .EQ. 0) THEN
      CALL naga_say(rank,'naga_setup','Flowcase 2 selected: Vortex ring')
    ENDIF
    CALL naga_case_torus(info)
    IF (info .NE. 0) THEN
      CALL naga_say(rank,'Naga_setup','Failed to setup flowcase 2.')
      GOTO 9999
    ENDIF
  CASE(3)
    !-------------------------------------------------------------------------
    ! Sphere
    !-------------------------------------------------------------------------
    IF (rank .EQ. 0) THEN
      CALL naga_say(rank,'naga_setup','Flowcase 3: sphere')
    ENDIF
    CALL naga_case_zero(info)
    IF (info .NE. 0) THEN
      CALL naga_say(rank,'Naga_setup','Failed to setup flowcase 3.')
      GOTO 9999
    ENDIF
    CALL naga_structure_sphere(info)
    IF (info .NE. 0) THEN
      CALL naga_say(rank,'Naga_setup','Failed to setup mask of static sphere.')
      GOTO 9999
    ENDIF
  CASE(31)
    !-------------------------------------------------------------------------
    ! Concentration
    !-------------------------------------------------------------------------
    IF (rank .EQ. 0) THEN
      CALL naga_say(rank,'naga_setup','Flowcase 31 selected: concentration')
    ENDIF
    CALL naga_case_zero(info)
    CALL naga_structure_sphere(info)
    CALL naga_concentration_sin(info)
    IF (info .NE. 0) THEN
      CALL naga_say(rank,'Naga_setup','Failed to setup flowcase 31.')
      GOTO 9999
    ENDIF
  CASE(32)
    !-------------------------------------------------------------------------
    ! Sphere with disruption of y component of free stream velocity (ploumhans)
    !-------------------------------------------------------------------------
    IF (rank .EQ. 0) THEN
      CALL naga_say(rank,'naga_setup','Flowcase 32: sphere with disruption')
    ENDIF
    CALL naga_case_zero(info)
    IF (info .NE. 0) THEN
      CALL naga_say(rank,'Naga_setup','Failed to setup flowcase 32.')
      GOTO 9999
    ENDIF
    CALL naga_structure_sphere(info)
    IF (info .NE. 0) THEN
      CALL naga_say(rank,'Naga_setup','Failed to setup mask of static sphere.')
      GOTO 9999
    ENDIF
  CASE(5)
    !-------------------------------------------------------------------------
    ! STL input
    !-------------------------------------------------------------------------
    IF (rank .EQ. 0) THEN
      CALL naga_say(rank,'naga_setup','Flowcase 5 selected: STL input')
    ENDIF
    CALL naga_case_zero(info)
    CALL naga_init_stepfunc1(info)
    IF (info .NE. 0) THEN
      CALL naga_say(rank,'Naga_setup','Failed to initialise step function 1.')
      GOTO 9999
    ENDIF
    ALLOCATE(stl(1))
    ALLOCATE(stlt(1))
    CALL naga_stl_read(info)
    IF (info .NE. 0) THEN
      CALL naga_say(rank,'Naga_setup','Failed to setup flowcase 5.')
      GOTO 9999
    ENDIF
    CALL naga_stl_init(info)
    IF (info .NE. 0) THEN
      CALL naga_say(rank,'Naga_setup','Failed to setup flowcase 5.')
      GOTO 9999
    ENDIF
    DEALLOCATE(stlt)
  CASE(6)
    !-------------------------------------------------------------------------
    ! Taylor Green vortices
    !-------------------------------------------------------------------------
    IF (rank .EQ. 0) THEN
      CALL naga_say(rank,'naga_setup','Flowcase 6 selected: Taylor-Green vortices')
    ENDIF
    CALL naga_case_inittaylorgreen(info)
    IF (info .NE. 0) THEN
      CALL naga_say(rank,'Naga_setup','Failed to setup flowcase 6.')
      GOTO 9999
    ENDIF
  CASE(301)
    !-------------------------------------------------------------------------
    ! Test flow case
    !-------------------------------------------------------------------------
    validation = .TRUE.
    IF (rank .EQ. 0) THEN
      CALL naga_say(rank,'naga_setup','Flowcase 301 selected: test of periodic Greens function')
    ENDIF
    CALL naga_case_testperiodic(info)
    IF (info .NE. 0) THEN
      CALL naga_say(rank,'Naga_setup','Failed to setup flowcase 301.')
      GOTO 9999
    ENDIF
  CASE(302)
    !-------------------------------------------------------------------------
    ! Vortex ring test case for stream function
    !-------------------------------------------------------------------------
    validation = .TRUE.
    IF (rank .EQ. 0) THEN
      CALL naga_say(rank,'naga_setup','Flowcase 302 selected: Vortex ring validation')
    ENDIF
    CALL naga_case_testring(info)
    IF (info .NE. 0) THEN
      CALL naga_say(rank,'Naga_setup','Failed to setup flowcase 302.')
      GOTO 9999
    ENDIF
  CASE(303)
    !-------------------------------------------------------------------------
    ! Test periodic reprojection
    !-------------------------------------------------------------------------
    validation = .TRUE.
    IF (rank .EQ. 0) THEN
      CALL naga_say(rank,'naga_setup','Flowcase 303 selected: test periodic reprojection')
    ENDIF
    !CALL naga_case_inittaylorgreen(info)
    CALL naga_case_zero(info)
    IF (info .NE. 0) THEN
      CALL naga_say(rank,'Naga_setup','Failed to setup flowcase 303.')
      GOTO 9999
    ENDIF
    CALL naga_case_testreprojection(info)
    IF (info .NE. 0) THEN
      CALL naga_say(rank,'Naga_setup','Failed to add reprojection test field.')
      GOTO 9999
    ENDIF
  CASE(304)
    !-------------------------------------------------------------------------
    ! Vortex freespace reprojection
    !-------------------------------------------------------------------------
    validation = .TRUE.
    IF (rank .EQ. 0) THEN
      CALL naga_say(rank,'naga_setup','Flowcase 304 selected: test freespace reprojection')
    ENDIF
    CALL naga_case_testring(info)
    !CALL naga_case_zero(info)
    IF (info .NE. 0) THEN
      CALL naga_say(rank,'Naga_setup','Failed to setup flowcase 304.')
      GOTO 9999
    ENDIF
    CALL naga_case_testreprojection(info)
    IF (info .NE. 0) THEN
      CALL naga_say(rank,'Naga_setup','Failed to add reprojection test field.')
      GOTO 9999
    ENDIF
  END SELECT
END IF


!----------------------------------------------------------------------------
! Initialize Poisson plan
!----------------------------------------------------------------------------
ALLOCATE(ppmpoisson(nlevels,maxpatches))
DO ilevel=1,nlevels
  ipatch=1
  CALL ppm_poisson_init(ptcset(ilevel,ipatch)%topoid,&
       &ptcset(ilevel,ipatch)%meshid,ppmpoisson(ilevel,ipatch),&
       &wf(ilevel,ipatch)%fld,uf(ilevel,ipatch)%fld,&
       &poissongreen,poissonbc,poissonoperation,poissonderivative,info)
  IF (info .NE. 0) THEN
    CALL naga_say(rank,'Naga_setup','Failed to initialise poisson routine.')
    GOTO 9999
  ENDIF
END DO


!----------------------------------------------------------------------------
! Create particles
!----------------------------------------------------------------------------
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


!----------------------------------------------------------------------------
! Return
!----------------------------------------------------------------------------
 9999 CONTINUE
CALL substop('naga_setup',t0,info)
RETURN


END SUBROUTINE naga_setup

END MODULE naga_mod_setup

