!-------------------------------------------------------------------------------
! Subroutine : Subroutine : naga_setup_topologies.f90
!-------------------------------------------------------------------------------
! _ __
! / |/ /__ ____ ____ _
! / / _ `/ _ `/ _ `/
! /_/|_/\_,_/\_, /\_,_/
! /___/
! A PPM Vortex-In-Cell client Vortex-In-Cell client - 2011 (c)
! Naga is distributed under the GNU Lesser General Public License version 3
! Johannes Tophøj Rasmussen, johannes.tophoj@gmail.com
! Technical University of Denmark
! Department of Mechanical Engineering
!
!
! This routines creates the topologies of the patches and obtains the size of
! the decomposed parallel domains
!-------------------------------------------------------------------------------
MODULE naga_mod_setup_topologies
  IMPLICIT NONE
  INTERFACE naga_setup_topologies
    MODULE PROCEDURE naga_setup_topologies
  END INTERFACE
  CONTAINS
SUBROUTINE naga_setup_topologies(psetting,decomposition,ghostwidth,info)
USE naga_mod_globals
USE naga_mod_say
USE ppm_module_mktopo
USE ppm_module_topo_get
IMPLICIT NONE
!-----------------------------------------------------------------------------
! Arguments
!-----------------------------------------------------------------------------
TYPE(patch_setup),DIMENSION(:,:),POINTER :: psetting
INTEGER, INTENT(IN) :: decomposition
INTEGER, DIMENSION(:), INTENT(IN) :: ghostwidth
INTEGER, INTENT(INOUT) :: info
!-----------------------------------------------------------------------------
! Local variables
!-----------------------------------------------------------------------------
REAL(MK) :: t0
INTEGER :: ilevel
INTEGER :: meshid
INTEGER,SAVE :: topoid !must have the SAVE
INTEGER, DIMENSION(ndim*2) :: bcdef !boundary conditions
INTEGER :: assigning
INTEGER :: ipatch !1 until supported by topo
REAL(MK),DIMENSION(:,:),POINTER :: tmpxp=>NULL()
INTEGER, DIMENSION(ndim) :: maxsnx !max extent of all subs
INTEGER,DIMENSION(:),POINTER :: dummynm=>NULL()
INTEGER,DIMENSION(:),POINTER :: dummyisublist=>NULL()
INTEGER :: dummynsublist
!-----------------------------------------------------------------------------
! Initialise routine
!-----------------------------------------------------------------------------
CALL substart('naga_setup_topologies',t0,info)
!-----------------------------------------------------------------------------
! Loop over all resolution levels
!-----------------------------------------------------------------------------
DO ilevel=1,nlevels
  ipatch = 1
  !---------------------------------------------------------------------------
  ! Get ready to create new topology:
  !---------------------------------------------------------------------------
  topoid = 0
  meshid = -1
  IF (domainbc .EQ. 0) THEN
    bcdef = ppm_param_bcdef_freespace
  ELSE IF (domainbc .EQ.1) THEN
    bcdef = ppm_param_bcdef_periodic
  ELSE
    CALL naga_say(rank, 'naga_setup_topologies','Invalid boundary condition.')
    GOTO 9999
  ENDIF
  assigning = ppm_param_assign_internal
  !---------------------------------------------------------------------------
  ! Create new topology
  ! This requires:
  ! - an index for the topology
  ! - an index for the mesh definition
  ! - a pointer to the particle positions
  ! - the number of particles. (= 0) as the domain decomposition is not guided
  ! - the decomposition type
  ! - assignement of the domain to the computational nodes (?)
  ! * minimum of physical domain
  ! * maximum of physical domain
  ! - boundary conditions for the domain
  ! !@!should be modified to freespace/semiperiodic
  ! - computational cost of the subdomain
  ! !@should be stored in a patch like structure
  ! - starting index of subdomain in the global mesh frame (istart)
  ! !@should be stored in a patch like structure
  ! - number of grid points on sub mesh
  ! !@should be stored in a patch like structure
  ! * number of grid points on the global mesh
  !---------------------------------------------------------------------------
  CALL ppm_mktopo(topoid,meshid,tmpxp,0,decomposition,assigning,&
  & psetting(ilevel,ipatch)%min, psetting(ilevel,ipatch)%max, &
  & bcdef,ghostwidth, psetting(ilevel,ipatch)%scost,&
  & psetting(ilevel,ipatch)%gnx,info)
  IF (info .NE. 0) THEN
    CALL naga_say(rank, 'naga_setup_topologies','Failed to create topology.')
    GOTO 9999
  ENDIF
  DO ipatch=1,npatches(ilevel)
    psetting(ilevel,ipatch)%topoid = topoid
    psetting(ilevel,ipatch)%meshid = meshid
    psetting(ilevel,ipatch)%gst(1) = ghostwidth(1)
    psetting(ilevel,ipatch)%gst(2) = ghostwidth(2)
    psetting(ilevel,ipatch)%gst(3) = ghostwidth(3)
  ENDDO
  !-------------------------------------------------------------------------
  ! Get additional mesh information
  !-------------------------------------------------------------------------
  ipatch = 1!@
  NULLIFY(dummynm)
  NULLIFY(psetting(ilevel,ipatch)%sistr)
  NULLIFY(psetting(ilevel,ipatch)%snx)
  NULLIFY(dummyisublist)
  CALL ppm_topo_get_meshinfo(psetting(ilevel,ipatch)%topoid, &
                           & psetting(ilevel,ipatch)%meshid,&
                           & dummynm,&
                           & psetting(ilevel,ipatch)%sistr,&
                           & psetting(ilevel,ipatch)%snx,&
                           & maxsnx, &
                           & dummyisublist, &
                           & dummynsublist, &
                           & info)
ENDDO
!----------------------------------------------------------------------------
! Return
!----------------------------------------------------------------------------
 9999 CONTINUE
CALL substop('naga_setup_topologies',t0,info)
RETURN
END SUBROUTINE naga_setup_topologies
END MODULE naga_mod_setup_topologies
