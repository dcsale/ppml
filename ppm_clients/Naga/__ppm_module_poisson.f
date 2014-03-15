      !-------------------------------------------------------------------------
      ! Subroutine : ppm_module_poisson
      !-------------------------------------------------------------------------
      ! Copyright (c) 2010 CSE Lab (ETH Zurich), MOSAIC Group (ETH Zurich),
      ! Center for Fluid Dynamics (DTU)
      !
      ! Notes on the init routine:
      ! Arguments:
      ! topoid
      ! meshid
      ! ppmpoisson
      ! fieldin
      ! fieldout
      ! greens function: custom array, custom function, integer to build-in
      ! info - return variable
      ! periodic/freespace flag
      ! real/fourier greens function (optional, only for custom array)
      ! derivatives (optional): none,curl,gradient,laplace, - spectral/FD
      ! possibly the option to de- and reallocate arrays of slabs/pencils
      !
      ! The fields returned from the subroutine have been ghosted/extrapolated
      ! as necessary.
      !
      ! Add finalize routine
      ! Allow custom kernels/greens functions
      ! Add switch to allocate and deallocate all arrays to save memory
      !
      ! This routine needs to be modified to allow cell centred data
      ! These routines should all be renamed - potentially not just Poisson
      ! equations can be solved.
      !
      ! IDEA for GREAT improvement: instead of
      ! -> map to slabs, fftXY, map to pencil, fftZ, Greens ... and back
      ! considerable memory, mapping and computational resources may be saved by
      ! -> map to pencils, fftX, map to slabs, fftYZ, Greens ... and back
      ! There by the mapping after fftX only needs to transfer 1/4th of the data
      ! the pencil domain requires 1/4th memory and fftX only needs treat the
      ! 1/4th of the data that is non-zero. This goes for the return as well.
      !
      ! ANOTHER IDEA: Make a flag for doing both vorticity reprojection on the
      ! input vorticity field and output the corresponding velocity.
      !-------------------------------------------------------------------------
      MODULE ppm_module_poisson
      !!! This module provides routines for solving the Poisson equation via
      !!! a Greens function. The Greens function is convolved with the RHS
      !!! spectrally using the FFTW library. The routines also do Helmholtz
      !!! reprojection to make fields solenoidal.
      !!!
      !!! Usage:
      !!! First a ppm_poisson_plan must initialised by calling ppm_poisson_init.
      !!! This is mainly done to initialise the FFTW library. An optional 'derive'
      !!! argument toggles various curl operations on the solution to the Poisson
      !!! equation.
      !!!
      !!! Then ppm_poisson_solve does the actual execution of the initialised FFTW
      !!! plans, convolutions, derivations, etc. The initialised Fourier
      !!! transformations and arrays may be used for some other operations;
      !!! presently just Helmholtz reprojection.
      !!!
      !!! So far no routine exists to remove any initialised ppm_poisson_plan
      USE ppm_module_fft
      USE ppm_module_substart
      USE ppm_module_substop
      USE ppm_module_write
      USE ppm_module_data,ONLY:ppm_rank,ppm_kind_single,ppm_kind_double,&
                          &ppm_t_topo,ppm_t_equi_mesh,&
                          &ppm_param_assign_internal,ppm_param_bcdef_periodic,&
                          &ppm_param_bcdef_freespace,&
                          &ppm_param_decomp_xy_slab,ppm_param_decomp_zpencil
      !-------------------------------------------------------------------------
      ! PPM Poisson type
      !-------------------------------------------------------------------------
      !!! Type containing the FFTW plan and its settings
      TYPE ppm_poisson_plan
         INTEGER :: derivatives
         INTEGER :: operation
         INTEGER :: bc
         INTEGER :: green
         REAL(ppm_kind_double) :: normkx
         REAL(ppm_kind_double) :: normky
         REAL(ppm_kind_double) :: normkz
         INTEGER, DIMENSION(3) :: nmfft
         !!!Size of the FFT
         INTEGER :: topoidxy
         !!!Topology id of xy topology
         INTEGER :: meshidxy
         !!!mesh id for xy mesh
         REAL(ppm_kind_double), DIMENSION(:), POINTER :: costxy=>NULL()
         !!!sub cost for xy topology
         INTEGER , DIMENSION(:,:),POINTER :: istartxy=>NULL()
         !!!sub index start for xy topology
         INTEGER, DIMENSION(3) :: nmxy
         !!!global number of grid points for xy topology
         INTEGER , DIMENSION(:,:),POINTER :: ndataxy=>NULL()
         !!!sub no. grid cells for xy topology
         INTEGER,DIMENSION(:),POINTER :: isublistxy=>NULL()
         !!!list of sub domains on current CPU on the xy topology
         INTEGER :: nsublistxy
         !!!number of sub domains for current CPU on the xy topology
         REAL(ppm_kind_double),DIMENSION(:,:,:,:,:),POINTER :: fldxyr=>NULL()
         !!!real slab field (xy topology)
         COMPLEX(ppm_kind_double),DIMENSION(:,:,:,:,:),POINTER:: fldxyc=>NULL()
         !!!complex slab field (xy topology)
         INTEGER :: topoidxyc
         !!!Topology id of complex xy topology
         INTEGER :: meshidxyc
         !!!mesh id for comple xy mesh
         REAL(ppm_kind_double), DIMENSION(:), POINTER :: costxyc=>NULL()
         !!!sub cost for xy topology
         INTEGER , DIMENSION(:,:),POINTER :: istartxyc=>NULL()
         !!!sub index start for complex xy topology
         INTEGER, DIMENSION(3) :: nmxyc
         !!!global number of grid points for complex xy topology
         INTEGER , DIMENSION(:,:),POINTER :: ndataxyc=>NULL()
         !!!sub no. grid cells for complex xy topology
         INTEGER,DIMENSION(:),POINTER :: isublistxyc=>NULL()
         !!!list of sub domains on current CPU on the complex xy topology
         INTEGER :: nsublistxyc
         !!!number of sub domains for current CPU on the complex xy topology
         TYPE(ppm_fft_plan) :: planfxy
         !!!fft plans for r2c complex xy slab
         TYPE(ppm_fft_plan) :: planbxy
         !!!fft plans for c2r complex xy slab
         TYPE(ppm_fft_plan) :: planfz
         !!!fft plans for forward c2c z pencil
         TYPE(ppm_fft_plan) :: planbz
         !!!fft plans for backward c2c z pencil
         INTEGER :: topoidz
         !!!Topology id of z topology
         INTEGER :: meshidz
         !!!mesh id for z mesh
         REAL(ppm_kind_double), DIMENSION(:), POINTER :: costz =>NULL()
         !!!sub cost for z topology
         INTEGER , DIMENSION(:,:),POINTER :: istartz=>NULL()
         !!!sub index start for z topology
         INTEGER, DIMENSION(3) :: nmz
         !!!global number of grid points for z topology
         INTEGER , DIMENSION(:,:),POINTER :: ndataz=>NULL()
         !!!sub no. grid cells for z topology
         INTEGER,DIMENSION(:),POINTER :: isublistz=>NULL()
         !!!list of sub domains on current CPU on the z topology
         INTEGER :: nsublistz
         !!!number of sub domains for current CPU on the z topology
         COMPLEX(ppm_kind_double),DIMENSION(:,:,:,:,:),POINTER:: fldzc1=>NULL()
         !!!complex pencil field 1 (z topology)
         COMPLEX(ppm_kind_double),DIMENSION(:,:,:,:,:),POINTER:: fldzc2=>NULL()
         !!!complex pencil field 2 (z topology)
         INTEGER :: topoidfr
         !!!Topology id of freespace topology
         INTEGER :: meshidfr
         !!!mesh id for freespace mesh
         REAL(ppm_kind_double), DIMENSION(:), POINTER :: costfr=>NULL()
         !!!sub cost
         INTEGER , DIMENSION(:,:),POINTER :: istartfr=>NULL()
         !!!sub index start
         INTEGER, DIMENSION(3) :: nmfr
         !!!global number of grid points
         INTEGER , DIMENSION(:,:),POINTER :: ndatafr=>NULL()
         !!!sub no. grid cells
         INTEGER,DIMENSION(:),POINTER :: isublistfr=>NULL()
         !!!list of sub domains on current CPU on the freespace topology
         INTEGER :: nsublistfr
         !!!number of sub domains for current CPU on the freespace topology
         !These field are only allocated as necessary:
         REAL(ppm_kind_double),DIMENSION(:,:,:,:),POINTER :: fldfrs=>NULL()
         !!!dummy array for the right hand side, for free space, scalar fields
         REAL(ppm_kind_double),DIMENSION(:,:,:,:,:),POINTER :: fldfrv=>NULL()
         !!!dummy array for the right hand side, for free space, vector fields
         REAL(ppm_kind_double),DIMENSION(:,:,:,:),POINTER :: fldgrnr=>NULL()
         !!!real Greens field, z-pencils, scalar
         COMPLEX(ppm_kind_double),DIMENSION(:,:,:,:),POINTER :: fldgrnc=>NULL()
         !!!complex Greens field, z-pencils, scalar
         REAL(ppm_kind_double),DIMENSION(:,:,:,:,:),POINTER :: drv_vr=>NULL()
         !!!dummy array for the right hand side, for free space, vector fields
      END TYPE ppm_poisson_plan
      INTEGER,PARAMETER :: ppm_poisson_opr_none = 1
      INTEGER,PARAMETER :: ppm_poisson_opr_vel = 2
      INTEGER,PARAMETER :: ppm_poisson_opr_repr = 3
      INTEGER,PARAMETER :: ppm_poisson_grn_pois_fre = 11
      INTEGER,PARAMETER :: ppm_poisson_grn_pois_per = 12
      INTEGER,PARAMETER :: ppm_poisson_grn_pois_blob2 = 13
      INTEGER,PARAMETER :: ppm_poisson_grn_pois_blob4 = 14
      INTEGER,PARAMETER :: ppm_poisson_grn_pois_blob6 = 15
      INTEGER,PARAMETER :: ppm_poisson_grn_pois_blob8 = 16
      INTEGER,PARAMETER :: ppm_poisson_grn_pois_blob10 = 17
      INTEGER,PARAMETER :: ppm_poisson_drv_sp = 21
      INTEGER,PARAMETER :: ppm_poisson_drv_fd2 = 22
      INTEGER,PARAMETER :: ppm_poisson_drv_fd4 = 23
      INTEGER,PARAMETER :: ppm_poisson_bc_per = 31
      INTEGER,PARAMETER :: ppm_poisson_bc_fre = 32
      !internal flags:
      INTEGER,PARAMETER :: ppm_poisson_curl = 91
      INTEGER,PARAMETER :: ppm_poisson_divergence = 92
      INTEGER,PARAMETER :: ppm_poisson_subtract = 93
      INTERFACE ppm_poisson_init
         MODULE PROCEDURE ppm_poisson_init
      END INTERFACE
      INTERFACE ppm_poisson_solve
         MODULE PROCEDURE ppm_poisson_solve
      END INTERFACE
      INTERFACE ppm_poisson_finalize
         MODULE PROCEDURE ppm_poisson_finalize
      END INTERFACE
      INTERFACE ppm_poisson_fd
         MODULE PROCEDURE ppm_poisson_fd
      END INTERFACE
      INTERFACE ppm_poisson_extrapolateghost
         MODULE PROCEDURE ppm_poisson_extrapolateghost_vr
      END INTERFACE
      CONTAINS
      !-------------------------------------------------------------------------
      ! Subroutine : ppm_poisson_init.f90
      !-------------------------------------------------------------------------
      ! Copyright (c) 2010 CSE Lab (ETH Zurich), MOSAIC Group (ETH Zurich),
      ! Center for Fluid Dynamics (DTU)
      !
      !-------------------------------------------------------------------------
      SUBROUTINE ppm_poisson_init(topoid,meshid,ppmpoisson,fieldin,fieldout,green,bc,&
                           & operation,derive,info,sigma)
      !!! Routine to initialise Greens function solution of the Poisson
      !!! equation. green is the flag defining which Greens function to use:
      !!! * ppm_poisson_grn_pois_per - Poisson equation, periodic boundaries
      !!! * ppm_poisson_grn_pois_fre - Poisson equation, freespace boundaries
      !!! * ppm_poisson_grn_reprojec - Do vorticity reprojection to kill divergence
      !!!
      !!! Eventually the routine should be overloaded to accept custom Greens
      !!! functions such that more general convolutions can be performed.
      !!! green should be expanded to include more buildin Greens functions.
      !!!
      !!! The routine should eventually accept an optional flag to toggle
      !!! deallocation of work arrays between calls to ppm_poisson_solve.
      USE ppm_module_mktopo
      USE ppm_module_topo_get
      USE ppm_module_mesh_define
      USE ppm_module_map_field
      USE ppm_module_map_field_global
      USE ppm_module_map
      IMPLICIT NONE
      !-------------------------------------------------------------------------
      ! Arguments
      !-------------------------------------------------------------------------
      INTEGER, INTENT(IN) :: topoid
      !!! Topology ID
      INTEGER, INTENT(IN) :: meshid
      !!! Mesh ID
      TYPE(ppm_poisson_plan),INTENT(INOUT) :: ppmpoisson
      !!! The PPM Poisson plan type (inspired by the FFTW plan)
      REAL(ppm_kind_double),DIMENSION(:,:,:,:,:),POINTER :: fieldin
      !!! Input data field. RHS to the Poisson equation/field to be convolved
      !@ strictly speaking fieldin is not being used in the init routine
      REAL(ppm_kind_double),DIMENSION(:,:,:,:,:),POINTER :: fieldout
      !!! Output data field
      INTEGER, INTENT(IN) :: green
      !!!flag to select build-in Greens functions:
      !!!ppm_poisson_grn_pois_per - Poisson equation, periodic boundaries
      !!!ppm_poisson_grn_pois_fre - Poisson equation, freespace boundaries
      !!!Mollified, high order kernels for freespace Poisson:
      !!!ppm_poisson_grn_pois_blob2 - 2nd order
      !!!ppm_poisson_grn_pois_blob4 - 4th order
      !!!ppm_poisson_grn_pois_blob6 - 6th order
      !!!ppm_poisson_grn_pois_blob8 - 8th order
      !!!ppm_poisson_grn_pois_blob10 - 10th order
      !!!
      !!!Eventually this should also accept custom Greens function
      INTEGER,INTENT(IN) :: bc
      !!!boundary condition for the convolution. Can be on of the following:
      !!!ppm_poisson_grn_bc_per, ppm_poisson_grn_bc_fre.
      !!!One could argue that this is redundant in the build-in operation
      INTEGER,INTENT(IN) :: operation
      !!!the functionality of the call:
      !!!ppm_poisson_opr_none, ppm_poisson_opr_vel, ppm_poisson_opr_repr
      INTEGER,INTENT(IN) :: derive
      !!!flag to toggle various derivatives of the solution (not to be used with
      !!!green=ppm_poisson_grn_reprojec):
      !!! * ppm_poisson_drv_sp
      !!! * ppm_poisson_drv_fd2
      !!! * ppm_poisson_drv_fd4
      INTEGER, INTENT(OUT) :: info
      INTEGER,INTENT(IN),OPTIONAL :: sigma
      !!!Mollification width relative to the mesh spacing
      !-------------------------------------------------------------------------
      ! Local variables
      !-------------------------------------------------------------------------
      REAL(ppm_kind_double) :: t0
      REAL(ppm_kind_double),DIMENSION(:,:),POINTER :: xp=>NULL() !particle positions
      TYPE(ppm_t_topo),POINTER :: topology=>NULL()
      TYPE(ppm_t_equi_mesh) :: mesh
      INTEGER ,DIMENSION(3) :: indl,indu
      INTEGER,PARAMETER :: MK = ppm_kind_double
      REAL(ppm_kind_double),PARAMETER :: PI=ACOS(-1.0_MK)
      REAL(ppm_kind_double) :: normfac
      !!!factor for the Greens function, including FFT normalization
      INTEGER :: i,j,k
      INTEGER :: kx,ky,kz
      INTEGER :: isubl,isub
      INTEGER,DIMENSION(3*2) :: bcdef
      INTEGER :: assigning
      INTEGER :: decomposition
      INTEGER,SAVE :: ttopoid
      INTEGER :: tmeshid
      REAL(ppm_kind_double) :: dx,dy,dz
      REAL(ppm_kind_double) :: Lx2,Ly2,Lz2
      REAL(ppm_kind_double) :: rho
      REAL(ppm_kind_double) :: sigmawidth
      REAL(ppm_kind_double) :: gzero
      REAL(ppm_kind_double),DIMENSION(3) :: tmpmin,tmpmax
      INTEGER, DIMENSION(3) :: maxndataxy,maxndataz
      INTEGER, DIMENSION(: ), POINTER :: dummynmxy,dummynmz
      !-------------------------------------------------------------------------
      ! Initialise routine
      !-------------------------------------------------------------------------
      CALL substart('ppm_poisson_init',t0,info)
      !-------------------------------------------------------------------------
      ! Investigate optional arguments, setup routine accordingly
      ! !@TODO: Also check if the input/output and derivatives match
      !-------------------------------------------------------------------------
      IF (operation .EQ. ppm_poisson_opr_none) THEN
        ppmpoisson%operation = ppm_poisson_opr_none
      ELSE IF (operation .EQ. ppm_poisson_opr_vel) THEN
        ppmpoisson%operation = ppm_poisson_opr_vel
      ELSE IF (operation .EQ. ppm_poisson_opr_repr) THEN
        ppmpoisson%operation = ppm_poisson_opr_repr
      !!ELSE IF (operation .EQ. ppm_poisson_opr_velrepr) THEN
      !! ppmpoisson%operation = ppm_poisson_opr_velrepr
      ELSE
        CALL ppm_write(ppm_rank,'ppm_poisson_init','Undefined operation input.',isub)
        info = -1
        GOTO 9999
      ENDIF
      !-------------------------------------------------------------------------
      ! Set the greens function
      !-------------------------------------------------------------------------
      ppmpoisson%green = green
      !-------------------------------------------------------------------------
      ! Nullify pointers from the ppmpoisson plans and the fftplans
      !-------------------------------------------------------------------------
      NULLIFY(xp)
      NULLIFY(ppmpoisson%costxy)
      NULLIFY(ppmpoisson%istartxy)
      NULLIFY(ppmpoisson%ndataxy)
      NULLIFY(ppmpoisson%istartxyc)
      NULLIFY(ppmpoisson%ndataxyc)
      NULLIFY(ppmpoisson%costz)
      NULLIFY(ppmpoisson%istartz)
      NULLIFY(ppmpoisson%ndataz)
      NULLIFY(ppmpoisson%planfxy%plan)
      NULLIFY(ppmpoisson%planbxy%plan)
      NULLIFY(ppmpoisson%planfz%plan)
      NULLIFY(ppmpoisson%planbz%plan)
      NULLIFY(dummynmxy)
      NULLIFY(dummynmz)
      !-------------------------------------------------------------------------
      ! Get topology and mesh values of input/output
      !-------------------------------------------------------------------------
      CALL ppm_topo_get(topoid,topology,info)
      IF (info .NE. 0) THEN
        CALL ppm_write(ppm_rank,'ppm_poisson_init','Failed to get topology.',isub)
        GOTO 9999
      ENDIF
      mesh = topology%mesh(meshid)
      !-------------------------------------------------------------------------
      ! Setup mesh sizes for intermediate meshes/topologies
      !-------------------------------------------------------------------------
      IF (bc .EQ. ppm_poisson_bc_per) THEN
        !size of real slabs
        ppmpoisson%nmxy (1) = mesh%nm(1)
        ppmpoisson%nmxy (2) = mesh%nm(2)
        ppmpoisson%nmxy (3) = mesh%nm(3)
        !size of complex slabs
        ppmpoisson%nmxyc(1) = (mesh%nm(1)-1)/2+1
        !!ppmpoisson%nmxyc(1) = mesh%nm(1)
        ppmpoisson%nmxyc(2) = mesh%nm(2)
        ppmpoisson%nmxyc(3) = mesh%nm(3)
        !size of complex pencils
        ppmpoisson%nmz (1) = (ppmpoisson%nmxyc(1))
        ppmpoisson%nmz (2) = (ppmpoisson%nmxyc(2))
        ppmpoisson%nmz (3) = (ppmpoisson%nmxyc(3))
        !size of the fft
        ppmpoisson%nmfft(1) = mesh%nm(1)-1
        ppmpoisson%nmfft(2) = mesh%nm(2)-1
        ppmpoisson%nmfft(3) = mesh%nm(3)-1
        !Inverse of the size of the domain squared
        Lx2 = 1.0_MK/(topology%max_physd(1)-topology%min_physd(1))**2
        Ly2 = 1.0_MK/(topology%max_physd(2)-topology%min_physd(2))**2
        Lz2 = 1.0_MK/(topology%max_physd(3)-topology%min_physd(3))**2
      ELSE IF (bc .EQ. ppm_poisson_bc_fre) THEN !vertex
        !size of real slabs
        ppmpoisson%nmxy (1) = mesh%nm(1)*2
        ppmpoisson%nmxy (2) = mesh%nm(2)*2
        ppmpoisson%nmxy (3) = mesh%nm(3)*2
        !size of complex slabs
        ppmpoisson%nmxyc(1) = (mesh%nm(1)*2)/2+1
        !!ppmpoisson%nmxyc(1) = mesh%nm(1)*2
        ppmpoisson%nmxyc(2) = mesh%nm(2)*2
        ppmpoisson%nmxyc(3) = mesh%nm(3)*2
        !size of complex pencils
        ppmpoisson%nmz (1) = (ppmpoisson%nmxyc(1))
        ppmpoisson%nmz (2) = (ppmpoisson%nmxyc(2))
        ppmpoisson%nmz (3) = (ppmpoisson%nmxyc(3))
        !size of the fft
        ppmpoisson%nmfft(1) = mesh%nm(1)*2
        ppmpoisson%nmfft(2) = mesh%nm(2)*2
        ppmpoisson%nmfft(3) = mesh%nm(3)*2
        !Determine the grid spacing !vertex
        dx = (topology%max_physd(1)-topology%min_physd(1))/REAL(mesh%nm(1)-1,MK) !vertex
        dy = (topology%max_physd(2)-topology%min_physd(2))/REAL(mesh%nm(2)-1,MK) !vertex
        dz = (topology%max_physd(3)-topology%min_physd(3))/REAL(mesh%nm(3)-1,MK) !vertex
      ENDIF
      !-------------------------------------------------------------------------
      ! Register derivation type.
      ! Create temporary derivation arrays if necessary.
      !-------------------------------------------------------------------------
      IF (( derive .EQ. ppm_poisson_drv_fd2 &
        & .OR. derive .EQ. ppm_poisson_drv_fd4)) THEN
        ppmpoisson%derivatives = derive
        IF (operation .EQ. ppm_poisson_opr_vel) THEN
          ALLOCATE(ppmpoisson%drv_vr(LBOUND(fieldout,1):UBOUND(fieldout,1),&
                                   & LBOUND(fieldout,2):UBOUND(fieldout,2),&
                                   & LBOUND(fieldout,3):UBOUND(fieldout,3),&
                                   & LBOUND(fieldout,4):UBOUND(fieldout,4),&
                                   & LBOUND(fieldout,5):UBOUND(fieldout,5)))
        ENDIF
      ELSE IF (derive .EQ. ppm_poisson_drv_sp) THEN
        ppmpoisson%derivatives = ppm_poisson_drv_sp
      ELSE
        CALL ppm_write(ppm_rank,'ppm_poisson_init','Undefined derivation input.',isub)
        info = -1
        GOTO 9999
      ENDIF
      !-------------------------------------------------------------------------
      ! Create spectral scaling components always. Just in case some
      ! reprojection comes up.
      ! The conditionals need to be for not just the Poisson equation
      !-------------------------------------------------------------------------
      IF (bc .EQ. ppm_poisson_bc_per) THEN
        ppmpoisson%bc = ppm_poisson_bc_per
        ppmpoisson%normkx = &
        & 2.0_MK*PI/(topology%max_physd(1)-topology%min_physd(1))
        ppmpoisson%normky = &
        & 2.0_MK*PI/(topology%max_physd(2)-topology%min_physd(2))
        ppmpoisson%normkz = &
        & 2.0_MK*PI/(topology%max_physd(3)-topology%min_physd(3))
      ELSE IF (bc .EQ. ppm_poisson_bc_fre) THEN
        ppmpoisson%bc = ppm_poisson_bc_fre
        ppmpoisson%normkx = &
        !!& 2.0_MK*PI/((topology%max_physd(1)-topology%min_physd(1))*2.0_MK)
        & 2.0_MK*PI/(dx*REAL(2*mesh%nm(1),MK))
        ppmpoisson%normky = &
        !!& 2.0_MK*PI/((topology%max_physd(2)-topology%min_physd(2))*2.0_MK)
        & 2.0_MK*PI/(dy*REAL(2*mesh%nm(2),MK))
        ppmpoisson%normkz = &
        !!& 2.0_MK*PI/((topology%max_physd(3)-topology%min_physd(3))*2.0_MK)
        & 2.0_MK*PI/(dz*REAL(2*mesh%nm(3),MK))
      ELSE
        CALL ppm_write(ppm_rank,'ppm_poisson_init','Undefined boundary condition input.',isub)
        info = -1
        GOTO 9999
      ENDIF
      !-------------------------------------------------------------------------
      ! Create new slab topology
      !-------------------------------------------------------------------------
      ttopoid = 0
      tmeshid = -1
      decomposition = ppm_param_decomp_xy_slab
      assigning = ppm_param_assign_internal
      IF (bc .EQ. ppm_poisson_bc_per) THEN
        bcdef = ppm_param_bcdef_periodic
      ELSE IF (bc .EQ. ppm_poisson_bc_fre) THEN
        bcdef = ppm_param_bcdef_freespace
      ENDIF
      tmpmin = topology%min_physd
      tmpmax = topology%max_physd
      CALL ppm_mktopo(ttopoid,tmeshid,xp,0,&
      & decomposition,assigning,&
      & tmpmin,tmpmax,bcdef,&
      & (/0,0,0/),ppmpoisson%costxy,&
      & ppmpoisson%nmxy,info)
      IF (info .NE. 0) THEN
        CALL ppm_write(ppm_rank,'ppm_poisson_init','Failed to create xy-topology.',isub)
        GOTO 9999
      ENDIF
      ppmpoisson%topoidxy = ttopoid
      ppmpoisson%meshidxy = tmeshid
      !-------------------------------------------------------------------------
      ! Get additional xy-mesh information
      !-------------------------------------------------------------------------
      CALL ppm_topo_get_meshinfo(ppmpoisson%topoidxy,ppmpoisson%meshidxy, &
         & dummynmxy,ppmpoisson%istartxy,ppmpoisson%ndataxy,maxndataxy, &
         & ppmpoisson%isublistxy,ppmpoisson%nsublistxy,info)
      !-------------------------------------------------------------------------
      ! Create complex slab mesh
      !-------------------------------------------------------------------------
      ttopoid = ppmpoisson%topoidxy
      tmeshid = -1
      CALL ppm_mesh_define(ttopoid,tmeshid,&
      & ppmpoisson%nmxyc,ppmpoisson%istartxyc,ppmpoisson%ndataxyc,info)
      IF (info .NE. 0) THEN
        CALL ppm_write(ppm_rank,'ppm_poisson_init','Failed to create complex xy mesh definition.',isub)
        GOTO 9999
      ENDIF
      ppmpoisson%meshidxyc = tmeshid
      !-------------------------------------------------------------------------
      ! Create new pencil topology
      !-------------------------------------------------------------------------
      ttopoid = 0
      tmeshid = -1
      IF (bc .EQ. ppm_poisson_bc_per) THEN
        bcdef = ppm_param_bcdef_periodic
      ELSE IF (bc .EQ. ppm_poisson_bc_fre) THEN
        bcdef = ppm_param_bcdef_freespace
      ENDIF
      assigning = ppm_param_assign_internal
      decomposition = ppm_param_decomp_zpencil
      CALL ppm_mktopo(ttopoid,tmeshid,xp,0,&
      & decomposition,assigning,&
      & tmpmin,tmpmax,bcdef,&
      & (/0,0,0/),ppmpoisson%costz,&
      & ppmpoisson%nmz,info)
      IF (info .NE. 0) THEN
        CALL ppm_write(ppm_rank,'ppm_poisson_init','Failed to create z-topology.',isub)
        GOTO 9999
      ENDIF
      ppmpoisson%topoidz = ttopoid
      ppmpoisson%meshidz = tmeshid
      !-------------------------------------------------------------------------
      ! Get additional z-mesh information
      !-------------------------------------------------------------------------
      CALL ppm_topo_get_meshinfo(ppmpoisson%topoidz,ppmpoisson%meshidz, &
         & dummynmz,ppmpoisson%istartz,ppmpoisson%ndataz,maxndataz, &
         & ppmpoisson%isublistz,ppmpoisson%nsublistz,info)
      !-------------------------------------------------------------------------
      ! Set and get minimum and maximum indicies
      !-------------------------------------------------------------------------
      indl(1) = 1
      indl(2) = 1
      indl(3) = 1
      !-------------------------------------------------------------------------
      ! Allocate real xy slabs
      !-------------------------------------------------------------------------
      ALLOCATE(ppmpoisson%fldxyr(3,&
      & indl(1):maxndataxy(1),indl(2):maxndataxy(2),indl(3):maxndataxy(3),&
      & 1:ppmpoisson%nsublistxy),stat=info)
      !-------------------------------------------------------------------------
      ! Set and get minimum and maximum indicies of COMPLEX xy slabs
      !-------------------------------------------------------------------------
      indl(1) = 1
      indl(2) = 1
      indl(3) = 1
      indu(1) = 0
      indu(2) = 0
      indu(3) = 0
      DO isub=1,ppmpoisson%nsublistxy
        isubl = ppmpoisson%isublistxy(isub)
        indu(1) = MAX(indu(1),ppmpoisson%ndataxyc(1,isubl))
        indu(2) = MAX(indu(2),ppmpoisson%ndataxyc(2,isubl))
        indu(3) = MAX(indu(3),ppmpoisson%ndataxyc(3,isubl))
      ENDDO
      !-------------------------------------------------------------------------
      ! Allocate complex xy slabs
      !-------------------------------------------------------------------------
      ALLOCATE(ppmpoisson%fldxyc(3,&
      & indl(1):indu(1),indl(2):indu(2),indl(3):indu(3),&
      & 1:ppmpoisson%nsublistxy),stat=info)
      !-------------------------------------------------------------------------
      ! Allocate two complex z pencils + Greens fcn array !@check return vars.
      !-------------------------------------------------------------------------
      ALLOCATE(ppmpoisson%fldzc1(3,&
      & indl(1):maxndataz(1),indl(2):maxndataz(2),indl(3):maxndataz(3),&
      & 1:ppmpoisson%nsublistz),stat=info)
      ALLOCATE(ppmpoisson%fldzc2(3,&
      & indl(1):maxndataz(1),indl(2):maxndataz(2),indl(3):maxndataz(3),&
      & 1:ppmpoisson%nsublistz),stat=info)
      !-------------------------------------------------------------------------
      ! The complex Greens function is always kept on the z-pencil topology
      !-------------------------------------------------------------------------
      IF (bc .EQ. ppm_poisson_bc_per) THEN
        ALLOCATE(ppmpoisson%fldgrnr(&
        & indl(1):maxndataz(1),indl(2):maxndataz(2),indl(3):maxndataz(3),&
        & 1:ppmpoisson%nsublistz),stat=info)
      ELSE IF (bc .EQ. ppm_poisson_bc_fre) THEN
        ALLOCATE(ppmpoisson%fldgrnc(&
        & indl(1):maxndataz(1),indl(2):maxndataz(2),indl(3):maxndataz(3),&
        & 1:ppmpoisson%nsublistz),stat=info)
      ENDIF
      !-------------------------------------------------------------------------
      ! Set up xy FFT plans
      ! The inverse plan takes the returning topology since it has the full size
      !-------------------------------------------------------------------------
      CALL ppm_fft_forward_2d(ppmpoisson%topoidxy,ppmpoisson%meshidxy,&
      & ppmpoisson%planfxy,ppmpoisson%fldxyr,&
      & ppmpoisson%fldxyc,info)
      CALL ppm_fft_backward_2d(ppmpoisson%topoidxy,ppmpoisson%meshidxy,&
      & ppmpoisson%planbxy,ppmpoisson%fldxyc,&
      & ppmpoisson%fldxyr,info)
      !-------------------------------------------------------------------------
      ! Set up z FFT plans
      !-------------------------------------------------------------------------
      CALL ppm_fft_forward_1d(ppmpoisson%topoidz,ppmpoisson%meshidz,&
      & ppmpoisson%planfz,ppmpoisson%fldzc1,&
      & ppmpoisson%fldzc2,info)
      CALL ppm_fft_backward_1d(ppmpoisson%topoidz,ppmpoisson%meshidz,&
      & ppmpoisson%planbz,ppmpoisson%fldzc2,&
      & ppmpoisson%fldzc1,info)
      !-------------------------------------------------------------------------
      ! Compute Greens function. Analytic, periodic
      !
      ! (d2_/dx2 + d2_/dy2 + d2_/dz2)psi = -omega =>
      ! -4*pi2(kx2 + ky2 + kz2)PSI = -OMEGA =>
      ! PSI = 1/(4*pi2)*1/(kx2 + ky2 + kz2)OMEGA
      ! Consider allocating the periodic Greens function in a complex array
      ! (for simplicity)
      !-------------------------------------------------------------------------
      IF (green .EQ. ppm_poisson_grn_pois_per) THEN
        ! Scaling the spectral coefficients...
        ! one minus due to (i*k)^2 and another due to the Poisson equation
        normfac = 1.0_MK/(4.0_MK*PI*PI * &
                !and normalisation of FFTs (full domain) !vertex
                & REAL(ppmpoisson%nmfft(1),MK)* &
                & REAL(ppmpoisson%nmfft(2),MK)* &
                & REAL(ppmpoisson%nmfft(3),MK))
        DO isub=1,ppmpoisson%nsublistz
          isubl=ppmpoisson%isublistz(isub)
          DO k=1,ppmpoisson%ndataz(3,isubl)
            DO j=1,ppmpoisson%ndataz(2,isubl)
              DO i=1,ppmpoisson%ndataz(1,isubl)
                kx = i-1 + (ppmpoisson%istartz(1,isubl)-1)
                ky = j-1 + (ppmpoisson%istartz(2,isubl)-1)
                kz = k-1 + (ppmpoisson%istartz(3,isubl)-1)
                !This is a nasty way to do this but it is only done once so...:
                IF (kx .GT. (ppmpoisson%nmfft(1)/2)) kx = kx-(ppmpoisson%nmfft(1))
                IF (ky .GT. (ppmpoisson%nmfft(2)/2)) ky = ky-(ppmpoisson%nmfft(2))
                IF (kz .GT. (ppmpoisson%nmfft(3)/2)) kz = kz-(ppmpoisson%nmfft(3))
                ppmpoisson%fldgrnr(i,j,k,isub) = &
                      & normfac/(REAL(kx*kx,ppm_kind_double)*Lx2 &
                      & + REAL(ky*ky,ppm_kind_double)*Ly2 &
                      & + REAL(kz*kz,ppm_kind_double)*Lz2)
                !Take care of singularity
                !This is nasty as well
                IF ((kx*kx+ky*ky+kz*kz) .EQ. 0) THEN
                  ppmpoisson%fldgrnr(i,j,k,isub) = 0.0_MK
                ENDIF
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      !-------------------------------------------------------------------------
      ! Compute real free space Greens function. Analytic
      ! The Greens function is initialised temporarily in the real xy slabs,
      ! then FFTed in those directions, mapped to z-pencils, FFTed in z and
      ! finally copied to ppmpoisson%fldgrnc. The real xy slabs have already
      ! been setup for FFTs etc so they offer a convenient container for the
      ! FFTing the Greens function instead of setting up the whole apparatus
      ! for this one-time affair.
      ! These loops must run over the padded(extended) domain thus %ndataxy
      ! \nabla \Psi = -\omega
      ! The Greens function is 1/4piR (R is distance to origo) and includes the
      ! minus of the RHS
      !-------------------------------------------------------------------------
      ELSE !IF (green .EQ. ppm_poisson_grn_pois_fre) THEN
        !-----------------------------------------------------------------------
        ! First initialise the real Greens function
        !@alternatively this could come from as input
        !-----------------------------------------------------------------------
        !there should NOT be a minus here since THIS Greens function takes
        !the minus of the Poisson equation into account
        normfac = 1.0_MK/(4.0_MK*PI* &
        !remembering FFT normalization of ALL points: !vertex
        !this includes 3 FFTs: 1) FFT of Green 2) FFT of omega 3) IFFT of phi
                & REAL(ppmpoisson%nmfft(1),MK)* &
                & REAL(ppmpoisson%nmfft(2),MK)* &
                & REAL(ppmpoisson%nmfft(3),MK))*dx*dy*dz
        !-------------------------------------------------------------------------
        ! Determine mollification width
        !-------------------------------------------------------------------------
        IF (PRESENT(sigma)) THEN
          sigmawidth = dx*REAL(sigma,MK)
        ELSE
          sigmawidth = dx
        ENDIF
        DO isub=1,ppmpoisson%nsublistxy
          isubl=ppmpoisson%isublistxy(isub)
          DO k=1,ppmpoisson%ndataxy(3,isubl)
            DO j=1,ppmpoisson%ndataxy(2,isubl)
              DO i=1,ppmpoisson%ndataxy(1,isubl)
                !kx,ky,kz implies that they are spectral coordinates. they are not
                kx = i-1 + (ppmpoisson%istartxy(1,isubl)-1)
                ky = j-1 + (ppmpoisson%istartxy(2,isubl)-1)
                kz = k-1 + (ppmpoisson%istartxy(3,isubl)-1)
                !This is a nasty way to do this but it is only done once so...:
                IF (kx .GT. (ppmpoisson%nmfft(1))/2) kx = kx-(ppmpoisson%nmfft(1))
                IF (ky .GT. (ppmpoisson%nmfft(2))/2) ky = ky-(ppmpoisson%nmfft(2))
                IF (kz .GT. (ppmpoisson%nmfft(3))/2) kz = kz-(ppmpoisson%nmfft(3))
                !Hoping for compiler optimization here. Cross your fingers:
                !----FREESPACE GREENS FUNCTION----
                IF (green .EQ. ppm_poisson_grn_pois_fre) THEN
                  ppmpoisson%fldxyr(1,i,j,k,isub) = &
                          & normfac/(SQRT( REAL(kx*kx,ppm_kind_double)*dx*dx+ &
                          & REAL(ky*ky,ppm_kind_double)*dy*dy+ &
                          & REAL(kz*kz,ppm_kind_double)*dz*dz))
                  ppmpoisson%fldxyr(2,i,j,k,isub) = 0.0_MK
                  ppmpoisson%fldxyr(3,i,j,k,isub) = 0.0_MK
                  !Take care of singularity (This is nasty as well)
                  gzero = 1.0_MK*normfac*4.0_MK*PI
                !----BLOB----2nd ORDER---:
                ELSE IF (green .EQ. ppm_poisson_grn_pois_blob2) THEN
                  rho = SQRT( REAL(kx*kx,ppm_kind_double)*dx*dx+ &
                            & REAL(ky*ky,ppm_kind_double)*dy*dy+ &
                            & REAL(kz*kz,ppm_kind_double)*dz*dz)/sigmawidth
                  ppmpoisson%fldxyr(1,i,j,k,isub) = &
                          & normfac/(rho)*ERF(rho/SQRT(2.0_MK))/sigmawidth
                  ppmpoisson%fldxyr(2,i,j,k,isub) = 0.0_MK
                  ppmpoisson%fldxyr(3,i,j,k,isub) = 0.0_MK
                  !Take care of singularity (This is nasty as well)
                  gzero = 4.0_MK*PI*normfac* &
                        & 1.0_MK/4.0_MK*SQRT(2.0_MK) &
                        & /(PI**(3.0_MK/2.0_MK)*sigmawidth)
                !----BLOB----4th ORDER---cell-normalised:
                ELSE IF (green .EQ. ppm_poisson_grn_pois_blob4) THEN
                  rho = SQRT( REAL(kx*kx,ppm_kind_double)*dx*dx+ &
                            & REAL(ky*ky,ppm_kind_double)*dy*dy+ &
                            & REAL(kz*kz,ppm_kind_double)*dz*dz)/sigmawidth
                  ppmpoisson%fldxyr(1,i,j,k,isub) = &
                          & normfac/(rho*sigmawidth)* &
                          & ( &
                          & ERF(rho/SQRT(2.0_MK)) + &
                          & 1.0_MK/SQRT(PI*2.0_MK) * rho &
                          & * EXP(-rho**2/2.0_MK) &
                          & )
                  ppmpoisson%fldxyr(2,i,j,k,isub) = 0.0_MK
                  ppmpoisson%fldxyr(3,i,j,k,isub) = 0.0_MK
                  !Take care of singularity (This is nasty as well)
                  gzero = 4.0_MK*PI*normfac* &
                        & 3.0_MK/8.0_MK*SQRT(2.0_MK) &
                        & /(PI**(3.0_MK/2.0_MK)*sigmawidth)
                !----BLOB----6th ORDER---cell-normalised:
                ELSE IF (green .EQ. ppm_poisson_grn_pois_blob6) THEN
                  rho = SQRT( REAL(kx*kx,ppm_kind_double)*dx*dx+ &
                            & REAL(ky*ky,ppm_kind_double)*dy*dy+ &
                            & REAL(kz*kz,ppm_kind_double)*dz*dz)/sigmawidth
                  ppmpoisson%fldxyr(1,i,j,k,isub) = &
                          & normfac/(rho*sigmawidth)* &
                          & ( &
                          & ERF(rho/SQRT(2.0_MK)) + &
                          & 1.0_MK/SQRT(PI*2.0_MK) * &
                          & (7.0_MK/4.0_MK *rho - 1.0_MK/4.0_MK *rho**3) &
                          & * EXP(-rho**2/2.0_MK) &
                          & )
                  ppmpoisson%fldxyr(2,i,j,k,isub) = 0.0_MK
                  ppmpoisson%fldxyr(3,i,j,k,isub) = 0.0_MK
                  !Take care of singularity (This is nasty as well)
                  gzero = 4.0_MK*PI*normfac* &
                        & 15.0_MK/32.0_MK*SQRT(2.0_MK) &
                        & /(PI**(3.0_MK/2.0_MK)*sigmawidth)
                !----BLOB----8th ORDER---cell-normalised:
                ELSE IF (green .EQ. ppm_poisson_grn_pois_blob8) THEN
                  rho = SQRT( REAL(kx*kx,ppm_kind_double)*dx*dx+ &
                            & REAL(ky*ky,ppm_kind_double)*dy*dy+ &
                            & REAL(kz*kz,ppm_kind_double)*dz*dz)/sigmawidth
                  ppmpoisson%fldxyr(1,i,j,k,isub) = &
                          & normfac/(rho*sigmawidth)* &
                          & ( &
                          & ERF(rho/SQRT(2.0_MK)) + &
                          & 1.0_MK/SQRT(PI*2.0_MK) * &
                          & (19.0_MK/8.0_MK *rho - 2.0_MK/3.0_MK *rho**3 &
                          & + 1.0_MK/24.0_MK *rho**5) &
                          & * EXP(-rho**2/2.0_MK) &
                          & )
                  ppmpoisson%fldxyr(2,i,j,k,isub) = 0.0_MK
                  ppmpoisson%fldxyr(3,i,j,k,isub) = 0.0_MK
                  !Take care of singularity (This is nasty as well)
                  gzero = 4.0_MK*PI*normfac* &
                        & 35.0_MK/64.0_MK*SQRT(2.0_MK) &
                        & /(PI**(3.0_MK/2.0_MK)*sigmawidth)
                !----BLOB----10th ORDER---cell-normalised:
                ELSE IF (green .EQ. ppm_poisson_grn_pois_blob10) THEN
                  rho = SQRT( REAL(kx*kx,ppm_kind_double)*dx*dx+ &
                            & REAL(ky*ky,ppm_kind_double)*dy*dy+ &
                            & REAL(kz*kz,ppm_kind_double)*dz*dz)/sigmawidth
                  ppmpoisson%fldxyr(1,i,j,k,isub) = &
                          & normfac/(rho*sigmawidth)* &
                          & ( &
                          & ERF(rho/SQRT(2.0_MK)) + &
                          & 1.0_MK/SQRT(PI*2.0_MK) * &
                          & (187.0_MK/64.0_MK *rho - &
                          & 233.0_MK/192_MK *rho**3 + &
                          & 29_MK/192.0_MK *rho**5 - &
                          & 1.0_MK/192_MK *rho**7) &
                          & * EXP(-rho**2/2.0_MK) &
                          & )
                  ppmpoisson%fldxyr(2,i,j,k,isub) = 0.0_MK
                  ppmpoisson%fldxyr(3,i,j,k,isub) = 0.0_MK
                  !Take care of singularity (This is nasty as well)
                  gzero = 4.0_MK*PI*normfac* &
                        & 315.0_MK/512.0_MK*SQRT(2.0_MK) &
                        & /(PI**(3.0_MK/2.0_MK)*sigmawidth)
                ENDIF
                !And as promised: Take care of singularity (This is nasty as well)
                IF ((kx*kx+ky*ky+kz*kz) .EQ. 0) THEN
                  !Professor Nutbutter is out.
                  ppmpoisson%fldxyr(:,i,j,k,isub) = gzero
                ENDIF
              ENDDO
            ENDDO
          ENDDO
        ENDDO
        !-------------------------------------------------------------------------
        ! FOURIER TRANSFORM AND MAP GALORE
        ! This part should be used both for freespace and a custom Greens function
        !-------------------------------------------------------------------------
        !-----------------------------------------------------------------------
        ! Do slab FFT (XY)
        !-----------------------------------------------------------------------
        CALL ppm_fft_execute_2d(ppmpoisson%topoidxy,&
        & ppmpoisson%meshidxy, ppmpoisson%planfxy, &
        & ppmpoisson%fldxyr, ppmpoisson%fldxyc, &
        & info)
        !-----------------------------------------------------------------------
        ! Map to the pencils (Z)
        !-----------------------------------------------------------------------
        !Initialise
        CALL ppm_map_field_global(&
        & ppmpoisson%topoidxy, &
        & ppmpoisson%topoidz, &
        & ppmpoisson%meshidxyc, &
        & ppmpoisson%meshidz,info)
        IF (info .NE. 0) THEN
          CALL ppm_write(ppm_rank, 'ppm_poisson_solve','Failed to initialise field mapping.',isub)
          GOTO 9999
        ENDIF
        !Push the data
        CALL ppm_map_field_push(&
        & ppmpoisson%topoidxy, &
        & ppmpoisson%meshidxyc,ppmpoisson%fldxyc,3,info)
        IF (info .NE. 0) THEN
          CALL ppm_write(ppm_rank, 'ppm_poisson_solve','Failed to push vector field.',isub)
          GOTO 9999
        ENDIF
        !Send
        CALL ppm_map_field_send(info)
        IF (info .NE. 0) THEN
          CALL ppm_write(ppm_rank, 'ppm_poisson_solve','Failed to send field.',isub)
          GOTO 9999
        ENDIF
        !Retrieve
        CALL ppm_map_field_pop(&
        & ppmpoisson%topoidz, &
        & ppmpoisson%meshidz,ppmpoisson%fldzc1, &
        & 3,(/0,0,0/),info)
        IF (info .NE. 0) THEN
          CALL ppm_write(ppm_rank, 'ppm_poisson_solve','Failed to pop vector field.',isub)
          GOTO 9999
        ENDIF
        !-----------------------------------------------------------------------
        ! Do pencil FFT (Z)
        !-----------------------------------------------------------------------
        CALL ppm_fft_execute_1d(ppmpoisson%topoidz,&
        & ppmpoisson%meshidz, ppmpoisson%planfz, &
        & ppmpoisson%fldzc1, ppmpoisson%fldzc2, &
        & info)
        !-----------------------------------------------------------------------
        ! Copy first component of the Fourier transformed vector to fldgrnc
        !-----------------------------------------------------------------------
        DO isub=1,ppmpoisson%nsublistz
          isubl=ppmpoisson%isublistz(isub)
          DO k=1,ppmpoisson%ndataz(3,isubl)
            DO j=1,ppmpoisson%ndataz(2,isubl)
              DO i=1,ppmpoisson%ndataz(1,isubl)
                ppmpoisson%fldgrnc(i,j,k,isub) = ppmpoisson%fldzc2(1,i,j,k,isub)
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      END IF
      !-------------------------------------------------------------------------
      ! Or alternatively FFT real Green from input
      !-------------------------------------------------------------------------
      !-------------------------------------------------------------------------
      ! Deallocate fields? !@
      !-------------------------------------------------------------------------
      !-------------------------------------------------------------------------
      ! Return
      !-------------------------------------------------------------------------
 9999 CONTINUE
      CALL substop('ppm_poisson_init',t0,info)
      RETURN
      END SUBROUTINE ppm_poisson_init
      !-------------------------------------------------------------------------
      ! Subroutine : ppm_poisson_solve.f90
      !-------------------------------------------------------------------------
      ! Copyright (c) 2010 CSE Lab (ETH Zurich), MOSAIC Group (ETH Zurich),
      ! Center for Fluid Dynamics (DTU)
      !
      !-------------------------------------------------------------------------
      SUBROUTINE ppm_poisson_solve(topoid,meshid,ppmpoisson,fieldin,fieldout,gstw,info,&
                         & tmpoperation,tmpderivative)
      !!! Routine to perform the Greens function solution of the Poisson
      !!! equation. All settings are defined in ppm_poisson_initdef and stored
      !!! in the ppmpoisson plan.
      !!!
      !!! The tmpoperation argument allows the use of a
      !!! different Greens function or operation than initialised. This is
      !!! particularly useful for Helmholtz reprojection
      !!! (ppm_poisson_grn_reprojec).
      !!!
      !!! [NOTE]
      !!! fieldin and fieldout must NOT be the same fields. In-place FFTs have
      !!! not been implemented.
      USE ppm_module_map_field
      USE ppm_module_map_field_global
      USE ppm_module_map
      !USE ppm_module_topo_get !@PERFREE
      IMPLICIT NONE
      include 'mpif.h'
      !-------------------------------------------------------------------------
      ! Arguments
      !-------------------------------------------------------------------------
      INTEGER, INTENT(IN) :: topoid
      !!! Topology ID
      INTEGER, INTENT(IN) :: meshid
      !!! Mesh ID
      TYPE(ppm_poisson_plan),INTENT(INOUT) :: ppmpoisson
      !!! The PPM Poisson plan
      REAL(ppm_kind_double),DIMENSION(:,:,:,:,:),POINTER :: fieldin
      !!! Input data field
      REAL(ppm_kind_double),DIMENSION(:,:,:,:,:),POINTER :: fieldout
      !!! Output data field
      INTEGER,DIMENSION(3),INTENT(IN) :: gstw
      !!! Ghost layer width
      INTEGER, INTENT(OUT) :: info
      !!! Return status, 0 upon succes
      INTEGER,OPTIONAL,INTENT(IN) :: tmpoperation
      !!! Temporary operation (useful for reprojection)
      INTEGER,OPTIONAL,INTENT(IN) :: tmpderivative
      !!! Temporary operation (useful for reprojection)
      !-------------------------------------------------------------------------
      ! Local variables
      !-------------------------------------------------------------------------
      INTEGER,PARAMETER :: MK = ppm_kind_double
      REAL(ppm_kind_double) :: t0
      INTEGER :: isub,isubl
      INTEGER :: i,j,k
      INTEGER :: info2
      INTEGER :: presentoperation
      INTEGER :: presentderivatives
      COMPLEX(ppm_kind_double) :: divomega
      INTEGER :: gi,gj,gk
      COMPLEX(ppm_kind_double) :: kx,ky,kz
      COMPLEX(ppm_kind_double) :: phix,phiy,phiz
      REAL(ppm_kind_double) :: normfac
      !-------------------------------------------------------------------------
      ! Initialise routine
      !-------------------------------------------------------------------------
      CALL substart('ppm_poisson_solve',t0,info)
      !-------------------------------------------------------------------------
      ! Check if we run a different/temporary operation
      !-------------------------------------------------------------------------
      IF (PRESENT(tmpoperation)) THEN
         presentoperation = tmpoperation
      ELSE
         presentoperation = ppmpoisson%operation
      ENDIF
      IF (PRESENT(tmpderivative)) THEN
         presentderivatives = tmpderivative
      ELSE
         presentderivatives = ppmpoisson%derivatives
      ENDIF
      !Vorticity reprojection is not fully implemented using finite differences
      IF ((presentderivatives .NE. ppm_poisson_drv_sp) .AND. &
        & (presentoperation .EQ. ppm_poisson_opr_repr)) THEN
         CALL ppm_write(ppm_rank,'ppm_poisson_solve','Warning: Vorticity reprojection is performed spectrally.',info2)
         presentderivatives = ppm_poisson_drv_sp
      ENDIF
      !-------------------------------------------------------------------------
      ! Perhaps check if ghostlayer suffices for a given fd stencil
      !-------------------------------------------------------------------------
      !-------------------------------------------------------------------------
      ! Perhaps allocate (and deallocate) arrays if a flag has been set
      !-------------------------------------------------------------------------
      !-----------------------------------------------------------------------
      ! Set the real xy slabs 0 (part of the 0 padding) for free-space
      !@ free-space calculations and reprojection may cause problems !why?
      ! this should be valid for all freespace as repr_fd is only for freespace
      !-----------------------------------------------------------------------
      !!IF ((presentoperation .EQ. ppm_poisson_grn_repr .AND. &
        !!& ppmpoisson%operation .EQ. ppm_poisson_grn_pois_fre) .OR. &
        !!& presentoperation .EQ. ppm_poisson_grn_repr_fd) THEN
        !!ppmpoisson%fldxyr = 0.0_MK
      !!ENDIF
      IF (ppmpoisson%bc .EQ. ppm_poisson_grn_pois_fre) THEN
        ppmpoisson%fldxyr = 0.0_MK
      ENDIF
      !!IF (ppmpoisson%bc .EQ. ppm_poisson_grn_repr) THEN
        !!ppmpoisson%fldxyr = 0.0_MK
      !!ENDIF
      !-----------------------------------------------------------------------
      ! For discontinuous vorticity at boundaries, calculate vort. divergence
      ! I.e. do vorticity reprojection in a freespace domain using FD
      !-----------------------------------------------------------------------
      IF (presentoperation .EQ. ppm_poisson_opr_repr .AND. &
        & (presentderivatives .EQ. ppm_poisson_drv_fd2 .OR. &
        & presentderivatives .EQ. ppm_poisson_drv_fd4)) THEN! .AND. &
        !& ppmpoisson%bc .EQ. ppm_poisson_grn_pois_fre) THEN
        fieldout = 0.0_MK
        CALL ppm_poisson_fd(topoid,meshid,fieldin,fieldout,&
                           & ppm_poisson_divergence, presentderivatives ,info)
      ENDIF
      !-----------------------------------------------------------------------
      ! Map data globally to the slabs (XY)
      ! This is where the vorticity is extended and padded with 0 for free-space
      !-----------------------------------------------------------------------
      !Initialise
      CALL ppm_map_field_global(&
      & topoid, &
      & ppmpoisson%topoidxy, &
      & meshid, &
      & ppmpoisson%meshidxy,info)
      IF (info .NE. 0) THEN
        CALL ppm_write(ppm_rank,'ppm_poisson_solve','Failed to initialise field mapping.',info2)
        GOTO 9999
      ENDIF
      !Push the data
      ! (in case of freespace reprojection push fieldout)
      IF (presentoperation .EQ. ppm_poisson_opr_repr .AND. &
        & (presentderivatives .EQ. ppm_poisson_drv_fd2 .OR. &
        & presentderivatives .EQ. ppm_poisson_drv_fd4)) THEN! .AND. &
        !& ppmpoisson%bc .EQ. ppm_poisson_grn_pois_fre) THEN
      !!IF (presentoperation .EQ. ppm_poisson_grn_repr_fd) THEN
        CALL ppm_map_field_push(&
        & topoid, &
        & meshid,fieldout,3,info)
        IF (info .NE. 0) THEN
          CALL ppm_write(ppm_rank, 'ppm_poisson_solve','Failed to push vector field.',info2)
          GOTO 9999
        ENDIF
      ELSE
        CALL ppm_map_field_push(&
        & topoid, &
        & meshid,fieldin,3,info)
        IF (info .NE. 0) THEN
          CALL ppm_write(ppm_rank, 'ppm_poisson_solve','Failed to push vector field.',info2)
          GOTO 9999
        ENDIF
      ENDIF
      !Send
      CALL ppm_map_field_send(info)
      IF (info .NE. 0) THEN
        CALL ppm_write(ppm_rank, 'ppm_poisson_solve','Failed to send field.',info2)
        GOTO 9999
      ENDIF
      !Retrieve
      CALL ppm_map_field_pop(&
      & ppmpoisson%topoidxy, &
      & ppmpoisson%meshidxy,ppmpoisson%fldxyr, &
      & 3,(/0,0,0/),info)
      IF (info .NE. 0) THEN
        CALL ppm_write(ppm_rank, 'ppm_poisson_solve','Failed to pop vector field.',info2)
        GOTO 9999
      ENDIF
      !-----------------------------------------------------------------------
      ! Do slab FFT (XY) - use the xy topology as its extent has not been halved
      !-----------------------------------------------------------------------
      CALL ppm_fft_execute_2d(ppmpoisson%topoidxy,&
      & ppmpoisson%meshidxy, ppmpoisson%planfxy, &
      & ppmpoisson%fldxyr, ppmpoisson%fldxyc, &
      & info)
      !-----------------------------------------------------------------------
      ! Map to the pencils (Z)
      !-----------------------------------------------------------------------
      !Initialise
      CALL ppm_map_field_global(&
      & ppmpoisson%topoidxy, &
      & ppmpoisson%topoidz, &
      & ppmpoisson%meshidxyc, &
      & ppmpoisson%meshidz,info)
      IF (info .NE. 0) THEN
        CALL ppm_write(ppm_rank, 'ppm_poisson_solve','Failed to initialise field mapping.',info2)
        GOTO 9999
      ENDIF
      !Push the data
      CALL ppm_map_field_push(&
      & ppmpoisson%topoidxy, &
      & ppmpoisson%meshidxyc,ppmpoisson%fldxyc,3,info)
      IF (info .NE. 0) THEN
        CALL ppm_write(ppm_rank, 'ppm_poisson_solve','Failed to push vector field.',info2)
        GOTO 9999
      ENDIF
      !Send
      CALL ppm_map_field_send(info)
      IF (info .NE. 0) THEN
        CALL ppm_write(ppm_rank, 'ppm_poisson_solve','Failed to send field.',info2)
        GOTO 9999
      ENDIF
      !Retrieve
      CALL ppm_map_field_pop(&
      & ppmpoisson%topoidz, &
      & ppmpoisson%meshidz,ppmpoisson%fldzc1, &
      & 3,(/0,0,0/),info)
      IF (info .NE. 0) THEN
        CALL ppm_write(ppm_rank, 'ppm_poisson_solve','Failed to pop vector field.',info2)
        GOTO 9999
      ENDIF
      !-----------------------------------------------------------------------
      ! Do pencil FFT (Z)
      !-----------------------------------------------------------------------
      CALL ppm_fft_execute_1d(ppmpoisson%topoidz,&
      & ppmpoisson%meshidz, ppmpoisson%planfz, &
      & ppmpoisson%fldzc1, ppmpoisson%fldzc2, &
      & info)
      !-----------------------------------------------------------------------
      ! HERE CHOOSE FIRST THE OPERATION
      !-----------------------------------------------------------------------
      ! Do vorticity reprojection if requested
      !-----------------------------------------------------------------------
      IF (presentoperation .EQ. ppm_poisson_opr_repr) THEN
        normfac = 1.0_MK/ REAL((ppmpoisson%nmfft(1))* & !vertex
                             & (ppmpoisson%nmfft(2))* &
                             & (ppmpoisson%nmfft(3)),MK)
        DO isub=1,ppmpoisson%nsublistz
          isubl=ppmpoisson%isublistz(isub)
          DO k=1,ppmpoisson%ndataz(3,isubl)
            gk = k - 1 + (ppmpoisson%istartz(3,isubl)-1)
            IF (gk .GT. (ppmpoisson%nmfft(3))/2) gk = gk-(ppmpoisson%nmfft(3))
            kz = CMPLX(0.0_MK,REAL(gk,MK),MK)*ppmpoisson%normkz
            DO j=1,ppmpoisson%ndataz(2,isubl)
              gj = j - 1 + (ppmpoisson%istartz(2,isubl)-1)
              IF (gj .GT. (ppmpoisson%nmfft(2))/2) gj = gj-(ppmpoisson%nmfft(2))
              ky = CMPLX(0.0_MK,REAL(gj,MK),MK)*ppmpoisson%normky
              DO i=1,ppmpoisson%ndataz(1,isubl)
                gi = i - 1 + (ppmpoisson%istartz(1,isubl)-1)
                IF (gi .GT. (ppmpoisson%nmfft(1))/2) gi = gi-(ppmpoisson%nmfft(1))
                kx = CMPLX(0.0_MK,REAL(gi,MK),MK)*ppmpoisson%normkx
                !IF (presentderivatives .EQ. ppm_poisson_drv_sp) THEN
                ! !this looks wacky
                ! ppmpoisson%fldzc2(1,i,j,k,isub) = &
                ! & -ppmpoisson%fldzc2(1,i,j,k,isub)*kx * &
                ! & ppmpoisson%fldgrnc( i,j,k,isub)
                ! ppmpoisson%fldzc2(2,i,j,k,isub) = &
                ! & -ppmpoisson%fldzc2(2,i,j,k,isub)*ky * &
                ! & ppmpoisson%fldgrnc( i,j,k,isub)
                ! ppmpoisson%fldzc2(3,i,j,k,isub) = &
                ! & -ppmpoisson%fldzc2(3,i,j,k,isub)*kz * &
                ! & ppmpoisson%fldgrnc( i,j,k,isub)
                !ELSE IF (ppmpoisson%bc .EQ. ppm_poisson_bc_per .AND. &
                IF (ppmpoisson%bc .EQ. ppm_poisson_bc_per .AND. &
                         presentderivatives .EQ. ppm_poisson_drv_sp) THEN
                  !compute spectral divergence....
                  divomega = (ppmpoisson%fldzc2(1,i,j,k,isub) * kx + &
                               & ppmpoisson%fldzc2(2,i,j,k,isub) * ky + &
                               & ppmpoisson%fldzc2(3,i,j,k,isub) * kz) * &
                               & ppmpoisson%fldgrnr( i,j,k,isub)
                  !...and subtract its gradient
                  ppmpoisson%fldzc2(1,i,j,k,isub) = &
                    & (ppmpoisson%fldzc2(1,i,j,k,isub)*normfac + divomega *kx)
                  ppmpoisson%fldzc2(2,i,j,k,isub) = &
                    & (ppmpoisson%fldzc2(2,i,j,k,isub)*normfac + divomega *ky)
                  ppmpoisson%fldzc2(3,i,j,k,isub) = &
                    & (ppmpoisson%fldzc2(3,i,j,k,isub)*normfac + divomega *kz)
                ELSE IF (ppmpoisson%bc .EQ. ppm_poisson_bc_fre .AND. &
                         presentderivatives .EQ. ppm_poisson_drv_sp) THEN
                  !compute spectral divergence and solve poisson....
                  divomega = (ppmpoisson%fldzc2(1,i,j,k,isub) * kx + &
                               & ppmpoisson%fldzc2(2,i,j,k,isub) * ky + &
                               & ppmpoisson%fldzc2(3,i,j,k,isub) * kz) * &
                               & ppmpoisson%fldgrnc( i,j,k,isub)
                  !...and subtract its gradient
                  ppmpoisson%fldzc2(1,i,j,k,isub) = &
                    & (ppmpoisson%fldzc2(1,i,j,k,isub)*normfac + divomega *kx)
                  ppmpoisson%fldzc2(2,i,j,k,isub) = &
                    & (ppmpoisson%fldzc2(2,i,j,k,isub)*normfac + divomega *ky)
                  ppmpoisson%fldzc2(3,i,j,k,isub) = &
                    & (ppmpoisson%fldzc2(3,i,j,k,isub)*normfac + divomega *kz)
                ENDIF
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      !-----------------------------------------------------------------------
      ! Or solve poisson equation for all components if requested
      !-----------------------------------------------------------------------
      ELSEIF (presentoperation .EQ. ppm_poisson_opr_none .OR. &
        & presentoperation .EQ. ppm_poisson_opr_vel) THEN
        !& presentoperation .EQ. ppm_poisson_opr_velk .OR. &
        !& presentoperation .EQ. ppm_poisson_opr_velrepr) THEN
        !-----------------------------------------------------------------------
        ! Apply the periodic Greens function
        !-----------------------------------------------------------------------
        IF (ppmpoisson%green .EQ. ppm_poisson_grn_pois_per) THEN
          DO isub=1,ppmpoisson%nsublistz
            isubl=ppmpoisson%isublistz(isub)
            DO k=1,ppmpoisson%ndataz(3,isubl)
              DO j=1,ppmpoisson%ndataz(2,isubl)
                DO i=1,ppmpoisson%ndataz(1,isubl)
                  ppmpoisson%fldzc2(1,i,j,k,isub) = ppmpoisson%fldgrnr( i,j,k,isub)*&
                                                  & ppmpoisson%fldzc2(1,i,j,k,isub)
                  ppmpoisson%fldzc2(2,i,j,k,isub) = ppmpoisson%fldgrnr( i,j,k,isub)*&
                                                  & ppmpoisson%fldzc2(2,i,j,k,isub)
                  ppmpoisson%fldzc2(3,i,j,k,isub) = ppmpoisson%fldgrnr( i,j,k,isub)*&
                                                  & ppmpoisson%fldzc2(3,i,j,k,isub)
                ENDDO
              ENDDO
            ENDDO
          ENDDO
        !-----------------------------------------------------------------------
        ! Apply the free-space Greens function
        !-----------------------------------------------------------------------
        ELSE IF (ppmpoisson%green .EQ. ppm_poisson_grn_pois_fre .OR. &
               & ppmpoisson%green .EQ. ppm_poisson_grn_pois_blob2 .OR. &
               & ppmpoisson%green .EQ. ppm_poisson_grn_pois_blob4 .OR. &
               & ppmpoisson%green .EQ. ppm_poisson_grn_pois_blob6 .OR. &
               & ppmpoisson%green .EQ. ppm_poisson_grn_pois_blob8 .OR. &
               & ppmpoisson%green .EQ. ppm_poisson_grn_pois_blob10) THEN
          DO isub=1,ppmpoisson%nsublistz
            isubl=ppmpoisson%isublistz(isub)
            DO k=1,ppmpoisson%ndataz(3,isubl)
              DO j=1,ppmpoisson%ndataz(2,isubl)
                DO i=1,ppmpoisson%ndataz(1,isubl)
                  ppmpoisson%fldzc2(1,i,j,k,isub) = ppmpoisson%fldgrnc( i,j,k,isub)*&
                                                  & ppmpoisson%fldzc2(1,i,j,k,isub)
                  ppmpoisson%fldzc2(2,i,j,k,isub) = ppmpoisson%fldgrnc( i,j,k,isub)*&
                                                  & ppmpoisson%fldzc2(2,i,j,k,isub)
                  ppmpoisson%fldzc2(3,i,j,k,isub) = ppmpoisson%fldgrnc( i,j,k,isub)*&
                                                  & ppmpoisson%fldzc2(3,i,j,k,isub)
                ENDDO
              ENDDO
            ENDDO
          ENDDO
        ENDIF
      ENDIF
      !-----------------------------------------------------------------------
      ! Spectral derivatives
      ! normkx, etc contains 2pi/Lx
      !-----------------------------------------------------------------------
      IF (presentoperation .EQ. ppm_poisson_opr_vel .AND. &
        & presentderivatives .EQ. ppm_poisson_drv_sp) THEN
        normfac = 1.0_MK/ REAL((ppmpoisson%nmfft(1))* & !vertex
                             & (ppmpoisson%nmfft(2))* &
                             & (ppmpoisson%nmfft(3)),MK)
        DO isub=1,ppmpoisson%nsublistz
          isubl=ppmpoisson%isublistz(isub)
          DO k=1,ppmpoisson%ndataz(3,isubl)
            gk = k - 1 + (ppmpoisson%istartz(3,isubl)-1)
            IF (gk .GT. (ppmpoisson%nmfft(3)/2)) gk = gk-(ppmpoisson%nmfft(3))
            kz = CMPLX(0.0_MK,REAL(gk,MK),MK)*ppmpoisson%normkz
            DO j=1,ppmpoisson%ndataz(2,isubl)
              gj = j - 1 + (ppmpoisson%istartz(2,isubl)-1)
              IF (gj .GT. (ppmpoisson%nmfft(2)/2)) gj = gj-(ppmpoisson%nmfft(2))
              ky = CMPLX(0.0_MK,REAL(gj,MK),MK)*ppmpoisson%normky
              DO i=1,ppmpoisson%ndataz(1,isubl)
                gi = i - 1 + (ppmpoisson%istartz(1,isubl)-1)
                IF (gi .GT. (ppmpoisson%nmfft(1)/2)) gi = gi-(ppmpoisson%nmfft(1))
                kx = CMPLX(0.0_MK,REAL(gi,MK),MK)*ppmpoisson%normkx
                phix = ppmpoisson%fldzc2(1,i,j,k,isub)
                phiy = ppmpoisson%fldzc2(2,i,j,k,isub)
                phiz = ppmpoisson%fldzc2(3,i,j,k,isub)
                ppmpoisson%fldzc2(1,i,j,k,isub) = (ky*phiz-kz*phiy)
                ppmpoisson%fldzc2(2,i,j,k,isub) = (kz*phix-kx*phiz)
                ppmpoisson%fldzc2(3,i,j,k,isub) = (kx*phiy-ky*phix)
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ENDIF
      !-----------------------------------------------------------------------
      ! IFFT pencil (Z)
      !-----------------------------------------------------------------------
      CALL ppm_fft_execute_1d(ppmpoisson%topoidz,&
      & ppmpoisson%meshidz, ppmpoisson%planbz, &
      & ppmpoisson%fldzc2, ppmpoisson%fldzc1, &
      & info)
      !-----------------------------------------------------------------------
      ! Map back to slabs (XY)
      !-----------------------------------------------------------------------
      !Initialise
      CALL ppm_map_field_global(&
      & ppmpoisson%topoidz, &
      & ppmpoisson%topoidxy, &
      & ppmpoisson%meshidz, &
      & ppmpoisson%meshidxyc,info)
      IF (info .NE. 0) THEN
        CALL ppm_write(ppm_rank, 'ppm_poisson_solve','Failed to initialise field mapping.',info2)
        GOTO 9999
      ENDIF
      !Push the data
      CALL ppm_map_field_push(&
      & ppmpoisson%topoidz, &
      & ppmpoisson%meshidz,ppmpoisson%fldzc1,3,info)
      IF (info .NE. 0) THEN
        CALL ppm_write(ppm_rank, 'ppm_poisson_solve','Failed to push vector field.',info2)
        GOTO 9999
      ENDIF
      !Send
      CALL ppm_map_field_send(info)
      IF (info .NE. 0) THEN
        CALL ppm_write(ppm_rank, 'ppm_poisson_solve','Failed to send field.',info2)
        GOTO 9999
      ENDIF
      !Retrieve
      CALL ppm_map_field_pop(&
      & ppmpoisson%topoidxy, &
      & ppmpoisson%meshidxyc,ppmpoisson%fldxyc, &
      & 3,(/0,0,0/),info)
      IF (info .NE. 0) THEN
        CALL ppm_write(ppm_rank, 'ppm_poisson_solve','Failed to pop vector field.',info2)
        GOTO 9999
      ENDIF
      !-----------------------------------------------------------------------
      ! IFFT (XY) use the non-reduced topology
      !-----------------------------------------------------------------------
      CALL ppm_fft_execute_2d(ppmpoisson%topoidxy,&
      & ppmpoisson%meshidxy, ppmpoisson%planbxy, &
      & ppmpoisson%fldxyc, ppmpoisson%fldxyr, &
      & info)
      !-----------------------------------------------------------------------
      ! Map back to standard topology (XYZ)
      !-----------------------------------------------------------------------
      !Initialise
      CALL ppm_map_field_global(&
      & ppmpoisson%topoidxy, &
      & topoid, &
      & ppmpoisson%meshidxy, &
      & meshid,info)
      IF (info .NE. 0) THEN
        CALL ppm_write(ppm_rank, 'ppm_poisson_solve','Failed to initialise field mapping.',info2)
        GOTO 9999
      ENDIF
      !Push the data
      CALL ppm_map_field_push(&
      & ppmpoisson%topoidxy, &
      & ppmpoisson%meshidxy,ppmpoisson%fldxyr,3,info)
      IF (info .NE. 0) THEN
        CALL ppm_write(ppm_rank, 'ppm_poisson_solve','Failed to push vector field.',info2)
        GOTO 9999
      ENDIF
      !Send
      CALL ppm_map_field_send(info)
      IF (info .NE. 0) THEN
        CALL ppm_write(ppm_rank, 'ppm_poisson_solve','Failed to send field.',info2)
        GOTO 9999
      ENDIF
      !-------------------------------------------------------------------------
      ! FINAL RETRIEVE - Here we do different things depending on the task
      ! i.e. the receiver varies
      !-------------------------------------------------------------------------
      ! * Map to temporary array to accomodate FD curl
      !-------------------------------------------------------------------------
      IF ((presentderivatives .EQ. ppm_poisson_drv_fd2 .OR. &
        & presentderivatives .EQ. ppm_poisson_drv_fd4) .AND. &
        & presentoperation .EQ. ppm_poisson_opr_vel) THEN
        CALL ppm_map_field_pop(&
        & topoid, &
        & meshid,ppmpoisson%drv_vr, &
        & 3,gstw,info)
        IF (info .NE. 0) THEN
          CALL ppm_write(ppm_rank, 'ppm_poisson_solve','Failed to pop vector field.',info2)
          GOTO 9999
        ENDIF
        !-------------------------------------------------------------------------
        ! Ghost the temporary array for derivatives (drv_vr)
        !-------------------------------------------------------------------------
        CALL ppm_map_field_ghost_get(topoid,meshid,gstw,info)
        IF (info .NE. 0) THEN
           CALL ppm_write(ppm_rank, 'ppm_poisson_solve','Failed to initialise ghosts.',info2)
           GOTO 9999
        ENDIF
        CALL ppm_map_field_push(topoid,meshid,ppmpoisson%drv_vr,3,info)
        IF (info .NE. 0) THEN
           CALL ppm_write(ppm_rank, 'ppm_poisson_solve','Failed to push ghosts.',info2)
           GOTO 9999
        ENDIF
        CALL ppm_map_field_send(info)
        IF (info .NE. 0) THEN
           CALL ppm_write(ppm_rank, 'ppm_poisson_solve','Failed to send ghosts.',info2)
           GOTO 9999
        ENDIF
        CALL ppm_map_field_pop(topoid,meshid,ppmpoisson%drv_vr,3,gstw,info)
        IF (info .NE. 0) THEN
           CALL ppm_write(ppm_rank, 'ppm_poisson_solve','Failed to pop ghosts.',info2)
           GOTO 9999
        ENDIF
      !-------------------------------------------------------------------------
      ! * For reprojection of continuous fields: Map back to the input array
      !-------------------------------------------------------------------------
      ELSE IF (presentoperation .EQ. ppm_poisson_opr_repr .AND. &
               presentderivatives .EQ. ppm_poisson_drv_sp) THEN
        CALL ppm_map_field_pop(&
        & topoid, &
        & meshid,fieldin, &
        & 3,gstw,info)
        IF (info .NE. 0) THEN
          CALL ppm_write(ppm_rank, 'ppm_poisson_solve','Failed to pop vector field.',info2)
          GOTO 9999
        ENDIF
      !-------------------------------------------------------------------------
      ! * For reprojection of dis-continuous fields: Map back to the output array
      !-------------------------------------------------------------------------
      ELSE IF (presentoperation .EQ. ppm_poisson_opr_repr .AND. &
            & (presentderivatives .EQ. ppm_poisson_drv_fd2 .OR. &
            & presentderivatives .EQ. ppm_poisson_drv_fd4)) THEN
        CALL ppm_map_field_pop(&
        & topoid, &
        & meshid,fieldout, &
        & 3,gstw,info)
        IF (info .NE. 0) THEN
          CALL ppm_write(ppm_rank, 'ppm_poisson_solve','Failed to pop vector field.',info2)
          GOTO 9999
        ENDIF
        !-----------------------------------------------------------------------
        ! Now subtract the correction (fieldout) from the vorticity (fieldin)
        !-----------------------------------------------------------------------
        CALL ppm_poisson_fd(topoid,meshid,fieldin,fieldout,&
                      & ppm_poisson_subtract, ppm_poisson_subtract,info)
      !-------------------------------------------------------------------------
      ! * Or map directly to the output array (spectral or no derivatives)
      !-------------------------------------------------------------------------
      ELSE IF (presentoperation .EQ. ppm_poisson_opr_none .OR. &
            & (presentoperation .EQ. ppm_poisson_opr_vel .AND. &
             & presentderivatives .EQ. ppm_poisson_drv_sp)) THEN
        CALL ppm_map_field_pop(&
        & topoid, &
        & meshid,fieldout, &
        & 3,gstw,info)
        IF (info .NE. 0) THEN
          CALL ppm_write(ppm_rank, 'ppm_poisson_solve','Failed to pop vector field.',info2)
          GOTO 9999
        ENDIF
      ELSE
        info = -1
        CALL ppm_write(ppm_rank, 'ppm_poisson_solve','Unsupported numerical setup.',info2)
      ENDIF
      !-------------------------------------------------------------------------
      ! In case of finite difference curl:
      !-------------------------------------------------------------------------
      ! Treat ghost layer to support FD stencils
      !-------------------------------------------------------------------------
      IF (presentoperation .EQ. ppm_poisson_opr_vel .AND. &
        & (presentderivatives .EQ. ppm_poisson_drv_fd2 .OR. &
        & presentderivatives .EQ. ppm_poisson_drv_fd4)) THEN
         CALL ppm_poisson_extrapolateghost(topoid,meshid,ppmpoisson%drv_vr,&
                                          & 2,4,gstw,info)
         CALL ppm_poisson_fd(topoid,meshid,ppmpoisson%drv_vr,fieldout,&
                           & ppm_poisson_curl,presentderivatives ,info)
      ENDIF
      !-------------------------------------------------------------------------
      ! Finally ghost the velocity/stream function field before returning it
      ! Also extrapolate if freespace
      !-------------------------------------------------------------------------
      CALL ppm_map_field_ghost_get(topoid,meshid,gstw,info)
      CALL ppm_map_field_push(topoid,meshid,fieldout,3,info)
      CALL ppm_map_field_send(info)
      CALL ppm_map_field_pop(topoid,meshid,fieldout,3,gstw,info)
      IF (ppmpoisson%bc .EQ. ppm_poisson_bc_fre) THEN
        CALL ppm_poisson_extrapolateghost(topoid,meshid,fieldout,2,4,gstw,info)
      ENDIF
      !-------------------------------------------------------------------------
      ! Perhaps allocate (and deallocate) arrays !@
      !-------------------------------------------------------------------------
      !-------------------------------------------------------------------------
      ! Return
      !-------------------------------------------------------------------------
 9999 CONTINUE
      CALL substop('ppm_poisson_solve',t0,info)
      RETURN
      END SUBROUTINE ppm_poisson_solve
      !-------------------------------------------------------------------------
      ! Subroutine : ppm_poisson_finalize.f90
      !-------------------------------------------------------------------------
      ! Copyright (c) 2010 CSE Lab (ETH Zurich), MOSAIC Group (ETH Zurich),
      ! Center for Fluid Dynamics (DTU)
      !
      !-------------------------------------------------------------------------
      SUBROUTINE ppm_poisson_finalize(ppmpoisson,info)
      !!! Routine to deallocate the data initialised for Greens function
      !!! solution of the Poisson
      IMPLICIT NONE
      !-------------------------------------------------------------------------
      ! Arguments
      !-------------------------------------------------------------------------
      TYPE(ppm_poisson_plan),INTENT(INOUT) :: ppmpoisson
      !!! The PPM Poisson plan type (inspired by the FFTW plan)
      INTEGER, INTENT(OUT) :: info
      !-------------------------------------------------------------------------
      ! Local variables
      !-------------------------------------------------------------------------
      REAL(ppm_kind_double) :: t0
      INTEGER :: stat,info2
      !-------------------------------------------------------------------------
      ! Initialise routine
      !-------------------------------------------------------------------------
      CALL substart('ppm_poisson_finalize',t0,info)
      IF(ASSOCIATED(ppmpoisson%drv_vr)) THEN
        DEALLOCATE(ppmpoisson%drv_vr,stat=info)
        IF (info .NE. 0) THEN
          CALL ppm_write(ppm_rank,'ppm_poisson_finalize','Failed deallocate drv_vr.',info2)
          GOTO 9999
        ENDIF
      ENDIF
      IF(ASSOCIATED(ppmpoisson%fldxyr)) THEN
        DEALLOCATE(ppmpoisson%fldxyr,stat=info)
        IF (info .NE. 0) THEN
          CALL ppm_write(ppm_rank,'ppm_poisson_finalize','Failed deallocate fldxyr.',info2)
          GOTO 9999
        ENDIF
      ENDIF
      IF(ASSOCIATED(ppmpoisson%fldxyc)) THEN
        DEALLOCATE(ppmpoisson%fldxyc,stat=info)
        IF (info .NE. 0) THEN
          CALL ppm_write(ppm_rank,'ppm_poisson_finalize','Failed deallocate fldxyc.',info2)
          GOTO 9999
        ENDIF
      ENDIF
      IF(ASSOCIATED(ppmpoisson%fldzc1)) THEN
        DEALLOCATE(ppmpoisson%fldzc1,stat=info)
        IF (info .NE. 0) THEN
          CALL ppm_write(ppm_rank,'ppm_poisson_finalize','Failed deallocate fldzc1.',info2)
          GOTO 9999
        ENDIF
      ENDIF
      IF(ASSOCIATED(ppmpoisson%fldzc2)) THEN
        DEALLOCATE(ppmpoisson%fldzc2,stat=info)
        IF (info .NE. 0) THEN
          CALL ppm_write(ppm_rank,'ppm_poisson_finalize','Failed deallocate fldzc2.',info2)
          GOTO 9999
        ENDIF
      ENDIF
      IF(ASSOCIATED(ppmpoisson%fldgrnr)) THEN
        DEALLOCATE(ppmpoisson%fldgrnr,stat=info)
        IF (info .NE. 0) THEN
          CALL ppm_write(ppm_rank,'ppm_poisson_finalize','Failed deallocate fldgrnr.',info2)
          GOTO 9999
        ENDIF
      ENDIF
      IF(ASSOCIATED(ppmpoisson%fldgrnc)) THEN
        DEALLOCATE(ppmpoisson%fldgrnc,stat=info)
        IF (info .NE. 0) THEN
          CALL ppm_write(ppm_rank,'ppm_poisson_finalize','Failed deallocate fldgrnc.',info2)
          GOTO 9999
        ENDIF
      ENDIF
      !the following topologies should be removed
      ! ppmpoisson%topoidxy
      ! ppmpoisson%meshidxy
      ! ppmpoisson%meshidxyc
      ! ppmpoisson%topoidz
      ! ppmpoisson%meshidz
      !Eventually the following fft plans must be deallocated
      ! ppmpoisson%planfxy
      ! ppmpoisson%planbxy
      ! ppmpoisson%planfz
      ! ppmpoisson%planbz
      !-------------------------------------------------------------------------
      ! Return
      !-------------------------------------------------------------------------
 9999 CONTINUE
      CALL substop('ppm_poisson_finalize',t0,info)
      RETURN
      END SUBROUTINE ppm_poisson_finalize
      !-------------------------------------------------------------------------
      ! Subroutine : ppm_poisson_fd.f90
      !-------------------------------------------------------------------------
      ! Copyright (c) 2010 CSE Lab (ETH Zurich), MOSAIC Group (ETH Zurich),
      ! Center for Fluid Dynamics (DTU)
      !
      !-------------------------------------------------------------------------
      SUBROUTINE ppm_poisson_fd(topoid,meshid,fieldin,fieldout,dtype,order,info)
      !!! This routine computes the finite difference gradients, curl etc of
      !!! fieldin and outputs to fieldout. Both in and out fields are on
      !!! the mesh with meshid belonging to the topology with id topoid.
      !!! The finite difference to be carried out is determined by dtype which
      !!! must be one of the following:
      !!! * ppm_poisson_drv_curl_fd2
      !!! * ppm_poisson_drv_grad_fd2 (not implemented yet)
      !!! * ppm_poisson_drv_lapl_fd2 (not implemented yet)
      !!! * ppm_poisson_drv_div_fd2
      !!! * ppm_poisson_drv_curl_fd4
      !!! * ppm_poisson_drv_grad_fd4 (not implemented yet)
      !!! * ppm_poisson_drv_lapl_fd4 (not implemented yet)
      !!! * ppm_poisson_drv_div_fd4
      !!!
      !!! [NOTE] fieldin and fieldout must NOT be the same array. A check
      !!! should be added.
      !@ TODO: Somewhere check if fieldin is equal to fieldout and give a warning
      USE ppm_module_topo_get
      IMPLICIT NONE
      !-------------------------------------------------------------------------
      ! Arguments
      !-------------------------------------------------------------------------
      INTEGER, INTENT(IN) :: topoid
      !!! ID of the topology
      INTEGER, INTENT(IN) :: meshid
      !!! Mesh ID
      REAL(ppm_kind_double),DIMENSION(:,:,:,:,:),POINTER :: fieldin
      !!! Input field data
      REAL(ppm_kind_double),DIMENSION(:,:,:,:,:),POINTER :: fieldout
      !!! Output field data
      INTEGER, INTENT(IN) :: dtype
      !!! Derivation type. Can be one of the types:
      !!! * ppm_poisson_curl
      !!! * ppm_poisson_divergence
      !!! * ppm_poisson_subtract
      INTEGER, INTENT(IN) :: order
      !!! Order of the stencil:
      !!! * ppm_poisson_drv_fd2
      !!! * ppm_poisson_drv_fd4
      INTEGER, INTENT(OUT) :: info
      !!! Return status, 0 upon succes
      !-------------------------------------------------------------------------
      ! Local variables
      !-------------------------------------------------------------------------
      INTEGER,PARAMETER :: MK = ppm_kind_double
      REAL(ppm_kind_double) :: t0
      TYPE(ppm_t_topo),POINTER :: topology
      TYPE(ppm_t_equi_mesh) :: mesh
      REAL(ppm_kind_double) :: dx,dy,dz
      REAL(ppm_kind_double) :: facx,facy,facz
      REAL(ppm_kind_double) :: divergence
      INTEGER :: isub,isubl
      INTEGER :: i,j,k
      !-------------------------------------------------------------------------
      ! Initialise routine
      !-------------------------------------------------------------------------
      CALL substart('ppm_poisson_fd',t0,info)
      !-------------------------------------------------------------------------
      ! Get topology and mesh values
      !-------------------------------------------------------------------------
      CALL ppm_topo_get(topoid,topology,info)
      IF (info .NE. 0) THEN
        CALL ppm_write(ppm_rank,'ppm_poisson_init_predef','Failed to get topology.',isub)
        GOTO 9999
      ENDIF
      mesh = topology%mesh(meshid)
      dx = (topology%max_physd(1)-topology%min_physd(1))/REAL(mesh%nm(1)-1) !vertex
      dy = (topology%max_physd(2)-topology%min_physd(2))/REAL(mesh%nm(2)-1)
      dz = (topology%max_physd(3)-topology%min_physd(3))/REAL(mesh%nm(3)-1)
      !-----------------------------------------------------------------------
      ! Do the finite difference calculation
      !-----------------------------------------------------------------------
      !-----------------------------------------------------------------------
      ! Curl, 2nd order FD
      !-----------------------------------------------------------------------
      IF (dtype .EQ. ppm_poisson_curl .AND. &
        & order .EQ. ppm_poisson_drv_fd2) THEN
      write(*,*) 'position 1'
        facx = 1.0_MK/(2.0_MK*dx)
        facy = 1.0_MK/(2.0_MK*dy)
        facz = 1.0_MK/(2.0_MK*dz)
        DO isub=1,topology%nsublist
          isubl=topology%isublist(isub)
          DO k=1,mesh%nnodes(3,isubl)
            DO j=1,mesh%nnodes(2,isubl)
              DO i=1,mesh%nnodes(1,isubl)
                fieldout(1,i,j,k,isub) = &
                & facy*(fieldin(3,i ,j+1,k ,isub)- &
                       & fieldin(3,i ,j-1,k ,isub)) &
                & -facz*(fieldin(2,i ,j ,k+1,isub)- &
                       & fieldin(2,i ,j ,k-1,isub))
                fieldout(2,i,j,k,isub) = &
                & facz*(fieldin(1,i ,j ,k+1,isub)- &
                       & fieldin(1,i ,j ,k-1,isub)) &
                & -facx*(fieldin(3,i+1,j ,k ,isub)- &
                       & fieldin(3,i-1,j ,k ,isub))
                fieldout(3,i,j,k,isub) = &
                & facx*(fieldin(2,i+1,j ,k ,isub)- &
                       & fieldin(2,i-1,j ,k ,isub)) &
                & -facy*(fieldin(1,i ,j+1,k ,isub)- &
                       & fieldin(1,i ,j-1,k ,isub))
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      !-----------------------------------------------------------------------
      ! Curl, 4th order FD
      !-----------------------------------------------------------------------
      ELSE IF (dtype .EQ. ppm_poisson_curl .AND. &
             & order .EQ. ppm_poisson_drv_fd4) THEN
      write(*,*) 'position 2'
        facx = 1.0_MK/(12.0_MK*dx)
        facy = 1.0_MK/(12.0_MK*dy)
        facz = 1.0_MK/(12.0_MK*dz)
        DO isub=1,topology%nsublist
          isubl=topology%isublist(isub)
          DO k=1,mesh%nnodes(3,isubl)
            DO j=1,mesh%nnodes(2,isubl)
              DO i=1,mesh%nnodes(1,isubl)
                fieldout(1,i,j,k,isub) = &
                & facy*( -fieldin(3,i ,j+2,k ,isub) &
                       & +8.0_MK*fieldin(3,i ,j+1,k ,isub) &
                       & -8.0_MK*fieldin(3,i ,j-1,k ,isub) &
                       & +fieldin(3,i ,j-2,k ,isub)) &
                & -facz*( -fieldin(2,i ,j ,k+2,isub) &
                       & +8.0_MK*fieldin(2,i ,j ,k+1,isub) &
                       & -8.0_MK*fieldin(2,i ,j ,k-1,isub) &
                       & +fieldin(2,i ,j ,k-2,isub))
                fieldout(2,i,j,k,isub) = &
                & facz*( -fieldin(1,i ,j ,k+2,isub) &
                       & +8.0_MK*fieldin(1,i ,j ,k+1,isub) &
                       & -8.0_MK*fieldin(1,i ,j ,k-1,isub) &
                       & +fieldin(1,i ,j ,k-2,isub)) &
                & -facx*( -fieldin(3,i+2,j ,k ,isub) &
                       & +8.0_MK*fieldin(3,i+1,j ,k ,isub) &
                       & -8.0_MK*fieldin(3,i-1,j ,k ,isub) &
                       & +fieldin(3,i-2,j ,k ,isub))
                fieldout(3,i,j,k,isub) = &
                & facx*( -fieldin(2,i+2,j ,k ,isub) &
                       & +8.0_MK*fieldin(2,i+1,j ,k ,isub) &
                       & -8.0_MK*fieldin(2,i-1,j ,k ,isub) &
                       & +fieldin(2,i-2,j ,k ,isub)) &
                & -facy*( -fieldin(1,i ,j+2,k ,isub) &
                       & +8.0_MK*fieldin(1,i ,j+1,k ,isub) &
                       & -8.0_MK*fieldin(1,i ,j-1,k ,isub) &
                       & +fieldin(1,i ,j-2,k ,isub))
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      !-----------------------------------------------------------------------
      ! Divergence, 2nd order FD
      !-----------------------------------------------------------------------
      ELSE IF (dtype .EQ. ppm_poisson_divergence .AND. &
             & order .EQ. ppm_poisson_drv_fd2) THEN
      write(*,*) 'position 3'
        facx = 1.0_MK/(2.0_MK*dx)
        facy = 1.0_MK/(2.0_MK*dy)
        facz = 1.0_MK/(2.0_MK*dz)
        DO isub=1,topology%nsublist
          isubl=topology%isublist(isub)
          DO k=1+topology%subs_bc(5,isubl)*1, &
                         & mesh%nnodes(3,isubl)-topology%subs_bc(6,isubl)*1
            DO j=1+topology%subs_bc(3,isubl)*1, &
                           & mesh%nnodes(2,isubl)-topology%subs_bc(4,isubl)*1
              DO i=1+topology%subs_bc(1,isubl)*1, &
                             & mesh%nnodes(1,isubl)-topology%subs_bc(2,isubl)*1
                divergence = &
                & facx*(fieldin(1,i+1,j ,k ,isub)- &
                       & fieldin(1,i-1,j ,k ,isub)) &
                & +facy*(fieldin(2,i ,j+1,k ,isub)- &
                       & fieldin(2,i ,j-1,k ,isub)) &
                & +facz*(fieldin(3,i ,j ,k+1,isub)- &
                       & fieldin(3,i ,j ,k-1,isub))
                fieldout(1,i,j,k,isub) = divergence
                fieldout(2,i,j,k,isub) = divergence
                fieldout(3,i,j,k,isub) = divergence
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      !-----------------------------------------------------------------------
      ! Divergence, 4th order FD
      !-----------------------------------------------------------------------
      ELSE IF (dtype .EQ. ppm_poisson_divergence .AND. &
             & order .EQ. ppm_poisson_drv_fd4) THEN
      write(*,*) 'position 4'
        facx = 1.0_MK/(12.0_MK*dx)
        facy = 1.0_MK/(12.0_MK*dy)
        facz = 1.0_MK/(12.0_MK*dz)
        DO isub=1,topology%nsublist
          isubl=topology%isublist(isub)
          DO k=1+topology%subs_bc(5,isubl)*2, &
                         & mesh%nnodes(3,isubl)-topology%subs_bc(6,isubl)*2
            DO j=1+topology%subs_bc(3,isubl)*2, &
                           & mesh%nnodes(2,isubl)-topology%subs_bc(4,isubl)*2
              DO i=1+topology%subs_bc(1,isubl)*2, &
                             & mesh%nnodes(1,isubl)-topology%subs_bc(2,isubl)*2
                divergence = &
                & facx*( -fieldin(1,i+2,j ,k ,isub) &
                       & +8.0_MK*fieldin(1,i+1,j ,k ,isub) &
                       & -8.0_MK*fieldin(1,i-1,j ,k ,isub) &
                       & +fieldin(1,i-2,j ,k ,isub)) &
                & +facy*( -fieldin(2,i ,j+2,k ,isub) &
                       & +8.0_MK*fieldin(2,i ,j+1,k ,isub) &
                       & -8.0_MK*fieldin(2,i ,j-1,k ,isub) &
                       & +fieldin(2,i ,j-2,k ,isub)) &
                & +facz*( -fieldin(3,i ,j ,k+2,isub) &
                       & +8.0_MK*fieldin(3,i ,j ,k+1,isub) &
                       & -8.0_MK*fieldin(3,i ,j ,k-1,isub) &
                       & +fieldin(3,i ,j ,k-2,isub))
                fieldout(1,i,j,k,isub) = divergence
                fieldout(2,i,j,k,isub) = divergence
                fieldout(3,i,j,k,isub) = divergence
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      !-----------------------------------------------------------------------
      ! Subtract fieldout from fieldin
      !-----------------------------------------------------------------------
      ELSE IF (dtype .EQ. ppm_poisson_subtract) THEN
      write(*,*) 'position 5'
        DO isub=1,topology%nsublist
          isubl=topology%isublist(isub)
          DO k=1,mesh%nnodes(3,isubl)
            DO j=1,mesh%nnodes(2,isubl)
              DO i=1,mesh%nnodes(1,isubl)
                fieldin(1,i,j,k,isub) = fieldin(1,i,j,k,isub) - &
                                        fieldout(1,i,j,k,isub)
                fieldin(2,i,j,k,isub) = fieldin(2,i,j,k,isub) - &
                                        fieldout(2,i,j,k,isub)
                fieldin(3,i,j,k,isub) = fieldin(3,i,j,k,isub) - &
                                        fieldout(3,i,j,k,isub)
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ENDIF
      !-------------------------------------------------------------------------
      ! Return
      !-------------------------------------------------------------------------
 9999 CONTINUE
      CALL substop('ppm_poisson_fd',t0,info)
      RETURN
      END SUBROUTINE ppm_poisson_fd
      !-------------------------------------------------------------------------
      ! Subroutine : ppm_poisson_extrapolateghost.f90
      !-------------------------------------------------------------------------
      ! Copyright (c) 2010 CSE Lab (ETH Zurich), MOSAIC Group (ETH Zurich),
      ! Center for Fluid Dynamics (DTU)
      !
      !-------------------------------------------------------------------------
      SUBROUTINE ppm_poisson_extrapolateghost_vr(topoid,meshid,field,nextra,nbase,gstw,info)
      !!! This routine extrapolates field values of field with topology and mesh
      !!! id topoid,meshid respectively into the ghost layer of width gstw.
      !!! nbase points are used to extrapolate into nextra points.
      !!!
      !!! [NOTE]
      !!! Presently extrapolation can only be done to 1 or two points into the
      !!! ghostlayer always based on 4 points (fourth order spatial convergence)
      !!! A general nbase,nextra extrapolation can be implemented vi solution
      !!! of a small linear system of equations. This has not been done.
      !!! This routine is in need of loop unrollling in particular for
      !!! typical choices of nbase,nextra pairs. Extrapolation is necessary for
      !!! e.g. freespace FD curl of stream function
      USE ppm_module_topo_get
      IMPLICIT NONE
      !-------------------------------------------------------------------------
      ! Arguments
      !-------------------------------------------------------------------------
      INTEGER, INTENT(IN) :: topoid
      !!! Topology id of the field
      INTEGER, INTENT(IN) :: meshid
      !!! Mesh id of the field
      REAL(ppm_kind_double),DIMENSION(:,:,:,:,:),POINTER :: field
      !!! Field to extrapolate ghost layer
      INTEGER, INTENT(IN) :: nextra
      !!! Number of points to extrapolate into the ghostlayer
      INTEGER, INTENT(IN) :: nbase
      !!! Number of points to base the extrapolation on
      INTEGER,DIMENSION(3),INTENT(IN) :: gstw
      !!! Width of the ghotslayer
      INTEGER, INTENT(OUT) :: info
      !!! Return state
      !-------------------------------------------------------------------------
      ! Local variables
      !-------------------------------------------------------------------------
      INTEGER,PARAMETER :: MK = ppm_kind_double
      REAL(ppm_kind_double) :: t0
      TYPE(ppm_t_topo),POINTER :: topology
      TYPE(ppm_t_equi_mesh) :: mesh
      INTEGER :: isub,isubl
      INTEGER :: i,j,k,iextra,ibase
      REAL(ppm_kind_double),DIMENSION(:,:),POINTER :: coeff
      REAL(ppm_kind_double),DIMENSION(3) :: tmpbuf
      !-------------------------------------------------------------------------
      ! Initialise routine
      !-------------------------------------------------------------------------
      CALL substart('ppm_poisson_extrapolateghost',t0,info)
      !-------------------------------------------------------------------------
      ! Compare the number of points to extrapolate to the ghost layer width
      !-------------------------------------------------------------------------
      IF (nextra .GT. gstw(1) .OR. &
        & nextra .GT. gstw(2) .OR. &
        & nextra .GT. gstw(3)) THEN
         CALL ppm_write(ppm_rank,'ppm_poisson_extrapolateghost',&
         & 'The points to extrapolate exceeds the ghost layer.',info)
         info = -1
         GOTO 9999
      ENDIF
      !-------------------------------------------------------------------------
      ! Determine weights
      !-------------------------------------------------------------------------
      ALLOCATE(coeff(nbase,nextra))
      IF (nbase .EQ. 4) THEN
        IF (nextra .GE. 1) THEN
          coeff(:,1) = (/4.0_MK,-6.0_MK,4.0_MK,-1.0_MK/)
        ENDIF
        IF (nextra .GE. 2) THEN
          coeff(:,2) = (/10.0_MK,-20.0_MK,15.0_MK,-4.0_MK/)
        ENDIF
        IF (nextra .GE. 3) THEN
          CALL ppm_write(ppm_rank,'ppm_poisson_extrapolateghost',&
          & 'Extrapolation to more than two points has not been implemented.',info)
          info = -1
          GOTO 9999
        ENDIF
      ELSE
        CALL ppm_write(ppm_rank,'ppm_poisson_extrapolateghost',&
        & 'Only extrapolation based on 4 points has been implemented.',info)
        info = -1
        GOTO 9999
      ENDIF
      !-------------------------------------------------------------------------
      ! Get topology
      !-------------------------------------------------------------------------
      CALL ppm_topo_get(topoid,topology,info)
      IF (info .NE. 0) THEN
         CALL ppm_write(ppm_rank,'ppm_poisson_extrapolateghost',&
         & 'Failed to get topology.',isub)
         GOTO 9999
      ENDIF
      mesh = topology%mesh(meshid)
      !-------------------------------------------------------------------------
      ! Extrapolate field into ghost layer
      ! The indicies of subs_bc represent:
      ! west,east(x),south,north(y),bottom,top(z)
      !@ Some more unrolling here would be nice
      !-------------------------------------------------------------------------
      DO isub=1,topology%nsublist
         isubl=topology%isublist(isub)
         !West (-x)
         IF (topology%subs_bc(1,isubl) .EQ. 1) THEN
            DO k=1-gstw(3),mesh%nnodes(3,isubl)+gstw(3)
               DO j=1-gstw(2),mesh%nnodes(2,isubl)+gstw(2)
                  !DO i=1-gstw(1),mesh%nnodes(1,isubl)+gstw(1)
                  i = 1
                  DO iextra=1,nextra
                     tmpbuf = 0.0_MK
                     DO ibase=0,nbase-1
                        tmpbuf(1) = tmpbuf(1) + &
                         & coeff(ibase+1,iextra)*field(1,i+ibase,j,k,isub)
                        tmpbuf(2) = tmpbuf(2) + &
                         & coeff(ibase+1,iextra)*field(2,i+ibase,j,k,isub)
                        tmpbuf(3) = tmpbuf(3) + &
                         & coeff(ibase+1,iextra)*field(3,i+ibase,j,k,isub)
                     END DO !ibase
                     field(1,i-iextra,j,k,isub) = tmpbuf(1)
                     field(2,i-iextra,j,k,isub) = tmpbuf(2)
                     field(3,i-iextra,j,k,isub) = tmpbuf(3)
                  END DO !iextra
               ENDDO !j
            ENDDO !k
         ENDIF
         !East (+x)
         IF (topology%subs_bc(2,isubl) .EQ. 1) THEN
            DO k=1-gstw(3),mesh%nnodes(3,isubl)+gstw(3)
               DO j=1-gstw(2),mesh%nnodes(2,isubl)+gstw(2)
                  !DO i=1-gstw(1),mesh%nnodes(1,isubl)+gstw(1)
                  i = mesh%nnodes(1,isubl)
                  DO iextra=1,nextra
                     tmpbuf = 0.0_MK
                     DO ibase=0,nbase-1
                        tmpbuf(1) = tmpbuf(1) + &
                         & coeff(ibase+1,iextra)*field(1,i-ibase,j,k,isub)
                        tmpbuf(2) = tmpbuf(2) + &
                         & coeff(ibase+1,iextra)*field(2,i-ibase,j,k,isub)
                        tmpbuf(3) = tmpbuf(3) + &
                         & coeff(ibase+1,iextra)*field(3,i-ibase,j,k,isub)
                     END DO !ibase
                     field(1,i+iextra,j,k,isub) = tmpbuf(1)
                     field(2,i+iextra,j,k,isub) = tmpbuf(2)
                     field(3,i+iextra,j,k,isub) = tmpbuf(3)
                  END DO !iextra
               ENDDO !j
            ENDDO !k
         ENDIF
         !South (-y)
         IF (topology%subs_bc(3,isubl) .EQ. 1) THEN
            DO k=1-gstw(3),mesh%nnodes(3,isubl)+gstw(3)
               !DO j=1-gstw(2),mesh%nnodes(2,isubl)+gstw(2)
               DO i=1-gstw(1),mesh%nnodes(1,isubl)+gstw(1)
                  j = 1
                  DO iextra=1,nextra
                     tmpbuf = 0.0_MK
                     DO ibase=0,nbase-1
                        tmpbuf(1) = tmpbuf(1) + &
                         & coeff(ibase+1,iextra)*field(1,i,j+ibase,k,isub)
                        tmpbuf(2) = tmpbuf(2) + &
                         & coeff(ibase+1,iextra)*field(2,i,j+ibase,k,isub)
                        tmpbuf(3) = tmpbuf(3) + &
                         & coeff(ibase+1,iextra)*field(3,i,j+ibase,k,isub)
                     END DO !ibase
                     field(1,i,j-iextra,k,isub) = tmpbuf(1)
                     field(2,i,j-iextra,k,isub) = tmpbuf(2)
                     field(3,i,j-iextra,k,isub) = tmpbuf(3)
                  END DO !iextra
               ENDDO !i
            ENDDO !k
         ENDIF
         !North (+y)
         IF (topology%subs_bc(4,isubl) .EQ. 1) THEN
            DO k=1-gstw(3),mesh%nnodes(3,isubl)+gstw(3)
               !DO j=1-gstw(2),mesh%nnodes(2,isubl)+gstw(2)
               DO i=1-gstw(1),mesh%nnodes(1,isubl)+gstw(1)
                  j = mesh%nnodes(2,isubl)
                  DO iextra=1,nextra
                     tmpbuf = 0.0_MK
                     DO ibase=0,nbase-1
                        tmpbuf(1) = tmpbuf(1) + &
                         & coeff(ibase+1,iextra)*field(1,i,j-ibase,k,isub)
                        tmpbuf(2) = tmpbuf(2) + &
                         & coeff(ibase+1,iextra)*field(2,i,j-ibase,k,isub)
                        tmpbuf(3) = tmpbuf(3) + &
                         & coeff(ibase+1,iextra)*field(3,i,j-ibase,k,isub)
                     END DO !ibase
                     field(1,i,j+iextra,k,isub) = tmpbuf(1)
                     field(2,i,j+iextra,k,isub) = tmpbuf(2)
                     field(3,i,j+iextra,k,isub) = tmpbuf(3)
                  END DO !iextra
               ENDDO !i
            ENDDO !k
         ENDIF
         !Bottom (-z)
         IF (topology%subs_bc(5,isubl) .EQ. 1) THEN
            !DO k=1-gstw(3),mesh%nnodes(3,isubl)+gstw(3)
            DO j=1-gstw(2),mesh%nnodes(2,isubl)+gstw(2)
               DO i=1-gstw(1),mesh%nnodes(1,isubl)+gstw(1)
                  k = 1
                  DO iextra=1,nextra
                     tmpbuf = 0.0_MK
                     DO ibase=0,nbase-1
                        tmpbuf(1) = tmpbuf(1) + &
                         & coeff(ibase+1,iextra)*field(1,i,j,k+ibase,isub)
                        tmpbuf(2) = tmpbuf(2) + &
                         & coeff(ibase+1,iextra)*field(2,i,j,k+ibase,isub)
                        tmpbuf(3) = tmpbuf(3) + &
                         & coeff(ibase+1,iextra)*field(3,i,j,k+ibase,isub)
                     END DO !ibase
                     field(1,i,j,k-iextra,isub) = tmpbuf(1)
                     field(2,i,j,k-iextra,isub) = tmpbuf(2)
                     field(3,i,j,k-iextra,isub) = tmpbuf(3)
                  END DO !iextra
               ENDDO !i
            ENDDO !j
         ENDIF
         !Top (+z)
         IF (topology%subs_bc(6,isubl) .EQ. 1) THEN
            !DO k=1-gstw(3),mesh%nnodes(3,isubl)+gstw(3)
            DO j=1-gstw(2),mesh%nnodes(2,isubl)+gstw(2)
               DO i=1-gstw(1),mesh%nnodes(1,isubl)+gstw(1)
                  k = mesh%nnodes(3,isubl)
                  DO iextra=1,nextra
                     tmpbuf = 0.0_MK
                     DO ibase=0,nbase-1
                        tmpbuf(1) = tmpbuf(1) + &
                         & coeff(ibase+1,iextra)*field(1,i,j,k-ibase,isub)
                        tmpbuf(2) = tmpbuf(2) + &
                         & coeff(ibase+1,iextra)*field(2,i,j,k-ibase,isub)
                        tmpbuf(3) = tmpbuf(3) + &
                         & coeff(ibase+1,iextra)*field(3,i,j,k-ibase,isub)
                     END DO !ibase
                     field(1,i,j,k+iextra,isub) = tmpbuf(1)
                     field(2,i,j,k+iextra,isub) = tmpbuf(2)
                     field(3,i,j,k+iextra,isub) = tmpbuf(3)
                  END DO !iextra
               ENDDO !i
            ENDDO !j
         ENDIF
      ENDDO !isub
 9999 CONTINUE
      CALL substop('ppm_poisson_extrapolateghost',t0,info)
      RETURN
      END SUBROUTINE ppm_poisson_extrapolateghost_vr
      END MODULE ppm_module_poisson
