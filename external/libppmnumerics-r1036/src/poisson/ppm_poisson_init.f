      !-------------------------------------------------------------------------
      !  Subroutine   : ppm_poisson_init.f90
      !-------------------------------------------------------------------------
      ! Copyright (c) 2010 CSE Lab (ETH Zurich), MOSAIC Group (ETH Zurich),
      !                    Center for Fluid Dynamics (DTU)
      !
      !-------------------------------------------------------------------------
      SUBROUTINE __ROUTINE(topoid,meshid,ppmpoisson,fieldin,fieldout,green,info&
                                &,bc,derive)
      !!! Routine to initialise Greens function solution of the Poisson
      !!! equation. green is the flag defining which Greens function to use:
      !!! * ppm_poisson_grn_pois_per - Poisson equation, periodic boundaries
      !!! * ppm_poisson_grn_pois_fre - Poisson equation, freespace boundaries (not implemented)
      !!! * ppm_poisson_grn_reprojec - Do vorticity reprojection to kill divergence
      !!! Eventually the routine should be overloaded to accept custom Greens
      !!! functions such that more general convolutions can be performed.
      !!! green should be expanded to include more buildin Greens functions.
      !!!
      !!! The routine should accept an optional flag to toggle deallocation of
      !!! work arrays between calls to ppm_poisson_solve
      !!!
      !!! [NOTE]
      !!! When meshid > 1 works with the mapping the memory consumption can be
      !!! halved. Sketches have been implemented.
      !!! This routine should also be renamed to _init


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
      INTEGER, INTENT(IN)                                         :: topoid
      !!! Topology ID
      INTEGER, INTENT(IN)                                         :: meshid
      !!! Mesh ID
      TYPE(ppm_poisson_plan),INTENT(INOUT)                        :: ppmpoisson
      !!! The PPM Poisson plan type (inspired by the FFTW plan)
      REAL(__PREC),DIMENSION(:,:,:,:,:),POINTER                   :: fieldin
      !!! Input data field
      !@ strictly speaking fieldin is not being used in the init routine
      REAL(__PREC),DIMENSION(:,:,:,:,:),POINTER                   :: fieldout
      !!! Output data field
      INTEGER, INTENT(IN)                                         :: green
      !!!flag to select build-in Greens functions:
      !!!ppm_poisson_grn_pois_per - Poisson equation, periodic boundaries
      !!!ppm_poisson_grn_pois_fre - Poisson equation, freespace boundaries
      !!!ppm_poisson_grn_reprojec - Do vorticity reprojection to kill divergence
      !!!
      !!!Eventually this should also accept custom Greens function
      INTEGER, INTENT(OUT)                                        :: info
      INTEGER,INTENT(IN),OPTIONAL                                 :: bc
      !!!boundary condition for the convolution. Can be on of the following:
      !!!ppm_poisson_grn_pois_per, ppm_poisson_grn_pois_fre.
      !!!One could argue that this is redundant in the build-in case
      !@ Isnt this redundant because of green?
      INTEGER,INTENT(IN),OPTIONAL                                 :: derive
      !!!flag to toggle various derivatives of the solution (not to be used with
      !!!green=ppm_poisson_grn_reprojec):
      !!! * ppm_poisson_drv_none
      !!! * ppm_poisson_drv_curl_sp  (only for periodic BC)
      !!! * ppm_poisson_drv_grad_sp  (not implemented)
      !!! * ppm_poisson_drv_lapl_sp  (not implemented)
      !!! * ppm_poisson_drv_div_sp   (not implemented)
      !!! * ppm_poisson_drv_curl_fd2
      !!! * ppm_poisson_drv_grad_fd2 (not implemented)
      !!! * ppm_poisson_drv_lapl_fd2 (not implemented)
      !!! * ppm_poisson_drv_div_fd2  (not implemented)
      !!! * ppm_poisson_drv_curl_fd4
      !!! * ppm_poisson_drv_grad_fd4 (not implemented)
      !!! * ppm_poisson_drv_lapl_fd4 (not implemented)
      !!! * ppm_poisson_drv_div_fd4  (not implemented)
      !!!
      !!! curl=curl, grad=gradient,lapl=laplace operator,div=divergence
      !!! sp=spectrally, fd2=2nd order finite differences, fd4=4th order FD
      !!!
      !!!The spectral derivatives can only be computed in this routine. Since 
      !!!the flag exists finite difference derivatives have also been included.

      !-------------------------------------------------------------------------
      ! Local variables
      !-------------------------------------------------------------------------
      REAL(__PREC)                            :: t0
      REAL(__PREC),DIMENSION(:,:),POINTER     :: xp      !particle positions
      TYPE(ppm_t_topo),POINTER                :: topology,topologyxy,topologyz
      !!TYPE(ppm_t_topo),POINTER                :: topologyxyc
      TYPE(ppm_t_equi_mesh)                   :: mesh
      INTEGER ,DIMENSION(__DIM)               :: indl,indu
      INTEGER,PARAMETER                       :: MK = __PREC
      REAL(__PREC),PARAMETER                  :: PI=ACOS(-1.0_MK) !@ use ppm pi
      !factor for the Greens function, including FFT normalization
      REAL(__PREC)                            :: normfac
      INTEGER                                 :: i,j,k
      INTEGER                                 :: kx,ky,kz
      INTEGER                                 :: isubl,isub
      INTEGER,DIMENSION(__DIM*2)              :: bcdef
      INTEGER                                 :: assigning
      INTEGER                                 :: decomposition
      INTEGER,SAVE                            :: ttopoid
      INTEGER                                 :: tmeshid
      INTEGER , DIMENSION(:  ), POINTER       :: isublist => NULL()
      INTEGER                                 :: nsublist
      REAL(__PREC)                            :: dx,dy,dz
      REAL(__PREC)                            :: Lx2,Ly2,Lz2

      REAL(__PREC),DIMENSION(__DIM)           :: tmpmin,tmpmax

      !-------------------------------------------------------------------------
      ! Initialise routine
      !-------------------------------------------------------------------------
      CALL substart('ppm_poisson_init',t0,info)

      !-------------------------------------------------------------------------
      ! Investigate optional arguments, setup routine accordingly
      ! !@ Also check if the input/output and derivatives match
      !-------------------------------------------------------------------------
      IF (green .EQ. ppm_poisson_grn_pois_per) THEN
        ppmpoisson%case  = ppm_poisson_grn_pois_per
      ELSE IF (green .EQ. ppm_poisson_grn_pois_fre) THEN
        ppmpoisson%case  = ppm_poisson_grn_pois_fre
      ELSE IF (green .EQ. ppm_poisson_grn_reprojec) THEN
        ppmpoisson%case  = ppm_poisson_grn_reprojec
      ENDIF

      !-------------------------------------------------------------------------
      ! Get topology and mesh values
      !-------------------------------------------------------------------------
      CALL ppm_topo_get(topoid,topology,info)
      IF (info .NE. 0) THEN
        CALL ppm_write(ppm_rank,'ppm_poisson_init','Failed to get topology.',isub)
        GOTO 9999
      ENDIF
      !nsubs = topology%nsublist
      mesh  = topology%mesh(meshid)

      !-------------------------------------------------------------------------
      ! Setup mesh sizes for intermediate meshes/topologies
      !-------------------------------------------------------------------------
      IF (green .EQ. ppm_poisson_grn_pois_per) THEN
        !Copy size of global mesh
        ppmpoisson%nmxy (1) =  mesh%nm(1)
        ppmpoisson%nmxy (2) =  mesh%nm(2)
        ppmpoisson%nmxy (3) =  mesh%nm(3)
        !ppmpoisson%nmxyc(1) = (mesh%nm(1)-1)/2+1+1
        ppmpoisson%nmxyc(1) =  mesh%nm(1)
        ppmpoisson%nmxyc(2) =  mesh%nm(2)
        ppmpoisson%nmxyc(3) =  mesh%nm(3)
        ppmpoisson%nmz  (1) = (ppmpoisson%nmxyc(1))
        ppmpoisson%nmz  (2) = (ppmpoisson%nmxyc(2))
        ppmpoisson%nmz  (3) = (ppmpoisson%nmxyc(3))
        !Inverse of the size of the domain squared
        Lx2 = 1.0_MK/(topology%max_physd(1)-topology%min_physd(1))**2
        Ly2 = 1.0_MK/(topology%max_physd(2)-topology%min_physd(2))**2
        Lz2 = 1.0_MK/(topology%max_physd(3)-topology%min_physd(3))**2
      ELSE IF (green .EQ. ppm_poisson_grn_pois_fre) THEN !vertex
        !Copy size of global mesh
        ppmpoisson%nmxy (1) =  mesh%nm(1)*2
        ppmpoisson%nmxy (2) =  mesh%nm(2)*2
        ppmpoisson%nmxy (3) =  mesh%nm(3)*2
        !ppmpoisson%nmxyc(1) = (mesh%nm(1)-1)/2+1
        !ppmpoisson%nmxyc(1) = (mesh%nm(1)*2-1)/2+1
        ppmpoisson%nmxyc(1) =  mesh%nm(1)*2
        ppmpoisson%nmxyc(2) =  mesh%nm(2)*2
        ppmpoisson%nmxyc(3) =  mesh%nm(3)*2
        ppmpoisson%nmz  (1) = (ppmpoisson%nmxyc(1))
        ppmpoisson%nmz  (2) = (ppmpoisson%nmxyc(2))
        ppmpoisson%nmz  (3) = (ppmpoisson%nmxyc(3))
        !Determine the grid spacing !vertex
        dx = (topology%max_physd(1)-topology%min_physd(1))/(mesh%nm(1)-1)
        dy = (topology%max_physd(2)-topology%min_physd(2))/(mesh%nm(2)-1)
        dz = (topology%max_physd(3)-topology%min_physd(3))/(mesh%nm(3)-1)
      ENDIF

      !-------------------------------------------------------------------------
      ! Create temporary derivation arrays if necessary
      ! or spectral scaling coefficients
      !-------------------------------------------------------------------------
      IF (PRESENT(derive)) THEN
        IF ((      derive .EQ. ppm_poisson_drv_curl_fd2  &
          &  .OR.  derive .EQ. ppm_poisson_drv_curl_fd4)) THEN
          ppmpoisson%derivatives = derive
          ALLOCATE(ppmpoisson%drv_vr(LBOUND(fieldout,1):UBOUND(fieldout,1),&
                                   & LBOUND(fieldout,2):UBOUND(fieldout,2),&
                                   & LBOUND(fieldout,3):UBOUND(fieldout,3),&
                                   & LBOUND(fieldout,4):UBOUND(fieldout,4),&
                                   & LBOUND(fieldout,5):UBOUND(fieldout,5)))
        ELSE IF (derive .EQ. ppm_poisson_drv_curl_sp) THEN
          ppmpoisson%derivatives = ppm_poisson_drv_curl_sp
          IF (ppmpoisson%case .EQ. ppm_poisson_grn_pois_fre) THEN
            CALL ppm_write(ppm_rank,'ppm_poisson_init', &
            & 'WARNING: Spectral curl is not fully implemented in Freespace.',isub)
          ENDIF
        ELSE IF (derive .EQ. ppm_poisson_drv_none) THEN
          ppmpoisson%derivatives = ppm_poisson_drv_none
        ELSE
          CALL ppm_write(ppm_rank,'ppm_poisson_init','Undefined derivation input.',isub)
          GOTO 9999
        ENDIF
      ELSE
        ppmpoisson%derivatives = ppm_poisson_drv_none
      ENDIF
      !-------------------------------------------------------------------------
      ! Create spectral scaling components always. Just in case some 
      ! reprojection comes up
      ! The conditionals need to be for not just the Poisson equation
      !-------------------------------------------------------------------------
      IF (ppmpoisson%case .EQ. ppm_poisson_grn_pois_per) THEN
        ppmpoisson%normkx = &
        & 2.0_MK*PI/(topology%max_physd(1)-topology%min_physd(1))
        ppmpoisson%normky = &
        & 2.0_MK*PI/(topology%max_physd(2)-topology%min_physd(2))
        ppmpoisson%normkz = &
        & 2.0_MK*PI/(topology%max_physd(3)-topology%min_physd(3))
      ELSE IF (ppmpoisson%case .EQ. ppm_poisson_grn_pois_fre) THEN
        ppmpoisson%normkx = &
        & 2.0_MK*PI/((topology%max_physd(1)-topology%min_physd(1))*2.0_MK)
        ppmpoisson%normky = &
        & 2.0_MK*PI/((topology%max_physd(2)-topology%min_physd(2))*2.0_MK)
        ppmpoisson%normkz = &
        & 2.0_MK*PI/((topology%max_physd(3)-topology%min_physd(3))*2.0_MK)
      ENDIF



      !-------------------------------------------------------------------------
      ! Create new slab topology
      !-------------------------------------------------------------------------
      ttopoid = 0
      tmeshid = -1
      decomposition       = ppm_param_decomp_xy_slab
      assigning           = ppm_param_assign_internal
      IF (ppmpoisson%case .EQ. ppm_poisson_grn_pois_per) THEN
        bcdef               = ppm_param_bcdef_periodic
      ELSE IF (ppmpoisson%case .EQ. ppm_poisson_grn_pois_fre) THEN
        bcdef               = ppm_param_bcdef_freespace
      ENDIF
      tmpmin              = topology%min_physd
      tmpmax              = topology%max_physd

      CALL ppm_mktopo(ttopoid,tmeshid,xp,0,&
      & decomposition,assigning,&
      !& topology%min_physd,topology%max_physd,bcdef,&
      & tmpmin,tmpmax,bcdef,&
      & __ZEROSI,ppmpoisson%costxy,&
      & ppmpoisson%nmxy,info)
      CALL ppm_meshinfo(ttopoid,tmeshid,ppmpoisson%nmxy,ppmpoisson%istartxy,&
      &                 ppmpoisson%ndataxy,ppmpoisson%maxndataxy,&
      &                 isublist,nsublist,info)

      IF (info .NE. 0) THEN
        CALL ppm_write(ppm_rank,'ppm_poisson_init','Failed to create xy-topology.',isub)
        GOTO 9999
      ENDIF
      ppmpoisson%topoidxy = ttopoid
      ppmpoisson%meshidxy = tmeshid

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
      IF (ppmpoisson%case .EQ. ppm_poisson_grn_pois_per) THEN
        bcdef               = ppm_param_bcdef_periodic
      ELSE IF (ppmpoisson%case .EQ. ppm_poisson_grn_pois_fre) THEN
        bcdef               = ppm_param_bcdef_freespace
      ENDIF
      assigning       = ppm_param_assign_internal
      decomposition   = ppm_param_decomp_zpencil

      CALL ppm_mktopo(ttopoid,tmeshid,xp,0,&
      & decomposition,assigning,&
      & tmpmin,tmpmax,bcdef,&
      & __ZEROSI,ppmpoisson%costz,&
      & ppmpoisson%nmz,info)
      CALL ppm_meshinfo(ttopoid,tmeshid,ppmpoisson%nmz,ppmpoisson%istartz,&
      &                 ppmpoisson%ndataz,ppmpoisson%maxndataz,&
      &                 isublist,nsublist,info)
      IF (info .NE. 0) THEN
        CALL ppm_write(ppm_rank,'ppm_poisson_init','Failed to create z-topology.',isub)
        GOTO 9999
      ENDIF
      ppmpoisson%topoidz = ttopoid
      ppmpoisson%meshidz = tmeshid

      !-------------------------------------------------------------------------
      ! Get the new topologies (slabs and pencils)
      !-------------------------------------------------------------------------
      CALL ppm_topo_get(ppmpoisson%topoidxy,topologyxy,info)
      IF (info .NE. 0) THEN
        CALL ppm_write(ppm_rank,'ppm_poisson_init','Failed to get xy topology.',isub)
        GOTO 9999
      ENDIF

      !!CALL ppm_topo_get(ppmpoisson%topoidxyc,topologyxyc,info)
      !!IF (info .NE. 0) THEN
        !!CALL ppm_write(ppm_rank,'ppm_poisson_init','Failed to get complex xy topology.',isub)
        !!GOTO 9999
      !!ENDIF

      CALL ppm_topo_get(ppmpoisson%topoidz,topologyz,info)
      IF (info .NE. 0) THEN
        CALL ppm_write(ppm_rank,'ppm_poisson_init','Failed to get z topology.',isub)
        GOTO 9999
      ENDIF


      !-------------------------------------------------------------------------
      ! Copy data to ppm_poisson type
      !-------------------------------------------------------------------------
      ppmpoisson%nsublistxy  = topologyxy %nsublist
      !!ppmpoisson%nsublistxyc = topologyxyc%nsublist
      ppmpoisson%nsublistz   = topologyz  %nsublist

      ALLOCATE(ppmpoisson%isublistxy (ppmpoisson%nsublistxy))
      !!ALLOCATE(ppmpoisson%isublistxyc(ppmpoisson%nsublistxyc))
      ALLOCATE(ppmpoisson%isublistz  (ppmpoisson%nsublistz))

      DO isub=1,ppmpoisson%nsublistxy
        ppmpoisson%isublistxy(isub)  = topologyxy%isublist(isub)
      ENDDO
      !!DO isub=1,ppmpoisson%nsublistxyc
        !!ppmpoisson%isublistxyc(isub) = topologyxyc%isublist(isub)
      !!ENDDO
      DO isub=1,ppmpoisson%nsublistz
        ppmpoisson%isublistz(isub)   = topologyz%isublist(isub)
      ENDDO

      !-------------------------------------------------------------------------
      ! Set and get minimum and maximum indicies of xy slabs
      !-------------------------------------------------------------------------
      indl(1) = 1
      indl(2) = 1
      indl(3) = 1
      indu(1) = 0
      indu(2) = 0
      indu(3) = 0
      DO isub=1,ppmpoisson%nsublistxy
        isubl = ppmpoisson%isublistxy(isub)
        indu(1) = MAX(indu(1),ppmpoisson%ndataxy(1,isubl))
        indu(2) = MAX(indu(2),ppmpoisson%ndataxy(2,isubl))
        indu(3) = MAX(indu(3),ppmpoisson%ndataxy(3,isubl))
      ENDDO

      !-------------------------------------------------------------------------
      ! Allocate real xy slabs
      !-------------------------------------------------------------------------
      ALLOCATE(ppmpoisson%fldxyr(__DIM,&
      & indl(1):indu(1),indl(2):indu(2),indl(3):indu(3),&
      & 1:ppmpoisson%nsublistxy),stat=info)

      !!ALLOCATE(ppmpoisson%fldxyc(__DIM,&
      !!& indl(1):indu(1),indl(2):indu(2),indl(3):indu(3),&
      !!& 1:ppmpoisson%nsublistxy),stat=info)

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
      ALLOCATE(ppmpoisson%fldxyc(__DIM,&
      & indl(1):indu(1),indl(2):indu(2),indl(3):indu(3),&
      & 1:ppmpoisson%nsublistxy),stat=info)

      !-------------------------------------------------------------------------
      ! Set and get minimum and maximum indicies of z pencils
      !-------------------------------------------------------------------------
      indl(1) = 1
      indl(2) = 1
      indl(3) = 1
      indu(1) = 0
      indu(2) = 0
      indu(3) = 0
      DO isub=1,ppmpoisson%nsublistz
        isubl = ppmpoisson%isublistz(isub)
        indu(1) = MAX(indu(1),ppmpoisson%ndataz(1,isubl))
        indu(2) = MAX(indu(2),ppmpoisson%ndataz(2,isubl))
        indu(3) = MAX(indu(3),ppmpoisson%ndataz(3,isubl))
      ENDDO

      !-------------------------------------------------------------------------
      ! Allocate two complex z pencils + Greens fcn array !@check return vars.
      !-------------------------------------------------------------------------
      ALLOCATE(ppmpoisson%fldzc1(__DIM,&
      & indl(1):indu(1),indl(2):indu(2),indl(3):indu(3),&
      & 1:ppmpoisson%nsublistz),stat=info)

      ALLOCATE(ppmpoisson%fldzc2(__DIM,&
      & indl(1):indu(1),indl(2):indu(2),indl(3):indu(3),&
      & 1:ppmpoisson%nsublistz),stat=info)

      !-------------------------------------------------------------------------
      ! The complex Greens function is always kept on the z-pencil topology
      !-------------------------------------------------------------------------
      IF (green .EQ. ppm_poisson_grn_pois_per) THEN
        ALLOCATE(ppmpoisson%fldgrnr(&
        & indl(1):indu(1),indl(2):indu(2),indl(3):indu(3),&
        & 1:ppmpoisson%nsublistz),stat=info)
      ELSE IF (green .EQ. ppm_poisson_grn_pois_fre) THEN
        ALLOCATE(ppmpoisson%fldgrnc(&
        & indl(1):indu(1),indl(2):indu(2),indl(3):indu(3),&
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
      ! (d2_/dx2 + d2_/dy2 + d2_/dz2)psi = -omega     =>
      ! -4*pi2(kx2 + ky2 + kz2)PSI       = -OMEGA     =>
      ! PSI                              = 1/(4*pi2)*1/(kx2 + ky2 + kz2)OMEGA
      !-------------------------------------------------------------------------
      IF (green .EQ. ppm_poisson_grn_pois_per) THEN
        ! Scaling the spectral coefficients... 
        ! one minus due to (i*k)^2 and another due to the Poisson equation
        normfac = 1.0_MK/(4.0_MK*PI*PI * &
                !and normalisation of FFTs (full domain) !vertex
                & REAL((mesh%nm(1)-1)* &
                &      (mesh%nm(2)-1)* &
                &      (mesh%nm(3)-1),MK))
        DO isub=1,ppmpoisson%nsublistz
          isubl=ppmpoisson%isublistz(isub)
          DO k=1,ppmpoisson%ndataz(3,isubl)
            DO j=1,ppmpoisson%ndataz(2,isubl)
              DO i=1,ppmpoisson%ndataz(1,isubl)
                kx = i-1 + (ppmpoisson%istartz(1,isubl)-1)
                ky = j-1 + (ppmpoisson%istartz(2,isubl)-1)
                kz = k-1 + (ppmpoisson%istartz(3,isubl)-1)
                !The integer division depends on mesh%nm to include N+1 points:
                !This is a nasty way to do this but it is only done once so...:
                !(we subtract nm '-1' as not all points are included (periodic))
                IF (kx .GT. (ppmpoisson%nmz(1)-1)/2) kx = kx-(ppmpoisson%nmz(1)-1)
                IF (ky .GT. (ppmpoisson%nmz(2)-1)/2) ky = ky-(ppmpoisson%nmz(2)-1)
                IF (kz .GT. (ppmpoisson%nmz(3)-1)/2) kz = kz-(ppmpoisson%nmz(3)-1)
                ppmpoisson%fldgrnr(i,j,k,isub) = &
                      & normfac/(REAL(kx*kx,__PREC)*Lx2 &
                      &        + REAL(ky*ky,__PREC)*Ly2 &
                      &        + REAL(kz*kz,__PREC)*Lz2)
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
      ELSE IF (green .EQ. ppm_poisson_grn_pois_fre) THEN
        !-----------------------------------------------------------------------
        ! First initialise the real Greens function
        !@alternatively this could come from as input
        !-----------------------------------------------------------------------
        !there should NOT be a minus here since THIS Greens function takes
        !the minus of the Poisson equation into account
        normfac = 1.0_MK/(4.0_MK*PI* &
        !remembering FFT normalization of ALL points: !vertex
        & REAL((ppmpoisson%nmxy(1))*(ppmpoisson%nmxy(2))*(ppmpoisson%nmxy(3)),MK))*dx*dy*dz
        DO isub=1,ppmpoisson%nsublistxy
          isubl=ppmpoisson%isublistxy(isub)
          DO k=1,ppmpoisson%ndataxy(3,isubl)
            DO j=1,ppmpoisson%ndataxy(2,isubl)
              DO i=1,ppmpoisson%ndataxy(1,isubl)
                kx = i-1 + (ppmpoisson%istartxy(1,isubl)-1)
                ky = j-1 + (ppmpoisson%istartxy(2,isubl)-1)
                kz = k-1 + (ppmpoisson%istartxy(3,isubl)-1)
                !The integer division depends on mesh%nm to include N+1 points:
                !This is a nasty way to do this but it is only done once so...:
                !(we subtract all nmxy as all points are included in fresspace)
                IF (kx .GT. (ppmpoisson%nmxy(1)-2)/2) kx = kx-(ppmpoisson%nmxy(1))
                IF (ky .GT. (ppmpoisson%nmxy(2)-2)/2) ky = ky-(ppmpoisson%nmxy(2))
                IF (kz .GT. (ppmpoisson%nmxy(3)-2)/2) kz = kz-(ppmpoisson%nmxy(3))
                ppmpoisson%fldxyr(1,i,j,k,isub) = &
                        & normfac/(SQRT( REAL(kx*kx,__PREC)*dx*dx+ &
                        &                REAL(ky*ky,__PREC)*dy*dy+ &
                        &                REAL(kz*kz,__PREC)*dz*dz))
                ppmpoisson%fldxyr(2,i,j,k,isub) = 0.0_MK
                ppmpoisson%fldxyr(3,i,j,k,isub) = 0.0_MK
                !Take care of singularity
                !This is nasty as well
                IF ((kx*kx+ky*ky+kz*kz) .EQ. 0) THEN
                  !Simply zero
                  !ppmpoisson%fldxyr(1,i,j,k,isub) = 0.0_MK
                  !Simply one (H&E style)
                  !ppmpoisson%fldxyr(1,i,j,k,isub) = 1.0_MK*normfac
                  !Equal to the neighbour point (first order)
                  !ppmpoisson%fldxyr(1,i,j,k,isub) = normfac/(dx)
                  !Linear extrapolation into zero (2nd order)
                  ppmpoisson%fldxyr(1,i,j,k,isub) =   2.0_MK*normfac/(       dx) &
                                                  & - 1.0_MK*normfac/(2.0_MK*dx)
                  !extrapolation 3rd order
                  !ppmpoisson%fldxyr(1,i,j,k,isub) =   3.0_MK*normfac/(       dx) &
                                                  !& - 3.0_MK*normfac/(2.0_MK*dx) &
                                                  !& + 1.0_MK*normfac/(3.0_MK*dx)
                  !extrapolation 4th order
                  !ppmpoisson%fldxyr(1,i,j,k,isub) =   4.0_MK*normfac/(       dx) &
                                                  !& - 6.0_MK*normfac/(2.0_MK*dx) &
                                                  !& + 4.0_MK*normfac/(3.0_MK*dx) &
                                                  !& - 1.0_MK*normfac/(4.0_MK*dx)
                  !extrapolation 6th order
                  !ppmpoisson%fldxyr(1,i,j,k,isub) =   6.0_MK*normfac/(       dx) &
                                                  !& -15.0_MK*normfac/(2.0_MK*dx) &
                                                  !& +20.0_MK*normfac/(3.0_MK*dx) &
                                                  !& -15.0_MK*normfac/(4.0_MK*dx) &
                                                  !& + 6.0_MK*normfac/(5.0_MK*dx) &
                                                  !& - 1.0_MK*normfac/(6.0_MK*dx)
                  !one thousand!!!
                  !ppmpoisson%fldxyr(1,i,j,k,isub) = normfac*1000.0_MK
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
        & ppmpoisson%meshidxyc,ppmpoisson%fldxyc,__NCOM,info)
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
        & __NCOM,__ZEROSI,info)
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

      END SUBROUTINE __ROUTINE

