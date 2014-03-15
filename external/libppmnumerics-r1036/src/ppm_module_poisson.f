      !-------------------------------------------------------------------------
      !  Subroutine   : ppm_module_poisson
      !-------------------------------------------------------------------------
      ! Copyright (c) 2010 CSE Lab (ETH Zurich), MOSAIC Group (ETH Zurich),
      !                    Center for Fluid Dynamics (DTU)
      !
      ! Notes on the init routine:
      !   Arguments:
      !     topoid
      !     meshid
      !     ppmpoisson
      !     fieldin
      !     fieldout
      !     greens function: custom array, custom function, integer to build-in
      !     info - return variable
      !     periodic/freespace flag
      !     real/fourier greens function (optional, only for custom array)
      !     derivatives (optional): none,curl,gradient,laplace, - spectral/FD
      !     possibly the option to de- and reallocate arrays of slabs/pencils
      !
      ! The fields returned from the subroutine have been ghosted/extrapolated
      ! as necessary.
      !
      ! Add finalize routine
      ! This routine needs to be modified to allow cell centred data
      ! These routines should all be renamed - potentially not just Poisson 
      ! equations can be solved.
      !-------------------------------------------------------------------------
#define __SINGLE 0
#define __DOUBLE 1

#define __DIM 3

      MODULE ppm_module_poisson

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
         INTEGER                                              :: derivatives
         INTEGER                                              :: case

         REAL(ppm_kind_double)                                :: normkx
         REAL(ppm_kind_double)                                :: normky
         REAL(ppm_kind_double)                                :: normkz


         INTEGER                                              :: topoidxy
         !!!Topology id of xy topology
         INTEGER                                              :: meshidxy
         !!!mesh id for xy mesh
         REAL(ppm_kind_double), DIMENSION(:),  POINTER        :: costxy=>NULL()
         !!!sub cost for xy topology
         INTEGER , DIMENSION(:,:),POINTER                     :: istartxy=>NULL()
         !!!sub index start for xy topology
         INTEGER, DIMENSION(__DIM)                            :: nmxy
         !!!global number of grid points for xy topology
         INTEGER , DIMENSION(:,:),POINTER                     :: ndataxy=>NULL()
         !!!sub no. grid cells for xy topology
         INTEGER , DIMENSION(2)                               :: maxndataxy
         INTEGER,DIMENSION(:),POINTER                         :: isublistxy
         !!!list of sub domains on current CPU on the xy topology
         INTEGER                                              :: nsublistxy
         !!!number of sub domains for current CPU on the xy topology
         REAL(ppm_kind_double),DIMENSION(:,:,:,:,:),POINTER   :: fldxyr=>NULL()
         !!!real slab field (xy topology)
         COMPLEX(ppm_kind_double),DIMENSION(:,:,:,:,:),POINTER:: fldxyc=>NULL()
         !!!complex slab field (xy topology)

         INTEGER                                              :: topoidxyc
         !!!Topology id of complex xy topology
         INTEGER                                              :: meshidxyc
         !!!mesh id for comple xy mesh
         REAL(ppm_kind_double), DIMENSION(:),  POINTER        :: costxyc=>NULL()
         !!!sub cost for xy topology
         INTEGER , DIMENSION(:,:),POINTER                     :: istartxyc=>NULL()
         !!!sub index start for complex xy topology
         INTEGER, DIMENSION(__DIM)                            :: nmxyc
         !!!global number of grid points for complex xy topology
         INTEGER , DIMENSION(:,:),POINTER                     :: ndataxyc=>NULL()
         !!!sub no. grid cells for complex xy topology
         INTEGER,DIMENSION(:),POINTER                         :: isublistxyc
         !!!list of sub domains on current CPU on the complex xy topology
         INTEGER                                              :: nsublistxyc
         !!!number of sub domains for current CPU on the complex xy topology
         TYPE(ppm_fft_plan)                                   :: planfxy
         !!!fft plans for r2c complex xy slab
         TYPE(ppm_fft_plan)                                   :: planbxy
         !!!fft plans for c2r complex xy slab

         TYPE(ppm_fft_plan)                                   :: planfz
         !!!fft plans for forward c2c z pencil
         TYPE(ppm_fft_plan)                                   :: planbz
         !!!fft plans for backward c2c z pencil
         INTEGER                                              :: topoidz
         !!!Topology id of z topology
         INTEGER                                              :: meshidz
         !!!mesh id for z mesh
         REAL(ppm_kind_double), DIMENSION(:),  POINTER        :: costz =>NULL()
         !!!sub cost for z topology
         INTEGER , DIMENSION(:,:),POINTER                     :: istartz=>NULL()
         !!!sub index start for z topology
         INTEGER, DIMENSION(__DIM)                            :: nmz
         !!!global number of grid points for z topology
         INTEGER , DIMENSION(:,:),POINTER                     :: ndataz=>NULL()
         !!!sub no. grid cells for z topology
         INTEGER , DIMENSION(3)                               :: maxndataz
         INTEGER,DIMENSION(:),POINTER                         :: isublistz
         !!!list of sub domains on current CPU on the z topology
         INTEGER                                              :: nsublistz
         !!!number of sub domains for current CPU on the z topology
         COMPLEX(ppm_kind_double),DIMENSION(:,:,:,:,:),POINTER:: fldzc1=>NULL()
         !!!complex pencil field 1 (z topology)
         COMPLEX(ppm_kind_double),DIMENSION(:,:,:,:,:),POINTER:: fldzc2=>NULL()
         !!!complex pencil field 2 (z topology)

         INTEGER                                              :: topoidfr
         !!!Topology id of freespace topology
         INTEGER                                              :: meshidfr
         !!!mesh id for freespace mesh
         REAL(ppm_kind_double), DIMENSION(:),  POINTER        :: costfr=>NULL()
         !!!sub cost
         INTEGER , DIMENSION(:,:),POINTER                     :: istartfr=>NULL()
         !!!sub index start
         INTEGER, DIMENSION(__DIM)                            :: nmfr
         !!!global number of grid points
         INTEGER , DIMENSION(:,:),POINTER                     :: ndatafr=>NULL()
         !!!sub no. grid cells
         INTEGER,DIMENSION(:),POINTER                         :: isublistfr
         !!!list of sub domains on current CPU on the freespace topology
         INTEGER                                              :: nsublistfr
         !!!number of sub domains for current CPU on the freespace topology

         !These field are only allocated as necessary:
         REAL(ppm_kind_double),DIMENSION(:,:,:,:),POINTER     :: fldfrs=>NULL()
         !!!dummy array for the right hand side, for free space, scalar fields
         REAL(ppm_kind_double),DIMENSION(:,:,:,:,:),POINTER   :: fldfrv=>NULL()
         !!!dummy array for the right hand side, for free space, vector fields
         REAL(ppm_kind_double),DIMENSION(:,:,:,:),POINTER     :: fldgrnr=>NULL()
         !!!real Greens field, z-pencils, scalar
         COMPLEX(ppm_kind_double),DIMENSION(:,:,:,:),POINTER  :: fldgrnc=>NULL()
         !!!complex Greens field, z-pencils, scalar
         REAL(ppm_kind_double),DIMENSION(:,:,:,:,:),POINTER   :: drv_vr=>NULL()
         !!!dummy array for the right hand side, for free space, vector fields

         ! COMPLEX(ppm_kind_double),DIMENSION(:,:,:,:,:),POINTER:: fldgrnc3=>NULL()
         ! !!!complex Greens 3 component vector field intended for spectral
         ! !!!derivatives, z-pencils, scalar
      END TYPE ppm_poisson_plan

      INTEGER,PARAMETER                          :: ppm_poisson_drv_none    =0
      INTEGER,PARAMETER                          :: ppm_poisson_drv_curl_sp =1
      INTEGER,PARAMETER                          :: ppm_poisson_drv_grad_sp =2
      INTEGER,PARAMETER                          :: ppm_poisson_drv_lapl_sp =3
      INTEGER,PARAMETER                          :: ppm_poisson_drv_div_sp  =4
      INTEGER,PARAMETER                          :: ppm_poisson_drv_curl_fd2=11
      INTEGER,PARAMETER                          :: ppm_poisson_drv_grad_fd2=12
      INTEGER,PARAMETER                          :: ppm_poisson_drv_lapl_fd2=13
      INTEGER,PARAMETER                          :: ppm_poisson_drv_div_fd2 =14
      INTEGER,PARAMETER                          :: ppm_poisson_drv_curl_fd4=21
      INTEGER,PARAMETER                          :: ppm_poisson_drv_grad_fd4=22
      INTEGER,PARAMETER                          :: ppm_poisson_drv_lapl_fd4=23
      INTEGER,PARAMETER                          :: ppm_poisson_drv_div_fd4 =24

      INTEGER,PARAMETER                          :: ppm_poisson_grn_pois_per =1
      INTEGER,PARAMETER                          :: ppm_poisson_grn_pois_fre =2
      INTEGER,PARAMETER                          :: ppm_poisson_grn_reprojec =3

      INTERFACE ppm_poisson_init
         MODULE PROCEDURE ppm_poisson_init
      END INTERFACE

      INTERFACE ppm_poisson_solve
         MODULE PROCEDURE ppm_poisson_solve
      END INTERFACE

      INTERFACE ppm_poisson_fd
         MODULE PROCEDURE ppm_poisson_fd
      END INTERFACE

      INTERFACE ppm_poisson_extrapolateghost
         MODULE PROCEDURE ppm_poisson_extrapolateghost_vr
      END INTERFACE

      CONTAINS
#define __KIND __SINGLE

#define __PREC ppm_kind_double
        !!#define __CMPLXDEF DCMPLX
#define __DIM  3
#define __ZEROSI (/0,0,0/)
#define __NCOM  3
#define __ROUTINE ppm_poisson_init
#include "poisson/ppm_poisson_init.f"
#undef __ROUTINE
#define __ROUTINE ppm_poisson_solve
#include "poisson/ppm_poisson_solve.f"
#undef __ROUTINE
#define __ROUTINE ppm_poisson_fd
#include "poisson/ppm_poisson_fd.f"
#undef __ROUTINE
#define __ROUTINE ppm_poisson_extrapolateghost_vr
#include "poisson/ppm_poisson_extrapolateghost.f"
#undef __ROUTINE
#undef __ZEROSI
#undef __DIM
#undef __NCOM

#undef __KIND
#define __KIND __DOUBLE


#undef __KIND
      END MODULE ppm_module_poisson


