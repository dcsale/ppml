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
      !
      !  Changelog: 
      !  5/16/2014 - Danny Sale - University of Washington - dsale@uw.edu - sale.danny@gmail.com
      !              * Modify include directories for Poisson routines
      !
      !-------------------------------------------------------------------------
#define __SINGLE 0
#define __DOUBLE 1

#define __DIM 3


#define __SIGMA  1.0_MK

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
         INTEGER                                              :: derivatives
         INTEGER                                              :: operation
         INTEGER                                              :: bc
         INTEGER                                              :: green

         REAL(ppm_kind_double)                                :: normkx
         REAL(ppm_kind_double)                                :: normky
         REAL(ppm_kind_double)                                :: normkz

         INTEGER, DIMENSION(__DIM)                            :: nmfft
         !!!Size of the FFT

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
         INTEGER,DIMENSION(:),POINTER                         :: isublistxy=>NULL()
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
         INTEGER,DIMENSION(:),POINTER                         :: isublistxyc=>NULL()
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
         INTEGER,DIMENSION(:),POINTER                         :: isublistz=>NULL()
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
         INTEGER,DIMENSION(:),POINTER                         :: isublistfr=>NULL()
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
      END TYPE ppm_poisson_plan


      INTEGER,PARAMETER                          :: ppm_poisson_opr_none     = 1
      INTEGER,PARAMETER                          :: ppm_poisson_opr_vel      = 2
      INTEGER,PARAMETER                          :: ppm_poisson_opr_repr     = 3

      INTEGER,PARAMETER                          :: ppm_poisson_grn_pois_fre    = 11
      INTEGER,PARAMETER                          :: ppm_poisson_grn_pois_per    = 12
      INTEGER,PARAMETER                          :: ppm_poisson_grn_pois_blob2  = 13
      INTEGER,PARAMETER                          :: ppm_poisson_grn_pois_blob4  = 14
      INTEGER,PARAMETER                          :: ppm_poisson_grn_pois_blob6  = 15
      INTEGER,PARAMETER                          :: ppm_poisson_grn_pois_blob8  = 16
      INTEGER,PARAMETER                          :: ppm_poisson_grn_pois_blob10 = 17

      INTEGER,PARAMETER                          :: ppm_poisson_drv_sp          = 21
      INTEGER,PARAMETER                          :: ppm_poisson_drv_fd2         = 22
      INTEGER,PARAMETER                          :: ppm_poisson_drv_fd4         = 23

      INTEGER,PARAMETER                          :: ppm_poisson_bc_per          = 31
      INTEGER,PARAMETER                          :: ppm_poisson_bc_fre          = 32

      !internal flags:
      INTEGER,PARAMETER                          :: ppm_poisson_curl            = 91
      INTEGER,PARAMETER                          :: ppm_poisson_divergence      = 92
      INTEGER,PARAMETER                          :: ppm_poisson_subtract        = 93


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
#define __KIND __SINGLE

#define __PREC ppm_kind_double
#define __DIM  3
#define __ZEROSI (/0,0,0/)
#define __NCOM  3
#define __ROUTINE ppm_poisson_init
#include "ppm_poisson_init.f"
#undef __ROUTINE
#define __ROUTINE ppm_poisson_solve
#include "ppm_poisson_solve.f"
#undef __ROUTINE
#define __ROUTINE ppm_poisson_finalize
#include "ppm_poisson_finalize.f"
#undef __ROUTINE
#define __ROUTINE ppm_poisson_fd
#include "ppm_poisson_fd.f"
#undef __ROUTINE
#define __ROUTINE ppm_poisson_extrapolateghost_vr
#include "ppm_poisson_extrapolateghost.f"
#undef __ROUTINE
#undef __ZEROSI
#undef __DIM
#undef __NCOM

#undef __KIND
#define __KIND __DOUBLE


#undef __KIND
      END MODULE ppm_module_poisson


