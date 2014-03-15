      !-------------------------------------------------------------------------
      !  Subroutine   : ppm_module_fft
      !-------------------------------------------------------------------------
      ! Copyright (c) 2010 CSE Lab (ETH Zurich), MOSAIC Group (ETH Zurich),
      !                    Center for Fluid Dynamics (DTU)
      !
      ! This modules contains routines for creating FFTW plans and executing them
      ! Plan routines are required for 
      !  z   means fft in x (1D)
      !  xy  means fft in xy (2D)
      !  xyz means fft in xyz (3D)
      !  sca means scalar
      !  vec means vector
      !  f   means forward fft
      !  b   means backward fft
      !  c   means complex
      !  r   means real
      !  s   means single precision
      !  d   means double precision
      !
      !  e.g. 3d_vec_fr2c_xy_d: 3D vector array forward transform from
      !                         real to complex and ffts in x and y directions
      !                         (2d). For double data
      !
      !  The following variants have been implemented
      !   3d_vec_fc2c_z_d
      !   3d_vec_bc2c_z_d
      !   3d_vec_fr2c_xy_d
      !   3d_vec_bc2r_xy_d
      !
      !  a normalization routine exists for debugging purposes
      !
      !  sofar all fftw calls are to double precision routines!
      !
      !  All IFFT routines should be called with the 'realest' topoid/meshid
      !
      !  The routines respects the periodic N+1 points periodic BC cf topo%bcdef
      !  but also does full domain FFTs for freespace BC
      !  It does not work on mixed periodic/freespace BC in the XY direction
      !  topo%bcdef(1) is assumed to be in x, (2) to be in y, (3) z, (4) x...
      !-------------------------------------------------------------------------
#define __SINGLE 0
#define __DOUBLE 1

      MODULE ppm_module_fft
      !-------------------------------------------------------------------------
      ! 1D transforms
      !-------------------------------------------------------------------------
      INTERFACE ppm_fft_forward_1d
         MODULE PROCEDURE ppm_fft_plan_3d_vec_fc2c_z_s
         MODULE PROCEDURE ppm_fft_plan_3d_vec_fc2c_z_d
      END INTERFACE

      INTERFACE ppm_fft_backward_1d
         MODULE PROCEDURE ppm_fft_plan_3d_vec_bc2c_z_s
         MODULE PROCEDURE ppm_fft_plan_3d_vec_bc2c_z_d
      END INTERFACE

      INTERFACE ppm_fft_execute_1d
         MODULE PROCEDURE ppm_fft_exec_3d_vec_c2c_z_s
         MODULE PROCEDURE ppm_fft_exec_3d_vec_c2c_z_d
      END INTERFACE

      !-------------------------------------------------------------------------
      ! 2D transforms
      !-------------------------------------------------------------------------
      INTERFACE ppm_fft_forward_2d
         MODULE PROCEDURE ppm_fft_plan_3d_vec_fr2c_xy_s
         MODULE PROCEDURE ppm_fft_plan_3d_vec_fr2c_xy_d
      END INTERFACE

      INTERFACE ppm_fft_backward_2d
         MODULE PROCEDURE ppm_fft_plan_3d_vec_bc2r_xy_s
         MODULE PROCEDURE ppm_fft_plan_3d_vec_bc2r_xy_d
      END INTERFACE

      INTERFACE ppm_fft_execute_2d
         MODULE PROCEDURE ppm_fft_exec_3d_vec_fr2c_xy_s
         MODULE PROCEDURE ppm_fft_exec_3d_vec_bc2r_xy_s
         MODULE PROCEDURE ppm_fft_exec_3d_vec_fr2c_xy_d
         MODULE PROCEDURE ppm_fft_exec_3d_vec_bc2r_xy_d
      END INTERFACE

      !-------------------------------------------------------------------------
      ! Normalization
      !-------------------------------------------------------------------------
      INTERFACE ppm_fft_normalize
         MODULE PROCEDURE ppm_fft_normalize_rs
         MODULE PROCEDURE ppm_fft_normalize_cs
         MODULE PROCEDURE ppm_fft_normalize_rd
         MODULE PROCEDURE ppm_fft_normalize_cd
      END INTERFACE

      !-------------------------------------------------------------------------
      ! PPM FFT plan type
      !-------------------------------------------------------------------------
      !!! Type containing the FFTW plan and its settings
      TYPE ppm_fft_plan
         !!!array of plan pointers, index for subs
         INTEGER*8,DIMENSION(:),POINTER  :: plan => NULL()
         !!!the dimension of the FFT (1D/2D/3D)
         INTEGER                         :: rank
         !!!number of points along each direction of the piece to be transformed
         !!!index is for rank and subs
         INTEGER,DIMENSION(:,:),POINTER  :: nx
         !!!the direction of the transform (forward/backward)
         INTEGER                         :: sign
         !!!the method to setup the optimal plan
         INTEGER                         :: flag
         !!!the number of components to transform
         INTEGER                         :: howmany
         !!!the size of the input array, index is for rank
         INTEGER,DIMENSION(:),POINTER    :: inembed
         !!!the size of the output array, index is for rank
         INTEGER,DIMENSION(:),POINTER    :: onembed
         !!!istride tells how same component data points are spaced in memory
         INTEGER                         :: istride
         INTEGER                         :: ostride
         !!!idist tells how multiple arrays are spaced. I.e. a memory offset
         INTEGER                         :: idist
         INTEGER                         :: odist
      END TYPE ppm_fft_plan

      CONTAINS
#define __KIND __SINGLE

      !FORWARD TRANSFORMS - PLAN
#include "fft/ppm_fft_plan_3d_vec_fc2c_z.f"
#include "fft/ppm_fft_plan_3d_vec_fr2c_xy.f"
      !BACKWARD TRANSFORMS - PLAN
#include "fft/ppm_fft_plan_3d_vec_bc2c_z.f"
#include "fft/ppm_fft_plan_3d_vec_bc2r_xy.f"
      !EXECUTION OF PLANS
#include "fft/ppm_fft_exec_3d_vec_c2c_z.f"
#include "fft/ppm_fft_exec_3d_vec_fr2c_xy.f"
#include "fft/ppm_fft_exec_3d_vec_bc2r_xy.f"
      !NORMALIZATION
#include "fft/ppm_fft_normalize_r.f"
#include "fft/ppm_fft_normalize_c.f"

#undef __KIND
#define __KIND __DOUBLE

      !FORWARD TRANSFORMS - PLAN
#include "fft/ppm_fft_plan_3d_vec_fc2c_z.f"
#include "fft/ppm_fft_plan_3d_vec_fr2c_xy.f"
      !BACKWARD TRANSFORMS - PLAN
#include "fft/ppm_fft_plan_3d_vec_bc2c_z.f"
#include "fft/ppm_fft_plan_3d_vec_bc2r_xy.f"
      !EXECUTION OF PLANS
#include "fft/ppm_fft_exec_3d_vec_c2c_z.f"
#include "fft/ppm_fft_exec_3d_vec_fr2c_xy.f"
#include "fft/ppm_fft_exec_3d_vec_bc2r_xy.f"
      !NORMALIZATION
#include "fft/ppm_fft_normalize_r.f"
#include "fft/ppm_fft_normalize_c.f"

#undef __KIND
      END MODULE ppm_module_fft


