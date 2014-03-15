      !-------------------------------------------------------------------------
      !  Subroutine   :               ppm_comp_pp_mk_table
      !-------------------------------------------------------------------------
      ! Copyright (c) 2010 CSE Lab (ETH Zurich), MOSAIC Group (ETH Zurich), 
      !                    Center for Fluid Dynamics (DTU)
      !
      !
      ! This file is part of the Parallel Particle Mesh Library (PPM).
      !
      ! PPM is free software: you can redistribute it and/or modify
      ! it under the terms of the GNU Lesser General Public License 
      ! as published by the Free Software Foundation, either 
      ! version 3 of the License, or (at your option) any later 
      ! version.
      !
      ! PPM is distributed in the hope that it will be useful,
      ! but WITHOUT ANY WARRANTY; without even the implied warranty of
      ! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
      ! GNU General Public License for more details.
      !
      ! You should have received a copy of the GNU General Public License
      ! and the GNU Lesser General Public License along with PPM. If not,
      ! see <http://www.gnu.org/licenses/>.
      !
      ! Parallel Particle Mesh Library (PPM)
      ! ETH Zurich
      ! CH-8092 Zurich, Switzerland
      !-------------------------------------------------------------------------

#if   __KERNEL == __INTERNAL
#if   __KIND == __SINGLE_PRECISION
      SUBROUTINE ppm_comp_pp_mk_table_si(kernel,kpar,ntable,cutoff2,   &
     &    ktab,dxtableinv,info)
#elif __KIND == __DOUBLE_PRECISION
      SUBROUTINE ppm_comp_pp_mk_table_di(kernel,kpar,ntable,cutoff2,   &
     &    ktab,dxtableinv,info)
#elif __KIND == __SINGLE_PRECISION_COMPLEX
      SUBROUTINE ppm_comp_pp_mk_table_sci(kernel,kpar,ntable,cutoff2,   &
     &    ktab,dxtableinv,info)
#elif __KIND == __DOUBLE_PRECISION_COMPLEX
      SUBROUTINE ppm_comp_pp_mk_table_dci(kernel,kpar,ntable,cutoff2,   &
     &    ktab,dxtableinv,info)
#endif
#elif __KERNEL == __USER_FUNCTION
#if   __KIND == __SINGLE_PRECISION
      SUBROUTINE ppm_comp_pp_mk_table_su(kernel,kpar,ntable,cutoff2,   &
     &    ktab,dxtableinv,info)
#elif __KIND == __DOUBLE_PRECISION
      SUBROUTINE ppm_comp_pp_mk_table_du(kernel,kpar,ntable,cutoff2,   &
     &    ktab,dxtableinv,info)
#elif __KIND == __SINGLE_PRECISION_COMPLEX
      SUBROUTINE ppm_comp_pp_mk_table_scu(kernel,kpar,ntable,cutoff2,   &
     &    ktab,dxtableinv,info)
#elif __KIND == __DOUBLE_PRECISION_COMPLEX
      SUBROUTINE ppm_comp_pp_mk_table_dcu(kernel,kpar,ntable,cutoff2,   &
     &    ktab,dxtableinv,info)
#endif
#endif
      !!! Creates a lookup table for a PP interaction kernel.
      !!! The result is returned in the variables ktab with
      !!! indices running from 0 to ntable-1. The table
      !!! contains kernel values eta as a function of the
      !!! _squared_ distance x**2.
      !-------------------------------------------------------------------------
      !  Modules
      !-------------------------------------------------------------------------
      USE ppm_module_data
      USE ppm_module_substart
      USE ppm_module_substop
      USE ppm_module_error
      USE ppm_module_alloc
      IMPLICIT NONE
#if   __KIND == __SINGLE_PRECISION | __KIND == __SINGLE_PRECISION_COMPLEX
      INTEGER, PARAMETER :: MK = ppm_kind_single
#elif __KIND == __DOUBLE_PRECISION | __KIND == __DOUBLE_PRECISION_COMPLEX
      INTEGER, PARAMETER :: MK = ppm_kind_double
#endif
      !-------------------------------------------------------------------------
      !  Includes
      !-------------------------------------------------------------------------
#include "ppm_define.h"
      !-------------------------------------------------------------------------
      !  Arguments
      !-------------------------------------------------------------------------
      INTEGER                    , INTENT(IN   ) :: ntable
      !!! Number of table entries to be created.
      REAL(MK)                   , INTENT(IN   ) :: cutoff2
      !!! Square of the cutoff used for PP interactions. The kernel will be
      !!! tabulated up to this value.
#if   __KIND == __SINGLE_PRECISION | __KIND == __DOUBLE_PRECISION
      REAL(MK)     , DIMENSION(:), POINTER       :: ktab
      !!! Tabulated kernel values 0...ntable-1.
      !!! Overloaded types: single,double.
#elif __KIND == __SINGLE_PRECISION_COMPLEX | __KIND == __DOUBLE_PRECISION_COMPLEX
      COMPLEX(MK)  , DIMENSION(:), POINTER       :: ktab
      !!! Tabulated kernel values 0...ntable-1.
      !!! Overloaded types: single complex,double complex.
#endif
#if   __KERNEL == __INTERNAL
      INTEGER                    , INTENT(IN   ) :: kernel
      !!! kernel for which to compute the correction. To use ppm-internal
      !!! kernels, specify one of:
      !!!
      !!! ---------------------------------------
      !!!    ppm_param_kerel_laplace2d_2p
      !!!      (2nd order Laplacian,
      !!!      polynomial in 2D)
      !!!    ppm_param_kerel_laplace3d_2p
      !!!      (2nd order Laplacian,
      !!!      polynomial in 3D)
      !!! ---------------------------------------
      !!!
      !!! To use your own kernel function, pass the function pointer here. Your
      !!! function should take one argument and return one value.
#if   __KIND == __SINGLE_PRECISION | __KIND == __DOUBLE_PRECISION
      REAL(MK)   , DIMENSION(:)  , INTENT(IN   ) :: kpar
#elif __KIND == __SINGLE_PRECISION_COMPLEX | __KIND == __DOUBLE_PRECISION_COMPLEX
      COMPLEX(MK), DIMENSION(:)  , INTENT(IN   ) :: kpar
#endif
      !!! Kernel parameters. See documentation or ppm_comp_pp_kernels.inc for
      !!! description. Type can be single, double, single complex or double
      !!! complex. Omit this argument when using a lookup table.
#endif
      REAL(MK)                   , INTENT(  OUT) :: dxtableinv
      !!! inverse of the dx (kernel evaluation locations) used to create
      !!! the table. To use the table, look up the value
      !!! `ktab(INT(dxtableinv*x))`.
      INTEGER                    , INTENT(  OUT) :: info
      !!! Returns status, 0 upon success
      !-------------------------------------------------------------------------
      !  Local variables
      !-------------------------------------------------------------------------
      INTEGER                  :: i,atable,iopt,idx
      INTEGER, DIMENSION(1)    :: ldl,ldu
      REAL(MK), PARAMETER      :: Rmin = 0.0_MK
      REAL(MK)                 :: dxtable,Rmax,factor,factor2
      REAL(MK)                 :: t0,dij,dij2,dij4,dij5,dx,dy,dz
#if   __KIND == __SINGLE_PRECISION | __KIND == __DOUBLE_PRECISION
      REAL(MK)                 :: eta
#elif __KIND == __SINGLE_PRECISION_COMPLEX | __KIND == __DOUBLE_PRECISION_COMPLEX
      COMPLEX(MK)              :: eta
#endif
      !-------------------------------------------------------------------------
      !  Externals
      !-------------------------------------------------------------------------
#if   __KERNEL == __USER_FUNCTION
#if   __KIND == __SINGLE_PRECISION | __KIND == __DOUBLE_PRECISION
      REAL(MK)   , DIMENSION(:)  , INTENT(IN   ) :: kpar
      INTERFACE
          FUNCTION kernel(x,kpar)
              IMPLICIT NONE
#if   __KIND == __SINGLE_PRECISION | __KIND == __SINGLE_PRECISION_COMPLEX
              INTEGER, PARAMETER :: MK = KIND(1.0E0)
#elif __KIND == __DOUBLE_PRECISION | __KIND == __DOUBLE_PRECISION_COMPLEX
              INTEGER, PARAMETER :: MK = KIND(1.0D0)
#endif
              REAL(MK), INTENT(IN) :: x
              REAL(MK), DIMENSION(:), INTENT(IN) :: kpar
              REAL(MK) :: kernel
          END FUNCTION kernel
      END INTERFACE
#elif __KIND == __SINGLE_PRECISION_COMPLEX | __KIND == __DOUBLE_PRECISION_COMPLEX
      COMPLEX(MK), DIMENSION(:)  , INTENT(IN   ) :: kpar
      INTERFACE
          FUNCTION kernel(x,kpar)
              IMPLICIT NONE
#if   __KIND == __SINGLE_PRECISION | __KIND == __SINGLE_PRECISION_COMPLEX
              INTEGER, PARAMETER :: MK = KIND(1.0E0)
#elif __KIND == __DOUBLE_PRECISION | __KIND == __DOUBLE_PRECISION_COMPLEX
              INTEGER, PARAMETER :: MK = KIND(1.0D0)
#endif
              REAL(MK), INTENT(IN) :: x
              COMPLEX(MK), DIMENSION(:), INTENT(IN) :: kpar
              COMPLEX(MK) :: kernel
          END FUNCTION kernel
      END INTERFACE
#endif
#endif

      !-------------------------------------------------------------------------
      !  Initialise
      !-------------------------------------------------------------------------
      CALL substart('ppm_comp_pp_mk_table',t0,info)

      !-------------------------------------------------------------------------
      !  Check arguments.
      !-------------------------------------------------------------------------
      IF (ppm_debug .GT. 0) THEN
          IF (ntable .LT. 1) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_comp_pp_mk_table',  &
     &            'ntable must be >0',__LINE__,info)
              GOTO 9999
          ENDIF
          IF (cutoff2 .LT. 0.0_MK) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_comp_pp_mk_table',  &
     &            'cutoff2 must be >=0',__LINE__,info)
              GOTO 9999
          ENDIF
      ENDIF

      !-------------------------------------------------------------------------
      !  Compute Rmax
      !-------------------------------------------------------------------------
      Rmax = cutoff2

      !-------------------------------------------------------------------------
      !  Allocate memory for table lookup
      !-------------------------------------------------------------------------
      atable = ntable-1

      iopt = ppm_param_alloc_fit
      ldl(1) = 0
      ldu(1) = atable
      CALL ppm_alloc(ktab,ldl,ldu,iopt,info)
      IF (info .NE. ppm_param_success) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_comp_pp_mk_table',     &
     &        'kernel lookup table KTAB',__LINE__,info)
          GOTO 9999
      ENDIF

      !-------------------------------------------------------------------------
      !  Compute values
      !-------------------------------------------------------------------------
      dxtable    = (Rmax - Rmin)/REAL(atable,MK)
      dxtableinv = 1.0_MK/dxtable

      ! This hack is needed since some kernels use dx,dy,dz explicitly. DO
      ! NOT MAKE TABLES FOR THOSE KERNELS!
      dx = 0.0_MK
      dy = 0.0_MK
      dz = 0.0_MK
      DO i=0,atable
         dij = Rmin + REAL(i,MK)*dxtable
#include "ppm_comp_pp_kernels.inc"
         ktab(i) = eta
      ENDDO

      !-------------------------------------------------------------------------
      !  Return
      !-------------------------------------------------------------------------
 9999 CONTINUE
      CALL substop('ppm_comp_pp_mk_table',t0,info)
      RETURN
#if   __KERNEL == __INTERNAL
#if   __KIND == __SINGLE_PRECISION
      END SUBROUTINE ppm_comp_pp_mk_table_si
#elif __KIND == __DOUBLE_PRECISION
      END SUBROUTINE ppm_comp_pp_mk_table_di
#elif __KIND == __SINGLE_PRECISION_COMPLEX
      END SUBROUTINE ppm_comp_pp_mk_table_sci
#elif __KIND == __DOUBLE_PRECISION_COMPLEX
      END SUBROUTINE ppm_comp_pp_mk_table_dci
#endif

#elif __KERNEL == __USER_FUNCTION
#if   __KIND == __SINGLE_PRECISION
      END SUBROUTINE ppm_comp_pp_mk_table_su
#elif __KIND == __DOUBLE_PRECISION
      END SUBROUTINE ppm_comp_pp_mk_table_du
#elif __KIND == __SINGLE_PRECISION_COMPLEX
      END SUBROUTINE ppm_comp_pp_mk_table_scu
#elif __KIND == __DOUBLE_PRECISION_COMPLEX
      END SUBROUTINE ppm_comp_pp_mk_table_dcu
#endif
#endif
