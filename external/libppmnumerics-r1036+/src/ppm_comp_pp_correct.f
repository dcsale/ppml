      !-------------------------------------------------------------------------
      !  Subroutine   :               ppm_comp_pp_correct
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
      SUBROUTINE ppm_comp_pp_correct_si(dh,cutoff,kernel,eps,kpar,targt,  &
     &    kappa,info)
#elif __KIND == __DOUBLE_PRECISION
      SUBROUTINE ppm_comp_pp_correct_di(dh,cutoff,kernel,eps,kpar,targt,  &
     &    kappa,info)
#elif __KIND == __SINGLE_PRECISION_COMPLEX
      SUBROUTINE ppm_comp_pp_correct_sci(dh,cutoff,kernel,eps,kpar,targt,  &
     &    kappa,info)
#elif __KIND == __DOUBLE_PRECISION_COMPLEX
      SUBROUTINE ppm_comp_pp_correct_dci(dh,cutoff,kernel,eps,kpar,targt,  &
     &    kappa,info)
#endif
#elif __KERNEL == __USER_FUNCTION
#if   __KIND == __SINGLE_PRECISION
      SUBROUTINE ppm_comp_pp_correct_su(dh,cutoff,kernel,eps,kpar,targt,  &
     &    kappa,info)
#elif __KIND == __DOUBLE_PRECISION
      SUBROUTINE ppm_comp_pp_correct_du(dh,cutoff,kernel,eps,kpar,targt,  &
     &    kappa,info)
#elif __KIND == __SINGLE_PRECISION_COMPLEX
      SUBROUTINE ppm_comp_pp_correct_scu(dh,cutoff,kernel,eps,kpar,targt,  &
     &    kappa,info)
#elif __KIND == __DOUBLE_PRECISION_COMPLEX
      SUBROUTINE ppm_comp_pp_correct_dcu(dh,cutoff,kernel,eps,kpar,targt,  &
     &    kappa,info)
#endif
#elif __KERNEL == __LOOKUP_TABLE
#if   __KIND == __SINGLE_PRECISION
      SUBROUTINE ppm_comp_pp_correct_st(dh,cutoff,kernel,eps,kpar,targt,  &
     &    kappa,info)
#elif __KIND == __DOUBLE_PRECISION
      SUBROUTINE ppm_comp_pp_correct_dt(dh,cutoff,kernel,eps,kpar,targt,  &
     &    kappa,info)
#elif __KIND == __SINGLE_PRECISION_COMPLEX
      SUBROUTINE ppm_comp_pp_correct_sct(dh,cutoff,kernel,eps,kpar,targt,  &
     &    kappa,info)
#elif __KIND == __DOUBLE_PRECISION_COMPLEX
      SUBROUTINE ppm_comp_pp_correct_dct(dh,cutoff,kernel,eps,kpar,targt,  &
     &    kappa,info)
#endif
#endif
      !!! Computes the discretisation correction factor
      !!! needed to ensure that the discrete second-order
      !!! moment of the kernel function is equal to some
      !!! value. The kernel is evaluated on a regular mesh
      !!! here, so this factor can only be used if the
      !!! particles are on a regular grid (e.g. after
      !!! remeshing).
      !-------------------------------------------------------------------------
      !  Modules
      !-------------------------------------------------------------------------
      USE ppm_module_data
      USE ppm_module_substart
      USE ppm_module_substop
      USE ppm_module_error
      USE ppm_module_write
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
      REAL(MK), DIMENSION(:), INTENT(IN   ) :: dh
      !!! grid spacings in all directions 1..ppm_dim
      REAL(MK)              , INTENT(IN   ) :: cutoff
      !!! PP interaction cutoff used. If no cutoff is used, pass a negative
      !!! value here.
      REAL(MK)              , INTENT(IN   ) :: eps
      !!! Kernel size (width)
#if   __KIND == __SINGLE_PRECISION | __KIND == __DOUBLE_PRECISION
      REAL(MK)              , INTENT(IN   ) :: targt
      !!! target value for the second order kernel moment. Overloaded types:
      !!! single,double
      REAL(MK)              , INTENT(  OUT) :: kappa
      !!! computed discretisation correction factor. Every kernel evaluation has
      !!! to be multiplied with this in order to correct for discretisation
      !!! effects. Overloaded types: single,double
#elif __KIND == __SINGLE_PRECISION_COMPLEX | __KIND == __DOUBLE_PRECISION_COMPLEX
      COMPLEX(MK)           , INTENT(IN   ) :: targt
      !!! target value for the second order kernel moment. Overloaded types:
      !!! single complex, double complex
      COMPLEX(MK)           , INTENT(  OUT) :: kappa
      !!! computed discretisation correction factor. Every kernel evaluation has
      !!! to be multiplied with this in order to correct for discretisation
      !!! effects. Overloaded types: single complex, double complex
#endif
#if   __KERNEL == __INTERNAL
      INTEGER               , INTENT(IN   ) :: kernel
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
      REAL(MK)   , DIMENSION(:), INTENT(IN   ) :: kpar
      !!! Kernel parameters. See documentation or ppm_comp_pp_kernels.inc for
      !!! description. Type can be single or complex.
#elif __KIND == __SINGLE_PRECISION_COMPLEX | __KIND == __DOUBLE_PRECISION_COMPLEX
      COMPLEX(MK), DIMENSION(:), INTENT(IN   ) :: kpar
      !!! Kernel parameters. See documentation or ppm_comp_pp_kernels.inc for
      !!! description. Type can be single complex or double complex.
#endif
#elif __KERNEL == __LOOKUP_TABLE
      REAL(MK)                   , INTENT(IN   ) :: kpar
      !!! Kernel parameters. See documentation or ppm_comp_pp_kernels.inc for
      !!! description. Pass dxtableinv (the inverse of the table spacing) as
      !!! a scalar here.
#if   __KIND == __SINGLE_PRECISION | __KIND == __DOUBLE_PRECISION
      REAL(MK)   , DIMENSION(:  ), INTENT(IN   ) :: kernel
#elif __KIND == __SINGLE_PRECISION_COMPLEX | __KIND == __DOUBLE_PRECISION_COMPLEX
      COMPLEX(MK), DIMENSION(:  ), INTENT(IN   ) :: kernel
#endif
      !!! Lookup table with tabulated kernel values.
      !!! Such a table can be created using <<ppm_comp_pp_mk_table>>.
#endif
      INTEGER               , INTENT(  OUT) :: info
      !!! Returns 0 upon success
      !-------------------------------------------------------------------------
      !  Local variables
      !-------------------------------------------------------------------------
      REAL(MK)                 :: t0,dij2,dij4,dij5
      REAL(MK)                 :: hi,dx,dy,dz,dij,factor,factor2
#if   __KIND == __SINGLE_PRECISION | __KIND == __DOUBLE_PRECISION
      REAL(MK)                 :: eta,summ
#elif __KIND == __SINGLE_PRECISION_COMPLEX | __KIND == __DOUBLE_PRECISION_COMPLEX
      COMPLEX(MK)              :: eta,summ
#endif
      INTEGER                  :: lhx,lhy,lhz,i,j,k,idx
      CHARACTER(LEN=ppm_char)  :: cbuf
      !-------------------------------------------------------------------------
      !  Externals
      !-------------------------------------------------------------------------
#if   __KERNEL == __USER_FUNCTION
#if   __KIND == __SINGLE_PRECISION | __KIND == __DOUBLE_PRECISION
      REAL(MK)   , DIMENSION(:), INTENT(IN   ) :: kpar
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
      COMPLEX(MK), DIMENSION(:), INTENT(IN   ) :: kpar
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
      CALL substart('ppm_comp_pp_correct',t0,info)
      kappa = 1.0_MK

      !-------------------------------------------------------------------------
      !  Check arguments.
      !-------------------------------------------------------------------------
      IF (ppm_debug .GT. 0) THEN
          IF (.NOT. ppm_initialized) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_ppm_noinit,'ppm_comp_pp_correct',  &
     &            'Please call ppm_init first!',__LINE__,info)
              GOTO 9999
          ENDIF
          IF (eps .LT. 0.0_MK) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_comp_pp_correct',  &
     &            'eps must be > 0',__LINE__,info)
              GOTO 9999
          ENDIF
          IF (SIZE(dh,1) .LT. ppm_dim) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_comp_pp_correct',  &
     &            'dh must be at least of length ppm_dim',__LINE__,info)
              GOTO 9999
          ENDIF
          DO i=1,ppm_dim
              IF (dh(i) .LE. 0.0_MK) THEN
                  info = ppm_error_error
                  CALL ppm_error(ppm_err_argument,'ppm_comp_pp_correct',  &
     &                'dh must be >0 in all directions',__LINE__,info)
                  GOTO 9999
              ENDIF
          ENDDO
#if   __KERNEL == __LOOKUP_TABLE
          IF (kpar .LT. 0.0_MK) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_comp_pp_correct',  &
     &            'kpar (dxtableinv) must be >=0',__LINE__,info)
              GOTO 9999
          ENDIF
#endif
      ENDIF

      !-------------------------------------------------------------------------
      !  Integration bounds
      !-------------------------------------------------------------------------
      IF (cutoff .GE. 0.0_MK) THEN
          hi  = cutoff
      ELSE
          hi  = 50.0_MK
      ENDIF

      !-------------------------------------------------------------------------
      !  Summation bounds
      !-------------------------------------------------------------------------
      lhx = INT(hi/dh(1))
      lhy = INT(hi/dh(2))
      IF (ppm_dim .GT. 2) lhz = INT(hi/dh(3))

      !-------------------------------------------------------------------------
      !  Discrete sum
      !-------------------------------------------------------------------------
      summ = 0.0_MK
      IF (ppm_dim .GT. 2) THEN
          DO i=-lhx,lhx
              dx = REAL(i,MK)*dh(1)
              DO j=-lhy,lhy
                  dy = REAL(j,MK)*dh(2)
                  DO k=-lhz,lhz
                      dz = REAL(k,MK)*dh(3)
                      dij  = dx*dx + dy*dy + dz*dz
#include "ppm_comp_pp_kernels.inc"
                      ! second order moment
                      eta = dij*eta
                      summ = summ + eta
                  ENDDO
              ENDDO
          ENDDO
          summ = summ*dh(1)*dh(2)*dh(3)
      ELSE
          DO i=-lhx,lhx
              dx = REAL(i,MK)*dh(1)
              DO j=-lhy,lhy
                  dy = REAL(j,MK)*dh(2)
                  dij  = dx*dx + dy*dy
#include "ppm_comp_pp_kernels.inc"
                  ! second order moment
                  eta = dij*eta
                  summ = summ + eta
              ENDDO
          ENDDO
          summ = summ*dh(1)*dh(2)
      ENDIF

      !-------------------------------------------------------------------------
      !  Correction factor
      !-------------------------------------------------------------------------
#if __KERNEL != __LOOKUP_TABLE
      kappa = (REAL(ppm_dim,MK)*targt*eps)/summ
#endif

      !-------------------------------------------------------------------------
      !  Debug output
      !-------------------------------------------------------------------------
      IF (ppm_debug .GT. 1) THEN
          WRITE(cbuf,'(A,F16.10)') 'Discretization correction factor: ',kappa
          CALL ppm_write(ppm_rank,'ppm_comp_pp_correct',cbuf,info)
      ENDIF

      !-------------------------------------------------------------------------
      !  Return
      !-------------------------------------------------------------------------
 9999 CONTINUE
      CALL substop('ppm_comp_pp_correct',t0,info)
      RETURN
#if   __KERNEL == __INTERNAL
#if   __KIND == __SINGLE_PRECISION
      END SUBROUTINE ppm_comp_pp_correct_si
#elif __KIND == __DOUBLE_PRECISION
      END SUBROUTINE ppm_comp_pp_correct_di
#elif __KIND == __SINGLE_PRECISION_COMPLEX
      END SUBROUTINE ppm_comp_pp_correct_sci
#elif __KIND == __DOUBLE_PRECISION_COMPLEX
      END SUBROUTINE ppm_comp_pp_correct_dci
#endif

#elif __KERNEL == __USER_FUNCTION
#if   __KIND == __SINGLE_PRECISION
      END SUBROUTINE ppm_comp_pp_correct_su
#elif __KIND == __DOUBLE_PRECISION
      END SUBROUTINE ppm_comp_pp_correct_du
#elif __KIND == __SINGLE_PRECISION_COMPLEX
      END SUBROUTINE ppm_comp_pp_correct_scu
#elif __KIND == __DOUBLE_PRECISION_COMPLEX
      END SUBROUTINE ppm_comp_pp_correct_dcu
#endif

#elif __KERNEL == __LOOKUP_TABLE
#if   __KIND == __SINGLE_PRECISION
      END SUBROUTINE ppm_comp_pp_correct_st
#elif __KIND == __DOUBLE_PRECISION
      END SUBROUTINE ppm_comp_pp_correct_dt
#elif __KIND == __SINGLE_PRECISION_COMPLEX
      END SUBROUTINE ppm_comp_pp_correct_sct
#elif __KIND == __DOUBLE_PRECISION_COMPLEX
      END SUBROUTINE ppm_comp_pp_correct_dct
#endif
#endif
