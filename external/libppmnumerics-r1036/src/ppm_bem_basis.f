      !-------------------------------------------------------------------------
      !  Subroutine   :                  ppm_bem_basis
      !-------------------------------------------------------------------------
      !
      !  Purpose      : Returns the coefficients of the basis functions for
      !                 a given (s,t) pair in the triangle
      !
      !  Input        : qp      (F) : coordinates where the basis functions
      !                               has to be evaluated
      !                 basis   (I) : set of basis functions. One of the
      !                               following values:
      !                                ppm_bem_basis_const  :
      !                                  constant interpolation, 1 coefficient
      !                                ppm_bem_basis_linear :
      !                                  linear interpolation, 3 coefficients
      !                                ppm_bem_basis_quad   :
      !                                  quadratic interpolation, 6 coefficients
      !
      !  Input/output : lb      (F) : coefficients of the basis functions
      !
      !  Output       : info    (I) : return status. 0 on success.
      !
      !  Remarks      : lb has to be big enough to hold the coefficients
      !                 for a given set of basis functions.
      !
      !                  t
      !                  3 (0,1)
      !                  |\
      !                  | \
      !                  4  5
      !                  |   \
      !                  |    \
      !                  1--6--2 (1,0) s
      !                 (0,0)
      !
      !                 The order of the coefficients in the lb array is
      !                 mentioned in the figure above. E.g. for the linear
      !                 case, lb(1) contains the coefficient of the value
      !                 at 1, lb(2) contains the coefficient of the value
      !                 at 2, lb(3) contains the coefficient of the value
      !                 at 3.
      !
      !  References   :
      !
      !  Revisions    :
      !-------------------------------------------------------------------------
      !  $Log: ppm_bem_basis.f,v $
      !  Revision 1.1.1.1  2007/07/13 10:18:54  ivos
      !  CBL version of the PPM library
      !
      !  Revision 1.6  2006/09/04 18:34:40  pchatela
      !  Fixes and cleanups to make it compile for
      !  release-1-0
      !
      !  Revision 1.5  2004/10/01 16:08:55  ivos
      !  Replaced REAL(ppm_kind_double) :: t0 with REAL(MK) t0.
      !
      !  Revision 1.4  2004/07/27 07:06:48  ivos
      !  Renamed parameters after moving their definition to ppm_param.h.
      !
      !  Revision 1.3  2004/07/26 11:46:52  ivos
      !  Fixes to make it compile.
      !
      !  Revision 1.2  2004/07/26 07:46:35  ivos
      !  Changed to use single-interface modules. Updated all USE statements.
      !
      !  Revision 1.1  2004/07/16 08:32:46  oingo
      !  Initial release. Not tested
      !
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
#if   __KIND == __SINGLE_PRECISION
      SUBROUTINE ppm_bem_basis_s(qp,lb,basis,info)
#elif __KIND == __DOUBLE_PRECISION
      SUBROUTINE ppm_bem_basis_d(qp,lb,basis,info)
#endif
      !-------------------------------------------------------------------------
      !  Includes
      !-------------------------------------------------------------------------
#include "ppm_define.h"
      !-------------------------------------------------------------------------
      !  Modules
      !-------------------------------------------------------------------------
      USE ppm_module_data
      USE ppm_module_numerics_data
      USE ppm_module_error
      USE ppm_module_substart
      USE ppm_module_substop
      IMPLICIT NONE
      !-------------------------------------------------------------------------
      !  Type kind
      !-------------------------------------------------------------------------
#if   __KIND == __SINGLE_PRECISION
      INTEGER, PARAMETER :: MK = ppm_kind_single
#elif __KIND == __DOUBLE_PRECISION
      INTEGER, PARAMETER :: MK = ppm_kind_double
#endif
      !-------------------------------------------------------------------------
      !  Arguments     
      !-------------------------------------------------------------------------
      REAL(MK), DIMENSION(:), INTENT(IN   ) :: qp
      REAL(MK), DIMENSION(:), POINTER       :: lb
      INTEGER               , INTENT(IN   ) :: basis
      INTEGER               , INTENT(  OUT) :: info
      !-------------------------------------------------------------------------
      !  Local variables 
      !-------------------------------------------------------------------------
      REAL(MK)               :: u
      REAL(MK)               :: t0
      !-------------------------------------------------------------------------
      !  Externals 
      !-------------------------------------------------------------------------

      !-------------------------------------------------------------------------
      !  Initialise 
      !-------------------------------------------------------------------------
      CALL substart('ppm_bem_basis',t0,info)
      !-------------------------------------------------------------------------
      !  Check arguments
      !-------------------------------------------------------------------------
      IF (ppm_debug .GT. 0) THEN
         IF ((basis .LT. 1) .OR. (basis .GT. 3)) THEN
            info = ppm_error_fatal
            CALL ppm_error(ppm_err_argument,'ppm_bem_basis',      &
     &           'Invalid value for BASIS',__LINE__,info)
            GOTO 9999
         ENDIF
         IF ((SIZE(lb,1) .LT. 1) .AND. (basis .EQ.    &
     &        ppm_param_bem_basis_const)) THEN
            info = ppm_error_fatal
            CALL ppm_error(ppm_err_argument,'ppm_bem_basis',  &
     &           'LB must have a size of at least 1',__LINE__,info)
            GOTO 9999
         ENDIF
         IF ((SIZE(lb,1) .LT. 3) .AND. (basis .EQ.    &
     &        ppm_param_bem_basis_linear)) THEN
            info = ppm_error_fatal
            CALL ppm_error(ppm_err_argument,'ppm_bem_basis',  &
     &           'LB must have a size of at least 3',__LINE__,info)
            GOTO 9999
         ENDIF
         IF ((SIZE(lb,1) .LT. 6) .AND. (basis .EQ.    &
     &        ppm_param_bem_basis_quad)) THEN
            info = ppm_error_fatal
            CALL ppm_error(ppm_err_argument,'ppm_bem_basis',  &
     &           'LB must have a size of at least 6',__LINE__,info)
            GOTO 9999
         ENDIF
      ENDIF
      !-------------------------------------------------------------------------
      !  Calculate the coefficients
      !-------------------------------------------------------------------------
      IF (basis .EQ. ppm_param_bem_basis_const) THEN
         lb(1) = 1.0
      ELSEIF (basis .EQ. ppm_param_bem_basis_linear) THEN
         lb(1) = 1.0 - qp(1) - qp(2)
         lb(2) = qp(1)
         lb(3) = qp(2)
      ELSE
         u = 1.0 - qp(1) - qp(2)
         lb(1) = u * (2.0*u - 1.0)
         lb(2) = qp(1) * (2.0*qp(1) - 1.0)
         lb(3) = qp(2) * (2.0*qp(2) - 1.0)
         lb(4) = 4.0 * qp(2) * u
         lb(5) = 4.0 * qp(1) * qp(2)
         lb(6) = 4.0 * qp(1) * u
      ENDIF
      !-------------------------------------------------------------------------
      !  Return 
      !-------------------------------------------------------------------------
 9999 CONTINUE
      CALL substop('ppm_bem_basis',t0,info)
      RETURN

#if   __KIND == __SINGLE_PRECISION
      END SUBROUTINE ppm_bem_basis_s
#elif __KIND == __DOUBLE_PRECISION
      END SUBROUTINE ppm_bem_basis_d
#endif
