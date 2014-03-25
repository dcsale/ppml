      !-------------------------------------------------------------------------
      !  Subroutine   :                 ppm_util_gmres_solveupper
      !-------------------------------------------------------------------------
      !
      !  Purpose      : Solve Ux=b
      !
      !  Input        : U          (F) upper triangular matrix
      !                 b          (F) Right-hand side
      !                 n          (I) size
      !
      !  Output       : info       (I) return status.
      !                 x          (F) Solution
      !
      !  Remarks      :
      !
      !  References   : 
      !
      !
      !  Revisions    :
      !-------------------------------------------------------------------------
      !  $Log: ppm_util_gmres_solveupper.f,v $
      !  Revision 1.1.1.1  2007/07/13 10:19:01  ivos
      !  CBL version of the PPM library
      !
      !  Revision 1.1  2006/05/11 10:27:00  pchatela
      !  Initial insertion
      !
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
      SUBROUTINE ppm_util_gmres_solveupper_s(U,b,x,n,info)
#elif __KIND == __DOUBLE_PRECISION
      SUBROUTINE ppm_util_gmres_solveupper_d(U,b,x,n,info)
#endif
      !-------------------------------------------------------------------------
      !  Includes
      !-------------------------------------------------------------------------
#include "ppm_define.h"

      !-------------------------------------------------------------------------
      !  Modules
      !-------------------------------------------------------------------------
      USE ppm_module_data
      USE ppm_module_substart
      USE ppm_module_substop
      USE ppm_module_error
      USE ppm_module_alloc
      IMPLICIT NONE

#if   __KIND == __SINGLE_PRECISION
      INTEGER, PARAMETER :: MK = ppm_kind_single
#elif __KIND == __DOUBLE_PRECISION
      INTEGER, PARAMETER :: MK = ppm_kind_double
#endif

      !-------------------------------------------------------------------------
      !  Arguments
      !-------------------------------------------------------------------------
      REAL(MK), DIMENSION(:,:), INTENT(IN)  :: U
      REAL(MK), DIMENSION(:),   POINTER     :: x
      REAL(MK), DIMENSION(:),   INTENT(IN)  :: b
      INTEGER               , INTENT(IN)    :: n

      INTEGER               , INTENT(OUT)   :: info
      !-------------------------------------------------------------------------
      !  Local variables
      !-------------------------------------------------------------------------
      REAL(MK)                              :: t0, lmyeps
      REAL(KIND(1.0D0))                     :: sum
      INTEGER                               :: i,j

      !-----------------------------------------------------------------------
      !  call substart
      !-----------------------------------------------------------------------
      CALL substart('ppm_util_gmres_solveupper',t0,info)

      !-----------------------------------------------------------------------
      !  check input arguments
      !-----------------------------------------------------------------------
      IF (ppm_debug .GT. 0) THEN
         IF (n.LT.0) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_util_gmres_solveupper',  &
     &            'passed a negative size n',__LINE__,info)
              GOTO 9999
          ENDIF
      ENDIF

#if   __KIND == __SINGLE_PRECISION
      lmyeps = ppm_myepss
#elif __KIND == __DOUBLE_PRECISION
      lmyeps = ppm_myepsd
#endif
      
      !-----------------------------------------------------------------------
      !  Elimination
      !-----------------------------------------------------------------------
      DO i=n,1,-1
         sum = b(i)
         DO j = i+1,n
            sum = sum - U(i,j)*x(j)
         END DO
         IF (ABS(U(i,i)).GT.lmyeps) THEN
            x(i) = sum/U(i,i)
         ELSE
            info = ppm_error_error
            CALL ppm_error(ppm_err_argument,'ppm_util_gmres_solveupper',  &
     &            'a pivot of U is null',__LINE__,info)
            GOTO 9999
         END IF
      END DO
      
      !-------------------------------------------------------------------------
      !  Return
      !-------------------------------------------------------------------------
 9999 CONTINUE
      CALL substop('ppm_util_gmres_solveupper',t0,info)
      RETURN

#if   __KIND == __SINGLE_PRECISION
END SUBROUTINE ppm_util_gmres_solveupper_s
#elif __KIND == __DOUBLE_PRECISION
END SUBROUTINE ppm_util_gmres_solveupper_d
#endif
