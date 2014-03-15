      !-------------------------------------------------------------------------
      !  Subroutine   :                ppm_gmm_reinitialize
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
#if    __DIM == __2D
#if    __KIND == __SINGLE_PRECISION
      SUBROUTINE ppm_gmm_reinitialize_2ds(fdata,tol,width,   &
     &    order,info,thresh,chi,MaxIter)
#elif  __KIND == __DOUBLE_PRECISION
      SUBROUTINE ppm_gmm_reinitialize_2dd(fdata,tol,width,   &
     &    order,info,thresh,chi,MaxIter)
#endif 
#elif  __DIM == __3D
#if    __KIND == __SINGLE_PRECISION
      SUBROUTINE ppm_gmm_reinitialize_3ds(fdata,tol,width,   &
     &    order,info,thresh,chi,MaxIter)
#elif  __KIND == __DOUBLE_PRECISION
      SUBROUTINE ppm_gmm_reinitialize_3dd(fdata,tol,width,   &
     &    order,info,thresh,chi,MaxIter)
#endif 
#endif
      !!! This routine re-initializes a signed distance level function with the
      !!! zero level representing the interface. Shift accordingly if other
      !!! level is to be used. The group marching method is used.
      !!!
      !!! === References ===
      !!!
      !!! S. Kim. An O(N) level set method for Eikonal equations. SIAM J. Sci.
      !!! Comput. 22(6):2178-2193, 2001.
      !-------------------------------------------------------------------------
      !  Includes
      !-------------------------------------------------------------------------
#include "ppm_define.h"
      !-------------------------------------------------------------------------
      !  Modules 
      !-------------------------------------------------------------------------
      USE ppm_module_data
      
      USE ppm_module_data_gmm
      USE ppm_module_gmm_init
      USE ppm_module_gmm_kickoff
      USE ppm_module_gmm_march
      USE ppm_module_gmm_finalize
      USE ppm_module_substart
      USE ppm_module_substop
      USE ppm_module_error
      USE ppm_module_typedef
      IMPLICIT NONE
#if    __KIND == __SINGLE_PRECISION | __KIND == __SINGLE_PRECISION_COMPLEX
      INTEGER, PARAMETER :: MK = ppm_kind_single
#else
      INTEGER, PARAMETER :: MK = ppm_kind_double
#endif
      !-------------------------------------------------------------------------
      !  Arguments     
      !-------------------------------------------------------------------------
#if   __DIM == __2D
      REAL(MK), DIMENSION(:,:,:  )   , POINTER             :: fdata
      !!! Field data. Rank 3 (for 2D scalar fields).
      !!! Indices: (i,j,[k],isub). On input: old level function values. The
      !!! interface is at level zero. A ghostsize of 1 is needed on all sides
      !!! which must be filled with the old level function value on input!!
      !!! On output: reinitialized signed distance function using the
      !!! interpolation method of Chopp. Points far from the interface will have
      !!! the value HUGE.
      REAL(MK), DIMENSION(:,:,:,:)   , INTENT(IN), OPTIONAL:: chi
      !!! Rank 4 (2d) field specifying the positions of the grid nodes.
      !!! 1st index: 1..ppm_dim, then i,j,[k],isub.  OPTIONAL. Uniform grid is
      !!! assumed if absent. Ghostlayers of size >=1 must be pre-filled.
#elif __DIM == __3D
      REAL(MK), DIMENSION(:,:,:,:)   , POINTER             :: fdata
      !!! Field data. Rank 4 (for 3D scalar fields).
      !!! Indices: (i,j,[k],isub). On input: old level function values. The
      !!! interface is at level zero. A ghostsize of 1 is needed on all sides
      !!! which must be filled with the old level function value on input!!
      !!! On output: reinitialized signed distance function using the
      !!! interpolation method of Chopp. Points far from the interface will have
      !!! the value HUGE.
      REAL(MK), DIMENSION(:,:,:,:,:) , INTENT(IN), OPTIONAL:: chi
      !!! Rank 5 (3d) field specifying the positions
      !!! of the grid nodes. 1st index: 1..ppm_dim, then i,j,[k],isub.
      !!! OPTIONAL. Uniform grid is assumed if absent. Ghostlayers of size >=1
      !!! must be pre-filled.
#endif
      INTEGER                        , INTENT(IN   )       :: order
      !!! Order of the method to be used. One of
      !!!
      !!! *ppm_param_order_1
      !!! *ppm_param_order_2
      !!! *ppm_param_order_3
      REAL(MK)                       , INTENT(IN   )       :: tol
      !!! Relative tolerance for the determined distance to the interface.
      !!! 1E-3 is a good choice. The tolerance is in multiples of grid spacings.
      REAL(MK)                       , INTENT(IN   )       :: width
      !!! Width of the narrow band to be produced on each side of the interface.
      INTEGER                        , INTENT(  OUT)       :: info
      REAL(MK) , OPTIONAL            , INTENT(IN   )       :: thresh
      !!! OPTIONAL. Threshold for interface detection. If this is not specified,
      !!! it is set to MAXVAL(ABS(fdata)).
      INTEGER  , OPTIONAL            , INTENT(IN   )       :: MaxIter
      !!! OPTIONAL argument specifying the maximum number of allowed iterations.
      !!! This can be useful since a cyclic dependency in the GMM algorithms
      !!! could cause infinite loops. In each iteration at least one point is
      !!! computed.
      !-------------------------------------------------------------------------
      !  Local variables 
      !-------------------------------------------------------------------------
      INTEGER                     :: xhi,i,isub,Nminit,MaxIt
      REAL(MK)                    :: t0,th
      LOGICAL                     :: lok  
      TYPE(ppm_t_topo), POINTER   :: topo
      !-------------------------------------------------------------------------
      !  Initialise 
      !-------------------------------------------------------------------------
      CALL substart('ppm_gmm_reinitialize',t0,info)
      topo => ppm_topo(gmm_topoid)%t
      !-------------------------------------------------------------------------
      !  Check arguments
      !-------------------------------------------------------------------------
      IF (ppm_debug .GT. 0) THEN
          IF (.NOT. ppm_initialized) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_ppm_noinit,'ppm_gmm_reinitialize',  &
     &            'Please call ppm_init first!',__LINE__,info)
              GOTO 9999
          ENDIF
          IF (tol .LE. 0.0_MK) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_gmm_reinitialize',  &
     &            'tolerance must be >0!',__LINE__,info)
              GOTO 9999
          ENDIF
          IF (width .LE. 0.0_MK) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_gmm_reinitialize',  &
     &            'width must be >0!',__LINE__,info)
              GOTO 9999
          ENDIF
#if   __DIM == __3D
          IF (SIZE(fdata,4) .LT. topo%nsublist) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_gmm_reinitialize',  &
     &            'field data for some subs is missing',__LINE__,info)
              GOTO 9999
          ENDIF
          IF (UBOUND(fdata,1) .LT. maxxhi) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_gmm_reinitialize',  &
     &            'x dimension of field data does not match mesh',__LINE__,info)
              GOTO 9999
          ENDIF
          IF (UBOUND(fdata,2) .LT. maxyhi) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_gmm_reinitialize',  &
     &            'y dimension of field data does not match mesh',__LINE__,info)
              GOTO 9999
          ENDIF
          IF (UBOUND(fdata,3) .LT. maxzhi) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_gmm_reinitialize',  &
     &            'z dimension of field data does not match mesh',__LINE__,info)
              GOTO 9999
          ENDIF
#elif __DIM == __2D
          IF (SIZE(fdata,3) .LT. topo%nsublist) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_gmm_reinitialize',  &
     &            'field data for some subs is missing',__LINE__,info)
              GOTO 9999
          ENDIF
          IF (UBOUND(fdata,1) .LT. maxxhi) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_gmm_reinitialize',  &
     &            'x dimension of field data does not match mesh',__LINE__,info)
              GOTO 9999
          ENDIF
          IF (UBOUND(fdata,2) .LT. maxyhi) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_gmm_reinitialize',  &
     &            'y dimension of field data does not match mesh',__LINE__,info)
              GOTO 9999
          ENDIF
#endif
      ENDIF          ! ppm_debug for argument check
      !-------------------------------------------------------------------------
      !  Kickoff GMM. We assume that the largest occuring values are
      !  outside of the band.
      !-------------------------------------------------------------------------
      IF (PRESENT(thresh)) THEN
          th = thresh
      ELSE
          th = 0.99_MK*MAXVAL(ABS(fdata))
      ENDIF
      IF (PRESENT(chi)) THEN
          CALL ppm_gmm_kickoff(fdata,tol,th,info,chi=chi)
      ELSE
          CALL ppm_gmm_kickoff(fdata,tol,th,info)
      ENDIF
      IF (info .NE. ppm_param_success) THEN
          info = ppm_error_error
          CALL ppm_error(ppm_err_sub_failed,'ppm_gmm_reinitialize',  &
     &        'Starting GMM failed.',__LINE__,info)
          GOTO 9999
      ENDIF
      !-------------------------------------------------------------------------
      !  Check if the maximum number of iterations was specified
      !-------------------------------------------------------------------------
      IF (PRESENT(MaxIter)) THEN
          MaxIt = MaxIter
      ELSE
          MaxIt = HUGE(MaxIt)
      ENDIF
      !-------------------------------------------------------------------------
      !  Marching GMM
      !-------------------------------------------------------------------------
      IF (PRESENT(chi)) THEN
          CALL ppm_gmm_march(width,order,fdata,1.0_MK,MaxIt,info,chi=chi)
      ELSE
          CALL ppm_gmm_march(width,order,fdata,1.0_MK,MaxIt,info)
      ENDIF
      IF (info .NE. ppm_param_success) THEN
          info = ppm_error_error
          CALL ppm_error(ppm_err_sub_failed,'ppm_gmm_reinitialize',  &
     &        'Marching GMM failed.',__LINE__,info)
          GOTO 9999
      ENDIF
      !-------------------------------------------------------------------------
      !  Return 
      !-------------------------------------------------------------------------
 9999 CONTINUE
      CALL substop('ppm_gmm_reinitialize',t0,info)
      RETURN
#if    __DIM == __2D
#if    __KIND == __SINGLE_PRECISION
      END SUBROUTINE ppm_gmm_reinitialize_2ds
#elif  __KIND == __DOUBLE_PRECISION
      END SUBROUTINE ppm_gmm_reinitialize_2dd
#endif 
#elif  __DIM == __3D
#if    __KIND == __SINGLE_PRECISION
      END SUBROUTINE ppm_gmm_reinitialize_3ds
#elif  __KIND == __DOUBLE_PRECISION
      END SUBROUTINE ppm_gmm_reinitialize_3dd
#endif 
#endif
