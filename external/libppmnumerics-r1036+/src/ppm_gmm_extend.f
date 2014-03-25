      !-------------------------------------------------------------------------
      !  Subroutine   :                  ppm_gmm_extend
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
#if    __KICKOFF == __YES
#if    __DIM == __2D
#if    __TYPE == __SFIELD
#if    __KIND == __SINGLE_PRECISION
      SUBROUTINE ppm_gmm_extend_2d_ksca_s(ivalue,fdata,udata,tol,  &
     &    width,order,info,chi,MaxIter)
#elif  __KIND == __DOUBLE_PRECISION
      SUBROUTINE ppm_gmm_extend_2d_ksca_d(ivalue,fdata,udata,tol,  &
     &    width,order,info,chi,MaxIter)
#endif 
#elif  __TYPE == __VFIELD
#if    __KIND == __SINGLE_PRECISION
      SUBROUTINE ppm_gmm_extend_2d_kvec_s(ivalue,fdata,udata,lda,tol, &
     &    width,order,info,chi,MaxIter)
#elif  __KIND == __DOUBLE_PRECISION
      SUBROUTINE ppm_gmm_extend_2d_kvec_d(ivalue,fdata,udata,lda,tol, &
     &    width,order,info,chi,MaxIter)
#endif 
#endif 
#elif  __DIM == __3D
#if    __TYPE == __SFIELD
#if    __KIND == __SINGLE_PRECISION
      SUBROUTINE ppm_gmm_extend_3d_ksca_s(ivalue,fdata,udata,tol,    &
     &    width,order,info,chi,MaxIter)
#elif  __KIND == __DOUBLE_PRECISION
      SUBROUTINE ppm_gmm_extend_3d_ksca_d(ivalue,fdata,udata,tol,    &
     &    width,order,info,chi,MaxIter)
#endif 
#elif  __TYPE == __VFIELD
#if    __KIND == __SINGLE_PRECISION
      SUBROUTINE ppm_gmm_extend_3d_kvec_s(ivalue,fdata,udata,lda,    &
     &    tol,width,order,info,chi,MaxIter)
#elif  __KIND == __DOUBLE_PRECISION
      SUBROUTINE ppm_gmm_extend_3d_kvec_d(ivalue,fdata,udata,lda,    &
     &    tol,width,order,info,chi,MaxIter)
#endif 
#endif 
#endif
#elif  __KICKOFF == __NO
#if    __DIM == __2D
#if    __TYPE == __SFIELD
#if    __KIND == __SINGLE_PRECISION
      SUBROUTINE ppm_gmm_extend_2d_tsca_s(ivalue,fdata,udata,tol,  &
     &    width,order,info,chi,MaxIter)
#elif  __KIND == __DOUBLE_PRECISION
      SUBROUTINE ppm_gmm_extend_2d_tsca_d(ivalue,fdata,udata,tol,  &
     &    width,order,info,chi,MaxIter)
#endif 
#elif  __TYPE == __VFIELD
#if    __KIND == __SINGLE_PRECISION
      SUBROUTINE ppm_gmm_extend_2d_tvec_s(ivalue,fdata,udata,lda,tol, &
     &    width,order,info,chi,MaxIter)
#elif  __KIND == __DOUBLE_PRECISION
      SUBROUTINE ppm_gmm_extend_2d_tvec_d(ivalue,fdata,udata,lda,tol, &
     &    width,order,info,chi,MaxIter)
#endif 
#endif 
#elif  __DIM == __3D
#if    __TYPE == __SFIELD
#if    __KIND == __SINGLE_PRECISION
      SUBROUTINE ppm_gmm_extend_3d_tsca_s(ivalue,fdata,udata,tol,    &
     &    width,order,info,chi,MaxIter)
#elif  __KIND == __DOUBLE_PRECISION
      SUBROUTINE ppm_gmm_extend_3d_tsca_d(ivalue,fdata,udata,tol,    &
     &    width,order,info,chi,MaxIter)
#endif 
#elif  __TYPE == __VFIELD
#if    __KIND == __SINGLE_PRECISION
      SUBROUTINE ppm_gmm_extend_3d_tvec_s(ivalue,fdata,udata,lda,    &
     &    tol,width,order,info,chi,MaxIter)
#elif  __KIND == __DOUBLE_PRECISION
      SUBROUTINE ppm_gmm_extend_3d_tvec_d(ivalue,fdata,udata,lda,    &
     &    tol,width,order,info,chi,MaxIter)
#endif 
#endif 
#endif
#endif
      !!! This routine extends a function defined on the interface to the whole
      !!! band on which the level function is defined. The extension is done
      !!! such that the gradient of the function is perpendicular to the
      !!! gradient of the level function. ppm_gmm_init must be called BEFORE
      !!! this routine is invoked.
      !-------------------------------------------------------------------------
      !  Includes
      !-------------------------------------------------------------------------
#include "ppm_define.h"
      !-------------------------------------------------------------------------
      !  Modules 
      !-------------------------------------------------------------------------
      USE ppm_module_data
      USE ppm_module_numerics_data
      USE ppm_module_data_gmm
      USE ppm_module_gmm_init
      USE ppm_module_gmm_cpt
      USE ppm_module_gmm_march
      USE ppm_module_gmm_finalize
      USE ppm_module_substart
      USE ppm_module_substop
      USE ppm_module_error
      USE ppm_module_alloc
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
#if   __KICKOFF == __YES
#if   __KIND == __SINGLE_PRECISION
      INTERFACE
          FUNCTION ivalue(x,y,info)
              REAL(KIND(1.0E0)), INTENT(IN)            :: x,y
              REAL(KIND(1.0E0))                        :: ivalue
              INTEGER, INTENT(OUT)                     :: info
          END FUNCTION ivalue
      END INTERFACE
#elif __KIND == __DOUBLE_PRECISION
      INTERFACE
          FUNCTION ivalue(x,y,info)
              REAL(KIND(1.0D0)), INTENT(IN)            :: x,y
              REAL(KIND(1.0D0))                        :: ivalue
              INTEGER, INTENT(OUT)                     :: info
          END FUNCTION ivalue
      END INTERFACE
#endif
#else
      REAL(MK)                       , INTENT(IN   )   :: ivalue
      !!! A (F) scalar value defining a cutoff. Points closer to the interface
      !!! than this cutoff will be kept to initialize the marching. The cpt is
      !!! skipped in this case.
#endif
#if   __TYPE == __SFIELD
      REAL(MK), DIMENSION(:,:,:)     , POINTER         :: udata
#elif __TYPE == __VFIELD
      REAL(MK), DIMENSION(:,:,:,:)   , POINTER         :: udata
      INTEGER                        , INTENT(IN   )   :: lda
      !!! Only present for vector udata. Gives the length of the leading
      !!! dimension of udata. All elements will be equally extended.
#endif
      REAL(MK), DIMENSION(:,:,:)     , POINTER         :: fdata
      !!! Level function data. Needs to be defined in a narrow band of
      !!! specific width around the interface. The zero level is interpreted
      !!! to be the interface. Always a scalar field. THIS NEEDS TO BE
      !!! PROPERLY ALLOCATED ON INPUT, INCLUDING GHOST LAYERS OF SIZE order.
      REAL(MK), DIMENSION(:,:,:,:)   , INTENT(IN),OPTIONAL:: chi
      !!! rank 4 (2d) field specifying the positions of the
      !!! grid nodes. 1st index: 1..ppm_dim, then i,j,[k],isub. OPTIONAL.
      !!! Uniform grid is assumed if absent. Ghostlayers of size >=1 must
      !!! be pre-filled.
#elif __DIM == __3D
#if   __KICKOFF == __YES
#if   __KIND == __SINGLE_PRECISION
      INTERFACE
          FUNCTION ivalue(x,y,z,info)
              REAL(KIND(1.0E0)), INTENT(IN)            :: x,y,z
              REAL(KIND(1.0E0))                        :: ivalue
              INTEGER, INTENT(OUT)                     :: info
          END FUNCTION ivalue
      END INTERFACE
#elif __KIND == __DOUBLE_PRECISION
      INTERFACE
          FUNCTION ivalue(x,y,z,info)
              REAL(KIND(1.0D0)), INTENT(IN)            :: x,y,z
              REAL(KIND(1.0D0))                        :: ivalue
              INTEGER, INTENT(OUT)                     :: info
          END FUNCTION ivalue
      END INTERFACE
#endif
      !!! Function pointer to the function computing the value of the function
      !!! on the interface: (F) ivalue(x,y[,z]) The function may assume that
      !!! the point (x,y[,z]) is on the interface.
#else
      REAL(MK)                       , INTENT(IN   )   :: ivalue
      !!! A (F) scalar value defining a cutoff. Points closer to the interface
      !!! than this cutoff will be kept to initialize the marching. The cpt is
      !!! skipped in this case.
#endif
#if   __TYPE == __SFIELD
      REAL(MK), DIMENSION(:,:,:,:)   , POINTER         :: udata
      !!! Field of the function to be extended, defined in the same narrow
      !!! band as the level function. Values outside this band are set to
      !!! HUGE. This has to be allocated to proper size incl. ghost layers
      !!! of size order! If ivalue is a function pointer, udata will be
      !!! completely replaced. If ivalue is a scalar, the points closer to
      !!! the interface than this scalar are kept unchanged. Can be a vector
      !!! or a scalar field. If vector, lda must be given.
#elif __TYPE == __VFIELD
      REAL(MK), DIMENSION(:,:,:,:,:) , POINTER         :: udata
      !!! Field of the function to be extended, defined in the same narrow
      !!! band as the level function. Values outside this band are set to
      !!! HUGE. This has to be allocated to proper size incl. ghost layers
      !!! of size order! If ivalue is a function pointer, udata will be
      !!! completely replaced. If ivalue is a scalar, the points closer to
      !!! the interface than this scalar are kept unchanged. Can be a vector
      !!! or a scalar field. If vector, lda must be given.
      INTEGER                        , INTENT(IN   )   :: lda
      !!! Only present for vector udata. Gives the length of the leading
      !!! dimension of udata. All elements will be equally extended.
#endif
      REAL(MK), DIMENSION(:,:,:,:)   , POINTER         :: fdata
      !!! Level function data. Needs to be defined in a narrow band of
      !!! specific width around the interface. The zero level is interpreted
      !!! to be the interface. Always a scalar field. THIS NEEDS TO BE
      !!! PROPERLY ALLOCATED ON INPUT, INCLUDING GHOST LAYERS OF SIZE order.
      REAL(MK), DIMENSION(:,:,:,:,:) , INTENT(IN),OPTIONAL:: chi
      !!! rank 5 (3d) field specifying the positions of the
      !!! grid nodes. 1st index: 1..ppm_dim, then i,j,[k],isub. OPTIONAL.
      !!! Uniform grid is assumed if absent. Ghostlayers of size >=1 must
      !!! be pre-filled.
#endif
      INTEGER                        , INTENT(IN   )   :: order
      !!! Desired order of the method. One of:
      !!!
      !!! *ppm_param_order_1
      !!! *ppm_param_order_2
      !!! *ppm_param_order_3
      REAL(MK)                       , INTENT(IN   )   :: tol
      !!! Relative tolerance for the determined distance to the interface.
      !!! 1E-3 is a good choice. The tolerance is in multiples of grid spacings
      REAL(MK)                       , INTENT(IN   )   :: width
      !!! Width of the narrow band to be produced on each side of the interface.
      INTEGER                        , INTENT(  OUT)   :: info
      !!! Return status. 0 upon success
      INTEGER , OPTIONAL             , INTENT(IN   )   :: MaxIter
      !!! OPTIONAL argument specifying the maximum number of allowed
      !!! iterations. This can be useful since a cyclic dependency in the
      !!! GMM algorithms could cause infinite loops. In each iteration at least
      !!! one point is computed.
      !-------------------------------------------------------------------------
      !  Local variables 
      !-------------------------------------------------------------------------
      INTEGER                               :: isub,npts,MaxIt
      INTEGER                               :: i,j,k,p
      INTEGER                               :: iopt,jsub,ida
      INTEGER , DIMENSION(4)                :: ldl,ldu
      REAL(MK)                              :: t0,x,y,z,big
      LOGICAL                               :: lok
      REAL(MK), DIMENSION(:,:), POINTER     :: closest
      TYPE(ppm_t_topo),         POINTER     :: topo
      TYPE(ppm_t_equi_mesh),    POINTER     :: mesh
#if   __TYPE == __VFIELD
#if   __DIM == __2D
      REAL(MK), DIMENSION(:,:,:  ), POINTER :: ext_wrk
#elif __DIM == __3D
      REAL(MK), DIMENSION(:,:,:,:), POINTER :: ext_wrk
#endif
#endif
      !-------------------------------------------------------------------------
      !  Initialise 
      !-------------------------------------------------------------------------
      CALL substart('ppm_gmm_extend',t0,info)
      big = HUGE(big)
      topo => ppm_topo(gmm_topoid)%t
      mesh => topo%mesh(gmm_meshid)
      !-------------------------------------------------------------------------
      !  Set pointers
      !-------------------------------------------------------------------------
#if   __KIND == __SINGLE_PRECISION
      closest => gmm_clos2
#elif __KIND == __DOUBLE_PRECISION
      closest => gmm_clod2
#endif
#if   __TYPE == __VFIELD
#if   __DIM == __2D
#if   __KIND == __SINGLE_PRECISION
      ext_wrk => ext_wrk_2ds
#elif __KIND == __DOUBLE_PRECISION
      ext_wrk => ext_wrk_2dd
#endif
#elif __DIM == __3D
#if   __KIND == __SINGLE_PRECISION
      ext_wrk => ext_wrk_3ds
#elif __KIND == __DOUBLE_PRECISION
      ext_wrk => ext_wrk_3dd
#endif
#endif
#endif
      !-------------------------------------------------------------------------
      !  Check arguments
      !-------------------------------------------------------------------------
      IF (ppm_debug .GT. 0) THEN
          IF (.NOT. ppm_initialized) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_ppm_noinit,'ppm_gmm_extend',  &
     &            'Please call ppm_init first!',__LINE__,info)
              GOTO 9999
          ENDIF
          IF (tol .LE. 0.0_MK) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_gmm_extend',  &
     &            'tolerance must be >0!',__LINE__,info)
              GOTO 9999
          ENDIF
          IF (width .LE. 0.0_MK) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_gmm_extend',  &
     &            'width must be >0!',__LINE__,info)
              GOTO 9999
          ENDIF
#if   __TYPE == __VFIELD
          IF (lda .LT. 1) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_gmm_extend',  &
     &            'lda must be >=1 for vector data',__LINE__,info)
              GOTO 9999
          ENDIF
#endif
      ENDIF          ! ppm_debug for argument check

#if   __KICKOFF == __YES
      !-------------------------------------------------------------------------
      !  Closest point transform on fdata
      !-------------------------------------------------------------------------
      IF (PRESENT(chi)) THEN
          CALL ppm_gmm_cpt(fdata,tol,npts,iptstmp2,closest,info,chi)
      ELSE
          CALL ppm_gmm_cpt(fdata,tol,npts,iptstmp2,closest,info)
      ENDIF
      IF (info .NE. ppm_param_success) THEN
          info = ppm_error_error
          CALL ppm_error(ppm_err_sub_failed,'ppm_gmm_extend',  &
     &        'Closest point transform failed.',__LINE__,info)
          GOTO 9999
      ENDIF
      !-------------------------------------------------------------------------
      !  Allocate udata
      !-------------------------------------------------------------------------
!     iopt = ppm_param_alloc_fit
!     ldl(1) = 1 - ghostsize(1)
!     ldl(2) = 1 - ghostsize(2)
!     ldl(3) = 1 - ghostsize(3)
!     ldl(4) = 1
!     ldu(1) = maxxhi + ghostsize(1)
!     ldu(2) = maxyhi + ghostsize(2)
!     ldu(3) = maxzhi + ghostsize(3)
!     ldu(4) = topo%nsublist
!     CALL ppm_alloc(udata,ldl,ldu,iopt,info)
!     IF (info .NE. ppm_param_success) THEN
!         info = ppm_error_fatal
!         CALL ppm_error(ppm_err_alloc,'ppm_gmm_extend',      &
!    &        'function data UDATA',__LINE__,info)
!         GOTO 9999
!     ENDIF
      udata = HUGE(x)
      !-------------------------------------------------------------------------
      !  Assign each point the value of its closest point on the
      !  interface.
      !-------------------------------------------------------------------------
#if   __DIM == __3D
      DO p=1,npts
          i = iptstmp2(1,p)
          j = iptstmp2(2,p)
          k = iptstmp2(3,p)
          isub = iptstmp2(4,p)
          x = closest(1,p)
          y = closest(2,p)
          z = closest(3,p)
#if   __TYPE == __SFIELD
          udata(i,j,k,isub) = ivalue(x,y,z,info)
#elif __TYPE == __VFIELD
          udata(1,i,j,k,isub) = ivalue(x,y,z,info)
          DO ida=2,lda
              udata(ida,i,j,k,isub) = udata(1,i,j,k,isub)
          ENDDO
#endif
      ENDDO
#elif __DIM == __2D
      DO p=1,npts
          i = iptstmp2(1,p)
          j = iptstmp2(2,p)
          isub = iptstmp2(3,p)
          x = closest(1,p)
          y = closest(2,p)
#if   __TYPE == __SFIELD
          udata(i,j,isub) = ivalue(x,y,info)
#elif __TYPE == __VFIELD
          udata(1,i,j,isub) = ivalue(x,y,info)
          DO ida=2,lda
              udata(ida,i,j,isub) = udata(1,i,j,isub)
          ENDDO
#endif
      ENDDO
#endif
      !-------------------------------------------------------------------------
      !  Deallocate
      !-------------------------------------------------------------------------
      iopt = ppm_param_dealloc
      CALL ppm_alloc(closest,ldu,iopt,info)
      IF (info .NE. ppm_param_success) THEN
          info = ppm_error_error
          CALL ppm_error(ppm_err_dealloc,'ppm_gmm_extend',      &
     &        'closest points locations CLOSEST',__LINE__,info)
      ENDIF
#if   __KIND == __SINGLE_PRECISION
      NULLIFY(gmm_clos2)
#elif __KIND == __DOUBLE_PRECISION
      NULLIFY(gmm_clod2)
#endif
      CALL ppm_alloc(iptstmp2,ldu,iopt,info)
      IF (info .NE. ppm_param_success) THEN
          info = ppm_error_error
          CALL ppm_error(ppm_err_dealloc,'ppm_gmm_extend',      &
     &        'close mesh points IPTS',__LINE__,info)
      ENDIF
#elif __KICKOFF == __NO
      !-------------------------------------------------------------------------
      !  Add ghost layers to field if needed
      !-------------------------------------------------------------------------
!     iopt = ppm_param_alloc_grow_preserve
!     ldl(1) = 1 - ghostsize(1)
!     ldl(2) = 1 - ghostsize(2)
!     ldl(3) = 1 - ghostsize(3)
!     ldl(4) = 1
!     ldu(1) = maxxhi + ghostsize(1)
!     ldu(2) = maxyhi + ghostsize(2)
!     ldu(3) = maxzhi + ghostsize(3)
!     ldu(4) = topo%nsublist
!     CALL ppm_alloc(udata,ldu,iopt,info)
!     IF (info .NE. ppm_param_success) THEN
!         info = ppm_error_fatal
!         CALL ppm_error(ppm_err_alloc,'ppm_gmm_extend',      &
!    &        'function data UDATA',__LINE__,info)
!         GOTO 9999
!     ENDIF
      !-------------------------------------------------------------------------
      !  Nuke points farther from the interface than ivalue
      !-------------------------------------------------------------------------
#if   __DIM == __3D
      DO isub=1,topo%nsublist
          jsub = topo%isublist(isub)
          DO k=1,mesh%nnodes(3,jsub)
              DO j=1,mesh%nnodes(2,jsub)
                  DO i=1,mesh%nnodes(1,jsub)
                      IF (ABS(fdata(i,j,k,isub)) .GT. ivalue) THEN
#if   __TYPE == __VFIELD
                          DO ida=1,lda
                              udata(ida,i,j,k,isub) = big
                          ENDDO
#elif __TYPE == __SFIELD
                          udata(i,j,k,isub) = big
#endif
                      ENDIF
                  ENDDO
              ENDDO
          ENDDO
      ENDDO
#elif __DIM == __2D
      DO isub=1,topo%nsublist
          jsub = topo%isublist(isub)
          DO j=1,mesh%nnodes(2,jsub)
              DO i=1,mesh%nnodes(1,jsub)
                  IF (ABS(fdata(i,j,isub)) .GT. ivalue) THEN
#if   __TYPE == __VFIELD
                      DO ida=1,lda
                          udata(ida,i,j,isub) = big
                      ENDDO
#elif __TYPE == __SFIELD
                      udata(i,j,isub) = big
#endif
                  ENDIF
              ENDDO
          ENDDO
      ENDDO
#endif
#endif
      !-------------------------------------------------------------------------
      !  Check the maximum number of allowed iterations
      !-------------------------------------------------------------------------
      IF (PRESENT(MaxIter)) THEN
          MaxIt = MaxIter
      ELSE
          MaxIt = HUGE(MaxIt)
      ENDIF
      !-------------------------------------------------------------------------
      !  Marching udata
      !-------------------------------------------------------------------------
#if   __TYPE == __VFIELD
      !-------------------------------------------------------------------------
      !  Allocate work array
      !-------------------------------------------------------------------------
      iopt = ppm_param_alloc_fit
      ldl(1) = LBOUND(udata,2)
      ldl(2) = LBOUND(udata,3)
#if   __DIM == __3D
      ldl(3) = LBOUND(udata,4)
#endif
      ldl(4) = 1
      ldu(1) = UBOUND(udata,2)
      ldu(2) = UBOUND(udata,3)
#if   __DIM == __3D
      ldu(3) = UBOUND(udata,4)
#endif
      ldu(4) = topo%nsublist
      CALL ppm_alloc(ext_wrk,ldl,ldu,iopt,info)
      IF (info .NE. ppm_param_success) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_gmm_extend',      &
     &        'work memory EXT_WRK',__LINE__,info)
          GOTO 9999
      ENDIF
      !-------------------------------------------------------------------------
      !  Do each vector dimension separately using the work memory
      !-------------------------------------------------------------------------
      DO ida=1,lda
#if   __DIM == __3D
          DO isub=1,ldu(4)
              DO k=ldl(3),ldu(3)
                  DO j=ldl(2),ldu(2)
                      DO i=ldl(1),ldu(1)
                          ext_wrk(i,j,k,isub) = udata(ida,i,j,k,isub)
                      ENDDO
                  ENDDO
              ENDDO
          ENDDO
#elif __DIM == __2D
          DO isub=1,ldu(4)
              DO j=ldl(2),ldu(2)
                  DO i=ldl(1),ldu(1)
                      ext_wrk(i,j,isub) = udata(ida,i,j,isub)
                  ENDDO
              ENDDO
          ENDDO
#endif
          !---------------------------------------------------------------------
          !  March
          !---------------------------------------------------------------------
          IF (PRESENT(chi)) THEN
              CALL ppm_gmm_march(width,order,fdata,0.0_MK,MaxIt,info,    &
     &            udata=ext_wrk,chi=chi)
          ELSE
              CALL ppm_gmm_march(width,order,fdata,0.0_MK,MaxIt,info,    &
     &            udata=ext_wrk)
          ENDIF
          IF (info .NE. ppm_param_success) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_sub_failed,'ppm_gmm_extend',  &
     &            'Marching GMM failed.',__LINE__,info)
              GOTO 9999
          ENDIF

          !---------------------------------------------------------------------
          !  The ghostsize might have changed in the marching, depending on
          !  the order of the marching algorithm.
          !---------------------------------------------------------------------
          ldl(1) = MAX(ldl(1),LBOUND(ext_wrk,1))
          ldl(2) = MAX(ldl(2),LBOUND(ext_wrk,2))
#if   __DIM == __3D
          ldl(3) = MAX(ldl(3),LBOUND(ext_wrk,3))
#endif
          ldu(1) = MIN(ldu(1),UBOUND(ext_wrk,1))
          ldu(2) = MIN(ldu(2),UBOUND(ext_wrk,2))
#if   __DIM == __3D
          ldu(3) = MIN(ldu(3),UBOUND(ext_wrk,3))
#endif
          !---------------------------------------------------------------------
          !  Copy results back
          !---------------------------------------------------------------------
#if   __DIM == __3D
          DO isub=1,ldu(4)
              DO k=ldl(3),ldu(3)
                  DO j=ldl(2),ldu(2)
                      DO i=ldl(1),ldu(1)
                          udata(ida,i,j,k,isub) = ext_wrk(i,j,k,isub)
                      ENDDO
                  ENDDO
              ENDDO
          ENDDO
#elif __DIM == __2D
          DO isub=1,ldu(4)
              DO j=ldl(2),ldu(2)
                  DO i=ldl(1),ldu(1)
                      udata(ida,i,j,isub) = ext_wrk(i,j,isub)
                  ENDDO
              ENDDO
          ENDDO
#endif
      ENDDO              ! ida
      !-------------------------------------------------------------------------
      !  Free work memory
      !-------------------------------------------------------------------------
      iopt = ppm_param_dealloc
      CALL ppm_alloc(ext_wrk,ldl,ldu,iopt,info)
      IF (info .NE. ppm_param_success) THEN
          info = ppm_error_error
          CALL ppm_error(ppm_err_dealloc,'ppm_gmm_extend',      &
     &        'work memory EXT_WRK',__LINE__,info)
      ENDIF
#if   __TYPE == __VFIELD
#if   __DIM == __2D
#if   __KIND == __SINGLE_PRECISION
      NULLIFY(ext_wrk_2ds)
#elif __KIND == __DOUBLE_PRECISION
      NULLIFY(ext_wrk_2dd)
#endif
#elif __DIM == __3D
#if   __KIND == __SINGLE_PRECISION
      NULLIFY(ext_wrk_3ds)
#elif __KIND == __DOUBLE_PRECISION
      NULLIFY(ext_wrk_3dd)
#endif
#endif
#endif
#elif __TYPE == __SFIELD
      !-------------------------------------------------------------------------
      !  Do scalar marching directly
      !-------------------------------------------------------------------------
      IF (PRESENT(chi)) THEN
          CALL ppm_gmm_march(width,order,fdata,0.0_MK,MaxIt,info,     &
     &        udata=udata,chi=chi)
      ELSE
          CALL ppm_gmm_march(width,order,fdata,0.0_MK,MaxIt,info,udata=udata)
      ENDIF
      IF (info .NE. ppm_param_success) THEN
          info = ppm_error_error
          CALL ppm_error(ppm_err_sub_failed,'ppm_gmm_extend',  &
     &        'Marching GMM failed.',__LINE__,info)
          GOTO 9999
      ENDIF
#endif
      !-------------------------------------------------------------------------
      !  Return 
      !-------------------------------------------------------------------------
 9999 CONTINUE
      CALL substop('ppm_gmm_extend',t0,info)
      RETURN
#if    __KICKOFF == __YES
#if    __DIM == __2D
#if    __TYPE == __SFIELD
#if    __KIND == __SINGLE_PRECISION
      END SUBROUTINE ppm_gmm_extend_2d_ksca_s
#elif  __KIND == __DOUBLE_PRECISION
      END SUBROUTINE ppm_gmm_extend_2d_ksca_d
#endif 
#elif  __TYPE == __VFIELD
#if    __KIND == __SINGLE_PRECISION
      END SUBROUTINE ppm_gmm_extend_2d_kvec_s
#elif  __KIND == __DOUBLE_PRECISION
      END SUBROUTINE ppm_gmm_extend_2d_kvec_d
#endif 
#endif 
#elif  __DIM == __3D
#if    __TYPE == __SFIELD
#if    __KIND == __SINGLE_PRECISION
      END SUBROUTINE ppm_gmm_extend_3d_ksca_s
#elif  __KIND == __DOUBLE_PRECISION
      END SUBROUTINE ppm_gmm_extend_3d_ksca_d
#endif 
#elif  __TYPE == __VFIELD
#if    __KIND == __SINGLE_PRECISION
      END SUBROUTINE ppm_gmm_extend_3d_kvec_s
#elif  __KIND == __DOUBLE_PRECISION
      END SUBROUTINE ppm_gmm_extend_3d_kvec_d
#endif 
#endif 
#endif
#elif  __KICKOFF == __NO
#if    __DIM == __2D
#if    __TYPE == __SFIELD
#if    __KIND == __SINGLE_PRECISION
      END SUBROUTINE ppm_gmm_extend_2d_tsca_s
#elif  __KIND == __DOUBLE_PRECISION
      END SUBROUTINE ppm_gmm_extend_2d_tsca_d
#endif 
#elif  __TYPE == __VFIELD
#if    __KIND == __SINGLE_PRECISION
      END SUBROUTINE ppm_gmm_extend_2d_tvec_s
#elif  __KIND == __DOUBLE_PRECISION
      END SUBROUTINE ppm_gmm_extend_2d_tvec_d
#endif 
#endif 
#elif  __DIM == __3D
#if    __TYPE == __SFIELD
#if    __KIND == __SINGLE_PRECISION
      END SUBROUTINE ppm_gmm_extend_3d_tsca_s
#elif  __KIND == __DOUBLE_PRECISION
      END SUBROUTINE ppm_gmm_extend_3d_tsca_d
#endif 
#elif  __TYPE == __VFIELD
#if    __KIND == __SINGLE_PRECISION
      END SUBROUTINE ppm_gmm_extend_3d_tvec_s
#elif  __KIND == __DOUBLE_PRECISION
      END SUBROUTINE ppm_gmm_extend_3d_tvec_d
#endif 
#endif 
#endif
#endif
