      !-------------------------------------------------------------------------
      !  Subroutine   :                 ppm_gmm_extend_bkwd
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
      SUBROUTINE ppm_gmm_extend_bkwd_2ds(fdta,dta,width,order,npos,TM,   &
     &    rhscst,dxinv,dyinv,dzinv,ghostsize,info,speed,chi)
#elif  __KIND == __DOUBLE_PRECISION
      SUBROUTINE ppm_gmm_extend_bkwd_2dd(fdta,dta,width,order,npos,TM,   &
     &    rhscst,dxinv,dyinv,dzinv,ghostsize,info,speed,chi)
#endif 
#elif  __DIM == __3D
#if    __KIND == __SINGLE_PRECISION
      SUBROUTINE ppm_gmm_extend_bkwd_3ds(fdta,dta,width,order,npos,TM,   &
     &    rhscst,dxinv,dyinv,dzinv,ghostsize,info,speed,chi)
#elif  __KIND == __DOUBLE_PRECISION
      SUBROUTINE ppm_gmm_extend_bkwd_3dd(fdta,dta,width,order,npos,TM,   &
     &    rhscst,dxinv,dyinv,dzinv,ghostsize,info,speed,chi)
#endif 
#endif
      !!! This routine performs the backward marching step of the GMM. See
      !!! ppm_gmm_march for details.
      !!!
      !!! === References ===
      !!!
      !!! Chopp:2001, Kim:2001b
      !-------------------------------------------------------------------------
      !  Modules 
      !-------------------------------------------------------------------------
      USE ppm_module_data
      USE ppm_module_numerics_data
      USE ppm_module_data_gmm
      USE ppm_module_substart
      USE ppm_module_substop
      USE ppm_module_error
      USE ppm_module_alloc
      USE ppm_module_typedef
      USE ppm_module_map_field
      IMPLICIT NONE
#if    __KIND == __SINGLE_PRECISION
      INTEGER, PARAMETER :: MK = ppm_kind_single
#else
      INTEGER, PARAMETER :: MK = ppm_kind_double
#endif
      !-------------------------------------------------------------------------
      !  Includes
      !-------------------------------------------------------------------------
#include "ppm_define.h"
#ifdef __MPI
      INCLUDE 'mpif.h'
#endif
      !-------------------------------------------------------------------------
      !  Arguments     
      !-------------------------------------------------------------------------
#if   __DIM == __2D
      REAL(MK), DIMENSION(:,:,:)     , POINTER          :: fdta
#elif __DIM == __3D
      REAL(MK), DIMENSION(:,:,:,:)   , POINTER          :: fdta
#endif
      !!! pointer to level function. Needs to be defined in a band
      !!! (width+order*dx).
#if   __DIM == __2D
      REAL(MK), DIMENSION(:,:,:)     , POINTER          :: dta
#elif __DIM == __3D
      REAL(MK), DIMENSION(:,:,:,:)   , POINTER          :: dta
#endif
      !!! pointer to value function.
#if   __DIM == __2D
      REAL(MK), DIMENSION(:,:,:)     , INTENT(IN), OPTIONAL :: speed
#elif __DIM == __3D
      REAL(MK), DIMENSION(:,:,:,:)   , INTENT(IN), OPTIONAL :: speed
#endif
      !!! rank 4 (3d) or rank 3 (2d) field of front speeds.
      !!! OPTIONAL to override rhscst.
#if   __DIM == __2D
      REAL(MK), DIMENSION(:,:,:,:)   , INTENT(IN), OPTIONAL :: chi
#elif __DIM == __3D
      REAL(MK), DIMENSION(:,:,:,:,:) , INTENT(IN), OPTIONAL :: chi
#endif
      !!! rank 5 (3d) or rank 4 (2d) field specifying the positions
      !!! of the grid nodes. 1st index: 1..ppm_dim, then i,j,[k],isub.
      !!! OPTIONAL. Uniform grid is assumed if absent.
      REAL(MK)                       , INTENT(IN   )    :: width
      !!! Width of the narrow band to be produced on each side of
      !!! the interface.
      REAL(MK)                       , INTENT(IN   )    :: rhscst
      !!! constant value for the right hand side of grad u * grad f = c.
      !!! If speed is present, this argument will be ignored.
      REAL(MK)                       , INTENT(IN   )    :: TM
      !!! Current threshold for wave front location.
      REAL(MK)                       , INTENT(IN   )    :: dxinv
      !!! inverse of the x grid spacing.
      REAL(MK)                       , INTENT(IN   )    :: dyinv
      !!! inverse of the y grid spacing.
      REAL(MK)                       , INTENT(IN   )    :: dzinv
      !!! inverse of the z grid spacing. (Not used in 2D version).
      INTEGER                        , INTENT(IN   )    :: order
      !!! Order of the method to be used. One of
      !!!
      !!! *ppm_param_order_1
      !!! *ppm_param_order_2
      !!! *ppm_param_order_3
      INTEGER, DIMENSION(3)          , INTENT(IN   )    :: ghostsize
      !!! Size of the ghostlayer on all sides.
      INTEGER                        , INTENT(INOUT)    :: npos
      !!! Current number of points in the close set.
      INTEGER                        , INTENT(  OUT)    :: info
      !!! Return status, 0 upon success
      !-------------------------------------------------------------------------
      !  Local variables 
      !-------------------------------------------------------------------------
      INTEGER                          :: i,j,k,p,xhi,yhi,zhi,ii,jj,kk
      INTEGER                          :: jsub,isub
      INTEGER                          :: i1,i2,i3
      INTEGER, DIMENSION(-3:3)         :: sx,sy,sz
      REAL(MK)                         :: t0,onethird,onetwelfth
      REAL(MK)                         :: valijk,det,hsave,fdta0
      REAL(MK)                         :: lmyeps,ainv,big,absfdta0
      REAL(MK)                         :: dxihalf,dxitwelve,agag
      REAL(MK)                         :: dyihalf,dyitwelve
      REAL(MK)                         :: dzihalf,dzitwelve
      REAL(MK), DIMENSION(3)           :: coefs,gphi,gpp
      REAL(MK), DIMENSION(3,3)         :: jac,ji
      REAL(MK), DIMENSION(-3:3,ppm_dim):: phi,psi
      REAL(MK), DIMENSION(ppm_dim)     :: alpha,beta
      REAL(MK), DIMENSION(2)           :: roots
      TYPE(ppm_t_topo),      POINTER   :: topo
      TYPE(ppm_t_equi_mesh), POINTER   :: mesh
      !-------------------------------------------------------------------------
      !  Externals 
      !-------------------------------------------------------------------------
      
      !-------------------------------------------------------------------------
      !  Initialise 
      !-------------------------------------------------------------------------
      CALL substart('ppm_gmm_extend_bkwd',t0,info)
      topo => ppm_topo(gmm_topoid)%t
      mesh => topo%mesh(gmm_meshid)
      phi      = 0.0_MK
      psi      = 0.0_MK
      big      = HUGE(big)
      hsave    = 0.9_MK*big
      onethird = 1.0_MK/3.0_MK
      onetwelfth = 1.0_MK/12.0_MK
#if   __KIND == __SINGLE_PRECISION
      lmyeps   = ppm_myepss
#else
      lmyeps   = ppm_myepsd
#endif
      dxihalf  = 0.5_MK*dxinv
      dyihalf  = 0.5_MK*dyinv
      dxitwelve  = onetwelfth*dxinv
      dyitwelve  = onetwelfth*dyinv
#if   __DIM == __3D
      dzihalf  = 0.5_MK*dzinv
      dzitwelve  = onetwelfth*dzinv
#endif 
      !-------------------------------------------------------------------------
      !  Check arguments
      !-------------------------------------------------------------------------
      IF (ppm_debug .GT. 0) THEN
          IF (width .LT. 0.0_MK) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_ppm_noinit,'ppm_gmm_extend_bkwd',  &
     &            'width must be positive!',__LINE__,info)
              GOTO 9999
          ENDIF
          IF ((order.NE.ppm_param_order_1).AND.(order.NE.ppm_param_order_2)  &
     &         .AND.(order.NE.ppm_param_order_3)) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_ppm_noinit,'ppm_gmm_extend_bkwd',  &
     &            'order must be 1, 2, or 3!',__LINE__,info)
              GOTO 9999
          ENDIF
      ENDIF          ! ppm_debug for argument check
#if   __DIM == __3D
      !-------------------------------------------------------------------------
      !  Reverse order: recompute neighbors of points in ggm_ipos
      !-------------------------------------------------------------------------
      DO p=npos,1,-1
          ii   = gmm_ipos(1,p)
          jj   = gmm_ipos(2,p)
          kk   = gmm_ipos(3,p)
          jsub = gmm_ipos(4,p)
          isub = topo%isublist(jsub)
          xhi  = mesh%nnodes(1,isub)
          yhi  = mesh%nnodes(2,isub)
          zhi  = mesh%nnodes(3,isub)
          fdta0= fdta(ii,jj,kk,jsub)
          absfdta0 = fdta0
          IF (absfdta0 .LT. 0.0_MK) absfdta0 = -absfdta0
          !---------------------------------------------------------------------
          !  GMM update condition (see Kim:2001a)
          !---------------------------------------------------------------------
          IF (.NOT.(absfdta0.GT.TM)) THEN
              !-----------------------------------------------------------------
              !  Recompute non-accepted close neighbors
              !-----------------------------------------------------------------
              i = ii - 1
              j = jj
              k = kk
              IF (i.GT.0) THEN
                  IF ((gmm_state3d(i,j,k,jsub) .NE.      &
     &                ppm_gmm_param_accepted) .AND.      &
     &                (ABS(fdta(i,j,k,jsub)).GT.absfdta0)) THEN
                      !---------------------------------------------------------
                      !  Update point i,j,k
                      !---------------------------------------------------------
#include "ppm_gmm_slvextn.inc"
                      IF (ABS(valijk) .LT. hsave) dta(i,j,k,jsub) = valijk 
                  ENDIF
              ENDIF
              i = ii + 1
              j = jj
              k = kk
              IF (i.LE.xhi) THEN
                  IF ((gmm_state3d(i,j,k,jsub) .NE.      &
     &                ppm_gmm_param_accepted) .AND.      &
     &                (ABS(fdta(i,j,k,jsub)).GT.absfdta0)) THEN
                      !---------------------------------------------------------
                      !  Update point i,j,k
                      !---------------------------------------------------------
#include "ppm_gmm_slvextn.inc"
                      IF (ABS(valijk) .LT. hsave) dta(i,j,k,jsub) = valijk 
                  ENDIF
              ENDIF
              i = ii
              j = jj - 1 
              k = kk
              IF (j.GT.0) THEN
                  IF ((gmm_state3d(i,j,k,jsub) .NE.      &
     &                ppm_gmm_param_accepted) .AND.      &
     &                (ABS(fdta(i,j,k,jsub)).GT.absfdta0)) THEN
                      !---------------------------------------------------------
                      !  Update point i,j,k
                      !---------------------------------------------------------
#include "ppm_gmm_slvextn.inc"
                      IF (ABS(valijk) .LT. hsave) dta(i,j,k,jsub) = valijk 
                  ENDIF
              ENDIF
              i = ii
              j = jj + 1 
              k = kk
              IF (j.LE.yhi) THEN
                  IF ((gmm_state3d(i,j,k,jsub) .NE.      &
     &                ppm_gmm_param_accepted) .AND.      &
     &                (ABS(fdta(i,j,k,jsub)).GT.absfdta0)) THEN
                      !---------------------------------------------------------
                      !  Update point i,j,k
                      !---------------------------------------------------------
#include "ppm_gmm_slvextn.inc"
                      IF (ABS(valijk) .LT. hsave) dta(i,j,k,jsub) = valijk 
                  ENDIF
              ENDIF
              i = ii
              j = jj 
              k = kk - 1
              IF (k.GT.0) THEN
                  IF ((gmm_state3d(i,j,k,jsub) .NE.      &
     &                ppm_gmm_param_accepted) .AND.      &
     &                (ABS(fdta(i,j,k,jsub)).GT.absfdta0)) THEN
                      !---------------------------------------------------------
                      !  Update point i,j,k
                      !---------------------------------------------------------
#include "ppm_gmm_slvextn.inc"
                      IF (ABS(valijk) .LT. hsave) dta(i,j,k,jsub) = valijk 
                  ENDIF
              ENDIF
              i = ii
              j = jj 
              k = kk + 1
              IF (k.LE.zhi) THEN
                  IF ((gmm_state3d(i,j,k,jsub) .NE.      &
     &                ppm_gmm_param_accepted) .AND.      &
     &                (ABS(fdta(i,j,k,jsub)).GT.absfdta0)) THEN
                      !---------------------------------------------------------
                      !  Update point i,j,k
                      !---------------------------------------------------------
#include "ppm_gmm_slvextn.inc"
                      IF (ABS(valijk) .LT. hsave) dta(i,j,k,jsub) = valijk 
                  ENDIF
              ENDIF
          ENDIF           ! TT .LE. TM
      ENDDO         ! p=npos,1,-1
      !-------------------------------------------------------------------------
      !  Update ghost layers for dta
      !-------------------------------------------------------------------------
      CALL ppm_map_field_push(gmm_topoid,gmm_meshid,dta,info)
      IF (info .NE. ppm_param_success) THEN
          info = ppm_error_error
          CALL ppm_error(ppm_err_sub_failed,'ppm_gmm_extend',  &
     &        'pushing field data failed',__LINE__,info)
          GOTO 9999
      ENDIF
      CALL ppm_map_field_send(info)
      IF (info .NE. ppm_param_success) THEN
          info = ppm_error_error
          CALL ppm_error(ppm_err_sub_failed,'ppm_gmm_extend',  &
     &        'sending ghosts failed',__LINE__,info)
          GOTO 9999
      ENDIF
      CALL ppm_map_field_pop(gmm_topoid,gmm_meshid,dta,ghostsize,info)
      IF (info .NE. ppm_param_success) THEN
          info = ppm_error_error
          CALL ppm_error(ppm_err_sub_failed,'ppm_gmm_extend',  &
     &        'popping field data failed',__LINE__,info)
          GOTO 9999
      ENDIF

#elif __DIM == __2D
      !-------------------------------------------------------------------------
      !  Reverse order: recompute neighbors of points in ggm_ipos
      !-------------------------------------------------------------------------
      DO p=npos,1,-1
          ii   = gmm_ipos(1,p)
          jj   = gmm_ipos(2,p)
          jsub = gmm_ipos(3,p)
          isub = topo%isublist(jsub)
          xhi  = mesh%nnodes(1,isub)
          yhi  = mesh%nnodes(2,isub)
          fdta0= fdta(ii,jj,jsub)
          !---------------------------------------------------------------------
          !  GMM update condition (see Kim:2001a)
          !---------------------------------------------------------------------
          IF (.NOT.(ABS(fdta0).GT.TM)) THEN
              !-----------------------------------------------------------------
              !  Recompute non-accepted close neighbors
              !-----------------------------------------------------------------
              i = ii - 1
              j = jj
              IF (i.GT.0) THEN
                  IF ((gmm_state2d(i,j,jsub) .NE.      &
     &                ppm_gmm_param_accepted) .AND.      &
     &                (ABS(fdta(i,j,jsub)).GT.ABS(fdta0))) THEN
                      !---------------------------------------------------------
                      !  Update point i,j
                      !---------------------------------------------------------
#include "ppm_gmm_slvextn.inc"
                      IF (valijk .LT. hsave) dta(i,j,jsub) = valijk 
                  ENDIF
              ENDIF
              i = ii + 1
              j = jj
              IF (i.LE.xhi) THEN
                  IF ((gmm_state2d(i,j,jsub) .NE.      &
     &                ppm_gmm_param_accepted) .AND.      &
     &                (ABS(fdta(i,j,jsub)).GT.ABS(fdta0))) THEN
                      !---------------------------------------------------------
                      !  Update point i,j
                      !---------------------------------------------------------
#include "ppm_gmm_slvextn.inc"
                      IF (valijk .LT. hsave) dta(i,j,jsub) = valijk 
                  ENDIF
              ENDIF
              i = ii
              j = jj - 1 
              IF (j.GT.0) THEN
                  IF ((gmm_state2d(i,j,jsub) .NE.      &
     &                ppm_gmm_param_accepted) .AND.      &
     &                (ABS(fdta(i,j,jsub)).GT.ABS(fdta0))) THEN
                      !---------------------------------------------------------
                      !  Update point i,j
                      !---------------------------------------------------------
#include "ppm_gmm_slvextn.inc"
                      IF (valijk .LT. hsave) dta(i,j,jsub) = valijk 
                  ENDIF
              ENDIF
              i = ii
              j = jj + 1 
              IF (j.LE.yhi) THEN
                  IF ((gmm_state2d(i,j,jsub) .NE.      &
     &                ppm_gmm_param_accepted) .AND.      &
     &                (ABS(fdta(i,j,jsub)).GT.ABS(fdta0))) THEN
                      !---------------------------------------------------------
                      !  Update point i,j
                      !---------------------------------------------------------
#include "ppm_gmm_slvextn.inc"
                      IF (valijk .LT. hsave) dta(i,j,jsub) = valijk 
                  ENDIF
              ENDIF
          ENDIF           ! TT .LE. TM
      ENDDO         ! p=npos,1,-1
      !-------------------------------------------------------------------------
      !  Update ghost layers for dta
      !-------------------------------------------------------------------------
      CALL ppm_map_field_push(gmm_topoid,gmm_meshid,dta,info)
      IF (info .NE. ppm_param_success) THEN
          info = ppm_error_error
          CALL ppm_error(ppm_err_sub_failed,'ppm_gmm_extend',  &
     &        'pushing field data failed',__LINE__,info)
          GOTO 9999
      ENDIF
      CALL ppm_map_field_send(info)
      IF (info .NE. ppm_param_success) THEN
          info = ppm_error_error
          CALL ppm_error(ppm_err_sub_failed,'ppm_gmm_extend',  &
     &        'sending ghosts failed',__LINE__,info)
          GOTO 9999
      ENDIF
      CALL ppm_map_field_pop(gmm_topoid,gmm_meshid,dta,ghostsize,info)
      IF (info .NE. ppm_param_success) THEN
          info = ppm_error_error
          CALL ppm_error(ppm_err_sub_failed,'ppm_gmm_extend',  &
     &        'popping field data failed',__LINE__,info)
          GOTO 9999
      ENDIF
#endif 
      !-------------------------------------------------------------------------
      !  Return 
      !-------------------------------------------------------------------------
 9999 CONTINUE
      CALL substop('ppm_gmm_extend_bkwd',t0,info)
      RETURN

#if    __DIM == __2D
#if    __KIND == __SINGLE_PRECISION
      END SUBROUTINE ppm_gmm_extend_bkwd_2ds
#elif  __KIND == __DOUBLE_PRECISION
      END SUBROUTINE ppm_gmm_extend_bkwd_2dd
#endif 
#elif  __DIM == __3D
#if    __KIND == __SINGLE_PRECISION
      END SUBROUTINE ppm_gmm_extend_bkwd_3ds
#elif  __KIND == __DOUBLE_PRECISION
      END SUBROUTINE ppm_gmm_extend_bkwd_3dd
#endif 
#endif
