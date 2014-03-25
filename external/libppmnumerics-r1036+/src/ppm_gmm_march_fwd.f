      !-------------------------------------------------------------------------
      !  Subroutine   :                 ppm_gmm_march_fwd
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
      SUBROUTINE ppm_gmm_march_fwd_2ds(fdta,width,order,npos,TM,    &
     &    rhscst,dxinv,dyinv,dzinv,ghostsize,info,speed,chi)
#elif  __KIND == __DOUBLE_PRECISION
      SUBROUTINE ppm_gmm_march_fwd_2dd(fdta,width,order,npos,TM,    &
     &    rhscst,dxinv,dyinv,dzinv,ghostsize,info,speed,chi)
#endif 
#elif  __DIM == __3D
#if    __KIND == __SINGLE_PRECISION
      SUBROUTINE ppm_gmm_march_fwd_3ds(fdta,width,order,npos,TM,    &
     &    rhscst,dxinv,dyinv,dzinv,ghostsize,info,speed,chi)
#elif  __KIND == __DOUBLE_PRECISION
      SUBROUTINE ppm_gmm_march_fwd_3dd(fdta,width,order,npos,TM,    &
     &    rhscst,dxinv,dyinv,dzinv,ghostsize,info,speed,chi)
#endif 
#endif
      !!! This routine performs the forward marching step of the GMM. See
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
      REAL(MK), DIMENSION(:,:,:)     , POINTER              :: fdta
      !!! pointer to level function.
      REAL(MK), DIMENSION(:,:,:)     , INTENT(IN), OPTIONAL :: speed
      !!! rank 4 (3d) or rank 3 (2d) field of front speeds.
      !!! OPTIONAL to override rhscst.
      REAL(MK), DIMENSION(:,:,:,:)   , INTENT(IN), OPTIONAL :: chi
      !!! rank 5 (3d) or rank 4 (2d) field specifying the positions of the
      !!! grid nodes. 1st index: 1..ppm_dim, then i,j,[k],isub.
      !!! OPTIONAL. Uniform grid is assumed if absent.
#elif __DIM == __3D
      REAL(MK), DIMENSION(:,:,:,:)   , POINTER              :: fdta
      !!! pointer to level function.
      REAL(MK), DIMENSION(:,:,:,:)   , INTENT(IN), OPTIONAL :: speed
      !!! rank 4 (3d) or rank 3 (2d) field of front speeds.
      !!! OPTIONAL to override rhscst.
      REAL(MK), DIMENSION(:,:,:,:,:) , INTENT(IN), OPTIONAL :: chi
      !!! rank 5 (3d) or rank 4 (2d) field specifying the positions of the
      !!! grid nodes. 1st index: 1..ppm_dim, then i,j,[k],isub.
      !!! OPTIONAL. Uniform grid is assumed if absent.
#endif
      REAL(MK)                       , INTENT(IN   )    :: width
      !!! Width of the narrow band to be produced on each side of the interface.
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
      INTEGER, DIMENSION(3)          , INTENT(IN   )    :: ghostsize
      !!! Size of the ghostlayer on all sides.
      INTEGER                        , INTENT(IN   )    :: order
      !!! Order of the method to be used. One of
      !!!
      !!! *ppm_param_order_1
      !!! *ppm_param_order_2
      !!! *ppm_param_order_3
      INTEGER                        , INTENT(INOUT)    :: npos
      !!! Current number of points in the close set.
      INTEGER                        , INTENT(  OUT)    :: info
      !!! Return status, 0 upon success
      !-------------------------------------------------------------------------
      !  Local variables 
      !-------------------------------------------------------------------------
      INTEGER                          :: i,j,k,p,xhi,yhi,zhi,ii,jj,kk
      INTEGER                          :: npos0,jsub,isub,iopt
      INTEGER                          :: i1,i2,i3
      INTEGER, DIMENSION(-3:3)         :: sx,sy,sz
      INTEGER, DIMENSION(4)            :: ldu
      REAL(MK)                         :: t0,onethird,onetwelfth
      REAL(MK)                         :: dxihalf,dyihalf,dzihalf
      REAL(MK)                         :: dxitwelve,dyitwelve,dzitwelve
      REAL(MK)                         :: valijk,det,hsave,fdta0
      REAL(MK)                         :: lmyeps,ainv,big,absfdta0
      REAL(MK), DIMENSION(3)           :: coefs
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
      CALL substart('ppm_gmm_march_fwd',t0,info)
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
              CALL ppm_error(ppm_err_ppm_noinit,'ppm_gmm_march_fwd',  &
     &            'width must be positive!',__LINE__,info)
              GOTO 9999
          ENDIF
          IF ((order.NE.ppm_param_order_1).AND.(order.NE.ppm_param_order_2)  &
     &         .AND.(order.NE.ppm_param_order_3)) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_ppm_noinit,'ppm_gmm_march_fwd',  &
     &            'order must be 1, 2, or 3!',__LINE__,info)
              GOTO 9999
          ENDIF
      ENDIF          ! ppm_debug for argument check

#if   __DIM == __3D
      npos0 = npos
      !-------------------------------------------------------------------------
      !  Forward order: update neighbors of points in ggm_ipos
      !-------------------------------------------------------------------------
      DO p=1,npos
          IF (p .GT. npos0) EXIT
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
              !  Compute non-accepted neighbors
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
#include "ppm_gmm_slvupwnd.inc"
                      IF (valijk .LT. hsave) THEN
                          IF (ABS(valijk).LT.ABS(fdta(i,j,k,jsub))) THEN
                              IF ((valijk*fdta(i,j,k,jsub).GE.0.0_MK) .OR.   &
     &                            (fdta(i,j,k,jsub).GT.hsave)) THEN 
                                  fdta(i,j,k,jsub) = valijk 
                              ENDIF
                          ENDIF
                          !-----------------------------------------------------
                          !  Keep in or add to close set
                          !-----------------------------------------------------
                          IF ((gmm_state3d(i,j,k,jsub) .EQ.     &
     &                        ppm_gmm_param_far) .AND.          &
     &                        (ABS(fdta(i,j,k,jsub)) .LT. width)) THEN
                              gmm_state3d(i,j,k,jsub) =         &
     &                            ppm_gmm_param_close
#include "ppm_gmm_add_to_list.inc"
                          ENDIF
                      ENDIF
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
#include "ppm_gmm_slvupwnd.inc"
                      IF (valijk .LT. hsave) THEN
                          IF (ABS(valijk).LT.ABS(fdta(i,j,k,jsub))) THEN
                              IF ((valijk*fdta(i,j,k,jsub).GE.0.0_MK) .OR.   &
     &                            (fdta(i,j,k,jsub).GT.hsave)) THEN 
                                  fdta(i,j,k,jsub) = valijk 
                              ENDIF
                          ENDIF
                          !-----------------------------------------------------
                          !  Keep in or add to close set
                          !-----------------------------------------------------
                          IF ((gmm_state3d(i,j,k,jsub) .EQ.     &
     &                        ppm_gmm_param_far) .AND.          &
     &                        (ABS(fdta(i,j,k,jsub)) .LT. width)) THEN
                              gmm_state3d(i,j,k,jsub) =         &
     &                            ppm_gmm_param_close
#include "ppm_gmm_add_to_list.inc"
                          ENDIF
                      ENDIF
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
#include "ppm_gmm_slvupwnd.inc"
                      IF (valijk .LT. hsave) THEN
                          IF (ABS(valijk).LT.ABS(fdta(i,j,k,jsub))) THEN
                              IF ((valijk*fdta(i,j,k,jsub).GE.0.0_MK) .OR.   &
     &                            (fdta(i,j,k,jsub).GT.hsave)) THEN 
                                  fdta(i,j,k,jsub) = valijk 
                              ENDIF
                          ENDIF
                          !-----------------------------------------------------
                          !  Keep in or add to close set
                          !-----------------------------------------------------
                          IF ((gmm_state3d(i,j,k,jsub) .EQ.     &
     &                        ppm_gmm_param_far) .AND.          &
     &                        (ABS(fdta(i,j,k,jsub)) .LT. width)) THEN
                              gmm_state3d(i,j,k,jsub) =         &
     &                            ppm_gmm_param_close
#include "ppm_gmm_add_to_list.inc"
                          ENDIF
                      ENDIF
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
#include "ppm_gmm_slvupwnd.inc"
                      IF (valijk .LT. hsave) THEN
                          IF (ABS(valijk).LT.ABS(fdta(i,j,k,jsub))) THEN
                              IF ((valijk*fdta(i,j,k,jsub).GE.0.0_MK) .OR.   &
     &                            (fdta(i,j,k,jsub).GT.hsave)) THEN 
                                  fdta(i,j,k,jsub) = valijk 
                              ENDIF
                          ENDIF
                          !-----------------------------------------------------
                          !  Keep in or add to close set
                          !-----------------------------------------------------
                          IF ((gmm_state3d(i,j,k,jsub) .EQ.     &
     &                        ppm_gmm_param_far) .AND.          &
     &                        (ABS(fdta(i,j,k,jsub)) .LT. width)) THEN
                              gmm_state3d(i,j,k,jsub) =         &
     &                            ppm_gmm_param_close
#include "ppm_gmm_add_to_list.inc"
                          ENDIF
                      ENDIF
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
#include "ppm_gmm_slvupwnd.inc"
                      IF (valijk .LT. hsave) THEN
                          IF (ABS(valijk).LT.ABS(fdta(i,j,k,jsub))) THEN
                              IF ((valijk*fdta(i,j,k,jsub).GE.0.0_MK) .OR.   &
     &                            (fdta(i,j,k,jsub).GT.hsave)) THEN 
                                  fdta(i,j,k,jsub) = valijk 
                              ENDIF
                          ENDIF
                          !-----------------------------------------------------
                          !  Keep in or add to close set
                          !-----------------------------------------------------
                          IF ((gmm_state3d(i,j,k,jsub) .EQ.     &
     &                        ppm_gmm_param_far) .AND.          &
     &                        (ABS(fdta(i,j,k,jsub)) .LT. width)) THEN
                              gmm_state3d(i,j,k,jsub) =         &
     &                            ppm_gmm_param_close
#include "ppm_gmm_add_to_list.inc"
                          ENDIF
                      ENDIF
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
#include "ppm_gmm_slvupwnd.inc"
                      IF (valijk .LT. hsave) THEN
                          IF (ABS(valijk).LT.ABS(fdta(i,j,k,jsub))) THEN
                              IF ((valijk*fdta(i,j,k,jsub).GE.0.0_MK) .OR.   &
     &                            (fdta(i,j,k,jsub).GT.hsave)) THEN 
                                  fdta(i,j,k,jsub) = valijk 
                              ENDIF
                          ENDIF
                          !-----------------------------------------------------
                          !  Keep in or add to close set
                          !-----------------------------------------------------
                          IF ((gmm_state3d(i,j,k,jsub) .EQ.     &
     &                        ppm_gmm_param_far) .AND.          &
     &                        (ABS(fdta(i,j,k,jsub)) .LT. width)) THEN
                              gmm_state3d(i,j,k,jsub) =         &
     &                            ppm_gmm_param_close
#include "ppm_gmm_add_to_list.inc"
                          ENDIF
                      ENDIF
                  ENDIF
              ENDIF
              !-----------------------------------------------------------------
              !  Accept this point and remove it from the list
              !-----------------------------------------------------------------
              gmm_state3d(ii,jj,kk,jsub) = ppm_gmm_param_accepted
              gmm_ipos(1,p) = gmm_ipos(1,npos0)
              gmm_ipos(2,p) = gmm_ipos(2,npos0)
              gmm_ipos(3,p) = gmm_ipos(3,npos0)
              gmm_ipos(4,p) = gmm_ipos(4,npos0)
              npos0 = npos0 - 1
          ENDIF           ! TT .LE. TM
      ENDDO          ! p=1,npos
      npos = npos0
      !       WRITE(cbuf,'(A,I2.2,A)') 'state_',Marchit,'.out'
      !       OPEN(40,FILE=cbuf,STATUS='REPLACE',ACTION='WRITE')
      !       WRITE(cbuf,'(A,I2.2,A)') 'value_',Marchit,'.out'
      !       OPEN(30,FILE=cbuf,STATUS='REPLACE',ACTION='WRITE')
      !       DO kk=1,zhi
      !           DO jj=1,yhi
      !               DO ii=1,xhi
      !                   WRITE(40,'(I3)') gmm_state3d(ii,jj,kk,jsub)
      !                   IF (fdta(ii,jj,kk,jsub) .GT. hsave) THEN
      !                       WRITE(30,'(E20.8)') 0.0_MK
      !                   ELSE
      !                       WRITE(30,'(E20.8)') fdta(ii,jj,kk,jsub)
      !                   ENDIF
      !               ENDDO
      !           ENDDO
      !       ENDDO
      !       CLOSE(30)
      !       CLOSE(40)
      !-------------------------------------------------------------------------
      !  Update ghost layers for both fdta AND gmm_state3d 
      !-------------------------------------------------------------------------
      CALL ppm_map_field_push(gmm_topoid,gmm_meshid,fdta,info)
      IF (info .NE. ppm_param_success) THEN
          info = ppm_error_error
          CALL ppm_error(ppm_err_sub_failed,'ppm_gmm_march_fwd',  &
     &        'pushing field data failed',__LINE__,info)
          GOTO 9999
      ENDIF
      CALL ppm_map_field_push(gmm_topoid,gmm_meshid,gmm_state3d,info)
      IF (info .NE. ppm_param_success) THEN
          info = ppm_error_error
          CALL ppm_error(ppm_err_sub_failed,'ppm_gmm_march_fwd',  &
     &        'pushing status data failed',__LINE__,info)
          GOTO 9999
      ENDIF
      CALL ppm_map_field_send(info)
      IF (info .NE. ppm_param_success) THEN
          info = ppm_error_error
          CALL ppm_error(ppm_err_sub_failed,'ppm_gmm_march_fwd',  &
     &        'sending ghosts failed',__LINE__,info)
          GOTO 9999
      ENDIF
      CALL ppm_map_field_pop(gmm_topoid,gmm_meshid,gmm_state3d,ghostsize,info)
      IF (info .NE. ppm_param_success) THEN
          info = ppm_error_error
          CALL ppm_error(ppm_err_sub_failed,'ppm_gmm_march_fwd',  &
     &        'popping status data failed',__LINE__,info)
          GOTO 9999
      ENDIF
      CALL ppm_map_field_pop(gmm_topoid,gmm_meshid,fdta,ghostsize,info)
      IF (info .NE. ppm_param_success) THEN
          info = ppm_error_error
          CALL ppm_error(ppm_err_sub_failed,'ppm_gmm_march_fwd',  &
     &        'popping field data failed',__LINE__,info)
          GOTO 9999
      ENDIF
#elif __DIM == __2D
      npos0 = npos
      !-------------------------------------------------------------------------
      !  Forward order: update neighbors of points in ggm_ipos
      !-------------------------------------------------------------------------
      DO p=1,npos
          IF (p .GT. npos0) EXIT
          ii   = gmm_ipos(1,p)
          jj   = gmm_ipos(2,p)
          jsub = gmm_ipos(3,p)
          isub = topo%isublist(jsub)
          xhi  = mesh%nnodes(1,isub)
          yhi  = mesh%nnodes(2,isub)
          fdta0= fdta(ii,jj,jsub)
          absfdta0 = fdta0
          IF (absfdta0 .LT. 0.0_MK) absfdta0 = -absfdta0
          !---------------------------------------------------------------------
          !  GMM update condition (see Kim:2001a)
          !---------------------------------------------------------------------
          IF (.NOT.(absfdta0.GT.TM)) THEN
              !-----------------------------------------------------------------
              !  Compute non-accepted neighbors
              !-----------------------------------------------------------------
              i = ii - 1
              j = jj 
              IF (i.GT.0) THEN
                  IF ((gmm_state2d(i,j,jsub) .NE.      &
     &                ppm_gmm_param_accepted) .AND.      &
     &                (ABS(fdta(i,j,jsub)).GT.absfdta0)) THEN
                      !---------------------------------------------------------
                      !  Update point i,j
                      !---------------------------------------------------------
#include "ppm_gmm_slvupwnd.inc"
                      IF (valijk .LT. hsave) THEN
                          IF (ABS(valijk).LT.ABS(fdta(i,j,jsub))) THEN
                              IF ((valijk*fdta(i,j,jsub).GE.0.0_MK) .OR.   &
     &                            (fdta(i,j,jsub).GT.hsave)) THEN 
                                  fdta(i,j,jsub) = valijk 
                              ENDIF
                          ENDIF
                          !-----------------------------------------------------
                          !  Keep in or add to close set
                          !-----------------------------------------------------
                          IF ((gmm_state2d(i,j,jsub) .EQ.     &
     &                        ppm_gmm_param_far) .AND.          &
     &                        (ABS(fdta(i,j,jsub)) .LT. width)) THEN
                              gmm_state2d(i,j,jsub) =         &
     &                            ppm_gmm_param_close
#include "ppm_gmm_add_to_list.inc"
                          ENDIF
                      ENDIF
                  ENDIF
              ENDIF
              i = ii + 1
              j = jj 
              IF (i.LE.xhi) THEN
                  IF ((gmm_state2d(i,j,jsub) .NE.      &
     &                ppm_gmm_param_accepted) .AND.      &
     &                (ABS(fdta(i,j,jsub)).GT.absfdta0)) THEN
                      !---------------------------------------------------------
                      !  Update point i,j
                      !---------------------------------------------------------
#include "ppm_gmm_slvupwnd.inc"
                      IF (valijk .LT. hsave) THEN
                          IF (ABS(valijk).LT.ABS(fdta(i,j,jsub))) THEN
                              IF ((valijk*fdta(i,j,jsub).GE.0.0_MK) .OR.   &
     &                            (fdta(i,j,jsub).GT.hsave)) THEN 
                                  fdta(i,j,jsub) = valijk 
                              ENDIF
                          ENDIF
                          !-----------------------------------------------------
                          !  Keep in or add to close set
                          !-----------------------------------------------------
                          IF ((gmm_state2d(i,j,jsub) .EQ.     &
     &                        ppm_gmm_param_far) .AND.          &
     &                        (ABS(fdta(i,j,jsub)) .LT. width)) THEN
                              gmm_state2d(i,j,jsub) =         &
     &                            ppm_gmm_param_close
#include "ppm_gmm_add_to_list.inc"
                          ENDIF
                      ENDIF
                  ENDIF
              ENDIF
              i = ii 
              j = jj - 1
              IF (j.GT.0) THEN
                  IF ((gmm_state2d(i,j,jsub) .NE.      &
     &                ppm_gmm_param_accepted) .AND.      &
     &                (ABS(fdta(i,j,jsub)).GT.absfdta0)) THEN
                      !---------------------------------------------------------
                      !  Update point i,j
                      !---------------------------------------------------------
#include "ppm_gmm_slvupwnd.inc"
                      IF (valijk .LT. hsave) THEN
                          IF (ABS(valijk).LT.ABS(fdta(i,j,jsub))) THEN
                              IF ((valijk*fdta(i,j,jsub).GE.0.0_MK) .OR.   &
     &                            (fdta(i,j,jsub).GT.hsave)) THEN 
                                  fdta(i,j,jsub) = valijk 
                              ENDIF
                          ENDIF
                          !-----------------------------------------------------
                          !  Keep in or add to close set
                          !-----------------------------------------------------
                          IF ((gmm_state2d(i,j,jsub) .EQ.     &
     &                        ppm_gmm_param_far) .AND.          &
     &                        (ABS(fdta(i,j,jsub)) .LT. width)) THEN
                              gmm_state2d(i,j,jsub) =         &
     &                            ppm_gmm_param_close
#include "ppm_gmm_add_to_list.inc"
                          ENDIF
                      ENDIF
                  ENDIF
              ENDIF
              i = ii 
              j = jj + 1
              IF (j.LE.yhi) THEN
                  IF ((gmm_state2d(i,j,jsub) .NE.      &
     &                ppm_gmm_param_accepted) .AND.      &
     &                (ABS(fdta(i,j,jsub)).GT.absfdta0)) THEN
                      !---------------------------------------------------------
                      !  Update point i,j
                      !---------------------------------------------------------
#include "ppm_gmm_slvupwnd.inc"
                      IF (valijk .LT. hsave) THEN
                          IF (ABS(valijk).LT.ABS(fdta(i,j,jsub))) THEN
                              IF ((valijk*fdta(i,j,jsub).GE.0.0_MK) .OR.   &
     &                            (fdta(i,j,jsub).GT.hsave)) THEN 
                                  fdta(i,j,jsub) = valijk 
                              ENDIF
                          ENDIF
                          !-----------------------------------------------------
                          !  Keep in or add to close set
                          !-----------------------------------------------------
                          IF ((gmm_state2d(i,j,jsub) .EQ.     &
     &                        ppm_gmm_param_far) .AND.          &
     &                        (ABS(fdta(i,j,jsub)) .LT. width)) THEN
                              gmm_state2d(i,j,jsub) =         &
     &                            ppm_gmm_param_close
#include "ppm_gmm_add_to_list.inc"
                          ENDIF
                      ENDIF
                  ENDIF
              ENDIF
              !-----------------------------------------------------------------
              !  Accept this point and remove it from the list
              !-----------------------------------------------------------------
              gmm_state2d(ii,jj,jsub) = ppm_gmm_param_accepted
              gmm_ipos(1,p) = gmm_ipos(1,npos0)
              gmm_ipos(2,p) = gmm_ipos(2,npos0)
              gmm_ipos(3,p) = gmm_ipos(3,npos0)
              npos0 = npos0 - 1
          ENDIF           ! TT .LE. TM
      ENDDO          ! p=1,npos
      npos = npos0
      !-------------------------------------------------------------------------
      !  Update ghost layers for both fdta AND gmm_state2d
      !-------------------------------------------------------------------------
      CALL ppm_map_field_push(gmm_topoid,gmm_meshid,fdta,info)
      IF (info .NE. ppm_param_success) THEN
          info = ppm_error_error
          CALL ppm_error(ppm_err_sub_failed,'ppm_gmm_march_fwd',  &
     &        'pushing field data failed',__LINE__,info)
          GOTO 9999
      ENDIF
      CALL ppm_map_field_push(gmm_topoid,gmm_meshid,gmm_state2d,info)
      IF (info .NE. ppm_param_success) THEN
          info = ppm_error_error
          CALL ppm_error(ppm_err_sub_failed,'ppm_gmm_march_fwd',  &
     &        'pushing status data failed',__LINE__,info)
          GOTO 9999
      ENDIF
      CALL ppm_map_field_send(info)
      IF (info .NE. ppm_param_success) THEN
          info = ppm_error_error
          CALL ppm_error(ppm_err_sub_failed,'ppm_gmm_march_fwd',  &
     &        'sending ghosts failed',__LINE__,info)
          GOTO 9999
      ENDIF
      CALL ppm_map_field_pop(gmm_topoid,gmm_meshid,gmm_state2d,ghostsize,info)
      IF (info .NE. ppm_param_success) THEN
          info = ppm_error_error
          CALL ppm_error(ppm_err_sub_failed,'ppm_gmm_march_fwd',  &
     &        'popping status data failed',__LINE__,info)
          GOTO 9999
      ENDIF
      CALL ppm_map_field_pop(gmm_topoid,gmm_meshid,fdta,ghostsize,info)
      IF (info .NE. ppm_param_success) THEN
          info = ppm_error_error
          CALL ppm_error(ppm_err_sub_failed,'ppm_gmm_march_fwd',  &
     &        'popping field data failed',__LINE__,info)
          GOTO 9999
      ENDIF
#endif 
      !-------------------------------------------------------------------------
      !  Return 
      !-------------------------------------------------------------------------
 9999 CONTINUE
      CALL substop('ppm_gmm_march_fwd',t0,info)
      RETURN

#if    __DIM == __2D
#if    __KIND == __SINGLE_PRECISION
      END SUBROUTINE ppm_gmm_march_fwd_2ds
#elif  __KIND == __DOUBLE_PRECISION
      END SUBROUTINE ppm_gmm_march_fwd_2dd
#endif 
#elif  __DIM == __3D
#if    __KIND == __SINGLE_PRECISION
      END SUBROUTINE ppm_gmm_march_fwd_3ds
#elif  __KIND == __DOUBLE_PRECISION
      END SUBROUTINE ppm_gmm_march_fwd_3dd
#endif 
#endif
