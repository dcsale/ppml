      !-------------------------------------------------------------------------
      !  Subroutine   :                   ppm_gmm_march
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
      SUBROUTINE ppm_gmm_march_2ds(width,order,fdata,rhscst,MaxIter,info, &
     &    speed,udata,chi)
#elif  __KIND == __DOUBLE_PRECISION
      SUBROUTINE ppm_gmm_march_2dd(width,order,fdata,rhscst,MaxIter,info, &
     &    speed,udata,chi)
#endif
#elif  __DIM == __3D
#if    __KIND == __SINGLE_PRECISION
      SUBROUTINE ppm_gmm_march_3ds(width,order,fdata,rhscst,MaxIter,info, &
     &    speed,udata,chi)
#elif  __KIND == __DOUBLE_PRECISION
      SUBROUTINE ppm_gmm_march_3dd(width,order,fdata,rhscst,MaxIter,info, &
     &    speed,udata,chi)
#endif
#endif
      !!! This routine performs the group marching to
      !!! extend the infprmation from the interface to a
      !!! narrow band. This is used for level set
      !!! reinitialization or function extension. The
      !!! interface is assumed to be at the zero level.
      !!! Use postprocess shifting if this is not the case.
      !!!
      !!! [NOTE]
      !!! This routine uses the o(N) group marching method
      !!! of Kim, rather than o(NlogN) fast marching. If
      !!! combines this with the high order equations by
      !!! Chopp to get 2nd or 3rd order accurate
      !!! extrapolations in linear time.
      !!! speed should maybe be passed as speed**2. This
      !!! would save the multiplication in gmm_slvupwind,
      !!! but would require SQRTs in the initial
      !!! conversion from distance to travel time.
      !!! When recomputing the close neighbors one could
      !!! detect local minima and directly remove (i.e.
      !!! accept) them to avoid unnecessary double
      !!! computations.
      !!! To be more memory efficient, we could maybe
      !!! allocate gmm_state only for the sub we are
      !!! currently working on instead of for the full local
      !!! field. This will however only give cashback if
      !!! there are lots of subs per proc, which is currently
      !!! not the case for mesh topologies.
      !!! 
      !!! === References ===
      !!!
      !!! D.L. Chopp. Some improvements on the fast marching method. SIAM J. 
      !!! Sci. Comput. 23(1):230-244, 2001.
      !!! S. Kim. An O(N) level set method for Eikonal equations. SIAM J. Sci. 
      !!! Comput. 22(6):2178-2193, 2001.
      !-------------------------------------------------------------------------
      !  Modules
      !-------------------------------------------------------------------------
      USE ppm_module_data
      USE ppm_module_numerics_data
      USE ppm_module_data_gmm
      USE ppm_module_gmm_march_fwd
      USE ppm_module_gmm_march_bkwd
      USE ppm_module_gmm_extend_fwd
      USE ppm_module_gmm_extend_bkwd
      USE ppm_module_substart
      USE ppm_module_substop
      USE ppm_module_error
      USE ppm_module_write
      USE ppm_module_alloc
      USE ppm_module_map_field
      USE ppm_module_map_field_ghost
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
      REAL(MK), DIMENSION(:,:,:)     , INTENT(IN), OPTIONAL :: speed
#elif __DIM == __3D
      REAL(MK), DIMENSION(:,:,:,:)   , INTENT(IN), OPTIONAL :: speed
#endif
      !!! rank 4 (3d) or rank 3 (2d) field of front speeds.
      !!! OPTIONAL to override rhscst.
#if   __DIM == __2D
      REAL(MK), DIMENSION(:,:,:)     , POINTER              :: fdata
#elif __DIM == __3D
      REAL(MK), DIMENSION(:,:,:,:)   , POINTER              :: fdata
#endif
      !!! Level data. Either rank 3 (for 2D scalar fields), or rank
      !!! 4 (for 3D scalar fields). Indices: (i,k,[k],isub). On input:
      !!! signed distance function after initialization as returned by
      !!! ppm_gmm_kickoff. Uninitialized values must be set to HUGE.
      !!! On output: approximation of the extended signed distance
      !!! function in the whole band. Points far from the interface will have
      !!! the value HUGE.
#if   __DIM == __2D
      REAL(MK), DIMENSION(:,:,:)     , POINTER, OPTIONAL    :: udata
#elif __DIM == __3D
      REAL(MK), DIMENSION(:,:,:,:)   , POINTER, OPTIONAL    :: udata
#endif
      !!! Function data. OPTIONAL. If present, udata will be extended in the
      !!! whole band such that its gradient is orthogonal to the gradient
      !!! of fdata. fdata in this case is a pure input and needs to be
      !!! available in the whole band. Needs to be a scalar field.
#if   __DIM == __2D
      REAL(MK), DIMENSION(:,:,:,:)   , INTENT(IN), OPTIONAL :: chi
#elif __DIM == __3D
      REAL(MK), DIMENSION(:,:,:,:,:) , INTENT(IN), OPTIONAL :: chi
#endif
      !!! rank 5 (3d) or rank 4 (2d) field specifying the positions of the
      !!! grid nodes. 1st index: 1..ppm_dim, then i,j,[k],isub.
      !!! OPTIONAL. Uniform grid is assumed if absent. Ghostlayers of size
      !!! >=1 must be pre-filled.
      REAL(MK)                       , INTENT(IN   )        :: width
      !!! Width of the narrow band to be produced on each side of the interface.
      REAL(MK)                       , INTENT(IN   )        :: rhscst
      !!! constant value for the right hand side of grad u * grad f = c. If
      !!! speed is present, this argument will be ignored.
      INTEGER                        , INTENT(IN   )        :: order
      !!! Order of the method to be used. One of
      !!!
      !!! *ppm_param_order_1
      !!! *ppm_param_order_2
      !!! *ppm_param_order_3
      INTEGER                        , INTENT(IN   )        :: MaxIter
      !!! argument specifying the maximum number of allowed iterations.
      !!! This can be useful since a cyclic dependency in the GMM algorithms
      !!! could cause infinite loops.
      INTEGER                        , INTENT(  OUT)        :: info
      !!! Return status, 0 upon success
      !-------------------------------------------------------------------------
      !  Local variables
      !-------------------------------------------------------------------------
      INTEGER                          :: i,j,k,p,iopt,xhi,yhi,zhi
      INTEGER                          :: npos,npos0,jsub,isub
      INTEGER                          :: nposg
      INTEGER                          :: nord
      INTEGER                          :: Marchit
      INTEGER, DIMENSION(3)            :: ghostsize
      INTEGER, DIMENSION(4)            :: ldl,ldu
      REAL(MK),DIMENSION(2)            :: coefs
      REAL(MK)                         :: t0,deltaT,sqrtdiminv,smin
      REAL(MK)                         :: TM,dx,dy,dz,dxinv,dyinv,dzinv
      REAL(MK)                         :: hsave,mindx
      REAL(MK)                         :: ainv,big
      TYPE(ppm_t_topo)     , POINTER   :: topo
      TYPE(ppm_t_equi_mesh), POINTER   :: mesh
#ifdef __MPI
      REAL(MK)                         :: allsmin
#endif
      LOGICAL                          :: lok
      CHARACTER(LEN=ppm_char)          :: cbuf
#if   __DIM == __2D
      REAL(MK), DIMENSION(:,:,:    ), POINTER :: fdta
      REAL(MK), DIMENSION(:,:,:    ), POINTER :: dta
#elif __DIM == __3D
      REAL(MK), DIMENSION(:,:,:,:  ), POINTER :: fdta
      REAL(MK), DIMENSION(:,:,:,:  ), POINTER :: dta
#endif
      !-------------------------------------------------------------------------
      !  Externals
      !-------------------------------------------------------------------------

      !-------------------------------------------------------------------------
      !  Initialise
      !-------------------------------------------------------------------------
      CALL substart('ppm_gmm_march',t0,info)
      big        = HUGE(big)
      hsave      = 0.9_MK*big
#if   __DIM == __2D
      sqrtdiminv = 1.0_MK/SQRT(2.0_MK)
#else
      sqrtdiminv = 1.0_MK/SQRT(3.0_MK)
#endif
      topo => ppm_topo(gmm_topoid)%t
      mesh => topo%mesh(gmm_meshid)
      !-------------------------------------------------------------------------
      !  Check arguments
      !-------------------------------------------------------------------------
      IF (ppm_debug .GT. 0) THEN
          IF (.NOT. ppm_initialized) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_ppm_noinit,'ppm_gmm_march',  &
     &            'Please call ppm_init first!',__LINE__,info)
              GOTO 9999
          ENDIF
          IF (width .LT. 0.0_MK) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_ppm_noinit,'ppm_gmm_march',  &
     &            'width must be positive!',__LINE__,info)
              GOTO 9999
          ENDIF
          IF ((order.NE.ppm_param_order_1).AND.(order.NE.ppm_param_order_2)  &
     &         .AND.(order.NE.ppm_param_order_3)) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_ppm_noinit,'ppm_gmm_march',  &
     &            'order must be 1, 2, or 3!',__LINE__,info)
              GOTO 9999
          ENDIF
#if   __DIM == __3D
          IF (SIZE(fdata,4) .LT. topo%nsublist) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_gmm_march',  &
     &            'field data for some subs is missing',__LINE__,info)
              GOTO 9999
          ENDIF
          IF (UBOUND(fdata,1) .LT. maxxhi) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_gmm_march',  &
     &            'x dimension of field data does not match mesh',__LINE__,info)
              GOTO 9999
          ENDIF
          IF (UBOUND(fdata,2) .LT. maxyhi) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_gmm_march',  &
     &            'y dimension of field data does not match mesh',__LINE__,info)
              GOTO 9999
          ENDIF
          IF (UBOUND(fdata,3) .LT. maxzhi) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_gmm_march',  &
     &            'z dimension of field data does not match mesh',__LINE__,info)
              GOTO 9999
          ENDIF
#elif __DIM == __2D
          IF (SIZE(fdata,3) .LT. topo%nsublist) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_gmm_march',  &
     &            'field data for some subs is missing',__LINE__,info)
              GOTO 9999
          ENDIF
          IF (UBOUND(fdata,1) .LT. maxxhi) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_gmm_march',  &
     &            'x dimension of field data does not match mesh',__LINE__,info)
              GOTO 9999
          ENDIF
          IF (UBOUND(fdata,2) .LT. maxyhi) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_gmm_march',  &
     &            'y dimension of field data does not match mesh',__LINE__,info)
              GOTO 9999
          ENDIF
#endif
      ENDIF          ! ppm_debug for argument check
      !-------------------------------------------------------------------------
      !  Determine ghostsize needed. Only change if larger than what it
      !  already is. Do not shrink as this would freak out the client...
      !-------------------------------------------------------------------------
#if   __DIM == __3D
      i = MIN(1-LBOUND(fdata,1),1-LBOUND(fdata,2),1-LBOUND(fdata,3))
#else
      i = MIN(1-LBOUND(fdata,1),1-LBOUND(fdata,2))
#endif
      IF (order .EQ. ppm_param_order_1) THEN
          ghostsize = MAX(1,i)
          nord = 1
       !  IF (ppm_debug .GT. 0) THEN
              WRITE(cbuf,'(A)') 'Using first order upwinding'
              CALL ppm_write(ppm_rank,'ppm_gmm_march',cbuf,info)
       !  ENDIF
      ELSEIF (order .EQ. ppm_param_order_2) THEN
          ghostsize = MAX(2,i)
          nord = 2
       !  IF (ppm_debug .GT. 0) THEN
              WRITE(cbuf,'(A)') 'Using second order upwinding'
              CALL ppm_write(ppm_rank,'ppm_gmm_march',cbuf,info)
       !  ENDIF
      ELSE
          ghostsize = MAX(3,i)
          nord = 3
       !  IF (ppm_debug .GT. 0) THEN
              WRITE(cbuf,'(A)') 'Using third order upwinding'
              CALL ppm_write(ppm_rank,'ppm_gmm_march',cbuf,info)
       !  ENDIF
      ENDIF
      !-------------------------------------------------------------------------
      !  Set the pointers
      !-------------------------------------------------------------------------
      IF (PRESENT(udata)) THEN
          !---------------------------------------------------------------------
          !  March udata and use fdata as level set
          !---------------------------------------------------------------------
          dta => udata
          fdta => fdata
      ELSE
          !---------------------------------------------------------------------
          !  March the level set itself
          !---------------------------------------------------------------------
          dta => fdata
          fdta => fdata
      ENDIF
      !-------------------------------------------------------------------------
      !  Find mesh spacing
      !-------------------------------------------------------------------------
      IF (ppm_kind .EQ. ppm_kind_single) THEN
          dx = (topo%max_physs(1)-topo%min_physs(1))/   &
     &        REAL(mesh%Nm(1)-1,ppm_kind_single)
          dy = (topo%max_physs(2)-topo%min_physs(2))/   &
     &        REAL(mesh%Nm(2)-1,ppm_kind_single)
          IF (ppm_dim .GT. 2) THEN
              dz = (topo%max_physs(3)-topo%min_physs(3))/ &
     &            REAL(mesh%Nm(3)-1,     &
     &            ppm_kind_single)
          ENDIF
      ELSE
          dx = (topo%max_physs(1)-topo%min_physs(1))/   &
     &        REAL(mesh%Nm(1)-1,ppm_kind_double)
          dy = (topo%max_physs(2)-topo%min_physs(2))/   &
     &        REAL(mesh%Nm(2)-1,ppm_kind_double)
          IF (ppm_dim .GT. 2) THEN
              dz = (topo%max_physs(3)-topo%min_physs(3))/ &
     &            REAL(mesh%Nm(3)-1,     &
     &            ppm_kind_double)
          ENDIF
      ENDIF
      dxinv = 1.0_MK/dx
      dyinv = 1.0_MK/dy
#if   __DIM == __3D
      dzinv = 1.0_MK/dz
      mindx = MIN(dx,dy,dz)
#else
      dzinv = 0.0_MK
      mindx = MIN(dx,dy)
#endif
      !-------------------------------------------------------------------------
      !  Initialize the ghost layers.
      !-------------------------------------------------------------------------
      CALL ppm_map_field_ghost_get(gmm_topoid,gmm_meshid,ghostsize,info)
      IF (info .NE. ppm_param_success) THEN
          info = ppm_error_error
          CALL ppm_error(ppm_err_sub_failed,'ppm_gmm_march',  &
     &        'ghost get mapping failed',__LINE__,info)
          GOTO 9999
      ENDIF
      CALL ppm_map_field_push(gmm_topoid,gmm_meshid,fdta,info)
      IF (info .NE. ppm_param_success) THEN
          info = ppm_error_error
          CALL ppm_error(ppm_err_sub_failed,'ppm_gmm_march',  &
     &        'pushing level data failed',__LINE__,info)
          GOTO 9999
      ENDIF
      IF (PRESENT(udata)) THEN
          CALL ppm_map_field_push(gmm_topoid,gmm_meshid,dta,info)
          IF (info .NE. ppm_param_success) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_sub_failed,'ppm_gmm_march',  &
     &        'pushing field data failed',__LINE__,info)
              GOTO 9999
          ENDIF
      ENDIF
      CALL ppm_map_field_send(info)
      IF (info .NE. ppm_param_success) THEN
          info = ppm_error_error
          CALL ppm_error(ppm_err_sub_failed,'ppm_gmm_march',  &
     &        'sending ghosts failed',__LINE__,info)
          GOTO 9999
      ENDIF
      IF (PRESENT(udata)) THEN
          CALL ppm_map_field_pop(gmm_topoid,gmm_meshid,dta,ghostsize,info)
          IF (info .NE. ppm_param_success) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_sub_failed,'ppm_gmm_march',  &
     &        'popping field data failed',__LINE__,info)
              GOTO 9999
          ENDIF
      ENDIF
      CALL ppm_map_field_pop(gmm_topoid,gmm_meshid,fdta,ghostsize,info)
      IF (info .NE. ppm_param_success) THEN
          info = ppm_error_error
          CALL ppm_error(ppm_err_sub_failed,'ppm_gmm_march',  &
     &        'popping level data failed',__LINE__,info)
          GOTO 9999
      ENDIF
      ! because it has not been mapped
      IF (.NOT.(PRESENT(udata))) dta => fdta
      !-------------------------------------------------------------------------
      !  Allocate and initialize status array
      !-------------------------------------------------------------------------
      iopt = ppm_param_alloc_grow
#if   __DIM == __3D
      ldl(1) = 1-ghostsize(1)
      ldl(2) = 1-ghostsize(2)
      ldl(3) = 1-ghostsize(3)
      ldl(4) = 1
      ldu(1) = maxxhi+ghostsize(1)
      ldu(2) = maxyhi+ghostsize(2)
      ldu(3) = maxzhi+ghostsize(3)
      ldu(4) = topo%nsublist
      CALL ppm_alloc(gmm_state3d,ldl,ldu,iopt,info)
#else
      ldl(1) = 1-ghostsize(1)
      ldl(2) = 1-ghostsize(2)
      ldl(3) = 1
      ldu(1) = maxxhi+ghostsize(1)
      ldu(2) = maxyhi+ghostsize(2)
      ldu(3) = topo%nsublist
      CALL ppm_alloc(gmm_state2d,ldl,ldu,iopt,info)
#endif
      IF (info .NE. 0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_gmm_march',      &
     &        'point status array GMM_STATE',__LINE__,info)
          GOTO 9999
      ENDIF
#if   __DIM == __3D
      gmm_state3d = ppm_gmm_param_far
      !-------------------------------------------------------------------------
      !  Convert distance to travel time and find minimum slowness (=
      !  max speed)
      !-------------------------------------------------------------------------
      DO jsub=1,topo%nsublist
          IF (PRESENT(speed)) THEN
              isub = topo%isublist(jsub)
              xhi = mesh%nnodes(1,isub)
              yhi = mesh%nnodes(2,isub)
              zhi = mesh%nnodes(3,isub)
              smin  = -big
              DO k=1,zhi
                  DO j=1,yhi
                      DO i=1,xhi
                          IF (speed(i,j,k,jsub) .GT. smin)    &
     &                        smin = speed(i,j,k,jsub)
                          IF (dta(i,j,k,jsub) .LT. hsave) THEN
                              dta(i,j,k,jsub)=dta(i,j,k,jsub)/speed(i,j,k,jsub)
                          ENDIF
                      ENDDO
                  ENDDO
              ENDDO
              smin = 1.0_MK/smin
#ifdef __MPI
              !-----------------------------------------------------------------
              !  Ensure consistent computations on all processors
              !-----------------------------------------------------------------
#if   __KIND == __SINGLE_PRECISION
              CALL MPI_Allreduce(smin,allsmin,1,MPI_REAL,MPI_MIN,   &
     &            ppm_comm,info)
#else
              CALL MPI_Allreduce(smin,allsmin,1,MPI_DOUBLE_PRECISION,MPI_MIN, &
     &            ppm_comm,info)
#endif
              IF (info .NE. 0) THEN
                  info = ppm_error_fatal
                  CALL ppm_error(ppm_err_mpi_fail,'ppm_gmm_march',      &
     &                'MPI_ALLREDUCE of smin',__LINE__,info)
                  GOTO 9999
              ENDIF
              smin = allsmin
#endif
          ELSE
              smin = 1.0_MK
          ENDIF
      ENDDO        ! jsub
      !-------------------------------------------------------------------------
      !  Initialize travel limit (see Kim:2001a)
      !-------------------------------------------------------------------------
      deltaT = 0.99_MK*sqrtdiminv*mindx*smin
      !IF (.NOT.PRESENT(udata)) deltaT = (10.0_MK**REAL(-nord+1,MK))*deltaT
      !-------------------------------------------------------------------------
      !  Make a list of all initialized points and find the minimum
      !  travel time on the boundary of the band.
      !-------------------------------------------------------------------------
      npos0 = 0
      DO jsub=1,topo%nsublist
          isub = topo%isublist(jsub)
          xhi  = mesh%nnodes(1,isub)
          yhi  = mesh%nnodes(2,isub)
          zhi  = mesh%nnodes(3,isub)
          TM   = big
          DO k=1,zhi
              DO j=1,yhi
                  DO i=1,xhi
                      coefs(1) = dta(i,j,k,jsub)
                      coefs(2) = fdta(i,j,k,jsub)
                      IF ((coefs(1).LT.hsave).AND.(gmm_state3d(i,j,k,jsub).EQ.&
     &                    ppm_gmm_param_far)) THEN
                          IF ((dta(i-1,j,k,jsub) .GT. hsave) .OR.   &
     &                        (dta(i+1,j,k,jsub) .GT. hsave) .OR.   &
     &                        (dta(i,j-1,k,jsub) .GT. hsave) .OR.   &
     &                        (dta(i,j+1,k,jsub) .GT. hsave) .OR.   &
     &                        (dta(i,j,k-1,jsub) .GT. hsave) .OR.   &
     &                        (dta(i,j,k+1,jsub) .GT. hsave)) THEN
                              !-------------------------------------------------
                              !  Add points on the surface to close set and
                              !  set TM to min distance to interface
                              !-------------------------------------------------
                              IF (ABS(coefs(2)) .LT. TM) TM = ABS(coefs(2))
                              gmm_state3d(i,j,k,jsub) = ppm_gmm_param_close
#include "ppm_gmm_add_to_list.inc"
                          ELSE
                              gmm_state3d(i,j,k,jsub) = ppm_gmm_param_accepted
                          ENDIF
                      ENDIF
                  ENDDO       ! mesh cell i
              ENDDO       ! mesh cell j
          ENDDO       ! mesh cell k
      ENDDO       ! jsub
      npos = npos0
      !-------------------------------------------------------------------------
      !  Initialize ghost layers for state array
      !-------------------------------------------------------------------------
      CALL ppm_map_field_push(gmm_topoid,gmm_meshid,gmm_state3d,info)
      IF (info .NE. ppm_param_success) THEN
          info = ppm_error_error
          CALL ppm_error(ppm_err_sub_failed,'ppm_gmm_march',  &
     &        'pushing status data failed',__LINE__,info)
          GOTO 9999
      ENDIF
      CALL ppm_map_field_send(info)
      IF (info .NE. ppm_param_success) THEN
          info = ppm_error_error
          CALL ppm_error(ppm_err_sub_failed,'ppm_gmm_march',  &
     &        'sending status ghosts failed',__LINE__,info)
          GOTO 9999
      ENDIF
      CALL ppm_map_field_pop(gmm_topoid,gmm_meshid,gmm_state3d,ghostsize,info)
      IF (info .NE. ppm_param_success) THEN
          info = ppm_error_error
          CALL ppm_error(ppm_err_sub_failed,'ppm_gmm_march',  &
     &        'popping status data failed',__LINE__,info)
          GOTO 9999
      ENDIF

      !-------------------------------------------------------------------------
      !  Initialize the global gamma size
      !-------------------------------------------------------------------------
#ifdef __MPI
      CALL MPI_AllReduce(npos,nposg,1,MPI_INTEGER,MPI_SUM,ppm_comm,info)
#else
      nposg = npos
#endif
      !-------------------------------------------------------------------------
      !  March
      !-------------------------------------------------------------------------
      Marchit = 0
      DO WHILE ((nposg .GT. 0) .AND. (Marchit .LT. MaxIter))
          Marchit = Marchit + 1
          TM = TM + deltaT
          IF (ppm_debug .GT. 0) THEN
              WRITE(cbuf,'(A,I4,A,I9)') 'Iteration ',Marchit,   &
     &            ': Number of points in list: ',npos
              CALL ppm_write(ppm_rank,'ppm_gmm_march',cbuf,info)
              IF (npos.LE.5) THEN
                   DO i=1,npos
                       WRITE(cbuf,'(A,I9,I9,I9)') 'Point ', gmm_ipos(1,i), gmm_ipos(2,i), gmm_ipos(3,i)
                       CALL ppm_write(ppm_rank,'ppm_gmm_march',cbuf,info)
                   END DO
              END IF
          ENDIF
          !---------------------------------------------------------------------
          !  Hard-core debuging
          !---------------------------------------------------------------------
!         WRITE(cbuf,'(A,I4.4,A)') 'state_',Marchit,'.out'
!         OPEN(50,FILE=cbuf,STATUS='REPLACE',ACTION='WRITE')
!         WRITE(cbuf,'(A,I4.4,A)') 'level_',Marchit,'.out'
!         OPEN(40,FILE=cbuf,STATUS='REPLACE',ACTION='WRITE')
!         DO k=1,zhi
!             DO j=1,yhi
!                 DO i=1,xhi
!                     IF (ABS(fdta(i,j,k,1)) .LT. hsave) THEN
!                         WRITE(40,*) fdta(i,j,k,1)
!                     ELSE
!                         WRITE(40,*) 0.0_MK
!                     ENDIF
!                     WRITE(50,*) gmm_state3d(i,j,k,1)
!                 ENDDO
!             ENDDO
!         ENDDO
!         CLOSE(40)
!         CLOSE(50)
          !---------------------------------------------------------------------
          !  Reverse order: recompute neighbors of points in ggm_ipos to
          !  fix stability.
          !---------------------------------------------------------------------
          IF (PRESENT(udata)) THEN
              DO i=1,nord
                  IF (PRESENT(speed)) THEN
                      CALL ppm_gmm_extend_bkwd(fdta,dta,width,order,npos,TM,  &
     &                    rhscst,dxinv,dyinv,dzinv,ghostsize,info,speed)
                  ELSE
                      CALL ppm_gmm_extend_bkwd(fdta,dta,width,order,npos,TM,  &
     &                    rhscst,dxinv,dyinv,dzinv,ghostsize,info)
                  ENDIF
                  IF (info .NE. 0) THEN
                      info = ppm_error_error
                      CALL ppm_error(ppm_err_sub_failed,'ppm_gmm_march',      &
     &                    'Backward extension failed.',__LINE__,info)
                      GOTO 8888
                  ENDIF
              ENDDO
          ELSE
              IF (PRESENT(speed)) THEN
                  CALL ppm_gmm_march_bkwd(fdta,width,order,npos,TM,  &
     &                rhscst,dxinv,dyinv,dzinv,ghostsize,info,speed)
              ELSE
                  CALL ppm_gmm_march_bkwd(fdta,width,order,npos,TM,  &
     &                rhscst,dxinv,dyinv,dzinv,ghostsize,info)
              ENDIF
              IF (info .NE. 0) THEN
                  info = ppm_error_error
                  CALL ppm_error(ppm_err_sub_failed,'ppm_gmm_march',      &
     &                'Backward marching failed.',__LINE__,info)
                  GOTO 8888
              ENDIF
          ENDIF
          !---------------------------------------------------------------------
          !  Forward order: update neighbors of points in ggm_ipos and
          !  advance front
          !---------------------------------------------------------------------
          IF (PRESENT(udata)) THEN
              IF (PRESENT(speed)) THEN
                  CALL ppm_gmm_extend_fwd(fdta,dta,width,order,npos,TM,  &
     &                rhscst,dxinv,dyinv,dzinv,ghostsize,info,speed)
              ELSE
                  CALL ppm_gmm_extend_fwd(fdta,dta,width,order,npos,TM,  &
     &                rhscst,dxinv,dyinv,dzinv,ghostsize,info)
              ENDIF
          ELSE
              IF (PRESENT(speed)) THEN
                  CALL ppm_gmm_march_fwd(fdta,width,order,npos,TM,  &
     &                rhscst,dxinv,dyinv,dzinv,ghostsize,info,speed)
              ELSE
                  CALL ppm_gmm_march_fwd(fdta,width,order,npos,TM,  &
     &                rhscst,dxinv,dyinv,dzinv,ghostsize,info)
              ENDIF
          ENDIF
          IF (info .NE. 0) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_sub_failed,'ppm_gmm_march',      &
     &            'Forward marching failed.',__LINE__,info)
              GOTO 8888
          ENDIF
          !---------------------------------------------------------------------
          !  Update the global gamma size
          !---------------------------------------------------------------------
#ifdef __MPI
          CALL MPI_AllReduce(npos,nposg,1,MPI_INTEGER,MPI_SUM,ppm_comm,info)
#else
          nposg = npos
#endif
      ENDDO          ! while nposg.GT.0
#elif __DIM == __2D
      gmm_state2d = ppm_gmm_param_far
      !-------------------------------------------------------------------------
      !  Convert distance to travel time and find minimum slowness (=
      !  max speed)
      !-------------------------------------------------------------------------
      DO jsub=1,topo%nsublist
          IF (PRESENT(speed)) THEN
              isub = topo%isublist(jsub)
              xhi = mesh%nnodes(1,isub)
              yhi = mesh%nnodes(2,isub)
              smin  = -big
              DO j=1,yhi
                  DO i=1,xhi
                      IF (speed(i,j,jsub) .GT. smin)    &
     &                    smin = speed(i,j,jsub)
                      IF (dta(i,j,jsub) .LT. hsave) THEN
                          dta(i,j,jsub) =               &
     &                        dta(i,j,jsub)/speed(i,j,jsub)
                      ENDIF
                  ENDDO
              ENDDO
              smin = 1.0_MK/smin
#ifdef __MPI
              !-----------------------------------------------------------------
              !  Ensure consistent computations on all processors
              !-----------------------------------------------------------------
#if   __KIND == __SINGLE_PRECISION
              CALL MPI_Allreduce(smin,allsmin,1,MPI_REAL,MPI_MIN,   &
     &            ppm_comm,info)
#else
              CALL MPI_Allreduce(smin,allsmin,1,MPI_DOUBLE_PRECISION,MPI_MIN, &
     &            ppm_comm,info)
#endif
              IF (info .NE. 0) THEN
                  info = ppm_error_fatal
                  CALL ppm_error(ppm_err_mpi_fail,'ppm_gmm_march',      &
     &                'MPI_ALLREDUCE of smin',__LINE__,info)
                  GOTO 9999
              ENDIF
              smin = allsmin
#endif
          ELSE
              smin = 1.0_MK
          ENDIF
      ENDDO        ! jsub
      !-------------------------------------------------------------------------
      !  Initialize travel limit (see Kim:2001a)
      !-------------------------------------------------------------------------
      deltaT = 0.99_MK*sqrtdiminv*mindx*smin
      !IF (.NOT.PRESENT(udata)) deltaT = (10.0_MK**REAL(-nord+1,MK))*deltaT
      !-------------------------------------------------------------------------
      !  Make a list of all initialized points and find the minimum
      !  travel time on the boundary of the band.
      !-------------------------------------------------------------------------
      npos0 = 0
      DO jsub=1,topo%nsublist
          isub = topo%isublist(jsub)
          xhi  = mesh%nnodes(1,isub)
          yhi  = mesh%nnodes(2,isub)
          TM   = big
          DO j=1,yhi
              DO i=1,xhi
                  coefs(1) = dta(i,j,jsub)
                  coefs(2) = fdta(i,j,jsub)
                  IF ((coefs(1).LT.hsave).AND.(gmm_state2d(i,j,jsub).EQ.&
     &                ppm_gmm_param_far)) THEN
                      IF ((dta(i-1,j,jsub) .GT. hsave) .OR.   &
     &                    (dta(i+1,j,jsub) .GT. hsave) .OR.   &
     &                    (dta(i,j-1,jsub) .GT. hsave) .OR.   &
     &                    (dta(i,j+1,jsub) .GT. hsave)) THEN
                          !-----------------------------------------------------
                          !  Add points on the surface to close set and
                          !  set TM to min distance to interface
                          !-----------------------------------------------------
                          IF (ABS(coefs(2)) .LT. TM) TM = ABS(coefs(2))
                          gmm_state2d(i,j,jsub) = ppm_gmm_param_close
#include "ppm_gmm_add_to_list.inc"
                      ELSE
                          gmm_state2d(i,j,jsub) = ppm_gmm_param_accepted
                      ENDIF
                  ENDIF
              ENDDO       ! mesh cell i
          ENDDO       ! mesh cell j
      ENDDO       ! jsub
      npos = npos0
      !-------------------------------------------------------------------------
      !  Initialize ghost layers for state array
      !-------------------------------------------------------------------------
      CALL ppm_map_field_push(gmm_topoid,gmm_meshid,gmm_state2d,info)
      IF (info .NE. ppm_param_success) THEN
          info = ppm_error_error
          CALL ppm_error(ppm_err_sub_failed,'ppm_gmm_march',  &
     &        'pushing status data failed',__LINE__,info)
          GOTO 9999
      ENDIF
      CALL ppm_map_field_send(info)
      IF (info .NE. ppm_param_success) THEN
          info = ppm_error_error
          CALL ppm_error(ppm_err_sub_failed,'ppm_gmm_march',  &
     &        'sending status ghosts failed',__LINE__,info)
          GOTO 9999
      ENDIF
      CALL ppm_map_field_pop(gmm_topoid,gmm_meshid,gmm_state2d,ghostsize,info)
      IF (info .NE. ppm_param_success) THEN
          info = ppm_error_error
          CALL ppm_error(ppm_err_sub_failed,'ppm_gmm_march',  &
     &        'popping status data failed',__LINE__,info)
          GOTO 9999
      ENDIF

      !-------------------------------------------------------------------------
      !  Initialize the global gamma size
      !-------------------------------------------------------------------------
#ifdef __MPI
      CALL MPI_AllReduce(npos,nposg,1,MPI_INTEGER,MPI_SUM,ppm_comm,info)
#else
      nposg = npos
#endif
      !-------------------------------------------------------------------------
      !  March
      !-------------------------------------------------------------------------
      Marchit = 0
      DO WHILE ((nposg .GT. 0) .AND. (Marchit .LT. MaxIter))
          Marchit = Marchit + 1
          TM = TM + deltaT
          IF (ppm_debug .GT. 0) THEN
              WRITE(cbuf,'(A,I4,A,I9)') 'Iteration ',Marchit,   &
     &            ': Number of points in list: ',npos
              CALL ppm_write(ppm_rank,'ppm_gmm_march',cbuf,info)
          ENDIF
          !---------------------------------------------------------------------
          !  Reverse order: recompute neighbors of points in ggm_ipos to
          !  fix stability
          !---------------------------------------------------------------------
          IF (PRESENT(udata)) THEN
              DO i=1,nord
                  IF (PRESENT(speed)) THEN
                      IF (PRESENT(chi)) THEN
                          CALL ppm_gmm_extend_bkwd(fdta,dta,width,order,npos,&
     &                        TM,rhscst,dxinv,dyinv,dzinv,ghostsize,info,    &
     &                        speed,chi)
                      ELSE
                          CALL ppm_gmm_extend_bkwd(fdta,dta,width,order,npos,&
     &                        TM,rhscst,dxinv,dyinv,dzinv,ghostsize,info,speed)
                      ENDIF
                  ELSE
                      IF (PRESENT(chi)) THEN
                          CALL ppm_gmm_extend_bkwd(fdta,dta,width,order,npos, &
     &                        TM,rhscst,dxinv,dyinv,dzinv,ghostsize,info,     &
     &                        chi=chi)
                      ELSE
                          CALL ppm_gmm_extend_bkwd(fdta,dta,width,order,npos, &
     &                        TM,rhscst,dxinv,dyinv,dzinv,ghostsize,info)
                      ENDIF
                  ENDIF
                  IF (info .NE. 0) THEN
                      info = ppm_error_error
                      CALL ppm_error(ppm_err_sub_failed,'ppm_gmm_march',      &
     &                    'Backward extension failed.',__LINE__,info)
                      GOTO 8888
                  ENDIF
              ENDDO
          ELSE
              IF (PRESENT(speed)) THEN
                  IF (PRESENT(chi)) THEN
                      CALL ppm_gmm_march_bkwd(fdta,width,order,npos,TM,  &
     &                    rhscst,dxinv,dyinv,dzinv,ghostsize,info,speed,chi)
                  ELSE
                      CALL ppm_gmm_march_bkwd(fdta,width,order,npos,TM,  &
     &                    rhscst,dxinv,dyinv,dzinv,ghostsize,info,speed)
                  ENDIF
              ELSE
                  IF (PRESENT(chi)) THEN
                      CALL ppm_gmm_march_bkwd(fdta,width,order,npos,TM,  &
     &                    rhscst,dxinv,dyinv,dzinv,ghostsize,info,chi=chi)
                  ELSE
                      CALL ppm_gmm_march_bkwd(fdta,width,order,npos,TM,  &
     &                    rhscst,dxinv,dyinv,dzinv,ghostsize,info)
                  ENDIF
              ENDIF
          ENDIF
          IF (info .NE. 0) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_sub_failed,'ppm_gmm_march',      &
     &            'Backward marching failed.',__LINE__,info)
              GOTO 8888
          ENDIF
          !---------------------------------------------------------------------
          !  Forward order: update neighbors of points in ggm_ipos and
          !  advance front
          !---------------------------------------------------------------------
          IF (PRESENT(udata)) THEN
              IF (PRESENT(speed)) THEN
                  IF (PRESENT(chi)) THEN
                      CALL ppm_gmm_extend_fwd(fdta,dta,width,order,npos,TM,  &
     &                    rhscst,dxinv,dyinv,dzinv,ghostsize,info,speed,chi)
                  ELSE
                      CALL ppm_gmm_extend_fwd(fdta,dta,width,order,npos,TM,  &
     &                    rhscst,dxinv,dyinv,dzinv,ghostsize,info,speed)
                  ENDIF
              ELSE
                  IF (PRESENT(chi)) THEN
                      CALL ppm_gmm_extend_fwd(fdta,dta,width,order,npos,TM,  &
     &                    rhscst,dxinv,dyinv,dzinv,ghostsize,info,chi=chi)
                  ELSE
                      CALL ppm_gmm_extend_fwd(fdta,dta,width,order,npos,TM,  &
     &                    rhscst,dxinv,dyinv,dzinv,ghostsize,info)
                  ENDIF
              ENDIF
          ELSE
              IF (PRESENT(speed)) THEN
                  IF (PRESENT(chi)) THEN
                      CALL ppm_gmm_march_fwd(fdta,width,order,npos,TM,  &
     &                    rhscst,dxinv,dyinv,dzinv,ghostsize,info,speed,chi)
                  ELSE
                      CALL ppm_gmm_march_fwd(fdta,width,order,npos,TM,  &
     &                    rhscst,dxinv,dyinv,dzinv,ghostsize,info,speed)
                  ENDIF
              ELSE
                  IF (PRESENT(chi)) THEN
                      CALL ppm_gmm_march_fwd(fdta,width,order,npos,TM,  &
     &                    rhscst,dxinv,dyinv,dzinv,ghostsize,info,chi=chi)
                  ELSE
                      CALL ppm_gmm_march_fwd(fdta,width,order,npos,TM,  &
     &                    rhscst,dxinv,dyinv,dzinv,ghostsize,info)
                  ENDIF
              ENDIF
          ENDIF
          IF (info .NE. 0) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_sub_failed,'ppm_gmm_march',      &
     &            'Forward marching failed.',__LINE__,info)
              GOTO 8888
          ENDIF
          !---------------------------------------------------------------------
          !  Update the global gamma size
          !---------------------------------------------------------------------
#ifdef __MPI
          CALL MPI_AllReduce(npos,nposg,1,MPI_INTEGER,MPI_SUM,ppm_comm,info)
#else
          nposg = npos
#endif
      ENDDO          ! while nposg.GT.0
#endif
      !-------------------------------------------------------------------------
      !  Deallocate status array
      !-------------------------------------------------------------------------
 8888 iopt = ppm_param_dealloc
#if   __DIM == __3D
      CALL ppm_alloc(gmm_state3d,ldu,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_error
          CALL ppm_error(ppm_err_dealloc,'ppm_gmm_march',      &
     &        'point status array GMM_STATE3D',__LINE__,info)
      ENDIF
#elif __DIM == __2D
      CALL ppm_alloc(gmm_state2d,ldu,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_error
          CALL ppm_error(ppm_err_dealloc,'ppm_gmm_march',      &
     &        'point status array GMM_STATE2D',__LINE__,info)
      ENDIF
#endif
      !-------------------------------------------------------------------------
      !  In case one of the mappings reallocated the data, we might want to
      !  pass back out the new data and save it from destruction upon
      !  subroutine termination.
      !-------------------------------------------------------------------------
      fdata => fdta
      IF (PRESENT(udata)) udata => dta
      NULLIFY(dta)
      NULLIFY(fdta)
      !-------------------------------------------------------------------------
      !  Return
      !-------------------------------------------------------------------------
 9999 CONTINUE
      CALL substop('ppm_gmm_march',t0,info)
      RETURN

#if    __DIM == __2D
#if    __KIND == __SINGLE_PRECISION
      END SUBROUTINE ppm_gmm_march_2ds
#elif  __KIND == __DOUBLE_PRECISION
      END SUBROUTINE ppm_gmm_march_2dd
#endif
#elif  __DIM == __3D
#if    __KIND == __SINGLE_PRECISION
      END SUBROUTINE ppm_gmm_march_3ds
#elif  __KIND == __DOUBLE_PRECISION
      END SUBROUTINE ppm_gmm_march_3dd
#endif
#endif
