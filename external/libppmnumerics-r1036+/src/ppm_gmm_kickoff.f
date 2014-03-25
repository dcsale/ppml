      !-------------------------------------------------------------------------
      !  Subroutine   :                 ppm_gmm_kickoff
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
      SUBROUTINE ppm_gmm_kickoff_2ds(fdata,tol,thresh,info,       &
     &    npts,ipts,closest,chi)
#elif  __KIND == __DOUBLE_PRECISION
      SUBROUTINE ppm_gmm_kickoff_2dd(fdata,tol,thresh,info,       &
     &    npts,ipts,closest,chi)
#endif 
#elif  __DIM == __3D
#if    __KIND == __SINGLE_PRECISION
      SUBROUTINE ppm_gmm_kickoff_3ds(fdata,tol,thresh,info,        &
     &    npts,ipts,closest,chi)
#elif  __KIND == __DOUBLE_PRECISION
      SUBROUTINE ppm_gmm_kickoff_3dd(fdata,tol,thresh,info,        &
     &    npts,ipts,closest,chi)
#endif 
#endif
      !!! This routine starts the gmm by computing second order accurate
      !!! approximations of the distance to the interface near it, i.e. on
      !!! neighboring grid points. Every sign reversal is interpreted as the
      !!! location of an interface! To use interfaces represented by a non-zero
      !!! level set, apply shifting as pre-processing.
      !!!
      !!! [NOTE]
      !!! ======================================================================
      !!! This routine creates second-order accurate data near the interface
      !!! and is needed to kick off a reinitialize a level set or extend a
      !!! function to a narrow band.
      !!!
      !!! The order of the method is limited by the order of the finite
      !!! differences used in computing the rhs of the interpolation system.
      !!! Use higher order FD to get higher order initialization.
      !!!
      !!! Tests have shown that storing shifted indices (i.e. ip1 = i+1) is
      !!! faster than using i+1 in the array index directly. We thus use
      !!! this technique here.
      !!!
      !!! Maybe we shoud actually allocate ipos and copy the stuff so it will
      !!! survive ppm_gmm_finalize?? This will be easy: just change the pointer
      !!! assignment at the end of the routine to a physical copy operation.
      !!! ======================================================================
      !!!
      !!! === References ====
      !!!
      !!! D.L. Chopp. Some improvements on the fast marching method. SIAM J.
      !!! Sci. Comput. 23(1):230-244, 2001.
      !-------------------------------------------------------------------------
      !  Includes
      !-------------------------------------------------------------------------
#include "ppm_define.h"
      !-------------------------------------------------------------------------
      !  Modules 
      !-------------------------------------------------------------------------
      USE ppm_module_data
      
      USE ppm_module_data_gmm
      USE ppm_module_substart
      USE ppm_module_substop
      USE ppm_module_error
      USE ppm_module_alloc
      USE ppm_module_typedef
      USE ppm_module_write
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
      REAL(MK), DIMENSION(:,:,:)       , POINTER              :: fdata
#elif __DIM == __3D
      REAL(MK), DIMENSION(:,:,:,:)     , POINTER              :: fdata
#endif
      !!! Level set data. Either rank 3 (for 2D scalar fields), or rank
      !!! 4 (for 3D scalar fields). Indices: (i,j,[k],isub). On input: old
      !!! level function values. The interface is at level zero.
      !!! A ghostsize of 1 is needed on all sides which must be filled with the
      !!! old level function value on input!! On output: 2nd order approximation
      !!! of the signed distance function using the interpolation method of
      !!! Chopp Points far from the interface will have the value HUGE.
#if   __DIM == __2D
      REAL(MK), DIMENSION(:,:,:,:)     , INTENT(IN), OPTIONAL :: chi
#elif __DIM == __3D
      REAL(MK), DIMENSION(:,:,:,:,:)   , INTENT(IN), OPTIONAL :: chi
#endif
      !!! rank 5 (3d) or rank 4 (2d) field specifying the positions of the grid
      !!! nodes. 1st index: 1..ppm_dim, then i,j,[k],isub. OPTIONAL. Uniform
      !!! grid is assumed if absent. Ghostlayers of size >=1 must be pre-filled.

      REAL(MK)                         , INTENT(IN   )        :: tol
      !!! Relative tolerance for the determined distance to the interface.
      !!! 1E-3 is a good choice. The tolerance is in multiples of grid spacings.
      REAL(MK)                         , INTENT(IN   )        :: thresh
      !!! Threshold for values in narrow band on input. Any zero-crossing from
      !!! a value larger than this will not be considered as interface. Set
      !!! this to the absolute value of whatever the field outside the band is
      !!! initialized to.
      REAL(MK), DIMENSION(:,:)         , POINTER,    OPTIONAL    :: closest
      !!! coordinates of the closest points on the interface from the mesh
      !!! points in ipts. Not unique! OPTIONAL. Only computed and returned if
      !!! present.
      INTEGER , DIMENSION(:,:)         , POINTER,    OPTIONAL    :: ipts
      !!! indices of mesh points adjacent to the interface. 1st index:
      !!! i,j,(k),jsub (local sub ID); 2nd: 1..npts. OPTIONAL has to be
      !!! present if closest is present! Not unique. The same point will occur
      !!! multiple times.
      INTEGER                          , INTENT(INOUT), OPTIONAL :: npts
      !!! ON input: if .LT. 0 the level is not recomputed and remains untouched.
      !!! Only the closest points are returned. On out: total number of points
      !!! adjacent to the interface. OPTIONAL. Has to be present if closest or
      !!! ipts are present!
      INTEGER                          , INTENT(  OUT)        :: info
      !!! Return status, 0 upon success.
      !-------------------------------------------------------------------------
      !  Local variables 
      !-------------------------------------------------------------------------
      INTEGER                     :: i,j,k,xhi,yhi,zhi,ip1,jp1,kp1,npos
      INTEGER                     :: isub,jsub,sneg,nit,ind,f,iopt,ii
      INTEGER                     :: m,n,p,kpp,jpn,pp2,np2,ll,ntotal,ltotal
      INTEGER ,DIMENSION(2)       :: ldu
      REAL(MK)                    :: t0,dx,dy,dz,nrm2,big,hsave
      REAL(MK)                    :: facx,facy,facz,facxy,facxz,facyz,facxyz
      REAL(MK)                    :: err,tol2,pxk,sprod,x,y,z
      REAL(MK),DIMENSION(ppm_dim) :: xk,x0,gradpxk,delta1,delta2,x0mxk,xkhalf
      REAL(MK),DIMENSION(:),POINTER :: phi
      LOGICAL                     :: lok
      INTEGER , DIMENSION(8)      :: neg
      REAL(MK), DIMENSION(0:3)    :: xv,yv,gxv,gyv
      CHARACTER(LEN=ppm_char)     :: mesg,cbuf
      TYPE(ppm_t_topo),      POINTER   :: topo
      TYPE(ppm_t_equi_mesh), POINTER   :: mesh
#if   __DIM == __2D
      REAL(MK), DIMENSION(16,16)  :: A
      INTEGER , DIMENSION(16)     :: Aind
      REAL(MK), DIMENSION(16)     :: rhs,coef
      REAL(MK), DIMENSION(4,4)    :: psi
      REAL(MK), DIMENSION(4)      :: ch
      REAL(MK), DIMENSION(2)      :: xr,xi
      INTEGER , PARAMETER         :: Nel = 16
#elif __DIM == __3D
      REAL(MK), DIMENSION(64,64)  :: A
      INTEGER , DIMENSION(64)     :: Aind
      REAL(MK), DIMENSION(64)     :: rhs,coef
      REAL(MK), DIMENSION(4,4,4)  :: psi
      REAL(MK), DIMENSION(0:3)    :: zv,gzv
      REAL(MK), DIMENSION(8)      :: ch
      REAL(MK), DIMENSION(3)      :: xr,xi
      INTEGER , PARAMETER         :: Nel = 64
#endif
      !-------------------------------------------------------------------------
      !  DATA
      !-------------------------------------------------------------------------
#if   __DIM == __3D
      !-------------------------------------------------------------------------
      !  LU factorization of the system matrix for the 3D tri-cubic
      !  interpolation polynomial. Created using ppm_gmm_create_3dmatrix.m.
      !-------------------------------------------------------------------------
#include "ppm_gmm_3dmatrix_data.inc"
      !-------------------------------------------------------------------------
      !  Permutation index vector of the LU factorization of the matrix.
      !  Created using ppm_gmm_create_3dmatrix.m.
      !-------------------------------------------------------------------------
#include "ppm_gmm_3dpermut_data.inc"
#elif __DIM == __2D
      !-------------------------------------------------------------------------
      !  LU factorization of the system matrix for the 2D bi-cubic
      !  interpolation polynomial. Created using ppm_gmm_create_2dmatrix.m.
      !-------------------------------------------------------------------------
#include "ppm_gmm_2dmatrix_data.inc"
      !-------------------------------------------------------------------------
      !  Permutation index vector of the LU factorization of the matrix.
      !  Created using ppm_gmm_create_2dmatrix.m.
      !-------------------------------------------------------------------------
#include "ppm_gmm_2dpermut_data.inc"
#endif   
      !-------------------------------------------------------------------------
      !  Initialise 
      !-------------------------------------------------------------------------
      CALL substart('ppm_gmm_kickoff',t0,info)
      ntotal = 0
      big    = HUGE(big)
      hsave  = 0.9_MK*big
      topo => ppm_topo(gmm_topoid)%t
      mesh => topo%mesh(gmm_meshid)
      !-------------------------------------------------------------------------
      !  Check arguments
      !-------------------------------------------------------------------------
      IF (ppm_debug .GT. 0) THEN
          IF (.NOT. ppm_initialized) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_ppm_noinit,'ppm_gmm_kickoff',  &
     &            'Please call ppm_init first!',__LINE__,info)
              GOTO 9999
          ENDIF
          IF (gmm_lsiz .LT. 0) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_ppm_noinit,'ppm_gmm_kickoff',  &
     &            'Please call gmm_init first!',__LINE__,info)
              GOTO 9999
          ENDIF
          IF (tol .LE. 0.0_MK) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_gmm_kickoff',  &
     &            'tolerance must be >0!',__LINE__,info)
              GOTO 9999
          ENDIF
          IF (PRESENT(closest)) THEN
              IF (.NOT.PRESENT(ipts)) THEN
                  info = ppm_error_error
                  CALL ppm_error(ppm_err_argument,'ppm_gmm_kickoff',  &
     &                'ipts must be present if closest is',__LINE__,info)
                  GOTO 9999
              ENDIF
              IF (.NOT.PRESENT(npts)) THEN
                  info = ppm_error_error
                  CALL ppm_error(ppm_err_argument,'ppm_gmm_kickoff',  &
     &                'npts must be present if closest is',__LINE__,info)
                  GOTO 9999
              ENDIF
          ENDIF
          IF (PRESENT(ipts)) THEN
              IF (.NOT.PRESENT(npts)) THEN
                  info = ppm_error_error
                  CALL ppm_error(ppm_err_argument,'ppm_gmm_kickoff',  &
     &                'npts must be present if ipts is',__LINE__,info)
                  GOTO 9999
              ENDIF
          ENDIF
          IF (.NOT.ASSOCIATED(fdata)) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_gmm_kickoff',  &
     &            'field data is not allocated!',__LINE__,info)
              GOTO 9999
          ENDIF
#if   __DIM == __3D
          IF (SIZE(fdata,4) .LT. topo%nsublist) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_gmm_kickoff',  &
     &            'field data for some subs is missing',__LINE__,info)
              GOTO 9999
          ENDIF
          IF (UBOUND(fdata,1) .LT. maxxhi) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_gmm_kickoff',  &
     &            'x dimension of field data does not match mesh',__LINE__,info)
              GOTO 9999
          ENDIF
          IF (UBOUND(fdata,2) .LT. maxyhi) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_gmm_kickoff',  &
     &            'y dimension of field data does not match mesh',__LINE__,info)
              GOTO 9999
          ENDIF
          IF (UBOUND(fdata,3) .LT. maxzhi) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_gmm_kickoff',  &
     &            'z dimension of field data does not match mesh',__LINE__,info)
              GOTO 9999
          ENDIF
#elif __DIM == __2D
          IF (SIZE(fdata,3) .LT. topo%nsublist) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_gmm_kickoff',  &
     &            'field data for some subs is missing',__LINE__,info)
              GOTO 9999
          ENDIF
          IF (UBOUND(fdata,1) .LT. maxxhi) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_gmm_kickoff',  &
     &            'x dimension of field data does not match mesh',__LINE__,info)
              GOTO 9999
          ENDIF
          IF (UBOUND(fdata,2) .LT. maxyhi) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_gmm_kickoff',  &
     &            'y dimension of field data does not match mesh',__LINE__,info)
              GOTO 9999
          ENDIF
#endif
      ENDIF          ! ppm_debug for argument check
      !-------------------------------------------------------------------------
      !  Set the pointers to work memory (in module data_gmm)
      !-------------------------------------------------------------------------
#if    __KIND == __SINGLE_PRECISION
      phi => gmm_phis
#elif  __KIND == __DOUBLE_PRECISION
      phi => gmm_phid
#endif 
      !-------------------------------------------------------------------------
      !  Set constants
      !-------------------------------------------------------------------------
      xv(0) = 1.0_MK
      yv(0) = 1.0_MK
      gxv(0) = 0.0_MK
      gxv(1) = 1.0_MK
      gyv(0) = 0.0_MK
      gyv(1) = 1.0_MK
#if   __DIM == __3D
      zv(0) = 1.0_MK
      gzv(0) = 0.0_MK
      gzv(1) = 1.0_MK
#endif
      !-------------------------------------------------------------------------
      !  Nuke OPTIONAL arrays
      !-------------------------------------------------------------------------
      IF (PRESENT(closest)) THEN
          iopt = ppm_param_dealloc
          CALL ppm_alloc(closest,ldu,iopt,info)
          IF (info .NE. 0) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_dealloc,'ppm_gmm_kickoff',      &
     &            'closest points CLOSEST',__LINE__,info)
          ENDIF
      ENDIF
      IF (PRESENT(ipts)) THEN
          iopt = ppm_param_dealloc
          CALL ppm_alloc(ipts,ldu,iopt,info)
          IF (info .NE. 0) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_dealloc,'ppm_gmm_kickoff',      &
     &            'point locations IPTS',__LINE__,info)
          ENDIF
      ENDIF
      ltotal = 0
      !-------------------------------------------------------------------------
      !  Find mesh spacing
      !-------------------------------------------------------------------------
      IF (ppm_kind .EQ. ppm_kind_single) THEN
          dx = (topo%max_physs(1)-topo%min_physs(1))/   &
     &        REAL(mesh%Nm(1)-1,ppm_kind_single)
          dy = (topo%max_physs(2)-topo%min_physs(2))/  &
     &        REAL(mesh%Nm(2)-1,ppm_kind_single)
          IF (ppm_dim .GT. 2) THEN
              dz = (topo%max_physs(3)-topo%min_physs(3))/ &
     &            REAL(mesh%Nm(3)-1,     &
     &            ppm_kind_single)
          ENDIF
      ELSE
          dx = (topo%max_physs(1)-topo%min_physs(1))/   &
     &        REAL(mesh%Nm(1)-1,ppm_kind_double)
          dy = (topo%max_physs(2)-topo%min_physs(2))/  &
     &        REAL(mesh%Nm(2)-1,ppm_kind_double)
          IF (ppm_dim .GT. 2) THEN
              dz = (topo%max_physs(3)-topo%min_physs(3))/ &
     &            REAL(mesh%Nm(3)-1,     &
     &            ppm_kind_double)
          ENDIF
      ENDIF
      !-------------------------------------------------------------------------
      !  Loop over all mesh blocks 
      !-------------------------------------------------------------------------
      DO jsub=1,topo%nsublist
          npos= 0
          isub = topo%isublist(jsub)
          xhi = mesh%nnodes(1,isub)-1
          yhi = mesh%nnodes(2,isub)-1
          facx= 0.5_MK
          facy= 0.5_MK
          facxy= 0.25_MK
#if   __DIM == __3D
          zhi = mesh%nnodes(3,isub)-1
          facz= 0.5_MK
          facxz= 0.25_MK
          facyz= 0.25_MK
          facxyz= 0.125_MK
          !---------------------------------------------------------------------
          !  Find mesh cells through which the zero level passes and
          !  that are close to the interface (to exclude the boundaries
          !  of the narrow band)
          !---------------------------------------------------------------------
          DO k=1,zhi
              kp1 = k+1
              DO j=1,yhi
                  jp1 = j+1
                  DO i=1,xhi
                      ! storing the shifted index is 50% faster then
                      ! doing the additions. According to test program.
                      ip1 = i+1
                      neg(1:8) = 0
                      lok = .TRUE.
                      IF (fdata(i  ,j  ,k  ,jsub) .LT. 0.0_MK) neg(1) = 1
                      IF (ABS(fdata(i,j,k,jsub)) .GT. thresh) lok = .FALSE.
                      IF (fdata(ip1,j  ,k  ,jsub) .LT. 0.0_MK) neg(2) = 1
                      IF (ABS(fdata(ip1,j,k,jsub)) .GT. thresh) lok = .FALSE.
                      IF (fdata(i  ,jp1,k  ,jsub) .LT. 0.0_MK) neg(3) = 1
                      IF (ABS(fdata(i,jp1,k,jsub)) .GT. thresh) lok = .FALSE.
                      IF (fdata(ip1,jp1,k  ,jsub) .LT. 0.0_MK) neg(4) = 1
                      IF (ABS(fdata(ip1,jp1,k,jsub)) .GT. thresh) lok = .FALSE.
                      IF (fdata(i  ,j  ,kp1,jsub) .LT. 0.0_MK) neg(5) = 1
                      IF (ABS(fdata(i,j,kp1,jsub)) .GT. thresh) lok = .FALSE.
                      IF (fdata(ip1,j  ,kp1,jsub) .LT. 0.0_MK) neg(6) = 1
                      IF (ABS(fdata(ip1,j,kp1,jsub)) .GT. thresh) lok = .FALSE.
                      IF (fdata(i  ,jp1,kp1,jsub) .LT. 0.0_MK) neg(7) = 1
                      IF (ABS(fdata(i,jp1,kp1,jsub)) .GT. thresh) lok = .FALSE.
                      IF (fdata(ip1,jp1,kp1,jsub) .LT. 0.0_MK) neg(8) = 1
                      IF (ABS(fdata(ip1,jp1,kp1,jsub)).GT.thresh) lok = .FALSE.
                      sneg = neg(1)+neg(2)+neg(3)+neg(4)+neg(5)+neg(6)+  &
     &                       neg(7)+neg(8)
                      IF ((sneg.GT.0).AND.(sneg.LT.8).AND.(lok)) THEN
                          !-----------------------------------------------------
                          !  Old value of the level function
                          !-----------------------------------------------------
                          DO p=-1,2
                              pp2 = p+2
                              kpp = k+p
                              DO n=-1,2
                                  np2 = n+2
                                  jpn = j+n
                                  DO m=-1,2
                                      psi(m+2,np2,pp2) = fdata(i+m,jpn,kpp,jsub)
                                  ENDDO
                              ENDDO
                          ENDDO
                          !-----------------------------------------------------
                          !  Right-hand-side for tri-cubic interpolation
                          !  polynomial. The order of the method is
                          !  limited by the order of the finite
                          !  differences used here.
                          !-----------------------------------------------------
                          ! value at (i,j,k)
                          rhs(1) = psi(2,2,2)
                          ! value at (i+1,j,k)
                          rhs(2) = psi(3,2,2)
                          ! value at (i,j+1,k)
                          rhs(3) = psi(2,3,2)
                          ! value at (i+1,j+1,k)
                          rhs(4) = psi(3,3,2)
                          ! value at (i,j,k+1)
                          rhs(5) = psi(2,2,3)
                          ! value at (i+1,j,k+1)
                          rhs(6) = psi(3,2,3)
                          ! value at (i,j+1,k+1)
                          rhs(7) = psi(2,3,3)
                          ! value at (i+1,j+1,k+1)
                          rhs(8) = psi(3,3,3)
                          !-----------------------------------------------------
                          ! Dx at (i,j,k)
                          rhs(9) = facx*(psi(3,2,2)-psi(1,2,2))
                          ! Dx at (i+1,j,k)
                          rhs(10) = facx*(psi(4,2,2)-psi(2,2,2))
                          ! Dx at (i,j+1,k)
                          rhs(11) = facx*(psi(3,3,2)-psi(1,3,2))
                          ! Dx at (i+1,j+1,k)
                          rhs(12) = facx*(psi(4,3,2)-psi(2,3,2))
                          ! Dx at (i,j,k+1)
                          rhs(13) = facx*(psi(3,2,3)-psi(1,2,3))
                          ! Dx at (i+1,j,k+1)
                          rhs(14) = facx*(psi(4,2,3)-psi(2,2,3))
                          ! Dx at (i,j+1,k+1)
                          rhs(15) = facx*(psi(3,3,3)-psi(1,3,3))
                          ! Dx at (i+1,j+1,k+1)
                          rhs(16) = facx*(psi(4,3,3)-psi(2,3,3))
                          !-----------------------------------------------------
                          ! Dy at (i,j,k)
                          rhs(17) = facy*(psi(2,3,2)-psi(2,1,2))
                          ! Dy at (i+1,j,k)
                          rhs(18) = facy*(psi(3,3,2)-psi(3,1,2))
                          ! Dy at (i,j+1,k)
                          rhs(19) = facy*(psi(2,4,2)-psi(2,2,2))
                          ! Dy at (i+1,j+1,k)
                          rhs(20) = facy*(psi(3,4,2)-psi(3,2,2))
                          ! Dy at (i,j,k+1)
                          rhs(21) = facy*(psi(2,3,3)-psi(2,1,3))
                          ! Dy at (i+1,j,k+1)
                          rhs(22) = facy*(psi(3,3,3)-psi(3,1,3))
                          ! Dy at (i,j+1,k+1)
                          rhs(23) = facy*(psi(2,4,3)-psi(2,2,3))
                          ! Dy at (i+1,j+1,k+1)
                          rhs(24) = facy*(psi(3,4,3)-psi(3,2,3))
                          !-----------------------------------------------------
                          ! Dz at (i,j,k)
                          rhs(25) = facz*(psi(2,2,3)-psi(2,2,1))
                          ! Dz at (i+1,j,k)
                          rhs(26) = facz*(psi(3,2,3)-psi(3,2,1))
                          ! Dz at (i,j+1,k)
                          rhs(27) = facz*(psi(2,3,3)-psi(2,3,1))
                          ! Dz at (i+1,j+1,k)
                          rhs(28) = facz*(psi(3,3,3)-psi(3,3,1))
                          ! Dz at (i,j,k+1)
                          rhs(29) = facz*(psi(2,2,4)-psi(2,2,2))
                          ! Dz at (i+1,j,k+1)
                          rhs(30) = facz*(psi(3,2,4)-psi(3,2,2))
                          ! Dz at (i,j+1,k+1)
                          rhs(31) = facz*(psi(2,3,4)-psi(2,3,2))
                          ! Dz at (i+1,j+1,k+1)
                          rhs(32) = facz*(psi(3,3,4)-psi(3,3,2))
                          !-----------------------------------------------------
                          ! DxDy at (i,j,k)
                          rhs(33) = facxy*(psi(3,3,2)-psi(1,3,2)-   &
     &                                     psi(3,1,2)+psi(1,1,2))
                          ! DxDy at (i+1,j,k)
                          rhs(34) = facxy*(psi(4,3,2)-psi(2,3,2)-   &
     &                                     psi(4,1,2)+psi(2,1,2))
                          ! DxDy at (i,j+1,k)
                          rhs(35) = facxy*(psi(3,4,2)-psi(1,4,2)-   &
     &                                     psi(3,2,2)+psi(1,2,2))
                          ! DxDy at (i+1,j+1,k)
                          rhs(36) = facxy*(psi(4,4,2)-psi(2,4,2)-   &
     &                                     psi(4,2,2)+psi(2,2,2))
                          ! DxDy at (i,j,k+1)
                          rhs(37) = facxy*(psi(3,3,3)-psi(1,3,3)-   &
     &                                     psi(3,1,3)+psi(1,1,3))
                          ! DxDy at (i+1,j,k+1)
                          rhs(38) = facxy*(psi(4,3,3)-psi(2,3,3)-   &
     &                                     psi(4,1,3)+psi(2,1,3))
                          ! DxDy at (i,j+1,k+1)
                          rhs(39) = facxy*(psi(3,4,3)-psi(1,4,3)-   &
     &                                     psi(3,2,3)+psi(1,2,3))
                          ! DxDy at (i+1,j+1,k+1)
                          rhs(40) = facxy*(psi(4,4,3)-psi(2,4,3)-   &
     &                                     psi(4,2,3)+psi(2,2,3))
                          !-----------------------------------------------------
                          ! DxDz at (i,j,k)
                          rhs(41) = facxz*(psi(3,2,3)-psi(1,2,3)-   &
     &                                     psi(3,2,1)+psi(1,2,1))
                          ! DxDz at (i+1,j,k)
                          rhs(42) = facxz*(psi(4,2,3)-psi(2,2,3)-   &
     &                                     psi(4,2,1)+psi(2,2,1))
                          ! DxDz at (i,j+1,k)
                          rhs(43) = facxz*(psi(3,3,3)-psi(1,3,3)-   &
     &                                     psi(3,3,1)+psi(1,3,1))
                          ! DxDz at (i+1,j+1,k)
                          rhs(44) = facxz*(psi(4,3,3)-psi(2,3,3)-   &
     &                                     psi(4,3,1)+psi(2,3,1))
                          ! DxDz at (i,j,k+1)
                          rhs(45) = facxz*(psi(3,2,4)-psi(1,2,4)-   &
     &                                     psi(3,2,2)+psi(1,2,2))
                          ! DxDz at (i+1,j,k+1)
                          rhs(46) = facxz*(psi(4,2,4)-psi(2,2,4)-   &
     &                                     psi(4,2,2)+psi(2,2,2))
                          ! DxDz at (i,j+1,k+1)
                          rhs(47) = facxz*(psi(3,3,4)-psi(1,3,4)-   &
     &                                     psi(3,3,2)+psi(1,3,2))
                          ! DxDz at (i+1,j+1,k+1)
                          rhs(48) = facxz*(psi(4,3,4)-psi(2,3,4)-   &
     &                                     psi(4,3,2)+psi(2,3,2))
                          !-----------------------------------------------------
                          ! DyDz at (i,j,k)
                          rhs(49) = facyz*(psi(2,3,3)-psi(2,1,3)-   &
     &                                     psi(2,3,1)+psi(2,1,1))
                          ! DyDz at (i+1,j,k)
                          rhs(50) = facyz*(psi(3,3,3)-psi(3,1,3)-   &
     &                                     psi(3,3,1)+psi(3,1,1))
                          ! DyDz at (i,j+1,k)
                          rhs(51) = facyz*(psi(2,4,3)-psi(2,2,3)-   &
     &                                     psi(2,4,1)+psi(2,2,1))
                          ! DyDz at (i+1,j+1,k)
                          rhs(52) = facyz*(psi(3,4,3)-psi(3,2,3)-   &
     &                                     psi(3,4,1)+psi(3,2,1))
                          ! DyDz at (i,j,k+1)
                          rhs(53) = facyz*(psi(2,3,4)-psi(2,1,4)-   &
     &                                     psi(2,3,2)+psi(2,1,2))
                          ! DyDz at (i+1,j,k+1)
                          rhs(54) = facyz*(psi(3,3,4)-psi(3,1,4)-   &
     &                                     psi(3,3,2)+psi(3,1,2))
                          ! DyDz at (i,j+1,k+1)
                          rhs(55) = facyz*(psi(2,4,4)-psi(2,2,4)-   &
     &                                     psi(2,4,2)+psi(2,2,2))
                          ! DyDz at (i+1,j+1,k+1)
                          rhs(56) = facyz*(psi(3,4,4)-psi(3,2,4)-   &
     &                                     psi(3,4,2)+psi(3,2,2))
                          !-----------------------------------------------------
                          ! DxDyDz at (i,j,k)
                          rhs(57) = facxyz*(psi(3,3,3)-psi(3,3,1)-  &
     &                        psi(3,1,3)+psi(3,1,1)-psi(1,3,3)+     &
     &                        psi(1,3,1)+psi(1,1,3)-psi(1,1,1))
                          ! DxDyDz at (i+1,j,k)
                          rhs(58) = facxyz*(psi(4,3,3)-psi(4,3,1)-  &
     &                        psi(4,1,3)+psi(4,1,1)-psi(2,3,3)+     &
     &                        psi(2,3,1)+psi(2,1,3)-psi(2,1,1))
                          ! DxDyDz at (i,j+1,k)
                          rhs(59) = facxyz*(psi(3,4,3)-psi(3,4,1)-  &
     &                        psi(3,2,3)+psi(3,2,1)-psi(1,4,3)+     &
     &                        psi(1,4,1)+psi(1,2,3)-psi(1,2,1))
                          ! DxDyDz at (i+1,j+1,k)
                          rhs(60) = facxyz*(psi(4,4,3)-psi(4,4,1)-  &
     &                        psi(4,2,3)+psi(4,2,1)-psi(2,4,3)+     &
     &                        psi(2,4,1)+psi(2,2,3)-psi(2,2,1))
                          ! DxDyDz at (i,j,k+1)
                          rhs(61) = facxyz*(psi(3,3,4)-psi(3,3,2)-  &
     &                        psi(3,1,4)+psi(3,1,2)-psi(1,3,4)+     &
     &                        psi(1,3,2)+psi(1,1,4)-psi(1,1,2))
                          ! DxDyDz at (i+1,j,k+1)
                          rhs(62) = facxyz*(psi(4,3,4)-psi(4,3,2)-  &
     &                        psi(4,1,4)+psi(4,1,2)-psi(2,3,4)+     &
     &                        psi(2,3,2)+psi(2,1,4)-psi(2,1,2))
                          ! DxDyDz at (i,j+1,k+1)
                          rhs(63) = facxyz*(psi(3,4,4)-psi(3,4,2)-  &
     &                        psi(3,2,4)+psi(3,2,2)-psi(1,4,4)+     &
     &                        psi(1,4,2)+psi(1,2,4)-psi(1,2,2))
                          ! DxDyDz at (i+1,j+1,k+1)
                          rhs(64) = facxyz*(psi(4,4,4)-psi(4,4,2)-  &
     &                        psi(4,2,4)+psi(4,2,2)-psi(2,4,4)+     &
     &                        psi(2,4,2)+psi(2,2,4)-psi(2,2,2))
                          !-----------------------------------------------------
                          !  Check if needed
                          !-----------------------------------------------------
                          IF (ppm_debug .GT. 1) THEN
                              DO m=1,64
                                  IF ((rhs(m) .NE. rhs(m)) .OR.    &
     &                                (ABS(rhs(m)) .GT. hsave)) THEN
                                      !-----------------------------------------
                                      !  rhs(m) is Inf or NaN
                                      !-----------------------------------------
                                      CALL ppm_write(ppm_rank,   &
     &                                    'ppm_gmm_kickoff',     &
     &          'WARNING: RHS not defined! Check support of input field!',info)
                                  ENDIF
                              ENDDO
                          ENDIF
                          !-----------------------------------------------------
                          !  Solve for polynomial coefficients. coef will
                          !  contain the result. 
                          !-----------------------------------------------------
#include "ppm_gmm_lubksb.inc"
! CHECK SOLUTION
!                         U = 0.0_MK
!                         L = 0.0_MK
!                         LU = 0.0_MK
!                         DO ii=1,64
!                             DO jj=ii,64
!                                 U(ii,jj) = A(ii,jj)
!                             ENDDO
!                         ENDDO
!                         DO ii=1,64
!                             DO jj=1,(ii-1)
!                                 L(ii,jj) = A(ii,jj)
!                             ENDDO
!                         ENDDO
!                         DO ii=1,64
!                             L(ii,ii) = L(ii,ii) + 1.0_MK
!                         ENDDO
!                         LU = MATMUL(L,U)
!                         gaga = MATMUL(LU,coef)
!                         nrm2 = 0.0_MK
!                         DO ii=1,64
!                             ll = Aind(ii)
!                             nrm2 = nrm2 + (gaga(ii)-rhs(ll))**2
!                         ENDDO
!                    !     IF (nrm2 .NE. nrm2) THEN
!                             PRINT*,'ERROR: ',nrm2
!                    !         PRINT*,'rhs: ',rhs
!                    !         PRINT*,'coefs: ',coef
!                    !     ENDIF
! Wanna see them interpolants in MATLAB??? Uncomment this!
!                     WRITE(cbuf,'(A,I3.3,A)') 'interp_',npos/8,'.out'
!                     OPEN(30,FILE=cbuf,STATUS='REPLACE',ACTION='WRITE')
!                     WRITE(30,'(3E20.8)') topo%min_subs(1,isub,gmm_topoid)+(i-1)*dx, &
!    &                               topo%min_subs(2,isub,gmm_topoid)+(j-1)*dy, &
!    &                               topo%min_subs(3,isub,gmm_topoid)+(k-1)*dz
!                     WRITE(30,'(3E20.8)') topo%min_subs(1,isub,gmm_topoid)+(i)*dx, &
!    &                               topo%min_subs(2,isub,gmm_topoid)+(j-1)*dy, &
!    &                               topo%min_subs(3,isub,gmm_topoid)+(k-1)*dz
!                     WRITE(30,'(3E20.8)') topo%min_subs(1,isub,gmm_topoid)+(i-1)*dx, &
!    &                               topo%min_subs(2,isub,gmm_topoid)+(j)*dy, &
!    &                               topo%min_subs(3,isub,gmm_topoid)+(k-1)*dz
!                     WRITE(30,'(3E20.8)') topo%min_subs(1,isub,gmm_topoid)+(i)*dx, &
!    &                               topo%min_subs(2,isub,gmm_topoid)+(j)*dy, &
!    &                               topo%min_subs(3,isub,gmm_topoid)+(k-1)*dz
!                     WRITE(30,'(3E20.8)') topo%min_subs(1,isub,gmm_topoid)+(i-1)*dx, &
!    &                               topo%min_subs(2,isub,gmm_topoid)+(j-1)*dy, &
!    &                               topo%min_subs(3,isub,gmm_topoid)+(k)*dz
!                     WRITE(30,'(3E20.8)') topo%min_subs(1,isub,gmm_topoid)+(i)*dx, &
!    &                               topo%min_subs(2,isub,gmm_topoid)+(j-1)*dy, &
!    &                               topo%min_subs(3,isub,gmm_topoid)+(k)*dz
!                     WRITE(30,'(3E20.8)') topo%min_subs(1,isub,gmm_topoid)+(i-1)*dx, &
!    &                               topo%min_subs(2,isub,gmm_topoid)+(j)*dy, &
!    &                               topo%min_subs(3,isub,gmm_topoid)+(k)*dz
!                     WRITE(30,'(3E20.8)') topo%min_subs(1,isub,gmm_topoid)+(i)*dx, &
!    &                               topo%min_subs(2,isub,gmm_topoid)+(j)*dy, &
!    &                               (3,isub,gmm_topoid)+(k)*dz
!                     WRITE(30,'(2E20.8)') rhs(1),rhs(2)
!                     WRITE(30,'(2E20.8)') rhs(3),rhs(4)
!                     WRITE(30,'(2E20.8)') rhs(5),rhs(6)
!                     WRITE(30,'(2E20.8)') rhs(7),rhs(8)
!                     DO ii=1,64
!                         WRITE(30,'(E20.8)') coef(ii)
!                     ENDDO
!                     CLOSE(30)
                          !-----------------------------------------------------
                          !  Check that sparse structure has at least 8
                          !  spaces to hold the new data
                          !-----------------------------------------------------
                          IF (gmm_lsiz .LT. npos+8) THEN
                              gmm_lsiz = gmm_lsiz + incr
                              iopt = ppm_param_alloc_grow_preserve
                              ldu(1) = gmm_lsiz
                              CALL ppm_alloc(phi,ldu,iopt,info)
                              IF (info .NE. 0) THEN
                                  info = ppm_error_fatal
                                  CALL ppm_error(ppm_err_alloc,   &
     &                                'ppm_gmm_kickoff',      &
     &                                'sparse values GMM_PHI',__LINE__,info)
                                  GOTO 9999
                              ENDIF
                              ldu(1) = 4
                              ldu(2) = gmm_lsiz
                              CALL ppm_alloc(gmm_ipos,ldu,iopt,info)
                              IF (info .NE. 0) THEN
                                  info = ppm_error_fatal
                                  CALL ppm_error(ppm_err_alloc,   &
     &                                'ppm_gmm_kickoff',      &
     &                                'sparse positions GMM_IPOS',__LINE__,info)
                                  GOTO 9999
                              ENDIF
                          ENDIF
                          IF (PRESENT(closest)) THEN
                              IF (ltotal .LT. ntotal+8) THEN
                                  ltotal = ntotal + incr
                                  iopt   = ppm_param_alloc_grow_preserve
                                  ldu(1) = 3
                                  ldu(2) = ltotal
                                  CALL ppm_alloc(closest,ldu,iopt,info)
                                  IF (info .NE. 0) THEN
                                      info = ppm_error_fatal
                                      CALL ppm_error(ppm_err_alloc,   &
     &                                    'ppm_gmm_kickoff',          &
     &                                    'closest points CLOSEST',   &
     &                                    __LINE__,info)
                                      GOTO 9999
                                  ENDIF
                              ENDIF
                          ENDIF
                          !-----------------------------------------------------
                          !  Newton iteration for each point bounding
                          !  the mesh cell to find the new level
                          !  function value and the closest point on the
                          !  interface.
                          !-----------------------------------------------------
                          IF (PRESENT(closest).AND..NOT.PRESENT(chi)) THEN
#if   __KIND == __SINGLE_PRECISION
                              x = topo%min_subs(1,isub)+(i-1)*dx
                              y = topo%min_subs(2,isub)+(j-1)*dy
                              z = topo%min_subs(3,isub)+(k-1)*dz
#elif __KIND == __DOUBLE_PRECISION
                              x = topo%min_subd(1,isub)+(i-1)*dx
                              y = topo%min_subd(2,isub)+(j-1)*dy
                              z = topo%min_subd(3,isub)+(k-1)*dz
#endif
                          ENDIF
                          x0    = 0.0_MK
                          IF (PRESENT(chi)) THEN
                              x = chi(1,i,j,k,jsub)
                              y = chi(2,i,j,k,jsub)
                              z = chi(3,i,j,k,jsub)
                          ENDIF
#include "ppm_gmm_quasinewton.inc"
                          IF (neg(1).EQ.1) sprod = -sprod
                          npos                   = npos + 1
                          ntotal                 = ntotal + 1
                          gmm_ipos(1,npos)       = i
                          gmm_ipos(2,npos)       = j
                          gmm_ipos(3,npos)       = k
                          gmm_ipos(4,npos)       = jsub
                          phi(npos)              = sprod
                          ! the closest point on the interface is xk
                          IF (PRESENT(closest)) THEN
                              IF (PRESENT(chi)) THEN
                                  closest(1,ntotal)  = xr(1)
                                  closest(2,ntotal)  = xr(2)
                                  closest(3,ntotal)  = xr(3)
                              ELSE
                                  closest(1,ntotal)  = x + dx*xk(1)
                                  closest(2,ntotal)  = y + dy*xk(2)
                                  closest(3,ntotal)  = z + dz*xk(3)
                              ENDIF
                          ENDIF
                          !-----------------------------------------------------
                          x0(1) = 1.0_MK
                          x0(2) = 0.0_MK
                          x0(3) = 0.0_MK
                          IF (PRESENT(chi)) THEN
                              x = chi(1,ip1,j,k,jsub)
                              y = chi(2,ip1,j,k,jsub)
                              z = chi(3,ip1,j,k,jsub)
                          ENDIF
#include "ppm_gmm_quasinewton.inc"
                          IF (neg(2).EQ.1) sprod = -sprod
                          npos                   = npos + 1
                          ntotal                 = ntotal + 1
                          gmm_ipos(1,npos)       = ip1
                          gmm_ipos(2,npos)       = j
                          gmm_ipos(3,npos)       = k
                          gmm_ipos(4,npos)       = jsub
                          phi(npos)              = sprod
                          ! the closest point on the interface is xk
                          IF (PRESENT(closest)) THEN
                              IF (PRESENT(chi)) THEN
                                  closest(1,ntotal)  = xr(1)
                                  closest(2,ntotal)  = xr(2)
                                  closest(3,ntotal)  = xr(3)
                              ELSE
                                  closest(1,ntotal)  = x + dx*xk(1)
                                  closest(2,ntotal)  = y + dy*xk(2)
                                  closest(3,ntotal)  = z + dz*xk(3)
                              ENDIF
                          ENDIF
                          !-----------------------------------------------------
                          x0(1) = 0.0_MK
                          x0(2) = 1.0_MK
                          x0(3) = 0.0_MK
                          IF (PRESENT(chi)) THEN
                              x = chi(1,i,jp1,k,jsub)
                              y = chi(2,i,jp1,k,jsub)
                              z = chi(3,i,jp1,k,jsub)
                          ENDIF
#include "ppm_gmm_quasinewton.inc"
                          IF (neg(3).EQ.1) sprod = -sprod
                          npos                   = npos + 1
                          ntotal                 = ntotal + 1
                          gmm_ipos(1,npos)       = i
                          gmm_ipos(2,npos)       = jp1
                          gmm_ipos(3,npos)       = k
                          gmm_ipos(4,npos)       = jsub
                          phi(npos)              = sprod
                          IF (PRESENT(closest)) THEN
                              IF (PRESENT(chi)) THEN
                                  closest(1,ntotal)  = xr(1)
                                  closest(2,ntotal)  = xr(2)
                                  closest(3,ntotal)  = xr(3)
                              ELSE
                                  closest(1,ntotal)  = x + dx*xk(1)
                                  closest(2,ntotal)  = y + dy*xk(2)
                                  closest(3,ntotal)  = z + dz*xk(3)
                              ENDIF
                          ENDIF
                          !-----------------------------------------------------
                          x0(1) = 1.0_MK
                          x0(2) = 1.0_MK
                          x0(3) = 0.0_MK
                          IF (PRESENT(chi)) THEN
                              x = chi(1,ip1,jp1,k,jsub)
                              y = chi(2,ip1,jp1,k,jsub)
                              z = chi(3,ip1,jp1,k,jsub)
                          ENDIF
#include "ppm_gmm_quasinewton.inc"
                          IF (neg(4).EQ.1) sprod = -sprod
                          npos                   = npos + 1
                          ntotal                 = ntotal + 1
                          gmm_ipos(1,npos)       = ip1
                          gmm_ipos(2,npos)       = jp1
                          gmm_ipos(3,npos)       = k
                          gmm_ipos(4,npos)       = jsub
                          phi(npos)              = sprod
                          ! the closest point on the interface is xk
                          IF (PRESENT(closest)) THEN
                              IF (PRESENT(chi)) THEN
                                  closest(1,ntotal)  = xr(1)
                                  closest(2,ntotal)  = xr(2)
                                  closest(3,ntotal)  = xr(3)
                              ELSE
                                  closest(1,ntotal)  = x + dx*xk(1)
                                  closest(2,ntotal)  = y + dy*xk(2)
                                  closest(3,ntotal)  = z + dz*xk(3)
                              ENDIF
                          ENDIF
                          !-----------------------------------------------------
                          x0(1) = 0.0_MK
                          x0(2) = 0.0_MK
                          x0(3) = 1.0_MK
                          IF (PRESENT(chi)) THEN
                              x = chi(1,i,j,kp1,jsub)
                              y = chi(2,i,j,kp1,jsub)
                              z = chi(3,i,j,kp1,jsub)
                          ENDIF
#include "ppm_gmm_quasinewton.inc"
                          IF (neg(5).EQ.1) sprod = -sprod
                          npos                   = npos + 1
                          ntotal                 = ntotal + 1
                          gmm_ipos(1,npos)       = i
                          gmm_ipos(2,npos)       = j
                          gmm_ipos(3,npos)       = kp1
                          gmm_ipos(4,npos)       = jsub
                          phi(npos)              = sprod
                          ! the closest point on the interface is xk
                          IF (PRESENT(closest)) THEN
                              IF (PRESENT(chi)) THEN
                                  closest(1,ntotal)  = xr(1)
                                  closest(2,ntotal)  = xr(2)
                                  closest(3,ntotal)  = xr(3)
                              ELSE
                                  closest(1,ntotal)  = x + dx*xk(1)
                                  closest(2,ntotal)  = y + dy*xk(2)
                                  closest(3,ntotal)  = z + dz*xk(3)
                              ENDIF
                          ENDIF
                          !-----------------------------------------------------
                          x0(1) = 1.0_MK
                          x0(2) = 0.0_MK
                          x0(3) = 1.0_MK
                          IF (PRESENT(chi)) THEN
                              x = chi(1,ip1,j,kp1,jsub)
                              y = chi(2,ip1,j,kp1,jsub)
                              z = chi(3,ip1,j,kp1,jsub)
                          ENDIF
#include "ppm_gmm_quasinewton.inc"
                          IF (neg(6).EQ.1) sprod = -sprod
                          npos                   = npos + 1
                          ntotal                 = ntotal + 1
                          gmm_ipos(1,npos)       = ip1
                          gmm_ipos(2,npos)       = j
                          gmm_ipos(3,npos)       = kp1
                          gmm_ipos(4,npos)       = jsub
                          phi(npos)              = sprod
                          ! the closest point on the interface is xk
                          IF (PRESENT(closest)) THEN
                              IF (PRESENT(chi)) THEN
                                  closest(1,ntotal)  = xr(1)
                                  closest(2,ntotal)  = xr(2)
                                  closest(3,ntotal)  = xr(3)
                              ELSE
                                  closest(1,ntotal)  = x + dx*xk(1)
                                  closest(2,ntotal)  = y + dy*xk(2)
                                  closest(3,ntotal)  = z + dz*xk(3)
                              ENDIF
                          ENDIF
                          !-----------------------------------------------------
                          x0(1) = 0.0_MK
                          x0(2) = 1.0_MK
                          x0(3) = 1.0_MK
                          IF (PRESENT(chi)) THEN
                              x = chi(1,i,jp1,kp1,jsub)
                              y = chi(2,i,jp1,kp1,jsub)
                              z = chi(3,i,jp1,kp1,jsub)
                          ENDIF
#include "ppm_gmm_quasinewton.inc"
                          IF (neg(7).EQ.1) sprod = -sprod
                          npos                   = npos + 1
                          ntotal                 = ntotal + 1
                          gmm_ipos(1,npos)       = i
                          gmm_ipos(2,npos)       = jp1
                          gmm_ipos(3,npos)       = kp1
                          gmm_ipos(4,npos)       = jsub
                          phi(npos)              = sprod
                          ! the closest point on the interface is xk
                          IF (PRESENT(closest)) THEN
                              IF (PRESENT(chi)) THEN
                                  closest(1,ntotal)  = xr(1)
                                  closest(2,ntotal)  = xr(2)
                                  closest(3,ntotal)  = xr(3)
                              ELSE
                                  closest(1,ntotal)  = x + dx*xk(1)
                                  closest(2,ntotal)  = y + dy*xk(2)
                                  closest(3,ntotal)  = z + dz*xk(3)
                              ENDIF
                          ENDIF
                          !-----------------------------------------------------
                          x0(1) = 1.0_MK
                          x0(2) = 1.0_MK
                          x0(3) = 1.0_MK
                          IF (PRESENT(chi)) THEN
                              x = chi(1,ip1,jp1,kp1,jsub)
                              y = chi(2,ip1,jp1,kp1,jsub)
                              z = chi(3,ip1,jp1,kp1,jsub)
                          ENDIF
#include "ppm_gmm_quasinewton.inc"
                          IF (neg(8).EQ.1) sprod = -sprod
                          npos                   = npos + 1
                          ntotal                 = ntotal + 1
                          gmm_ipos(1,npos)       = ip1
                          gmm_ipos(2,npos)       = jp1
                          gmm_ipos(3,npos)       = kp1
                          gmm_ipos(4,npos)       = jsub
                          phi(npos)              = sprod
                          ! the closest point on the interface is xk
                          IF (PRESENT(closest)) THEN
                              IF (PRESENT(chi)) THEN
                                  closest(1,ntotal)  = xr(1)
                                  closest(2,ntotal)  = xr(2)
                                  closest(3,ntotal)  = xr(3)
                              ELSE
                                  closest(1,ntotal)  = x + dx*xk(1)
                                  closest(2,ntotal)  = y + dy*xk(2)
                                  closest(3,ntotal)  = z + dz*xk(3)
                              ENDIF
                          ENDIF
                          !-----------------------------------------------------
                      ENDIF          ! cell near interface
                  ENDDO       ! mesh cell i
              ENDDO       ! mesh cell j 
          ENDDO       ! mesh cell k
#elif __DIM == __2D
          !---------------------------------------------------------------------
          !  Find mesh cells through which the zero level passes and
          !  that are close to the interface (to exclude the boundaries
          !  of the narrow band)
          !---------------------------------------------------------------------
          DO j=1,yhi
              jp1 = j+1
              DO i=1,xhi
                  ip1 = i+1
                  neg(1:8) = 0
                  lok = .TRUE.
                  IF (fdata(i  ,j  ,jsub) .LT. 0.0_MK) neg(1) = 1
                  IF (ABS(fdata(i  ,j  ,jsub)) .GT. thresh) lok = .FALSE.
                  IF (fdata(ip1,j  ,jsub) .LT. 0.0_MK) neg(2) = 1
                  IF (ABS(fdata(ip1,j  ,jsub)) .GT. thresh) lok = .FALSE.
                  IF (fdata(i  ,jp1,jsub) .LT. 0.0_MK) neg(3) = 1
                  IF (ABS(fdata(i  ,jp1,jsub)) .GT. thresh) lok = .FALSE.
                  IF (fdata(ip1,jp1,jsub) .LT. 0.0_MK) neg(4) = 1
                  IF (ABS(fdata(ip1,jp1,jsub)) .GT. thresh) lok = .FALSE.
                  sneg = neg(1)+neg(2)+neg(3)+neg(4)
                  IF ((sneg.GT.0).AND.(sneg.LT.4).AND.(lok)) THEN
                      !---------------------------------------------------------
                      !  Old value of the level function
                      !---------------------------------------------------------
                      DO n=-1,2
                          np2 = n+2
                          jpn = j+n
                          DO m=-1,2
                              psi(m+2,np2) = fdata(i+m,jpn,jsub)
                          ENDDO
                      ENDDO
                      !---------------------------------------------------------
                      !  Right-hand-side for bi-cubic interpolation
                      !  polynomial. The order of the method is
                      !  limited by the order of the finite
                      !  differences used here.
                      !---------------------------------------------------------
                      ! value at (i,j)
                      rhs(1) = psi(2,2)
                      ! value at (i+1,j)
                      rhs(2) = psi(3,2)
                      ! value at (i,j+1)
                      rhs(3) = psi(2,3)
                      ! value at (i+1,j+1)
                      rhs(4) = psi(3,3)
                      !---------------------------------------------------------
                      ! Dx at (i,j)
                      rhs(5) = facx*(psi(3,2)-psi(1,2))
                      ! Dx at (i+1,j)
                      rhs(6) = facx*(psi(4,2)-psi(2,2))
                      ! Dx at (i,j+1)
                      rhs(7) = facx*(psi(3,3)-psi(1,3))
                      ! Dx at (i+1,j+1)
                      rhs(8) = facx*(psi(4,3)-psi(2,3))
                      !---------------------------------------------------------
                      ! Dy at (i,j)
                      rhs(9) = facy*(psi(2,3)-psi(2,1))
                      ! Dy at (i+1,j)
                      rhs(10) = facy*(psi(3,3)-psi(3,1))
                      ! Dy at (i,j+1)
                      rhs(11) = facy*(psi(2,4)-psi(2,2))
                      ! Dy at (i+1,j+1)
                      rhs(12) = facy*(psi(3,4)-psi(3,2))
                      !---------------------------------------------------------
                      ! DxDy at (i,j)
                      rhs(13) = facxy*(psi(3,3)-psi(1,3)-psi(3,1)+psi(1,1))
                      ! DxDy at (i+1,j)
                      rhs(14) = facxy*(psi(4,3)-psi(2,3)-psi(4,1)+psi(2,1))
                      ! DxDy at (i,j+1)
                      rhs(15) = facxy*(psi(3,4)-psi(1,4)-psi(3,2)+psi(1,2))
                      ! DxDy at (i+1,j+1)
                      rhs(16) = facxy*(psi(4,4)-psi(2,4)-psi(4,2)+psi(2,2))
                      !---------------------------------------------------------
                      !  Check if needed
                      !---------------------------------------------------------
                      IF (ppm_debug .GT. 1) THEN
                          DO m=1,16
                              IF ((rhs(m) .NE. rhs(m)) .OR.    &
     &                            (ABS(rhs(m)) .GT. hsave)) THEN
                                  !---------------------------------------------
                                  !  rhs(m) is Inf or NaN or too big
                                  !---------------------------------------------
                                  CALL ppm_write(ppm_rank,   &
     &                                'ppm_gmm_kickoff',     &
     &          'WARNING: RHS not defined! Check support of input field!',info)
                              ENDIF
                          ENDDO
                      ENDIF
                      !---------------------------------------------------------
                      !  Solve for polynomial coefficients. coef will
                      !  contain the result. 
                      !---------------------------------------------------------
#include "ppm_gmm_lubksb.inc"
! CHECK SOLUTION
!                         U = 0.0_MK
!                         L = 0.0_MK
!                         LU = 0.0_MK
!                         DO ii=1,16
!                             DO jj=ii,16
!                                 U(ii,jj) = A(ii,jj)
!                             ENDDO
!                         ENDDO
!                         DO ii=1,16
!                             DO jj=1,(ii-1)
!                                 L(ii,jj) = A(ii,jj)
!                             ENDDO
!                         ENDDO
!                         DO ii=1,16
!                             L(ii,ii) = L(ii,ii) + 1.0_MK
!                         ENDDO
!                         LU = MATMUL(L,U)
!                         gaga = MATMUL(LU,coef)
!                         nrm2 = 0.0_MK
!                         DO ii=1,16
!                             ll = Aind(ii)
!                             nrm2 = nrm2 + (gaga(ii)-rhs(ll))**2
!                         ENDDO
!                         PRINT*,'ERROR: ',nrm2
! Wanna see them interpolants in MATLAB??? Uncomment this!
!                     WRITE(cbuf,'(A,I3.3,A)') 'interp_',npos/4,'.out'
!                     OPEN(30,FILE=cbuf,STATUS='REPLACE',ACTION='WRITE')
!                     WRITE(30,'(2E20.8)') topo%min_subd(1,isub,gmm_topoid)+(i-1)*dx, &
!    &                               topo%min_subd(2,isub,gmm_topoid)+(j-1)*dy
!                     WRITE(30,'(2E20.8)') topo%min_subd(1,isub,gmm_topoid)+(i)*dx, &
!    &                               topo%min_subd(2,isub,gmm_topoid)+(j-1)*dy
!                     WRITE(30,'(2E20.8)') topo%min_subd(1,isub,gmm_topoid)+(i-1)*dx, &
!    &                               topo%min_subd(2,isub,gmm_topoid)+(j)*dy
!                     WRITE(30,'(2E20.8)') topo%min_subd(1,isub,gmm_topoid)+(i)*dx, &
!    &                               topo%min_subd(2,isub,gmm_topoid)+(j)*dy
!                     WRITE(30,'(2E20.8)') rhs(1),rhs(2)
!                     WRITE(30,'(2E20.8)') rhs(3),rhs(4)
!                     WRITE(30,'(2E20.8)') coef(1:2)
!                     WRITE(30,'(2E20.8)') coef(3:4)
!                     WRITE(30,'(2E20.8)') coef(5:6)
!                     WRITE(30,'(2E20.8)') coef(7:8)
!                     WRITE(30,'(2E20.8)') coef(9:10)
!                     WRITE(30,'(2E20.8)') coef(11:12)
!                     WRITE(30,'(2E20.8)') coef(13:14)
!                     WRITE(30,'(2E20.8)') coef(15:16)
!                     CLOSE(30)
                      !---------------------------------------------------------
                      !  Check that sparse structure has at least 4
                      !  spaces to hold the new data
                      !---------------------------------------------------------
                      IF (gmm_lsiz .LT. npos+4) THEN
                          gmm_lsiz = gmm_lsiz + incr
                          iopt = ppm_param_alloc_grow_preserve
                          ldu(1) = gmm_lsiz
                          CALL ppm_alloc(phi,ldu,iopt,info)
                          IF (info .NE. 0) THEN
                              info = ppm_error_fatal
                              CALL ppm_error(ppm_err_alloc,   &
     &                            'ppm_gmm_kickoff',      &
     &                            'sparse values GMM_PHI',__LINE__,info)
                              GOTO 9999
                          ENDIF
                          ldu(1) = 3
                          ldu(2) = gmm_lsiz
                          CALL ppm_alloc(gmm_ipos,ldu,iopt,info)
                          IF (info .NE. 0) THEN
                              info = ppm_error_fatal
                              CALL ppm_error(ppm_err_alloc,   &
     &                            'ppm_gmm_kickoff',      &
     &                            'sparse positions GMM_IPOS',__LINE__,info)
                              GOTO 9999
                          ENDIF
                      ENDIF
                      IF (PRESENT(closest)) THEN
                          IF (ltotal .LT. ntotal+4) THEN
                              ltotal = ntotal + incr
                              iopt   = ppm_param_alloc_grow_preserve
                              ldu(1) = 2
                              ldu(2) = ltotal
                              CALL ppm_alloc(closest,ldu,iopt,info)
                              IF (info .NE. 0) THEN
                                  info = ppm_error_fatal
                                  CALL ppm_error(ppm_err_alloc,   &
     &                                'ppm_gmm_kickoff',          &
     &                                'closest points CLOSEST',   &
     &                                __LINE__,info)
                                  GOTO 9999
                              ENDIF
                          ENDIF
                      ENDIF
                      !---------------------------------------------------------
                      !  Newton iteration for each point bounding
                      !  the mesh cell to find the new level
                      !  function value and the closest point on the
                      !  interface.
                      !---------------------------------------------------------
                      IF (PRESENT(closest).AND..NOT.PRESENT(chi)) THEN
#if   __KIND == __SINGLE_PRECISION
                          x = topo%min_subs(1,isub)+(i-1)*dx
                          y = topo%min_subs(2,isub)+(j-1)*dy
#elif __KIND == __DOUBLE_PRECISION
                          x = topo%min_subd(1,isub)+(i-1)*dx
                          y = topo%min_subd(2,isub)+(j-1)*dy
#endif
                      ENDIF
                      x0    = 0.0_MK
                      IF (PRESENT(chi)) THEN
                          x = chi(1,i,j,jsub)
                          y = chi(2,i,j,jsub)
                      ENDIF
#include "ppm_gmm_quasinewton.inc"
                      ! BTW: the closest point on the interface is xk
                      IF (neg(1).EQ.1) sprod = -sprod
                      npos                   = npos + 1
                      ntotal                 = ntotal + 1
                      gmm_ipos(1,npos)       = i
                      gmm_ipos(2,npos)       = j
                      gmm_ipos(3,npos)       = jsub
                      phi(npos)              = sprod
                      ! the closest point on the interface is xk
                      IF (PRESENT(closest)) THEN
                          IF (PRESENT(chi)) THEN
                              closest(1,ntotal)  = xr(1)
                              closest(2,ntotal)  = xr(2)
                          ELSE
                              closest(1,ntotal)  = x + dx*xk(1)
                              closest(2,ntotal)  = y + dy*xk(2)
                          ENDIF
                      ENDIF
                      !---------------------------------------------------------
                      x0(1) = 1.0_MK
                      x0(2) = 0.0_MK
                      IF (PRESENT(chi)) THEN
                          x = chi(1,ip1,j,jsub)
                          y = chi(2,ip1,j,jsub)
                      ENDIF
#include "ppm_gmm_quasinewton.inc"
                      IF (neg(2).EQ.1) sprod = -sprod
                      npos                   = npos + 1
                      ntotal                 = ntotal + 1
                      gmm_ipos(1,npos)       = ip1
                      gmm_ipos(2,npos)       = j
                      gmm_ipos(3,npos)       = jsub
                      phi(npos)              = sprod
                      ! the closest point on the interface is xk
                      IF (PRESENT(closest)) THEN
                          IF (PRESENT(chi)) THEN
                              closest(1,ntotal)  = xr(1)
                              closest(2,ntotal)  = xr(2)
                          ELSE
                              closest(1,ntotal)  = x + dx*xk(1)
                              closest(2,ntotal)  = y + dy*xk(2)
                          ENDIF
                      ENDIF
                      !---------------------------------------------------------
                      x0(1) = 0.0_MK
                      x0(2) = 1.0_MK
                      IF (PRESENT(chi)) THEN
                          x = chi(1,i,jp1,jsub)
                          y = chi(2,i,jp1,jsub)
                      ENDIF
#include "ppm_gmm_quasinewton.inc"
                      IF (neg(3).EQ.1) sprod = -sprod
                      npos                   = npos + 1
                      ntotal                 = ntotal + 1
                      gmm_ipos(1,npos)       = i
                      gmm_ipos(2,npos)       = jp1
                      gmm_ipos(3,npos)       = jsub
                      phi(npos)              = sprod
                      ! the closest point on the interface is xk
                      IF (PRESENT(closest)) THEN
                          IF (PRESENT(chi)) THEN
                              closest(1,ntotal)  = xr(1)
                              closest(2,ntotal)  = xr(2)
                          ELSE
                              closest(1,ntotal)  = x + dx*xk(1)
                              closest(2,ntotal)  = y + dy*xk(2)
                          ENDIF
                      ENDIF
                      !---------------------------------------------------------
                      x0(1) = 1.0_MK
                      x0(2) = 1.0_MK
                      IF (PRESENT(chi)) THEN
                          x = chi(1,ip1,jp1,jsub)
                          y = chi(2,ip1,jp1,jsub)
                      ENDIF
#include "ppm_gmm_quasinewton.inc"
                      IF (neg(4).EQ.1) sprod = -sprod
                      npos                   = npos + 1
                      ntotal                 = ntotal + 1
                      gmm_ipos(1,npos)       = ip1
                      gmm_ipos(2,npos)       = jp1
                      gmm_ipos(3,npos)       = jsub
                      phi(npos)              = sprod
                      ! the closest point on the interface is xk
                      IF (PRESENT(closest)) THEN
                          IF (PRESENT(chi)) THEN
                              closest(1,ntotal)  = xr(1)
                              closest(2,ntotal)  = xr(2)
                          ELSE
                              closest(1,ntotal)  = x + dx*xk(1)
                              closest(2,ntotal)  = y + dy*xk(2)
                          ENDIF
                      ENDIF
                      !---------------------------------------------------------
                  ENDIF        ! mesh cell near interface
              ENDDO        ! mesh cell i
          ENDDO         ! mesh cell j
#endif
          !---------------------------------------------------------------------
          !  Check if we need to return the level function or not
          !---------------------------------------------------------------------
          lok = .TRUE.
          IF (PRESENT(npts)) THEN
              IF (npts .LT. 0) lok = .FALSE.
          ENDIF

          IF (lok) THEN
#if    __DIM == __3D
              !-----------------------------------------------------------------
              !  Flush field data incl. ghost layer
              !-----------------------------------------------------------------
              DO k=0,zhi+2   ! +2 because zhi=ndata-1 (see above)
                  DO j=0,yhi+2
                      DO i=0,xhi+2
                          fdata(i,j,k,jsub) = big
                      ENDDO
                  ENDDO
              ENDDO
              !-----------------------------------------------------------------
              !  Copy back from sparse structure
              !-----------------------------------------------------------------
              DO i=1,npos
                  m = gmm_ipos(1,i)
                  n = gmm_ipos(2,i)
                  p = gmm_ipos(3,i)
                  IF (ABS(phi(i)).LT.ABS(fdata(m,n,p,jsub)))    &
     &                fdata(m,n,p,jsub)=phi(i)
              ENDDO
#elif __DIM == __2D
              !-----------------------------------------------------------------
              !  Flush field data incl. ghost layer
              !-----------------------------------------------------------------
              DO j=0,yhi+2
                  DO i=0,xhi+2
                      fdata(i,j,jsub) = big
                  ENDDO
              ENDDO
              !-----------------------------------------------------------------
              !  Copy back from sparse structure
              !-----------------------------------------------------------------
              DO i=1,npos
                  m = gmm_ipos(1,i)
                  n = gmm_ipos(2,i)
                  IF (ABS(phi(i)).LT.ABS(fdata(m,n,jsub)))     &
     &                fdata(m,n,jsub)=phi(i)
              ENDDO
#endif
          ENDIF
          !---------------------------------------------------------------------
          !  Pass stuff back out if requested
          !---------------------------------------------------------------------
          IF (PRESENT(ipts)) THEN
              iopt = ppm_param_alloc_grow_preserve
              ldu(1) = ppm_dim+1
              ldu(2) = ntotal
              CALL ppm_alloc(ipts,ldu,iopt,info)
              IF (info .NE. ppm_param_success) THEN
                  info = ppm_error_fatal
                  CALL ppm_error(ppm_err_alloc,'ppm_gmm_kickoff',      &
     &                'point locations IPTS',__LINE__,info)
                  GOTO 9999
              ENDIF
              DO i=1,npos
                  j = i+ntotal-npos
#if   __DIM == __3D
                  ipts(1,j) = gmm_ipos(1,i)
                  ipts(2,j) = gmm_ipos(2,i)
                  ipts(3,j) = gmm_ipos(3,i)
                  ipts(4,j) = jsub
#elif __DIM == __2D
                  ipts(1,j) = gmm_ipos(1,i)
                  ipts(2,j) = gmm_ipos(2,i)
                  ipts(3,j) = jsub
#endif
              ENDDO
          ENDIF
      ENDDO             ! jsub
      !-------------------------------------------------------------------------
      !  Shrink to actual size to save memory
      !-------------------------------------------------------------------------
      IF (PRESENT(closest)) THEN
          iopt = ppm_param_alloc_fit_preserve
          ldu(1) = ppm_dim
          ldu(2) = ntotal
          CALL ppm_alloc(closest,ldu,iopt,info)
          IF (info .NE. ppm_param_success) THEN
              info = ppm_error_fatal
              CALL ppm_error(ppm_err_alloc,'ppm_gmm_kickoff',      &
     &            'closest points CLOSEST',__LINE__,info)
              GOTO 9999
          ENDIF
      ENDIF
      IF (PRESENT(npts)) npts = ntotal
      !-------------------------------------------------------------------------
      !  Nullify pointer to work memory and restore module pointers
      !-------------------------------------------------------------------------
#if    __KIND == __SINGLE_PRECISION
      gmm_phis => phi
#elif  __KIND == __DOUBLE_PRECISION
      gmm_phid => phi
#endif 
      NULLIFY(phi)
      !-------------------------------------------------------------------------
      !  Return 
      !-------------------------------------------------------------------------
 9999 CONTINUE
      CALL substop('ppm_gmm_kickoff',t0,info)
      RETURN

#if    __DIM == __2D
#if    __KIND == __SINGLE_PRECISION
      END SUBROUTINE ppm_gmm_kickoff_2ds
#elif  __KIND == __DOUBLE_PRECISION
      END SUBROUTINE ppm_gmm_kickoff_2dd
#endif 
#elif  __DIM == __3D
#if    __KIND == __SINGLE_PRECISION
      END SUBROUTINE ppm_gmm_kickoff_3ds
#elif  __KIND == __DOUBLE_PRECISION
      END SUBROUTINE ppm_gmm_kickoff_3dd
#endif 
#endif
