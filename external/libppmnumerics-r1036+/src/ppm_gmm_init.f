      !-------------------------------------------------------------------------
      !  Subroutine   :                   ppm_gmm_init
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
      SUBROUTINE ppm_gmm_init(field_topoid,meshid,Nest,prec,info)
      !!! This routine initializes the ppm_gmm module and allocates all data
      !!! structures.
      !-------------------------------------------------------------------------
      !  Includes
      !-------------------------------------------------------------------------
#include "ppm_define.h"
      !-------------------------------------------------------------------------
      !  Modules 
      !-------------------------------------------------------------------------
      USE ppm_module_data
      USE ppm_module_data_gmm
      
      USE ppm_module_error
      USE ppm_module_typedef
      USE ppm_module_check_id
      USE ppm_module_alloc
      USE ppm_module_substart
      USE ppm_module_substop

      IMPLICIT NONE
      !-------------------------------------------------------------------------
      !  Arguments     
      !-------------------------------------------------------------------------
      INTEGER, INTENT(IN   ) :: Nest
      !!! Estimated number of grid points adjacent to the zero level. This is
      !!! only used to estimate the size of the memory needed and will be grown
      !!! automatically if too small.
      INTEGER, INTENT(IN   ) :: prec
      !!! Precision of the field data to be reinitialized. One of:
      !!!
      !!! *ppm_kind_single
      !!! *ppm_kind_double
      INTEGER, INTENT(IN   ) :: meshid
      !!! Mesh ID (user numbering) for which a GMM should be initialized.
      INTEGER, INTENT(IN   ) :: field_topoid
      !!! Topo ID of the field
      INTEGER, INTENT(  OUT) :: info
      !!! Return status. 0 upon success
      !-------------------------------------------------------------------------
      !  Local variables 
      !-------------------------------------------------------------------------
      INTEGER, DIMENSION(3)            :: ldu
      INTEGER                          :: iopt,i,isub
      LOGICAL                          :: lok
      REAL(ppm_kind_double)            :: t0
      TYPE(ppm_t_topo),      POINTER   :: topo
      TYPE(ppm_t_equi_mesh), POINTER   :: mesh
      !-------------------------------------------------------------------------
      !  Externals 
      !-------------------------------------------------------------------------
      
      !-------------------------------------------------------------------------
      !  Initialise 
      !-------------------------------------------------------------------------
      CALL substart('ppm_gmm_init',t0,info)
      topo => ppm_topo(field_topoid)%t
      mesh => topo%mesh(gmm_meshid)
      !-------------------------------------------------------------------------
      !  Check arguments
      !-------------------------------------------------------------------------
      IF (ppm_debug .GT. 0) THEN
          IF (.NOT. ppm_initialized) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_ppm_noinit,'ppm_gmm_init',  &
     &            'Please call ppm_init first!',__LINE__,info)
              GOTO 9999
          ENDIF
          CALL ppm_check_meshid(field_topoid,meshid,lok,info)
          IF (.NOT. lok) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_gmm_init',  &
     &            'meshid out of range',__LINE__,info)
              GOTO 9999
          ENDIF
          IF (Nest .LT. 1) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_gmm_init',  &
     &            'Nest must be >0!',__LINE__,info)
              GOTO 9999
          ENDIF
          IF ((prec.NE.ppm_kind_single).AND.(prec.NE.ppm_kind_double)) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_gmm_init',  &
     &            'Illegal precision specifiec!',__LINE__,info)
              GOTO 9999
          ENDIF
      ENDIF
      !-------------------------------------------------------------------------
      !  Nullify the pointers to allow proper inquiry of ASSOCIATED status.
      !-------------------------------------------------------------------------
      NULLIFY(gmm_ipos)
      NULLIFY(gmm_phis)
      NULLIFY(gmm_phid)
      NULLIFY(gmm_state2d)
      NULLIFY(gmm_state3d)
      !-------------------------------------------------------------------------
      !  Allocate sparse work space structure
      !-------------------------------------------------------------------------
      gmm_lsiz = Nest
      iopt = ppm_param_alloc_fit
      ldu(1) = gmm_lsiz
      IF (prec .EQ. ppm_kind_double) THEN
          CALL ppm_alloc(gmm_phid,ldu,iopt,info)
      ELSE
          CALL ppm_alloc(gmm_phis,ldu,iopt,info)
      ENDIF
      IF (info .NE. 0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_gmm_init',     &
     &        'sparse data values GMM_PHI',__LINE__,info)
          GOTO 9999
      ENDIF
      ldu(1) = ppm_dim + 1
      ldu(2) = gmm_lsiz
      CALL ppm_alloc(gmm_ipos,ldu,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_gmm_init',     &
     &        'sparse data locations GMM_IPOS',__LINE__,info)
          GOTO 9999
      ENDIF
      !-------------------------------------------------------------------------
      !  Translate meshid to internal numbering and store it
      !-------------------------------------------------------------------------

      gmm_meshid = mesh%ID
      !-------------------------------------------------------------------------
      !  Determine max extent of mesh in any sub
      !-------------------------------------------------------------------------
      maxxhi = 0
      maxyhi = 0
      maxzhi = 0
      DO i=1,topo%nsublist
          isub = topo%isublist(i)
          IF (mesh%nnodes(1,isub).GT.maxxhi) &
     &        maxxhi = mesh%nnodes(1,isub)
          IF (mesh%nnodes(2,isub).GT.maxyhi) &
     &        maxyhi = mesh%nnodes(2,isub)
          IF (ppm_dim .GT. 2) THEN 
             IF (mesh%nnodes(3,isub).GT.maxzhi)&
     &           maxzhi = mesh%nnodes(3,isub)
          ENDIF
      ENDDO
      !-------------------------------------------------------------------------
      !  Memory increment step size 
      !-------------------------------------------------------------------------
      incr = MAX(maxxhi,maxyhi,maxzhi)
      incr = 10*incr*incr
      !-------------------------------------------------------------------------
      !  Return 
      !-------------------------------------------------------------------------
 9999 CONTINUE
      CALL substop('ppm_gmm_init',t0,info)
      RETURN

      END SUBROUTINE ppm_gmm_init
