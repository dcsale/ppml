      !-------------------------------------------------------------------------
      !  Subroutine   :                ppm_fdsolver_map_3d
      !-------------------------------------------------------------------------
      !
      !  Purpose      : maps the data from topology topo_ids(1)) / mesh 
      !               (mesh_ids(1)) to topology (topo_ids(2))/mesh (mesh_ids(2))
      !                 
      !
      !  Input        :  
      !                 topo_ids(2)        (I)     
      !                                         first: current topology   
      !                                         second: destination topology   
      !                 mesh_ids(2)        (I)
      !                                         first: current mesh of the data
      !                                         second: destination mesh   
      !                                
      !
      !  Input/output : data_f(:,:,:,:,:)  (F)  data to be mapped
      !                                
      !                                
      !
      !  Output       : 
      !                 info          (I) return status. =0 if no error.
      !
      !  Remarks      : 
      !                                                  
      !  References   :
      !
      !  Revisions    :
      !-------------------------------------------------------------------------
      !  $Log: ppm_fdsolver_map_3d.f,v $
      !  Revision 1.1.1.1  2007/07/13 10:18:54  ivos
      !  CBL version of the PPM library
      !
      !  Revision 1.11  2006/09/04 18:34:44  pchatela
      !  Fixes and cleanups to make it compile for
      !  release-1-0
      !
      !  Revision 1.10  2005/05/03 13:33:10  hiebers
      !  Bugfix: call ppm_check_meshid with internal topo ID
      !
      !  Revision 1.9  2004/10/01 16:08:58  ivos
      !  Replaced REAL(ppm_kind_double) :: t0 with REAL(MK) t0.
      !
      !  Revision 1.8  2004/08/31 15:14:26  hiebers
      !  added argument checks for topo and mesh ids
      !
      !  Revision 1.7  2004/08/27 09:45:47  hiebers
      !  added ghostsize as an argument
      !
      !  Revision 1.4  2004/08/18 16:34:07  hiebers
      !  added scalar vector version
      !
      !  Revision 1.3  2004/07/26 13:49:17  ivos
      !  Removed Routines sections from the header comment.
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
#if   __DIM == __SFIELD
#if   __KIND == __SINGLE_PRECISION
      SUBROUTINE ppm_fdsolver_map_3d_sca_s(data_fv, topo_ids, mesh_ids,&
     &  ghostsize, info)
#elif __KIND == __DOUBLE_PRECISION
      SUBROUTINE ppm_fdsolver_map_3d_sca_d(data_fv, topo_ids, mesh_ids,&
     &  ghostsize, info)
#elif __KIND == __COMPLEX
      SUBROUTINE ppm_fdsolver_map_3d_sca_c(data_fv, topo_ids, mesh_ids,&
     &  ghostsize, info)
#elif __KIND == __DOUBLE_COMPLEX
      SUBROUTINE ppm_fdsolver_map_3d_sca_cc(data_fv, topo_ids, mesh_ids,&
     &  ghostsize, info)
#endif
#endif
#if   __DIM == __VFIELD
#if   __KIND == __SINGLE_PRECISION
      SUBROUTINE ppm_fdsolver_map_3d_vec_s(data_fv, lda, topo_ids, mesh_ids, &
     &  ghostsize, info)
#elif __KIND == __DOUBLE_PRECISION
      SUBROUTINE ppm_fdsolver_map_3d_vec_d(data_fv, lda, topo_ids, mesh_ids, &
     &  ghostsize, info)
#elif __KIND == __COMPLEX
      SUBROUTINE ppm_fdsolver_map_3d_vec_c(data_fv, lda, topo_ids, mesh_ids, &
     &  ghostsize, info)
#elif __KIND == __DOUBLE_COMPLEX
      SUBROUTINE ppm_fdsolver_map_3d_vec_cc(data_fv, lda, topo_ids, mesh_ids, &
     & ghostsize, info)
#endif
#endif
      !-------------------------------------------------------------------------
      !  Modules
      !-------------------------------------------------------------------------
      USE ppm_module_data
      USE ppm_module_map
      USE ppm_module_check_id
      USE ppm_module_error
      USE ppm_module_write
      USE ppm_module_substart
      USE ppm_module_substop
      USE ppm_module_map_field

      IMPLICIT NONE
#if   __KIND == __SINGLE_PRECISION | __KIND == __COMPLEX 
      INTEGER, PARAMETER :: MK = ppm_kind_single
#elif __KIND == __DOUBLE_PRECISION | __KIND == __DOUBLE_COMPLEX 
      INTEGER, PARAMETER :: MK = ppm_kind_double
#endif
      !-------------------------------------------------------------------------
      !  Arguments     
      !-------------------------------------------------------------------------
      ! POINTER to data
#if   __DIM == __SFIELD
#if   __KIND == __SINGLE_PRECISION | __KIND == __DOUBLE_PRECISION
      REAL(MK), DIMENSION(:,:,:,:),          POINTER   :: data_fv
#elif __KIND == __COMPLEX | __KIND == __DOUBLE_COMPLEX
      COMPLEX(MK), DIMENSION(:,:,:,:),       POINTER   :: data_fv
#endif
#endif
#if   __DIM == __VFIELD
#if   __KIND == __SINGLE_PRECISION | __KIND == __DOUBLE_PRECISION
      REAL(MK), DIMENSION(:,:,:,:,:),        POINTER   :: data_fv
#elif __KIND == __COMPLEX | __KIND == __DOUBLE_COMPLEX
      COMPLEX(MK), DIMENSION(:,:,:,:,:),     POINTER   :: data_fv
#endif
      INTEGER,                          INTENT(IN)     :: lda
#endif   
      INTEGER, DIMENSION(2),            INTENT(IN)     :: topo_ids
      INTEGER, DIMENSION(2),            INTENT(IN)     :: mesh_ids
      INTEGER, DIMENSION(3),            INTENT(IN)     :: ghostsize
      INTEGER              ,            INTENT(  OUT)  :: info
      !-------------------------------------------------------------------------
      !  Local variables 
      !-------------------------------------------------------------------------
      ! timer
      REAL(MK)                                :: t0
      ! counters
      INTEGER                                 :: k, i, j
      !Size of the data_in 
      INTEGER                                 :: from_topo, to_topo
      INTEGER                                 :: from_mesh, to_mesh
#if   __DIM == __SFIELD
      INTEGER, PARAMETER                      :: lda = 1
#endif   
      INTEGER                                 :: maptype
      LOGICAL                                 :: valid
      CHARACTER(LEN=ppm_char)                 :: mesg
      !-------------------------------------------------------------------------
      !  Externals 
      !-------------------------------------------------------------------------

      !-------------------------------------------------------------------------
      !  Initialise
      !-------------------------------------------------------------------------
      CALL substart('ppm_fdsolver_map_3d',t0,info)
      !-------------------------------------------------------------------------
      !  Check arguments
      !-------------------------------------------------------------------------
      IF (ppm_debug .GT. 0) THEN
          IF (topo_ids(1).GE.0) THEN
              CALL ppm_check_topoid(topo_ids(1),valid,info)
              IF (.NOT. valid) THEN
                  info = ppm_error_error
                  CALL ppm_error(ppm_err_argument,'ppm_fdsolver_map_2d',  &
     &                'Topology ID (from_topo) is invalid!',__LINE__,info)
                  GOTO 9999
              ENDIF
          ENDIF
          IF (topo_ids(2).GE.0) THEN
              CALL ppm_check_topoid(topo_ids(2),valid,info)
              IF (.NOT. valid) THEN
                  info = ppm_error_error
                  CALL ppm_error(ppm_err_argument,'ppm_fdsolver_map_2d',  &
     &                'Topology ID (to_topo) is invalid!',__LINE__,info)
                  GOTO 9999
              ENDIF           
          ENDIF
          IF (mesh_ids(1) .GT. 0) THEN
                CALL ppm_check_meshid(topo_ids(1),mesh_ids(1),valid,info)
                IF (.NOT. valid) THEN
                    info = ppm_error_error
                    CALL ppm_error(ppm_err_argument,'ppm_fdsolver_map_2d',  &
     &                  'Mesh ID (from_mesh) invalid!',__LINE__,info)
                    GOTO 9999
                ENDIF
          ENDIF
          IF (mesh_ids(2) .GT. 0) THEN
                CALL ppm_check_meshid(topo_ids(2),mesh_ids(2),valid,info)
                IF (.NOT. valid) THEN
                    info = ppm_error_error
                    CALL ppm_error(ppm_err_argument,'ppm_fdsolver_map_2d',  &
     &                  'Mesh ID (to_mesh) invalid!',__LINE__,info)
                    GOTO 9999
                ENDIF
          ENDIF
      ENDIF
      !-------------------------------------------------------------------------
      ! Define source and destination 
      !-------------------------------------------------------------------------
      from_topo = topo_ids(1)
      to_topo   = topo_ids(2)
      from_mesh = mesh_ids(1)
      to_mesh   = mesh_ids(2)
      IF (ppm_debug .GT. 0) THEN
         WRITE(mesg,'(A,I5,A,I5)' )'Mapping from topo ',from_topo,', mesh ',from_mesh 
         CALL ppm_write(ppm_rank,'ppm_fdsolver_map',mesg,j)      
         WRITE(mesg,'(A,I5,A,I5)' )'          to topo ',to_topo,', mesh ',to_mesh 
         CALL ppm_write(ppm_rank,'ppm_fdsolver_map',mesg,j)      
      ENDIF
      !-------------------------------------------------------------------------
      !  Map field 
      !-------------------------------------------------------------------------
#if   __DIM == __SFIELD
      CALL ppm_map_field_global(from_topo,to_topo,from_mesh,to_mesh,info)
      CALL ppm_map_field_push(from_topo,from_mesh,DATA_fv,lda,info)
      CALL ppm_map_field_send(info)
      CALL ppm_map_field_pop(to_topo,to_mesh,DATA_fv,lda,ghostsize,info)
#elif __DIM == __VFIELD
      CALL ppm_map_field_global(from_topo,to_topo,from_mesh,to_mesh,info)
      CALL ppm_map_field_push(from_topo,from_mesh,DATA_fv,lda,info)
      CALL ppm_map_field_send(info)
      CALL ppm_map_field_pop(to_topo,to_mesh,DATA_fv,lda,ghostsize,info)
#endif
      !-------------------------------------------------------------------------
      !  Return
      !-------------------------------------------------------------------------
 9999 CONTINUE
      CALL substop('ppm_fdsolver_map_3d',t0,info)
      RETURN

#if   __DIM == __SFIELD
#if   __KIND == __SINGLE_PRECISION
      END SUBROUTINE ppm_fdsolver_map_3d_sca_s
#elif __KIND == __DOUBLE_PRECISION
      END SUBROUTINE ppm_fdsolver_map_3d_sca_d
#elif __KIND == __COMPLEX
      END SUBROUTINE ppm_fdsolver_map_3d_sca_c
#elif __KIND == __DOUBLE_COMPLEX
      END SUBROUTINE ppm_fdsolver_map_3d_sca_cc
#endif
#endif
#if   __DIM == __VFIELD
#if   __KIND == __SINGLE_PRECISION
      END SUBROUTINE ppm_fdsolver_map_3d_vec_s
#elif __KIND == __DOUBLE_PRECISION
      END SUBROUTINE ppm_fdsolver_map_3d_vec_d
#elif __KIND == __COMPLEX
      END SUBROUTINE ppm_fdsolver_map_3d_vec_c
#elif __KIND == __DOUBLE_COMPLEX
      END SUBROUTINE ppm_fdsolver_map_3d_vec_cc
#endif
#endif
