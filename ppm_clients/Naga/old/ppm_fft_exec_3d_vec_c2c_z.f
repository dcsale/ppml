      !-------------------------------------------------------------------------
      !  Subroutine   : ppm_fft_exec_3d_vec_c2c_z
      !-------------------------------------------------------------------------
      ! Copyright (c) 2010 CSE Lab (ETH Zurich), MOSAIC Group (ETH Zurich),
      !                    Center for Fluid Dynamics (DTU)
      !
      ! FFTW execute wrapper for 3d arrays, 1d complex to complex 
      ! (forward and backward) FFT in the z direction
      ! The routine does not work with fields that include ghost layers
      !-------------------------------------------------------------------------
#if __KIND == __SINGLE
#define __ROUTINE ppm_fft_exec_3d_vec_c2c_z_s
#define __PREC ppm_kind_single
#elif __KIND == __DOUBLE
#define __ROUTINE ppm_fft_exec_3d_vec_c2c_z_d
#define __PREC ppm_kind_double
#endif
      SUBROUTINE __ROUTINE(topoid,meshid,ppmplan,infield,outfield,info)
      !!! FFTW execute wrapper for 3d arrays, 1d complex to complex 
      !!! (forward and backward) FFT in the z direction
      !!! Before calling this routine a ppm_fft_plan_ routine must be called
      !!! The routine does not work with fields that include ghost layers
      USE ppm_module_substart
      USE ppm_module_substop
      USE ppm_module_typedef
      USE ppm_module_topo_get
      USE ppm_module_write
      USE ppm_module_data,ONLY:ppm_rank,ppm_kind_single,ppm_kind_double

      IMPLICIT NONE

      INCLUDE 'fftw3.f'

      ! if debug check if dimensions are 2a 3b 5c 7d 11e 13f
      !-------------------------------------------------------------------------
      ! Arguments
      !-------------------------------------------------------------------------
      !!!topology identifier of target
      INTEGER,INTENT(IN)                                             :: topoid
      !!!id of the mesh
      INTEGER,INTENT(IN)                                             :: meshid
      !!!ppm fft plan type
      TYPE(ppm_fft_plan),INTENT(INOUT)                               :: ppmplan
      !!!input field to fourier transform
      !COMPLEX(__PREC),DIMENSION(:,:,:,:,:),POINTER,INTENT(INOUT)     :: infield
      COMPLEX(__PREC),DIMENSION(:,:,:,:,:),POINTER                   :: infield
      !!!output field for the result of the fourier transform
      !COMPLEX(__PREC),DIMENSION(:,:,:,:,:),POINTER,INTENT(INOUT)     :: outfield
      COMPLEX(__PREC),DIMENSION(:,:,:,:,:),POINTER                   :: outfield
      !!!Returns status, 0 upon success
      INTEGER,INTENT(OUT)                                            :: info
      !in time perhaps an argument for alternate directions

      !-------------------------------------------------------------------------
      ! Local variables
      !-------------------------------------------------------------------------
      REAL(__PREC)                  :: t0
      INTEGER                       :: i,j
      INTEGER                       :: isub,isubl
      INTEGER                       :: nsubs
      INTEGER,DIMENSION(:),POINTER  :: isublist
      TYPE(ppm_t_topo),POINTER      :: topology
      TYPE(ppm_t_equi_mesh)         :: mesh

      !-------------------------------------------------------------------------
      ! Initialise routine
      !-------------------------------------------------------------------------
      CALL substart('ppm_fft_exec',t0,info)

      !-------------------------------------------------------------------------
      ! Get topology and mesh values
      !-------------------------------------------------------------------------
      CALL ppm_topo_get(topoid,topology,info)
      IF (info .NE. 0) THEN
         CALL ppm_write(ppm_rank,'ppm_fft_exec','Failed to get topology.',isub)
         GOTO 9999
      ENDIF
      nsubs = topology%nsublist
      ALLOCATE(isublist(nsubs))
      DO isub=1,nsubs
        isublist(isub) = topology%isublist(isub)
      ENDDO
      mesh  = topology%mesh(meshid)

      !outfield = 0.0_ppm_kind_double !@
      !-------------------------------------------------------------------------
      ! Execute plan
      !-------------------------------------------------------------------------
      DO isub=1,nsubs
         isubl=isublist(isub)
         DO j=1,mesh%nnodes(2,isubl) !@ add '-1' to exclude n+1 points
            DO i=1,mesh%nnodes(1,isubl) !@ add '-1' to exclude n+1 points
               CALL dfftw_execute_dft(ppmplan%plan(isub),&
               & infield(1,i,j,1,isub),outfield(1,i,j,1,isub))
            END DO
         END DO
      END DO

      !-------------------------------------------------------------------------
      ! Copy periodic - this is only for the periodic 'N+1' vertex points
      !-------------------------------------------------------------------------
      IF (topology%bcdef(3) .EQ. ppm_param_bcdef_periodic) THEN !vertex
         IF (ppmplan%sign .EQ. FFTW_BACKWARD) THEN
            DO isub=1,nsubs
               isubl=isublist(isub)
               DO j=1,mesh%nnodes(2,isubl) !@ add '-1' to exclude n+1 points
                  DO i=1,mesh%nnodes(1,isubl) !@ add '-1' to exclude n+1 points
                     outfield(1,i,j,mesh%nnodes(3,isubl),isub) = outfield(1,i,j,1,isub)
                     outfield(2,i,j,mesh%nnodes(3,isubl),isub) = outfield(2,i,j,1,isub)
                     outfield(3,i,j,mesh%nnodes(3,isubl),isub) = outfield(3,i,j,1,isub)
                  END DO
               END DO
            END DO
         ENDIF
      ENDIF

      !-------------------------------------------------------------------------
      ! Return
      !-------------------------------------------------------------------------
 9999 CONTINUE
      CALL substop('ppm_fft_exec',t0,info)
      RETURN

      END SUBROUTINE __ROUTINE
#undef __ROUTINE
#undef __PREC
