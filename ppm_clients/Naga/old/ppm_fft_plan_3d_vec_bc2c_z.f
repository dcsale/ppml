      !-------------------------------------------------------------------------
      !  Subroutine   : ppm_fft_plan_3d_vec_bc2c_z
      !-------------------------------------------------------------------------
      ! Copyright (c) 2010 CSE Lab (ETH Zurich), MOSAIC Group (ETH Zurich),
      !                    Center for Fluid Dynamics (DTU)
      !
      ! FFTW plan wrapper for 3d arrays, 1d complex to complex
      ! (backward) FFT in the z direction
      ! The routine does not work with fields that include ghost layers
      !-------------------------------------------------------------------------
#if __KIND == __SINGLE
#define __ROUTINE ppm_fft_plan_3d_vec_bc2c_z_s
#define __PREC ppm_kind_single
#elif __KIND == __DOUBLE
#define __ROUTINE ppm_fft_plan_3d_vec_bc2c_z_d
#define __PREC ppm_kind_double
#endif
      SUBROUTINE __ROUTINE(topoid,meshid,ppmplan,infield,outfield,info)
      !!! FFTW plan wrapper for 3d arrays, 1d complex to complex
      !!! (backward) FFT in the z direction
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
      INTEGER                       :: isub,isubl
      INTEGER                       :: nsubs
      INTEGER,DIMENSION(:),POINTER  :: isublist
      TYPE(ppm_t_topo),POINTER      :: topology
      TYPE(ppm_t_equi_mesh)         :: mesh

      !-------------------------------------------------------------------------
      ! Initialise routine
      !-------------------------------------------------------------------------
      CALL substart('ppm_fft_plan',t0,info)

      !-------------------------------------------------------------------------
      ! Get topology and mesh values
      !-------------------------------------------------------------------------
      CALL ppm_topo_get(topoid,topology,info)
      IF (info .NE. 0) THEN
         CALL ppm_write(ppm_rank,'ppm_fft_plan','Failed to get topology.',isub)
         GOTO 9999
      ENDIF
      nsubs = topology%nsublist
      ALLOCATE(isublist(nsubs))
      DO isub=1,nsubs
        isublist(isub) = topology%isublist(isub)
      ENDDO
      mesh  = topology%mesh(meshid)

      !-------------------------------------------------------------------------
      ! Setup parameters for this particular routine
      !-------------------------------------------------------------------------
      !the dimension of the FFT (1D/2D/3D)
      ppmplan%rank=1
      !the number of points along each direction of the piece to be transformed
      ALLOCATE(ppmplan%nx(ppmplan%rank,nsubs))
      !the direction of the transform
      ppmplan%sign=FFTW_BACKWARD
      !the method to setup the optimal plan
      ppmplan%flag=FFTW_MEASURE
      !the number of components to transform - 3 component vector
      ppmplan%howmany=3
      !the size of the input array - full size (assuming LBOUND=1 thus UBOUND)
      ALLOCATE(ppmplan%inembed(ppmplan%rank))
      ppmplan%inembed(1) = UBOUND(infield,4)
      !the size of the output array - full size (assuming LBOUND=1 thus UBOUND)
      ALLOCATE(ppmplan%onembed(ppmplan%rank))
      ppmplan%onembed(1) = UBOUND(outfield,4)
      !istride tells how the same componenet data points are spaced in memory
      !e.g. z values recur every x-dim*y-dim*component values
      ppmplan%istride = UBOUND(infield,2) *UBOUND(infield,3)*3
      ppmplan%ostride = UBOUND(outfield,2)*UBOUND(outfield,3)*3
      !idist tells how multiple arrays are spaced in memory. I.e. a memory 
      !offset. e.g. vector components (idist=1) or scalar 2D arrays in 
      !3D array(idist=NxNy)
      ppmplan%idist = 1
      ppmplan%odist = 1

      !-------------------------------------------------------------------------
      ! Allocate plan array
      !-------------------------------------------------------------------------
      IF(ASSOCIATED(ppmplan%plan)) THEN
         DEALLOCATE(ppmplan%plan,stat=info)
         IF (info .NE. 0) THEN
            CALL ppm_write(ppm_rank,'ppm_fft_plan','Failed to deallocate plan-array.',isub)
            GOTO 9999
         ENDIF
      END IF
      ALLOCATE(ppmplan%plan(nsubs))

      DO isub=1,nsubs
         isubl=isublist(isub)
         !@ maybe the -1 needs to be removed when doing cell data
         !we subtract the -1 to avoid the periodic vertex point
         IF (topology%bcdef(3) .EQ. ppm_param_bcdef_periodic) THEN !vertex
            ppmplan%nx(1,isub) = mesh%nm(3)-1
         ELSE
            ppmplan%nx(1,isub) = mesh%nm(3)
         ENDIF

         CALL dfftw_plan_many_dft(ppmplan%plan(isub),ppmplan%rank,&
         & ppmplan%nx(:,isub),ppmplan%howmany,infield(1,1,1,1,isub),&
         & ppmplan%inembed(1),ppmplan%istride,ppmplan%idist,&
         & outfield(1,1,1,1,isub),ppmplan%onembed(1),ppmplan%ostride,&
         & ppmplan%odist,ppmplan%sign,ppmplan%flag)
      END DO


      !-------------------------------------------------------------------------
      ! Return
      !-------------------------------------------------------------------------
 9999 CONTINUE
      CALL substop('ppm_fft_plan',t0,info)
      RETURN

      END SUBROUTINE __ROUTINE 
#undef __ROUTINE
#undef __PREC
