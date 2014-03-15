      !-------------------------------------------------------------------------
      ! Subroutine : ppm_module_fft
      !-------------------------------------------------------------------------
      ! Copyright (c) 2010 CSE Lab (ETH Zurich), MOSAIC Group (ETH Zurich),
      ! Center for Fluid Dynamics (DTU)
      !
      ! This modules contains routines for creating FFTW plans and executing them
      ! Plan routines are required for
      ! z means fft in x (1D)
      ! xy means fft in xy (2D)
      ! xyz means fft in xyz (3D)
      ! sca means scalar
      ! vec means vector
      ! f means forward fft
      ! b means backward fft
      ! c means complex
      ! r means real
      ! s means single precision
      ! d means double precision
      !
      ! e.g. 3d_vec_fr2c_xy_d: 3D vector array forward transform from
      ! real to complex and ffts in x and y directions
      ! (2d). For double data
      !
      ! The following variants have been implemented
      ! 3d_vec_fc2c_z_d
      ! 3d_vec_bc2c_z_d
      ! 3d_vec_fr2c_xy_d
      ! 3d_vec_bc2r_xy_d
      !
      ! a normalization routine exists for debugging purposes
      !
      ! sofar all fftw calls are to double precision routines!
      !
      ! All IFFT routines should be called with the 'realest' topoid/meshid
      !
      ! The routines respects the periodic N+1 points periodic BC cf topo%bcdef
      ! but also does full domain FFTs for freespace BC
      ! It does not work on mixed periodic/freespace BC in the XY direction
      ! topo%bcdef(1) is assumed to be in x, (2) to be in y, (3) z, (4) x...
      !-------------------------------------------------------------------------
      MODULE ppm_module_fft
      !-------------------------------------------------------------------------
      ! 1D transforms
      !-------------------------------------------------------------------------
      INTERFACE ppm_fft_forward_1d
         MODULE PROCEDURE ppm_fft_plan_3d_vec_fc2c_z_s
         MODULE PROCEDURE ppm_fft_plan_3d_vec_fc2c_z_d
      END INTERFACE
      INTERFACE ppm_fft_backward_1d
         MODULE PROCEDURE ppm_fft_plan_3d_vec_bc2c_z_s
         MODULE PROCEDURE ppm_fft_plan_3d_vec_bc2c_z_d
      END INTERFACE
      INTERFACE ppm_fft_execute_1d
         MODULE PROCEDURE ppm_fft_exec_3d_vec_c2c_z_s
         MODULE PROCEDURE ppm_fft_exec_3d_vec_c2c_z_d
      END INTERFACE
      !-------------------------------------------------------------------------
      ! 2D transforms
      !-------------------------------------------------------------------------
      INTERFACE ppm_fft_forward_2d
         MODULE PROCEDURE ppm_fft_plan_3d_vec_fr2c_xy_s
         MODULE PROCEDURE ppm_fft_plan_3d_vec_fr2c_xy_d
      END INTERFACE
      INTERFACE ppm_fft_backward_2d
         MODULE PROCEDURE ppm_fft_plan_3d_vec_bc2r_xy_s
         MODULE PROCEDURE ppm_fft_plan_3d_vec_bc2r_xy_d
         !MODULE PROCEDURE ppm_fft_plan_3d_vec_bc2c_xy_s
         !MODULE PROCEDURE ppm_fft_plan_3d_vec_bc2c_xy_d
      END INTERFACE
      INTERFACE ppm_fft_execute_2d
         MODULE PROCEDURE ppm_fft_exec_3d_vec_fr2c_xy_s
         MODULE PROCEDURE ppm_fft_exec_3d_vec_bc2r_xy_s
         MODULE PROCEDURE ppm_fft_exec_3d_vec_fr2c_xy_d
         MODULE PROCEDURE ppm_fft_exec_3d_vec_bc2r_xy_d
         !MODULE PROCEDURE ppm_fft_exec_3d_vec_bc2c_xy_s
         !MODULE PROCEDURE ppm_fft_exec_3d_vec_bc2c_xy_d
      END INTERFACE
      !-------------------------------------------------------------------------
      ! Normalization
      !-------------------------------------------------------------------------
      INTERFACE ppm_fft_normalize
         MODULE PROCEDURE ppm_fft_normalize_rs
         MODULE PROCEDURE ppm_fft_normalize_cs
         MODULE PROCEDURE ppm_fft_normalize_rd
         MODULE PROCEDURE ppm_fft_normalize_cd
      END INTERFACE
      !-------------------------------------------------------------------------
      ! PPM FFT plan type
      !-------------------------------------------------------------------------
      !!! Type containing the FFTW plan and its settings
      TYPE ppm_fft_plan
         !!!array of plan pointers, index for subs
         INTEGER*8,DIMENSION(:),POINTER :: plan => NULL()
         !!!the dimension of the FFT (1D/2D/3D)
         INTEGER :: rank
         !!!number of points along each direction of the piece to be transformed
         !!!index is for rank and subs
         INTEGER,DIMENSION(:,:),POINTER :: nx => NULL()
         !!!the direction of the transform (forward/backward)
         INTEGER :: sign
         !!!the method to setup the optimal plan
         INTEGER :: flag
         !!!the number of components to transform
         INTEGER :: howmany
         !!!the size of the input array, index is for rank
         INTEGER,DIMENSION(:),POINTER :: inembed => NULL()
         !!!the size of the output array, index is for rank
         INTEGER,DIMENSION(:),POINTER :: onembed => NULL()
         !!!istride tells how same component data points are spaced in memory
         INTEGER :: istride
         INTEGER :: ostride
         !!!idist tells how multiple arrays are spaced. I.e. a memory offset
         INTEGER :: idist
         INTEGER :: odist
      END TYPE ppm_fft_plan
      CONTAINS
      !FORWARD TRANSFORMS - PLAN
      !-------------------------------------------------------------------------
      ! Subroutine : ppm_fft_plan_3d_vec_fc2c_z
      !-------------------------------------------------------------------------
      ! Copyright (c) 2010 CSE Lab (ETH Zurich), MOSAIC Group (ETH Zurich),
      ! Center for Fluid Dynamics (DTU)
      !
      ! FFTW plan wrapper for 3d arrays, 1d complex to complex
      ! (forward) FFT in the xy directions
      ! The routine does not work with fields that include ghost layers
      !-------------------------------------------------------------------------
      SUBROUTINE ppm_fft_plan_3d_vec_fc2c_z_s(topoid,meshid,ppmplan,infield,outfield,info)
      !!! FFTW plan wrapper for 3d arrays, 1d complex to complex
      !!! (forward) FFT in the xy directions
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
      INTEGER,INTENT(IN) :: topoid
      !!!id of the mesh
      INTEGER,INTENT(IN) :: meshid
      !!!ppm fft plan type
      TYPE(ppm_fft_plan),INTENT(INOUT) :: ppmplan
      !!!input field to fourier transform
      !COMPLEX(ppm_kind_single),DIMENSION(:,:,:,:,:),POINTER,INTENT(INOUT) :: infield
      COMPLEX(ppm_kind_single),DIMENSION(:,:,:,:,:),POINTER :: infield
      !!!output field for the result of the fourier transform
      !COMPLEX(ppm_kind_single),DIMENSION(:,:,:,:,:),POINTER,INTENT(INOUT) :: outfield
      COMPLEX(ppm_kind_single),DIMENSION(:,:,:,:,:),POINTER :: outfield
      !!!Returns status, 0 upon success
      INTEGER,INTENT(OUT) :: info
      !in time perhaps an argument for alternate directions
      !-------------------------------------------------------------------------
      ! Local variables
      !-------------------------------------------------------------------------
      REAL(ppm_kind_single) :: t0
      INTEGER :: isub,isubl
      INTEGER :: nsubs
      INTEGER,DIMENSION(:),POINTER :: isublist
      TYPE(ppm_t_topo),POINTER :: topology
      TYPE(ppm_t_equi_mesh) :: mesh
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
      mesh = topology%mesh(meshid)
      !-------------------------------------------------------------------------
      ! Setup parameters for this particular routine
      !-------------------------------------------------------------------------
      !the dimension of the FFT (1D/2D/3D)
      ppmplan%rank=1
      !the number of points along each direction of the piece to be transformed
      ALLOCATE(ppmplan%nx(ppmplan%rank,nsubs))
      !the direction of the transform
      ppmplan%sign=FFTW_FORWARD
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
      !e.g. for 2/3 component vector istride = 2/3 or for scalar istride = 1
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
      END SUBROUTINE ppm_fft_plan_3d_vec_fc2c_z_s
      !-------------------------------------------------------------------------
      ! Subroutine : ppm_fft_plan_3d_vec_fr2c_xy
      !-------------------------------------------------------------------------
      ! Copyright (c) 2010 CSE Lab (ETH Zurich), MOSAIC Group (ETH Zurich),
      ! Center for Fluid Dynamics (DTU)
      !
      ! FFTW plan wrapper for 3d arrays, 2d real to complex
      ! (forward) FFT in the xy directions
      ! The routine does not work with fields that include ghost layers
      !-------------------------------------------------------------------------
      SUBROUTINE ppm_fft_plan_3d_vec_fr2c_xy_s(topoid,meshid,ppmplan,infield,outfield,info)
      !!! FFTW plan wrapper for 3d arrays, 2d real to complex
      !!! (forward) FFT in the xy directions
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
      INTEGER,INTENT(IN) :: topoid
      !!!id of the mesh
      INTEGER,INTENT(IN) :: meshid
      !!!ppm fft plan type
      TYPE(ppm_fft_plan),INTENT(INOUT) :: ppmplan
      !!!input field to fourier transform
      !REAL(ppm_kind_single),DIMENSION(:,:,:,:,:),POINTER,INTENT(INOUT) :: infield
      REAL(ppm_kind_single),DIMENSION(:,:,:,:,:),POINTER :: infield
      !!!output field for the result of the fourier transform
      !COMPLEX(ppm_kind_single),DIMENSION(:,:,:,:,:),POINTER,INTENT(INOUT) :: outfield
      COMPLEX(ppm_kind_single),DIMENSION(:,:,:,:,:),POINTER :: outfield
      !!!Returns status, 0 upon success
      INTEGER,INTENT(OUT) :: info
      !in time perhaps an argument for alternate directions
      !-------------------------------------------------------------------------
      ! Local variables
      !-------------------------------------------------------------------------
      REAL(ppm_kind_single) :: t0
      INTEGER :: isub,isubl
      INTEGER :: nsubs
      INTEGER,DIMENSION(:),POINTER :: isublist
      TYPE(ppm_t_topo),POINTER :: topology
      TYPE(ppm_t_equi_mesh) :: mesh
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
      mesh = topology%mesh(meshid)
      !-------------------------------------------------------------------------
      ! Setup parameters for this particular routine
      !-------------------------------------------------------------------------
      !the dimension of the FFT (1D/2D/3D)
      ppmplan%rank=2
      !the number of points along each direction of the piece to be transformed
      ALLOCATE(ppmplan%nx(ppmplan%rank,nsubs))
      !the direction of the transform
      ppmplan%sign=FFTW_FORWARD
      !the method to setup the optimal plan
      ppmplan%flag=FFTW_MEASURE
      !the number of components to transform - 3 component vector
      ppmplan%howmany=3
      !the size of the input array - full size (assuming LBOUND=1 thus UBOUND)
      ALLOCATE(ppmplan%inembed(ppmplan%rank))
      ppmplan%inembed(1) = UBOUND(infield,2)
      ppmplan%inembed(2) = UBOUND(infield,3)
      !the size of the output array - full size (assuming LBOUND=1 thus UBOUND)
      ALLOCATE(ppmplan%onembed(ppmplan%rank))
      ppmplan%onembed(1) = UBOUND(outfield,2)
      ppmplan%onembed(2) = UBOUND(outfield,3)
      !istride tells how the same componenet data points are spaced in memory
      !e.g. for 2/3 component vector istride = 2/3 or for scalar istride = 1
      ppmplan%istride = 3
      ppmplan%ostride = 3
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
         IF (topology%bcdef(1) .EQ. ppm_param_bcdef_periodic) THEN !vertex
            ppmplan%nx(1,isub) = mesh%nm(1)-1
            ppmplan%nx(2,isub) = mesh%nm(2)-1
         ELSE
            ppmplan%nx(1,isub) = mesh%nm(1)
            ppmplan%nx(2,isub) = mesh%nm(2)
         ENDIF
         CALL dfftw_plan_many_dft_r2c(ppmplan%plan(isub),ppmplan%rank,&
         & ppmplan%nx(:,isub),ppmplan%howmany,infield(1,1,1,1,isub),&
         & ppmplan%inembed(1),ppmplan%istride,ppmplan%idist,&
         & outfield(1,1,1,1,isub),ppmplan%onembed(1),ppmplan%ostride,&
         & ppmplan%odist,ppmplan%flag)
      END DO
      !-------------------------------------------------------------------------
      ! Return
      !-------------------------------------------------------------------------
 9999 CONTINUE
      CALL substop('ppm_fft_plan',t0,info)
      RETURN
      END SUBROUTINE ppm_fft_plan_3d_vec_fr2c_xy_s
      !BACKWARD TRANSFORMS - PLAN
      !-------------------------------------------------------------------------
      ! Subroutine : ppm_fft_plan_3d_vec_bc2c_z
      !-------------------------------------------------------------------------
      ! Copyright (c) 2010 CSE Lab (ETH Zurich), MOSAIC Group (ETH Zurich),
      ! Center for Fluid Dynamics (DTU)
      !
      ! FFTW plan wrapper for 3d arrays, 1d complex to complex
      ! (backward) FFT in the z direction
      ! The routine does not work with fields that include ghost layers
      !-------------------------------------------------------------------------
      SUBROUTINE ppm_fft_plan_3d_vec_bc2c_z_s(topoid,meshid,ppmplan,infield,outfield,info)
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
      INTEGER,INTENT(IN) :: topoid
      !!!id of the mesh
      INTEGER,INTENT(IN) :: meshid
      !!!ppm fft plan type
      TYPE(ppm_fft_plan),INTENT(INOUT) :: ppmplan
      !!!input field to fourier transform
      !COMPLEX(ppm_kind_single),DIMENSION(:,:,:,:,:),POINTER,INTENT(INOUT) :: infield
      COMPLEX(ppm_kind_single),DIMENSION(:,:,:,:,:),POINTER :: infield
      !!!output field for the result of the fourier transform
      !COMPLEX(ppm_kind_single),DIMENSION(:,:,:,:,:),POINTER,INTENT(INOUT) :: outfield
      COMPLEX(ppm_kind_single),DIMENSION(:,:,:,:,:),POINTER :: outfield
      !!!Returns status, 0 upon success
      INTEGER,INTENT(OUT) :: info
      !in time perhaps an argument for alternate directions
      !-------------------------------------------------------------------------
      ! Local variables
      !-------------------------------------------------------------------------
      REAL(ppm_kind_single) :: t0
      INTEGER :: isub,isubl
      INTEGER :: nsubs
      INTEGER,DIMENSION(:),POINTER :: isublist
      TYPE(ppm_t_topo),POINTER :: topology
      TYPE(ppm_t_equi_mesh) :: mesh
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
      mesh = topology%mesh(meshid)
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
      END SUBROUTINE ppm_fft_plan_3d_vec_bc2c_z_s
      !-------------------------------------------------------------------------
      ! Subroutine : ppm_fft_plan_3d_vec_bc2r_xy
      !-------------------------------------------------------------------------
      ! Copyright (c) 2010 CSE Lab (ETH Zurich), MOSAIC Group (ETH Zurich),
      ! Center for Fluid Dynamics (DTU)
      !
      ! FFTW plan wrapper for 3d arrays, 2d complex to real
      ! (backward) FFT in the xy directions
      ! The routine does not work with fields that include ghost layers
      !
      ! For the C2R inverse transform FFTW always destroys the input data
      !-------------------------------------------------------------------------
      SUBROUTINE ppm_fft_plan_3d_vec_bc2r_xy_s(topoid,meshid,ppmplan,infield,outfield,info)
      !!! FFTW plan wrapper for 3d arrays, 2d complex to real
      !!! (backward) FFT in the xy directions
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
      INTEGER,INTENT(IN) :: topoid
      !!!id of the mesh
      INTEGER,INTENT(IN) :: meshid
      !!!ppm fft plan type
      TYPE(ppm_fft_plan),INTENT(INOUT) :: ppmplan
      !!!input field to fourier transform
      !COMPLEX(ppm_kind_single),DIMENSION(:,:,:,:,:),POINTER,INTENT(INOUT) :: infield
      COMPLEX(ppm_kind_single),DIMENSION(:,:,:,:,:),POINTER :: infield
      !!!output field for the result of the fourier transform
      !REAL(ppm_kind_single),DIMENSION(:,:,:,:,:),POINTER,INTENT(INOUT) :: outfield
      REAL(ppm_kind_single),DIMENSION(:,:,:,:,:),POINTER :: outfield
      !!!Returns status, 0 upon success
      INTEGER,INTENT(OUT) :: info
      !in time perhaps an argument for alternate directions
      !-------------------------------------------------------------------------
      ! Local variables
      !-------------------------------------------------------------------------
      REAL(ppm_kind_single) :: t0
      INTEGER :: isub,isubl
      INTEGER :: nsubs
      INTEGER,DIMENSION(:),POINTER :: isublist
      TYPE(ppm_t_topo),POINTER :: topology
      TYPE(ppm_t_equi_mesh) :: mesh
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
      mesh = topology%mesh(meshid)
      !-------------------------------------------------------------------------
      ! Setup parameters for this particular routine
      !-------------------------------------------------------------------------
      !the dimension of the FFT (1D/2D/3D)
      ppmplan%rank=2
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
      ppmplan%inembed(1) = UBOUND(infield,2)
      ppmplan%inembed(2) = UBOUND(infield,3)
      !the size of the output array - full size (assuming LBOUND=1 thus UBOUND)
      ALLOCATE(ppmplan%onembed(ppmplan%rank))
      ppmplan%onembed(1) = UBOUND(outfield,2)
      ppmplan%onembed(2) = UBOUND(outfield,3)
      !istride tells how the same componenet data points are spaced in memory
      !e.g. for 2/3 component vector istride = 2/3 or for scalar istride = 1
      ppmplan%istride = 3
      ppmplan%ostride = 3
      !idist tells how multiple arrays are spaced in memory. I.e. a memory
      !offset. e.g. vector components (idist=1) or scalar 2D arrays
      !in 3D array(idist=NxNy)
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
         !the number of points in each direction of the piece to be transformed
         !we subtract the -1 to avoid the periodic vertex point
         IF (topology%bcdef(1) .EQ. ppm_param_bcdef_periodic) THEN !vertex
            ppmplan%nx(1,isub) = mesh%nm(1)-1
            ppmplan%nx(2,isub) = mesh%nm(2)-1
         ELSE
            ppmplan%nx(1,isub) = mesh%nm(1)
            ppmplan%nx(2,isub) = mesh%nm(2)
         ENDIF
         CALL dfftw_plan_many_dft_c2r(ppmplan%plan(isub),ppmplan%rank,&
         & ppmplan%nx,ppmplan%howmany,infield(1,1,1,1,isub),ppmplan%inembed(1),&
         & ppmplan%istride,ppmplan%idist,outfield(1,1,1,1,isub),&
         & ppmplan%onembed(1),ppmplan%ostride,ppmplan%odist,ppmplan%flag)
      END DO
      !-------------------------------------------------------------------------
      ! Return
      !-------------------------------------------------------------------------
 9999 CONTINUE
      CALL substop('ppm_fft_plan',t0,info)
      RETURN
      END SUBROUTINE ppm_fft_plan_3d_vec_bc2r_xy_s
!#include "ppm_fft_plan_3d_vec_bc2c_xy.f"
      !EXECUTION OF PLANS
      !-------------------------------------------------------------------------
      ! Subroutine : ppm_fft_exec_3d_vec_c2c_z
      !-------------------------------------------------------------------------
      ! Copyright (c) 2010 CSE Lab (ETH Zurich), MOSAIC Group (ETH Zurich),
      ! Center for Fluid Dynamics (DTU)
      !
      ! FFTW execute wrapper for 3d arrays, 1d complex to complex
      ! (forward and backward) FFT in the z direction
      ! The routine does not work with fields that include ghost layers
      !-------------------------------------------------------------------------
      SUBROUTINE ppm_fft_exec_3d_vec_c2c_z_s(topoid,meshid,ppmplan,infield,outfield,info)
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
      INTEGER,INTENT(IN) :: topoid
      !!!id of the mesh
      INTEGER,INTENT(IN) :: meshid
      !!!ppm fft plan type
      TYPE(ppm_fft_plan),INTENT(INOUT) :: ppmplan
      !!!input field to fourier transform
      !COMPLEX(ppm_kind_single),DIMENSION(:,:,:,:,:),POINTER,INTENT(INOUT) :: infield
      COMPLEX(ppm_kind_single),DIMENSION(:,:,:,:,:),POINTER :: infield
      !!!output field for the result of the fourier transform
      !COMPLEX(ppm_kind_single),DIMENSION(:,:,:,:,:),POINTER,INTENT(INOUT) :: outfield
      COMPLEX(ppm_kind_single),DIMENSION(:,:,:,:,:),POINTER :: outfield
      !!!Returns status, 0 upon success
      INTEGER,INTENT(OUT) :: info
      !in time perhaps an argument for alternate directions
      !-------------------------------------------------------------------------
      ! Local variables
      !-------------------------------------------------------------------------
      REAL(ppm_kind_single) :: t0
      INTEGER :: i,j
      INTEGER :: isub,isubl
      INTEGER :: nsubs
      INTEGER,DIMENSION(:),POINTER :: isublist
      TYPE(ppm_t_topo),POINTER :: topology
      TYPE(ppm_t_equi_mesh) :: mesh
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
      mesh = topology%mesh(meshid)
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
      END SUBROUTINE ppm_fft_exec_3d_vec_c2c_z_s
      !-------------------------------------------------------------------------
      ! Subroutine : ppm_fft_exec_3d_vec_fr2c_xy
      !-------------------------------------------------------------------------
      ! Copyright (c) 2010 CSE Lab (ETH Zurich), MOSAIC Group (ETH Zurich),
      ! Center for Fluid Dynamics (DTU)
      !
      ! FFTW execute wrapper for 3d arrays, 2d real to complex
      ! (forward) FFT in the xy directions
      ! The routine does not work with fields that include ghost layers
      !-------------------------------------------------------------------------
      SUBROUTINE ppm_fft_exec_3d_vec_fr2c_xy_s(topoid,meshid,ppmplan,infield,outfield,info)
      !!! FFTW execute wrapper for 3d arrays, 2d real to complex
      !!! (forward) FFT in the xy directions
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
      INTEGER,INTENT(IN) :: topoid
      !!!id of the mesh
      INTEGER,INTENT(IN) :: meshid
      !!!ppm fft plan type
      TYPE(ppm_fft_plan),INTENT(INOUT) :: ppmplan
      !!!input field to fourier transform
      !REAL(ppm_kind_single),DIMENSION(:,:,:,:,:),POINTER,INTENT(INOUT) :: infield
      REAL(ppm_kind_single),DIMENSION(:,:,:,:,:),POINTER :: infield
      !!!output field for the result of the fourier transform
      !COMPLEX(ppm_kind_single),DIMENSION(:,:,:,:,:),POINTER,INTENT(INOUT) :: outfield
      COMPLEX(ppm_kind_single),DIMENSION(:,:,:,:,:),POINTER :: outfield
      !!!Returns status, 0 upon success
      INTEGER,INTENT(OUT) :: info
      !in time perhaps an argument for alternate directions
      !-------------------------------------------------------------------------
      ! Local variables
      !-------------------------------------------------------------------------
      REAL(ppm_kind_single) :: t0
      INTEGER :: k
      INTEGER :: isub,isubl
      INTEGER :: nsubs
      INTEGER,DIMENSION(:),POINTER :: isublist
      TYPE(ppm_t_topo),POINTER :: topology
      TYPE(ppm_t_equi_mesh) :: mesh
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
      mesh = topology%mesh(meshid)
      !-------------------------------------------------------------------------
      ! Execute plan
      !-------------------------------------------------------------------------
      DO isub=1,nsubs
         isubl=isublist(isub)
         DO k=1,mesh%nnodes(3,isubl) !@ add '-1' to exclude n+1 slabs
            CALL dfftw_execute_dft_r2c(ppmplan%plan(isub),&
            & infield(1,1,1,k,isub),outfield(1,1,1,k,isub))
         END DO
      END DO
      !-------------------------------------------------------------------------
      ! Return
      !-------------------------------------------------------------------------
 9999 CONTINUE
      CALL substop('ppm_fft_exec',t0,info)
      RETURN
      END SUBROUTINE ppm_fft_exec_3d_vec_fr2c_xy_s
      !-------------------------------------------------------------------------
      ! Subroutine : ppm_fft_exec_3d_vec_bc2r_xy
      !-------------------------------------------------------------------------
      ! Copyright (c) 2010 CSE Lab (ETH Zurich), MOSAIC Group (ETH Zurich),
      ! Center for Fluid Dynamics (DTU)
      !
      ! FFTW execute wrapper for 3d arrays, 2d complex to real
      ! (backward) FFT in the xy directions
      ! The routine does not work with fields that include ghost layers
      !-------------------------------------------------------------------------
      SUBROUTINE ppm_fft_exec_3d_vec_bc2r_xy_s(topoid,meshid,ppmplan,infield,outfield,info)
      !!! FFTW execute wrapper for 3d arrays, 2d complex to real
      !!! (backward) FFT in the xy directions
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
      INTEGER,INTENT(IN) :: topoid
      !!!id of the mesh
      INTEGER,INTENT(IN) :: meshid
      !!!ppm fft plan type
      TYPE(ppm_fft_plan),INTENT(INOUT) :: ppmplan
      !!!input field to fourier transform
      !COMPLEX(ppm_kind_single),DIMENSION(:,:,:,:,:),POINTER,INTENT(INOUT) :: infield
      COMPLEX(ppm_kind_single),DIMENSION(:,:,:,:,:),POINTER :: infield
      !!!output field for the result of the fourier transform
      !REAL(ppm_kind_single),DIMENSION(:,:,:,:,:),POINTER,INTENT(INOUT) :: outfield
      REAL(ppm_kind_single),DIMENSION(:,:,:,:,:),POINTER :: outfield
      !!!Returns status, 0 upon success
      INTEGER,INTENT(OUT) :: info
      !in time perhaps an argument for alternate directions
      !-------------------------------------------------------------------------
      ! Local variables
      !-------------------------------------------------------------------------
      REAL(ppm_kind_single) :: t0
      INTEGER :: i,k
      INTEGER :: isub,isubl
      INTEGER :: nsubs
      INTEGER,DIMENSION(:),POINTER :: isublist
      TYPE(ppm_t_topo),POINTER :: topology
      TYPE(ppm_t_equi_mesh) :: mesh
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
      mesh = topology%mesh(meshid)
      !-------------------------------------------------------------------------
      ! Execute plan
      !-------------------------------------------------------------------------
      DO isub=1,nsubs
         isubl=isublist(isub)
         DO k=1,mesh%nnodes(3,isubl) !@ add '-1' to exclude n+1 slabs
            CALL dfftw_execute_dft_c2r(ppmplan%plan(isub),&
            & infield(1,1,1,k,isub),outfield(1,1,1,k,isub))
            IF (topology%bcdef(1) .EQ. ppm_param_bcdef_periodic) THEN !vertex
               !-------------------------------------------------------------------
               ! Copy periodic layer back - only for periodic 'N' vertex points
               !-------------------------------------------------------------------
               DO i=1,mesh%nnodes(1,isubl)
                  outfield(1,i,mesh%nnodes(2,isubl),k,isub) = outfield(1,i,1,k,isub)
                  outfield(2,i,mesh%nnodes(2,isubl),k,isub) = outfield(2,i,1,k,isub)
                  outfield(3,i,mesh%nnodes(2,isubl),k,isub) = outfield(3,i,1,k,isub)
               END DO
               DO i=1,mesh%nnodes(2,isubl)
                  outfield(1,mesh%nnodes(1,isubl),i,k,isub) = outfield(1,1,i,k,isub)
                  outfield(2,mesh%nnodes(1,isubl),i,k,isub) = outfield(2,1,i,k,isub)
                  outfield(3,mesh%nnodes(1,isubl),i,k,isub) = outfield(3,1,i,k,isub)
               END DO
            END IF
         END DO
      END DO
      !-------------------------------------------------------------------------
      ! Return
      !-------------------------------------------------------------------------
 9999 CONTINUE
      CALL substop('ppm_fft_exec',t0,info)
      RETURN
      END SUBROUTINE ppm_fft_exec_3d_vec_bc2r_xy_s
!#include "ppm_fft_exec_3d_vec_bc2c_xy.f"
      !NORMALIZATION
      !-------------------------------------------------------------------------
      ! Subroutine : ppm_fft_normalize_r
      !-------------------------------------------------------------------------
      ! Copyright (c) 2010 CSE Lab (ETH Zurich), MOSAIC Group (ETH Zurich),
      ! Center for Fluid Dynamics (DTU)
      !
      ! The routine does not work with fields that include ghost layers.
      ! Fields are noralized in the 1/N - fashion that scales fields to the
      ! proper level after being FFTed and iFFTed
      ! The routines are not meant for production code as the normalization can
      ! be done for all dimensions in one pass - which may and should be done
      ! in a routine (looping through the data anyway) following the FFTs.
      !-------------------------------------------------------------------------
      SUBROUTINE ppm_fft_normalize_rs(topoid,meshid,ppmplan,infield,info)
      !!! The routine does not work with fields that include ghost layers.
      !!! Fields are noralized in the 1/N - fashion that scales fields to the
      !!! proper level after being FFTed and iFFTed
      !!! The routines are not meant for production code as the normalization
      !!! can be done for all dimensions in one pass - which may and should be
      !!! done in a routine (looping through the data anyway) following the FFTs
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
      INTEGER,INTENT(IN) :: topoid
      !!!id of the mesh
      INTEGER,INTENT(IN) :: meshid
      !!!ppm fft plan type
      TYPE(ppm_fft_plan),INTENT(INOUT) :: ppmplan
      !!!input field to normalize
      !REAL(ppm_kind_single),DIMENSION(:,:,:,:,:),POINTER,INTENT(INOUT) :: infield
      REAL(ppm_kind_single),DIMENSION(:,:,:,:,:),POINTER :: infield
      !!!Returns status, 0 upon success
      INTEGER,INTENT(OUT) :: info
      !in time perhaps an argument for alternate directions
      !-------------------------------------------------------------------------
      ! Local variables
      !-------------------------------------------------------------------------
      !timer
      INTEGER,PARAMETER :: MK = ppm_kind_double
      REAL(ppm_kind_single) :: t0
      !normalization factor
      REAL(MK) :: fac
      INTEGER :: i,j,k
      INTEGER :: isub,isubl
      INTEGER :: nsubs
      INTEGER,DIMENSION(:),POINTER :: isublist
      TYPE(ppm_t_topo),POINTER :: topology
      TYPE(ppm_t_equi_mesh) :: mesh
      !-------------------------------------------------------------------------
      ! Initialise routine
      !-------------------------------------------------------------------------
      CALL substart('ppm_fft_normalize',t0,info)
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
      mesh = topology%mesh(meshid)
      DO isub=1,nsubs
         isubl=isublist(isub)
         !determine normalization factor
         IF (ppmplan%rank .EQ. 1) THEN
            IF (topology%bcdef(3) .EQ. ppm_param_bcdef_periodic) THEN !vertex
               fac = 1.0_MK/REAL((mesh%nm(3)-1),MK)
            ELSE
               fac = 1.0_MK/REAL((mesh%nm(3)),MK)
            ENDIF
         ELSE IF (ppmplan%rank .EQ. 2) THEN
            IF (topology%bcdef(1) .EQ. ppm_param_bcdef_periodic) THEN !vertex
               fac = 1.0_MK/REAL((mesh%nm(1)-1)*(mesh%nm(2)-1),MK)
            ELSE
               fac = 1.0_MK/REAL((mesh%nm(1) )*(mesh%nm(2) ),MK)
            ENDIF
         ENDIF
         DO k=1,mesh%nnodes(3,isubl)
            DO j=1,mesh%nnodes(2,isubl)
               DO i=1,mesh%nnodes(1,isubl)
                  infield(1,i,j,k,isub) = fac*infield(1,i,j,k,isub)
                  infield(2,i,j,k,isub) = fac*infield(2,i,j,k,isub)
                  infield(3,i,j,k,isub) = fac*infield(3,i,j,k,isub)
               END DO
            END DO
         END DO
      END DO
      !-------------------------------------------------------------------------
      ! Return
      !-------------------------------------------------------------------------
 9999 CONTINUE
      CALL substop('ppm_fft_normalize',t0,info)
      RETURN
      END SUBROUTINE ppm_fft_normalize_rs
      !-------------------------------------------------------------------------
      ! Subroutine : ppm_fft_normalize_c
      !-------------------------------------------------------------------------
      ! Copyright (c) 2010 CSE Lab (ETH Zurich), MOSAIC Group (ETH Zurich),
      ! Center for Fluid Dynamics (DTU)
      !
      ! The routine does not work with fields that include ghost layers.
      ! Fields are noralized in the 1/N - fashion that scales fields to the
      ! proper level after being FFTed and iFFTed
      ! The routines are not meant for production code as the normalization can
      ! be done for all dimensions in one pass - which may and should be done
      ! in a routine (looping through the data anyway) following the FFTs.
      !-------------------------------------------------------------------------
      SUBROUTINE ppm_fft_normalize_cs(topoid,meshid,ppmplan,infield,info)
      !!! The routine does not work with fields that include ghost layers.
      !!! Fields are noralized in the 1/N - fashion that scales fields to the
      !!! proper level after being FFTed and iFFTed
      !!! The routines are not meant for production code as the normalization
      !!! can be done for all dimensions in one pass - which may and should be
      !!! done in a routine (looping through the data anyway) following the FFTs
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
      INTEGER,INTENT(IN) :: topoid
      !!!id of the mesh
      INTEGER,INTENT(IN) :: meshid
      !!!ppm fft plan type
      TYPE(ppm_fft_plan),INTENT(INOUT) :: ppmplan
      !!!input field to normalize
      COMPLEX(ppm_kind_single),DIMENSION(:,:,:,:,:),POINTER :: infield
      !!!Returns status, 0 upon success
      INTEGER,INTENT(OUT) :: info
      !in time perhaps an argument for alternate directions
      !-------------------------------------------------------------------------
      ! Local variables
      !-------------------------------------------------------------------------
      !timer
      INTEGER,PARAMETER :: MK = ppm_kind_double
      REAL(ppm_kind_single) :: t0
      !normalization factor
      REAL(MK) :: fac
      INTEGER :: i,j,k
      INTEGER :: isub,isubl
      INTEGER :: nsubs
      INTEGER,DIMENSION(:),POINTER :: isublist
      TYPE(ppm_t_topo),POINTER :: topology
      TYPE(ppm_t_equi_mesh) :: mesh
      !-------------------------------------------------------------------------
      ! Initialise routine
      !-------------------------------------------------------------------------
      CALL substart('ppm_fft_normalize',t0,info)
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
      mesh = topology%mesh(meshid)
      DO isub=1,nsubs
         isubl=isublist(isub)
         !determine normalization factor
         IF (ppmplan%rank .EQ. 1) THEN
            IF (topology%bcdef(3) .EQ. ppm_param_bcdef_periodic) THEN !vertex
               fac = 1.0_MK/REAL((mesh%nm(3)-1),MK)
            ELSE
               fac = 1.0_MK/REAL((mesh%nm(3)),MK)
            ENDIF
         ELSE IF (ppmplan%rank .EQ. 2) THEN
            IF (topology%bcdef(1) .EQ. ppm_param_bcdef_periodic) THEN !vertex
               fac = 1.0_MK/REAL((mesh%nm(1)-1)*(mesh%nm(2)-1),MK)
            ELSE
               fac = 1.0_MK/REAL((mesh%nm(1) )*(mesh%nm(2) ),MK)
            ENDIF
         ENDIF
         DO k=1,mesh%nnodes(3,isubl)
            DO j=1,mesh%nnodes(2,isubl)
               DO i=1,mesh%nnodes(1,isubl)
                  infield(1,i,j,k,isub) = fac*infield(1,i,j,k,isub)
                  infield(2,i,j,k,isub) = fac*infield(2,i,j,k,isub)
                  infield(3,i,j,k,isub) = fac*infield(3,i,j,k,isub)
               END DO
            END DO
         END DO
      END DO
      !-------------------------------------------------------------------------
      ! Return
      !-------------------------------------------------------------------------
 9999 CONTINUE
      CALL substop('ppm_fft_normalize',t0,info)
      RETURN
      END SUBROUTINE ppm_fft_normalize_cs
      !FORWARD TRANSFORMS - PLAN
      !-------------------------------------------------------------------------
      ! Subroutine : ppm_fft_plan_3d_vec_fc2c_z
      !-------------------------------------------------------------------------
      ! Copyright (c) 2010 CSE Lab (ETH Zurich), MOSAIC Group (ETH Zurich),
      ! Center for Fluid Dynamics (DTU)
      !
      ! FFTW plan wrapper for 3d arrays, 1d complex to complex
      ! (forward) FFT in the xy directions
      ! The routine does not work with fields that include ghost layers
      !-------------------------------------------------------------------------
      SUBROUTINE ppm_fft_plan_3d_vec_fc2c_z_d(topoid,meshid,ppmplan,infield,outfield,info)
      !!! FFTW plan wrapper for 3d arrays, 1d complex to complex
      !!! (forward) FFT in the xy directions
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
      INTEGER,INTENT(IN) :: topoid
      !!!id of the mesh
      INTEGER,INTENT(IN) :: meshid
      !!!ppm fft plan type
      TYPE(ppm_fft_plan),INTENT(INOUT) :: ppmplan
      !!!input field to fourier transform
      !COMPLEX(ppm_kind_double),DIMENSION(:,:,:,:,:),POINTER,INTENT(INOUT) :: infield
      COMPLEX(ppm_kind_double),DIMENSION(:,:,:,:,:),POINTER :: infield
      !!!output field for the result of the fourier transform
      !COMPLEX(ppm_kind_double),DIMENSION(:,:,:,:,:),POINTER,INTENT(INOUT) :: outfield
      COMPLEX(ppm_kind_double),DIMENSION(:,:,:,:,:),POINTER :: outfield
      !!!Returns status, 0 upon success
      INTEGER,INTENT(OUT) :: info
      !in time perhaps an argument for alternate directions
      !-------------------------------------------------------------------------
      ! Local variables
      !-------------------------------------------------------------------------
      REAL(ppm_kind_double) :: t0
      INTEGER :: isub,isubl
      INTEGER :: nsubs
      INTEGER,DIMENSION(:),POINTER :: isublist
      TYPE(ppm_t_topo),POINTER :: topology
      TYPE(ppm_t_equi_mesh) :: mesh
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
      mesh = topology%mesh(meshid)
      !-------------------------------------------------------------------------
      ! Setup parameters for this particular routine
      !-------------------------------------------------------------------------
      !the dimension of the FFT (1D/2D/3D)
      ppmplan%rank=1
      !the number of points along each direction of the piece to be transformed
      ALLOCATE(ppmplan%nx(ppmplan%rank,nsubs))
      !the direction of the transform
      ppmplan%sign=FFTW_FORWARD
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
      !e.g. for 2/3 component vector istride = 2/3 or for scalar istride = 1
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
      END SUBROUTINE ppm_fft_plan_3d_vec_fc2c_z_d
      !-------------------------------------------------------------------------
      ! Subroutine : ppm_fft_plan_3d_vec_fr2c_xy
      !-------------------------------------------------------------------------
      ! Copyright (c) 2010 CSE Lab (ETH Zurich), MOSAIC Group (ETH Zurich),
      ! Center for Fluid Dynamics (DTU)
      !
      ! FFTW plan wrapper for 3d arrays, 2d real to complex
      ! (forward) FFT in the xy directions
      ! The routine does not work with fields that include ghost layers
      !-------------------------------------------------------------------------
      SUBROUTINE ppm_fft_plan_3d_vec_fr2c_xy_d(topoid,meshid,ppmplan,infield,outfield,info)
      !!! FFTW plan wrapper for 3d arrays, 2d real to complex
      !!! (forward) FFT in the xy directions
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
      INTEGER,INTENT(IN) :: topoid
      !!!id of the mesh
      INTEGER,INTENT(IN) :: meshid
      !!!ppm fft plan type
      TYPE(ppm_fft_plan),INTENT(INOUT) :: ppmplan
      !!!input field to fourier transform
      !REAL(ppm_kind_double),DIMENSION(:,:,:,:,:),POINTER,INTENT(INOUT) :: infield
      REAL(ppm_kind_double),DIMENSION(:,:,:,:,:),POINTER :: infield
      !!!output field for the result of the fourier transform
      !COMPLEX(ppm_kind_double),DIMENSION(:,:,:,:,:),POINTER,INTENT(INOUT) :: outfield
      COMPLEX(ppm_kind_double),DIMENSION(:,:,:,:,:),POINTER :: outfield
      !!!Returns status, 0 upon success
      INTEGER,INTENT(OUT) :: info
      !in time perhaps an argument for alternate directions
      !-------------------------------------------------------------------------
      ! Local variables
      !-------------------------------------------------------------------------
      REAL(ppm_kind_double) :: t0
      INTEGER :: isub,isubl
      INTEGER :: nsubs
      INTEGER,DIMENSION(:),POINTER :: isublist
      TYPE(ppm_t_topo),POINTER :: topology
      TYPE(ppm_t_equi_mesh) :: mesh
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
      mesh = topology%mesh(meshid)
      !-------------------------------------------------------------------------
      ! Setup parameters for this particular routine
      !-------------------------------------------------------------------------
      !the dimension of the FFT (1D/2D/3D)
      ppmplan%rank=2
      !the number of points along each direction of the piece to be transformed
      ALLOCATE(ppmplan%nx(ppmplan%rank,nsubs))
      !the direction of the transform
      ppmplan%sign=FFTW_FORWARD
      !the method to setup the optimal plan
      ppmplan%flag=FFTW_MEASURE
      !the number of components to transform - 3 component vector
      ppmplan%howmany=3
      !the size of the input array - full size (assuming LBOUND=1 thus UBOUND)
      ALLOCATE(ppmplan%inembed(ppmplan%rank))
      ppmplan%inembed(1) = UBOUND(infield,2)
      ppmplan%inembed(2) = UBOUND(infield,3)
      !the size of the output array - full size (assuming LBOUND=1 thus UBOUND)
      ALLOCATE(ppmplan%onembed(ppmplan%rank))
      ppmplan%onembed(1) = UBOUND(outfield,2)
      ppmplan%onembed(2) = UBOUND(outfield,3)
      !istride tells how the same componenet data points are spaced in memory
      !e.g. for 2/3 component vector istride = 2/3 or for scalar istride = 1
      ppmplan%istride = 3
      ppmplan%ostride = 3
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
         IF (topology%bcdef(1) .EQ. ppm_param_bcdef_periodic) THEN !vertex
            ppmplan%nx(1,isub) = mesh%nm(1)-1
            ppmplan%nx(2,isub) = mesh%nm(2)-1
         ELSE
            ppmplan%nx(1,isub) = mesh%nm(1)
            ppmplan%nx(2,isub) = mesh%nm(2)
         ENDIF
         CALL dfftw_plan_many_dft_r2c(ppmplan%plan(isub),ppmplan%rank,&
         & ppmplan%nx(:,isub),ppmplan%howmany,infield(1,1,1,1,isub),&
         & ppmplan%inembed(1),ppmplan%istride,ppmplan%idist,&
         & outfield(1,1,1,1,isub),ppmplan%onembed(1),ppmplan%ostride,&
         & ppmplan%odist,ppmplan%flag)
      END DO
      !-------------------------------------------------------------------------
      ! Return
      !-------------------------------------------------------------------------
 9999 CONTINUE
      CALL substop('ppm_fft_plan',t0,info)
      RETURN
      END SUBROUTINE ppm_fft_plan_3d_vec_fr2c_xy_d
      !BACKWARD TRANSFORMS - PLAN
      !-------------------------------------------------------------------------
      ! Subroutine : ppm_fft_plan_3d_vec_bc2c_z
      !-------------------------------------------------------------------------
      ! Copyright (c) 2010 CSE Lab (ETH Zurich), MOSAIC Group (ETH Zurich),
      ! Center for Fluid Dynamics (DTU)
      !
      ! FFTW plan wrapper for 3d arrays, 1d complex to complex
      ! (backward) FFT in the z direction
      ! The routine does not work with fields that include ghost layers
      !-------------------------------------------------------------------------
      SUBROUTINE ppm_fft_plan_3d_vec_bc2c_z_d(topoid,meshid,ppmplan,infield,outfield,info)
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
      INTEGER,INTENT(IN) :: topoid
      !!!id of the mesh
      INTEGER,INTENT(IN) :: meshid
      !!!ppm fft plan type
      TYPE(ppm_fft_plan),INTENT(INOUT) :: ppmplan
      !!!input field to fourier transform
      !COMPLEX(ppm_kind_double),DIMENSION(:,:,:,:,:),POINTER,INTENT(INOUT) :: infield
      COMPLEX(ppm_kind_double),DIMENSION(:,:,:,:,:),POINTER :: infield
      !!!output field for the result of the fourier transform
      !COMPLEX(ppm_kind_double),DIMENSION(:,:,:,:,:),POINTER,INTENT(INOUT) :: outfield
      COMPLEX(ppm_kind_double),DIMENSION(:,:,:,:,:),POINTER :: outfield
      !!!Returns status, 0 upon success
      INTEGER,INTENT(OUT) :: info
      !in time perhaps an argument for alternate directions
      !-------------------------------------------------------------------------
      ! Local variables
      !-------------------------------------------------------------------------
      REAL(ppm_kind_double) :: t0
      INTEGER :: isub,isubl
      INTEGER :: nsubs
      INTEGER,DIMENSION(:),POINTER :: isublist
      TYPE(ppm_t_topo),POINTER :: topology
      TYPE(ppm_t_equi_mesh) :: mesh
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
      mesh = topology%mesh(meshid)
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
      END SUBROUTINE ppm_fft_plan_3d_vec_bc2c_z_d
      !-------------------------------------------------------------------------
      ! Subroutine : ppm_fft_plan_3d_vec_bc2r_xy
      !-------------------------------------------------------------------------
      ! Copyright (c) 2010 CSE Lab (ETH Zurich), MOSAIC Group (ETH Zurich),
      ! Center for Fluid Dynamics (DTU)
      !
      ! FFTW plan wrapper for 3d arrays, 2d complex to real
      ! (backward) FFT in the xy directions
      ! The routine does not work with fields that include ghost layers
      !
      ! For the C2R inverse transform FFTW always destroys the input data
      !-------------------------------------------------------------------------
      SUBROUTINE ppm_fft_plan_3d_vec_bc2r_xy_d(topoid,meshid,ppmplan,infield,outfield,info)
      !!! FFTW plan wrapper for 3d arrays, 2d complex to real
      !!! (backward) FFT in the xy directions
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
      INTEGER,INTENT(IN) :: topoid
      !!!id of the mesh
      INTEGER,INTENT(IN) :: meshid
      !!!ppm fft plan type
      TYPE(ppm_fft_plan),INTENT(INOUT) :: ppmplan
      !!!input field to fourier transform
      !COMPLEX(ppm_kind_double),DIMENSION(:,:,:,:,:),POINTER,INTENT(INOUT) :: infield
      COMPLEX(ppm_kind_double),DIMENSION(:,:,:,:,:),POINTER :: infield
      !!!output field for the result of the fourier transform
      !REAL(ppm_kind_double),DIMENSION(:,:,:,:,:),POINTER,INTENT(INOUT) :: outfield
      REAL(ppm_kind_double),DIMENSION(:,:,:,:,:),POINTER :: outfield
      !!!Returns status, 0 upon success
      INTEGER,INTENT(OUT) :: info
      !in time perhaps an argument for alternate directions
      !-------------------------------------------------------------------------
      ! Local variables
      !-------------------------------------------------------------------------
      REAL(ppm_kind_double) :: t0
      INTEGER :: isub,isubl
      INTEGER :: nsubs
      INTEGER,DIMENSION(:),POINTER :: isublist
      TYPE(ppm_t_topo),POINTER :: topology
      TYPE(ppm_t_equi_mesh) :: mesh
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
      mesh = topology%mesh(meshid)
      !-------------------------------------------------------------------------
      ! Setup parameters for this particular routine
      !-------------------------------------------------------------------------
      !the dimension of the FFT (1D/2D/3D)
      ppmplan%rank=2
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
      ppmplan%inembed(1) = UBOUND(infield,2)
      ppmplan%inembed(2) = UBOUND(infield,3)
      !the size of the output array - full size (assuming LBOUND=1 thus UBOUND)
      ALLOCATE(ppmplan%onembed(ppmplan%rank))
      ppmplan%onembed(1) = UBOUND(outfield,2)
      ppmplan%onembed(2) = UBOUND(outfield,3)
      !istride tells how the same componenet data points are spaced in memory
      !e.g. for 2/3 component vector istride = 2/3 or for scalar istride = 1
      ppmplan%istride = 3
      ppmplan%ostride = 3
      !idist tells how multiple arrays are spaced in memory. I.e. a memory
      !offset. e.g. vector components (idist=1) or scalar 2D arrays
      !in 3D array(idist=NxNy)
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
         !the number of points in each direction of the piece to be transformed
         !we subtract the -1 to avoid the periodic vertex point
         IF (topology%bcdef(1) .EQ. ppm_param_bcdef_periodic) THEN !vertex
            ppmplan%nx(1,isub) = mesh%nm(1)-1
            ppmplan%nx(2,isub) = mesh%nm(2)-1
         ELSE
            ppmplan%nx(1,isub) = mesh%nm(1)
            ppmplan%nx(2,isub) = mesh%nm(2)
         ENDIF
         CALL dfftw_plan_many_dft_c2r(ppmplan%plan(isub),ppmplan%rank,&
         & ppmplan%nx,ppmplan%howmany,infield(1,1,1,1,isub),ppmplan%inembed(1),&
         & ppmplan%istride,ppmplan%idist,outfield(1,1,1,1,isub),&
         & ppmplan%onembed(1),ppmplan%ostride,ppmplan%odist,ppmplan%flag)
      END DO
      !-------------------------------------------------------------------------
      ! Return
      !-------------------------------------------------------------------------
 9999 CONTINUE
      CALL substop('ppm_fft_plan',t0,info)
      RETURN
      END SUBROUTINE ppm_fft_plan_3d_vec_bc2r_xy_d
!#include "ppm_fft_plan_3d_vec_bc2c_xy.f"
      !EXECUTION OF PLANS
      !-------------------------------------------------------------------------
      ! Subroutine : ppm_fft_exec_3d_vec_c2c_z
      !-------------------------------------------------------------------------
      ! Copyright (c) 2010 CSE Lab (ETH Zurich), MOSAIC Group (ETH Zurich),
      ! Center for Fluid Dynamics (DTU)
      !
      ! FFTW execute wrapper for 3d arrays, 1d complex to complex
      ! (forward and backward) FFT in the z direction
      ! The routine does not work with fields that include ghost layers
      !-------------------------------------------------------------------------
      SUBROUTINE ppm_fft_exec_3d_vec_c2c_z_d(topoid,meshid,ppmplan,infield,outfield,info)
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
      INTEGER,INTENT(IN) :: topoid
      !!!id of the mesh
      INTEGER,INTENT(IN) :: meshid
      !!!ppm fft plan type
      TYPE(ppm_fft_plan),INTENT(INOUT) :: ppmplan
      !!!input field to fourier transform
      !COMPLEX(ppm_kind_double),DIMENSION(:,:,:,:,:),POINTER,INTENT(INOUT) :: infield
      COMPLEX(ppm_kind_double),DIMENSION(:,:,:,:,:),POINTER :: infield
      !!!output field for the result of the fourier transform
      !COMPLEX(ppm_kind_double),DIMENSION(:,:,:,:,:),POINTER,INTENT(INOUT) :: outfield
      COMPLEX(ppm_kind_double),DIMENSION(:,:,:,:,:),POINTER :: outfield
      !!!Returns status, 0 upon success
      INTEGER,INTENT(OUT) :: info
      !in time perhaps an argument for alternate directions
      !-------------------------------------------------------------------------
      ! Local variables
      !-------------------------------------------------------------------------
      REAL(ppm_kind_double) :: t0
      INTEGER :: i,j
      INTEGER :: isub,isubl
      INTEGER :: nsubs
      INTEGER,DIMENSION(:),POINTER :: isublist
      TYPE(ppm_t_topo),POINTER :: topology
      TYPE(ppm_t_equi_mesh) :: mesh
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
      mesh = topology%mesh(meshid)
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
      END SUBROUTINE ppm_fft_exec_3d_vec_c2c_z_d
      !-------------------------------------------------------------------------
      ! Subroutine : ppm_fft_exec_3d_vec_fr2c_xy
      !-------------------------------------------------------------------------
      ! Copyright (c) 2010 CSE Lab (ETH Zurich), MOSAIC Group (ETH Zurich),
      ! Center for Fluid Dynamics (DTU)
      !
      ! FFTW execute wrapper for 3d arrays, 2d real to complex
      ! (forward) FFT in the xy directions
      ! The routine does not work with fields that include ghost layers
      !-------------------------------------------------------------------------
      SUBROUTINE ppm_fft_exec_3d_vec_fr2c_xy_d(topoid,meshid,ppmplan,infield,outfield,info)
      !!! FFTW execute wrapper for 3d arrays, 2d real to complex
      !!! (forward) FFT in the xy directions
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
      INTEGER,INTENT(IN) :: topoid
      !!!id of the mesh
      INTEGER,INTENT(IN) :: meshid
      !!!ppm fft plan type
      TYPE(ppm_fft_plan),INTENT(INOUT) :: ppmplan
      !!!input field to fourier transform
      !REAL(ppm_kind_double),DIMENSION(:,:,:,:,:),POINTER,INTENT(INOUT) :: infield
      REAL(ppm_kind_double),DIMENSION(:,:,:,:,:),POINTER :: infield
      !!!output field for the result of the fourier transform
      !COMPLEX(ppm_kind_double),DIMENSION(:,:,:,:,:),POINTER,INTENT(INOUT) :: outfield
      COMPLEX(ppm_kind_double),DIMENSION(:,:,:,:,:),POINTER :: outfield
      !!!Returns status, 0 upon success
      INTEGER,INTENT(OUT) :: info
      !in time perhaps an argument for alternate directions
      !-------------------------------------------------------------------------
      ! Local variables
      !-------------------------------------------------------------------------
      REAL(ppm_kind_double) :: t0
      INTEGER :: k
      INTEGER :: isub,isubl
      INTEGER :: nsubs
      INTEGER,DIMENSION(:),POINTER :: isublist
      TYPE(ppm_t_topo),POINTER :: topology
      TYPE(ppm_t_equi_mesh) :: mesh
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
      mesh = topology%mesh(meshid)
      !-------------------------------------------------------------------------
      ! Execute plan
      !-------------------------------------------------------------------------
      DO isub=1,nsubs
         isubl=isublist(isub)
         DO k=1,mesh%nnodes(3,isubl) !@ add '-1' to exclude n+1 slabs
            CALL dfftw_execute_dft_r2c(ppmplan%plan(isub),&
            & infield(1,1,1,k,isub),outfield(1,1,1,k,isub))
         END DO
      END DO
      !-------------------------------------------------------------------------
      ! Return
      !-------------------------------------------------------------------------
 9999 CONTINUE
      CALL substop('ppm_fft_exec',t0,info)
      RETURN
      END SUBROUTINE ppm_fft_exec_3d_vec_fr2c_xy_d
      !-------------------------------------------------------------------------
      ! Subroutine : ppm_fft_exec_3d_vec_bc2r_xy
      !-------------------------------------------------------------------------
      ! Copyright (c) 2010 CSE Lab (ETH Zurich), MOSAIC Group (ETH Zurich),
      ! Center for Fluid Dynamics (DTU)
      !
      ! FFTW execute wrapper for 3d arrays, 2d complex to real
      ! (backward) FFT in the xy directions
      ! The routine does not work with fields that include ghost layers
      !-------------------------------------------------------------------------
      SUBROUTINE ppm_fft_exec_3d_vec_bc2r_xy_d(topoid,meshid,ppmplan,infield,outfield,info)
      !!! FFTW execute wrapper for 3d arrays, 2d complex to real
      !!! (backward) FFT in the xy directions
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
      INTEGER,INTENT(IN) :: topoid
      !!!id of the mesh
      INTEGER,INTENT(IN) :: meshid
      !!!ppm fft plan type
      TYPE(ppm_fft_plan),INTENT(INOUT) :: ppmplan
      !!!input field to fourier transform
      !COMPLEX(ppm_kind_double),DIMENSION(:,:,:,:,:),POINTER,INTENT(INOUT) :: infield
      COMPLEX(ppm_kind_double),DIMENSION(:,:,:,:,:),POINTER :: infield
      !!!output field for the result of the fourier transform
      !REAL(ppm_kind_double),DIMENSION(:,:,:,:,:),POINTER,INTENT(INOUT) :: outfield
      REAL(ppm_kind_double),DIMENSION(:,:,:,:,:),POINTER :: outfield
      !!!Returns status, 0 upon success
      INTEGER,INTENT(OUT) :: info
      !in time perhaps an argument for alternate directions
      !-------------------------------------------------------------------------
      ! Local variables
      !-------------------------------------------------------------------------
      REAL(ppm_kind_double) :: t0
      INTEGER :: i,k
      INTEGER :: isub,isubl
      INTEGER :: nsubs
      INTEGER,DIMENSION(:),POINTER :: isublist
      TYPE(ppm_t_topo),POINTER :: topology
      TYPE(ppm_t_equi_mesh) :: mesh
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
      mesh = topology%mesh(meshid)
      !-------------------------------------------------------------------------
      ! Execute plan
      !-------------------------------------------------------------------------
      DO isub=1,nsubs
         isubl=isublist(isub)
         DO k=1,mesh%nnodes(3,isubl) !@ add '-1' to exclude n+1 slabs
            CALL dfftw_execute_dft_c2r(ppmplan%plan(isub),&
            & infield(1,1,1,k,isub),outfield(1,1,1,k,isub))
            IF (topology%bcdef(1) .EQ. ppm_param_bcdef_periodic) THEN !vertex
               !-------------------------------------------------------------------
               ! Copy periodic layer back - only for periodic 'N' vertex points
               !-------------------------------------------------------------------
               DO i=1,mesh%nnodes(1,isubl)
                  outfield(1,i,mesh%nnodes(2,isubl),k,isub) = outfield(1,i,1,k,isub)
                  outfield(2,i,mesh%nnodes(2,isubl),k,isub) = outfield(2,i,1,k,isub)
                  outfield(3,i,mesh%nnodes(2,isubl),k,isub) = outfield(3,i,1,k,isub)
               END DO
               DO i=1,mesh%nnodes(2,isubl)
                  outfield(1,mesh%nnodes(1,isubl),i,k,isub) = outfield(1,1,i,k,isub)
                  outfield(2,mesh%nnodes(1,isubl),i,k,isub) = outfield(2,1,i,k,isub)
                  outfield(3,mesh%nnodes(1,isubl),i,k,isub) = outfield(3,1,i,k,isub)
               END DO
            END IF
         END DO
      END DO
      !-------------------------------------------------------------------------
      ! Return
      !-------------------------------------------------------------------------
 9999 CONTINUE
      CALL substop('ppm_fft_exec',t0,info)
      RETURN
      END SUBROUTINE ppm_fft_exec_3d_vec_bc2r_xy_d
!#include "ppm_fft_exec_3d_vec_bc2c_xy.f"
      !NORMALIZATION
      !-------------------------------------------------------------------------
      ! Subroutine : ppm_fft_normalize_r
      !-------------------------------------------------------------------------
      ! Copyright (c) 2010 CSE Lab (ETH Zurich), MOSAIC Group (ETH Zurich),
      ! Center for Fluid Dynamics (DTU)
      !
      ! The routine does not work with fields that include ghost layers.
      ! Fields are noralized in the 1/N - fashion that scales fields to the
      ! proper level after being FFTed and iFFTed
      ! The routines are not meant for production code as the normalization can
      ! be done for all dimensions in one pass - which may and should be done
      ! in a routine (looping through the data anyway) following the FFTs.
      !-------------------------------------------------------------------------
      SUBROUTINE ppm_fft_normalize_rd(topoid,meshid,ppmplan,infield,info)
      !!! The routine does not work with fields that include ghost layers.
      !!! Fields are noralized in the 1/N - fashion that scales fields to the
      !!! proper level after being FFTed and iFFTed
      !!! The routines are not meant for production code as the normalization
      !!! can be done for all dimensions in one pass - which may and should be
      !!! done in a routine (looping through the data anyway) following the FFTs
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
      INTEGER,INTENT(IN) :: topoid
      !!!id of the mesh
      INTEGER,INTENT(IN) :: meshid
      !!!ppm fft plan type
      TYPE(ppm_fft_plan),INTENT(INOUT) :: ppmplan
      !!!input field to normalize
      !REAL(ppm_kind_double),DIMENSION(:,:,:,:,:),POINTER,INTENT(INOUT) :: infield
      REAL(ppm_kind_double),DIMENSION(:,:,:,:,:),POINTER :: infield
      !!!Returns status, 0 upon success
      INTEGER,INTENT(OUT) :: info
      !in time perhaps an argument for alternate directions
      !-------------------------------------------------------------------------
      ! Local variables
      !-------------------------------------------------------------------------
      !timer
      INTEGER,PARAMETER :: MK = ppm_kind_double
      REAL(ppm_kind_double) :: t0
      !normalization factor
      REAL(MK) :: fac
      INTEGER :: i,j,k
      INTEGER :: isub,isubl
      INTEGER :: nsubs
      INTEGER,DIMENSION(:),POINTER :: isublist
      TYPE(ppm_t_topo),POINTER :: topology
      TYPE(ppm_t_equi_mesh) :: mesh
      !-------------------------------------------------------------------------
      ! Initialise routine
      !-------------------------------------------------------------------------
      CALL substart('ppm_fft_normalize',t0,info)
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
      mesh = topology%mesh(meshid)
      DO isub=1,nsubs
         isubl=isublist(isub)
         !determine normalization factor
         IF (ppmplan%rank .EQ. 1) THEN
            IF (topology%bcdef(3) .EQ. ppm_param_bcdef_periodic) THEN !vertex
               fac = 1.0_MK/REAL((mesh%nm(3)-1),MK)
            ELSE
               fac = 1.0_MK/REAL((mesh%nm(3)),MK)
            ENDIF
         ELSE IF (ppmplan%rank .EQ. 2) THEN
            IF (topology%bcdef(1) .EQ. ppm_param_bcdef_periodic) THEN !vertex
               fac = 1.0_MK/REAL((mesh%nm(1)-1)*(mesh%nm(2)-1),MK)
            ELSE
               fac = 1.0_MK/REAL((mesh%nm(1) )*(mesh%nm(2) ),MK)
            ENDIF
         ENDIF
         DO k=1,mesh%nnodes(3,isubl)
            DO j=1,mesh%nnodes(2,isubl)
               DO i=1,mesh%nnodes(1,isubl)
                  infield(1,i,j,k,isub) = fac*infield(1,i,j,k,isub)
                  infield(2,i,j,k,isub) = fac*infield(2,i,j,k,isub)
                  infield(3,i,j,k,isub) = fac*infield(3,i,j,k,isub)
               END DO
            END DO
         END DO
      END DO
      !-------------------------------------------------------------------------
      ! Return
      !-------------------------------------------------------------------------
 9999 CONTINUE
      CALL substop('ppm_fft_normalize',t0,info)
      RETURN
      END SUBROUTINE ppm_fft_normalize_rd
      !-------------------------------------------------------------------------
      ! Subroutine : ppm_fft_normalize_c
      !-------------------------------------------------------------------------
      ! Copyright (c) 2010 CSE Lab (ETH Zurich), MOSAIC Group (ETH Zurich),
      ! Center for Fluid Dynamics (DTU)
      !
      ! The routine does not work with fields that include ghost layers.
      ! Fields are noralized in the 1/N - fashion that scales fields to the
      ! proper level after being FFTed and iFFTed
      ! The routines are not meant for production code as the normalization can
      ! be done for all dimensions in one pass - which may and should be done
      ! in a routine (looping through the data anyway) following the FFTs.
      !-------------------------------------------------------------------------
      SUBROUTINE ppm_fft_normalize_cd(topoid,meshid,ppmplan,infield,info)
      !!! The routine does not work with fields that include ghost layers.
      !!! Fields are noralized in the 1/N - fashion that scales fields to the
      !!! proper level after being FFTed and iFFTed
      !!! The routines are not meant for production code as the normalization
      !!! can be done for all dimensions in one pass - which may and should be
      !!! done in a routine (looping through the data anyway) following the FFTs
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
      INTEGER,INTENT(IN) :: topoid
      !!!id of the mesh
      INTEGER,INTENT(IN) :: meshid
      !!!ppm fft plan type
      TYPE(ppm_fft_plan),INTENT(INOUT) :: ppmplan
      !!!input field to normalize
      COMPLEX(ppm_kind_double),DIMENSION(:,:,:,:,:),POINTER :: infield
      !!!Returns status, 0 upon success
      INTEGER,INTENT(OUT) :: info
      !in time perhaps an argument for alternate directions
      !-------------------------------------------------------------------------
      ! Local variables
      !-------------------------------------------------------------------------
      !timer
      INTEGER,PARAMETER :: MK = ppm_kind_double
      REAL(ppm_kind_double) :: t0
      !normalization factor
      REAL(MK) :: fac
      INTEGER :: i,j,k
      INTEGER :: isub,isubl
      INTEGER :: nsubs
      INTEGER,DIMENSION(:),POINTER :: isublist
      TYPE(ppm_t_topo),POINTER :: topology
      TYPE(ppm_t_equi_mesh) :: mesh
      !-------------------------------------------------------------------------
      ! Initialise routine
      !-------------------------------------------------------------------------
      CALL substart('ppm_fft_normalize',t0,info)
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
      mesh = topology%mesh(meshid)
      DO isub=1,nsubs
         isubl=isublist(isub)
         !determine normalization factor
         IF (ppmplan%rank .EQ. 1) THEN
            IF (topology%bcdef(3) .EQ. ppm_param_bcdef_periodic) THEN !vertex
               fac = 1.0_MK/REAL((mesh%nm(3)-1),MK)
            ELSE
               fac = 1.0_MK/REAL((mesh%nm(3)),MK)
            ENDIF
         ELSE IF (ppmplan%rank .EQ. 2) THEN
            IF (topology%bcdef(1) .EQ. ppm_param_bcdef_periodic) THEN !vertex
               fac = 1.0_MK/REAL((mesh%nm(1)-1)*(mesh%nm(2)-1),MK)
            ELSE
               fac = 1.0_MK/REAL((mesh%nm(1) )*(mesh%nm(2) ),MK)
            ENDIF
         ENDIF
         DO k=1,mesh%nnodes(3,isubl)
            DO j=1,mesh%nnodes(2,isubl)
               DO i=1,mesh%nnodes(1,isubl)
                  infield(1,i,j,k,isub) = fac*infield(1,i,j,k,isub)
                  infield(2,i,j,k,isub) = fac*infield(2,i,j,k,isub)
                  infield(3,i,j,k,isub) = fac*infield(3,i,j,k,isub)
               END DO
            END DO
         END DO
      END DO
      !-------------------------------------------------------------------------
      ! Return
      !-------------------------------------------------------------------------
 9999 CONTINUE
      CALL substop('ppm_fft_normalize',t0,info)
      RETURN
      END SUBROUTINE ppm_fft_normalize_cd
      END MODULE ppm_module_fft
