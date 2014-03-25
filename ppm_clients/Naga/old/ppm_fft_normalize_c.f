      !-------------------------------------------------------------------------
      !  Subroutine   : ppm_fft_normalize_c
      !-------------------------------------------------------------------------
      ! Copyright (c) 2010 CSE Lab (ETH Zurich), MOSAIC Group (ETH Zurich),
      !                    Center for Fluid Dynamics (DTU)
      !
      ! The routine does not work with fields that include ghost layers.
      ! Fields are noralized in the 1/N - fashion that scales fields to the
      ! proper level after being FFTed and iFFTed
      ! The routines are not meant for production code as the normalization can
      ! be done for all dimensions in one pass - which may and should be done
      ! in a routine (looping through the data anyway) following the FFTs.
      !-------------------------------------------------------------------------
#if __KIND == __SINGLE
#define __ROUTINE ppm_fft_normalize_cs
#define __PREC ppm_kind_single
#elif __KIND == __DOUBLE
#define __ROUTINE ppm_fft_normalize_cd
#define __PREC ppm_kind_double
#endif
      SUBROUTINE __ROUTINE(topoid,meshid,ppmplan,infield,info)
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
      INTEGER,INTENT(IN)                                             :: topoid
      !!!id of the mesh
      INTEGER,INTENT(IN)                                             :: meshid
      !!!ppm fft plan type
      TYPE(ppm_fft_plan),INTENT(INOUT)                               :: ppmplan
      !!!input field to normalize
      COMPLEX(__PREC),DIMENSION(:,:,:,:,:),POINTER                   :: infield
      !!!Returns status, 0 upon success
      INTEGER,INTENT(OUT)                                            :: info
      !in time perhaps an argument for alternate directions

      !-------------------------------------------------------------------------
      ! Local variables
      !-------------------------------------------------------------------------
      !timer
      INTEGER,PARAMETER             :: MK = ppm_kind_double
      REAL(__PREC)                  :: t0
      !normalization factor
      REAL(MK)                      :: fac
      INTEGER                       :: i,j,k
      INTEGER                       :: isub,isubl
      INTEGER                       :: nsubs
      INTEGER,DIMENSION(:),POINTER  :: isublist
      TYPE(ppm_t_topo),POINTER      :: topology
      TYPE(ppm_t_equi_mesh)         :: mesh

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
      mesh  = topology%mesh(meshid)

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
               fac = 1.0_MK/REAL((mesh%nm(1)  )*(mesh%nm(2)  ),MK)
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

      END SUBROUTINE __ROUTINE

#undef __ROUTINE
#undef __PREC
