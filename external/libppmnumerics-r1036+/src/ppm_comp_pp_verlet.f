      !-------------------------------------------------------------------------
      !  Subroutine   :                 ppm_comp_pp_verlet
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

#if   __KERNEL == __INTERNAL
#if   __KIND == __SINGLE_PRECISION
      SUBROUTINE ppm_comp_pp_verlet_si(xp,Np,pdata,lda,lsymm,kernel,kpar,  &
     &    nvlist,vlist,dpd,info)
#elif __KIND == __DOUBLE_PRECISION
      SUBROUTINE ppm_comp_pp_verlet_di(xp,Np,pdata,lda,lsymm,kernel,kpar,  &
     &    nvlist,vlist,dpd,info)
#elif __KIND == __SINGLE_PRECISION_COMPLEX
      SUBROUTINE ppm_comp_pp_verlet_sci(xp,Np,pdata,lda,lsymm,kernel,kpar, &
     &    nvlist,vlist,dpd,info)
#elif __KIND == __DOUBLE_PRECISION_COMPLEX
      SUBROUTINE ppm_comp_pp_verlet_dci(xp,Np,pdata,lda,lsymm,kernel,kpar, &
     &    nvlist,vlist,dpd,info)
#endif
#elif __KERNEL == __USER_FUNCTION
#if   __KIND == __SINGLE_PRECISION
      SUBROUTINE ppm_comp_pp_verlet_su(xp,Np,pdata,lda,lsymm,kernel,kpar,  &
     &    nvlist,vlist,dpd,info)
#elif __KIND == __DOUBLE_PRECISION
      SUBROUTINE ppm_comp_pp_verlet_du(xp,Np,pdata,lda,lsymm,kernel,kpar,  &
     &    nvlist,vlist,dpd,info)
#elif __KIND == __SINGLE_PRECISION_COMPLEX
      SUBROUTINE ppm_comp_pp_verlet_scu(xp,Np,pdata,lda,lsymm,kernel,kpar, &
     &    nvlist,vlist,dpd,info)
#elif __KIND == __DOUBLE_PRECISION_COMPLEX
      SUBROUTINE ppm_comp_pp_verlet_dcu(xp,Np,pdata,lda,lsymm,kernel,kpar, &
     &    nvlist,vlist,dpd,info)
#endif
#elif __KERNEL == __LOOKUP_TABLE
#if   __KIND == __SINGLE_PRECISION
      SUBROUTINE ppm_comp_pp_verlet_st(xp,Np,pdata,lda,lsymm,kernel,kpar,  &
     &    nvlist,vlist,dpd,info)
#elif __KIND == __DOUBLE_PRECISION
      SUBROUTINE ppm_comp_pp_verlet_dt(xp,Np,pdata,lda,lsymm,kernel,kpar,  &
     &    nvlist,vlist,dpd,info)
#elif __KIND == __SINGLE_PRECISION_COMPLEX
      SUBROUTINE ppm_comp_pp_verlet_sct(xp,Np,pdata,lda,lsymm,kernel,kpar, &
     &    nvlist,vlist,dpd,info)
#elif __KIND == __DOUBLE_PRECISION_COMPLEX
      SUBROUTINE ppm_comp_pp_verlet_dct(xp,Np,pdata,lda,lsymm,kernel,kpar, &
     &    nvlist,vlist,dpd,info)
#endif
#endif
      !!! Subroutine which computes kernel interactions by direct
      !!! particle-particle interactions using Verlet lists.
      !!!
      !!! [NOTE]
      !!! ======================================================================
      !!! The loops for lda.EQ.1 and lda.EQ.2 are explicitly
      !!! unfolded to allow vectorization in those cases. For
      !!! lda.GT.2 this routine will not vectorize.
      !!!
      !!! Both nvlist and vlist can be created using
      !!! ppm_neighlist_vlist. This should be done by the
      !!! user program since there could be several sets of
      !!! Verlet lists. This routine can then be called with
      !!! the appropriate one. Do not forget to DEALLOCATE
      !!! nvlist and vlist when they are not needed any more.
      !!!
      !!! dpd needs to be allocated to proper size before
      !!! calling this routine. Also, this routine is not
      !!! resetting dpd to zero before doing the PP
      !!! interactions. This allows contributions from
      !!! different kernels to be accumulated. If needed,
      !!! set it to zero before calling this routine the
      !!! first time.
      !!! ======================================================================
      !!!
      !-------------------------------------------------------------------------
      !  Modules
      !-------------------------------------------------------------------------
      USE ppm_module_data
      USE ppm_module_substart
      USE ppm_module_substop
      USE ppm_module_error
      IMPLICIT NONE
#if   __KIND == __SINGLE_PRECISION | __KIND == __SINGLE_PRECISION_COMPLEX
      INTEGER, PARAMETER :: MK = ppm_kind_single
#elif __KIND == __DOUBLE_PRECISION | __KIND == __DOUBLE_PRECISION_COMPLEX
      INTEGER, PARAMETER :: MK = ppm_kind_double
#endif
      !-------------------------------------------------------------------------
      !  Includes
      !-------------------------------------------------------------------------
#include "ppm_define.h"
      !-------------------------------------------------------------------------
      !  Arguments     
      !-------------------------------------------------------------------------
      REAL(MK)   , DIMENSION(:,:), INTENT(IN   ) :: xp
      !!! particle coordinates
#if   __KIND == __SINGLE_PRECISION | __KIND == __DOUBLE_PRECISION
      REAL(MK)   , DIMENSION(:,:), INTENT(IN   ) :: pdata
      !!! particle data (strengths).
      !!! Overloaded types: single,double
      REAL(MK)   , DIMENSION(:,:), INTENT(  OUT) :: dpd
      !!! Change of particle data (strengths).
      !!! Overloaded types: single,double
#elif __KIND == __SINGLE_PRECISION_COMPLEX | __KIND == __DOUBLE_PRECISION_COMPLEX
      COMPLEX(MK), DIMENSION(:,:), INTENT(IN   ) :: pdata
      !!! particle data (strengths).
      !!! Overloaded types: single complex,double complex
      COMPLEX(MK), DIMENSION(:,:), INTENT(  OUT) :: dpd
      !!! Change of particle data (strengths).
      !!! Overloaded types: single complex,double complex.
#endif
      LOGICAL                    , INTENT(IN   ) :: lsymm
      !!! using symmetry or not?
      INTEGER                    , INTENT(IN   ) :: Np
      !!! number of particles on local proc.
      INTEGER                    , INTENT(IN   ) :: lda
      !!! leading dimension of pdata.
#if   __KERNEL == __INTERNAL
      INTEGER                    , INTENT(IN   ) :: kernel
      !!! kernel for which to compute the correction. To use ppm-internal
      !!! kernels, specify one of:
      !!!
      !!! ---------------------------------------
      !!!    ppm_param_kerel_laplace2d_2p
      !!!      (2nd order Laplacian,
      !!!      polynomial in 2D)
      !!!    ppm_param_kerel_laplace3d_2p
      !!!      (2nd order Laplacian,
      !!!      polynomial in 3D)
      !!! ---------------------------------------
      !!!
      !!! To use your own kernel function, pass the function pointer here. Your
      !!! function should take one argument and return one value.
#if   __KIND == __SINGLE_PRECISION | __KIND == __DOUBLE_PRECISION
      REAL(MK)   , DIMENSION(:)  , INTENT(IN   ) :: kpar
      !!! Kernel parameters. See documentation or ppm_comp_pp_kernels.inc for
      !!! description. Type can be single or complex.
#elif __KIND == __SINGLE_PRECISION_COMPLEX | __KIND == __DOUBLE_PRECISION_COMPLEX
      COMPLEX(MK), DIMENSION(:)  , INTENT(IN   ) :: kpar
      !!! Kernel parameters. See documentation or ppm_comp_pp_kernels.inc for
      !!! description. Type can be single complex or double complex.
#endif
#elif __KERNEL == __LOOKUP_TABLE
      REAL(MK)                   , INTENT(IN   ) :: kpar
      !!! Kernel parameters. See documentation or ppm_comp_pp_kernels.inc for
      !!! description. Pass dxtableinv (the inverse of the table spacing) as
      !!! a scalar here.
#if   __KIND == __SINGLE_PRECISION | __KIND == __DOUBLE_PRECISION
      REAL(MK)   , DIMENSION(:  ), INTENT(IN   ) :: kernel
#elif __KIND == __SINGLE_PRECISION_COMPLEX | __KIND == __DOUBLE_PRECISION_COMPLEX
      COMPLEX(MK), DIMENSION(:  ), INTENT(IN   ) :: kernel
#endif
      !!! Lookup table with tabulated kernel values.
      !!! Such a table can be created using <<ppm_comp_pp_mk_table>>.
#endif
      INTEGER    , DIMENSION(:,:), INTENT(IN   ) :: vlist
      !!! Verlet lists of all particles. 1st index: particles to interact with
      !!! (1..nvlist), 2nd: particle (1..Np).
      INTEGER    , DIMENSION(  :), INTENT(IN   ) :: nvlist
      !!! length of the Verlet lists of all particles
      INTEGER                    , INTENT(  OUT) :: info
      !!! Returns status, 0 upon success
      !-------------------------------------------------------------------------
      !  Local variables 
      !-------------------------------------------------------------------------
      INTEGER                                    :: jpart,ip,jp,ispec,idx
      REAL(MK)                                   :: dx,dy,dz
      REAL(MK)                                   :: factor,factor2
      REAL(MK)                                   :: dij,dij2,dij4,dij5
#if   __KIND == __SINGLE_PRECISION | __KIND == __DOUBLE_PRECISION
      REAL(MK)                                   :: summ,summ2
      REAL(MK)                                   :: eta,dm
#elif __KIND == __SINGLE_PRECISION_COMPLEX | __KIND == __DOUBLE_PRECISION_COMPLEX
      COMPLEX(MK)                                :: summ,summ2
      COMPLEX(MK)                                :: eta,dm
#endif
      ! timer 
      REAL(MK)                                   :: t0
      CHARACTER(LEN=ppm_char)                    :: mesg
      !-------------------------------------------------------------------------
      !  Externals 
      !-------------------------------------------------------------------------
#if   __KERNEL == __USER_FUNCTION
#if   __KIND == __SINGLE_PRECISION | __KIND == __DOUBLE_PRECISION
      REAL(MK)   , DIMENSION(:)  , INTENT(IN   ) :: kpar
      INTERFACE
          FUNCTION kernel(x,kpar)
              IMPLICIT NONE
#if   __KIND == __SINGLE_PRECISION | __KIND == __SINGLE_PRECISION_COMPLEX
              INTEGER, PARAMETER :: MK = KIND(1.0E0)
#elif __KIND == __DOUBLE_PRECISION | __KIND == __DOUBLE_PRECISION_COMPLEX
              INTEGER, PARAMETER :: MK = KIND(1.0D0)
#endif
              REAL(MK), INTENT(IN) :: x
              REAL(MK), DIMENSION(:), INTENT(IN) :: kpar
              REAL(MK) :: kernel
          END FUNCTION kernel
      END INTERFACE
#elif __KIND == __SINGLE_PRECISION_COMPLEX | __KIND == __DOUBLE_PRECISION_COMPLEX
      COMPLEX(MK), DIMENSION(:)  , INTENT(IN   ) :: kpar
      INTERFACE
          FUNCTION kernel(x,kpar)
              IMPLICIT NONE
#if   __KIND == __SINGLE_PRECISION | __KIND == __SINGLE_PRECISION_COMPLEX
              INTEGER, PARAMETER :: MK = KIND(1.0E0)
#elif __KIND == __DOUBLE_PRECISION | __KIND == __DOUBLE_PRECISION_COMPLEX
              INTEGER, PARAMETER :: MK = KIND(1.0D0)
#endif
              REAL(MK), INTENT(IN) :: x
              COMPLEX(MK), DIMENSION(:), INTENT(IN) :: kpar
              COMPLEX(MK) :: kernel
          END FUNCTION kernel
      END INTERFACE
#endif
#endif

      !-------------------------------------------------------------------------
      !  Initialise
      !-------------------------------------------------------------------------
      CALL substart('ppm_comp_pp_verlet',t0,info)

      !-------------------------------------------------------------------------
      !  Check arguments.
      !-------------------------------------------------------------------------
      IF (ppm_debug .GT. 0) THEN
          IF (.NOT. ppm_initialized) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_ppm_noinit,'ppm_comp_pp_verlet',  &
     &            'Please call ppm_init first!',__LINE__,info)
              GOTO 9999
          ENDIF
          IF (Np .LT. 0) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_comp_pp_verlet',  &
     &            'Np must be >=0',__LINE__,info)
              GOTO 9999
          ENDIF
          IF (lda .LT. 1) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_comp_pp_verlet',  &
     &            'lda must be >0',__LINE__,info)
              GOTO 9999
          ENDIF
#if   __KERNEL == __LOOKUP_TABLE
          IF (kpar .LT. 0.0_MK) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_comp_pp_verlet',  &
     &            'kpar (dxtableinv) must be >=0',__LINE__,info)
              GOTO 9999
          ENDIF
#endif
      ENDIF

      !-------------------------------------------------------------------------
      !  PARTICLE-PARTICLE INTERACTIONS using symmetry
      !-------------------------------------------------------------------------
      IF (lsymm) THEN
          IF (lda .EQ. 1) THEN
              DO ip=1,Np
                  summ = 0.0_MK
#ifdef __SXF90
!CDIR NODEP
#endif
                  DO jpart=1,nvlist(ip)
                      jp = vlist(jpart,ip)
                      !---------------------------------------------------------
                      !  Calculate the square of the distance between the two 
                      !  particles. It will always be .LE. (cutoff+skin)**2 by 
                      !  construction of the Verlet list.
                      !---------------------------------------------------------
                      dx  = xp(1,ip) - xp(1,jp)
                      dy  = xp(2,ip) - xp(2,jp)
                      IF (ppm_dim .GT. 2) THEN
                          dz  = xp(3,ip) - xp(3,jp)
                          dij = (dx*dx) + (dy*dy) + (dz*dz)
                      ELSE
                          dz = 0.0_MK
                          dij = (dx*dx) + (dy*dy)
                      ENDIF
                      !---------------------------------------------------------
                      !  Particle ip interacts with jp and vice versa
                      !---------------------------------------------------------
#include "ppm_comp_pp_kernels.inc"
                      dm        = eta*(pdata(1,jp) - pdata(1,ip))
                      summ      = summ      + dm
                      dpd(1,jp) = dpd(1,jp) - dm
                  ENDDO
                  dpd(1,ip) = dpd(1,ip) + summ
              ENDDO
          ELSEIF (lda .EQ. 2) THEN
              DO ip=1,Np
                  summ = 0.0_MK
                  summ2 = 0.0_MK
#ifdef __SXF90
!CDIR NODEP
#endif
                  DO jpart=1,nvlist(ip)
                      jp = vlist(jpart,ip)
                      !---------------------------------------------------------
                      !  Calculate the square of the distance between the two 
                      !  particles. It will always be .LE. (cutoff+skin)**2 by 
                      !  construction of the Verlet list.
                      !---------------------------------------------------------
                      dx  = xp(1,ip) - xp(1,jp)
                      dy  = xp(2,ip) - xp(2,jp)
                      IF (ppm_dim .GT. 2) THEN
                          dz  = xp(3,ip) - xp(3,jp)
                          dij = (dx*dx) + (dy*dy) + (dz*dz)
                      ELSE
                          dz = 0.0_MK
                          dij = (dx*dx) + (dy*dy)
                      ENDIF
                      !---------------------------------------------------------
                      !  Particle ip interacts with jp and vice versa
                      !---------------------------------------------------------
#include "ppm_comp_pp_kernels.inc"
                      dm        = eta*(pdata(1,jp) - pdata(1,ip))
                      summ      = summ      + dm
                      dpd(1,jp) = dpd(1,jp) - dm
                      dm        = eta*(pdata(2,jp) - pdata(2,ip))
                      summ2     = summ2     + dm
                      dpd(2,jp) = dpd(2,jp) - dm
                  ENDDO
                  dpd(1,ip) = dpd(1,ip) + summ
                  dpd(2,ip) = dpd(2,ip) + summ2
              ENDDO
          ELSE
              DO ip=1,Np
                  DO jpart=1,nvlist(ip)
                      jp = vlist(jpart,ip)
                      !---------------------------------------------------------
                      !  Calculate the square of the distance between the two 
                      !  particles. It will always be .LE. (cutoff+skin)**2 by 
                      !  construction of the Verlet list.
                      !---------------------------------------------------------
                      dx  = xp(1,ip) - xp(1,jp)
                      dy  = xp(2,ip) - xp(2,jp)
                      IF (ppm_dim .GT. 2) THEN
                          dz  = xp(3,ip) - xp(3,jp)
                          dij = (dx*dx) + (dy*dy) + (dz*dz)
                      ELSE
                          dz = 0.0_MK
                          dij = (dx*dx) + (dy*dy)
                      ENDIF
                      !---------------------------------------------------------
                      !  Particle ip interacts with jp and vice versa
                      !---------------------------------------------------------
#include "ppm_comp_pp_kernels.inc"
                      DO ispec=1,lda
                          dm   = eta*(pdata(ispec,jp) - pdata(ispec,ip))
                          dpd(ispec,ip) = dpd(ispec,ip) + dm
                          dpd(ispec,jp) = dpd(ispec,jp) - dm
                      ENDDO
                  ENDDO
              ENDDO
          ENDIF

      !-------------------------------------------------------------------------
      !  PARTICLE-PARTICLE INTERACTIONS not using symmetry
      !-------------------------------------------------------------------------
      ELSE
          IF (lda .EQ. 1) THEN
              DO ip=1,Np
#ifdef __SXF90
!CDIR NODEP
#endif
                  DO jpart=1,nvlist(ip)
                      jp = vlist(jpart,ip)
                      !---------------------------------------------------------
                      !  Calculate the square of the distance between the two 
                      !  particles. It will always be .LE. (cutoff+skin)**2 by 
                      !  construction of the Verlet list.
                      !---------------------------------------------------------
                      dx  = xp(1,ip) - xp(1,jp)
                      dy  = xp(2,ip) - xp(2,jp)
                      IF (ppm_dim .GT. 2) THEN
                          dz  = xp(3,ip) - xp(3,jp)
                          dij = (dx*dx) + (dy*dy) + (dz*dz)
                      ELSE
                          dz = 0.0_MK
                          dij = (dx*dx) + (dy*dy)
                      ENDIF
                      !---------------------------------------------------------
                      !  Particle ip interacts with jp
                      !---------------------------------------------------------
#include "ppm_comp_pp_kernels.inc"
                      dm = eta*(pdata(1,jp) - pdata(1,ip))
                      dpd(1,ip) = dpd(1,ip) + dm
                  ENDDO
              ENDDO
          ELSEIF (lda .EQ. 2) THEN
              DO ip=1,Np
#ifdef __SXF90
!CDIR NODEP
#endif
                  DO jpart=1,nvlist(ip)
                      jp = vlist(jpart,ip)
                      !---------------------------------------------------------
                      !  Calculate the square of the distance between the two 
                      !  particles. It will always be .LE. (cutoff+skin)**2 by 
                      !  construction of the Verlet list.
                      !---------------------------------------------------------
                      dx  = xp(1,ip) - xp(1,jp)
                      dy  = xp(2,ip) - xp(2,jp)
                      IF (ppm_dim .GT. 2) THEN
                          dz  = xp(3,ip) - xp(3,jp)
                          dij = (dx*dx) + (dy*dy) + (dz*dz)
                      ELSE
                          dz = 0.0_MK
                          dij = (dx*dx) + (dy*dy)
                      ENDIF
                      !---------------------------------------------------------
                      !  Particle ip interacts with jp
                      !---------------------------------------------------------
#include "ppm_comp_pp_kernels.inc"
                      dm = eta*(pdata(1,jp) - pdata(1,ip))
                      dpd(1,ip) = dpd(1,ip) + dm
                      dm = eta*(pdata(2,jp) - pdata(2,ip))
                      dpd(2,ip) = dpd(2,ip) + dm
                  ENDDO
              ENDDO
          ELSE
              DO ip=1,Np
                  DO jpart=1,nvlist(ip)
                      jp = vlist(jpart,ip)
                      !---------------------------------------------------------
                      !  Calculate the square of the distance between the two 
                      !  particles. It will always be .LE. (cutoff+skin)**2 by 
                      !  construction of the Verlet list.
                      !---------------------------------------------------------
                      dx  = xp(1,ip) - xp(1,jp)
                      dy  = xp(2,ip) - xp(2,jp)
                      IF (ppm_dim .GT. 2) THEN
                          dz  = xp(3,ip) - xp(3,jp)
                          dij = (dx*dx) + (dy*dy) + (dz*dz)
                      ELSE
                          dz = 0.0_MK
                          dij = (dx*dx) + (dy*dy)
                      ENDIF
                      !---------------------------------------------------------
                      !  Particle ip interacts with jp
                      !---------------------------------------------------------
#include "ppm_comp_pp_kernels.inc"
                      DO ispec=1,lda
                          dm = eta*(pdata(ispec,jp) - pdata(ispec,ip))
                          dpd(ispec,ip) = dpd(ispec,ip) + dm
                      ENDDO
                  ENDDO
              ENDDO
          ENDIF
      ENDIF

      !-------------------------------------------------------------------------
      !  Return
      !-------------------------------------------------------------------------
 9999 CONTINUE
      CALL substop('ppm_comp_pp_verlet',t0,info)
      RETURN
#if   __KERNEL == __INTERNAL
#if   __KIND == __SINGLE_PRECISION
      END SUBROUTINE ppm_comp_pp_verlet_si
#elif __KIND == __DOUBLE_PRECISION
      END SUBROUTINE ppm_comp_pp_verlet_di
#elif __KIND == __SINGLE_PRECISION_COMPLEX
      END SUBROUTINE ppm_comp_pp_verlet_sci
#elif __KIND == __DOUBLE_PRECISION_COMPLEX
      END SUBROUTINE ppm_comp_pp_verlet_dci
#endif

#elif __KERNEL == __USER_FUNCTION
#if   __KIND == __SINGLE_PRECISION
      END SUBROUTINE ppm_comp_pp_verlet_su
#elif __KIND == __DOUBLE_PRECISION
      END SUBROUTINE ppm_comp_pp_verlet_du
#elif __KIND == __SINGLE_PRECISION_COMPLEX
      END SUBROUTINE ppm_comp_pp_verlet_scu
#elif __KIND == __DOUBLE_PRECISION_COMPLEX
      END SUBROUTINE ppm_comp_pp_verlet_dcu
#endif

#elif __KERNEL == __LOOKUP_TABLE
#if   __KIND == __SINGLE_PRECISION
      END SUBROUTINE ppm_comp_pp_verlet_st
#elif __KIND == __DOUBLE_PRECISION
      END SUBROUTINE ppm_comp_pp_verlet_dt
#elif __KIND == __SINGLE_PRECISION_COMPLEX
      END SUBROUTINE ppm_comp_pp_verlet_sct
#elif __KIND == __DOUBLE_PRECISION_COMPLEX
      END SUBROUTINE ppm_comp_pp_verlet_dct
#endif
#endif
