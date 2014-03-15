      !-------------------------------------------------------------------------
      !  Module       :            ppm_module_ode_step
      !-------------------------------------------------------------------------
      !
      !  Purpose      : procedure module for ppm_ode_step
      !
      !  Remarks      : 
      !
      !  References   : 
      !
      !  Revisions    :
      !-------------------------------------------------------------------------
      !  $Log: ppm_module_ode_step.f,v $
      !  Revision 1.1.1.1  2007/07/13 10:19:00  ivos
      !  CBL version of the PPM library
      !
      !  Revision 1.5  2004/07/26 14:58:05  michaebe
      !  renamed the preprocessing defines vector and scalar
      !
      !  Revision 1.4  2004/07/26 12:00:24  ivos
      !  Fixes to make it compile.
      !
      !  Revision 1.3  2004/07/26 11:48:28  michaebe
      !  added vector and scalar defines
      !
      !  Revision 1.2  2004/07/26 08:14:15  michaebe
      !  Added overloading for scalar lda.
      !
      !  Revision 1.1  2004/07/26 07:45:50  michaebe
      !  Procedure modules created in the course of atomization.
      !
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

      !-------------------------------------------------------------------------
      !  Define data types
      !-------------------------------------------------------------------------
#define __SINGLE_PRECISION 1
#define __DOUBLE_PRECISION 2
#define __SCA 3
#define __VEC 4

      MODULE ppm_module_ode_step

        !-----------------------------------------------------
        !  Interface
        !-----------------------------------------------------

        INTERFACE ppm_ode_step
           MODULE PROCEDURE ppm_ode_step_ss
           MODULE PROCEDURE ppm_ode_step_ds
           MODULE PROCEDURE ppm_ode_step_sv
           MODULE PROCEDURE ppm_ode_step_dv
        END INTERFACE

#ifdef __F2003
        !-----------------------------------------------------
        !  Abstract RHS Interface
        !-----------------------------------------------------

        ABSTRACT INTERFACE

           ! single precision

           ! scalar
           FUNCTION rhsfunc_ss(topoid,xp,up,dup,lda,npart,ipack,&
                &lpack,rpack,info)
             INTEGER                          , INTENT(IN)  :: topoid
             INTEGER                          , INTENT(IN)  :: lda,npart
             INTEGER                          , INTENT(OUT) :: info
             REAL(KIND(1.0E0)), DIMENSION(:,:), POINTER     :: xp
             REAL(KIND(1.0E0)), DIMENSION(:),   POINTER     :: up
             REAL(KIND(1.0E0)), DIMENSION(:),   POINTER     :: dup
             INTEGER,  DIMENSION(:,:), INTENT(IN), OPTIONAL :: ipack
             LOGICAL,  DIMENSION(:,:), INTENT(IN), OPTIONAL :: lpack
             REAL(kind(1.0E0)), DIMENSION(:,:), INTENT(IN), OPTIONAL :: rpack
             INTEGER                                     :: rhsfunc
           END FUNCTION rhsfunc

           ! vector
           FUNCTION rhsfunc_sv(topoid,xp,up,dup,lda,npart,ipack,&
                &lpack,rpack,info)
             INTEGER                          , INTENT(IN)  :: topoid
             INTEGER                          , INTENT(IN)  :: lda,npart
             INTEGER                          , INTENT(OUT) :: info
             REAL(KIND(1.0E0)), DIMENSION(:,:), POINTER     :: xp,up
             REAL(KIND(1.0E0)), DIMENSION(:,:), POINTER     :: dup
             INTEGER,  DIMENSION(:,:), INTENT(IN), OPTIONAL :: ipack
             LOGICAL,  DIMENSION(:,:), INTENT(IN), OPTIONAL :: lpack
             REAL(kind(1.0E0)), DIMENSION(:,:), INTENT(IN), OPTIONAL :: rpack
             INTEGER                                     :: rhsfunc
           END FUNCTION rhsfunc

           ! double precision

           ! scalar
           FUNCTION rhsfunc_ds(topoid,xp,up,dup,lda,npart,ipack,&
                &lpack,rpack,info)
             INTEGER                          , INTENT(IN)  :: topoid
             INTEGER                          , INTENT(IN)  :: lda,npart
             INTEGER                          , INTENT(OUT) :: info
             REAL(KIND(1.0D0)), DIMENSION(:,:), POINTER     :: xp
             REAL(KIND(1.0D0)), DIMENSION(:),   POINTER     :: up
             REAL(KIND(1.0D0)), DIMENSION(:),   POINTER     :: dup
             INTEGER,  DIMENSION(:,:), INTENT(IN), OPTIONAL :: ipack
             LOGICAL,  DIMENSION(:,:), INTENT(IN), OPTIONAL :: lpack
             REAL(kind(1.0D0)), DIMENSION(:,:), INTENT(IN), OPTIONAL :: rpack
             INTEGER                                     :: rhsfunc
           END FUNCTION rhsfunc

           ! vector
           FUNCTION rhsfunc_dv(topoid,xp,up,dup,lda,npart,ipack,&
                &lpack,rpack,info)
             INTEGER                          , INTENT(IN)  :: topoid
             INTEGER                          , INTENT(IN)  :: lda,npart
             INTEGER                          , INTENT(OUT) :: info
             REAL(KIND(1.0D0)), DIMENSION(:,:), POINTER     :: xp,up
             REAL(KIND(1.0D0)), DIMENSION(:,:), POINTER     :: dup
             INTEGER,  DIMENSION(:,:), INTENT(IN), OPTIONAL :: ipack
             LOGICAL,  DIMENSION(:,:), INTENT(IN), OPTIONAL :: lpack
             REAL(kind(1.0D0)), DIMENSION(:,:), INTENT(IN), OPTIONAL :: rpack
             INTEGER                                     :: rhsfunc
           END FUNCTION rhsfunc

        END INTERFACE

#endif

      CONTAINS
#define __MODE __SCA
#define __KIND __SINGLE_PRECISION
#include "ppm_ode_step.f"
#undef  __KIND
#define __KIND __DOUBLE_PRECISION
#include "ppm_ode_step.f"
#undef  __KIND
#undef  __MODE
#define __MODE __VEC
#define __KIND __SINGLE_PRECISION
#include "ppm_ode_step.f"
#undef  __KIND
#define __KIND __DOUBLE_PRECISION
#include "ppm_ode_step.f"
#undef  __KIND
#undef  __MODE

      END MODULE ppm_module_ode_step


        
