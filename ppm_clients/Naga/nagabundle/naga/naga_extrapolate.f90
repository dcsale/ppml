#ifndef __skipheaderextrapolate
!------------------------------------------------------------------------------
! Subroutine :  naga_extrapolate.f90 
!------------------------------------------------------------------------------ 
!      _  __
!     / |/ /__ ____ ____ _
!    /    / _ `/ _ `/ _ `/
!   /_/|_/\_,_/\_, /\_,_/
!             /___/
!   A PPM Vortex-In-Cell client - 2011 (c)
!   Naga is distributed under the GNU Lesser General Public License version 3
!   Johannes Tophøj Rasmussen, johannes.tophoj@gmail.com
!   Technical University of Denmark
!   Department of Mechanical Engineering
!
!
! This routine extrapolates field values into the ghost layer
!------------------------------------------------------------------------------
MODULE naga_mod_extrapolate

INTERFACE naga_extrapolate
    MODULE PROCEDURE naga_extrapolate_s
    MODULE PROCEDURE naga_extrapolate_v
  END INTERFACE

  CONTAINS

#define __skipheaderextrapolate

#define __scalar
#define __ROUTINE naga_extrapolate_s
#include "naga_extrapolate.f90"
#undef __ROUTINE
#undef __scalar
#define __vector
#define __ROUTINE naga_extrapolate_v
#include "naga_extrapolate.f90"
#undef __ROUTINE
#undef __vector

END MODULE naga_mod_extrapolate

#else

SUBROUTINE __ROUTINE(pset,field,nextra,nbase,info)

USE naga_mod_globals
USE naga_mod_say
USE ppm_module_topo_get


IMPLICIT NONE
!-------------------------------------------------------------------------
! Arguments
!-------------------------------------------------------------------------
TYPE(patch_setup)                                           :: pset
!!!patch setting - dereferenced
!TYPE(patch_field_v)                                         :: field
#ifdef __scalar
REAL(MK),DIMENSION(:,:,:,:),POINTER                         :: field
#endif
#ifdef __vector
REAL(MK),DIMENSION(:,:,:,:,:),POINTER                       :: field
#endif
!!! Field to extrapolate ghost layer
INTEGER, INTENT(IN)                                         :: nextra
!!! Number of points to extrapolate into the ghostlayer
INTEGER, INTENT(IN)                                         :: nbase
!!! Number of points to base the extrapolation on
INTEGER, INTENT(OUT)                                        :: info
!!! Return state

!-------------------------------------------------------------------------
! Local variables
!-------------------------------------------------------------------------
REAL(MK)                             :: t0
TYPE(ppm_t_topo),POINTER             :: topology =>NULL()
INTEGER                              :: isub,isubl
INTEGER                              :: i,j,k,iextra,ibase
INTEGER                              :: ilevel,ipatch
REAL(MK)    ,DIMENSION(:,:),POINTER  :: coeff =>NULL()
#ifdef __scalar
REAL(MK)    ,DIMENSION(1)            :: tmpbuf
#endif
#ifdef __vector
REAL(MK)    ,DIMENSION(ncom)         :: tmpbuf
#endif


!-------------------------------------------------------------------------
! Initialise routine
!-------------------------------------------------------------------------
CALL substart('naga_extrapolate',t0,info)


!-------------------------------------------------------------------------
! Compare the number of points to extrapolate to the ghost layer width
!-------------------------------------------------------------------------
IF (nextra .GT. gstw(1) .OR. &
  & nextra .GT. gstw(2) .OR. &
  & nextra .GT. gstw(3)) THEN
  CALL naga_say(rank,'naga_extrapolate',&
  & 'The points to extrapolate exceeds the ghost layer.',info)
  info = -1
  GOTO 9999
ENDIF

!-------------------------------------------------------------------------
! Determine weights
!-------------------------------------------------------------------------
ALLOCATE(coeff(nbase,nextra))
IF (nbase .EQ. 4) THEN
  IF (nextra .GE. 1) THEN
    coeff(:,1) = (/4.0_MK,-6.0_MK,4.0_MK,-1.0_MK/)
  ENDIF
  IF (nextra .GE. 2) THEN
    coeff(:,2) = (/10.0_MK,-20.0_MK,15.0_MK,-4.0_MK/)
  ENDIF
  IF (nextra .GE. 3) THEN
    CALL naga_say(rank,'naga_extrapolate',&
    & 'Extrapolation to more than two points has not been implemented.',info)
    info = -1
    GOTO 9999
  ENDIF
ELSE
  CALL naga_say(rank,'naga_extrapolate',&
  & 'Only extrapolation based on 4 points has been implemented.',info)
  info = -1
  GOTO 9999
ENDIF


!-------------------------------------------------------------------------
! Extrapolate field into ghost layer
! The indicies of subs_bc represent: 
! west,east(x),south,north(y),bottom,top(z)
!@ Some more unrolling here would be nice
!-------------------------------------------------------------------------
DO ilevel=1,nlevels
  DO ipatch=1,npatches(ilevel)
    !-------------------------------------------------------------------------
    ! Get topology
    !-------------------------------------------------------------------------
    CALL ppm_topo_get(pset%topoid,topology,info)
    IF (info .NE. 0) THEN
      CALL naga_say(rank,'naga_extrapolate', 'Failed to get topology.',isub)
      GOTO 9999
    ENDIF

    DO isub=1,topos(ilevel)%nsublist
      isubl=topos(ilevel)%isublist(isub)
      !West (-x)
      IF (topology%subs_bc(1,isubl) .EQ. 1) THEN
        DO k=1-gstw(3),pset%snx(3,isubl)+gstw(3)
          DO j=1-gstw(2),pset%snx(2,isubl)+gstw(2)
            i = 1
            DO iextra=1,nextra
              tmpbuf = 0.0_MK
#ifdef __vector
              DO ibase=0,nbase-1
                tmpbuf(1) = tmpbuf(1) + &
                & coeff(ibase+1,iextra)*field(1,i+ibase,j,k,isub)
                tmpbuf(2) = tmpbuf(2) + &
                & coeff(ibase+1,iextra)*field(2,i+ibase,j,k,isub)
                tmpbuf(3) = tmpbuf(3) + &
                & coeff(ibase+1,iextra)*field(3,i+ibase,j,k,isub)
              END DO !ibase
              field(1,i-iextra,j,k,isub) = tmpbuf(1)
              field(2,i-iextra,j,k,isub) = tmpbuf(2)
              field(3,i-iextra,j,k,isub) = tmpbuf(3)
#endif
#ifdef __scalar
              DO ibase=0,nbase-1
                tmpbuf(1) = tmpbuf(1) + &
                & coeff(ibase+1,iextra)*field(i+ibase,j,k,isub)
              END DO !ibase
              field(i-iextra,j,k,isub) = tmpbuf(1)
#endif
            END DO !iextra
          ENDDO !j
        ENDDO !k
      ENDIF
      !East (+x)
      IF (topology%subs_bc(2,isubl) .EQ. 1) THEN
        DO k=1-gstw(3),pset%snx(3,isubl)+gstw(3)
          DO j=1-gstw(2),pset%snx(2,isubl)+gstw(2)
            i = pset%snx(1,isubl)
            DO iextra=1,nextra
              tmpbuf = 0.0_MK
#ifdef __vector
              DO ibase=0,nbase-1
                tmpbuf(1) = tmpbuf(1) + &
                & coeff(ibase+1,iextra)*field(1,i-ibase,j,k,isub)
                tmpbuf(2) = tmpbuf(2) + &
                & coeff(ibase+1,iextra)*field(2,i-ibase,j,k,isub)
                tmpbuf(3) = tmpbuf(3) + &
                & coeff(ibase+1,iextra)*field(3,i-ibase,j,k,isub)
              END DO !ibase
              field(1,i+iextra,j,k,isub) = tmpbuf(1)
              field(2,i+iextra,j,k,isub) = tmpbuf(2)
              field(3,i+iextra,j,k,isub) = tmpbuf(3)
#endif
#ifdef __scalar
              DO ibase=0,nbase-1
                tmpbuf(1) = tmpbuf(1) + &
                & coeff(ibase+1,iextra)*field(i-ibase,j,k,isub)
              END DO !ibase
              field(i+iextra,j,k,isub) = tmpbuf(1)
#endif
            END DO !iextra
          ENDDO !j
        ENDDO !k
      ENDIF
      !South (-y)
      IF (topology%subs_bc(3,isubl) .EQ. 1) THEN
        DO k=1-gstw(3),pset%snx(3,isubl)+gstw(3)
          DO i=1-gstw(1),pset%snx(1,isubl)+gstw(1)
            j = 1
            DO iextra=1,nextra
              tmpbuf = 0.0_MK
#ifdef __vector
              DO ibase=0,nbase-1
                tmpbuf(1) = tmpbuf(1) + &
                & coeff(ibase+1,iextra)*field(1,i,j+ibase,k,isub)
                tmpbuf(2) = tmpbuf(2) + &
                & coeff(ibase+1,iextra)*field(2,i,j+ibase,k,isub)
                tmpbuf(3) = tmpbuf(3) + &
                & coeff(ibase+1,iextra)*field(3,i,j+ibase,k,isub)
              END DO !ibase
              field(1,i,j-iextra,k,isub) = tmpbuf(1)
              field(2,i,j-iextra,k,isub) = tmpbuf(2)
              field(3,i,j-iextra,k,isub) = tmpbuf(3)
#endif
#ifdef __scalar
              DO ibase=0,nbase-1
                tmpbuf(1) = tmpbuf(1) + &
                & coeff(ibase+1,iextra)*field(i,j+ibase,k,isub)
              END DO !ibase
              field(i,j-iextra,k,isub) = tmpbuf(1)
#endif
            END DO !iextra
          ENDDO !i
        ENDDO !k
      ENDIF
      !North (+y)
      IF (topology%subs_bc(4,isubl) .EQ. 1) THEN
        DO k=1-gstw(3),pset%snx(3,isubl)+gstw(3)
          DO i=1-gstw(1),pset%snx(1,isubl)+gstw(1)
            j = pset%snx(2,isubl)
            DO iextra=1,nextra
              tmpbuf = 0.0_MK
#ifdef __vector
              DO ibase=0,nbase-1
                tmpbuf(1) = tmpbuf(1) + &
                & coeff(ibase+1,iextra)*field(1,i,j-ibase,k,isub)
                tmpbuf(2) = tmpbuf(2) + &
                & coeff(ibase+1,iextra)*field(2,i,j-ibase,k,isub)
                tmpbuf(3) = tmpbuf(3) + &
                & coeff(ibase+1,iextra)*field(3,i,j-ibase,k,isub)
              END DO !ibase
              field(1,i,j+iextra,k,isub) = tmpbuf(1)
              field(2,i,j+iextra,k,isub) = tmpbuf(2)
              field(3,i,j+iextra,k,isub) = tmpbuf(3)
#endif
#ifdef __scalar
              DO ibase=0,nbase-1
                tmpbuf(1) = tmpbuf(1) + &
                & coeff(ibase+1,iextra)*field(i,j-ibase,k,isub)
              END DO !ibase
              field(i,j+iextra,k,isub) = tmpbuf(1)
#endif
            END DO !iextra
          ENDDO !i
        ENDDO !k
      ENDIF
      !Bottom (-z)
      IF (topology%subs_bc(5,isubl) .EQ. 1) THEN
        DO j=1-gstw(2),pset%snx(2,isubl)+gstw(2)
          DO i=1-gstw(1),pset%snx(1,isubl)+gstw(1)
            k = 1
            DO iextra=1,nextra
              tmpbuf = 0.0_MK
#ifdef __vector
              DO ibase=0,nbase-1
                tmpbuf(1) = tmpbuf(1) + &
                & coeff(ibase+1,iextra)*field(1,i,j,k+ibase,isub)
                tmpbuf(2) = tmpbuf(2) + &
                & coeff(ibase+1,iextra)*field(2,i,j,k+ibase,isub)
                tmpbuf(3) = tmpbuf(3) + &
                & coeff(ibase+1,iextra)*field(3,i,j,k+ibase,isub)
              END DO !ibase
              field(1,i,j,k-iextra,isub) = tmpbuf(1)
              field(2,i,j,k-iextra,isub) = tmpbuf(2)
              field(3,i,j,k-iextra,isub) = tmpbuf(3)
#endif
#ifdef __scalar
              DO ibase=0,nbase-1
                tmpbuf(1) = tmpbuf(1) + &
                & coeff(ibase+1,iextra)*field(i,j,k+ibase,isub)
              END DO !ibase
              field(i,j,k-iextra,isub) = tmpbuf(1)
#endif
            END DO !iextra
          ENDDO !i
        ENDDO !j
      ENDIF
      !Top (+z)
      IF (topology%subs_bc(6,isubl) .EQ. 1) THEN
        DO j=1-gstw(2),pset%snx(2,isubl)+gstw(2)
          DO i=1-gstw(1),pset%snx(1,isubl)+gstw(1)
            k = pset%snx(3,isubl)
            DO iextra=1,nextra
              tmpbuf = 0.0_MK
#ifdef __vector
              DO ibase=0,nbase-1
                tmpbuf(1) = tmpbuf(1) + &
                & coeff(ibase+1,iextra)*field(1,i,j,k-ibase,isub)
                tmpbuf(2) = tmpbuf(2) + &
                & coeff(ibase+1,iextra)*field(2,i,j,k-ibase,isub)
                tmpbuf(3) = tmpbuf(3) + &
                & coeff(ibase+1,iextra)*field(3,i,j,k-ibase,isub)
              END DO !ibase
              field(1,i,j,k+iextra,isub) = tmpbuf(1)
              field(2,i,j,k+iextra,isub) = tmpbuf(2)
              field(3,i,j,k+iextra,isub) = tmpbuf(3)
#endif
#ifdef __scalar
              DO ibase=0,nbase-1
                tmpbuf(1) = tmpbuf(1) + &
                & coeff(ibase+1,iextra)*field(i,j,k-ibase,isub)
              END DO !ibase
              field(i,j,k+iextra,isub) = tmpbuf(1)
#endif
            END DO !iextra
          ENDDO !i
        ENDDO !j
      ENDIF
    ENDDO !isub
  ENDDO !ipatch
ENDDO !ilevel


9999 CONTINUE
CALL substop('ppm_poisson_extrapolateghost',t0,info)
RETURN

END SUBROUTINE __ROUTINE

#endif
