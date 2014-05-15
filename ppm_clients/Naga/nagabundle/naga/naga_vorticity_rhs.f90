!-------------------------------------------------------------------------------
! Subroutine :  naga_vorticity_rhs.f90
!-------------------------------------------------------------------------------
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
! This routines computes the rhs of the vorticity equation. Stretching diffusion
! 2nd or 4th order finite differences are used to evaluate the equation.
!-------------------------------------------------------------------------------
#define __WITHSTRETCHING

MODULE naga_mod_vorticity_rhs

IMPLICIT NONE

INTERFACE naga_vorticity_rhs
  MODULE PROCEDURE naga_vorticity_rhs
END INTERFACE

CONTAINS

SUBROUTINE naga_vorticity_rhs(info)

USE naga_mod_globals
USE naga_mod_say

IMPLICIT NONE

!-------------------------------------------------------------------------
! Arguments
!-------------------------------------------------------------------------
INTEGER, INTENT(OUT)          :: info


!-------------------------------------------------------------------------
! Local variables
!-------------------------------------------------------------------------
REAL(MK)                      :: t0
INTEGER                       :: ilevel,ipatch
REAL(MK)                      :: dx,dy,dz
REAL(MK)                      :: facx1,facy1,facz1
REAL(MK)                      :: facx2,facy2,facz2
REAL(MK)                      :: fac1,fac2,fac3,fac4,fac5,fac6,fac7,fac8,fac9
INTEGER                       :: isub,isubl
INTEGER                       :: i,j,k


!----------------------------------------------------------------------------
! Initialise routine
!----------------------------------------------------------------------------
CALL substart('naga_vorticity_rhs',t0,info)
info = 0

DO ilevel=1,nlevels
  DO ipatch=1,npatches(ilevel)
    dx = ptcset(ilevel,ipatch)%dx(1)
    dy = ptcset(ilevel,ipatch)%dx(2)
    dz = ptcset(ilevel,ipatch)%dx(3)

    !-----------------------------------------------------------------------
    ! No rhs
    !-----------------------------------------------------------------------
    IF (rhsscheme .EQ. 0) THEN
      dwf(ilevel,ipatch)%fld = 0.0_MK
    ENDIF


    !-----------------------------------------------------------------------
    ! Take curl, 2nd order FD
    !-----------------------------------------------------------------------
    IF (rhsscheme .EQ. 1) THEN
      facx1 = 1.0_MK/(2.0_MK*dx)
      facy1 = 1.0_MK/(2.0_MK*dy)
      facz1 = 1.0_MK/(2.0_MK*dz)
      facx2 = nu*1.0_MK/(dx**2)
      facy2 = nu*1.0_MK/(dy**2)
      facz2 = nu*1.0_MK/(dz**2)

      DO isub=1,topos(ilevel)%nsublist
        isubl=topos(ilevel)%isublist(isub)
        DO k=1,ptcset(ilevel,ipatch)%snx(3,isubl)
          DO j=1,ptcset(ilevel,ipatch)%snx(2,isubl)
            DO i=1,ptcset(ilevel,ipatch)%snx(1,isubl)
              dwf(ilevel,ipatch)%fld(1,i,j,k,isub) =              &
              &  facx2*(wf(ilevel,ipatch)%fld(1,i+1,j  ,k  ,isub) &
              &        +wf(ilevel,ipatch)%fld(1,i-1,j  ,k  ,isub))&
              & +facy2*(wf(ilevel,ipatch)%fld(1,i  ,j+1,k  ,isub) &
              &        +wf(ilevel,ipatch)%fld(1,i  ,j-1,k  ,isub))&
              & +facz2*(wf(ilevel,ipatch)%fld(1,i  ,j  ,k+1,isub) &
              &        +wf(ilevel,ipatch)%fld(1,i  ,j  ,k-1,isub))&
              &        -wf(ilevel,ipatch)%fld(1,i  ,j  ,k  ,isub) &
              & *2.0_MK*(facx2+facy2+facz2)                       &
              #ifdef __WITHSTRETCHING
              & +        wf(ilevel,ipatch)%fld(1,i  ,j  ,k  ,isub)  &
              & * facx1*(uf(ilevel,ipatch)%fld(1,i+1,j  ,k  ,isub)- &
                       & uf(ilevel,ipatch)%fld(1,i-1,j  ,k  ,isub)) &
              & +        wf(ilevel,ipatch)%fld(2,i  ,j  ,k  ,isub)  &
              & * facy1*(uf(ilevel,ipatch)%fld(1,i  ,j+1,k  ,isub)- &
                       & uf(ilevel,ipatch)%fld(1,i  ,j-1,k  ,isub)) &
              & +        wf(ilevel,ipatch)%fld(3,i  ,j  ,k  ,isub)  &
              & * facz1*(uf(ilevel,ipatch)%fld(1,i  ,j  ,k+1,isub)- &
                       & uf(ilevel,ipatch)%fld(1,i  ,j  ,k-1,isub))
              #else
              & + 0.0_MK
              #endif
              dwf(ilevel,ipatch)%fld(2,i,j,k,isub) =              &
              &  facx2*(wf(ilevel,ipatch)%fld(2,i+1,j  ,k  ,isub) &
              &        +wf(ilevel,ipatch)%fld(2,i-1,j  ,k  ,isub))&
              & +facy2*(wf(ilevel,ipatch)%fld(2,i  ,j+1,k  ,isub) &
              &        +wf(ilevel,ipatch)%fld(2,i  ,j-1,k  ,isub))&
              & +facz2*(wf(ilevel,ipatch)%fld(2,i  ,j  ,k+1,isub) &
              &        +wf(ilevel,ipatch)%fld(2,i  ,j  ,k-1,isub))&
              &        -wf(ilevel,ipatch)%fld(2,i  ,j  ,k  ,isub) &
              & *2.0_MK*(facx2+facy2+facz2)                       &
              #ifdef __WITHSTRETCHING
              & +        wf(ilevel,ipatch)%fld(1,i  ,j  ,k  ,isub)  &
              & * facx1*(uf(ilevel,ipatch)%fld(2,i+1,j  ,k  ,isub)- &
                       & uf(ilevel,ipatch)%fld(2,i-1,j  ,k  ,isub)) &
              & +        wf(ilevel,ipatch)%fld(2,i  ,j  ,k  ,isub)  &
              & * facy1*(uf(ilevel,ipatch)%fld(2,i  ,j+1,k  ,isub)- &
                       & uf(ilevel,ipatch)%fld(2,i  ,j-1,k  ,isub)) &
              & +        wf(ilevel,ipatch)%fld(3,i  ,j  ,k  ,isub)  &
              & * facz1*(uf(ilevel,ipatch)%fld(2,i  ,j  ,k+1,isub)- &
                       & uf(ilevel,ipatch)%fld(2,i  ,j  ,k-1,isub))
              #else
              & + 0.0_MK
              #endif
              dwf(ilevel,ipatch)%fld(3,i,j,k,isub) =              &
              &  facx2*(wf(ilevel,ipatch)%fld(3,i+1,j  ,k  ,isub) &
              &        +wf(ilevel,ipatch)%fld(3,i-1,j  ,k  ,isub))&
              & +facy2*(wf(ilevel,ipatch)%fld(3,i  ,j+1,k  ,isub) &
              &        +wf(ilevel,ipatch)%fld(3,i  ,j-1,k  ,isub))&
              & +facz2*(wf(ilevel,ipatch)%fld(3,i  ,j  ,k+1,isub) &
              &        +wf(ilevel,ipatch)%fld(3,i  ,j  ,k-1,isub))&
              &        -wf(ilevel,ipatch)%fld(3,i  ,j  ,k  ,isub) &
              & *2.0_MK*(facx2+facy2+facz2)                       &
              #ifdef __WITHSTRETCHING
              & +        wf(ilevel,ipatch)%fld(1,i  ,j  ,k  ,isub)  &
              & * facx1*(uf(ilevel,ipatch)%fld(3,i+1,j  ,k  ,isub)- &
                       & uf(ilevel,ipatch)%fld(3,i-1,j  ,k  ,isub)) &
              & +        wf(ilevel,ipatch)%fld(2,i  ,j  ,k  ,isub)  &
              & * facy1*(uf(ilevel,ipatch)%fld(3,i  ,j+1,k  ,isub)- &
                       & uf(ilevel,ipatch)%fld(3,i  ,j-1,k  ,isub)) &
              & +        wf(ilevel,ipatch)%fld(3,i  ,j  ,k  ,isub)  &
              & * facz1*(uf(ilevel,ipatch)%fld(3,i  ,j  ,k+1,isub)- &
                       & uf(ilevel,ipatch)%fld(3,i  ,j  ,k-1,isub))
              #else
              & + 0.0_MK
              #endif

            ENDDO
          ENDDO
        ENDDO
      ENDDO !isub


    !-----------------------------------------------------------------------
    ! Take curl, standard 4th order FD - not div(omega u)
    !-----------------------------------------------------------------------
    ELSE IF (rhsscheme .EQ. 2) THEN
      fac1 = 1.0_MK/(dx**2)*nu/12.0_MK
      fac2 = 1.0_MK/(dy**2)*nu/12.0_MK
      fac3 = 1.0_MK/(dz**2)*nu/12.0_MK

      facx1 = 1.0_MK/dx/12.0_MK
      facy1 = 1.0_MK/dy/12.0_MK
      facz1 = 1.0_MK/dz/12.0_MK

      DO isub=1,topos(ilevel)%nsublist
        isubl=topos(ilevel)%isublist(isub)
        DO k=1,ptcset(ilevel,ipatch)%snx(3,isubl)
          DO j=1,ptcset(ilevel,ipatch)%snx(2,isubl)
            DO i=1,ptcset(ilevel,ipatch)%snx(1,isubl)
              dwf(ilevel,ipatch)%fld(1,i,j,k,isub) = (&
                &        -fac1*(wf(ilevel,ipatch)%fld(1,i+2,j,k,isub)  &
                &              +wf(ilevel,ipatch)%fld(1,i-2,j,k,isub)) &
                &+16.0_MK*fac1*(wf(ilevel,ipatch)%fld(1,i+1,j,k,isub)  &
                &              +wf(ilevel,ipatch)%fld(1,i-1,j,k,isub)) &
                &        -fac2*(wf(ilevel,ipatch)%fld(1,i,j+2,k,isub)  &
                &              +wf(ilevel,ipatch)%fld(1,i,j-2,k,isub)) &
                &+16.0_MK*fac2*(wf(ilevel,ipatch)%fld(1,i,j+1,k,isub)  &
                &              +wf(ilevel,ipatch)%fld(1,i,j-1,k,isub)) &
                &        -fac3*(wf(ilevel,ipatch)%fld(1,i,j,k+2,isub)  &
                &              +wf(ilevel,ipatch)%fld(1,i,j,k-2,isub)) &
                &+16.0_MK*fac3*(wf(ilevel,ipatch)%fld(1,i,j,k+1,isub)  &
                &              +wf(ilevel,ipatch)%fld(1,i,j,k-1,isub)) &
                &-30.0_MK*(fac1+fac2+fac3)*wf(ilevel,ipatch)%fld(1,i,j,k,isub)) +&
              #ifdef __WITHSTRETCHING
                &               wf(ilevel,ipatch)%fld(1,i  ,j  ,k  ,isub)  *&
                &(       -facx1*uf(ilevel,ipatch)%fld(1,i+2,j  ,k  ,isub)   &
                & +8.0_MK*facx1*uf(ilevel,ipatch)%fld(1,i+1,j  ,k  ,isub)   &
                & -8.0_MK*facx1*uf(ilevel,ipatch)%fld(1,i-1,j  ,k  ,isub)   &
                &        +facx1*uf(ilevel,ipatch)%fld(1,i-2,j  ,k  ,isub)) +&
                &               wf(ilevel,ipatch)%fld(2,i  ,j  ,k  ,isub)  *&
                &(       -facy1*uf(ilevel,ipatch)%fld(1,i  ,j+2,k  ,isub)   &
                & +8.0_MK*facy1*uf(ilevel,ipatch)%fld(1,i  ,j+1,k  ,isub)   &
                & -8.0_MK*facy1*uf(ilevel,ipatch)%fld(1,i  ,j-1,k  ,isub)   &
                &        +facy1*uf(ilevel,ipatch)%fld(1,i  ,j-2,k  ,isub)) +&
                &               wf(ilevel,ipatch)%fld(3,i  ,j  ,k  ,isub)  *&
                &(       -facz1*uf(ilevel,ipatch)%fld(1,i  ,j  ,k+2,isub)   &
                & +8.0_MK*facz1*uf(ilevel,ipatch)%fld(1,i  ,j  ,k+1,isub)   &
                & -8.0_MK*facz1*uf(ilevel,ipatch)%fld(1,i  ,j  ,k-1,isub)   &
                &        +facz1*uf(ilevel,ipatch)%fld(1,i  ,j  ,k-2,isub))
              #else
                & 0.0_MK
              #endif

              dwf(ilevel,ipatch)%fld(2,i,j,k,isub) = (&
                &       -fac1*(wf(ilevel,ipatch)%fld(2,i+2,j,k,isub)  &
                &             +wf(ilevel,ipatch)%fld(2,i-2,j,k,isub)) &
                &+16.0_MK*fac1*(wf(ilevel,ipatch)%fld(2,i+1,j,k,isub) &
                &             +wf(ilevel,ipatch)%fld(2,i-1,j,k,isub)) &
                &       -fac2*(wf(ilevel,ipatch)%fld(2,i,j+2,k,isub)  &
                &             +wf(ilevel,ipatch)%fld(2,i,j-2,k,isub)) &
                &+16.0_MK*fac2*(wf(ilevel,ipatch)%fld(2,i,j+1,k,isub) &
                &             +wf(ilevel,ipatch)%fld(2,i,j-1,k,isub)) &
                &       -fac3*(wf(ilevel,ipatch)%fld(2,i,j,k+2,isub)  &
                &             +wf(ilevel,ipatch)%fld(2,i,j,k-2,isub)) &
                &+16.0_MK*fac3*(wf(ilevel,ipatch)%fld(2,i,j,k+1,isub) &
                &             +wf(ilevel,ipatch)%fld(2,i,j,k-1,isub)) &
                &-30.0_MK*(fac1+fac2+fac3)*wf(ilevel,ipatch)%fld(2,i,j,k,isub)) +&
              #ifdef __WITHSTRETCHING
                &               wf(ilevel,ipatch)%fld(1,i  ,j  ,k  ,isub)  *&
                &(       -facx1*uf(ilevel,ipatch)%fld(2,i+2,j  ,k  ,isub)   &
                & +8.0_MK*facx1*uf(ilevel,ipatch)%fld(2,i+1,j  ,k  ,isub)   &
                & -8.0_MK*facx1*uf(ilevel,ipatch)%fld(2,i-1,j  ,k  ,isub)   &
                &        +facx1*uf(ilevel,ipatch)%fld(2,i-2,j  ,k  ,isub)) +&
                &               wf(ilevel,ipatch)%fld(2,i  ,j  ,k  ,isub)  *&
                &(       -facy1*uf(ilevel,ipatch)%fld(2,i  ,j+2,k  ,isub)   &
                & +8.0_MK*facy1*uf(ilevel,ipatch)%fld(2,i  ,j+1,k  ,isub)   &
                & -8.0_MK*facy1*uf(ilevel,ipatch)%fld(2,i  ,j-1,k  ,isub)   &
                &        +facy1*uf(ilevel,ipatch)%fld(2,i  ,j-2,k  ,isub)) +&
                &               wf(ilevel,ipatch)%fld(3,i  ,j  ,k  ,isub)  *&
                &(       -facz1*uf(ilevel,ipatch)%fld(2,i  ,j  ,k+2,isub)   &
                & +8.0_MK*facz1*uf(ilevel,ipatch)%fld(2,i  ,j  ,k+1,isub)   &
                & -8.0_MK*facz1*uf(ilevel,ipatch)%fld(2,i  ,j  ,k-1,isub)   &
                &        +facz1*uf(ilevel,ipatch)%fld(2,i  ,j  ,k-2,isub))
              #else
                & 0.0_MK
              #endif

              dwf(ilevel,ipatch)%fld(3,i,j,k,isub) = (&
                &       -fac1*(wf(ilevel,ipatch)%fld(3,i+2,j,k,isub)  &
                &             +wf(ilevel,ipatch)%fld(3,i-2,j,k,isub)) &
                &+16.0_MK*fac1*(wf(ilevel,ipatch)%fld(3,i+1,j,k,isub) &
                &             +wf(ilevel,ipatch)%fld(3,i-1,j,k,isub)) &
                &       -fac2*(wf(ilevel,ipatch)%fld(3,i,j+2,k,isub)  &
                &             +wf(ilevel,ipatch)%fld(3,i,j-2,k,isub)) &
                &+16.0_MK*fac2*(wf(ilevel,ipatch)%fld(3,i,j+1,k,isub) &
                &             +wf(ilevel,ipatch)%fld(3,i,j-1,k,isub)) &
                &       -fac3*(wf(ilevel,ipatch)%fld(3,i,j,k+2,isub)  &
                &             +wf(ilevel,ipatch)%fld(3,i,j,k-2,isub)) &
                &+16.0_MK*fac3*(wf(ilevel,ipatch)%fld(3,i,j,k+1,isub) &
                &             +wf(ilevel,ipatch)%fld(3,i,j,k-1,isub)) &
                &-30.0_MK*(fac1+fac2+fac3)*wf(ilevel,ipatch)%fld(3,i,j,k,isub)) +&
              #ifdef __WITHSTRETCHING
                &               wf(ilevel,ipatch)%fld(1,i  ,j  ,k  ,isub)  *&
                &(       -facx1*uf(ilevel,ipatch)%fld(3,i+2,j  ,k  ,isub)   &
                & +8.0_MK*facx1*uf(ilevel,ipatch)%fld(3,i+1,j  ,k  ,isub)   &
                & -8.0_MK*facx1*uf(ilevel,ipatch)%fld(3,i-1,j  ,k  ,isub)   &
                &        +facx1*uf(ilevel,ipatch)%fld(3,i-2,j  ,k  ,isub)) +&
                &               wf(ilevel,ipatch)%fld(2,i  ,j  ,k  ,isub)  *&
                &(       -facy1*uf(ilevel,ipatch)%fld(3,i  ,j+2,k  ,isub)   &
                & +8.0_MK*facy1*uf(ilevel,ipatch)%fld(3,i  ,j+1,k  ,isub)   &
                & -8.0_MK*facy1*uf(ilevel,ipatch)%fld(3,i  ,j-1,k  ,isub)   &
                &        +facy1*uf(ilevel,ipatch)%fld(3,i  ,j-2,k  ,isub)) +&
                &               wf(ilevel,ipatch)%fld(3,i  ,j  ,k  ,isub)  *&
                &(       -facz1*uf(ilevel,ipatch)%fld(3,i  ,j  ,k+2,isub)   &
                & +8.0_MK*facz1*uf(ilevel,ipatch)%fld(3,i  ,j  ,k+1,isub)   &
                & -8.0_MK*facz1*uf(ilevel,ipatch)%fld(3,i  ,j  ,k-1,isub)   &
                &        +facz1*uf(ilevel,ipatch)%fld(3,i  ,j  ,k-2,isub))
              #else
                & 0.0_MK
              #endif

            ENDDO
          ENDDO
        ENDDO
      ENDDO !isub


    !-----------------------------------------------------------------------
    ! Take curl, 4th order FD conservative - untested in this implementation
    !-----------------------------------------------------------------------
    ELSE IF (rhsscheme .EQ. 3) THEN
      fac1 = 1.0_MK/dx**2*nu/12.0_MK
      fac2 = 1.0_MK/dy**2*nu/12.0_MK
      fac3 = 1.0_MK/dz**2*nu/12.0_MK

      fac4 = 8.0_MK/dx/12.0_MK
      fac5 = 8.0_MK/dy/12.0_MK
      fac6 = 8.0_MK/dz/12.0_MK

      fac7 = 1.0_MK/dx/12.0_MK
      fac8 = 1.0_MK/dy/12.0_MK
      fac9 = 1.0_MK/dz/12.0_MK

      DO isub=1,topos(ilevel)%nsublist
        isubl=topos(ilevel)%isublist(isub)
        DO k=1,ptcset(ilevel,ipatch)%snx(3,isubl)
          DO j=1,ptcset(ilevel,ipatch)%snx(2,isubl)
            DO i=1,ptcset(ilevel,ipatch)%snx(1,isubl)
              dwf(ilevel,ipatch)%fld(1,i,j,k,isub) = (&
                &16.0_MK*fac1*(wf(ilevel,ipatch)%fld(1,i+1,j,k,isub)  &
                &             +wf(ilevel,ipatch)%fld(1,i-1,j,k,isub))-&
                &        fac1*(wf(ilevel,ipatch)%fld(1,i+2,j,k,isub)  &
                &             +wf(ilevel,ipatch)%fld(1,i-2,j,k,isub))+&
                &16.0_MK*fac2*(wf(ilevel,ipatch)%fld(1,i,j+1,k,isub)  &
                &             +wf(ilevel,ipatch)%fld(1,i,j-1,k,isub))-&
                &        fac2*(wf(ilevel,ipatch)%fld(1,i,j+2,k,isub)  &
                &             +wf(ilevel,ipatch)%fld(1,i,j-2,k,isub))+&
                &16.0_MK*fac3*(wf(ilevel,ipatch)%fld(1,i,j,k+1,isub)  &
                &             +wf(ilevel,ipatch)%fld(1,i,j,k-1,isub))-&
                &        fac3*(wf(ilevel,ipatch)%fld(1,i,j,k+2,isub)  &
                &             +wf(ilevel,ipatch)%fld(1,i,j,k-2,isub))-&
                &30.0_MK*(fac1+fac2+fac3)*wf(ilevel,ipatch)%fld(1,i,j,k,isub)) +&
                & fac4*(wf(ilevel,ipatch)%fld(1,i+1,j,k,isub) &
                &      *uf(ilevel,ipatch)%fld(1,i+1,j,k,isub)-&
                &       wf(ilevel,ipatch)%fld(1,i-1,j,k,isub) &
                &      *uf(ilevel,ipatch)%fld(1,i-1,j,k,isub))+&
                & fac5*(wf(ilevel,ipatch)%fld(2,i,j+1,k,isub) &
                &      *uf(ilevel,ipatch)%fld(1,i,j+1,k,isub)-&
                &       wf(ilevel,ipatch)%fld(2,i,j-1,k,isub) &
                &      *uf(ilevel,ipatch)%fld(1,i,j-1,k,isub))+&
                & fac6*(wf(ilevel,ipatch)%fld(3,i,j,k+1,isub) &
                &      *uf(ilevel,ipatch)%fld(1,i,j,k+1,isub)-&
                &       wf(ilevel,ipatch)%fld(3,i,j,k-1,isub) &
                &      *uf(ilevel,ipatch)%fld(1,i,j,k-1,isub))-&
                & fac7*(wf(ilevel,ipatch)%fld(1,i+2,j,k,isub) &
                &      *uf(ilevel,ipatch)%fld(1,i+2,j,k,isub)-&
                &       wf(ilevel,ipatch)%fld(1,i-2,j,k,isub) &
                &      *uf(ilevel,ipatch)%fld(1,i-2,j,k,isub))-&
                & fac8*(wf(ilevel,ipatch)%fld(2,i,j+2,k,isub) &
                &      *uf(ilevel,ipatch)%fld(1,i,j+2,k,isub)-&
                &       wf(ilevel,ipatch)%fld(2,i,j-2,k,isub) &
                &      *uf(ilevel,ipatch)%fld(1,i,j-2,k,isub))-&
                & fac9*(wf(ilevel,ipatch)%fld(3,i,j,k+2,isub) &
                &      *uf(ilevel,ipatch)%fld(1,i,j,k+2,isub)-&
                &       wf(ilevel,ipatch)%fld(3,i,j,k-2,isub) &
                &      *uf(ilevel,ipatch)%fld(1,i,j,k-2,isub))

              dwf(ilevel,ipatch)%fld(2,i,j,k,isub) = (&
                &16.0_MK*fac1*(wf(ilevel,ipatch)%fld(2,i+1,j,k,isub)  &
                &             +wf(ilevel,ipatch)%fld(2,i-1,j,k,isub))-&
                &        fac1*(wf(ilevel,ipatch)%fld(2,i+2,j,k,isub)  &
                &             +wf(ilevel,ipatch)%fld(2,i-2,j,k,isub))+&
                &16.0_MK*fac2*(wf(ilevel,ipatch)%fld(2,i,j+1,k,isub)  &
                &             +wf(ilevel,ipatch)%fld(2,i,j-1,k,isub))-&
                &        fac2*(wf(ilevel,ipatch)%fld(2,i,j+2,k,isub)  &
                &             +wf(ilevel,ipatch)%fld(2,i,j-2,k,isub))+&
                &16.0_MK*fac3*(wf(ilevel,ipatch)%fld(2,i,j,k+1,isub)  &
                &             +wf(ilevel,ipatch)%fld(2,i,j,k-1,isub))-&
                &        fac3*(wf(ilevel,ipatch)%fld(2,i,j,k+2,isub)  &
                &             +wf(ilevel,ipatch)%fld(2,i,j,k-2,isub))-&
                &30.0_MK*(fac1+fac2+fac3)*wf(ilevel,ipatch)%fld(2,i,j,k,isub)) +&
                & fac4*(wf(ilevel,ipatch)%fld(1,i+1,j,k,isub) &
                &      *uf(ilevel,ipatch)%fld(2,i+1,j,k,isub)-&
                &       wf(ilevel,ipatch)%fld(1,i-1,j,k,isub) &
                &      *uf(ilevel,ipatch)%fld(2,i-1,j,k,isub))+&
                & fac5*(wf(ilevel,ipatch)%fld(2,i,j+1,k,isub) &
                &      *uf(ilevel,ipatch)%fld(2,i,j+1,k,isub)-&
                &       wf(ilevel,ipatch)%fld(2,i,j-1,k,isub) &
                &      *uf(ilevel,ipatch)%fld(2,i,j-1,k,isub))+&
                & fac6*(wf(ilevel,ipatch)%fld(3,i,j,k+1,isub) &
                &      *uf(ilevel,ipatch)%fld(2,i,j,k+1,isub)-&
                &       wf(ilevel,ipatch)%fld(3,i,j,k-1,isub) &
                &      *uf(ilevel,ipatch)%fld(2,i,j,k-1,isub))-&
                & fac7*(wf(ilevel,ipatch)%fld(1,i+2,j,k,isub) &
                &      *uf(ilevel,ipatch)%fld(2,i+2,j,k,isub)-&
                &       wf(ilevel,ipatch)%fld(1,i-2,j,k,isub) &
                &      *uf(ilevel,ipatch)%fld(2,i-2,j,k,isub))-&
                & fac8*(wf(ilevel,ipatch)%fld(2,i,j+2,k,isub) &
                &      *uf(ilevel,ipatch)%fld(2,i,j+2,k,isub)-&
                &       wf(ilevel,ipatch)%fld(2,i,j-2,k,isub) &
                &      *uf(ilevel,ipatch)%fld(2,i,j-2,k,isub))-&
                & fac9*(wf(ilevel,ipatch)%fld(3,i,j,k+2,isub) &
                &      *uf(ilevel,ipatch)%fld(2,i,j,k+2,isub)-&
                &       wf(ilevel,ipatch)%fld(3,i,j,k-2,isub) &
                &      *uf(ilevel,ipatch)%fld(2,i,j,k-2,isub))

              dwf(ilevel,ipatch)%fld(3,i,j,k,isub) = (&
                &16.0_MK*fac1*(wf(ilevel,ipatch)%fld(3,i+1,j,k,isub)  &
                &             +wf(ilevel,ipatch)%fld(3,i-1,j,k,isub))-&
                &        fac1*(wf(ilevel,ipatch)%fld(3,i+2,j,k,isub)  &
                &             +wf(ilevel,ipatch)%fld(3,i-2,j,k,isub))+&
                &16.0_MK*fac2*(wf(ilevel,ipatch)%fld(3,i,j+1,k,isub)  &
                &             +wf(ilevel,ipatch)%fld(3,i,j-1,k,isub))-&
                &        fac2*(wf(ilevel,ipatch)%fld(3,i,j+2,k,isub)  &
                &             +wf(ilevel,ipatch)%fld(3,i,j-2,k,isub))+&
                &16.0_MK*fac3*(wf(ilevel,ipatch)%fld(3,i,j,k+1,isub)  &
                &             +wf(ilevel,ipatch)%fld(3,i,j,k-1,isub))-&
                &        fac3*(wf(ilevel,ipatch)%fld(3,i,j,k+2,isub)  &
                &             +wf(ilevel,ipatch)%fld(3,i,j,k-2,isub))-&
                &30.0_MK*(fac1+fac2+fac3)*wf(ilevel,ipatch)%fld(3,i,j,k,isub)) +&
                & fac4*(wf(ilevel,ipatch)%fld(1,i+1,j,k,isub) &
                &      *uf(ilevel,ipatch)%fld(3,i+1,j,k,isub)-&
                &       wf(ilevel,ipatch)%fld(1,i-1,j,k,isub) &
                &      *uf(ilevel,ipatch)%fld(3,i-1,j,k,isub))+&
                & fac5*(wf(ilevel,ipatch)%fld(2,i,j+1,k,isub) &
                &      *uf(ilevel,ipatch)%fld(3,i,j+1,k,isub)-&
                &       wf(ilevel,ipatch)%fld(2,i,j-1,k,isub) &
                &      *uf(ilevel,ipatch)%fld(3,i,j-1,k,isub))+&
                & fac6*(wf(ilevel,ipatch)%fld(3,i,j,k+1,isub) &
                &      *uf(ilevel,ipatch)%fld(3,i,j,k+1,isub)-&
                &       wf(ilevel,ipatch)%fld(3,i,j,k-1,isub) &
                &      *uf(ilevel,ipatch)%fld(3,i,j,k-1,isub))-&
                & fac7*(wf(ilevel,ipatch)%fld(1,i+2,j,k,isub) &
                &      *uf(ilevel,ipatch)%fld(3,i+2,j,k,isub)-&
                &       wf(ilevel,ipatch)%fld(1,i-2,j,k,isub) &
                &      *uf(ilevel,ipatch)%fld(3,i-2,j,k,isub))-&
                & fac8*(wf(ilevel,ipatch)%fld(2,i,j+2,k,isub) &
                &      *uf(ilevel,ipatch)%fld(3,i,j+2,k,isub)-&
                &       wf(ilevel,ipatch)%fld(2,i,j-2,k,isub) &
                &      *uf(ilevel,ipatch)%fld(3,i,j-2,k,isub))-&
                & fac9*(wf(ilevel,ipatch)%fld(3,i,j,k+2,isub) &
                &      *uf(ilevel,ipatch)%fld(3,i,j,k+2,isub)-&
                &       wf(ilevel,ipatch)%fld(3,i,j,k-2,isub) &
                &      *uf(ilevel,ipatch)%fld(3,i,j,k-2,isub))

            ENDDO
          ENDDO
        ENDDO
      ENDDO !isub

    ENDIF
  ENDDO
ENDDO

IF ((penalization .GT. 0) .AND. clearinteriorrhs) THEN
  DO ilevel=1,nlevels
    DO ipatch=1,npatches(ilevel)
      DO isub=1,topos(ilevel)%nsublist
        isubl=topos(ilevel)%isublist(isub)
        DO k=1,ptcset(ilevel,ipatch)%snx(3,isubl)
          DO j=1,ptcset(ilevel,ipatch)%snx(2,isubl)
            DO i=1,ptcset(ilevel,ipatch)%snx(1,isubl)
              dwf(ilevel,ipatch)%fld(1,i,j,k,isub) = &
                & dwf(ilevel,ipatch)%fld(1,i,j,k,isub) * &
                & (1.0_MK-chif(ilevel,ipatch)%fld(i,j,k,isub))
              dwf(ilevel,ipatch)%fld(2,i,j,k,isub) = &
                & dwf(ilevel,ipatch)%fld(2,i,j,k,isub) * &
                & (1.0_MK-chif(ilevel,ipatch)%fld(i,j,k,isub))
              dwf(ilevel,ipatch)%fld(3,i,j,k,isub) = &
                & dwf(ilevel,ipatch)%fld(3,i,j,k,isub) * &
                & (1.0_MK-chif(ilevel,ipatch)%fld(i,j,k,isub))
            ENDDO
          ENDDO
        ENDDO
      ENDDO !isub
    ENDDO
  ENDDO
ENDIF


!----------------------------------------------------------------------------
! Return
!----------------------------------------------------------------------------
 9999 CONTINUE
CALL substop('naga_vorticity_rhs',t0,info)
RETURN


END SUBROUTINE naga_vorticity_rhs

END MODULE naga_mod_vorticity_rhs
