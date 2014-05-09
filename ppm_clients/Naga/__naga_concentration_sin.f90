!------------------------------------------------------------------------------
! Subroutine : naga_concentration_sin.f90
!------------------------------------------------------------------------------
! _ __
! / |/ /__ ____ ____ _
! / / _ `/ _ `/ _ `/
! /_/|_/\_,_/\_, /\_,_/
! /___/
! A PPM Vortex-In-Cell client - 2011 (c)
! Naga is distributed under the GNU Lesser General Public License version 3
! Johannes Tophøj Rasmussen, johannes.tophoj@gmail.com
! Technical University of Denmark
! Department of Mechanical Engineering
!
!
! This routines sets up an initial scalar field (concentration, heat, etc.).
! It is sinusoidal in the z-direction
!------------------------------------------------------------------------------
MODULE naga_mod_concentration_sin
IMPLICIT NONE
INTERFACE naga_concentration_sin
  MODULE PROCEDURE naga_concentration_sin
END INTERFACE
CONTAINS
SUBROUTINE naga_concentration_sin(info)
USE naga_mod_globals
IMPLICIT NONE
!-------------------------------------------------------------------------
! Arguments
!-------------------------------------------------------------------------
INTEGER, INTENT(OUT) :: info
!-------------------------------------------------------------------------
! Local variables
!-------------------------------------------------------------------------
REAL(MK),PARAMETER :: PI=ACOS(-1.0_mk)
REAL(MK) :: t0
REAL(MK) :: px, py, pz
INTEGER :: i,j,k
INTEGER :: ilevel,ipatch,isubl,isub
REAL(MK) :: centerx,centery,centerz
REAL(MK) :: dx,dy,dz
!----------------------------------------------------------------------------
! Initialise routine
!----------------------------------------------------------------------------
CALL substart('naga_concentration_sin',t0,info)
!----------------------------------------------------------------------------
! Determine center of the base patch
!----------------------------------------------------------------------------
centerx = (ptcset(1,1)%max(1) + ptcset(1,1)%min(1))*0.5_MK
centery = (ptcset(1,1)%max(2) + ptcset(1,1)%min(2))*0.5_MK
centerz = (ptcset(1,1)%max(3) + ptcset(1,1)%min(3))*0.9_MK
!----------------------------------------------------------------------------
! Loop through the sub domain points
!----------------------------------------------------------------------------
DO ilevel=1,nlevels
  DO ipatch=1,npatches(ilevel)
    DO isub=1,topos(ilevel)%nsublist
      dx = ptcset(ilevel,ipatch)%dx(1)
      dy = ptcset(ilevel,ipatch)%dx(2)
      dz = ptcset(ilevel,ipatch)%dx(3)
      isubl=topos(ilevel)%isublist(isub)
      DO k=1,ptcset(ilevel,ipatch)%snx(3,isubl)
        pz = topos(ilevel)%min_subd(3,isubl) + REAL(k-1,MK)*dx-centerz
        DO j=1,ptcset(ilevel,ipatch)%snx(2,isubl)
          py = topos(ilevel)%min_subd(2,isubl) + REAL(j-1,MK)*dy-centery
          DO i=1,ptcset(ilevel,ipatch)%snx(1,isubl)
            px = topos(ilevel)%min_subd(1,isubl) + REAL(i-1,MK)*dz-centerx
            cf(ilevel,ipatch)%fld(i,j,k,isub) = SIN(2*PI*pz)
          END DO !i
        END DO !j
      END DO !k
    END DO !isub
  END DO !ipatch
END DO !ilevel
!----------------------------------------------------------------------------
! Return
!----------------------------------------------------------------------------
 9999 CONTINUE
CALL substop('naga_concentration_sin',t0,info)
RETURN
END SUBROUTINE naga_concentration_sin
END MODULE naga_mod_concentration_sin
