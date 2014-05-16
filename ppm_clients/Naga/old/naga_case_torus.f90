!------------------------------------------------------------------------------
! Subroutine :  naga_case_torus.f90
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
! This routines sets up the initial vorticity field of a vortex ring with 
! sinusoidal core. Parameters torus_radius1 (size), torus_radius2 (thickness)
! The ring is centred in the domain and is initialised with periodic images
!------------------------------------------------------------------------------

MODULE naga_mod_case_torus

IMPLICIT NONE

INTERFACE naga_case_torus
  MODULE PROCEDURE naga_case_torus
END INTERFACE

CONTAINS

SUBROUTINE naga_case_torus(info)

USE naga_mod_globals

IMPLICIT NONE
!-------------------------------------------------------------------------
! Arguments
!-------------------------------------------------------------------------
INTEGER, INTENT(OUT) :: info


!-------------------------------------------------------------------------
! Local variables
!-------------------------------------------------------------------------
REAL(MK),PARAMETER   :: PI=ACOS(-1.0_mk)
REAL(MK)             :: t0
REAL(MK)             :: px, py, pz
INTEGER              :: i,j,k
INTEGER              :: ilevel,ipatch,isubl,isub
REAL(MK)             :: centerx,centery,centerz
REAL(MK)             :: planx,rho,phi
REAL(MK)             :: strength
REAL(MK)             :: dx,dy,dz


!----------------------------------------------------------------------------
! Initialise routine
!----------------------------------------------------------------------------
CALL substart('naga_case_torus',t0,info)


!----------------------------------------------------------------------------
! Determine center of the base patch
!----------------------------------------------------------------------------
centerx = (ptcset(1,1)%max(1) + ptcset(1,1)%min(1))*0.5_MK
centery = (ptcset(1,1)%max(2) + ptcset(1,1)%min(2))*0.5_MK
centerz = (ptcset(1,1)%max(3) + ptcset(1,1)%min(3))*0.5_MK


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

            phi   = ATAN2(pz,px)
            planx = SQRT(px*px+pz*pz)-torus_radius1
            rho   = SQRT(planx*planx+py*py)

            strength = exp(-rho*rho/(2*torus_radius2**2)) &
                     & /SQRT(2*PI*torus_radius2**2)

            wf(ilevel,ipatch)%fld(1,i,j,k,isub) = SIN(phi)*strength
            wf(ilevel,ipatch)%fld(2,i,j,k,isub) = COS(phi)*strength
            wf(ilevel,ipatch)%fld(3,i,j,k,isub) = 0.0_MK

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
CALL substop('naga_case_torus',t0,info)


RETURN

END SUBROUTINE naga_case_torus

END MODULE naga_mod_case_torus

