!------------------------------------------------------------------------------
! Subroutine : naga_structure_sphere.f90
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
! This routines sets up the mask corresponding to a sphere centred around 0,0,0
! and sets ubar 0
!------------------------------------------------------------------------------
MODULE naga_mod_structure_sphere
IMPLICIT NONE
INTERFACE naga_structure_sphere
  MODULE PROCEDURE naga_structure_sphere
END INTERFACE
CONTAINS
SUBROUTINE naga_structure_sphere(info)
USE naga_mod_globals
USE naga_mod_stepfunction
IMPLICIT NONE
!-------------------------------------------------------------------------
! Arguments
!-------------------------------------------------------------------------
INTEGER, INTENT(OUT) :: info
!-------------------------------------------------------------------------
! Local variables
!-------------------------------------------------------------------------
REAL(MK) :: t0
REAL(MK) :: px, py, pz
INTEGER :: i,j,k
INTEGER :: ilevel,ipatch,isubl,isub
REAL(MK) :: centerx,centery,centerz
REAL(MK) :: rho
REAL(MK) :: radius1,radius3
REAL(MK) :: dx,dy,dz
REAL(MK) :: norm
!----------------------------------------------------------------------------
! Initialise routine
!----------------------------------------------------------------------------
CALL substart('naga_structure_sphere',t0,info)
!@ get from input file:
!!radius1 = -0.15_MK*sphere_radius
!!radius3 = 0.15_MK*sphere_radius
!----------------------------------------------------------------------------
! Determine centre of the sphere
!----------------------------------------------------------------------------
centerx = 0.0_MK
centery = 0.0_MK
centerz = 0.0_MK
!----------------------------------------------------------------------------
! Loop through the sub domain points
!----------------------------------------------------------------------------
DO ilevel=1,nlevels
  DO ipatch=1,npatches(ilevel)
    DO isub=1,topos(ilevel)%nsublist
      dx = ptcset(ilevel,ipatch)%dx(1)
      dy = ptcset(ilevel,ipatch)%dx(2)
      dz = ptcset(ilevel,ipatch)%dx(3)
      norm = 1.0_MK/SQRT(dx**2+dy**2+dz**2)
      isubl=topos(ilevel)%isublist(isub)
      DO k=1-gstw(1),ptcset(ilevel,ipatch)%snx(3,isubl)+gstw(1)
        pz = topos(ilevel)%min_subd(3,isubl) + REAL(k-1,MK)*dz-centerz
        DO j=1-gstw(2),ptcset(ilevel,ipatch)%snx(2,isubl)+gstw(2)
          py = topos(ilevel)%min_subd(2,isubl) + REAL(j-1,MK)*dy-centery
          DO i=1-gstw(3),ptcset(ilevel,ipatch)%snx(1,isubl)+gstw(3)
            px = topos(ilevel)%min_subd(1,isubl) + REAL(i-1,MK)*dx-centerx
            rho = (SQRT(px*px+py*py+pz*pz) - &
                  & sphere_radius)*norm
            chif(ilevel,ipatch)%fld(i,j,k,isub) = naga_stepfunc1(rho)
            ubarf(ilevel,ipatch)%fld(1,i,j,k,isub) = 0.0_MK
            ubarf(ilevel,ipatch)%fld(2,i,j,k,isub) = 0.0_MK
            ubarf(ilevel,ipatch)%fld(3,i,j,k,isub) = 0.0_MK
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
CALL substop('naga_structure_sphere',t0,info)
RETURN
END SUBROUTINE naga_structure_sphere
END MODULE naga_mod_structure_sphere
