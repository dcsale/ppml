!------------------------------------------------------------------------------
! Subroutine :  naga_stl_init.f90
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
! This routines initialises the solid mask from triangulated surface STL files.
! Calculates intersections between (a line from the point to a 
! reference point) and (the triangles). Calculates barycentric coordinates and
! determines whether the line intersects the triangle. If so the direction of
! interesection is stored and after checking all triangles it can be determined
! if the point is interior or exterior.
!------------------------------------------------------------------------------

SUBROUTINE naga_stl_init(info)

  USE naga_mod_globals
  USE naga_mod_stepfunction
  USE naga_mod_say

  IMPLICIT NONE

  INTEGER                  :: info


  INTEGER                  :: i,j,k,tr,isub,isubl
  INTEGER                  :: inout,intersections
  REAL(Mk), DIMENSION(3)   :: vecT,vecP,vecS,vece
  REAL(Mk), DIMENSION(3)   :: vecX
  REAL(Mk)                 :: Pdotu,Pdotv,a,b,c
  REAL(Mk)                 :: epsilon
  INTEGER                  :: ibmin,ibmax,jbmin,jbmax,kbmin,kbmax
  INTEGER                  :: imin,imax,jmin,jmax,kmin,kmax
  REAL(Mk)                 :: minimumdist,diag
  INTEGER                  :: npenextra
  REAL(Mk)                 :: udotx,vdotx,dist
  INTEGER                  :: ilevel,ipatch
  REAL(Mk)                 :: dx,dy,dz


  info = 0


  DO ilevel=1,nlevels
    DO ipatch=1,npatches(ilevel)
      dx = ptcset(ilevel,ipatch)%dx(1)
      dy = ptcset(ilevel,ipatch)%dx(2)
      dz = ptcset(ilevel,ipatch)%dx(3)
      epsilon = 1.0_MK / (sqrt(dx**2 + dy**2 + dz**2))
      chif(ilevel,ipatch)%fld = 0.0_MK
      diag = (ptcset(ilevel,ipatch)%max(1) - ptcset(ilevel,ipatch)%min(1))**2+ &
           & (ptcset(ilevel,ipatch)%max(2) - ptcset(ilevel,ipatch)%min(2))**2+ &
           & (ptcset(ilevel,ipatch)%max(3) - ptcset(ilevel,ipatch)%min(3))**2
      DO isub=1,topos(ilevel)%nsublist
        isubl=topos(ilevel)%isublist(isub)

        !---------------------------------------------------------------------------
        ! determine boundaries for the cells to be evaluated
        ! the boundaries are expanded by the number of cells equal to 50% of 
        ! step1_interval cell diagonals
        ! ijkbmin will only exceed ijkbmax when a part of the CV is in the subdomain
        ! so that DO loops only initiate when CVs are present in the subdomain:
        !---------------------------------------------------------------------------
        npenextra = ceiling((step1_interval*0.5_Mk + step1_offset)/epsilon)
        ibmin = NINT((stl(1)%bndminx-ptcset(1,1)%min(1))/dx)+1 - &
              & (ptcset(ilevel,ipatch)%sistr(1,isubl)-1)-npenextra
        jbmin = NINT((stl(1)%bndminy-ptcset(1,1)%min(2))/dy)+1 - &
              & (ptcset(ilevel,ipatch)%sistr(2,isubl)-1)-npenextra
        kbmin = NINT((stl(1)%bndminz-ptcset(1,1)%min(3))/dz)+1 - &
              & (ptcset(ilevel,ipatch)%sistr(3,isubl)-1)-npenextra
        ibmax = NINT((stl(1)%bndmaxx-ptcset(1,1)%min(1))/dx)+1 - &
              & (ptcset(ilevel,ipatch)%sistr(1,isubl)-1)+npenextra
        jbmax = NINT((stl(1)%bndmaxy-ptcset(1,1)%min(2))/dy)+1 - &
              & (ptcset(ilevel,ipatch)%sistr(2,isubl)-1)+npenextra
        kbmax = NINT((stl(1)%bndmaxz-ptcset(1,1)%min(3))/dz)+1 - &
              & (ptcset(ilevel,ipatch)%sistr(3,isubl)-1)+npenextra
        imin=max(ibmin,1-gstw(1)) 
        jmin=max(jbmin,1-gstw(2))
        kmin=max(kbmin,1-gstw(3))
        imax=min(ibmax,ptcset(ilevel,ipatch)%snx(1,isubl)+gstw(1))
        jmax=min(jbmax,ptcset(ilevel,ipatch)%snx(2,isubl)+gstw(2))
        kmax=min(kbmax,ptcset(ilevel,ipatch)%snx(3,isubl)+gstw(3))

        IF (stlopt_inout_direction .EQ. 1) THEN
#define __z    1
#define __y    3
#define __x    2
#include "naga_stl_inout.f90"
        ELSE IF (stlopt_inout_direction .EQ. 2) THEN
#define __z    2
#define __y    1
#define __x    3
#include "naga_stl_inout.f90"
        ELSE IF (stlopt_inout_direction .EQ. 3) THEN
#define __z    3
#define __y    2
#define __x    1
#include "naga_stl_inout.f90"
        ENDIF
#undef __y
#undef __x
#undef __z

      END DO !isub
    END DO !ipatch
  END DO !ilevel

END SUBROUTINE naga_stl_init

