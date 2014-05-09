!------------------------------------------------------------------------------
! Subroutine : naga_case_testring.f90
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
! This routines sets up the initial vorticity field
! this is for freespace validation of poisson solution by using the
! bump vorticity ring
! (cf. Particle Methods in Bluff Body Aerodynamics, Rasmussen 2011)
!------------------------------------------------------------------------------
MODULE naga_mod_case_testring
IMPLICIT NONE
INTERFACE naga_case_testring
  MODULE PROCEDURE naga_case_testring
END INTERFACE
CONTAINS
SUBROUTINE naga_case_testring(info)
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
REAL(MK) :: centerx,centery,centerz,centerz2
REAL(MK) :: planx,rho,phi
REAL(MK) :: strength1,strength2
REAL(MK) :: strvelrho,strvelz
REAL(MK) :: r0,lrad,sigma
REAL(MK) :: dx,dy,dz
REAL(MK) :: prefactor
!----------------------------------------------------------------------------
! Initialise routine
!----------------------------------------------------------------------------
CALL substart('naga_case_testring',t0,info)
!----------------------------------------------------------------------------
! Determine center of the base patch
!----------------------------------------------------------------------------
centerx = (ptcset(1,1)%max(1) + ptcset(1,1)%min(1))*0.5_MK
centery = (ptcset(1,1)%max(2) + ptcset(1,1)%min(2))*0.5_MK
centerz = (ptcset(1,1)%max(3) + ptcset(1,1)%min(3))*0.5_MK
r0 = torus_radius1
sigma = torus_radius2
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
      DO k=1-gstw(3),ptcset(ilevel,ipatch)%snx(3,isubl)+gstw(3)
        pz = topos(ilevel)%min_subd(3,isubl) + REAL(k-1,MK)*dz-centerz
        DO j=1-gstw(2),ptcset(ilevel,ipatch)%snx(2,isubl)+gstw(2)
          py = topos(ilevel)%min_subd(2,isubl) + REAL(j-1,MK)*dy-centery
          DO i=1-gstw(1),ptcset(ilevel,ipatch)%snx(1,isubl)+gstw(1)
            px = topos(ilevel)%min_subd(1,isubl) + REAL(i-1,MK)*dx-centerx
            phi = ATAN2(py,px)
            rho = SQRT(px*px+py*py)
            r0 = 1.0_MK
            lrad = 10.0_MK
            prefactor = 1.0E5_MK
            prefactor = 1.0E0_MK
            IF (sqrt((rho-r0)**2+pz**2) .LT. r0-1.0E-12_MK .AND. &
               & sqrt((rho-r0)**2+pz**2) .GT. 1.0E-12_MK) THEN
              strength1 = prefactor * 1.0_MK * exp(-2.0_MK * lrad / (r0) / &
                & (1.0_MK - ((rho ** 2 - 2 * rho * r0 + r0 ** 2 + pz ** 2) / r0 ** 2)))
            strength2 = &
              & prefactor * -(-16.0_MK * lrad / (r0 ** 5) / (1.0_MK - ((rho ** 2 - 2.0_MK * &
              & rho * r0 + r0 ** 2 + pz ** 2) / r0 ** 2)) ** 3 * (pz ** 2) * exp(-2.0_MK * &
              & lrad / (r0) / (1.0_MK - ((rho ** 2 - 2.0_MK * rho * r0 + r0 ** 2 + pz ** &
              & 2) / r0 ** 2))) - 4.0_MK * lrad / (r0 ** 3) / (1.0_MK - ((rho ** 2 - &
              & 2.0_MK * rho * r0 + r0 ** 2 + pz ** 2) / r0 ** 2)) ** 2 * exp(-2.0_MK * &
              & lrad / (r0) / (1.0_MK - ((rho ** 2 - 2.0_MK * rho * r0 + r0 ** 2 + pz ** &
              & 2) / r0 ** 2))) + 16.0_MK * lrad ** 2 / (r0 ** 6) / (1.0_MK - ((rho ** 2 - &
              & 2.0_MK * rho * r0 + r0 ** 2 + pz ** 2) / r0 ** 2)) ** 4 * (pz ** 2) * &
              & exp(-2.0_MK * lrad / (r0) / (1.0_MK - ((rho ** 2 - 2.0_MK * rho * r0 + r0 &
              & ** 2 + pz ** 2) / r0 ** 2))) - 1.0_MK / (rho ** 2) * (1.0_MK * exp(-2.0_MK &
              & * lrad / (r0) / (1.0_MK - ((rho ** 2 - 2.0_MK * rho * r0 + r0 ** 2 + pz ** &
              & 2) / r0 ** 2))) - 2.0_MK * (rho) * lrad / (r0 ** 3) / (1.0_MK - ((rho ** 2 &
              & - 2.0_MK * rho * r0 + r0 ** 2 + pz ** 2) / r0 ** 2)) ** 2 * (2.0_MK * rho &
              & - 2.0_MK * r0) * exp(-2.0_MK * lrad / (r0) / (1.0_MK - ((rho ** 2 - 2.0_MK &
              & * rho * r0 + r0 ** 2 + pz ** 2) / r0 ** 2)))) + 1.0_MK / (rho) * (-4.0_MK &
              & * lrad / (r0 ** 3) / (1.0_MK - ((rho ** 2 - 2.0_MK * rho * r0 + r0 ** 2 + &
              & pz ** 2) / r0 ** 2)) ** 2 * (2.0_MK * rho - 2.0_MK * r0) * exp(-2.0_MK * &
              & lrad / (r0) / (1.0_MK - ((rho ** 2 - 2.0_MK * rho * r0 + r0 ** 2 + pz ** &
              & 2) / r0 ** 2))) - 4.0_MK * (rho) * lrad / (r0 ** 5) / (1.0_MK - ((rho ** 2 &
              & - 2.0_MK * rho * r0 + r0 ** 2 + pz ** 2) / r0 ** 2)) ** 3 * ((2.0_MK * rho &
              & - 2.0_MK * r0) ** 2) * exp(-2.0_MK * lrad / (r0) / (1.0_MK - ((rho ** 2 - &
              & 2.0_MK * rho * r0 + r0 ** 2 + pz ** 2) / r0 ** 2))) - 4.0_MK * (rho) * &
              & lrad / (r0 ** 3) / (1.0_MK - ((rho ** 2 - 2.0_MK * rho * r0 + r0 ** 2 + pz &
              & ** 2) / r0 ** 2)) ** 2 * exp(-2.0_MK * lrad / (r0) / (1.0_MK - ((rho ** 2 &
              & - 2.0_MK * rho * r0 + r0 ** 2 + pz ** 2) / r0 ** 2))) + 4.0_MK * (rho) * &
              & lrad ** 2 / (r0 ** 6) / (1.0_MK - ((rho ** 2 - 2.0_MK * rho * r0 + r0 ** 2 &
              & + pz ** 2) / r0 ** 2)) ** 4 * ((2.0_MK * rho - 2.0_MK * r0) ** 2) * &
              & exp(-2.0_MK * lrad / (r0) / (1.0_MK - ((rho ** 2 - 2.0_MK * rho * r0 + r0 &
              & ** 2 + pz ** 2) / r0 ** 2)))))
            strvelrho = &
              & prefactor * 4.0_MK * lrad / (r0 ** 3) / (1.0_MK - ((rho ** 2 - 2.0_MK * &
              & rho * r0 + r0 ** 2 + pz ** 2) / r0 ** 2)) ** 2 * (pz) * exp(-2.0_MK * lrad &
              & / (r0) / (1.0_MK - ((rho ** 2 - 2.0_MK * rho * r0 + r0 ** 2 + pz ** 2) / &
              & r0 ** 2)))
            strvelz = &
              & prefactor * 1.0_MK / (rho) * (1.0_MK * exp(-2.0_MK * lrad / (r0) / (1.0_MK &
              & - ((rho ** 2 - 2.0_MK * rho * r0 + r0 ** 2 + pz ** 2) / r0 ** 2))) - &
              & 2.0_MK * (rho) * lrad / (r0 ** 3) / (1.0_MK - ((rho ** 2 - 2.0_MK * rho * &
              & r0 + r0 ** 2 + pz ** 2) / r0 ** 2)) ** 2 * (2.0_MK * rho - 2.0_MK * r0) * &
              & exp(-2.0_MK * lrad / (r0) / (1.0_MK - ((rho ** 2 - 2.0_MK * rho * r0 + r0 &
              & ** 2 + pz ** 2) / r0 ** 2))))
            ELSE
              strength1 = 0.0_MK
              strength2 = 0.0_MK
              strvelrho = 0.0_MK
              strvelz = 0.0_MK
            ENDIF
            IF (velocityscheme .EQ. 0) THEN
              psif(ilevel,ipatch)%fld(1,i,j,k,isub) = -SIN(phi)*strength1!comment if K
              psif(ilevel,ipatch)%fld(2,i,j,k,isub) = COS(phi)*strength1!comment if K
              psif(ilevel,ipatch)%fld(3,i,j,k,isub) = 0.0_MK!comment if K
            ELSE!comment if K
              psif(ilevel,ipatch)%fld(1,i,j,k,isub) = strvelrho*px/rho
              psif(ilevel,ipatch)%fld(2,i,j,k,isub) = strvelrho*py/rho
              psif(ilevel,ipatch)%fld(3,i,j,k,isub) = strvelz
            ENDIF
            wf(ilevel,ipatch)%fld(1,i,j,k,isub) = -SIN(phi)*strength2
            wf(ilevel,ipatch)%fld(2,i,j,k,isub) = COS(phi)*strength2
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
CALL substop('naga_case_testring',t0,info)
RETURN
END SUBROUTINE naga_case_testring
END MODULE naga_mod_case_testring
