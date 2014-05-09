!------------------------------------------------------------------------------
! Subroutine : naga_case_testperiodic.f90
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
! to test the periodic Green's function. The Poisson initialisation must be
! without derivatives and some additional output (comparison/error) must be
! added e.g. in Naga.f90. Also rk1 and one timestep should be used,
! no reprojection etc etc
!------------------------------------------------------------------------------
MODULE naga_mod_case_testperiodic
IMPLICIT NONE
INTERFACE naga_case_testperiodic
  MODULE PROCEDURE naga_case_testperiodic
END INTERFACE
CONTAINS
SUBROUTINE naga_case_testperiodic(info)
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
REAL(MK) :: normx,normy,normz
REAL(MK) :: dx,dy,dz
REAL(MK) :: tmp,k1,k2 !@tmp
!----------------------------------------------------------------------------
! Initialise routine
!----------------------------------------------------------------------------
CALL substart('naga_case_testperiodic',t0,info)
!----------------------------------------------------------------------------
! Determine normalization according to the base patch (level 1)
!----------------------------------------------------------------------------
normx = 2.0_MK*PI/(ptcset(1,1)%max(1) - ptcset(1,1)%min(1))*1.0_MK
normy = 2.0_MK*PI/(ptcset(1,1)%max(2) - ptcset(1,1)%min(2))*1.0_MK
normz = 2.0_MK*PI/(ptcset(1,1)%max(3) - ptcset(1,1)%min(3))*1.0_MK
!----------------------------------------------------------------------------
! Loop through the sub domain points
!----------------------------------------------------------------------------
DO ilevel=1,nlevels
  DO ipatch=1,npatches(ilevel)
    DO isub=1,topos(ilevel)%nsublist
      dx = ptcset(ilevel,ipatch)%dx(1)
      dy = ptcset(ilevel,ipatch)%dx(2)
      dz = ptcset(ilevel,ipatch)%dx(3)
      psif(ilevel,ipatch)%fld = 0.0_MK
      wf(ilevel,ipatch)%fld = 0.0_MK
      isubl=topos(ilevel)%isublist(isub)
      DO k=1-gstw(3),ptcset(ilevel,ipatch)%snx(3,isubl)+gstw(3)
        pz = topos(ilevel)%min_subd(3,isubl) + REAL(k-1,MK)*dz
        DO j=1-gstw(2),ptcset(ilevel,ipatch)%snx(2,isubl)+gstw(2)
          py = topos(ilevel)%min_subd(2,isubl) + REAL(j-1,MK)*dy
          DO i=1-gstw(1),ptcset(ilevel,ipatch)%snx(1,isubl)+gstw(1)
            px = topos(ilevel)%min_subd(1,isubl) + REAL(i-1,MK)*dx
            !----------Some old test that show spectral convergence. Only for that
            !vorticity
            wf(ilevel,ipatch)%fld(1,i,j,k,isub) = -exp(sin(normx * px)) *&
            & sin(normy * py) * sin(normz * pz) * (sin(normx * px) * normx ** 2&
            & - cos(normx * px) ** 2 * normx ** 2 + normy ** 2 + normz ** 2)
            wf(ilevel,ipatch)%fld(2,i,j,k,isub) = -exp(sin(normx * px)) *&
            & sin(normy * py) * sin(normz * pz) * (sin(normx * px) * normx ** 2&
            & - cos(normx * px) ** 2 * normx ** 2 + normy ** 2 + normz ** 2)
            wf(ilevel,ipatch)%fld(3,i,j,k,isub) = -exp(sin(normx * px)) *&
            & sin(normy * py) * sin(normz * pz) * (sin(normx * px) * normx ** 2&
            & - cos(normx * px) ** 2 * normx ** 2 + normy ** 2 + normz ** 2)
            !stream function
            psif(ilevel,ipatch)%fld(1,i,j,k,isub) = -exp(sin(normx * px)) &
            &* sin(normy * py) * sin(normz * pz)
            psif(ilevel,ipatch)%fld(2,i,j,k,isub) = -exp(sin(normx * px)) &
            &* sin(normy * py) * sin(normz * pz)
            psif(ilevel,ipatch)%fld(3,i,j,k,isub) = -exp(sin(normx * px)) &
            &* sin(normy * py) * sin(normz * pz)
            IF (VELOCITYSCHEME .NE. 0) THEN
              !velocity
              psif(ilevel,ipatch)%fld(1,i,j,k,isub) = &
                & exp(sin(px * normx)) * (-cos(py * normy) * normy * sin(pz * normz) &
                & + sin(py * normy) * cos(pz * normz) * normz)
              psif(ilevel,ipatch)%fld(2,i,j,k,isub) = &
                & exp(sin(px * normx)) * sin(py * normy) * (-cos(pz * normz) * normz &
                & + cos(px * normx) * normx * sin(pz * normz))
              psif(ilevel,ipatch)%fld(3,i,j,k,isub) = &
                & -exp(sin(px * normx)) * sin(pz * normz) * (cos(px * normx) * normx &
                & * sin(py * normy) - cos(py * normy) * normy)
            ENDIF
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
CALL substop('naga_case_testperiodic',t0,info)
RETURN
END SUBROUTINE naga_case_testperiodic
END MODULE naga_mod_case_testperiodic
