!------------------------------------------------------------------------------
! Subroutine : naga_case_testreprojection.f90
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
! This routine adds to the vorticity field (which must have been initialised)
! the gradient of a scalar field. The gradient is compact rho < 1 and rho=0 at
! the domain centre.
! The unperturbed vorticity field is saved in psif
!------------------------------------------------------------------------------
MODULE naga_mod_case_testreprojection
IMPLICIT NONE
INTERFACE naga_case_testreprojection
  MODULE PROCEDURE naga_case_testreprojection
END INTERFACE
CONTAINS
SUBROUTINE naga_case_testreprojection(info)
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
CALL substart('naga_case_testreprojection',t0,info)
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
      DO k=1-gstw(3),ptcset(ilevel,ipatch)%snx(3,isubl)+gstw(3)
        pz = topos(ilevel)%min_subd(3,isubl) + REAL(k-1,MK)*dz - centerz
        DO j=1-gstw(2),ptcset(ilevel,ipatch)%snx(2,isubl)+gstw(2)
          py = topos(ilevel)%min_subd(2,isubl) + REAL(j-1,MK)*dy - centery
          DO i=1-gstw(1),ptcset(ilevel,ipatch)%snx(1,isubl)+gstw(1)
            px = topos(ilevel)%min_subd(1,isubl) + REAL(i-1,MK)*dx - centerx
            psif(ilevel,ipatch)%fld(1,i,j,k,isub) = &
              & wf(ilevel,ipatch)%fld(1,i,j,k,isub)
            psif(ilevel,ipatch)%fld(2,i,j,k,isub) = &
              & wf(ilevel,ipatch)%fld(2,i,j,k,isub)
            psif(ilevel,ipatch)%fld(3,i,j,k,isub) = &
              & wf(ilevel,ipatch)%fld(3,i,j,k,isub)
            IF (px**2+py**2+pz**2 .LT. 1.0_MK) THEN
              wf(ilevel,ipatch)%fld(1,i,j,k,isub) = &
                & wf(ilevel,ipatch)%fld(1,i,j,k,isub) - & !@ changed + to -
                & 1.0E4_MK*(-20.0_MK / (1.0_MK - px ** 2 - py ** 2 - pz ** 2) ** 2 * px &
                & * exp(-10.0_MK / (1.0_MK - px ** 2 - py ** 2 - pz ** 2)))
              wf(ilevel,ipatch)%fld(2,i,j,k,isub) = &
                & wf(ilevel,ipatch)%fld(2,i,j,k,isub) - & !@ changed + to -
                & 1.0E4_MK*(-20.0_MK / (1.0_MK - px ** 2 - py ** 2 - pz ** 2) ** 2 * py &
                & * exp(-10.0_MK / (1.0_MK - px ** 2 - py ** 2 - pz ** 2)))
              wf(ilevel,ipatch)%fld(3,i,j,k,isub) = &
                & wf(ilevel,ipatch)%fld(3,i,j,k,isub) - & !@ changed + to -
                & 1.0E4_MK*(-20.0_MK / (1.0_MK - px ** 2 - py ** 2 - pz ** 2) ** 2 * pz &
                & * exp(-10.0_MK / (1.0_MK - px ** 2 - py ** 2 - pz ** 2)))
              !!!#
              !!psif(ilevel,ipatch)%fld(1,i,j,k,isub) = &
                !!& -exp(-10.0_MK / (1.0_MK - px ** 2 - py ** 2 - pz ** 2))
              !!psif(ilevel,ipatch)%fld(2,i,j,k,isub) = &
                !!& -exp(-10.0_MK / (1.0_MK - px ** 2 - py ** 2 - pz ** 2))
              !!psif(ilevel,ipatch)%fld(3,i,j,k,isub) = &
                !!& -exp(-10.0_MK / (1.0_MK - px ** 2 - py ** 2 - pz ** 2))
              !!!#
            !!ELSE
              !!!#
              !!psif(ilevel,ipatch)%fld(1,i,j,k,isub) = 0.0_MK
              !!psif(ilevel,ipatch)%fld(2,i,j,k,isub) = 0.0_MK
              !!psif(ilevel,ipatch)%fld(3,i,j,k,isub) = 0.0_MK
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
CALL substop('naga_case_testreprojection',t0,info)
RETURN
END SUBROUTINE naga_case_testreprojection
END MODULE naga_mod_case_testreprojection
