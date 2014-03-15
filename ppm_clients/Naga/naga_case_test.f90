!------------------------------------------------------------------------------
! Subroutine :  naga_case_test.f90
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
! This routines sets up the initial vorticity field
! corresponding to test. Use the routine to initialise various vorticity fields
!------------------------------------------------------------------------------

MODULE naga_mod_case_test

IMPLICIT NONE

INTERFACE naga_case_test
  MODULE PROCEDURE naga_case_test
END INTERFACE

CONTAINS

SUBROUTINE naga_case_test(info)

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
REAL(MK)             :: normx,normy,normz
REAL(MK)             :: dx,dy,dz
REAL(MK)             :: tmp,k1,k2 !@tmp

!----------------------------------------------------------------------------
! Initialise routine
!----------------------------------------------------------------------------
CALL substart('naga_case_test',t0,info)


!----------------------------------------------------------------------------
! Determine normalization according to the base patch (level 1)
!----------------------------------------------------------------------------
normx = 2.0_MK*PI/(ptcset(1,1)%max(1) - ptcset(1,1)%min(1))*1.0_MK!obs
normy = 2.0_MK*PI/(ptcset(1,1)%max(2) - ptcset(1,1)%min(2))*1.0_MK!obs
normz = 2.0_MK*PI/(ptcset(1,1)%max(3) - ptcset(1,1)%min(3))*1.0_MK!obs
!normx = 2.0_MK/(ptcset(1,1)%max(1) - ptcset(1,1)%min(1)) !for 3rd ord pol
!normy = 2.0_MK/(ptcset(1,1)%max(2) - ptcset(1,1)%min(2))
!normz = 2.0_MK/(ptcset(1,1)%max(3) - ptcset(1,1)%min(3))
normx = 2.0_MK*PI/(ptcset(1,1)%max(1) - ptcset(1,1)%min(1))*1.0_MK!obs
normy = 2.0_MK*PI/(ptcset(1,1)%max(2) - ptcset(1,1)%min(2))*1.0_MK!obs
normz = 2.0_MK*PI/(ptcset(1,1)%max(3) - ptcset(1,1)%min(3))*1.0_MK!obs

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
      wf(ilevel,ipatch)%fld   = 0.0_MK
      isubl=topos(ilevel)%isublist(isub)
      DO k=1-gstw(3),ptcset(ilevel,ipatch)%snx(3,isubl)+gstw(3)
        pz = topos(ilevel)%min_subd(3,isubl) + REAL(k-1,MK)*dz
        DO j=1-gstw(2),ptcset(ilevel,ipatch)%snx(2,isubl)+gstw(2)
          py = topos(ilevel)%min_subd(2,isubl) + REAL(j-1,MK)*dy
          DO i=1-gstw(1),ptcset(ilevel,ipatch)%snx(1,isubl)+gstw(1)
            px = topos(ilevel)%min_subd(1,isubl) + REAL(i-1,MK)*dx
            !test first and then shift to interior of domain
            !or just do convolution

            !--------------------------------------------------------------------
            ! THESE ARE TESTS FOR PERIODIC DOMAINS
            !--------------------------------------------------------------------
            !--------Taylor Green vortices
            !stream function reference
            !psif(ilevel,ipatch)%fld(1,i,j,k,isub) =  sin(py * normy) * cos(pz * normz)
            !psif(ilevel,ipatch)%fld(2,i,j,k,isub) =  sin(pz * normz) * cos(px * normx)
            !psif(ilevel,ipatch)%fld(3,i,j,k,isub) =  sin(px * normx) * cos(py * normy)

            !vorticity
            !wf(ilevel,ipatch)%fld(1,i,j,k,isub) = &
              !& -sin(py * normy) * cos(pz * normz) * normy ** 2 &
              !& -sin(py * normy) * cos(pz * normz) * normz ** 2
            !wf(ilevel,ipatch)%fld(2,i,j,k,isub) = &
              !& -sin(pz * normz) * cos(px * normx) * normz ** 2 &
              !& -sin(pz * normz) * cos(px * normx) * normx ** 2
            !wf(ilevel,ipatch)%fld(3,i,j,k,isub) = &
              !& -sin(px * normx) * cos(py * normy) * normx ** 2 &
              !& -sin(px * normx) * cos(py * normy) * normy ** 2

            !velocity reference
            !psif(ilevel,ipatch)%fld(1,i,j,k,isub) = &
              !& -sin(px * normx) * sin(py * normy) * normy &
              !& -cos(pz * normz) * cos(px * normx) * normz
            !psif(ilevel,ipatch)%fld(2,i,j,k,isub) = &
              !& -sin(py * normy) * sin(pz * normz) * normz &
              !& -cos(px * normx) * cos(py * normy) * normx
            !psif(ilevel,ipatch)%fld(3,i,j,k,isub) = &
              !& -sin(pz * normz) * sin(px * normx) * normx &
              !& -cos(py * normy) * cos(pz * normz) * normy



            !!-------Custom vortex (TG vorticies with a little extra)
            !!MAYBE THIS IS BEST TO TEST DIRECTLY ON STREAM FUNCTION
            !!stream function
            !psif(ilevel,ipatch)%fld(1,i,j,k,isub) = &
            !    & py ** 2 * normy ** 2 * sin(py * normy) * cos(pz * normz)
            !psif(ilevel,ipatch)%fld(2,i,j,k,isub) = &
            !    & pz ** 2 * normz ** 2 * sin(pz * normz) * cos(px * normx)
            !psif(ilevel,ipatch)%fld(3,i,j,k,isub) = &
            !    & px ** 2 * normx ** 2 * sin(px * normx) * cos(py * normy)

            !!vorticity
            !wf(ilevel,ipatch)%fld(1,i,j,k,isub) = &
            !& 2.0_MK * sin(py * normy) * normy ** 2 * cos(pz * normz) &
            !& + 4.0_MK * py * normy ** 3 * cos(py * normy) * cos(pz * normz) &
            !& - py ** 2 * normy ** 4 * sin(py * normy) * cos(pz * normz) &
            !& - py ** 2 * normy ** 2 * sin(py * normy) * cos(pz * normz) * normz ** 2
            !wf(ilevel,ipatch)%fld(2,i,j,k,isub) = &
            !& -pz ** 2 * normz ** 2 * sin(pz * normz) * cos(px * normx) * normx ** 2 &
            !& + 2.0_MK * sin(pz * normz) * normz ** 2 * cos(px * normx) &
            !& + 4.0_MK * pz * normz ** 3 * cos(pz * normz) * cos(px * normx) &
            !& - pz ** 2 * normz ** 4 * sin(pz * normz) * cos(px * normx)
            !wf(ilevel,ipatch)%fld(3,i,j,k,isub) = &
            !& 2.0_MK * sin(px * normx) * normx ** 2 * cos(py * normy) &
            !& + 4.0_MK * px * normx ** 3 * cos(px * normx) * cos(py * normy) &
            !& - px ** 2 * normx ** 4 * sin(px * normx) * cos(py * normy) &
            !& - px ** 2 * normx ** 2 * sin(px * normx) * cos(py * normy) * normy ** 2

            !!velocity
            !psif(ilevel,ipatch)%fld(1,i,j,k,isub) = &
            !& -px ** 2 * normx ** 2 * sin(px * normx) * sin(py * normy) * normy &
            !& - 2.0_MK * pz * normz ** 2 * sin(pz * normz) * cos(px * normx) &
            !& - pz ** 2 * normz ** 3 * cos(pz * normz) * cos(px * normx)
            !psif(ilevel,ipatch)%fld(2,i,j,k,isub) = &
            !& -py ** 2 * normy ** 2 * sin(py * normy) * sin(pz * normz) * normz &
            !& - 2.0_MK * px * normx ** 2 * sin(px * normx) * cos(py * normy) &
            !& - px ** 2 * normx ** 3 * cos(px * normx) * cos(py * normy)
            !psif(ilevel,ipatch)%fld(3,i,j,k,isub) = &
            !& -pz ** 2 * normz ** 2 * sin(pz * normz) * sin(px * normx) * normx &
            !& - 2.0_MK * py * normy ** 2 * sin(py * normy) * cos(pz * normz) &
            !& - py ** 2 * normy ** 3 * cos(py * normy) * cos(pz * normz)

            !-------------Third order polynomials instead of sine
            !psif(ilevel,ipatch)%fld(1,i,j,k,isub) = py * normy * &
            !& (py * normy - 1.0_MK) * (py * normy - 2_MK) * pz * normz * &
            !& (pz * normz - 1.0_MK) * (pz * normz - 2_MK)
            !psif(ilevel,ipatch)%fld(2,i,j,k,isub) = pz * normz * &
            !& (pz * normz - 1.0_MK) * (pz * normz - 2_MK) * px * normx * &
            !& (px * normx - 1.0_MK) * (px * normx - 2_MK)
            !psif(ilevel,ipatch)%fld(3,i,j,k,isub) = px * normx * &
            !& (px * normx - 1.0_MK) * (px * normx - 2_MK) * py * normy * &
            !& (py * normy - 1.0_MK) * (py * normy - 2_MK)

            !wf(ilevel,ipatch)%fld(1,i,j,k,isub) = &
            !& 2.0_MK * normy ** 2 * (py * normy - 2.0_MK) * pz * normz &
            !& * (pz * normz - 1.0_MK) * (pz * normz - 2.0_MK) + 2.0_MK * &
            !& normy ** 2 * (py * normy - 1.0_MK) * pz * normz * (pz * normz &
            !& - 1.0_MK) * (pz * normz - 2.0_MK) + 2.0_MK * py * normy ** 3 * &
            !& pz * normz * (pz * normz - 1.0_MK) * (pz * normz - 2.0_MK) + &
            !& 2.0_MK * py * normy * (py * normy - 1.0_MK) * (py * normy - &
            !& 2.0_MK) * normz ** 2 * (pz * normz - 2.0_MK) + 2.0_MK * py * &
            !& normy * (py * normy - 1.0_MK) * (py * normy - 2.0_MK) * &
            !& normz ** 2 * (pz * normz - 1.0_MK) + 2.0_MK * py * normy * (py &
            !& * normy - 1.0_MK) * (py * normy - 2.0_MK) * pz * normz ** 3
            !wf(ilevel,ipatch)%fld(2,i,j,k,isub) = &
            !& 2.0_MK * pz * normz * (pz * normz - 1.0_MK) * (pz * normz &
            !& - 2.0_MK) * normx ** 2 * (px * normx - 2.0_MK) + 2.0_MK * pz * &
            !& normz * (pz * normz - 1.0_MK) * (pz * normz - 2.0_MK) * &
            !& normx ** 2 * (px * normx - 1.0_MK) + 2.0_MK * pz * normz * (pz &
            !& * normz - 1.0_MK) * (pz * normz - 2.0_MK) * px * normx ** 3 + &
            !& 2.0_MK * normz ** 2 * (pz * normz - 2.0_MK) * px * normx * (px &
            !& * normx - 1.0_MK) * (px * normx - 2.0_MK) + 2.0_MK * normz ** 2 &
            !& * (pz * normz - 1.0_MK) * px * normx * (px * normx - 1.0_MK) * &
            !& (px * normx - 2.0_MK) + 2.0_MK * pz * normz ** 3 * px * normx * &
            !& (px * normx - 1.0_MK) * (px * normx - 2.0_MK)
            !wf(ilevel,ipatch)%fld(3,i,j,k,isub) = &
            !& 2.0_MK * normx ** 2 * (px * normx - 2.0_MK) * py * normy &
            !& * (py * normy - 1.0_MK) * (py * normy - 2.0_MK) + 2.0_MK * &
            !& normx ** 2 * (px * normx - 1.0_MK) * py * normy * (py * normy - &
            !& 1.0_MK) * (py * normy - 2.0_MK) + 2.0_MK * px * normx ** 3 * py &
            !& * normy * (py * normy - 1.0_MK) * (py * normy - 2.0_MK) + 2.0_MK &
            !& * px * normx * (px * normx - 1.0_MK) * (px * normx - 2.0_MK) * &
            !& normy ** 2 * (py * normy - 2.0_MK) + 2.0_MK * px * normx * (px * &
            !& normx - 1.0_MK) * (px * normx - 2.0_MK) * normy ** 2 * (py * &
            !& normy - 1.0_MK) + 2.0_MK * px * normx * (px * normx - 1.0_MK) * &
            !& (px * normx - 2.0_MK) * py * normy ** 3

            !psif(ilevel,ipatch)%fld(1,i,j,k,isub) = px * normx * (px * normx &
            !& - 1.0_MK) * (px * normx - 2.0_MK) * normy * (py * normy - &
            !& 1.0_MK) * (py * normy - 2.0_MK) + px * normx * (px * normx &
            !& - 1.0_MK) * (px * normx - 2.0_MK) * py * normy ** 2 * (py * &
            !& normy - 2.0_MK) + px * normx * (px * normx - 1.0_MK) * (px * &
            !& normx - 2.0_MK) * py * normy ** 2 * (py * normy - 1.0_MK) - &
            !& normz * (pz * normz - 1.0_MK) * (pz * normz - 2.0_MK) * px * &
            !& normx * (px * normx - 1.0_MK) * (px * normx - 2.0_MK) - pz * &
            !& normz ** 2 * (pz * normz - 2.0_MK) * px * normx * (px * normx &
            !& - 1.0_MK) * (px * normx - 2.0_MK) - pz * normz ** 2 * (pz * &
            !& normz - 1.0_MK) * px * normx * (px * normx - 1.0_MK) * (px * &
            !& normx - 2.0_MK)
            !psif(ilevel,ipatch)%fld(2,i,j,k,isub) = py * normy * (py * normy &
            !& - 1.0_MK) * (py * normy - 2.0_MK) * normz * (pz * normz - &
            !& 1.0_MK) * (pz * normz - 2.0_MK) + py * normy * (py * normy &
            !& - 1.0_MK) * (py * normy - 2.0_MK) * pz * normz ** 2 * (pz * &
            !& normz - 2.0_MK) + py * normy * (py * normy - 1.0_MK) * (py * &
            !& normy - 2.0_MK) * pz * normz ** 2 * (pz * normz - 1.0_MK) - &
            !& normx * (px * normx - 1.0_MK) * (px * normx - 2.0_MK) * py * &
            !& normy * (py * normy - 1.0_MK) * (py * normy - 2.0_MK) - px * &
            !& normx ** 2 * (px * normx - 2.0_MK) * py * normy * (py * normy &
            !& - 1.0_MK) * (py * normy - 2.0_MK) - px * normx ** 2 * (px * &
            !& normx - 1.0_MK) * py * normy * (py * normy - 1.0_MK) * (py * &
            !& normy - 2.0_MK)
            !psif(ilevel,ipatch)%fld(3,i,j,k,isub) = pz * normz * (pz * normz &
            !& - 1.0_MK) * (pz * normz - 2.0_MK) * normx * (px * normx - &
            !& 1.0_MK) * (px * normx - 2.0_MK) + pz * normz * (pz * normz - &
            !& 1.0_MK) * (pz * normz - 2.0_MK) * px * normx ** 2 * (px * normx &
            !& - 2.0_MK) + pz * normz * (pz * normz - 1.0_MK) * (pz * normz - &
            !& 2.0_MK) * px * normx ** 2 * (px * normx - 1.0_MK) - normy * (py &
            !& * normy - 1.0_MK) * (py * normy - 2.0_MK) * pz * normz * (pz * &
            !& normz - 1.0_MK) * (pz * normz - 2.0_MK) - py * normy ** 2 * (py &
            !& * normy - 2.0_MK) * pz * normz * (pz * normz - 1.0_MK) * (pz * &
            !& normz - 2.0_MK) - py * normy ** 2 * (py * normy - 1.0_MK) * pz &
            !& * normz * (pz * normz - 1.0_MK) * (pz * normz - 2.0_MK)

            !--------------------------------------------------------------------
            ! THESE ARE TESTS FOR FREE-SPACE DOMAINS
            !--------------------------------------------------------------------
            !!-------------Isolated TG vortex set
            !IF (px .GT. 0.0_MK .AND. px .LT. 1.0_MK .AND. &
                !py .GT. 0.0_MK .AND. py .LT. 1.0_MK .AND. &
                !pz .GT. 0.0_MK .AND. pz .LT. 1.0_MK) THEN
                !psif(ilevel,ipatch)%fld(1,i,j,k,isub) = &
                !& cos(PI/2.0_MK+px * normx) * sin(py * normy) * sin(pz * normz)
                !psif(ilevel,ipatch)%fld(2,i,j,k,isub) = &
                !& sin(px * normx) * cos(PI/2.0_MK+py * normy) * sin(pz * normz)
                !psif(ilevel,ipatch)%fld(3,i,j,k,isub) = &
                !& -0.2D1 * sin(px * normx) * sin(py * normy) * cos(PI/2.0_MK+pz * normz)

                !wf(ilevel,ipatch)%fld(1,i,j,k,isub) = &
                !& -(normx**2+normy**2+normz**2)*cos(PI/2.0_MK+px * normx) * sin(py * normy) * sin(pz * normz)
                !wf(ilevel,ipatch)%fld(2,i,j,k,isub) = &
                !& -(normx**2+normy**2+normz**2)*sin(px * normx) * cos(PI/2.0_MK+py * normy) * sin(pz * normz)
                !wf(ilevel,ipatch)%fld(3,i,j,k,isub) = &
                !& (normx**2+normy**2+normz**2)*0.2D1 * sin(px * normx) * sin(py * normy) * cos(PI/2.0_MK+pz * normz)
            !ENDIF

            !!-------------Isolated sines - not suitable for freespace (very discontinuous)
            !IF (px .GT. 0.0_MK .AND. px .LT. 1.0_MK .AND. &
                !py .GT. 0.0_MK .AND. py .LT. 1.0_MK .AND. &
                !pz .GT. 0.0_MK .AND. pz .LT. 1.0_MK) THEN
                !!stream function
                !psif(ilevel,ipatch)%fld(1,i,j,k,isub) = sin(py * normy)
                !psif(ilevel,ipatch)%fld(2,i,j,k,isub) = sin(pz * normz)
                !psif(ilevel,ipatch)%fld(3,i,j,k,isub) = sin(px * normx)

                !!vorticity
                !wf(ilevel,ipatch)%fld(1,i,j,k,isub) = -sin(py * normy) * normy**2
                !wf(ilevel,ipatch)%fld(2,i,j,k,isub) = -sin(pz * normz) * normz**2
                !wf(ilevel,ipatch)%fld(3,i,j,k,isub) = -sin(px * normx) * normx**2

                !!velocity
                !psif(ilevel,ipatch)%fld(1,i,j,k,isub) = -cos(pz * normz) * normz
                !psif(ilevel,ipatch)%fld(2,i,j,k,isub) = -cos(px * normx) * normx
                !psif(ilevel,ipatch)%fld(3,i,j,k,isub) = -cos(py * normy) * normy

            !ENDIF

            !!-------- 2D Perlman
            !tmp = (px-1.5_MK)**2+(py-1.5_MK)**2;
            !IF (tmp .LT. 1.0_MK) THEN
              !wf(ilevel,ipatch)%fld(3,i,j,k,isub) = (1.0_MK-tmp)**7
            !ELSE
              !wf(ilevel,ipatch)%fld(3,i,j,k,isub) = 0.0_MK
            !ENDIF
            !wf(ilevel,ipatch)%fld(1,i,j,k,isub) = 0.0_MK
            !wf(ilevel,ipatch)%fld(2,i,j,k,isub) = 0.0_MK


            !IF (tmp .LE. 1.0_MK) THEN
              !tmp = (-1.0_MK/(16.0_MK*tmp)*(1.0_MK-(1.0_MK-tmp)**8))
            !ELSE IF (tmp .GT. 1.0_MK) THEN
              !tmp = (-1.0_MK/(16.0_MK*tmp));
            !ENDIF

            !psif(ilevel,ipatch)%fld(1,i,j,k,isub) = -(py - 1.5_MK)*tmp
            !psif(ilevel,ipatch)%fld(2,i,j,k,isub) = -(px - 1.5_MK)*tmp
            !psif(ilevel,ipatch)%fld(3,i,j,k,isub) = 0.0_MK

            !psif(ilevel,ipatch)%fld(1,i,j,k,isub) = 0.0_MK
            !psif(ilevel,ipatch)%fld(2,i,j,k,isub) = 0.0_MK
            !k2 = 0.0_MK
            !k1 = k2-(-0.084933035714286)

            !IF (tmp .LE. 1.0_MK) THEN
              !psif(ilevel,ipatch)%fld(3,i,j,k,isub) = &
              !& (-(-1.0_MK/16.0_MK*tmp**16 + 4.0_MK/7.0_MK*tmp**14 - 7.0_MK/3.0_MK*tmp**12 + &
              !& 28.0_MK/5.0_MK*tmp**10 - 35.0_MK/4.0_MK*tmp**8 + 28.0_MK/3.0_MK*tmp**6 - 7.0_MK*tmp**4 &
              !& + 4.0_MK*tmp**2)/16.0_MK + k1)
            !ELSE IF (tmp .GT. 1.0_MK) THEN
              !psif(ilevel,ipatch)%fld(3,i,j,k,isub) = &
              !& (-log(tmp)/16.0_MK + k2)
            !ENDIF


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
CALL substop('naga_case_test',t0,info)


RETURN

END SUBROUTINE naga_case_test

END MODULE naga_mod_case_test

