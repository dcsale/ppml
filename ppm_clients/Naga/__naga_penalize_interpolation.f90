!------------------------------------------------------------------------------
! Subroutine : naga_penalize_interpolation.f90
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
! This routines performs interpolation penalization, cf Rasmussen:2011 doing
! a replacement of the vorticity field in and around the solid using FD curls.
! The order of the FD stencil is set in the variable 'penalizationscheme'
!------------------------------------------------------------------------------
MODULE naga_mod_penalize_interpolation
IMPLICIT NONE
INTERFACE naga_penalize_interpolation
  MODULE PROCEDURE naga_penalize_interpolation
END INTERFACE
CONTAINS
SUBROUTINE naga_penalize_interpolation(info)
USE naga_mod_globals
USE naga_mod_say
IMPLICIT NONE
!-------------------------------------------------------------------------
! Arguments
!-------------------------------------------------------------------------
INTEGER, INTENT(OUT) :: info
!-------------------------------------------------------------------------
! Local variables
!-------------------------------------------------------------------------
REAL(MK) :: t0
INTEGER :: ilevel,ipatch
REAL(MK) :: dx,dy,dz
REAL(MK) :: facx,facy,facz
INTEGER :: isub,isubl
INTEGER :: i,j,k
!----------------------------------------------------------------------------
! Initialise routine
!----------------------------------------------------------------------------
CALL substart('naga_penalize_interpolation',t0,info)
info = 0
!----------------------------------------------------------------------------
! Penalize patches
!----------------------------------------------------------------------------
DO ilevel=1,nlevels
  DO ipatch=1,npatches(ilevel)
    dx = ptcset(ilevel,ipatch)%dx(1)
    dy = ptcset(ilevel,ipatch)%dx(2)
    dz = ptcset(ilevel,ipatch)%dx(3)
    !-----------------------------------------------------------------------
    ! Do the actual interpolation and store in velocity array
    ! The interpolation is done as a correction, therefore the interpolated
    ! velocity field is subtracted the unpenalized velocity field and the
    ! curl of this difference is to be added to the vorticity field
    ! delta u = chi*(ubar-u)
    ! omega = omega + curl(delta u)
    !-----------------------------------------------------------------------
    DO isub=1,topos(ilevel)%nsublist
      isubl=topos(ilevel)%isublist(isub)
      DO k=1-gstw(1),ptcset(ilevel,ipatch)%snx(3,isubl)+gstw(1)
        DO j=1-gstw(2),ptcset(ilevel,ipatch)%snx(2,isubl)+gstw(2)
          DO i=1-gstw(3),ptcset(ilevel,ipatch)%snx(1,isubl)+gstw(3)
            uf(ilevel,ipatch)%fld(1,i,j,k,isub) = &
             & ubarf(ilevel,ipatch)%fld(1,i,j,k,isub) &
             & *chif(ilevel,ipatch)%fld(i,j,k,isub) &
             & +uf (ilevel,ipatch)%fld(1,i,j,k,isub) &
             & *(1.0_MK - chif(ilevel,ipatch)%fld(i,j,k,isub))
            uf(ilevel,ipatch)%fld(2,i,j,k,isub) = &
             & ubarf(ilevel,ipatch)%fld(2,i,j,k,isub) &
             & *chif(ilevel,ipatch)%fld(i,j,k,isub) &
             & +uf (ilevel,ipatch)%fld(2,i,j,k,isub) &
             & *(1.0_MK - chif(ilevel,ipatch)%fld(i,j,k,isub))
            uf(ilevel,ipatch)%fld(3,i,j,k,isub) = &
             & ubarf(ilevel,ipatch)%fld(3,i,j,k,isub) &
             & *chif(ilevel,ipatch)%fld(i,j,k,isub) &
             & +uf (ilevel,ipatch)%fld(3,i,j,k,isub) &
             & *(1.0_MK - chif(ilevel,ipatch)%fld(i,j,k,isub))
          ENDDO
        ENDDO
      ENDDO
    ENDDO !isub
    !-----------------------------------------------------------------------
    ! Take curl, 2nd order FD
    !-----------------------------------------------------------------------
    IF (penalizationscheme .EQ. 1) THEN
      facx = 1.0_MK/(2.0_MK*dx)
      facy = 1.0_MK/(2.0_MK*dy)
      facz = 1.0_MK/(2.0_MK*dz)
      DO isub=1,topos(ilevel)%nsublist
        isubl=topos(ilevel)%isublist(isub)
        DO k=1,ptcset(ilevel,ipatch)%snx(3,isubl)
          DO j=1,ptcset(ilevel,ipatch)%snx(2,isubl)
            DO i=1,ptcset(ilevel,ipatch)%snx(1,isubl)
              wf(ilevel,ipatch)%fld(1,i,j,k,isub) = &
              & facy*(uf(ilevel,ipatch)%fld(3,i ,j+1,k ,isub)- &
                     & uf(ilevel,ipatch)%fld(3,i ,j-1,k ,isub)) &
              & -facz*(uf(ilevel,ipatch)%fld(2,i ,j ,k+1,isub)- &
                     & uf(ilevel,ipatch)%fld(2,i ,j ,k-1,isub))
              wf(ilevel,ipatch)%fld(2,i,j,k,isub) = &
              & facz*(uf(ilevel,ipatch)%fld(1,i ,j ,k+1,isub)- &
                     & uf(ilevel,ipatch)%fld(1,i ,j ,k-1,isub)) &
              & -facx*(uf(ilevel,ipatch)%fld(3,i+1,j ,k ,isub)- &
                     & uf(ilevel,ipatch)%fld(3,i-1,j ,k ,isub))
              wf(ilevel,ipatch)%fld(3,i,j,k,isub) = &
              & facx*(uf(ilevel,ipatch)%fld(2,i+1,j ,k ,isub)- &
                     & uf(ilevel,ipatch)%fld(2,i-1,j ,k ,isub)) &
              & -facy*(uf(ilevel,ipatch)%fld(1,i ,j+1,k ,isub)- &
                     & uf(ilevel,ipatch)%fld(1,i ,j-1,k ,isub))
            ENDDO
          ENDDO
        ENDDO
      ENDDO !isub
    !-----------------------------------------------------------------------
    ! Take curl, 4th order FD
    !-----------------------------------------------------------------------
    ELSE IF (penalizationscheme .EQ. 2) THEN
      facx = 1.0_MK/(12.0_MK*dx)
      facy = 1.0_MK/(12.0_MK*dy)
      facz = 1.0_MK/(12.0_MK*dz)
      DO isub=1,topos(ilevel)%nsublist
        isubl=topos(ilevel)%isublist(isub)
        DO k=1,ptcset(ilevel,ipatch)%snx(3,isubl)
          DO j=1,ptcset(ilevel,ipatch)%snx(2,isubl)
            DO i=1,ptcset(ilevel,ipatch)%snx(1,isubl)
              wf(ilevel,ipatch)%fld(1,i,j,k,isub) = &
              & facy*( -uf(ilevel,ipatch)%fld(3,i ,j+2,k ,isub) &
                     & +8.0_MK*uf(ilevel,ipatch)%fld(3,i ,j+1,k ,isub) &
                     & -8.0_MK*uf(ilevel,ipatch)%fld(3,i ,j-1,k ,isub) &
                     & +uf(ilevel,ipatch)%fld(3,i ,j-2,k ,isub)) &
              & -facz*( -uf(ilevel,ipatch)%fld(2,i ,j ,k+2,isub) &
                     & +8.0_MK*uf(ilevel,ipatch)%fld(2,i ,j ,k+1,isub) &
                     & -8.0_MK*uf(ilevel,ipatch)%fld(2,i ,j ,k-1,isub) &
                     & +uf(ilevel,ipatch)%fld(2,i ,j ,k-2,isub))
              wf(ilevel,ipatch)%fld(2,i,j,k,isub) = &
              & facz*( -uf(ilevel,ipatch)%fld(1,i ,j ,k+2,isub) &
                     & +8.0_MK*uf(ilevel,ipatch)%fld(1,i ,j ,k+1,isub) &
                     & -8.0_MK*uf(ilevel,ipatch)%fld(1,i ,j ,k-1,isub) &
                     & +uf(ilevel,ipatch)%fld(1,i ,j ,k-2,isub)) &
              & -facx*( -uf(ilevel,ipatch)%fld(3,i+2,j ,k ,isub) &
                     & +8.0_MK*uf(ilevel,ipatch)%fld(3,i+1,j ,k ,isub) &
                     & -8.0_MK*uf(ilevel,ipatch)%fld(3,i-1,j ,k ,isub) &
                     & +uf(ilevel,ipatch)%fld(3,i-2,j ,k ,isub))
              wf(ilevel,ipatch)%fld(3,i,j,k,isub) = &
              & facx*( -uf(ilevel,ipatch)%fld(2,i+2,j ,k ,isub) &
                     & +8.0_MK*uf(ilevel,ipatch)%fld(2,i+1,j ,k ,isub) &
                     & -8.0_MK*uf(ilevel,ipatch)%fld(2,i-1,j ,k ,isub) &
                     & +uf(ilevel,ipatch)%fld(2,i-2,j ,k ,isub)) &
              & -facy*( -uf(ilevel,ipatch)%fld(1,i ,j+2,k ,isub) &
                     & +8.0_MK*uf(ilevel,ipatch)%fld(1,i ,j+1,k ,isub) &
                     & -8.0_MK*uf(ilevel,ipatch)%fld(1,i ,j-1,k ,isub) &
                     & +uf(ilevel,ipatch)%fld(1,i ,j-2,k ,isub))
            ENDDO
          ENDDO
        ENDDO
      ENDDO !isub
    ENDIF
  ENDDO
ENDDO
!----------------------------------------------------------------------------
! Return
!----------------------------------------------------------------------------
 9999 CONTINUE
CALL substop('naga_penalize_interpolation',t0,info)
RETURN
END SUBROUTINE naga_penalize_interpolation
END MODULE naga_mod_penalize_interpolation
