!------------------------------------------------------------------------------
! Subroutine :  naga_strainrate.f90
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
! This routines computes the strain rate to output in diagnostics or to use
! in LES models. LES models have not been (fully) implemented.
!------------------------------------------------------------------------------

MODULE naga_mod_strainrate

IMPLICIT NONE

INTERFACE naga_strainrate
  MODULE PROCEDURE naga_strainrate
END INTERFACE

CONTAINS

SUBROUTINE naga_strainrate(info)

USE naga_mod_globals
USE naga_mod_say

IMPLICIT NONE

!-------------------------------------------------------------------------
! Arguments
!-------------------------------------------------------------------------
INTEGER, INTENT(INOUT)                                  :: info


!-------------------------------------------------------------------------
! Local variables
!-------------------------------------------------------------------------
REAL(MK)                        :: t0
INTEGER                         :: ilevel,ipatch
INTEGER                         :: isub,isubl
INTEGER                         :: i,j,k
REAL(MK)                        :: dx,dy,dz
REAL(MK)                        :: facx,facy,facz
INTEGER                         :: imax,jmax,kmax


!------------------------------------------------------------------------------
! Variables for strainrate
!------------------------------------------------------------------------------
REAL(MK)                        :: maxstrain,pmaxstrain !strain in point,patch


!----------------------------------------------------------------------------
! Initialise routine
!----------------------------------------------------------------------------
CALL substart('naga_strainrate',t0,info)


!----------------------------------------------------------------------------
! Set variables zero. p/l/g corresponds to patch/local/global
!----------------------------------------------------------------------------
lmaxstrain    = 0.0_MK


!----------------------------------------------------------------------------
! Compute maximum strain rate of sub without storing the field
!----------------------------------------------------------------------------
IF (lesmodel .EQ. 0) THEN
  DO ilevel=1,nlevels
    DO ipatch=1,npatches(ilevel)
      pmaxstrain    = 0.0_MK

      dx = ptcset(ilevel,ipatch)%dx(1)
      dy = ptcset(ilevel,ipatch)%dx(2)
      dz = ptcset(ilevel,ipatch)%dx(3)

      facx = 1.0_MK/(dx*12.0_MK)
      facy = 1.0_MK/(dy*12.0_MK)
      facz = 1.0_MK/(dz*12.0_MK)

      DO isub=1,topos(ilevel)%nsublist
        isubl=topos(ilevel)%isublist(isub)

        !determine max indicies depending on periodicity (N+1)
        IF (domainbc .EQ. 0) THEN
          IF (topos(ilevel)%subs_bc(2,isubl) .EQ. 0) THEN
            imax = ptcset(ilevel,ipatch)%snx(1,isubl)-1
          ELSE
            imax = ptcset(ilevel,ipatch)%snx(1,isubl)
          ENDIF
          IF (topos(ilevel)%subs_bc(4,isubl) .EQ. 0) THEN
            jmax = ptcset(ilevel,ipatch)%snx(2,isubl)-1
          ELSE
            jmax = ptcset(ilevel,ipatch)%snx(2,isubl)
          ENDIF
          IF (topos(ilevel)%subs_bc(6,isubl) .EQ. 0) THEN
            kmax = ptcset(ilevel,ipatch)%snx(3,isubl)-1
          ELSE
            kmax = ptcset(ilevel,ipatch)%snx(3,isubl)
          ENDIF
        ELSE
          imax = ptcset(ilevel,ipatch)%snx(1,isubl)-1
          jmax = ptcset(ilevel,ipatch)%snx(2,isubl)-1
          kmax = ptcset(ilevel,ipatch)%snx(3,isubl)-1
        ENDIF
        DO k=1,kmax
          DO j=1,jmax
            DO i=1,imax
              !-------------------------------------------------------------------
              ! Compute maximum strainrate O(4) FD
              !-------------------------------------------------------------------
              pmaxstrain    = MAX(pmaxstrain, &  !du/dx
              &(-       facx*uf(ilevel,ipatch)%fld(1,i+2,j  ,k  ,isub)   &
              & +8.0_MK*facx*uf(ilevel,ipatch)%fld(1,i+1,j  ,k  ,isub)   &
              & -8.0_MK*facx*uf(ilevel,ipatch)%fld(1,i-1,j  ,k  ,isub)   &
              & +       facx*uf(ilevel,ipatch)%fld(1,i-2,j  ,k  ,isub)))
              pmaxstrain    = MAX(pmaxstrain, &  !dv/dx
              &(-       facx*uf(ilevel,ipatch)%fld(2,i+2,j  ,k  ,isub)   &
              & +8.0_MK*facx*uf(ilevel,ipatch)%fld(2,i+1,j  ,k  ,isub)   &
              & -8.0_MK*facx*uf(ilevel,ipatch)%fld(2,i-1,j  ,k  ,isub)   &
              & +       facx*uf(ilevel,ipatch)%fld(2,i-2,j  ,k  ,isub)))
              pmaxstrain    = MAX(pmaxstrain, &  !dw/dx
              &(-       facx*uf(ilevel,ipatch)%fld(3,i+2,j  ,k  ,isub)   &
              & +8.0_MK*facx*uf(ilevel,ipatch)%fld(3,i+1,j  ,k  ,isub)   &
              & -8.0_MK*facx*uf(ilevel,ipatch)%fld(3,i-1,j  ,k  ,isub)   &
              & +       facx*uf(ilevel,ipatch)%fld(3,i-2,j  ,k  ,isub)))
              pmaxstrain    = MAX(pmaxstrain, &  !du/dy
              &(-       facy*uf(ilevel,ipatch)%fld(1,i  ,j+2,k  ,isub)   &
              & +8.0_MK*facy*uf(ilevel,ipatch)%fld(1,i  ,j+1,k  ,isub)   &
              & -8.0_MK*facy*uf(ilevel,ipatch)%fld(1,i  ,j-1,k  ,isub)   &
              & +       facy*uf(ilevel,ipatch)%fld(1,i  ,j-2,k  ,isub)))
              pmaxstrain    = MAX(pmaxstrain, &  !dv/dy
              &(-       facy*uf(ilevel,ipatch)%fld(2,i  ,j+2,k  ,isub)   &
              & +8.0_MK*facy*uf(ilevel,ipatch)%fld(2,i  ,j+1,k  ,isub)   &
              & -8.0_MK*facy*uf(ilevel,ipatch)%fld(2,i  ,j-1,k  ,isub)   &
              & +       facy*uf(ilevel,ipatch)%fld(2,i  ,j-2,k  ,isub)))
              pmaxstrain    = MAX(pmaxstrain, &  !dw/dy
              &(-       facy*uf(ilevel,ipatch)%fld(1,i  ,j+2,k  ,isub)   &
              & +8.0_MK*facy*uf(ilevel,ipatch)%fld(1,i  ,j+1,k  ,isub)   &
              & -8.0_MK*facy*uf(ilevel,ipatch)%fld(1,i  ,j-1,k  ,isub)   &
              & +       facy*uf(ilevel,ipatch)%fld(1,i  ,j-2,k  ,isub)))
              pmaxstrain    = MAX(pmaxstrain, &  !du/dz
              &(-       facz*uf(ilevel,ipatch)%fld(2,i  ,j  ,k+2,isub)   &
              & +8.0_MK*facz*uf(ilevel,ipatch)%fld(2,i  ,j  ,k+1,isub)   &
              & -8.0_MK*facz*uf(ilevel,ipatch)%fld(2,i  ,j  ,k-1,isub)   &
              & +       facz*uf(ilevel,ipatch)%fld(2,i  ,j  ,k-2,isub)))
              pmaxstrain    = MAX(pmaxstrain, &  !dv/dz
              &(-       facz*uf(ilevel,ipatch)%fld(3,i  ,j  ,k+2,isub)   &
              & +8.0_MK*facz*uf(ilevel,ipatch)%fld(3,i  ,j  ,k+1,isub)   &
              & -8.0_MK*facz*uf(ilevel,ipatch)%fld(3,i  ,j  ,k-1,isub)   &
              & +       facz*uf(ilevel,ipatch)%fld(3,i  ,j  ,k-2,isub)))
              pmaxstrain    = MAX(pmaxstrain, &  !dw/dz
              &(-       facz*uf(ilevel,ipatch)%fld(3,i  ,j  ,k+2,isub)   &
              & +8.0_MK*facz*uf(ilevel,ipatch)%fld(3,i  ,j  ,k+1,isub)   &
              & -8.0_MK*facz*uf(ilevel,ipatch)%fld(3,i  ,j  ,k-1,isub)   &
              & +       facz*uf(ilevel,ipatch)%fld(3,i  ,j  ,k-2,isub)))
            ENDDO
          ENDDO
        ENDDO
      ENDDO
      lmaxstrain    = MAX(lmaxstrain ,pmaxstrain)
    ENDDO
  ENDDO

!----------------------------------------------------------------------------
! Compute maximum strain rate and store the field
!----------------------------------------------------------------------------
ELSE IF(lesmodel .EQ. 1) THEN
  DO ilevel=1,nlevels
    DO ipatch=1,npatches(ilevel)
      pmaxstrain    = 0.0_MK

      dx = ptcset(ilevel,ipatch)%dx(1)
      dy = ptcset(ilevel,ipatch)%dx(2)
      dz = ptcset(ilevel,ipatch)%dx(3)

      facx = 1.0_MK/(dx*12.0_MK)
      facy = 1.0_MK/(dy*12.0_MK)
      facz = 1.0_MK/(dz*12.0_MK)

      DO isub=1,topos(ilevel)%nsublist
        isubl=topos(ilevel)%isublist(isub)

        !determine max indicies depending on periodicity (N+1)
        IF (domainbc .EQ. 0) THEN
          IF (topos(ilevel)%subs_bc(2,isubl) .EQ. 0) THEN
            imax = ptcset(ilevel,ipatch)%snx(1,isubl)-1
          ELSE
            imax = ptcset(ilevel,ipatch)%snx(1,isubl)
          ENDIF
          IF (topos(ilevel)%subs_bc(4,isubl) .EQ. 0) THEN
            jmax = ptcset(ilevel,ipatch)%snx(2,isubl)-1
          ELSE
            jmax = ptcset(ilevel,ipatch)%snx(2,isubl)
          ENDIF
          IF (topos(ilevel)%subs_bc(6,isubl) .EQ. 0) THEN
            kmax = ptcset(ilevel,ipatch)%snx(3,isubl)-1
          ELSE
            kmax = ptcset(ilevel,ipatch)%snx(3,isubl)
          ENDIF
        ELSE
          imax = ptcset(ilevel,ipatch)%snx(1,isubl)-1
          jmax = ptcset(ilevel,ipatch)%snx(2,isubl)-1
          kmax = ptcset(ilevel,ipatch)%snx(3,isubl)-1
        ENDIF
        DO k=1,kmax
          DO j=1,jmax
            DO i=1,imax
              maxstrain    = 0.0_MK
              !-------------------------------------------------------------------
              ! Compute maximum strainrate O(4) FD
              !-------------------------------------------------------------------
              maxstrain    = MAX(maxstrain, &  !du/dx
              &(-       facx*uf(ilevel,ipatch)%fld(1,i+2,j  ,k  ,isub)   &
              & +8.0_MK*facx*uf(ilevel,ipatch)%fld(1,i+1,j  ,k  ,isub)   &
              & -8.0_MK*facx*uf(ilevel,ipatch)%fld(1,i-1,j  ,k  ,isub)   &
              & +       facx*uf(ilevel,ipatch)%fld(1,i-2,j  ,k  ,isub)))
              maxstrain    = MAX(maxstrain, &  !dv/dx
              &(-       facx*uf(ilevel,ipatch)%fld(2,i+2,j  ,k  ,isub)   &
              & +8.0_MK*facx*uf(ilevel,ipatch)%fld(2,i+1,j  ,k  ,isub)   &
              & -8.0_MK*facx*uf(ilevel,ipatch)%fld(2,i-1,j  ,k  ,isub)   &
              & +       facx*uf(ilevel,ipatch)%fld(2,i-2,j  ,k  ,isub)))
              maxstrain    = MAX(maxstrain, &  !dw/dx
              &(-       facx*uf(ilevel,ipatch)%fld(3,i+2,j  ,k  ,isub)   &
              & +8.0_MK*facx*uf(ilevel,ipatch)%fld(3,i+1,j  ,k  ,isub)   &
              & -8.0_MK*facx*uf(ilevel,ipatch)%fld(3,i-1,j  ,k  ,isub)   &
              & +       facx*uf(ilevel,ipatch)%fld(3,i-2,j  ,k  ,isub)))
              maxstrain    = MAX(maxstrain, &  !du/dy
              &(-       facy*uf(ilevel,ipatch)%fld(1,i  ,j+2,k  ,isub)   &
              & +8.0_MK*facy*uf(ilevel,ipatch)%fld(1,i  ,j+1,k  ,isub)   &
              & -8.0_MK*facy*uf(ilevel,ipatch)%fld(1,i  ,j-1,k  ,isub)   &
              & +       facy*uf(ilevel,ipatch)%fld(1,i  ,j-2,k  ,isub)))
              maxstrain    = MAX(maxstrain, &  !dv/dy
              &(-       facy*uf(ilevel,ipatch)%fld(2,i  ,j+2,k  ,isub)   &
              & +8.0_MK*facy*uf(ilevel,ipatch)%fld(2,i  ,j+1,k  ,isub)   &
              & -8.0_MK*facy*uf(ilevel,ipatch)%fld(2,i  ,j-1,k  ,isub)   &
              & +       facy*uf(ilevel,ipatch)%fld(2,i  ,j-2,k  ,isub)))
              maxstrain    = MAX(maxstrain, &  !dw/dy
              &(-       facy*uf(ilevel,ipatch)%fld(1,i  ,j+2,k  ,isub)   &
              & +8.0_MK*facy*uf(ilevel,ipatch)%fld(1,i  ,j+1,k  ,isub)   &
              & -8.0_MK*facy*uf(ilevel,ipatch)%fld(1,i  ,j-1,k  ,isub)   &
              & +       facy*uf(ilevel,ipatch)%fld(1,i  ,j-2,k  ,isub)))
              maxstrain    = MAX(maxstrain, &  !du/dz
              &(-       facz*uf(ilevel,ipatch)%fld(2,i  ,j  ,k+2,isub)   &
              & +8.0_MK*facz*uf(ilevel,ipatch)%fld(2,i  ,j  ,k+1,isub)   &
              & -8.0_MK*facz*uf(ilevel,ipatch)%fld(2,i  ,j  ,k-1,isub)   &
              & +       facz*uf(ilevel,ipatch)%fld(2,i  ,j  ,k-2,isub)))
              maxstrain    = MAX(maxstrain, &  !dv/dz
              &(-       facz*uf(ilevel,ipatch)%fld(3,i  ,j  ,k+2,isub)   &
              & +8.0_MK*facz*uf(ilevel,ipatch)%fld(3,i  ,j  ,k+1,isub)   &
              & -8.0_MK*facz*uf(ilevel,ipatch)%fld(3,i  ,j  ,k-1,isub)   &
              & +       facz*uf(ilevel,ipatch)%fld(3,i  ,j  ,k-2,isub)))
              maxstrain    = MAX(maxstrain, &  !dw/dz
              &(-       facz*uf(ilevel,ipatch)%fld(3,i  ,j  ,k+2,isub)   &
              & +8.0_MK*facz*uf(ilevel,ipatch)%fld(3,i  ,j  ,k+1,isub)   &
              & -8.0_MK*facz*uf(ilevel,ipatch)%fld(3,i  ,j  ,k-1,isub)   &
              & +       facz*uf(ilevel,ipatch)%fld(3,i  ,j  ,k-2,isub)))
              ! Store strain rate in field sf
              sf(ilevel,ipatch)%fld(i,j,k,isub) = maxstrain
              ! Get max value for the sub/patch
              pmaxstrain    = MAX(pmaxstrain,maxstrain)
            ENDDO
          ENDDO
        ENDDO
      ENDDO
      lmaxstrain    = MAX(lmaxstrain ,pmaxstrain)
    ENDDO
  ENDDO
ENDIF

!----------------------------------------------------------------------------
! Return
!----------------------------------------------------------------------------
 9999 CONTINUE
CALL substop('naga_strainrate',t0,info)
RETURN


END SUBROUTINE naga_strainrate

END MODULE naga_mod_strainrate

