!------------------------------------------------------------------------------
! Subroutine : naga_diagnostics.f90
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
! This routines computes various diagnostics and outputs them to a file
! File handle 10 is used
!
! time step
! time
! +Maximum absolute velocity
! +kinetic energy
! +Sum (integral) of vorticity (3 components)
! +Sum (integral) of magnitude of vorticity
! +Enstrophy
! Fourier number
! CFL number
! solid force
! solid force moment
! +Divergence of velocity : Maximum absolute, integral of absolute
! +Divergence of vorticity: Maximum absolute, integral of absolute
!------------------------------------------------------------------------------
MODULE naga_mod_diagnostics
IMPLICIT NONE
INTERFACE naga_diagnostics
  MODULE PROCEDURE naga_diagnostics
END INTERFACE
CONTAINS
SUBROUTINE naga_diagnostics(info)
USE naga_mod_output
USE naga_mod_globals
USE naga_mod_say
IMPLICIT NONE
!-------------------------------------------------------------------------
! Arguments
!-------------------------------------------------------------------------
INTEGER, INTENT(INOUT) :: info
!-------------------------------------------------------------------------
! Local variables
!-------------------------------------------------------------------------
REAL(MK) :: t0
INTEGER :: ilevel,ipatch
INTEGER :: isub,isubl
INTEGER :: i,j,k
INTEGER :: ios
REAL(MK) :: dx,dy,dz
REAL(MK) :: px,py,pz
REAL(MK) :: facx,facy,facz
CHARACTER(256) :: filename
REAL(MK) :: tmp,tmp2,tmp3
LOGICAL :: firsttime
INTEGER :: imax,jmax,kmax
!-------------------------------------------------------------------------
! Variables for diagnostics
!-------------------------------------------------------------------------
REAL(MK),DIMENSION(ncom) :: ptotvort,ltotvort,gtotvort !sum vort
REAL(MK),DIMENSION(ncom) :: pkinenrg,lkinenrg,gkinenrg !sum kin energy
REAL(MK),DIMENSION(ncom) :: pforce,lforce,gforce !sum force
REAL(MK),DIMENSION(ncom) :: pmoment,lmoment,gmoment !sum moment
REAL(MK) :: ptotabsvort,ltotabsvort,gtotabsvort!sum abs(vort)
REAL(MK) :: penstrophy,lenstrophy,genstrophy !enstrophy
REAL(MK) :: pmaxabsvel,lmaxabsvel,gmaxabsvel !max(abs(vel)
REAL(MK) :: pmaxdivvel,lmaxdivvel,gmaxdivvel !max(div(vel))
REAL(MK) :: pmaxdivvort,lmaxdivvort,gmaxdivvort!max(div(vort))
REAL(MK) :: pintdivvel,lintdivvel,gintdivvel !int(div(vel))
REAL(MK) :: pintdivvort,lintdivvort,gintdivvort!int(div(vort))
REAL(MK) :: gmaxstrain !max(strain)
REAL(MK) :: pinterot,linterot,ginterot !curl(vel)-vrt
REAL(MK) :: lfourier,gfourier
REAL(MK) :: lcourant,gcourant
REAL(MK) :: gkinenrgdt
REAL(MK) :: gvisceff
!----------------------------------------------------------------------------
! Initialise routine
!----------------------------------------------------------------------------
CALL substart('naga_diagnostics',t0,info)
!----------------------------------------------------------------------------
! Set variables zero. p/l/g corresponds to patch/local/global
!----------------------------------------------------------------------------
ltotvort = 0.0_MK
gtotvort = 0.0_MK
lkinenrg = 0.0_MK
gkinenrg = 0.0_MK
lforce = 0.0_MK
gforce = 0.0_MK
lmoment = 0.0_MK
gmoment = 0.0_MK
ltotabsvort = 0.0_MK
gtotabsvort = 0.0_MK
lenstrophy = 0.0_MK
genstrophy = 0.0_MK
lmaxabsvel = 0.0_MK
gmaxabsvel = 0.0_MK
lfourier = 0.0_MK
lcourant = 0.0_MK
lintdivvel = 0.0_MK
lintdivvort = 0.0_MK
lmaxdivvel = 0.0_MK
lmaxdivvort = 0.0_MK
!!lmaxstrain = 0.0_MK
linterot = 0.0_MK
!----------------------------------------------------------------------------
! Do local sums
!----------------------------------------------------------------------------
DO ilevel=1,nlevels
  DO ipatch=1,npatches(ilevel)
    ptotvort = 0.0_MK
    pkinenrg = 0.0_MK
    ptotabsvort = 0.0_MK
    penstrophy = 0.0_MK
    pmaxabsvel = 0.0_MK
    pintdivvel = 0.0_MK
    !! pintdivvort = 0.0_MK
    pmaxdivvel = 0.0_MK
    !! pmaxdivvort = 0.0_MK
    !!pmaxstrain = 0.0_MK
    pinterot = 0.0_MK
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
            ! Do local sums
            !-------------------------------------------------------------------
            ! Max absolute velocity - to be post processed
            ! and total of kinetic energy components - to be post processed
            !-------------------------------------------------------------------
            tmp = uf(ilevel,ipatch)%fld(1,i,j,k,isub)**2
            tmp2 = uf(ilevel,ipatch)%fld(2,i,j,k,isub)**2
            tmp3 = uf(ilevel,ipatch)%fld(3,i,j,k,isub)**2
            pmaxabsvel = MAX(tmp+tmp2+tmp3,pmaxabsvel)
            pkinenrg(1) = pkinenrg(1) + tmp
            pkinenrg(2) = pkinenrg(2) + tmp2
            pkinenrg(3) = pkinenrg(3) + tmp3
            !-------------------------------------------------------------------
            ! Sum of vorticity
            !-------------------------------------------------------------------
            ptotvort(1) = ptotvort(1) + wf(ilevel,ipatch)%fld(1,i,j,k,isub)
            ptotvort(2) = ptotvort(2) + wf(ilevel,ipatch)%fld(2,i,j,k,isub)
            ptotvort(3) = ptotvort(3) + wf(ilevel,ipatch)%fld(3,i,j,k,isub)
            !-------------------------------------------------------------------
            ! Sum of absolute of vorticity
            ! and enstrophy
            !-------------------------------------------------------------------
            tmp = &
                 & (wf(ilevel,ipatch)%fld(1,i,j,k,isub)**2 + &
                 & wf(ilevel,ipatch)%fld(2,i,j,k,isub)**2 + &
                 & wf(ilevel,ipatch)%fld(3,i,j,k,isub)**2)
            penstrophy = penstrophy + tmp
            ptotabsvort = ptotabsvort + SQRT(tmp)
            !-------------------------------------------------------------------
            ! Divergence of velocity (4th order FD)
            !-------------------------------------------------------------------
            tmp = &
            &(- facx*uf(ilevel,ipatch)%fld(1,i+2,j ,k ,isub) &
            & +8.0_MK*facx*uf(ilevel,ipatch)%fld(1,i+1,j ,k ,isub) &
            & -8.0_MK*facx*uf(ilevel,ipatch)%fld(1,i-1,j ,k ,isub) &
            & + facx*uf(ilevel,ipatch)%fld(1,i-2,j ,k ,isub)) +&
            &(- facy*uf(ilevel,ipatch)%fld(2,i ,j+2,k ,isub) &
            & +8.0_MK*facy*uf(ilevel,ipatch)%fld(2,i ,j+1,k ,isub) &
            & -8.0_MK*facy*uf(ilevel,ipatch)%fld(2,i ,j-1,k ,isub) &
            & + facy*uf(ilevel,ipatch)%fld(2,i ,j-2,k ,isub)) +&
            &(- facz*uf(ilevel,ipatch)%fld(3,i ,j ,k+2,isub) &
            & +8.0_MK*facz*uf(ilevel,ipatch)%fld(3,i ,j ,k+1,isub) &
            & -8.0_MK*facz*uf(ilevel,ipatch)%fld(3,i ,j ,k-1,isub) &
            & + facz*uf(ilevel,ipatch)%fld(3,i ,j ,k-2,isub))
            pintdivvel = pintdivvel + ABS(tmp)
            pmaxdivvel = MAX(pmaxdivvel,ABS(tmp))
            !-------------------------------------------------------------------
            ! Divergence of vorticity (4th order FD)
            !@ consider only doing vorticity divergence 2 points inside the
            ! domain such that ghosting is not necessary...
            !-------------------------------------------------------------------
            !! tmp = &
            !! &(- facx*wf(ilevel,ipatch)%fld(1,i+2,j ,k ,isub) &
            !! & +8.0_MK*facx*wf(ilevel,ipatch)%fld(1,i+1,j ,k ,isub) &
            !! & -8.0_MK*facx*wf(ilevel,ipatch)%fld(1,i-1,j ,k ,isub) &
            !! & + facx*wf(ilevel,ipatch)%fld(1,i-2,j ,k ,isub)) +&
            !! &(- facy*wf(ilevel,ipatch)%fld(2,i ,j+2,k ,isub) &
            !! & +8.0_MK*facy*wf(ilevel,ipatch)%fld(2,i ,j+1,k ,isub) &
            !! & -8.0_MK*facy*wf(ilevel,ipatch)%fld(2,i ,j-1,k ,isub) &
            !! & + facy*wf(ilevel,ipatch)%fld(2,i ,j-2,k ,isub)) +&
            !! &(- facz*wf(ilevel,ipatch)%fld(3,i ,j ,k+2,isub) &
            !! & +8.0_MK*facz*wf(ilevel,ipatch)%fld(3,i ,j ,k+1,isub) &
            !! & -8.0_MK*facz*wf(ilevel,ipatch)%fld(3,i ,j ,k-1,isub) &
            !! & + facz*wf(ilevel,ipatch)%fld(3,i ,j ,k-2,isub))
            !! pintdivvort = pintdivvort + ABS(tmp)
            !! pmaxdivvort = MAX(pmaxdivvort,ABS(tmp))
            !-------------------------------------------------------------------
            ! Compute maximum strainrate O(4) FD
            !-------------------------------------------------------------------
            !!pmaxstrain = MAX(pmaxstrain, & !du/dx
            !!&(- facx*uf(ilevel,ipatch)%fld(1,i+2,j ,k ,isub) &
            !!& +8.0_MK*facx*uf(ilevel,ipatch)%fld(1,i+1,j ,k ,isub) &
            !!& -8.0_MK*facx*uf(ilevel,ipatch)%fld(1,i-1,j ,k ,isub) &
            !!& + facx*uf(ilevel,ipatch)%fld(1,i-2,j ,k ,isub)))
            !!pmaxstrain = MAX(pmaxstrain, & !dv/dx
            !!&(- facx*uf(ilevel,ipatch)%fld(2,i+2,j ,k ,isub) &
            !!& +8.0_MK*facx*uf(ilevel,ipatch)%fld(2,i+1,j ,k ,isub) &
            !!& -8.0_MK*facx*uf(ilevel,ipatch)%fld(2,i-1,j ,k ,isub) &
            !!& + facx*uf(ilevel,ipatch)%fld(2,i-2,j ,k ,isub)))
            !!pmaxstrain = MAX(pmaxstrain, & !dw/dx
            !!&(- facx*uf(ilevel,ipatch)%fld(3,i+2,j ,k ,isub) &
            !!& +8.0_MK*facx*uf(ilevel,ipatch)%fld(3,i+1,j ,k ,isub) &
            !!& -8.0_MK*facx*uf(ilevel,ipatch)%fld(3,i-1,j ,k ,isub) &
            !!& + facx*uf(ilevel,ipatch)%fld(3,i-2,j ,k ,isub)))
            !!pmaxstrain = MAX(pmaxstrain, & !du/dy
            !!&(- facy*uf(ilevel,ipatch)%fld(1,i ,j+2,k ,isub) &
            !!& +8.0_MK*facy*uf(ilevel,ipatch)%fld(1,i ,j+1,k ,isub) &
            !!& -8.0_MK*facy*uf(ilevel,ipatch)%fld(1,i ,j-1,k ,isub) &
            !!& + facy*uf(ilevel,ipatch)%fld(1,i ,j-2,k ,isub)))
            !!pmaxstrain = MAX(pmaxstrain, & !dv/dy
            !!&(- facy*uf(ilevel,ipatch)%fld(2,i ,j+2,k ,isub) &
            !!& +8.0_MK*facy*uf(ilevel,ipatch)%fld(2,i ,j+1,k ,isub) &
            !!& -8.0_MK*facy*uf(ilevel,ipatch)%fld(2,i ,j-1,k ,isub) &
            !!& + facy*uf(ilevel,ipatch)%fld(2,i ,j-2,k ,isub)))
            !!pmaxstrain = MAX(pmaxstrain, & !dw/dy
            !!&(- facy*uf(ilevel,ipatch)%fld(1,i ,j+2,k ,isub) &
            !!& +8.0_MK*facy*uf(ilevel,ipatch)%fld(1,i ,j+1,k ,isub) &
            !!& -8.0_MK*facy*uf(ilevel,ipatch)%fld(1,i ,j-1,k ,isub) &
            !!& + facy*uf(ilevel,ipatch)%fld(1,i ,j-2,k ,isub)))
            !!pmaxstrain = MAX(pmaxstrain, & !du/dz
            !!&(- facz*uf(ilevel,ipatch)%fld(2,i ,j ,k+2,isub) &
            !!& +8.0_MK*facz*uf(ilevel,ipatch)%fld(2,i ,j ,k+1,isub) &
            !!& -8.0_MK*facz*uf(ilevel,ipatch)%fld(2,i ,j ,k-1,isub) &
            !!& + facz*uf(ilevel,ipatch)%fld(2,i ,j ,k-2,isub)))
            !!pmaxstrain = MAX(pmaxstrain, & !dv/dz
            !!&(- facz*uf(ilevel,ipatch)%fld(3,i ,j ,k+2,isub) &
            !!& +8.0_MK*facz*uf(ilevel,ipatch)%fld(3,i ,j ,k+1,isub) &
            !!& -8.0_MK*facz*uf(ilevel,ipatch)%fld(3,i ,j ,k-1,isub) &
            !!& + facz*uf(ilevel,ipatch)%fld(3,i ,j ,k-2,isub)))
            !!pmaxstrain = MAX(pmaxstrain, & !dw/dz
            !!&(- facz*uf(ilevel,ipatch)%fld(3,i ,j ,k+2,isub) &
            !!& +8.0_MK*facz*uf(ilevel,ipatch)%fld(3,i ,j ,k+1,isub) &
            !!& -8.0_MK*facz*uf(ilevel,ipatch)%fld(3,i ,j ,k-1,isub) &
            !!& + facz*uf(ilevel,ipatch)%fld(3,i ,j ,k-2,isub)))
            !-------------------------------------------------------------------
            ! Compute omega-curl(velocity)
            !-------------------------------------------------------------------
            pinterot = pinterot &
            & + wf(ilevel,ipatch)%fld(1,i ,j ,k ,isub) - &
            &((- facy*uf(ilevel,ipatch)%fld(3,i ,j+2,k ,isub) &
            & +8.0_MK*facy*uf(ilevel,ipatch)%fld(3,i ,j+1,k ,isub) &
            & -8.0_MK*facy*uf(ilevel,ipatch)%fld(3,i ,j-1,k ,isub) &
            & + facy*uf(ilevel,ipatch)%fld(3,i ,j-2,k ,isub)) -&
            & (- facz*uf(ilevel,ipatch)%fld(2,i ,j ,k+2,isub) &
            & +8.0_MK*facz*uf(ilevel,ipatch)%fld(2,i ,j ,k+1,isub) &
            & -8.0_MK*facz*uf(ilevel,ipatch)%fld(2,i ,j ,k-1,isub) &
            & + facz*uf(ilevel,ipatch)%fld(2,i ,j ,k-2,isub))) &
            & + wf(ilevel,ipatch)%fld(2,i ,j ,k ,isub) - &
            &((- facz*uf(ilevel,ipatch)%fld(1,i ,j ,k+2,isub) &
            & +8.0_MK*facz*uf(ilevel,ipatch)%fld(1,i ,j ,k+1,isub) &
            & -8.0_MK*facz*uf(ilevel,ipatch)%fld(1,i ,j ,k-1,isub) &
            & + facz*uf(ilevel,ipatch)%fld(1,i ,j ,k-2,isub)) -&
            & (- facx*uf(ilevel,ipatch)%fld(3,i+2,j ,k ,isub) &
            & +8.0_MK*facx*uf(ilevel,ipatch)%fld(3,i+1,j ,k ,isub) &
            & -8.0_MK*facx*uf(ilevel,ipatch)%fld(3,i-1,j ,k ,isub) &
            & + facx*uf(ilevel,ipatch)%fld(3,i-2,j ,k ,isub))) &
            & + wf(ilevel,ipatch)%fld(3,i ,j ,k ,isub) - &
            &((- facx*uf(ilevel,ipatch)%fld(2,i+2,j ,k ,isub) &
            & +8.0_MK*facx*uf(ilevel,ipatch)%fld(2,i+1,j ,k ,isub) &
            & -8.0_MK*facx*uf(ilevel,ipatch)%fld(2,i-1,j ,k ,isub) &
            & + facx*uf(ilevel,ipatch)%fld(2,i-2,j ,k ,isub)) -&
            & (- facy*uf(ilevel,ipatch)%fld(1,i ,j+2,k ,isub) &
            & +8.0_MK*facy*uf(ilevel,ipatch)%fld(1,i ,j+1,k ,isub) &
            & -8.0_MK*facy*uf(ilevel,ipatch)%fld(1,i ,j-1,k ,isub) &
            & + facy*uf(ilevel,ipatch)%fld(1,i ,j-2,k ,isub)))
          ENDDO
        ENDDO
      ENDDO
    ENDDO
    !----------------------------------------------------------------------------
    ! Postprocessing
    !----------------------------------------------------------------------------
    pmaxabsvel = SQRT(pmaxabsvel)
    !penstrophy = penstrophy
    pkinenrg = 0.5_MK*pkinenrg
    !----------------------------------------------------------------------------
    ! Add patch values
    !----------------------------------------------------------------------------
    ltotvort = ltotvort + ptotvort *dx*dy*dz
    lkinenrg = lkinenrg + pkinenrg *dx*dy*dz
    ltotabsvort = ltotabsvort + ptotabsvort *dx*dy*dz
    lenstrophy = lenstrophy + penstrophy *dx*dy*dz
    lmaxabsvel = MAX(lmaxabsvel,pmaxabsvel)
    lintdivvel = lintdivvel + pintdivvel *dx*dy*dz
    !! lintdivvort = lintdivvort + pintdivvort *dx*dy*dz
    lmaxdivvel = MAX(lmaxdivvel ,pmaxdivvel)
    !! lmaxdivvort = MAX(lmaxdivvort,pmaxdivvort)
    !! lmaxstrain = MAX(lmaxstrain ,pmaxstrain)
    linterot = linterot + pinterot *dx*dy*dz
    !----------------------------------------------------------------------------
    ! Misc
    !----------------------------------------------------------------------------
    lfourier = MAX(lfourier,nu*dtime/MIN(dx**2,MIN(dy**2,dz**2)))
    lcourant = MAX(lcourant,pmaxabsvel*dtime/MIN(dx,MIN(dy,dz)))
  ENDDO
ENDDO
!----------------------------------------------------------------------------
! Do local sums with a FD-stencil-width offset into the domain
!----------------------------------------------------------------------------
DO ilevel=1,nlevels
  DO ipatch=1,npatches(ilevel)
    pintdivvort = 0.0_MK
    pmaxdivvort = 0.0_MK
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
        imax = ptcset(ilevel,ipatch)%snx(1,isubl)-1-gstw(1)
        jmax = ptcset(ilevel,ipatch)%snx(2,isubl)-1-gstw(2)
        kmax = ptcset(ilevel,ipatch)%snx(3,isubl)-1-gstw(3)
      ENDIF
      DO k=1+gstw(1),kmax
        DO j=1+gstw(2),jmax
          DO i=1+gstw(3),imax
            !-------------------------------------------------------------------
            ! Divergence of vorticity (4th order FD)
            !-------------------------------------------------------------------
            tmp = &
            &(- facx*wf(ilevel,ipatch)%fld(1,i+2,j ,k ,isub) &
            & +8.0_MK*facx*wf(ilevel,ipatch)%fld(1,i+1,j ,k ,isub) &
            & -8.0_MK*facx*wf(ilevel,ipatch)%fld(1,i-1,j ,k ,isub) &
            & + facx*wf(ilevel,ipatch)%fld(1,i-2,j ,k ,isub)) +&
            &(- facy*wf(ilevel,ipatch)%fld(2,i ,j+2,k ,isub) &
            & +8.0_MK*facy*wf(ilevel,ipatch)%fld(2,i ,j+1,k ,isub) &
            & -8.0_MK*facy*wf(ilevel,ipatch)%fld(2,i ,j-1,k ,isub) &
            & + facy*wf(ilevel,ipatch)%fld(2,i ,j-2,k ,isub)) +&
            &(- facz*wf(ilevel,ipatch)%fld(3,i ,j ,k+2,isub) &
            & +8.0_MK*facz*wf(ilevel,ipatch)%fld(3,i ,j ,k+1,isub) &
            & -8.0_MK*facz*wf(ilevel,ipatch)%fld(3,i ,j ,k-1,isub) &
            & + facz*wf(ilevel,ipatch)%fld(3,i ,j ,k-2,isub))
            pintdivvort = pintdivvort + ABS(tmp)
            pmaxdivvort = MAX(pmaxdivvort,ABS(tmp))
          ENDDO
        ENDDO
      ENDDO
    ENDDO
    lintdivvort = lintdivvort + pintdivvort *dx*dy*dz
    lmaxdivvort = MAX(lmaxdivvort,pmaxdivvort)
  ENDDO
ENDDO
!-------------------------------------------------------------------
! Penalization forces
!-------------------------------------------------------------------
IF (penalization .GT. 0) THEN
  DO ilevel=1,nlevels
    DO ipatch=1,npatches(ilevel)
      pforce = 0.0_MK
      pmoment = 0.0_MK
      dx = ptcset(ilevel,ipatch)%dx(1)
      dy = ptcset(ilevel,ipatch)%dx(2)
      dz = ptcset(ilevel,ipatch)%dx(3)
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
          pz = topos(ilevel)%min_subd(3,isubl) + REAL(k-1,MK)*dz
          DO j=1,jmax
            py = topos(ilevel)%min_subd(2,isubl) + REAL(j-1,MK)*dy
            DO i=1,imax
              px = topos(ilevel)%min_subd(1,isubl) + REAL(i-1,MK)*dx
              tmp = -(ubarf(ilevel,ipatch)%fld(1,i,j,k,isub) &
                        & -uf (ilevel,ipatch)%fld(1,i,j,k,isub)) &
                        & * chif(ilevel,ipatch)%fld( i,j,k,isub)
              tmp2 = -(ubarf(ilevel,ipatch)%fld(2,i,j,k,isub) &
                        & -uf (ilevel,ipatch)%fld(2,i,j,k,isub)) &
                        & * chif(ilevel,ipatch)%fld( i,j,k,isub)
              tmp3 = -(ubarf(ilevel,ipatch)%fld(3,i,j,k,isub) &
                        & -uf (ilevel,ipatch)%fld(3,i,j,k,isub)) &
                        & * chif(ilevel,ipatch)%fld( i,j,k,isub)
              !-------------------------------------------------------------------
              ! Integrate linear force
              !-------------------------------------------------------------------
              pforce(1) = pforce(1) + tmp
              pforce(2) = pforce(2) + tmp2
              pforce(3) = pforce(3) + tmp3
              !-------------------------------------------------------------------
              ! Integrate force moment
              !-------------------------------------------------------------------
              pmoment(1) = pmoment(1) + py*tmp3 - pz*tmp2
              pmoment(2) = pmoment(2) + pz*tmp - px*tmp3
              pmoment(3) = pmoment(3) + px*tmp2 - py*tmp
            ENDDO
          ENDDO
        ENDDO
      ENDDO
      !multiply with cell volume
      IF (penalization .EQ. 1) THEN
        !remembering to take time step average (/dtime)
        lforce = lforce + pforce *dx*dy*dz /dtime
        lmoment = lmoment + pmoment *dx*dy*dz /dtime
      ELSE
        lforce = lforce + pforce *dx*dy*dz
        lmoment = lmoment + pmoment *dx*dy*dz
      ENDIF
    ENDDO
  ENDDO
ENDIF
!----------------------------------------------------------------------------
! Check if a forced vti dump is required. And clean up files
!----------------------------------------------------------------------------
forcedumpvrt = .FALSE.
forcedumpvel = .FALSE.
forcedumpchi = .FALSE.
forcedumpcnc = .FALSE.
forcedumpdwp = .FALSE.
forceabort = .FALSE.
IF (rank .EQ. 0) THEN
  WRITE(filename,'(A)') 'DUMPVRT'
  INQUIRE(FILE=filename,EXIST=forcedumpvrt)
  IF (forcedumpvrt) THEN
    OPEN(10,FILE=filename)
    CLOSE(10,STATUS='delete')
  ENDIF
  WRITE(filename,'(A)') 'DUMPVEL'
  INQUIRE(FILE=filename,EXIST=forcedumpvel)
  IF (forcedumpvel) THEN
    OPEN(10,FILE=filename)
    CLOSE(10,STATUS='delete')
  ENDIF
  WRITE(filename,'(A)') 'DUMPMSK'
  INQUIRE(FILE=filename,EXIST=forcedumpchi)
  IF (forcedumpchi) THEN
    OPEN(10,FILE=filename)
    CLOSE(10,STATUS='delete')
  ENDIF
  WRITE(filename,'(A)') 'DUMPCNC'
  INQUIRE(FILE=filename,EXIST=forcedumpcnc)
  IF (forcedumpcnc) THEN
    OPEN(10,FILE=filename)
    CLOSE(10,STATUS='delete')
  ENDIF
  WRITE(filename,'(A)') 'DUMPRHS'
  INQUIRE(FILE=filename,EXIST=forcedumpdwp)
  IF (forcedumpdwp) THEN
    OPEN(10,FILE=filename)
    CLOSE(10,STATUS='delete')
  ENDIF
  WRITE(filename,'(A)') 'ABORT'
  INQUIRE(FILE=filename,EXIST=forceabort)
ENDIF
!----------------------------------------------------------------------------
! Get global values
!----------------------------------------------------------------------------
CALL MPI_Reduce(ltotvort ,gtotvort ,3, mpi_prec,MPI_SUM,0,comm,info)
CALL MPI_Reduce(lkinenrg ,gkinenrg ,3, mpi_prec,MPI_SUM,0,comm,info)
CALL MPI_Reduce(lforce ,gforce ,3, mpi_prec,MPI_SUM,0,comm,info)
CALL MPI_Reduce(lmoment ,gmoment ,3, mpi_prec,MPI_SUM,0,comm,info)
CALL MPI_Reduce(ltotabsvort,gtotabsvort,1, mpi_prec,MPI_SUM,0,comm,info)
CALL MPI_Reduce(lenstrophy ,genstrophy ,1, mpi_prec,MPI_SUM,0,comm,info)
CALL MPI_Reduce(lmaxabsvel ,gmaxabsvel ,1, mpi_prec,MPI_MAX,0,comm,info)
CALL MPI_Reduce(lintdivvel ,gintdivvel ,1, mpi_prec,MPI_SUM,0,comm,info)
CALL MPI_Reduce(lintdivvort,gintdivvort,1, mpi_prec,MPI_SUM,0,comm,info)
CALL MPI_Reduce(lmaxdivvel ,gmaxdivvel ,1, mpi_prec,MPI_MAX,0,comm,info)
CALL MPI_Reduce(lmaxdivvort,gmaxdivvort,1, mpi_prec,MPI_MAX,0,comm,info)
CALL MPI_Reduce(lmaxstrain ,gmaxstrain ,1, mpi_prec,MPI_MAX,0,comm,info)
CALL MPI_Reduce(lfourier ,gfourier ,1, mpi_prec,MPI_MAX,0,comm,info)
CALL MPI_Reduce(lcourant ,gcourant ,1, mpi_prec,MPI_MAX,0,comm,info)
CALL MPI_Reduce(linterot ,ginterot ,1, mpi_prec,MPI_SUM,0,comm,info)
CALL MPI_Bcast(forcedumpvrt,1,MPI_LOGICAL,0,comm,info)
CALL MPI_Bcast(forcedumpvel,1,MPI_LOGICAL,0,comm,info)
CALL MPI_Bcast(forcedumpchi,1,MPI_LOGICAL,0,comm,info)
CALL MPI_Bcast(forcedumpcnc,1,MPI_LOGICAL,0,comm,info)
CALL MPI_Bcast(forcedumpdwp,1,MPI_LOGICAL,0,comm,info)
CALL MPI_Bcast(forceabort ,1,MPI_LOGICAL,0,comm,info)
!----------------------------------------------------------------------------
! Postprocessing
!----------------------------------------------------------------------------
IF (rank .EQ. 0) THEN
  IF (itime .EQ. 0) THEN
    gkinenrgdt = 0.0_MK
    gvisceff = 1.0_MK
    kinenrgold = gkinenrg(1)+gkinenrg(2)+gkinenrg(3)
  ELSE
    gkinenrgdt = (gkinenrg(1)+gkinenrg(2)+gkinenrg(3) - kinenrgold)/dtime
    gvisceff = -gkinenrgdt/(genstrophy*nu) !normalised by molec. kin. visc.
    kinenrgold = gkinenrg(1)+gkinenrg(2)+gkinenrg(3)
  ENDIF
ENDIF
!----------------------------------------------------------------------------
! Adaptive time-stepping goes here
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
! If the penalization parameter is adapted to the time step calculate it here
! ... and multiply the integrated solid force with the penalizationparameter
! but only in the case of explicit penalization
!----------------------------------------------------------------------------
IF (penalization .EQ. 2 .AND. penalizationadapt) THEN
  penalizationstrength = penalizationparam/dtime
ELSE
  penalizationstrength = penalizationparam
ENDIF
IF (penalization .EQ. 2) THEN
  gforce = gforce *penalizationstrength
  gmoment = gmoment *penalizationstrength
ENDIF
!----------------------------------------------------------------------------
! Maybe add warning here if fourier number etc is too high
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
! Write to file (only root)
!----------------------------------------------------------------------------
IF (rank .EQ. 0) THEN
  WRITE(filename,'(A,A)') runtag(1:iruntag), '.diag'
  OPEN(10,FILE=filename,IOSTAT=ios,POSITION='append',STATUS='unknown')
  IF (ios .NE. 0) THEN
    CALL naga_say(rank, 'naga_diagnostics','Failed to open output file.')
    GOTO 9999
  ENDIF
  !initialise output char arrays
  IF (itime .EQ. 0) THEN
    firsttime = .TRUE.
    CALL naga_output(.TRUE.,info)
  ELSE
    firsttime = .FALSE.
    CALL naga_output(.FALSE.)
  ENDIF
  !add output
  call naga_output('itime' ,itime) !1
  call naga_output('time' ,time)
  call naga_output('dtime' ,dtime)
  call naga_output('max_abs_vel' ,gmaxabsvel)
  call naga_output('kinetic_energy_x' ,gkinenrg(1)) !5
  call naga_output('kinetic_energy_y' ,gkinenrg(2))
  call naga_output('kinetic_energy_z' ,gkinenrg(3))
  call naga_output('tot_kin_enrg/dt' ,gkinenrgdt)
  call naga_output('eff_visc(norm.)' ,gvisceff)
  CALL naga_output('int_vorticity_x' ,gtotvort(1)) !10
  CALL naga_output('int_vorticity_y' ,gtotvort(2))
  CALL naga_output('int_vorticity_z' ,gtotvort(3))
  CALL naga_output('int_abs_vorticity' ,gtotabsvort)
  CALL naga_output('enstrophy' ,genstrophy)
  CALL naga_output('Fourier_number' ,gfourier) !15
  CALL naga_output('CFL_number' ,gcourant)
  CALL naga_output('int_omega-curl(vel)' ,ginterot)
  CALL naga_output('max_vel_divergence' ,gmaxdivvel)
  CALL naga_output('int_vel_divergence' ,gintdivvel)
  CALL naga_output('max_vort_divergence' ,gmaxdivvort)
  CALL naga_output('int_vort_divergence' ,gintdivvort) !20
  CALL naga_output('max_strain_rate' ,gmaxstrain)
  CALL naga_output('penalization_param' ,penalizationstrength)
  CALL naga_output('aerodyn_force_x' ,gforce(1))
  CALL naga_output('aerodyn_force_y' ,gforce(2))
  CALL naga_output('aerodyn_force_z' ,gforce(3))
  CALL naga_output('aerodyn_moment_x' ,gmoment(1))
  CALL naga_output('aerodyn_moment_y' ,gmoment(2)) !25
  CALL naga_output('aerodyn_moment_z' ,gmoment(3))
  !write data to file
  CALL naga_output(10)
  CLOSE(10)
END IF
!----------------------------------------------------------------------------
! Return
!----------------------------------------------------------------------------
 9999 CONTINUE
CALL substop('naga_diagnostics',t0,info)
RETURN
END SUBROUTINE naga_diagnostics
END MODULE naga_mod_diagnostics
