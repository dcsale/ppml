!------------------------------------------------------------------------------
! Subroutine :  naga_validation.f90
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
! This routine compares the numerical solution with an analytic reference and
! outputs error statistics to stdout
!------------------------------------------------------------------------------

MODULE naga_mod_validation

IMPLICIT NONE

INTERFACE naga_validation
  MODULE PROCEDURE naga_validation
END INTERFACE

CONTAINS

SUBROUTINE naga_validation(info)

USE naga_mod_globals
USE naga_mod_say
USE naga_mod_dump_vtk
USE naga_mod_save_fields

IMPLICIT NONE
!-------------------------------------------------------------------------
! Arguments
!-------------------------------------------------------------------------
INTEGER, INTENT(OUT) :: info


!-------------------------------------------------------------------------
! Local variables
!-------------------------------------------------------------------------
REAL(MK)             :: t0
INTEGER              :: i,j,k
INTEGER              :: ilevel,ipatch,isubl,isub
CHARACTER(LEN=256)   :: msg
REAL(MK)             ::  rmserr, rmsref, maxerr, maxref,npoints
REAL(MK)             :: grmserr,grmsref,gmaxerr,gmaxref


!----------------------------------------------------------------------------
! Initialise routine
!----------------------------------------------------------------------------
CALL substart('naga_validation',t0,info)


!----------------------------------------------------------------------------
! Loop through the sub domain points
!----------------------------------------------------------------------------
IF (.NOT. (flowcase .EQ. 303 .OR. flowcase .EQ. 304)) THEN
  DO ilevel=1,nlevels
    DO ipatch=1,npatches(ilevel)
      rmserr = 0.0_MK
      rmsref = 0.0_MK
      maxerr = 0.0_MK
      maxref = 0.0_MK
      DO isub=1,topos(ilevel)%nsublist
        wf(ilevel,ipatch)%fld = 0.0_MK
        isubl=topos(ilevel)%isublist(isub)
        DO k=1,ptcset(ilevel,ipatch)%snx(3,isubl)
          DO j=1,ptcset(ilevel,ipatch)%snx(2,isubl)
            DO i=1,ptcset(ilevel,ipatch)%snx(1,isubl)
              wf    (ilevel,ipatch)%fld(1,i,j,k,isub) = &
              & uf  (ilevel,ipatch)%fld(1,i,j,k,isub) - &
              & psif(ilevel,ipatch)%fld(1,i,j,k,isub)
              wf    (ilevel,ipatch)%fld(2,i,j,k,isub) = &
              & uf  (ilevel,ipatch)%fld(2,i,j,k,isub) - &
              & psif(ilevel,ipatch)%fld(2,i,j,k,isub)
              wf    (ilevel,ipatch)%fld(3,i,j,k,isub) = &
              & uf  (ilevel,ipatch)%fld(3,i,j,k,isub) - &
              & psif(ilevel,ipatch)%fld(3,i,j,k,isub)
              rmserr = rmserr + wf(ilevel,ipatch)%fld(1,i,j,k,isub)**2
              rmserr = rmserr + wf(ilevel,ipatch)%fld(2,i,j,k,isub)**2
              rmserr = rmserr + wf(ilevel,ipatch)%fld(3,i,j,k,isub)**2
              rmsref = rmsref + psif(ilevel,ipatch)%fld(1,i,j,k,isub)**2
              rmsref = rmsref + psif(ilevel,ipatch)%fld(2,i,j,k,isub)**2
              rmsref = rmsref + psif(ilevel,ipatch)%fld(3,i,j,k,isub)**2
              maxerr = MAX(maxerr,abs(wf(ilevel,ipatch)%fld(1,i,j,k,isub)))
              maxerr = MAX(maxerr,abs(wf(ilevel,ipatch)%fld(2,i,j,k,isub)))
              maxerr = MAX(maxerr,abs(wf(ilevel,ipatch)%fld(3,i,j,k,isub)))
              maxref = MAX(maxref,abs(psif(ilevel,ipatch)%fld(1,i,j,k,isub)))
              maxref = MAX(maxref,abs(psif(ilevel,ipatch)%fld(2,i,j,k,isub)))
              maxref = MAX(maxref,abs(psif(ilevel,ipatch)%fld(3,i,j,k,isub)))
            ENDDO
          ENDDO
        ENDDO
      ENDDO
      IF (domainbc .EQ. 0) THEN
        npoints = REAL((ptcset(ilevel,ipatch)%gnx(1)* &
                 & ptcset(ilevel,ipatch)%gnx(2)* &
                 & ptcset(ilevel,ipatch)%gnx(3)),MK)
      ELSE
        npoints = REAL(((ptcset(ilevel,ipatch)%gnx(1)-1)* &
                 & (ptcset(ilevel,ipatch)%gnx(2)-1)* &
                 & (ptcset(ilevel,ipatch)%gnx(3)-1)),MK)
      ENDIF

      !MPI away
      CALL MPI_Reduce(rmserr,grmserr,1, mpi_prec,MPI_SUM,0,comm,info)
      CALL MPI_Reduce(rmsref,grmsref,1, mpi_prec,MPI_SUM,0,comm,info)
      CALL MPI_Reduce(maxerr,gmaxerr,1, mpi_prec,MPI_MAX,0,comm,info)
      CALL MPI_Reduce(maxref,gmaxref,1, mpi_prec,MPI_MAX,0,comm,info)

      IF (rank .EQ. 0) THEN
        IF (domainbc .EQ. 0) THEN
          WRITE(*,*) 'errorNX ', ptcset(ilevel,ipatch)%gnx(1)
        ELSE
          WRITE(*,*) 'errorNX ', ptcset(ilevel,ipatch)%gnx(1)-1
        ENDIF
        WRITE(*,*) 'errorMAX ',gmaxerr
        WRITE(*,*) 'errorRMS ',SQRT(grmserr/npoints)
        WRITE(*,*) 'errorREL ',SQRT(grmserr)/SQRT(grmsref)
        WRITE(*,*) 'errorL2  ',SQRT(grmserr)
        WRITE(*,*) 'errorREF ',gmaxref
      ENDIF
    ENDDO
  ENDDO

  DO ilevel=1,nlevels
    DO ipatch=1,npatches(ilevel)
      msg = 'valdiff'
      CALL naga_dump_vtk(ptcset(ilevel,ipatch)%topoid, ptcset(ilevel,ipatch)%meshid, wf(ilevel,ipatch)%fld,msg,0,info)
      IF (info .NE. 0) THEN
        CALL naga_say(rank, 'naga_validation','Failed to save field_w as vtk.')
        GOTO 9999
      ENDIF

      msg = 'valcomp'
      CALL naga_dump_vtk(ptcset(ilevel,ipatch)%topoid, ptcset(ilevel,ipatch)%meshid, uf(ilevel,ipatch)%fld,msg,0,info)
      IF (info .NE. 0) THEN
        CALL naga_say(rank, 'naga_validation','Failed to save field_w as vtk.')
        GOTO 9999
      ENDIF

      msg = 'valref'
      CALL naga_dump_vtk(ptcset(ilevel,ipatch)%topoid, ptcset(ilevel,ipatch)%meshid, psif(ilevel,ipatch)%fld,msg,0,info)
      IF (info .NE. 0) THEN
        CALL naga_say(rank, 'naga_validation','Failed to save field_w as vtk.')
        GOTO 9999
      ENDIF
    ENDDO
  ENDDO
ENDIF

IF (flowcase .EQ. 303 .OR. flowcase .EQ. 304) THEN
  DO ilevel=1,nlevels
    DO ipatch=1,npatches(ilevel)
      rmserr = 0.0_MK
      rmsref = 0.0_MK
      maxerr = 0.0_MK
      maxref = 0.0_MK
      DO isub=1,topos(ilevel)%nsublist
        isubl=topos(ilevel)%isublist(isub)
        DO k=1,ptcset(ilevel,ipatch)%snx(3,isubl)
          DO j=1,ptcset(ilevel,ipatch)%snx(2,isubl)
            DO i=1,ptcset(ilevel,ipatch)%snx(1,isubl)
              uf    (ilevel,ipatch)%fld(1,i,j,k,isub) = &
              & wf  (ilevel,ipatch)%fld(1,i,j,k,isub) - &
              & psif(ilevel,ipatch)%fld(1,i,j,k,isub)
              uf    (ilevel,ipatch)%fld(2,i,j,k,isub) = &
              & wf  (ilevel,ipatch)%fld(2,i,j,k,isub) - &
              & psif(ilevel,ipatch)%fld(2,i,j,k,isub)
              uf    (ilevel,ipatch)%fld(3,i,j,k,isub) = &
              & wf  (ilevel,ipatch)%fld(3,i,j,k,isub) - &
              & psif(ilevel,ipatch)%fld(3,i,j,k,isub)
              rmserr = rmserr + uf(ilevel,ipatch)%fld(1,i,j,k,isub)**2
              rmserr = rmserr + uf(ilevel,ipatch)%fld(2,i,j,k,isub)**2
              rmserr = rmserr + uf(ilevel,ipatch)%fld(3,i,j,k,isub)**2
              rmsref = rmsref + psif(ilevel,ipatch)%fld(1,i,j,k,isub)**2
              rmsref = rmsref + psif(ilevel,ipatch)%fld(2,i,j,k,isub)**2
              rmsref = rmsref + psif(ilevel,ipatch)%fld(3,i,j,k,isub)**2
              maxerr = MAX(maxerr,abs(uf(ilevel,ipatch)%fld(1,i,j,k,isub)))
              maxerr = MAX(maxerr,abs(uf(ilevel,ipatch)%fld(2,i,j,k,isub)))
              maxerr = MAX(maxerr,abs(uf(ilevel,ipatch)%fld(3,i,j,k,isub)))
              maxref = MAX(maxref,abs(psif(ilevel,ipatch)%fld(1,i,j,k,isub)))
              maxref = MAX(maxref,abs(psif(ilevel,ipatch)%fld(2,i,j,k,isub)))
              maxref = MAX(maxref,abs(psif(ilevel,ipatch)%fld(3,i,j,k,isub)))
            ENDDO
          ENDDO
        ENDDO
      ENDDO
      IF (domainbc .EQ. 0) THEN
        npoints = (ptcset(ilevel,ipatch)%gnx(1)* &
                 & ptcset(ilevel,ipatch)%gnx(2)* &
                 & ptcset(ilevel,ipatch)%gnx(3))
      ELSE
        npoints = ((ptcset(ilevel,ipatch)%gnx(1)-1)* &
                 & (ptcset(ilevel,ipatch)%gnx(2)-1)* &
                 & (ptcset(ilevel,ipatch)%gnx(3)-1))
      ENDIF

      !MPI away
      CALL MPI_Reduce(rmserr,grmserr,1, mpi_prec,MPI_SUM,0,comm,info)
      CALL MPI_Reduce(rmsref,grmsref,1, mpi_prec,MPI_SUM,0,comm,info)
      CALL MPI_Reduce(maxerr,gmaxerr,1, mpi_prec,MPI_MAX,0,comm,info)
      CALL MPI_Reduce(maxref,gmaxref,1, mpi_prec,MPI_MAX,0,comm,info)

      IF (rank .EQ. 0) THEN
        IF (domainbc .EQ. 0) THEN
          WRITE(*,*) 'errorNX ', ptcset(ilevel,ipatch)%gnx(1)
        ELSE
          WRITE(*,*) 'errorNX ', ptcset(ilevel,ipatch)%gnx(1)-1
        ENDIF
        WRITE(*,*) 'errorMAX ',maxerr
        WRITE(*,*) 'errorRMS ',SQRT(rmserr/npoints)
        WRITE(*,*) 'errorREL ',SQRT(rmserr)/SQRT(rmsref)
        WRITE(*,*) 'errorL2  ',SQRT(rmserr)
        WRITE(*,*) 'errorREF ',maxref
      ENDIF
    ENDDO
  ENDDO

  DO ilevel=1,nlevels
    DO ipatch=1,npatches(ilevel)
      msg = 'vrtref'
      CALL naga_dump_vtk(ptcset(ilevel,ipatch)%topoid, ptcset(ilevel,ipatch)%meshid, psif(ilevel,ipatch)%fld,msg,0,info)
      IF (info .NE. 0) THEN
        CALL naga_say(rank, 'naga_validation','Failed to save field_psi as vtk.')
        GOTO 9999
      ENDIF

      msg = 'vrtrepr'
      CALL naga_dump_vtk(ptcset(ilevel,ipatch)%topoid, ptcset(ilevel,ipatch)%meshid, wf(ilevel,ipatch)%fld,msg,0,info)
      IF (info .NE. 0) THEN
        CALL naga_say(rank, 'naga_validation','Failed to save field_w as vtk.')
        GOTO 9999
      ENDIF

      msg = 'vrtdiff'
      CALL naga_dump_vtk(ptcset(ilevel,ipatch)%topoid, ptcset(ilevel,ipatch)%meshid, uf(ilevel,ipatch)%fld,msg,0,info)
      IF (info .NE. 0) THEN
        CALL naga_say(rank, 'naga_validation','Failed to save field_u as vtk.')
        GOTO 9999
      ENDIF
    ENDDO
  ENDDO
END IF


!----------------------------------------------------------------------------
! Return
!----------------------------------------------------------------------------
 9999 CONTINUE
CALL substop('naga_validation',t0,info)
RETURN


END SUBROUTINE naga_validation

END MODULE naga_mod_validation

