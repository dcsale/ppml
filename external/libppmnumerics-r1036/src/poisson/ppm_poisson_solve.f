      !-------------------------------------------------------------------------
      !  Subroutine   : ppm_poisson_solve.f90
      !-------------------------------------------------------------------------
      ! Copyright (c) 2010 CSE Lab (ETH Zurich), MOSAIC Group (ETH Zurich),
      !                    Center for Fluid Dynamics (DTU)
      !
      !-------------------------------------------------------------------------
      !@mapping calls can be cleaned up
      !-------------------------------------------------------------------------
      SUBROUTINE __ROUTINE(topoid,meshid,ppmpoisson,fieldin,fieldout,gstw,info,&
                         & tmpcase)
      !!! Routine to perform the Greens function solution of the Poisson
      !!! equation. All settings are defined in ppm_poisson_initdef and stored 
      !!! in the ppmpoisson plan. The tmpcase argument allows the use of a
      !!! different Greens function or operation than initialised. This is 
      !!! particularly useful for vorticity reprojection 
      !!! (ppm_poisson_grn_reprojec).
      !!! 
      !!! The code will eventually get a needed cleanup overhaul.

      USE ppm_module_map_field
      USE ppm_module_map_field_global
      USE ppm_module_map

      IMPLICIT NONE
#ifdef __MPI
      include 'mpif.h'
#endif
      !-------------------------------------------------------------------------
      ! Arguments
      !-------------------------------------------------------------------------
      INTEGER, INTENT(IN)                                         :: topoid
      !!! Topology ID
      INTEGER, INTENT(IN)                                         :: meshid
      !!! Mesh ID
      TYPE(ppm_poisson_plan),INTENT(INOUT)                        :: ppmpoisson
      !!! The PPM Poisson plan
      REAL(__PREC),DIMENSION(:,:,:,:,:),POINTER                   :: fieldin
      !!! Input data field
      REAL(__PREC),DIMENSION(:,:,:,:,:),POINTER                   :: fieldout
      !!! Output data field
      INTEGER,DIMENSION(__DIM),INTENT(IN)                         :: gstw
      !!! Ghost layer width
      INTEGER, INTENT(OUT)                                        :: info
      !!! Return status, 0 upon succes
      INTEGER,OPTIONAL,INTENT(IN)                                 :: tmpcase
      !!! Temporary operation (useful for ppm_poisson_grn_reprojec)

      !-------------------------------------------------------------------------
      ! Local variables
      !-------------------------------------------------------------------------
      INTEGER,PARAMETER                 :: MK = __PREC
      REAL(__PREC)                      :: t0
      INTEGER                           :: isub,isubl
      INTEGER                           :: i,j,k
      INTEGER                           :: info2
      INTEGER                           :: presentcase
      COMPLEX(__PREC)                   :: helmholtzpic
      INTEGER                           :: gi,gj,gk
      COMPLEX(__PREC)                   :: kx,ky,kz
      COMPLEX(__PREC)                   :: phix,phiy,phiz
      REAL(__PREC)                      :: normfac

      !-------------------------------------------------------------------------
      ! Initialise routine
      !-------------------------------------------------------------------------
      CALL substart('ppm_poisson_solve',t0,info)

      !-------------------------------------------------------------------------
      ! Check if we run a different/temporary case
      !-------------------------------------------------------------------------
      IF (PRESENT(tmpcase)) THEN
         presentcase = tmpcase
      ELSE
         presentcase = ppmpoisson%case
      ENDIF

      !-------------------------------------------------------------------------
      !@ Perhaps check if ghostlayer suffices for a given fd stencil
      !-------------------------------------------------------------------------


      !-------------------------------------------------------------------------
      ! Perhaps allocate (and deallocate) arrays
      !-------------------------------------------------------------------------


      !-----------------------------------------------------------------------
      ! Set the real xy slabs 0 (part of the 0 padding) for free-space
      !@ free-space calculations and reprojection may cause problems !why?
      !-----------------------------------------------------------------------
      IF (presentcase .EQ. ppm_poisson_grn_pois_fre) THEN
        ppmpoisson%fldxyr = 0.0_MK
      ENDIF
      !-----------------------------------------------------------------------
      ! Map data globally to the slabs (XY)
      ! This is where the vorticity is extended and padded with 0 for free-space
      !-----------------------------------------------------------------------
      !Initialise
      CALL ppm_map_field_global(&
      & topoid, &
      & ppmpoisson%topoidxy, &
      & meshid, &
      & ppmpoisson%meshidxy,info)
      IF (info .NE. 0) THEN
        CALL ppm_write(ppm_rank,'ppm_poisson_solve','Failed to initialise field mapping.',info2)
        GOTO 9999
      ENDIF

      !Push the data
      CALL ppm_map_field_push(&
      & topoid, &
      & meshid,fieldin,__NCOM,info)
      IF (info .NE. 0) THEN
        CALL ppm_write(ppm_rank, 'ppm_poisson_solve','Failed to push vector field.',info2)
        GOTO 9999
      ENDIF

      !Send
      CALL ppm_map_field_send(info)
      IF (info .NE. 0) THEN
        CALL ppm_write(ppm_rank, 'ppm_poisson_solve','Failed to send field.',info2)
        GOTO 9999
      ENDIF

      !Retrieve
      CALL ppm_map_field_pop(&
      & ppmpoisson%topoidxy, &
      & ppmpoisson%meshidxy,ppmpoisson%fldxyr, &
      & __NCOM,__ZEROSI,info)
      IF (info .NE. 0) THEN
        CALL ppm_write(ppm_rank, 'ppm_poisson_solve','Failed to pop vector field.',info2)
        GOTO 9999
      ENDIF

      !-----------------------------------------------------------------------
      ! Do slab FFT (XY) - use the non-reduced topology !@what does this mean
      !-----------------------------------------------------------------------
      CALL ppm_fft_execute_2d(ppmpoisson%topoidxy,&
      & ppmpoisson%meshidxy, ppmpoisson%planfxy, &
      & ppmpoisson%fldxyr, ppmpoisson%fldxyc, &
      & info)


      !-----------------------------------------------------------------------
      ! Map to the pencils (Z)
      !-----------------------------------------------------------------------
      !Initialise
      CALL ppm_map_field_global(&
      & ppmpoisson%topoidxy, &
      & ppmpoisson%topoidz, &
      & ppmpoisson%meshidxyc, & !@to be used when meshid>1 works
      !& ppmpoisson%meshidxy, &
      & ppmpoisson%meshidz,info)
      IF (info .NE. 0) THEN
        CALL ppm_write(ppm_rank, 'ppm_poisson_solve','Failed to initialise field mapping.',info2)
        GOTO 9999
      ENDIF

      !Push the data
      CALL ppm_map_field_push(&
      & ppmpoisson%topoidxy, &
      & ppmpoisson%meshidxyc,ppmpoisson%fldxyc,__NCOM,info)
      !& ppmpoisson%meshidxyc,ppmpoisson%fldxyc,__NCOM,info)!@to be used when meshid>1 works
      IF (info .NE. 0) THEN
        CALL ppm_write(ppm_rank, 'ppm_poisson_solve','Failed to push vector field.',info2)
        GOTO 9999
      ENDIF

      !Send
      CALL ppm_map_field_send(info)
      IF (info .NE. 0) THEN
        CALL ppm_write(ppm_rank, 'ppm_poisson_solve','Failed to send field.',info2)
        GOTO 9999
      ENDIF

      !Retrieve
      CALL ppm_map_field_pop(&
      & ppmpoisson%topoidz, &
      & ppmpoisson%meshidz,ppmpoisson%fldzc1, &
      & __NCOM,__ZEROSI,info)
      IF (info .NE. 0) THEN
        CALL ppm_write(ppm_rank, 'ppm_poisson_solve','Failed to pop vector field.',info2)
        GOTO 9999
      ENDIF

      !-----------------------------------------------------------------------
      ! Do pencil FFT (Z)
      !-----------------------------------------------------------------------
      CALL ppm_fft_execute_1d(ppmpoisson%topoidz,&
      & ppmpoisson%meshidz, ppmpoisson%planfz, &
      & ppmpoisson%fldzc1, ppmpoisson%fldzc2, &
      & info)


      !-----------------------------------------------------------------------
      ! Apply the periodic Greens function
      !-----------------------------------------------------------------------
      IF (presentcase .EQ. ppm_poisson_grn_pois_per) THEN
        DO isub=1,ppmpoisson%nsublistz
          isubl=ppmpoisson%isublistz(isub)
          DO k=1,ppmpoisson%ndataz(3,isubl)
            DO j=1,ppmpoisson%ndataz(2,isubl)
              DO i=1,ppmpoisson%ndataz(1,isubl)
                ppmpoisson%fldzc2(1,i,j,k,isub) = ppmpoisson%fldgrnr( i,j,k,isub)*&
                                                & ppmpoisson%fldzc2(1,i,j,k,isub)
                ppmpoisson%fldzc2(2,i,j,k,isub) = ppmpoisson%fldgrnr( i,j,k,isub)*&
                                                & ppmpoisson%fldzc2(2,i,j,k,isub)
                ppmpoisson%fldzc2(3,i,j,k,isub) = ppmpoisson%fldgrnr( i,j,k,isub)*&
                                                & ppmpoisson%fldzc2(3,i,j,k,isub)
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      !-----------------------------------------------------------------------
      ! Apply the free-space Greens function
      !-----------------------------------------------------------------------
      ELSE IF (presentcase .EQ. ppm_poisson_grn_pois_fre) THEN
        DO isub=1,ppmpoisson%nsublistz
          isubl=ppmpoisson%isublistz(isub)
          DO k=1,ppmpoisson%ndataz(3,isubl)
            DO j=1,ppmpoisson%ndataz(2,isubl)
              DO i=1,ppmpoisson%ndataz(1,isubl)
                ppmpoisson%fldzc2(1,i,j,k,isub) = ppmpoisson%fldgrnc( i,j,k,isub)*&
                                                & ppmpoisson%fldzc2(1,i,j,k,isub)
                ppmpoisson%fldzc2(2,i,j,k,isub) = ppmpoisson%fldgrnc( i,j,k,isub)*&
                                                & ppmpoisson%fldzc2(2,i,j,k,isub)
                ppmpoisson%fldzc2(3,i,j,k,isub) = ppmpoisson%fldgrnc( i,j,k,isub)*&
                                                & ppmpoisson%fldzc2(3,i,j,k,isub)
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      !-----------------------------------------------------------------------
      ! Vorticity re-projection
      !-----------------------------------------------------------------------
      ELSE IF (presentcase .EQ. ppm_poisson_grn_reprojec) THEN
        !remembering to normalize the FFT
        IF (ppmpoisson%case .EQ. ppm_poisson_grn_pois_fre) THEN
          normfac = 1.0_MK/ REAL((ppmpoisson%nmz(1))* & !vertex
                               & (ppmpoisson%nmz(2))* &
                               & (ppmpoisson%nmz(3)),MK)
        ELSE IF (ppmpoisson%case .EQ. ppm_poisson_grn_pois_per) THEN
          normfac = 1.0_MK/ REAL((ppmpoisson%nmz(1)-1)* & !vertex
                               & (ppmpoisson%nmz(2)-1)* &
                               & (ppmpoisson%nmz(3)-1),MK)
        ENDIF
        DO isub=1,ppmpoisson%nsublistz
          isubl=ppmpoisson%isublistz(isub)
          DO k=1,ppmpoisson%ndataz(3,isubl)
            gk = k - 1 + (ppmpoisson%istartz(3,isubl)-1)
            IF (gk .GT. (ppmpoisson%nmz(3)-1)/2) gk = gk-(ppmpoisson%nmz(3)-1)
            kz = REAL(gk,MK)*ppmpoisson%normkz
            DO j=1,ppmpoisson%ndataz(2,isubl)
              gj = j - 1 + (ppmpoisson%istartz(2,isubl)-1)
              IF (gj .GT. (ppmpoisson%nmz(2)-1)/2) gj = gj-(ppmpoisson%nmz(2)-1)
              ky = REAL(gj,MK)*ppmpoisson%normky
              DO i=1,ppmpoisson%ndataz(1,isubl)
                gi = i - 1 + (ppmpoisson%istartz(1,isubl)-1)
                IF (gi .GT. (ppmpoisson%nmz(1)-1)/2) gi = gi-(ppmpoisson%nmz(1)-1)
                kx = REAL(gi,MK)*ppmpoisson%normkx

                !IF (gi .EQ. 0 .AND. gj .EQ. 0 .AND. gk .EQ. 0) THEN
                  !wdotk = 0.0_mk
                !ELSE
                  !wdotk = (ppmpoisson%fldzc2(1,i,j,k,isub) * kx +  &
                        !&  ppmpoisson%fldzc2(2,i,j,k,isub) * ky +  &
                        !&  ppmpoisson%fldzc2(3,i,j,k,isub) * kz) / &
                        !&  (kx*kx+ky*ky+kz*kz)
                !ENDIF

                !ppmpoisson%fldzc2(1,i,j,k,isub) = &
                  !& (ppmpoisson%fldzc2(1,i,j,k,isub) - wdotk*kx)*normfac
                !ppmpoisson%fldzc2(2,i,j,k,isub) = &
                  !& (ppmpoisson%fldzc2(2,i,j,k,isub) - wdotk*ky)*normfac
                !ppmpoisson%fldzc2(3,i,j,k,isub) = &
                  !& (ppmpoisson%fldzc2(3,i,j,k,isub) - wdotk*kz)*normfac

                IF (ppmpoisson%case .EQ. ppm_poisson_grn_pois_fre) THEN
                  helmholtzpic = (ppmpoisson%fldzc2(1,i,j,k,isub) * kx +  &
                               &  ppmpoisson%fldzc2(2,i,j,k,isub) * ky +  &
                               &  ppmpoisson%fldzc2(3,i,j,k,isub) * kz) * &
                               &  ppmpoisson%fldgrnc( i,j,k,isub)
                      ppmpoisson%fldzc2(1,i,j,k,isub) = & !@formatting
                   & (ppmpoisson%fldzc2(1,i,j,k,isub)*normfac - helmholtzpic*kx)
                      ppmpoisson%fldzc2(2,i,j,k,isub) = &
                   & (ppmpoisson%fldzc2(2,i,j,k,isub)*normfac - helmholtzpic*ky)
                      ppmpoisson%fldzc2(3,i,j,k,isub) = &
                   & (ppmpoisson%fldzc2(3,i,j,k,isub)*normfac - helmholtzpic*kz)
                ELSE IF (ppmpoisson%case .EQ. ppm_poisson_grn_pois_per) THEN
                  helmholtzpic = (ppmpoisson%fldzc2(1,i,j,k,isub) * kx +  &
                               &  ppmpoisson%fldzc2(2,i,j,k,isub) * ky +  &
                               &  ppmpoisson%fldzc2(3,i,j,k,isub) * kz) * &
                               &  ppmpoisson%fldgrnr( i,j,k,isub)
                      ppmpoisson%fldzc2(1,i,j,k,isub) = & !@formatting
                   & (ppmpoisson%fldzc2(1,i,j,k,isub)*normfac - helmholtzpic*kx)
                      ppmpoisson%fldzc2(2,i,j,k,isub) = &
                   & (ppmpoisson%fldzc2(2,i,j,k,isub)*normfac - helmholtzpic*ky)
                      ppmpoisson%fldzc2(3,i,j,k,isub) = &
                   & (ppmpoisson%fldzc2(3,i,j,k,isub)*normfac - helmholtzpic*kz)
                ENDIF

              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ENDIF
      !-----------------------------------------------------------------------
      ! Spectral derivatives
      ! normkx, etc contains 2pi/Lx
      !-----------------------------------------------------------------------
      IF (ppmpoisson%derivatives .EQ. ppm_poisson_drv_curl_sp .AND.&
         & (presentcase .EQ. ppm_poisson_grn_pois_per .OR. &
            presentcase .EQ. ppm_poisson_grn_pois_fre)) THEN
        DO isub=1,ppmpoisson%nsublistz
          isubl=ppmpoisson%isublistz(isub)
          DO k=1,ppmpoisson%ndataz(3,isubl)
            gk = k - 1 + (ppmpoisson%istartz(3,isubl)-1)
            IF (gk .GT. (ppmpoisson%nmz(3)-1)/2) gk = gk-(ppmpoisson%nmz(3)-1)
            kz = CMPLX(0.0_MK,gk*ppmpoisson%normkz,MK)
            DO j=1,ppmpoisson%ndataz(2,isubl)
              gj = j - 1 + (ppmpoisson%istartz(2,isubl)-1)
              IF (gj .GT. (ppmpoisson%nmz(2)-1)/2) gj = gj-(ppmpoisson%nmz(2)-1)
              ky = CMPLX(0.0_MK,gj*ppmpoisson%normky,MK)
              DO i=1,ppmpoisson%ndataz(1,isubl)
                gi = i - 1 + (ppmpoisson%istartz(1,isubl)-1)
                IF (gi .GT. (ppmpoisson%nmz(1)-1)/2) gi = gi-(ppmpoisson%nmz(1)-1)
                kx = CMPLX(0.0_MK,gi*ppmpoisson%normkx,MK)

                phix = ppmpoisson%fldzc2(1,i,j,k,isub)
                phiy = ppmpoisson%fldzc2(2,i,j,k,isub)
                phiz = ppmpoisson%fldzc2(3,i,j,k,isub)

                ppmpoisson%fldzc2(1,i,j,k,isub) = (ky*phiz-kz*phiy)
                ppmpoisson%fldzc2(2,i,j,k,isub) = (kz*phix-kx*phiz)
                ppmpoisson%fldzc2(3,i,j,k,isub) = (kx*phiy-ky*phix)
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ENDIF


      !-----------------------------------------------------------------------
      ! IFFT pencil (Z)
      !-----------------------------------------------------------------------
      CALL ppm_fft_execute_1d(ppmpoisson%topoidz,&
      & ppmpoisson%meshidz, ppmpoisson%planbz, &
      & ppmpoisson%fldzc2, ppmpoisson%fldzc1, &
      & info)


      !-----------------------------------------------------------------------
      ! Map back to slabs (XY)
      !-----------------------------------------------------------------------
      !Initialise
      CALL ppm_map_field_global(&
      & ppmpoisson%topoidz, &
      & ppmpoisson%topoidxy, &
      & ppmpoisson%meshidz, &
      & ppmpoisson%meshidxyc,info)
      IF (info .NE. 0) THEN
        CALL ppm_write(ppm_rank, 'ppm_poisson_solve','Failed to initialise field mapping.',info2)
        GOTO 9999
      ENDIF

      !Push the data
      CALL ppm_map_field_push(&
      & ppmpoisson%topoidz, &
      & ppmpoisson%meshidz,ppmpoisson%fldzc1,__NCOM,info)
      IF (info .NE. 0) THEN
        CALL ppm_write(ppm_rank, 'ppm_poisson_solve','Failed to push vector field.',info2)
        GOTO 9999
      ENDIF

      !Send
      CALL ppm_map_field_send(info)
      IF (info .NE. 0) THEN
        CALL ppm_write(ppm_rank, 'ppm_poisson_solve','Failed to send field.',info2)
        GOTO 9999
      ENDIF

      !Retrieve
      CALL ppm_map_field_pop(&
      & ppmpoisson%topoidxy, &
      & ppmpoisson%meshidxyc,ppmpoisson%fldxyc, &
      & __NCOM,__ZEROSI,info)
      IF (info .NE. 0) THEN
        CALL ppm_write(ppm_rank, 'ppm_poisson_solve','Failed to pop vector field.',info2)
        GOTO 9999
      ENDIF

      !-----------------------------------------------------------------------
      ! IFFT (XY) use the non-reduced topology
      !-----------------------------------------------------------------------
      CALL ppm_fft_execute_2d(ppmpoisson%topoidxy,&
      & ppmpoisson%meshidxy, ppmpoisson%planbxy, &
      & ppmpoisson%fldxyc, ppmpoisson%fldxyr, &
      & info)


      !-----------------------------------------------------------------------
      ! Map back to standard topology (XYZ)
      !-----------------------------------------------------------------------
      !Initialise
      CALL ppm_map_field_global(&
      & ppmpoisson%topoidxy, &
      & topoid, &
      & ppmpoisson%meshidxy, &
      & meshid,info)
      IF (info .NE. 0) THEN
        CALL ppm_write(ppm_rank, 'ppm_poisson_solve','Failed to initialise field mapping.',info2)
        GOTO 9999
      ENDIF


      !Push the data
      CALL ppm_map_field_push(&
      & ppmpoisson%topoidxy, &
      & ppmpoisson%meshidxy,ppmpoisson%fldxyr,__NCOM,info)
      IF (info .NE. 0) THEN
        CALL ppm_write(ppm_rank, 'ppm_poisson_solve','Failed to push vector field.',info2)
        GOTO 9999
      ENDIF

      !Send
      CALL ppm_map_field_send(info)
      IF (info .NE. 0) THEN
        CALL ppm_write(ppm_rank, 'ppm_poisson_solve','Failed to send field.',info2)
        GOTO 9999
      ENDIF

      !-------------------------------------------------------------------------
      ! FINAL RETRIEVE - Here we do different things depending on the task
      ! i.e. the receiver varies
      !-------------------------------------------------------------------------
      IF ((ppmpoisson%derivatives .EQ. ppm_poisson_drv_curl_fd2 .OR. &
           ppmpoisson%derivatives .EQ. ppm_poisson_drv_curl_fd4) .AND. &
         & (presentcase .EQ. ppm_poisson_grn_pois_per .OR. &  !@these may be unnecessary - perhaps just the derive value. Or maybe not: in case of vorticity reprojection we could get lost
            presentcase .EQ. ppm_poisson_grn_pois_fre)) THEN
        CALL ppm_map_field_pop(&
        & topoid, &
        & meshid,ppmpoisson%drv_vr, &
        & __NCOM,gstw,info)
        IF (info .NE. 0) THEN
          CALL ppm_write(ppm_rank, 'ppm_poisson_solve','Failed to pop vector field.',info2)
          GOTO 9999
        ENDIF
        !-------------------------------------------------------------------------
        ! Ghost the temporary array for derivatives (drv_vr)
        !-------------------------------------------------------------------------
        CALL ppm_map_field_ghost_get(topoid,meshid,gstw,info)
        IF (info .NE. 0) THEN
           CALL ppm_write(ppm_rank, 'ppm_poisson_solve','Failed to initialise ghosts.',info2)
           GOTO 9999
        ENDIF
        CALL ppm_map_field_push(topoid,meshid,ppmpoisson%drv_vr,__NCOM,info)
        IF (info .NE. 0) THEN
           CALL ppm_write(ppm_rank, 'ppm_poisson_solve','Failed to push ghosts.',info2)
           GOTO 9999
        ENDIF
        CALL ppm_map_field_send(info)
        IF (info .NE. 0) THEN
           CALL ppm_write(ppm_rank, 'ppm_poisson_solve','Failed to send ghosts.',info2)
           GOTO 9999
        ENDIF
        CALL ppm_map_field_pop(topoid,meshid,ppmpoisson%drv_vr,__NCOM,gstw,info)
        IF (info .NE. 0) THEN
           CALL ppm_write(ppm_rank, 'ppm_poisson_solve','Failed to pop ghosts.',info2)
           GOTO 9999
        ENDIF

      ELSE
        CALL ppm_map_field_pop(&
        & topoid, &
        & meshid,fieldout, &
        & __NCOM,gstw,info)
      ENDIF
      IF (info .NE. 0) THEN
        CALL ppm_write(ppm_rank, 'ppm_poisson_solve','Failed to pop vector field.',info2)
        GOTO 9999
      ENDIF


      !-------------------------------------------------------------------------
      ! Treat ghost layer to make FD stencils work
      !-------------------------------------------------------------------------
      IF (ppmpoisson%derivatives .EQ. ppm_poisson_drv_curl_fd2 .AND.&
         & (presentcase .EQ. ppm_poisson_grn_pois_fre)) THEN
         CALL ppm_poisson_extrapolateghost(topoid,meshid,ppmpoisson%drv_vr,&
                                     & 2,4,gstw,info)
      ENDIF
      IF (ppmpoisson%derivatives .EQ. ppm_poisson_drv_curl_fd4 .AND.&
         & (presentcase .EQ. ppm_poisson_grn_pois_fre)) THEN
         CALL ppm_poisson_extrapolateghost(topoid,meshid,ppmpoisson%drv_vr,&
                                     & 2,4,gstw,info)
      ENDIF

      !-------------------------------------------------------------------------
      ! Optionally do derivatives
      ! Perhaps make ppm_poisson_fd take _none as argument. Then maybe no
      ! if-statement is required
      !-------------------------------------------------------------------------
      IF (presentcase .NE. ppm_poisson_grn_reprojec) THEN
        IF (ppmpoisson%derivatives .EQ. ppm_poisson_drv_curl_fd2) THEN
           CALL ppm_poisson_fd(topoid,meshid,ppmpoisson%drv_vr,fieldout,&
                             & ppm_poisson_drv_curl_fd2,info)
        ENDIF
        IF (ppmpoisson%derivatives .EQ. ppm_poisson_drv_curl_fd4) THEN
           CALL ppm_poisson_fd(topoid,meshid,ppmpoisson%drv_vr,fieldout,&
                             & ppm_poisson_drv_curl_fd4,info)
        ENDIF
      ENDIF

      !-------------------------------------------------------------------------
      ! Finally ghost the velocity/stream function field before returning it
      ! Also extrapolate if freespace
      !-------------------------------------------------------------------------
      CALL ppm_map_field_ghost_get(topoid,meshid,gstw,info)
      CALL ppm_map_field_push(topoid,meshid,fieldout,__NCOM,info)
      CALL ppm_map_field_send(info)
      CALL ppm_map_field_pop(topoid,meshid,fieldout,__NCOM,gstw,info)
      IF (presentcase .EQ. ppm_poisson_grn_pois_fre) THEN
         CALL ppm_poisson_extrapolateghost(topoid,meshid,fieldout,&
                                     & 2,4,gstw,info)
      ENDIF

      !-------------------------------------------------------------------------
      ! Perhaps allocate (and deallocate) arrays !@
      !-------------------------------------------------------------------------

      !-------------------------------------------------------------------------
      ! Return
      !-------------------------------------------------------------------------
 9999 CONTINUE
      CALL substop('ppm_poisson_solve',t0,info)
      RETURN

      END SUBROUTINE __ROUTINE

