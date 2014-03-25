      !-------------------------------------------------------------------------
      !  Subroutine   : ppm_poisson_finalize.f90
      !-------------------------------------------------------------------------
      ! Copyright (c) 2010 CSE Lab (ETH Zurich), MOSAIC Group (ETH Zurich),
      !                    Center for Fluid Dynamics (DTU)
      !
      !-------------------------------------------------------------------------
      SUBROUTINE __ROUTINE(ppmpoisson,info)
      !!! Routine to deallocate the data initialised for Greens function 
      !!! solution of the Poisson


      IMPLICIT NONE
      !-------------------------------------------------------------------------
      ! Arguments
      !-------------------------------------------------------------------------
      TYPE(ppm_poisson_plan),INTENT(INOUT)                        :: ppmpoisson
      !!! The PPM Poisson plan type (inspired by the FFTW plan)
      INTEGER, INTENT(OUT)                                        :: info

      !-------------------------------------------------------------------------
      ! Local variables
      !-------------------------------------------------------------------------
      REAL(__PREC)                            :: t0
      INTEGER                                 :: stat,info2

      !-------------------------------------------------------------------------
      ! Initialise routine
      !-------------------------------------------------------------------------
      CALL substart('ppm_poisson_finalize',t0,info)


      IF(ASSOCIATED(ppmpoisson%drv_vr)) THEN
        DEALLOCATE(ppmpoisson%drv_vr,stat=info)
        IF (info .NE. 0) THEN
          CALL ppm_write(ppm_rank,'ppm_poisson_finalize','Failed deallocate drv_vr.',info2)
          GOTO 9999
        ENDIF
      ENDIF

      IF(ASSOCIATED(ppmpoisson%fldxyr)) THEN
        DEALLOCATE(ppmpoisson%fldxyr,stat=info)
        IF (info .NE. 0) THEN
          CALL ppm_write(ppm_rank,'ppm_poisson_finalize','Failed deallocate fldxyr.',info2)
          GOTO 9999
        ENDIF
      ENDIF

      IF(ASSOCIATED(ppmpoisson%fldxyc)) THEN
        DEALLOCATE(ppmpoisson%fldxyc,stat=info)
        IF (info .NE. 0) THEN
          CALL ppm_write(ppm_rank,'ppm_poisson_finalize','Failed deallocate fldxyc.',info2)
          GOTO 9999
        ENDIF
      ENDIF

      IF(ASSOCIATED(ppmpoisson%fldzc1)) THEN
        DEALLOCATE(ppmpoisson%fldzc1,stat=info)
        IF (info .NE. 0) THEN
          CALL ppm_write(ppm_rank,'ppm_poisson_finalize','Failed deallocate fldzc1.',info2)
          GOTO 9999
        ENDIF
      ENDIF

      IF(ASSOCIATED(ppmpoisson%fldzc2)) THEN
        DEALLOCATE(ppmpoisson%fldzc2,stat=info)
        IF (info .NE. 0) THEN
          CALL ppm_write(ppm_rank,'ppm_poisson_finalize','Failed deallocate fldzc2.',info2)
          GOTO 9999
        ENDIF
      ENDIF

      IF(ASSOCIATED(ppmpoisson%fldgrnr)) THEN
        DEALLOCATE(ppmpoisson%fldgrnr,stat=info)
        IF (info .NE. 0) THEN
          CALL ppm_write(ppm_rank,'ppm_poisson_finalize','Failed deallocate fldgrnr.',info2)
          GOTO 9999
        ENDIF
      ENDIF

      IF(ASSOCIATED(ppmpoisson%fldgrnc)) THEN
        DEALLOCATE(ppmpoisson%fldgrnc,stat=info)
        IF (info .NE. 0) THEN
          CALL ppm_write(ppm_rank,'ppm_poisson_finalize','Failed deallocate fldgrnc.',info2)
          GOTO 9999
        ENDIF
      ENDIF

      !the following topologies should be removed
      ! ppmpoisson%topoidxy
      ! ppmpoisson%meshidxy
      ! ppmpoisson%meshidxyc
      ! ppmpoisson%topoidz
      ! ppmpoisson%meshidz

      !Eventually the following fft plans must be deallocated
      ! ppmpoisson%planfxy
      ! ppmpoisson%planbxy
      ! ppmpoisson%planfz
      ! ppmpoisson%planbz



      !-------------------------------------------------------------------------
      ! Return
      !-------------------------------------------------------------------------
 9999 CONTINUE
      CALL substop('ppm_poisson_finalize',t0,info)
      RETURN

      END SUBROUTINE __ROUTINE

