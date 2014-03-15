      !-------------------------------------------------------------------------
      !  Subroutine   : ppm_poisson_fd.f90
      !-------------------------------------------------------------------------
      ! Copyright (c) 2010 CSE Lab (ETH Zurich), MOSAIC Group (ETH Zurich),
      !                    Center for Fluid Dynamics (DTU)
      !
      !-------------------------------------------------------------------------
      SUBROUTINE __ROUTINE(topoid,meshid,fieldin,fieldout,dtype,info)
      !!! This routine computes the finite difference gradients, curl etc of 
      !!! fieldin and outputs to fieldout. Both in and out fields are on
      !!! the mesh with meshid belonging to the topology with id topoid.
      !!! The finite difference to be carried out is determined by dtype which
      !!! must be one of the following:
      !!! * ppm_poisson_drv_curl_fd2
      !!! * ppm_poisson_drv_grad_fd2 (not implemented yet)
      !!! * ppm_poisson_drv_lapl_fd2 (not implemented yet)
      !!! * ppm_poisson_drv_div_fd2  (not implemented yet)
      !!! * ppm_poisson_drv_curl_fd4
      !!! * ppm_poisson_drv_grad_fd4 (not implemented yet)
      !!! * ppm_poisson_drv_lapl_fd4 (not implemented yet)
      !!! * ppm_poisson_drv_div_fd4  (not implemented yet)
      !!!
      !!! [NOTE] fieldin and fieldout must NOT be the same array. A check
      !!! should be added.
      !@ TODO: Somewhere check if fieldin is equal to fieldout and give a warning

      USE ppm_module_topo_get

      IMPLICIT NONE
      !-------------------------------------------------------------------------
      ! Arguments
      !-------------------------------------------------------------------------
      INTEGER, INTENT(IN)                                         :: topoid
      !!! ID of the topology
      INTEGER, INTENT(IN)                                         :: meshid
      !!! Mesh ID
      REAL(__PREC),DIMENSION(:,:,:,:,:),POINTER                   :: fieldin
      !!! Input field data
      REAL(__PREC),DIMENSION(:,:,:,:,:),POINTER                   :: fieldout
      !!! Output field data
      INTEGER, INTENT(IN)                                         :: dtype
      !!! Derivation type. Can be one of the types:
      !!! * ppm_poisson_drv_curl_fd2
      !!! * ppm_poisson_drv_grad_fd2 (not implemented yet)
      !!! * ppm_poisson_drv_lapl_fd2 (not implemented yet)
      !!! * ppm_poisson_drv_div_fd2  (not implemented yet)
      !!! * ppm_poisson_drv_curl_fd4
      !!! * ppm_poisson_drv_grad_fd4 (not implemented yet)
      !!! * ppm_poisson_drv_lapl_fd4 (not implemented yet)
      !!! * ppm_poisson_drv_div_fd4  (not implemented yet)
      INTEGER, INTENT(OUT)                                        :: info
      !!! Return status, 0 upon succes

      !-------------------------------------------------------------------------
      ! Local variables
      !-------------------------------------------------------------------------
      INTEGER,PARAMETER                 :: MK = __PREC
      REAL(__PREC)                      :: t0
      TYPE(ppm_t_topo),POINTER          :: topology
      TYPE(ppm_t_equi_mesh)             :: mesh
      REAL(__PREC)                      :: dx,dy,dz
      REAL(__PREC)                      :: facx,facy,facz
      INTEGER                           :: isub,isubl
      INTEGER                           :: i,j,k

      !-------------------------------------------------------------------------
      ! Initialise routine
      !-------------------------------------------------------------------------
      CALL substart('ppm_poisson_fd',t0,info)

      !-------------------------------------------------------------------------
      ! Get topology and mesh values
      !-------------------------------------------------------------------------
      CALL ppm_topo_get(topoid,topology,info)
      IF (info .NE. 0) THEN
        CALL ppm_write(ppm_rank,'ppm_poisson_init_predef','Failed to get topology.',isub)
        GOTO 9999
      ENDIF
      mesh  = topology%mesh(meshid)

      dx = (topology%max_physd(1)-topology%min_physd(1))/REAL(mesh%nm(1)-1) !vertex
      dy = (topology%max_physd(2)-topology%min_physd(2))/REAL(mesh%nm(2)-1)
      dz = (topology%max_physd(3)-topology%min_physd(3))/REAL(mesh%nm(3)-1)

      !-----------------------------------------------------------------------
      ! Do the finite difference calculation
      !-----------------------------------------------------------------------
      !-----------------------------------------------------------------------
      ! Curl, 2nd order FD
      !-----------------------------------------------------------------------
      IF (dtype .EQ. ppm_poisson_drv_curl_fd2) THEN
        facx = 1.0_MK/(2.0_MK*dx)
        facy = 1.0_MK/(2.0_MK*dy)
        facz = 1.0_MK/(2.0_MK*dz)

        DO isub=1,topology%nsublist
          isubl=topology%isublist(isub)
          DO k=1,mesh%nnodes(3,isubl)
            DO j=1,mesh%nnodes(2,isubl)
              DO i=1,mesh%nnodes(1,isubl)
                fieldout(1,i,j,k,isub) = &
                &  facy*(fieldin(3,i  ,j+1,k  ,isub)- &
                       & fieldin(3,i  ,j-1,k  ,isub)) &
                & -facz*(fieldin(2,i  ,j  ,k+1,isub)- &
                       & fieldin(2,i  ,j  ,k-1,isub))
                fieldout(2,i,j,k,isub) = &
                &  facz*(fieldin(1,i  ,j  ,k+1,isub)- &
                       & fieldin(1,i  ,j  ,k-1,isub)) &
                & -facx*(fieldin(3,i+1,j  ,k  ,isub)- &
                       & fieldin(3,i-1,j  ,k  ,isub))
                fieldout(3,i,j,k,isub) = &
                &  facx*(fieldin(2,i+1,j  ,k  ,isub)- &
                       & fieldin(2,i-1,j  ,k  ,isub)) &
                & -facy*(fieldin(1,i  ,j+1,k  ,isub)- &
                       & fieldin(1,i  ,j-1,k  ,isub))
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      !-----------------------------------------------------------------------
      ! Curl, 4th order FD - untested
      !-----------------------------------------------------------------------
      ELSE IF (dtype .EQ. ppm_poisson_drv_curl_fd4) THEN
        facx = 1.0_MK/(12.0_MK*dx)
        facy = 1.0_MK/(12.0_MK*dy)
        facz = 1.0_MK/(12.0_MK*dz)

        DO isub=1,topology%nsublist
          isubl=topology%isublist(isub)
          DO k=1,mesh%nnodes(3,isubl)
            DO j=1,mesh%nnodes(2,isubl)
              DO i=1,mesh%nnodes(1,isubl)
                fieldout(1,i,j,k,isub) = &
                &  facy*(         -fieldin(3,i  ,j+2,k  ,isub)  &
                       &   +8.0_MK*fieldin(3,i  ,j+1,k  ,isub)  &
                       &   -8.0_MK*fieldin(3,i  ,j-1,k  ,isub)  &
                       &          +fieldin(3,i  ,j-2,k  ,isub)) &
                & -facz*(         -fieldin(2,i  ,j  ,k+2,isub)  &
                       &   +8.0_MK*fieldin(2,i  ,j  ,k+1,isub)  &
                       &   -8.0_MK*fieldin(2,i  ,j  ,k-1,isub)  &
                       &          +fieldin(2,i  ,j  ,k-2,isub))
                fieldout(2,i,j,k,isub) = &
                &  facz*(         -fieldin(1,i  ,j  ,k+2,isub)  &
                       &   +8.0_MK*fieldin(1,i  ,j  ,k+1,isub)  &
                       &   -8.0_MK*fieldin(1,i  ,j  ,k-1,isub)  &
                       &          +fieldin(1,i  ,j  ,k-2,isub)) &
                & -facx*(         -fieldin(3,i+2,j  ,k  ,isub)  &
                       &   +8.0_MK*fieldin(3,i+1,j  ,k  ,isub)  &
                       &   -8.0_MK*fieldin(3,i-1,j  ,k  ,isub)  &
                       &          +fieldin(3,i-2,j  ,k  ,isub))
                fieldout(3,i,j,k,isub) = &
                &  facx*(         -fieldin(2,i+2,j  ,k  ,isub)  &
                       &   +8.0_MK*fieldin(2,i+1,j  ,k  ,isub)  &
                       &   -8.0_MK*fieldin(2,i-1,j  ,k  ,isub)  &
                       &          +fieldin(2,i-2,j  ,k  ,isub)) &
                & -facy*(         -fieldin(1,i  ,j+2,k  ,isub)  &
                       &   +8.0_MK*fieldin(1,i  ,j+1,k  ,isub)  &
                       &   -8.0_MK*fieldin(1,i  ,j-1,k  ,isub)  &
                       &          +fieldin(1,i  ,j-2,k  ,isub))
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ENDIF

      !-------------------------------------------------------------------------
      ! Return
      !-------------------------------------------------------------------------
 9999 CONTINUE
      CALL substop('ppm_poisson_fd',t0,info)
      RETURN

      END SUBROUTINE __ROUTINE

