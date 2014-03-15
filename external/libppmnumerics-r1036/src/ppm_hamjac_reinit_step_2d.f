      !-------------------------------------------------------------------------
      !     Subroutine   :            ppm_hamjac_reinit_step_2d
      !-------------------------------------------------------------------------
      !     
      !     Purpose      : Solve Hamilton-Jacobi for Gowas reinit
      !      
      !     Input        : 
      !                    
      !     Input/Output : 
      !                    
      !     Output       : 
      !      
      !     Remarks      : 
      !                    
      !     
      !     References   :
      !     
      !     Revisions    :
      !-------------------------------------------------------------------------
      !     $Log: ppm_hamjac_reinit_step_2d.f,v $
      !     Revision 1.1.1.1  2006/07/25 15:18:19  menahel
      !     initial import
      !
      !     Revision 1.1  2005/07/25 00:34:04  ivos
      !     Initial check-in.
      !
      !-------------------------------------------------------------------------
      !     Parallel Particle Mesh Library (PPM)
      !     Institute of Computational Science
      !     ETH Zentrum, Hirschengraben 84
      !     CH-8092 Zurich, Switzerland
      !-------------------------------------------------------------------------


#if   __MODE == __SCA
#if   __KIND == __SINGLE_PRECISION
      SUBROUTINE ppm_hamjac_reinit_step_2ds (phi, tphi, trgt, res, &
           &                          topo_id, mesh_id, ghostsize, info)
#elif __KIND == __DOUBLE_PRECISION
      SUBROUTINE ppm_hamjac_reinit_step_2dd (phi, tphi, trgt, res, topo_id, mesh_id, ghostsize, info)
#endif
#elif __MODE == __VEC
#error VECTOR NOT IMPLEMENTED       
#endif

        USE ppm_module_data
        
        USE ppm_module_error
        USE ppm_module_substart
        USE ppm_module_substop
        USE ppm_module_typedef
        
        IMPLICIT NONE
        
#if    __KIND == __SINGLE_PRECISION
        INTEGER, PARAMETER :: MK = ppm_kind_single
#elif  __KIND == __DOUBLE_PRECISION       
        INTEGER, PARAMETER :: MK = ppm_kind_double
#endif

        !-----------------------------------------------------
        !  Arguments
        !-----------------------------------------------------
        REAL(MK), DIMENSION(:,:,:  ), POINTER :: phi, tphi
        INTEGER, INTENT(in)                   :: topo_id, mesh_id
        INTEGER, DIMENSION(2), INTENT(in)     :: ghostsize
        INTEGER, INTENT(inout)                :: info
        REAL(mk),INTENT(out)                  :: res
        REAL(mk), INTENT(in)                  :: trgt

        !-----------------------------------------------------
        !  Aliases
        !-----------------------------------------------------
        INTEGER, DIMENSION(:), POINTER        :: isublist
        INTEGER                               :: nsublist
        INTEGER, DIMENSION(:,:), POINTER      :: ndata
        INTEGER                               :: topoid, meshid
        REAL(MK), DIMENSION(:), POINTER       :: min_phys, max_phys
        TYPE(ppm_t_topo),      POINTER        :: topo
        TYPE(ppm_t_equi_mesh), POINTER        :: mesh
        
        !-----------------------------------------------------
        !  standard stuff
        !-----------------------------------------------------
        INTEGER                               :: isub,isubl,i,j,k
        REAL(MK)                              :: len_phys(2)
        !-----------------------------------------------------
        !  WENO stuff
        !-----------------------------------------------------
        REAL(mk) :: oneg(2), opos(2), wenoeps, wenotau, pbs
        REAL(mk) :: laps(-1:1,2), rpos(2), rneg(2), dx(2), dxi(2)
        REAL(mk) :: phip(2), phin(2), phimid(2), rms, dphi_dt
        INTEGER  :: ilap
        INTEGER, PARAMETER, DIMENSION(3,3) :: offs &
             & = RESHAPE((/2,1,0,1,0,-1,0,-1,-2/),(/3,3/))
        REAL(mk) :: t0

        
        CALL substart('ppm_hamjac_reinit_step_2d',t0,info)
        
        !-----------------------------------------------------
        !  Get the mesh data
        !-----------------------------------------------------
        topo => ppm_topo(topo_id)%t
        mesh => topo%mesh(mesh_id)
        meshid = mesh%ID
        nsublist = topo%nsublist
        ndata    => mesh%nnodes
        isublist => topo%isublist
#if    __KIND == __SINGLE_PRECISION
        min_phys => topo%min_physs
        max_phys => topo%max_physs
#elif  __KIND == __DOUBLE_PRECISION       
        min_phys => topo%min_physd
        max_phys => topo%max_physd
#endif

        len_phys(1) = max_phys(1) - min_phys(1)
        len_phys(2) = max_phys(2) - min_phys(2)
        dx(1)       = len_phys(1)/REAL(mesh%Nm(1)-1,mk)
        dx(2)       = len_phys(2)/REAL(mesh%Nm(2)-1,mk)
        dxi(1)      = 1.0_mk/dx(1)
        dxi(2)      = 1.0_mk/dx(2)
        wenoeps = 1.0e-6_mk
        wenotau = 0.25_mk*MINVAL(dx)

        rms = -HUGE(rms)

        DO isub=1,nsublist
           isubl = isublist(isub)
           DO j=1,ndata(2,isubl)
              DO i=1,ndata(1,isubl)
                 
                 !DO ilap=1,3
                 !   laps(2-ilap,1) = phi(i+offs(1,ilap),j,isub)   &
                 !        & -2.0_mk * phi(i+offs(2,ilap),j,isub) &
                 !        &       + phi(i+offs(3,ilap),j,isub)
                 !   laps(2-ilap,2) = phi(i,j+offs(1,ilap),isub)   &
                 !        & -2.0_mk * phi(i,j+offs(2,ilap),isub) &
                 !        &       + phi(i,j+offs(3,ilap),isub)
                 !END DO
                 laps(1,1) = phi(i+2,j,isub)   &
                         & -2.0_mk * phi(i+1,j,isub) &
                         &       + phi(i,j,isub)
                 laps(1,2) = phi(i,j+2,isub)   &
                         & -2.0_mk * phi(i,j+1,isub) &
                         &       + phi(i,j,isub)
                 laps(0,1) = phi(i+1,j,isub)   &
                         & -2.0_mk * phi(i,j,isub) &
                         &       + phi(i-1,j,isub)
                 laps(0,2) = phi(i,j+1,isub)   &
                      & -2.0_mk * phi(i,j,isub) &
                      &       + phi(i,j-1,isub)
                 laps(-1,1) = phi(i,j,isub)   &
                      & -2.0_mk * phi(i-1,j,isub) &
                      &       + phi(i-2,j,isub)
                 laps(-1,2) = phi(i,j,isub)   &
                      & -2.0_mk * phi(i,j-1,isub) &
                      &       + phi(i,j-2,isub)
                 
                 

                 rpos(1) = (wenoeps + laps( 1,1)**2)/(wenoeps + laps(0,1)**2)
                 rneg(1) = (wenoeps + laps(-1,1)**2)/(wenoeps + laps(0,1)**2)
                 rpos(2) = (wenoeps + laps( 1,2)**2)/(wenoeps + laps(0,2)**2)
                 rneg(2) = (wenoeps + laps(-1,2)**2)/(wenoeps + laps(0,2)**2)

                 opos(1) = 1.0_mk/(1.0_mk+2.0_mk*rpos(1)**2)
                 opos(2) = 1.0_mk/(1.0_mk+2.0_mk*rpos(2)**2)
                 oneg(1) = 1.0_mk/(1.0_mk+2.0_mk*rneg(1)**2)
                 oneg(2) = 1.0_mk/(1.0_mk+2.0_mk*rneg(2)**2)

                 phimid(1) = phi(i+1,j,isub)-phi(i-1,j,isub)
                 phimid(2) = phi(i,j+1,isub)-phi(i,j-1,isub)

                 phip(1) = 0.5_mk*(phimid(1) - &
                      & opos(1)*( &
                      &         phi(i+2,j,isub) - &
                      & 3.0_mk*(phi(i+1,j,isub) - phi(i  ,j,isub)) - &
                      &         phi(i-1,j,isub)))*dxi(1)
                 phip(2) = 0.5_mk*(phimid(2) - &
                      & opos(2)*( &
                      &         phi(i,j+2,isub) - &
                      & 3.0_mk*(phi(i,j+1,isub) - phi(i  ,j,isub)) - &
                      &         phi(i,j-1,isub)))*dxi(2)
                 phin(1) = 0.5_mk*(phimid(1) - &
                      & oneg(1)*( &
                      &         phi(i+1,j,isub) - &
                      & 3.0_mk*(phi(i  ,j,isub) - phi(i-1,j,isub)) - &
                      &         phi(i-2,j,isub)))*dxi(1)
                 phin(2) = 0.5_mk*(phimid(2) - &
                      & oneg(2)*( &
                      &         phi(i,j+1,isub) - &
                      & 3.0_mk*(phi(i,j  ,isub) - phi(i,j-1,isub)) - &
                      &         phi(i,j-2,isub)))*dxi(2)

                 !--- collect
                 IF(phi(i,j,isub).GT.0.0_mk) THEN
                    pbs = SQRT( &
                         & MAX(-MIN(phip(1),0.0_mk),MAX(phin(1),0.0_mk))**2+&
                         & MAX(-MIN(phip(2),0.0_mk),MAX(phin(2),0.0_mk))**2)&
                         & - trgt
                 ELSEIF(phi(i,j,isub).LT.0.0_mk) THEN
                    pbs = SQRT( &
                         & MAX(MAX(phip(1),0.0_mk),-MIN(phin(1),0.0_mk))**2+&
                         & MAX(MAX(phip(2),0.0_mk),-MIN(phin(2),0.0_mk))**2)&
                         & - trgt
                 ELSE
                    pbs = 0.0_mk
                 END IF
                 dphi_dt =  pbs * phi(i,j,isub) / &
                      & SQRT(phi(i,j,isub)**2+0.25_mk*SUM(phimid**2)) 
                 tphi(i,j,isub) = phi(i,j,isub) - wenotau * dphi_dt

                 rms = MAX(rms,ABS(dphi_dt))

              END DO
              
           END DO
           
        END DO
           


        res = rms

        CALL substop('ppm_hamjac_reinit_step_2d',t0,info)

#if   __KIND == __SINGLE_PRECISION
      END SUBROUTINE ppm_hamjac_reinit_step_2ds 
#elif __KIND == __DOUBLE_PRECISION
      END SUBROUTINE ppm_hamjac_reinit_step_2dd 
#endif

      
                    


                    
           
           



        
        
        
        
