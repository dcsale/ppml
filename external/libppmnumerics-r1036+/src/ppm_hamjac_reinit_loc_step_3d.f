      !-------------------------------------------------------------------------
      !     Subroutine   :           ppm_hamjac_reinit_step_3d
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
      !     $Log: ppm_hamjac_reinit_loc_step_3d.f,v $
      !     Revision 1.1.1.1  2006/07/25 15:18:19  menahel
      !     initial import
      !
      !     Revision 1.2  2005/08/12 14:38:01  ivos
      !     bugfix: index bounds in loop corrected.
      !
      !     Revision 1.1  2005/07/25 00:34:05  ivos
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
      SUBROUTINE ppm_hamjac_reinit_loc_step_3ds(phi,tphi,iloc,np,trgt,res, &
     &                                     topo_id,mesh_id,ghostsize,info)
#elif __KIND == __DOUBLE_PRECISION
      SUBROUTINE ppm_hamjac_reinit_loc_step_3dd(phi,tphi,iloc,np,trgt,res, &
     &                                     topo_id,mesh_id,ghostsize,info)
#endif
#elif __MODE == __VEC
#if   __KIND == __SINGLE_PRECISION
      SUBROUTINE ppm_hamjac_reinit_loc_step_3dsV(phi,idx,tphi,iloc,np,trgt,&
     &                                     res,topo_id,mesh_id,ghostsize,info)
#elif __KIND == __DOUBLE_PRECISION
      SUBROUTINE ppm_hamjac_reinit_loc_step_3ddV(phi,idx,tphi,iloc,np,trgt,&
     &                                     res,topo_id,mesh_id,ghostsize,info)
#endif
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
#if __MODE == __SCA
        REAL(MK), DIMENSION(:,:,:,:), POINTER :: phi
#elif __MODE == __VEC
        REAL(MK), DIMENSION(:,:,:,:,:), POINTER :: phi
#endif        
        REAL(MK), DIMENSION(:,:,:,:), POINTER :: tphi
        INTEGER, INTENT(in)                   :: topo_id, mesh_id
        INTEGER, DIMENSION(3), INTENT(in)     :: ghostsize
        INTEGER, INTENT(inout)                :: info
        REAL(mk), INTENT(out)                 :: res
#if __MODE == __VEC
        INTEGER, INTENT(in)                   :: idx
#endif        
        REAL(mk), INTENT(in)                  :: trgt
        INTEGER, DIMENSION(:,:), INTENT(in)   :: iloc
        INTEGER                               :: np, p
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
        REAL(MK)                              :: len_phys(3)
        !-----------------------------------------------------
        !  WENO stuff
        !-----------------------------------------------------
        REAL(mk) :: oneg(3), opos(3), wenoeps, wenotau, pbs
        REAL(mk) :: laps(-1:1,3), rpos(3), rneg(3), dx(3), dxi(3)
        REAL(mk) :: phip(3), phin(3), phimid(3), rms, dphi_dt
        INTEGER  :: ilap
        INTEGER, PARAMETER, DIMENSION(3,3) :: offs &
             & = RESHAPE((/2,1,0,1,0,-1,0,-1,-2/),(/3,3/))
        REAL(mk) :: t0

        
        CALL substart('ppm_hamjac_reinit_loc_step_3d',t0,info)
        
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
        len_phys(3) = max_phys(3) - min_phys(3)
        dx(1)       = len_phys(1)/REAL(mesh%Nm(1)-1,mk)
        dx(2)       = len_phys(2)/REAL(mesh%Nm(2)-1,mk)
        dx(3)       = len_phys(3)/REAL(mesh%Nm(3)-1,mk)
        dxi(1)      = 1.0_mk/dx(1)
        dxi(2)      = 1.0_mk/dx(2)
        dxi(3)      = 1.0_mk/dx(3)
        wenoeps = 1.0e-6_mk
        wenotau = 0.25_mk*MINVAL(dx)

        rms = -HUGE(rms)

           DO p=1,np
              isub = iloc(4,p)
              i = iloc(1,p)
              j = iloc(2,p)
              k = iloc(3,p)
                    ! hack
#if __MODE == __SCA                    
!                    IF(phi(i+1,j,k,isub).EQ.phi(i-1,j,k,isub).AND. &
!                         & phi(i,j+1,k,isub).EQ.phi(i,j-1,k,isub).AND. &
!                         & phi(i,j,k+1,isub).EQ.phi(i,j,k-1,isub).AND.ABS(phi(i&
!                         &,j,k,isub)).LT.14.0_mk*dx(1)) CYCLE
#endif                    
#if __MODE == __SCA
                    phimid(1) = phi(i+1,j,k,isub)-phi(i-1,j,k,isub)
                    phimid(2) = phi(i,j+1,k,isub)-phi(i,j-1,k,isub)
                    phimid(3) = phi(i,j,k+1,isub)-phi(i,j,k-1,isub)
#else
                    phimid(1) = phi(idx,i+1,j,k,isub)-phi(idx,i-1,j,k,isub)
                    phimid(2) = phi(idx,i,j+1,k,isub)-phi(idx,i,j-1,k,isub)
                    phimid(3) = phi(idx,i,j,k+1,isub)-phi(idx,i,j,k-1,isub)
#endif
                    
#if __MODE == __SCA
                       laps(2-3,1) = phi(i+offs(1,3),j,k,isub)   &
                            & -2.0_mk * phi(i+offs(2,3),j,k,isub) &
                            &       + phi(i+offs(3,3),j,k,isub)
                       laps(2-3,2) = phi(i,j+offs(1,3),k,isub)   &
                            & -2.0_mk * phi(i,j+offs(2,3),k,isub) &
                            &       + phi(i,j+offs(3,3),k,isub)
                       laps(2-3,3) = phi(i,j,k+offs(1,3),isub)   &
                            & -2.0_mk * phi(i,j,k+offs(2,3),isub) &
                            &       + phi(i,j,k+offs(3,3),isub)
                       laps(2-2,1) = phi(i+offs(1,2),j,k,isub)   &
                            & -2.0_mk * phi(i+offs(2,2),j,k,isub) &
                            &       + phi(i+offs(3,2),j,k,isub)
                       laps(2-2,2) = phi(i,j+offs(1,2),k,isub)   &
                            & -2.0_mk * phi(i,j+offs(2,2),k,isub) &
                            &       + phi(i,j+offs(3,2),k,isub)
                       laps(2-2,3) = phi(i,j,k+offs(1,2),isub)   &
                            & -2.0_mk * phi(i,j,k+offs(2,2),isub) &
                            &       + phi(i,j,k+offs(3,2),isub)
                       laps(2-1,1) = phi(i+offs(1,1),j,k,isub)   &
                            & -2.0_mk * phi(i+offs(2,1),j,k,isub) &
                            &       + phi(i+offs(3,1),j,k,isub)
                       laps(2-1,2) = phi(i,j+offs(1,1),k,isub)   &
                            & -2.0_mk * phi(i,j+offs(2,1),k,isub) &
                            &       + phi(i,j+offs(3,1),k,isub)
                       laps(2-1,3) = phi(i,j,k+offs(1,1),isub)   &
                            & -2.0_mk * phi(i,j,k+offs(2,1),isub) &
                            &       + phi(i,j,k+offs(3,1),isub)
#elif __MODE == __VEC
                    DO ilap=1,3
                       laps(2-ilap,1) = phi(idx,i+offs(1,ilap),j,k,isub)   &
                            & -2.0_mk * phi(idx,i+offs(2,ilap),j,k,isub) &
                            &       + phi(idx,i+offs(3,ilap),j,k,isub)
                       laps(2-ilap,2) = phi(idx,i,j+offs(1,ilap),k,isub)   &
                            & -2.0_mk * phi(idx,i,j+offs(2,ilap),k,isub) &
                            &       + phi(idx,i,j+offs(3,ilap),k,isub)
                       laps(2-ilap,3) = phi(idx,i,j,k+offs(1,ilap),isub)   &
                            & -2.0_mk * phi(idx,i,j,k+offs(2,ilap),isub) &
                            &       + phi(idx,i,j,k+offs(3,ilap),isub)
                    END DO
#endif                       

                    rpos(1) = (wenoeps + laps( 1,1)**2)/(wenoeps + laps(0,1)**2)
                    rneg(1) = (wenoeps + laps(-1,1)**2)/(wenoeps + laps(0,1)**2)
                    rpos(2) = (wenoeps + laps( 1,2)**2)/(wenoeps + laps(0,2)**2)
                    rneg(2) = (wenoeps + laps(-1,2)**2)/(wenoeps + laps(0,2)**2)
                    rpos(3) = (wenoeps + laps( 1,3)**2)/(wenoeps + laps(0,3)**2)
                    rneg(3) = (wenoeps + laps(-1,3)**2)/(wenoeps + laps(0,3)**2)

                    opos(1) = 1.0_mk/(1.0_mk+2.0_mk*rpos(1)**2)
                    opos(2) = 1.0_mk/(1.0_mk+2.0_mk*rpos(2)**2)
                    opos(3) = 1.0_mk/(1.0_mk+2.0_mk*rpos(3)**2)
                    oneg(1) = 1.0_mk/(1.0_mk+2.0_mk*rneg(1)**2)
                    oneg(2) = 1.0_mk/(1.0_mk+2.0_mk*rneg(2)**2)
                    oneg(3) = 1.0_mk/(1.0_mk+2.0_mk*rneg(3)**2)

#if __MODE == __SCA
                    phip(1) = 0.5_mk*(phimid(1) - &
                         & opos(1)*( &
                         &         phi(i+2,j,k,isub) - &
                         & 3.0_mk*(phi(i+1,j,k,isub) - phi(i  ,j,k,isub)) - &
                         &         phi(i-1,j,k,isub)))*dxi(1)
                    phip(2) = 0.5_mk*(phimid(2) - &
                         & opos(2)*( &
                         &         phi(i,j+2,k,isub) - &
                         & 3.0_mk*(phi(i,j+1,k,isub) - phi(i  ,j,k,isub)) - &
                         &         phi(i,j-1,k,isub)))*dxi(2)
                    phip(3) = 0.5_mk*(phimid(3) - &
                         & opos(3)*( &
                         &         phi(i,j,k+2,isub) - &
                         & 3.0_mk*(phi(i,j,k+1,isub) - phi(i  ,j,k,isub)) - &
                         &         phi(i,j,k-1,isub)))*dxi(3)
                    phin(1) = 0.5_mk*(phimid(1) - &
                         & oneg(1)*( &
                         &         phi(i+1,j,k,isub) - &
                         & 3.0_mk*(phi(i  ,j,k,isub) - phi(i-1,j,k,isub)) - &
                         &         phi(i-2,j,k,isub)))*dxi(1)
                    phin(2) = 0.5_mk*(phimid(2) - &
                         & oneg(2)*( &
                         &         phi(i,j+1,k,isub) - &
                         & 3.0_mk*(phi(i,j  ,k,isub) - phi(i,j-1,k,isub)) - &
                         &         phi(i,j-2,k,isub)))*dxi(2)
                    phin(3) = 0.5_mk*(phimid(3) - &
                         & oneg(3)*( &
                         &         phi(i,j,k+1,isub) - &
                         & 3.0_mk*(phi(i,j,k  ,isub) - phi(i,j,k-1,isub)) - &
                         &         phi(i,j,k-2,isub)))*dxi(3)
#else
                    phip(1) = 0.5_mk*(phimid(1) - &
                         & opos(1)*( &
                         &         phi(idx,i+2,j,k,isub) - &
                         & 3.0_mk*(phi(idx,i+1,j,k,isub)-phi(idx,i,j,k,isub))-&
                         &         phi(idx,i-1,j,k,isub)))*dxi(1)
                    phip(2) = 0.5_mk*(phimid(2) - &
                         & opos(2)*( &
                         &         phi(idx,i,j+2,k,isub) - &
                         & 3.0_mk*(phi(idx,i,j+1,k,isub)-phi(idx,i,j,k,isub))-&
                         &         phi(idx,i,j-1,k,isub)))*dxi(2)
                    phip(3) = 0.5_mk*(phimid(3) - &
                         & opos(3)*( &
                         &         phi(idx,i,j,k+2,isub) - &
                         & 3.0_mk*(phi(idx,i,j,k+1,isub)-phi(idx,i,j,k,isub))-&
                         &         phi(idx,i,j,k-1,isub)))*dxi(3)
                    phin(1) = 0.5_mk*(phimid(1) - &
                         & oneg(1)*( &
                         &         phi(idx,i+1,j,k,isub) - &
                         & 3.0_mk*(phi(idx,i  ,j,k,isub)-phi(idx,i-1,j,k,isub))-&
                         &         phi(idx,i-2,j,k,isub)))*dxi(1)
                    phin(2) = 0.5_mk*(phimid(2) - &
                         & oneg(2)*( &
                         &         phi(idx,i,j+1,k,isub) - &
                         & 3.0_mk*(phi(idx,i,j  ,k,isub)-phi(idx,i,j-1,k,isub))-&
                         &         phi(idx,i,j-2,k,isub)))*dxi(2)
                    phin(3) = 0.5_mk*(phimid(3) - &
                         & oneg(3)*( &
                         &         phi(idx,i,j,k+1,isub) - &
                         & 3.0_mk*(phi(idx,i,j,k  ,isub)-phi(idx,i,j,k-1,isub))-&
                         &         phi(idx,i,j,k-2,isub)))*dxi(3)
#endif

#if __MODE == __SCA                    
                    !--- collect
                    IF(phi(i,j,k,isub).GT.0.0_mk) THEN
                       pbs = SQRT( &
                            & MAX(-MIN(phip(1),0.0_mk),MAX(phin(1),0.0_mk))**2+&
                            & MAX(-MIN(phip(2),0.0_mk),MAX(phin(2),0.0_mk))**2+&
                            & MAX(-MIN(phip(3),0.0_mk),MAX(phin(3),0.0_mk))**2)&
                            & - trgt
                    ELSEIF(phi(i,j,k,isub).LT.0.0_mk) THEN
                       pbs = SQRT( &
                            & MAX(MAX(phip(1),0.0_mk),-MIN(phin(1),0.0_mk))**2+&
                            & MAX(MAX(phip(2),0.0_mk),-MIN(phin(2),0.0_mk))**2+&
                            & MAX(MAX(phip(3),0.0_mk),-MIN(phin(3),0.0_mk))**2)&
                            & - trgt
                    ELSE
                       pbs = 0.0_mk
                    END IF
                    dphi_dt =  pbs * phi(i,j,k,isub) / &
                         & SQRT(phi(i,j,k,isub)**2+0.25_mk*SUM(phimid**2)) 
                    tphi(i,j,k,isub) = phi(i,j,k,isub) - wenotau * dphi_dt
                    
#else
                    !--- collect
                    IF(phi(idx,i,j,k,isub).GT.0.0_mk) THEN
                       pbs = SQRT( &
                            & MAX(-MIN(phip(1),0.0_mk),MAX(phin(1),0.0_mk))**2+&
                            & MAX(-MIN(phip(2),0.0_mk),MAX(phin(2),0.0_mk))**2+&
                            & MAX(-MIN(phip(3),0.0_mk),MAX(phin(3),0.0_mk))**2)&
                            & - trgt
                    ELSEIF(phi(idx,i,j,k,isub).LT.0.0_mk) THEN
                       pbs = SQRT( &
                            & MAX(MAX(phip(1),0.0_mk),-MIN(phin(1),0.0_mk))**2+&
                            & MAX(MAX(phip(2),0.0_mk),-MIN(phin(2),0.0_mk))**2+&
                            & MAX(MAX(phip(3),0.0_mk),-MIN(phin(3),0.0_mk))**2)&
                            & - trgt
                    ELSE
                       pbs = 0.0_mk
                    END IF
                    dphi_dt =  pbs * phi(idx,i,j,k,isub) / &
                         & SQRT(phi(idx,i,j,k,isub)**2+0.25_mk*SUM(phimid**2)) 
                    tphi(i,j,k,isub) = phi(idx,i,j,k,isub) - wenotau * dphi_dt
                    
#endif                    

                    rms = MAX(rms,ABS(dphi_dt))



        END DO

        res = rms

        CALL substop('ppm_hamjac_reinit_loc_step_3d',t0,info)
#if   __MODE == __SCA
#if   __KIND == __SINGLE_PRECISION
      END SUBROUTINE ppm_hamjac_reinit_loc_step_3ds 
#elif __KIND == __DOUBLE_PRECISION
      END SUBROUTINE ppm_hamjac_reinit_loc_step_3dd 
#endif
#elif __MODE == __VEC
#if   __KIND == __SINGLE_PRECISION
      END SUBROUTINE ppm_hamjac_reinit_loc_step_3dsV 
#elif __KIND == __DOUBLE_PRECISION
      END SUBROUTINE ppm_hamjac_reinit_loc_step_3ddV 
#endif
#endif      

      
                    


                    
           
           



        
        
        
        
