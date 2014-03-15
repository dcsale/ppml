      !-------------------------------------------------------------------------
      !     Subroutine   :           ppm_hamjac_reinit_step_ref_3d
      !-------------------------------------------------------------------------
      !     
      !     Purpose      : Solve Hamilton-Jacobi for Gowas reinit in ref
      !                    spc
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
      !     $Log: ppm_hamjac_reinit_step_ref_3d.f,v $
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
      SUBROUTINE ppm_hamjac_reinit_step_ref_3ds (phi, chi, tphi, trgt, res, &
     &                                topo_id, mesh_id, ghostsize, info)
#elif __KIND == __DOUBLE_PRECISION
      SUBROUTINE ppm_hamjac_reinit_step_ref_3dd (phi, chi, tphi, trgt, res, &
     &                                topo_id, mesh_id, ghostsize, info)
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
        REAL(MK), DIMENSION(:,:,:,:  ), POINTER :: phi, tphi
        REAL(MK), DIMENSION(:,:,:,:,:), POINTER :: chi
        INTEGER, INTENT(in)                   :: topo_id, mesh_id
        INTEGER, DIMENSION(3), INTENT(in)     :: ghostsize
        INTEGER, INTENT(inout)                :: info
        real(mk),INTENT(out)                  :: res
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
        REAL(MK)                              :: len_phys(3)
        !-----------------------------------------------------
        !  WENO stuff
        !-----------------------------------------------------
        REAL(mk) :: oneg(3), opos(3), wenoeps, wenotau, pbs
        REAL(mk) :: laps(-1:1,3), rpos(3), rneg(3), dx(3), dxi(3)
        REAL(mk) :: phip(3), phin(3), phimid(3), rms, dphi_dt
        REAL(mk) :: dxihalf, dyihalf, dzihalf, phinx(3), phipx(3)
        REAL(mk) :: dxitwelve, dyitwelve, dzitwelve
        INTEGER  :: ilap, order, jsub
        REAL(mk) :: ji(3,3), jac(3,3)
        INTEGER, PARAMETER, DIMENSION(3,3) :: offs &
             & = RESHAPE((/2,1,0,1,0,-1,0,-1,-2/),(/3,3/))
        REAL(mk) :: t0

        
        CALL substart('ppm_hamjac_step_3d',t0,info)
        
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
        dxihalf     = 0.5_mk*dxi(1)
        dyihalf     = 0.5_mk*dxi(2)
        dzihalf     = 0.5_mk*dxi(3)
        dxitwelve   = dxi(1)/12.0_mk
        dyitwelve   = dxi(2)/12.0_mk
        dzitwelve   = dxi(3)/12.0_mk
        
        wenoeps = 1.0e-6_mk
        wenotau = 0.5_mk*MINVAL(dx)

        rms = -HUGE(rms)

        DO isub=1,nsublist
           isubl = isublist(isub)
           DO k=1,ndata(3,isubl)
              DO j=1,ndata(2,isubl)
                 DO i=1,ndata(1,isubl)

                    DO ilap=1,3
                       laps(2-ilap,1) = phi(i+offs(1,ilap),j,k,isub)   &
                            & -2.0_mk * phi(i+offs(2,ilap),j,k,isub) &
                            &       + phi(i+offs(3,ilap),j,k,isub)
                       laps(2-ilap,2) = phi(i,j+offs(1,ilap),k,isub)   &
                            & -2.0_mk * phi(i,j+offs(2,ilap),k,isub) &
                            &       + phi(i,j+offs(3,ilap),k,isub)
                       laps(2-ilap,3) = phi(i,j,k+offs(1,ilap),isub)   &
                            & -2.0_mk * phi(i,j,k+offs(2,ilap),isub) &
                            &       + phi(i,j,k+offs(3,ilap),isub)
                    END DO

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

                    phimid(1) = phi(i+1,j,k,isub)-phi(i-1,j,k,isub)
                    phimid(2) = phi(i,j+1,k,isub)-phi(i,j-1,k,isub)
                    phimid(3) = phi(i,j,k+1,isub)-phi(i,j,k-1,isub)

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

                    jsub = isub
#include "ppm_gmm_jacobian.inc"

                    phinx(1) = jac(1,1)*phin(1)+jac(2,1)*phin(2)+jac(3,1)*phin(3)
                    phinx(2) = jac(1,2)*phin(1)+jac(2,2)*phin(2)+jac(3,2)*phin(3)
                    phinx(3) = jac(1,3)*phin(1)+jac(2,3)*phin(2)+jac(3,3)*phin(3)
                    phipx(1) = jac(1,1)*phip(1)+jac(2,1)*phip(2)+jac(3,1)*phip(3)
                    phipx(2) = jac(1,2)*phip(1)+jac(2,2)*phip(2)+jac(3,2)*phip(3)
                    phipx(3) = jac(1,3)*phip(1)+jac(2,3)*phip(2)+jac(3,3)*phip(3)
                    
                    
                    !--- collect
                    IF(phi(i,j,k,isub).GT.0.0_mk) THEN
                       pbs = SQRT( &
                            & MAX(-MIN(phipx(1),0.0_mk),MAX(phinx(1),0.0_mk))**2+&
                            & MAX(-MIN(phipx(2),0.0_mk),MAX(phinx(2),0.0_mk))**2+&
                            & MAX(-MIN(phipx(3),0.0_mk),MAX(phinx(3),0.0_mk))**2)&
                            & - trgt
                    ELSEIF(phi(i,j,k,isub).LT.0.0_mk) THEN
                       pbs = SQRT( &
                            & MAX(MAX(phipx(1),0.0_mk),-MIN(phinx(1),0.0_mk))**2+&
                            & MAX(MAX(phipx(2),0.0_mk),-MIN(phinx(2),0.0_mk))**2+&
                            & MAX(MAX(phipx(3),0.0_mk),-MIN(phinx(3),0.0_mk))**2)&
                            & - trgt
                    ELSE
                       pbs = 0.0_mk
                    END IF
                    dphi_dt =  pbs * phi(i,j,k,isub) / &
                         & SQRT(phi(i,j,k,isub)**2+0.25_mk*SUM(phimid**2)) 
                    tphi(i,j,k,isub) = phi(i,j,k,isub) - wenotau * dphi_dt

                    rms = MAX(rms,ABS(dphi_dt))

                 END DO

              END DO

           END DO

        END DO

        res = rms

        CALL substop('ppm_hamjac_step_3d',t0,info)

#if   __KIND == __SINGLE_PRECISION
      END SUBROUTINE ppm_hamjac_reinit_step_ref_3ds 
#elif __KIND == __DOUBLE_PRECISION
      END SUBROUTINE ppm_hamjac_reinit_step_ref_3dd 
#endif

      
                    


                    
           
           



        
        
        
        
