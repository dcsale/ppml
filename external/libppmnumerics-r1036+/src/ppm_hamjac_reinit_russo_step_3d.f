      !-------------------------------------------------------------------------
      !     Subroutine   :            ppm_hamjac_reinit_step_3d
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
      !     $Log: ppm_hamjac_reinit_step_3d.f,v $
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
      SUBROUTINE ppm_hamjac_reinit_russo_step_3ds (phi, phi0, phirhs,phigrad,res, s2didx,&
           &                          topo_id, mesh_id, ghostsize, info)
#elif __KIND == __DOUBLE_PRECISION
      SUBROUTINE ppm_hamjac_reinit_russo_step_3dd (phi, phi0, phirhs,phigrad,res, s2didx,& 
	       &						  topo_id, mesh_id, ghostsize, info)
#endif
#elif __MODE == __VEC
#error VECTOR NOT IMPLEMENTED       
#endif

        USE ppm_module_data
        
        USE ppm_module_error
        USE ppm_module_substart
        USE ppm_module_substop
        
        IMPLICIT NONE
        
#if    __KIND == __SINGLE_PRECISION
        INTEGER, PARAMETER :: MK = ppm_kind_single
#elif  __KIND == __DOUBLE_PRECISION       
        INTEGER, PARAMETER :: MK = ppm_kind_double
#endif

        !-----------------------------------------------------
        !  Arguments
        !-----------------------------------------------------
        REAL(MK), DIMENSION(:,:,:,:  ), POINTER :: phi, phi0,phirhs
		REAL(MK), DIMENSION(:,:,:,:,:), POINTER :: phigrad
        INTEGER, INTENT(in)                   :: topo_id, mesh_id
        INTEGER, DIMENSION(3), INTENT(in)     :: ghostsize
        INTEGER, INTENT(inout)                :: info
        REAL(mk),INTENT(out)                  :: res
		INTEGER, INTENT(in)					  :: s2didx

        !-----------------------------------------------------
        !  Aliases
        !-----------------------------------------------------
        INTEGER, DIMENSION(:), POINTER        :: isublist
        INTEGER                               :: nsublist
        INTEGER, DIMENSION(:,:), POINTER      :: ndata
        INTEGER                               :: topoid, meshid
        REAL(mk), DIMENSION(:), POINTER       :: min_phys, max_phys
        TYPE(ppm_t_topo),      POINTER        :: topo
        TYPE(ppm_t_equi_mesh), POINTER        :: mesh
        
        !-----------------------------------------------------
        !  standard stuff
        !-----------------------------------------------------
        INTEGER                               :: isub,isubl,i,j,k,kmin,kmax
        REAL(mk)                              :: len_phys(3)
		REAL(mk)							  :: delta_phis(3,3), delta(3), D
        !-----------------------------------------------------
        !  WENO stuff
        !-----------------------------------------------------
        REAL(mk) :: oneg(3), opos(3), closeeps, wenoeps, pbs
        REAL(mk) :: laps(-1:1,3), rpos(3), rneg(3), dx(3), dxi(3)
        REAL(mk) :: phip(3), phin(3), phimid(3), rms, dphi_dt
        INTEGER  :: ilap
        INTEGER, PARAMETER, DIMENSION(3,3) :: offs &
             & = RESHAPE((/2,1,0,1,0,-1,0,-1,-2/),(/3,3/))
        REAL(mk) :: t0

        
        CALL substart('ppm_hamjac_reinit_russo_step_3d',t0,info)
        
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
		dxi(3)		= 1.0_mk/dx(3)
		
        closeeps = 1.0e-6_mk
		wenoeps = 1.0e-6_mk		

        rms = -HUGE(rms)

		!-----------------------------------------------------
		!  compute right hand side
		!-----------------------------------------------------
		
        DO isub=1,nsublist
           isubl = isublist(isub)
		   IF(s2didx.LT.0) THEN
			  kmin = 1 + ghostsize(3)
			  kmax = ndata(3,isubl) - ghostsize(3)
		   ELSE
			  kmin = s2didx
			  kmax = s2didx
		   END IF		   
		   DO k=kmin,kmax
           DO j=1,ndata(2,isubl)
              DO i=1,ndata(1,isubl)

				IF( phi0(i,j,k,isub)*phi0(i+1,j,k,isub).LT.0.0_mk.OR.&
				  & phi0(i,j,k,isub)*phi0(i-1,j,k,isub).LT.0.0_mk.OR.&
				  & phi0(i,j,k,isub)*phi0(i,j+1,k,isub).LT.0.0_mk.OR.&
				  & phi0(i,j,k,isub)*phi0(i,j-1,k,isub).LT.0.0_mk.OR.&
				  & phi0(i,j,k,isub)*phi0(i,j,k+1,isub).LT.0.0_mk.OR.&
				  & phi0(i,j,k,isub)*phi0(i,j,k-1,isub).LT.0.0_mk) THEN
				  ! near interface: apply russo smereka (2000) technique to 
				  ! prevent interface from moving
					
					delta_phis(1,1) = 0.5_mk*(phi0(i+1,j,k,isub)-phi0(i-1,j,k,isub))**2
					delta_phis(1,2) = 0.5_mk*(phi0(i,j+1,k,isub)-phi0(i,j-1,k,isub))**2
					delta_phis(1,3) = 0.5_mk*(phi0(i,j,k+1,isub)-phi0(i,j,k-1,isub))**2
					
					delta_phis(2,1) = (phi0(i+1,j,k,isub)-phi0(i,j,k,isub))**2
					delta_phis(2,2) = (phi0(i,j+1,k,isub)-phi0(i,j,k,isub))**2
					delta_phis(2,3) = (phi0(i,j,k+1,isub)-phi0(i,j,k,isub))**2
					
					delta_phis(3,1) = (phi0(i,j,k,isub)-phi0(i-1,j,k,isub))**2
					delta_phis(3,2) = (phi0(i,j,k,isub)-phi0(i,j-1,k,isub))**2
					delta_phis(3,3) = (phi0(i,j,k,isub)-phi0(i,j,k-1,isub))**2					
					
					delta(1) = MAX(MAX(delta_phis(1,1),delta_phis(2,1)), &
							     & MAX(delta_phis(3,1),closeeps))
					delta(2) = MAX(MAX(delta_phis(1,2),delta_phis(2,2)), &
							     & MAX(delta_phis(3,2),closeeps))
				    delta(3) = MAX(MAX(delta_phis(1,3),delta_phis(2,3)), &
								 & MAX(delta_phis(3,3),closeeps))
					
					!here assume dx(1)==dx(2) need to do the dx(1)~=dx(2) case
					D = dx(1)*phi0(i,j,k,isub)/(SQRT(delta(1)+delta(2)+delta(3)))
					
					IF(phi0(i,j,k,isub).GT.0.0_mk) THEN
						phirhs(i,j,k,isub) = -dxi(1)*(ABS(phi(i,j,k,isub))-D)

					ELSE IF(phi0(i,j,k,isub).LT.0.0_mk) THEN
						phirhs(i,j,k,isub) = -dxi(1)*(-1.0_mk*ABS(phi(i,j,k,isub))-D)
					ELSE
						phirhs(i,j,k,isub) = 0.0_mk
					END IF
						! gradient not yet correctly implemented 
					    phigrad(1,i,j,k,isub) =  dxi(1)*D
					    phigrad(2,i,j,k,isub) =  dxi(2)*D
					    phigrad(3,i,j,k,isub) =  dxi(3)*D

				ELSE
					! away from interface: normal flux calculation

                    phimid(1) = phi(i+1,j,k,isub)-phi(i-1,j,k,isub)
                    phimid(2) = phi(i,j+1,k,isub)-phi(i,j-1,k,isub)
                    phimid(3) = phi(i,j,k+1,isub)-phi(i,j,k-1,isub)
					
					laps(-1,1) = phi(i+2,j,k,isub)   &
							 & -2.0_mk * phi(i+1,j,k,isub) &
							 &       + phi(i,j,k,isub)
					 laps(-1,2) = phi(i,j+2,k,isub)   &
							 & -2.0_mk * phi(i,j+1,k,isub) &
							 &       + phi(i,j,k,isub)
					 laps(-1,3) = phi(i,j,k+2,isub)   &
							 & -2.0_mk * phi(i,j,k+1,isub) &
							 &       + phi(i,j,k,isub)
							 
					 laps(0,1) = phi(i+1,j,k,isub)   &
							 & -2.0_mk * phi(i,j,k,isub) &
							 &       + phi(i-1,j,k,isub)
					 laps(0,2) = phi(i,j+1,k,isub)   &
						  & -2.0_mk * phi(i,j,k,isub) &
						  &       + phi(i,j-1,k,isub)
					 laps(0,3) = phi(i,j,k+1,isub)   &
						  & -2.0_mk * phi(i,j,k,isub) &
						  &       + phi(i,j,k-1,isub)	  
						
					 laps(1,1) = phi(i,j,k,isub)   &
						  & -2.0_mk * phi(i-1,j,k,isub) &
						  &       + phi(i-2,j,k,isub)
					 laps(1,2) = phi(i,j,k,isub)   &
						  & -2.0_mk * phi(i,j-1,k,isub) &
						  &       + phi(i,j-2,k,isub)
					 laps(1,3) = phi(i,j,k,isub)   &
						  & -2.0_mk * phi(i,j,k-1,isub) &
						  &       + phi(i,j,k-2,isub)
						  

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
                    
                    !--- collect
					IF(phi0(i,j,k,isub).GT.0.0_mk) THEN
                       pbs = SQRT( &
                            & MAX(-MIN(phip(1),0.0_mk),MAX(phin(1),0.0_mk))**2+&
                            & MAX(-MIN(phip(2),0.0_mk),MAX(phin(2),0.0_mk))**2+&
                            & MAX(-MIN(phip(3),0.0_mk),MAX(phin(3),0.0_mk))**2)&
                            & - 1.0_mk
					   
					   phirhs(i,j,k,isub) =  -1.0_mk*pbs

					   					
					ELSE IF(phi0(i,j,k,isub).LT.0.0_mk) THEN
                       pbs = SQRT( &
                            & MAX(MAX(phip(1),0.0_mk),-MIN(phin(1),0.0_mk))**2+&
                            & MAX(MAX(phip(2),0.0_mk),-MIN(phin(2),0.0_mk))**2+&
                            & MAX(MAX(phip(3),0.0_mk),-MIN(phin(3),0.0_mk))**2)&
                            & - 1.0_mk

					   phirhs(i,j,k,isub) =  +1.0_mk*pbs
					   					   
					ELSE
					   pbs = 0.0_mk
					   phirhs(i,j,k,isub) = 0.0_mk				
					END IF	
					   ! gradient not yet correctly implemented 	
					   phigrad(1,i,j,k,isub) = MAX(-MIN(phip(1),0.0_mk),MAX(phin(1),0.0_mk))
					   phigrad(2,i,j,k,isub) = MAX(-MIN(phip(2),0.0_mk),MAX(phin(2),0.0_mk))
					   phigrad(3,i,j,k,isub) = MAX(-MIN(phip(3),0.0_mk),MAX(phin(3),0.0_mk))			


				END IF	 
				! simple euler for now
                !tphi(i,j,k,isub) = phi(i,j,k,isub) + reinit_tau * flux

                rms = MAX(rms,ABS(phirhs(i,j,k,isub)))

              END DO
              
           END DO
		   END DO
           
        END DO
		
		
		IF(s2didx.GT.0) THEN
        DO isub=1,nsublist		
           isubl = isublist(isub)
		   DO k=1,ndata(3,isubl)		   
			   DO j=1,ndata(2,isubl)
				  DO i=1,ndata(1,isubl)		
					phirhs(i,j,k,isub) = phirhs(i,j,s2didx,isub)
					phigrad(1,i,j,k,isub) = phigrad(1,i,j,s2didx,isub)		
					phigrad(2,i,j,k,isub) = phigrad(2,i,j,s2didx,isub)		
					phigrad(3,i,j,k,isub) = phigrad(3,i,j,s2didx,isub)												
				  END DO
				END DO
			END DO
		END DO
		END IF
        res = rms

        CALL substop('ppm_hamjac_reinit_russo_step_3d',t0,info)

#if   __KIND == __SINGLE_PRECISION
      END SUBROUTINE ppm_hamjac_reinit_russo_step_3ds 
#elif __KIND == __DOUBLE_PRECISION
      END SUBROUTINE ppm_hamjac_reinit_russo_step_3dd 
#endif

      
                    


                    
           
           



        
        
        
        
