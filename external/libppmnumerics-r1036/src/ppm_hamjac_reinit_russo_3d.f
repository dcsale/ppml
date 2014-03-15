      !-------------------------------------------------------------------------
      !     Subroutine   :                 ppm_hamjac_reinit_3d
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
      !     $Log: ppm_hamjac_reinit_3d.f,v $
      !     Revision 1.1.1.1  2006/07/25 15:18:19  menahel
      !     initial import
      !
      !     Revision 1.1  2005/07/25 00:34:02  ivos
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
      SUBROUTINE ppm_hamjac_reinit_russo_3ds (phi, phigrad, tol, maxstep, &
           &                     topo_id, mesh_id, ghostsize, info, indx)
#elif __KIND == __DOUBLE_PRECISION
      SUBROUTINE ppm_hamjac_reinit_russo_3dd (phi, phigrad, tol, maxstep, &
           &                     topo_id, mesh_id, ghostsize, info, indx)
#endif
#elif __MODE == __VEC
#error VECTOR NOT IMPLEMENTED       
#endif

        USE ppm_module_data
        
        USE ppm_module_error
        USE ppm_module_write
        USE ppm_module_substart
        USE ppm_module_alloc
        USE ppm_module_substop
        USE ppm_module_map
        USE ppm_module_map_field
        USE ppm_module_map_field_ghost
        IMPLICIT NONE

#if    __KIND == __SINGLE_PRECISION
        INTEGER, PARAMETER :: MK = ppm_kind_single
#elif  __KIND == __DOUBLE_PRECISION       
        INTEGER, PARAMETER :: MK = ppm_kind_double
#endif

#ifdef __MPI
      INCLUDE  'mpif.h'
#endif
        !-----------------------------------------------------
        !  Arguments
        !-----------------------------------------------------
        REAL(MK), DIMENSION(:,:,:,:  ), POINTER :: phi
        REAL(MK), DIMENSION(:,:,:,:,:), POINTER :: phigrad
        INTEGER, INTENT(in)                   :: topo_id, mesh_id
        INTEGER, DIMENSION(3), INTENT(in)     :: ghostsize
        INTEGER, INTENT(inout)                :: info
        INTEGER, INTENT(in)                   :: maxstep
        REAL(mk), INTENT(in)                  :: tol
        INTEGER,INTENT(IN),OPTIONAL              :: indx
        !-----------------------------------------------------
        !  Aliases
        !-----------------------------------------------------
        INTEGER, DIMENSION(:), POINTER        :: isublist
        REAL(mk), DIMENSION(:,:,:,:  ), POINTER :: phi0,phirhs,phiprev
        INTEGER                               :: nsublist
        INTEGER, DIMENSION(:,:), POINTER      :: ndata
        INTEGER                               :: topoid,meshid
        REAL(mk), DIMENSION(:), POINTER       :: min_phys, max_phys
        INTEGER, DIMENSION(6)                 :: orgbcdef
        INTEGER                               :: s2didx, mpi_prec
        TYPE(ppm_t_topo),      POINTER        :: topo
        TYPE(ppm_t_equi_mesh), POINTER        :: mesh
        !-----------------------------------------------------
        !  standard stuff
        !-----------------------------------------------------
        INTEGER                               :: isub,isubl,i,j,k,maptype,istep,iopt
        INTEGER                               :: ldl(4), ldu(4), ndata_max(3)
        REAL(mk)                              :: len_phys(3), dx(3)
        REAL(mk) :: t0, res, gres,tau
        CHARACTER(len=256)                    :: msg

        CALL substart('ppm_hamjac_reinit_russo_3d',t0,info)

        IF(PRESENT(indx)) THEN
            s2didx=indx
        ELSE
            s2didx= -1
        END IF

#ifdef __MPI        
        IF (ppm_kind.EQ.ppm_kind_single) THEN
           MPI_PREC = MPI_REAL
        ELSE
           MPI_PREC = MPI_DOUBLE_PRECISION
        ENDIF
#endif
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
        
        ! timestep
        tau = 0.25_mk*MINVAL(dx)
        
        !-----------------------------------------------------
        !  allocate temporary storage
        !-----------------------------------------------------
        ldl(1:3) = 1 - ghostsize(1:3); ldl(4) = 1
        ndata_max(1) = MAXVAL(ndata(1,1:nsublist))
        ndata_max(2) = MAXVAL(ndata(2,1:nsublist))
        ndata_max(3) = MAXVAL(ndata(3,1:nsublist))
        ldu(1)   = ndata_max(1) + ghostsize(1)
        ldu(2)   = ndata_max(2) + ghostsize(2)
        ldu(3)   = ndata_max(3) + ghostsize(3)
        ldu(4)   = nsublist
        iopt     = ppm_param_alloc_fit
        CALL ppm_alloc(phi0,ldl,ldu,iopt,info)
        CALL ppm_alloc(phirhs,ldl,ldu,iopt,info)
        CALL ppm_alloc(phiprev,ldl,ldu,iopt,info)
        IF(info.NE.0) THEN
           info = ppm_error_fatal
           CALL ppm_error(ppm_err_alloc,'ppm_hamjac_reinit_russo_3d', &
                &        'temp storage for hamjac',__LINE__,info)
           GOTO 9999
        END IF
        
        ! fill phi0 with initial condition
        DO isub=1,nsublist
           isubl = isublist(isub)
           DO k=1,ndata(3,isubl)
           DO j=1,ndata(2,isubl)
              DO i=1,ndata(1,isubl)
                phi0(i,j,k,isub) = phi(i,j,k,isub)
                phiprev(i,j,k,isub) = phi(i,j,k,isub)
                phigrad(1:3,i,j,k,isub) = 0.0_mk
              END DO
            END DO
            END DO
        END DO
        
        CALL ppm_map_field_ghost_get(topo_id,mesh_id,ghostsize,info)
        CALL ppm_map_field_push(topo_id,mesh_id,phi0,info)
        CALL ppm_map_field_send(info)
        CALL ppm_map_field_pop(topo_id,mesh_id,phi0,ghostsize,info)
                
        !-----------------------------------------------------
        ! start time integration loop: using TVD RK3
        !-----------------------------------------------------
        
        DO istep=1,maxstep
           !--- map the ghosts
           CALL ppm_map_field_ghost_get(topo_id,mesh_id,ghostsize,info)
           CALL ppm_map_field_push(topo_id,mesh_id,phi,info)
           CALL ppm_map_field_send(info)
           CALL ppm_map_field_pop(topo_id,mesh_id,phi,ghostsize,info)
           
           ! TVDRK3 substep A
           
           CALL ppm_hamjac_reinit_russo_step(phi,phi0,phirhs,phigrad,res,s2didx,topo_id,mesh_id&
                &,                  ghostsize,info)
                
           DO isub=1,nsublist
              isubl = isublist(isub)              
              DO k=1,ndata(3,isubl); DO j=1,ndata(2,isubl);DO i=1,ndata(1,isubl)
                 phi(i,j,k,isub) = phi(i,j,k,isub) + tau*phirhs(i,j,k,isub)
              END DO; END DO; END DO
           END DO
           
           !--- map the ghosts
           CALL ppm_map_field_ghost_get(topo_id,mesh_id,ghostsize,info)
           CALL ppm_map_field_push(topo_id,mesh_id,phi,info)
           CALL ppm_map_field_send(info)
           CALL ppm_map_field_pop(topo_id,mesh_id,phi,ghostsize,info)
           
           ! TVDRK3 substep B
           
           CALL ppm_hamjac_reinit_russo_step(phi,phi0,phirhs,phigrad,res,s2didx,topo_id,mesh_id&
                &,                  ghostsize,info)
                
           DO isub=1,nsublist
              isubl = isublist(isub)
              DO k=1,ndata(3,isubl); DO j=1,ndata(2,isubl);DO i=1,ndata(1,isubl)
                 phi(i,j,k,isub) = phi(i,j,k,isub) + tau*phirhs(i,j,k,isub)
                 phi(i,j,k,isub) = 0.25_mk*phiprev(i,j,k,isub) + 0.75_mk*phi(i,j,k,isub)
              END DO; END DO; END DO
           END DO
           
           !--- map the ghosts
           CALL ppm_map_field_ghost_get(topo_id,mesh_id,ghostsize,info)
           CALL ppm_map_field_push(topo_id,mesh_id,phi,info)
           CALL ppm_map_field_send(info)
           CALL ppm_map_field_pop(topo_id,mesh_id,phi,ghostsize,info)
           
           ! TVDRK3 substep C
           
           CALL ppm_hamjac_reinit_russo_step(phi,phi0,phirhs,phigrad,res,s2didx,topo_id,mesh_id&
                &,                  ghostsize,info)
                
           DO isub=1,nsublist
              isubl = isublist(isub)
              DO k=1,ndata(3,isubl); DO j=1,ndata(2,isubl);DO i=1,ndata(1,isubl)
                 phi(i,j,k,isub) = phi(i,j,k,isub) + tau*phirhs(i,j,k,isub)
                 phi(i,j,k,isub) = (2.0_mk*phiprev(i,j,k,isub) + phi(i,j,k,isub))/3.0_mk
                 phiprev(i,j,k,isub) = phi(i,j,k,isub)
              END DO; END DO; END DO
           END DO           
#ifdef __MPI
           CALL MPI_AllReduce(res,gres,1,mpi_prec,MPI_MIN,0,ppm_comm,info)
           res = gres
#endif   
           IF(res.LT.tol) GOTO 666
           
           IF(ppm_rank.EQ.0.AND.MOD(istep,5).EQ.0) THEN
                WRITE(msg,*) 'iteration #',istep,' res=',res
                CALL ppm_write(ppm_Rank,'ppm_hamjac',msg,info)
           END IF
           
           !IF(res.LT.tol) GOTO 666 ! does not work in parallel with ghosting, because all 
                                    ! processors need to be in loop for ghosting to work
        END DO

        !info = ppm_error_warning
        !CALL ppm_error(ppm_err_converge,'ppm_hamjac_reinit_russo_3d', &
        !     &         'failed to reach target residual',__LINE__,info)

666     CONTINUE
       !--- map the ghosts for last time
       CALL ppm_map_field_ghost_get(topo_id,mesh_id,ghostsize,info)
       CALL ppm_map_field_push(topo_id,mesh_id,phi,info)
       CALL ppm_map_field_send(info)
       CALL ppm_map_field_pop(topo_id,mesh_id,phi,ghostsize,info)
   
       IF(ppm_rank.EQ.0) THEN
            WRITE(msg,*) 'ended after iteration #',istep,' res=',res
            CALL ppm_write(ppm_Rank,'ppm_hamjac',msg,info)
       END IF

        iopt = ppm_param_dealloc
        CALL ppm_alloc(phi0,ldl,ldu,iopt,info)
        CALL ppm_alloc(phirhs,ldl,ldu,iopt,info)
        CALL ppm_alloc(phiprev,ldl,ldu,iopt,info)
        IF(info.NE.0) THEN
           info = ppm_error_error
           CALL ppm_error(ppm_err_dealloc,'ppm_hamjac_reinit_russo_3d', &
                &        'temp storage for hamjac not freed',__LINE__,info)
           GOTO 9999
        END IF


9999    CONTINUE
        
        CALL substop('ppm_hamjac_reinit_russo_3d',t0,info)
        
        
#if   __KIND == __SINGLE_PRECISION
      END SUBROUTINE ppm_hamjac_reinit_russo_3ds 
#elif __KIND == __DOUBLE_PRECISION
      END SUBROUTINE ppm_hamjac_reinit_russo_3dd 
#endif


        
           

        
        

