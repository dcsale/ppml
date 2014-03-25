      !-------------------------------------------------------------------------
      !     Subroutine   :                 ppm_hamjac_reinit_ref_3d
      !-------------------------------------------------------------------------
      !     
      !     Purpose      : Solve Hamilton-Jacobi for Gowas reinit in ref spc
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
      !     $Log: ppm_hamjac_reinit_ref_3d.f,v $
      !     Revision 1.1.1.1  2006/07/25 15:18:19  menahel
      !     initial import
      !
      !     Revision 1.1  2005/07/25 00:34:03  ivos
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
      SUBROUTINE ppm_hamjac_reinit_ref_3ds (phi, chi, trgt, tol, maxstep, &
           &                     topo_id, mesh_id, ghostsize, info)
#elif __KIND == __DOUBLE_PRECISION
      SUBROUTINE ppm_hamjac_reinit_ref_3dd (phi, chi, trgt, tol, maxstep, &
           &                     topo_id, mesh_id, ghostsize, info)
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
        USE ppm_module_typedef
        USE ppm_module_map_field
        USE ppm_module_map_field_ghost
        IMPLICIT NONE

#if    __KIND == __SINGLE_PRECISION
        INTEGER, PARAMETER :: MK = ppm_kind_single
#elif  __KIND == __DOUBLE_PRECISION       
        INTEGER, PARAMETER :: MK = ppm_kind_double
#endif

        !-----------------------------------------------------
        !  Arguments
        !-----------------------------------------------------
        REAL(MK), DIMENSION(:,:,:,:), POINTER :: phi
        REAL(mk), DIMENSION(:,:,:,:,:), POINTER :: chi
        INTEGER, INTENT(in)                   :: topo_id, mesh_id
        INTEGER, DIMENSION(3), INTENT(in)     :: ghostsize
        INTEGER, INTENT(inout)                :: info
        INTEGER, INTENT(in)                   :: maxstep
        REAL(mk), INTENT(in)                  :: tol, trgt

        !-----------------------------------------------------
        !  Aliases
        !-----------------------------------------------------
        INTEGER, DIMENSION(:), POINTER        :: isublist
        REAL(mk), DIMENSION(:,:,:,:), POINTER :: tphi
        INTEGER                               :: nsublist
        INTEGER, DIMENSION(:,:), POINTER      :: ndata
        INTEGER                               :: topoid,meshid
        REAL(MK), DIMENSION(:), POINTER       :: min_phys, max_phys
        TYPE(ppm_t_topo),      POINTER        :: topo
        TYPE(ppm_t_equi_mesh), POINTER        :: mesh
        
        !-----------------------------------------------------
        !  standard stuff
        !-----------------------------------------------------
        INTEGER                               :: isub,isubl,i,j,k,maptype,istep,iopt
        INTEGER                               :: ldl(4), ldu(4), ndata_max(3)
        REAL(mk)                              :: len_phys(3)
        REAL(mk) :: t0, res
        CHARACTER(LEN=ppm_char)               :: cbuf

        CALL substart('ppm_hamjac_reinit_ref_3d',t0,info)
        
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


        !-----------------------------------------------------
        !  RATIONALE Thu May 26 20:51:19 PDT 2005:
        !  loop ghostmap doit. easy.
        !-----------------------------------------------------


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
        CALL ppm_alloc(tphi,ldl,ldu,iopt,info)
        IF(info.NE.0) THEN
           info = ppm_error_fatal
           CALL ppm_error(ppm_err_alloc,'ppm_hamjac_reinit_ref_3d', &
                &        'temp storage for hamjac',__LINE__,info)
           GOTO 9999
        END IF

        !--- map the map
        CALL ppm_map_field_ghost_get(topo_id,mesh_id,ghostsize,info)
        CALL ppm_map_field_push(topo_id,mesh_id,chi,3,info)
        CALL ppm_map_field_send(info)
        CALL ppm_map_field_pop(topo_id,mesh_id,chi,3,ghostsize,info)
           

        !--- COMMENT Thu May 26 21:05:23 PDT 2005:  simple euler here, do TVD
        DO istep=1,maxstep
           !--- map the gowas
           CALL ppm_map_field_ghost_get(topo_id,mesh_id,ghostsize,info)
           CALL ppm_map_field_push(topo_id,mesh_id,phi,info)
           CALL ppm_map_field_send(info)
           CALL ppm_map_field_pop(topo_id,mesh_id,phi,ghostsize,info)

           CALL ppm_hamjac_reinit_step_ref(phi,chi,tphi,trgt,res,topo_id,   &
     &         mesh_id,ghostsize,info)

           ! IF (ppm_debug .GT. 0) THEN
           !TODO Uncomment WRITE statement
           !WRITE(cbuf,'(A,I4,A,ES)') 'Iteration ',istep,' Residual: ', res
           CALL ppm_write(ppm_rank,'ppm_hamjac_reinit_ref_3d',cbuf,info)
           ! ENDIF

           DO isub=1,nsublist
              isubl = isublist(isub)
              DO k=1,ndata(3,isubl);DO j=1,ndata(2,isubl);DO i=1,ndata(1,isubl)
                 phi(i,j,k,isub) = tphi(i,j,k,isub)
              END DO; END DO; END DO
           END DO
           IF(res.LT.tol) GOTO 666
        END DO

        info = ppm_error_warning
        CALL ppm_error(ppm_err_converge,'ppm_hamjac_reinit_ref_3d', &
             &         'failed to reach target residual',__LINE__,info)
        info = ppm_param_success

666     CONTINUE

        iopt = ppm_param_dealloc
        CALL ppm_alloc(tphi,ldl,ldu,iopt,info)
        IF(info.NE.0) THEN
           info = ppm_error_error
           CALL ppm_error(ppm_err_dealloc,'ppm_hamjac_reinit_ref_3d', &
                &        'temp storage for hamjac not freed',__LINE__,info)
           GOTO 9999
        END IF


9999    CONTINUE

#if   __KIND == __SINGLE_PRECISION
      END SUBROUTINE ppm_hamjac_reinit_ref_3ds 
#elif __KIND == __DOUBLE_PRECISION
      END SUBROUTINE ppm_hamjac_reinit_ref_3dd 
#endif


        
           

        
        

