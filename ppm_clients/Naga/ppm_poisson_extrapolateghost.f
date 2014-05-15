      !-------------------------------------------------------------------------
      !  Subroutine   : ppm_poisson_extrapolateghost.f90
      !-------------------------------------------------------------------------
      ! Copyright (c) 2010 CSE Lab (ETH Zurich), MOSAIC Group (ETH Zurich),
      !                    Center for Fluid Dynamics (DTU)
      !
      !-------------------------------------------------------------------------
      SUBROUTINE __ROUTINE(topoid,meshid,field,nextra,nbase,gstw,info)
      !!! This routine extrapolates field values of field with topology and mesh
      !!! id topoid,meshid respectively into the ghost layer of width gstw.
      !!! nbase points are used to extrapolate into nextra points.
      !!!
      !!! [NOTE]
      !!! Presently extrapolation can only be done to 1 or two points into the
      !!! ghostlayer always based on 4 points (fourth order spatial convergence)
      !!! A general nbase,nextra extrapolation can be implemented vi solution
      !!! of a small linear system of equations. This has not been done.
      !!! This routine is in need of loop unrollling in particular for
      !!! typical choices of nbase,nextra pairs. Extrapolation is necessary for
      !!! e.g. freespace FD curl of stream function

      USE ppm_module_topo_get


      IMPLICIT NONE
      !-------------------------------------------------------------------------
      ! Arguments
      !-------------------------------------------------------------------------
      INTEGER, INTENT(IN)                                         :: topoid
      !!! Topology id of the field
      INTEGER, INTENT(IN)                                         :: meshid
      !!! Mesh id of the field
      REAL(__PREC),DIMENSION(:,:,:,:,:),POINTER                   :: field
      !!! Field to extrapolate ghost layer
      INTEGER, INTENT(IN)                                         :: nextra
      !!! Number of points to extrapolate into the ghostlayer
      INTEGER, INTENT(IN)                                         :: nbase
      !!! Number of points to base the extrapolation on
      INTEGER,DIMENSION(__DIM),INTENT(IN)                         :: gstw
      !!! Width of the ghotslayer
      INTEGER, INTENT(OUT)                                        :: info
      !!! Return state

      !-------------------------------------------------------------------------
      ! Local variables
      !-------------------------------------------------------------------------
      INTEGER,PARAMETER                    :: MK = __PREC
      REAL(__PREC)                         :: t0
      TYPE(ppm_t_topo),POINTER             :: topology
      TYPE(ppm_t_equi_mesh)                :: mesh
      INTEGER                              :: isub,isubl
      INTEGER                              :: i,j,k,iextra,ibase
      REAL(__PREC),DIMENSION(:,:),POINTER  :: coeff
      REAL(__PREC),DIMENSION(__DIM)        :: tmpbuf


      !-------------------------------------------------------------------------
      ! Initialise routine
      !-------------------------------------------------------------------------
      CALL substart('ppm_poisson_extrapolateghost',t0,info)


      !-------------------------------------------------------------------------
      ! Compare the number of points to extrapolate to the ghost layer width
      !-------------------------------------------------------------------------
      IF (nextra .GT. gstw(1) .OR. &
        & nextra .GT. gstw(2) .OR. &
        & nextra .GT. gstw(3)) THEN
         CALL ppm_write(ppm_rank,'ppm_poisson_extrapolateghost',&
         & 'The points to extrapolate exceeds the ghost layer.',info)
         info = -1
         GOTO 9999
      ENDIF

      !-------------------------------------------------------------------------
      ! Determine weights
      !-------------------------------------------------------------------------
      ALLOCATE(coeff(nbase,nextra))
      IF (nbase .EQ. 4) THEN
        IF (nextra .GE. 1) THEN
          coeff(:,1) = (/4.0_MK,-6.0_MK,4.0_MK,-1.0_MK/)
        ENDIF
        IF (nextra .GE. 2) THEN
          coeff(:,2) = (/10.0_MK,-20.0_MK,15.0_MK,-4.0_MK/)
        ENDIF
        IF (nextra .GE. 3) THEN
          CALL ppm_write(ppm_rank,'ppm_poisson_extrapolateghost',&
          & 'Extrapolation to more than two points has not been implemented.',info)
          info = -1
          GOTO 9999
        ENDIF
      ELSE
        CALL ppm_write(ppm_rank,'ppm_poisson_extrapolateghost',&
        & 'Only extrapolation based on 4 points has been implemented.',info)
        info = -1
        GOTO 9999
      ENDIF


      !-------------------------------------------------------------------------
      ! Get topology
      !-------------------------------------------------------------------------
      CALL ppm_topo_get(topoid,topology,info)
      IF (info .NE. 0) THEN
         CALL ppm_write(ppm_rank,'ppm_poisson_extrapolateghost',&
         & 'Failed to get topology.',isub)
         GOTO 9999
      ENDIF
      mesh = topology%mesh(meshid)

      !-------------------------------------------------------------------------
      ! Extrapolate field into ghost layer
      ! The indicies of subs_bc represent: 
      ! west,east(x),south,north(y),bottom,top(z)
      !@ Some more unrolling here would be nice
      !-------------------------------------------------------------------------
      DO isub=1,topology%nsublist
         isubl=topology%isublist(isub)
         !West (-x)
         IF (topology%subs_bc(1,isubl) .EQ. 1) THEN
            DO k=1-gstw(3),mesh%nnodes(3,isubl)+gstw(3)
               DO j=1-gstw(2),mesh%nnodes(2,isubl)+gstw(2)
                  !DO i=1-gstw(1),mesh%nnodes(1,isubl)+gstw(1)
                  i = 1
                  DO iextra=1,nextra
                     tmpbuf = 0.0_MK
                     DO ibase=0,nbase-1
                        tmpbuf(1) = tmpbuf(1) + &
                         & coeff(ibase+1,iextra)*field(1,i+ibase,j,k,isub)
                        tmpbuf(2) = tmpbuf(2) + &
                         & coeff(ibase+1,iextra)*field(2,i+ibase,j,k,isub)
                        tmpbuf(3) = tmpbuf(3) + &
                         & coeff(ibase+1,iextra)*field(3,i+ibase,j,k,isub)
                     END DO !ibase
                     field(1,i-iextra,j,k,isub) = tmpbuf(1)
                     field(2,i-iextra,j,k,isub) = tmpbuf(2)
                     field(3,i-iextra,j,k,isub) = tmpbuf(3)
                  END DO !iextra
               ENDDO !j
            ENDDO !k
         ENDIF
         !East (+x)
         IF (topology%subs_bc(2,isubl) .EQ. 1) THEN
            DO k=1-gstw(3),mesh%nnodes(3,isubl)+gstw(3)
               DO j=1-gstw(2),mesh%nnodes(2,isubl)+gstw(2)
                  !DO i=1-gstw(1),mesh%nnodes(1,isubl)+gstw(1)
                  i = mesh%nnodes(1,isubl)
                  DO iextra=1,nextra
                     tmpbuf = 0.0_MK
                     DO ibase=0,nbase-1
                        tmpbuf(1) = tmpbuf(1) + &
                         & coeff(ibase+1,iextra)*field(1,i-ibase,j,k,isub)
                        tmpbuf(2) = tmpbuf(2) + &
                         & coeff(ibase+1,iextra)*field(2,i-ibase,j,k,isub)
                        tmpbuf(3) = tmpbuf(3) + &
                         & coeff(ibase+1,iextra)*field(3,i-ibase,j,k,isub)
                     END DO !ibase
                     field(1,i+iextra,j,k,isub) = tmpbuf(1)
                     field(2,i+iextra,j,k,isub) = tmpbuf(2)
                     field(3,i+iextra,j,k,isub) = tmpbuf(3)
                  END DO !iextra
               ENDDO !j
            ENDDO !k
         ENDIF
         !South (-y)
         IF (topology%subs_bc(3,isubl) .EQ. 1) THEN
            DO k=1-gstw(3),mesh%nnodes(3,isubl)+gstw(3)
               !DO j=1-gstw(2),mesh%nnodes(2,isubl)+gstw(2)
               DO i=1-gstw(1),mesh%nnodes(1,isubl)+gstw(1)
                  j = 1
                  DO iextra=1,nextra
                     tmpbuf = 0.0_MK
                     DO ibase=0,nbase-1
                        tmpbuf(1) = tmpbuf(1) + &
                         & coeff(ibase+1,iextra)*field(1,i,j+ibase,k,isub)
                        tmpbuf(2) = tmpbuf(2) + &
                         & coeff(ibase+1,iextra)*field(2,i,j+ibase,k,isub)
                        tmpbuf(3) = tmpbuf(3) + &
                         & coeff(ibase+1,iextra)*field(3,i,j+ibase,k,isub)
                     END DO !ibase
                     field(1,i,j-iextra,k,isub) = tmpbuf(1)
                     field(2,i,j-iextra,k,isub) = tmpbuf(2)
                     field(3,i,j-iextra,k,isub) = tmpbuf(3)
                  END DO !iextra
               ENDDO !i
            ENDDO !k
         ENDIF
         !North (+y)
         IF (topology%subs_bc(4,isubl) .EQ. 1) THEN
            DO k=1-gstw(3),mesh%nnodes(3,isubl)+gstw(3)
               !DO j=1-gstw(2),mesh%nnodes(2,isubl)+gstw(2)
               DO i=1-gstw(1),mesh%nnodes(1,isubl)+gstw(1)
                  j = mesh%nnodes(2,isubl)
                  DO iextra=1,nextra
                     tmpbuf = 0.0_MK
                     DO ibase=0,nbase-1
                        tmpbuf(1) = tmpbuf(1) + &
                         & coeff(ibase+1,iextra)*field(1,i,j-ibase,k,isub)
                        tmpbuf(2) = tmpbuf(2) + &
                         & coeff(ibase+1,iextra)*field(2,i,j-ibase,k,isub)
                        tmpbuf(3) = tmpbuf(3) + &
                         & coeff(ibase+1,iextra)*field(3,i,j-ibase,k,isub)
                     END DO !ibase
                     field(1,i,j+iextra,k,isub) = tmpbuf(1)
                     field(2,i,j+iextra,k,isub) = tmpbuf(2)
                     field(3,i,j+iextra,k,isub) = tmpbuf(3)
                  END DO !iextra
               ENDDO !i
            ENDDO !k
         ENDIF
         !Bottom (-z)
         IF (topology%subs_bc(5,isubl) .EQ. 1) THEN
            !DO k=1-gstw(3),mesh%nnodes(3,isubl)+gstw(3)
            DO j=1-gstw(2),mesh%nnodes(2,isubl)+gstw(2)
               DO i=1-gstw(1),mesh%nnodes(1,isubl)+gstw(1)
                  k = 1
                  DO iextra=1,nextra
                     tmpbuf = 0.0_MK
                     DO ibase=0,nbase-1
                        tmpbuf(1) = tmpbuf(1) + &
                         & coeff(ibase+1,iextra)*field(1,i,j,k+ibase,isub)
                        tmpbuf(2) = tmpbuf(2) + &
                         & coeff(ibase+1,iextra)*field(2,i,j,k+ibase,isub)
                        tmpbuf(3) = tmpbuf(3) + &
                         & coeff(ibase+1,iextra)*field(3,i,j,k+ibase,isub)
                     END DO !ibase
                     field(1,i,j,k-iextra,isub) = tmpbuf(1)
                     field(2,i,j,k-iextra,isub) = tmpbuf(2)
                     field(3,i,j,k-iextra,isub) = tmpbuf(3)
                  END DO !iextra
               ENDDO !i
            ENDDO !j
         ENDIF
         !Top (+z)
         IF (topology%subs_bc(6,isubl) .EQ. 1) THEN
            !DO k=1-gstw(3),mesh%nnodes(3,isubl)+gstw(3)
            DO j=1-gstw(2),mesh%nnodes(2,isubl)+gstw(2)
               DO i=1-gstw(1),mesh%nnodes(1,isubl)+gstw(1)
                  k = mesh%nnodes(3,isubl)
                  DO iextra=1,nextra
                     tmpbuf = 0.0_MK
                     DO ibase=0,nbase-1
                        tmpbuf(1) = tmpbuf(1) + &
                         & coeff(ibase+1,iextra)*field(1,i,j,k-ibase,isub)
                        tmpbuf(2) = tmpbuf(2) + &
                         & coeff(ibase+1,iextra)*field(2,i,j,k-ibase,isub)
                        tmpbuf(3) = tmpbuf(3) + &
                         & coeff(ibase+1,iextra)*field(3,i,j,k-ibase,isub)
                     END DO !ibase
                     field(1,i,j,k+iextra,isub) = tmpbuf(1)
                     field(2,i,j,k+iextra,isub) = tmpbuf(2)
                     field(3,i,j,k+iextra,isub) = tmpbuf(3)
                  END DO !iextra
               ENDDO !i
            ENDDO !j
         ENDIF
      ENDDO !isub


 9999 CONTINUE
      CALL substop('ppm_poisson_extrapolateghost',t0,info)
      RETURN

      END SUBROUTINE __ROUTINE

