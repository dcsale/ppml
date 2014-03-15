SUBROUTINE lj_newtopo(info)
!======================================================================!
! creates new topology, maps particles, creates/maps ghosts, and builds
! cell lists
!======================================================================!
  USE ppm_module_core
  USE ppm_module_map
  USE lj_module_global
  IMPLICIT NONE
  INCLUDE 'mpif.h'
  ! arguments
  INTEGER,INTENT(OUT) :: info
  ! local variables
  INTEGER :: partindex,cellindex,istart,iend
  INTEGER :: icell,jcell,kcell,ip,isub
  CHARACTER(LEN=256) :: myformat
  !====================================================================!
  ! set up a topology
  ALLOCATE(bcdef(2*ndim),STAT=info)
  IF (info .NE. 0) THEN
     WRITE(*,*)'ERROR in lj_newtopo: Allocation failed.'
     GOTO 9999
  ENDIF
  decomp = ppm_param_decomp_bisection
  assig = ppm_param_assign_internal
  bcdef = ppm_param_bcdef_periodic ! make all six boundaries periodic
  ghostsize = cutoff + skin
  topo_id = 0
  CALL ppm_mktopo(topo_id,decomp,assig,mindom,maxdom,bcdef,ghostsize, &
       cost,info,minsub,maxsub,nsubs,sub2proc)
  IF (info .NE. 0) THEN
     WRITE(*,*)'ERROR in lj_newtopo: ppm_mktopo failed.'
     GOTO 9999
  ENDIF
  !====================================================================!
  ! global mapping of the data onto the topology
  CALL ppm_map_part_global(topo_id,xp,Npart,info)
  IF (info .NE. 0) THEN
     WRITE(*,*)'ERROR in lj_newtopo: ppm_map_part (global) failed.'
     GOTO 9999
  ENDIF
  CALL ppm_map_part_push(vp,ndim,Npart,info)
  IF (info .NE. 0) THEN
     WRITE(*,*)'ERROR in lj_newtopo: ppm_map_part (push) failed.'
     GOTO 9999
  ENDIF
  CALL ppm_map_part_send(Npart,Mpart,info)
  IF (info .NE. 0) THEN
     WRITE(*,*)'ERROR in lj_newtopo: ppm_map_part (send) failed.'
     GOTO 9999
  ENDIF
  CALL ppm_map_part_pop(vp,ndim,Npart,Mpart,info)
  IF (info .NE. 0) THEN
     WRITE(*,*)'ERROR in lj_newtopo: ppm_map_part (pop) failed.'
     GOTO 9999
  ENDIF
  CALL ppm_map_part_pop(xp,ndim,Npart,Mpart,info)
  IF (info .NE. 0) THEN
     WRITE(*,*)'ERROR in lj_newtopo: ppm_map_part (pop) failed.'
     GOTO 9999
  ENDIF
  Npart = Mpart
  !====================================================================!
  ! get the ghosts
  CALL ppm_map_part_ghost_get(topo_id,xp,ndim,Npart,isymm,ghostsize,info)
  IF (info .NE. 0) THEN
     WRITE(*,*)'ERROR in lj_newtopo: ppm_map_part_ghost (get) failed.'
     GOTO 9999
  ENDIF
  CALL ppm_map_part_push(vp,ndim,Npart,info)
  IF (info .NE. 0) THEN
     WRITE(*,*)'ERROR in lj_newtopo: ppm_map_part_ghost (get) failed.'
     GOTO 9999
  ENDIF
  CALL ppm_map_part_send(Npart,Mpart,info)
  IF (info .NE. 0) THEN
     WRITE(*,*)'ERROR in lj_newtopo: ppm_map_part_ghost (send) failed.'
     GOTO 9999
  ENDIF
  CALL ppm_map_part_pop(vp,ndim,Npart,Mpart,info)
  IF (info .NE. 0) THEN
     WRITE(*,*)'ERROR in lj_newtopo: ppm_map_part_ghost (pop) failed.'
     GOTO 9999
  ENDIF
  CALL ppm_map_part_pop(xp,ndim,Npart,Mpart,info)
  IF (info .NE. 0) THEN
     WRITE(*,*)'ERROR in lj_newtopo: ppm_map_part_ghost (pop) failed.'
     GOTO 9999
  ENDIF
  WRITE(myformat,'(A,I1,A,I1,A,I1,A)'),'(A,I',CEILING(LOG10(REAL(nproc))),',A,I', &
       INT(LOG10(REAL(MAX(1,Npart))))+1,',A,I',INT(LOG10(REAL(MAX(1,Mpart-Npart))))+1,',A)'
  WRITE(*,myformat)'[',rank, &
       '] (lj_newtopo) : number of real (ghost) particles  = ',Npart, &
       ' (',Mpart-Npart,')'
  !====================================================================!
  ! build cell lists
  lsymm = .FALSE.
  IF (isymm .GT. 0) lsymm = .TRUE.
  cellsize = ghostsize
  CALL ppm_neighlist_clist(topo_id,xp,Mpart,cellsize,lsymm,clist,info)
  IF (info .NE. 0) THEN
     WRITE(*,*)'ERROR in lj_newtopo: ppm_neighlist_clist failed.'
     GOTO 9999
  ENDIF
9999 CONTINUE ! jump here upon error
END SUBROUTINE lj_newtopo
