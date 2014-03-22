SUBROUTINE lj_remap(info)
!======================================================================!
! partial mapping, creation of ghosts, creation of cell lists
!======================================================================!
  USE ppm_module_core
  USE ppm_module_map
  USE lj_module_global
  IMPLICIT NONE
  INCLUDE 'mpif.h'
  ! arguments
  INTEGER,INTENT(OUT)                   :: info
  ! local variables

  !====================================================================!
  ! initialize some general variables
  info = 0 ! change if error occurs

  !====================================================================!
  ! partial mapping of the data onto the topology
  CALL ppm_map_part_partial(topo_id,xp,Npart,info)
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

  !====================================================================!
  ! build cell lists
  CALL ppm_neighlist_clist(topo_id,xp,Mpart,cellsize,lsymm,clist,info)
  IF (info .NE. 0) THEN
     WRITE(*,*)'ERROR in lj_newtopo: ppm_neighlist_clist failed.'
     GOTO 9999
  ENDIF

  DEALLOCATE(accp)
  IF (lsymm) THEN
     ALLOCATE(accp(ndim,Mpart))
     IF (info .NE. 0) THEN
        WRITE(*,*)'ERROR in lj_init: Allocation failed.'
        GOTO 9999
     ENDIF
  ELSE
     ALLOCATE(accp(ndim,Npart))
     IF (info .NE. 0) THEN
        WRITE(*,*)'ERROR in lj_init: Allocation failed.'
        GOTO 9999
     ENDIF
  ENDIF

9999  CONTINUE ! jump here upon error

END SUBROUTINE lj_remap
