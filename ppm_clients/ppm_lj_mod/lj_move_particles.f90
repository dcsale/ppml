SUBROUTINE lj_move_particles(info)
!======================================================================!
! changes particle positions and velocities following the forces 
! resulting from the Lennard-Jones potential
!======================================================================!
  USE lj_module_global
  USE ppm_module_core
  USE ppm_module_user_numerics
  USE ppm_module_map
  IMPLICIT NONE
  ! arguments
  INTEGER, INTENT(OUT)                          :: info
  ! local variables
  INTEGER                                       :: ip
  REAL(prec)                                    :: cutoff2

  !====================================================================!
  ! initialize
  info = 0
  comp_accel_arg(1) = sigma
  comp_accel_arg(2) = epsil
  comp_accel_arg(3) = mass
  cutoff2 = cutoff**2
  IF (ndim .EQ. 2) THEN
     IF (lsymm) THEN
        DO ip=1,Mpart
           accp(1,ip) = 0._prec
           accp(2,ip) = 0._prec
        ENDDO
     ELSE
        DO ip=1,Npart
           accp(1,ip) = 0._prec
           accp(2,ip) = 0._prec
        ENDDO
     ENDIF
  ELSE
     IF (lsymm) THEN
        DO ip=1,Mpart
           accp(1,ip) = 0._prec
           accp(2,ip) = 0._prec
           accp(3,ip) = 0._prec
        ENDDO
     ELSE
        DO ip=1,Npart
           accp(1,ip) = 0._prec
           accp(2,ip) = 0._prec
           accp(3,ip) = 0._prec
        ENDDO
     ENDIF
  ENDIF

  !====================================================================!
  ! compute particle interactions (=> acceleration) using cell lists
  CALL ppm_comp_pp_cell(topo_id,xp,Mpart,xp,ndim,lsymm,comp_accel,  &
       comp_accel_arg,clist,cutoff2,accp,info)
  IF (info .NE. 0) THEN
     WRITE(*,*)'ERROR in lj_move_particles: ppm_comp_pp_cell failed.'
     GOTO 9999
  ENDIF
  
  IF (lsymm) THEN
  !====================================================================!
  ! put contributions to acceleration from being ghost onto 
  ! acceleration of the real particle
     CALL ppm_map_part_ghost_put(topo_id,info)
     IF (info .NE. 0) THEN
        WRITE(*,*)'ERROR in lj_newtopo: ppm_map_part_ghost (put) failed.'
        GOTO 9999
     ENDIF
     CALL ppm_map_part_push(accp,ndim,Npart,info)
     IF (info .NE. 0) THEN
        WRITE(*,*)'ERROR in lj_newtopo: ppm_map_part_ghost (push) failed.'
        GOTO 9999
     ENDIF
     CALL ppm_map_part_send(Npart,Mpart,info)
     IF (info .NE. 0) THEN
        WRITE(*,*)'ERROR in lj_newtopo: ppm_map_part_ghost (send) failed.'
        GOTO 9999
     ENDIF
     CALL ppm_map_part_pop(accp,ndim,Npart,Mpart,info)
     IF (info .NE. 0) THEN
        WRITE(*,*)'ERROR in lj_newtopo: ppm_map_part_ghost (pop) failed.'
        GOTO 9999
     ENDIF
  ENDIF

  !====================================================================!
  ! change particle positions and velocities (time integration)
  CALL lj_change_pos_vel(info)
  IF (info .NE. 0) THEN
     WRITE(*,*)'ERROR in lj_move_particles: lj_change_pos_vel failed.'
     GOTO 9999
  ENDIF

9999 CONTINUE ! jump here upon error

CONTAINS
  
  FUNCTION comp_accel(d2,args)
    ! returns the acceleration divided by the distance vector(!).
    ! args(1) is sigma, args(2) is epsilon of the LJ potential,
    ! args(3) is the mass of a particle.
    ! d2 is the squared distance between the two particles under 
    ! consideration.]
    USE lj_module_global
    IMPLICIT NONE
    REAL(prec)                            :: comp_accel
    REAL(prec),             INTENT(IN)    :: d2
    REAL(prec),DIMENSION(:),INTENT(IN)    :: args
    REAL(prec)                            :: fraction6

    fraction6 = (args(1)**2/d2)**3
    comp_accel = -24._prec*args(2)*(fraction6*(2._prec*fraction6 - 1._prec))/(d2*args(3))

  END FUNCTION comp_accel

END SUBROUTINE lj_move_particles
