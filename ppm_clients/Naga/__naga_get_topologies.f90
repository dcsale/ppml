!-------------------------------------------------------------------------------
! Subroutine : naga_get_topologies.f90
!-------------------------------------------------------------------------------
! _ __
! / |/ /__ ____ ____ _
! / / _ `/ _ `/ _ `/
! /_/|_/\_,_/\_, /\_,_/
! /___/
! A PPM Vortex-In-Cell client - 2011 (c)
! Naga is distributed under the GNU Lesser General Public License version 3
! Johannes Tophøj Rasmussen, johannes.tophoj@gmail.com
! Technical University of Denmark
! Department of Mechanical Engineering
!
!
! This routines gets topologies for all levels and stores them in a vector
!-------------------------------------------------------------------------------
MODULE naga_mod_get_topologies
IMPLICIT NONE
INTERFACE naga_get_topologies
  MODULE PROCEDURE naga_get_topologies
END INTERFACE
CONTAINS
SUBROUTINE naga_get_topologies(psetting,topoarray,info)
USE naga_mod_globals
USE naga_mod_say
USE ppm_module_topo_get
IMPLICIT NONE
!-----------------------------------------------------------------------------
! Arguments
!-----------------------------------------------------------------------------
TYPE(patch_setup),DIMENSION(:,:),POINTER :: psetting
TYPE(ppm_t_topo),DIMENSION(:), POINTER :: topoarray
INTEGER, INTENT(INOUT) :: info
!-----------------------------------------------------------------------------
! Local variables
!-----------------------------------------------------------------------------
REAL(MK) :: t0
INTEGER :: ilevel
INTEGER :: ipatch !1 until supported by topo
TYPE(ppm_t_topo),POINTER :: topo => NULL() !temporary topology
!-----------------------------------------------------------------------------
! Initialise routine
!-----------------------------------------------------------------------------
CALL substart('naga_get_topologies',t0,info)
!-----------------------------------------------------------------------------
! Allocate topology vector for explicit storage of the topology data
!-----------------------------------------------------------------------------
ALLOCATE(topoarray(nlevels))
!-----------------------------------------------------------------------------
! Loop over all resolution levels
!-----------------------------------------------------------------------------
DO ilevel=1,nlevels
  ipatch = 1
  !----------------------------------------------------------------------------
  ! Get the topology (struct) from the topoid and save it to topoarray vector
  !----------------------------------------------------------------------------
  CALL ppm_topo_get(psetting(ilevel,ipatch)%topoid,topo,info)
  IF (info .NE. 0) THEN
    CALL naga_say(rank, 'naga_get_topologies','Failed to get topology.')
    GOTO 9999
  ENDIF
  topoarray(ilevel) = topo
ENDDO
!----------------------------------------------------------------------------
! Return
!----------------------------------------------------------------------------
 9999 CONTINUE
CALL substop('naga_get_topologies',t0,info)
RETURN
END SUBROUTINE naga_get_topologies
END MODULE naga_mod_get_topologies
