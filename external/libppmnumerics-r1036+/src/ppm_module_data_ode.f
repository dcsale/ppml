
!-------------------------------------------------------------------------
      !  Module       :                 ppm_module_data_ode
      !-------------------------------------------------------------------------
      !
      !  Purpose      : ode module
      !
      !  Remarks      : 
      !
      !  References   : 
      !
      !  Revisions    :
      !-------------------------------------------------------------------------
      !  $Log: ppm_module_data_ode.f,v $
      !  Revision 1.1.1.1  2007/07/13 10:18:57  ivos
      !  CBL version of the PPM library
      !
      !  Revision 1.7  2005/08/05 09:26:58  ivos
      !  Added STS stuff (credits: MB)
      !
      !  Revision 1.6  2004/10/01 15:15:40  hiebers
      !  added Runge Kutta 4th order
      !
      !  Revision 1.5  2004/08/13 15:32:30  michaebe
      !  modified scheme information for midrk2
      !
      !  Revision 1.4  2004/08/12 12:53:40  michaebe
      !  inserted a hack for some compilers
      !
      !  Revision 1.3  2004/07/27 09:11:24  michaebe
      !  removed user accessible params and moved them to ppm_param.h
      !
      !  Revision 1.2  2004/07/26 11:39:23  michaebe
      !  removed PRIVATE statement as invalid after atomization
      !
      !  Revision 1.1  2004/07/26 11:32:31  michaebe
      !  Data module for the ode
      !
      !
      !-------------------------------------------------------------------------
      ! Copyright (c) 2010 CSE Lab (ETH Zurich), MOSAIC Group (ETH Zurich), 
      !                    Center for Fluid Dynamics (DTU)
      !
      !
      ! This file is part of the Parallel Particle Mesh Library (PPM).
      !
      ! PPM is free software: you can redistribute it and/or modify
      ! it under the terms of the GNU Lesser General Public License 
      ! as published by the Free Software Foundation, either 
      ! version 3 of the License, or (at your option) any later 
      ! version.
      !
      ! PPM is distributed in the hope that it will be useful,
      ! but WITHOUT ANY WARRANTY; without even the implied warranty of
      ! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
      ! GNU General Public License for more details.
      !
      ! You should have received a copy of the GNU General Public License
      ! and the GNU Lesser General Public License along with PPM. If not,
      ! see <http://www.gnu.org/licenses/>.
      !
      ! Parallel Particle Mesh Library (PPM)
      ! ETH Zurich
      ! CH-8092 Zurich, Switzerland

      !-------------------------------------------------------------------------

      MODULE ppm_module_data_ode
        !-----------------------------------------------------------------------
        !  Includes
        !-----------------------------------------------------------------------


        !-----------------------------------------------------------------------
        !  Time
        !-----------------------------------------------------------------------
        REAL(KIND(1.0D0)) :: t0

        !-----------------------------------------------------------------------
        ! implemented schemes
        !-----------------------------------------------------------------------
        ! have been moved to ppm_param.h
        
        !-----------------------------------------------------------------------
        ! scheme stuff
        ! _o : order
        ! _m : memory needs
        ! _s : number of stages
        ! _k : suitable kick off scheme
        !-----------------------------------------------------------------------
        INTEGER, DIMENSION(7)       :: ppm_ode_scheme_o, &
     &                                 ppm_ode_scheme_m, &
     &                                 ppm_ode_scheme_s, &
     &                                 ppm_ode_scheme_k
        DATA ppm_ode_scheme_o /1,2,2,4,2,3,1/
        DATA ppm_ode_scheme_m /1,2,2,4,2,1,0/
        DATA ppm_ode_scheme_s /1,2,2,4,2,3,999999/
        DATA ppm_ode_scheme_k /1,2,2,4,2,6,7/
        
        
        !-----------------------------------------------------------------------
        ! what scheme for which mode
        !-----------------------------------------------------------------------
        INTEGER, DIMENSION(:), POINTER :: ppm_ode_ischeme
        INTEGER, DIMENSION(:), POINTER :: ppm_ode_kscheme
        !-----------------------------------------------------------------------
        ! use an adaptive timestep
        !-----------------------------------------------------------------------
        LOGICAL, DIMENSION(:), POINTER :: ppm_ode_adaptive
        !-----------------------------------------------------------------------
        ! number of stages a mode uses
        !-----------------------------------------------------------------------
        INTEGER, DIMENSION(:), POINTER :: ppm_ode_stages
        !-----------------------------------------------------------------------
        ! state of a mode
        !-----------------------------------------------------------------------
        INTEGER, DIMENSION(:), POINTER :: ppm_ode_state
        !-----------------------------------------------------------------------
        ! number of sent stages
        !-----------------------------------------------------------------------
        INTEGER, DIMENSION(:), POINTER :: ppm_ode_sent
        !-----------------------------------------------------------------------
        ! size of the buffer
        !-----------------------------------------------------------------------
        INTEGER, DIMENSION(:), POINTER :: ppm_ode_bfrsize
        !-----------------------------------------------------------------------
        ! id lists
        !-----------------------------------------------------------------------
        INTEGER, DIMENSION(:), POINTER :: ppm_user_mid
        INTEGER, DIMENSION(:), POINTER :: ppm_internal_mid
        

        !-----------------------------------------------------------------------
        ! some stages for ppm_ode_state
        !-----------------------------------------------------------------------
        INTEGER, PARAMETER :: ppm_ode_state_finished  = 3
        INTEGER, PARAMETER :: ppm_ode_state_kickoff   = 2
        INTEGER, PARAMETER :: ppm_ode_state_running   = 1
        INTEGER, PARAMETER :: ppm_ode_state_inited    = 0
        
        
        !-----------------------------------------------------------------------
        ! number of modes
        !-----------------------------------------------------------------------
        INTEGER            :: ppm_max_mid
        INTEGER            :: ppm_max_mid_allocd

        !-----------------------------------------------------------------------
        ! topology ID to work on
        !-----------------------------------------------------------------------
        INTEGER            :: ppm_ode_topoid

      CONTAINS
        SUBROUTINE ppm_module_data_ode_activate

        END SUBROUTINE ppm_module_data_ode_activate

      END MODULE ppm_module_data_ode
