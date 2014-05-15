!------------------------------------------------------------------------------
! Subroutine :  naga_stl.f90
!------------------------------------------------------------------------------
!      _  __
!     / |/ /__ ____ ____ _
!    /    / _ `/ _ `/ _ `/
!   /_/|_/\_,_/\_, /\_,_/
!             /___/
!   A PPM Vortex-In-Cell client - 2011 (c)
!   Naga is distributed under the GNU Lesser General Public License version 3
!   Johannes Tophøj Rasmussen, johannes.tophoj@gmail.com
!   Technical University of Denmark
!   Department of Mechanical Engineering
!
!
! This module contains routines for initialising solids from STL input
!------------------------------------------------------------------------------

MODULE naga_mod_stl

USE naga_mod_globals

IMPLICIT NONE


INTERFACE naga_stl_read
  MODULE PROCEDURE naga_stl_read
END INTERFACE

INTERFACE naga_stl_init
  MODULE PROCEDURE naga_stl_init
END INTERFACE

CONTAINS

#include "naga_stl_read.f90"
#include "naga_stl_init.f90"

END MODULE naga_mod_stl


