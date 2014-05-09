!------------------------------------------------------------------------------
! Naga.f90
!------------------------------------------------------------------------------
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
! Copyright 2011
! Naga is free software: you can redistribute it and/or modify
! it under the terms of the GNU Lesser General Public License
! as published by the Free Software Foundation, either
! version 3 of the License, or (at your option) any later
! version. See <http:
!------------------------------------------------------------------------------
! Naga is vortex-in-cell (VIC) implementation using the PPM Library. The name
! originates from a D&D Creature that may live in water.
! This software is copyrighted but released under the
! GNU Lesser General Public License. If you use this software please give
! proper credits by citing the following articles:
! * I.F. Sbalzarini, J.H. Walther, M. Bergdorf, S.E. Hieber, E.M. Kotsalis,
! P. Koumoutsakos, PPM – a highly efficient parallel particle-mesh library for
! the simulation of continuum systems, J. Comput. Phys. 215 (2006) 566–588.
! Hejlesen: kernel article
!
! If you find it relevant you may cite the following paper on the VIC:
! * J. T. Rasmussen, G.-H. Cottet, and J. H. Walther. A multiresolution remeshed
! particle vortex method using patches. J. Comput. Phys., 230(17):6742–6755, 2011.
!
! For more information on particle methods, PPM and VIC consult the dissertation:
! Particle Methods in Bluff Body Aerodynamics, J.T. Rasmussen, DTU, 2012
! <http:
!
! Naga is designed towards multi-resolution but this will not be fully
! implemented until staggered meshes become available in PPM. Hence the
! ipatch/ilevel loops throughout the code
!
! - cpp statements are used throughout to expand implementations for several
! data types.
!
! The STL-initialisation routine may incorrectly determine inside/outside
! values of the solid field if the STL-coordinates align with the grid points.
! Shift the STL file by a random value (by a very low value) using the CTRL
! file parameters. Always check the solid field before running a simulation!
!------------------------------------------------------------------------------
PROGRAM naga
USE naga_mod_globals
USE naga_mod_say
USE naga_mod_init
USE naga_mod_setup
USE naga_mod_integrate_rk1
USE naga_mod_integrate_rk2
USE naga_mod_validation
USE ppm_module_poisson
IMPLICIT NONE
!-------------------------------------------------------------------------
! Local variables
!-------------------------------------------------------------------------
INTEGER :: info !return variable
INTEGER :: ilevel,ipatch
info = 0
!-------------------------------------------------------------------------
! Initialise MPI, PPM, read input files, set defaults
!-------------------------------------------------------------------------
CALL naga_init(info)
IF (info .NE. 0) THEN
  CALL naga_say(rank,'Naga','Failed to initialize.')
  GOTO 9999
ENDIF
!-------------------------------------------------------------------------
! Say hello to the world and introduce ourselves
!-------------------------------------------------------------------------
IF (rank .EQ. 0) THEN
  write(*,*) "     _  __"
  write(*,*) "    / |/ /__ ____ ____ _"
  write(*,*) "   /    / _ `/ _ `/ _ `/"
  write(*,*) "  /_/|_/\_,_/\_, /\_,_/"
  write(*,*) "            /___/         says hello,"
  write(*,*) "                          today I feel like driving Monster Trucks"
  write(*,*) ""
ENDIF
!-------------------------------------------------------------------------
! Setup the simulation: Allocate and initialise arrays. All job specific
! initialisation
!-------------------------------------------------------------------------
CALL naga_setup(info)
IF (info .NE. 0) THEN
  CALL naga_say(rank,'Naga','Failed to setup the simulation.')
  GOTO 9999
ENDIF
!-------------------------------------------------------------------------
! Do the time integration: Particle advection, particle strength update
!-------------------------------------------------------------------------
IF (timeintscheme .EQ. 1) THEN
  CALL naga_integrate_rk1(info)
  IF (info .NE. 0) THEN
    CALL naga_say(rank,'Naga','Error performing the rk1 time integration.')
    GOTO 9999
  ENDIF
ELSE IF (timeintscheme .EQ. 2) THEN
  CALL naga_integrate_rk2(info)
  IF (info .NE. 0) THEN
    CALL naga_say(rank,'Naga','Error performing the rk2 time integration.')
    GOTO 9999
  ENDIF
ENDIF
!-------------------------------------------------------------------------
! When doing validation call the comparison routine after the simulation
!-------------------------------------------------------------------------
IF (validation) THEN
  CALL naga_validation(info)
  IF (info .NE. 0) THEN
    CALL naga_say(rank,'Naga','Error calling the validation routine.')
    GOTO 9999
  ENDIF
ENDIF
!-------------------------------------------------------------------------
! Finalize: Deallocate arrays etc. More to come here
!-------------------------------------------------------------------------
! not sure this is needed anymore in v1.2.2
!DO ilevel=1,nlevels
! ipatch=1
! CALL ppm_poisson_finalize(ppmpoisson(ilevel,ipatch),info)
!END DO
!-------------------------------------------------------------------------
! Return·
!-------------------------------------------------------------------------
 9999 CONTINUE
CALL MPI_finalize(info)
IF (info .NE. 0) THEN
  CALL naga_say(rank, 'naga_setup','Failed to finalize MPI.')
ENDIF
END PROGRAM naga
