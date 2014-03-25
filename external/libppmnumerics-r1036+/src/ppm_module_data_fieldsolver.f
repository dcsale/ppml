      !-------------------------------------------------------------------------
      ! Module         :            ppm_module_data_fieldsolver
      !-------------------------------------------------------------------------
      !
      ! Purpose       :  data module of fieldsolver mostly containing the FFT
      !                        plans
      !               
      !
      ! Remarks       :
      !
      ! References    : 
      !
      ! Revisions     :
      !-------------------------------------------------------------------------
      !  $Log: ppm_module_data_fieldsolver.f,v $
      !  Revision 1.1.1.1  2007/07/13 10:18:57  ivos
      !  CBL version of the PPM library
      !
      !  Revision 1.3  2005/02/17 17:47:00  hiebers
      !  Reimplementation
      !
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

MODULE ppm_module_data_fieldsolver

   USE ppm_module_data,ONLY:ppm_kind_single,ppm_kind_double
   PRIVATE :: ppm_kind_single,ppm_kind_double

      ! FFTW Plans
      INTEGER*8 Plan_fd_s,     Plan_fd_d 
      INTEGER*8 Plan_slab_fd_s,Plan_slab_fd_d 
      INTEGER*8 Plan_fd_c_y,   Plan_fd_cc_y
      INTEGER*8 Plan_fd_c_z,   Plan_fd_cc_z 
      INTEGER*8 Plan_bd_s,     Plan_bd_d
      INTEGER*8 Plan_slab_bd_s,Plan_slab_bd_d 
      INTEGER*8 Plan_bd_c_y,   Plan_bd_cc_y
      INTEGER*8 Plan_bd_c_z,   Plan_bd_cc_z 


      ! MATHKEISAN variables for MathKeisan FFTs
      ! working storage
      REAL(ppm_kind_single), DIMENSION(:),POINTER :: table_fd_s,   table_bd_s
      REAL(ppm_kind_double), DIMENSION(:),POINTER :: table_fd_d,   table_bd_d
      REAL(ppm_kind_single), DIMENSION(:),POINTER :: table_fd_c_y, table_bd_c_y
      REAL(ppm_kind_double), DIMENSION(:),POINTER :: table_fd_cc_y,table_bd_cc_y
      REAL(ppm_kind_single), DIMENSION(:),POINTER :: table_fd_c_z, table_bd_c_z
      REAL(ppm_kind_double), DIMENSION(:),POINTER :: table_fd_cc_z,table_bd_cc_z

      ! the size of the working storage
      INTEGER, DIMENSION(1)              :: lda_table, lda_table_y, lda_table_z


END MODULE ppm_module_data_fieldsolver
