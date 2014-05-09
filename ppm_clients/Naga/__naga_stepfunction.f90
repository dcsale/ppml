!------------------------------------------------------------------------------
! Subroutine : naga_stl_init.f90
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
!
! Contains a step function for initialising a mollified solid mask.
! Contains an initialisation routine for the step function
! The step function is changes from 1 to 0 over a step1_interval wide interval
! via a polynomial function which is linear in the step1_linearfraction center
! of the interval and second order polynomial at the ends. When called the
! distances are typically normalised by the cell diagonal.
!-------------------------------------------------------------------------------
MODULE naga_mod_stepfunction
USE naga_mod_globals
IMPLICIT NONE
!----------------------------------------------------------------------------
! Step function variables
!----------------------------------------------------------------------------
REAL(MK) :: step1A,step1B,step1C,step1E !step function bounds
REAL(MK) :: step1vara,step1varb, & !polynomial coeff.
                             & step1varc,step1vard, &
                             & step1vare,step1varj, &
                             & step1vark,step1varl
INTERFACE naga_stepfunc1
  MODULE PROCEDURE naga_stepfunc1
END INTERFACE
INTERFACE naga_init_stepfunc1
  MODULE PROCEDURE naga_init_stepfunc1
END INTERFACE
CONTAINS
!------------------------------------------------------------------------------
! Initialise step function 1
!------------------------------------------------------------------------------
SUBROUTINE naga_init_stepfunc1(info)
USE naga_mod_globals
IMPLICIT NONE
!-------------------------------------------------------------------------
! Arguments
!-------------------------------------------------------------------------
INTEGER, INTENT(OUT) :: info
!-------------------------------------------------------------------------
! Local variables
!-------------------------------------------------------------------------
REAL(MK) :: t0
!----------------------------------------------------------------------------
! Initialise routine
!----------------------------------------------------------------------------
CALL substart('naga_add_uinfinity',t0,info)
step1A = step1_interval*0.5_mk + step1_offset
step1B = step1_interval*0.5_mk * step1_linearfraction + step1_offset
step1C = -step1_interval*0.5_mk * step1_linearfraction + step1_offset
step1E = -step1_interval*0.5_mk + step1_offset
step1vara = 1.0_mk/ (step1C * step1B - step1C * step1A - step1B ** 2 + &
          & step1A ** 2 + step1E * step1B - step1E * step1A)
step1varb = -2.0_mk / (step1C * step1B - step1C * step1A - step1B ** 2 + &
          & step1A ** 2 + step1E * step1B - step1E * step1A) * step1A
step1varc = 1.0_mk/ (step1C * step1B - step1C * step1A - step1B ** 2 + &
            & step1A ** 2 + step1E * step1B - step1E * step1A) * step1A ** 2
step1vard = 1.0_mk/ (step1A - step1C - step1E + step1B) / (-step1C + step1E)
step1vare = -2.0_mk / (step1A - step1C - step1E + step1B) / (-step1C + &
          & step1E) * step1E
step1varj = (-step1C * step1A + step1E * step1A + step1C ** 2 - step1C * &
          & step1B + step1E * step1B) / (step1A - step1C - step1E + step1B) &
          & / (-step1C + step1E)
step1vark = -2.0_mk / (step1A - step1C - step1E + step1B)
step1varl = (step1A + step1B) / (step1A - step1C - step1E + step1B)
!----------------------------------------------------------------------------
! Return
!----------------------------------------------------------------------------
 9999 CONTINUE
CALL substop('naga_add_uinfinity',t0,info)
RETURN
END SUBROUTINE naga_init_stepfunc1
!------------------------------------------------------------------------------
! Step function 1
!------------------------------------------------------------------------------
FUNCTION naga_stepfunc1(xin)
USE naga_mod_globals
IMPLICIT NONE
!-------------------------------------------------------------------------
! Arguments
!-------------------------------------------------------------------------
REAL(MK), INTENT(IN) :: xin
!-------------------------------------------------------------------------
! Local variables
!-------------------------------------------------------------------------
REAL(MK) :: x
REAL(MK) :: naga_stepfunc1
x = xin
IF (x .GE. step1A) THEN
  naga_stepfunc1 = 0.0_MK
  ELSEIF (x .GE. step1B) THEN
  naga_stepfunc1 = step1vara*x**2 + step1varb*x + step1varc
  ELSEIF (x .GE. step1C) THEN
  naga_stepfunc1 = step1vark*x + step1varl
  ELSEIF (x .GE. step1E) THEN
  naga_stepfunc1 = step1vard*x**2 + step1vare*x + step1varj
ELSE
  naga_stepfunc1 = 1.0_MK
END IF
END FUNCTION naga_stepfunc1
END MODULE naga_mod_stepfunction
