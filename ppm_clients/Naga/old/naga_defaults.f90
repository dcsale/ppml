!------------------------------------------------------------------------------
! Subroutine :  naga_defaults.f90
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
! This routines sets default values
!------------------------------------------------------------------------------

MODULE naga_mod_defaults

IMPLICIT NONE

INTERFACE naga_defaults
  MODULE PROCEDURE naga_defaults
END INTERFACE

CONTAINS

SUBROUTINE naga_defaults(info)

USE naga_mod_globals

IMPLICIT NONE
!-------------------------------------------------------------------------
! Arguments
!-------------------------------------------------------------------------
INTEGER, INTENT(OUT) :: info


!-------------------------------------------------------------------------
! Local variables
!-------------------------------------------------------------------------
REAL(MK)             :: t0


!----------------------------------------------------------------------------
! Initialise routine
!----------------------------------------------------------------------------
CALL substart('naga_defaults',t0,info)


!----------------------------------------------------------------------------
! Input
!----------------------------------------------------------------------------
ctrlfile         = 'CTRL'
flowcase         = 0


!----------------------------------------------------------------------------
! Output
!----------------------------------------------------------------------------
idumpvrt        = 0
idumpdvrt       = 0
idumpvel        = 0
idumpchi        = 0
idumpconc       = 0
runtag          = 'nagarun'
debug           = 0
validationfield = 0
ppm_log_unit    = 99
naga_log_unit   = 88

iruntag = LEN_TRIM(runtag)


!----------------------------------------------------------------------------
! Physics
!----------------------------------------------------------------------------
uinfinity            = (/0.0,0.0,0.0/)
nu                   = 0.01
penalization         = 0
penalizationparam    = 1.0
penalizationstrength = 1.0
penalizationadapt    = .TRUE.
concentration        = 0
lesmodel             = 0


!----------------------------------------------------------------------------
! Time
!----------------------------------------------------------------------------
itime              = 0
dtime              = 0.0_MK
endtime            = 0.0_MK
maxitime           = 0


!----------------------------------------------------------------------------
! Numerics
!----------------------------------------------------------------------------
timeintscheme      = 1
penalizationscheme = 1
clearinteriorrhs   = .FALSE.
rhsscheme          = 2
velocityscheme     = 0
velocityscheme     = 1
gstw               = 2
restarted          = .FALSE.
validation         = .FALSE.
iremesh            = 1
ireproject         = 0


!----------------------------------------------------------------------------
! Step function input
!----------------------------------------------------------------------------
step1_interval       = 5.0_MK
step1_linearfraction = 0.0_MK
step1_offset         = 0.0_MK


!----------------------------------------------------------------------------
! Input
!----------------------------------------------------------------------------


!----------------------------------------------------------------------------
! Geometry
!----------------------------------------------------------------------------
nlevels     = 0
maxpatches  = 0
domainbc    = 1

!----------------------------------------------------------------------------
! Specific for flowcases
!----------------------------------------------------------------------------
! torus
torus_radius1 = 0.25_MK
torus_radius2 = 0.05_MK
! sphere
sphere_radius = 0.25_MK
! STL
stlopt_filename        = ''
stlopt_scale           = 1.0_MK
stlopt_translate       = (/0.0_MK,0.0_MK,0.0_MK/)
stlopt_check_bounding  = .FALSE.
stlopt_inout_direction = 0

!----------------------------------------------------------------------------
! Return
!----------------------------------------------------------------------------
 9999 CONTINUE
CALL substop('naga_defaults',t0,info)
RETURN


END SUBROUTINE naga_defaults

END MODULE naga_mod_defaults

