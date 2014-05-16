!------------------------------------------------------------------------------
! Subroutine :  naga_globals.f90
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
! This routines declares all global data
!------------------------------------------------------------------------------

MODULE naga_mod_globals

USE ppm_module_data
USE ppm_module_substart
USE ppm_module_substop
USE ppm_module_fft,     ONLY:ppm_fft_plan
USE ppm_module_poisson, ONLY:ppm_poisson_plan
!!USE naga_mod_stl      , ONLY:stl_data,stl_data_temporary

include 'mpif.h'

!----------------------------------------------------------------------------
! Arithmetic precision
!----------------------------------------------------------------------------
INTEGER, PARAMETER          :: MK = ppm_kind_double
INTEGER, PARAMETER          :: mpi_prec = MPI_DOUBLE_PRECISION


!----------------------------------------------------------------------------
! Data types
!----------------------------------------------------------------------------
INTEGER, PARAMETER          :: MAXCHAR = ppm_char


!----------------------------------------------------------------------------
! Geometry
!----------------------------------------------------------------------------
INTEGER,PARAMETER           :: ndim = 3         !number of physical dimensions
INTEGER,PARAMETER           :: ncom = 3         !number of vector components
!width of ghost layer. for multi-resolution (patches) this needs [lvl,patch] 
!dimensions
INTEGER , DIMENSION(3)      :: gstw


!----------------------------------------------------------------------------
! Time
!----------------------------------------------------------------------------
INTEGER                     :: itime            !time step
REAL(MK)                    :: time             !time


!----------------------------------------------------------------------------
! DEFINE PATCH: The multiresolution data format
! Split into:
! (1) a patch_setup containing all geometry and ppm bindings
! (2) the patch field data
! (3) the patch particle data
! Data field indicies are:
! - direction/component (only for the vector field)
! - x-index
! - y-index
! - z-index
! - index for the parallel sub division
! The index for parallel sub division (isub) is handled when declaring 3Darray
!----------------------------------------------------------------------------
TYPE patch_setup
  INTEGER                                :: lvl       !resolution level
  INTEGER                                :: par       !patch parent
  INTEGER, DIMENSION(ndim)               :: buf       !width of buffer
  REAL(MK),DIMENSION(ndim)               :: min       !minimum physical extent
  REAL(MK),DIMENSION(ndim)               :: max       !maximum physical extent
  INTEGER, DIMENSION(ndim)               :: gnx       !number of grid points
  REAL(MK),DIMENSION(ndim)               :: dx        !resolution
  INTEGER                                :: topoid    !topology id of level
  INTEGER                                :: meshid    !mesh id of level
  INTEGER, DIMENSION(ndim)               :: gst       !number of ghost points
  INTEGER                                :: np        !number of particles (local?)
  REAL(MK), DIMENSION(:),  POINTER       :: scost => NULL()  !sub cost
  INTEGER , DIMENSION(:,:),POINTER       :: sistr => NULL()  !sub index start
  INTEGER , DIMENSION(:,:),POINTER       :: snx   => NULL()  !sub no. grid cells
END TYPE patch_setup

!the scalar particle data
TYPE patch_part_s
  REAL(MK),DIMENSION(:), POINTER           :: val => NULL()
END TYPE patch_part_s
!the vector particle data
TYPE patch_part_v
  REAL(MK),DIMENSION(:,:), POINTER         :: val => NULL()
END TYPE patch_part_v

!the data scalar field
TYPE patch_field_s
  REAL(MK),DIMENSION(:,:,:,:), POINTER     :: fld => NULL()
END TYPE patch_field_s
!the data vector field
TYPE patch_field_v
  REAL(MK),DIMENSION(:,:,:,:,:), POINTER   :: fld => NULL()
END TYPE patch_field_v
!the data scalar complex field
TYPE patch_field_sc
  COMPLEX(MK),DIMENSION(:,:,:,:), POINTER   :: fld => NULL()
END TYPE patch_field_sc
!the data vector complex field
TYPE patch_field_vc
  COMPLEX(MK),DIMENSION(:,:,:,:,:), POINTER :: fld => NULL()
END TYPE patch_field_vc


!----------------------------------------------------------------------------
! Patch data
!----------------------------------------------------------------------------
!number of res. levels
INTEGER                                      :: nlevels
!patches pr. level (level)
INTEGER,DIMENSION(:),POINTER                 :: npatches => NULL()
!maximum number of patches pr. level
INTEGER                                      :: maxpatches
!settings for all patch (level,patch)
TYPE(patch_setup),DIMENSION(:,:),POINTER     :: ptcset => NULL()
!patch field vorticity (level,patch)
TYPE(patch_field_v),DIMENSION(:,:),POINTER   :: wf => NULL()
!patch field stream function (level,patch)
!@ Eventually we can remove psif. It is only being used for validation
TYPE(patch_field_v),DIMENSION(:,:),POINTER   :: psif => NULL()
!patch field vorticity update (level,patch)
TYPE(patch_field_v),DIMENSION(:,:),POINTER   :: dwf => NULL()
!patch field velocity (level,patch)
TYPE(patch_field_v),DIMENSION(:,:),POINTER   :: uf => NULL()
!patch field solid mask(level,patch)
TYPE(patch_field_s),DIMENSION(:,:),POINTER   :: chif => NULL()
!patch field concentration(level,patch)
TYPE(patch_field_s),DIMENSION(:,:),POINTER   :: cf => NULL()
!patch field strain rate(level,patch)
TYPE(patch_field_s),DIMENSION(:,:),POINTER   :: sf => NULL()
!patch field solid velocity (level,patch)
TYPE(patch_field_v),DIMENSION(:,:),POINTER   :: ubarf => NULL()

!---- particles and integration states
!patch particle positions(level,patch)
TYPE(patch_part_v),DIMENSION(:,:),POINTER    :: xp => NULL()
!patch particle vorticity (level,patch)
TYPE(patch_part_v),DIMENSION(:,:),POINTER    :: wp => NULL()
!patch particle velocity (level,patch)
TYPE(patch_part_v),DIMENSION(:,:),POINTER    :: up => NULL()
!patch particle vorticity rate of change (level,patch)
TYPE(patch_part_v),DIMENSION(:,:),POINTER    :: dwp => NULL()
!patch particle concentration (level,patch)
TYPE(patch_part_s),DIMENSION(:,:),POINTER    :: cp => NULL()

!patch particle positions(level,patch)
TYPE(patch_part_v),DIMENSION(:,:),POINTER    :: xp0 => NULL()
!patch particle vorticity (level,patch)
TYPE(patch_part_v),DIMENSION(:,:),POINTER    :: wp0 => NULL()
!patch particle velocity (level,patch)
TYPE(patch_part_v),DIMENSION(:,:),POINTER    :: up0 => NULL()
!patch particle vorticity rate of change (level,patch)
TYPE(patch_part_v),DIMENSION(:,:),POINTER    :: dwp0 => NULL()


!----------------------------------------------------------------------------
! PPM data
!----------------------------------------------------------------------------
!vector array of topologies (lvl)
TYPE(ppm_t_topo),DIMENSION(:), POINTER :: topos => NULL()
!!TYPE(ppm_t_topo),DIMENSION(:), POINTER :: toposxy => NULL() !@?
!!TYPE(ppm_t_topo),DIMENSION(:), POINTER :: toposz => NULL() !@?

TYPE(ppm_poisson_plan),DIMENSION(:,:), POINTER :: ppmpoisson => NULL()


!----------------------------------------------------------------------------
! MPI
!----------------------------------------------------------------------------
INTEGER                     :: comm             !MPI communicator
INTEGER                     :: rank             !MPI rank
INTEGER                     :: Nproc            !number of MPI nodes


!----------------------------------------------------------------------------
! Numerics
!----------------------------------------------------------------------------
LOGICAL                     :: restarted        !true if restarted job
LOGICAL                     :: validation       !true if restarted job


!----------------------------------------------------------------------------
! Input
!----------------------------------------------------------------------------
CHARACTER(LEN=MAXCHAR)      :: ctrlfile
CHARACTER(LEN=MAXCHAR)      :: layoutfile       !name of the patch layout file!@
CHARACTER(LEN=MAXCHAR)      :: runtag           !name of the job
REAL(MK)                    :: nu               !kinematic viscosity
INTEGER                     :: iruntag          !length of job name
INTEGER                     :: ilayoutfile      !length of layout file!@
INTEGER                     :: flowcase         !determines simulation to run
INTEGER                     :: idumpvrt         !interval to dump vorticity
INTEGER                     :: idumpdvrt        !interval to dump vorticity RHS
INTEGER                     :: idumpvel         !interval to dump velocity
INTEGER                     :: idumpchi         !interval to dump velocity
INTEGER                     :: idumpconc        !interval to dump concentration
REAL(MK)                    :: dtime            !time step size
REAL(MK)                    :: endtime          !endtime
INTEGER                     :: iremesh          !interval to remesh
INTEGER                     :: ireproject       !interval to reproject (diver.)
REAL(MK),DIMENSION(ncom)    :: uinfinity        !number of physical dimensions
INTEGER                     :: domainbc         !boundary condition for domain
INTEGER                     :: maxitime         !maximum number of time steps
INTEGER                     :: timeintscheme    !time integration scheme
INTEGER                     :: penalizationscheme!penalization FD scheme
LOGICAL                     :: clearinteriorrhs !toggle RHS in solid interior
INTEGER                     :: rhsscheme        !vorticity rhs FD scheme
INTEGER                     :: velocityscheme   !strm func. derivative scheme
INTEGER                     :: freespacekernel  !greens function for poisson eq
INTEGER                     :: lesmodel         !toggles scheme for LES model
INTEGER                     :: penalization     !toggles penalization
REAL(MK)                    :: penalizationparam!input penalization parameter
REAL(MK)                    :: penalizationstrength!parameter for penalization
LOGICAL                     :: penalizationadapt!make pen.param. depend on 1/dt
INTEGER                     :: concentration    !toggles concentration
LOGICAL                     :: forcedumpvrt     !toggle spontaneous dump of vrt
LOGICAL                     :: forcedumpvel     !toggle spontaneous dump of vel
LOGICAL                     :: forcedumpchi     !toggle spontaneous dump of mask
LOGICAL                     :: forcedumpcnc     !toggle spontaneous dump of conc
LOGICAL                     :: forcedumpdwp     !toggle spontaneous dump of dw
LOGICAL                     :: forceabort       !toggle job abortion


!----------------------------------------------------------------------------
! Flowcase input
!----------------------------------------------------------------------------
REAL(MK)                    :: torus_radius1    !ring radius
REAL(MK)                    :: torus_radius2    !ring thickness
REAL(MK)                    :: sphere_radius    !sphere radius


!----------------------------------------------------------------------------
! STL input
!----------------------------------------------------------------------------
TYPE stl_data_temporary !temporary auxillary data of the STL triangles
  CHARACTER(LEN=MAXCHAR)                  :: filename
  INTEGER                                 :: lenfilename
  REAL(MK), DIMENSION(:), ALLOCATABLE     :: tri_denom,tri_udotv
  REAL(MK), DIMENSION(:,:), ALLOCATABLE   :: tri_norm,tri_base
  REAL(MK), DIMENSION(:), ALLOCATABLE     :: tri_udotu,tri_vdotv,tri_wdotw
  REAL(MK), DIMENSION(:), ALLOCATABLE     :: tri_udotu2d,tri_vdotv2d,tri_udotv2d
  REAL(MK), DIMENSION(:), ALLOCATABLE     :: tri_denom2d
  REAL(MK), DIMENSION(:,:), ALLOCATABLE   :: tri_vecu,tri_vecv,tri_vecw
  INTEGER                                 :: tri_count
END TYPE stl_data_temporary
TYPE stl_data
  REAL(MK)                                :: bndminx,bndminy,bndminz
  REAL(MK)                                :: bndmaxx,bndmaxy,bndmaxz
END TYPE stl_data

CHARACTER(LEN=MAXCHAR)      :: stlopt_filename  !name of STL file
REAL(MK)                    :: stlopt_scale     !scaling of the STL points
REAL(MK),DIMENSION(ndim)    :: stlopt_translate !translation of the STL points
LOGICAL                     :: stlopt_check_bounding !check if domain fits STL
LOGICAL                     :: stlopt_check_intersections!check normal vectors
INTEGER                     :: stlopt_inout_direction!do additional input check
TYPE(stl_data_temporary),POINTER,DIMENSION(:) :: stlt => NULL()
TYPE(stl_data),POINTER,DIMENSION(:)           :: stl => NULL()


!----------------------------------------------------------------------------
! Step function input
!----------------------------------------------------------------------------
REAL(MK)                    :: step1_interval
REAL(MK)                    :: step1_linearfraction
REAL(MK)                    :: step1_offset


!----------------------------------------------------------------------------
! Output
!----------------------------------------------------------------------------
INTEGER                     :: debug
INTEGER                     :: validationfield  !toggles validation
INTEGER                     :: ppm_log_unit
INTEGER                     :: naga_log_unit


!----------------------------------------------------------------------------
! Stored values
!----------------------------------------------------------------------------
REAL(MK)                    :: kinenrgold
REAL(MK)                    :: lmaxstrain

END MODULE naga_mod_globals

