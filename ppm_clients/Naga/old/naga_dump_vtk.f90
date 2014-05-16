#ifndef __skipheaderdumpvtk
!------------------------------------------------------------------------------
! Subroutine :  naga_dump_vtk.f90
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
! These routines dumps scalar and vector fields as vtk files
! Both cell centered and vertex centered
!------------------------------------------------------------------------------
MODULE naga_mod_dump_vtk

IMPLICIT NONE

!---------------------------------------------------------------------------
! precision_tvk set the precision for vtk output
! my experience is that this can either be 4 or 8
!---------------------------------------------------------------------------
INTEGER, PARAMETER       :: precision_vtk = 8

INTERFACE naga_dump_vtk
  MODULE PROCEDURE naga_dump_vtk_s
  MODULE PROCEDURE naga_dump_vtk_v
END INTERFACE
INTERFACE naga_dump_vtk_cell
  MODULE PROCEDURE naga_dump_vtk_cell_s
  MODULE PROCEDURE naga_dump_vtk_cell_v
END INTERFACE

CONTAINS

#define __skipheaderdumpvtk
#define __dumpvtkscalar
#include "naga_dump_vtk.f90"
#undef __dumpvtkscalar
#define __dumpvtkvector
#include "naga_dump_vtk.f90"
#undef __dumpvtkvector
#define __dumpvtkcellscalar
#include "naga_dump_vtk.f90"
#undef __dumpvtkcellscalar
#define __dumpvtkcellvector
#include "naga_dump_vtk.f90"
#undef __dumpvtkcellvector

END MODULE naga_mod_dump_vtk

#else

#ifdef __dumpvtkscalar
SUBROUTINE naga_dump_vtk_s(topoid, meshid, field, tagname, iter,info)
#endif
#ifdef __dumpvtkvector
SUBROUTINE naga_dump_vtk_v(topoid, meshid, field, tagname, iter,info)
#endif
#ifdef __dumpvtkcellscalar
SUBROUTINE naga_dump_vtk_cell_s(topoid, meshid, field, tagname, iter,info)
#endif
#ifdef __dumpvtkcellvector
SUBROUTINE naga_dump_vtk_cell_v(topoid, meshid, field, tagname, iter,info)
#endif

USE naga_mod_globals
USE naga_mod_say
USE ppm_module_topo_get

IMPLICIT NONE

!-----------------------------------------------------------------------------
! arguments
!-----------------------------------------------------------------------------
INTEGER, INTENT(IN)                      :: topoid
INTEGER, INTENT(IN)                      :: meshid
#ifdef __dumpvtkscalar
REAL(MK),DIMENSION(:,:,:,:), POINTER     :: field
#endif
#ifdef __dumpvtkvector
REAL(MK),DIMENSION(:,:,:,:,:), POINTER     :: field
#endif
#ifdef __dumpvtkcellscalar
REAL(MK),DIMENSION(:,:,:,:), POINTER     :: field
#endif
#ifdef __dumpvtkcellvector
REAL(MK),DIMENSION(:,:,:,:,:), POINTER     :: field
#endif
CHARACTER(len=256), INTENT(IN)           :: tagname
INTEGER, INTENT(IN)                      :: iter
INTEGER, INTENT(INOUT)                   :: info


!-----------------------------------------------------------------------------
! Local variables
!-----------------------------------------------------------------------------
TYPE(ppm_t_topo),POINTER :: topology =>NULL()
TYPE(ppm_t_equi_mesh)    :: mesh
CHARACTER(len=256)       :: vtkfile,datatype
CHARACTER(len=256)       :: asciiline
INTEGER                  :: iasciiline
INTEGER                  :: i,j,k
INTEGER                  :: imin,jmin,kmin,imax,jmax,kmax
INTEGER                  :: iming,jming,kming,imaxg,jmaxg,kmaxg
INTEGER                  :: ioff0, nbytes0
INTEGER                  :: isubl,isub
INTEGER                  :: itagname,idatatype
REAL(MK)                 :: dx,dy,dz


!-----------------------------------------------------------------------------
! Get topology and mesh from topoid and meshid
! Compute dx,dy,dz
!-----------------------------------------------------------------------------
CALL ppm_topo_get(topoid,topology,info)
IF (info .NE. 0) THEN
  CALL naga_say(rank,'naga_dump_vtk','Failed to get topology.')
  GOTO 9999
ENDIF
mesh  = topology%mesh(meshid)

dx = (topology%max_physd(1)-topology%min_physd(1))/REAL(mesh%nm(1)-1,MK)
dy = (topology%max_physd(2)-topology%min_physd(2))/REAL(mesh%nm(2)-1,MK)
dz = (topology%max_physd(3)-topology%min_physd(3))/REAL(mesh%nm(3)-1,MK)


!-----------------------------------------------------------------------------
! Calculate data size (and offset in case of multiple fields in the file),
! Precision is set trough the precision_vtk parameter
!-----------------------------------------------------------------------------
IF (precision_vtk .EQ. 1) THEN
  datatype  = 'Float8'
  idatatype = LEN_TRIM(datatype)
  ELSEIF (precision_vtk .EQ. 2) THEN
  datatype  = 'Float16'
  idatatype = LEN_TRIM(datatype)
  ELSEIF (precision_vtk .EQ. 4) THEN
  datatype  = 'Float32'
  idatatype = LEN_TRIM(datatype)
  ELSEIF (precision_vtk .EQ. 8) THEN
  datatype  = 'Float64'
  idatatype = LEN_TRIM(datatype)
END IF


!-----------------------------------------------------------------------------
! Ensure that a proper handle has been specified
!-----------------------------------------------------------------------------
itagname = LEN_TRIM(tagname)
IF (itagname.EQ.0) THEN
  CALL naga_say(rank,'naga_dump_vtk_s','Error: Name for field must be specified')
  info = -1
  GOTO 9999
END IF


!-------------------------------------------------------------------------
! If RANK 0 output the parallel pvti file
!-------------------------------------------------------------------------
IF (rank .EQ. 0) THEN
  !-----------------------------------------------------------------------
  ! Global bounds in accordance with the ImageData format (vti)
  ! min is (0,0,0) and max is number of points (-1 since min is 0)
  !-----------------------------------------------------------------------
  iming = 0
  jming = 0
  kming = 0
  imaxg = mesh%nm(1)-1
  jmaxg = mesh%nm(2)-1
  kmaxg = mesh%nm(3)-1

  WRITE(vtkfile,'(A,A,A,A,I5.5,A)')&
  & runtag(1:iruntag),'_',tagname(1:itagname),'I',iter,'.pvti'

  OPEN(17,file=vtkfile)

  WRITE(17,'(A)') '<?xml version="1.0"?>'
  WRITE(17,'(2A)') '<VTKFile type="PImageData" version="0.1" ',&
  &'byte_order="LittleEndian">'
  WRITE(17,'(A,6(I8),A,3(E20.12),A,3(E20.12),A)') &
  &'  <PImageData WholeExtent="', iming,imaxg,jming,jmaxg,kming,kmaxg, &
  & '" Ghostlevel="0" Origin="',topology%min_physd(1),&
  & topology%min_physd(2),topology%min_physd(3), &
  & '" Spacing="',dx,dy,dz,'">'
#ifdef __dumpvtkscalar
  WRITE(17,'(4A)') '    <PPointData Scalars="', tagname(1:itagname),'">'
  WRITE(17,'(5A,I1,A)') '      <PDataArray type="',datatype(1:idatatype),&
  &'" Name="', tagname(1:itagname),'" NumberOfComponents="',1,'"/>'
  WRITE(17,'(A)') '    </PPointData>'
#endif
#ifdef __dumpvtkvector
  WRITE(17,'(3A)') '    <PPointData Vectors="', tagname(1:itagname),'">'
  WRITE(17,'(5A,I1,A)') '      <PDataArray type="',datatype(1:idatatype),&
  &'" Name="', tagname(1:itagname),'" NumberOfComponents="',ncom,'"/>'
  WRITE(17,'(A)') '    </PPointData>'
#endif
#ifdef __dumpvtkcellscalar
  WRITE(17,'(4A)') '    <PCellData Scalars="', tagname(1:itagname),'">'
  WRITE(17,'(5A,I1,A)') '      <PDataArray type="',datatype(1:idatatype),&
  &'" Name="', tagname(1:itagname),'" NumberOfComponents="',1,'"/>'
  WRITE(17,'(A)') '    </PCellData>'
#endif
#ifdef __dumpvtkcellvector
  WRITE(17,'(3A)') '    <PCellData Vectors="', tagname(1:itagname),'">'
  WRITE(17,'(5A,I1,A)') '      <PDataArray type="',datatype(1:idatatype),&
  &'" Name="', tagname(1:itagname),'" NumberOfComponents="',ncom,'"/>'
  WRITE(17,'(A)') '    </PCellData>'
#endif

  DO isub=1,topology%nsubs
    WRITE(17,'(A,A,6(I8),A,A,A,A,A,I5.5,A,I5.5,A)') '    <Piece ', &
    & 'Extent="',mesh%istart(1,isub)-1,&
    & mesh%istart(1,isub)-1+&
    &   mesh%nnodes(1,isub)-1,&
    & mesh%istart(2,isub)-1,&
    & mesh%istart(2,isub)-1+&
    &   mesh%nnodes(2,isub)-1,&
    & mesh%istart(3,isub)-1,&
    & mesh%istart(3,isub)-1+&
    &   mesh%nnodes(3,isub)-1,&
    & '" Source="', runtag(1:iruntag), &
    & '_',tagname(1:itagname),'S', isub,'I',iter, '.vti"/>'
  END DO

  WRITE(17,'(A)') '  </PImageData>'
  WRITE(17,'(A)') '</VTKFile>'
  CLOSE(17)

END IF
! End of outputting the parallel pvti file
!-----------------------------------------------------------------------------


!-----------------------------------------------------------------------------
! Loop through subdomains
!-----------------------------------------------------------------------------
DO isub=1,topology%nsublist
  isubl=topology%isublist(isub)

  !-----------------------------------------------------------------------
  ! Local bounds in accordance with the ImageData format (vti)
  !-----------------------------------------------------------------------
  ! -1 ensures the right starting index (due to the c++ 0,0,0 indexing)
  imin = mesh%istart(1,isubl)-1
  jmin = mesh%istart(2,isubl)-1
  kmin = mesh%istart(3,isubl)-1
  ! The first -1 is due to the c++ 0,0,0 indexing. The second is to get
  ! the correct end index (e.g. 3 points from 0 to 3-1 [0 1 2])
  ! The ppm decomposition entails an additional layer of points in all
  ! directions. These additional points exist on neighbouring subdomains.
  ! this corresponds well with the pvti requirements of a redundant
  ! layer of points. Therefore sub_ngrid is used
  imax = mesh%istart(1,isubl)-1+&
       & mesh%nnodes(1,isubl)-1
  jmax = mesh%istart(2,isubl)-1+&
       & mesh%nnodes(2,isubl)-1
  kmax = mesh%istart(3,isubl)-1+&
       & mesh%nnodes(3,isubl)-1

  !-----------------------------------------------------------------------
  ! Determine file name
  !-----------------------------------------------------------------------
  WRITE(vtkfile,'(A,A,A,A,I5.5,A,I5.5,A,I5.5,A)')&
  & runtag(1:iruntag),'_',tagname(1:itagname),'S',isubl,'I',iter,'.vti'

  !-----------------------------------------------------------------------
  ! Determine data chunk size (and offset which is 0 so far)
  !-----------------------------------------------------------------------
#ifdef __dumpvtkscalar
  nbytes0 = (mesh%nnodes(1,isubl))*&
          & (mesh%nnodes(2,isubl))*&
          & (mesh%nnodes(3,isubl)) &
          & * precision_vtk
#endif
#ifdef __dumpvtkvector
  nbytes0 = ncom*(mesh%nnodes(1,isubl))*&
          & (mesh%nnodes(2,isubl))*&
          & (mesh%nnodes(3,isubl))*&
          & precision_vtk
#endif
#ifdef __dumpvtkcellscalar
  nbytes0 = (mesh%nnodes(1,isubl)+1-1)*&
          & (mesh%nnodes(2,isubl)+1-1)*&
          & (mesh%nnodes(3,isubl)+1-1) &
          & * precision_vtk
#endif
#ifdef __dumpvtkcellvector
  nbytes0 = ncom*(mesh%nnodes(1,isubl)-1)*&
          & (mesh%nnodes(2,isubl)-1)*&
          & (mesh%nnodes(3,isubl)-1) * precision_vtk
#endif
  ioff0   = 0


  !-----------------------------------------------------------------------
  ! Write header
  ! - this is done to a formatted string and then written to file via
  !   unformatted output. This makes the code more portable
  !-----------------------------------------------------------------------
  ! Warning! ACCESS='STREAM' is Fortran 95 syntax but necessary!
  !-----------------------------------------------------------------------
  OPEN(17,FILE=vtkfile,FORM='UNFORMATTED',ACCESS='STREAM')

  WRITE(asciiline,'(A,A)') '<?xml version="1.0"?>',CHAR(10)
  iasciiline = LEN_TRIM(asciiline);  WRITE(17) asciiline(1:iasciiline)

  WRITE(asciiline,'(2A,A)') '<VTKFile type="ImageData" version="0.1" ',&
  &'byte_order="LittleEndian">',CHAR(10)
  iasciiline = LEN_TRIM(asciiline);  WRITE(17) asciiline(1:iasciiline)

  WRITE(asciiline,'(A,A,6(I8),A,3(E20.12),A,3(E20.12),A,A)') '  <ImageData',&
  &  ' WholeExtent="',imin,imax,jmin,jmax,kmin,kmax, '" Origin="',& 
  & topology%min_physd(1),topology%min_physd(2),&
  & topology%min_physd(3), '" Spacing="',dx,dy,dz,'">',CHAR(10)
  iasciiline = LEN_TRIM(asciiline);  WRITE(17) asciiline(1:iasciiline)

  WRITE(asciiline,'(A,A,6(I8),A,A)') '    <Piece ',&
  &  'Extent="',imin,imax,jmin,jmax,kmin,kmax,'">',CHAR(10)
  iasciiline = LEN_TRIM(asciiline);  WRITE(17) asciiline(1:iasciiline)

#ifdef __dumpvtkscalar
  WRITE(asciiline,'(3A,A)') '      <PointData Scalars="',tagname(1:itagname),'">'&
    & ,CHAR(10)
  iasciiline = LEN_TRIM(asciiline);  WRITE(17) asciiline(1:iasciiline)

  WRITE(asciiline,'(A,A,A,A,A,A,I8,A,A)') '        <DataArray type="', &
  & datatype(1:idatatype),'" Name="',tagname(1:itagname), &
  & '" NumberOfComponents="1" Format="appended" ', &
  & 'offset="',ioff0,'" />',CHAR(10)
  iasciiline = LEN_TRIM(asciiline);  WRITE(17) asciiline(1:iasciiline)

  WRITE(asciiline,'(A,A)') '      </PointData>',CHAR(10)
  iasciiline = LEN_TRIM(asciiline);  WRITE(17) asciiline(1:iasciiline)

  WRITE(asciiline,'(A,A)') '      <CellData>',CHAR(10)
  iasciiline = LEN_TRIM(asciiline);  WRITE(17) asciiline(1:iasciiline)

  WRITE(asciiline,'(A,A)') '      </CellData>',CHAR(10)
  iasciiline = LEN_TRIM(asciiline);  WRITE(17) asciiline(1:iasciiline)

#endif
#ifdef __dumpvtkvector
  WRITE(asciiline,'(3A,A)') '      <PointData Vectors="',tagname(1:itagname),'">'&
    & ,CHAR(10)
  iasciiline = LEN_TRIM(asciiline);  WRITE(17) asciiline(1:iasciiline)

  WRITE(asciiline,'(A,A,A,A,A,A,I8,A,A)') '        <DataArray type="', &
  & datatype(1:idatatype),'" Name="',tagname(1:itagname), &
  & '" NumberOfComponents="3" Format="appended" ', &
  & 'offset="',ioff0,'" />',CHAR(10)
  iasciiline = LEN_TRIM(asciiline);  WRITE(17) asciiline(1:iasciiline)

  WRITE(asciiline,'(A,A)') '      </PointData>',CHAR(10)
  iasciiline = LEN_TRIM(asciiline);  WRITE(17) asciiline(1:iasciiline)

  WRITE(asciiline,'(A,A)') '      <CellData>',CHAR(10)
  iasciiline = LEN_TRIM(asciiline);  WRITE(17) asciiline(1:iasciiline)

  WRITE(asciiline,'(A,A)') '      </CellData>',CHAR(10)
  iasciiline = LEN_TRIM(asciiline);  WRITE(17) asciiline(1:iasciiline)

#endif
#ifdef __dumpvtkcellscalar
  WRITE(asciiline,'(A,A)') '      <PointData>',CHAR(10)
  iasciiline = LEN_TRIM(asciiline);  WRITE(17) asciiline(1:iasciiline)

  WRITE(asciiline,'(A,A)') '      </PointData>',CHAR(10)
  iasciiline = LEN_TRIM(asciiline);  WRITE(17) asciiline(1:iasciiline)

  WRITE(asciiline,'(3A,A)') '      <CellData Scalars="',tagname(1:itagname),'">'&
    & ,CHAR(10)
  iasciiline = LEN_TRIM(asciiline);  WRITE(17) asciiline(1:iasciiline)

  WRITE(asciiline,'(A,A,A,A,A,A,I8,A,A)') '        <DataArray type="', &
  & datatype(1:idatatype),'" Name="',tagname(1:itagname), &
  & '" NumberOfComponents="1" Format="appended" ', &
  & 'offset="',ioff0,'" />',CHAR(10)
  iasciiline = LEN_TRIM(asciiline);  WRITE(17) asciiline(1:iasciiline)

  WRITE(asciiline,'(A,A)') '      </CellData>',CHAR(10)
  iasciiline = LEN_TRIM(asciiline);  WRITE(17) asciiline(1:iasciiline)

#endif
#ifdef __dumpvtkcellvector
  WRITE(asciiline,'(A,A)') '      <PointData>',CHAR(10)
  iasciiline = LEN_TRIM(asciiline);  WRITE(17) asciiline(1:iasciiline)

  WRITE(asciiline,'(A,A)') '      </PointData>',CHAR(10)
  iasciiline = LEN_TRIM(asciiline);  WRITE(17) asciiline(1:iasciiline)

  WRITE(asciiline,'(3A,A)') '      <CellData Vectors="',tagname(1:itagname),'">'&
    & ,CHAR(10)
  iasciiline = LEN_TRIM(asciiline);  WRITE(17) asciiline(1:iasciiline)

  WRITE(asciiline,'(A,A,A,A,A,A,I8,A,A)') '        <DataArray type="', &
  & datatype(1:idatatype),'" Name="',tagname(1:itagname), &
  & '" NumberOfComponents="3" Format="appended" ', &
  & 'offset="',ioff0,'" />',CHAR(10)
  iasciiline = LEN_TRIM(asciiline);  WRITE(17) asciiline(1:iasciiline)

  WRITE(asciiline,'(A,A)') '      </CellData>',CHAR(10)
  iasciiline = LEN_TRIM(asciiline);  WRITE(17) asciiline(1:iasciiline)

#endif
  WRITE(asciiline,'(A,A)') '    </Piece>',CHAR(10)
  iasciiline = LEN_TRIM(asciiline);  WRITE(17) asciiline(1:iasciiline)

  WRITE(asciiline,'(A,A)') '  </ImageData>',CHAR(10)
  iasciiline = LEN_TRIM(asciiline);  WRITE(17) asciiline(1:iasciiline)

  WRITE(asciiline,'(A,A)') '  <AppendedData encoding="raw">',CHAR(10)
  iasciiline = LEN_TRIM(asciiline);  WRITE(17) asciiline(1:iasciiline)

  WRITE(17) '_'

  !-----------------------------------------------------------------------
  ! Finished XML file - Write data, and close
  ! output the redundant point layer which matches pvti point requirements
  !-----------------------------------------------------------------------
  WRITE(17) nbytes0
#ifdef __dumpvtkscalar
  DO k=1,mesh%nnodes(3,isubl)
    DO j=1,mesh%nnodes(2,isubl)
      DO i=1,mesh%nnodes(1,isubl)
        WRITE(17) REAL(field(i,j,k,isub),precision_vtk)
#endif
#ifdef __dumpvtkvector
  DO k=1,mesh%nnodes(3,isubl)
    DO j=1,mesh%nnodes(2,isubl)
      DO i=1,mesh%nnodes(1,isubl)
        WRITE(17) REAL(field(1,i,j,k,isub), &
        & precision_vtk), &
        & REAL(field(2,i,j,k,isub), &
        & precision_vtk), &
        & REAL(field(3,i,j,k,isub), &
        & precision_vtk)
#endif
#ifdef __dumpvtkcellscalar
  DO k=1,mesh%nnodes(3,isubl)-1
    DO j=1,mesh%nnodes(2,isubl)-1
      DO i=1,mesh%nnodes(1,isubl)-1
        WRITE(17) REAL(field(i,j,k,isub),precision_vtk)
#endif
#ifdef __dumpvtkcellvector
  DO k=1,mesh%nnodes(3,isubl)-1
    DO j=1,mesh%nnodes(2,isubl)-1
      DO i=1,mesh%nnodes(1,isubl)-1
        WRITE(17) REAL(field(1,i,j,k,isub), &
        & precision_vtk), &
        & REAL(field(2,i,j,k,isub), &
        & precision_vtk), &
        & REAL(field(3,i,j,k,isub), &
        & precision_vtk)
#endif
      END DO
    END DO
  END DO

  WRITE(asciiline,'(A,A,A)') CHAR(10),'  </AppendedData>',CHAR(10)
  iasciiline = LEN_TRIM(asciiline);  WRITE(17) asciiline(1:iasciiline)

  WRITE(asciiline,'(A,A)') '</VTKFile>',CHAR(10)
  CLOSE(17)

END DO !ilevel


9999 CONTINUE
RETURN


#ifdef __dumpvtkscalar
END SUBROUTINE naga_dump_vtk_s
#endif
#ifdef __dumpvtkvector
END SUBROUTINE naga_dump_vtk_v
#endif
#ifdef __dumpvtkcellscalar
END SUBROUTINE naga_dump_vtk_cell_s
#endif
#ifdef __dumpvtkcellvector
END SUBROUTINE naga_dump_vtk_cell_v
#endif

#endif
