      !-------------------------------------------------------------------------
      ! Module : client_module_realloc
      !-------------------------------------------------------------------------
      !
      ! Purpose : This module contains subroutines to dynamically
      ! reallocate arrays.
      !
      ! Remarks :
      !
      ! References :
      !-------------------------------------------------------------------------
      ! Parallel Particle Mesh Library (PPM)
      ! ETH Zurich
      ! CH-8092 Zurich, Switzerland
      !-------------------------------------------------------------------------
      MODULE client_module_realloc
      IMPLICIT NONE
      !-------------------------------------------------------------------------
      ! Reallocate a 1d or 2d array to a new length,preserving its
      ! contents.
      !-------------------------------------------------------------------------
      INTERFACE reallocate
         MODULE PROCEDURE reallocate_rv,reallocate_rm,reallocate_iv, &
              reallocate_im,reallocate_hv
      END INTERFACE
      CONTAINS
      !-------------------------------------------------------------------------
      ! The source
      !-------------------------------------------------------------------------
      ! for real vectors
      function reallocate_rv(p,n)
        USE client_global
        IMPLICIT NONE
        REAL(MK),DIMENSION(:),POINTER :: p,reallocate_rv
        INTEGER,INTENT(IN) :: n
        INTEGER :: nold,ierr
        ALLOCATE(reallocate_rv(n),STAT=ierr)
        IF (ierr .NE. 0) THEN
           CALL pwrite(rank,'client_module_realloc','Not enough memory',ierr)
           RETURN
        ENDIF
        IF (.NOT. ASSOCIATED(p)) RETURN
        nold = SIZE(p)
        reallocate_rv(1:MIN(nold,n))=p(1:MIN(nold,n))
        DEALLOCATE(p)
      END function reallocate_rv
      ! for integer vectors
      function reallocate_iv(p,n)
        USE client_global
        IMPLICIT NONE
        INTEGER,DIMENSION(:),POINTER :: p,reallocate_iv
        INTEGER,INTENT(IN) :: n
        INTEGER :: nold,ierr
        ALLOCATE(reallocate_iv(n),STAT=ierr)
        IF (ierr .NE. 0) THEN
           CALL pwrite(rank,'client_module_realloc','Not enough memory',ierr)
           RETURN
        ENDIF
        IF (.NOT. ASSOCIATED(p)) RETURN
        nold = SIZE(p)
        reallocate_iv(1:MIN(nold,n))=p(1:MIN(nold,n))
        DEALLOCATE(p)
      END function reallocate_iv
      ! for character vectors
      function reallocate_hv(p,n)
        USE client_global
        IMPLICIT NONE
        CHARACTER(1),DIMENSION(:),POINTER :: p,reallocate_hv
        INTEGER,INTENT(IN) :: n
        INTEGER :: nold,ierr
        ALLOCATE(reallocate_hv(n),STAT=ierr)
        IF (ierr .NE. 0) THEN
           CALL pwrite(rank,'client_module_realloc','Not enough memory',ierr)
           RETURN
        ENDIF
        IF (.NOT. ASSOCIATED(p)) RETURN
        nold = SIZE(p)
        reallocate_hv(1:MIN(nold,n))=p(1:MIN(nold,n))
        DEALLOCATE(p)
      END function reallocate_hv
      ! for real matrices
      function reallocate_rm(p,n1,n2)
        USE client_global
        IMPLICIT NONE
        REAL(MK),DIMENSION(:,:),POINTER :: p,reallocate_rm
        INTEGER,INTENT(IN) :: n1,n2
        INTEGER :: nold,mold,ierr
        ALLOCATE(reallocate_rm(n1,n2),STAT=ierr)
        IF (ierr .NE. 0) THEN
           CALL pwrite(rank,'client_module_realloc','Not enough memory',ierr)
           RETURN
        ENDIF
        IF (.NOT. ASSOCIATED(p)) RETURN
        nold = SIZE(p,1)
        mold = SIZE(p,2)
        reallocate_rm(1:MIN(nold,n1),1:MIN(mold,n2)) = &
     & p(1:MIN(nold,n1),1:MIN(mold,n2))
        DEALLOCATE(p)
      END function reallocate_rm
      ! for integer matrices
      function reallocate_im(p,n1,n2)
        USE client_global
        IMPLICIT NONE
        INTEGER,DIMENSION(:,:),POINTER :: p,reallocate_im
        INTEGER,INTENT(IN) :: n1,n2
        INTEGER :: nold,mold,ierr
        ALLOCATE(reallocate_im(n1,n2),STAT=ierr)
        IF (ierr .NE. 0) THEN
           CALL pwrite(rank,'client_module_realloc','Not enough memory',ierr)
           RETURN
        ENDIF
        IF (.NOT. ASSOCIATED(p)) RETURN
        nold = SIZE(p,1)
        mold = SIZE(p,2)
        reallocate_im(1:MIN(nold,n1),1:MIN(mold,n2)) = &
     & p(1:MIN(nold,n1),1:MIN(mold,n2))
        DEALLOCATE(p)
      END function reallocate_im
      END MODULE client_module_realloc
