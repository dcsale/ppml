      !-------------------------------------------------------------------------
      !  Include file for adding a point to the list gmm_ipos and 
      !  growing the list if needed.
      !
      !  INPUT:   INTEGER :: i,j,(k) -- point to add
      !           INTEGER :: jsub -- local sub index of point to add
      !  OUTPUT:  
      !
      !-------------------------------------------------------------------------
      !  $Log: ppm_gmm_add_to_list.inc,v $
      !  Revision 1.1.1.1  2007/07/13 10:18:55  ivos
      !  CBL version of the PPM library
      !
      !  Revision 1.3  2005/05/10 04:48:44  ivos
      !  Split marching and extension routines for faster compilation,
      !  Sharked extension routines, moved all initialization to gmm_init, and
      !  code cosmetics.
      !
      !  Revision 1.2  2005/04/27 01:06:09  ivos
      !  Convergence tests completed, cleaned up code, optmized code (Shark),
      !  and changed structure to allow faster compilation.
      !
      !  Revision 1.1  2005/03/10 01:38:42  ivos
      !  Initial check-in.
      !
      !-------------------------------------------------------------------------
      !  Parallel Particle Mesh Library (PPM)
      !  Institute of Computational Science
      !  ETH Zentrum, Hirschengraben 84
      !  CH-8092 Zurich, Switzerland
      !-------------------------------------------------------------------------
#if   __DIM == __3D
      npos0 = npos0 + 1
      IF (gmm_lsiz .LT. npos0) THEN
          gmm_lsiz = gmm_lsiz + incr
          iopt = ppm_param_alloc_grow_preserve
          ldu(1) = 4
          ldu(2) = gmm_lsiz
          CALL ppm_alloc(gmm_ipos,ldu,iopt,info)
          IF (info .NE. 0) THEN
              info = ppm_error_fatal
              CALL ppm_error(ppm_err_alloc,   &
     &            'ppm_gmm_march','sparse positions GMM_IPOS',      &
     &            __LINE__,info)
              GOTO 9999
          ENDIF
      ENDIF
      ! append to gmm_ipos
      gmm_ipos(1,npos0) = i
      gmm_ipos(2,npos0) = j
      gmm_ipos(3,npos0) = k
      gmm_ipos(4,npos0) = jsub
#elif __DIM == __2D
      npos0 = npos0 + 1
      IF (gmm_lsiz .LT. npos0) THEN
          gmm_lsiz = gmm_lsiz + incr
          iopt = ppm_param_alloc_grow_preserve
          ldu(1) = 3
          ldu(2) = gmm_lsiz
          CALL ppm_alloc(gmm_ipos,ldu,iopt,info)
          IF (info .NE. 0) THEN
              info = ppm_error_fatal
              CALL ppm_error(ppm_err_alloc,   &
     &            'ppm_gmm_march','sparse positions GMM_IPOS',      &
     &            __LINE__,info)
              GOTO 9999
          ENDIF
      ENDIF
      ! append to gmm_ipos
      gmm_ipos(1,npos0) = i
      gmm_ipos(2,npos0) = j
      gmm_ipos(3,npos0) = jsub
#endif
