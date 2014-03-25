      !-------------------------------------------------------------------------
      !  Include file for reading the data from fdta (level function).
      !
      !  INPUT:   INTEGER :: i,j,k -- Point to solve for
      !           INTEGER :: order -- Desired order of FD scheme
      !  OUTPUT:  REAL(MK)(-order:order,1:ppm_dim) :: phi -- level
      !                      function at shifted locations in each dir.
      !
      !-------------------------------------------------------------------------
      !  $Log: ppm_gmm_getfdta.inc,v $
      !  Revision 1.1.1.1  2007/07/13 10:18:55  ivos
      !  CBL version of the PPM library
      !
      !  Revision 1.1  2005/05/10 04:41:15  ivos
      !  Newly created during modularization of ppm_gmm_march. Marching
      !  and orthogonal extendion are now in separate routines for faster
      !  compilation.
      !
      !-------------------------------------------------------------------------
      !  Parallel Particle Mesh Library (PPM)
      !  Institute of Computational Science
      !  ETH Zentrum, Hirschengraben 84
      !  CH-8092 Zurich, Switzerland
      !-------------------------------------------------------------------------

      phi = 0.0_MK
#if   __DIM == __3D
      IF (order .EQ. ppm_param_order_3) THEN
          phi(-3,1) = fdta(i-3,j,k,jsub)
          phi(-2,1) = fdta(i-2,j,k,jsub)
          phi(-1,1) = fdta(i-1,j,k,jsub)
          phi( 0,1) = fdta(i  ,j,k,jsub)
          phi( 1,1) = fdta(i+1,j,k,jsub)
          phi( 2,1) = fdta(i+2,j,k,jsub)
          phi( 3,1) = fdta(i+3,j,k,jsub)

          phi(-3,2) = fdta(i,j-3,k,jsub)
          phi(-2,2) = fdta(i,j-2,k,jsub)
          phi(-1,2) = fdta(i,j-1,k,jsub)
          phi( 0,2) = fdta(i,j  ,k,jsub)
          phi( 1,2) = fdta(i,j+1,k,jsub)
          phi( 2,2) = fdta(i,j+2,k,jsub)
          phi( 3,2) = fdta(i,j+3,k,jsub)

          phi(-3,3) = fdta(i,j,k-3,jsub)
          phi(-2,3) = fdta(i,j,k-2,jsub)
          phi(-1,3) = fdta(i,j,k-1,jsub)
          phi( 0,3) = fdta(i,j,k  ,jsub)
          phi( 1,3) = fdta(i,j,k+1,jsub)
          phi( 2,3) = fdta(i,j,k+2,jsub)
          phi( 3,3) = fdta(i,j,k+3,jsub)
      ELSEIF (order .EQ. ppm_param_order_2) THEN
          phi(-2,1) = fdta(i-2,j,k,jsub)
          phi(-1,1) = fdta(i-1,j,k,jsub)
          phi( 0,1) = fdta(i  ,j,k,jsub)
          phi( 1,1) = fdta(i+1,j,k,jsub)
          phi( 2,1) = fdta(i+2,j,k,jsub)

          phi(-2,2) = fdta(i,j-2,k,jsub)
          phi(-1,2) = fdta(i,j-1,k,jsub)
          phi( 0,2) = fdta(i,j  ,k,jsub)
          phi( 1,2) = fdta(i,j+1,k,jsub)
          phi( 2,2) = fdta(i,j+2,k,jsub)

          phi(-2,3) = fdta(i,j,k-2,jsub)
          phi(-1,3) = fdta(i,j,k-1,jsub)
          phi( 0,3) = fdta(i,j,k  ,jsub)
          phi( 1,3) = fdta(i,j,k+1,jsub)
          phi( 2,3) = fdta(i,j,k+2,jsub)
      ELSE
          phi(-1,1) = fdta(i-1,j,k,jsub)
          phi( 0,1) = fdta(i  ,j,k,jsub)
          phi( 1,1) = fdta(i+1,j,k,jsub)

          phi(-1,2) = fdta(i,j-1,k,jsub)
          phi( 0,2) = fdta(i,j  ,k,jsub)
          phi( 1,2) = fdta(i,j+1,k,jsub)

          phi(-1,3) = fdta(i,j,k-1,jsub)
          phi( 0,3) = fdta(i,j,k  ,jsub)
          phi( 1,3) = fdta(i,j,k+1,jsub)
      ENDIF
#elif __DIM == __2D
      IF (order .EQ. ppm_param_order_3) THEN
          phi(-3,1) = fdta(i-3,j,jsub)
          phi(-2,1) = fdta(i-2,j,jsub)
          phi(-1,1) = fdta(i-1,j,jsub)
          phi( 0,1) = fdta(i  ,j,jsub)
          phi( 1,1) = fdta(i+1,j,jsub)
          phi( 2,1) = fdta(i+2,j,jsub)
          phi( 3,1) = fdta(i+3,j,jsub)

          phi(-3,2) = fdta(i,j-3,jsub)
          phi(-2,2) = fdta(i,j-2,jsub)
          phi(-1,2) = fdta(i,j-1,jsub)
          phi( 0,2) = fdta(i,j  ,jsub)
          phi( 1,2) = fdta(i,j+1,jsub)
          phi( 2,2) = fdta(i,j+2,jsub)
          phi( 3,2) = fdta(i,j+3,jsub)
      ELSEIF (order .EQ. ppm_param_order_2) THEN
          phi(-2,1) = fdta(i-2,j,jsub)
          phi(-1,1) = fdta(i-1,j,jsub)
          phi( 0,1) = fdta(i  ,j,jsub)
          phi( 1,1) = fdta(i+1,j,jsub)
          phi( 2,1) = fdta(i+2,j,jsub)

          phi(-2,2) = fdta(i,j-2,jsub)
          phi(-1,2) = fdta(i,j-1,jsub)
          phi( 0,2) = fdta(i,j  ,jsub)
          phi( 1,2) = fdta(i,j+1,jsub)
          phi( 2,2) = fdta(i,j+2,jsub)
      ELSE
          phi(-1,1) = fdta(i-1,j,jsub)
          phi( 0,1) = fdta(i  ,j,jsub)
          phi( 1,1) = fdta(i+1,j,jsub)

          phi(-1,2) = fdta(i,j-1,jsub)
          phi( 0,2) = fdta(i,j  ,jsub)
          phi( 1,2) = fdta(i,j+1,jsub)
      ENDIF
#endif
