      !-------------------------------------------------------------------------
      !  Include file for reading the function values from dta.
      !
      !  INPUT:   INTEGER :: i,j,k -- Point to solve for
      !           INTEGER :: order -- Desired order of FD scheme
      !  OUTPUT:  REAL(MK)(-order:order,1:ppm_dim) :: psi -- function
      !                      values at shifted locations in all directions.
      !
      !-------------------------------------------------------------------------
      !  $Log: ppm_gmm_getdta.inc,v $
      !  Revision 1.1.1.1  2007/07/13 10:18:55  ivos
      !  CBL version of the PPM library
      !
      !  Revision 1.1  2005/05/10 04:41:14  ivos
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

      psi = 0.0_MK
#if   __DIM == __3D
      IF (order .EQ. ppm_param_order_3) THEN
          psi(-3,1) = dta(i-3,j,k,jsub)
          psi(-2,1) = dta(i-2,j,k,jsub)
          psi(-1,1) = dta(i-1,j,k,jsub)
          psi( 0,1) = dta(i  ,j,k,jsub)
          psi( 1,1) = dta(i+1,j,k,jsub)
          psi( 2,1) = dta(i+2,j,k,jsub)
          psi( 3,1) = dta(i+3,j,k,jsub)

          psi(-3,2) = dta(i,j-3,k,jsub)
          psi(-2,2) = dta(i,j-2,k,jsub)
          psi(-1,2) = dta(i,j-1,k,jsub)
          psi( 0,2) = dta(i,j  ,k,jsub)
          psi( 1,2) = dta(i,j+1,k,jsub)
          psi( 2,2) = dta(i,j+2,k,jsub)
          psi( 3,2) = dta(i,j+3,k,jsub)

          psi(-3,3) = dta(i,j,k-3,jsub)
          psi(-2,3) = dta(i,j,k-2,jsub)
          psi(-1,3) = dta(i,j,k-1,jsub)
          psi( 0,3) = dta(i,j,k  ,jsub)
          psi( 1,3) = dta(i,j,k+1,jsub)
          psi( 2,3) = dta(i,j,k+2,jsub)
          psi( 3,3) = dta(i,j,k+3,jsub)
      ELSEIF (order .EQ. ppm_param_order_2) THEN
          psi(-2,1) = dta(i-2,j,k,jsub)
          psi(-1,1) = dta(i-1,j,k,jsub)
          psi( 0,1) = dta(i  ,j,k,jsub)
          psi( 1,1) = dta(i+1,j,k,jsub)
          psi( 2,1) = dta(i+2,j,k,jsub)

          psi(-2,2) = dta(i,j-2,k,jsub)
          psi(-1,2) = dta(i,j-1,k,jsub)
          psi( 0,2) = dta(i,j  ,k,jsub)
          psi( 1,2) = dta(i,j+1,k,jsub)
          psi( 2,2) = dta(i,j+2,k,jsub)

          psi(-2,3) = dta(i,j,k-2,jsub)
          psi(-1,3) = dta(i,j,k-1,jsub)
          psi( 0,3) = dta(i,j,k  ,jsub)
          psi( 1,3) = dta(i,j,k+1,jsub)
          psi( 2,3) = dta(i,j,k+2,jsub)
      ELSE
          psi(-1,1) = dta(i-1,j,k,jsub)
          psi( 0,1) = dta(i  ,j,k,jsub)
          psi( 1,1) = dta(i+1,j,k,jsub)

          psi(-1,2) = dta(i,j-1,k,jsub)
          psi( 0,2) = dta(i,j  ,k,jsub)
          psi( 1,2) = dta(i,j+1,k,jsub)

          psi(-1,3) = dta(i,j,k-1,jsub)
          psi( 0,3) = dta(i,j,k  ,jsub)
          psi( 1,3) = dta(i,j,k+1,jsub)
      ENDIF
#elif __DIM == __2D
      IF (order .EQ. ppm_param_order_3) THEN
          psi(-3,1) = dta(i-3,j,jsub)
          psi(-2,1) = dta(i-2,j,jsub)
          psi(-1,1) = dta(i-1,j,jsub)
          psi( 0,1) = dta(i  ,j,jsub)
          psi( 1,1) = dta(i+1,j,jsub)
          psi( 2,1) = dta(i+2,j,jsub)
          psi( 3,1) = dta(i+3,j,jsub)

          psi(-3,2) = dta(i,j-3,jsub)
          psi(-2,2) = dta(i,j-2,jsub)
          psi(-1,2) = dta(i,j-1,jsub)
          psi( 0,2) = dta(i,j  ,jsub)
          psi( 1,2) = dta(i,j+1,jsub)
          psi( 2,2) = dta(i,j+2,jsub)
          psi( 3,2) = dta(i,j+3,jsub)
      ELSEIF (order .EQ. ppm_param_order_2) THEN
          psi(-2,1) = dta(i-2,j,jsub)
          psi(-1,1) = dta(i-1,j,jsub)
          psi( 0,1) = dta(i  ,j,jsub)
          psi( 1,1) = dta(i+1,j,jsub)
          psi( 2,1) = dta(i+2,j,jsub)

          psi(-2,2) = dta(i,j-2,jsub)
          psi(-1,2) = dta(i,j-1,jsub)
          psi( 0,2) = dta(i,j  ,jsub)
          psi( 1,2) = dta(i,j+1,jsub)
          psi( 2,2) = dta(i,j+2,jsub)
      ELSE
          psi(-1,1) = dta(i-1,j,jsub)
          psi( 0,1) = dta(i  ,j,jsub)
          psi( 1,1) = dta(i+1,j,jsub)

          psi(-1,2) = dta(i,j-1,jsub)
          psi( 0,2) = dta(i,j  ,jsub)
          psi( 1,2) = dta(i,j+1,jsub)
      ENDIF
#endif
