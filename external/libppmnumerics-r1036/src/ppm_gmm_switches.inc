      !-------------------------------------------------------------------------
      !  Include file for determining the switches for the higher order 
      !  differences.
      !
      !  INPUT:   INTEGER :: i,j,k -- Point to solve for
      !           INTEGER :: order -- Desired order of FD scheme
      !  OUTPUT:  INTEGER(-order:order) :: sx,sy,sz -- gmm_state of the 
      !                                    neighbors in each direction
      !
      !-------------------------------------------------------------------------
      !  $Log: ppm_gmm_switches.inc,v $
      !  Revision 1.1.1.1  2007/07/13 10:18:55  ivos
      !  CBL version of the PPM library
      !
      !  Revision 1.2  2005/05/10 04:48:47  ivos
      !  Split marching and extension routines for faster compilation,
      !  Sharked extension routines, moved all initialization to gmm_init, and
      !  code cosmetics.
      !
      !  Revision 1.1  2005/04/27 01:08:39  ivos
      !  Initial commit, but tested.
      !
      !-------------------------------------------------------------------------
      !  Parallel Particle Mesh Library (PPM)
      !  Institute of Computational Science
      !  ETH Zentrum, Hirschengraben 84
      !  CH-8092 Zurich, Switzerland
      !-------------------------------------------------------------------------

#if   __DIM == __3D
      IF (order .EQ. ppm_param_order_3) THEN
          sx     = 0
          sx(-3) = gmm_state3d(i-3,j,k,jsub)
          sx(-2) = gmm_state3d(i-2,j,k,jsub)
          sx(-1) = gmm_state3d(i-1,j,k,jsub)
          sx( 0) = gmm_state3d(i  ,j,k,jsub)
          sx( 1) = gmm_state3d(i+1,j,k,jsub)
          sx( 2) = gmm_state3d(i+2,j,k,jsub)
          sx( 3) = gmm_state3d(i+3,j,k,jsub)

          sy     = 0
          sy(-3) = gmm_state3d(i,j-3,k,jsub)
          sy(-2) = gmm_state3d(i,j-2,k,jsub)
          sy(-1) = gmm_state3d(i,j-1,k,jsub)
          sy( 0) = gmm_state3d(i,j  ,k,jsub)
          sy( 1) = gmm_state3d(i,j+1,k,jsub)
          sy( 2) = gmm_state3d(i,j+2,k,jsub)
          sy( 3) = gmm_state3d(i,j+3,k,jsub)

          sz     = 0
          sz(-3) = gmm_state3d(i,j,k-3,jsub)
          sz(-2) = gmm_state3d(i,j,k-2,jsub)
          sz(-1) = gmm_state3d(i,j,k-1,jsub)
          sz( 0) = gmm_state3d(i,j,k  ,jsub)
          sz( 1) = gmm_state3d(i,j,k+1,jsub)
          sz( 2) = gmm_state3d(i,j,k+2,jsub)
          sz( 3) = gmm_state3d(i,j,k+3,jsub)
      ELSEIF (order .EQ. ppm_param_order_2) THEN
          sx     = 0
          sx(-2) = gmm_state3d(i-2,j,k,jsub)
          sx(-1) = gmm_state3d(i-1,j,k,jsub)
          sx( 0) = gmm_state3d(i  ,j,k,jsub)
          sx( 1) = gmm_state3d(i+1,j,k,jsub)
          sx( 2) = gmm_state3d(i+2,j,k,jsub)

          sy     = 0
          sy(-2) = gmm_state3d(i,j-2,k,jsub)
          sy(-1) = gmm_state3d(i,j-1,k,jsub)
          sy( 0) = gmm_state3d(i,j  ,k,jsub)
          sy( 1) = gmm_state3d(i,j+1,k,jsub)
          sy( 2) = gmm_state3d(i,j+2,k,jsub)

          sz     = 0
          sz(-2) = gmm_state3d(i,j,k-2,jsub)
          sz(-1) = gmm_state3d(i,j,k-1,jsub)
          sz( 0) = gmm_state3d(i,j,k  ,jsub)
          sz( 1) = gmm_state3d(i,j,k+1,jsub)
          sz( 2) = gmm_state3d(i,j,k+2,jsub)
      ELSE
          sx     = 0
          sx(-1) = gmm_state3d(i-1,j,k,jsub)
          sx( 0) = gmm_state3d(i  ,j,k,jsub)
          sx( 1) = gmm_state3d(i+1,j,k,jsub)

          sy     = 0
          sy(-1) = gmm_state3d(i,j-1,k,jsub)
          sy( 0) = gmm_state3d(i,j  ,k,jsub)
          sy( 1) = gmm_state3d(i,j+1,k,jsub)

          sz     = 0
          sz(-1) = gmm_state3d(i,j,k-1,jsub)
          sz( 0) = gmm_state3d(i,j,k  ,jsub)
          sz( 1) = gmm_state3d(i,j,k+1,jsub)
      ENDIF
#elif __DIM == __2D
      IF (order .EQ. ppm_param_order_3) THEN
          sx     = 0
          sx(-3) = gmm_state2d(i-3,j,jsub)
          sx(-2) = gmm_state2d(i-2,j,jsub)
          sx(-1) = gmm_state2d(i-1,j,jsub)
          sx( 0) = gmm_state2d(i  ,j,jsub)
          sx( 1) = gmm_state2d(i+1,j,jsub)
          sx( 2) = gmm_state2d(i+2,j,jsub)
          sx( 3) = gmm_state2d(i+3,j,jsub)

          sy     = 0
          sy(-3) = gmm_state2d(i,j-3,jsub)
          sy(-2) = gmm_state2d(i,j-2,jsub)
          sy(-1) = gmm_state2d(i,j-1,jsub)
          sy( 0) = gmm_state2d(i,j  ,jsub)
          sy( 1) = gmm_state2d(i,j+1,jsub)
          sy( 2) = gmm_state2d(i,j+2,jsub)
          sy( 3) = gmm_state2d(i,j+3,jsub)
      ELSEIF (order .EQ. ppm_param_order_2) THEN
          sx     = 0
          sx(-2) = gmm_state2d(i-2,j,jsub)
          sx(-1) = gmm_state2d(i-1,j,jsub)
          sx( 0) = gmm_state2d(i  ,j,jsub)
          sx( 1) = gmm_state2d(i+1,j,jsub)
          sx( 2) = gmm_state2d(i+2,j,jsub)

          sy     = 0
          sy(-2) = gmm_state2d(i,j-2,jsub)
          sy(-1) = gmm_state2d(i,j-1,jsub)
          sy( 0) = gmm_state2d(i,j  ,jsub)
          sy( 1) = gmm_state2d(i,j+1,jsub)
          sy( 2) = gmm_state2d(i,j+2,jsub)
      ELSE
          sx     = 0
          sx(-1) = gmm_state2d(i-1,j,jsub)
          sx( 0) = gmm_state2d(i  ,j,jsub)
          sx( 1) = gmm_state2d(i+1,j,jsub)

          sy     = 0
          sy(-1) = gmm_state2d(i,j-1,jsub)
          sy( 0) = gmm_state2d(i,j  ,jsub)
          sy( 1) = gmm_state2d(i,j+1,jsub)
      ENDIF
#endif
