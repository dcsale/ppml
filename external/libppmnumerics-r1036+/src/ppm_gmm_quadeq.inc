      !-------------------------------------------------------------------------
      !  Include file for solving a quadratic equation.
      !
      !  INPUT:   REAL(MK), DIMENSION(3) :: coefs -- coefficients:
      !                  coefs(3)*x**2 + coefs(2)*x + coefs(1) = 0
      !  OUTPUT:  REAL(MK), DIMENSION(2) :: roots -- the two roots.
      !
      !-------------------------------------------------------------------------
      !  $Log: ppm_gmm_quadeq.inc,v $
      !  Revision 1.1.1.1  2007/07/13 10:18:55  ivos
      !  CBL version of the PPM library
      !
      !  Revision 1.3  2005/04/27 01:06:12  ivos
      !  Convergence tests completed, cleaned up code, optmized code (Shark),
      !  and changed structure to allow faster compilation.
      !
      !  Revision 1.2  2005/03/10 01:48:27  ivos
      !  Removed debug output.
      !
      !  Revision 1.1  2005/03/10 01:37:20  ivos
      !  Initial check-in. BEWARE: Not tested in parallel yet!
      !
      !-------------------------------------------------------------------------
      !  Parallel Particle Mesh Library (PPM)
      !  Institute of Computational Science
      !  ETH Zentrum, Hirschengraben 84
      !  CH-8092 Zurich, Switzerland
      !-------------------------------------------------------------------------
      IF (ABS(coefs(3)) .GT. lmyeps) THEN
          ainv = 1.0_MK/coefs(3)
          det  = ainv*ainv
          det  = 0.25_MK*coefs(2)*coefs(2)*det
          det  = det - coefs(1)*ainv
          IF (det .GT. 0.0_MK) THEN
              det = SQRT(det)
              roots(1) = -0.5_MK*coefs(2)*ainv
              roots(2) = roots(1) + det
              roots(1) = roots(1) - det
          ELSE
              roots(1) = big
              roots(2) = big
          ENDIF
      ELSE
          ! linear equation
          IF (ABS(coefs(2)) .GT. lmyeps) THEN
              roots(1) = 0.0_MK
              roots(2) = -coefs(1)/coefs(2)
          ELSE
              roots(1) = big
              roots(2) = big
          ENDIF
      ENDIF
