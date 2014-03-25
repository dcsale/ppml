      !-------------------------------------------------------------------------
      !  Include file for determining and solving the quadratic equation
      !  using upwind finite differences.
      !
      !  INPUT:   INTEGER :: i,j,k -- Point to solve for
      !           INTEGER :: order -- Desired order of FD scheme
      !           REAL(MK), POINTER :: dta -- data
      !           REAL(MK), POINTER :: fdta -- level function
      !  OUTPUT:  REAL(MK):: valijk -- Computed value for point i,j,k
      !
      !-------------------------------------------------------------------------
      !  $Log: ppm_gmm_slvupwnd.inc,v $
      !  Revision 1.1.1.1  2007/07/13 10:18:55  ivos
      !  CBL version of the PPM library
      !
      !  Revision 1.7  2005/07/25 00:31:09  ivos
      !  bugfix: index errors in jacobian fixed.
      !
      !  Revision 1.6  2005/07/14 19:58:16  ivos
      !  Added OPTIONAL argument chi for mesh node positions in distorted
      !  (mapped) meshes. For use with AGM for example.
      !
      !  Revision 1.5  2005/05/10 04:48:47  ivos
      !  Split marching and extension routines for faster compilation,
      !  Sharked extension routines, moved all initialization to gmm_init, and
      !  code cosmetics.
      !
      !  Revision 1.4  2005/04/27 01:06:13  ivos
      !  Convergence tests completed, cleaned up code, optmized code (Shark),
      !  and changed structure to allow faster compilation.
      !
      !  Revision 1.3  2005/03/16 06:20:09  ivos
      !  Several bugfixes. 1st order version is now tested. Moved all large
      !  data to the module.
      !
      !  Revision 1.2  2005/03/11 04:17:01  ivos
      !  Added possiblity to for velocity extension and to get back
      !  the closest point transform.
      !
      !  Revision 1.1  2005/03/10 01:37:18  ivos
      !  Initial check-in. BEWARE: Not tested in parallel yet!
      !
      !-------------------------------------------------------------------------
      !  Parallel Particle Mesh Library (PPM)
      !  Institute of Computational Science
      !  ETH Zentrum, Hirschengraben 84
      !  CH-8092 Zurich, Switzerland
      !-------------------------------------------------------------------------

      !-------------------------------------------------------------------------
      !  Determine switches
      !-------------------------------------------------------------------------
#include "ppm_gmm_switches.inc"

      !-------------------------------------------------------------------------
      !  Read the level set data
      !-------------------------------------------------------------------------
#include "ppm_gmm_getfdta.inc"

      !-------------------------------------------------------------------------
      !  Compute upwind finite differences of appropriate
      !  order. 
      !     alpha(1:3) x,y(,z) coefs of the unknown term in fdta
      !     beta(1:3) x,y(,z) coefs of the constant term in fdta
      !-------------------------------------------------------------------------
      alpha = 0.0_MK
      beta  = 0.0_MK

      !-------------------------------------------------------------------------
      !  --- X DIRECTION
      !-------------------------------------------------------------------------
      ! Use upwind direction
      IF ((ABS(phi(-1,1)) .LT. ABS(phi(0,1))) .OR. (phi(1,1) .GT. hsave)) THEN
          i1 = -1
          i2 = -2
          i3 = -3
      ELSE
          i1 = 1
          i2 = 2
          i3 = 3
      ENDIF
      IF (sx(i1) .GT. 0) THEN
          !---------------------------------------------------------------------
          !  First order in x direction
          !---------------------------------------------------------------------
          alpha(1) = -1.0_MK
          beta(1)  = phi(i1,1)
          IF (sx(i2) .GT. 0) THEN
              !-----------------------------------------------------------------
              !  Second order in x direction
              !-----------------------------------------------------------------
              alpha(1) = alpha(1) - 0.5_MK
              beta(1)  = beta(1) + phi(i1,1) - 0.5_MK*phi(i2,1)
              IF (sx(i3) .GT. 0) THEN
                  !-------------------------------------------------------------
                  !  Third order in negative x direction
                  !-------------------------------------------------------------
                  alpha(1) = alpha(1) - onethird
                  beta(1)  = beta(1)+phi(i1,1)-phi(i2,1)+onethird*phi(i3,1)
              ENDIF
          ENDIF
      ENDIF

      IF (i1 .LT. 0) THEN 
          alpha(1) = -alpha(1)*dxinv
          beta(1)  = -beta(1)*dxinv
      ELSE
          alpha(1) = alpha(1)*dxinv
          beta(1)  = beta(1)*dxinv
      ENDIF

      !-------------------------------------------------------------------------
      !  --- Y DIRECTION
      !-------------------------------------------------------------------------
      ! Use upwind direction
      IF ((ABS(phi(-1,2)) .LT. ABS(phi(0,2))) .OR. (phi(1,2) .GT. hsave)) THEN
          i1 = -1
          i2 = -2
          i3 = -3
      ELSE
          i1 = 1
          i2 = 2
          i3 = 3
      ENDIF
      IF (sy(i1) .GT. 0) THEN
          !---------------------------------------------------------------------
          !  First order in y direction
          !---------------------------------------------------------------------
          alpha(2) = -1.0_MK
          beta(2)  = phi(i1,2)
          IF (sy(i2) .GT. 0) THEN
              !-----------------------------------------------------------------
              !  Second order in y direction
              !-----------------------------------------------------------------
              alpha(2) = alpha(2) - 0.5_MK
              beta(2)  = beta(2) + phi(i1,2) - 0.5_MK*phi(i2,2)
              IF (sy(i3) .GT. 0) THEN
                  !-------------------------------------------------------------
                  !  Third order in y direction
                  !-------------------------------------------------------------
                  alpha(2) = alpha(2) - onethird
                  beta(2)  = beta(2)+phi(i1,2)-phi(i2,2)+onethird*phi(i3,2)
              ENDIF
          ENDIF
      ENDIF

      IF (i1 .LT. 0) THEN
          alpha(2) = -alpha(2)*dyinv
          beta(2)  = -beta(2)*dyinv
      ELSE
          alpha(2) = alpha(2)*dyinv
          beta(2)  = beta(2)*dyinv
      ENDIF

#if    __DIM == __3D
      !-------------------------------------------------------------------------
      !  --- Z DIRECTION
      !-------------------------------------------------------------------------
      ! Use upwind direction
      IF ((ABS(phi(-1,3)) .LT. ABS(phi(0,3))) .OR. (phi(1,3) .GT. hsave)) THEN
          i1 = -1
          i2 = -2
          i3 = -3
      ELSE
          i1 = 1
          i2 = 2
          i3 = 3
      ENDIF
      IF (sz(i1) .GT. 0) THEN
          !---------------------------------------------------------------------
          !  First order in z direction
          !---------------------------------------------------------------------
          alpha(3) = -1.0_MK
          beta(3)  = phi(i1,3)
          IF (sz(i2) .GT. 0) THEN
              !-----------------------------------------------------------------
              !  Second order in z direction
              !-----------------------------------------------------------------
              alpha(3) = alpha(3) - 0.5_MK
              beta(3)  = beta(3) + phi(i1,3) - 0.5_MK*phi(i2,3)
              IF (sz(i3) .GT. 0) THEN
                  !-------------------------------------------------------------
                  !  Third order in z direction
                  !-------------------------------------------------------------
                  alpha(3) = alpha(3) - onethird
                  beta(3)  = beta(3)+phi(i1,3)-phi(i2,3)+onethird*phi(i3,3)
              ENDIF
          ENDIF
      ENDIF

      IF (i1 .LT. 0) THEN
          alpha(3) = -alpha(3)*dzinv
          beta(3)  = -beta(3)*dzinv
      ELSE
          alpha(3) = alpha(3)*dzinv
          beta(3)  = beta(3)*dzinv
      ENDIF
#endif

      !-------------------------------------------------------------------------
      !  Compute coefficients for the quadratic equation for
      !  the value at (i,j,k)
      !-------------------------------------------------------------------------
#if   __DIM == __3D
      IF (PRESENT(chi)) THEN
          !---------------------------------------------------------------------
          !  If we have a non-uniform mesh, compute the Jacobian
          !---------------------------------------------------------------------
#include "ppm_gmm_jacobian.inc"
          !---------------------------------------------------------------------
          !  Determine the quadratic equation for the unknown node
          !---------------------------------------------------------------------
          ! x**0 term
          coefs(1) =                                                        &
     &       (jac(1,1)*beta(1)+jac(2,1)*beta(2)+jac(3,1)*beta(3))**2.0_MK+  &
     &       (jac(1,2)*beta(1)+jac(2,2)*beta(2)+jac(3,2)*beta(3))**2.0_MK+  &
     &       (jac(1,3)*beta(1)+jac(2,3)*beta(2)+jac(3,3)*beta(3))**2.0_MK
          IF(PRESENT(speed)) THEN
              coefs(1) = coefs(1) - 1.0_MK/(speed(i,j,k,jsub)*speed(i,j,k,jsub))
          ELSE
              coefs(1) = coefs(1) - rhscst
          ENDIF

          ! x**2 term
          coefs(3) =                                                          &
     &       (jac(1,1)*alpha(1)+jac(2,1)*alpha(2)+jac(3,1)*alpha(3))**2.0_MK+ &
     &       (jac(1,2)*alpha(1)+jac(2,2)*alpha(2)+jac(3,2)*alpha(3))**2.0_MK+ &
     &       (jac(1,3)*alpha(1)+jac(2,3)*alpha(2)+jac(3,3)*alpha(3))**2.0_MK

          ! x**1 term
          coefs(2) =                                                     &
             (jac(1,1)*beta(1) +jac(2,1)*beta(2) +jac(3,1)*beta(3) ) *   &
     &       (jac(1,1)*alpha(1)+jac(2,1)*alpha(2)+jac(3,1)*alpha(3))
          coefs(2) = coefs(2)                                        +   &
             (jac(1,2)*beta(1) +jac(2,2)*beta(2) +jac(3,2)*beta(3) ) *   &
     &       (jac(1,2)*alpha(1)+jac(2,2)*alpha(2)+jac(3,2)*alpha(3))
          coefs(2) = coefs(2)                                        +   &
             (jac(1,3)*beta(1) +jac(2,3)*beta(2) +jac(3,3)*beta(3) ) *   &
     &       (jac(1,3)*alpha(1)+jac(2,3)*alpha(2)+jac(3,3)*alpha(3))
          coefs(2) = 2.0_MK*coefs(2)
      ELSE
          !---------------------------------------------------------------------
          !  Uniform mesh
          !---------------------------------------------------------------------
          ! x**0 term
          coefs(1) = (beta(1)*beta(1))+(beta(2)*beta(2))+(beta(3)*beta(3))
          IF(PRESENT(speed)) THEN
              coefs(1) = coefs(1) - 1.0_MK/(speed(i,j,k,jsub)*speed(i,j,k,jsub))
          ELSE
              coefs(1) = coefs(1) - rhscst
          ENDIF

          ! x**2 term
          coefs(3) = (alpha(1)*alpha(1))+(alpha(2)*alpha(2))+(alpha(3)*alpha(3))

          ! x**1 term
          coefs(2) = (alpha(1)*beta(1))+(alpha(2)*beta(2))+(alpha(3)*beta(3))
          coefs(2) = 2.0_MK*coefs(2)
      ENDIF
#elif __DIM == __2D
      IF (PRESENT(chi)) THEN
          !---------------------------------------------------------------------
          !  If we have a non-uniform mesh, compute the Jacobian
          !---------------------------------------------------------------------
#include "ppm_gmm_jacobian.inc"
          !---------------------------------------------------------------------
          !  Determine the quadratic equation for the unknown node
          !---------------------------------------------------------------------
          ! x**0 term
          coefs(1) =                                       &
     &       (jac(1,1)*beta(1)+jac(2,1)*beta(2))**2.0_MK+  &
     &       (jac(1,2)*beta(1)+jac(2,2)*beta(2))**2.0_MK
          IF(PRESENT(speed)) THEN
              coefs(1) = coefs(1) - 1.0_MK/(speed(i,j,jsub)*speed(i,j,jsub))
          ELSE
              coefs(1) = coefs(1) - rhscst
          ENDIF

          ! x**2 term
          coefs(3) =                                        &
     &       (jac(1,1)*alpha(1)+jac(2,1)*alpha(2))**2.0_MK+ &
     &       (jac(1,2)*alpha(1)+jac(2,2)*alpha(2))**2.0_MK

          ! x**1 term
          coefs(2) =                                   &
             (jac(1,1)*beta(1) +jac(2,1)*beta(2) ) *   &
     &       (jac(1,1)*alpha(1)+jac(2,1)*alpha(2))
          coefs(2) = coefs(2)                      +   &
             (jac(1,2)*beta(1) +jac(2,2)*beta(2) ) *   &
     &       (jac(1,2)*alpha(1)+jac(2,2)*alpha(2))
          coefs(2) = 2.0_MK*coefs(2)
      ELSE
          !---------------------------------------------------------------------
          !  Uniform mesh
          !---------------------------------------------------------------------
          ! x**0 term
          coefs(1) = (beta(1)*beta(1))+(beta(2)*beta(2))
          IF(PRESENT(speed)) THEN
              coefs(1) = coefs(1) - 1.0_MK/(speed(i,j,jsub)*speed(i,j,jsub))
          ELSE
              coefs(1) = coefs(1) - rhscst
          ENDIF

          ! x**2 term
          coefs(3) = (alpha(1)*alpha(1))+(alpha(2)*alpha(2))

          ! x**1 term
          coefs(2) = (alpha(1)*beta(1))+(alpha(2)*beta(2))
          coefs(2) = 2.0_MK*coefs(2)
      ENDIF
#endif

      !-------------------------------------------------------------------------
      !  Solve quadratic equation for the unknown node value
      !  Cannot use ppm_util_quadeq_real because this equation can have
      !  complex roots. Result in roots(1:2)
      !-------------------------------------------------------------------------
#include "ppm_gmm_quadeq.inc"

      !-------------------------------------------------------------------------
      !  Update value using the larger of the two roots
      !-------------------------------------------------------------------------
      IF (ABS(roots(1)) .GT. ABS(roots(2))) THEN
          valijk = roots(1)
      ELSE
          valijk = roots(2)
      ENDIF
