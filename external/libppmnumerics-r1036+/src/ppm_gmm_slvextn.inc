      !-------------------------------------------------------------------------
      !  Include file for determining and solving the linear equation
      !  for orthogonal extension using upwind differences in the unknown
      !  funciton and centered differences for the level function.
      !
      !  INPUT:   INTEGER :: i,j,k -- Point to solve for
      !           INTEGER :: order -- Desired order of FD scheme
      !           REAL(MK), POINTER :: dta -- data
      !           REAL(MK), POINTER :: fdta -- level function
      !  OUTPUT:  REAL(MK):: valijk -- Computed value for point i,j,k
      !
      !-------------------------------------------------------------------------
      !  $Log: ppm_gmm_slvextn.inc,v $
      !  Revision 1.1.1.1  2007/07/13 10:18:55  ivos
      !  CBL version of the PPM library
      !
      !  Revision 1.4  2005/07/25 00:31:10  ivos
      !  bugfix: index errors in jacobian fixed.
      !
      !  Revision 1.3  2005/07/14 19:58:16  ivos
      !  Added OPTIONAL argument chi for mesh node positions in distorted
      !  (mapped) meshes. For use with AGM for example.
      !
      !  Revision 1.2  2005/05/24 23:24:44  ivos
      !  Added checks for uninitialized values in the centered differences
      !  to avoid float overflows in underresolved regions. This is a
      !  technical fix as the result will still be wrong in such cases...
      !
      !  Revision 1.1  2005/05/10 04:41:16  ivos
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

      !-------------------------------------------------------------------------
      !  Determine switches
      !-------------------------------------------------------------------------
#include "ppm_gmm_switches.inc"

      !-------------------------------------------------------------------------
      !  Read level function
      !-------------------------------------------------------------------------
#include "ppm_gmm_getfdta.inc"

      !-------------------------------------------------------------------------
      !  Read value function
      !-------------------------------------------------------------------------
#include "ppm_gmm_getdta.inc"

      !-------------------------------------------------------------------------
      !  Compute upwind finite differences of appropriate
      !  order. 
      !     alpha(1:3) x,y(,z) coefs of the unknown term in dta
      !     beta(1:3) x,y(,z) coefs of the constant term in dta
      !     gphi(1:3) x,y(,z) gradient of the level function 
      !-------------------------------------------------------------------------
      alpha = 0.0_MK
      beta  = 0.0_MK
      gphi  = 0.0_MK

      !-------------------------------------------------------------------------
      !  --- X DIRECTION
      !-------------------------------------------------------------------------
      ! Use upwind direction in level function
      IF ((ABS(phi(-1,1)) .LT. ABS(phi(0,1))) .OR. (psi(1,1) .GT. hsave)) THEN
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
          beta(1)  = psi(i1,1)
          IF (sx(i2) .GT. 0) THEN
              !-----------------------------------------------------------------
              !  Second order in x direction
              !-----------------------------------------------------------------
              alpha(1) = alpha(1) - 0.5_MK
              beta(1)  = beta(1) + psi(i1,1) - 0.5_MK*psi(i2,1)
              IF (sx(i3) .GT. 0) THEN
                  !-------------------------------------------------------------
                  !  Third order in negative x direction
                  !-------------------------------------------------------------
                  alpha(1) = alpha(1) - onethird
                  beta(1)  = beta(1)+psi(i1,1)-psi(i2,1)+onethird*psi(i3,1)
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
      IF ((ABS(phi(-1,2)) .LT. ABS(phi(0,2))) .OR. (psi(1,2) .GT. hsave)) THEN
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
          beta(2)  = psi(i1,2)
          IF (sy(i2) .GT. 0) THEN
              !-----------------------------------------------------------------
              !  Second order in y direction
              !-----------------------------------------------------------------
              alpha(2) = alpha(2) - 0.5_MK
              beta(2)  = beta(2) + psi(i1,2) - 0.5_MK*psi(i2,2)
              IF (sy(i3) .GT. 0) THEN
                  !-------------------------------------------------------------
                  !  Third order in y direction
                  !-------------------------------------------------------------
                  alpha(2) = alpha(2) - onethird
                  beta(2)  = beta(2)+psi(i1,2)-psi(i2,2)+onethird*psi(i3,2)
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
      IF ((ABS(phi(-1,3)) .LT. ABS(phi(0,3))) .OR. (psi(1,3) .GT. hsave)) THEN
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
          beta(3)  = psi(i1,3)
          IF (sz(i2) .GT. 0) THEN
              !-----------------------------------------------------------------
              !  Second order in z direction
              !-----------------------------------------------------------------
              alpha(3) = alpha(3) - 0.5_MK
              beta(3)  = beta(3) + psi(i1,3) - 0.5_MK*psi(i2,3)
              IF (sz(i3) .GT. 0) THEN
                  !-------------------------------------------------------------
                  !  Third order in z direction
                  !-------------------------------------------------------------
                  alpha(3) = alpha(3) - onethird
                  beta(3)  = beta(3)+psi(i1,3)-psi(i2,3)+onethird*psi(i3,3)
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
      !  Compute the gradient of the level function using centered
      !  differences. 
      !-------------------------------------------------------------------------
      IF (order .EQ. ppm_param_order_1) THEN
          IF (ABS(phi(-1,1)).LT.hsave .AND. ABS(phi(1,1)).LT.hsave) THEN
              gphi(1) = phi(1,1)-phi(-1,1)
              gphi(1) = gphi(1)*dxihalf
          ENDIF
          IF (ABS(phi(-1,2)).LT.hsave .AND. ABS(phi(1,2)).LT.hsave) THEN
              gphi(2) = phi(1,2)-phi(-1,2)
              gphi(2) = gphi(2)*dyihalf
          ENDIF
#if   __DIM == __3D
          IF (ABS(phi(-1,3)).LT.hsave .AND. ABS(phi(1,3)).LT.hsave) THEN
              gphi(3) = phi(1,3)-phi(-1,3)
              gphi(3) = gphi(3)*dzihalf
          ENDIF
#endif
      ELSEIF (order .EQ. ppm_param_order_2) THEN
          IF (ABS(phi(-1,1)).LT.hsave .AND. ABS(phi(1,1)).LT.hsave) THEN
              gphi(1) = phi(1,1)-phi(-1,1)
              gphi(1) = gphi(1)*dxihalf
          ENDIF
          IF (ABS(phi(-1,2)).LT.hsave .AND. ABS(phi(1,2)).LT.hsave) THEN
              gphi(2) = phi(1,2)-phi(-1,2)
              gphi(2) = gphi(2)*dyihalf
          ENDIF
#if   __DIM == __3D
          IF (ABS(phi(-1,3)).LT.hsave .AND. ABS(phi(1,3)).LT.hsave) THEN
              gphi(3) = phi(1,3)-phi(-1,3)
              gphi(3) = gphi(3)*dzihalf
          ENDIF
#endif
      ELSE
          IF (ABS(phi(-2,1)).LT.hsave .AND. ABS(phi(-1,1)).LT.hsave .AND.  &
     &        ABS(phi(1,1)).LT.hsave .AND. ABS(phi(2,1)).LT.hsave) THEN
              gphi(1) = phi(-2,1)-8.0_MK*phi(-1,1)+8.0_MK*phi(1,1)-phi(2,1)
              gphi(1) = gphi(1)*dxitwelve
          ENDIF
          IF (ABS(phi(-2,2)).LT.hsave .AND. ABS(phi(-1,2)).LT.hsave .AND.  &
     &        ABS(phi(1,2)).LT.hsave .AND. ABS(phi(2,2)).LT.hsave) THEN
              gphi(2) = phi(-2,2)-8.0_MK*phi(-1,2)+8.0_MK*phi(1,2)-phi(2,2)
              gphi(2) = gphi(2)*dyitwelve
          ENDIF
#if   __DIM == __3D
          IF (ABS(phi(-2,3)).LT.hsave .AND. ABS(phi(-1,3)).LT.hsave .AND.  &
     &        ABS(phi(1,3)).LT.hsave .AND. ABS(phi(2,3)).LT.hsave) THEN
              gphi(3) = phi(-2,3)-8.0_MK*phi(-1,3)+8.0_MK*phi(1,3)-phi(2,3)
              gphi(3) = gphi(3)*dzitwelve
          ENDIF
#endif
      ENDIF

      !-------------------------------------------------------------------------
      !  Compute coefficients for the linear equation for
      !  the value at (i,j,k)
      !-------------------------------------------------------------------------
#if   __DIM == __3D
      IF (PRESENT(chi)) THEN
          !---------------------------------------------------------------------
          !  If we have a non-uniform mesh, compute the Jacobian
          !---------------------------------------------------------------------
#include "ppm_gmm_jacobian.inc"
          !---------------------------------------------------------------------
          !  Grad phi in physical space
          !---------------------------------------------------------------------
          gpp(1) = jac(1,1)*gphi(1) + jac(2,1)*gphi(2) + jac(3,1)*gphi(3)
          gpp(2) = jac(1,2)*gphi(1) + jac(2,2)*gphi(2) + jac(3,2)*gphi(3)
          gpp(3) = jac(1,3)*gphi(1) + jac(2,3)*gphi(2) + jac(3,3)*gphi(3)
          !---------------------------------------------------------------------
          !  Determine the linear equation for the unknown node
          !---------------------------------------------------------------------
          ! x**0 term
          coefs(1) =                                                        &
     &        gpp(1)*(jac(1,1)*beta(1)+jac(2,1)*beta(2)+jac(3,1)*beta(3)) + &
     &        gpp(2)*(jac(1,2)*beta(1)+jac(2,2)*beta(2)+jac(3,2)*beta(3)) + &
     &        gpp(3)*(jac(1,3)*beta(1)+jac(2,3)*beta(2)+jac(3,3)*beta(3))
          IF(PRESENT(speed)) THEN
              coefs(1) = coefs(1) -                              &
     &            1.0_MK/(speed(i,j,k,jsub)*speed(i,j,k,jsub))
          ELSE
              coefs(1) = coefs(1) - rhscst
          ENDIF

          ! x**1 term
          coefs(3) =                                                        &
     &        gpp(1)*(jac(1,1)*alpha(1)+jac(2,1)*alpha(2)+jac(3,1)*alpha(3))+ &
     &        gpp(2)*(jac(1,2)*alpha(1)+jac(2,2)*alpha(2)+jac(3,2)*alpha(3))+ &
     &        gpp(3)*(jac(1,3)*alpha(1)+jac(2,3)*alpha(2)+jac(3,3)*alpha(3))
      ELSE
          !---------------------------------------------------------------------
          !  Uniform mesh
          !---------------------------------------------------------------------
          ! x**0 term
          coefs(1) = (beta(1)*gphi(1))+(beta(2)*gphi(2))+(beta(3)*gphi(3))
          IF(PRESENT(speed)) THEN
              coefs(1) = coefs(1) -                              &
     &            1.0_MK/(speed(i,j,k,jsub)*speed(i,j,k,jsub))
          ELSE
              coefs(1) = coefs(1) - rhscst
          ENDIF

          ! x**1 term
          coefs(3) = (alpha(1)*gphi(1))+(alpha(2)*gphi(2))+(alpha(3)*gphi(3))
      ENDIF
#elif __DIM == __2D
      IF (PRESENT(chi)) THEN
          !---------------------------------------------------------------------
          !  If we have a non-uniform mesh, compute the Jacobian
          !---------------------------------------------------------------------
#include "ppm_gmm_jacobian.inc"
          !---------------------------------------------------------------------
          !  Grad phi in physical space
          !---------------------------------------------------------------------
          gpp(1) = jac(1,1)*gphi(1) + jac(2,1)*gphi(2)
          gpp(2) = jac(1,2)*gphi(1) + jac(2,2)*gphi(2)
          !---------------------------------------------------------------------
          !  Determine the linear equation for the unknown node
          !---------------------------------------------------------------------
          ! x**0 term
          coefs(1) =                                                        &
     &        gpp(1)*(jac(1,1)*beta(1)+jac(2,1)*beta(2)) + &
     &        gpp(2)*(jac(1,2)*beta(1)+jac(2,2)*beta(2))
          IF(PRESENT(speed)) THEN
              coefs(1) = coefs(1) -                              &
     &            1.0_MK/(speed(i,j,jsub)*speed(i,j,jsub))
          ELSE
              coefs(1) = coefs(1) - rhscst
          ENDIF

          ! x**1 term
          coefs(3) =                                                        &
     &        gpp(1)*(jac(1,1)*alpha(1)+jac(2,1)*alpha(2))+ &
     &        gpp(2)*(jac(1,2)*alpha(1)+jac(2,2)*alpha(2))
      ELSE
          !---------------------------------------------------------------------
          !  Uniform mesh
          !---------------------------------------------------------------------
          ! x**0 term
          coefs(1) = (beta(1)*gphi(1))+(beta(2)*gphi(2))
          IF(PRESENT(speed)) THEN
              coefs(1) = coefs(1) -                              &
     &            1.0_MK/(speed(i,j,jsub)*speed(i,j,jsub))
          ELSE
              coefs(1) = coefs(1) - rhscst
          ENDIF

          ! x**1 term
          coefs(3) = (alpha(1)*gphi(1))+(alpha(2)*gphi(2))
      ENDIF
#endif
      !-------------------------------------------------------------------------
      !  Solve for the unknown node value
      !-------------------------------------------------------------------------
      IF (ABS(coefs(3)) .GT. lmyeps) THEN
          valijk = -coefs(1)/coefs(3)
      ELSE
          valijk = big
      ENDIF

