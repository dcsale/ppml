      !-------------------------------------------------------------------------
      !  Include file for the backsubstitution to solve the linear system
      !  for the interpolation polynomial coefficients based on the LU
      !  factorization of the system matrix.
      !
      !  INPUT:   REAL(MK), DIMENSION(:) :: rhs -- right-hand side
      !           REAL(MK), DIMENSION(:,:) :: A -- LU factorization of
      !                                               system matrix
      !           INTEGER :: Nel -- number of elements in rhs
      !           INTEGER, DIMENSION(:) :: Aind -- permutation index of LU
      !                                            factorization
      !  OUTPUT:  REAL(MK), DIMENSION(:) :: coef -- solution vector
      !
      !  VARIABLES NEEDED: INTEGER :: m,n,ll
      !
      !-------------------------------------------------------------------------
      !  $Log: ppm_gmm_lubksb.inc,v $
      !  Revision 1.1.1.1  2007/07/13 10:18:55  ivos
      !  CBL version of the PPM library
      !
      !  Revision 1.2  2005/03/11 04:18:59  ivos
      !  Cosmetics.
      !
      !  Revision 1.1  2005/03/10 01:37:21  ivos
      !  Initial check-in. BEWARE: Not tested in parallel yet!
      !
      !-------------------------------------------------------------------------
      !  Parallel Particle Mesh Library (PPM)
      !  Institute of Computational Science
      !  ETH Zentrum, Hirschengraben 84
      !  CH-8092 Zurich, Switzerland
      !-------------------------------------------------------------------------
      DO n=1,Nel
          ll = Aind(n)
          coef(n) = rhs(ll)
          DO m=1,n-1
              coef(n) = coef(n) - A(n,m)*coef(m)
          ENDDO
      ENDDO
      coef(Nel) = coef(Nel)/A(Nel,Nel)
      DO n=Nel-1,1,-1
          DO m=n+1,Nel
              coef(n) = coef(n) - A(n,m)*coef(m)
          ENDDO
          coef(n) = coef(n)/A(n,n)
      ENDDO
