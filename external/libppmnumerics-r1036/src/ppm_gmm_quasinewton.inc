      !-------------------------------------------------------------------------
      !  Include file for the quasi-Newton iteration in the initialization
      !  of the high-order marching method.
      !
      !  INPUT:   REAL(MK), DIMENSION(2 or 3) :: x0 -- start value
      !  OUTPUT:  REAL(MK) :: sprod -- new distance to interface
      !
      !-------------------------------------------------------------------------
      !  $Log: ppm_gmm_quasinewton.inc,v $
      !  Revision 1.1.1.1  2007/07/13 10:18:55  ivos
      !  CBL version of the PPM library
      !
      !  Revision 1.6  2005/07/14 19:58:15  ivos
      !  Added OPTIONAL argument chi for mesh node positions in distorted
      !  (mapped) meshes. For use with AGM for example.
      !
      !  Revision 1.5  2005/05/24 23:21:55  ivos
      !  Commented out the convergence warning as I think this can naturally
      !  occur in under-resolved regions...
      !
      !  Revision 1.4  2005/05/10 04:48:46  ivos
      !  Split marching and extension routines for faster compilation,
      !  Sharked extension routines, moved all initialization to gmm_init, and
      !  code cosmetics.
      !
      !  Revision 1.3  2005/04/21 04:53:27  ivos
      !  If convergence fails, an error message is now produced, rather
      !  than returning HUGE. The latter could go undetected and produce
      !  funny results downstream...
      !
      !  Revision 1.2  2005/03/11 04:20:13  ivos
      !  Removed debug output.
      !
      !  Revision 1.1  2005/03/10 01:37:19  ivos
      !  Initial check-in. BEWARE: Not tested in parallel yet!
      !
      !-------------------------------------------------------------------------
      !  Parallel Particle Mesh Library (PPM)
      !  Institute of Computational Science
      !  ETH Zentrum, Hirschengraben 84
      !  CH-8092 Zurich, Switzerland
      !-------------------------------------------------------------------------
#if   __DIM == __3D
      xk   = x0
      err  = big
      tol2 = tol*tol
      nit  = 0
      DO WHILE ((err.GT.tol2).AND.(nit.LT.40))
          ! interpolation polynomial and its gradient at xk
          ! This is a stupid way to evaluate the 
          ! polynomial, but I could not find a better
          ! one immediately...
          xv(1) = xk(1)
          xv(2) = xk(1)*xk(1)
          xv(3) = xv(2)*xk(1)
          yv(1) = xk(2)
          yv(2) = xk(2)*xk(2)
          yv(3) = yv(2)*xk(2)
          zv(1) = xk(3)
          zv(2) = xk(3)*xk(3)
          zv(3) = zv(2)*xk(3)
          gxv(2) = 2.0_MK*xk(1)
          gxv(3) = 3.0_MK*xv(2)
          gyv(2) = 2.0_MK*xk(2)
          gyv(3) = 3.0_MK*yv(2)
          gzv(2) = 2.0_MK*xk(3)
          gzv(3) = 3.0_MK*zv(2)
          pxk = 0.0_MK
          gradpxk = 0.0_MK
          DO p=0,3
              DO n=0,3
                  DO m=0,3
                      ind = 16*p+4*n+m+1
                      pxk = pxk+coef(ind)*xv(m)*yv(n)*zv(p)
                      gradpxk(1) = gradpxk(1) +   &
     &                    (coef(ind)*gxv(m)*yv(n)*zv(p))
                      gradpxk(2) = gradpxk(2) +   &
     &                    (coef(ind)*xv(m)*gyv(n)*zv(p))
                      gradpxk(3) = gradpxk(3) +   &
     &                    (coef(ind)*xv(m)*yv(n)*gzv(p))
                  ENDDO
              ENDDO
          ENDDO
          ! equation (3.4) in Chopp:2001
          nrm2      =         gradpxk(1)*gradpxk(1)
          nrm2      = nrm2 + (gradpxk(2)*gradpxk(2))
          nrm2      = nrm2 + (gradpxk(3)*gradpxk(3))
          nrm2      = 1.0_MK/nrm2
          delta1(1) = -pxk*gradpxk(1)*nrm2
          delta1(2) = -pxk*gradpxk(2)*nrm2
          delta1(3) = -pxk*gradpxk(3)*nrm2
          err       =        delta1(1)*delta1(1)
          err       = err + (delta1(2)*delta1(2))
          err       = err + (delta1(3)*delta1(3))
          ! equation (3.5) in Chopp:2001
          xkhalf(1) = xk(1) + delta1(1)
          xkhalf(2) = xk(2) + delta1(2)
          xkhalf(3) = xk(3) + delta1(3)
          ! equation (3.6) in Chopp:2001
          x0mxk(1)  = x0(1)-xk(1)
          x0mxk(2)  = x0(2)-xk(2)
          x0mxk(3)  = x0(3)-xk(3)
          sprod     =          x0mxk(1)*gradpxk(1)
          sprod     = sprod + (x0mxk(2)*gradpxk(2))
          sprod     = sprod + (x0mxk(3)*gradpxk(3))
          nrm2      = nrm2*sprod
          delta2(1) = x0mxk(1) - nrm2*gradpxk(1)
          delta2(2) = x0mxk(2) - nrm2*gradpxk(2)
          delta2(3) = x0mxk(3) - nrm2*gradpxk(3)
          err       = err + (delta2(1)*delta2(1))
          err       = err + (delta2(2)*delta2(2))
          err       = err + (delta2(3)*delta2(3))
          ! equation (3.7) in Chopp:2001
          xk(1)     = xkhalf(1) + delta2(1)
          xk(2)     = xkhalf(2) + delta2(2)
          xk(3)     = xkhalf(3) + delta2(3)
          nit       = nit + 1
      ENDDO
      IF (nit .GT. 39) THEN
          ! did not converge
    !     info = ppm_error_warning
    !     CALL ppm_error(ppm_err_converge,   &
    !&        'ppm_gmm_kickoff',               &
    !&        'WARNING: Tolerance not reached!'&
    !&        ,__LINE__,info)
    !     info = ppm_param_success
          sprod = big
      ELSE
          ! distance to the interface in real
          ! coordinates 
          IF (PRESENT(chi)) THEN
              !-----------------------------------------------------------------
              !  Map the point to physical space using trilinear
              !  interpolation of the map (AGM uses linear basis functions)
              !-----------------------------------------------------------------
              xi(1) = 1.0_MK-xk(1)
              xi(2) = 1.0_MK-xk(2)
              xi(3) = 1.0_MK-xk(3)
              ch(1) = chi(1,i  ,j  ,k  ,jsub)
              ch(2) = chi(1,ip1,j  ,k  ,jsub)
              ch(3) = chi(1,i  ,jp1,k  ,jsub)
              ch(4) = chi(1,i  ,j  ,kp1,jsub)
              ch(5) = chi(1,ip1,j  ,kp1,jsub)
              ch(6) = chi(1,i  ,jp1,kp1,jsub)
              ch(7) = chi(1,ip1,jp1,k  ,jsub)
              ch(8) = chi(1,ip1,jp1,kp1,jsub)
              xr(1) =         ch(1)*xi(1)*xi(2)*xi(3)
              xr(1) = xr(1) + ch(2)*xk(1)*xi(2)*xi(3)
              xr(1) = xr(1) + ch(3)*xi(1)*xk(2)*xi(3)
              xr(1) = xr(1) + ch(4)*xi(1)*xi(2)*xk(3)
              xr(1) = xr(1) + ch(5)*xk(1)*xi(2)*xk(3)
              xr(1) = xr(1) + ch(6)*xi(1)*xk(2)*xk(3)
              xr(1) = xr(1) + ch(7)*xk(1)*xk(2)*xi(3)
              xr(1) = xr(1) + ch(8)*xk(1)*xk(2)*xk(3)
              ch(1) = chi(2,i  ,j  ,k  ,jsub)
              ch(2) = chi(2,ip1,j  ,k  ,jsub)
              ch(3) = chi(2,i  ,jp1,k  ,jsub)
              ch(4) = chi(2,i  ,j  ,kp1,jsub)
              ch(5) = chi(2,ip1,j  ,kp1,jsub)
              ch(6) = chi(2,i  ,jp1,kp1,jsub)
              ch(7) = chi(2,ip1,jp1,k  ,jsub)
              ch(8) = chi(2,ip1,jp1,kp1,jsub)
              xr(2) =         ch(1)*xi(1)*xi(2)*xi(3)
              xr(2) = xr(2) + ch(2)*xk(1)*xi(2)*xi(3)
              xr(2) = xr(2) + ch(3)*xi(1)*xk(2)*xi(3)
              xr(2) = xr(2) + ch(4)*xi(1)*xi(2)*xk(3)
              xr(2) = xr(2) + ch(5)*xk(1)*xi(2)*xk(3)
              xr(2) = xr(2) + ch(6)*xi(1)*xk(2)*xk(3)
              xr(2) = xr(2) + ch(7)*xk(1)*xk(2)*xi(3)
              xr(2) = xr(2) + ch(8)*xk(1)*xk(2)*xk(3)
              ch(1) = chi(3,i  ,j  ,k  ,jsub)
              ch(2) = chi(3,ip1,j  ,k  ,jsub)
              ch(3) = chi(3,i  ,jp1,k  ,jsub)
              ch(4) = chi(3,i  ,j  ,kp1,jsub)
              ch(5) = chi(3,ip1,j  ,kp1,jsub)
              ch(6) = chi(3,i  ,jp1,kp1,jsub)
              ch(7) = chi(3,ip1,jp1,k  ,jsub)
              ch(8) = chi(3,ip1,jp1,kp1,jsub)
              xr(3) =         ch(1)*xi(1)*xi(2)*xi(3)
              xr(3) = xr(3) + ch(2)*xk(1)*xi(2)*xi(3)
              xr(3) = xr(3) + ch(3)*xi(1)*xk(2)*xi(3)
              xr(3) = xr(3) + ch(4)*xi(1)*xi(2)*xk(3)
              xr(3) = xr(3) + ch(5)*xk(1)*xi(2)*xk(3)
              xr(3) = xr(3) + ch(6)*xi(1)*xk(2)*xk(3)
              xr(3) = xr(3) + ch(7)*xk(1)*xk(2)*xi(3)
              xr(3) = xr(3) + ch(8)*xk(1)*xk(2)*xk(3)
              !-----------------------------------------------------------------
              !  Compute the distance vector in physical space
              !-----------------------------------------------------------------
              x0mxk(1) = x-xr(1)
              x0mxk(2) = y-xr(2)
              x0mxk(3) = z-xr(3)
          ELSE
              !-----------------------------------------------------------------
              !  Undo scaling and compute distance vector
              !-----------------------------------------------------------------
              x0mxk(1)  = dx*(x0(1)-xk(1))
              x0mxk(2)  = dy*(x0(2)-xk(2))
              x0mxk(3)  = dz*(x0(3)-xk(3))
          ENDIF
          sprod     =          x0mxk(1)*x0mxk(1)
          sprod     = sprod + (x0mxk(2)*x0mxk(2))
          sprod     = sprod + (x0mxk(3)*x0mxk(3))
          sprod     = SQRT(sprod)
      ENDIF
#elif  __DIM == __2D
      xk   = x0
      err  = big
      tol2 = tol*tol
      nit  = 0
      DO WHILE ((err.GT.tol2).AND.(nit.LT.40))
          ! interpolation polynomial and its gradient at xk
          ! This is a stupid way to evaluate the 
          ! polynomial, but I could not find a better
          ! one immediately...
          xv(1) = xk(1)
          xv(2) = xk(1)*xk(1)
          xv(3) = xv(2)*xk(1)
          yv(1) = xk(2)
          yv(2) = xk(2)*xk(2)
          yv(3) = yv(2)*xk(2)
          gxv(2) = 2.0_MK*xk(1)
          gxv(3) = 3.0_MK*xv(2)
          gyv(2) = 2.0_MK*xk(2)
          gyv(3) = 3.0_MK*yv(2)
          pxk = 0.0_MK
          gradpxk = 0.0_MK
          DO n=0,3
              DO m=0,3
                  ind = 4*n+m+1
                  pxk = pxk+coef(ind)*xv(m)*yv(n)
                  gradpxk(1) = gradpxk(1) +   &
     &                (coef(ind)*gxv(m)*yv(n))
                  gradpxk(2) = gradpxk(2) +   &
     &                (coef(ind)*xv(m)*gyv(n))
              ENDDO
          ENDDO
          ! equation (3.4) in Chopp:2001
          nrm2      =         gradpxk(1)*gradpxk(1)
          nrm2      = nrm2 + (gradpxk(2)*gradpxk(2))
          nrm2      = 1.0_MK/nrm2
          delta1(1) = -pxk*gradpxk(1)*nrm2
          delta1(2) = -pxk*gradpxk(2)*nrm2
          err       =        delta1(1)*delta1(1)
          err       = err + (delta1(2)*delta1(2))
          ! equation (3.5) in Chopp:2001
          xkhalf(1) = xk(1) + delta1(1)
          xkhalf(2) = xk(2) + delta1(2)
          ! equation (3.6) in Chopp:2001
          x0mxk(1)  = x0(1)-xk(1)
          x0mxk(2)  = x0(2)-xk(2)
          sprod     =          x0mxk(1)*gradpxk(1)
          sprod     = sprod + (x0mxk(2)*gradpxk(2))
          nrm2      = nrm2*sprod
          delta2(1) = x0mxk(1) - nrm2*gradpxk(1)
          delta2(2) = x0mxk(2) - nrm2*gradpxk(2)
          err       = err + (delta2(1)*delta2(1))
          err       = err + (delta2(2)*delta2(2))
          ! equation (3.7) in Chopp:2001
          xk(1)     = xkhalf(1) + delta2(1)
          xk(2)     = xkhalf(2) + delta2(2)
          nit       = nit + 1
      ENDDO
      IF (nit .GT. 39) THEN
          ! did not converge
    !     info = ppm_error_warning
    !     CALL ppm_error(ppm_err_converge,   &
    !&        'ppm_gmm_kickoff',               &
    !&        'WARNING: Tolerance not reached!'&
    !&        ,__LINE__,info)
    !     info = ppm_param_success
          sprod = big
      ELSE
          IF (PRESENT(chi)) THEN
              !-----------------------------------------------------------------
              !  Map the point to physical space using bilinear
              !  interpolation of the map (AGM uses linear basis functions)
              !-----------------------------------------------------------------
              xi(1) = 1.0_MK-xk(1)
              xi(2) = 1.0_MK-xk(2)
              ch(1) = chi(1,i  ,j  ,jsub)
              ch(2) = chi(1,ip1,j  ,jsub)
              ch(3) = chi(1,i  ,jp1,jsub)
              ch(4) = chi(1,ip1,jp1,jsub)
              xr(1) =         ch(1)*xi(1)*xi(2)
              xr(1) = xr(1) + ch(2)*xk(1)*xi(2)
              xr(1) = xr(1) + ch(3)*xi(1)*xk(2)
              xr(1) = xr(1) + ch(4)*xk(1)*xk(2)
              ch(1) = chi(2,i  ,j  ,jsub)
              ch(2) = chi(2,ip1,j  ,jsub)
              ch(3) = chi(2,i  ,jp1,jsub)
              ch(4) = chi(2,ip1,jp1,jsub)
              xr(2) =         ch(1)*xi(1)*xi(2)
              xr(2) = xr(2) + ch(2)*xk(1)*xi(2)
              xr(2) = xr(2) + ch(3)*xi(1)*xk(2)
              xr(2) = xr(2) + ch(4)*xk(1)*xk(2)
              !-----------------------------------------------------------------
              !  Compute the distance vector in physical space
              !-----------------------------------------------------------------
              x0mxk(1) = x-xr(1)
              x0mxk(2) = y-xr(2)
          ELSE
              !-----------------------------------------------------------------
              !  Undo scaling and compute distance vector
              !-----------------------------------------------------------------
              x0mxk(1)  = dx*(x0(1)-xk(1))
              x0mxk(2)  = dy*(x0(2)-xk(2))
          ENDIF
          sprod     =          x0mxk(1)*x0mxk(1)
          sprod     = sprod + (x0mxk(2)*x0mxk(2))
          sprod     = SQRT(sprod)
      ENDIF
#endif
