      !-------------------------------------------------------------------------
      !  Include file for computing the jacobian of a mesh map
      !  The indexing is in mathematical ordering with the first index
      !  the row number and the second one the column number. This way
      !  multiplication with the transpose of the 
      !  Jacobian will be stride-1 in memory.
      !
      !  INPUT:   INTEGER :: i,j,k,jsub -- Point to compute at
      !           INTEGER :: order -- Desired order of FD scheme
      !           REAL(MK)(field) :: chi -- mesh node positions
      !  OUTPUT:  REAL(MK)(3,3)   :: jac -- The Jacobian
      !                                     [grad = transp(jac)*gref]
      !
      !-------------------------------------------------------------------------
      !  $Log: ppm_gmm_jacobian.inc,v $
      !  Revision 1.1.1.1  2007/07/13 10:18:55  ivos
      !  CBL version of the PPM library
      !
      !  Revision 1.3  2005/09/06 14:32:07  ivos
      !  Now returns the Jacobian and NOT its transpose. Ordering is (a)ij=a(i,j).
      !
      !  Revision 1.2  2005/07/25 00:31:09  ivos
      !  bugfix: index errors in jacobian fixed.
      !
      !  Revision 1.1  2005/07/14 19:59:59  ivos
      !  Initial implementation.
      !
      !-------------------------------------------------------------------------
      !  Parallel Particle Mesh Library (PPM)
      !  Institute of Computational Science
      !  ETH Zentrum, Hirschengraben 84
      !  CH-8092 Zurich, Switzerland
      !-------------------------------------------------------------------------

#if   __DIM == __3D
      !-------------------------------------------------------------------------
      !  Compute the inverse Jacobian using centered FD of order 2 or 4
      !-------------------------------------------------------------------------
      IF (order .NE. ppm_param_order_3) THEN
          ji(1,1) = dxihalf*(chi(1,i+1,j,k,jsub)-chi(1,i-1,j,k,jsub))
          ji(1,2) = dyihalf*(chi(1,i,j+1,k,jsub)-chi(1,i,j-1,k,jsub))
          ji(1,3) = dzihalf*(chi(1,i,j,k+1,jsub)-chi(1,i,j,k-1,jsub))
          ji(2,1) = dxihalf*(chi(2,i+1,j,k,jsub)-chi(2,i-1,j,k,jsub))
          ji(2,2) = dyihalf*(chi(2,i,j+1,k,jsub)-chi(2,i,j-1,k,jsub))
          ji(2,3) = dzihalf*(chi(2,i,j,k+1,jsub)-chi(2,i,j,k-1,jsub))
          ji(3,1) = dxihalf*(chi(3,i+1,j,k,jsub)-chi(3,i-1,j,k,jsub))
          ji(3,2) = dyihalf*(chi(3,i,j+1,k,jsub)-chi(3,i,j-1,k,jsub))
          ji(3,3) = dzihalf*(chi(3,i,j,k+1,jsub)-chi(3,i,j,k-1,jsub))
      ELSE
          ji(1,1) =                                                   &
     &    dxitwelve*(chi(1,i-2,j,k,jsub)-8.0_MK*chi(1,i-1,j,k,jsub)   &
     &       +8.0_MK*chi(1,i+1,j,k,jsub)-       chi(1,i+2,j,k,jsub))
          ji(1,2) =                                                   &
     &    dyitwelve*(chi(1,i,j-2,k,jsub)-8.0_MK*chi(1,i,j-1,k,jsub)   &
     &       +8.0_MK*chi(1,i,j+1,k,jsub)-       chi(1,i,j+2,k,jsub))
          ji(1,3) =                                                   &
     &    dzitwelve*(chi(1,i,j,k-2,jsub)-8.0_MK*chi(1,i,j,k-1,jsub)   &
     &       +8.0_MK*chi(1,i,j,k+1,jsub)-       chi(1,i,j,k+2,jsub))
          ji(2,1) =                                                   &
     &    dxitwelve*(chi(2,i-2,j,k,jsub)-8.0_MK*chi(2,i-1,j,k,jsub)   &
     &       +8.0_MK*chi(2,i+1,j,k,jsub)-       chi(2,i+2,j,k,jsub))
          ji(2,2) =                                                   &
     &    dyitwelve*(chi(2,i,j-2,k,jsub)-8.0_MK*chi(2,i,j-1,k,jsub)   &
     &       +8.0_MK*chi(2,i,j+1,k,jsub)-       chi(2,i,j+2,k,jsub))
          ji(2,3) =                                                   &
     &    dzitwelve*(chi(2,i,j,k-2,jsub)-8.0_MK*chi(2,i,j,k-1,jsub)   &
     &       +8.0_MK*chi(2,i,j,k+1,jsub)-       chi(2,i,j,k+2,jsub))
          ji(3,1) =                                                   &
     &    dxitwelve*(chi(3,i-2,j,k,jsub)-8.0_MK*chi(3,i-1,j,k,jsub)   &
     &       +8.0_MK*chi(3,i+1,j,k,jsub)-       chi(3,i+2,j,k,jsub))
          ji(3,2) =                                                   &
     &    dyitwelve*(chi(3,i,j-2,k,jsub)-8.0_MK*chi(3,i,j-1,k,jsub)   &
     &       +8.0_MK*chi(3,i,j+1,k,jsub)-       chi(3,i,j+2,k,jsub))
          ji(3,3) =                                                   &
     &    dzitwelve*(chi(3,i,j,k-2,jsub)-8.0_MK*chi(3,i,j,k-1,jsub)   &
     &       +8.0_MK*chi(3,i,j,k+1,jsub)-       chi(3,i,j,k+2,jsub))
      ENDIF
      !-------------------------------------------------------------------------
      !  Invert (grad = Transp(jac)*gref)
      !-------------------------------------------------------------------------
      det = ji(3,1)*ji(1,2)*ji(2,3) - ji(3,1)*ji(1,3)*ji(2,2) -     &
     &      ji(2,1)*ji(1,2)*ji(3,3) + ji(2,1)*ji(1,3)*ji(3,2) +     &
     &      ji(1,1)*ji(2,2)*ji(3,3) - ji(1,1)*ji(2,3)*ji(3,2)
      det = 1.0_MK/det
      jac(1,1) = ( ji(2,2)*ji(3,3) - ji(2,3)*ji(3,2))*det
      jac(2,1) = ( ji(3,1)*ji(2,3) - ji(2,1)*ji(3,3))*det
      jac(3,1) = (-ji(3,1)*ji(2,2) + ji(2,1)*ji(3,2))*det
      jac(1,2) = (-ji(1,2)*ji(3,3) + ji(1,3)*ji(3,2))*det
      jac(2,2) = (-ji(3,1)*ji(1,3) + ji(1,1)*ji(3,3))*det
      jac(3,2) = ( ji(3,1)*ji(1,2) - ji(1,1)*ji(3,2))*det
      jac(1,3) = ( ji(1,2)*ji(2,3) - ji(1,3)*ji(2,2))*det
      jac(2,3) = ( ji(2,1)*ji(1,3) - ji(1,1)*ji(2,3))*det
      jac(3,3) = (-ji(2,1)*ji(1,2) + ji(1,1)*ji(2,2))*det
#elif __DIM == __2D
      !-------------------------------------------------------------------------
      !  Compute the inverse Jacobian using centered FD of order 2 or 4
      !-------------------------------------------------------------------------
      IF (order .NE. ppm_param_order_3) THEN
          ji(1,1) = dxihalf*(chi(1,i+1,j,jsub)-chi(1,i-1,j,jsub))
          ji(1,2) = dyihalf*(chi(1,i,j+1,jsub)-chi(1,i,j-1,jsub))
          ji(2,1) = dxihalf*(chi(2,i+1,j,jsub)-chi(2,i-1,j,jsub))
          ji(2,2) = dyihalf*(chi(2,i,j+1,jsub)-chi(2,i,j-1,jsub))
      ELSE
          ji(1,1) =                                               &
     &    dxitwelve*(chi(1,i-2,j,jsub)-8.0_MK*chi(1,i-1,j,jsub)   &
     &       +8.0_MK*chi(1,i+1,j,jsub)-       chi(1,i+2,j,jsub))
          ji(1,2) =                                               &
     &    dyitwelve*(chi(1,i,j-2,jsub)-8.0_MK*chi(1,i,j-1,jsub)   &
     &       +8.0_MK*chi(1,i,j+1,jsub)-       chi(1,i,j+2,jsub))
          ji(2,1) =                                               &
     &    dxitwelve*(chi(2,i-2,j,jsub)-8.0_MK*chi(2,i-1,j,jsub)   &
     &       +8.0_MK*chi(2,i+1,j,jsub)-       chi(2,i+2,j,jsub))
          ji(2,2) =                                               &
     &    dyitwelve*(chi(2,i,j-2,jsub)-8.0_MK*chi(2,i,j-1,jsub)   &
     &       +8.0_MK*chi(2,i,j+1,jsub)-       chi(2,i,j+2,jsub))
      ENDIF
      !-------------------------------------------------------------------------
      !  Invert
      !-------------------------------------------------------------------------
      det = ji(1,1)*ji(2,2) - ji(1,2)*ji(2,1)
      det = 1.0_MK/det
      jac(1,1) = ji(2,2)*det
      jac(2,1) = -ji(2,1)*det
      jac(1,2) = -ji(1,2)*det
      jac(2,2) = ji(1,1)*det
#endif

