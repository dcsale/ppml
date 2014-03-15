      !-------------------------------------------------------------------------
      !  Module       :               ppm_module_data_gmm
      !-------------------------------------------------------------------------
      !
      !  Purpose      : This module contains all data structures and
      !                 definitions that are PRIVATE to the group 
      !                 marching method routines.
      !                
      !  Remarks      :
      !
      !  References   :
      !
      !  Revisions    :
      !-------------------------------------------------------------------------
      !  $Log: ppm_module_data_gmm.f,v $
      !  Revision 1.1.1.1  2007/07/13 10:18:57  ivos
      !  CBL version of the PPM library
      !
      !  Revision 1.5  2006/07/03 12:56:56  ivos
      !  Added comments to explain all the variables.
      !
      !  Revision 1.4  2005/04/27 01:06:14  ivos
      !  Convergence tests completed, cleaned up code, optmized code (Shark),
      !  and changed structure to allow faster compilation.
      !
      !  Revision 1.3  2005/04/21 04:48:25  ivos
      !  Cleaned interfaces and removed unnecessary overloaded versions.
      !
      !  Revision 1.2  2005/03/16 06:20:10  ivos
      !  Several bugfixes. 1st order version is now tested. Moved all large
      !  data to the module.
      !
      !  Revision 1.1  2005/03/10 01:37:17  ivos
      !  Initial check-in. BEWARE: Not tested in parallel yet!
      !
      !-------------------------------------------------------------------------
      !  Perallel Particle Mesh Library (PPM)
      !  Institute of Computational Science
      !  ETH Zentrum, Hirschengraben 84
      !  CH-8092 Zurich, Switzerland
      !-------------------------------------------------------------------------

      MODULE ppm_module_data_gmm

         !----------------------------------------------------------------------
         !  Modules
         !----------------------------------------------------------------------
         USE ppm_module_data, ONLY: ppm_kind_single, ppm_kind_double 
         PRIVATE :: ppm_kind_single, ppm_kind_double
         
         !----------------------------------------------------------------------
         !  Parameters
         !----------------------------------------------------------------------
         ! FLAGS for the three different types of points in GMM.
         ! It is important that far is 0 and far<close<accepted
         ! Otherwise the computation of the switches in gmm_slvupwnd needs to
         ! be adjusted !
         INTEGER, PARAMETER      :: ppm_gmm_param_far      = 0
         INTEGER, PARAMETER      :: ppm_gmm_param_close    = 1
         INTEGER, PARAMETER      :: ppm_gmm_param_accepted = 2

         !----------------------------------------------------------------------
         !  Global values
         !----------------------------------------------------------------------
         ! GMM step size and max extent of mesh in all directions
         INTEGER                  :: incr,maxxhi,maxyhi,maxzhi
         ! topoID and meshID on which the GMM operates
         INTEGER                  :: gmm_topoid,gmm_meshid
         
         !----------------------------------------------------------------------
         !  Sparse matrix structure as workspace
         !----------------------------------------------------------------------
         ! locations (mesh indices) of the data points. Sparse structure.
         INTEGER              , DIMENSION(:,:)  , POINTER   :: gmm_ipos
         ! values at these locations. Sparse structure.
         REAL(ppm_kind_double), DIMENSION(:  )  , POINTER   :: gmm_phid
         REAL(ppm_kind_single), DIMENSION(:  )  , POINTER   :: gmm_phis
         ! length of the sparse structure (i.e. number of points in ipos)
         INTEGER                                            :: gmm_lsiz = -1
         ! positions of the closest points on the interface (for CPT)
         REAL(ppm_kind_double), DIMENSION(:,:)  , POINTER   :: gmm_clod
         REAL(ppm_kind_single), DIMENSION(:,:)  , POINTER   :: gmm_clos
         REAL(ppm_kind_double), DIMENSION(:,:)  , POINTER   :: gmm_clod2
         REAL(ppm_kind_single), DIMENSION(:,:)  , POINTER   :: gmm_clos2
         ! states of the mesh points (PARAMETER flags above)
         INTEGER              , DIMENSION(:,:,:,:), POINTER :: gmm_state3d
         INTEGER              , DIMENSION(:,:,:), POINTER   :: gmm_state2d
         ! list of mesh points closest to the interface
         INTEGER              , DIMENSION(:,:)  , POINTER   :: iptstmp
         INTEGER              , DIMENSION(:,:)  , POINTER   :: iptstmp2
         ! index and sort key for ranking of points
         INTEGER              , DIMENSION(:  )  , POINTER   :: idx,key

      END MODULE ppm_module_data_gmm
