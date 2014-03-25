      !-------------------------------------------------------------------------
      !  Module       :                 ppm_module_gmm
      !-------------------------------------------------------------------------
      !
      !  Purpose      : This module contains all user-callable routines
      !                 needed for the group marching method.
      !                
      !  Remarks      :
      !
      !  References   :
      !
      !  Revisions    :
      !-------------------------------------------------------------------------
      !  $Log: ppm_module_gmm.f,v $
      !  Revision 1.1.1.1  2007/07/13 10:18:58  ivos
      !  CBL version of the PPM library
      !
      !  Revision 1.3  2005/03/12 04:08:36  ivos
      !  Misc bug fixes.
      !
      !  Revision 1.2  2005/03/11 04:17:59  ivos
      !  Added ppm_gmm_extend and ppm_gmm_reinitialize.
      !
      !  Revision 1.1  2005/03/10 01:37:16  ivos
      !  Initial check-in. BEWARE: Not tested in parallel yet!
      !
      !-------------------------------------------------------------------------
      !  Perallel Particle Mesh Library (PPM)
      !  Institute of Computational Science
      !  ETH Zentrum, Hirschengraben 84
      !  CH-8092 Zurich, Switzerland
      !-------------------------------------------------------------------------

      MODULE ppm_module_gmm

         !----------------------------------------------------------------------
         !  PPM modules
         !----------------------------------------------------------------------
         USE ppm_module_gmm_init
         USE ppm_module_gmm_kickoff
         USE ppm_module_gmm_march
         USE ppm_module_gmm_cpt
         USE ppm_module_gmm_reinitialize
         USE ppm_module_gmm_extend
         USE ppm_module_gmm_finalize
         
      END MODULE ppm_module_gmm
