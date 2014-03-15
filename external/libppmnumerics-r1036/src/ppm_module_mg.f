      !-------------------------------------------------------------------------
      !  Module       :                ppm_module_mg
      !-------------------------------------------------------------------------
      !
      !  Purpose      : This module contains the user-callable functions of
      !                 the mg solver.
      !                
      !  Remarks      :
      !
      !  References   :
      !
      !  Revisions    :
      !-------------------------------------------------------------------------
      !  $Log: ppm_module_mg.f,v $
      !  Revision 1.1.1.1  2007/07/13 10:19:00  ivos
      !  CBL version of the PPM library
      !
      !  Revision 1.1  2004/09/22 18:39:26  kotsalie
      !  MG new version
      !
      !-------------------------------------------------------------------------
      !  Perallel Particle Mesh Library (PPM)
      !  Institute of Computational Science
      !  ETH Zentrum, Hirschengraben 84
      !  CH-8092 Zurich, Switzerland
      !-------------------------------------------------------------------------

      MODULE ppm_module_mg

         !----------------------------------------------------------------------
         !  PPM modules
         !----------------------------------------------------------------------
         USE ppm_module_mg_init
         USE ppm_module_mg_solv
         USE ppm_module_mg_finalize
         USE ppm_module_mg_alloc
         USE ppm_module_mg_core
         USE ppm_module_mg_prolong
         USE ppm_module_mg_res
         USE ppm_module_mg_restrict
         USE ppm_module_mg_smooth
         
      END MODULE ppm_module_mg
