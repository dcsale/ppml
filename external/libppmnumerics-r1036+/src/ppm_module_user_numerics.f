      !-------------------------------------------------------------------------
      !  Module       :                 ppm_module_user_numerics
      !-------------------------------------------------------------------------
      !
      !  Purpose      : This is the global user module. It contains all
      !                 user-callable routines of the entire ppm numerics 
      !                 library.
      !                
      !  Remarks      :
      !
      !  References   :
      !
      !  Revisions    :
      !-------------------------------------------------------------------------
      !  $Log:$
      !-------------------------------------------------------------------------
      !  Perallel Particle Mesh Library (PPM)
      !  Institute of Computational Science
      !  ETH Zentrum, Hirschengraben 84
      !  CH-8092 Zurich, Switzerland
      !-------------------------------------------------------------------------
 
      MODULE ppm_module_user_numerics

         !----------------------------------------------------------------------
         ! PPM numerics routines
         !----------------------------------------------------------------------
         USE ppm_module_comp_part
         USE ppm_module_bem
         USE ppm_module_fieldsolver
         USE ppm_module_ode
         USE ppm_module_mg
         USE ppm_module_fmm
         USE ppm_module_gmm
         USE ppm_module_hamjac

      END MODULE ppm_module_user_numerics
