      !-------------------------------------------------------------------------
      !  Module       :                     ppm numerics header file
      !-------------------------------------------------------------------------
      !
      !-------------------------------------------------------------------------
      !  Parallel Particle Mesh Library (PPM)
      !  Institute of Computational Science
      !  ETH Zentrum, Hirschengraben 84
      !  CH-8092 Zurich, Switzerland
      !-------------------------------------------------------------------------

      !-------------------------------------------------------------------------
      !  Define quadrature rules for the unit triangle (BEM)
      !-------------------------------------------------------------------------
      INTEGER, PARAMETER :: ppm_param_bem_quadrule_center = 1
      INTEGER, PARAMETER :: ppm_param_bem_quadrule_nodes  = 2
      INTEGER, PARAMETER :: ppm_param_bem_quadrule_edges  = 3
      INTEGER, PARAMETER :: ppm_param_bem_quadrule_cne    = 4
      INTEGER, PARAMETER :: ppm_param_bem_quadrule_stroud = 5
      INTEGER, PARAMETER :: ppm_param_bem_quadrule_hammer3= 6
      INTEGER, PARAMETER :: ppm_param_bem_quadrule_hammer4= 7
      INTEGER, PARAMETER :: ppm_param_bem_quadrule_hammer7= 8
      INTEGER, PARAMETER :: ppm_param_bem_quadrule_hammer12= 9

      !-------------------------------------------------------------------------
      !  Define basis functions for the unit triangle (BEM)
      !-------------------------------------------------------------------------
      INTEGER, PARAMETER :: ppm_param_bem_basis_const     = 1
      INTEGER, PARAMETER :: ppm_param_bem_basis_linear    = 2
      INTEGER, PARAMETER :: ppm_param_bem_basis_quad      = 3

      !-------------------------------------------------------------------------
      !  Define field solvers
      !-------------------------------------------------------------------------
      INTEGER, PARAMETER :: ppm_param_solver_mg           = 0
      INTEGER, PARAMETER :: ppm_param_solver_sor          = 1
      INTEGER, PARAMETER :: ppm_param_solver_fft          = 2
      
      !------------------------------------------------------------------------
      !  Define equation ppm solvers 
      !------------------------------------------------------------------------
      INTEGER, PARAMETER :: ppm_param_eq_poisson          = 1 

      !------------------------------------------------------------------------
      !  Define order for finite difference
      !------------------------------------------------------------------------ 
      INTEGER, PARAMETER :: ppm_param_order_1             = 1
      INTEGER, PARAMETER :: ppm_param_order_2             = 2
      INTEGER, PARAMETER :: ppm_param_order_3             = 3
      INTEGER, PARAMETER :: ppm_param_order_4             = 4

      !------------------------------------------------------------------------
      !  Define smoother
      !------------------------------------------------------------------------
      INTEGER, PARAMETER :: ppm_param_smooth_rbsor        = 1 

      !-------------------------------------------------------------------------
      !  ODE parameters
      !-------------------------------------------------------------------------
      INTEGER, PARAMETER :: ppm_param_ode_scheme_eulerf   = 1
      INTEGER, PARAMETER :: ppm_param_ode_scheme_tvdrk2   = 2
      INTEGER, PARAMETER :: ppm_param_ode_scheme_midrk2   = 3
      INTEGER, PARAMETER :: ppm_param_ode_scheme_rk4      = 4
      INTEGER, PARAMETER :: ppm_param_ode_scheme_trapez   = 5
      INTEGER, PARAMETER :: ppm_param_ode_scheme_tvdrk3   = 6
      INTEGER, PARAMETER :: ppm_param_ode_scheme_sts      = 7

      
