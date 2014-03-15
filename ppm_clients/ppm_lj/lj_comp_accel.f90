FUNCTION lj_comp_accel(arg)
  USE lj_module_global
  IMPLICIT NONE
  REAL(prec),POINTER                    :: lj_comp_accel
  REAL(prec),DIMENSION(3),INTENT(IN)    :: arg

  comp_accel = 0._prec

END FUNCTION lj_comp_accel
