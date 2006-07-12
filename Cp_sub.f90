SUBROUTINE Cp_sub(Ti,Si,cp,mm)
     
  USE mean_param

  IMPLICIT NONE

  TYPE(gr_d), INTENT(IN)::mm
  DOUBLE PRECISION, DIMENSION(0:mm%M+1), INTENT(IN)::Ti, Si
  TYPE(constant), INTENT(IN OUT)::cp  

  DOUBLE PRECISION::cp_w
  DOUBLE PRECISION, DIMENSION(0:mm%M+1)::tt
  INTEGER::x
  
  cp%g = 0.
  cp%i = 0.
  tt = Ti-273.15

  DO x = 1, mm%M+1

     cp_w = 4217.4 - 3.720283*tt(x) + 0.1412855*tt(x)**2.0 - 2.654387D-03*tt(x)**3.0 + &
          2.093236D-05*tt(x)**4.0
     cp%g(x) = cp_w + Si(x)*(-7.6444 + 0.107276*tt(x) - 1.3839D-03*tt(x)**2.0) + &
          Si(x)**(3.0/2.0)*(0.17709 - 4.0772D-03*tt(x) + 5.3539D-05*tt(x)**2.0)

  END DO

  cp%g(0) = cp%g(1)
  cp%i(0) = cp%g(1)

  DO x = 1, mm%M
     cp%i(x) = cp%g(x) + (cp%g(x+1)-cp%g(x))*(mm%d_i(x)-mm%d_g(x))/mm%g_space(x)
  END DO

END SUBROUTINE Cp_sub

