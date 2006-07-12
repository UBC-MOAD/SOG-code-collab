SUBROUTINE horizontal_adv(mm,Gx,P,density,factor,dt)

  USE mean_param
  
  IMPLICIT NONE

  TYPE(gr_d),INTENT(IN)::mm      
  DOUBLE PRECISION,DIMENSION(mm%M),INTENT(IN OUT)::Gx
  DOUBLE PRECISION,DIMENSION(mm%M),INTENT(IN)::P  !P_q or -P_f
  DOUBLE PRECISION,DIMENSION(mm%M),INTENT(IN)::density !density%new(1:mm%M)
  DOUBLE PRECISION,DIMENSION(mm%M),INTENT(IN)::factor !Cp(1:mm%M) or 1./S%new(1:mm%M)
  DOUBLE PRECISION, INTENT(IN)::dt

  Gx(1:mm%M) = Gx(1:mm%M) +dt*P(1:mm%M)/density(1:mm%M)/factor(1:mm%M)

END SUBROUTINE horizontal_adv
