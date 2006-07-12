SUBROUTINE U_H(mm,Hx,Gx,Gx_o,Gx_o_o,Gx_c,Gx_co,Gx_co_o,Bx_o,Bx_o_o,tstep,vel)                    

      USE mean_param
      USE IMEX_constants

      IMPLICIT NONE

      TYPE(gr_d),INTENT(IN)::mm                 
      DOUBLE PRECISION, DIMENSION(mm%M), INTENT(OUT)::Hx  !zonal vector     
      DOUBLE PRECISION, DIMENSION(mm%M), INTENT(IN)::Gx, Gx_o, Gx_o_o, Gx_c, Gx_co, Gx_co_o
      TYPE(trivector),INTENT(IN)::Bx_o,Bx_o_o
      INTEGER,INTENT(IN)::tstep
      TYPE(prop),INTENT(IN)::vel   !U
 
      INTEGER::index

      Hx = 0.

         
      DO index = 1, mm%M
         Hx(index) = vel%old(index) + (1.0-a_IMEX1)*(Gx_o(index) + Gx_co(index) + &
              Bx_o%A(index)*vel%old(index-1) + &
              Bx_o%B(index)*vel%old(index) + Bx_o%C(index)*vel%old(index+1)) + a_IMEX1*(Gx(index) + &
              Gx_c(index))
      END DO

END SUBROUTINE U_H




