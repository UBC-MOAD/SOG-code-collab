SUBROUTINE scalar_H(mm,Hx,Gx,Gx_o,Gx_o_o,Bx_o,Bx_o_o,tstep,scalar)     

      USE mean_param
      USE surface_forcing
      USE IMEX_constants

      IMPLICIT NONE

      TYPE(gr_d),INTENT(IN)::mm
      DOUBLE PRECISION, DIMENSION(mm%M), INTENT(IN OUT)::Hx  !scalar vector
      DOUBLE PRECISION, DIMENSION(mm%M), INTENT(IN)::Gx, Gx_o, Gx_o_o
      TYPE(trivector),INTENT(IN)::Bx_o,Bx_o_o
      INTEGER,INTENT(IN)::tstep
      TYPE(prop),INTENT(IN)::scalar  !eg. T
      
         
      INTEGER::index

      Hx = 0.

      DO index = 1, mm%M
         Hx(index) = scalar%old(index) + (1.0-a_IMEX1)*(Gx_o(index) + Bx_o%A(index)*scalar%old(index-1)+&
              Bx_o%B(index)*scalar%old(index) + Bx_o%C(index)*scalar%old(index+1)) + a_IMEX1*Gx(index)
      END DO

    END SUBROUTINE scalar_H






