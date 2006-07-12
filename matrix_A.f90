SUBROUTINE matrix_A(Ax,Bx,tstep)

      USE mean_param
      USE surface_forcing
      USE IMEX_constants

      IMPLICIT NONE

      TYPE(trivector), INTENT(OUT)::Ax !integration matrix      
      TYPE(trivector), INTENT(IN)::Bx !integration matrix
      INTEGER, INTENT(IN)::tstep  !time_step

      Ax%A = 0.
      Ax%B = 0.
      AX%C = 0.

         Ax%A = -a_IMEX1*Bx%A
         Ax%B = 1.0 - a_IMEX1*Bx%B
         Ax%C = -a_IMEX1*Bx%C

END SUBROUTINE matrix_A





