SUBROUTINE alpha_sub(Ti, Si, alp, d)

      USE mean_param

      IMPLICIT NONE

      TYPE(constant), INTENT(IN OUT)::alp
      TYPE(gr_d), INTENT(IN)::d
      DOUBLE PRECISION, DIMENSION(0:d%M+1), INTENT(IN)::Ti,Si !current T and S

      INTEGER:: i, k
      DOUBLE PRECISION::div
      
      alp%g = 0.0
      alp%i = 0.0

      DO k = 1, d%M+1
         IF (Ti(k) >= alp%data(1)%TT) THEN
               alp%g(k) = alp%data(1)%value
               alp%g(k) = alp%g(k) + alp%data(1)%dS*(Si(k)-alp%data(1)%SS)
         ELSE 
            DO i = 2,13
               IF (alp%data(i)%TT <= Ti(k)) THEN
                  alp%g(k) = alp%data(i)%value +(alp%data(i-1)%value- &
                  alp%data(i)%value)*&
                  (Ti(k)-alp%data(i)%TT)/(alp%data(i-1)%TT-alp%data(i)%TT)
                  div = alp%data(i)%dS + (Ti(k)-alp%data(i)%TT)*(alp%data(i)%dS-alp%data(i-1)%dS)/&
                       (alp%data(i)%TT-alp%data(i-1)%TT)
                  alp%g(k) = alp%g(k) + div*(Si(k)-alp%data(i)%SS)
                  EXIT
               ELSE IF (alp%data(13)%TT > Ti(k)) THEN
                  alp%g(k) = alp%data(13)%value
                  alp%g(k) = alp%g(k) + alp%data(13)%dS*(Si(k)-alp%data(13)%SS)
                  EXIT
               END IF
            END DO
         END IF
      END DO

      alp%g(0) = alp%g(1)      
      alp%i(0) = alp%g(1)

      DO k = 1, d%M
         alp%i(k) = alp%g(k) + (alp%g(k+1)-alp%g(k))*(d%d_i(k)-d%d_g(k))/&
                    d%g_space(k)
      END DO


END SUBROUTINE alpha_sub      
