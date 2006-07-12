SUBROUTINE define_adv(mm,P_q,P_f,fraction,A_q,W_q,A_f,W_f,h)

  USE mean_param

  IMPLICIT NONE

  TYPE(gr_d),INTENT(IN)::mm  !grid
  DOUBLE PRECISION, DIMENSION(mm%M),INTENT(OUT)::P_q,P_f
  DOUBLE PRECISION,INTENT(IN)::fraction  !P_q_fraction(month)
  DOUBLE PRECISION, INTENT(IN)::A_q,W_q,A_f,W_f
  TYPE(height), INTENT(IN)::h !h

  INTEGER::index
  DOUBLE PRECISION,DIMENSION(mm%M)::Heavyside

  Heavyside(1:h%g) = 1.
  Heavyside(h%g+1:mm%M) = 0.
  P_q = 0.
  P_f = 0.

print*,'def adv pq',fraction

  DO index = 1, mm%M
     IF (mm%d_g(index) < 100.) THEN
        P_q(index) = fraction*12.*A_q/h%new*Heavyside(index)
        P_f(index) = fraction*12.*A_f/h%new*Heavyside(index)
     ELSE
        P_q(index) = fraction*12.*A_q/h%new*Heavyside(index)+W_q/100.
        P_f(index) = fraction*12.*A_f/h%new*Heavyside(index)+W_f/100.
     END IF
  END DO

 ! PRINT "(A)","P_q"
 ! PRINT *,P_q
 ! PRINT "(A)","P_f"
 ! PRINT *,P_f

END SUBROUTINE define_adv
