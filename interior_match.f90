! $Id$
! $Source$

SUBROUTINE interior_match(mm, hh, X, X_min)

      USE mean_param
      USE surface_forcing

      IMPLICIT NONE

      TYPE(gr_d), INTENT(IN)::mm !
      TYPE(diff), INTENT(IN OUT)::X
      TYPE(height), INTENT(IN)::hh
      DOUBLE PRECISION, INTENT(IN)::X_min

      DOUBLE PRECISION::R

      !!!!!Case (1) h%g-1 == h%i and n = h%g - 1 or h%i
      !!!!!Case (2) h%g == h%i and n = h%g or h%i

         R = (hh%new - mm%d_i(hh%i-1))/mm%i_space(hh%i) ! R below (D5b)
         X%div = (1.0 - R)*((X%total(hh%i-1)-X%total(hh%i))/mm%i_space(hh%i))&
         + R*((X%total(hh%i)-X%total(hh%i+1))/mm%i_space(hh%i+1)) ! (D5b)

         X%h = X%total(hh%i) + X%div*(mm%d_i(hh%i)-hh%new) ! (D5a)


! note neither of these should happen
         IF (X%h < 0.) THEN
            X%h = 0.
            X%div = 0.
         ELSE IF (X%div > 0.) THEN
            X%div = 0.
         END IF

END SUBROUTINE interior_match

