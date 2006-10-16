! $Id$
! $Source$

SUBROUTINE initialize

      USE declarations

      IMPLICIT NONE

      time = dble(t_o)
      dummy_time = 0.
      day = day_o   !initialize julian day
     ! month = month_o   !KC june 28, 2004
       year = year_o

      day_time = time
      time_step = 1

      j_gamma = 0  !h_gamma = 0 or wt_r contribution to gamma%t vanishes

      neg_count = 0  ! counts number of times scheme biology becomes negative

      f_ratio = 0.
      null_vector = 0.
end subroutine initialize
