SUBROUTINE interpolate(Jday, Q96, Qtot)  !interpolates to the day

  USE mean_param

  IMPLICIT NONE

  INTEGER, INTENT(IN)::Jday
  TYPE(Large1996_data), DIMENSION(14), INTENT(IN)::Q96
  DOUBLE PRECISION, INTENT(OUT)::Qtot

  INTEGER::j, end

  end = 0

  DO j = 1, 14
     IF (Q96(j)%day >= DBLE(Jday)) THEN
        end = j
        EXIT
     END IF
  END DO
  
  IF (end == 0) THEN
     PRINT "(A)","Problem in interpolate.f90: end = 0."
     STOP
  ELSE IF (end == 1) THEN
     Qtot = Q96(end)%data
  ELSE 
     Qtot = Q96(end)%data - (Q96(end)%day - DBLE(Jday))*(Q96(end)%data-Q96(end-1)%data)/&
          (Q96(end)%day-Q96(end-1)%day)
  END IF

END SUBROUTINE interpolate
