SUBROUTINE y_jday_t(L_data)

  USE mean_param
  USE surface_forcing

  IMPLICIT NONE

  EXTERNAL Julian_day
  TYPE(papmd), INTENT(IN OUT)::L_data  !Large_data(xx)

  INTEGER::mdh,m,dh,d,h, k, marker

  marker = 0
  !time is in Greenwich mean time  hour + long
  L_data%y = L_data%ymdh/1000000
  mdh = L_data%ymdh-L_data%y*1000000
  m = mdh/10000            !month
  dh = mdh-m*10000
  d = dh/100               !day
  h = dh-100*d             !hour
 ! L_data%t = DBLE(h*3600)  !day_time (s)

  L_data%y = 1900+L_data%y  !year
  L_data%leap = 0           !is leap year? NO
  DO k = 1,size(leap_year)
     IF (L_data%y == leap_year(k)) THEN
        L_data%leap = 1     !is leap year? YES
        EXIT
     END IF
  END DO

  L_data%t =DBLE(h)-Lon/15.  !Lon in degrees, h in hours
  IF (L_data%t <= 0.) THEN
     L_data%t = 3600*(24. + L_data%t)
     marker = 1
  ELSE
     L_data%t = 3600*L_data%t
  END IF

  CALL Julian_day(L_data%leap,m,d,L_data%jday,1)  !day and month ==> Julian day

  IF (marker == 1) THEN
     L_data%jday = L_data%jday-1
     IF (L_data%jday == 0) THEN
        L_data%y = L_data%y-1
        L_data%leap = 0
        DO k = 1,size(leap_year)
           IF (L_data%y == leap_year(k)) THEN
              L_data%leap = 1
              EXIT
           END IF
        END DO
        L_data%jday = 365 + L_data%leap
     END IF
  END IF

  !!Change units

  L_data%Uten = L_data%Uten/100. !m/s
  L_data%SST = L_data%SST + 273.15  !K
  L_data%Ta = L_data%Ta + 273.15    !K
  L_data%Td = L_data%Td + 273.15    !K
  
END SUBROUTINE y_jday_t
