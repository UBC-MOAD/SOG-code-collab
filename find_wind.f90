SUBROUTINE find_wind(year,day,time,ecmapp,unow,vnow,wind)

      USE declarations
      USE surface_forcing

      IMPLICIT NONE

      INTEGER, INTENT(IN)::year,day
      DOUBLE PRECISION, INTENT (IN):: time
      INTEGER, INTENT(IN OUT)::ecmapp      
      TYPE(wind_ecmwf), DIMENSION(wind_n), INTENT(IN)::wind
      DOUBLE PRECISION, INTENT(OUT):: unow,vnow


      INTEGER :: deltat,requiredtime,done,wind_end
      DOUBLE PRECISION :: withinhour,otime

      count = ecmapp
      deltat = wind(2)%time-wind(1)%time


      if (deltat <> 1) then
           write (*,*) 'find wind oh no'
           stop
        endif

!      wind_end= 30802

      otime = time - int(time/(365*86400))*365*86400
      otime = otime - int(otime/86400)*86400
      requiredtime = int(otime/3600)
      withinhour = (otime-requiredtime*3600)/(deltat*3600)
     
      done = 0

!print*,day,count,ecmapp,'fw1'

if (day==1 .AND. year==2002) then
ecmapp=8761
count=ecmapp
else if (day==1 .AND. year==2003) then
ecmapp=17521
count=ecmapp
else if (day==1 .AND. year==2004) then
ecmapp=26281
count=ecmapp
else if (day==1 .AND. year==2005) then
ecmapp=35065
count=ecmapp
else if (day==1 .AND. year==2006) then
ecmapp=43825
count=ecmapp
endif

      do while (done == 0)
         if (wind(count)%Jday > day) then
            count = count -1
         else if (wind(count)%Jday < day) then
            count = count + 1 
        else
            if (wind(count)%time > requiredtime) then
               count = count - 1
               write (*,*) 'time too far',count
            else if (wind(count)%time < requiredtime) then
               count = count + 1
               write (*,*) 'time too small',count,day,wind(count)%Jday
            else
               vnow = wind(count)%meridional*(1-withinhour)+wind(count+1)%meridional*withinhour
               unow = wind(count)%zonal*(1-withinhour)+wind(count+1)%zonal*withinhour
!print*,'unow,vnow,wind(count)%meridional,wind(count)%zonal, otime,withinhour,requiredtime,time'
!print*,unow,vnow,wind(count)%meridional,wind(count)%zonal,otime,withinhour,requiredtime,time
!pause

              done = 1
               ecmapp = count
               write (*,*) 'time', count

            end if
         end if
         if (count < 0 .OR. count > wind_n) then

            write (*,*) count, ' find wind'
            write (*,*) count,year,day,time,wind(count)%year,wind(count)%Jday,wind(count)%time
            stop
         endif
      end do

    END SUBROUTINE find_wind



