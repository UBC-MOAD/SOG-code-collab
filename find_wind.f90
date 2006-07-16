! $Id$
! $Source

subroutine find_wind(year, day, time, ecmapp, wind_n, unow, vnow, wind)
! Find the unow and vnow components of the wind at the specified 
! year/day/time by interpolation in the wind data structure.


  use mean_param  ! *** Changed from declarations to resolve conflicts
  use surface_forcing

  implicit none

  ! Arguments:
  integer, intent(in) :: year, day, wind_n
  double precision, intent(in) :: time
  integer, intent(in out) :: ecmapp      
  double precision, intent(out) :: unow, vnow
  type(wind_ecmwf), dimension(wind_n), intent(in) :: wind

  ! Local variables:
  integer :: count, deltat, requiredtime, done
  double precision :: withinhour, otime

  count = ecmapp
  deltat = wind(2)%time - wind(1)%time
  if (deltat /= 1) then
     write (*,*) '!!! Error: find_wind: deltat != 1'
     stop
  endif

  otime = time - int(time / (365 * 86400)) *365 * 86400
  otime = otime - int(otime / 86400) * 86400
  requiredtime = int(otime / 3600)
  withinhour = (otime - requiredtime * 3600) / (deltat * 3600)

  ! Set counter to frist day of year in wind data
  ! *** This is not robust and needs to be fixed
  if (day == 1 .AND. year == 2002) then
     ecmapp = 8761
     count = ecmapp
  else if (day == 1 .AND. year == 2003) then
     ecmapp = 17521
     count = ecmapp
  else if (day == 1 .AND. year == 2004) then
     ecmapp = 26288
     count = ecmapp
  else if (day == 1 .AND. year == 2005) then
     ecmapp = 35065
     count = ecmapp
  else if (day == 1 .AND. year == 2006) then
     ecmapp = 43825
     count = ecmapp
  endif

  ! Search and interpolation loop
  done = 0
  do while (done == 0)
     if (wind(count)%Jday > day) then
        count = count - 1
     else if (wind(count)%Jday < day) then
        count = count + 1 
     else
        if (wind(count)%time > requiredtime) then
           count = count - 1
           write (*,*) '!!! Warning: find_wind: ', &
                'time beyond end of wind data', count
        else if (wind(count)%time < requiredtime) then
           count = count + 1
           write (*,*) '!!! Warning: find_wind: ', &
                'time before start of wind data', count, day, wind(count)%Jday
        else
           ! Interpolation
           vnow = wind(count)%meridional * (1 - withinhour) &
                + wind(count+1)%meridional * withinhour
           unow = wind(count)%zonal * (1 - withinhour) &
                + wind(count + 1)%zonal * withinhour
           done = 1
           ecmapp = count
           write (*,*) 'find_wind: time', count
        end if
     end if

     if (count < 0 .OR. count > wind_n) then
        write (*,*) count, '!!! Error: find_wind: Interpolation failed'
        write (*,*) count, year, day, time, wind(count)%year, &
             wind(count)%Jday, wind(count)%time
        stop
     endif
  end do

end subroutine find_wind
