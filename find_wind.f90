! $Id$
! $Source$

module find_wind_mod
  ! Module contains:
  !  subroutine find_wind(...)

  use mean_param

  implicit none

contains

  subroutine find_wind(year, day, time, ecmapp, wind_n, wind, unow, vnow)
    ! Find the unow and vnow components of the wind at the specified 
    ! year/day/time by interpolation in the wind data structure.

    ! Arguments:
    !  Inputs:
    integer, intent(in) :: year           ! year for which wind is req'd
    integer, intent(in) :: day            ! year-day for which wind is req'd
    double precision, intent(in) :: time  ! time for which wind is req'd [s]
    integer, intent(in out) :: ecmapp     ! wind struct index of year's day 1 
    integer, intent(in) :: wind_n         ! # of lines in the wind struct
    type(wind_ecmwf), dimension(wind_n), intent(in) :: wind ! wind data struct
    !  Outputs:
    double precision, intent(out) :: unow ! u component of wind
    double precision, intent(out) :: vnow ! v component of wind

    ! Local variables:
    integer :: i                  ! wind data structure index
    integer :: hr                 ! integer part of hour in which wind is req'd
    double precision :: hr_frac   ! fraction part of hr in which wind is req'd
    double precision :: time_part ! temp storage for time calculations
    logical :: done               ! flag for search & interpolation loop

    ! Check that the 1st 2 elements of the wind struct are 1 hr apart
    ! This confirms that wind data structure probably contains hourly data
    ! *** This really only needs to be done during code initialization
    ! *** (maybe when wind data is read) rather than here at every time step
    if (wind(2)%time - wind(1)%time /= 1.) then
       ! *** Change this to write to stderr
       write (*,*) '!!! Error: find_wind: wind_dt != 1'
       stop
    endif

    ! Extract the hour from the time
    ! Get rid of whole years
    time_part = time - int(time / (365 * 86400)) * 365 * 86400
    ! Get rid of whole days (24 hr = 86400 s)
    time_part = time_part - int(time_part / 86400) * 86400
    ! Separate the hour into integer and fractional parts
    hr = int(time_part / 3600)
    hr_frac = (time_part - hr * 3600) / 3600

    ! Set index to first day of correct year in wind data structure.
    ! This is only executed if a year boundary have been crossed in
    ! the time step loop.
    ! *** This is not robust and needs to be fixed
    i = ecmapp
    if (day == 1 .AND. year == 2002) then
       ecmapp = 8761
       i = ecmapp
    else if (day == 1 .AND. year == 2003) then
       ecmapp = 17521
       i = ecmapp
    else if (day == 1 .AND. year == 2004) then
       ecmapp = 26288
       i = ecmapp
    else if (day == 1 .AND. year == 2005) then
       ecmapp = 35065
       i = ecmapp
    else if (day == 1 .AND. year == 2006) then
       ecmapp = 43825
       i = ecmapp
    endif

    ! Search and interpolation loop
    done = .false.
    do while (.not. done)
       ! Find the matching day
       if (wind(i)%Jday > day) then
          i = i - 1
       else if (wind(i)%Jday < day) then
          i = i + 1 
       else
          ! Find the matching hour
          if (wind(i)%time > hr) then
             i = i - 1
          else if (wind(i)%time < hr) then
             i = i + 1
          else
             ! Interpolate for the fractional part of the hour
             vnow = wind(i)%meridional * (1 - hr_frac) &
                  + wind(i+1)%meridional * hr_frac
             unow = wind(i)%zonal * (1 - hr_frac) &
                  + wind(i + 1)%zonal * hr_frac
             done = .true.
             ecmapp = i
          end if
       end if

       if (i < 0 .OR. i > wind_n) then
          ! *** Change this to write to stderr
          write (*,*) '!!! Error: find_wind: Search failed', &
               '\n i = ', i,          &
               '\n year = ', year,    &
               '\n day = ', day,      &
               '\n time = ', time,    &
               '\n wind_n = ', wind_n
          stop
       endif
    end do

  end subroutine find_wind

end module find_wind_mod
