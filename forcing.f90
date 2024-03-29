module forcing
  ! Calculate forcing values and allow them to be shifted in time, or
  ! scale them, or set them constant.

  ! Data is coming from historic data files. We have ctd data for some
  ! 70s, 80s, 90s.  Metar data goes back to 50s and river data goes
  ! back to 1910s.

  ! Read in NY years of data at a time. Wind data is in a single array
  ! with 1 hour per line, cf, atemp and humid are in array of day and
  ! time, and the river flows are one value per day.

  use precision_defs, only: dp, sp

  implicit none

  private
  public :: &
       ! Types:
       vary_quantity, vary_forcing, &
       ! Variables:
       UseRiverTemp, vary, &
       Qinter, Einter, RiverTemp, &
       cf_value, atemp_value, humid_value, &
       unow, vnow, &
       ! Subroutines:
       init_forcing,   &
       read_variation, &
       read_forcing,   &
       get_forcing,    &
       dalloc_forcing_variables

  ! Public type definitions:
  !
  ! Facilitate time shifting, scaling, and setting to constant value
  ! the forcing quantities
  type :: vary_quantity
     logical :: enabled ! is the parameter (wind, cf, rivers, T) to be varied
     real(kind=dp) :: shift ! use wind/rivers shifted in time (in years)
     real(kind=dp) :: fraction ! multiply wind/rivers strength by this amount
     logical :: fixed ! use a fixed value
     ! (if this is set to true, shift and fraction are not used and value is)
     real(kind=dp) :: value ! use this value if fixed is set true
     real(kind=dp) :: addition ! add this value (after fraction is multiplied)
  end type vary_quantity

  type :: vary_forcing
     type(vary_quantity) :: wind, cf, rivers, temperature
  end type vary_forcing

  ! Public Variable Declarations:
  type(vary_forcing) :: vary
  ! Do we add the cooling/warming effect of the Major River?
  logical :: UseRiverTemp
  real(kind=dp) :: Qinter, Einter, &! Current flow of Major and Minor Rivers
                   RiverTemp  ! Current temperature of Major River
  ! current values of cloud fraction, air temperature and humidity
  real(kind=sp) :: cf_value, atemp_value, humid_value
  real(kind=dp) :: unow, vnow ! Current wind components

  ! Private module variable declarations:
  integer ::     &
       NY,       &  ! Number of years of forcing data
       wind_n,   &  ! Size of hourly wind speed forcing data arrays
       met_n,    &  ! Size of daily meteo forcing data arrays
       river_n,  &  ! Size of daily river flow forcing data arrays
       startyear    ! Starting year (startyear + NY-1 more years of data)
  real(kind=dp), dimension(:), allocatable :: &
       wind_eastnorth,  &  ! North-east wind speed forcing data array
       wind_northwest,  &  ! North-west wind speed forcing data array
       Qriver,          &  ! Major river flow forcing data array
       Eriver              ! Minor river flow forcing data array

  real(kind=sp), dimension(:,:), allocatable :: &
       cf,              &  ! Cloud fraction forcing data array
       atemp,           &  ! Air temperature forcing data array
       humid               ! Relative humidity forcing data array

contains

  subroutine init_forcing
    ! Allocate memory for forcing data arrays,
    ! and read forcing data values.
    use input_processor, only: getpari
    implicit none


    NY = getpari('years of forcing data')
    wind_n = NY * 366 * 24
    met_n =  NY * 366
    river_n = NY * 366
    call alloc_forcing_variables(wind_n, met_n, river_n)
    call read_forcing
  end subroutine init_forcing


  subroutine read_variation
    ! Read variation parameters from the infile.
    use io_unit_defs, only: stderr
    implicit none

    ! Read wind variation parameters
    call read_vary ("vary%wind", vary%wind)
    if (vary%wind%enabled) then
       if  (vary%wind%fixed) then
          write (stderr, *) 'Fixed wind magnitude has not been implemented'
          call exit(1)
       else
          if (vary%wind%shift /= dint(vary%wind%shift)) then
             write (stderr, *) 'Non-full year shifts have not been ', &
                  'implemented for the wind'
             call exit(1)
          endif
          if (vary%wind%addition /= 0) then
             write (stderr, *) "addtion not implemented for wind"
             call exit(1)
          endif
       endif
    endif
    ! Read cloud fraction variation parameters
    call read_vary ("vary%cf", vary%cf)
    if (vary%cf%enabled .and. .not. vary%cf%fixed) then
       if (vary%cf%shift /= int(vary%cf%shift)) then
          write (stderr, *) "non-integer cf shifts not implemented"
          call exit(1)
       endif
       if (vary%cf%addition /= 0) then
          write (stderr, *) "addition not impleemented for cf"
          call exit(1)
       endif
    endif
    ! Read river flow variation parameters
    call read_vary ("vary%rivers",vary%rivers)
    if (vary%rivers%enabled .and. .not.vary%rivers%fixed) then
       if (vary%rivers%addition /= 0 ) then
          write (stderr, *) "addition not implemented for rivers"
          call exit(1)
       endif
    endif
    ! Read air temperature variation parameters
    call read_vary ("vary%temperature", vary%temperature)
    if (vary%temperature%enabled) then
       if (vary%temperature%fixed) then
          write (stderr, *) "Fixed temperature not enabled"
          call exit (1)
       else
          if (vary%temperature%fraction /= 1) then
             write (stderr, *) "Multiplication of temperature not enabled, ", &
                  "use addition"
             call exit(1)
          endif
       endif
    endif
  end subroutine read_variation


  subroutine read_forcing
    ! Read forcing data and apply variations (time shift, scale,
    ! etc.), if requested.
    use input_processor, only: getpari, getpars, getparl
    use io_unit_defs, only: forcing_data, stdout, stderr
    use unit_conversions, only: CtoK
    use numerics, only: &
         initDatetime   ! Date/time of initial conditions
    implicit none
    ! Local variables:
    !

    ! Should average data be used
    character*80 :: use_average_forcing_data
    ! File name for minor and major river data files
    character*80 :: minor_river_file, major_river_file, river_nut_file
    ! Indices
    integer :: ic, jc, j
    ! Data file values
    integer :: para, stn, year, day, month
    real(kind=sp) :: hour
    ! Parameters for integration of minor river flow
    integer :: integ_days  ! number of days to integrate over
    real(kind=dp) :: integ_minor_river(river_n)  ! integrated flow
    ! Years and days for data shifting
    integer :: &
         wind_startyear,    &
         temp_startyear,    &
         cf_startyear,      &
         rivers_startyear,  &
         rivers_startmonth, &
         rivers_startday
    logical ::  found_data


    ! Start year for forcing data from initialization date/time for run
    startyear = initDatetime%yr

    ! Wind
    ! shift the start year if requested
    if (vary%wind%enabled) then
       wind_startyear = startyear + int(vary%wind%shift)
    else
       wind_startyear = startyear
    endif

    ! check if using average data
    use_average_forcing_data = getpars("use average/hist forcing")
    flush(stdout)
    if (use_average_forcing_data .eq. "yes"  &
        .or. use_average_forcing_data .eq. "fill") then
       open(unit=forcing_data, file=getpars("average/hist wind"))
       flush(stdout)
       do jc = 1, wind_n/NY
          read(forcing_data,*) day, month, year, hour, &
                  wind_eastnorth(jc), wind_northwest(jc)
          wind_eastnorth(jc) = 0.94* wind_eastnorth(jc)  ! tune the wind down
          wind_northwest(jc) = 0.94* wind_northwest(jc)
       enddo
       do jc = wind_n/NY+1, wind_n
          wind_eastnorth(jc) = wind_eastnorth(jc-wind_n/NY)
          wind_northwest(jc) = wind_northwest(jc-wind_n/NY)
       enddo
       close(forcing_data)
    elseif (use_average_forcing_data .eq. "histfill") then
       open(unit=forcing_data, file=getpars("average/hist wind"))
       flush(stdout)
       do jc = 1, wind_n
          read(forcing_data,*) day, month, year, hour, &
                  wind_eastnorth(jc), wind_northwest(jc)
       enddo
       close(forcing_data)
    endif

    ! standard data
    if (use_average_forcing_data .eq. "fill"          &
       .or. use_average_forcing_data .eq. "histfill"  &
       .or. use_average_forcing_data .eq. "no") then
       open(unit=forcing_data, file=getpars("wind"))
       flush(stdout)
       found_data = .false.
       do while (.not. found_data)
          read(forcing_data, *, end=777) day, month, year, hour, &
               wind_eastnorth(1), wind_northwest(1)
          ! Assign hour and month to themselves to prevent g95 from
          ! throwing "set but never used" warnings
          hour = hour
          month = month
          if (year == wind_startyear) then
             do jc = 2, wind_n
                read(forcing_data,*,end=777) day, month, year, hour, &
                     wind_eastnorth(jc), wind_northwest(jc)
             enddo
             found_data = .true.
          endif
       enddo
777    if (use_average_forcing_data .eq. "no" .and. .not. found_data) then
          write (stderr, *) "End of File on Wind Data"
          call exit(1)
       endif
       close(forcing_data)
    endif

    ! Air Temperature
    ! shift the start year if requested
    if (vary%temperature%enabled) then
       temp_startyear = startyear + int(vary%temperature%shift)
    else
       temp_startyear = startyear
    endif

    ! Average data
    if (use_average_forcing_data .eq. "yes" &
        .or. use_average_forcing_data .eq. "fill") then
       open(unit=forcing_data, file=getpars("average/hist air temp"))
       flush(stdout)
       do jc = 1, met_n/NY
          read(forcing_data,*) stn, year, month, day, para, &
               (atemp(jc,j), j=1,24)
          do j= 1, 24
             atemp(jc,j) = CtoK(atemp(jc,j) / 10.)
          enddo
       enddo
       do jc = met_n/NY+1, met_n
          atemp(jc,:) = atemp(jc-met_n/NY,:)
       enddo
       close(forcing_data)
    elseif (use_average_forcing_data .eq. "histfill") then
       open(unit=forcing_data, file=getpars("average/hist air temp"))
       flush(stdout)
       do jc = 1, met_n
          read(forcing_data,*) stn, year, month, day, para, &
               (atemp(jc,j), j=1,24)
          do j= 1, 24
             atemp(jc,j) = CtoK(atemp(jc,j) / 10.)
          enddo
       enddo
       close(forcing_data)
    endif

    ! Standard data
    if (use_average_forcing_data .eq. "fill" &
        .or. use_average_forcing_data .eq. "histfill" &
        .or. use_average_forcing_data .eq. "no") then
       open(unit=forcing_data, file=getpars("air temp"))
       flush(stdout)
       found_data = .false.
       do while (.not. found_data)
          read(forcing_data, *, end=780) stn, year, month, day, para, &
               (atemp(1,j), j=1,24)
          do j= 1, 24
             atemp(1,j) = CtoK(atemp(1,j) / 10.)
          enddo
          if (year == temp_startyear) then
             do jc = 2, met_n
                read (forcing_data, *, end=780) stn, year, month, day, para, &
                  (atemp(jc,j),j=1,24)
                do j= 1, 24
                   atemp(jc,j) = CtoK(atemp(jc,j) / 10.)
                enddo
             enddo
             found_data = .true.
          endif
       enddo
780    if (use_average_forcing_data .eq. "no" .and. .not. found_data) then
          write (stderr, *) "End of File on Air Temp Data"
          call exit(1)
       endif
       close(forcing_data)
    endif

    ! Cloud Fraction
    ! shift the start year if requested
    if (vary%cf%enabled .and. .not.(vary%cf%fixed) ) then
       cf_startyear = startyear + int(vary%cf%shift)
    else
       cf_startyear = startyear
    endif

    ! Average data
    if (use_average_forcing_data .eq. "yes" &
        .or. use_average_forcing_data .eq. "fill") then
       open(unit=forcing_data, file=getpars("average/hist cloud"))
       flush(stdout)
       do jc = 1, met_n/NY
          read(forcing_data,*) stn, year, month, day, para, (cf(jc,j), j=1,24)
          ! later on basically use integers by taking floor
          cf(jc,:) = floor(cf(jc,:) + 0.5)
       enddo
       do jc = met_n/NY+1, met_n
          cf(jc,:) = cf(jc-met_n/NY,:)
       enddo
       close(forcing_data)
    elseif  (use_average_forcing_data .eq. "histfill") then
       open(unit=forcing_data, file=getpars("average/hist cloud"))
       flush(stdout)
       do jc = 1, met_n
          read(forcing_data,*) stn, year, month, day, para, (cf(jc,j), j=1,24)
       enddo
       close(forcing_data)
    endif

    ! standard data
    if (use_average_forcing_data .eq. "fill" &
         .or. use_average_forcing_data .eq. "histfill" &
         .or. use_average_forcing_data .eq. "no") then
       open(unit=forcing_data, file=getpars("cloud"))
       flush(stdout)
       found_data = .false.
       do while (.not. found_data)
          read (forcing_data, *, end=779) stn, year, month, day, para, &
               (cf(1,j), j=1,24)
          ! Assign stn and para to themselves to prevent g95 from
          ! throwing "set but never used" warnings
          stn = stn
          para = para
          if (year == cf_startyear) then
             do jc = 2, met_n
                read(forcing_data, *, end=779) stn, year, month, day, para, &
                     (cf(jc,j), j=1,24)
             enddo
             found_data = .true.
          endif
       enddo
779    if (use_average_forcing_data .eq. "no" .and. .not. found_data) then
          write (stderr, *) "End of File on Cloud Data"
          call exit(1)
       endif
       close (forcing_data)
    endif

    ! Humidity
    ! Average data
    if (use_average_forcing_data .eq. "yes" &
        .or. use_average_forcing_data .eq. "fill") then
       open(unit=forcing_data, file=getpars("average/hist humidity"))
       flush(stdout)
       do jc = 1, met_n/NY
          read(forcing_data,*) stn, year, month, day, para, &
               (humid(jc,j), j=1,24)
       enddo
       do jc = met_n/NY+1, met_n
          humid(jc,:) = humid(jc-met_n/NY,:)
       enddo
       close(forcing_data)
    elseif (use_average_forcing_data .eq. "histfill") then
       open(unit=forcing_data, file=getpars("average/hist humidity"))
       flush(stdout)
       do jc = 1, met_n
          read(forcing_data,*) stn, year, month, day, para, &
               (humid(jc,j), j=1,24)
       enddo
       close(forcing_data)
    endif

    ! Standard data
    if (use_average_forcing_data .eq. "fill" &
        .or. use_average_forcing_data .eq. "histfill" &
        .or. use_average_forcing_data .eq. "no") then
       open(unit=forcing_data, file=getpars("humidity"))
       flush(stdout)
       found_data = .false.
       do while (.not. found_data)
          read (forcing_data, *, end=781) stn, year, month, day, para, &
               (humid(1,j), j=1,24)
          if (year == startyear) then
             do jc = 2, met_n
                read(forcing_data, *, end=781)stn, year, month, day, para, &
                     (humid(jc,j), j=1,24)
             enddo
             found_data = .true.
          endif
       enddo
781    if (use_average_forcing_data .eq. "no" .and. .not. found_data) then
          write (stderr, *) "End of File on Humidity Data"
          call exit(1)
       endif
       close(forcing_data)
    endif

    ! Major river
    ! shift the start year if requested
    if (vary%rivers%enabled .and. .not. vary%rivers%fixed) then
       call shift_time (1, 1, startyear, vary%rivers%shift, &
            rivers_startday, rivers_startmonth, rivers_startyear)
    else
       rivers_startday = 1
       rivers_startmonth = 1
       rivers_startyear = startyear
    endif

    ! check if using average data
    if (use_average_forcing_data .eq. "yes" &
        .or. use_average_forcing_data .eq. "fill") then
       major_river_file = getpars("average/hist major river")
       flush(stdout)
       open(unit=forcing_data, file=major_river_file)
       do jc = 1, rivers_startday-1
          read(forcing_data, *) year, month, day, Qriver(jc)
       enddo
       do jc = rivers_startday, river_n/NY
          read(forcing_data,*) year, month, day, Qriver(jc+1-rivers_startday)
       enddo
       close(forcing_data)
       if (rivers_startday .ne. 1) then
          open(unit=forcing_data, file=major_river_file)
          do jc = 1, rivers_startday-1
             read(forcing_data,*) year, month, day, &
                  Qriver(river_n/NY+1-rivers_startday+jc)
          enddo
          close(forcing_data)
       endif
       do jc = river_n/NY+1, river_n
          Qriver(jc) = Qriver(jc-river_n/NY)
       enddo
    elseif  (use_average_forcing_data .eq. "histfill") then
       major_river_file = getpars("average/hist major river")
       flush(stdout)
       open(unit=forcing_data, file=major_river_file)
       do jc = 1, rivers_startday-1
          read(forcing_data, *) year, month, day, Qriver(jc)
       enddo
       do jc = rivers_startday, river_n
          read(forcing_data,*) year, month, day, Qriver(jc+1-rivers_startday)
       enddo
       close(forcing_data)
       if (rivers_startday .ne. 1) then
          open(unit=forcing_data, file=major_river_file)
          do jc = 1, rivers_startday-1
             read(forcing_data,*) year, month, day, &
                  Qriver(river_n+1-rivers_startday+jc)
          enddo
          close(forcing_data)
       endif
    endif

    ! standard data
    if (use_average_forcing_data .eq. "fill" &
        .or. use_average_forcing_data .eq. "histfill" &
        .or. use_average_forcing_data .eq. "no") then
       open(forcing_data, file=getpars("major river"))
       flush(stdout)
       found_data = .false.
       do while(.not. found_data)
          read(forcing_data, *, end=778) year, month, day, Qriver(1)
          if(year == rivers_startyear .and. month == rivers_startmonth &
               .and. day == rivers_startday) then
             do jc = 2, river_n
                read(forcing_data,*,end=778) year, month, day, Qriver(jc)
             enddo
             found_data = .true.
          endif
       enddo
778    if (use_average_forcing_data .eq. "no" .and. .not. found_data) then
          write (stderr, *) "End of File on Major River Data"
          stop
       endif
       close(forcing_data)
    endif
    ! Include cooling/heating effect of major river?
    UseRiverTemp = getparl("use river temp")

    ! Minor river
    ! check if using average data
    if (use_average_forcing_data .eq. "yes" &
        .or. use_average_forcing_data .eq. "fill") then
       minor_river_file = getpars("average/hist minor river")
       flush(stdout)
       if (minor_river_file /= "N/A") then
          open(unit=forcing_data, file=minor_river_file)
          ! Use same shift as major river
          do jc = 1, rivers_startday-1
             read(forcing_data, *) year, month, day, Eriver(jc)
          enddo
          do jc = rivers_startday, river_n/NY
             read(forcing_data,*) year, month, day, Eriver(jc+1-rivers_startday)
          enddo
          close(forcing_data)
          if (rivers_startday .ne. 1) then
             open(unit=forcing_data, file=minor_river_file)
             do jc = 1, rivers_startday-1
                read(forcing_data,*) year, month, day, &
                     Eriver(river_n/NY+1-rivers_startday+jc)
             enddo
             close(forcing_data)
          endif
          do jc = river_n/NY+1, river_n
             Eriver(jc) = Eriver(jc-river_n/NY)
          enddo
       else ! no average minor river
          do jc = 1, river_n
             Eriver(jc) = 0.0d0
          enddo
       endif
    elseif  (use_average_forcing_data .eq. "histfill") then
       minor_river_file = getpars("average/hist minor river")
       flush(stdout)
       if (minor_river_file /= "N/A") then
          open(unit=forcing_data, file=minor_river_file)
          ! Use same shift as major river
          do jc = 1, rivers_startday-1
             read(forcing_data, *) year, month, day, Eriver(jc)
          enddo
          do jc = rivers_startday, river_n
             read(forcing_data,*) year, month, day, Eriver(jc+1-rivers_startday)
          enddo
          close(forcing_data)
          if (rivers_startday .ne. 1) then
             open(unit=forcing_data, file=minor_river_file)
             do jc = 1, rivers_startday-1
                read(forcing_data,*) year, month, day, &
                     Eriver(river_n+1-rivers_startday+jc)
             enddo
             close(forcing_data)
          endif
       else ! no average minor river
          do jc = 1, river_n
             Eriver(jc) = 0.0d0
          enddo
       endif
    endif

    ! standard data
    if (use_average_forcing_data .eq. "fill" &
        .or. use_average_forcing_data .eq. "histfill" &
        .or. use_average_forcing_data .eq. "no") then
       minor_river_file = getpars("minor river")
       flush(stdout)
       if(minor_river_file /= "N/A") then
          open(forcing_data, file=minor_river_file)
          found_data = .false.
          do while (.not. found_data)
             read (forcing_data, *, end=175) year, month, day, Eriver(1)
             ! Use same shift as major river
             if (year == rivers_startyear .and. month == rivers_startmonth &
                  .and. day == rivers_startday) then
                do jc = 2, river_n
                   read (12,*, end=175) year, month, day, Eriver(jc)
                enddo
                found_data = .true.
             endif
          enddo
          ! Read minor river flow data from an alternative file. One
          ! example use for this is historical Strait of Georgia runs for
          ! which there is no Englishman River data available, so the
          ! Nanaimo River is used instead.
175       if(.not. found_data) then
             close(forcing_data)
             minor_river_file = getpars("alt minor river")
             flush(stdout)
             open(forcing_data, file=minor_river_file)
             do while (.not. found_data)
                read (forcing_data, *, end=782) year, month, day, Eriver(1)
                if (year == rivers_startyear .and. month == rivers_startmonth &
                     .and. day == rivers_startday) then
                   do jc = 2, river_n
                      read (12,*) year, month, day, Eriver(jc)
                   enddo
                   found_data = .true.
                endif
             enddo
782          if (use_average_forcing_data .eq. "no" &
                 .and. .not. found_data) then
                write (stderr, *) "End of File on Minor River Data and ", &
                                  "Alt Minor River Data"
                call exit(1)
             endif
           endif
          close(forcing_data)
       else
          do jc = 1, river_n
             Eriver(jc) = 0.0d0
          enddo
       endif
    endif

    ! Number of days to integrate the minor river data over.  This is
    ! used because the minor river can be a proxy for regional
    ! rainfall events. In the case of the Strait of Georgia, the
    ! Englishman River is used to represent all of the fresh water
    ! inputs other than the Fraser. Some of those rivers flow into
    ! lakes, or dammed reservoirs, so there is a lag between rainfall
    ! events and their effect at the model location. Integration of
    ! the minor river flow data accounts for these effects. The number
    ! of days over which to integrate needs to be determined by data
    ! fitting.
    if (minor_river_file /= "N/A") then
       integ_days = getpari("minor river integ days")
       ! Integrate minor river data over the specified number of
       ! days.
       do ic = 1, integ_days
          integ_minor_river(ic) = sum(Eriver(1:ic)) / dble(ic)
       enddo
       do ic = integ_days + 1, river_n
          integ_minor_river(ic) = sum(Eriver(ic - integ_days:ic)) / integ_days
       enddo
       Eriver(1:river_n) = integ_minor_river(1:river_n)
    endif

    ! read River Nutrients file (just a stub here)

    river_nut_file = getpars("river nutrients file")
end subroutine read_forcing


  subroutine get_forcing (year, day, day_time)

    use unit_conversions, only: CtoK

    ! given the year, day and day_time find the physical forcing values
    ! taking into account any shifts or other changes set by vary

    integer, intent(in) :: year ! current year
    integer, intent(in) :: day ! current year day
    real(kind=dp), intent(in) :: day_time ! current time of day in seconds

    ! arrays we are interpolating across to get current values:
    ! Wind data is in a single array with 1 hour per line, cf, atemp and
    ! humid are in array of day and time
    ! the river flows are one value per day.

    ! local variable

    integer:: accul_day ! accumulated day number starting from startyear
    integer:: j ! index to hour, 1 am has j=2
    integer:: index_day ! index into wind data for the end of the day before
    real(kind=dp) :: hr_frac ! fraction of the hour
    real(kind=dp) :: dday ! decimal yearday


    accul_day = accum_day(year, day)

    ! RIVERS

    if (vary%rivers%enabled .and. vary%rivers%fixed) then
       Qinter = SNGL(vary%rivers%value)
       Einter = 0.
    else
       Qinter = (day_time * Qriver(accul_day) &
            + (86400 - day_time) * Qriver(accul_day - 1)) / 86400.
       Einter = (day_time * Eriver(accul_day) &
            + (86400 - day_time) * Eriver(accul_day-1)) / 86400.
       if (vary%rivers%enabled) then
          Qinter = Qinter * SNGL(vary%rivers%fraction)
          Einter = Einter * SNGL(vary%rivers%fraction)
       endif
    endif

    ! RIVER temperature (hard coded values here are for the Fraser)
    ! fit May 2010 by SEA, labbook page 65
    ! /ocean/sallen/allen/research/sog/sog-forcing/rivers/fit_river_temp.m

    dday = day+day_time/86400.

    if (UseRiverTemp) then
       if (dday < 52.8 .or. dday > 334.4) then
          RiverTemp = CtoK(2.5)
       elseif (dday < 232.9) then
          RiverTemp = CtoK(2.5 + (dday - 52.8) * (19.3 - 2.5)/(232.9 - 52.8))
       else
          RiverTemp = CtoK(19.3 + (dday - 232.9) * (2.5 - 19.3)/(334.4 - 232.9))
       endif
    endif


    ! MET DATA (cloud fraction, air temperature and humidity)

    ! index into met data array for the appropriate hr
    j = floor(day_time/3600.0) + 1

    if (vary%cf%enabled) then
       if (vary%cf%fixed) then
          cf_value = sngl(vary%cf%value)
       else
          cf_value = sngl(cf(accul_day,j)*vary%cf%fraction)
       endif
    else
       cf_value = cf(accul_day,j)
    endif

    ! TEMPERATURE

    if (vary%temperature%enabled .and. .not.vary%temperature%fixed) then
       ! multiply water temperature increment but 1.33 to get observed
       ! Masson and Cummins 2004 air temperature increment
       atemp_value = atemp(accul_day,j) + sngl(vary%temperature%addition)*1.33
    else
       atemp_value = atemp(accul_day,j)
    endif

    humid_value = humid(accul_day,j)

    ! WIND

    ! wind for the previous day, last value is at
    index_day = (accul_day-1) * 24

    ! frac hour is
    hr_frac = day_time/3600. - (j-1)
    ! index the wind at the start of the hour is
    j = index_day + j

    vnow = wind_northwest(j) * (1 - hr_frac) &
         + wind_northwest(j+1) * hr_frac
    unow = wind_eastnorth(j) * (1 - hr_frac) &
         + wind_eastnorth(j+1) * hr_frac

    if (vary%wind%enabled) then
       vnow = vnow * vary%wind%fraction
       unow = unow * vary%wind%fraction
    endif
  end subroutine get_forcing


  subroutine read_vary (quantity_string, quantity)
    ! Read forcing data variation parameter values.
    use input_processor, only: getparl, getpard
    implicit none
    ! Arguments:
    ! the leading part of the variable quantity to be read in (eg vary%wind)
    character, intent(in) :: quantity_string*(*)
    type(vary_quantity), intent(out) :: quantity

    quantity%enabled = getparl(quantity_string//"%enabled")
    if (quantity%enabled) then
       quantity%fixed = getparl(quantity_string//"%fixed")
       if(quantity%fixed) then
          quantity%value = getpard(quantity_string//"%value")
       else
          quantity%shift = getpard(quantity_string//"%shift")
          quantity%fraction = getpard(quantity_string//"%fraction")
          quantity%addition = getpard(quantity_string//"%addition")
       endif
    endif
  end subroutine read_vary


  integer function accum_day (year, day) result(day_met)
    ! Calculate the accumulated day from Jan 1 of startyear. Based on
    ! the fact that NY years of data is read in.
    use io_unit_defs, only: stderr
    implicit none
    ! Arguments:
    integer, intent(in) :: day, year ! year day and year
    ! Local Variables:
    integer :: i, & ! counter
         numY, & ! number of years beyond startyear we are at
         nleap   ! number of leap years we have modelled

    if (year == startyear) then
       day_met = day
    else
       numY = year - startyear
       if (numY+1 .gt. NY) then
          write (stderr,*) "Run is exceeding data input"
          stop
       endif
       nleap = 0
       do i=1,numY
          if (leapyear(year-i)) nleap = nleap + 1
       enddo
       day_met = day + numY*365 + nleap
    endif

  end function accum_day


  subroutine shift_time (day, month, year, shift, &
       shift_day, shift_month, shift_year)

    ! subroutine to shift by part of a year

    integer, intent(in) :: day, month, year ! current day, month, and year
    real(kind=dp), intent(in) :: shift ! shift by fraction of a year

    ! shifted day and year
    integer, intent(out) :: shift_day, shift_month, shift_year

    ! Local variables
    integer :: yearday
    integer, dimension(12) :: month_perp, month_leap, month_index

    month_perp = (/ 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 /)
    month_leap = (/ 31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 /)
    shift_day   = day + floor(365 * shift)
    shift_month = month
    shift_year  = year

    if (leapyear(shift_year)) then
       month_index = month_leap
    else
       month_index = month_perp
    endif

    if (month == 1) then
       yearday = day
    else
       yearday = sum(month_index(1: month - 1)) + day
    endif

    do while (shift_day .gt. month_index(shift_month) .or. shift_day .lt. 0)
       if ((leapyear(shift_year) .and. shift_day .gt. 366 - yearday) .or. &
            shift_day .gt. 365 - yearday) then
          if (leapyear(shift_year)) then
             shift_day = shift_day - (366 - yearday)
          else
             shift_day = shift_day - (365 - yearday)
          endif
          shift_year  = shift_year + 1
          shift_month = 1
          yearday     = 1
          if (leapyear(shift_year)) then
             month_index = month_leap
          else
             month_index = month_perp
          endif
       elseif (shift_day .lt. -1 * yearday) then
          shift_day   = shift_day + yearday
          shift_year  = shift_year - 1
          shift_month = 12
          if (leapyear(shift_year)) then
             yearday = 366
             month_index = month_leap
          else
             yearday = 365
             month_index = month_perp
          endif
       else
          if (shift_day .gt. 0) then
             shift_day = shift_day - month_index(shift_month)
             if (shift_month == 12) then
                shift_month = 1
             else
                shift_month = shift_month + 1
             endif
          else
             if (shift_day .gt. -1 * month_index(month)) then
                shift_day = shift_day + month_index(shift_month)
             else
                shift_day = shift_day + month_index(shift_month)
                if (shift_month == 1) then
                   shift_month = 12
                else
                   shift_month = shift_month - 1
                endif
             endif
          endif
       endif
    enddo
  end subroutine shift_time


  logical function leapyear(year)
    ! Determine if year is a leap-year
    use io_unit_defs, only: stderr
    implicit none
    ! Argument:
    integer, intent(in) :: year

    leapyear = .FALSE.
    if (year .gt. 1953 .and. year .lt. 2025) then
       if (year .eq. 1956 .or. year .eq. 1960 .or. year .eq. 1964 .or.  &
            year .eq. 1968 .or. year .eq. 1972 .or. year .eq. 1976 .or.  &
            year .eq. 1980 .or. year .eq. 1984 .or. year .eq. 1988 .or.  &
            year .eq. 1992 .or. year .eq. 1996 .or. year .eq. 2000 .or.  &
            year .eq. 2004 .or. year .eq. 2008 .or. year .eq. 2012 .or.  &
            year .eq. 2016 .or. year .eq. 2020 .or. year .eq. 2024) then
          leapyear = .TRUE.
       endif
    else
       write (stderr, *) 'Out of range for leap-year calculation. ', &
                         'Add more years to function leap year'
       call exit(1)
    endif
  end function leapyear


  subroutine alloc_forcing_variables(wind_n, met_n, river_n)
    ! Allocate memory for forcing variable arrays.
    use malloc, only: alloc_check
    implicit none
    ! Arguments:
    integer, intent(in) :: &
         wind_n,   &  ! Size of hourly wind speed forcing data arrays
         met_n,    &  ! Size of daily meteo forcing data arrays
         river_n      ! Size of daily river flow forcing data arrays
    ! Local variables:
    integer           :: allocstat  ! Allocation return status
    character(len=80) :: msg        ! Allocation failure message prefix

    msg = "Wind forcing data arrays"
    allocate(wind_eastnorth(wind_n), wind_northwest(wind_n), stat=allocstat)
    call alloc_check(allocstat, msg)
    msg = "Meteorological forcing data arrays"
    allocate(cf(met_n, 24), atemp(met_n, 24), humid(met_n, 24), stat=allocstat)
    call alloc_check(allocstat, msg)
    msg = "River flow forcing data arrays"
    allocate(Qriver(river_n), Eriver(river_n), stat=allocstat)
    call alloc_check(allocstat, msg)
  end subroutine alloc_forcing_variables


  subroutine dalloc_forcing_variables
    ! Deallocate memory for forcing variable arrays.
    use malloc, only: dalloc_check
    implicit none
    ! Local variables:
    integer           :: dallocstat  ! Deallocation return status
    character(len=80) :: msg         ! Deallocation failure message prefix

    msg = "Wind forcing data arrays"
    deallocate(wind_eastnorth, wind_northwest, stat=dallocstat)
    call dalloc_check(dallocstat, msg)
    msg = "Meteorological forcing data arrays"
    deallocate(cf, atemp, humid, stat=dallocstat)
    call dalloc_check(dallocstat, msg)
    msg = "River flow forcing data arrays"
    deallocate(Qriver, Eriver, stat=dallocstat)
    call dalloc_check(dallocstat, msg)
  end subroutine dalloc_forcing_variables

end module forcing
