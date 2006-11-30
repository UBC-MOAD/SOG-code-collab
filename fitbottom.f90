! $Id$
! $Source$

module fitbottom
  ! Constants, functions, and subroutines for calculation of boundary
  ! condition values at bottom of grid.
  !
  ! Public Subroutines:
  !
  ! bot__bound_time (day, day_time, Tbot, Sbot, Nobot, Sibot, Pmbot, Pnbot)
  !  -- Calculate the values at the bottom of the grid for those
  !     quantities that we have data for from an annual fit.
  !
  ! bot_bound_uniform (M, Z, NH, Detritus)
  !  -- Set the values at the bottom of the grid for those quantities
  ! that we don't have time series data for.

  use precision_defs, only: dp

  implicit none

  private
  public :: bot_bound_time, bot_bound_uniform

  ! Private module variable declarations:
  ! 
  ! Number of different quantities
  integer, parameter :: NQ = 6
  ! Quantity names
  character(len=12), dimension(NQ), parameter :: quantity &
       = ['salinity', 'temperature', 'chl fluor', 'nitrate', 'silicon', &
          'plank ratio']
  ! Fit coefficients (see fitbottom.py and fitted.m in
  ! sog-forcing/bottom for details of how these coefficient values
  ! were derived.
  real(kind=dp), dimension(5, NQ), parameter :: c &
       = reshape([ &
          ! Salinity (constant, seasonal and biseasonal components)
          29.62605865, 0.10374454, -0.03562458, -0.14156091, -0.06348989, &
          ! Temperature (constant, seasonal and biseasonal components)
          9.34995751, -0.27245442, -0.90151197, 0.17933473, -0.05590939, &
          ! Phytoplankton from fluor (constant and seasonal)
          0.58316137, -0.11206845, 0.26241523, 0., 0., &
          ! Nitrate (constant) !*** not from fit... what SEA thinks
          30.0, 0., 0., 0., 0., & 
          ! Silicon (constant) !*** not from fit... what SEA thinks
          50.0, 0., 0., 0., 0., &
          ! Ratio of pico/micro plankton (constant and biseasonal)
          1.25948868, 0., 0., 0.97697686, 0.46289294 &
         ], [5, NQ])

contains

  function bottom_value (arg, qty) result(value)
    ! Calculate the value of the specified quantitity at the bottom
    ! grid boundary for the given value of arg.
    use precision_defs, only: dp
    use io_unit_defs, only: stdout
    implicit none
    ! Arguments:
    real(kind=dp), intent(in):: arg
    character(len=*), intent(in):: qty
    ! Result:
    real(kind=dp) :: value
    ! Local variable:
    integer:: index

    if (qty == quantity(1)) then
       index = 1
    elseif (qty == quantity(2)) then
       index = 2
    elseif (qty == quantity(3)) then
       index = 3
    elseif (qty == quantity(4)) then
       index = 4
    elseif (qty == quantity(5)) then
       index = 5
    elseif (qty == quantity(6)) then
       index = 6
    else
       write (stdout,*) 'bottom_value in fitbottom.f90:', &
            'Unexpected quantity: ', qty
       call exit(1)
    endif

    value = c(1,index) + c(2,index)*cos(arg) + c(3,index)*sin(arg) &
         + c(4,index)*cos(2.*arg) + c(5,index)*sin(2.*arg)
  end function bottom_value


  subroutine bot_bound_time (day, day_time, &
       Tbot, Sbot, Nobot, Sibot, Pmbot, Pnbot)
    ! Calculate the values at the bottom of the grid for those
    ! quantities that we have data for from an annual fit.
    use precision_defs, only: dp
    use unit_conversions, only: CtoK
    use fundamental_constants, only: pi
    implicit none
    ! Arguments:
    integer, intent(in) :: day
    real(kind=dp), intent(in) :: day_time
    real(kind=dp), intent(out) :: Tbot, Sbot, Nobot, Sibot, Pmbot, Pnbot
    ! Local variables:
    real(kind=dp) :: arg, chl, ratio
    
    arg = 2 * pi * (day + day_time / 86400.) / 365.25
    
    Tbot = bottom_value(arg, 'temperature')          ! in Celcius
    Tbot = CtoK(Tbot)
    Sbot = bottom_value(arg, 'salinity')
    Nobot = bottom_value(arg, 'nitrate')
    Sibot = bottom_value(arg, 'silicon')
    
    chl = bottom_value(arg, 'chl fluor')
    ratio = bottom_value(arg, 'plank ratio')
    
    Pmbot = chl / (ratio + 1)
    Pnbot = ratio * Pmbot
  end subroutine bot_bound_time
  

  subroutine bot_bound_uniform (M, Z, NH, D_DON, D_PON, D_refr, D_bSi)
    ! Set the values at the bottom of the grid for those
    ! quantities that we don't have time series data for.
    use precision_defs, only: dp
    implicit none
    ! Arguments:
    integer, intent(in) :: M  ! Number of grid points
    real(kind=dp), intent(inout), dimension(0:) :: &
         Z,      &  ! microzooplantkon
         NH,     &  ! Ammonium concentration profile array
         D_DON,  &  ! Dissolved organic nitrogen concentration profile array
         D_PON,  &  ! Particulate organic nitrogen concentration profile array
         D_refr, &  ! Refractory nitrogen concentration profile array
         D_bSi      ! Biogenic silicon concentration profile array

    Z(M+1) = min(Z(M), 5.0d0)
    NH(M+1) = min(NH(M), 0.5d0)
    D_DON(M+1) = D_DON(M)
    D_PON(M+1) = D_PON(M)
    D_refr(M+1) = D_refr(M)
    D_bSi(M+1) = D_bSi(M)
  end subroutine bot_bound_uniform
  
end module fitbottom
