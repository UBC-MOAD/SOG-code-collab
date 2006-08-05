! $Id$
! $Source$

module unit_conversions
  ! A collection of functions that do unit conversions:
  !  CtoK -- degrees Celcius to Kelvins
  !  KtoC -- Kelvins to degrees Celcius

  implicit none

  ! *** We need interfaces for these functions to overload them for
  ! *** various types as appropriate (e.g. both reals and doubles)

contains

  real function CtoK(degC) result(K)
    ! Convert temperature from Celcius to Kelvin
    implicit none
    real, intent(in)  :: degC

    K = degC + 273.15
  end function CtoK


  double precision function KtoC(K) result(degC)
    ! Convert temperature from Kelvin to Celcius
    implicit none
    double precision, intent(in)  :: K
    
    degC = K - 273.15
  end function KtoC


end module unit_conversions
