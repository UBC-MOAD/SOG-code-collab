! $Id$
! $Source$

module fundamental_constants
  ! Parameter value declarations of fundamental constants widely used
  ! in other modules of the SOG code.
  !
  ! Public Parameters:
  !
  !   f -- Coriolis factor
  !
  !   g -- Acceleration due to gravity [m/s^2]
  !
  !   latitude -- Station S3 latitude [deg]
  !
  !   pi -- Ratio of circumference to diameter of a circle [-]

  use precision_defs, only: dp
  implicit none

  private
  public :: &
       ! Parameter values:
       f,        &  !  Coriolis factor
       g,        &  ! Acceleration due to gravity [m/s^2]
       latitude, &  ! Station S3 latitude [deg]
       pi,       &  ! Ratio of circumference to diameter of a circle [-]
       ! Subroutine:
       init_constants

  ! Public parameter declarations:
  real(kind=dp) :: &
       f  ! Coriolis factor (would be a parameter but for a pgf90 bug)
  real(kind=dp), parameter :: &
       g = 9.80665, &                  ! Acceleration due to gravity [m/s^2]
       latitude = 49. + 7.517 / 60., & ! Station S3 latitude [deg]
       pi = 3.141592653589793

contains

  subroutine init_constants()
    ! Coriolis factor
    ! *** This must be calculated because pgf90 will not accept an
    ! *** intrinsic in parameter statement
    f = 2. * (2. * pi / 86400.) * sin(pi * latitude / 180.)
  end subroutine init_constants
  
end module fundamental_constants
