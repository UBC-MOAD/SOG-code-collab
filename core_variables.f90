! $Id$
! $Source$

module core_variables
  ! Type declarations, variables, and subroutines related to the core
  ! variables that the code calculates.
  ! 
  ! Public Variables:
  !
  !   U -- Velocity component in the u (cross-strait, 35 deg) direction [m/s]
  !
  !   V -- Velocity component in the v (along-strait, 305 deg) direction [m/s]
  !
  !   T -- Water column temperature [K]
  !
  !   S -- Water column salinity [-]
  !
  !   P%micro -- Micro phytoplankton (diatoms) biomass [uM N]
  !
  !   P%nano -- Nano phytoplankton (flagellates) biomass [uM N]
  !
  !   N%O -- Nitrate concentration [uM N]
  !
  !   N%H -- Ammonium concentration [uM N]
  !
  !   Si -- Silicon concentration [uM]
  !
  ! Public Subroutines:
  !
  !   alloc_core_variables -- Allocate memory for core variables arrays.
  !
  !   dalloc_core_variables -- De-allocate memory for core variables arrays.

  use precision_defs, only: dp
  implicit none

  private
  public :: &
       ! Variables:
       T,  &  ! Temperature profile arrays
       S,  &  ! Salinity profile arrays
       P,  &  ! Micro & nano phytoplankton profile arrays
       N,  &  ! Nitrate & ammonium concentation profile arrays
       Si, &  ! Silicon concentration profile arrays
       ! Subroutines:
       alloc_core_variables, dalloc_core_variables

  ! Private module type definitions:
  !
  ! Velocities, temperature, and salinity
  type :: profiles
     real(kind=dp), dimension(:), allocatable :: &
          new, &  ! Profile of quantity at current time setp
          old, &  ! Profile of quantity at previous time step
          grad_i  ! Profile of gradient of quantity at grid layer interfaces
  end type profiles
  !
  ! Nitrogen compounds
  type :: nitrogen
     real(kind=dp), dimension(:), allocatable :: &
          O, &  ! N%O is nitrate (NO3) concentration profile
          H     ! H%H is ammonium (NH4) concentration profile
  end type nitrogen
  !
  ! Plankton
  type :: plankton
     real(kind=dp), dimension(:), allocatable :: &
          micro, &  ! P%micro is micro phytoplankton (diatoms) biomass profile
          nano      ! P%nano is nano phytoplankton (flagellate) biomass profile
  end type plankton


  ! Public variable declarations:
  type(profiles) :: &
       T, &  ! Temperature profile arrays
       S     ! Salinity profile arrays
  type(plankton) :: &
       P  ! Micro & nano phytoplankton profile arrays
  type(nitrogen) :: &
       N  ! Nitrate & ammonium profile arrays
  real(kind=dp), dimension(:), allocatable :: &
       Si ! Silicon concentration profile array

contains

  subroutine alloc_core_variables(M)
    ! Allocate memory for core variables arrays.
    use malloc, only: alloc_check
    implicit none
    ! Argument:
    integer, intent(in) :: M  ! Number of grid points
    ! Local variables:
    integer           :: allocstat  ! Allocation return status
    character(len=80) :: msg        ! Allocation failure message prefix

    msg = "Temperature profile arrays"
    allocate(T%new(0:M+1), T%old(0:M+1), T%grad_i(1:M), &
         stat=allocstat)
    call alloc_check(allocstat, msg)
    msg = "Salinity profile arrays"
    allocate(S%new(0:M+1), S%old(0:M+1), S%grad_i(1:M), &
         stat=allocstat)
    call alloc_check(allocstat, msg)
    msg = "Micro phytoplankton biomass profile arrays"
    allocate(P%micro(0:M+1), &
         stat=allocstat)
    call alloc_check(allocstat, msg)
    msg = "Nano phytoplankton biomass profile arrays"
    allocate(P%nano(0:M+1), &
         stat=allocstat)
    call alloc_check(allocstat, msg)
    msg = "Nitrate profile arrays"
    allocate(N%O(0:M+1), &
         stat=allocstat)
    call alloc_check(allocstat, msg)
    msg = "Ammonium profile arrays"
    allocate(N%H(0:M+1), &
         stat=allocstat)
    call alloc_check(allocstat, msg)
    msg = "Silicon concentration profile arrays"
    allocate(Si(0:M+1), &
         stat=allocstat)
    call alloc_check(allocstat, msg)
  end subroutine alloc_core_variables


  subroutine dalloc_core_variables
    ! Deallocate memory for core variables arrays.
    use malloc, only: dalloc_check
    implicit none
    ! Local variables:
    integer           :: dallocstat  ! Deallocation return status
    character(len=80) :: msg        ! Deallocation failure message prefix

    msg = "Temperature profile arrays"
    deallocate(T%new, T%old, T%grad_i, &
         stat=dallocstat)
    call dalloc_check(dallocstat, msg)
    msg = "Salinity profile arrays"
    deallocate(S%new, S%old, S%grad_i, &
         stat=dallocstat)
    call dalloc_check(dallocstat, msg)
    msg = "Micro phytoplankton biomass profile arrays"
    deallocate(P%micro, &
         stat=dallocstat)
    call dalloc_check(dallocstat, msg)
    msg = "Nano phytoplankton biomass profile arrays"
    deallocate(P%nano, &
         stat=dallocstat)
    call dalloc_check(dallocstat, msg)
    msg = "Nitrate profile arrays"
    deallocate(N%O, &
         stat=dallocstat)
    call dalloc_check(dallocstat, msg)
    msg = "Ammonium profile arrays"
    deallocate(N%H, &
         stat=dallocstat)
    call dalloc_check(dallocstat, msg)
    msg = "Silicon concentration profile arrays"
    deallocate(Si, &
         stat=dallocstat)
    call dalloc_check(dallocstat, msg)
  end subroutine dalloc_core_variables

end module core_variables
