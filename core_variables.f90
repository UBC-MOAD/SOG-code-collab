! $Id$
! $Source$

module core_variables
  ! Type declarations, variables, and subroutines related to the core
  ! variables that the code calculates.
  ! 
  ! Public Variables:
  !
  ! U -- Profile of velocity component in the u (cross-strait) direction [m/s]
  !
  ! T -- Profile of water column temperature [K]
  !
  ! Public Subroutines:
  !
  ! alloc_core_variables -- Allocate memory for core variables arrays.
  !
  ! dalloc_core_variables -- De-allocate memory for core variables arrays.

  use precision_defs, only: dp
  implicit none

  private
  public :: &
       ! Variables:
       T, &  ! Temperature profile arrays
       ! Subroutines:
       alloc_core_variables, dalloc_core_variables

  ! Core variables type definition:
  type :: quantity
     real(kind=dp), dimension(:), pointer :: &
          new, &  ! Profile of quantity at current time setp
          old, &  ! Profile of quantity at previous time step
          grad_i  ! Profile of gradient of quantity at grid layer interfaces
  end type quantity

  ! Public variable declarations:
  type(quantity) :: &
       T  ! Temperature profile arrays

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

    msg = "Water column temperature profile arrays"
    allocate(T%new(0:M+1), T%old(0:M+1), T%grad_i(1:M), &
         stat=allocstat)
    call alloc_check(allocstat, msg)
  end subroutine alloc_core_variables


  subroutine dalloc_core_variables
    ! Allocate memory for core variables arrays.
    use malloc, only: dalloc_check
    implicit none
    ! Local variables:
    integer           :: dallocstat  ! Allocation return status
    character(len=80) :: msg        ! Allocation failure message prefix

    msg = "Water column temperature profile arrays"
    deallocate(T%new, T%old, T%grad_i, &
         stat=dallocstat)
    call dalloc_check(dallocstat, msg)
  end subroutine dalloc_core_variables

end module core_variables
