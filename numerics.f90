! $Id$
! $Source$

module numerics
  ! Type definitions, variable & parameter value declarations, and
  ! subroutines related to the numerical simulation aspects of the
  ! SOG code.
  !
  ! Public Type Definitions:
  !
  !   tridiag -- Tridiagonal matrix arrays.
  !
  ! Public Variables:
  !
  !   startDatetime -- Date/time of initial conditions
  !
  !   year -- Year counter
  !
  !   day --  Year-day counter
  !
  !   day_time -- Day-sec counter
  !
  !   time -- Time counter through run; seconds since
  !           startDatetime%day_sec
  !
  !   dt -- Time step [s]
  !
  !   steps --  Number of time steps in the main time loop
  !
  !   max_iter -- Maximum number of iterations allowed for implicit
  !               solver loop.
  !
  !   new_weight, prev_weight -- Weighting factors for blending of
  !                              values from current and previous
  !                              implicit solver iterations.
  !
  !   vel_scale, h_scale -- Scale factors used to normalize the
  !                         velocity components, and mixing layer
  !                         depth values to calculate the convergence
  !                         metrics for the implicit solver loop.
  !
  ! Public Subroutines:
  !
  !   check_negative -- Check for negative values in a vector.  Output
  !                     a message if a negative value is found, and
  !                     optionally stop execution; default is to stop.
  !
  !   dalloc_numerics_variables -- Deallocate memory for numerics
  !                                variables.

  use precision_defs, only: dp
  use datetime, only: datetime_
  implicit none

  private
  public :: &
       ! Type definitions:
       tridiag, &
       ! Variables:
       startDatetime, &   ! Date/time of initial conditions
       year, &   ! Year counter
       day,  &   ! Year-day counter
       day_time, &  ! Day-sec counter
       time, &  ! Time counter through run; seconds since
                ! startDatetime%day_sec
       dt, &      ! Time step [s]
       steps, &   ! Number of time steps in the main time loop
       max_iter, &  ! Maximum number of iterations allowed for
                    ! implicit solver loop.
       hprev, &  ! Values of h%new, u%new and v%new from previous
       Uprev, &  ! iteration for blending with new value to stabilize
       Vprev, &  ! convergence of the implicit solver loop
       new_weight, prev_weight, &  ! Weighting factors for convergence
                                   ! stabilization blending.
       vel_scale, h_scale,  &  ! Scale factors used to normalize the
                               ! velocity components, and mixing layer
                               ! depth values to calculate the
                               ! convergence metrics for the implicit
                               ! solver loop.
       del, &  ! Convergence metrics for implicit solver loop.
       ! Subroutines:
       init_numerics, check_negative, dalloc_numerics_variables

  ! Type Definitions:
  !
  ! Public:
  !
  ! Tridiagnonal matrix vectors:
  type :: tridiag
     real(kind=dp), dimension(:), allocatable :: &
          sub,  &  ! Sub-diagonal vector of a tridiagonal matrix
          diag, &  ! Diagonal vector of a tridiagonal matrix
          sup      ! Super-diagonal vector of a tridiagonal matrix
  end type tridiag
  !
  ! Private:
  !
  ! Convergence metrics for implicit solver loop.
  type :: converg_metrics
     real(kind=dp) :: &
          U, &  ! Cross-strait (35 deg) velocity component arrays
          V, &  ! Along-strait (305 deg) velocity component arrays
          h     ! Mising layer depth
  end type converg_metrics

  ! Variable Declarations:
  !
  ! Public:
  type(datetime_) :: &
       startDatetime   ! Date/time of initial conditions
  ! *** These should be replaced with a datetime structure
  integer :: &
       year, &   ! Year counter
       day   ! Year-day counter
  real(kind=dp) :: &
       day_time, &  ! Day-sec counter
       time  ! Time counter through run; seconds since
             ! startDatetime%day_sec
  real(kind=dp) :: &
       dt      ! Time step [s]
  integer :: &
       steps  ! Number of time steps in the main time loop
  integer :: &
       max_iter  ! Maximum number of iterations allowed for implicit
                 ! solver loop.
  real(kind=dp) :: &
       hprev  ! Value of h%new from previous iteration for blending
              ! with new value to stabilize convergence of the
              ! implicit solver loop.
  real(kind=dp), dimension (:), allocatable :: &
       Uprev, Vprev  ! Values of u%new and v%new from previous
                     ! iteration for blending with new value to
                     ! stabilize convergence of the implicit solver
                     ! loop.
  real(kind=dp)  :: &
       new_weight,  &  ! Weighting factors for blending of
       prev_weight, &  ! values from current and previous implicit
                       ! solver iterations.  Goal of blending is to
                       ! stabilize convergence.  It is applied to
                       ! mixing layer depth (h), and velocity
                       ! components (U & V).
       vel_scale,   &  ! Scale factors used to normalize the velocity 
       h_scale         ! components, and mixing layer depth values to
                       ! calculate the convergence metrics for the
                       ! implicit solver loop.
  type(converg_metrics) :: &
       del  ! Convergence metrics for implicit solver loop.
  !
  ! Private:
  real(kind=dp) :: &
       t_f  ! Run duration [s]  *** should be replaced with endDatetime

contains

  subroutine init_numerics(M)
    ! Allocate memory for numerics arrays, and read parameter values
    ! from the infile.

    ! Subroutines from other modules:
    use datetime, only: calendar_date, clock_time
    
    implicit none

    ! Argument:
    integer, intent(in) :: M  ! Number of grid points

    ! Allocate memory for numerics quantity arrays
    call alloc_numerics_variables(M)
    ! Read numerics parameter values from the infile.
    call read_numerics_params()
    ! Calculate the month number and month day for output file headers
    call calendar_date(startDatetime)
    call clock_time(startDatetime)
    ! Initialize the time counters
    year = startDatetime%yr
    day = startDatetime%yr_day
    day_time = dble(startDatetime%day_sec)
    time = dble(startDatetime%day_sec)
    ! Calculate the number of time steps for the run (note that int()
    ! rounds down)
    steps = 1 + int((t_f - dble(startDatetime%day_sec)) / dt)
  end subroutine init_numerics


  subroutine read_numerics_params()
    ! Read the fresh water parameter values from the infile.
    use input_processor, only: getpari, getpard
    implicit none

    ! Initial conditions date/time for the run
    startDatetime%yr = getpari("year_o")
    startDatetime%yr_day = getpari("yr_day_o")
    startDatetime%day_sec = getpari("t_o")
    ! Run duration, and time step
    t_f = getpard("run_dur")
    dt = getpard("dt")
    ! Maximum number of iterations allowed for implicit solver loop.
    max_iter = getpari('max_iter')     
  end subroutine read_numerics_params


  subroutine check_negative(lbound, vector, msg, day, time, fatal)
    ! Check for negative values in a vector.  Output a message if a
    ! negative value is found, and optionally stop execution; default
    ! is to stop.
    use precision_defs, only: dp
    use io_unit_defs, only: stdout
    implicit none
    ! Arguments:
    integer, intent(in) :: &
         lbound  ! Lower bound of vector
    real(kind=dp), dimension(lbound:), intent(inout) :: &
         vector  ! Vector of values to check; inout so -ves can be set to zero
    character(len=*) :: &
         msg  ! Vector specific part of message to print
    integer, intent(in) :: &
         day  ! Year-day of current time step
    real(kind=dp), intent(in) :: &
         time  ! Time of current time step [s since start of run]
    logical, intent(in), optional :: &
         fatal  ! Is a negative value a fatal error? [Default = .true.]
    ! Local variable
    logical :: &
         die    ! Stop execution if a -ve value is found?
    integer :: &
         j          ! Indices of most negative values in vector

    ! Handle optional fatal argument
    if (present(fatal)) then
       die = fatal
    else
       die = .true.  ! Default
    endif
    ! Check for negative value(s)
    if (any(vector < 0.0d0)) then
       ! Build message depending on whether or not -ve values are fatal
       if (die) then
          write(stdout, *) "Error: Negative value(s) in " // msg, &
               " at day = ", day, " time = ", time
       else
          write(stdout, *) "Warning: Negative value(s) in " // msg &
               // " set to zero at day = ", day, " time = ", time
       endif
       ! Find the index of the most -ve value
       do j = lbound, size(vector) + lbound - 1
          if (vector(j) < 0.0d0) then
             ! Output the message
             write(stdout, *) " j = ", j
          endif
       enddo
       ! Stop execution if -ve values are fatal
       if (die) then
          call exit(1)
       else
          ! If not a fatal errror, set the -ve value to zero
          where (vector < 0.0d0) vector = 0.0d0
       endif
    endif
  end subroutine check_negative


  subroutine alloc_numerics_variables(M)
    ! Allocate memory for numerics quantity arrays.
    use malloc, only: alloc_check
    implicit none
    ! Argument:
    integer, intent(in) :: M  ! Number of grid points
    ! Local variables:
    integer           :: allocstat  ! Allocation return status
    character(len=80) :: msg        ! Allocation failure message prefix

    msg = "Velocity component profile arrays from previous iteration"
    allocate(Uprev(0:M+1), Vprev(0:M+1), &
         stat=allocstat)
    call alloc_check(allocstat, msg)
  end subroutine alloc_numerics_variables


  subroutine dalloc_numerics_variables()
    ! Deallocate memory for numerics quantity arrays.
    use malloc, only: dalloc_check
    implicit none
    ! Local variables:
    integer           :: dallocstat  ! Deallocation return status
    character(len=80) :: msg         ! Deallocation failure message prefix

    msg = "Velocity component profile arrays from previous iteration"
    deallocate(Uprev, Vprev, &
         stat=dallocstat)
    call dalloc_check(dallocstat, msg)
  end subroutine dalloc_numerics_variables

end module numerics
       
