module grid_mod
  ! Type definitions, variables, functions, and subroutines related to
  ! the grid.
  !
  ! Public Type:
  !
  !   grid_ -- Grid parameters:
  !              D -- depth of bottom of grid [m]
  !              M -- number of grid points
  !              d_g -- array of depths of grid points [m]
  !              d_i -- array of depths of grid cell interfaces [m]
  !              g_space -- array of depths of grid point spacings [m]
  !              i_space -- array of depths of grid cell interface spacings [m]
  !
  ! Public Variables:
  !
  !   grid   -- Grid point and interface depth and spacing arrays.
  !
  ! Public Functions:
  !
  !   depth_average(qty_g, d1, d2) -- Return the average value of the
  !                                   grid point stored quantity qty_g
  !                                   between depths d1, and d2.
  !
  !   full_depth_average(qty_g) -- Return the average value of the
  !                                grid point stored quantity qty_g
  !                                over the entire depth of the grid.
  !                                This is a special case wrapper for
  !                                the depth_average() function for
  !                                quantities defined as qty_g(1:M).
  !
  !   gradient_g(qty_g) -- Return the values of the gradients at the
  !                        grid point depths of a quantity stored
  !                        at the grid point depths.
  !
  !   gradient_i(qty_g) -- Return the values of the gradients at the
  !                        grid layer interface depths of a quantity stored
  !                        at the grid point depths.
  !
  !   interp_i(qty_g) -- Return the interpolated values of a quantity at
  !                      the grid interface depths from its values at the
  !                      grid point depths.
  !   interp_g(qty_i) -- Return the interpolated values of a quantity at
  !                      the grid pionts depths from its values at the
  !                      grid interface depths.
  !
  ! Public Subroutines:
  !
  !   dalloc_grid() -- Deallocate memory from grid arrays.
  !
  !   init_grid(M, grid) -- Allocate memory for grid arrays and
  !                         initialize grid interpolation factor arrays.
  !
  !   interp_value(value, lb, search_array, interp_array, interp_val, j_below)
  !      -- Calculate the value of the quantity in the interp_array that
  !         corresponds to the the specified value in the search_array,
  !         and the index of the grid point or interface immediately below
  !         the calculated value.  lb is the lower bound of
  !         search_array & interp_array.  This can be used to find
  !         the value of a quantity at a specified depth:
  !           call interp_value(d, lb, d_g, qty_g, qty_g(d), j_below) or
  !           call interp_value(d, lb, d_i, qty_i, qty_i(d), j_below)
  !         or to find the shallowest depth at which a specified quantity
  !         value occurs:
  !           call interp_value(q, lb, qty_g, d_g, d(q), j_below) or
  !           call interp_value(q, lb, qty_i, d_i, d(q), j_below)

  use precision_defs, only: dp
  implicit none

  private
  public :: &
       ! Types:
       grid_, &
       ! Variables:
       grid, &
       ! Functions:
       depth_average, full_depth_average, gradient_g, gradient_i, interp_i, &
       interp_g, &
       ! Subroutines:
       dalloc_grid, init_grid, interp_value

  ! Public type definition:
  !
  ! Grid parameters
  type :: grid_
     integer :: &
          M  ! Number of grid points
     real(kind=dp) :: &
          D  ! Depth of bottom of grid [m]
     real(kind=dp), dimension(:), allocatable :: &
          d_g,     &  ! Grid point depths [m]
          d_i,     &  ! Grid interface depths [m]
          g_space, &  ! Grid point spacing [m]
          i_space     ! Grid interface spacing [m]
  end type grid_

  ! Public Variable Declarations:
  type(grid_) :: &
       grid  ! Grid point and interface depth and spacing arrays

  ! Private Type Definition:
  !
  ! Interpolation factors
  type :: interp_factor
     real(kind=dp), dimension(:), allocatable :: &
          factor
  end type interp_factor

  ! Private module variable declarations:
  real(kind=dp) :: &
       lambda  ! Grid spacing parameter (<0 concentrates resolution
               ! near surface, =0 produced uniform grid, >0
               ! concentrates resolution near bottom of grid).  See
               ! Large, et al (1994), App. D.
  real(kind=dp), dimension(:), allocatable :: &
       xsi_i, &  ! Interface indices
       xsi_g     ! Grid point indices
  type(interp_factor) :: &
       above_g, &  ! Fraction of grid point spacing that interface j
                   ! is above grid point j+1.
       above_i, &  ! Fraction of interface spacing that gridpoint ***
       below_g, &  ! Fraction of grid point spacing that interface j
                   ! is below grid point j.
       below_i     ! Fraction of interface spacing that gridpoint  ***

contains

  subroutine init_grid
    ! Initialize grid.
    use precision_defs, only: dp, sp
    use io_unit_defs, only: stdout
    implicit none
    ! Local variables:
    integer :: j  ! Loop index over depth of grid

    ! Read the grid parameters
    call read_grid_params
    ! Allocate memory for grid arrays
    call alloc_grid

    ! Calculate the grid point and interface indices per Large, et al
    ! (1994) Appendix D (expressions are defined in the paragraph
    ! following the one containing eq'n D1)
    do j = 0, grid%M
       xsi_i(j) = dble(j) / dble(grid%M)
    enddo
    xsi_g(0) = 0.
    do j = 1, grid%M
       xsi_g(j) = (dble(j) - 0.5) / dble(grid%M)
    enddo
    xsi_g(grid%M + 1) = xsi_i(grid%M)

    ! Calculate the grid point and interface depths (Large, et al
    ! (1994) eq'n D1).  Note that abs(x) < epsilon(x) is a
    ! real-number-robust test for x == 0.
    if (abs(lambda) < epsilon(lambda)) then
       ! Uniform grid (lambda == 0)
       grid%d_g = grid%D * xsi_g
       grid%d_i = grid%D * xsi_i
    else
       ! Non-uniform grid (lambda < 0 concentrates resolution near
       ! surface, lambda > 0 concentrates resolution near D)
       ! *** Note that implementation of non-uniform grid is not
       ! *** consistent throughout the code
       write(stdout, *) "init_grid: Non-uniform grid is not ", &
            "fully implemented, lambda = ", lambda
       call exit(1)
       grid%d_g = (grid%D / lambda) * log(1. * xsi_g * (1. - exp(lambda)))
       grid%d_i = (grid%D / lambda) * log(1. * xsi_i * (1. - exp(lambda)))
    endif

    ! Calculate the grid point and interface spacings
    grid%g_space = grid%d_g(1:) - grid%d_g(0:grid%M)
    grid%i_space = grid%d_i(1:) - grid%d_i(0:grid%M-1)

    ! Calculate the grid interface interpolation factors
    above_g%factor(1:grid%M) = abs(grid%d_i(0:grid%M-1) - grid%d_g(1:grid%M)) &
         / grid%g_space(0:grid%M-1)
    above_g%factor(grid%M+1) = 0.
    below_g%factor(0) = 0.
    below_g%factor(1:grid%M) = abs(grid%d_g(1:grid%M) - grid%d_i(1:grid%M)) &
         / grid%g_space(1:grid%M)
    ! Calculate the grid point interpolation factors
    ! proportional to the other length
    above_i%factor(1:grid%M) = abs(grid%d_g(1:grid%M) - grid%d_i(1:grid%M)) &
         / grid%i_space(1:grid%M)
    below_i%factor(1:grid%M) = abs(grid%d_i(0:grid%M-1) - grid%d_g(1:grid%M)) &
         / grid%i_space(1:grid%M)
  end subroutine init_grid


  subroutine read_grid_params
    ! Read the grid parameters from the input file
    use input_processor, only: getpari, getpard
    implicit none

    ! Read the model depths, the number of grid points, and the grid
    ! spacing factor values.
    grid%D = getpard("maxdepth")
    grid%M = getpari("gridsize")
    lambda = getpard("lambda")
  end subroutine read_grid_params


  function depth_average(qty_g, d1, d2) result(avg)
    ! Return the average value of the grid point stored quantity qty_g
    ! between depths d1, and d2.
    use io_unit_defs, only: stdout
    implicit none
    ! Arguments:
    real(kind=dp), dimension(0:), intent(in) :: &
         qty_g  ! Quantity values at grid point depths
    real(kind=dp), intent(in) :: &
         d1, &  ! Beginning of depth range to calculate average over
         d2     ! End of depth range to calculate average over
    ! Result:
    real(kind=dp) :: avg  ! Average value
    ! Local variables:
    real(kind=dp) :: &
         q_d1, &    ! Value of quantity at depth d1
         q_d2, &    ! Value of quantity at depth d2
         q_i_d1, &  ! Value of qty at depth of grid interface layer below d1
         q_i_d2     ! Value of qty at depth of grid interface layer above d2
    integer :: &
         j_below_d1, &  ! Index of grid point below depth d1
         j_below_d2, &  ! Index of grid point below depth d2
         j_above_d2, &  ! Index of grid point above depth d2
         j_junk         ! Unused j_below argument of interp_value()

    ! Validate the depth arguments
    if (d1 < grid%d_g(0) .or. d2 > grid%D) then
       write(stdout, *) "depth_average(): depth(s) out of range: ", &
            "(", d1, ", ", d2, ")"
       call exit(1)
    endif
    ! Find the values of the quantity at the depth limits, and the
    ! indices of the grid points immediately below those depths
    call interp_value(d1, 0, grid%d_g, qty_g, q_d1, j_below_d1)
    call interp_value(d2, 0, grid%d_g, qty_g, q_d2, j_below_d2)
    ! Handle the special cases of d1 = d2 and q_d1 = q_d2.  Note that
    ! abs(x) < epsilon(x) is a real-number-robust test for x == 0.
    if (abs(d2 - d1) < epsilon(d2 - d1) &
         .and. abs(q_d2 - q_d1) < epsilon(q_d2 - q_d1)) then
       avg = q_d1
       return
    endif
    ! Set the indices of the grid points immediately above d2
    j_above_d2 = j_below_d2 - 1
    ! Check for d1 and d2 close together
    if (d2 < grid%d_g(j_below_d1+1)) then
       if (d2 <= grid%d_g(j_below_d1)) then
          ! In this case d1 and d2 are between the same pair of grid points
          avg = (q_d1 + q_d2) / 2.
       else
          ! In this case d1 and d2 are within three grid points
          avg = (0.5 * (q_d1 + qty_g(j_below_d1)) * (grid%d_g(j_below_d1)-d1) &
               + 0.5 * (qty_g(j_below_d1) + q_d2) * (d2-grid%d_g(j_below_d1))) &
               / (d2 - d1)
       endif
    else
       ! Calculate the depth weighted average value of the quantity
       ! through at least one grid point and then to the interface
       ! below that
       call interp_value(grid%d_i(j_below_d1), 0, grid%d_g, qty_g, &
            q_i_d1, j_junk)
       avg = 0.5 * (q_d1 + qty_g(j_below_d1)) * (grid%d_g(j_below_d1)-d1) &
            + 0.5 * (qty_g(j_below_d1) + q_i_d1) * &
            (grid%d_i(j_below_d1) - grid%d_g(j_below_d1))
       ! Portion between whole grid layer interfaces (note, if the second
       ! index is less that the first the dot product gets 0 (so we don't
       ! need to test to this)
       avg = avg + dot_product(qty_g(j_below_d1+1:j_above_d2-1), &
            grid%i_space(j_below_d1+1:j_above_d2-1))
       ! Portion below a layer interface through one
       ! grid point and to d2
       call interp_value(grid%d_i(j_above_d2-1), 0, grid%d_g, qty_g, &
         q_i_d2, j_junk)
       avg = avg + &
            0.5 * (q_i_d2 + qty_g(j_above_d2)) * &
            (grid%d_g(j_above_d2) - grid%d_i(j_above_d2-1)) + &
            0.5 * (qty_g(j_above_d2) + q_d2) * (d2 - grid%d_g(j_above_d2))
       avg = avg / (d2 - d1)
    endif
  end function depth_average


  function full_depth_average(qty_g) result(avg)
    ! Return the average value of the grid point stored quantity qty_g
    ! over the entire depth of the grid.  This is a special case
    ! wrapper for the depth_average() function for quantities defined
    ! as qty_g(1:M).
    use io_unit_defs, only: stdout
    implicit none
    ! Arguments:
    real(kind=dp), dimension(1:), intent(in) :: &
         qty_g  ! Quantity values at grid point depths
    ! Result:
    real(kind=dp) :: &
         avg  ! Average value
    ! Local Variables:
    real(kind=dp), dimension(0:grid%M+1) :: &
         qty_g_extended  ! qty_g array extended to include values at
                         ! depths 0 and grid%M+1.  Those values are
                         ! copies of the ones at 1 and grid%M.
    real(kind=dp) :: &
         d1, &  ! Beginning of depth range to calculate average over;
                ! i.e. 0.0d0
         d2     ! End of depth range to calculate average over;
                ! i.e. grid%D

    d1 = 0.0d0
    d2 = grid%D
    qty_g_extended(0) = qty_g(1)
    qty_g_extended(1:grid%M) = qty_g
    qty_g_extended(grid%M+1) = qty_g(grid%M)
    avg = depth_average(qty_g_extended, d1, d2)
  end function full_depth_average


  function gradient_g(qty_g) result(grad_g)
    ! Return the values of the gradients at the grid point
    ! depths of a quantity stored at the grid point depths.
    ! *** This function does not work for an unevenly spaced grid
    implicit none
    ! Argument:
    ! Quantity values at grid point depths
    real(kind=dp), dimension(0:), intent(in) :: qty_g
    ! Result:
    ! Values of the gradients at grid point depths
    real(kind=dp), dimension(1:grid%M) :: grad_g

    grad_g = (qty_g(0:grid%M-1) - qty_g(2:grid%M+1)) &
         / (grid%g_space(0:grid%M-1) + grid%g_space(1:))
  end function gradient_g


  function gradient_i(qty_g) result(grad_i)
    ! Return the values of the gradients at the grid layer interface
    ! depths of a quantity stored at the grid point depths.
    ! *** This function may not work for an unevenly spaced grid
    implicit none
    ! Argument:
    ! Quantity values at grid point depths
    real(kind=dp), dimension(0:), intent(in) :: qty_g
    ! Result:
    ! Values of the gradients at grid interface depths
    real(kind=dp), dimension(1:grid%M) :: grad_i

    grad_i = (qty_g(1:grid%M) - qty_g(2:grid%M+1)) / grid%g_space(1:grid%M)
  end function gradient_i


  function interp_i(qty_g) result(qty_i)
    ! Return the interpolated values of a quantity at the grid interface
    ! depths from its values at the grid point depths.
    implicit none
    ! Arguments:
    ! Quantity values at grid point depths
    real(kind=dp), dimension(0:), intent(in) :: qty_g
    ! Result:
    ! Quantity values at grid interface depths
    real(kind=dp), dimension(0:grid%M) :: qty_i

    qty_i = qty_g(0:size(qty_g) - 2) * above_g%factor &
         + qty_g(1:) * below_g%factor
  end function interp_i


  function interp_g(qty_i) result(qty_g)
    ! Return the interpolated values of a quantity at the grid point depths
    ! depths from its values at the grid interfaces.
    implicit none
    ! Arguments:
    ! Quantity values at grid interfaces
    real(kind=dp), dimension(0:), intent(in) :: qty_i
    ! Result:
    ! Quantity values at grid point depths
    real(kind=dp), dimension(1:grid%M) :: qty_g

    qty_g = qty_i(0:size(qty_g)-1) * above_i%factor &
         + qty_i(1:) * below_i%factor
  end function interp_g


  subroutine interp_value(value, lb, search_array, interp_array, &  ! in
       interp_val, j_below)                                         ! out
    ! Calculate the value of the quantity in the interp_array that
    ! corresponds to the the specified value in the search_array,
    ! and the index of the grid point or interface immediately below
    ! the calculated value.  This can be used to find the value of a
    ! quantity at a specified depth:
    !   call interp_value(d, 0, d_g, qty_g, qty_g(d), j_below) or
    !   call interp_value(d, 0, d_i, qty_i, qty_i(d), j_below)
    ! or to find the shallowest depth at which a specified quantity
    ! value occurs:
    !   call interp_value(q, lb_qty_g, qty_g, d_g, d(q), j_below) or
    !   call interp_value(q, lb_qty_i, qty_i, d_i, d(q), j_below)

    use io_unit_defs, only: stdout

    implicit none

    ! Arguments:
    real(kind=dp), intent(in) :: &
         value  ! Depth or quantity value to find
    integer, intent(in) :: &
         lb  ! Lower bound of search_array & interp_array
    real(kind=dp), dimension(lb:), intent(in) :: &
         search_array, &  ! Array to find value in
         interp_array     ! Array to interpolate over
    real(kind=dp), intent(out) :: &
         interp_val  ! Interpolated value
    integer, intent(out) :: &
         j_below  ! Index of grid point or layer interface immediately
                  ! below interpolated value

    ! Local variables:
    integer :: &
         ub, &  !  Upper bound of search_array
         i
    integer, dimension(lb:size(search_array) + lb - 1) :: &
         indices  ! Array of grid point/interface indices used to find
                  ! the index of the search array element below the
                  ! search value
    logical, dimension(lb:size(search_array) + lb - 1) :: &
         mask  ! Mask array used as a faster alternative to looping
               ! over the search array

    ! Mask the search array for values greater than or equal to the
    ! specified value.  We need to handle the 1st and last values in
    ! the search array separately so that we catch cases where they
    ! are exactly equal to the specified value.
    mask(lb) = search_array(lb) >= (value + epsilon(value))
    ub = size(search_array) + lb - 1
    mask(lb + 1 : ub - 1) = search_array(lb + 1 : ub - 1) >= value
    mask(ub) = search_array(ub) >= (value - epsilon(value))
    ! Confirm that the specified value is within the range of the
    ! search array.
    if (any(mask) .and. .not. all(mask)) then
       ! Find the index of the search array element below the
       ! specified value
       indices = (/ (i, i = lb, size(mask) + lb - 1) /)
       j_below = minval(indices, mask)
       ! Interpolate the value of the quantity in interp_array that
       ! corresponds to the specified value in the search array
       interp_val = interp_array(j_below - 1) &
            + (interp_array(j_below) - interp_array(j_below - 1)) &
            * (value - search_array(j_below - 1)) &
            / (search_array(j_below) - search_array(j_below - 1))
    else
       write(stdout, *) "Warning: value = ", value, &
            " out of range in interp_value()"
       call exit(1)
    endif
  end subroutine interp_value


  subroutine alloc_grid
    ! Allocate memory for grid arrays.
    use malloc, only: alloc_check
    implicit none
    ! Local variables:
    integer           :: allocstat  ! Allocation return status
    character(len=80) :: msg        ! Allocation failure message prefix

    msg = "Grid point and interface index arrays"
    allocate(xsi_g(0:grid%M+1), xsi_i(0:grid%M), stat=allocstat)
    call alloc_check(allocstat, msg)
    msg = "Grid point and interface depth and spacing arrays"
    allocate(grid%d_g(0:grid%M+1), grid%d_i(0:grid%M), &
         grid%g_space(0:grid%M), grid%i_space(1:grid%M), stat=allocstat)
    call alloc_check(allocstat, msg)
    msg = "Grid interface interpolation factor arrays"
    allocate(above_g%factor(1:grid%M+1), below_g%factor(0:grid%M), &
         above_i%factor(1:grid%M), below_i%factor(1:grid%M), &
         stat=allocstat)
    call alloc_check(allocstat, msg)
  end subroutine alloc_grid


  subroutine dalloc_grid
    ! Deallocate memory for grid arrays.
    use malloc, only: dalloc_check
    implicit none
    ! Local variables:
    integer           :: dallocstat  ! Allocation return status
    character(len=80) :: msg         ! Allocation failure message prefix

    msg = "Grid point and interface index arrays"
    deallocate(xsi_g, xsi_i, stat=dallocstat)
    call dalloc_check(dallocstat, msg)
    msg = "Grid point and interface depth and spacing arrays"
    deallocate(grid%d_g, grid%d_i, grid%g_space, grid%i_space, &
         stat=dallocstat)
    call dalloc_check(dallocstat, msg)
    msg = "Grid interface interpolation factor arrays"
    deallocate(above_g%factor, below_g%factor, above_i%factor, &
         below_i%factor, stat=dallocstat)
    call dalloc_check(dallocstat, msg)
  end subroutine dalloc_grid

end module grid_mod

