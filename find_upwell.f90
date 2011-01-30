module find_upwell
  ! Subroutines that calculate the vertical profile of entrainment
  ! velocity (upwelling), and the vertical advection of quantities (T,
  ! S, etc.).
  !
  ! Public Subroutines:
  !
  ! init_find_upwell()
  !   -Read parameter values from the infile.
  !
  ! upwell_profile(grid, Qriver, upwell, wupwell)
  !   -Calculate the vertical profile of the entrainment velocity based
  !    on the strength of the River, and the maximum upwelling velocity.
  !
  ! vertical_advection(grid, dt, property, wupwell, gvector)
  !   -Calculate the vertical advection of a quantity based on the
  !    vertical profile of the entrainment velocity.
  use precision_defs, only: dp
  implicit none

  private
  public :: &
       ! subroutines:
       init_upwelling, w_upwell, &
       upwelling_profile, vertical_advection

  ! Variable Declarations:
  !
  ! Private:
  real(kind=dp) :: &
       d  ! depth in m, of 68% fresh water content
  real(kind=dp), dimension(:), allocatable :: &
       w_upwell  ! Vertical profile of the entrainment velocity

contains

  subroutine init_upwelling(M)
    ! Allocate memory for upwelling quantity arrays, and read
    ! parameter values from the infile.
    implicit none
    ! Argument:
    integer, intent(in) :: M  ! Number of grid points
  
    ! Allocate memory for upwelling quantity arrays
    call alloc_upwelling_variables(M)
    ! Read upwelling parameter values from the infile.
    call read_upwelling_params()
  end subroutine init_upwelling


  subroutine read_upwelling_params
    ! Read the fresh water depth parameter value from the infile.
    use input_processor, only: getpard
    implicit none
 
    ! Read values for d fit equation from infile
    ! fix d because fit is poor and not enough surface NO3 in winter
    d = getpard("d")
  end subroutine read_upwelling_params


  subroutine upwelling_profile(Qriver)
    ! Calculate the vertical profile of the entrainment velocity based
    ! on the river flow, and the maximum freshwater upwelling
    ! velocity.  The latter is a function of the river flow, and is
    ! calculated in the surface_flux subroutine.  The details of the
    ! model to do that are in the the document entrainment.pdf written
    ! in late June 2006.
    use precision_defs, only: dp
    use grid_mod, only: grid
    use freshwater, only: upwell  ! Maximum freshwater upwelling velocity

    implicit none

    ! Argument:
    real(kind=dp), intent(in) :: Qriver  ! River flow
    ! Local variables:
    real(kind=dp) :: d25  ! Depth of upwelling variation
    integer :: index      ! Grid step index

    ! Three models for the upwelling variation depth have been used:
    !   1. Use a value based on the current salinity profile. Found to
    !      be not sufficiently stable. See code version 1.8 for
    !      calculation of d from the salinity profile.
    !   2. Set the depth of 68% fresh water content to the fit found
    !      in entrainment.pdf.
!!$    d = 11.7 * (Qriver/2720.)**(-0.23)
    !   3. A value of d calculated according to entrainment.pdf for
    !      the specific estuary being modeled is read from the infile.
 
   ! Depth of upwelling variation is defined as 2.5*d (see entrainment.pdf)
    d25 = 2.5 * d
    ! Set upwelling velocity profile (on interfaces!)
    do index = 0, grid%M
       if (grid%d_g(index) < d25) then
          w_upwell(index) = upwell * (1. - ( (1.- grid%d_i(index) / d25)**2) )
       else
          w_upwell(index) = upwell
       endif
    enddo
  end subroutine upwelling_profile


  subroutine vertical_advection (grid, dt, qty, wupwell, gvector)
    ! Calculate the vertical advection of a quantity based on the
    ! vertical profile of the entrainment velocity.  Note that the
    ! vertical velocity decreases as the surface is approaches and so
    ! water and property is squeezed out the sides of the water column.

    use precision_defs, only: dp
    use grid_mod, only: grid_

    implicit none

    ! Arguments:
    type(grid_), intent(in) :: grid 
    real(kind=dp), intent(in) :: dt ! time step
    real(kind=dp), dimension (0:), intent(in) :: qty ! quantity to be advected
    real(kind=dp), intent(in), dimension(0:) :: wupwell
    real(kind=dp), intent(inout), dimension(1:) :: Gvector

    ! Local variables:
    integer :: index  ! counter through depth
    real(kind=dp) :: inbottom, outtop, squashing, valuelost

    outtop = 0.
    do index = 1, grid%M
       inbottom = wupwell(index) * qty(index+1) ! upwind scheme
       squashing = wupwell(index) - wupwell(index-1)
       if (-squashing*dt > grid%i_space(index)) then
          write (*,*) "Problem in find_upwell, upwelling too strong for dt"
          stop
       endif
       valuelost = ( qty(index) * & 
            (grid%i_space(index) - wupwell(index-1) * dt) + & 
            qty(index+1) * (wupwell(index) * dt) ) / &
            (grid%i_space(index) + squashing * dt)
       Gvector(index) = Gvector(index) + dt / grid%i_space(index) &
            * (inbottom - outtop - squashing * valuelost)
       outtop = inbottom    ! of previous cell
    enddo
  end subroutine vertical_advection


  subroutine alloc_upwelling_variables(M)
    ! Allocate memory for upwelling variable arrays.
    use malloc, only: alloc_check
    implicit none
    ! Argument:
    integer, intent(in) :: M  ! Number of grid points
    ! Local variables:
    integer           :: allocstat  ! Allocation return status
    character(len=80) :: msg        ! Allocation failure message prefix

    msg = "Vertical profile of entrainment velocity array"
    allocate(w_upwell(0:M), &
         stat=allocstat) 
    call alloc_check(allocstat, msg)
  end subroutine alloc_upwelling_variables


  subroutine dalloc_upwelling_variables()
    ! Deallocate memory for upwelling variable array.
    use malloc, only: dalloc_check
    implicit none
    ! Local variables:
    integer           :: dallocstat  ! Deallocation return status
    character(len=80) :: msg         ! Deallocation failure message prefix

    msg = "Vertical profile of entrainment velocity array"
    deallocate(w_upwell, &
         stat=dallocstat)
    call dalloc_check(dallocstat, msg)
  end subroutine dalloc_upwelling_variables

end module find_upwell
