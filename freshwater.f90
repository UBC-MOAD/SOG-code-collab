! $Id$
! $Source$

module freshwater
  ! Type definitions, variable & parameter value declarations, and
  ! subroutines related to the fresh water flow influences SOG code.
  !
  ! Public Variables:
  !
  !   Fw_surface -- Add all of the fresh water on the surface?
  !
  !   Ft -- Total fresh water flux
  !
  !   F_n -- Fresh water contribution to salinity flux
  !
  ! Public Subroutines:
  !
  !   init_freshwater -- Allocate memory for fresh water quantity
  !                      arrays, and read parameter values from the
  !                      infile.
  !
  !   freshwater_phys -- Calculate the strength of
  !                      upwelling/entrainment, the freshwater flux,
  !                      the surface turbulent kinematic salinity
  !                      flux, and the profile of fresh water
  !                      contribution to the salinity flux.
  !
  !   freshwater_bio -- Calculate the freshwater biological fluxes.

  use precision_defs, only: dp
  implicit none

  private
  public :: &
       ! Variables:
       Fw_surface, &  ! Add all of the fresh water on the surface?
       Ft,         &  ! Total fresh water flux
       F_n,        &  ! Fresh water contribution to salinity flux
       upwell,     &  ! Upwelling velocity from river flows
                      ! parameterization.
       ! Diagnostics:
       S_riv, &  ! Surface salinity prediction from fit
       ! Subroutines:
       init_freshwater, freshwater_phys, freshwater_bio, &
       dalloc_freshwater_variables

  ! Parameter Value Declarations:
  !
  ! Private:
  !
  ! Circulation strength of a scalar is equal to deep value minus
  ! river value
  real(kind=dp), parameter :: phys_circ_nitrate = 30.5d0 - 13.0d0
  real(kind=dp), parameter :: phys_circ_Pmicro = 0.0d0
  real(kind=dp), parameter :: phys_circ_Pnano = 0.0d0
  real(kind=dp), parameter :: phys_circ_Zoo = 0.0d0
  real(kind=dp), parameter :: phys_circ_silicon = 54.0d0 - 80.0d0

  ! Variable Declarations:
  !
  ! Public:
  logical :: &
       Fw_surface  ! Add all of the fresh water on the surface?
  real(kind=dp) :: &
       Ft,  &  ! Total fresh water flux
       upwell  ! Upwelling velocity from river flows parameterization
  real(kind=dp), dimension(:), allocatable :: &
       F_n  ! Fresh water contribution to salinity flux
  !
  ! Diagnostic:
  real(kind=dp) :: &
       S_riv  ! Surface salinity prediction from fit
  !
  ! Private:
  logical :: &
       use_Fw_nutrients  ! Include influence of Fw nutrients?
  real(kind=dp), dimension(:), allocatable :: &
       Fw  ! Fresh water flux profile
  real(kind=dp) :: &
       Fw_scale, &   ! Fresh water scale factor for river flows
       Fw_depth, &   ! Depth to distribute fresh water flux over
       upwell_const  ! Maximum upwelling velocity (tuning parameter)

contains

  subroutine init_freshwater(M)
    ! Allocate memory for fresh water quantity arrays, and read
    ! parameter values from the infile.
    implicit none
    ! Argument:
    integer, intent(in) :: M  ! Number of grid points

    ! Allocate memory for fresh water quantity arrays
    call alloc_freshwater_variables(M)
    ! Read fresh water parameter values from the infile.
    call read_freshwater_params()
  end subroutine init_freshwater


  subroutine read_freshwater_params()
    ! Read the fresh water parameter values from the infile.
    use input_processor, only: getpard, getparl
    implicit none

    ! Maximum upwelling velocity (tuning parameter)
    upwell_const = getpard("upwell_const")
    ! Fresh water scale factor for river flows
    Fw_scale = getpard('Fw_scale')     
    ! Add all fresh water on surface?
    Fw_surface = getparl('Fw_surface') 
    ! Depth to distribute freshwater flux over
    if (.not. Fw_surface) then
       Fw_depth = getpard('Fw_depth')  
    endif
    ! Include effect of Fw nutrients?
    use_Fw_nutrients = getparl('use_Fw_nutrients')
  end subroutine read_freshwater_params


  subroutine freshwater_phys(Qriver, Eriver, S_old, h)
    ! Calculate the strength of upwelling/entrainment, the freshwater
    ! flux, the surface turbulent kinematic salinity flux, and the
    ! profile of fresh water contribution to the salinity flux.

    ! Elements from other modules:
    !
    ! Type Definitions:
    use precision_defs, only: sp
    ! Variables:
    use grid_mod, only: &
         grid  ! Grid parameters and depth & spacing arrays
    use turbulence, only: &
         wbar  ! Turbulent kinematic flux profile arrays

    implicit none

    ! Arguments
    real(kind=sp), intent(in) :: &
         Qriver, &  ! Fraser River flow
         Eriver     ! Englishman River flow
    real(kind=dp), intent(in) :: &
         S_old,        &  ! Surface salinity
         h                ! Mixing layer depth

    ! Local variables
    real(kind=dp), parameter :: &
         Qmean = 2720.0d0 ! Mean fraser river flow from entrainment fit

    ! Parameterized fit of the surface salinity of the Strait of
    ! Georgia at station S3 based on the river flows.  (Derived by
    ! Kate Collins 16-Jun-2005)  This value is not directly used
    ! in the model but is used to make sure the tuned Ft value
    ! is correct.
    S_riv = 29.1166d0 - Qriver * (0.0019d0) - Eriver * (0.0392d0)
    
    ! Tuned fresh water flux value (to give, on average) the parameterized
    ! value above.  
    Ft = Fw_scale * (0.0019d0 * Qriver + 0.0392d0 * Eriver) * &
         (Qriver/ Qmean) ** (0.d0-0.41d0)
    
    ! The entrainment of deep water into the bottom of the
    ! grid is based on the parameterization derived by Susan Allen in
    ! Jun-2006 (See entrainment.pdf) exponent is 0.41 from entrainment
    ! reduced to 0.25 based on match to Olivier and more NO3 needed in winter
    upwell = upwell_const * (Qriver / Qmean) ** 0.d0

    ! Calculate the surface turbulent kinematic salinity flux (Large,
    ! et al (1994), eqn A2d).  Note that fresh water flux is added via
    ! Bf in calc_buoyancy().
    if (Fw_surface) then
       wbar%s(0) = Ft * S_old
    else
       wbar%s(0) = 0.0d0
    endif

    ! Calculate the nonturbulent fresh water flux profile, and its
    ! contribution to the salinity profile
    if (Fw_surface) then
       F_n = 0.0d0
    else
       Fw = Ft * exp(-grid%d_i / (Fw_depth * h))
       F_n = 29.626d0 * Fw
    endif
  end subroutine freshwater_phys


  subroutine freshwater_bio(qty, current_value, surf_flux, distrib_flux)
    ! Calculate the freshwater biological fluxes.

    use grid_mod, only: grid
    implicit none
    
    character (len=*), intent(in) :: qty
    real (kind=dp), dimension(0:), intent(in) :: current_value
    real (kind=dp), intent(out) :: surf_flux
    real (kind=dp), dimension(0:), intent(out) :: distrib_flux

    integer :: i

    if (use_Fw_nutrients) then
       if (Fw_surface) then
          distrib_flux = 0.
          if (qty.eq."nitrate") then
             ! if phyto are using the nitrate the grad is zero'd. 
             ! 4.0 = 2 x (Michalis-Menton K)
             surf_flux = Ft * phys_circ_nitrate * &
                  min(4.0, current_value(1))/ 4.0
          elseif (qty.eq."Pmicro") then
             surf_flux = Ft * (phys_circ_Pmicro)
          elseif (qty.eq."Pnano") then
             surf_flux = Ft * (phys_circ_Pnano)
          elseif (qty.eq."Zoo") then
             surf_flux = Ft * (phys_circ_Zoo)
          elseif (qty.eq."silicon") then
             surf_flux = Ft * phys_circ_silicon * &
                  min(4.0, current_value(1)) / 4.0
          else
             write (*,*) "problems in freshwater, river flux for ", qty, &
                  " is not defined."
             call exit(1)
          endif
       else ! distributed fresh water and fluxes
          surf_flux = 0.
          if (qty.eq."nitrate") then
             ! vector variation does not give the same results as the loop
             do i=0,grid%M
                ! if phyto are using the nitrate the grad is zero'd. 
                ! 4.0 = 2 x (Michalis-Menton K)
                distrib_flux(i) = Fw(i) * phys_circ_nitrate * &
                     min(4.0,current_value(i))/4.0
             enddo
          elseif (qty.eq."Pmicro") then
             distrib_flux = Fw * (phys_circ_Pmicro)
          elseif (qty.eq."Pnano") then
             distrib_flux = Fw * (phys_circ_Pnano)
          elseif (qty.eq."Zoo") then
             distrib_flux = Fw * (phys_circ_Zoo)
          elseif (qty.eq."silicon") then
             do i=0,grid%M
                distrib_flux(i) = Fw(i) * phys_circ_silicon * &
                     min(4.0, current_value(i))/4.0
             enddo
          else
             write (*,*) "problems in freshwater, river flux for ", qty, &
                  " is not defined."
             call exit(1)
          endif
       endif
    else
       distrib_flux = 0
       surf_flux = 0
    endif

  end subroutine freshwater_bio


  subroutine alloc_freshwater_variables(M)
    ! Allocate memory for fresh water variable arrays.
    use malloc, only: alloc_check
    implicit none
    ! Argument:
    integer, intent(in) :: M  ! Number of grid points
    ! Local variables:
    integer           :: allocstat  ! Allocation return status
    character(len=80) :: msg        ! Allocation failure message prefix

    msg = "Fresh water flux and salinity contribution profile arrays"
    allocate(Fw(0:M), F_n(0:M), &
         stat=allocstat)
    call alloc_check(allocstat, msg)
  end subroutine alloc_freshwater_variables


  subroutine dalloc_freshwater_variables()
    ! Deallocate memory for fresh water variable arrays.
    use malloc, only: dalloc_check
    implicit none
    ! Local variables:
    integer           :: dallocstat  ! Deallocation return status
    character(len=80) :: msg         ! Deallocation failure message prefix

    msg = "Fresh water flux and salinity contribution profile arrays"
    deallocate(Fw, F_n, &
         stat=dallocstat)
    call dalloc_check(dallocstat, msg)
  end subroutine dalloc_freshwater_variables

end module freshwater
