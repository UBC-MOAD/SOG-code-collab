module baroclinic_pressure
  ! Type definitions, variable & parameter value declarations, and
  ! subroutines related to the baroclinic pressure gradient
  ! calculation in the SOG code.
  !
  ! Public Variables:
  !
  !   dPdx_b -- Baroclinic pressure gradient x (cross-strait) component
  !             [kg/m^2 . s^2]
  !
  !   dPdy_b -- Baroclinic pressure gradient y (along-strait) component
  !             [kg/m^2 . s^2]
  !
  !   ut -- Integral of u (cross-strait) velocity component over time
  !
  !   vt -- Integral of v (along-strait) velocity component over time
  !
  ! Public Subroutines:
  !
  !   init_baroclinic_pressure -- Initialize baroclinic pressure gradient
  !                               calculation quantities.
  !
  !   baroclinic_P_gradient -- Calculate baroclinic pressure gradient
  !                            components.
  !
  !   new_to_old_vel_integrals -- Copy %new components of the velocity
  !                               component integrals to %old for use in
  !                               next time step.

  use precision_defs, only: dp
  implicit none

  private
  public :: &
       ! Variables:
       dPdx_b, &  ! Baroclinic pressure gradient x (cross-strait) component
       dPdy_b, &  ! Baroclinic pressure gradient y (along-strait) component
       ut,     &  ! Layer expansion (m) at eastern side
       vt,     &  ! layer expansion (m) at western side
       w_wind, &  ! wind-driven upwelling (m/s)
       ! Types (as required by new pg compiler)
       new_old_arrays, & ! type for ut, vt
       ! Subroutines:
       init_baroclinic_pressure, baroclinic_P_gradient, &
       new_to_old_vel_integrals, dalloc_baro_press_variables

  ! Type Definitions:
  !
  ! Private:
  !
  ! New/old array components:
  type :: new_old_arrays
     real(kind=dp), dimension(:), allocatable :: &
          new, &  ! Current time step values
          old     ! Previous time step values
  end type new_old_arrays

  ! Parameter Value Declarations:
  !
  ! Private:
  real(kind=dp) :: &
       Lx, &  ! Minor axis (cross-strait) of model domain [m]
       Ly     ! Major axis (along-strait) of model domain [m]
  logical :: openEnd ! Is this a fjord system? If so, we use an open ended basin system.

  ! Variable Declarations:
  !
  ! Public:
  real(kind=dp), dimension(:), allocatable :: &
       dPdx_b, &  ! Baroclinic pressure gradient x (cross-strait) component
       dPdy_b, &  ! Baroclinic pressure gradient y (along-strait) component
       w_wind     ! wind-driven upwelling
  type(new_old_arrays) :: &
       ut, &  ! Layer expansion (m) on eastern side
       vt     ! Layer expansion (m) on northern side
  !
  ! Private:
  real(kind=dp), dimension(:), allocatable :: &
       dze, &  ! Isopycnal depth profile array on eastern side
       dzw, &  ! Isopycnal depth profile array on western side
       dzn, &  ! Isopycnal depth profile array on northern side
       dzs     ! Isopycnal depth profile array on southern side

contains

  subroutine init_baroclinic_pressure(M)
    ! Initialize baroclinic pressure gradient calculation quantities.
    use input_processor, only: getpard,getparl
    implicit none
    ! Argument:
    integer :: M  ! Number of grid points

    ! Initializes x and y model domain axes
    Lx = getpard('Lx')
    Ly = getpard('Ly')
    openEnd =getparl('openEnd')
    ! Allocate memory for baroclinic pressure gradient calculation
    ! quantities
    call alloc_baro_press_variables(M)
    ! Initialize the baroclinic pressure gradient components, and the
    ! velocity component integrals used in their calculation
    dPdx_b = 0.0d0
    dPdy_b = 0.0d0
    w_wind = 0.0d0
    ut%new = 0.0d0
    vt%new = 0.0d0
  end subroutine init_baroclinic_pressure


  subroutine baroclinic_P_gradient(dt, U_new, V_new, rho_g)
    ! Calculate the baroclinic pressure gradient components.

    ! Elements from other modules:
    !
    ! Type Definitions:
    use precision_defs, only: dp
    ! Functions:
    use grid_mod, only: depth_average
    ! Parameter Values:
    use fundamental_constants, only: g
    ! Variable Declarations:
    use grid_mod, only: &
         grid, &       ! Grid parameters and depth & spacing arrays
         interp_value
    use turbulence, only: &
         wbar  ! Turbulent kinematic flux profile arrays; we need
               ! wbar%u(0) & wbar%v(0)

    implicit none

    ! Agruments:
    real(kind=dp), intent(in) :: &
         dt    ! Time step [s]
    real(kind=dp), dimension(0:), intent(inout) :: &
         U_new, &  ! U component (cross-strait, 35 deg) velocity profile [m/s]
         V_new, &  ! V component (along-strait, 305 deg) velocity profile [m/s]
         rho_g     ! Density profile at grid point depths [kg/m^3]

    ! Local parameter:
    !
    ! Surface layer depth of calculation of wind-driven upwelling
    ! effect on baroclinic pressure gradient in open-ended (fjord)
    ! domains like Rivers Inlet. The value of 15 m is based on Knight
    ! Inlet.
    real(kind=dp), parameter :: surface_layer_depth = 15.0d0

    ! Local variables:
    real(kind=dp) :: &
         decayscale, &  ! to keep the isopycnals less tilted than size
                        ! of domain
         unused         ! Unused value calculated by interp_value
    real(kind=dp), dimension(0:grid%M) :: kink  ! in m, kink the in isopycnals
    integer :: &
         gridbotsurf, &  ! grid point at bottom of surface layer
         yy              ! Loop index over depth of grid

    ! Remove barotropic mode from velocity field
    U_new = U_new - depth_average(U_new, 0.0d0, grid%D)
    V_new = V_new - depth_average(V_new, 0.0d0, grid%D)

    ! Calculate the distortion of the isopycnals at the grid layer
    ! interfaces that the velocity field produces.
    !
    ! Calculate integrals of u and v over time.  This is from a
    ! geometrical argument that:
    !
    !    u * grid%i_space * dt = ut * grid*i_space * 2 / (Lx / 2)
    !
    ! Thus, ut(j) is the amount in m that the u velocity deepens the
    ! isopycnal.  That is, the extra water goes into a triangular
    ! region of height ut, and length 1/2 width of the strait.
    !
    ! If the velocities get large enough the displacement of the
    ! isopycnals reaches the 40 m scale of the model.  To keep the
    ! isopycnals displacements reasonable, say 6 m, introduce a decay which
    ! particularly impacts high flow situations (by being proportional to
    ! the square) but only takes the 6 m displacements down slowly (90 days).
    !
    ! *** This statement probably needs to be modified if the grid
    ! *** depth is changed from 40 m.
    decayscale = sum(abs(ut%old(1:grid%M))) ** 2 / (90.0d0 * 86400.0d0) &
         / (6.0d0 ** 2)

    ! Calculate the added thickness of each layer at the "east" and
    ! "north" side of basin, respectively, in metres.
    ut%new(1:grid%M) = ut%old(1:grid%M) * (1.0d0 - decayscale * dt) &
         + U_new(1:grid%M) * dt * grid%i_space(1:grid%M) &
           * (2.0d0 / (Lx / 2.0d0))
    vt%new(1:grid%M) = vt%old(1:grid%M) * (1.0d0 - decayscale * dt) &
         + V_new(1:grid%M) * dt * grid%i_space(1:grid%M) &
           * (2.0d0 / (Ly / 2.0d0))

    ! Get grid point index at the depth of the bottom of the surface
    ! layer
    call interp_value(surface_layer_depth, 0, grid%d_g, grid%d_g, &
         unused, gridbotsurf)
    gridbotsurf = gridbotsurf - 1

    ! Calculate the depths of the distorted isopycnals in metres.
    call iso_distortion(ut%new, dze)
    call iso_distortion(-ut%new, dzw)
    call iso_distortion(vt%new, dzn)
    if (openend) then
       call iso_distortion_openend(-vt%new, dzs, gridbotsurf)
    else
       call iso_distortion(-vt%new, dzs)
    endif

    ! Make sure values are not below 40m
    dze = min(dze, grid%D)
    dzw = min(dzw, grid%D)
    dzn = min(dzn, grid%D)
    dzs = min(dzs, grid%D)

    ! Calculate the baroclinic pressure gradient
    call delta_p(dze, dzw, rho_g, dPdx_b)
    call delta_p(dzn, dzs, rho_g, dPdy_b)

    ! At this point we have calculated delta P.  To make it into the
    ! pressure gradient term we need to multiply by g, divide by rho
    ! and the distance between the sides and the center of the domain.
    !
    ! Also, to keep the flow purely baroclinic, take out the
    ! barotropic effect of the wind forcing
    dPdx_b = dPdx_b * g / (Lx * rho_g(1:grid%M)) - wbar%u(0) / grid%D
    dPdy_b = dPdy_b * g / (Ly * rho_g(1:grid%M)) - wbar%v(0) / grid%D


    if (openEnd) then
       ! Calculate the upwelling driven by the wind and its correction
       ! on ut. kink is the change in orientation of the isopycnals
       ! along the fjord
       kink(0) = 0.d0
       w_wind(0) = 0.d0
       do yy = 1, gridbotsurf
          kink(yy) = - dzn(yy) - dzs(yy) + 2 * grid%d_i(yy)
          kink(yy) = max(min(kink(yy), 1.d0), 0.d0)
          if (V_new(yy) < 0) then  ! velocity above the kinked isopycnal
             w_wind(yy) = w_wind(yy-1) &
                  - 2 * grid%i_space(yy) * V_new(yy) / Ly * kink(yy)
          else
             w_wind(yy) = w_wind(yy-1)
          endif
          vt%new(yy) = vt%new(yy) + (w_wind(yy) - w_wind(yy-1)) * dt
       enddo
       do yy = gridbotsurf + 1, grid%M
          w_wind(yy) = w_wind(yy-1)
       enddo
    endif
    !  Remove barotropic mode
    vt%new = vt%new - depth_average(vt%new, 0.d0, grid%D)
  end subroutine baroclinic_P_gradient


  subroutine iso_distortion(ut, dz)
    ! Calculate the depths of the distorted isopycnals in metres.

    ! Variables from other modules:
    use grid_mod, only: &
         grid  !  Grid parameters and depth & spacing arrays
    implicit none
    ! Arguments:
    real(kind=dp), intent(in), dimension(1:) :: &
         ut
    real(kind=dp), intent(out), dimension(1:) :: &
         dz
    ! Local Variables:
    real(kind=dp) :: &
         resid_u
    integer :: &
         yy

    if (ut(1) > -grid%i_space(1)) then
       dz(1) = ut(1) + grid%i_space(1)
       resid_u = 0.0d0
    else
       dz(1) = 0.0d0
       resid_u = ut(1) + grid%i_space(1)
    endif
    do yy = 2, grid%M
       if (ut(yy) + resid_u > -grid%i_space(yy)) then
          dz(yy) = dz(yy-1) + ut(yy) + grid%i_space(yy) + resid_u
          resid_u = 0.0d0
       else
          dz(yy) = dz(yy-1)
          resid_u = ut(yy) + grid%i_space(yy) + resid_u
       endif
    enddo
  end subroutine iso_distortion

  subroutine iso_distortion_openend(ut, dz, gridbotsurf)
    ! Calculate the depths of the distorted isopycnals in metres, assuming this end is open

    ! Variables from other modules:
    use grid_mod, only: &
         grid  !  Grid parameters and depth & spacing arrays
    implicit none
    ! Arguments:
    real(kind=dp), intent(in), dimension(1:) :: &
         ut
    real(kind=dp), intent(out), dimension(1:) :: &
         dz
    integer, intent(in) :: gridbotsurf
    ! Local Variables:
    real(kind=dp) :: &
         resid_u, &
         extra
    integer :: &
         yy

    if (ut(1) < -grid%i_space(1)) then
       dz(1) = 0.0d0
       resid_u = ut(1) + grid%i_space(1)
    elseif (ut(1) < 0.d0) then
       dz(1) = ut(1) + grid%i_space(1)
       resid_u = 0.0d0
    else
       dz(1) = grid%i_space(1)
       resid_u = ut(1)
    endif

    do yy = 2, gridbotsurf
       if (ut(yy) + resid_u < -grid%i_space(yy)) then
          dz(yy) = dz(yy-1)
          resid_u = ut(yy) + grid%i_space(yy) + resid_u
       elseif (dz(yy-1) + ut(yy) + resid_u < 0.25*(yy-1)) then  ! note explicit grid size!!!
          dz(yy) = dz(yy-1) + ut(yy) + grid%i_space(yy) + resid_u
          resid_u = 0.0d0
       else
          dz(yy) = dz(yy-1) + grid%i_space(yy)
          resid_u = ut(yy) + resid_u
       endif
    enddo

    extra = 0.d0
    if (resid_u > 0.d0) then  ! spread it over the bottom region
       extra = resid_u/(grid%M-gridbotsurf)
       resid_u = 0.d0
    endif

    do yy = gridbotsurf + 1, grid%M
       if (ut(yy) + resid_u + extra > -grid%i_space(yy)) then
          dz(yy) = dz(yy-1) + ut(yy) + grid%i_space(yy) + resid_u + extra
          resid_u = 0.0d0
       else
          dz(yy) = dz(yy-1)
          resid_u = ut(yy) + grid%i_space(yy) + resid_u + extra
       endif
    enddo


  end subroutine iso_distortion_openend


  subroutine delta_p (dzpos, dzneg, rho, dPres)
    use grid_mod, only: full_depth_average
    use grid_mod, only: &
         grid
    implicit none
    real(kind=dp), intent(in), dimension(1:) :: dzpos, dzneg
    real(kind=dp), intent(in), dimension(0:) :: rho
    real(kind=dp), intent(out), dimension(1:) :: dPres

    real(kind=dp) :: cpos, cneg
    integer :: ii, jj, yy

    cpos = 0.0d0
    cneg = 0.0d0
    ii = 1
    jj = 1

    ! calculated on grid points
    ! dPdx_b = P(eastern side) - P(western side)
    ! P(eastern side) = sum rho_g * eastern side thickness
    ! but note that this sum must go down to the depth grid%d_g(yy)
    ! and not below

    do yy = 1, grid%M
       ! Initialize
       if (yy == 1) then
          dPres(1) = 0.0d0
       else
          dPres(yy) = dPres(yy-1)
       endif

       ! First, add the pressure on the eastern side.
       do while (dzpos(ii) - grid%d_g(yy) < epsilon(dzpos(ii) - grid%d_g(yy))&
                .and. ii < grid%M)
          dPres(yy) = dPres(yy) + rho(ii) * (dzpos(ii) - cpos)
          cpos = dzpos(ii)
          ii = ii + 1
       enddo
       dPres(yy) = dPres(yy) + rho(ii) * (grid%d_g(yy) - cpos)
       cpos = grid%d_g(yy)

       ! Then, substract the pressure on the western side.
       do while (dzneg(jj) - grid%d_g(yy) < epsilon(dzneg(jj) - grid%d_g(yy))&
            .and. jj < grid%M)
          dPres(yy) = dPres(yy) - rho(jj) * (dzneg(jj) - cneg)
          cneg = dzneg(jj)
          jj = jj + 1
       enddo
       dPres(yy) = dPres(yy) - rho(jj) * (grid%d_g(yy) - cneg)
       cneg = grid%d_g(yy)
    enddo

    ! At this point we have pressure on the grid points, assuming that the
    ! pressure at the surface is zero.
    !
    ! If we assume that the average pressure should be zero, then we
    ! need to subtract the average of the calculated pressure.
    dPres = dPres - full_depth_average(dPres)
  end subroutine delta_p


  subroutine new_to_old_vel_integrals()
    ! Copy %new component of the velocity component integrals to %old
    ! component for use at next time step.

    ut%old = ut%new
    vt%old = vt%new
  end subroutine new_to_old_vel_integrals


  subroutine alloc_baro_press_variables(M)
    ! Allocate memory for baroclinic pressure gradient calculation
    ! quantities.
    use malloc, only: alloc_check
    implicit none
    ! Argument:
    integer, intent(in) :: M  ! Number of grid points
    ! Local variables:
    integer           :: allocstat  ! Allocation return status
    character(len=80) :: msg        ! Allocation failure message prefix

    msg = "Baroclinic pressure gradient profile arrays"
    allocate(dPdx_b(1:M), dPdy_b(1:M), &
         stat=allocstat)
    call alloc_check(allocstat, msg)
    msg = "Wind-driven upwelling array"
    allocate(w_wind(0:M), &
         stat=allocstat)
    call alloc_check(allocstat, msg)
    msg = "Velocity component integral profile arrays"
    allocate(ut%new(0:M+1), ut%old(0:M+1), vt%new(0:M+1), vt%old(0:M+1), &
         stat=allocstat)
    call alloc_check(allocstat, msg)
    msg = "Distorted isopicnal depth profile arrays"
    allocate(dze(1:M), dzw(1:M), dzn(1:M), dzs(1:M),  &
         stat=allocstat)
    call alloc_check(allocstat, msg)
  end subroutine alloc_baro_press_variables


  subroutine dalloc_baro_press_variables
    ! Deallocate memory for baroclinic pressure gradient calculation
    ! quantities.
    use malloc, only: dalloc_check
    implicit none
    ! Local variables:
    integer           :: dallocstat  ! Deallocation return status
    character(len=80) :: msg         ! Deallocation failure message prefix

    msg = "Baroclinic pressure gradient profile arrays"
    deallocate(dPdx_b, dPdy_b, &
         stat=dallocstat)
    call dalloc_check(dallocstat, msg)
    msg = "Wind-driven upwelling arrays"
    deallocate(w_wind, &
         stat=dallocstat)
    call dalloc_check(dallocstat, msg)

    msg = "Velocity component integral profile arrays"
    deallocate(ut%new, ut%old, vt%new, vt%old, &
         stat=dallocstat)
    call dalloc_check(dallocstat, msg)
    msg = "Distorted isopicnal depth profile arrays"
    deallocate(dze, dzw, dzn,dzs, &
         stat=dallocstat)
    call dalloc_check(dallocstat, msg)
  end subroutine dalloc_baro_press_variables

end module baroclinic_pressure
