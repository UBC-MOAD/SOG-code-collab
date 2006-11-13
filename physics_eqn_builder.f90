! $Id$
! $Source$

module physics_eqn_builder
  ! Type definitions, variable declarations, and
  ! subroutines related to the semi-implicit diffusion/advection PDEs
  ! (with Coriolis, and baroclinic pressure gradient terms) for the
  ! physics quantities in the SOG code.
  !
  ! Public Variables:
  !
  !   Bmatrix -- Tridiagonal matrix of diffusion coefficient values
  !              precursor to the LHS A matrix in the Aq = h matrix
  !              equation (see Large, et al (1994), App. D, pg 398)
  !
  !   Right-hand side vector arrays; precursor to the RHS h vector in
  !   the Aq = h matrix equation (see Large, et al (1994), App. D, pg 398):
  !
  !   U_RHS -- U velocity component (cross-strait) right-hand side arrays
  !
  !   V_RHS -- V velocity component (along-strait) right-hand side arrays
  !
  !   T_RHS -- Temperature right-hand side arrays
  !
  !   S_RHS -- Salinity right-hand side arrays
  !
  ! Public Subroutines:
  !
  !   build_physics_equations -- Build the right-hand side (RHS) arrays
  !                              for the diffusion/advection/Coriolis/
  !                              baroclinic pressure gradient equations for
  !                              the physics quantities.
  !
  !   new_to_old_phys_RHS -- Copy %new component of the physics *_RHS%diff_adv
  !                          arrays to the %old component for use by the
  !                          IMEX semi-impllicit PDE solver.
  !
  !   new_to_old_phys_Bmatrix -- Copy %new component of the Bmatrix%phys
  !                              arrays to the %old component for use by the
  !                              IMEX semi-impllicit PDE solver.
  !
  !   alloc_phys_RHS_variables -- Allocate memory for physics RHS arrays
  !
  !   dalloc_phys_RHS_variables -- Deallocate memory for physics RHS arrays

  use precision_defs, only: dp
  use numerics, only: tridiag  ! Tridiagonal matrix arrays type def'n
  implicit none

  private
  public :: &
       ! Variables:
       Bmatrix, &  ! Tridiagonal matrix of diffusion coefficient values
       U_RHS,   &  ! U velocity component (cross-strait) right-hand side arrays
       V_RHS,   &  ! V velocity component (along-strait) right-hand side arrays
       T_RHS,   &  ! Temperature right-hand side arrays
       S_RHS,   &  ! Salinity right-hand side arrays
       ! Subroutines:
       build_physics_equations, &
       new_to_old_phys_RHS, new_to_old_phys_Bmatrix, &
       alloc_phys_RHS_variables, dalloc_phys_RHS_variables

  ! Type Definitions:
  !
  ! Private to module:
  !
  ! New/old tridiagonal matrix components:
  type :: new_old
     type(tridiag) :: &
          new, &  ! Current time step values
          old     ! Previous time step values
  end type new_old
  !
  ! Diffusion coefficient matrix type:
  type :: diff_coeffs_matrix
     type(new_old) :: &
          vel, &  ! Velocity components diffusion coefficients matrix
          T,   &  ! Temperature diffusion coefficients matrix
          S       ! Salinity diffusion coefficients matrix
  end type diff_coeffs_matrix
  !
  ! New/old array components:
  type :: new_old_arrays
     real(kind=dp), dimension(:), allocatable :: &
          new, &  ! Current time step values
          old     ! Previous time step values
  end type new_old_arrays
  !
  ! Semi-implicit diffusion/advection/Coriolis/baroclinic pressure
  ! gradient equation right-hand side (RHS) arrays
  type :: RHS
     type(new_old_arrays) :: &
          diff_adv, &  ! Diffusion/advection component of RHS
          C_pg         ! Coriolis and baroclinic pressure grad component of RHS
  end type RHS

  ! Variable Declarations:
  !
  ! Public:
  type(diff_coeffs_matrix) :: &
       Bmatrix  ! Tridiagonal matrix of diffusion coefficient values
                ! Components:
                !   Bmatrix%vel
                !          %T
                !          %S
                !              %new
                !              %old
                !                  %sub
                !                  %diag
                !                  %sup
  type(RHS) :: &
       U_RHS, &  ! U velocity component (cross-strait) RHS arrays
       V_RHS, &  ! V velocity component (along-strait) RHS arrays
       T_RHS, &  ! Temperature RHS arrays
       S_RHS     ! Salinity RHS arrays
  !
  ! Private to module:

contains

  subroutine build_physics_equations(grid, dt, U, V, T, S)  ! in
    ! Build the terms for the diffusion/advection/Coriolis/baroclinic
    ! pressure gradient equations for the physics quantities.
    !
    ! This calculates the values of the precursor diffusion
    ! coefficients matrices (Bmatrix%vel%*, Bmatrix%T%*, and
    ! Bmatrix%S%*), the RHS diffusion/advection term vectors
    ! (*_RHS%diff_adv%new), and the RHS Coriolis and barolcinic
    ! pressure gradient term vectors (*_RHS%C_pg).

    ! Element from other modules:
    !
    ! Type Definitions:
    use precision_defs, only: dp
    use grid_mod, only: grid_
    ! Variable Declarations:
    use baroclinic_pressure, only: &
         dPdx_b, &  ! Cross-strait component of baroclinic pressure gradient
         dPdy_b     ! Along-strait component of baroclinic pressure gradient

    implicit none

    ! Arguments:
    type(grid_), intent(in) :: &
         grid  ! Grid arrays
    real(kind=dp), intent(in) :: &
         dt  ! Time step [s]
    real(kind=dp), dimension(0:), intent(in) :: &
         U, &  ! U velocity component(cross-strait) profile
         V, &  ! V velocity component(along-strait) profile
         T, &  ! Temperature profile
         S

    ! Calculate the precursor diffusion coefficients matrices for the
    ! physics model quantities.
    call calc_phys_Bmatrices(grid, dt)  ! in

    ! Initialize the RHS *%diff_adv%new arrays, and calculate the diffusive
    ! fluxes at the bottom and top of the grid
    call init_phys_RHS_fluxes(grid, dt, U, V, T, S)  ! in

!!$    ! Add vertical advection due to upwelling to velocity components,
!!$    ! temperature, and salinity RHS arrays
!!$    call calc_phys_upwell_advection(grid, dt, U, V, T, S)

    ! Calculate the Coriolis and baroclinic pressure gradient
    ! components of the RHS arrays for the velocity components
    call Coriolis_and_pg(dt, -U, dPdy_b, &  ! in
         V_RHS%C_pg%new)                  ! out
    call Coriolis_and_pg(dt, V, dPdx_b, &   ! in
         U_RHS%C_pg%new)                  ! out
  end subroutine build_physics_equations


  subroutine calc_phys_Bmatrices(grid, dt)  !  in
    ! Calculate the precursor diffusion coefficients matrices for the
    ! physics model quantities.

    ! Element from other modules:
    !
    ! Type Definitions:
    use precision_defs, only: dp
    use grid_mod, only: grid_
    ! Variable Declarations:
    use declarations, only: K  ! *** Should change to use turbulence
    ! Subroutines:
    use diffusion, only: diffusion_coeff
    
    implicit none

    ! Arguments:
    type(grid_), intent(in) :: &
         grid  ! Grid arrays
    real(kind=dp), intent(in) :: &
         dt  ! Time step [s]

    ! Velocity components
    call diffusion_coeff(grid, dt, K%u%all, & ! in
         Bmatrix%vel%new)                     ! out
    ! Temperature
    call diffusion_coeff(grid, dt, K%T%all, & ! in
         Bmatrix%T%new)                       ! out
    ! Salinity
    call diffusion_coeff(grid, dt, K%S%all, & ! in
         Bmatrix%S%new)                       ! out
  end subroutine calc_phys_Bmatrices


  subroutine init_phys_RHS_fluxes(grid, dt, U, V, T, S)  ! in
    ! Initialize the RHS *%diff_adv%new arrays, and calculate the diffusive
    ! fluxes at the bottom and top of the grid

    ! Element from other modules:
    !
    ! Type Definitions:
    use precision_defs, only: dp
    use grid_mod, only: grid_
    ! Variable Declarations:
    use declarations, only: K, w, gamma  ! *** Should change to use turbulence
    use declarations, only: Q_n, F_n  ! *** Should come from somewhere else
    ! Subroutines:
    use diffusion, only: diffusion_bot_surf_flux, diffusion_nonlocal_fluxes
    
    implicit none
    
    ! Arguments:
    type(grid_), intent(in) :: &
         grid  ! Grid arrays
    real(kind=dp), intent(in) :: &
         dt  ! Time step [s]
    real(kind=dp), dimension(0:), intent(in) :: &
         U, &  ! U velocity component(cross-strait) profile
         V, &  ! V velocity component(along-strait) profile
         T, &  ! Temperature profile
         S

    ! Velocity components
    call diffusion_bot_surf_flux (grid, dt, K%U%all, w%u(0), &  ! in
         U(grid%M+1), &                                         ! in
         U_RHS%diff_adv%new)                                    ! out
    call diffusion_bot_surf_flux (grid, dt, K%U%all, w%v(0), &  ! in
         V(grid%M+1),                                        &  ! in
         V_RHS%diff_adv%new)                                    ! out
    ! Temperature
    call diffusion_nonlocal_fluxes(grid, dt, K%T%all, gamma%t, &  ! in
         w%T(0), -Q_n, T(grid%M+1),                            &  ! in
         T_RHS%diff_adv%new)                                      ! out
    ! Salinity
    call diffusion_nonlocal_fluxes(grid, dt, K%S%all, gamma%s, &  ! in
         w%s(0), F_n, S(grid%M+1),                             &  ! in
         S_RHS%diff_adv%new)                                      ! out
  end subroutine init_phys_RHS_fluxes


!!$  subroutine calc_phys_upwell_advection(grid, dt, U, V, T, S)
!!$    ! Add vertical advection due to upwelling to velocity components,
!!$    ! temperature, and salinity RHS arrays
!!$    use precision_defs, only: dp
!!$    use grid_mod, only: grid_
!!$    use upwelling, only: upwell_profile, upwelling_advection
!!$    implicit none
!!$    ! Arguments:
!!$    type(grid_), intent(in) :: &
!!$         grid  ! Grid arrays
!!$    real(kind=dp), intent(in) :: &
!!$         dt  ! Time step [s]
!!$    real(kind=dp), dimension(0:), intent(in) :: &
!!$         U, &  ! U velocity component(cross-strait) profile
!!$         V, &  ! V velocity component(along-strait) profile
!!$         T, &  ! Temperature profile
!!$         S
!!$
!!$    ! Calculate profile of upwelling velocity
!!$    call upwell_profile(grid)
!!$    ! Apply upwelling advection
!!$    !
!!$    ! Velocity components
!!$    call upwelling_advection (grid, dt, U, &  ! in
!!$         U_RHS%diff_adv%new)                  ! out
!!$    call upwelling_advection (grid, dt, V, &  ! in
!!$         V_RHS%diff_adv%new)                  ! out
!!$    ! Temperature
!!$    call upwelling_advection (grid, dt, S, &  ! in
!!$         T_RHS%diff_adv%new)                  ! out
!!$    ! Salinity
!!$    call upwelling_advection (grid, dt, T, &  ! in
!!$         S_RHS%diff_adv%new)                  ! out
!!$  end subroutine calc_phys_upwell_advection


  subroutine Coriolis_and_pg(dt, vel, P_grad, RHS)
    ! Calculate the Coriolis and baroclinic pressure gradient
    ! components of the RHS vector for the specified velocity component.
    use fundamental_constants, only: f
    implicit none
    ! Arguments:
    real(kind=dp), intent(in) :: &
         dt  ! Time step [s]
    real(kind=dp), dimension(0:), intent(in) :: &
         vel  ! Velocity component profile
    real(kind=dp), dimension(1:), intent(in) :: &
         P_grad  ! Baroclinic pressure gradient component profile
    real(kind=dp), dimension(1:), intent(out) :: &
         RHS  ! Velocity component right-hand side array

    RHS = (f * vel(1:) - P_grad) * dt
  end subroutine Coriolis_and_pg


  subroutine new_to_old_phys_RHS()
    ! Copy %new component of the physics *_RHS%* arrays to the
    ! %old component for use by the IMEX semi-impllicit PDE solver.
    implicit none

    ! Diffusion/advection components
    U_RHS%diff_adv%old = U_RHS%diff_adv%new
    V_RHS%diff_adv%old = V_RHS%diff_adv%new
    T_RHS%diff_adv%old = T_RHS%diff_adv%new
    S_rhs%diff_adv%old = S_RHS%diff_adv%new
    ! Coriolis and baroclinic pressure gradient components
    U_RHS%C_pg%old = U_RHS%C_pg%new
    V_RHS%C_pg%old = V_RHS%C_pg%new
  end subroutine new_to_old_phys_RHS


  subroutine new_to_old_phys_Bmatrix()
    ! Copy %new component of the physics Bmatrix%* arrays to the
    ! %old component for use by the IMEX semi-impllicit PDE solver.
    implicit none

    Bmatrix%vel%old%sub = Bmatrix%vel%new%sub
    Bmatrix%vel%old%diag = Bmatrix%vel%new%diag
    Bmatrix%vel%old%sup = Bmatrix%vel%new%sup
    Bmatrix%T%old%sub = Bmatrix%T%new%sub
    Bmatrix%T%old%diag = Bmatrix%T%new%diag
    Bmatrix%T%old%sup = Bmatrix%T%new%sup
    Bmatrix%S%old%sub = Bmatrix%S%new%sub
    Bmatrix%S%old%diag = Bmatrix%S%new%diag
    Bmatrix%S%old%sup = Bmatrix%S%new%sup
  end subroutine new_to_old_phys_Bmatrix


  subroutine alloc_phys_RHS_variables(M)
    ! Allocate memory for arrays for right-hand sides of
    ! diffusion/advection equations for the physics model.
    use malloc, only: alloc_check
    implicit none
    ! Argument:
    integer, intent(in) :: M  ! Number of grid points
    ! Local variables:
    integer           :: allocstat  ! Allocation return status
    character(len=80) :: msg        ! Allocation failure message prefix

    msg = "Diffusion coefficients tridiagonal matrix arrays"
    allocate(Bmatrix%vel%new%sub(1:M), Bmatrix%vel%new%diag(1:M), &
         Bmatrix%vel%new%sup(1:M),                                &
         Bmatrix%vel%old%sub(1:M), Bmatrix%vel%old%diag(1:M),     &
         Bmatrix%vel%old%sup(1:M),                                &
         Bmatrix%T%new%sub(1:M), Bmatrix%T%new%diag(1:M),         &
         Bmatrix%T%new%sup(1:M),                                  &
         Bmatrix%T%old%sub(1:M), Bmatrix%T%old%diag(1:M),         &
         Bmatrix%T%old%sup(1:M),                                  &
         Bmatrix%S%new%sub(1:M), Bmatrix%S%new%diag(1:M),         &
         Bmatrix%S%new%sup(1:M),                                  &
         Bmatrix%S%old%sub(1:M), Bmatrix%S%old%diag(1:M),         &
         Bmatrix%S%old%sup(1:M),                                  &
         stat=allocstat)
    call alloc_check(allocstat, msg)
    msg = "U velocity component (cross-strait) RHS arrays"
    allocate(U_RHS%diff_adv%new(1:M), U_RHS%diff_adv%old(1:M), &
         U_RHS%C_pg%new(1:M), U_RHS%C_pg%old(1:M),             &
         stat=allocstat)
    call alloc_check(allocstat, msg)
    msg = "V velocity component (cross-strait) RHS arrays"
    allocate(V_RHS%diff_adv%new(1:M), V_RHS%diff_adv%old(1:M), &
         V_RHS%C_pg%new(1:M), V_RHS%C_pg%old(1:M),             &
         stat=allocstat)
    call alloc_check(allocstat, msg)
    msg = "Temperature RHS arrays"
    allocate(T_RHS%diff_adv%new(1:M), T_RHS%diff_adv%old(1:M), &
         stat=allocstat)
    call alloc_check(allocstat, msg)
    msg = "Salinity RHS arrays"
    allocate(S_RHS%diff_adv%new(1:M), S_RHS%diff_adv%old(1:M), &
         stat=allocstat)
    call alloc_check(allocstat, msg)
  end subroutine alloc_phys_RHS_variables


  subroutine dalloc_phys_RHS_variables
    ! Deallocate memory from arrays for right-hand sides of
    ! diffusion/advection equations for the physics model.
    use malloc, only: dalloc_check
    implicit none
    ! Local variables:
    integer           :: dallocstat  ! Allocation return status
    character(len=80) :: msg         ! Allocation failure message prefix

    msg = "Diffusion coefficients tridiagonal matrix arrays"
    deallocate(Bmatrix%vel%new%sub, Bmatrix%vel%new%diag, &
         Bmatrix%vel%new%sup,                             &
         Bmatrix%vel%old%sub, Bmatrix%vel%old%diag,       &
         Bmatrix%vel%old%sup,                             &
         Bmatrix%T%new%sub, Bmatrix%T%new%diag,           &
         Bmatrix%T%new%sup,                               &
         Bmatrix%T%old%sub, Bmatrix%T%old%diag,           &
         Bmatrix%T%old%sup,                               &
         Bmatrix%S%new%sub, Bmatrix%S%new%diag,           &
         Bmatrix%S%new%sup,                               &
         Bmatrix%S%old%sub, Bmatrix%S%old%diag,           &
         Bmatrix%S%old%sup,                               &
         stat=dallocstat)
    call dalloc_check(dallocstat, msg)
    msg = "U velocity component (cross-strait) RHS arrays"
    deallocate(U_RHS%diff_adv%new, U_RHS%diff_adv%old, &
         U_RHS%C_pg%new, U_RHS%C_pg%old,               &
         stat=dallocstat)
    call dalloc_check(dallocstat, msg)
    msg = "V velocity component (cross-strait) RHS arrays"
    deallocate(V_RHS%diff_adv%new, V_RHS%diff_adv%old, &
         V_RHS%C_pg%new, V_RHS%C_pg%old,               &
         stat=dallocstat)
    call dalloc_check(dallocstat, msg)
    msg = "Temperature RHS arrays"
    deallocate(T_RHS%diff_adv%new, T_RHS%diff_adv%old, &
         stat=dallocstat)
    call dalloc_check(dallocstat, msg)
    msg = "Salinity RHS arrays"
    deallocate(S_RHS%diff_adv%new, S_RHS%diff_adv%old, &
         stat=dallocstat)
    call dalloc_check(dallocstat, msg)
  end subroutine dalloc_phys_RHS_variables

end module physics_eqn_builder
       