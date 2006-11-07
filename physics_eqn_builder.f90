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
       nBmatrix, &  ! Tridiagonal matrix of diffusion coefficient values
       U_RHS,   &  ! U velocity component (cross-strait) right-hand side arrays
       V_RHS,   &  ! V velocity component (along-strait) right-hand side arrays
       T_RHS,   &  ! Temperature right-hand side arrays
       S_RHS,   &  ! Salinity right-hand side arrays
       ! Subroutines:
!!$       build_physics_equations, &
!!$       new_to_old_phys_RHS, new_to_old_phys_Bmatrix, &
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
          diff_adv  ! Diffusion/advection component of RHS
     real(kind=dp), dimension(:), allocatable :: &
          C_pg      ! Coriolis and baroclinic pressure grad component of RHS
  end type RHS

  ! Variable Declarations:
  !
  ! Public:
  type(diff_coeffs_matrix) :: &
       nBmatrix  ! Tridiagonal matrix of diffusion coefficient values
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
    allocate(nBmatrix%vel%new%sub(1:M), nBmatrix%vel%new%diag(1:M), &
         nBmatrix%vel%new%sup(1:M),                                &
         nBmatrix%vel%old%sub(1:M), nBmatrix%vel%old%diag(1:M),     &
         nBmatrix%vel%old%sup(1:M),                                &
         nBmatrix%T%new%sub(1:M), nBmatrix%T%new%diag(1:M),         &
         nBmatrix%T%new%sup(1:M),                                  &
         nBmatrix%T%old%sub(1:M), nBmatrix%T%old%diag(1:M),         &
         nBmatrix%T%old%sup(1:M),                                  &
         nBmatrix%S%new%sub(1:M), nBmatrix%S%new%diag(1:M),         &
         nBmatrix%S%new%sup(1:M),                                  &
         nBmatrix%S%old%sub(1:M), nBmatrix%S%old%diag(1:M),         &
         nBmatrix%S%old%sup(1:M),                                  &
         stat=allocstat)
    call alloc_check(allocstat, msg)
    msg = "U velocity component (cross-strait) RHS arrays"
    allocate(U_RHS%diff_adv%new(1:M), U_RHS%diff_adv%old(1:M), &
         U_RHS%C_pg(1:M),                                      &
         stat=allocstat)
    call alloc_check(allocstat, msg)
    msg = "V velocity component (cross-strait) RHS arrays"
    allocate(V_RHS%diff_adv%new(1:M), V_RHS%diff_adv%old(1:M), &
         V_RHS%C_pg(1:M),                                      &
         stat=allocstat)
    call alloc_check(allocstat, msg)
    msg = "Temperature RHS arrays"
    allocate(T_RHS%diff_adv%new(1:M), T_RHS%diff_adv%old(1:M), &
         T_RHS%C_pg(1:M),                                      &
         stat=allocstat)
    call alloc_check(allocstat, msg)
    msg = "Salinity RHS arrays"
    allocate(S_RHS%diff_adv%new(1:M), S_RHS%diff_adv%old(1:M), &
         S_RHS%C_pg(1:M),                                      &
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
    deallocate(nBmatrix%vel%new%sub, nBmatrix%vel%new%diag, &
         nBmatrix%vel%new%sup,                             &
         nBmatrix%vel%old%sub, nBmatrix%vel%old%diag,       &
         nBmatrix%vel%old%sup,                             &
         nBmatrix%T%new%sub, nBmatrix%T%new%diag,           &
         nBmatrix%T%new%sup,                               &
         nBmatrix%T%old%sub, nBmatrix%T%old%diag,           &
         nBmatrix%T%old%sup,                               &
         nBmatrix%S%new%sub, nBmatrix%S%new%diag,           &
         nBmatrix%S%new%sup,                               &
         nBmatrix%S%old%sub, nBmatrix%S%old%diag,           &
         nBmatrix%S%old%sup,                               &
         stat=dallocstat)
    call dalloc_check(dallocstat, msg)
    msg = "U velocity component (cross-strait) RHS arrays"
    deallocate(U_RHS%diff_adv%new, U_RHS%diff_adv%old, &
         U_RHS%C_pg,                                   &
         stat=dallocstat)
    call dalloc_check(dallocstat, msg)
    msg = "V velocity component (cross-strait) RHS arrays"
    deallocate(V_RHS%diff_adv%new, V_RHS%diff_adv%old, &
         V_RHS%C_pg,                                   &
         stat=dallocstat)
    call dalloc_check(dallocstat, msg)
    msg = "Temperature RHS arrays"
    deallocate(T_RHS%diff_adv%new, T_RHS%diff_adv%old, &
         T_RHS%C_pg,                                   &
         stat=dallocstat)
    call dalloc_check(dallocstat, msg)
    msg = "Salinity RHS arrays"
    deallocate(S_RHS%diff_adv%new, S_RHS%diff_adv%old, &
         S_RHS%C_pg,                                   &
         stat=dallocstat)
    call dalloc_check(dallocstat, msg)
  end subroutine dalloc_phys_RHS_variables

end module physics_eqn_builder
       
