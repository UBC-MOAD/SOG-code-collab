module IMEX_solver
  ! Type definitions, variable & parameter value declarations, and
  ! subroutines related to the IMEX semi-implicit PDE solver in the
  ! SOG code.
  !
  ! Public Subroutines:
  !
  !   init_IMEX_solver -- Initialize IMEX semi-implicit PDE solver variables.
  !
  !   dalloc_IMEX_variables -- Deallocate memory for arrays for tridiagonal
  !                            matrices and right-hand sides of
  !                            diffusion/advection PDEs.
  !
  !   solve_phys_eqns -- Solve the semi-implicit diffusion/advection
  !                      PDEs with Coriolis and baroclinic pressure
  !                      gradient terms for the physics quantities.
  !
  !   solve_bio_eqns -- Solve the semi-implicit diffusion/advection
  !                     PDEs for the biology quantities.

  ! Type definitions from other modules:
  use precision_defs, only: dp
  use numerics, only: tridiag  ! Tridiagonal matrix arrays type def'n
  
  implicit none

  private
  public :: &
       ! Parameter values:
       ! *** Temporary for refactoring
       a_IMEX1, &
       ! Variables:
       ! *** Temporary for refactoring
       Amatrix, Hvector, null_vector, &
       ! Types (as required by new pg compiler)
       RHS, & ! type for Hvector
       ! Subroutines:
       init_IMEX_solver, solve_phys_eqns, solve_bio_eqns, &
       dalloc_IMEX_variables, build_Amatrix, bio_Hvector, solve_tridiag

  ! Type Definitions:
  !
  ! Private to module:
  !
  ! LHS matrix types:
  type :: matrices
     type(tridiag) :: &
          vel, &  ! Velocity LHS matrix
          T,   &  ! Temperature LHS matrix
          S,   &  ! Salinity LHS matrix
          bio     ! Biology quantity LHS matrix
  end type matrices
  !
  ! Semi-implicit diffusion/advection PDE right-hand side (RHS) arrays
  type :: RHS
     real(kind=dp), dimension(:), allocatable :: &
          U,      &  ! U (cross-strait) velocity RHS array
          V,      &  ! V (along-strait) velocity RHS array
          T,      &  ! Temperature RHS array
          S,      &  ! Salinity RHS array
          Pmicro, &  ! Micro phytoplankton (diatoms) RHS array
          Pnano,  &  ! Nano phytoplankton (meso-rub) RHS array
          Ppico,  &  ! Pico phytoplankton (flagellates) RHS array
          Z,      &  ! Micro zooplankton RHS array
          NO,     &  ! Nitrate concentration RHS array
          NH,     &  ! Ammonium concentration RHS array
          Si,     &  ! Silicon concentration RHS array
          DIC,    &  ! Dissolved inorganic carbon RHS array
          Oxy,    &  ! Dissolved oxygen RHS array
          D_DOC,  &  ! Dissolved organic carbon detritus RHS array
          D_POC,  &  ! Particulate organic carbon detritus RHS array
          D_DON,  &  ! Dissolved organic nitrogen detritus RHS array
          D_PON,  &  ! Particulate organic nitrogen detritus RHS array
          D_refr, &  ! Refractory nitrogen detritus RHS array
          D_bSi      ! Biogenic silicon detritus RHS array
  end type RHS

  ! Parameter Value Declarations:
  !
  ! Private to module:
  !
  ! Implicit-explicit first order scheme:
  !   a_IMEX1 = 1.0 --> backward Euler in diffusion
  !           = 0.5 --> Crank-Nicolson in diffusion
  !           = 0.0 --> forward Euler in diffusion
  ! Both have forward Euler in explicit terms (i.e. Coriolis and
  ! boundary conditions)
  real(kind=dp), parameter :: a_IMEX1 = 1.0d0

  ! Variable Declarations:
  !
  ! Private to module:
  !
  type(matrices) :: &
       Amatrix  ! LHS matrices (A) for semi-implicit PDE matrix eqns
  type(RHS) :: &
       Hvector  ! RHS vectors (H) for semi-implicit PDE matrix eqns
  real(kind=dp), dimension(:), allocatable :: &
       null_vector  ! Vector of zeros; used for quantities that don't sink

contains

  subroutine init_IMEX_solver(M)
    ! Initialize IMEX semi-implicit PDE solver variables.
    implicit none
    ! Argument:
    integer :: M  ! Number of grid points

    ! Allocate memory for IMEX solver variables
    call alloc_IMEX_variables(M)
    ! Initialize the null vector of zeros used for non-sinking quantities
    null_vector = 0.
  end subroutine init_IMEX_solver


  subroutine solve_phys_eqns(M, day, time, U_old, V_old, T_old, S_old, &  ! in
       U_new, V_new, T_new, S_new)                                        ! out
    ! Solve the semi-implicit diffusion/advection PDEs with Coriolis
    ! and baroclinic pressure gradient terms for the physics quantities.

    ! Type definitions from other modules:
    use precision_defs, only: dp
    
    implicit none
    
    ! Arguments:
    integer, intent(in) :: &
         M  ! Number of grid points
    integer, intent(in) :: &
         day   ! Year-day of current time step
    real(kind=dp), intent(in) :: &
         time  ! Time of current time step [s since start of run]
    real(kind=dp), dimension(0:), intent(in) :: &
         U_old, &  ! U velocity component(cross-strait) profile; prev time step
         V_old, &  ! V velocity component(along-strait) profile; prev time step
         T_old, &  ! Temperature profile at previous time step
         S_old     ! Salinity profile at previous time step
    real(kind=dp), dimension(0:), intent(out) :: &
         U_new, &  ! U velocity component(cross-strait) profile; cur time step
         V_new, &  ! V velocity component(along-strait) profile; cur time step
         T_new, &  ! Temperature profile at current time step
         S_new     ! Salinity profile at current time step

    
    ! Build the RHS vectors (h) for the discretized semi-implicit PDE
    ! matrix equations Aq = h
    call build_phys_Hvectors(M, U_old, V_old, T_old, S_old)

    ! Build the LHS matrices (A) for the discretized semi-implicit PDE
    ! matrix equations Aq = h
    call build_phys_Amatrices()

    ! Solve the discretized semi-implicit PDE matrix equations Aq = h
    call solve_phys_tridiags(M, U_new, V_new, T_new, S_new)
  end subroutine solve_phys_eqns


  subroutine build_phys_Hvectors(M, U_old, V_old, T_old, S_old)
    ! Build the RHS vectors (h) for the discretized semi-implicit PDE
    ! matrix equations Aq = h

    ! Elements from other modules:
    !
    ! Type Definitions:
    use precision_defs, only: dp
    ! Variables:
    use physics_eqn_builder, only: &
         Bmatrix, &  ! Precursor diffusion coefficient matrices
         U_RHS, &     ! U velocity component (cross-strait) RHS arrays
         V_RHS, &     ! V velocity component (along-strait) RHS arrays
         T_RHS, &     ! Temperature RHS arrays
         S_RHS        ! Salinity RHS arrays

    implicit none
    
    ! Arguments:
    integer, intent(in) :: &
         M  ! Number of grid points
    real(kind=dp), dimension(0:), intent(in) :: &
         U_old, &  ! U velocity component(cross-strait) profile; prev time step
         V_old, &  ! V velocity component(along-strait) profile; prev time step
         T_old, &  ! Temperature profile at previous time step
         S_old     ! Salinity profile at previous time step

    ! U velocity component
    call phys_Hvector(M, U_old, U_RHS%diff_adv%new, U_RHS%diff_adv%old, & ! in
         U_RHS%C_pg%new, U_RHS%C_pg%old, Bmatrix%vel%old,              & ! in
         Hvector%U)                                                       ! out
    ! V velocity component
    call phys_Hvector(M, V_old, V_RHS%diff_adv%new, V_RHS%diff_adv%old, & ! in
         V_RHS%C_pg%new, V_RHS%C_pg%old, Bmatrix%vel%old,              & ! in
         Hvector%V)                                                       ! out
    ! Temperature
    ! Use null_vector because no Coriolis or pressure gradients terms
    call phys_Hvector(M, T_old, T_RHS%diff_adv%new, T_RHS%diff_adv%old, & ! in
         null_vector, null_vector, Bmatrix%T%old,                      & ! in
         Hvector%T)                                                       ! out
    ! Salinity
    ! Use null_vector because no Coriolis or pressure gradients terms
    call phys_Hvector(M, S_old, S_RHS%diff_adv%new, S_RHS%diff_adv%old, & ! in
         null_vector, null_vector, Bmatrix%S%old,                      & ! in
         Hvector%S)                                                       ! out
  end subroutine build_phys_Hvectors
  

  subroutine phys_Hvector(M, qty_old, flux, flux_old, C_pg, &
       C_pg_old, Bmatrix, Hvector)
    ! Calculate RHS (Hvector) of semi-implicit PDE for a physic quantity. 
    ! It includes the value of the quantity from the previous time step, and
    ! current time step non-local flux terms, bottom and surface 
    ! flux terms, and Coriolis and pressure gradient flux terms.  
    ! Depending on the IMEX scheme being used to solve the PDE values of
    ! the quanties above, and of the diffusion coefficient matrix from
    ! the previous time step may also be blended in.

    ! Elements from other modules:
    !
    ! Type Definitions:
    use precision_defs, only: dp
    use numerics, only: tridiag  ! Tridiagonal matrix arrays type def'n
    
    implicit none
    
    ! Arguments:
    integer :: &
         M  ! Number of grid points
    real(kind=dp), dimension(0:) :: &
         qty_old  ! Profile of the quantity at the previous time step
    real(kind=dp), dimension(1:), intent(in):: &
         flux,     &  ! Non-local flux term at current time step
         flux_old, &  ! Non-local flux term at previous time step
         C_pg,     &  ! Coriolis and pressure gradient flux term (cur)
         C_pg_old     ! Coriolis and pressure gradient flux term (prev)
    type(tridiag), intent(in) :: &
         Bmatrix  ! Diffussion coefficients matrix at previous time step
    real(kind=dp), dimension(1:), intent(out) :: &
         Hvector  ! RHS of discretized semi-implicit PDE matrix equation

    Hvector = qty_old(1:M)                                      &
         + (1.0 - a_IMEX1) * (flux_old + C_pg_old               &
                            + Bmatrix%sub * qty_old(0:M-1)  &
                            + Bmatrix%diag * qty_old(1:M)   &
                            + Bmatrix%sup * qty_old(2:M+1)) &
         + a_IMEX1 * (flux + C_pg)
  end subroutine phys_Hvector


  subroutine build_phys_Amatrices()
    ! Build the LHS matrices (A) for the discretized semi-implicit PDE
    ! matrix equations Aq = h for physics quantities.

    ! Elements from other modules:
    !
    ! Variables:
    use physics_eqn_builder, only: &
         Bmatrix  ! Precursor diffusion coefficient matrices

    implicit none

    ! Velocity components
    call build_Amatrix(Bmatrix%vel%new, Amatrix%vel)
    ! Temperature
    call build_Amatrix(Bmatrix%T%new, Amatrix%T)
    ! Salinity
    call build_Amatrix(Bmatrix%S%new, Amatrix%S)
  end subroutine build_phys_Amatrices


  subroutine solve_phys_tridiags(M, U_new, V_new, T_new, S_new)
    ! Solve the discretized semi-implicit PDE matrix equations Aq = h
    ! for the physics quantities.
    
    ! Elements from other modules:
    !
    ! Type Definitions:
    use precision_defs, only: dp
    
    implicit none
    
    ! Arguments:
    integer, intent(in) :: &
         M  ! Number of grid points
    real(kind=dp), dimension(0:), intent(inout) :: &
         U_new, &  ! U velocity component(cross-strait) profile; cur time step
         V_new, &  ! V velocity component(along-strait) profile; cur time step
         T_new, &  ! Temperature profile at current time step
         S_new     ! Salinity profile at current time step
    
    ! U velocity component
    call solve_tridiag(Amatrix%vel, Hvector%U, U_new(1:M))
    ! V velocity component
    call solve_tridiag(Amatrix%vel, Hvector%V, V_new(1:M))
    ! Temperature
    call solve_tridiag(Amatrix%T, Hvector%T, T_new(1:M))
    ! Salinity
    call solve_tridiag(Amatrix%S, Hvector%S, S_new(1:M))
  end subroutine solve_phys_tridiags


  subroutine solve_bio_eqns(M, Pmicro, Pnano, Ppico, Z, NO, NH, Si, &
       D_DOC, D_POC, D_DON, D_PON, D_refr, D_bSi, day, time)
    ! Solve the semi-implicit diffusion/advection PDEs for the
    ! biology quantities.

    ! Elements from other modules:
    !
    ! Type Definitions:
    use precision_defs, only: dp
    ! Variables:
    use biology_eqn_builder, only: &
         Bmatrix  ! Precursor diffusion coefficient matrices
    
    implicit none
    
    ! Arguments:
    integer, intent(in) :: &
         M  ! Number of grid points
    real(kind=dp), dimension(0:), intent(inout) :: &
         Pmicro, &  ! Micro phytoplankton
         Pnano,  &  ! Nano phytoplankton
         Ppico,  &  ! Pico phytoplankton
         Z,      &  ! Micro Zooplankton
         NO,     &  ! Nitrate
         NH,     &  ! Ammonium
         Si,     &  ! Silicon
         D_DOC,  &  ! Dissolved organic carbon detritus profile
         D_POC,  &  ! Particulate organic carbon detritus profile
         D_DON,  &  ! Dissolved organic nitrogen detritus profile
         D_PON,  &  ! Particulate organic nitrogen detritus profile
         D_refr, &  ! Refractory nitrogen detritus profile
         D_bSi      ! Biogenic silicon detritus profile
    integer, intent(in) :: &
         day   ! Year-day of current time step
    real(kind=dp), intent(in) :: &
         time  ! Time of current time step [s since start of run]

    
    ! Build the RHS vectors (h) for the discretized semi-implicit PDE
    ! matrix equations Aq = h
    call build_bio_Hvectors(M, Pmicro, Pnano, Ppico, Z, NO, NH, Si, &
       D_DOC, D_POC, D_DON, D_PON, D_refr, D_bSi)

    ! Build the LHS matrix (A) for the discretized semi-implicit PDE
    ! matrix equations Aq = h
    call build_Amatrix(Bmatrix%bio%new, Amatrix%bio)

    ! Solve the discretized semi-implicit PDE matrix equations Aq = h
    call solve_bio_tridiags(M, Pmicro, Pnano, Ppico, Z, NO, NH, Si, &
       D_DOC, D_POC, D_DON, D_PON, D_refr, D_bSi)

    ! Check for negative values in results, and print with a warning
    ! message if any are found
    call check_bio_negatives(Pmicro, Pnano, Ppico, Z, NO, NH, Si, &
       D_DOC, D_POC, D_DON, D_PON, D_refr, D_bSi, day, time, &
       fatal=.false.)
  end subroutine solve_bio_eqns


  subroutine build_bio_Hvectors(M, Pmicro, Pnano, Ppico, Z, NO, NH, Si, &
       D_DOC, D_POC, D_DON, D_PON, D_refr, D_bSi)
    ! Build the RHS vectors (h) for the discretized semi-implicit PDE
    ! matrix equations Aq = h

    ! Elements from other modules:
    !
    ! Type Definitions:
    use precision_defs, only: dp
    ! Variables:
    use biology_eqn_builder, only: &
         Bmatrix,    &  ! Precursor diffusion coefficient matrices
         Pmicro_RHS, &  ! Micro phytoplankton (diatoms) RHS arrays
         Pnano_RHS,  &  ! Nano phytoplankton (mesorub) RHS arrays         
         Ppico_RHS,  &  ! Pico phytoplankton (flagellates) RHS arrays
         Z_RHS,      &  ! Zooplankton (micro) RHS arrays
         NO_RHS,     &  ! Nitrate concentration RHS arrays
         NH_RHS,     &  ! Ammonium concentration RHS arrays
         Si_RHS,     &  ! Silicon concentration RHS arrays
         D_DOC_RHS,  &  ! Dissolved organic carbon detritus RHS arrays
         D_POC_RHS,  &  ! Particulate organic carbon detritus RHS arrays
         D_DON_RHS,  &  ! Dissolved organic nitrogen detritus RHS arrays
         D_PON_RHS,  &  ! Particulate organic nitrogen detritus RHS arrays
         D_refr_RHS, &  ! Refractory nitrogen detritus RHS arrays
         D_bSi_RHS      ! Biogenic silicon detritus RHS arrays
    
    implicit none
    
    ! Arguments:
    integer, intent(in) :: &
         M  ! Number of grid points
    real(kind=dp), dimension(0:), intent(inout) :: &
         Pmicro, &  ! Micro phytoplankton
         Pnano,  &  ! Nano phytoplankton
         Ppico,  &  ! Pico phytoplankton
         Z,      &  ! Microzooplankton
         NO,     &  ! Nitrate
         NH,     &  ! Ammonium
         Si,     &  ! Silicon
         D_DOC,  &  ! Dissolved organic carbon detritus profile
         D_POC,  &  ! Particulate organic carbon detritus profile
         D_DON,  &  ! Dissolved organic nitrogen detritus profile
         D_PON,  &  ! Particulate organic nitrogen detritus profile
         D_refr, &  ! Refractory nitrogen detritus profile
         D_bSi      ! Biogenic silicon detritus profile

    ! Micro phytoplankton (diatoms)
    call bio_Hvector(M, Pmicro, Pmicro_RHS%diff_adv%new, &
         Pmicro_RHS%diff_adv%old, Pmicro_RHS%bio, Pmicro_RHS%sink, &
         Bmatrix%bio%old, &
         Hvector%Pmicro)
    ! Nano phytoplankton (meso-rub); non-sinking
    call bio_Hvector(M, Pnano, Pnano_RHS%diff_adv%new, &
         Pnano_RHS%diff_adv%old, Pnano_RHS%bio, null_vector, &
         Bmatrix%bio%old, &
         Hvector%Pnano)
    ! Pico phytoplankton (flagellates); non-sinking
    call bio_Hvector(M, Ppico, Ppico_RHS%diff_adv%new, &
         Ppico_RHS%diff_adv%old, Ppico_RHS%bio, null_vector, &
         Bmatrix%bio%old, &
         Hvector%Ppico)
    ! Microzooplankton: non-sinking
    call bio_Hvector(M, Z, Z_RHS%diff_adv%new, &
         Z_RHS%diff_adv%old, Z_RHS%bio, null_vector, &
         Bmatrix%bio%old, &
         Hvector%Z)
    ! Nitrate; non-sinking
    call bio_Hvector(M, NO, NO_RHS%diff_adv%new, &
         NO_RHS%diff_adv%old, NO_RHS%bio, null_vector, &
         Bmatrix%bio%old, &
         Hvector%NO)
    ! Ammonium; non-sinking
    call bio_Hvector(M, NH, NH_RHS%diff_adv%new, &
         NH_RHS%diff_adv%old, NH_RHS%bio, null_vector, &
         Bmatrix%bio%old, &
         Hvector%NH)
    ! Silicon; non-sinking
    call bio_Hvector(M, Si, Si_RHS%diff_adv%new, &
         Si_RHS%diff_adv%old, Si_RHS%bio, null_vector, &
         Bmatrix%bio%old, &
         Hvector%Si)
    ! Dissolved organic carbon detritus; non-sinking
    call bio_Hvector(M, D_DOC, D_DOC_RHS%diff_adv%new, &
         D_DOC_RHS%diff_adv%old, D_DOC_RHS%bio, null_vector, &
         Bmatrix%bio%old, &
         Hvector%D_DOC)
!----- STILL NEED TO INVESTIGATE!!! ---------
    ! Particulate organic carbon detritus
    call bio_Hvector(M, D_POC, D_POC_RHS%diff_adv%new, &
         D_POC_RHS%diff_adv%old, D_POC_RHS%bio, D_POC_RHS%sink, & ! POC sinking
         Bmatrix%bio%old, &
         Hvector%D_POC)
!--------------------------------------------
    ! Dissolved organic nitrogen detritus; non-sinking
    call bio_Hvector(M, D_DON, D_DON_RHS%diff_adv%new, &
         D_DON_RHS%diff_adv%old, D_DON_RHS%bio, null_vector, &
         Bmatrix%bio%old, &
         Hvector%D_DON)
    ! Particulate organic nitrogen detritus
    call bio_Hvector(M, D_PON, D_PON_RHS%diff_adv%new, &
         D_PON_RHS%diff_adv%old, D_PON_RHS%bio, D_PON_RHS%sink, &
         Bmatrix%bio%old, &
         Hvector%D_PON)
    ! Refractory nitrogen detritus
    call bio_Hvector(M, D_refr, D_refr_RHS%diff_adv%new, &
         D_refr_RHS%diff_adv%old, D_refr_RHS%bio, D_refr_RHS%sink, &
         Bmatrix%bio%old, &
         Hvector%D_refr)
    ! Biogenic silicon detritus
    call bio_Hvector(M, D_bSi, D_bSi_RHS%diff_adv%new, &
         D_bSi_RHS%diff_adv%old, D_bSi_RHS%bio, D_bSi_RHS%sink, &
         Bmatrix%bio%old, &
         Hvector%D_bSi)
  end subroutine build_bio_Hvectors

  
  subroutine bio_Hvector(M, qty_old, diff_adv, diff_adv_old, bio, sink, &
       Bmatrix, Hvector)
    ! Calculate RHS (Hvector) of semi-implicit PDE for a biology quantity. 
    ! It includes the value of the quantity from the previous time step, and
    ! current time step upwelling flux terms, bottom and surface 
    ! flux terms, growth - mortality terms, and sinking terms.  
    ! Depending on the IMEX scheme being used to solve the PDE values of
    ! the quanties above, and of the diffusion coefficient matrix from
    ! the previous time step may also be blended in.

    ! Type definitions from other modules:
    use precision_defs, only: dp
    use numerics, only: tridiag  ! Tridiagonal matrix arrays type def'n
    
    implicit none
    
    ! Arguments:
    integer :: &
         M  ! Number of grid points
    real(kind=dp), dimension(0:) :: &
         qty_old  ! Profile of the quantity at the previous time step
    real(kind=dp), dimension(1:), intent(in):: &
         diff_adv,     &  ! Diffusion/advection term (current time step)
         diff_adv_old, &  ! Diffusion/advection term (previous time step)
         bio,          &  ! Growth - mortality term
         sink             ! Sinking term
    type(tridiag), intent(in) :: &
         Bmatrix  ! Diffussion coefficients matrix at previous time step
    real(kind=dp), dimension(1:), intent(out) :: &
         Hvector  ! RHS of discretized semi-implicit PDE matrix equation

    Hvector = qty_old(1:M) + bio + sink                           &
         + (1.0 - a_IMEX1) * (diff_adv_old                        &
                              + Bmatrix%sub * qty_old(0:M-1)  &
                              + Bmatrix%diag * qty_old(1:M)   &
                              + Bmatrix%sup * qty_old(2:M+1)) &
         + a_IMEX1 * diff_adv
    ! Note that abs(x) < epsilon(x) is a real-number-robust test for x
    ! == 0.
    where (abs(Hvector) < epsilon(Hvector)) Hvector = 0.
  end subroutine bio_Hvector


  subroutine solve_bio_tridiags(M, Pmicro, Pnano, Ppico, Z, NO, NH, Si, &
       D_DOC, D_POC, D_DON, D_PON, D_refr, D_bSi)
    ! Solve the discretized semi-implicit PDE matrix equations Aq = h.

    ! Type definitions from other modules:
    use precision_defs, only: dp

    implicit none

    ! Arguments:
    integer, intent(in) :: &
         M  ! Number of grid points
    real(kind=dp), dimension(0:), intent(inout) :: &
         Pmicro, &  ! Micro phytoplankton
         Pnano,  &  ! Nano phytoplankton
         Ppico,  &  ! Pico phytoplankton
         Z,      &  ! Micro zooplankton
         NO,     &  ! Nitrate
         NH,     &  ! Ammonium
         Si,     &  ! Silicon
         D_DOC,  &  ! Dissolved organic carbon detritus profile
         D_POC,  &  ! Particulate organic carbon detritus profile
         D_DON,  &  ! Dissolved organic nitrogen detritus profile
         D_PON,  &  ! Particulate organic nitrogen detritus profile
         D_refr, &  ! Refractory nitrogen detritus profile
         D_bSi      ! Biogenic silicon detritus profile

    ! Micro phytoplankton (diatoms)
    call solve_tridiag(Amatrix%bio, Hvector%Pmicro, Pmicro(1:M))
    ! Nano phytoplankton (meso-rub)
    call solve_tridiag(Amatrix%bio, Hvector%Pnano, Pnano(1:M))
    ! Pico phytoplankton (flagellates)
    call solve_tridiag(Amatrix%bio, Hvector%Ppico, Ppico(1:M))
    ! Microzooplantkon
    call solve_tridiag(Amatrix%bio, Hvector%Z, Z(1:M))
    ! Nitrate
    call solve_tridiag(Amatrix%bio, Hvector%NO, NO(1:M))
    ! Ammonium
    call solve_tridiag(Amatrix%bio, Hvector%NH, NH(1:M))
    ! Silicon
    call solve_tridiag(Amatrix%bio, Hvector%Si, Si(1:M))
    ! Dissolved organic carbon detritus
    call solve_tridiag(Amatrix%bio, Hvector%D_DOC, D_DOC(1:M))
    ! Particulate organic carbon detritus
    call solve_tridiag(Amatrix%bio, Hvector%D_POC, D_POC(1:M))
    ! Dissolved organic nitrogen detritus
    call solve_tridiag(Amatrix%bio, Hvector%D_DON, D_DON(1:M))
    ! Particulate organic nitrogen detritus
    call solve_tridiag(Amatrix%bio, Hvector%D_PON, D_PON(1:M))
    ! Refractory nitrogen detritus
    call solve_tridiag(Amatrix%bio, Hvector%D_refr, D_refr(1:M))
    ! Biogenic silicon detritus
    call solve_tridiag(Amatrix%bio, Hvector%D_bSi, D_bSi(1:M))
  end subroutine solve_bio_tridiags


  subroutine check_bio_negatives(Pmicro, Pnano, Ppico, Z, NO, NH, Si, &
       D_DOC, D_POC, D_DON, D_PON, D_refr, D_bSi, day, time, fatal)
    ! Check for negative values in results.

    ! Elements from other modules:
    !
    ! Type Definitions:
    use precision_defs, only: dp
    ! Subroutines:
    use numerics, only: check_negative
    
    implicit none
    
    ! Arguments:
    real(kind=dp), dimension(0:), intent(inout) :: &
         Pmicro, &  ! Micro phytoplankton
         Pnano,  &  ! Nano phytoplankton
         Ppico,  &  ! Pico phytoplankton
         Z,      &  ! Microzooplankton
         NO,     &  ! Nitrate
         NH,     &  ! Ammonium
         Si,     &  ! Silicon
         D_DOC,  &  ! Dissolved organic carbon detritus profile
         D_POC,  &  ! Particulate organic carbon detritus profile
         D_DON,  &  ! Dissolved organic nitrogen detritus profile
         D_PON,  &  ! Particulate organic nitrogen detritus profile
         D_refr, &  ! Refractory nitrogen detritus profile
         D_bSi      ! Biogenic silicon detritus profile
    integer, intent(in) :: &
         day   ! Year-day of current time step
    real(kind=dp), intent(in) :: &
         time  ! Time of current time step [s since start of run]
    logical, intent(in), optional :: &
         fatal

    ! Micro phytoplankton (diatoms)
    call check_negative(0, Pmicro, "P%micro after solve_tridiag()", &
         day, time, fatal)
    ! Nano phytoplankton (meso-rub)
    call check_negative(0, Pnano, "P%nano after solve_tridiag()", &
         day, time, fatal)
    ! Nano phytoplankton (flagellates)
    call check_negative(0, Ppico, "P%pico after solve_tridiag()", &
         day, time, fatal)
    ! Micro zooplankton 
    call check_negative (0, Z, "Z after solve_tridiag()", &
         day, time, fatal)
    ! Nitrate
    call check_negative(0, NO, "N%O after solve_tridiag()", &
         day, time, fatal)
    ! Ammonium
    call check_negative(0, NH, "N%H after solve_tridiag()", &
         day, time, fatal)
    ! Silicon
    call check_negative(0, Si, "Si after solve_tridiag()", &
         day, time, fatal)
    ! Dissolved organic carbon detritus
    call check_negative(0, D_DOC, "D%DOC after solve_tridiag()", &
         day, time, fatal)
    ! Particulate organic carbon detritus
    call check_negative(0, D_POC, "D%POC after solve_tridiag()", &
         day, time, fatal)
    ! Dissolved organic nitrogen detritus
    call check_negative(0, D_DON, "D%DON after solve_tridiag()", &
         day, time, fatal)
    ! Particulate organic nitrogen detritus
    call check_negative(0, D_PON, "D%PON after solve_tridiag()", &
         day, time, fatal)
    ! Refractory nitrogen detritus
    call check_negative(0, D_refr, "D%refr after solve_tridiag()", &
         day, time, fatal)
    ! Biogenic silicon detritus
    call check_negative(0, D_bSi, "D%bSi after solve_tridiag()", &
         day, time, fatal)
  end subroutine check_bio_negatives


  subroutine build_Amatrix(Bmatrix, Amatrix)
    ! Return the LHS matrix (A) for semi-implicit PDE matrix equation
    ! from the matrix of diffusion coefficients (B) for the quantity
    ! being solved.

    ! Type definitions from other modules:
    use precision_defs, only: dp
    use numerics, only: tridiag  ! Tridiagonal matrix arrays type def'n
    
    implicit none
    
    ! Arguments:
    type(tridiag), intent(in) :: &
         Bmatrix  ! Diffusion coefficient matrix sub-diagonal (prev timestep)
    type(tridiag), intent(inout) :: &  ! inout 'cause not all components set
         Amatrix  ! LHS matrix for discretized semi-implicit PDE matrix eqn

    Amatrix%sub = -a_IMEX1 * Bmatrix%sub
    Amatrix%diag = 1.0 - a_IMEX1 * Bmatrix%diag
    Amatrix%sup = -a_IMEX1 * Bmatrix%sup
  end subroutine build_Amatrix


  subroutine solve_tridiag(A, h, q)
    ! Solve the matrix equation Aq = h, where A is a tridiagonal matrix,
    ! h is the RHS side vector, and q is the vector of the quantity to
    ! be solved for.

    ! Elements from other modules:
    !
    ! Type Definitions:
    use precision_defs, only: dp
    use numerics, only: tridiag  ! Tridiagonal matrix arrays type def'n
    ! Variables:
    use io_unit_defs, only: stdout
    
    implicit none
    
    ! Arguments:
    type(tridiag), intent(in) :: &
         A  ! Tridiagonal matrix in matrix eq'n Aq = h
    real(kind=dp), dimension(1:), intent(in) :: &
         h     ! Right-hand-side vector of matrix eq'm Aq = h
    real(kind=dp), dimension(1:), intent(out) :: &
         q  ! Vector of quantity to be solved for in matrix eq'n Aq = h
    ! Local variables:
    real(kind=dp), dimension(size(h)) :: gamma
    real(kind=dp) :: beta
    integer :: j

    ! Confirm that the matrix equation is properly expressed.  Note
    ! that abs(x) < epsilon(x) is a real-number-robust test for x == 0.
    if (abs(A%diag(1)) < epsilon(A%diag(1))) then
       write(stdout, *) "tridiag: Malformed matrix, Amatrix(1) = 0"
       call exit(1)
    endif
    ! Decomposition and forward substitution
    beta = A%diag(1)
    q(1) = h(1) / beta
    do j = 2, size(h)
       gamma(j) = A%sup(j-1) / beta
       beta = A%diag(j) - A%sub(J) * gamma(j)
       if (abs(beta) < epsilon(beta)) then
          write(stdout, *) "tridiag: Solver failed, beta = 0 at j = ", j
          call exit(1)
       endif
       q(j) = (h(j) - A%sub(j) * q(j-1)) / beta
    enddo
    ! Back substitution
    do j = size(h) - 1, 1, -1
       q(j) = q(j) - gamma(j+1) * q(j+1)
    enddo
  end subroutine solve_tridiag


  subroutine alloc_IMEX_variables(M)
    ! Allocate memory for arrays for tridiagonal matrices and
    ! right-hand sides of diffusion/advection PDEs.
    use malloc, only: alloc_check
    implicit none
    ! Argument:
    integer, intent(in) :: M  ! Number of grid points
    ! Local variables:
    integer           :: allocstat  ! Allocation return status
    character(len=80) :: msg        ! Allocation failure message prefix

    msg = "LHS matrices (A) for semi-implicit PDE matrix eqns"
    allocate(Amatrix%vel%sub(1:M), Amatrix%vel%diag(1:M), &
         Amatrix%vel%sup(1:M),                             &
         Amatrix%T%sub(1:M), Amatrix%T%diag(1:M),         &
         Amatrix%T%sup(1:M),                               &
         Amatrix%S%sub(1:M), Amatrix%S%diag(1:M),         &
         Amatrix%S%sup(1:M),                               &
         Amatrix%bio%sub(1:M), Amatrix%bio%diag(1:M),     &
         Amatrix%bio%sup(1:M),                             &
         stat=allocstat)
    call alloc_check(allocstat, msg)
    msg = "RHS vectors (H) for semi-implicit PDE matrix eqns"
    allocate(Hvector%U(1:M), Hvector%V(1:M), Hvector%T(1:M),          &
         Hvector%S(1:M), Hvector%Pmicro(1:M), Hvector%Pnano(1:M),     &
         Hvector%Ppico(1:M), Hvector%Z(1:M), Hvector%NO(1:M),         &
         Hvector%NH(1:M), Hvector%Si(1:M), Hvector%DIC(1:M),          &
         Hvector%Oxy(1:M), Hvector%D_DOC(1:M), Hvector%D_POC(1:M),    &
         Hvector%D_DON(1:M), Hvector%D_PON(1:M), Hvector%D_refr(1:M), &
         Hvector%D_bSi(1:M),                                          &
         stat=allocstat)
    call alloc_check(allocstat, msg)
    msg = "IMEX solver null vector"
    allocate(null_vector(1:M), &
         stat=allocstat)
    call alloc_check(allocstat, msg)
  end subroutine alloc_IMEX_variables


  subroutine dalloc_IMEX_variables
    ! Deallocate memory for arrays for tridiagonal matrices and
    ! right-hand sides of diffusion/advection PDEs.
    use malloc, only: dalloc_check
    implicit none
    ! Local variables:
    integer           :: dallocstat  ! Deallocation return status
    character(len=80) :: msg         ! Deallocation failure message prefix

    msg = "LHS matrices (A) for semi-implicit PDE matrix eqns"
    deallocate(Amatrix%vel%sub, Amatrix%vel%diag, &
         Amatrix%vel%sup,                          &
         Amatrix%T%sub, Amatrix%T%diag,           &
         Amatrix%T%sup,                            &
         Amatrix%S%sub, Amatrix%S%diag,           &
         Amatrix%S%sup,                            &
         Amatrix%bio%sub, Amatrix%bio%diag,       &
         Amatrix%bio%sup,                          &
         stat=dallocstat)
    call dalloc_check(dallocstat, msg)
    msg = "RHS vectors (H) for semi-implicit PDE matrix eqns"
    deallocate(Hvector%U, Hvector%V, Hvector%T,        &
         Hvector%S, Hvector%Pmicro, Hvector%Pnano,     &
         Hvector%Ppico, Hvector%Z, Hvector%NO,         &
         Hvector%NH, Hvector%Si, Hvector%DIC,          &
         Hvector%Oxy, Hvector%D_DOC, Hvector%D_POC,    &
         Hvector%D_DON, Hvector%D_PON, Hvector%D_refr, &
         Hvector%D_bSi,                                &
         stat=dallocstat)
    call dalloc_check(dallocstat, msg)
    msg = "IMEX solver null vector"
    deallocate(null_vector, &
         stat=dallocstat)
    call dalloc_check(dallocstat, msg)
  end subroutine dalloc_IMEX_variables

end module IMEX_solver
