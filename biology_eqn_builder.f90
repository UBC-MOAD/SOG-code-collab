! $Id$
! $Source$

module biology_eqn_builder
  ! Type definitions, variable declarations, and subroutines related
  ! to building the semi-implicit diffusion/advection PDEs for the
  ! biology quantities in the SOG code.
  !
  ! Public Variables:
  !
  !   Bmatrix -- Tridiagonal matrix of diffusion coefficient values;
  !              precursor to the LHS A matrix in the Aq = h matrix
  !              equation (see Large, et al (1994), App. D, pg 398)
  !
  !   Right-hand side vector arrays; precursor to the RHS h vector in
  !   the Aq = h matrix equation (see Large, et al (1994), App. D, pg 398):
  !
  !   Pmicro_RHS -- Micro phytoplankton (diatoms) right-hand side arrays
  !
  !   Pnano_RHS -- Nano phytoplankton (flagellates) right-hand side arrays
  !
  !   Z_RHS -- Micro zooplankton right-hand side arrays
  !
  !   NO_RHS -- Nitrate concentration right-hand side arrays
  !
  !   NH_RHS -- Ammonium concentration right-hand side arrays
  !
  !   Si_RHS -- Silicon concentration right-hand side arrays
  !
  !   D_DON_RHS -- Dissolved organic nitrogen detritus right-hand side arrays
  !
  !   D_PON_RHS -- Particulate organic nitro detritus right-hand side arrays
  !
  !   D_refr_RHS -- Refractory nitrogen detritus right-hand side arrays
  !
  !   D_bSi_RHS -- Biogenic silicon detritus right-hand side arrays
  !
  ! Public Subroutines:
  !
  !   read_sink_params -- Read the sinking rate parameter values from
  !                       the infile.
  !
  !   build_biology_equations -- Build the right-hand side (RHS) arrays
  !                              for the diffusion/advection equations for
  !                              the biology quantities.
  !
  !   new_to_old_bio_RHS -- Copy %new component of the biology *_RHS%diff_adv
  !                         arrays to the %old component for use by the
  !                         IMEX semi-impllicit PDE solver.
  !
  !   new_to_old_bio_Bmatrix -- Copy %new component of the Bmatrix%bio
  !                             arrays to the %old component for use by the
  !                             IMEX semi-impllicit PDE solver.
  !
  !   alloc_bio_RHS_variables -- Allocate memory for biology RHS arrays
  !
  !   dalloc_bio_RHS_variables -- Deallocate memory for biology RHS arrays

  use precision_defs, only: dp
  use numerics, only: tridiag  ! Tridiagonal matrix arrays type def'n
  implicit none

  private
  public :: &
       ! Variables:
       Bmatrix,         &  ! Tridiagonal matrix of diffusion coefficient values
       Pmicro_RHS,      &  ! Micro phytoplankton (diatoms) RHS arrays
       Pnano_RHS,       &  ! Nano phytoplankton (flagellates) RHS arrays
       Z_RHS,           &  ! Zooplankton (micro) RHS arrays
       NO_RHS,          &  ! Nitrate concentration RHS arrays
       NH_RHS,          &  ! Ammonium concentration RHS arrays
       Si_RHS,          &  ! Silicon concentration RHS arrays
       D_DON_RHS,       &  ! Dissolved organic nitrogen detritus RHS arrays
       D_PON_RHS,       &  ! Particulate organic nitrogen detritus RHS arrays
       D_refr_RHS,      &  ! Refractory nitrogen detritus RHS arrays
       D_bSi_RHS,       &  ! Biogenic silicon detritus RHS arrays
       ! Subroutines:
       read_sink_params, build_biology_equations, &
       new_to_old_bio_RHS, new_to_old_bio_Bmatrix, &
       alloc_bio_RHS_variables, dalloc_bio_RHS_variables

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
          bio  ! Biology quantities diffusion coefficients matrix
  end type diff_coeffs_matrix
  !
  ! New/old array components:
  type :: new_old_arrays
     real(kind=dp), dimension(:), allocatable :: &
          new, &  ! Current time step values
          old     ! Previous time step values
  end type new_old_arrays
  !
  ! Semi-implicit diffusion/advection equation right-hand side (RHS) arrays
  type :: RHS
     type(new_old_arrays) :: &
          diff_adv  ! Diffusion/advection component of RHS
     real(kind=dp), dimension(:), allocatable :: &
          bio, &  ! Biology (growth - mortality) component of RHS
          sink    ! Sinking component of RHS
  end type RHS
  !
  ! Sinking velocities of biology quantities
  type :: sink_vels
     real(kind=dp) :: &
          Pmicro_min, &  ! Minimum sinking velocity of micro phytos [m/s]
          Pmicro_max, &  ! Maximum sinking velocity of micro phytos [m/s]
          D_PON,      &  ! Sinking velocity of PON detritus [m/s]
          D_refr,     &  ! Sinking velocity of refractorydetritus [m/s]
          D_bSi          ! Sinking velocity of biogenic silicon detritus [m/s]
  end type sink_vels

  ! Variable Declarations:
  !
  ! Public:
  type(diff_coeffs_matrix) :: &
       Bmatrix  ! Tridiagonal matrix of diffusion coefficient values
                ! Components:
                !   Bmatrix%bio
                !              %new
                !              %old
                !                  %sub
                !                  %diag
                !                  %sup
  type(RHS) :: &
       Pmicro_RHS, &  ! Micro phytoplankton (diatoms) RHS arrays
       Pnano_RHS,  &  ! Nano phytoplankton (flagellates) RHS arrays
       Z_RHS,      &  ! micro Zooplankton RHS arrays
       NO_RHS,     &  ! Nitrate concentration RHS arrays
       NH_RHS,     &  ! Ammonium concentration RHS arrays
       Si_RHS,     &  ! Silicon concentration RHS arrays
       D_DON_RHS,  &  ! Dissolved organic nitrogen detritus RHS arrays
       D_PON_RHS,  &  ! Particulate organic nitro detritus RHS arrays
       D_refr_RHS, &  ! Refractory nitrogen detritus RHS arrays
       D_bSi_RHS      ! Biogenic silicon detritus RHS arrays
  !
  ! Private to module:
  !
  ! Sinking velocities:
  type(sink_vels) :: &
       w_sink  ! Sinking velocities of biology quantities [m/s]

contains

  subroutine read_sink_params()
    ! Read the sinking rate parameter values from the infile.
    use input_processor, only: getpard
    implicit none

    ! Micro phytoplankton (diatoms) minimum sinking rate [m/d]
    w_sink%Pmicro_min = getpard("Micro min sink rate") / 86400.
    ! Micro phytoplankton (diatoms) maximum sinking rate [m/d]
    w_sink%Pmicro_max = getpard("Micro max sink rate") / 86400.
    ! Particulate organic nitrogen (PON) detritus sinking rate [m/s]
    w_sink%D_PON = getpard("PON sink rate")
    ! Refractory nitrogen detritus sinking rate [m/s]
    w_sink%D_refr = getpard("refr sink rate")
    ! Biogenic silicon detritus sinking rate [m/s]
    w_sink%D_bSi = getpard("bSi sink rate")
  end subroutine read_sink_params
  

  subroutine build_biology_equations(grid, dt, Pmicro, Pnano, Z, NO, NH, & ! in
       Si, D_DON, D_PON, D_refr, D_bSi, Ft, wupwell)
    ! Build the terms for the diffusion/advection equations for the
    ! biology quantities.
    !
    ! This calculates the values of the precursor diffusion
    ! coefficients matrix (Bmatrix%bio%*), the RHS diffusion/advection
    ! term vectors (*_RHS%diff_adv%new), and the RHS sinking term
    ! vectors (*_RHS%sink).
    use precision_defs, only: dp
    use grid_mod, only: grid_
    use declarations, only: micro  ! *** This definitely needs to go away
    use turbulence, only: K
    use diffusion, only: diffusion_coeff, diffusion_bot_surf_flux, &
         diffusion_nonlocal_fluxes
    use find_upwell, only: vertical_advection
    use freshwater, only: freshwater_bio
    implicit none
    ! Arguments:
    type(grid_), intent(in) :: &
         grid  ! Grid arrays
    real(kind=dp), intent(in) :: &
         dt  ! Time step [s]
    real(kind=dp), dimension(0:), intent(in) :: &
         Pmicro, &  ! Micro phytoplankton
         Pnano,  &  ! Nano phytoplankton
         Z,      &  ! Micro Zooplankton
         NO,     &  ! Nitrate
         NH,     &  ! Ammonium
         Si,     &  ! Silicon
         D_DON,  &  ! Dissolved organic nitrogen detritus profile
         D_PON,  &  ! Particulate organic nitrogen detritus profile
         D_refr, &  ! Refractory nitrogen detritus profile
         D_bSi      ! Biogenic silicon detritus profile
    real(kind=dp), intent(in) :: &
         Ft  ! Total fresh water flux
    real(kind=dp), dimension(0:), intent(in) :: &
         wupwell  ! Profile of vertical upwelling velocity [m/s]
    ! Local variables:
    real(kind=dp) :: &
         surf_flux ! surface nutrient flux when all the river water on surface
    real(kind=dp), dimension(0:grid%M):: distrib_flux!distributed nutrient flux
    real(kind=dp), dimension(0:grid%M):: null ! null vector
    real(kind=dp), dimension(1:grid%M):: unit ! unit vector
    ! sinking defined by nutrient status
    real(kind=dp), dimension(1:grid%M):: Pmicro_w_sink 

    null = 0.d0
    unit = 1.d0

    ! Calculate the strength of the diffusion coefficients for biology
    ! model quantities.  They diffuse like salinity.
    call diffusion_coeff(grid, dt, K%S, & ! in
         Bmatrix%bio%new)                ! out

    ! Initialize the RHS *%diff_adv%new arrays, and calculate the diffusive
    ! fluxes at the bottom and top of the grid
    call freshwater_bio ('Pmicro', Pmicro(0:grid%M),     &
         surf_flux, distrib_flux)
    call diffusion_nonlocal_fluxes (grid, dt, K%S, null, &  ! in
         surf_flux, distrib_flux, Pmicro(grid%M+1),      &  ! in
         Pmicro_RHS%diff_adv%new)                           ! out
    call freshwater_bio ('Pnano', Pnano(0:grid%M),       &
         surf_flux, distrib_flux)
    call diffusion_nonlocal_fluxes(grid, dt, K%S, null,  &  ! in
         surf_flux, distrib_flux, Pnano(grid%M+1),       &  ! in
         Pnano_RHS%diff_adv%new)                            ! out
    call freshwater_bio ('Zoo', Z(0:grid%M),             &
         surf_flux, distrib_flux)
    call diffusion_nonlocal_fluxes (grid, dt, K%S, null, &  ! in
         surf_flux, distrib_flux, Z(grid%M+1),           &  ! in
         Z_RHS%diff_adv%new)                                ! out
    call freshwater_bio ('nitrate', NO(0:grid%M),        &
         surf_flux, distrib_flux)
    call diffusion_nonlocal_fluxes (grid, dt, K%S, null, &  ! in
         surf_flux, distrib_flux, NO(grid%M+1),          &  ! in
         NO_RHS%diff_adv%new)                               ! out
    call diffusion_bot_surf_flux(grid, dt, K%S, 0.d0,    &  ! in
         NH(grid%M+1),                                   &  ! in
         NH_RHS%diff_adv%new)                               ! out
    call freshwater_bio ('silicon', Si(0:grid%M),        &
         surf_flux, distrib_flux)
    call diffusion_nonlocal_fluxes (grid, dt, K%S, null, &  ! in
         surf_flux, distrib_flux, Si(grid%M+1),          &  ! in
         Si_RHS%diff_adv%new)                               ! out
    call diffusion_bot_surf_flux(grid, dt, K%S, 0.d0,    &  ! in
         D_DON(grid%M+1),                                &  ! in
         D_DON_RHS%diff_adv%new)                            ! out
    call diffusion_bot_surf_flux(grid, dt, K%S, 0.d0,    &  ! in
         D_PON(grid%M+1),                                &  ! in
         D_PON_RHS%diff_adv%new)                            ! out
    D_refr_RHS%diff_adv%new = 0.
    call diffusion_bot_surf_flux(grid, dt, K%S, 0.d0,    &  ! in
         D_bSi(grid%M+1),                                &  ! in
         D_bSi_RHS%diff_adv%new)                            ! out

    ! Add vertical advection due to upwelling
    call vertical_advection(grid, dt, Pmicro, wupwell, &
         Pmicro_RHS%diff_adv%new)
    call vertical_advection(grid, dt, Pnano, wupwell, &
         Pnano_RHS%diff_adv%new)
    call vertical_advection(grid, dt, Z, wupwell, &
         Z_RHS%diff_adv%new)
    call vertical_advection(grid, dt, NO, wupwell, &
         NO_RHS%diff_adv%new)
    call vertical_advection(grid, dt, NH, wupwell, &
         NH_RHS%diff_adv%new)
    call vertical_advection(grid, dt, Si, wupwell, &
         Si_RHS%diff_adv%new)
    call vertical_advection(grid, dt, D_DON, wupwell, &
         D_DON_RHS%diff_adv%new)
    call vertical_advection(grid, dt, D_PON, wupwell, &
         D_PON_RHS%diff_adv%new)
    call vertical_advection(grid, dt, D_bSi, wupwell, &
         D_bSi_RHS%diff_adv%new)

    ! Calculate the sinking term for the quantities that sink
    ! to  calculate sinking at the interfaces used the values above
     Pmicro_w_sink = w_sink%Pmicro_min * micro%Nlimit &
          + w_sink%Pmicro_max * (1 - micro%Nlimit)
     call sinking_advection(grid, dt, Pmicro, Pmicro_w_sink, &
          Pmicro_RHS%sink)
     call sinking_advection(grid, dt, D_PON, w_sink%D_PON*unit, &
          D_PON_RHS%sink)
     call sinking_advection(grid, dt, D_refr, w_sink%D_refr*unit, &
          D_refr_RHS%sink)
     call sinking_advection(grid, dt, D_bSi, w_sink%D_bSi*unit, &
          D_bSi_RHS%sink)
  end subroutine build_biology_equations


  subroutine sinking_advection(grid, dt, qty, w_sink, RHS)
    ! Calculate the sinking term of the semi-implicit PDEs for the biology
    ! quantities that sink (some classes of plankton, and some of detritus).
    ! use upwind advection

    use precision_defs, only: dp
    use grid_mod, only: grid_
    implicit none
    ! Arguments:
    type(grid_), intent(in) :: &
         grid  ! Grid arrays
    real(kind=dp) :: &
         dt  ! Time step [s]
    real(kind=dp), dimension(0:), intent(in) :: &
         qty  ! Profile array of sinking quantity
    real(kind=dp), dimension(1:), intent(in) :: &
         w_sink  ! Sinking velocity [m/s] on interfaces
    real(kind=dp), dimension(1:), intent(out) :: &
         RHS  ! RHS term vector for semi-implicit diffusion/advection PDE

    ! local variables
    integer :: ii ! interface index
    integer :: ig ! grid index 
    real(kind=dp), dimension(0:grid%M) :: flux ! positive is upward

    flux(0) = 0. ! no flux in or out from above

    do ii = 1,grid%M
       if (w_sink(ii).lt.0) then
          write (*,*) "Youve got biology sinking upward!"
          stop
       else
          ig = ii
          flux(ii) = - w_sink(ii) * qty(ig)
       endif
    enddo

    do ig = 1,grid%M
       ii = ig
       RHS(ig) = - dt * (flux(ii-1) - flux(ii))/grid%i_space(ii)
    enddo

  end subroutine sinking_advection


  subroutine new_to_old_bio_RHS()
    ! Copy %new component of the biology *_RHS%diff_adv arrays to the
    ! %old component for use by the IMEX semi-impllicit PDE solver.
    implicit none

    Pmicro_RHS%diff_adv%old = Pmicro_RHS%diff_adv%new
    Pnano_RHS%diff_adv%old = Pnano_RHS%diff_adv%new
    Z_RHS%diff_adv%old = Z_RHS%diff_adv%new
    NO_rhs%diff_adv%old = NO_RHS%diff_adv%new
    NH_rhs%diff_adv%old = NH_RHS%diff_adv%new
    Si_rhs%diff_adv%old = Si_RHS%diff_adv%new
    D_DON_rhs%diff_adv%old = D_DON_RHS%diff_adv%new
    D_PON_rhs%diff_adv%old = D_PON_RHS%diff_adv%new
    D_refr_rhs%diff_adv%old = D_refr_RHS%diff_adv%new
    D_bSi_rhs%diff_adv%old = D_bSi_RHS%diff_adv%new
  end subroutine new_to_old_bio_RHS


  subroutine new_to_old_bio_Bmatrix()
    ! Copy %new component of the Bmatrix%bio arrays to the
    ! %old component for use by the IMEX semi-impllicit PDE solver.
    implicit none

    Bmatrix%bio%old%sub = Bmatrix%bio%new%sub
    Bmatrix%bio%old%diag = Bmatrix%bio%new%diag
    Bmatrix%bio%old%sup = Bmatrix%bio%new%sup
  end subroutine new_to_old_bio_Bmatrix


  subroutine alloc_bio_RHS_variables(M)
    ! Allocate memory for arrays for right-hand sides of
    ! diffusion/advection equations for the biology model.
    use malloc, only: alloc_check
    implicit none
    ! Argument:
    integer, intent(in) :: M  ! Number of grid points
    ! Local variables:
    integer           :: allocstat  ! Allocation return status
    character(len=80) :: msg        ! Allocation failure message prefix

    msg = "Diffusion coefficients tridiagonal matrix arrays"
    allocate(Bmatrix%bio%new%sub(1:M), Bmatrix%bio%new%diag(1:M), &
         Bmatrix%bio%new%sup(1:M),                                &
         Bmatrix%bio%old%sub(1:M), Bmatrix%bio%old%diag(1:M),     &
         Bmatrix%bio%old%sup(1:M),                                &
         stat=allocstat)
    call alloc_check(allocstat, msg)
    msg = "Micro phytoplankton RHS arrays"
    allocate(Pmicro_RHS%diff_adv%new(1:M), Pmicro_RHS%diff_adv%old(1:M), &
         Pmicro_RHS%bio(1:M), Pmicro_RHS%sink(1:M), &
         stat=allocstat)
    call alloc_check(allocstat, msg)
    msg = "Nano phytoplankton RHS arrays"
    allocate(Pnano_RHS%diff_adv%new(1:M), Pnano_RHS%diff_adv%old(1:M), &
         Pnano_RHS%bio(1:M), &
         stat=allocstat)
    call alloc_check(allocstat, msg)
    msg = "Micro zooplankton RHS arrays"
    allocate(Z_RHS%diff_adv%new(1:M), Z_RHS%diff_adv%old(1:M), &
         Z_RHS%bio(1:M), &
         stat=allocstat)
    call alloc_check(allocstat, msg)
    msg = "Nitrate concentration RHS arrays"
    allocate(NO_RHS%diff_adv%new(1:M), NO_RHS%diff_adv%old(1:M), &
         NO_RHS%bio(1:M), &
         stat=allocstat)
    call alloc_check(allocstat, msg)
    msg = "Ammonium concentration RHS arrays"
    allocate(NH_RHS%diff_adv%new(1:M), NH_RHS%diff_adv%old(1:M), &
         NH_RHS%bio(1:M), &
         stat=allocstat)
    call alloc_check(allocstat, msg)
    msg = "Silicon concentration RHS arrays"
    allocate(Si_RHS%diff_adv%new(1:M), Si_RHS%diff_adv%old(1:M), &
         Si_RHS%bio(1:M), &
         stat=allocstat)
    call alloc_check(allocstat, msg)
    msg = "Dissolved organic nitrogen detritus RHS arrays"
    allocate(D_DON_RHS%diff_adv%new(1:M), D_DON_RHS%diff_adv%old(1:M), &
         D_DON_RHS%bio(1:M), &
         stat=allocstat)
    call alloc_check(allocstat, msg)
    msg = "Particulate organic nitrogen detritus RHS arrays"
    allocate(D_PON_RHS%diff_adv%new(1:M), D_PON_RHS%diff_adv%old(1:M), &
         D_PON_RHS%bio(1:M), D_PON_RHS%sink(1:M), &
         stat=allocstat)
    call alloc_check(allocstat, msg)
    msg = "Refractory nitrogen detritus RHS arrays"
    allocate(D_refr_RHS%diff_adv%new(1:M), D_refr_RHS%diff_adv%old(1:M), &
         D_refr_RHS%bio(1:M), D_refr_RHS%sink(1:M), &
         stat=allocstat)
    call alloc_check(allocstat, msg)
    msg = "Biogenic silicon detritus RHS arrays"
    allocate(D_bSi_RHS%diff_adv%new(1:M), D_bSi_RHS%diff_adv%old(1:M), &
         D_bSi_RHS%bio(1:M), D_bSi_RHS%sink(1:M), &
         stat=allocstat)
    call alloc_check(allocstat, msg)
  end subroutine alloc_bio_RHS_variables


  subroutine dalloc_bio_RHS_variables
    ! Deallocate memory from arrays for right-hand sides of
    ! diffusion/advection equations for the biology model.
    use malloc, only: dalloc_check
    implicit none
    ! Local variables:
    integer           :: dallocstat  ! Allocation return status
    character(len=80) :: msg         ! Allocation failure message prefix

    msg = "Diffusion coefficients tridiagonal matrix arrays"
    deallocate(Bmatrix%bio%new%sub, Bmatrix%bio%new%diag, &
         Bmatrix%bio%new%sup,                             &
         Bmatrix%bio%old%sub, Bmatrix%bio%old%diag,       &
         Bmatrix%bio%old%sup,                             &
         stat=dallocstat)
    call dalloc_check(dallocstat, msg)
    msg = "Micro phytoplankton RHS arrays"
    deallocate(Pmicro_RHS%diff_adv%new, Pmicro_RHS%diff_adv%old, &
         Pmicro_RHS%bio, Pmicro_RHS%sink, &
         stat=dallocstat)
    call dalloc_check(dallocstat, msg)
    msg = "Nano phytoplankton RHS arrays"
    deallocate(Pnano_RHS%diff_adv%new, Pnano_RHS%diff_adv%old, &
         Pnano_RHS%bio, &
         stat=dallocstat)
    call dalloc_check(dallocstat, msg)
    msg = "Micro zooplankton RHS arrays"
    deallocate(Z_RHS%diff_adv%new, Z_RHS%diff_adv%old, &
         Z_RHS%bio, &
         stat=dallocstat)
    call dalloc_check(dallocstat, msg)
    msg = "Nitrate concentration RHS arrays"
    deallocate(NO_RHS%diff_adv%new, NO_RHS%diff_adv%old, &
         NO_RHS%bio, &
         stat=dallocstat)
    call dalloc_check(dallocstat, msg)
    msg = "Ammonium concentration RHS arrays"
    deallocate(NH_RHS%diff_adv%new, NH_RHS%diff_adv%old, &
         NH_RHS%bio, &
         stat=dallocstat)
    call dalloc_check(dallocstat, msg)
    msg = "Silicon concentration RHS arrays"
    deallocate(Si_RHS%diff_adv%new, Si_RHS%diff_adv%old, &
         Si_RHS%bio, &
         stat=dallocstat)
    call dalloc_check(dallocstat, msg)
    msg = "Dissolved organic nitrogen detritus RHS arrays"
    deallocate(D_DON_RHS%diff_adv%new, D_DON_RHS%diff_adv%old, &
         D_DON_RHS%bio, &
         stat=dallocstat)
    call dalloc_check(dallocstat, msg)
    msg = "Particulate organic nitrogen detritus RHS arrays"
    deallocate(D_PON_RHS%diff_adv%new, D_PON_RHS%diff_adv%old, &
         D_PON_RHS%bio, D_PON_RHS%sink, &
         stat=dallocstat)
    call dalloc_check(dallocstat, msg)
    msg = "Refractory nitrogen detritus RHS arrays"
    deallocate(D_refr_RHS%diff_adv%new, D_refr_RHS%diff_adv%old, &
         D_refr_RHS%bio, D_refr_RHS%sink, &
         stat=dallocstat)
    call dalloc_check(dallocstat, msg)
    msg = "Biogenic silicon detritus RHS arrays"
    deallocate(D_bSi_RHS%diff_adv%new, D_bSi_RHS%diff_adv%old, &
         D_bSi_RHS%bio, D_bSi_RHS%sink, &
         stat=dallocstat)
    call dalloc_check(dallocstat, msg)
  end subroutine dalloc_bio_RHS_variables

end module biology_eqn_builder
