! $Id$
! $Source$

program SOG      
  ! Coupled physical and biological model of the Strait of Georgia

  ! Things that we will use from other modules
  !
  ! Type definitions:
  use datetime, only: datetime_
  !
  ! Parameter values:
  use physics_model, only: f, g
  !
  ! Variables:
  use grid_mod, only: grid
  use core_variables, only: U, V, T, S, P, Z, N, Si, D
  use water_properties, only: rho, alpha, beta, Cp
  use physics_model, only: B, dPdx_b, dPdy_b
  ! *** Temporary until IMEX refactoring is completed
  use IMEX_solver, only: Hvector
  ! *** Temporary until physics equations refactoring is completed
  use physics_eqn_builder, only: nBmatrix
  !
  ! Subroutines and functions:
  use core_variables, only: alloc_core_variables, dalloc_core_variables
  use grid_mod, only: init_grid, dalloc_grid, gradient_g, gradient_i, &
       interp_i
  use physics_model, only: init_physics, double_diffusion, &
       baroclinic_P_gradient, new_to_old_physics, dalloc_physics_variables
  use water_properties, only: calc_rho_alpha_beta_Cp_profiles
  use physics_eqn_builder, only: build_physics_equations
  use air_sea_fluxes, only: wind_stress
  use biology_model, only: init_biology, calc_bio_rate
  use biology_eqn_builder, only: build_biology_equations, new_to_old_bio_RHS, &
       new_to_old_bio_Bmatrix
  use NPZD, only: dalloc_biology_variables
  use IMEX_solver, only: init_IMEX_solver, solve_bio_eqns, &
       dalloc_IMEX_variables
  use input_processor, only: init_input_processor, getpars, getpari, &
       getpard, getparl
  use timeseries_output, only: init_timeseries_output, write_std_timeseries, &
       timeseries_output_close
  use profiles_output, only: init_profiles_output, write_std_profiles, &
       profiles_output_close
  use user_output, only: write_user_timeseries, write_user_profiles
  use mixing_layer, only: find_mixing_layer_depth
  use find_upwell, only: upwell_profile, vertical_advection
  use diffusion, only: diffusion_coeff, diffusion_nonlocal_fluxes, &
       diffusion_bot_surf_flux
  use fitbottom, only: bot_bound_time, bot_bound_uniform
  use forcing, only: read_variation, read_forcing, get_forcing
  use precision_defs, only: dp
  use io_unit_defs, only: stripped_infile, stdout
  use unit_conversions, only: KtoC
  use datetime, only: os_datetime, calendar_date, clock_time, datetime_str

  ! Inherited modules
  ! *** Goal is to make these go away
  use declarations
  use surface_forcing, only: del_o, precision, step_guess, step_min
  use initial_sog, only: initial_mean
  ! Subroutine & function modules:
  ! (Wrapping subroutines and functions in modules provides compile-time
  !  checking of number and type of arguments - but not order!)
  ! *** These should eventually end up in refactored modules
  use Coriolis_and_pg_mod
  use define_flux_mod
  use tridiag_mod

  implicit none

  real(kind=dp) :: &
       upwell
  ! Upwelling constant (tuned parameter)
  !*** read in here, used by surface_flux_sog : eventually should be local
  ! to the surface_forcing module (not be be confused with current 
  ! surface_forcing module
  real(kind=dp) :: upwell_const, S_riv, sumS=0, sumSriv=0
  integer :: scount=0

  ! Internal wave breaking eddy viscosity for momentum and scalars
  ! (tuned parameters)
  real(kind=dp) :: nu_w_m, nu_w_s

  ! Interpolated river flows
  real(kind=sp) :: Qinter  ! Fraser River
  real(kind=sp) :: Einter  ! Englishman River
  ! Current time met data
  real(kind=sp) :: cf_value, atemp_value, humid_value
  ! Wind data
  real(kind=dp) unow, vnow

  integer :: wind_n, met_n, river_n ! length of various forcing files

  ! Initial month parameter read from run control input file
  ! and passed to new_year()
  integer :: month_o

  ! Iteration limit for inner loop that solves KPP turbulence model
  integer :: niter

  ! Variables for velocity component convergence metric
  double precision :: &
       &uprev, vprev, delu, delv

  ! Code identification string (maintained by CVS), for output file headers
  character*70 :: &
       codeId = "$Source$"

  ! Date/time structures for output file headers
  type(datetime_) :: runDatetime     ! Date/time of code run
  type(datetime_) :: startDatetime   ! Date/time of initial conditions


  ! ---------- Beginning of Initialization Section ----------
  ! Get the current date/time from operating system to timestamp the
  ! run with
  call os_datetime(runDatetime)

  ! Initialize the input processor, including preparing the infile
  ! (read from stdin) to be read by the init_* subroutines called
  ! below.
  call init_input_processor(datetime_str(runDatetime))

  ! Read the run start date/time, duration, and time step
  ! *** Not sure where these should go?
  year_o = getpari("year_o")
  month_o = getpari("month_o")
  day_o = getpari("yr_day_o")
  t_o = getpari("t_o")
  t_f = getpard("run_dur")
  dt = getpard("dt")
  ! Calculate the number of time steps for the run (note that int()
  ! rounds down)
  steps = 1 + int((t_f - t_o) / dt)
  ! Calculate the month number and month day for output file headers
  startDatetime%yr = year_o
  startDatetime%yr_day = day_o
  startDatetime%day_sec = t_o
  call calendar_date(startDatetime)
  call clock_time(startDatetime)

  ! Initialize time series writing code
  call init_timeseries_output(codeId, datetime_str(runDatetime), &
       startDatetime)
  ! Initialize profiles writing code
  call init_profiles_output(codeId, datetime_str(runDatetime), &
       startDatetime)

  ! Initialize the grid
  call init_grid

  ! Allocate memory for the profiles arrays of the core variables that
  ! we are calculating, i.e. U, V, T, S, etc.
  call alloc_core_variables(grid%M)

  ! Initialize the physics model
  call init_physics(grid%M)

  ! Initialize the biology model
  call init_biology(grid%M)

  ! Initialize the IMEX semi-implicit PDE solver
  ! *** This should eventually go into init_numerics()
  call init_IMEX_solver(grid%M)

  ! Allocate memory
  ! *** It would be nice if everything in allocate[13] could end up in 
  ! *** alloc_* subroutines private to various modules, and called by 
  ! *** their init_* subroutines.
  CALL allocate1(grid%M, alloc_stat) 
  DO xx = 1,12
     IF (alloc_stat(xx) /= 0) THEN
        PRINT "(A)","ALLOCATION failed.  KPP.f  xx:"
        PRINT *,xx,alloc_stat(xx)
        STOP
     END IF
  END DO

  ! Length of forcing data files (move to be read in)
  wind_n = 46056 - 8 ! with wind shifted to LST we lose 8 records
  met_n = 1918
  river_n = 1826
  call read_forcing (wind_n, met_n, river_n)
  call read_variation

  ! Read the physic model parameter values
  upwell_const = getpard("upwell_const")
  nu_w_m = getpard('nu_w_m')         ! Internal wave mixing momentum
  nu_w_s = getpard('nu_w_s')         ! Internal wave mixing scalar
  Fw_scale = getpard('Fw_scale')     ! Fresh water scale factor for river flows
  Fw_surface = getparl('Fw_surface') ! Add all fresh water on surface?
  if (.not. Fw_surface) then
     Fw_depth = getpard('Fw_depth')  ! Depth to distribute freshwater flux over
  endif

  CALL initialize ! initializes everything (biology too)

  ! Read the cruise id from stdin to use to build the file name for
  ! nutrient initial conditions data file
  cruise_id = getpars("cruise_id")
  CALL initial_mean(U%new, V%new, T%new, S%new, P%micro, P%nano, &
       Z, N%O, N%H, Si, D%DON, D%PON, D%refr, D%bSi, &
       h%new, dPdx_b, dPdy_b, grid, cruise_id)

  ! Read the iteration count limit for the physics model implicit
  ! solver
  ! *** This should go somewhere else - maybe init_physics() (yet to come)
  niter = getpari("niter")

  ! Initialize mixing, and mixed layer depths
  ! *** It would be cool to calculate these from the initial salinity
  ! *** profile read from the CTD file.
  IF(h%new < grid%d_g(1))THEN 
     h%new = grid%d_g(1)
  END IF

  ! Initialize the profiles of the water column properties
  ! Density (rho), thermal expansion (alpha) and saline contraction (beta)
  ! coefficients, and specific heat capacity (Cp).
  ! This calculates the values of the grid point depth arrays (*%g) of
  ! the 4 water properties.
  call calc_rho_alpha_beta_Cp_profiles(T%new, S%new)
  ! Interpolate the values of density and specific heat capacity at
  ! the grid layer interface depths
  rho%i = interp_i(rho%g)
  alpha%i = interp_i(alpha%g)
  beta%i = interp_i(beta%g)
  Cp%i = interp_i(Cp%g)
  ! Calculate the gradients of density, and thermal expansion and saline
  ! contraction coefficients at the grid layer interface depths
  rho%grad_i = gradient_i(rho%g)
  alpha%grad_i = gradient_i(alpha%g)
  beta%grad_i = gradient_i(beta%g)

  ! Close the input parameters file
  close(stripped_infile)
  ! ---------- End of Initialization Section ----------

  do time_step = 1, steps  !---------- Beginning of the time loop ----------
     ! Store %new components of various variables from time step just
     ! completed in %old their components for use in the next time
     ! step
     CALL define_sog(time_step)
     call new_to_old_physics()
     call new_to_old_bio_RHS()
     call new_to_old_bio_Bmatrix()

     ! Get forcing data
     call get_forcing (year, day, day_time, &
          Qinter, Einter, cf_value, atemp_value, humid_value, &
          unow, vnow)

     CALL irradiance_sog(cf_value, day_time, day, &
          I, I_par, grid, jmax_i, Q_sol, euph, Qinter, P%micro, P%nano)

     DO count = 1, niter !------ Beginning of the implicit solver loop ------
        ! Calculate gradient pofiles of the velocity field and water column
        ! temperature, and salinity at the grid layer interface depths
        U%grad_i = gradient_i(U%new)
        V%grad_i = gradient_i(V%new)
        T%grad_i = gradient_i(T%new)
        S%grad_i = gradient_i(S%new)

        ! *** I think this is finding the depths of the grid point
        ! *** and grid interface that bound the mixed layer depth
        ! d_g(j_max_g) > h ie, hh%g is grid point just below h
        CALL find_jmax_g(h, grid)
        ! d_i(j_max_i) > h ie h%i is just above h
        CALL find_jmax_i(h, grid)

        ! Calculate surface forcing components
        ! *** Confirm that all of these arguments are necessary
        CALL surface_flux_sog(grid%M, rho%g, w, wt_r, S%new(1),        &
             S%old(1), S_riv, T%new(0), j_gamma, I, Q_t,        &
             alpha%i(0), Cp%i(0), beta%i(0), unow, vnow, cf_value/10.,    &
             atemp_value, humid_value, Qinter, stress, &
             day, dt/h%new, h, upwell_const, upwell, Einter,       &
             u%new(1), dt, Fw_surface, Fw_scale, Ft, count) 
        call wind_stress (unow, vnow, rho%g(1), &
             w%u(0), w%v(0))

        ! Calculate nonturbulent heat flux profile
        Q_n = I / (Cp%i * rho%i)

        ! Calculate the nonturbulent fresh water flux profile, and its
        ! contribution to the salinity profile
        ! *** Move to a subroutine
        if (Fw_surface) then
           F_n = 0.
        else
           Fw = Ft * exp(-grid%d_i / (Fw_depth * h%old))
           F_n = S%new * Fw
        endif

        ! Store the surface buoyancy forcing value from the previous
        ! iteration so we can use it to blend with the new value to
        ! help the implicit solver converge more quickly
        Bf_old = Bf

        ! Calculate buoyancy profile, and surface buoyancy forcing
        CALL buoyancy(grid, T%new, S%new, h, I, F_n, w%b(0), &  ! in
             rho%g, alpha%g, beta%g, Cp%g, Fw_surface,       &  ! in
             B, Bf)                                             ! out

        ! Blend the values of the surface buoyancy forcing from current
        ! and previous iteration to help convergence
        Bf = (count * Bf_old + (niter - count) * Bf) / niter !!$

!!$        ! Calculate baroclinic pressure gradient components
!!$        ! *** This might be a better place to calculate these gradients
!!$        ! *** than at the end of the implicit loop.
!!$        call baroclinic_P_gradient(grid, dt, U%new, V%new, rho%g, &
!!$             stress%u%new, stress%v%new, dPdx_b, dPdy_b)

        CALL fun_constants(u_star, w_star, L_star, w, Bf, h%new)   !test conv

        CALL stability   !stable = 0 (unstable), stable = 1 (stable), stable = 2 (no forcing)  this is the stability of the water column.

        IF (u_star /= 0.)  THEN         !Wind stress /= 0.
           CALL ND_flux_profile(grid, L_star, phi)   ! define flux profiles aka (B1)
           CALL vel_scales(grid, omega, phi, u_star,L_star,h)
           ! calculates wx (13) as omega (not Omega's)
        ELSE IF (u_star == 0. .AND. Bf < 0.) THEN        !Convective unstable limit
           CALL convection_scales(grid,omega,h, w_star)   !test conv
           ! calculates wx (15) as omega (not Omega's)
        ELSE                !  No surface forcing or Bf > 0. 
           omega%m%value = 0.
           omega%s%value = 0.
           omega%m%div = 0.
           omega%s%div = 0.
           omega%m%h = 0.
           omega%m%h = 0.
        END IF

        CALL shear_diff(grid, U%grad_i, V%grad_i, rho, K%u%shear)  !test conv  !density instead of linear B  ! calculates ocean interior shear diffusion

        ! Calculate double diffusion mixing   
        call double_diffusion(grid%M, T%grad_i, S%grad_i, &  ! in
             alpha%i, beta%i, &                              ! in
             K%t%dd, K%s%dd)                                 ! out

        ! Define interior diffusivity K%x%total, nu_w_m and nu_w_s
        ! constant internal wave mixing
        K%u%total = 0.0
        K%s%total = 0.0
        K%t%total = 0.0          
        do xx = 1, grid%M
           K%u%total(xx) = K%u%shear(xx) + nu_w_m + K%s%dd(xx) !test conv
           K%s%total(xx) = K%u%shear(xx) + nu_w_s + K%s%dd(xx) !test conv
           K%t%total(xx) = K%u%shear(xx) + nu_w_s + K%t%dd(xx) !test conv
        enddo

        CALL interior_match(grid, h, K%t, nu_w_s)  ! calculate nu (D5)
        CALL interior_match(grid, h, K%u, nu_w_m)
        CALL interior_match(grid, h, K%s, nu_w_s)   
        !test conv

        IF (u_star /= 0. .OR. Bf < 0.) THEN
           CALL interior_match2(omega, L_star, u_star, h, grid) !test conv
        END IF

        !Define shape functions G_shape%x

        IF (u_star /= 0. .OR. Bf < 0.) THEN
           CALL shape_parameters(K%u,omega%m, h%new, a2%m, a3%m)
           CALL shape_parameters(K%s,omega%s, h%new, a2%s, a3%s)
           CALL shape_parameters(K%t,omega%s, h%new, a2%t, a3%t) !test conv
        ELSE
           a2%m = -2.0
           a2%s = -2.0
           a2%t = -2.0
           a3%m = 1.0
           a3%s = 1.0
           a3%t = 1.0
        END IF

        G_shape%m = 0.   
        G_shape%s = 0.
        G_shape%t = 0.

        DO index = 0, h%i  !Just in Surface Layer  (h%i-1)?
           IF (grid%d_i(index) > h%new) THEN
              G_shape%m(index) = 0.   
              G_shape%s(index) = 0.
              G_shape%t(index) = 0.
           ELSE ! not sigma = d_i/h%new, a1 = 1 
              G_shape%m(index) = grid%d_i(index)/h%new*(1.0 + grid%d_i(index)/h%new*(a2%m + &
                   grid%d_i(index)/h%new*a3%m)) ! (11)
              G_shape%s(index) = grid%d_i(index)/h%new*(1.0 + grid%d_i(index)/h%new*(a2%s + &
                   grid%d_i(index)/h%new*a3%s)) 
              G_shape%t(index) = grid%d_i(index)/h%new*(1.0 + grid%d_i(index)/h%new*(a2%t + &
                   grid%d_i(index)/h%new*a3%t))
           END IF
        END DO   !test conv

        !Define K%x%ML profiles (10)

        K%u%ML = 0.0
        K%s%ML = 0.0
        K%t%ML = 0.0

        DO xx = 0, h%i  !use only up to h%i-1
           IF (grid%d_i(xx) > h%new) THEN
              K%u%ML(xx) = 0.0
              K%s%ML(xx) = 0.0
              K%t%ML(xx) = 0.0
           ELSE
              K%u%ML(xx) = h%new*omega%m%value(xx)*G_shape%m(xx)
              K%s%ML(xx) = h%new*omega%s%value(xx)*G_shape%s(xx) 
              K%t%ML(xx) = h%new*omega%s%value(xx)*G_shape%t(xx)   !test conv
           END IF
           IF (K%u%ML(xx) < 0. .OR. K%s%ML(xx) < 0. .OR. K%t%ML(xx) < 0.) THEN
              PRINT "(A)","Diffusivities < 0. See KPP.f90:  time, Jday, h%new"
              PRINT *,time,day,h%new
              PRINT "(A)","K%u%ML(xx),K%s%ML(xx),K%t%ML(xx)"
              PRINT *,K%u%ML(xx),K%s%ML(xx),K%t%ML(xx)
              STOP
           END IF
        END DO

        ! K star, enhance diffusion at grid point just below or above h (see (D6)
        CALL modify_K(grid, h, K%u) !test conv  !!must modify G_shape!!
        CALL modify_K(grid, h, K%s) !test conv  !!If G_shape is used again, use modified value!!
        CALL modify_K(grid, h, K%t) !test conv

        !!Define flux profiles, w,  and Q_t (vertical heat flux)!!
        !!Don-t need for running of the model!!  (not sure what this means)
        IF (u_star /= 0. .OR. Bf /= 0.) THEN
           CALL def_gamma(L_star, grid, w,  wt_r, h, gamma, Bf, omega) 
           ! calculates non-local transport (20)
        ELSE
           gamma%t = 0.
           gamma%s = 0.
           gamma%m = 0.
        END IF

        !defines w%x, K%x%all, K%x%old, Bf%b, and F_n 
        CALL define_flux(grid, U%grad_i, V%grad_i, T%grad_i, S%grad_i, &
             alpha, beta)
! *** End of turbulence code???

        ! Build the terms of the semi-implicit diffusion/advection
        ! PDEs with Coriolis and baroclinic pressure gradient terms for
        ! the physics quantities.
        !
        ! This calculates the values of the precursor diffusion
        ! coefficients matrices (Bmatrix%vel%*, Bmatrix%T%*, and
        ! Bmatrix%S%*), the RHS diffusion/advection term vectors
        ! (*_RHS%diff_adv%new), and the RHS Coriolis and barolcinic
        ! pressure gradient term vectors (*_RHS%C_pg).
        call build_physics_equations(grid, dt, U%new, V%new, T%new, S%new) ! in
     
!!$        ! Calculate matrix B (changed to Amatrix later) 
!!$        call diffusion_coeff(grid, dt, K%t%all, &
!!$             Bmatrix%t)
!!$        call diffusion_coeff(grid, dt, K%s%all, &
!!$             Bmatrix%s)
!!$        call diffusion_coeff(grid, dt, K%u%all, &
!!$             Bmatrix%u)
! *** Physics equations refactoring bridge code
Bmatrix%u%A = nBmatrix%vel%new%sub
Bmatrix%u%B = nBmatrix%vel%new%diag
Bmatrix%u%C = nBmatrix%vel%new%sup
Bmatrix%T%A = nBmatrix%T%new%sub
Bmatrix%T%B = nBmatrix%T%new%diag
Bmatrix%T%C = nBmatrix%T%new%sup
Bmatrix%S%A = nBmatrix%S%new%sub
Bmatrix%S%B = nBmatrix%S%new%diag
Bmatrix%S%C = nBmatrix%S%new%sup

        ! Start calculation of RHS called Gvector (D9) (D10)
        call diffusion_nonlocal_fluxes(grid, dt, K%t%all, gamma%t, &   ! in
             w%t(0), -Q_n, T%new(grid%M+1), &                          ! in
             Gvector%t)                                                ! out
        call diffusion_nonlocal_fluxes(grid, dt, K%s%all, gamma%s, &   ! in
             w%s(0), F_n, S%new(grid%M+1), &                           ! in
             Gvector%s)                                                ! out
        call diffusion_bot_surf_flux (grid, dt, K%u%all, w%u(0), &    ! in
             U%new(grid%M+1), &                                       ! in
             Gvector%u)                                               ! out
        call diffusion_bot_surf_flux (grid, dt, K%u%all, w%v(0), &    ! in
             V%new(grid%M+1), &                                       ! in
             Gvector%v)                                               ! out

        ! Calculate profile of upwelling velocity
        call upwell_profile(grid, Qinter, upwell, wupwell)
        ! Upwell salinity, temperature, and u & v velocity components
        ! similarly to nitrates
        call vertical_advection (grid, dt, S%new, wupwell, &
             Gvector%s)
        call vertical_advection (grid, dt, T%new, wupwell, &
             Gvector%t)
        call vertical_advection (grid, dt, U%new, wupwell, &
             Gvector%u)
        call vertical_advection (grid, dt, V%new, wupwell, &
             Gvector%v)

        ! Calculate the Coriolis and baroclinic pressure gradient
        ! components of the G vector for each velocity component
        call Coriolis_and_pg(dt, V%new, dPdx_b, &
             Gvector_c%u)
        call Coriolis_and_pg(dt, -U%new, dPdy_b, &
             Gvector_c%v)      

        IF (time_step == 1 .AND. count  == 1) THEN
           Bmatrix_o%t%A = Bmatrix%t%A
           Bmatrix_o%t%B = Bmatrix%t%B
           Bmatrix_o%t%C = Bmatrix%t%C
           Bmatrix_o%s%A = Bmatrix%s%A
           Bmatrix_o%s%B = Bmatrix%s%B
           Bmatrix_o%s%C = Bmatrix%s%C         
           Bmatrix_o%u%A = Bmatrix%u%A
           Bmatrix_o%u%B = Bmatrix%u%B
           Bmatrix_o%u%C = Bmatrix%u%C

           Gvector_o%t = Gvector%t
           Gvector_o%s = Gvector%s
           Gvector_o%u = Gvector%u
           Gvector_o%v = Gvector%v

           Gvector_co%u = Gvector_c%u
           Gvector_co%v = Gvector_c%v
        END IF
!!$        ! Store %new components of RHS and Bmatrix variables in %old
!!$        ! their components for use by the IMEX solver.  Necessary for the
!!$        ! 1st time step because the values just calculated are a better
!!$        ! estimate than zero.
!!$        if (time_step == 1) then
!!$           call new_to_old_phys_RHS()
!!$           call new_to_old_phys_Bmatrix()
!!$        endif

!!!!! IMEX SCHEME !!!!   

!!$        ! Solve the semi-implicit diffusion/advection PDEs with
!!$        ! Coriolis and baroclinic pressure gradient terms for the
!!$        ! physics quantities.
!!$        call solve_phys_eqns(grid%M, U%new, V%new, T%new, S%new, &  ! in
!!$             day, time)                                             ! in

        CALL matrix_A (Amatrix%u, Bmatrix%u) ! changed from Bmatrix to Amatrix
        CALL matrix_A (Amatrix%t, Bmatrix%t) ! now have (D9)
        CALL matrix_A (Amatrix%s, Bmatrix%s)

        ! Add in Coriolis term (Gvector_c) and previous value to H vector (D7)
        call phys_Hvector(grid%M, U%old, Gvector%u, Gvector_o%u, &  ! in
             Gvector_c%u, Gvector_co%u, Bmatrix_o%u,             &  ! in
             Hvector%U)                                             ! out
        call phys_Hvector(grid%M, V%old, Gvector%v, Gvector_o%v, &  ! in
             Gvector_c%v, Gvector_co%v, Bmatrix_o%u,             &  ! in
             Hvector%V)                                             ! out
        ! Add Xt to H vector (D7)
        call phys_Hvector(grid%M, T%old, Gvector%t, Gvector_o%t, &  ! in
             null_vector, null_vector, Bmatrix_o%t,              &  ! in
             ! null_vector because no Coriolis or pressure gradients terms
             Hvector%T)
        call phys_Hvector(grid%M, S%old, Gvector%s, Gvector_o%s, &  ! in
             null_vector, null_vector, Bmatrix_o%s,              &  ! in
             ! null_vector because no Coriolis or pressure gradients terms
             Hvector%S)                                             ! out

        ! Solves tridiagonal system for physics quantities
        call tridiag(Amatrix%u%A, Amatrix%u%B, Amatrix%u%C, Hvector%U, &
             U%new(1:grid%M))
        call tridiag(Amatrix%u%A, Amatrix%u%B, Amatrix%u%C, Hvector%V, &
             V%new(1:grid%M))
        call tridiag(Amatrix%t%A, Amatrix%t%B, Amatrix%t%C, Hvector%T, &
             T%new(1:grid%M))
        call tridiag(Amatrix%s%A, Amatrix%s%B, Amatrix%s%C, Hvector%S, &
             S%new(1:grid%M))

        ! Update boundary conditions at surface, and bottom of grid
        U%new(0) = U%new(1)
        V%new(0) = V%new(1)
        S%new(0) = S%new(1)   
        T%new(0) = T%new(1)
        U%new(grid%M+1) = U%new(grid%M) 
        V%new(grid%M+1) = V%new(grid%M) 
        T%new(grid%M+1) = T%old(grid%M+1)
        S%new(grid%M+1) = S%old(grid%M+1)

        ! Update gradient pofiles of the water column
        ! temperature, and salinity at the grid layer interface depths
        T%grad_i = gradient_i(T%new)
        S%grad_i = gradient_i(S%new)

        ! Update the profiles of the water column properties
        ! Density (rho), thermal expansion (alpha) and saline
        ! contraction (beta) coefficients, and specific heat capacity (Cp)
        ! This calculates the values of the grid point depth arrays (*%g) of
        ! the 4 water properties.
        call calc_rho_alpha_beta_Cp_profiles(T%new, S%new)
        ! Interpolate the values of density and specific heat capacity at
        ! the grid layer interface depths
        rho%i = interp_i(rho%g)
        alpha%i = interp_i(alpha%g)
        beta%i = interp_i(beta%g)
        Cp%i = interp_i(Cp%g)
        ! Calculate the gradients of thermal expansion and saline
        ! contraction coefficients at the grid layer interface depths
        rho%grad_i = gradient_i(rho%g)
        alpha%grad_i = gradient_i(alpha%g)
        beta%grad_i = gradient_i(beta%g)

        ! Update buoyancy profile, surface buoyancy forcing, and local
        ! buoyancy frequency
        ! *** Local buoyancy freqency may be part of mixing layer
        ! *** depth calculation
        CALL buoyancy(grid, T%new, S%new, h, I, F_n, w%b(0), &  ! in
             rho%g, alpha%g, beta%g, Cp%g, Fw_surface,       &  ! in
             B, Bf)                                             ! out
        rho%grad_g = gradient_g(rho%g)
        DO xx = 1,grid%M
           N_2_g(xx) = -(g / rho%g(xx)) * rho%grad_g(xx)
        END DO

!!!Define turbulent velocity shear V_t_square for Ri_b  use last count values

!!!Calculate beta_t: need h%e and h%e%i and new w%b  w%i vertical turbulent flux of i


        w%b_err(0) = 0.

        DO xx = 1, grid%M           !uses K%old and T%new
           w%t(xx) = -K%t%all(xx)*(T%grad_i(xx) - gamma%t(xx))
           w%s(xx) = -K%s%all(xx)*(S%grad_i(xx) - gamma%s(xx))
           w%b(xx) = g * (alpha%i(xx) * w%t(xx) - beta%i(xx) * w%s(xx))   
           w%b_err(xx) = g * (alpha%grad_i(xx) * w%t(xx) &
                - beta%grad_i(xx) * w%s(xx))
        END DO


        ! beta_t is the ratio of the entrainment flux to the surface buoyancy flux -- set equal to -0.2 when convection is occuring

        if (w%b(0)>0) then
           beta_t = -0.2
        else
           beta_t = 0.
        endif

        V_t_square = 0.


        IF (beta_t < 0.) THEN  !omega, vertical velocity scale
           CALL def_v_t_sog(grid, h, N_2_g, omega%s%value, V_t_square, beta_t, L_star) !test conv  
           ! V_t_square, the turbulent velocity shear, (23)
        END IF

        ! Calculate the profile of bulk Richardson number (Large, etal
        ! (1994) eqn 21)
        CALL define_Ri_b_sog(grid, h, surface_h, U%new, V%new, rho%g, Ri_b, &
             V_t_square, N_2_g)
        ! Find the mixing layer depth by comparing Ri_b to Ri_c
        call find_mixing_layer_depth (grid, Ri_b, year, day, day_time, count, &
             h_i)
        ! Apply the Ekman depth criterion to the mixing layer depth
        ! when stable forcing exists
        ! *** This code could probably go into ML_height_sog.
        h_i = (h_i + h%new)/2.0 !use the average value!!!!!!!!##
        IF (stable == 1) THEN          !Stability criteria 
           h_Ekman = 0.7*u_star/f         ! see (24) and surrounding
           IF (h_Ekman < L_star) THEN
              IF (h_i > h_Ekman) THEN
                 h_i = h_Ekman
                 write (*,*) 'Ekman', h_i
              END IF
           ELSE
              IF (h_i > L_star) THEN
                 h_i = L_star
              END IF
           END IF
        END IF

        IF (h_i < grid%d_g(1)) THEN  !minimum mixing
           h_i = grid%d_g(1)
        END IF

        CALL find_jmax_g(h,grid) !***! can't be less than the grid

        del = ABS(h_i - h%new)/grid%i_space(1)
        h%new = h_i
        CALL find_jmax_g(h,grid)

        ! Calculate baroclinic pressure gradient components
        call baroclinic_P_gradient(grid, dt, U%new, V%new, rho%g, &
             w%u(0), w%v(0), dPdx_b, dPdy_b)

        ! Calculate convergence metric for velocity field
        delu = abs(U%new(1)-uprev)/0.01
        delv = abs(V%new(1)-vprev)/0.01
        uprev = U%new(1)
        vprev = V%new(1)
        if (del.lt.delu) del = delu
        if (del.lt.delv) del = delv

        ! Test for convergence of implicit solver
        if (count >= 2 .and. del < del_o) then
           exit
        endif
     enddo  !---------- End of the implicit solver loop ----------

     !---------- Biology Model ----------
     !
     ! Solve the biology model ODEs to advance the biology quantity values
     ! to the next time step, and calculate the growth - mortality terms
     ! (*_RHS%bio) of the semi-implicit diffusion/advection equations.
     call calc_bio_rate(time, day, dt, grid%M, precision, step_guess, step_min,  &
          T%new(0:grid%M), I_Par, P%micro, P%nano, Z, N%O, N%H, Si,              &
          D%DON, D%PON, D%refr, D%bSi)

     ! Build the rest of the terms of the semi-implicit diffusion/advection
     ! PDEs for the biology quantities.
     !
     ! This calculates the values of the precursor diffusion
     ! coefficients matrix (Bmatrix%bio%*), the RHS diffusion/advection
     ! term vectors (*_RHS%diff_adv%new), and the RHS sinking term
     ! vectors (*_RHS%sink).
     call build_biology_equations(grid, dt, P%micro, P%nano, Z, N%O, N%H, &! in
          Si, D%DON, D%PON, D%refr, D%bSi, Ft, wupwell)                    ! in

     ! Store %new components of RHS and Bmatrix variables in %old
     ! their components for use by the IMEX solver.  Necessary for the
     ! 1st time step because the values just calculated are a better
     ! estimate than zero.
     if (time_step == 1) then
        call new_to_old_bio_RHS()
        call new_to_old_bio_Bmatrix()
     endif

     ! Solve the semi-implicit diffusion/advection PDEs for the
     ! biology quantities.
     call solve_bio_eqns(grid%M, P%micro, P%nano, Z, N%O, N%H, Si, &
          D%DON, D%PON, D%refr, D%bSi, day, time)
     !
     !---------- End of Biology Model ----------

     ! Update boundary conditions at surface
     P%micro(0) = P%micro(1)
     P%nano(0) = P%nano(1)
     Z(0) = Z(1)
     N%O(0) = N%O(1)
     N%H(0) = N%H(1)
     Si(0) = Si(1)
     D%DON(0) = D%DON(1)
     D%PON(0) = D%PON(1)
     D%refr(0) = D%refr(1)
     D%bSi(0) = D%bSi(1)

     !--------bottom boundaries--------------------------
     ! Update boundary conditions at bottom of grid
     !
     ! For those variables that we have data, use the annual fit
     ! calculated from the data
     call bot_bound_time(day, day_time, &                             ! in
          T%new(grid%M+1), S%new(grid%M+1), N%O(grid%M+1), &          ! out
          Si(grid%M+1), P%micro(grid%M+1), P%nano(grid%M+1)) ! out
     ! For those variables that we have no data for, assume uniform at
     ! bottom of domain
     call bot_bound_uniform(grid%M, Z, N%H, D%DON, D%PON, D%refr, D%bSi)

     ! Write standard time series results
     ! !!! Please don't change this argument list without good reason. !!!
     ! !!! If it is changed, the change should be committed to CVS.    !!!
     ! !!! For exploratory, debugging, etc. output use                 !!!
     ! !!! write_user_timeseries() below.                              !!!
     call write_std_timeseries(time / 3600., grid,                    &
       ! Variables for standard physics model output
       count, h%new, U%new, V%new, T%new, S%new,                      &
       ! Variables for standard biology model output
       N%O , N%H, Si, P%micro, P%nano, Z, D%DON, D%PON, D%refr, D%bSi)

     ! Write user-specified time series results
     ! !!! Please don't add arguments to this call.           !!!
     ! !!! Instead put use statements in your local copy of   !!!
     ! !!! write_user_timeseries() in the user_output module. !!!
     call write_user_timeseries(time / 3600., grid)
     
     ! Write standard profiles results
     ! !!! Please don't change this argument list without good reason. !!!
     ! !!! If it is changed, the change should be committed to CVS.    !!!
     ! !!! For exploratory, debugging, etc. output use                 !!!
     ! !!! write_user_timeseries() below.                              !!!
     call write_std_profiles(codeId, datetime_str(runDatetime),       &
          datetime_str(startDatetime), year, day, day_time, dt, grid, &
          T%new, S%new, rho%g, P%micro, P%nano, Z, N%O, N%H, Si,      &
          D%DON, D%PON, D%refr, D%bSi, K%u%all, K%t%all, K%s%all,     &
          I_par, U%new, V%new)

     ! Write user-specified profiles results
     ! !!! Please don't add arguments to this call.           !!!
     ! !!! Instead put use statements in your local copy of   !!!
     ! !!! write_user_profiles() in the user_output module.   !!!
     call write_user_profiles(codeId, datetime_str(runDatetime),      &
          datetime_str(startDatetime), year, day, day_time, dt, grid)

     ! Increment time, calendar date and clock time
     call new_year(day_time, day, year, time, dt, month_o)
     scount = scount + 1
     !*** should compare to 1 m value... fix this to be grid independent
     sumS = sumS + 0.5*(S%new(2)+S%new(3))
     sumSriv = sumSriv + S_riv
  end do  !--------- End of time loop ----------

  write (stdout,*) "For Ft tuning"
  write (stdout,*) "Average SSS should be", sumSriv/float(scount)
  write (stdout,*) "Average SSS was", sumS/float(scount)

  ! Close output files
  call timeseries_output_close
  call profiles_output_close

  ! Deallocate memory
  call dalloc_grid
  call dalloc_core_variables
  call dalloc_physics_variables
  call dalloc_biology_variables
  call dalloc_IMEX_variables

end program SOG
