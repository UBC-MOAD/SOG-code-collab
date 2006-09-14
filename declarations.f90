! $Id$
! $Source$

module declarations
! *** What's in here?

  use precision_defs, only: dp
  use grid_mod, only: grid

  use mean_param

  implicit NONE

  ! Run parameters read from file specified as 'inputfile' in infile
  integer :: year_o         ! year in which run starts [yr]
  ! *** Why is month_o not declared here; rather it's local in SOG.f90
  integer :: day_o          ! year-day on which run starts [1-Jan = 1]
  real(kind=dp) :: t_o      ! start of run time on day_o [s]
                            ! *** t_o should match CTD profile time of day
  real(kind=dp) :: t_f      ! end of run time [s]
  real(kind=dp) :: dt       ! time step [s]
  character*4 :: cruise_id  ! four number code that labels the start cruise

  TYPE(constant)::alph, beta, dens
  TYPE(prop)::T, S, U, V, B, density, vt, ut
  type(prop):: Sil              ! silicon 
  TYPE(plankton):: P 
  TYPE(nutrient):: N
  TYPE(Knu)::K 
  TYPE(flux)::w
!!$,Bf  
  TYPE(height)::h, surface_h, h_m
  TYPE(MST)::gamma,G_shape
  TYPE(MS)::phi,omega
  TYPE(MSTscalar)::a2,a3  !shape constants
  TYPE(UVST)::Hvector
  type(UVST):: Gvector, Gvector_o     ! vert. adv. & nonlocal and surf/bot fluxes
  type(UVST):: Gvector_c, Gvector_co  ! Coriolis and pg forces
  type(UVST):: Gvector_ao             ! Sinking of particles
  type(UVST):: Gvector_ro             ! Contains effect of biol. model
  TYPE(UVSTmatrix)::Amatrix, Bmatrix, Bmatrix_o
  TYPE(okta)::cloud
  TYPE(windstress)::stress !wind stress (surface_flux.f90)
  TYPE(plankton2)::micro, nano, zmicro
  TYPE(losses)::waste
  TYPE(wind_ecmwf), DIMENSION(:), ALLOCATABLE::wind  
  TYPE(insol_daily), DIMENSION(:), ALLOCATABLE::insol   !2922 data points, about 3 years worth  some of 83, 84, 85, 86
  TYPE(entrain)::euph

  DOUBLE PRECISION, DIMENSION(365)::insol_avg, insol_actual_avg, insol_smooth, insol_actual_smooth
  !                                        avgdata_u, avgdata_v
  DOUBLE PRECISION, DIMENSION(0:366)::smooth_u, smooth_v
  REAL, DIMENSION(1919,24)::cf,atemp,humid
  REAL, DIMENSION(1919)::Qriver
  REAL, DIMENSION(1919)::Eriver
  DOUBLE PRECISION::begin_hour,end_hour  ! used in interpolate_dt
  DOUBLE PRECISION :: time, h_i, del, dummy_time, del_p

  ! Surface buoyancy forcing
  real(kind=dp) :: Bf

  ! Turbulent friction velocity
  real(kind=dp) :: u_star

  ! Fresh water flux quantities:
  ! *** Destined for a module (probably new surface_forcing) eventually
  logical :: Fw_surface  ! Add all of the fresh water on the surface?
  real(kind=dp) :: &
       Ft, &        ! Total fresh water flux
       Fw_scale, &  ! Fresh water scale factor for river flows
       Fw_depth     ! Depth to distribute fresh water flux over
  real(kind=dp), dimension(:), allocatable :: &
       Fw, &  ! Fresh water flux profile
       F_n    ! Fresh water contribution to salinity flux

!!!!For surface fluxes (surface_flux.f90): 
  DOUBLE PRECISION::U_ten,V_ten, & ! U and V velocities (m/s) at standard height (10m or 22m)
       T_atm,& !air temperature at standard height (17m) (K)
       Q_st,Q_atm, & !specific humidity of air in contact with salt water and at standard height (17m)
       Q_tot, F_tot, & !Total turbulent surface heat flux (W/m^2) and freshwater flux
       wt_r
  INTEGER::j_gamma !interface point corresponding to depth at which Radiation contributes to turbulent
  !surface heat flux
  DOUBLE PRECISION::Q_sol  !See irradiance.f90
  DOUBLE PRECISION::Br
  CHARACTER::ignored_input
  INTEGER, DIMENSION(20)::alloc_stat  !***
  INTEGER :: steps, time_step, index,xx, xx2, year, month, jmax_i,    &
       count, yy, yy2, gg, jmaxg, cloud_type, day, stable,  &
       water_type, neg_count
  REAL::D_test
  DOUBLE PRECISION :: w_star,L_star, day_time, Io, h_Ekman  !depth
  !friction vel, convective vel scale
  !monin_obukov length and julian day
  DOUBLE PRECISION::beta_t ! entrainment coefficient under convection
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE::V_t_square, &
       Ri_b, Q_n, N_2_g
  double precision, dimension(:), allocatable:: pbx,pby,dzx,dzy
  !I = intensity, Q_n is the
  !non_turbulent heat flux  (i.e. I)
  DOUBLE PRECISION, DIMENSION(:),  ALLOCATABLE:: I, I_par, &
       &Q_t

  DOUBLE PRECISION, DIMENSION(1)::microQ1_p,nanoQ1_p
  ! Temporary vectors U_i , U_p... to hold results from tridiagonal solver
  real(kind=dp), dimension(:), allocatable :: P1_p, Pnano1_p, &
       NO1_p, NH1_p, Sil1_p, &
       U_p, V_p, S_p, T_p
  real(kind=dp), dimension(:,:), allocatable :: Detritus1_p
  ! Surface biological flux == 0
  DOUBLE PRECISION::pflux_o

  !Variables for printing test functions !

  DOUBLE PRECISION :: avg_T, read_var 

  !Variable for odeint.f,  rkqs.f and derivs.f

  DOUBLE PRECISION::next_time  !see Read_data.f90 and input file 
  ! input/forcing.dat
  INTEGER::M2, M3, N_ok, N_bad      

  !Variable for wind and insol data
  DOUBLE PRECISION::last_wind_u, last_wind_v
  INTEGER::data_point_insol,data_point_ecmwf
  integer :: wind_n
  integer :: insol_n, smooth_x

  DOUBLE PRECISION::vapour_pressure, prain ! (atmosphere at 17m) and current value for &
  !precipitation from p_Knox (Large_param.f90) !
  DOUBLE PRECISION, DIMENSION(23)::p_Knox  !Piecewise linear function for precipitation
  !defined in coefficients.f90 !

  ! Vertical profile of the entrainment velocity
  real(kind=dp), dimension(:), allocatable :: wupwell

  !Detritus variables
  INTEGER::D_bins  !number of detrital compartments  see input/biology.dat
  TYPE(snow), DIMENSION(:), ALLOCATABLE::Detritus  

  DOUBLE PRECISION,DIMENSION(:),ALLOCATABLE::f_ratio

  ! An empty vector
  real(kind=dp), dimension(:), allocatable :: null_vector

  TYPE(bottom_fit), DIMENSION(1659):: ctd_bottom

end module declarations
