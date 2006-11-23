! $Id$
! $Source$

SUBROUTINE surface_flux_sog(mm,ro, wt_r, & 
                         salinity_n,salinity_o,S_riv,temp_o,j_gamma, I,Q_t,alp, Cp_o, &
                         bet, U_ten, V_ten, cf, atemp, humid, Qriver,&
                         stress,&
                         day,dtdz,&
!!$                         h,&
                         upwell_const,upwell,Eriver,u,dt, &
                         Fw_surface, Fw_scale, Ft, &
                         count)
  ! *** Check whether wt_r is needed
  
  use fundamental_constants, only: g
  use turbulence, only: wbar
      USE mean_param
  implicit none
  ! Arguments:
!!$      TYPE(height), INTENT(IN)::h
      INTEGER, INTENT(IN)::mm,day
      DOUBLE PRECISION,DIMENSION(0:mm+1), INTENT(IN)::ro ! density%new
      DOUBLE PRECISION, INTENT(IN):: alp, Cp_o, bet !alph%i(0), Cp%i(0), beta%i(0)
      DOUBLE PRECISION, INTENT(IN)::salinity_n, & ! current upper level Sal
           salinity_o, &  ! previous time step upper level Sal
           temp_o, & 
           U_ten, V_ten, &
           dtdz, &  ! time step divided by mixing layer depth
           u, dt
           !Q_tot,F_tot, Q_sol,  !U_ten, V_ten are unow, vnow
      REAL, INTENT(IN):: cf,atemp,humid, Qriver, Eriver
      INTEGER, INTENT(IN)::j_gamma
      DOUBLE PRECISION,DIMENSION(0:mm),INTENT(IN)::I
      TYPE(windstress), INTENT(IN OUT)::stress
      double precision, intent(out):: S_riv ! salinity goal
      logical, intent(in) :: Fw_surface
      real(kind=dp), intent(in) :: Fw_scale  ! Fresh water scale factor for river flows
      real(kind=dp), intent(out):: Ft  ! fresh water flux
      integer, intent(in) :: count ! iteration count used to stabilize Ft
      real(kind=dp), intent(in) :: upwell_const 
                                 ! upwelling constant, tuned parameter
      DOUBLE PRECISION, INTENT(OUT)::Q_t, &  !Q_t(0)
           wt_r, upwell

      DOUBLE PRECISION:: UU, rho_atm, Sa
      double precision:: r, Ce, sigma, lw_in, lw_out, lw_net
      double precision:: Cs, h_sens, h_latent, h_flux, Cp
      REAL:: epsilon_w,a,b,c,ea,es,cl,le

!           U_ten, V_ten used to find UU for drag coefficient calculation
!           U_ten, V_ten used to find wind stress (with drag coefficient)

!           fresh water flux is F_tot

!           Q_t depends on Q_tot and Q_sol

!           wbar%t is Q_t/(ro*cp)
!           wbar%s is Ft*sal/rf
!           wbar%b is g,alp,bet* wbar%t and wbar%s 

!           wt_r depends on I, ro, cp


         UU = SQRT(U_ten**2.0 + V_ten**2.0)   !note U_ten and V_ten at 22m height

!----------Salinity------------------------------------

     ! Parameterized fit of the surface salinity of the Strait of
     ! Georgia at station S3 based on the river flows.  (Derived by
     ! Kate Collins 16-Jun-2005)  This value is not directly used
     ! in the model but is used to make sure the tuned FT value
     ! is correct.
     S_riv = 29.1166 - Qriver * (0.0019) - Eriver * (0.0392)

     ! tuned fresh water flux value (to give, on average) the parameterized
     ! value above.  
     ! *** need to check linearity over river flows.
     Ft = Fw_scale * (0.0019 * Qriver + 0.0392 * Eriver)

     ! The entrainment of deep water into the bottom of the
     ! grid is based on the parameterization derived by Susan Allen in
     ! Jun-2006 (See entrainment.pdf)
     upwell = upwell_const * (Qriver/2720)**0.41


!net longwave radiation.
!-----------------------------------------------
r=0.03
Ce=9.37e-6
sigma=5.6697e-8  !stefan boltzmann constant
epsilon_w=0.96   !surface emissivity in the IR portion of the spectrum
! cloud fraction cf 
! air temperature atemp 
! water temperature temp_o 

lw_in=(1-r)*(1+0.17*cf**2)*Ce*atemp**2*sigma*atemp**4  
!downward radiation from atmosphere
!there are several different ways to calculate this

lw_out=-epsilon_w*sigma*temp_o**4                      
!upward emission of radiation from earth's surface, stull page 48

lw_net=lw_in + lw_out


!sensible heat flux
!-----------------------------------------------
Cs=1.3e-3 ! sensible heat transfer coefficient
!fix rho_atm at 1.25 kg/m^3
rho_atm = 1.25
Cp=1003 ! specific heat of air
! use UU (was U=metday(ii,8)) ! is this wind speed at the surface?

h_sens = Cs*rho_atm*Cp*UU*(atemp-temp_o) ! power

!latent heat flux
!-----------------------------------------------
!vapour pressure
a = 7.5
b = 237.3
c = 0.7858
ea = (humid/100)*exp(2.303*((a*(atemp-273.15)/((atemp-273.15) + b))+c));  ! mb
!saturated vapour pressure
es = exp(2.3026*((7.5*(temp_o-273.15)/((temp_o-273.15)+237.3))+0.7858))
!latent heat flux
CL=1.3e-3
LE=2.453e6
h_latent=(0.622/1013)*CL*rho_atm*LE*UU*(ea-es) ! latent heat
if (h_latent.lt.0) then
   h_latent=0
endif

!----Total heat flux (Q_t)
!-------------------------
      
h_flux = lw_net+h_sens+h_latent    

      Q_t = h_flux ! W/m2

      
!-----Surface fluxes---------------------  equations A2a-d
      ! Calculate surface flux components (see Large, etal (1994),
      ! eq'ns A2 & A3b)
      !
      ! Temperature (eq'n A2c)
      wbar%t(0) = -Q_t / (ro(0) * Cp_o)
      ! Salinity (eq'n A2d)
      ! Note that fresh water flux is added via Bf in buoyancy.f90
      ! *** Need to check the implications of wbar%s(0)=0 on def_gamma.f90
      if (Fw_surface) then
         wbar%s(0) = Ft * salinity_o
      else
         wbar%s(0) = 0.
      endif
      ! Buoyancy (eq'n A3b)
      wbar%b(0) = g * (alp * wbar%t(0) - bet * wbar%s(0))


   !!!Radiative contribution to surface heat flux!!! this is equal to zero because j_gamma=0 because it is never defined as anything

 ! need to think about this  eq.A4 in Large, used in def_gamma

      wt_r = -(I(0)/(ro(0)*Cp_o)-  &               
                 I(j_gamma)/(ro(j_gamma)*Cp_o))

!PRINT*,'wtr',wt_r, j_gamma
!pause

END SUBROUTINE surface_flux_sog


