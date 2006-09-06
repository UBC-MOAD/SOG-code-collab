! $Id$
! $Source$

SUBROUTINE surface_flux_sog(mm,ro,w, wt_r, & 
                         salinity_o,salinity_m,salinity_d,temp_o,j_gamma, I,Q_t,alp, Cp_o, &
                         bet, U_ten, V_ten, cf, atemp, humid, Qriver,&
                         stress,rho_fresh_o,day,dtdz,h,upwell_const,upwell,Eriver,u,dt)

      USE mean_param
      USE surface_forcing

      IMPLICIT NONE
      TYPE(height), INTENT(IN)::h
      INTEGER, INTENT(IN)::mm,day
      DOUBLE PRECISION,DIMENSION(0:mm+1), INTENT(IN)::ro ! density%new
      DOUBLE PRECISION, INTENT(IN):: alp, Cp_o, bet !alph%i(0), Cp%i(0), beta%i(0)
      DOUBLE PRECISION, INTENT(IN)::salinity_o, salinity_m, salinity_d,temp_o, & 
           U_ten, V_ten, &
           rho_fresh_o, dtdz, u, dt
           !Q_tot,F_tot, Q_sol,  !U_ten, V_ten are unow, vnow
      REAL, INTENT(IN):: cf,atemp,humid, Qriver, Eriver
      INTEGER, INTENT(IN)::j_gamma
      DOUBLE PRECISION,DIMENSION(0:mm),INTENT(IN)::I
      TYPE(windstress), INTENT(IN OUT)::stress
      TYPE(flux), INTENT(OUT)::w
      real(kind=dp), intent(in) :: upwell_const 
                                 ! upwelling constant, tuned parameter
      DOUBLE PRECISION, INTENT(OUT)::Q_t, &  !Q_t(0)
           wt_r, upwell

      DOUBLE PRECISION::Ft, C_D, UU, rho_atm, Sa, S_riv
      double precision:: r, Ce, sigma, lw_in, lw_out, lw_net
      double precision:: Cs, h_sens, h_latent, h_flux, Cp
      REAL:: epsilon_w,a,b,c,ea,es,cl,le

!           U_ten, V_ten used to find UU for drag coefficient calculation
!           U_ten, V_ten used to find wind stress (with drag coefficient)

!           fresh water flux is F_tot

!           Q_t depends on Q_tot and Q_sol

!           w%u is wind stress/ro
!           w%v is wind stress/ro
!           w%t is Q_t/(ro*cp)
!           w%s is Ft*sal/rf
!           w%b is g,alp,bet* w%t and w%s 

!           wt_r depends on I, ro, cp


         UU = SQRT(U_ten**2.0 + V_ten**2.0)   !note U_ten and V_ten at 22m height

         C_D = 1.0D-03*(2.70/UU + 0.142 + 0.0764*UU)



! other coefficients not used

!----------Wind Stress--------------------------------
      
      !fix rho_atm at 1.25 kg/m^3
      rho_atm = 1.25

      IF (UU /= 0.) THEN  ! UVten /= 0.         
         stress%u%new = U_ten/UU*C_D*rho_atm*(UU)**2.0 ! U_ten/UVten*C_D*rho_atm*(UVten)**2.0 !0.3  !N/m^2
         stress%v%new = V_ten/UU*C_D*rho_atm*(UU)**2.0 !V_ten/UVten*C_D*rho_atm*(UVten)**2.0  !N/m^2
!PRINT*,'UU,U_ten,day,C_D,ro(0)',UU,U_ten,day,C_D,ro(0)
      ELSE
         stress%u%new = 0.0    
         stress%v%new = 0.0    
      END IF 

!----------Salinity------------------------------------

     ! Parameterized fit of the surface salinity of the Strait of
     ! Georgia at station S3 based on the river flows.  (Derived by
     ! Kate Collins 16-Jun-2005)
     S_riv = 29.1166 - Qriver * (0.0019) - Eriver * (0.0392)
     ! Smear the change in salinity over a half day
     Sa = salinity_o + (S_riv - salinity_o) * dt / (0.5 * 24 * 3600) 
     ! Calculate the fresh water flux necessary to get Sa in the top
     ! layer
     Ft = (salinity_o - Sa) / (Sa * dtdz)
     ! ...but the rivers can't add salinity
     if (Ft < 0.) then
        Ft = 0.
     endif

     ! Calculate the entrainment of deep water into the bottom of the
     ! grid based on the parameterization derived by Susan Allen in
     ! Jun-2006 (See entrainment.pdf).
     upwell = upwell_const * (Qriver ** 0.45)


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
! use rho_atm (was rho_a=1.2256) 
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

      w%u(0) = -stress%u%new/ro(0) !
      w%v(0) = -stress%v%new/ro(0) !

      w%t(0) = -Q_t/(ro(0)*Cp_o)  


      w%s(0) = Ft*Sa
      w%b(0) = g*(alp*w%t(0)-bet*w%s(0))


   !!!Radiative contribution to surface heat flux!!! this is equal to zero because j_gamma=0 because it is never defined as anything

 ! need to think about this  eq.A4 in Large, used in def_gamma

      wt_r = -(I(0)/(ro(0)*Cp_o)-  &               
                 I(j_gamma)/(ro(j_gamma)*Cp_o))

!PRINT*,'wtr',wt_r, j_gamma
!pause

END SUBROUTINE surface_flux_sog


