SUBROUTINE surface_flux(mm,ro,w, wt_r, & 
                         time, salinity_o,temp_o,j_gamma, I,Q_t,alp, Cp_o, &
                         bet, U_ten, V_ten, UVten, T_atm, stable, Q_st, Q_atm, day, &
                         cloud_frac,stress,Q_sol, Q_tot,F_tot,vapour_pressure,pressure,prain,rho_fresh_o,&
                         year,stormday,day_time)

      USE mean_param
      USE surface_forcing

      IMPLICIT NONE
      
      INTEGER, INTENT(IN)::mm, day,year,stormday  !grid%M
      DOUBLE PRECISION,DIMENSION(0:mm+1), INTENT(IN)::ro ! density%new
      DOUBLE PRECISION, INTENT(IN)::time, day_time, &
           alp, Cp_o, bet !alph%i(0), Cp%i(0), beta%i(0)
      DOUBLE PRECISION, INTENT(IN)::salinity_o, temp_o, & !S%new(0), T%new(0)
           U_ten, V_ten,T_atm, Q_st, Q_atm, Q_tot,F_tot, Q_sol, &
           cloud_frac,pressure,& !cloud%type(cloud_type)%fraction and Large_data(xx)%P
           vapour_pressure,prain, rho_fresh_o 
      INTEGER, INTENT(IN)::j_gamma, stable
      DOUBLE PRECISION,DIMENSION(0:mm),INTENT(IN)::I
      DOUBLE PRECISION, INTENT(IN OUT)::UVten  !Large_data(xx)%Uten
      TYPE(wind), INTENT(IN OUT)::stress
      TYPE(flux), INTENT(OUT)::w
      DOUBLE PRECISION, INTENT(OUT)::Q_t, &  !Q_t(0)
           wt_r

      DOUBLE PRECISION::Ft, Fs,C_D, C_H, C_E, Qt_sens,&
                        Qt_latent, Qt_lw, evap, precip, Le, omega, UU, Cp_atm, rho_atm

      !IF ( 0.0 <= time .AND. time < T_stress) THEN
      !   stress%u%new = A_stress*SIN(PI*time/T_stress)*SIN(2*PI*time/T_stress)
      !   stress%v%new = A_stress*SIN(PI*time/T_stress)*COS(2*PI*time/T_stress)
      !ELSE IF (time >= T_stress) THEN 
      !   stress%u%new = 0.0   !N/m^2
      !   stress%v%new = 0.0   !N/m^2
      !END IF

      !!!Bulk coefficients!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      IF (year == 1980) THEN
         UU = UVten
      ELSE
         UU = SQRT(U_ten**2.0 + V_ten**2.0)   !note U_ten and V_ten at 22m height
      END IF
      C_D = 1.0D-03*(2.70/UU + 0.142 + 0.0764*UU)
      C_E = 1.0D-03*34.6*SQRT(C_D)

      IF (stable == 1) THEN  !Stable
         C_H = 1.0D-03*18.0*SQRT(C_D)
      ELSE IF (stable == 0) THEN !Unstable
         C_H = 1.0D-03*32.7*SQRT(C_D)
      END IF

!!!!!!!!!!!!!!!!Wind Stress!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      
      IF (year == 1980) THEN
         rho_atm = pressure*100./(R_v*T_atm*(1.0-Q_atm+Q_atm/0.62197)) !***!
      ELSE
         !fix rho_atm at 1.25 kg/m^3
         rho_atm = 1.25
      END IF
 
      IF (UU /= 0.) THEN  ! UVten /= 0.         
         stress%u%new = U_ten/UU*C_D*rho_atm*(UU)**2.0 ! U_ten/UVten*C_D*rho_atm*(UVten)**2.0 !0.3  !N/m^2
         stress%v%new = V_ten/UU*C_D*rho_atm*(UU)**2.0 !V_ten/UVten*C_D*rho_atm*(UVten)**2.0  !N/m^2
      ELSE
         stress%u%new = 0.0    
         stress%v%new = 0.0    
      END IF


!!!!!!!!!!!!!!!!!!!MODELLING A STORM ON DAY "STORMDAY"
      
!      IF (day == stormday .AND. year == 1967) THEN
!         stress%u%new = 0.0 !0.3  !0.05   
!         stress%v%new = 0.0        !!STORM  Large and Crawford [1996]
!                                                  
!         omega = -f      !angular frequency!!
!         IF (day_time >=  0. .AND. day_time <= 86400.0) THEN  !one day run
!            stress%u%new = -1.4*COS((omega + f)*(day_time ))*(SIN(PI*(day_time)/86400.0))**2.0
!            stress%v%new = 1.4*SIN((omega + f)*(day_time))*(SIN(PI*(day_time)/86400.0))**2.0 
!         END IF
!      END IF

!!!!!!!!!!!!!! END MODELLING A STORM !!!!!!!!!!!!!!!!

!!!!!!!!!!!Sensible heat flux!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!         

      Cp_atm = 1004.6*(1.0+0.8375*Q_atm)  !specific heat (See Gill)
      Qt_sens = 0.
      IF (UVten /= 0.) THEN         
         Qt_sens = rho_atm*C_H*Cp_atm*UU*(temp_o - T_atm)  !W/m^2
        ! Qt_sens = rho_atm*C_H*Cp_atm*UVten*(temp_o - T_atm)  !W/m^2
      END IF
!!!!!!!!!!!Latent heat and turbulent fresh water flux!!!!!!!!!!!!!!!!!!!!!!!!!

      evap = 0.
      IF (UVten /= 0.) THEN
     !    IF (Q_st <= Q_atm) THEN !Condensation (ignored)  surface assumed to be saturated
     !       evap = 0.
     !    ELSE
            evap = rho_atm*C_E*UU*( Q_atm -Q_st)
           ! evap = rho_atm*C_E*UVten*( Q_atm -Q_st)
     !    END IF
      END IF
      !precip = (23.2+9.1*COS(2*PI*(day/365.25 - 0.88)))*1.0D-06 !slightly different than Large 1996
      precip = (prain +9.5*COS(2*PI*(day/365.25 - 0.88)))*1.0D-06 !see Large 1996 (A.6)
      Le = 2.5008D+06 - 2.3D+03*(temp_o - 273.16)  !latent heat J/kg (not SST?)

      Qt_latent = Le*evap  !W/m^2
      IF (year == 1980) THEN
         Ft = precip + evap   !kg/m^2/s  !F_tot is from Large1996
      ELSE
         Ft = F_tot
      END IF
!!!!!!!!!!!Longwave radiative loss!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      
      Qt_lw = emiss*Stef_Boltz*T_atm**4.0*(0.39-0.05*SQRT(vapour_pressure))*(1.0 - &
              0.72*cloud_frac) + 4*emiss*Stef_Boltz*T_atm**3.0*(temp_o - T_atm)
              !emiss is negative according to Large 1996  !remove 2.2* fudge factor

!!!!!!!!!Total fluxes  heat (Q_t) and freshwater (Fs)!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      
      IF (year == 1980) THEN
         Q_t = Qt_lw+Qt_latent+Qt_sens !-100.0 W/m^2 or some function of time
         !                !Q_tot is Large1996 and Q_sol is smoothed daily integrated flux due to insolation 
      ELSE
         Q_t = Q_tot-Q_sol
      END IF
      Fs = 0.0  ! no ice conditions
      
!!!!!!!!Surface fluxes!!!!!!!!!!!!!!!

      w%u(0) = -stress%u%new/ro(0) !
      w%v(0) = -stress%v%new/ro(0) !
      w%t(0) = -Q_t/(ro(0)*Cp_o)  
      w%s(0) = Ft*salinity_o/rho_fresh_o + Fs*(salinity_o - &
                    S_ice)/rho_ice
      w%b(0) = g*(alp*w%t(0)-bet*w%s(0))

   !!!Radiative contribution to surface heat flux!!! Zero

      wt_r = -(I(0)/(ro(0)*Cp_o)-  &               
                 I(j_gamma)/(ro(j_gamma)*Cp_o))


END SUBROUTINE surface_flux



