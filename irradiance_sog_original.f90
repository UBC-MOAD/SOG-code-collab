SUBROUTINE irradiance_sog(clouds,cf,t_ime,Jday, In, Ipar, d, waters, I_k, Qs,&
     euphotic)
!clouds,cf,t_ime,Jday, In, Ipar, d, waters, w_type, I_k, Qs,& euphotic)

      USE mean_param
      USE surface_forcing

      IMPLICIT NONE

      TYPE(gr_d), INTENT(IN)::d
      TYPE(okta), INTENT(IN) :: clouds
      REAL, INTENT(IN) :: cf !cloud_fraction  (random variable)
      TYPE(light), INTENT(IN)::waters
 !     INTEGER, INTENT(IN)::w_type  
      INTEGER, INTENT(OUT)::I_k 
      DOUBLE PRECISION, INTENT(OUT)::Qs   !integrated daily solar contribution to
                                              !the heat flux 
      TYPE(entrain), INTENT(OUT)::euphotic !euphotic%depth, euphotic%i                                         
      DOUBLE PRECISION, DIMENSION(0:d%M), INTENT(OUT)::In, Ipar  
      DOUBLE PRECISION, INTENT(IN) :: t_ime  !day_time
      INTEGER, INTENT(IN)::Jday  ! julian day 

      INTEGER::k, check, j, of                            
      DOUBLE PRECISION::declination, hour, cos_Z, P_avg, P_a,P_b, day_length, hour_angle, &
                        sunrise, sunset, Qso, a, b
      DOUBLE PRECISION:: II, cos_Z_max, IImax      
      DOUBLE PRECISION, DIMENSION(0:d%M)::Iparmax
      DOUBLE PRECISION, DIMENSION(0:d%M+1)::plank

      check = 0
      j = t_ime/3600
      hour = (t_ime/3600.0-12.0)*15.  !degrees
      declination = 23.45*PI/180.0*SIN((284.0+DBLE(Jday))/365.25*2.0*PI)  !radians

      a = SIN(declination)*SIN(Lat) 
      b = COS(declination)*COS(Lat)
      cos_Z = a+b*COS(PI/180.0*hour)      !solar elevation
      hour_angle = ACOS(-(TAN(Lat)*TAN(declination)))  !radians
      day_length = hour_angle/15.0*2.0*180.0/PI !hours
      sunrise = 12.0 - 0.5*day_length  !hours
      sunset = 12.0 + 0.5*day_length   !hours
      cos_Z_max = COS(declination-Lat)  !zenith angle

      Qso = Q_o*(1.0+0.033*COS(DBLE(Jday)/365.25*2.0*PI))*(1.0-albedo) !*(1.0-insol)
      if (cf.ge.9) then
         of = cf-2
      elseif (cf.ge.3) then
         of = cf-1
      else 
         of = cf
      endif

      IImax = Qso*(clouds%type(of)%A + clouds%type(of)%B*cos_Z_max)*cos_Z_max  !!Jeffrey thesis page 124, Dobson and Smith

      IF (t_ime/3600.0 > sunrise .AND. t_ime/3600.0 < sunset) THEN    
          II = Qso*(clouds%type(of)%A + clouds%type(of)%B*cos_Z)*cos_Z   
      ELSE
         II = 0.
      END IF

! so II is the incoming radiation

      Qs = Qso*(((clouds%type(of)%A+clouds%type(of)%B*a)*a+clouds%type(of)%B*b**2.0/2.0)*&
           day_length+&
           (clouds%type(of)%A*b+2.0*clouds%type(of)%B*a*b)*&
           180.0/PI/15.0*(SIN(PI/180.0*(sunset-12.0)*15.0)-SIN(PI/180.0*(sunrise-12.0)*15.0)) + &
           clouds%type(of)%B*b**2.0*&
           180.0/PI/60.0*(SIN(PI/180.0*(sunset-12.0)*30.0)-SIN(PI/180.0*(sunrise-12.0)*30.0)))/24.0 

! Qs is the daily integrated value

!      PRINT "(A)","II"
!      PRINT *,II
!      PRINT "(A)","Qs"
!      PRINT *,Qs

      Ipar = 0.        
      In = 0.          

      Ipar(0) = II*(1.0-waters%R(Jday)) 
      In(0) =  II
!      write (*,*) II,of,cf
      Iparmax(0) = IImax*(1.0-waters%R(Jday))

     DO k = 1, d%M    

        P_avg = 0.

        Ipar(k) = II*(1.0-waters%R(Jday))*EXP(-d%d_i(k)/waters%mu_2(Jday)-&
                      P_avg/mu_3)         !attenuation of light below 10m

        Iparmax(k) = IImax*(1.0-waters%R(Jday))*EXP(-d%d_i(k)/waters%mu_2(Jday)-&
                      P_avg/mu_3)

         IF (Ipar(k) < 0.01*Ipar(0) .AND. check == 0) THEN            
            I_k = k                                          
            check = 1                 
         END IF               
         In(k) = II*waters%R(Jday)*EXP(-d%d_i(k)/waters%mu_1(Jday))+ &
                 Ipar(k)                 !!total attenuation of light = I/Io                                                
     END DO  

     ! splitting up the incoming depending on the penetration

     euphotic%i = 0
     euphotic%g = 0
     euphotic%depth = 0.
     DO k = 1,d%M
!        IF (Iparmax(k) <= Iparmax(0)*0.01) THEN
        IF (Ipar(k) <= 1.4) THEN  ! saturation light intensity for Thalasosira
           euphotic%i = k
           euphotic%g = k
           euphotic%depth = d%d_g(k)
           if (Ipar(k) == 0) euphotic%depth = 0.
           EXIT
        END IF
     END DO




END SUBROUTINE irradiance_sog


