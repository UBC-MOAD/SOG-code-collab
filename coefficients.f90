SUBROUTINE coefficients(alph, beta, Cp, dens, cloud,p_Knox)
      !(alph, beta, Cp, dens, cloud,water,p_Knox,P_q_fraction)
      USE mean_param

      IMPLICIT NONE

      TYPE(constant), INTENT(OUT)::alph, Cp, beta, dens
      TYPE(okta), INTENT(OUT)::cloud
!      TYPE(light), INTENT(OUT)::water
      DOUBLE PRECISION, DIMENSION(23)::p_Knox
!      DOUBLE PRECISION, DIMENSION(12)::P_q_fraction
      
 !!Define thermal expansion coefficient matrix 1/K!!      

      alph%data(1) = alpha(3413.0D-07, 304.15, 35.0, 8.D-07)
      alph%data(2) = alpha(3196.0D-07, 301.15, 35.0, 9.D-07)
      alph%data(3) = alpha(2970.0D-07, 298.15, 35.0, 11.D-07)
      alph%data(4) = alpha(2734.0D-07, 295.15, 35.0, 12.D-07)
      alph%data(5) = alpha(2489.0D-07, 292.15, 35.0, 14.D-07)
      alph%data(6) = alpha(2230.0D-07, 289.15, 35.0, 15.D-07)
      alph%data(7) = alpha(1958.0D-07, 286.15, 35.0, 17.D-07)
      alph%data(8) = alpha(1668.0D-07, 283.15, 35.0, 20.D-07)
      alph%data(9) = alpha(1357.0D-07, 280.15, 35.0, 23.D-07)
      alph%data(10)= alpha(1021.0D-07, 277.15, 35.0, 26.D-07)
      alph%data(11) = alpha(781.0D-07, 275.15, 35.0, 28.D-07)
      alph%data(12) = alpha(526.0D-07, 273.15, 35.0, 31.D-07)
      alph%data(13) = alpha(254.0D-07, 271.15, 35.0, 33.D-07)

  !!Define salinity expansion coefficient matrix 
  !! *actually (partial density/ partial S)_(P,T)!!

      beta%data(1) = alpha(0.749, 304.15, 35.0, 0.)
      beta%data(2) = alpha(0.752, 301.15, 35.0, 0.)
      beta%data(3) = alpha(0.756, 298.15, 35.0, 0.)
      beta%data(4) = alpha(0.760, 295.15, 35.0, 0.)
      beta%data(5) = alpha(0.764, 292.15, 35.0, 0.)
      beta%data(6) = alpha(0.769, 289.15, 35.0, 0.)
      beta%data(7) = alpha(0.775, 286.15, 35.0, 0.)
      beta%data(8) = alpha(0.781, 283.15, 35.0, 0.)
      beta%data(9) = alpha(0.788, 280.15, 35.0, 0.)
      beta%data(10) = alpha(0.796, 277.15, 35.0, 0.)
      beta%data(11) = alpha(0.801, 275.15, 35.0, 0.)
      beta%data(12) = alpha(0.808, 273.15, 35.0, 0.)
      beta%data(13) = alpha(0.814, 271.15, 35.0, 0.)
      
  !!Define heat capacity matrix J/kg/K!!

      Cp%data(1) = alpha(4002.0, 304.15, 35.0,-4.7)
      Cp%data(2) = alpha(4000.0, 301.15, 35.0,-4.8)
      Cp%data(3) = alpha(3998.0, 298.15, 35.0,-4.9)
      Cp%data(4) = alpha(3996.0, 295.15, 35.0,-4.9)
      Cp%data(5) = alpha(3993.0, 292.15, 35.0,-5.1)
      Cp%data(6) = alpha(3991.0, 289.15, 35.0,-5.2)
      Cp%data(7) = alpha(3988.0, 286.15, 35.0,-5.3)
      Cp%data(8) = alpha(3986.0, 283.15, 35.0,-5.5)
      Cp%data(9) = alpha(3985.0, 280.15, 35.0,-5.6)
      Cp%data(10) = alpha(3985.0, 277.15, 35.0,-5.8)
      Cp%data(11) = alpha(3985.0, 275.15, 35.0,-5.9)
      Cp%data(12) = alpha(3987.0, 273.15, 35.0,-6.1)
      Cp%data(13) = alpha(3989.0, 271.15, 35.0,-6.2)

!!Define surface dens%data matrix !!
      
      dens%data(1) = alpha(1021.384, 304.15, 35.0, 0.749)
      dens%data(2) = alpha(1022.397, 301.15, 35.0, 0.752)
      dens%data(3) = alpha(1023.343, 298.15, 35.0, 0.756)
      dens%data(4) = alpha(1024.219, 295.15, 35.0, 0.760)
      dens%data(5) = alpha(1025.022, 292.15, 35.0, 0.764)
      dens%data(6) = alpha(1025.748, 289.15, 35.0, 0.769) 
      dens%data(7) = alpha(1026.394, 286.15, 35.0, 0.775) 
      dens%data(8) = alpha(1026.952, 283.15, 35.0, 0.781)
      dens%data(9) = alpha(1027.419, 280.15, 35.0, 0.788)
      dens%data(10) = alpha(1027.786, 277.15, 35.0, 0.796)
      dens%data(11) = alpha(1027.972, 275.15, 35.0, 0.801)
      dens%data(12) = alpha(1028.106, 273.15, 35.0, 0.808)
      dens%data(13) = alpha(1028.187, 271.15, 35.0, 0.814)

      beta%data%value = beta%data%value/dens%data%value

!!!Define Okta Cloud Model!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! based on Dobson and Smith, table 5
      !cloud%type(0) = cloudy(0.400,0.386,  0) ! NJ: don't use at OSP.  For no cloud use 1 parameters 
                                                 !Large 1996 
!      cloud%type(0) = cloudy(0.517,0.317,0.51) 
!      cloud%type(1) = cloudy(0.517,0.317,0.51)  
!      cloud%type(2) = cloudy(0.474,0.381,0.65)
!      cloud%type(3) = cloudy(0.421,0.413 ,0.73 )
!      cloud%type(4) = cloudy(0.380,0.468, 0.82)
!      cloud%type(5) = cloudy(0.350, 0.457 ,0.91)
!      cloud%type(6) = cloudy(0.304, 0.438, 0.98)
!      cloud%type(7) = cloudy(0.230, 0.384 , 0.99)
!      cloud%type(8) = cloudy(0.106,0.285, 0.84)
!      cloud%type(9) = cloudy(0.134,0.295,0.75)

!!!KC-- new coefficients addded august, 2004. Check on the standard deviation values, but I don't think these are used anywhere, anyway

      cloud%type(0) = cloudy(0.4641,0.3304,0.0945) 
      cloud%type(1) = cloudy(0.4062,0.3799,0.0186)  
      cloud%type(2) = cloudy(0.4129,0.3420,0.0501)
      cloud%type(3) = cloudy(0.4263,0.3212,0.0743)
      cloud%type(4) = cloudy(0.4083,0.3060,0.0723)
      cloud%type(5) = cloudy(0.3360,0.3775,0.0294)
      cloud%type(6) = cloudy(0.3448,0.3128,0.0226)
      cloud%type(7) = cloudy(0.3232,0.3259,0.0019)
      cloud%type(8) = cloudy(0.2835,0.2949,0.0081)
      cloud%type(9) = cloudy(0.1482,0.3384,0.1345)

!!!!!!!!!!!!!!!!Jerlov Water Type!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!Jerlov water type:  1 == I, 2 == IA, 3 == IB, *4 == II*, 5 == III!!!!!

!      water%type(1) = light(0.58,0.35,23.0)
!      water%type(2) = light(0.62,0.60,20.0)
!      water%type(3) = light(0.67,1.0,17.0)
!      water%type(4) = light(0.77,1.5,14.0)
!      water%type(5) = light(0.78,1.4,7.9)

!!!!!!!!!!!!!!!!!Knox piecewise linear function for Precipitation!!!!!!!!!!!!!!!!
      
      p_Knox = (/ 21.772,23.741,22.185,19.151,22.108,18.602,19.082,20.504, &
           17.891,23.349,11.617,11.747,11.499,16.059,13.747,14.257, &
           15.112,17.472,18.455,22.199,17.954,25.277,23.985 /) 

!!!!!!!!!!!!!!!!!!Large 1996 Horizontal advection Fraction each month !!!!!!!!!!used in P_q

!      P_q_fraction = (/ 0.1, 0.1, 0., 0., 0., 0., 0., 0.05, 0.15, 0.1, 0.3, 0.2 /)
!      P_q_fraction = (/ 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0. /)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END SUBROUTINE coefficients
