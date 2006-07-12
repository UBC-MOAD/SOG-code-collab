SUBROUTINE define_adv_bio(mm,bio,Gx,dt,P_bio,Wz,dz)

  USE mean_param
  USE surface_forcing

  IMPLICIT NONE
  
  TYPE(gr_d),INTENT(IN)::mm      !grid

  DOUBLE PRECISION, DIMENSION(mm%M),INTENT(IN OUT)::Gx  !Gvector
  DOUBLE PRECISION, DIMENSION(0:mm%M+1),INTENT(IN)::bio  !N%O%new, Detritus(1)%D%new ...
  DOUBLE PRECISION, INTENT(IN)::dt,dz
  DOUBLE PRECISION, DIMENSION(mm%M),INTENT(OUT)::P_bio  !P_no,P_nh,P_di,P_na,P_mi,P_co,P_d1,P_d2
   
  INTEGER::index
  DOUBLE PRECISION, DIMENSION(mm%M+1), INTENT (IN)::Wz

  P_bio = 0.

  !!!for grid 3:  Nitrogen in = bio(4)*W(4)*dt
  !!!!                     out = bio(3)*W(3)*dt
  !!!! New density == (bio(3)*(ispace(3)-W(3)*dt)+bio(4)*W(4)*dt)/(ispace(3)-W(3)*dt+W(4)*dt)
  !!!! Nitrogen loss from same density water flux through the sides == (W(4)-W(3))*dt*(New density)

  DO index = 1,mm%M
        P_bio(index) = Wz(index+1)*bio(index+1) - Wz(index)*bio(index)-(Wz(index+1)-Wz(index))*&      !equation A.4 Jeffrey thesis
             (bio(index)*(mm%i_space(index)-Wz(index)*dt)+bio(index+1)*Wz(index+1)*dt)/&
             (mm%i_space(index)+(Wz(index+1)-Wz(index))*dt)
  END DO
!  write (*,*) 'salt in', P_bio(1)
  Gx(1:mm%M) = Gx(1:mm%M) + dt*P_bio(1:mm%M)/dz


END SUBROUTINE define_adv_bio
