SUBROUTINE define_Ri_b_sog(d, hh, surf_h, Bb, Uu, Vv, dens, Rib, Vt_sq, N2)

  USE mean_param
  USE surface_forcing

  TYPE(gr_d), INTENT(IN)::d
  TYPE(height), INTENT(IN)::hh 
  TYPE(height), INTENT(OUT)::surf_h
  TYPE(prop), INTENT(IN OUT)::Bb, Uu, Vv, dens
  DOUBLE PRECISION, DIMENSION(d%M), INTENT(IN)::Vt_sq
  DOUBLE PRECISION, DIMENSION(d%M), INTENT(OUT)::Rib
  DOUBLE PRECISION, DIMENSION(0:d%M+1), INTENT(IN)::N2
 
  DOUBLE PRECISION, DIMENSION(0:d%M+1)::test_vector, test_vector2
  INTEGER::yy

  surf_h%new = ep*hh%new
       
  CALL average(d,Uu,surf_h) !test conv  
  CALL average(d,Vv,surf_h) !test conv 
  CALL average(d,dens,surf_h)

  test_vector = (Uu%avg-Uu%new)**2.0+(Vv%avg-Vv%new)**2.0 
  test_vector2 =  -g/dens%new(0)*(dens%avg - dens%new) !Bb%avg - Bb%new  
         
  Rib = 0.0
  Ribmin = 1000.
           
  DO yy = 1, d%M !surf_h%g, d%M 
     IF (N2(yy) >= 0.) THEN
        IF (test_vector2(yy) > EPSILON(test_vector2(yy)) .OR. &
             test_vector2(yy) < -EPSILON(test_vector2(yy))) THEN
           Rib(yy) = test_vector2(yy)*d%d_g(yy)/&
                (test_vector(yy)+Vt_sq(yy) + 1.0D-30)
           if (Rib(yy).lt.Ribmin) then
               Ribmin = Rib(yy)
               ymin = yy
            endif
        END IF
     END IF
  END DO
  
END SUBROUTINE define_Ri_b_sog
