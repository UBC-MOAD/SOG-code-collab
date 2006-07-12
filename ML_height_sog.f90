SUBROUTINE ML_height_sog(mm,Rib, new_h, kmax)

      USE mean_param
      USE surface_forcing

      IMPLICIT NONE

      TYPE(gr_d), INTENT(IN)::mm 
      DOUBLE PRECISION, DIMENSION(mm%M), INTENT(IN)::Rib !Ri_b
      DOUBLE PRECISION, INTENT(OUT)::new_h
      INTEGER, INTENT(OUT)::kmax
      
      INTEGER::k

      kmax = mm%M             
      new_h =mm%d_g(kmax)    
!PRINT*,'1.new-h',new_h,kmax

      DO k = 1, mm%M-2
         IF  (Rib(k) > Ri_c) THEN
          !PRINT*,'Rib(k)',Rib
            kmax = k 
            EXIT
         END IF
      END DO

      IF (kmax == 1) THEN
         new_h = mm%d_g(kmax)*Ri_c/(Rib(kmax)+1.0D-30) !##
!PRINT*,'2.newh',new_h
      ELSE  IF (kmax >= mm%M - 2) THEN
         PRINT "(A)","**ERROR** Mixing too deep!  j_max_g:"
!PRINT*,'kmax',kmax
         kmax = mm%M-3
         new_h = mm%d_g(kmax) - mm%g_space(kmax)/4.0 !#
         PRINT "(A)","Set Kmax = mm%M-3"
!PRINT*,'3.new_h',new_h
      ELSE
!PRINT*,'yes'
         new_h = mm%d_g(kmax) +(Ri_c - Rib(kmax))/(Rib(kmax)-Rib(kmax-1) + &
             1.0D-30)*mm%g_space(kmax-1) ! interpolate

!PRINT*,'4.new-h',new_h,mm%d_g(kmax),Rib(kmax),Ri_c,Rib(kmax-1)
!pause
      END IF

END SUBROUTINE ML_height_sog

      
