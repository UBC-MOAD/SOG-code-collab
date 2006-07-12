SUBROUTINE find_new

  USE surface_forcing
  USE declarations
!  USE Copepod_mod

  IMPLICIT NONE

  EXTERNAL erfc, gammq, node_wt, rtbis, zbrac, zbrak, rtnewt, wt_dwt
!  INTEGER,INTENT(IN)::max_length
  DOUBLE PRECISION, DIMENSION(max_length)::wt_bin_size  !MAXVAL(Cevent%length)
  DOUBLE PRECISION::last_dwt
  REAL::t2
  DOUBLE PRECISION::node_wt_max,node_wt_min,accuracy, rtnewt
  DOUBLE PRECISION::node_wt,rtbis,node_wt_old,node_guess
  DOUBLE PRECISION::fun,dfun
  REAL::erfc,gammq
  LOGICAL::succes
  INTEGER::nb,stepb
  wt_bin_size = 1.0


  !test function node_wt


  DO yy = 1, M

     IF (P1_p(yy) < zero) THEN
        P1_p(yy) = zero
        neg_count = neg_count + 1
     END IF
     IF (NO1_p(yy) < zero) THEN
        NO1_p(yy) = zero
        neg_count = neg_count + 1
     END IF
     IF (NH1_p(yy) < zero) THEN
        NH1_p(yy) = zero
        neg_count = neg_count + 1
     END IF
     DO xx = 1,D_bins
        IF (Detritus1_p(xx,yy) < zero) THEN
           Detritus1_p(xx,yy) = zero
           neg_count = neg_count + 1
        END IF
        Detritus(xx)%D%new(yy) = Detritus1_p(xx,yy)
     END DO

     P%micro%new(yy) = P1_p(yy)
     N%O%new(yy) = NO1_p(yy)
     N%H%new(yy) = NH1_p(yy)
     
  END DO
  
  P%micro%new(0) = P%micro%new(1)
  N%O%new(0) = N%O%new(1)
  N%H%new(0) = N%H%new(1)
  
  P%micro%new(M+1) = P%micro%old(M+1)
  N%O%new(M+1) = N%O%old(M+1)
  N%H%new(M+1) = N%H%old(M+1)
 
 DO xx = 1,D_bins
     Detritus(xx)%D%new(0) = Detritus(xx)%D%new(1)
     Detritus(xx)%D%new(M+1) = Detritus(xx)%D%old(M+1)
 END DO
  
   
END SUBROUTINE find_new
