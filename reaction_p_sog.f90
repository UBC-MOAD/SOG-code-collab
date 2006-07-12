SUBROUTINE reaction_p_sog(mm,nvar,D_bins,PZ,PZ_det,Pmicro,NO,NH,Detritus,Gvector)

      USE mean_param

      IMPLICIT NONE

      TYPE(gr_d), INTENT(IN)::mm
      INTEGER, INTENT(IN)::nvar, D_bins  !
      DOUBLE PRECISION, DIMENSION(nvar), INTENT(IN)::PZ
      DOUBLE PRECISION, DIMENSION(D_bins*mm%M),INTENT(IN)::PZ_det
      DOUBLE PRECISION, DIMENSION(0:mm%M+1), INTENT(IN)::Pmicro,NO, NH
      TYPE(snow), DIMENSION(D_bins), INTENT(IN)::Detritus
      TYPE(UVST), INTENT(OUT)::Gvector

      INTEGER::j

      Gvector%p%micro = PZ(1:mm%M) - Pmicro(1:mm%M)
      Gvector%N%O = PZ(mm%M+1:2*mm%M) - NO(1:mm%M)
      Gvector%N%H = PZ(2*mm%M+1:3*mm%M) - NH(1:mm%M)

      DO j = 1,D_bins
         Gvector%d(j)%bin = PZ_det((j-1)*mm%M+1:j*mm%M) - Detritus(j)%D%old(1:mm%M)
      END DO

END SUBROUTINE reaction_p_sog
