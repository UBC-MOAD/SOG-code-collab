SUBROUTINE Coriolis_and_pg(mm,U,pg,Gx,step)

      USE mean_param
      USE surface_forcing

      IMPLICIT NONE

      TYPE(gr_d),INTENT(IN)::mm
      DOUBLE PRECISION, DIMENSION(0:mm%M+1),INTENT(IN)::U
      DOUBLE PRECISION, DIMENSION(mm%M), INTENT(IN)::pg
      DOUBLE PRECISION, DIMENSION(mm%M), INTENT(OUT)::Gx        
      DOUBLE PRECISION, INTENT(IN)::step  !dt     

      INTEGER::index

      Gx = 0.
      
      DO index = 1, mm%M
         Gx(index) = step*f*U(index) - step*pg(index)
!         write (*,*) index,Gx(index),-step*pg(index)
      END DO

END SUBROUTINE Coriolis_and_pg
