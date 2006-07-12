SUBROUTINE allocate3(STAT)

      USE mean_param
      USE surface_forcing
      USE declarations
!      USE Copepod_mod

      IMPLICIT NONE
          
      INTEGER, DIMENSION(20), INTENT(IN OUT)::STAT

      INTEGER::y_y, x_x


     ALLOCATE(cgraze(prey+d_prey,M))
      
      DO y_y = 1,D_bins
         ALLOCATE(Gvector%d(y_y)%bin(M),Gvector_o%d(y_y)%bin(M),Gvector_o_o%d(y_y)%bin(M),&
              Gvector_ao%d(y_y)%bin(M),Gvector_ao_o%d(y_y)%bin(M),Gvector_ro%d(y_y)%bin(M),&
              Gvector_ro_o%d(y_y)%bin(M),Hvector%d(y_y)%bin(M),&
              Detritus(y_y)%D%new(0:M+1),Detritus(y_y)%D%old(0:M+1),Detritus(y_y)%D%old_old(0:M+1))
      END DO

      
END SUBROUTINE allocate3


