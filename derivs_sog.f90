SUBROUTINE derivs_sog(x_var,nvar1,Y_var,DYDX_deriv,Temp)

      USE surface_forcing
      USE declarations
      USE reaction
      USE pdf

      IMPLICIT NONE

      INTEGER, INTENT(IN)::nvar1 !5*M+bin_tot + Csources*M+D_bins*M
      DOUBLE PRECISION, DIMENSION(nvar1), INTENT(IN)::Y_var  !INTENT(IN OUT)
      DOUBLE PRECISION, DIMENSION(nvar1), INTENT(OUT)::DYDX_deriv
      DOUBLE PRECISION, INTENT(IN)::x_var  !Y_var known at x_var

      INTEGER::j_j, k_k,bin_tot2,bin,isusan
      DOUBLE PRECISION, DIMENSION(nvar1)::Yplus
      DOUBLE PRECISION, DIMENSION(1:grid%M)::PO_deriv,pellet_prod!,N_urea_save
      DOUBLE PRECISION, DIMENSION(1:grid%M)::Temp,fratioavg
      DOUBLE PRECISION::bottom,bottom2
      DOUBLE PRECISION, DIMENSION(2*grid%M)::Resp

     
      bin_tot2 = M2
      bin = M2

      DO j_j = 1,nvar1  
         IF (Y_var(j_j) >= 0.) THEN 
            Yplus(j_j) = Y_var(j_j)
         ELSE
            Yplus(j_j) = 0.
         END IF
      END DO

      N%O_uptake%new = 0.
      N%H_uptake%new = 0.
      waste%medium = 0.
      N%remin = 0.


      bottom = 100.  !integrate rate of urea production over top 100 m
      bottom2 = 40.  

      CALL p_growth(Yplus,Yplus(1:M),nvar1,I_par,grid,N,micro,Temp,Resp) !microplankton
!      waste%medium = waste%medium + micro%M_z*Yplus(1:M)
      waste%medium = waste%medium + Resp(M+1:2*M)*Yplus(1:M)

      DYDX_deriv(1:nvar1) = 0.      

      !!!New quantity, bacterial 0xidation of NH to NO pool ==> NH^2
      N%bacteria(1:M) = N%r*Yplus(2*M+1:3*M)**2.0


      DO k_k = 1,D_bins-1
         N%remin(:) = N%remin(:) + Detritus(k_k)%r*Yplus(3*M+(k_k-1)*M+1:3*M+k_k*M)
      END DO
      IF (MINVAL(N%remin) < 0.) THEN
         PRINT "(A)","N%remin < 0. in derivs.f90"
         PRINT *,N%remin
         STOP
      END IF

      DO j_j = 1,M 
         IF (N%O_uptake%new(j_j)+N%H_uptake%new(j_j) > 0.) THEN
            f_ratio(j_j) = N%O_uptake%new(j_j)/(N%O_uptake%new(j_j)+N%H_uptake%new(j_j))

            !fratioavg(j_j) = f_ratio(j_j)*(N%H_uptake%new(j_j)+N%O_uptake%new(j_j))
         ELSE
            f_ratio(j_j) = 0.
         END IF
      END DO


      DO j_j = 1,nvar1
         IF(j_j  <= M .AND. Yplus(j_j) > 0.) THEN !Pmicro mort 
           DYDX_deriv(j_j) = (micro%growth%new(j_j)-Resp(j_j)-Resp(j_j+M))*Yplus(j_j)
!            DYDX_deriv(j_j) = (micro%growth%new(j_j)-Resp(j_j)-micro%M_z)*Yplus(j_j)
         ELSE IF (j_j > M .AND. j_j <= 2*M .AND. Yplus(j_j) > 0.) THEN  !NO
            DYDX_deriv(j_j) = -N%O_uptake%new(j_j-M) +N%bacteria(j_j-M)
         ELSE IF (j_j > 2*M .AND. j_j <= 3*M) THEN  !NH
            DYDX_deriv(j_j) = -N%H_uptake%new(j_j-2*M) + Resp(j_j-2*M) + waste%medium(j_j-2*M)*waste%m%destiny(0)+N%remin(j_j-2*M) - N%bacteria(j_j-2*M) 
         ELSE IF (j_j > 3*M .AND. j_j <= 3*M+D_bins*M) THEN !Detritus
            DO k_k = 1,D_bins-1
               IF (j_j > 3*M+(k_k-1)*M .AND. j_j <= 3*M+k_k*M) THEN
                     DYDX_deriv(j_j) = waste%medium(j_j-(3*M+(k_k-1)*M))*waste%m%destiny(k_k)-Detritus(k_k)%r*Yplus(j_j) 
               END IF
            END DO
        END IF
      END DO


      PO_deriv(1:M) = DYDX_deriv(1:M)

      CALL sum_g(grid,DYDX_deriv(M+1:2*M),bottom2,NO50_rate)  !in top 5 meters
      CALL sum_g(grid,PO_deriv(1:M),bottom2,PO50_rate)          !in top 5 meters

!      CALL sum_g(grid,DYDX_deriv(M+1:2*M),bottom,NO100_rate)  !in top 5 meters
!      CALL sum_g(grid,PO_deriv(1:M),bottom,PO100_rate)          !in top 5 meters

    END SUBROUTINE derivs_sog









