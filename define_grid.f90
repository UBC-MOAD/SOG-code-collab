SUBROUTINE define_grid(gridd, depth, spacing)
          USE mean_param
          
          IMPLICIT NONE

          TYPE(gr_d), INTENT(IN OUT)::gridd !g_space, i_space, d_i, d_g
          DOUBLE PRECISION, INTENT(IN)::depth, spacing !D, lambda
          DOUBLE PRECISION, DIMENSION(0:gridd%M)::xsi_i
          DOUBLE PRECISION, DIMENSION(0:gridd%M+1)::xsi_g
          INTEGER::i, j, k, l

     
          DO i = 0, gridd%M
             xsi_i(i) = DBLE(i)/DBLE(gridd%M)
          END DO
          !!!!KIND(0D0) = 8
          DO j = 1, gridd%M
             xsi_g(j) = (DBLE(j) - 0.5)/DBLE(gridd%M)
          END DO

          xsi_g(0) = xsi_i(0)
          xsi_g(gridd%M+1) = xsi_i(gridd%M)

          IF   (ABS(spacing) < EPSILON(spacing)) THEN
             gridd%d_i = depth*xsi_i
             gridd%d_g = depth*xsi_g
          ELSE 
             gridd%d_i = depth/spacing*LOG(1.0-xsi_i*(1.0-EXP(spacing)))
             gridd%d_g = depth/spacing*LOG(1.0-xsi_g*(1.0-EXP(spacing)))
          END IF

!PRINT*,gridd%d_i
!PRINT*,gridd%d_g
!pause

          DO k = 1, gridd%M
             gridd%i_space(k) = gridd%d_i(k)-gridd%d_i(k-1)
          END DO

          DO l = 0, gridd%M
             gridd%g_space(l) = gridd%d_g(l+1)-gridd%d_g(l)
          END DO

END SUBROUTINE define_grid




