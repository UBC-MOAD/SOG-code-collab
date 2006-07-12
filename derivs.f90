SUBROUTINE derivs(x_var,nvar1,Y_var, DYDX_deriv)

      USE surface_forcing
      USE declarations
      USE reaction
      USE pdf

      IMPLICIT NONE

      INTEGER, INTENT(IN)::nvar1 !5*M+bin_tot + Csources*M+D_bins*M
      DOUBLE PRECISION, DIMENSION(nvar1), INTENT(IN)::Y_var  !INTENT(IN OUT)
      DOUBLE PRECISION, DIMENSION(nvar1), INTENT(OUT)::DYDX_deriv
      DOUBLE PRECISION, INTENT(IN)::x_var  !Y_var known at x_var
     
      INTEGER::j_j, k_k,bin_tot2,bin
      DOUBLE PRECISION, DIMENSION(nvar1)::Yplus
      DOUBLE PRECISION, DIMENSION(2)::growth_o
      DOUBLE PRECISION::bottom,bottom2
      DOUBLE PRECISION, DIMENSION(1:grid%M)::PO_deriv,pellet_prod,N_urea_save

      bin_tot2 = Cevent(1)%length
      bin = MAXVAL(Cevent%length)
      !!odeint.f uses a fifth order Runge-Kutta

      DO j_j = 1,nvar1  !use only positive definite quantities in reaction equations
         IF (Y_var(j_j) >= 0.) THEN  !small
            Yplus(j_j) = Y_var(j_j)
         ELSE
            Yplus(j_j) = 0.
         END IF
         !Y_var(j_j) = Yplus(j_j)  !Try this
      END DO

      IF (MINVAL(Y_var) < 0.) THEN
     !    PRINT "(A)","Y_var < 0. in derivs.f90: Y_var"
     !    PRINT "(A)","P%micro"
     !    PRINT *,Y_var(1:M)
     !    PRINT "(A)","Z%micro"
     !    PRINT *,Y_var(M+1:2*M)
     !    PRINT "(A)","P%nano"
     !    PRINT *,Y_var(2*M+1:3*M)
     !    PRINT "(A)","NO"
     !    PRINT *,Y_var(3*M+1:4*M)
      !   PRINT "(A)","NH"
      !   PRINT *,Y_var(4*M+1:5*M)
      !   PRINT "(A)","Det1"
      !   PRINT *,Y_var(5*M+bin+M+1:5*M+bin+2*M)
      !   PRINT "(A)","Det2"
      !   PRINT *,Y_var(5*M+bin_tot+2*M+1:5*M+bin_tot+3*M)
      END IF

      N%O_uptake%new = 0.
      N%H_uptake%new = 0.
      N%urea%new = 0.
      waste%small = 0.
      waste%medium = 0.
      waste%large = 0.
      zmicro%growth%new = 0.
      N%remin = 0.
      urea_c = 0.
      urea_f = 0.
    !  NO_rate = 0.
    ! PO_rate = 0.

     ! PRINT "(A)","Yplus(nvar1-12:nvar1) before subroutines"
     ! PRINT *,Yplus(nvar1-12:nvar1)

      bottom = 200.  !integrate rate of urea production over top 100 m
      bottom2 = 50.

      CALL Copepod_growth(Yplus,nvar1,Zoo,species,Cevent,grid,Csources,C_types,waste,N%urea%new,bin,bin_tot)
      pellet_prod = waste%medium

      CALL sum_g(grid,N%urea%new(1:grid%M),h_m%new,urea_c)   !production of urea by copepods in ml
      urea_c = urea_c*3600.*24.  !daily production rate
      CALL Copepod_mortality(Yplus,nvar1,Zoo,species,Cevent,grid,Csources,C_types,waste,bin)
      CALL p_growth(Yplus,Yplus(1:M),nvar1,micro,I_par,grid,N,waste%medium,Zoo(1),Cevent(1))  !microplankton
      N_urea_save = Zoo(1)%eta*micro%mort%new+&
           Zoo(1)%a_Ex*Zoo(1)%molt_wt(0)**Zoo(1)%b_Ex*Cevent(1)%nauplii/80.*micro%M_z
     ! N%urea%new = N%urea%new+Zoo(1)%eta*micro%mort%new+&
     !      Zoo(1)%a_Ex*Zoo(1)%molt_wt(0)**Zoo(1)%b_Ex*Cevent(1)%nauplii/80.*micro%M_z     
      waste%medium = waste%medium + (1.-Zoo(1)%delta)*micro%mort%new +nano%M_z*Yplus(1:M)
      CALL p_growth(Yplus,Yplus(2*M+1:3*M),nvar1,nano,I_par,grid,N,waste%small,Zoo(1),Cevent(1)) 
      !nanoplankton
      waste%small = waste%small + nano%M_z*Yplus(2*M+1:3*M)   !natural mortality for nano
      CALL p_graze(Yplus,nvar1,grid,zmicro,micro,N,waste,Csources,bin_tot,Zoo(1),Cevent(1))
      waste%medium = waste%medium + nano%M_z*Yplus(M+1:2*M)
      CALL sum_g(grid,N%urea%new(1:grid%M),h_m%new,urea_f)
      urea_f = urea_f*3600.*24.
      urea_f = urea_f-urea_c 
      IF (urea_f < 0.) THEN
         PRINT "(A)","Error in derivs.f90.  urea_f < 0.; urea_f,urea_c,day,time,year"
         PRINT *,urea_f,urea_c,day,time,year
         STOP
      END IF
      !contribution of dead plankton excreted by background Mesozooplankton
      N%urea%new(1:M) =  N%urea%new(1:M) + Zoo(1)%eta*zmicro%mort%new(1:M)+N_urea_save
      !PRINT "(A)","Yplus(nvar1-12:nvar1) after subroutines"
      !PRINT *,Yplus(nvar1-12:nvar1)
     ! IF (I_par(1) > 0.) THEN 
     !    STOP
     ! END IF

      DYDX_deriv(1:nvar1) = 0.
      DO j_j = 1,prey+d_prey
      cgraze(j_j,:) = 0.
         DO k_k = 1,Csources  !grazing on prey j_j  due to all Copepod events!
            IF (species(k_k)%Ntot /= 0. .AND. Cevent(k_k)%on /= 0 ) THEN  !Copepods are present
               cgraze(j_j,:) = cgraze(j_j,:) + species(k_k)%graze(j_j,:)
            END IF
         END DO
      END DO

      !!!New quantity, bacterial 0xidation of NH to NO pool ==> NH^2
      N%bacteria(1:M) = N%r*Yplus(4*M+1:5*M)**2.0

      DO k_k = 1,D_bins-1
         N%remin(:) = N%remin(:) + Detritus(k_k)%r*Yplus(5*M+bin_tot+Csources*M+(k_k-1)*M+1:&
              5*M+bin_tot+Csources*M+k_k*M)
      END DO
      IF (MINVAL(N%remin) < 0.) THEN
         PRINT "(A)","N%remin < 0. in derivs.f90"
         PRINT *,N%remin
         STOP
      END IF
      IF (MINVAL(N%urea%new) < 0.) THEN
         PRINT "(A)","N%urea%new < 0. In derivs.f90"
         PRINT *,N%urea%new
         STOP
      END IF

     ! PRINT "(A)","MAXVAL(micro%mort%new),MAXVAL(zmicro%mort%new),nano%M_z,MAXVAL(Pnano)"
     ! PRINT *,MAXVAL(micro%mort%new),MAXVAL(zmicro%mort%new),nano%M_z,MAXVAL(Yplus(2*M+1:3*M))
     ! PRINT "(A)","Cevent%nauplii,micro%M_z,Zoo%molt_wt(1)"
     ! PRINT *,Cevent%nauplii,micro%M_z,Zoo%molt_wt(1)
     ! STOP
      DO j_j = 1,nvar1
         IF (j_j  <= M .AND. Yplus(j_j) > 0.) THEN  !Pmicro  (Diatoms)
            DYDX_deriv(j_j) = (micro%growth%new(j_j)-nano%M_z)*Yplus(j_j)-micro%mort%new(j_j) - &
                 cgraze(1,j_j)!note, diatom mort like nano + background mortality
         ELSE IF (j_j > M .AND. j_j <= 2*M .AND. Yplus(j_j) > 0.) THEN  !Zmicro mort is squared
            DYDX_deriv(j_j) = (zmicro%growth%new(j_j-M)-nano%M_z)*Yplus(j_j)-&
                 zmicro%mort%new(j_j-M) - cgraze(2,j_j-M)
            ! zmicro also have nano + background mortality
         ELSE IF(j_j > 2*M .AND. j_j <= 3*M .AND. Yplus(j_j) > 0.) THEN !Pnano mort 
          ! DYDX_deriv(j_j) = (nano%growth%new(j_j-2*M)-nano%mort%new(j_j-2*M))*&
          !      Yplus(j_j) - zmicro%graze(1,j_j-2*M)*Yplus(j_j-M) !*Yplus(j_j)/3.D-03 
           DYDX_deriv(j_j) = (nano%growth%new(j_j-2*M)-nano%M_z)*Yplus(j_j) - &
                zmicro%graze(1,j_j-2*M)*Yplus(j_j-M) !*Yplus(j_j)/3.D-03 
         ELSE IF (j_j > 3*M .AND. j_j <= 4*M .AND. Yplus(j_j) > 0.) THEN  !NO
            DYDX_deriv(j_j) = -N%O_uptake%new(j_j-3*M) +N%bacteria(j_j-3*M)
         ELSE IF (j_j > 4*M .AND. j_j <= 5*M) THEN  !NH
            DYDX_deriv(j_j) = -N%H_uptake%new(j_j-4*M)+N%urea%new(j_j-4*M)+waste%small(j_j-4*M)*&
                 waste%s%destiny(0)+waste%medium(j_j-4*M)*waste%m%destiny(0)+&
                 waste%large(j_j-4*M)*waste%l%destiny(0) + N%remin(j_j-4*M)  -N%bacteria(j_j-4*M)
         ELSE IF (j_j > 5*M .AND. j_j <= 5*M+bin_tot+Csources*M) THEN ! Copepods
            DO k_k = 1,Csources
               IF (Yplus(j_j) > 0. .AND. Cevent(k_k)%on /= 0) THEN !(Yplus(5*M +bin_tot2+k_k*M) /= 0.) THEN
                  IF (j_j == 5*M+bin_tot2+k_k*M+1)THEN  !beginning of yy=2
                     bin_tot2 = bin_tot2 + Cevent(k_k)%length
                  END IF
                  IF (j_j > 5*M+bin_tot2-Cevent(k_k)%length+(k_k-1)*M .AND. j_j <= 5*M+bin_tot2+k_k*M) THEN 
                     IF (j_j <= 5*M + bin_tot2 + (k_k-1)*M) THEN !species(k_k)%mature_pdf%wt
                        IF (species(k_k)%mature_pdf(j_j-(5*M+bin_tot2-&
                             Cevent(k_k)%length+(k_k-1)*M))%f > small) THEN
                           DYDX_deriv(j_j) = species(k_k)%ingest(j_j-(5*M+bin_tot2-&
                                Cevent(k_k)%length+(k_k-1)*M))*Yplus(j_j) - &
                                species(k_k)%Ex(j_j-(5*M+bin_tot2-Cevent(k_k)%length+(k_k-1)*M)) 
                        END IF
                     ELSE IF (j_j <= 5*M+bin_tot2+k_k*M) THEN  !species(k_k)%Z%new(j_j)
                        DYDX_deriv(j_j) = -species(k_k)%mort(j_j-(5*M+bin_tot2+(k_k-1)*(M+1)))
                     END IF
                  END IF
               END IF
            END DO
         ELSE IF (j_j > 5*M+bin_tot+Csources*M .AND. j_j <= 5*M+bin_tot+Csources*M+D_bins*M) THEN !Detritus
            DO k_k = 1,D_bins
               IF (j_j > 5*M+bin_tot+Csources*M+(k_k-1)*M .AND. j_j <= 5*M+bin_tot+Csources*M+k_k*M) THEN
                  IF (k_k < D_bins) THEN
                     DYDX_deriv(j_j) =  waste%small(j_j-(5*M+bin_tot+Csources*M+(k_k-1)*M))*&
                          waste%s%destiny(k_k) + &
                          waste%medium(j_j-(5*M+bin_tot+Csources*M+(k_k-1)*M))*waste%m%destiny(k_k) + &
                          waste%large(j_j-(5*M+bin_tot+Csources*M+(k_k-1)*M))*waste%l%destiny(k_k) - &
                          zmicro%graze(k_k+1,j_j-(5*M+bin_tot+Csources*M+(k_k-1)*M))*&
                          Yplus(j_j-(5*M+bin_tot+Csources*M+(k_k-1)*M)+M) -&
                          Detritus(k_k)%r*Yplus(j_j) 
                     IF (k_k == 2 .AND. Yplus(j_j) > 0.) THEN 
                        !Grazing by copepods only on PON  (ie. Detritus(2))
                        DYDX_deriv(j_j) = DYDX_deriv(j_j)-&
                             cgraze(3,j_j-(5*M+bin_tot+Csources*M+(k_k-1)*M))
                     END IF
                  ELSE !losses from model  (currently only mortality of copepods)
                     DYDX_deriv(j_j) = waste%small(j_j-(5*M+bin_tot+Csources*M+(k_k-1)*M))*&
                          waste%s%destiny(k_k) + &
                         waste%medium(j_j-(5*M+bin_tot+Csources*M+(k_k-1)*M))*waste%m%destiny(k_k) + &
                         waste%large(j_j-(5*M+bin_tot+Csources*M+(k_k-1)*M))*waste%l%destiny(k_k)
                     IF (DYDX_deriv(j_j) < 0.) THEN
                        PRINT "(A)","losses from the model are negative in derivs.f90: DYDX_deriv(j_j)"
                        PRINT *,DYDX_deriv(j_j),j_j
                        STOP
                     END IF
                  END IF
               END IF
            END DO
        ! ELSE IF (j_j == 5*M+Csources*(bin+M)+D_bins*M+1) THEN    !micro%Q
        !    DYDX_deriv(j_j) = micro%dlnQ_dt*Yplus(j_j)
        ! ELSE IF (j_j == 5*M+Csources*(bin+M)+D_bins*M+2) THEN    !nano%Q
        !    DYDX_deriv(j_j) = nano%dlnQ_dt*Yplus(j_j)
         END IF
      END DO

      
     ! IF (Cevent(1)%on /= 0) THEN
     !    IF (species(1)%mature_pdf(1)%f > 0. .OR. species(1)%mature_pdf(2)%f > 0.) THEN
     !       growth_o(1)= species(1)%ingest(1)/Yplus(5*M+1)**(Zoo(Cevent(1)%type)%b_Ex-1.)
     !       growth_o(2) = species(1)%ingest(2)/Yplus(5*M+2)**(Zoo(Cevent(1)%type)%b_Ex-1.)
     !       PRINT "(A)","wt1,wt2"
     !       PRINT *,Yplus(5*M+1),Yplus(5*M+2)
     !       PRINT "(A)","Copepod growth (weight independent portion for wt 1 and wt2"
     !       PRINT *,growth_o
     !       PRINT "(A)","Excretion (weight independent portion)  a_Ex"
     !       PRINT *,Zoo(Cevent(1)%type)%a_Ex
     !    END IF
     ! END IF

      
    !  IF (time >= 2984400.) THEN
    !     PRINT "(A)","time in derivs"
    !     PRINT *,time
    !     PRINT "(A)","DYDX_P%micro"
    !     PRINT *,DYDX_deriv(1:M)
    !     PRINT "(A)","DYDX_Z%micro"
    !     PRINT *,DYDX_deriv(M+1:2*M)
    !     PRINT "(A)","DYDX_P%nano"
    !     PRINT *,DYDX_deriv(2*M+1:3*M)
    !     PRINT "(A)","DYDX_NO"
    !     PRINT *,DYDX_deriv(3*M+1:4*M)
    !     PRINT "(A)","DYDX_NH"
    !     PRINT *,DYDX_deriv(4*M+1:5*M)
    !     PRINT "(A)","DYDX_wt"
    !     PRINT *,DYDX_deriv(5*M+1:5*M+bin_tot)
    !     PRINT "(A)","DYDX_Zcopepod"
    !     PRINT *,DYDX_deriv(5*M+bin_tot+1:6*M+bin_tot)
    !     PRINT "(A)","DYDX_D1"
    !     PRINT *,DYDX_deriv(6*M+bin_tot+1:7*M+bin_tot)
    !     PRINT "(A)","DYDX_D2"
    !     PRINT *,DYDX_deriv(7*M+bin_tot+1:8*M+bin_tot)
    !     PRINT "(A)","DYDX_D3"
    !     PRINT *,DYDX_deriv(8*M+bin_tot+1:9*M+bin_tot)
    !  END IF
    !     PRINT "(A)","DYDX_deriv(5*M+1:5*M+Cevent(1)%length) or dydx_wt"
    !     PRINT *,DYDX_deriv(5*M+1:5*M+Cevent(1)%length)
    !     PRINT "(A)","DYDX_deriv(5*M+Cevent(1)%length+1:5*M+Cevent(1)%length+M) or dydx_Znew"
    !     PRINT *,DYDX_deriv(5*M+Cevent(1)%length+1:5*M+Cevent(1)%length+M)
    !     PRINT "(A)","DYDX_deriv(5*M+Cevent(1)%length+3*M+1:5*M+Cevent(1)%length+4*M) or dydx_D3"
    !     PRINT *,DYDX_deriv(5*M+Cevent(1)%length+3*M+1:5*M+Cevent(1)%length+4*M)
    !  END IF
!      PRINT "(A)","nvar1"
!      PRINT *,nvar1
!      PRINT "(A)","Yplus(nvar1-12:nvar1)"
!      PRINT *,Yplus(nvar1-12:nvar1)
!      PRINT "(A)","Y_var(nvar1-12:nvar1)"
!      PRINT *,Y_var(nvar1-12:nvar1)
!      STOP
      PO_deriv(1:M) = DYDX_deriv(1:M)+DYDX_deriv(M+1:2*M)+DYDX_deriv(2*M+1:3*M)+&
           DYDX_deriv(7*M+bin_tot+1:8*M+bin_tot)

      CALL sum_g(grid,DYDX_deriv(3*M+1:4*M),bottom2,NO50_rate)  !in top 50 meters
      CALL sum_g(grid,PO_deriv(1:M),bottom2,PO50_rate)          !in top 50 meters
      CALL sum_g(grid,pellet_prod(1:M),bottom2,feacal50_rate)    !in top 50 meters

      CALL sum_g(grid,DYDX_deriv(3*M+1:4*M),bottom,NO100_rate)  !in top 50 meters
      CALL sum_g(grid,PO_deriv(1:M),bottom,PO100_rate)          !in top 50 meters
      CALL sum_g(grid,pellet_prod(1:M),bottom,feacal100_rate)    !in top 50 meters

      CALL sum_g(grid,pellet_prod(1:M),h_m%new,feacalml_rate)
      feacalml_rate = feacalml_rate*3600.*24.  !daily rate of production in ml

      DO j_j = 1,M 
         IF (N%O_uptake%new(j_j)+N%H_uptake%new(j_j) > 0.) THEN
            f_ratio(j_j) = N%O_uptake%new(j_j)/(N%O_uptake%new(j_j)+N%H_uptake%new(j_j))
      !      PRINT "(A)","f_ratio(j_j),j_j"
      !      PRINT *,f_ratio(j_j),j_j
         ELSE
            f_ratio(j_j) = 0.
         END IF
      END DO


    END SUBROUTINE derivs









