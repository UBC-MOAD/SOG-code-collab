! $Id$
! $Source$

subroutine allocate1(M, STAT)

      use mean_param
      use declarations
      use surface_forcing

      implicit none
     
      ! Arguments:
      integer, intent(in) :: M  ! Number of grid points
      integer, dimension(20), intent(out) :: STAT  ! Memory allocation status

      ALLOCATE(U%new(0:M+1),U%old(0:M+1),&
           V%new(0:M+1),V%old(0:M+1),&
           T%new(0:M+1),T%old(0:M+1),&
           S%new(0:M+1),S%old(0:M+1),&
           B%new(0:M+1),&
           P%micro%old(0:M+1),P%micro%new(0:M+1), &
           P%nano%new(0:M+1),P%nano%old(0:M+1), &
           N%O%new(0:M+1),N%O%old(0:M+1),&
           N%H%new(0:M+1),N%H%old(0:M+1),&
           Sil%new(0:M+1), Sil%old(0:M+1), &            ! silicon
           density%new(0:M+1),&
           STAT = alloc_stat(1))

      ALLOCATE(Detritus(D_bins),Hvector%d(D_bins),Gvector%d(D_bins),Gvector_o%d(D_bins),&
              Gvector_ro%d(D_bins),&
              Gvector_ao%d(D_bins),Detritus1_p(D_bins,M),&
              STAT = alloc_stat(3))

      ALLOCATE(alph%g(0:M+1), alph%i(0:M), alph%idiv(M), &
           beta%idiv(M), beta%g(0:M+1), beta%i(0:M),     &
           w%b(0:M),w%b_err(0:M),Bf%b(0:M),&
               Bf%b_err(0:M),w%t(0:M),w%s(0:M),w%u(0:M),&
               w%v(0:M),w%p%micro(0:M),Amatrix%u%A(M),Amatrix%u%B(M),&
               Amatrix%u%C(M),Amatrix%t%A(M),Amatrix%t%B(M),Amatrix%t%C(M),Amatrix%s%A(M),Amatrix%s%B(M),&
               Amatrix%s%C(M),Bmatrix%u%A(M),Bmatrix%u%B(M),&
               Bmatrix%u%C(M),Bmatrix%t%A(M),Bmatrix%t%B(M),Bmatrix%t%C(M),Bmatrix%s%A(M),Bmatrix%s%B(M),&
               Bmatrix%s%C(M),Bmatrix_o%u%A(M),Bmatrix_o%u%B(M),&
               Bmatrix_o%u%C(M),Bmatrix_o%t%A(M),Bmatrix_o%t%B(M),Bmatrix_o%t%C(M),Bmatrix_o%s%A(M),&
               Bmatrix_o%s%B(M),Bmatrix_o%s%C(M),&
               Amatrix%bio%A(M),Amatrix%bio%B(M),Amatrix%bio%C(M),&
               Bmatrix%bio%A(M),Bmatrix%bio%B(M),Bmatrix%bio%C(M),&
               Bmatrix%null%A(M),Bmatrix%null%B(M),Bmatrix%null%C(M),&
               Amatrix%null%A(M),Amatrix%null%B(M),&
               Bmatrix_o%bio%A(M),Bmatrix_o%bio%B(M),Bmatrix_o%bio%C(M),&
               Amatrix%no%A(M),Amatrix%no%B(M),Amatrix%no%C(M),&
               Bmatrix%no%A(M),Bmatrix%no%B(M),Bmatrix%no%C(M),&
               Bmatrix_o%no%A(M),Bmatrix_o%no%B(M),Bmatrix_o%no%C(M),&
               STAT = alloc_stat(4)) !**&

      ALLOCATE(&
           &T%div_g(M),&
           &U%div_g(M),V%div_g(M),&
           &S%div_g(M),&
           &B%div_g(M),&
           density%div_g(M),&
           T%div_i(M),&
           U%div_i(M), V%div_i(M), S%div_i(M), density%div_i(M),&
           K%u%shear(0:M),K%s%dd(0:M),K%t%dd(0:M),&
           K%u%total(0:M),K%s%total(0:M),K%t%total(0:M),&
           K%t%all(0:M),K%u%all(0:M),K%s%all(0:M),&
           STAT = alloc_stat(5))

      ALLOCATE(ref_T(0:M+1), &
           avg_12(0:M+1), tot_avg(0:M+1), Q_t(0:M), T_To(0:M+1), &
           I(0:M), I_par(0:M), Q_n(0:M), F_n(0:M), &
           phi%m%value(0:M), phi%s%value(0:M), Ri_b(M), &
           N_2_i(M), N_2_g(M), N_2_dens_g(M), &
           Q_test(M), V_t_square(M), omega%s%value(0:M), omega%m%value(0:M), &
           STAT = alloc_stat(6)) 
                     
      ALLOCATE(gamma%m(0:M),gamma%s(0:M),gamma%t(0:M), &
           STAT = alloc_stat(7))

      ALLOCATE(Hvector%s(M),Hvector%t(M),Hvector%u(M),Hvector%v(M),&
           Hvector%p%micro(M),Hvector%p%nano(M),&
               Gvector%s(M),Gvector%t(M),Gvector%u(M),Gvector%v(M),Gvector_o%s(M),Gvector_o%t(M),&
               Gvector_o%u(M),Gvector_o%v(M),&
               Gvector_c%u(M),Gvector_c%v(M),Gvector_co%u(M),&
               Gvector_co%v(M),&
               Gvector%p%micro(M),Gvector_o%p%micro(M),&
               Gvector%p%nano(M),Gvector_o%p%nano(M),&
               Gvector_ao%p%micro(M), &
               Gvector_ro%p%micro(M),&
               Gvector_ro%p%nano(M),&
               Hvector%n%o(M), Hvector%n%h(M), Hvector%sil(M), Gvector%n%o(M), Gvector%n%h(M), &
               Gvector%sil(M), &
               Gvector_o%n%o(M), Gvector_o%n%h(M), Gvector_o%sil(M), &
               Gvector_ro%n%o(M), Gvector_ro%n%h(M), Gvector_ro%sil(M), &
               null_vector(M),STAT = alloc_stat(8)) !**&
               
      ALLOCATE(U_p(M), V_p(M), S_p(M), T_p(M), P1_p(M), Pnano1_p(M), &
           NO1_p(M), NH1_p(M), Sil1_p(M), & 
           wind(wind_n), insol(insol_n), STAT = alloc_stat(9))

      ALLOCATE(micro%growth%light(M),micro%growth%new(M),nano%growth%light(M),nano%growth%new(M),&
           N%O_uptake%new(M),N%H_uptake%new(M),N%urea%new(M),N%remin(M),N%bacteria(M), &
           micro%mort%new(M),  nano%mort%new(M),&
           zmicro%growth%new(M), zmicro%mort%new(M),&
           zmicro%graze(zprey,M),zmicro%q(zprey-1), waste%small(M),waste%medium(M), &
           waste%large(M),waste%s%destiny(0:D_bins),waste%m%destiny(0:D_bins),waste%l%destiny(0:D_bins), &
           wupwell(1:M+1), &
           STAT = alloc_stat(10))

      ALLOCATE(f_ratio(M), &
           STAT = alloc_stat(11))

! sea add
      ALLOCATE (ut%new(M),vt%new(M),pbx(M),pby(M),&
                ut%old(M),vt%old(M),dzx(M),dzy(M), STAT=alloc_stat(12))

               
      
END SUBROUTINE allocate1



