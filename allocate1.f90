! $Id$
! $Source$

subroutine allocate1(M, STAT)

      use mean_param
      use declarations

      implicit none
     
      ! Arguments:
      integer, intent(in) :: M  ! Number of grid points
      integer, dimension(20), intent(out) :: STAT  ! Memory allocation status

      ALLOCATE(Fw(0:M), F_n(0:M), &
           Q_n(0:M), &
           I(0:M), I_par(0:M), &
           omega%s%value(0:M), omega%m%value(0:M), &
           Ri_b(0:M), N_2_g(M), V_t_square(M), &
           STAT = alloc_stat(6)) 
                     
      ALLOCATE(Gvector%s(M),Gvector%t(M),Gvector%u(M),Gvector%v(M),&
               STAT = alloc_stat(8))
               
      ALLOCATE(micro%growth%light(M), micro%growth%new(M), &
           nano%growth%light(M), nano%growth%new(M), &
           micro%Nlimit(M), nano%Nlimit(M), &
           wupwell(0:M), &
           STAT = alloc_stat(10))

      ALLOCATE(f_ratio(M), &
           STAT = alloc_stat(11))

END SUBROUTINE allocate1



