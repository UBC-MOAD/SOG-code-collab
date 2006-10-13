! $Id$
! $Source$

module initial_sog
  ! *** What's in here?
  implicit none

  ! *** These parameter values can probably be set in other, more
  ! *** appropriate modules
  DOUBLE PRECISION, PARAMETER:: &
       Uo = 0.0,     &  ! m/s
       Vo = 0.0,     &  ! m/s
       hm =  2.0, &     !28, & !75.0, & !Large1996
!!$       P_micro = 0.3D-3, &
! *** Parameter value setting of P_nano replaced by a variable version in 
! *** initial_mean below, so that initial value of flagellates biomass may
! *** be set to zero without recompiling
!       P_nano = 2.6D-3 * 0., & !V.flagella.01 add comm. 3.6D-3, &!2.6D-3 , & !7.5D-04 gN/m^3, winter estimate
       NHo = .5D-3 

contains

  subroutine initial_mean (Ui, Vi, T_new, S_new, Pmicro, Pnano, &
       NO, NH, Si, Detritus, &
       hi, ut, vt, &
       pbx, pby, d, D_bins, cruise_id)       
    ! *** What's it do?
    use precision_defs, only: dp
    use grid_mod, only: grid_
    use input_processor, only: getpars
    use mean_param, only: prop, snow
    implicit none
    ! Arguments:
    type(prop), intent(out) :: Ui, Vi
    real(kind=dp), dimension(0:), intent(out) :: &
         T_new,  &  ! Temperature profile
         S_new,  &  ! Salinity profile
         NO,     &  ! Salinity Nitrate
         NH,     &  ! Ammonium profile
         Si,     &  ! Silicon profile
         Pmicro, &  ! Micro phytoplankton profile
         Pnano      ! Nano phytoplankton profile
    integer, intent(in) :: D_bins
    type(snow), dimension(D_bins), intent(inout) :: Detritus
    real(kind=dp), intent(out) :: hi !h%new
    type(prop), intent(out) :: ut, vt
    type(grid_), intent(in) :: d
    real(kind=dp), dimension(d%M), intent(out) :: pbx, pby
    character*4  cruise_id           ! cruise_id

    ! Local variables:
    integer :: i, j      ! loop index
    ! File name to open
    character*80 :: fn           
    ! Place holders for reading CTD data file
    integer :: dum1
    real :: depth, dumc, dumt, dump, dumo
    ! *** P_nano variable replaces parameter version above,
    ! *** so that initial value of flagellates biomass may
    ! *** be set to zero without recompiling
    real(kind=dp) :: P_nano

       !V.flagella.01 add comm. 3.6D-3, &!2.6D-3 , & !7.5D-04 gN/m^3, winter estimate
       P_nano = 2.6D-3

    Ui%new(1) = Uo
    Vi%new(1) = Vo
!!$    Pmicro_new(1) = P_micro
    Pnano(1) = P_nano
    NH(1) = NHo

    !-----detritus loop added march 2006---------------------------------

    open(unit=44, file="input/initial_Detritus.dat", &
         status="OLD", action="READ")
    do i = 1, d%M + 1
       read(44, *) Detritus(1)%D%new(i), Detritus(2)%D%new(i)
    end do
    close(44)

    ! *** The values read here are all overwritten, but removing this code
    ! *** causes unclean diffs (profiles are visually the same though)
    open(unit=49, file="input/NH4.dat", status="OLD", &
         action="READ")
    do i = 1, d%M + 1
       read(49, *) NH(i)
    end do
    close(49)


    !---biology jan 10 2005---------------------------------

    DO i = 2, d%M+1   
       IF (d%d_g(i) <= hm) THEN  !Large1996  March 1960 initial profile        

!!$          Pmicro_new(i) = P_micro
          Pnano(i) = P_nano
          NH(i) = NHo
       ELSE
!!$          Pmicro_new(i) = 0.
          Pnano(i) = 0.
          NH(i) = 0.
       END IF
       Ui%new(i) = Uo
       Vi%new(i) = Vo

    END DO

!!$    Pmicro_new(d%M+1) = 0.
    Pnano(d%M+1) = 0.
    NH(d%M+1) = 0.
    Detritus(1)%D%new(d%M+1) = 0. ! (DON ==> Detritus(1), need some deep ocean value)
    Detritus(2)%D%new(d%M+1) =  0. !Detritus(2)%D%new(d%M) !PON needs a deep ocean value


    ! read in nutrients data
    write (fn,'("../sog-initial/Nuts_",a4,".txt")') cruise_id
    open(unit=66, file=fn, status="OLD", action="READ")

    do i = 0, d%M
       read(66, *) depth, NO(i), Si(i)
       if (depth.ne.i*0.5) then
          write (*,*) 'Expecting nutrients, NO3 and Silicon at 0.5 m intervals'
          stop
       endif
    enddo
    close (66)
    NO(d%M+1) = NO(d%M)
    Si(d%M+1) = Si(d%M)
    


    !---end biology------------------------------------------ 

    ! STRATOGEM CTD data from station S3 to initialize temperature,
    ! phytoplacnkton biomass, and salinity in water column
    fn = getpars("ctd_in")
    open(unit=46, file=fn, status="old")
    ! *** Maybe rework this so we can read orginal (not stripped) CTD
    ! *** data files?
    do i = 1, d%M + 1
       read(46, *) dum1, depth, T_new(i), dumc, Pmicro(i), &
            dumt, dump, dumo, S_new(i)  
       ! *** Maybe we should have a degC2degK function? 
       T_new(i) = T_new(i) + 273.15
       ! *** Does the next line actually do anything?
!!$       Pmicro_new(i) = Pmicro_new(i)
    enddo
    close(46)

    T_new(0) = T_new(1)  !Surface
    S_new(0) = S_new(1)  !Boundary
    Pmicro(0) = Pmicro(1)
    Pnano(0) = Pnano(1)
    NH(0) = NH(1)

    ! assuming dz = 0.5
    ! Interpolate CTD data values for grid points between measurements?
    do i = d%M + 1, 2, -1
       j = i / 2
       if (j * 2 == i) then
          T_new(i) = T_new(j)
          S_new(i) = S_new(j)
          Pmicro(i) = Pmicro(j) 
          Pnano(i) = Pnano(j)
       else
          ! *** This looks like a job for a arith_mean function
          T_new(i) = T_new(j) * 0.5 + T_new(j+1) * 0.5
          S_new(i) = S_new(j) * 0.5 + S_new(j+1) * 0.5
          Pmicro(i) = Pmicro(j) * 0.5 + Pmicro(j+1) * 0.5
          Pnano(i) = Pnano(j) * 0.5 + Pnano(j+1) * 0.5
       endif
    enddo

    hi = hm  !Large1996 !ho 

    !If the bottom fluxes are fixed, use the following 
    !tp reevaluate M+1 values:
    Vi%new(0) = Vi%new(1)  !Conditions
    Ui%new(0) = Ui%new(1)

    !            initialize ut and vt
    do i=1,d%M
       ut%new(i) = 0.
       vt%new(i) = 0.
       pbx(i) = 0.
       pby(i) = 0.
    enddo
  end subroutine initial_mean

end module initial_sog
