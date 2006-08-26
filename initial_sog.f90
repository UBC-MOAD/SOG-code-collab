! $Id$
! $Source$

module initial_sog

  use surface_forcing

  implicit none

  ! *** This appears to be a collection of parameters of the Large, et al KPP
  ! *** model.  Why are they declared here?
  DOUBLE PRECISION, PARAMETER::To = 5. +273.15, & !!Large1996, &!22.000 + 273.15000   
       So = 32.70, & ! Large199635.0 !32.8, &!PSU  !32.70 Large1996
       Uo = 0.0,     & ! m/s
       Vo = 0.0,     &   ! m/s
       ho = 0.0, & !20.0, &  !19.0,     & ! 
       !hm = 80.0, & !(28.0)Sep, Aug (19.0), Oct (38.0), Nov (58.0), Dec(80.0)
  hm =  2.0, & !28, & !75.0, & !Large1996
       hm2 = 2.0, & !Large et al 1994
       To2 = 6.0 + 273.15, & !Large et al 1994
       hs = 100.0, & !Large1996
       hs2 = 150.0, & !Large1996
       !Tm = 6.7, & !(12.8)Sep, Aug (12.9), Oct (10.8), Nov (8.1), Dec(6.7)
  Tm = 4.0 + 273.15, & !Large1996
  Zd = 200., & !Large1996 (m)
       Ss = 33.70, & !Large1996 
       Sd = 33.80, &  !Large1996
       dh_t = 20.0, & !(32.0)Sep, Aug (41.0), Oct (32.0), Nov (32.0), Dec(20.0)
       dT_t = 2.4, & !(6.2)Sep, Aug (6.7), Oct (4.9), Nov (2.8), Dec (2.4)
       h_t = 100.0, & !(60.0)Sep, Aug (60.0), Oct (70.0), Nov (90.0), Dec (100.0)
       Div_T_M = 0.0, & !!Bottom flux  1/s
       Div_S_M = 0.0, &  
       Div_U_M = 0.0, &
       Div_V_M = 0.0, &                             
       P_micro = 0.3D-3, &
! *** Parameter value setting of P_nano replaced by a variable version in 
! *** initial_mean below, so that initial value of flagellates biomass may
! *** be set to zero without recompiling
!       P_nano = 2.6D-3 * 0., & !V.flagella.01 add comm. 3.6D-3, &!2.6D-3 , & !7.5D-04 gN/m^3, winter estimate
       Z_micro = 1.6D-3, &
       Deto =  1.D-3, &
       NHo = .5D-3 

contains

  subroutine initial_mean (Ui, Vi, Ti, Si, Pi, NO, NH, Sil, Detritus, &
       hi, ut, vt, &
       pbx, pby, d, D_bins, cruise_id, flagellates)       

    ! Arguments:
    type(prop), intent(out) :: Ui, Vi, Ti, Si
    type(plankton), intent(out) :: Pi 
    double precision, dimension (0:), intent(out) :: NO  ! N%O%new,  nitrate
    double precision, dimension (0:), intent(out) :: NH  ! N%H%new,  ammonium
    double precision, dimension (0:), intent(out) :: Sil ! Sil%new,  silicon
    integer, intent(in) :: D_bins
    type(snow), dimension(D_bins), intent(inout) :: Detritus
    double precision, intent(out) :: hi !h%new
    type(prop), intent(out) :: ut, vt
    type(gr_d), intent(in) :: d
    double precision, dimension(d%M), intent(out) :: pbx, pby
    character*4  cruise_id           ! cruise_id
    ! *** Temporary flag to turn flagellates model on/off
    logical :: flagellates

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
    double precision :: P_nano

    ! *** Temporary code to allow flagellates biomass to be initialized to zero
    if (flagellates) then
       !V.flagella.01 add comm. 3.6D-3, &!2.6D-3 , & !7.5D-04 gN/m^3, winter estimate
       P_nano = 2.6D-3
    else
       P_nano = 0.
    endif

    Ui%new(1) = Uo
    Vi%new(1) = Vo
    Si%new(1) = So
    Ti%new(1) = To 
    Pi%micro%new(1) = P_micro
    Pi%nano%new(1) = P_nano !V.flagella.01
    NH(1) = NHo

    !-----detritus loop added march 2006---------------------------------

    open(unit=44, file="input/initial_Detritus.dat", &
         status="OLD", action="READ")
    do i = 1, d%M + 1
       read(44, *) Detritus(1)%D%new(i), Detritus(2)%D%new(i)
    end do
    close(44)

    open(unit=49, file="input/NH4.dat", status="OLD", &
         action="READ")
    do i = 1, d%M + 1
       read(49, *) NH(i)
    end do
    close(49)


    !---biology jan 10 2005---------------------------------

    DO i = 2, d%M+1   
       IF (d%d_g(i) <= hm) THEN  !Large1996  March 1960 initial profile        

          Pi%micro%new(i) = P_micro
          Pi%nano%new(i) = P_nano !V.flagella.01
          NH(i) = NHo
       ELSE       !IF (d%d_g(i) <= Zd) THEN
          Pi%micro%new(i) = 0.
          Pi%nano%new(i) = 0. !V.flagella.01
          NH(i) = 0.
       END IF
       Ui%new(i) = Uo
       Vi%new(i) = Vo

    END DO

    Pi%micro%new(d%M+1) = zero !0.
    Pi%nano%new(d%M+1) = zero !0. !V.flagella.01
    NH(d%M+1) = zero !0.
    Detritus(1)%D%new(d%M+1) = zero ! 0. (DON ==> Detritus(1), need some deep ocean value)
    Detritus(2)%D%new(d%M+1) =  zero !Detritus(2)%D%new(d%M) !PON needs a deep ocean value


    ! read in nutrients data
    write (fn,'("../sog-initial/Nuts_",a4,".txt")') cruise_id
    open(unit=66, file=fn, status="OLD", action="READ")

    do i = 0, d%M
       read(66, *) depth, NO(i), Sil(i)
       if (depth.ne.i*0.5) then
          write (*,*) 'Expecting nutrients, NO3 and Silicon at 0.5 m intervals'
          stop
       endif
    enddo
    close (66)
    NO(d%M+1) = NO(d%M)
    Sil(d%M+1) = Sil(d%M)
    


    !---end biology------------------------------------------ 

    ! STRATOGEM CTD data from station S3 to initialize temperature,
    ! phytoplacnkton biomass, and salinity in water column
    fn = getpars("ctd_in", 1)
    open(unit=46, file=fn, status="old")
    ! *** Maybe rework this so we can read orginal (not stripped) CTD
    ! *** data files?
    do i = 1, d%M + 1
       read(46, *) dum1, depth,Ti%new(i), dumc, Pi%micro%new(i), &
            dumt, dump, dumo, Si%new(i)  
       ! *** Maybe we should have a degC2degK function? 
       Ti%new(i) = Ti%new(i) + 273.15
       ! *** Does the next line actually do anything?
       Pi%micro%new(i) = Pi%micro%new(i)
    enddo
    close(46)

    Ti%new(0) = Ti%new(1)  !Surface
    Si%new(0) = Si%new(1)  !Boundary
    Pi%micro%new(0) = Pi%micro%new(1)
    Pi%nano%new(0) = Pi%nano%new(1) !V.flagella.02
    NH(0) = NH(1)

    ! assuming dz = 0.5
    ! Interpolate CTD data values for grid points between measurements?
    do i = d%M + 1, 2, -1
       j = i / 2
       if (j * 2 == i) then
          Ti%new(i) = Ti%new(j)
          Si%new(i) = Si%new(j)
          Pi%micro%new(i) = Pi%micro%new(j) 
          Pi%nano%new(i) = Pi%nano%new(j) !V.flagella.02
       else
          ! *** This looks like a job for a arith_mean function
          Ti%new(i) = Ti%new(j) * 0.5 + Ti%new(j+1) * 0.5
          Si%new(i) = Si%new(j) * 0.5 + Si%new(j+1) * 0.5
          Pi%micro%new(i) = Pi%micro%new(j) * 0.5 + Pi%micro%new(j+1) * 0.5
          Pi%nano%new(i) = Pi%nano%new(j) * 0.5 + Pi%nano%new(j+1) * 0.5 !V.flagella.02
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


  character*80 function getpars(name, flag)
    ! Part of an input processor originated by Joe Tam
    ! *** Move this into a new, real input processor

    ! Arguments:
    character*(*) :: name  ! Name of data item to read
    integer :: flag        ! Echo description and value read if non-zero

    ! Return value:
    character ::  s*80

    ! Local variables:
    character :: label*16, desc*50

    ! Read and validate data item
    read(5, *) label, s, desc
    if (name /= label) then
       print *, "!!! Error: getpars(): ", &
            "Expecting ", name, " but got ", label(1:len_trim(label))
       stop
    end if

    ! Echo description and value, if requested
    if (flag .ne. 0) then
       write(*,'(a50," = ",a)') desc(1:len_trim(desc)), s(1:len_trim(s))
       ! *** Value only returned if flag is non-zero ???
       getpars = s
       return
    end if
  end function getpars

end module initial_sog
