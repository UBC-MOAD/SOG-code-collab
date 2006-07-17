! $Id$
! $Source$

! *** This is a collection of subroutines that should maybe be broken out
! *** into separate files ???

subroutine read_sog
  ! Read in data from various files to initialize the run.
  ! Forcing data: wind, cloud fraction, air temperature, air humidity,
  !               Fraser River flow, Englishman River flow,
  !               salinity at bottom of model domain
  ! *** And some other stuff that I need to figure out.
  ! *** Need to get rid of explicit file name references here.

  use declarations
  use surface_forcing

  implicit none

  ! Local variables:
  integer :: ndays, ic, j, il, jl, iu, ju, para, stn, yr, check
  real :: idu,idl

  external y_jday_t, Julian_day  

  !  First read the wind
  yr=year_o   !this is so the actual year does not change to the last data year read when the data is read

  open(unit = 12, file = "input/wind/SH200123456_rot.dat", &
       status = "OLD", action = "READ")
  !OPEN(UNIT = 12,FILE = "input/wind/wtwice.dat",STATUS = "OLD", &
  !    ACTION = "READ")
  do xx = 1, wind_n   !one year of wind data
     read(12, *) wind(xx)%day, wind(xx)%month, wind(xx)%year, wind(xx)%time, &
          wind(xx)%zonal, wind(xx)%meridional

     DO yy = 1, SIZE(leap_year,1)
        IF (wind(xx)%year == leap_year(yy)) THEN
           wind(xx)%leap = 1
           EXIT
        ELSE IF (yy == SIZE(leap_year,1)) THEN
           wind(xx)%leap = 0
        END IF
     END DO
  END DO
  close(12)

  ! Find Julian day 
  ! *** But it's not really Julian day, just year-day
  CALL Julian_day(wind%leap, wind%month, wind%day, wind%Jday, wind_n) 


  !------second read the cloud fraction---------------------------------------

  open (12,file="input/met/yvr/cf200123456.dat",STATUS = "OLD", &
       ACTION = "READ")
  !  open (12,file="input/met/yvr/c9.dat",STATUS = "OLD", &
  !     ACTION = "READ")

  ndays = 1919

  do ic=1,ndays
     read (12,*) stn,year,month,day,para,(cf(ic,j),j=1,24)
  enddo



  do ic=1,ndays
     do j=1,24
        if (cf(ic,j).lt.-100) then
           idl = 0
           idu = 0
           il = ic
           jl = j
           do while (cf(il,jl).lt.0)
              call stepdown (il,jl)
              idl = idl + 1
              !                  write (*,*) il,jl,cf(il,jl),' cf'
           enddo
           iu = ic
           ju = j
           do while (cf(iu,ju).lt.0)
              call stepup(iu,ju)
              idu = idu + 1
           enddo
           cf(ic,j) = (idu*cf(il,jl)+idl*cf(iu,ju))/(idu+idl)
        endif
     enddo
  enddo

  close (12)


  !---third read the temperature-----------------------------------------

  open (12,file="input/met/yvr/temp200123456.dat",STATUS = "OLD", &
       ACTION = "READ")

  do ic=1,ndays
     read (12,*) stn,year,month,day,para,(atemp(ic,j),j=1,24)
     !         write (*,*) day
  enddo

  do ic=1,ndays
     do j=1,24

        if (atemp(ic,j).lt.-100) then
           idl = 0
           idu = 0
           il = ic
           jl = j
           do while (atemp(il,jl).lt.0)
              call stepdown (il,jl)
              idl = idl + 1
              !                  write (*,*) il,jl,atemp(il,jl)
           enddo
           iu = ic
           ju = j
           do while (atemp(iu,ju).lt.0)
              call stepup(iu,ju)
              idu = idu + 1
           enddo
           atemp(ic,j) = (idu*atemp(il,jl)+idl*atemp(iu,ju))/(idu+idl)
        endif
     enddo
  enddo
  do ic=1,ndays
     do j=1,24
        atemp(ic,j) = atemp(ic,j)/10.+273.15
        !            write (*,*) atemp(ic,j)
     enddo
  enddo
  close (12)

  !---fourth read the humidity--------------------------------------------

  open (12,file="input/met/yvr/hum200123456.dat",STATUS = "OLD", &
       ACTION = "READ")


  do ic=1,ndays
     read (12,*) stn,year,month,day,para,(humid(ic,j),j=1,24)
     !write (*,*) year,month,day,humid(ic,1)
  enddo

  do ic=1,ndays
     do j=1,24
        if (humid(ic,j).lt.-100) then
           idl = 0
           idu = 0
           il = ic
           jl = j
           do while (humid(il,jl).lt.0)
              call stepdown (il,jl)
              idl = idl + 1
              !write (*,*) il,jl,humid(il,jl)
           enddo
           iu = ic
           ju = j
           do while (humid(iu,ju).lt.0)
              call stepup(iu,ju)
              idu = idu + 1
           enddo
           humid(ic,j) = (idu*humid(il,jl)+idl*humid(iu,ju))/(idu+idl)
           !write (*,*) humid(ic,j)
        endif
     enddo
  enddo
  close (12)

  !--------Fraser River Flow---------------------------------------------

  open (12,file="input/rivers/fr200123456.dat",STATUS = "OLD", &
       ACTION = "READ")
  !open (12,file="input/rivers/fshift.dat",STATUS = "OLD", &
  !      ACTION = "READ")

  !      do ic=1,730
  do ic=1,1640
     read (12,*) year,month,day,Qriver(ic)
  enddo

  !--------read Englishman river -----------------------

  open (12,file="input/rivers/eng200123456.dat",STATUS = "OLD", &
       ACTION = "READ")
  do ic=1,1919
     read (12,*) year,month,day,Eriver(ic)
  enddo

  close (12) 

  !---Bottom salinity condition--------------------------

  OPEN(UNIT = 16,FILE = "input/CTD/bottom_200123456.dat", STATUS = "OLD", &
       ACTION = "READ")

  DO xx = 1, 1659
     READ(16,*)ctd_bottom(xx)%sal,ctd_bottom(xx)%temp,ctd_bottom(xx)%P,ctd_bottom(xx)%No,ctd_bottom(xx)%date
  END DO

  !---total light with depth (jerlov)--------------------------

  !OPEN(UNIT = 12,FILE = "input/jerlov_types2.dat", STATUS = "OLD", &
  !     ACTION = "READ")
  !    DO xx = 1, 82
  !       READ(12,*)water%type1(xx),water%type3(xx),water%type5(xx),water%type7(xx),water%type9(xx)
  !    END DO
  !--biology-----------------------

  OPEN(UNIT = 7,FILE = "input/biology_sog.dat", STATUS = "OLD", &  
       ACTION = "READ")
  OPEN(UNIT = 11,FILE = "input/nutrient.dat",STATUS = "OLD", &
       ACTION = "READ")

  READ (7,*)micro%R, micro%sigma, micro%gamma, micro%Rm, micro%M_z, micro%inhib!, &
  !       Csources, C_types,D_bins
  READ (11,*) micro%k,micro%kapa,micro%N_o,micro%N_x,micro%gamma_o,N%r


  !-------read file PAPMD.60 -----------------------------------------  


  OPEN(UNIT = 42, FILE = "input/PAPMD.60",STATUS = "OLD",&
       ACTION = "READ") 

  DO xx = 1,Large_data_size 
     !  DO xx = 1,wind_n
     READ(42, *)Large_data(xx)%ymdh, Large_data(xx)%Uten, Large_data(xx)%theta, Large_data(xx)%SST, &
          Large_data(xx)%Ta, Large_data(xx)%Tw, Large_data(xx)%Td, Large_data(xx)%cf, Large_data(xx)%P 
     CALL y_jday_t(Large_data(xx))  !find year and leap_year, Julian day, time (s) and change units 
  END DO

  !---Large Parameters--------------------------

  OPEN(UNIT = 14,FILE = "input/Large1996_FN.dat",STATUS = "OLD", &
       ACTION = "READ") 
  OPEN(UNIT = 15,FILE = "input/Large1996_QN.dat",STATUS = "OLD",&
       ACTION = "READ")
  DO xx = 1, 14 !data/Large1996_FN.dat
     ! EXIT
     READ(14, *)FN_1996(xx)%day, FN_1996(xx)%data ! nonturbulent freshwater flux
     READ(15, *)QN_1996(xx)%day, QN_1996(xx)%data ! heat flux in mixed layer (as opposed to at surface)
  END DO

  !-----------------------


  OPEN(UNIT = 43, FILE = "input/Detritus.dat",STATUS = "OLD",&
       ACTION = "READ")

  waste%m%destiny = 0.
  DO xx = 1, D_bins-1
     READ(43,*) waste%m%destiny(xx),Detritus(xx)%r,Detritus(xx)%v
  END DO
  waste%m%destiny(0) = 1.-waste%m%destiny(1)-waste%m%destiny(2)

  !-----------------------------------------------

END SUBROUTINE read_sog

!------------------------------------------------------------
subroutine stepdown (il, jl)

  ! Adjust array indices "downward" for interpolation of
  ! missing meteorological data

  implicit none

  integer, intent(inout) :: il, jl

  if (jl == 1) then
     il = il - 1
     jl = 24
  else
     il = il
     jl = jl - 1
  endif

end subroutine stepdown

subroutine stepup(iu, ju)

  ! Adjust array indices "upward" for interpolation of
  ! missing meteorological data

  implicit none

  integer, intent(inout) :: iu, ju

  if (ju == 24) then
     iu = iu + 1
     ju = 1
  else
     iu = iu
     ju = ju + 1
  endif

end subroutine stepup
