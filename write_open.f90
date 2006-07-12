SUBROUTINE write_open

  USE initial_sog

      character*80	str2	

      open (5,file="infile")

      str2=getpars("writefile1",1)
  
      open (293,file=str2,status="replace",action="write")


      str2=getpars("writefile2",1)
  
      open (295,file=str2,status="replace",action="write")



END SUBROUTINE write_open

!  OPEN(UNIT = 125, FILE = "output/NO_20res.dat",STATUS = "REPLACE", &
!       ACTION = "WRITE")
 ! OPEN(UNIT = 110, FILE = "output/nano_20res.dat",STATUS = "REPLACE", &
 !      ACTION = "WRITE")
!!  OPEN(UNIT = 126, FILE = "output/NH_20res.dat",STATUS = "REPLACE", &
!       ACTION = "WRITE")
!  OPEN(UNIT = 127, FILE = "output/pgrowth_20res.dat",STATUS = "REPLACE", &
!       ACTION = "WRITE")
!  OPEN(UNIT = 119, FILE = "output/NOup_20res.dat",STATUS = "REPLACE", &
!       ACTION = "WRITE")
!  OPEN(UNIT = 121, FILE = "output/NHup_20res.dat",STATUS = "REPLACE", &
!       ACTION = "WRITE")
!  OPEN(UNIT = 123, FILE = "output/fratio_20res.dat",STATUS = "REPLACE", &
!       ACTION = "WRITE")
!  OPEN(UNIT = 135, FILE = "output/PON_fecal_urea_20res.dat",STATUS = "REPLACE", &
!      ACTION = "WRITE")
!


!  OPEN(UNIT = 150, FILE = "output/physical_o_20res.dat",STATUS = "REPLACE", &
!       ACTION = "WRITE")
!  OPEN(UNIT = 152, FILE = "output/dhdt_ml_hourly_20res.dat",STATUS = "REPLACE", &
!       ACTION = "WRITE")
 ! OPEN(UNIT = 153, FILE = "output/Ks_ml_hourly_20res.dat",STATUS = "REPLACE", &
 !      ACTION = "WRITE")
 !! OPEN(UNIT = 154, FILE = "output/Ipar_ml_hourly_20res.dat",STATUS = "REPLACE", &
 !      ACTION = "WRITE")
!  OPEN(UNIT = 155, FILE = "output/time_hourly_20res.dat",STATUS = "REPLACE", &
!       ACTION = "WRITE")
!  OPEN(UNIT = 156, FILE = "output/phys_mar_20res.dat",STATUS = "REPLACE", &
!       ACTION = "WRITE")
!  OPEN(UNIT = 157, FILE = "output/phys_jun_20res.dat",STATUS = "REPLACE", &
!       ACTION = "WRITE") 
!  OPEN(UNIT = 158, FILE = "output/phys_sep_20res.dat",STATUS = "REPLACE", &
!       ACTION = "WRITE")
!  OPEN(UNIT = 159, FILE = "output/phys_dec_20res.dat",STATUS = "REPLACE", &
!       ACTION = "WRITE")
!   OPEN(UNIT = 160, FILE = "output/UV_20res.dat",STATUS = "REPLACE", &
!       ACTION = "WRITE")
!END SUBROUTINE write_open
