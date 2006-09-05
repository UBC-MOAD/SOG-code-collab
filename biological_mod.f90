module biological_mod

use declarations, only: D_bins, M2 ! hopefully can get these out of everything else

implicit none

private

public :: PZ_bins, & ! temporary
     derivs_sog, &
     define_PZ, & ! temporary
     reaction_p_sog, &
     init_biology

  TYPE :: bins
!     sequence
     INTEGER :: micro, nano, NO, NH, det, Quant
  END TYPE bins

  type(bins):: PZ_bins ! where quantities (eg. phyto, nitrate) are in PZ

  ! Size of the biology we are using (Quantities and Detritus)
  data PZ_bins%Quant /4/  
  data PZ_bins%micro /1/  ! Position of Diatoms (micro plankton)
  data PZ_bins%nano  /2/  ! Position of Flagellates (nano plankton)
  data PZ_bins%NO    /3/  ! Position of Nitrate
  data PZ_bins%NH    /4/  ! Position of Ammonium
  data PZ_bins%det   /5/  ! Start of detritus
  ! Number of detritus bins, dissolved, slow sink and fast sink
!  data D_bins  /3/ 


contains

  subroutine init_biology (M)

    integer, intent(in) :: M

    D_bins = 3
    M2 = (PZ_bins%Quant+D_bins) * M   !size of PZ in biology: 
! this subroutine should then allocate the PZ array etc
    
  end subroutine init_biology



  SUBROUTINE define_PZ(M, Pmicro, Pnano, NO, NH, Detritus, &
       PZ)

    ! This subroutine takes all the separate variables (microplankton,
    ! nanoplankton, nitrate, ammonium and detritus and loads them
    ! sequentially into the PZ vector for the ODE solver to use.

    USE mean_param, only : snow 


    IMPLICIT NONE

    integer, intent(in):: M
    DOUBLE PRECISION, DIMENSION(0:), INTENT(IN)::Pmicro, Pnano, NO, NH
    TYPE(snow), DIMENSION(D_bins), INTENT(IN)::Detritus
    DOUBLE PRECISION, DIMENSION (M2), INTENT (OUT) :: PZ

    ! Local Variables
    
    INTEGER :: bPZ, ePZ ! start position and end position in PZ array
    INTEGER :: j
    
    PZ = 0.

    ! Microplankton

    bPz = (PZ_bins%micro-1) * M + 1
    ePZ = PZ_bins%micro * M
    PZ(bPZ:ePZ) = Pmicro(1:M)
    
    ! Nanoplankton
    
    bPz = (PZ_bins%nano-1) * M + 1
    ePZ = PZ_bins%nano * M
    PZ(bPZ:ePZ) = Pnano(1:M)

    ! Nitrate

    bPz = (PZ_bins%NO-1) * M + 1
    ePZ = PZ_bins%NO * M
    PZ(bPZ:ePZ) = NO(1:M)

    ! Ammonium

    bPz = (PZ_bins%NH-1) * M + 1
    ePZ = PZ_bins%NH * M
    PZ(bPZ:ePZ) = NH(1:M)

    ! Detritus
    do j=1,D_bins
       bPz = (PZ_bins%Quant + (j-1) ) * M + 1
       ePz = (PZ_bins%Quant + j) * M
       PZ(bPz:ePz) = Detritus(j)%D%new(1:M)
    enddo

  END SUBROUTINE define_PZ

  SUBROUTINE reaction_p_sog(M, PZ, Pmicro, Pnano, NO, &
       NH, Detritus, Gvector)

    USE mean_param, only: snow, UVST
    
    IMPLICIT NONE
    
    integer, intent(in):: M
    DOUBLE PRECISION, DIMENSION(M2), INTENT(IN)::PZ
    DOUBLE PRECISION, DIMENSION(0:M+1), INTENT(IN)::Pmicro,Pnano, NO, NH
    TYPE(snow), DIMENSION(D_bins), INTENT(IN)::Detritus
    TYPE(UVST), INTENT(IN OUT)::Gvector  ! IN only 'cause not setting all

    ! Local variables
    INTEGER :: bPZ, ePZ ! start position and end position in PZ array
    INTEGER::j
    
    ! Micro plankton
    bPZ = (PZ_bins%micro-1) * M + 1
    ePZ = PZ_bins%micro * M
    Gvector%p%micro = PZ(bPZ:ePZ) - Pmicro(1:M)
    
    ! Nano plankton
    bPZ = (PZ_bins%nano-1) * M + 1
    ePZ = PZ_bins%nano * M
    Gvector%p%nano = PZ(bPZ:ePZ) - Pnano(1:M)
    
    ! Nitrate
    bPZ = (PZ_bins%NO-1) * M + 1
    ePZ = PZ_bins%NO * M
    Gvector%N%O = PZ(bPZ:ePZ) - NO(1:M)
    
    ! Ammonimum
    bPZ = (PZ_bins%NH-1) * M + 1
    ePZ = PZ_bins%NH * M
    Gvector%N%H = PZ(bPZ:ePZ) - NH(1:M)

    ! Detritus
    DO j = 1,D_bins
       bPz = (PZ_bins%Quant + (j-1) ) * M + 1
       ePz = (PZ_bins%Quant + j) * M
       Gvector%d(j)%bin = PZ(bPz:ePz) - Detritus(j)%D%old(1:M)
    END DO

  END SUBROUTINE reaction_p_sog

  SUBROUTINE p_growth(M, NO, NH, P, I_par, Temp, N, plank, Resp, Mort) 

! this subroutine calculates the growth, respiration and mortality
! (light limited or nutrient limited) of either phytoplankton class.
! All three are functions of temperature

  use mean_param, only: plankton2, nutrient
  use surface_forcing, only: small

  implicit none

  integer, intent(in) :: M
  ! Nitrate and ammonium concentrations
  double precision, dimension(M), intent(in) :: NO, NH 
  ! plankton concentraton (either Pmicro or Pnano)
  DOUBLE PRECISION, DIMENSION(M), INTENT(IN)::P !V.flagella.01 note: either Pmicro or Pnano
  
  DOUBLE PRECISION, DIMENSION(0:M), INTENT(IN):: I_par  ! light
  DOUBLE PRECISION, DIMENSION(0:M), INTENT(IN):: Temp  ! temperature


  ! in are the parameters of the growth equations, out are the growth values
  TYPE(plankton2), INTENT(IN OUT):: plank ! either micro or nano 
  
  ! nutrient uptake values (incremented)
  TYPE(nutrient), INTENT(IN OUT):: N
  
  DOUBLE PRECISION, DIMENSION(M), intent(out) :: Resp, Mort
    

  ! internal variables

  INTEGER::j ! counter through depth
  
  DOUBLE PRECISION, DIMENSION(M) :: Uc ! growth based on light
  double precision, dimension(M) :: Oup_cell ! NO uptake assuming full light
  double precision, dimension(M) :: Hup_cell ! NH uptake assuming full light
  DOUBLE PRECISION, DIMENSION(M) :: Rmax ! maximum growth rate for given Temp
  double precision :: temp_effect ! temperature limitation
  double precision :: NH_effect ! ammonium effect on nutrient uptake

!!!!!!!!!!!Define growth Due to I_par!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  plank%growth%light = 0.

  Oup_cell = 0.
  Hup_cell = 0.
    
! flagella_jun06/ Next loop is different from V.16.1 => Rmax(replaces plank%R in KPP) 
! is temp dependent in SOG 

  DO j = 1,M

     ! biological process impactedby Q10 temperature effect
     temp_effect = 1.88**(0.1*(Temp(j)-273.15-20.)) 
     
     ! maximum growth rate (before light/nutrient limitation)
     Rmax(j)=plank%R*temp_effect 

     Resp(j)=(plank%Rm)*temp_effect 
     Mort(j)=(plank%M_z)*temp_effect

     plank%growth%light(j) = Rmax(j)*(1.0-EXP(-plank%sigma*I_par(j)/Rmax(j)))
     
     Uc(j) = (1.0-plank%gamma)*plank%growth%light(j)* &
          (1/(1+plank%inhib*I_par(j)))
  END DO

!!!!!!!!!!!!!!!Define growth due to nutrients!!!!!!!!!!!!!!!!!!!!!!!!!

  DO j = 1,M

     NH_effect = exp(-plank%gamma_o * NH(j))
     IF (NO(j) > small) THEN
        Oup_cell(j) = Rmax(j) * NO(j) * plank%kapa * NH_effect / &
             (plank%k + NO(j) * plank%kapa * NH_effect + NH(j)) * &
             (NO(j) + NH(j))**plank%N_x / &
             (plank%N_o + (NO(j) + NH(j))**plank%N_x)
     ELSE
        Oup_cell(j) = 0.
     END IF
     IF (NH(j) > small) THEN
        Hup_cell(j) = Rmax(j) * NH(j) / &
             (plank%k + NO(j) * plank%kapa * NH_effect + NH(j))* &
             (NO(j) + NH(j))**plank%N_x / &
             (plank%N_o + (NO(j) + NH(j))**plank%N_x)
     ELSE
        Hup_cell(j) = 0.
     END IF
     
     IF (Oup_cell(j) < 0.) THEN
        PRINT "(A)","Oup_cell(j) < 0. in reaction.f90"
        PRINT *,Oup_cell(j)
        STOP
     END IF

     IF (Hup_cell(j) < 0.) THEN
        PRINT "(A)","Hup_cell(j) < 0. in reaction.f90"
        PRINT *,Hup_cell(j)
        STOP
     END IF

  END DO

  ! Choose light limitation or nutrient limitation

  DO j = 1,M
     
     IF (Uc(j) >= 0. .AND. Uc(j) < Oup_cell(j) + Hup_cell(j)) THEN 
        
        !LIGHT LIMITING
        plank%growth%new(j) = Uc(j) 
        
        ! split the nutrient uptake between NH and NO
        
        IF (Uc(j) <= Hup_cell(j)) THEN
           ! add to nutrient uptake so we combined the effects of
           ! different phyto classes
           N%H_uptake%new(j) = Uc(j) * P(j) + N%H_uptake%new(j)
           N%O_uptake%new(j) = 0.
        ELSE
           N%H_uptake%new(j) = Hup_cell(j) * P(j) + N%H_uptake%new(j)
           N%O_uptake%new(j) = (Uc(j) - Hup_cell(j)) * P(j) + &
                N%O_uptake%new(j)
        END IF
        
     ELSE IF (Uc(j) >= 0. .AND. Uc(j) >= Oup_cell(j) + Hup_cell(j)) THEN

        !NUTRIENT LIMITING
        plank%growth%new(j) = Oup_cell(j) + Hup_cell(j)  
        
        IF (plank%growth%new(j) < 0.) THEN
           plank%growth%new(j) = 0.
        ENDIF
        
        N%O_uptake%new(j) = Oup_cell(j) * P(j) + N%O_uptake%new(j)
        N%H_uptake%new(j) = Hup_cell(j) * P(j) + N%H_uptake%new(j)
        
     ELSE  !No nutrient uptake, no growth
        plank%growth%new(j) =  0
     END IF
     
  END DO
  
END SUBROUTINE p_growth

subroutine derivs_sog(time, M, M2a, PZ, dPZdt, Temp, I_par)
  ! Calculate the derivatives of the biological model for odeint to
  ! use to advance the biology to the next time step.

  use precision_defs, only: dp
  use mean_param, only: losses, plankton2, nutrient, snow
  use declarations, only: D_bins, N, micro, nano, &
       f_ratio, waste, Detritus

  implicit none

  ! Arguments:
  real(kind=dp), INTENT(IN):: time ! not used
  integer, intent(in) :: M ! grid size
  integer, intent(in) :: M2a ! need to get rid of this
  real(kind=dp), DIMENSION(M2), INTENT(IN):: PZ  ! values
  real(kind=dp), DIMENSION(M2), INTENT(OUT)::dPZdt ! derivatives
  real(kind=dp), dimension(0:M), INTENT(IN):: Temp ! temperature
  real(kind=dp), dimension(0:M), INTENT(IN):: I_par ! light

  ! Local variables

  integer :: ii ! counter through grid ie 1 to M
  integer :: jj ! counter through PZ
  integer :: kk ! counter through detritus

  real(kind=dp), DIMENSION(M):: Resp_micro, Mort_micro! respiration & mortality
  real(kind=dp), DIMENSION(M):: Resp_nano, Mort_nano  ! respiration & mortality

  real(kind=dp), dimension (M) :: Pmicro, Pnano ! micro/nano plankton conc.
  real(kind=dp), dimension (M) :: NO, NH ! nitrate and ammonium conc.
  real(kind=dp), dimension (D_bins, M) :: detr ! detritus

  ! Put PZ micro values into Pmicro variable, removing any negative values
  do ii = 1,M                        ! counter through grid
     jj = (PZ_bins%micro-1) * M + ii ! counter into PZ

     if (PZ(jj)>0) then
        Pmicro(ii) = PZ(jj)
     else
        Pmicro(ii) = 0
     endif
  enddo

! put PZ nano values into Pnano variable, removing any negative values

  do ii = 1,M                        ! counter through grid
     jj = (PZ_bins%nano-1) * M + ii ! counter into PZ

     if (PZ(jj)>0) then
        Pnano(ii) = PZ(jj)
     else
        Pnano(ii) = 0
     endif
  enddo


! put PZ nitrate values into NO variable, removing any negative values

  do ii = 1,M                        ! counter through grid
     jj = (PZ_bins%NO-1) * M + ii    ! counter into PZ
     
     if (PZ(jj)>0) then
        NO(ii) = PZ(jj)
     else
        NO(ii) = 0
     endif
  enddo

! put PZ ammonium values into NH variable, removing any negative values

  do ii = 1,M                        ! counter through grid
     jj = (PZ_bins%NH-1) * M + ii    ! counter into PZ
     
     if (PZ(jj)>0) then
        NH(ii) = PZ(jj)
     else
        NH(ii) = 0
     endif
  enddo

! put PZ detritus values into detr  variable, removing any negative values

  do kk = 1,D_bins
     do ii = 1,M                                 ! counter through grid
        jj = (PZ_bins%Quant + (kk-1)) * M + ii   ! counter through PZ

        if (PZ(jj)>0) then
           detr(kk,ii) = PZ(jj)
        else
           detr(kk,ii) = 0
        endif
     enddo
  enddo

! initialize transfer rates between the pools
  N%O_uptake%new = 0.
  N%H_uptake%new = 0.
  waste%medium = 0.
  waste%small = 0.
  N%remin = 0.

! phytoplankton growth: Nitrate and Ammonimum, conc. of micro plankton
! I_par is light, Temp is temperature 
! N ammonium and nitrate uptake rates (IN and OUT) incremented
! micro is the growth parameters for micro plankton (IN) and the growth rates 
! (OUT)
! Resp is respiration, Mort is mortality

  call p_growth(M, NO, NH, Pmicro, I_par, Temp, & ! in
       N, micro, &                                ! in and out
       Resp_micro, Mort_micro)                    ! out

! put microplankton mortality into the medium detritus flux
  waste%medium = waste%medium + Mort_micro*Pmicro

! phytoplankton growth: Nitrate and Ammonimum, conc. of nano plankton
! I_par is light, Temp is temperature 
! N ammonium and nitrate uptake rates (IN and OUT) incremented
! nano is the growth parameters for micro plankton (IN) and the growth rates 
! (OUT)
! Resp_nano is respiration, Mort_nano is mortality

  call p_growth(M, NO, NH, Pnano, I_par, Temp, & ! in
       N, nano, &                                ! in and out
       Resp_nano, Mort_nano)                     ! out

  waste%small = waste%small + Mort_nano*Pnano

!!!New quantity, bacterial 0xidation of NH to NO pool ==> NH^2
  N%bacteria(1:M) = N%r * NH**2

! remineralization of detritus groups 1 and 2 (not last one)
  do kk = 1,D_bins-1
     N%remin(:) = N%remin(:) + Detritus(kk)%r * detr(kk,:)
  enddo

  IF (MINVAL(N%remin) < 0.) THEN
     PRINT "(A)","N%remin < 0. in derivs.f90"
     PRINT *,N%remin
     STOP
  END IF

! calculate the f-ratio
  do jj = 1,M 
     IF (N%O_uptake%new(jj) + N%H_uptake%new(jj) > 0.) THEN
        f_ratio(jj) = N%O_uptake%new(jj) /  &
             (N%O_uptake%new(jj) + N%H_uptake%new(jj))
     ELSE
        f_ratio(jj) = 0.
     END IF
  END DO

! now put it all together into the derivatives

! initialize derivatives
  dPZdt(1:M2) = 0.      

! microplankton

  do ii = 1,M ! index for Pmicro
     jj = (PZ_bins%micro-1) * M + ii ! index for PZ, derivatives etc

     if (Pmicro(ii) > 0.) then 
        dPZdt(jj) = (micro%growth%new(ii) - Resp_micro(ii) &
             - Mort_micro(ii)) * Pmicro(ii)
! *** problem here, Resp and Mort are multiplied by Pmicro here but not in NH
! *** or detritus derivatives
! okay here, problem is below.
     endif
  enddo

! nanoplankton

  do ii = 1,M ! index for Pnano
     jj = (PZ_bins%nano-1) * M + ii ! index for PZ, derivatives etc

     if (Pnano(ii) > 0.) then 
        dPZdt(jj) = (nano%growth%new(ii) - Resp_nano(ii) &
             - Mort_nano(ii)) * Pnano(ii)
     endif
  enddo


! nitrate

  do ii = 1,M ! index for NO
     jj = (PZ_bins%NO-1) * M + ii ! index for PZ, derivatives etc

     if (NO(ii) > 0.) then  
        dPZdt(jj) = -N%O_uptake%new(ii) + N%bacteria(ii)
     endif
  enddo

! ammonium

  do ii = 1,M ! index for NH
     jj = (PZ_bins%NH-1) * M + ii ! index for PZ, derivatives etc
     
     if (NH(ii) > 0.) then  
        dPZdt(jj) = -N%H_uptake%new(ii) + Resp_micro(ii) + Resp_nano(ii) &
             + waste%medium(ii) * waste%m%destiny(0)                     &
             + waste%small(ii) * waste%s%destiny(0)                      &
             + N%remin(ii) - N%bacteria(ii)
     endif
  enddo

! detritus (not last bin)

  do ii = 1, M
     do kk = 1, D_bins-1
        jj = (PZ_bins%Quant + (kk-1)) * M + ii
        
        if (detr(kk,ii) > 0) then
           dPZdt(jj) = waste%medium(ii) * waste%m%destiny(kk)  &
                + waste%small(ii) * waste%s%destiny(kk) &
                - Detritus(kk)%r * detr(kk,ii) 
        end if
     end do
  end do

end subroutine derivs_sog

end module biological_mod
