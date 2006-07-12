SUBROUTINE write_biological_sog

  USE declarations
  USE surface_forcing

  IMPLICIT NONE

!new parameters

  EXTERNAL bin_wts
  DOUBLE PRECISION::ezone_NOup, ezone_NHup,  ezone_fratio,cohort_b,cohort_f,cohort_mwt,cohort_mwt_o,&
       mixlay_NOup,mixlay_NHup, mixlay_fratio,mixlay_up, mixlay_D2
  DOUBLE PRECISION,DIMENSION(M)::copepodml,NPPml,PONml,DNml,SPNml,tempT,fratioavg,Totalup,ngrowP,dgrowD,&
       D3rate
  DOUBLE PRECISION::mixlay_SPN
  INTEGER::jj,years,sizerun,jj_day,day_count,cohort_cnt

 100  FORMAT(1X, E15.8)
 200  FORMAT(1X, I10)
  
  jj_day = 0

  DO jj = 1,27
     IF (jj == 1 .AND. day < j_day(2)-7) THEN  !j_day(2) = 15
        jj_day = 1
        EXIT
     ELSE IF (jj == 27 .AND. day >= j_day(26)+7) THEN  !j_day(26) = 351
        jj_day = 27
        EXIT
     ELSE IF (day >= j_day(jj)-7 .AND. day < j_day(jj)+7) THEN
        jj_day = jj
        EXIT
     END IF
  END DO
  years = 1
  !years = 9
  sizerun = 27*years  !number of data points interested in saving
  day_count = 0

  DO jj = 1,years  !1960
     IF (jj == year-1966) THEN  !year - 1967
        day_count = (jj-1)*27+jj_day
        EXIT
     END IF
  END DO


  IF (year >= 1967   .AND. year < 1968 .AND. time_step /= steps .AND. jj_day /= 0) THEN
     o_cnt(day_count) = o_cnt(day_count)+1
     NO_o(day_count) = NO_o(day_count)+N%O%new(1)
     ngrow_o(day_count) = ngrow_o(day_count)+nano%growth%new(1)
     dgrow_o(day_count) = dgrow_o(day_count)+micro%growth%new(1)
    ! PRINT "(A)","NO_o(day_count),day_count,N%O%new(1),year"
    ! PRINT *,NO_o(day_count),day_count,N%O%new(1),year
     euph_new(day_count) = euph_new(day_count) + grid%d_g(euph%g)
     hm_new(day_count) = hm_new(day_count) + h_m%new
     euph_g(day_count) = euph%g + euph_g(day_count)
     hm_g(day_count) = h_m%g + hm_g(day_count)
     PONflux200(day_count) = PONflux200(day_count) + PONflux_200   !daily PON flux
     PONflux100(day_count) = PONflux100(day_count) + PONflux_100
     PONfluxml(day_count) = PONfluxml(day_count) + PONflux_ml  
     NOflux200(day_count) = NOflux200(day_count) + NOflux_200
     NOfluxml(day_count) = NOfluxml(day_count) + NOflux_ml
     NOflux100(day_count) = NOflux100(day_count) + NOflux_100     
     NHflux200(day_count) = NHflux200(day_count) + NHflux_200
     NHfluxml(day_count) = NHfluxml(day_count) + NHflux_ml
     NHflux100(day_count) = NHflux100(day_count) + NHflux_100
     feacalml(day_count) =  feacalml(day_count) + feacalml_rate  !daily production of feacal pellets in ml
     ureac(day_count) = ureac(day_count) + urea_c !daily production of urea in ml
     ureaf(day_count) = ureaf(day_count) + urea_f 
     migrateflux(day_count) = migrateflux(day_count) - species(1)%new_wt_out/dt*3600.*24.+&
          species(1)%new_wt_in/dt*3600.*24.  !gN/m^2/day into or out of upper 200m 
     species(1)%new_wt_out = 0.
     species(1)%new_wt_in = 0.
     IF (Cevent(1)%on /= 0.) THEN
        species_cnt2(day_count) = species_cnt2(day_count) + 1
        stage1_no(day_count) = stage1_no(day_count) + species(1)%Ntot*species(1)%stage_f(1)
        stage2_no(day_count) = stage2_no(day_count) + species(1)%Ntot*species(1)%stage_f(2)
        stage3_no(day_count) = stage3_no(day_count) + species(1)%Ntot*species(1)%stage_f(3)
        stage4_no(day_count) = stage4_no(day_count) + species(1)%Ntot*species(1)%stage_f(4)
        stage5_no(day_count) = stage5_no(day_count) + species(1)%Ntot*species(1)%stage_f(5)
        out_no(day_count) = out_no(day_count) + species(1)%n_out
        species_avgwt(day_count) = species_avgwt(day_count) + species(1)%avg_wt
     END IF

     !IF (year == 1967) THEN
      !  IF (Cevent(1)%on /= 0 .AND. day_check /= 0 .AND. day_time <= 43200.0 ) THEN
           !DO jj = 1,nobins+1
           !   IF (day == daybins(jj)) THEN
           
          ! CALL bin_wts  !find bin_cnt 
          ! bin_cnt_year(day_count-54,0:bin_no+2) = bin_cnt(0:bin_no+2)             
           !END IF
           !END DO
      !  END IF
     !END IF
     IF (Cevent(1)%on /= 0 .AND. day_check /= 0 .AND. day_time <= 43200.0 ) THEN  !only once a day when
        ! migration "may" occur !!!!!!!!!
        cohort_cnt = 0
        cohort_b = 0.
        cohort_f = 0.
        cohort_mwt = 0.
        cohort_mwt_o = 0.
        DO jj = 1,species(1)%current_arrival-1
           IF (species(1)%mature_pdf(jj)%f > 0. .AND. species(1)%mature_pdf(jj)%cnout > small2 &
                .AND. species(1)%mature_pdf(jj)%pdf_out <= 1.) THEN  !if some from jj have migrated
              cohort_f = cohort_f + species(1)%mature_pdf(jj)%cnout
              cohort_cnt = cohort_cnt +1
              cohort_b = cohort_b + species(1)%mature_pdf(jj)%b*species(1)%mature_pdf(jj)%cnout
              cohort_mwt = cohort_mwt + species(1)%mature_pdf(jj)%mwt_avg*species(1)%mature_pdf(jj)%cnout
              cohort_mwt_o = cohort_mwt_o+species(1)%mature_pdf(jj)%wt_m*species(1)%mature_pdf(jj)%cnout
              IF (species(1)%mature_pdf(jj)%b <= 0. .OR. species(1)%mature_pdf(jj)%mwt_avg <= 0.) THEN 
                 PRINT "(A)","in write_bio species(1)%mature_pdf(jj)%b <= 0..jj,day,time,species(1)%mature_pdf(jj)%cnout"
                 PRINT *,species(1)%mature_pdf(jj)%b,jj,day,time,species(1)%mature_pdf(jj)%cnout 
                 PRINT "(A)","or species(1)%mature_pdf(jj)%mwt_avg <= 0."
                 PRINT *,species(1)%mature_pdf(jj)%mwt_avg
                 STOP 
              END IF
           END IF
        END DO
        IF (cohort_cnt > 0 .AND. cohort_f > 0.) THEN
           cohort_b = cohort_b/cohort_f   !/DBLE(cohort_cnt)
           cohort_mwt = cohort_mwt/cohort_f  !/DBLE(cohort_cnt)
           cohort_mwt_o = cohort_mwt_o/cohort_f  !/DBLE(cohort_cnt)
        ELSE
           cohort_b = 0.
           cohort_mwt = 0.
           cohort_mwt_o = 0.
        END IF
        IF (cohort_cnt > 0 .AND. cohort_f > 0.) THEN 
           species_cnt(day_count) = species_cnt(day_count) + 1
           species_b(day_count) = species_b(day_count) + cohort_b
           species_mwt(day_count) = species_mwt(day_count) + cohort_mwt
           species_mwt_o(day_count) = species_mwt_o(day_count)+cohort_mwt_o
        END IF
     END IF
!for NPP only calculate integrated values i.e. per m^2 per day

! Mixed layer average quantities
     CALL average2(grid,P%nano%new(1:M),h_m%new,mixlay%nano)
     nano_ml(day_count) = nano_ml(day_count)+mixlay%nano
     CALL average2(grid,P%micro%new(1:M),h_m%new,mixlay%diatom)
     diatom_ml(day_count) = diatom_ml(day_count)+mixlay%diatom
     CALL average2(grid,Z%micro%new(1:M),h_m%new,mixlay%zmicro)
     zmicro_ml(day_count) = zmicro_ml(day_count)+mixlay%zmicro
     copepodml(1:M) = species(1)%Z%new(1:M)*species(1)%avg_wt
     CALL average2(grid,copepodml(1:M),h_m%new,mixlay%copepod)
     copepod_ml(day_count) = copepod_ml(day_count)+mixlay%copepod
     PONml(1:M) = Detritus(2)%D%new(1:M)+Detritus(1)%D%new(1:M)+ P%nano%new(1:M)+P%micro%new(1:M)+&
            Z%micro%new(1:M)
     SPNml(1:M) = Detritus(1)%D%new(1:M)+ P%nano%new(1:M)+Z%micro%new(1:M)
     CALL average2(grid,SPNml(1:M),h_m%new,mixlay_SPN)
     SPN_ml(day_count) = SPN_ml(day_count) + mixlay_SPN
     CALL average2(grid,PONml(1:M),h_m%new,mixlay%PN)
     PN_ml(day_count) = PN_ml(day_count)+mixlay%PN
     DNml(1:M) = Detritus(2)%D%new(1:M)+Detritus(1)%D%new(1:M)
     CALL average2(grid,DNml(1:M),h_m%new,mixlay%DN)
     DN_ml(day_count) = DN_ml(day_count)+mixlay%DN
     CALL average2(grid,N%O%new(1:M),h_m%new,mixlay%NO)
     NO_ml(day_count) = NO_ml(day_count)+mixlay%NO
     CALL average2(grid,N%H%new(1:M),h_m%new,mixlay%NH)
     NH_ml(day_count) = NH_ml(day_count)+mixlay%NH
     NPPml(1:M) = nano%growth%new(1:M)*P%nano%new(1:M)+&
                micro%growth%new(1:M)*P%micro%new(1:M)
     CALL average2(grid,NPPml(1:M),h_m%new,mixlay%NPP)
     NPP_ml(day_count) = NPP_ml(day_count)+mixlay%NPP*h_m%new  !integrated quantity
     ngrowP(1:M) = nano%growth%new(1:M)*P%nano%new(1:M)
     CALL average2(grid,ngrowP(1:M),h_m%new,mixlay%ngrow)
     mixlay%ngrow = mixlay%ngrow/mixlay%nano
     ngrow_ml(day_count) = ngrow_ml(day_count)+mixlay%ngrow
     dgrowD(1:M) = micro%growth%new(1:M)*P%micro%new(1:M)
     CALL average2(grid,dgrowD(1:M),h_m%new,mixlay%dgrow)
     mixlay%dgrow = mixlay%dgrow/mixlay%diatom
     dgrow_ml(day_count) = dgrow_ml(day_count)+mixlay%dgrow
     CALL average2(grid,zmicro%graze(1,1:M)*Z%micro%new(1:M),h_m%new,mixlay%zgrazen)    
     mixlay%zgrazen = mixlay%zgrazen/mixlay%nano !!! want graze*Zmicro/Pnano
     zgrazen_ml(day_count) = zgrazen_ml(day_count)+mixlay%zgrazen
     IF (mixlay%copepod > 0.) THEN
        CALL average2(grid,cgraze(1,1:M)+micro%mort%new(1:M),h_m%new,mixlay%cgrazed)
        mixlay%cgrazed = mixlay%cgrazed/mixlay%diatom !!! want (cgraze+background)/diatom
        cgrazed_ml(day_count) = cgrazed_ml(day_count)+mixlay%cgrazed
        CALL average2(grid,cgraze(2,1:M)+zmicro%mort%new(1:M),h_m%new,mixlay%cgrazez)
        mixlay%cgrazez = mixlay%cgrazez/mixlay%zmicro !!! want cgraze/zmicro
        cgrazez_ml(day_count) = cgrazez_ml(day_count) + mixlay%cgrazez
        CALL average2(grid,Detritus(2)%D%new(1:M),h_m%new,mixlay_D2)
        CALL average2(grid,cgraze(3,1:M),h_m%new,mixlay%cgrazed2)
        mixlay%cgrazed2 = mixlay%cgrazed2/mixlay_D2!!! want cgraze/D2
        cgrazed2_ml(day_count) = cgrazed2_ml(day_count) + mixlay%cgrazed2
     ELSE
        CALL average2(grid,micro%mort%new(1:M),h_m%new,mixlay%cgrazed)
        mixlay%cgrazed = mixlay%cgrazed/mixlay%diatom !!! want (cgraze+background)/diatom
        cgrazed_ml(day_count) = cgrazed_ml(day_count)+mixlay%cgrazed
        CALL average2(grid,zmicro%mort%new(1:M),h_m%new,mixlay%cgrazez)
        mixlay%cgrazez = mixlay%cgrazez/mixlay%zmicro !!! want cgraze/zmicro
        cgrazez_ml(day_count) = cgrazez_ml(day_count) + mixlay%cgrazez
     END IF
     CALL average2(grid, N%O_uptake%new(1:M),h_m%new,mixlay_NOup)
     NOup_ml(day_count) = NOup_ml(day_count) + mixlay_NOup
     CALL average2(grid, N%H_uptake%new(1:M),h_m%new,mixlay_NHup)
     NHup_ml(day_count) = NHup_ml(day_count) + mixlay_NHup
     D3rate(1:M) = (Detritus(3)%D%new(1:M)-Detritus(3)%D%old(1:M))/dt*3600.*24.  !daily rate of loss
     CALL sum_g(grid, D3rate(1:M),h_m%new,mixlay%out)
     D3_ml(day_count) = D3_ml(day_count) + mixlay%out  !average daily rate of loss frpm m.l.  gN m^-2 day^-1
    ! T2nano_ml(day_count) = T2nano_ml(day_count) + mixlay%ngrow
   !  TempT(1:M) = P%nano%new(1:M)*LOG(2.)/(nano%growth%new(1:M)+1.D-20)
   !  CALL average2(grid,TempT(1:M),h_m%new,doubleT%nano)
   !  doubleT%nano = doubleT%nano/mixlay%nano
   !  T2nano_ml(day_count) = T2nano_ml(day_count)+doubleT%nano
   !  TempT(1:M) = P%micro%new(1:M)*LOG(2.)/(micro%growth%new(1:M)+1.D-20)
   !  CALL average2(grid,TempT(1:M),h_m%new,doubleT%diatom)
   !  doubleT%diatom = doubleT%diatom/mixlay%diatom
   !  T2diatom_ml(day_count) = T2diatom_ml(day_count)+doubleT%diatom
     TempT(1:M) = Z%micro%new(1:M)*LOG(2.)/(zmicro%growth%new(1:M)+1.D-20)
     CALL average2(grid,TempT(1:M),h_m%new,doubleT%zmicro)
     doubleT%zmicro = doubleT%zmicro/mixlay%zmicro
     T2zmicro_ml(day_count) = T2zmicro_ml(day_count) + doubleT%zmicro
     TempT(1:M) = -P%nano%new(1:M)*LOG(0.5)/(nano%M_z+zmicro%graze(1,1:M)*Z%micro%new(1:M)/(P%nano%new(1:M)+1.D-20)+1.D-20)
     CALL average2(grid,TempT(1:M),h_m%new,halfT%nano)
     halfT%nano = halfT%nano/mixlay%nano
     Thalfnano_ml(day_count) = Thalfnano_ml(day_count)+halfT%nano
     TempT(1:M) = -P%micro%new(1:M)*LOG(0.5)/(nano%M_z+(micro%mort%new(1:M)+cgraze(1,1:M))/(P%micro%new(1:M)+1.D-20))
     CALL average2(grid,TempT(1:M),h_m%new,halfT%diatom)
     halfT%diatom = halfT%diatom/mixlay%diatom
     Thalfdiatom_ml(day_count) = Thalfdiatom_ml(day_count)+halfT%diatom
     TempT(1:M) = -Z%micro%new(1:M)*LOG(0.5)/(nano%M_z+(zmicro%mort%new(1:M)+cgraze(2,1:M))/(Z%micro%new(1:M)+1.D-20))
     CALL average2(grid,TempT(1:M),h_m%new,halfT%zmicro)
     halfT%zmicro = halfT%zmicro/mixlay%zmicro
     Thalfzmicro_ml(day_count) = Thalfzmicro_ml(day_count)+halfT%zmicro 
!!!ml fratio is <NOuptake>ml/<Totaluptake>ml and ezone fratio is <NOuptake/totaluptake>e
     Totalup(1:M) = N%H_uptake%new(1:M)+N%O_uptake%new(1:M)
     fratioavg(1:M) = f_ratio(1:M)*Totalup(1:M)
     CALL average2(grid, fratioavg(1:M),h_m%new,mixlay_fratio)
     CALL average2(grid, Totalup(1:M),h_m%new,mixlay_up)
     mixlay_fratio = mixlay_fratio/(mixlay_up+1.D-20)
     fratio_ml(day_count) = fratio_ml(day_count)+mixlay_fratio
!150 m average quantities 147.5 m
     CALL average2(grid,PONml(1:M),grid%d_g(g_150),depth150%PN)
     CALL average2(grid,DNml(1:M),grid%d_g(g_150),depth150%DN)
     CALL average2(grid,copepodml(1:M),grid%d_g(g_150),depth150%copepod)
     CALL average2(grid,N%H%new(1:M),grid%d_g(g_150),depth150%NH)
     CALL average2(grid,N%O%new(1:M),grid%d_g(g_150),depth150%NO)
     PN_150(day_count) = PN_150(day_count) +depth150%PN
     DN_150(day_count) = DN_150(day_count) +depth150%DN
     copepod_150(day_count) = copepod_150(day_count) +depth150%copepod
     NH_150(day_count) = NH_150(day_count) +depth150%NH
     NO_150(day_count) = NO_150(day_count) + depth150%NO
     
!100 m average quantities  actually 97.5 m
     CALL average2(grid,PONml(1:M),grid%d_g(g_100),depth100%PN)
     CALL average2(grid,DNml(1:M),grid%d_g(g_100),depth100%DN)
     PN_100(day_count) = PN_100(day_count) +depth100%PN
     DN_100(day_count) = DN_100(day_count) +depth100%DN

!80 m average quantities  77.5 m
     CALL average2(grid,NPPml(1:M),grid%d_g(g_80),depth80%NPP)
     NPP_80(day_count) = NPP_80(day_count)+depth80%NPP*grid%d_g(g_80)

!ezone average quantities
     CALL average2(grid,P%nano%new(1:M),grid%d_g(euph%g),ezone%nano)
     CALL average2(grid,P%micro%new(1:M),grid%d_g(euph%g),ezone%diatom)
     CALL average2(grid,PONml(1:M),grid%d_g(euph%g),ezone%PN)
     CALL average2(grid,DNml(1:M),grid%d_g(euph%g),ezone%DN)
     CALL average2(grid,N%O%new(1:M),grid%d_g(euph%g),ezone%NO)
     CALL average2(grid,N%H%new(1:M),grid%d_g(euph%g),ezone%NH)
     CALL average2(grid,NPPml(1:M),grid%d_g(euph%g),ezone%NPP)
     CALL average2(grid,ngrowP(1:M),grid%d_g(euph%g),ezone%ngrow)
     ezone%ngrow = ezone%ngrow/ezone%nano
     CALL average2(grid,dgrowD(1:M),grid%d_g(euph%g),ezone%dgrow)
     ezone%dgrow = ezone%dgrow/ezone%diatom
     CALL average2(grid,N%O_uptake%new(1:M),grid%d_g(euph%g),ezone_NOup)
     CALL average2(grid,N%H_uptake%new(1:M),grid%d_g(euph%g),ezone_NHup)
     CALL average2(grid,f_ratio(1:M),grid%d_g(euph%g),ezone_fratio)

     nano_e(day_count) = nano_e(day_count)+ezone%nano
     diatom_e(day_count) = diatom_e(day_count)+ezone%diatom
     PN_e(day_count) = PN_e(day_count)+ezone%PN
     DN_e(day_count) = DN_e(day_count)+ezone%DN
     NO_e(day_count) = NO_e(day_count)+ezone%NO
     NH_e(day_count) = NH_e(day_count)+ezone%NH
     NPP_e(day_count) = NPP_e(day_count)+ezone%NPP*grid%d_g(euph%g)
     ngrow_e(day_count) = ngrow_e(day_count)+ezone%ngrow
     dgrow_e(day_count) = dgrow_e(day_count)+ezone%dgrow
     NOup_e(day_count) = NOup_e(day_count) +ezone_NOup
     NHup_e(day_count) = NHup_e(day_count) +ezone_NHup
     fratio_e(day_count) = fratio_e(day_count) + ezone_fratio

!200 m total
     CALL sum_g(grid,D3rate(1:M),DBLE(grid%D),mixlay%out)
     CALL average2(grid,N%H%new(1:M),DBLE(grid%D),depth200%NH)
     CALL average2(grid,N%O%new(1:M),DBLE(grid%D),depth200%NO)
     CALL average2(grid,PONml(1:M),DBLE(grid%D),depth200%PN)
     D3_200(day_count) = D3_200(day_count) + mixlay%out
     NO_200(day_count) = NO_200(day_count) + depth200%NO
     NH_200(day_count) = NH_200(day_count) + depth200%NH 
     PN_200(day_count) = PN_200(day_count) +depth200%PN  
     !average daily rate of loss from the model gN m^-2 day^-1


  ELSE IF (time_step == steps) THEN  !time_step == steps
     DO xx = 1,sizerun  !27*years

      !  IF (cnt_wt(xx) <= 0) THEN
      !     molt_wt(xx) = 0.
      !  ELSE
      !     molt_wt(xx) = molt_wt(xx)/DBLE(cnt_wt(xx))
      !  END IF
      !  IF (cnt_avg_wt(xx) <= 0) THEN
      !     avg_wt(xx) = 0.
      !  ELSE
      !     avg_wt(xx) = avg_wt(xx)/DBLE(cnt_avg_wt(xx))
      !  END IF        
        IF (species_cnt(xx) <= 0) THEN          
           species_b(xx) =  0.
           species_mwt(xx) = 0.
           species_mwt_o(xx) = 0.
        ELSE
           species_b(xx) = species_b(xx)/DBLE(species_cnt(xx))
           species_mwt(xx) = species_mwt(xx)/DBLE(species_cnt(xx))
           species_mwt_o(xx) = species_mwt_o(xx)/DBLE(species_cnt(xx))
        END IF        
        IF (species_cnt2(xx) <= 0) THEN        
           stage1_no(xx) =  0.
           stage2_no(xx) =  0.
           stage3_no(xx) =  0.
           stage4_no(xx) =  0.
           stage5_no(xx) =  0.
           out_no(xx) =  0.
           species_avgwt(xx) =  0.
        ELSE      
           stage1_no(xx) = stage1_no(xx)/DBLE(species_cnt2(xx)) 
           stage2_no(xx) = stage2_no(xx)/DBLE(species_cnt2(xx)) 
           stage3_no(xx) = stage3_no(xx)/DBLE(species_cnt2(xx)) 
           stage4_no(xx) = stage4_no(xx)/DBLE(species_cnt2(xx)) 
           stage5_no(xx) = stage5_no(xx)/DBLE(species_cnt2(xx))  
           out_no(xx) = out_no(xx)/DBLE(species_cnt2(xx)) 
           species_avgwt(xx) = species_avgwt(xx)/DBLE(species_cnt2(xx)) 
        END IF
        IF (DBLE(o_cnt(xx)) <= 0.) THEN
           NO_o(xx) = 0.
           nano_ml(xx) = 0.
           diatom_ml(xx) = 0.
           zmicro_ml(xx) = 0.
           copepod_ml(xx) =0.
           NPP_ml(xx) = 0.
           PN_ml(xx) =0.
           DN_ml(xx) =0.
           NO_ml(xx) =0.
           NH_ml(xx) =0.
           SPN_ml(xx) = 0.
           T2nano_ml(xx) = 0.
           T2diatom_ml(xx) = 0.
           T2zmicro_ml(xx) = 0.
           Thalfnano_ml(xx) = 0.
           Thalfdiatom_ml(xx) = 0.
           Thalfzmicro_ml(xx) = 0.
           PONflux200(xx) = 0.
           PONflux100(xx) = 0.
           PONfluxml(xx) = 0.
           NOflux200(xx) = 0.
           NOfluxml(xx) = 0.
           NOflux100(xx) = 0.
           NHflux200(xx) = 0.
           NHfluxml(xx) = 0.
           NHflux100(xx) = 0.
           feacalml(xx) = 0.
           ureac(xx) = 0.
           ureaf(xx) = 0.
           migrateflux(xx) = 0.
           copepod_150(xx) = 0.
           NH_150(xx) = 0.
           NO_150(xx) = 0.
           PN_150(xx) = 0.
           DN_150(xx) = 0.
           PN_100(xx) = 0.
           DN_100(xx) = 0.
           nano_e(xx) = 0.
           diatom_e(xx) = 0.
           PN_e(xx) = 0.
           DN_e(xx) = 0.
           NO_e(xx) =0.
           NH_e(xx) =0.
           NOup_e(xx) =0.
           NHup_e(xx) =0.
           NOup_ml(xx) =0.
           NHup_ml(xx) =0.
           fratio_e(xx) = 0.
           fratio_ml(xx) = 0.
           NPP_e(xx) = 0.
           NPP_80(xx) = 0.
           ngrow_o(xx) = 0.
           dgrow_o(xx) = 0.
           ngrow_ml(xx) = 0.
           dgrow_ml(xx) = 0.
           zgrazen_ml(xx) = 0.
           cgrazed_ml(xx) = 0.
           cgrazez_ml(xx) = 0.
           cgrazed2_ml(xx) = 0.
           ngrow_e(xx) = 0.
           dgrow_e(xx) = 0.
           euph_new(xx) = 0.
           hm_new(xx) = 0.
           hm_g(xx) = 0
           euph_g(xx) = 0
           D3_ml(xx) = 0.
           D3_200(xx) = 0.
           NO_200(xx) = 0.
           NH_200(xx) = 0.
           PN_200(xx) = 0.  
          ! stage1_n(xx) = 0.
          ! stage2_n(xx) =  0. 
          ! stage3_n(xx) =  0.
          ! stage4_n(xx) =  0.
          ! stage5_n(xx) =  0.
          ! stage6_n(xx) =   0.
          ! nano_avg(xx) =  0.
          ! diatom_avg(xx) =  0. 
          ! zmicro_avg(xx) =  0.
          ! copepod_avg(xx) =  0.
          ! don_avg(xx) =  0.
          ! pon_avg(xx) =  0.
          ! out_avg(xx) = 0.
          ! NO_avg(xx) = 0. 
          ! NH_avg(xx) = 0.
          ! Ntot_avg(xx) =  0.
          ! urea_cop(xx) = 0.
          ! urea_fla(xx) = 0.
          ! wt_stage1(xx) = 0. 
           !wt_stage2(xx) = 0. 
          ! wt_stage3(xx) = 0.  
          ! wt_stage4(xx) = 0.
          ! wt_stage5(xx) = 0. 
          ! n_loss(xx) =  0.
          ! nano_gML(xx) =  0.
          ! micro_gML(xx) =  0.
          ! nano_g50(xx) =  0.
          ! micro_g50(xx) = 0. 
           !nano_go(xx)=  0.
           !micro_go(xx) = 0. 
          ! nano_zML(xx)=  0.
          ! micro_zML(xx) = 0. 
          ! nano_zo(xx)= 0.
          ! micro_zo(xx) = 0.
          ! NPPd(xx) =  0.
          ! NPPn(xx) = 0.
          ! NPPd_o(xx) = 0. 
          ! NPPn_o(xx) =  0.
          ! fratio(xx) = 0.
        ELSE
           NO_o(xx) = NO_o(xx)/DBLE(o_cnt(xx)) 
           nano_ml(xx) = nano_ml(xx)/DBLE(o_cnt(xx))
           diatom_ml(xx) = diatom_ml(xx)/DBLE(o_cnt(xx)) 
           zmicro_ml(xx) = zmicro_ml(xx)/DBLE(o_cnt(xx))
           copepod_ml(xx) = copepod_ml(xx)/DBLE(o_cnt(xx))
           NPP_ml(xx) = NPP_ml(xx)/DBLE(o_cnt(xx))
           PN_ml(xx) = PN_ml(xx)/DBLE(o_cnt(xx))
           SPN_ml(xx) = SPN_ml(xx)/DBLE(o_cnt(xx))
           DN_ml(xx) = DN_ml(xx)/DBLE(o_cnt(xx))
           NO_ml(xx) = NO_ml(xx)/DBLE(o_cnt(xx))
           NH_ml(xx) = NH_ml(xx)/DBLE(o_cnt(xx)) 
           T2nano_ml(xx) = T2nano_ml(xx)/DBLE(o_cnt(xx))
           T2diatom_ml(xx) = T2diatom_ml(xx)/DBLE(o_cnt(xx))
           T2zmicro_ml(xx) = T2zmicro_ml(xx)/DBLE(o_cnt(xx))
           Thalfnano_ml(xx) = Thalfnano_ml(xx)/DBLE(o_cnt(xx))
           Thalfdiatom_ml(xx) = Thalfdiatom_ml(xx)/DBLE(o_cnt(xx))
           Thalfzmicro_ml(xx) = Thalfzmicro_ml(xx)/DBLE(o_cnt(xx))
           PONflux200(xx) = PONflux200(xx)/DBLE(o_cnt(xx))
           PONflux100(xx) = PONflux100(xx)/DBLE(o_cnt(xx))
           PONfluxml(xx) = PONfluxml(xx)/DBLE(o_cnt(xx))
           NOflux200(xx) = NOflux200(xx)/DBLE(o_cnt(xx))
           NOfluxml(xx) = NOfluxml(xx)/DBLE(o_cnt(xx))
           NOflux100(xx) = NOflux100(xx)/DBLE(o_cnt(xx))
           NHflux200(xx) = NHflux200(xx)/DBLE(o_cnt(xx))
           NHfluxml(xx) = NHfluxml(xx)/DBLE(o_cnt(xx))
           NHflux100(xx) = NHflux100(xx)/DBLE(o_cnt(xx))
           feacalml(xx) = feacalml(xx)/DBLE(o_cnt(xx))
           ureac(xx) = ureac(xx)/DBLE(o_cnt(xx))
           ureaf(xx) = ureaf(xx)/DBLE(o_cnt(xx))
           migrateflux(xx) = migrateflux(xx)/DBLE(o_cnt(xx))
           copepod_150(xx) = copepod_150(xx)/DBLE(o_cnt(xx))
           PN_150(xx) = PN_150(xx)/DBLE(o_cnt(xx))
           PN_100(xx) = PN_100(xx)/DBLE(o_cnt(xx))
           DN_150(xx) = DN_150(xx)/DBLE(o_cnt(xx))
           DN_100(xx) = DN_100(xx)/DBLE(o_cnt(xx))
           PN_e(xx) = PN_e(xx)/DBLE(o_cnt(xx))
           DN_e(xx) = DN_e(xx)/DBLE(o_cnt(xx))
           NH_150(xx) = NH_150(xx)/DBLE(o_cnt(xx)) 
           NO_150(xx) = NO_150(xx)/DBLE(o_cnt(xx))  
           nano_e(xx) = nano_e(xx)/DBLE(o_cnt(xx)) 
           diatom_e(xx) = diatom_e(xx)/DBLE(o_cnt(xx)) 
           NO_e(xx) = NO_e(xx)/DBLE(o_cnt(xx)) 
           NH_e(xx) = NH_e(xx)/DBLE(o_cnt(xx))
           NOup_e(xx) =NOup_e(xx)/DBLE(o_cnt(xx)) 
           NHup_e(xx) = NHup_e(xx)/DBLE(o_cnt(xx)) 
           NOup_ml(xx) =NOup_ml(xx)/DBLE(o_cnt(xx))
           NHup_ml(xx) =NHup_ml(xx)/DBLE(o_cnt(xx))
           fratio_e(xx) = fratio_e(xx)/DBLE(o_cnt(xx)) 
           fratio_ml(xx) = fratio_ml(xx)/DBLE(o_cnt(xx))
           NPP_e(xx) = NPP_e(xx)/DBLE(o_cnt(xx)) 
           NPP_80(xx) = NPP_80(xx)/DBLE(o_cnt(xx)) 
           ngrow_o(xx) = ngrow_o(xx)/DBLE(o_cnt(xx))
           dgrow_o(xx) = dgrow_o(xx)/DBLE(o_cnt(xx))
           ngrow_ml(xx) = ngrow_ml(xx)/DBLE(o_cnt(xx))
           dgrow_ml(xx) = dgrow_ml(xx)/DBLE(o_cnt(xx))
           zgrazen_ml(xx) = zgrazen_ml(xx)/DBLE(o_cnt(xx))
           cgrazed_ml(xx) = cgrazed_ml(xx)/DBLE(o_cnt(xx))
           cgrazez_ml(xx) = cgrazez_ml(xx)/DBLE(o_cnt(xx))
           cgrazed2_ml(xx) = cgrazed2_ml(xx)/DBLE(o_cnt(xx))
           ngrow_e(xx) = ngrow_e(xx)/DBLE(o_cnt(xx))
           dgrow_e(xx) = dgrow_e(xx)/DBLE(o_cnt(xx)) 
           hm_new(xx) = hm_new(xx)/DBLE(o_cnt(xx))
           euph_new(xx) = euph_new(xx)/DBLE(o_cnt(xx))
           hm_g(xx) = INT(REAL(hm_g(xx))/REAL(o_cnt(xx)))
           euph_g(xx) = INT(REAL(euph_g(xx))/REAL(o_cnt(xx)))
           D3_ml(xx) = D3_ml(xx)/DBLE(o_cnt(xx))
           D3_200(xx) = D3_200(xx)/DBLE(o_cnt(xx))
           NO_200(xx) = NO_200(xx)/DBLE(o_cnt(xx))
           NH_200(xx) =NH_200(xx)/DBLE(o_cnt(xx))
           PN_200(xx) =PN_200(xx)/DBLE(o_cnt(xx))
        END IF
     END DO
     WRITE(125,100)REAL(NO_o(1:sizerun)),REAL(NO_ml(1:sizerun)),REAL(NO_e(1:sizerun)),REAL(NO_150(1:sizerun)),REAL(NO_200(1:sizerun))
     WRITE(110,100)REAL(nano_ml(1:sizerun)),REAL(nano_e(1:sizerun))
     WRITE(112,100)REAL(diatom_ml(1:sizerun)),REAL(diatom_e(1:sizerun))
     WRITE(114,100)REAL(zmicro_ml(1:sizerun))
     WRITE(116,100)REAL(copepod_ml(1:sizerun)),REAL(copepod_150(1:sizerun))
     WRITE(117,100)REAL(NPP_ml(1:sizerun)),REAL(NPP_80(1:sizerun)),REAL(NPP_e(1:sizerun))
     WRITE(122,100)REAL(PN_ml(1:sizerun)),REAL(PN_150(1:sizerun)),REAL(PN_e(1:sizerun)),REAL(PN_100(1:sizerun)),REAL(PN_200(1:sizerun))
     WRITE(122,100)REAL(SPN_ml(1:sizerun))  !suspended particulates
     WRITE(128,100)REAL(DN_ml(1:sizerun)),REAL(DN_150(1:sizerun)),REAL(DN_e(1:sizerun)),REAL(DN_100(1:sizerun))     
     WRITE(126,100)REAL(NH_ml(1:sizerun)),REAL(NH_e(1:sizerun)),REAL(NH_150(1:sizerun)),REAL(NH_200(1:sizerun))
     WRITE(119,100)REAL(NOup_e(1:sizerun)),REAL(NOup_ml(1:sizerun))
     WRITE(121,100)REAL(NHup_e(1:sizerun)),REAL(NHup_ml(1:sizerun))
     WRITE(123,100)REAL(fratio_ml(1:sizerun)),REAL(fratio_e(1:sizerun))
     WRITE(124,200)hm_g(1:sizerun),euph_g(1:sizerun)
     WRITE(124,100)REAL(hm_new(1:sizerun)),REAL(euph_new(1:sizerun))
     WRITE(127,100)REAL(ngrow_ml(1:sizerun)),REAL(ngrow_e(1:sizerun)),REAL(ngrow_o(1:sizerun))
     WRITE(127,100)REAL(dgrow_ml(1:sizerun)),REAL(dgrow_e(1:sizerun)),REAL(dgrow_o(1:sizerun))
     WRITE(129,100)REAL(species_b(1:sizerun)),REAL(species_mwt(1:sizerun)),REAL(species_mwt_o(1:sizerun))
     WRITE(129,100)REAL(species_avgwt(1:sizerun))
   !  WRITE(130,100)REAL(bin_logwt(1:bin_no))
  !   DO jj = 1,27
  !     WRITE(130,100) REAL(bin_cnt_year(jj,0:bin_no+2))
  !  END DO
     WRITE(131,100)REAL(stage1_no(1:sizerun)),REAL(stage2_no(1:sizerun)),REAL(stage3_no(1:sizerun))
     WRITE(131,100)REAL(stage4_no(1:sizerun)),REAL(stage5_no(1:sizerun)),REAL(out_no(1:sizerun))
     WRITE(133,100)REAL(T2nano_ml(1:sizerun)),REAL(T2diatom_ml(1:sizerun)),REAL(T2zmicro_ml(1:sizerun))
     WRITE(133,100)REAL(Thalfnano_ml(1:sizerun)),REAL(Thalfdiatom_ml(1:sizerun)),REAL(Thalfzmicro_ml(1:sizerun))
     WRITE(135,100)REAL(PONflux200(1:sizerun)),REAL(PONfluxml(1:sizerun)),REAL(PONflux100(1:sizerun)),REAL(feacalml(1:sizerun)),REAL(ureac(1:sizerun)),REAL(ureaf(1:sizerun)),REAL(D3_ml(1:sizerun)),REAL(D3_200(1:sizerun))
     WRITE(136,100)REAL(migrateflux(1:sizerun))
     WRITE(137,100)REAL(NOflux200(1:sizerun)),REAL(NOfluxml(1:sizerun)),REAL(NOflux100(1:sizerun)),REAL(NHflux200(1:sizerun)),REAL(NHfluxml(1:sizerun)),REAL(NHflux100(1:sizerun))
     WRITE(143,100)REAL(zgrazen_ml(1:sizerun)),REAL(cgrazed_ml(1:sizerun)),REAL(cgrazez_ml(1:sizerun)),REAL(cgrazed2_ml(1:sizerun))   !cgraze contains background term     
     o_cnt(1:sizerun) = 0.
     species_cnt(1:sizerun) = 0.
     species_cnt2(1:sizerun) = 0.
           NO_o(1:sizerun) = 0.
           nano_ml(1:sizerun) = 0.
           diatom_ml(1:sizerun) = 0.
           zmicro_ml(1:sizerun) = 0.
           copepod_ml(1:sizerun) =0.
           NPP_ml(1:sizerun) = 0.
           PN_ml(1:sizerun) =0.
           DN_ml(1:sizerun) =0.
           NO_ml(1:sizerun) =0.
           NH_ml(1:sizerun) =0.
           SPN_ml(1:sizerun) = 0.
           T2nano_ml(1:sizerun) = 0.
           T2diatom_ml(1:sizerun) = 0.
           T2zmicro_ml(1:sizerun) = 0.
           Thalfnano_ml(1:sizerun) = 0.
           Thalfdiatom_ml(1:sizerun) = 0.
           Thalfzmicro_ml(1:sizerun) = 0.
           PONflux200(1:sizerun) = 0.
           PONflux100(1:sizerun) = 0.
           PONfluxml(1:sizerun) = 0.
           NOflux200(1:sizerun) = 0.
           NOfluxml(1:sizerun) = 0.
           NOflux100(1:sizerun) = 0.
           NHflux200(1:sizerun) = 0.
           NHfluxml(1:sizerun) = 0.
           NHflux100(1:sizerun) = 0.
           feacalml(1:sizerun) = 0.
           ureac(1:sizerun) = 0.
           ureaf(1:sizerun) = 0.
           migrateflux(1:sizerun) = 0.
           copepod_150(1:sizerun) = 0.
           NH_150(1:sizerun) = 0.
           NO_150(1:sizerun) = 0.
           PN_150(1:sizerun) = 0.
           DN_150(1:sizerun) = 0.
           PN_100(1:sizerun) = 0.
           DN_100(1:sizerun) = 0.
           nano_e(1:sizerun) = 0.
           diatom_e(1:sizerun) = 0.
           PN_e(1:sizerun) = 0.
           DN_e(1:sizerun) = 0.
           NO_e(1:sizerun) =0.
           NH_e(1:sizerun) =0.
           NOup_e(1:sizerun) =0.
           NHup_e(1:sizerun) =0.
           NOup_ml(1:sizerun) =0.
           NHup_ml(1:sizerun) =0.
           fratio_e(1:sizerun) = 0.
           fratio_ml(1:sizerun) = 0.
           NPP_e(1:sizerun) = 0.
           NPP_80(1:sizerun) = 0.
           ngrow_o(1:sizerun) = 0.
           dgrow_o(1:sizerun) = 0.
           ngrow_ml(1:sizerun) = 0.
           dgrow_ml(1:sizerun) = 0.
           zgrazen_ml(1:sizerun) = 0.
           cgrazed_ml(1:sizerun) = 0.
           cgrazez_ml(1:sizerun) = 0.
           cgrazed2_ml(1:sizerun) = 0.
           ngrow_e(1:sizerun) = 0.
           dgrow_e(1:sizerun) = 0.
           euph_new(1:sizerun) = 0.
           hm_new(1:sizerun) = 0.
           hm_g(1:sizerun) = 0
           euph_g(1:sizerun) = 0
           D3_ml(1:sizerun) = 0.
           D3_200(1:sizerun) = 0.
           NO_200(1:sizerun) = 0.
           NH_200(1:sizerun) = 0.
           PN_200(1:sizerun) = 0.         
           species_b(1:sizerun) =  0.
           species_mwt(1:sizerun) = 0.
           species_mwt_o(1:sizerun) = 0.      
           stage1_no(1:sizerun) =  0.
           stage2_no(1:sizerun) =  0.
           stage3_no(1:sizerun) =  0.
           stage4_no(1:sizerun) =  0.
           stage5_no(1:sizerun) =  0.
           out_no(1:sizerun) =  0.
           species_avgwt(1:sizerun) =  0.    
          ! IF (ironday == 45) THEN
          !    ironday = 135
          ! ELSE IF (ironday == 135) THEN
          !    ironday = 225
          ! ELSE IF (ironday == 225) THEN
          !   ironday = 45
          ! END IF
          ! IF (ironday == 221) THEN
          !    ironday = 0
          ! END IF
!     PRINT "(A)","DONE"
!     PRINT "(A)","ironday"
!     PRINT *,ironday
   !  WRITE(147,100)REAL(molt_wt(1:sizerun))
   !  WRITE(149,100)REAL(avg_wt(1:sizerun))
   !  WRITE(151,100)REAL(urea_cop(1:sizerun))
   !  WRITE(153,100)REAL(urea_fla(1:sizerun))
   !  WRITE(155,100)REAL(PO_flux(1:sizerun))
   !  WRITE(157,100)REAL(NO_flux(1:sizerun))
   !  WRITE(159,100)REAL(fratio(1:sizerun))
   !  WRITE(129,100)REAL(stage1_n(1:sizerun)),REAL(stage2_n(1:sizerun)),REAL(stage3_n(1:sizerun)),REAL(stage4_n(1:sizerun)),REAL(stage5_n(1:sizerun)),REAL(stage6_n(1:sizerun))
     !WRITE(131,100)REAL(stage1_f),REAL(stage2_f),REAL(stage3_f),REAL(stage4_f),REAL(stage5_f)
    ! WRITE(133,100)REAL(Ntot_avg(1:sizerun))
    ! WRITE(135,100)REAL(nano_mar),REAL(diatom_mar),REAL(zmicro_mar),REAL(copepod_mar),REAL(NO_mar),REAL(NH_mar),REAL(don_mar),REAL(pon_mar)     
    ! WRITE(136,100)REAL(nano_jun),REAL(diatom_jun),REAL(zmicro_jun),REAL(copepod_jun),REAL(NO_jun),REAL(NH_jun),REAL(don_jun),REAL(pon_jun)
    ! WRITE(137,100)REAL(nano_sep),REAL(diatom_sep),REAL(zmicro_sep),REAL(copepod_sep),REAL(NO_sep),REAL(NH_sep),REAL(don_sep),REAL(pon_sep)
    ! WRITE(138,100)REAL(nano_dec),REAL(diatom_dec),REAL(zmicro_dec),REAL(copepod_dec),REAL(NO_dec),REAL(NH_dec),REAL(don_dec),REAL(pon_dec)
   !  WRITE(140,100)REAL(wt_stage1),REAL(wt_stage2),REAL(wt_stage3),REAL(wt_stage4),REAL(wt_stage5)
    ! WRITE(142,100)REAL(Ntot_avg(1:sizerun)), REAL(stage6_n(1:sizerun)),REAL( n_loss(1:sizerun))  !copeod_no.dat
    ! WRITE(144,100)REAL(nano_zML(1:sizerun)),REAL(nano_zo(1:sizerun)),REAL(micro_zML(1:sizerun)),REAL(micro_zo(1:sizerun))
    ! WRITE(145,100)REAL(NPPn(1:sizerun)),REAL(NPPn_o(1:sizerun)),REAL(NPPd(1:sizerun)),REAL(NPPnd_o(1:sizerun))
  END IF

!  IF (year == 1969 .AND. day /= ironday + 9 ) THEN
!     c_cnt(day) = c_cnt(day) + 1
!     c_biomass_150 = 0.
!     DO jj = 1,g_150 
!        IF (jj < g_150) THEN
!           c_biomass_150 = c_biomass_150 + species(1)%Z%new(jj)*species(1)%avg_wt*grid%i_space(jj) 
!        ELSE
!           c_biomass(day) = c_biomass(day) + (c_biomass_150+species(1)%Z%new(jj)*&
!                species(1)%avg_wt*grid%i_space(jj)/2.)/depth150%size + &
!                Zoo(1)%molt_wt(0)*Cevent(1)%nauplii/depth150%size*micro%M_z
!        END IF
!     END DO
!  ELSE IF (day == ironday + 9) THEN
 !    c_biomass(1:365) = c_biomass(1:365)/c_cnt(1:365)
!     WRITE(116,100)REAL(c_biomass(1:365))
!  END IF

END SUBROUTINE write_biological_sog






