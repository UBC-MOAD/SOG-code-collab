SUBROUTINE migrate_inout(Cevent,day,year,Csources,C_types,d,species,Zoo,euphotic_i,migrate,detritus_loss,&
     n_loss,timestep,check)

  USE mean_param
!  USE Copepod_mod

  IMPLICIT NONE


  INTEGER, INTENT(IN)::euphotic_i,check !interface level at about 1% Ipar (jmax_i)
  INTEGER, INTENT(IN)::Csources,C_types,day,year,timestep
  INTEGER, INTENT(IN OUT)::migrate
  TYPE(gr_d), INTENT(IN)::d !grid
  TYPE(event), DIMENSION(Csources),INTENT(IN OUT)::Cevent
  TYPE(copepod), DIMENSION(Csources),INTENT(IN OUT)::species
  TYPE(Cdata), DIMENSION(C_types),INTENT(IN)::Zoo
  DOUBLE PRECISION, DIMENSION(0:d%M+1),INTENT(IN OUT)::detritus_loss !Detritus(3)%D%new
  DOUBLE PRECISION, INTENT(IN OUT)::n_loss !total%n_loss
  INTEGER::ii

PRINT*,'migrate check',check
pause

  CALL Cevent_on(Cevent,day,year,Csources,species,Zoo,C_types,detritus_loss,d,n_loss,timestep)

  DO ii = 1,Csources
     IF (Cevent(ii)%on /= 0) THEN
        migrate = 1   !only use old values (not old_old) during timesteps when copepods migrate 
     END IF
     IF ((Cevent(ii)%on == 1 .OR. Cevent(ii)%on == 2) .AND. day /= Cevent(ii)%start .AND. timestep /= 1) THEN
        CALL bye_zoo(Cevent(ii),species(ii),day,Zoo(Cevent(ii)%type),detritus_loss,d,n_loss) !year
     END IF
     IF (Cevent(ii)%on == 1) THEN
        !PRINT "(A)","before new_zoo,species(1)%Ntot"
        !PRINT *,species(1)%Ntot
        CALL new_zoo(Cevent(ii),species(ii),day,euphotic_i,d,Zoo(Cevent(ii)%type),timestep)
     END IF
     species(ii)%Z%old = species(ii)%Z%new 
     species(ii)%mature_pdf%o_wt_avg = species(ii)%mature_pdf%wt_avg
  END DO


END SUBROUTINE migrate_inout

