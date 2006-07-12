F90 = pgf90
FC = pgf77
FFLAGS = -O3 
SRCS = mean_param.f90 declarations.f90 surface_forcing.f90 allocate1.f90 allocate2.f90 Julian_day.f90 y_jday_t.f90 read_sog.f90 find_wind.f90 smoothdata.f90 coefficients.f90 initialize.f90 define_grid.f90 define_sog.f90 initial_sog.f90 alpha_sub.f90 Cp_sub.f90 density_sub.f90 irradiance_sog.f90 define_adv.f90 buoyancy.f90 interpolate.f90 surface_flux_sog.f90 fun_constants.f90 stability.f90 ND_flux_profile.f90 vel_scales.f90 convection_scales.f90 shear_diff.f90 double_diff.f90 interior_match.f90 interior_match2.f90 shaper_parameters.f90 modify_K.f90 def_gamma.f90 define_flux.f90 diffusion.f90 horizontal_adv.f90 define_adv_bio.f90 Coriolis_and_pg.f90 pdf.f90 IMEX_constants.f90 matrix_A.f90 scalar_H.f90 U_H.f90 TRIDAG.f def_v_t_sog.f90 define_Ri_b_sog.f90 ML_height_sog.f90 define_hm_sog.f90 write_open.f90 write_physical_sog.f90 new_year.f90 find_upwell.f90 SOG.f90 define_PZ.f90 odeint.f90 reaction_p_sog.f90 P_H.f90 N_H.f90 find_new.f90 write_biological_sog.f90 allocate4.f90 fit.f polint.f gammq.f gser.f gcf.f gammln.f allocate3.f90 reaction.f90 derivs_sog.f90 rkqs.f rkck.f find_PON.f90 advection.f90

OBJS = mean_param.o declarations.o surface_forcing.o allocate1.o allocate2.o Julian_day.o y_jday_t.o read_sog.o find_wind.o smoothdata.o coefficients.o initialize.o define_grid.o define_sog.o initial_sog.o alpha_sub.o Cp_sub.o density_sub.o irradiance_sog.o define_adv.o buoyancy.o interpolate.o surface_flux_sog.o fun_constants.o stability.o ND_flux_profile.o vel_scales.o convection_scales.o shear_diff.o double_diff.o interior_match.o interior_match2.o shape_parameters.o modify_K.o def_gamma.o define_flux.o diffusion.o horizontal_adv.o define_adv_bio.o Coriolis_and_pg.o  pdf.o IMEX_constants.o matrix_A.o scalar_H.o U_H.o def_v_t_sog.o TRIDAG.o define_Ri_b_sog.o ML_height_sog.o define_hm_sog.o write_open.o write_physical_sog.o new_year.o find_upwell.o SOG.o define_PZ.o odeint.o reaction_p_sog.o P_H.o N_H.o find_new.o write_biological_sog.o allocate4.o fit.o polint.o gammq.o gser.o gcf.o gammln.o allocate3.o reaction.o derivs_sog.o rkqs.o rkck.o find_PON.o advection.o

SOG: Makefile $(OBJS)
	 $(F90) $(FFLAGS) -o SOG $(OBJS)
.SUFFIXES: 		# Delete built-in implicit suffix-based rules
.SUFFIXES: .f .c .f90 .o .mod

.f90.o:
	$(F90) -c $(FFLAGS) $<
.f.o: 
	$(F90) -c $(FFLAGS) $<

SOG.o: read_sog.o surface_forcing.o declarations.o mean_param.o initial_sog.o initialize.o define_adv.o coefficients.o surface_flux_sog.o shear_diff.o find_wind.o allocate1.o

shear_diff.o: surface_forcing.o 

clean:
	rm *.o *.mod
