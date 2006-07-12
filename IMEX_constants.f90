MODULE IMEX_constants

	IMPLICIT NONE

!!!See Ascher et al (1995), SIAM J. Numer Anal. 32(3), pp. 797-823!!!!

!!!!implicit-explicit first order scheme:  a_IMEX1 = 1.0 --> backward Euler in diffusion
!!!				                   = 0.5 --> Crank-Nicolson in diffusion
!!!                                                = 0.0 --> forward Euler in diffusion
!!!                                       Both have forward Euler in explicit  terms (i.e. Coriolis and 
!!!                                                                                 Boundary conditions)

	DOUBLE PRECISION, PARAMETER::a_IMEX1 = 1.0

!!!!implicit-explicit second order scheme: (a_IMEX2, b_IMEX2) = (0.5, 0) --> CNAB (Crank-Nicolson, Adams-
!!!!                                                                               Bashforth)
!!			                                      = (0.5, 1/8) --> MCNAB (modified Crank-Nicolson,
!!!!                                                             1/8 = 0.125              Adams-Bashforth)
!!!							      = (0, 1)--> CNLF (Crank-Nicolson, Leap frog)
!!!!                                                          = (1, 0)--> SBDF (second order backward 
!!!                                                                             differencing)
	DOUBLE PRECISION, PARAMETER::a_IMEX2 = 1.0, b_IMEX2 =0.0

END MODULE IMEX_constants









