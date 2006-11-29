! $Id$
! $Source$

module mean_param
  ! This should probably be renamed to something more descriptive like
  ! type_defs.
  ! It would also be good to re-organize somehow so that they types
  ! hierarchy is more evident.

  use precision_defs, only: dp

  implicit none

  ! Phytoplankton component of UVST type
  type :: phyto                     
     real(kind=dp), dimension(:), pointer :: micro, nano
  end type phyto

  ! Type for Gvectors and Hvector used in implicit solver
  type :: UVST              
     real(kind=dp), dimension(:), pointer :: u, v, s, t
  end type UVST

  ! Grow component of plankton2 type
  type :: grow               
     real(kind=dp), dimension(:), pointer :: light, new
  end type grow

  ! *** Plankton "behaviours" ??? type
  ! *** This will probably be refactored when Susan implements sinking
  type :: plankton2
     type(grow) :: growth
     real(kind=dp) :: sink_min, sink_max  ! Sinking velocity
     real(kind=dp), dimension(:), pointer:: Nlimit
  end type plankton2

  TYPE :: entrain           !Entrainment or mixed layer depth
     DOUBLE PRECISION :: depth
     INTEGER :: i,g
  END TYPE entrain

  TYPE :: height            !boundary layer depth 
     DOUBLE PRECISION :: old, new
     INTEGER :: i, g        !interface and grid index
     TYPE(entrain)::e, ml   !entrainment depth or mixed layer depth
  END TYPE height

  TYPE :: MST               ! Momentum, Salinity and Temp component vectors
     DOUBLE PRECISION, DIMENSION(:), POINTER::m,s,t    
  END TYPE MST

  TYPE :: boundary          !boundary conditions at h
     DOUBLE PRECISION :: h, div
     DOUBLE PRECISION, DIMENSION(:), POINTER::value
  END TYPE boundary

  TYPE :: MS                !Momemtum and scalar
     TYPE(boundary)::m, s
  END TYPE MS

  INTEGER::is_leap_year, was_leap_year

END MODULE mean_param
