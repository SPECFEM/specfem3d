!=====================================================================
!
!          S p e c f e m 3 D  B a s i n  V e r s i o n  1 . 3
!          --------------------------------------------------
!
!                 Dimitri Komatitsch and Jeroen Tromp
!    Seismological Laboratory - California Institute of Technology
!         (c) California Institute of Technology July 2005
!
!    A signed non-commercial agreement is required to use this program.
!   Please check http://www.gps.caltech.edu/research/jtromp for details.
!           Free for non-commercial academic research ONLY.
!      This program is distributed WITHOUT ANY WARRANTY whatsoever.
!      Do not redistribute this program without written permission.
!
!=====================================================================

  double precision function comp_source_time_function(t,hdur)

  implicit none

  include "constants.h"

  double precision t,hdur

  double precision, external :: erf

!!!!!!!! DK DK XXXXXXXXXXX YYYYYYYYYY UGLY for Carcione copper aniso
#ifdef CARCIONE_ANISO
  include "carcione_anisotropy.h"
  double precision a_source,t0_source
#endif

! quasi Heaviside, small Gaussian moment-rate tensor with hdur
  comp_source_time_function = 0.5d0*(1.0d0 + erf(t/hdur))

#ifdef CARCIONE_ANISO
! onset time is chosen as characteristic time plus 25 %
  a_source = (PI*SOURCE_DOMINANT_FREQ)**2
  t0_source = 1.25d0 / SOURCE_DOMINANT_FREQ

! Gaussian
  if(SOURCE_TIME_FUNCTION == 1) then
    comp_source_time_function = + FACTOR_SOURCE*exp(-a_source*(t-t0_source)**2)

! Ricker
  else if(SOURCE_TIME_FUNCTION == 2) then
    comp_source_time_function = - FACTOR_SOURCE*(1.d0-2.d0*a_source*(t-t0_source)**2)*exp(-a_source*(t-t0_source)**2)

  else
    stop 'incorrect type of source time function for Carcione copper'
  endif
#endif

  end function comp_source_time_function

