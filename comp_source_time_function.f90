!=====================================================================
!
!          S p e c f e m 3 D  G l o b e  V e r s i o n  3 . 5
!          --------------------------------------------------
!
!                 Dimitri Komatitsch and Jeroen Tromp
!    Seismological Laboratory - California Institute of Technology
!        (c) California Institute of Technology July 2004
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

  double precision hdur_gauss
  double precision, external :: erf

! Gaussian moment-rate tensor
! for Gaussian use 1.66667*hdur to get roughly a triangle with half-duration hdur
  hdur_gauss = hdur * 5. / 3.
  comp_source_time_function = 0.5d0*(1.0d0+erf(SOURCE_DECAY_RATE*t/hdur_gauss))

  end function comp_source_time_function

