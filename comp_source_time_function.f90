!=====================================================================
!
!          S p e c f e m 3 D  B a s i n  V e r s i o n  1 . 1
!          --------------------------------------------------
!
!                 Dimitri Komatitsch and Jeroen Tromp
!    Seismological Laboratory - California Institute of Technology
!         (c) California Institute of Technology October 2002
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

! source decay rate (also change in source spectrum if needed)
  double precision, parameter :: decay_rate = 2.628d0

  double precision t,hdur

  double precision, external :: erf

! Gaussian moment-rate tensor
  comp_source_time_function = 0.5d0*(1.0d0+erf(decay_rate*(t-hdur)/hdur))

  end function comp_source_time_function
