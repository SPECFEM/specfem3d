!=====================================================================
!
!          S p e c f e m 3 D  B a s i n  V e r s i o n  1 . 2
!          --------------------------------------------------
!
!                 Dimitri Komatitsch and Jeroen Tromp
!    Seismological Laboratory - California Institute of Technology
!         (c) California Institute of Technology July 2004
!
!    A signed non-commercial agreement is required to use this program.
!   Please check http://www.gps.caltech.edu/research/jtromp for details.
!           Free for non-commercial academic research ONLY.
!      This program is distributed WITHOUT ANY WARRANTY whatsoever.
!      Do not redistribute this program without written permission.
!
!=====================================================================

!=====================================================================
! 07/09/04 Last changed by Min Chen
! Users need to modify this subroutine to implement their own
! anisotropic models.
!=====================================================================

  subroutine aniso_model(idoubling,zmesh,rho,vp,vs,c11,c12,c13,c14,c15,c16, &
               c22,c23,c24,c25,c26,c33,c34,c35,c36,c44,c45,c46,c55,c56,c66)

  implicit none

  include "constants.h"

!------------------------------------------------------------------------------
! for anisotropy simulations in a halfspace model

! only related to body waves
! one-zeta term
  double precision, parameter :: FACTOR_CS1p_A = 0.01d0
  double precision, parameter :: FACTOR_CS1sv_A = 0.0d0
  double precision, parameter :: FACTOR_CS1sh_N = 0.d0
! three-zeta term
  double precision, parameter :: FACTOR_CS3_L = 0.0d0

! Relative to Love wave
! four-zeta term
  double precision, parameter :: FACTOR_N = 0.d0
  double precision, parameter :: FACTOR_E_N = 0.d0

! Relative to Rayleigh wave
! two-zeta term
  double precision, parameter :: FACTOR_A = 0.d0
  double precision, parameter :: FACTOR_C = 0.d0
  double precision, parameter :: FACTOR_F = 0.d0
  double precision, parameter :: FACTOR_H_F = 0.d0
  double precision, parameter :: FACTOR_B_A = 0.d0

! Relative to both Love wave and Rayleigh wave
! two-zeta term
  double precision, parameter :: FACTOR_L = 0.d0
  double precision, parameter :: FACTOR_G_L = 0.d0

!------------------------------------------------------------------------------

  integer idoubling

  double precision zmesh
  double precision rho,vp,vs
  double precision vpv,vph,vsv,vsh,eta_aniso
  double precision aa,cc,nn,ll,ff
  double precision A,C,F,AL,AN,Bc,Bs,Gc,Gs,Hc,Hs,Ec,Es,C1p,C1sv,C1sh,C3,S1p,S1sv,S1sh,S3
  double precision c11,c12,c13,c14,c15,c16,c22,c23,c24,c25,c26,c33,c34,c35,c36, &
                   c44,c45,c46,c55,c56,c66
  double precision d11,d12,d13,d14,d15,d16,d22,d23,d24,d25,d26,d33,d34,d35,d36, &
                   d44,d45,d46,d55,d56,d66

! implement the background model
  if(idoubling == IFLAG_HALFSPACE_MOHO) then
    vp=7.8d0
    vs=4.5d0
    rho=3.0d0
    vph = vp
    vpv = vp
    vsh = vs
    vsv = vs
    eta_aniso = 1.0d0

  else if(idoubling == IFLAG_MOHO_16km) then
    vp=7.8d0
    vs=4.5d0
    rho=3.0d0
    vph = vp
    vpv = vp
    vsh = vs
    vsv = vs
    eta_aniso = 1.0d0

  else if(zmesh >= DEPTH_5p5km_SOCAL) then
    vp=7.8d0
    vs=4.5d0
    rho=3.0d0
    vph = vp
    vpv = vp
    vsh = vs
    vsv = vs
    eta_aniso = 1.0d0

  else
    vp=7.8d0
    vs=4.5d0
    rho=3.0d0
    vph = vp
    vpv = vp
    vsh = vs
    vsv = vs
    eta_aniso = 1.0d0

  endif

! scale to standard units
  vp = vp * 1000.d0
  vs = vs * 1000.d0
  vph = vph * 1000.d0
  vpv = vpv * 1000.d0
  vsh = vsh * 1000.d0
  vsv = vsv * 1000.d0
  rho = rho * 1000.d0

  aa = rho*vph*vph
  cc = rho*vpv*vpv
  nn = rho*vsh*vsh
  ll = rho*vsv*vsv
  ff = eta_aniso*(aa - 2.*ll)

! Add anisotropic perturbation in the whole halfspace
! You can also add different perturbations to different layers
  A = aa*(1.0d0 + FACTOR_A)
  C = cc*(1.0d0 + FACTOR_C)
  AN = nn*(1.0d0 + FACTOR_N)
  AL = ll*(1.0d0 + FACTOR_L)
  F = ff*(1.0d0 + FACTOR_F)
  C1p = FACTOR_CS1p_A*aa
  S1p = 0.d0
  C1sv = FACTOR_CS1sv_A*aa
  S1sv = 0.d0
  C1sh = FACTOR_CS1sh_N*nn
  S1sh = 0.d0
  Gc = FACTOR_G_L*ll
  Gs = 0.d0
  Bc = FACTOR_B_A*aa
  Bs = 0.d0
  Hc = FACTOR_H_F*ff
  Hs = 0.d0
  C3 = FACTOR_CS3_L*ll
  S3 = 0.d0
  Ec = FACTOR_E_N*nn
  Es = 0.d0

! The mapping from the elastic coefficients to the elastic tensor elements
! in the local Cartesian coordinate system (classical geographic) used in the
! global code (1---South, 2---East, 3---up)
! Always keep the following part when you modify this subroutine
  d11 = A + Ec + Bc
  d12 = A - 2.*AN - Ec
  d13 = F + Hc
  d14 = S3 + 2.*S1sh + 2.*S1p
  d15 = 2.*C1p + C3
  d16 = -Bs/2. - Es
  d22 = A + Ec - Bc
  d23 = F - Hc
  d24 = 2.*S1p - S3
  d25 = 2.*C1p - 2.*C1sh - C3
  d26 = -Bs/2. + Es
  d33 = C
  d34 = 2.*(S1p - S1sv)
  d35 = 2.*(C1p - C1sv)
  d36 = -Hs
  d44 = AL - Gc
  d45 = -Gs
  d46 = C1sh - C3
  d55 = AL + Gc
  d56 = S3 - S1sh
  d66 = AN - Ec

! The mapping to the global Cartesian coordinate system used in the basin code
! (1---East, 2---North, 3---up)
  c11 = d22
  c12 = d12
  c13 = d23
  c14 = - d25
  c15 = d24
  c16 = - d26
  c22 = d11
  c23 = d13
  c24 = - d15
  c25 = d14
  c26 = - d16
  c33 = d33
  c34 = - d35
  c35 = d34
  c36 = - d36
  c44 = d55
  c45 = - d45
  c46 = d56
  c55 = d44
  c56 = - d46
  c66 = d66

  end subroutine aniso_model

