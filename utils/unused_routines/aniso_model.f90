!=====================================================================
!
!               S p e c f e m 3 D  V e r s i o n  2 . 1
!               ---------------------------------------
!
!          Main authors: Dimitri Komatitsch and Jeroen Tromp
!    Princeton University, USA and University of Pau / CNRS / INRIA
! (c) Princeton University / California Institute of Technology and University of Pau / CNRS / INRIA
!                            November 2010
!
! This program is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 2 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License along
! with this program; if not, write to the Free Software Foundation, Inc.,
! 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
!
!=====================================================================

!=====================================================================
! 07/09/04 Last changed by Min Chen
! Users need to modify this subroutine to implement their own
! anisotropic models.
!=====================================================================

  subroutine aniso_model(iflag_aniso,rho,vp,vs,c11,c12,c13,c14,c15,c16, &
               c22,c23,c24,c25,c26,c33,c34,c35,c36,c44,c45,c46,c55,c56,c66)

  implicit none

  include "constants.h"

! see for example:
!
! M. Chen & J. Tromp, 2006. Theoretical & numerical investigations
! of global and regional seismic wave propagation in weakly anisotropic earth models,
! GJI, 168, 1130-1152.

!------------------------------------------------------------------------------
! for anisotropy simulations in a halfspace model

! only related to body waves
! one-zeta term
  real(kind=CUSTOM_REAL), parameter :: FACTOR_CS1p_A = 0.2_CUSTOM_REAL
  real(kind=CUSTOM_REAL), parameter :: FACTOR_CS1sv_A = 0._CUSTOM_REAL
  real(kind=CUSTOM_REAL), parameter :: FACTOR_CS1sh_N = 0._CUSTOM_REAL
! three-zeta term
  real(kind=CUSTOM_REAL), parameter :: FACTOR_CS3_L = 0._CUSTOM_REAL

! Relative to Love wave
! four-zeta term
  real(kind=CUSTOM_REAL), parameter :: FACTOR_N = 0._CUSTOM_REAL
  real(kind=CUSTOM_REAL), parameter :: FACTOR_E_N = 0._CUSTOM_REAL

! Relative to Rayleigh wave
! two-zeta term
  real(kind=CUSTOM_REAL), parameter :: FACTOR_A = 0._CUSTOM_REAL
  real(kind=CUSTOM_REAL), parameter :: FACTOR_C = 0._CUSTOM_REAL
  real(kind=CUSTOM_REAL), parameter :: FACTOR_F = 0._CUSTOM_REAL
  real(kind=CUSTOM_REAL), parameter :: FACTOR_H_F = 0._CUSTOM_REAL
  real(kind=CUSTOM_REAL), parameter :: FACTOR_B_A = 0._CUSTOM_REAL

! Relative to both Love wave and Rayleigh wave
! two-zeta term
  real(kind=CUSTOM_REAL), parameter :: FACTOR_L = 0._CUSTOM_REAL
  real(kind=CUSTOM_REAL), parameter :: FACTOR_G_L = 0._CUSTOM_REAL

!------------------------------------------------------------------------------

  !integer idoubling
  integer iflag_aniso

  !real(kind=CUSTOM_REAL) zmesh
  real(kind=CUSTOM_REAL) rho,vp,vs
  real(kind=CUSTOM_REAL) c11,c12,c13,c14,c15,c16,c22,c23,c24,c25,c26,c33,c34,c35,c36, &
                   c44,c45,c46,c55,c56,c66

! local parameters
  real(kind=CUSTOM_REAL) vpv,vph,vsv,vsh,eta_aniso
  real(kind=CUSTOM_REAL) aa,cc,nn,ll,ff
  real(kind=CUSTOM_REAL) A,C,F,AL,AN,Bc,Bs,Gc,Gs,Hc,Hs,Ec,Es,C1p,C1sv,C1sh,C3,S1p,S1sv,S1sh,S3
  real(kind=CUSTOM_REAL) d11,d12,d13,d14,d15,d16,d22,d23,d24,d25,d26,d33,d34,d35,d36, &
                   d44,d45,d46,d55,d56,d66

! assumes vp,vs given in m/s, rho in kg/m**3
  vph = vp
  vpv = vp
  vsh = vs
  vsv = vs
  eta_aniso = 1.0_CUSTOM_REAL


! for definition, see for example:
!
! Dziewonski & Anderson, 1981. Preliminary reference earth model, PEPI, 25, 297-356.
! page 305:
  aa = rho*vph*vph
  cc = rho*vpv*vpv
  nn = rho*vsh*vsh
  ll = rho*vsv*vsv
  ff = eta_aniso*(aa - 2.*ll)

! Add anisotropic perturbation

! notation: see Chen & Tromp, 2006, appendix A, page 1151
!
! zeta-independant terms:
! A = \delta A
! C = \delta C
! AN = \delta N
! AL = \delta L
! F = \delta F
!
! zeta-dependant terms:
! C1p =  J_c
! C1sv = K_c
! C1sh = M_c
! S1p =  J_s
! S1sv = K_s
! S1sh = M_s
!
! two-zeta dependant terms:
! Gc = G_c
! Gs = G_s
! Bc = B_c
! Bs = B_s
! Hc = H_c
! Hs =  H_s
!
! three-zeta dependant terms:
! C3 = D_c
! S3 = D_s
!
! four-zeta dependant terms:
! Ec = E_c
! Es = E_s

! no anisotropic perturbation
  if( iflag_aniso <= 0 ) then
    ! zeta-independant
    A = aa
    C = cc
    AN = nn
    AL = ll
    F = ff

    ! zeta-dependant terms
    C1p = 0._CUSTOM_REAL
    C1sv = 0._CUSTOM_REAL
    C1sh = 0._CUSTOM_REAL
    S1p = 0._CUSTOM_REAL
    S1sv = 0._CUSTOM_REAL
    S1sh = 0._CUSTOM_REAL

    ! two-zeta dependant terms
    Gc = 0._CUSTOM_REAL
    Gs = 0._CUSTOM_REAL

    Bc = 0._CUSTOM_REAL
    Bs = 0._CUSTOM_REAL

    Hc = 0._CUSTOM_REAL
    Hs = 0._CUSTOM_REAL

    ! three-zeta dependant terms
    C3 = 0._CUSTOM_REAL
    S3 = 0._CUSTOM_REAL

    ! four-zeta dependant terms
    Ec = 0._CUSTOM_REAL
    Es = 0._CUSTOM_REAL
  endif

! perturbation model 1
  if( iflag_aniso == IANISOTROPY_MODEL1 ) then
    ! zeta-independant
    A = aa*(1.0_CUSTOM_REAL + FACTOR_A)
    C = cc*(1.0_CUSTOM_REAL + FACTOR_C)
    AN = nn*(1.0_CUSTOM_REAL + FACTOR_N)
    AL = ll*(1.0_CUSTOM_REAL + FACTOR_L)
    F = ff*(1.0_CUSTOM_REAL + FACTOR_F)

    ! zeta-dependant terms
    C1p = FACTOR_CS1p_A*aa
    C1sv = FACTOR_CS1sv_A*aa
    C1sh = FACTOR_CS1sh_N*nn
    S1p = 0._CUSTOM_REAL
    S1sv = 0._CUSTOM_REAL
    S1sh = 0._CUSTOM_REAL

    ! two-zeta dependant terms
    Gc = FACTOR_G_L*ll
    Bc = FACTOR_B_A*aa
    Hc = FACTOR_H_F*ff
    Gs = 0._CUSTOM_REAL
    Bs = 0._CUSTOM_REAL
    Hs = 0._CUSTOM_REAL

    ! three-zeta dependant terms
    C3 = FACTOR_CS3_L*ll
    S3 = 0._CUSTOM_REAL

    ! four-zeta dependant terms
    Ec = FACTOR_E_N*nn
    Es = 0._CUSTOM_REAL
  endif

! perturbation model 2
  if( iflag_aniso == IANISOTROPY_MODEL2 ) then
    ! zeta-independant
    A = aa*(1.0_CUSTOM_REAL + FACTOR_A + 0.1)
    C = cc*(1.0_CUSTOM_REAL + FACTOR_C + 0.1)
    AN = nn*(1.0_CUSTOM_REAL + FACTOR_N + 0.1)
    AL = ll*(1.0_CUSTOM_REAL + FACTOR_L + 0.1)
    F = ff*(1.0_CUSTOM_REAL + FACTOR_F + 0.1)

    ! zeta-dependant terms
    C1p = FACTOR_CS1p_A*aa
    C1sv = FACTOR_CS1sv_A*aa
    C1sh = FACTOR_CS1sh_N*nn
    S1p = 0._CUSTOM_REAL
    S1sv = 0._CUSTOM_REAL
    S1sh = 0._CUSTOM_REAL

    ! two-zeta dependant terms
    Gc = FACTOR_G_L*ll
    Bc = FACTOR_B_A*aa
    Hc = FACTOR_H_F*ff
    Gs = 0._CUSTOM_REAL
    Bs = 0._CUSTOM_REAL
    Hs = 0._CUSTOM_REAL

    ! three-zeta dependant terms
    C3 = FACTOR_CS3_L*ll
    S3 = 0._CUSTOM_REAL

    ! four-zeta dependant terms
    Ec = FACTOR_E_N*nn
    Es = 0._CUSTOM_REAL
  endif


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

! The mapping to the global Cartesian coordinate system used in the code
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

