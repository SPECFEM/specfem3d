!=====================================================================
!
!          S p e c f e m 3 D  G l o b e  V e r s i o n  4 . 0
!          --------------------------------------------------
!
!          Main authors: Dimitri Komatitsch and Jeroen Tromp
!    Seismological Laboratory, California Institute of Technology, USA
!             and University of Pau / CNRS / INRIA, France
! (c) California Institute of Technology and University of Pau / CNRS / INRIA
!                            February 2008
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

  subroutine read_aniso_inner_core_model

  implicit none

! one should add an MPI_BCAST in meshfem3D.f90 if one adds a read_aniso_inner_core_model subroutine

  end subroutine read_aniso_inner_core_model

!-----------------------------------

  subroutine aniso_inner_core_model(x,c11,c33,c12,c13,c44,REFERENCE_1D_MODEL)

  implicit none

  include "constants.h"

! given a normalized radius x, gives non-dimensionalized c11,c33,c12,c13,c44

  integer REFERENCE_1D_MODEL

  double precision x,c11,c33,c12,c13,c44

  double precision vp,vs,rho
  double precision vp0,vs0,rho0,A0
  double precision c66
  double precision scale_fac

  if(REFERENCE_1D_MODEL == REFERENCE_MODEL_IASP91) then
    vp=11.24094d0-4.09689d0*x*x
    vs=3.56454d0-3.45241d0*x*x
    rho=13.0885d0-8.8381d0*x*x

! values at center
    vp0=11.24094d0
    vs0=3.56454d0
    rho0=13.0885d0

  else if(REFERENCE_1D_MODEL == REFERENCE_MODEL_PREM) then
    vp=11.2622d0-6.3640d0*x*x
    vs=3.6678d0-4.4475d0*x*x
    rho=13.0885d0-8.8381d0*x*x

! values at center
    vp0=11.2622d0
    vs0=3.6678d0
    rho0=13.0885d0

  else
    stop 'unknown 1D reference Earth model in anisotropic inner core'
  endif

! elastic tensor for hexagonal symmetry in reduced notation:
!
!      c11 c12 c13  0   0        0
!      c12 c11 c13  0   0        0
!      c13 c13 c33  0   0        0
!       0   0   0  c44  0        0
!       0   0   0   0  c44       0
!       0   0   0   0   0  c66=(c11-c12)/2
!
!       in terms of the A, C, L, N and F of Love (1927):
!
!       c11 = A
!       c33 = C
!       c12 = A-2N
!       c13 = F
!       c44 = L
!       c66 = N
!
!       isotropic equivalent:
!
!       c11 = lambda+2mu
!       c33 = lambda+2mu
!       c12 = lambda
!       c13 = lambda
!       c44 = mu
!       c66 = mu

! non-dimensionalization of elastic parameters
  scale_fac=RHOAV*R_EARTH*R_EARTH*PI*GRAV*RHOAV

! Ishii et al. (2002):
!
! alpha = 3.490 % = (C-A)/A0    = (c33-c11)/A0
!  beta = 0.988 % = (L-N)/A0    = (c44-c66)/A0
! gamma = 0.881 % = (A-2N-F)/A0    = (c12-c13)/A0
! where A0 is A at the Earth's center
!
! assume c11 = lamda+2mu
!        c66 = (c11-c12)/2 = mu
!
! then   c33 = c11 + alpha*A0
!        c44 = c66 + beta*A0
!        c13 = c12 - gamma*A0
!
! Steinle-Neumann (2002):
!
!  r    T    rho    c11   c12  c13  c33  c44 KS   mu
! (km) (K) (Mg/m3) (GPa)
! 0    5735 13.09   1693 1253 1364 1813 154 1457 184
! 200  5729 13.08   1689 1251 1362 1809 154 1455 184
! 400  5711 13.05   1676 1243 1353 1795 151 1444 181
! 600  5682 13.01   1661 1232 1341 1779 150 1432 180
! 800  5642 12.95   1638 1214 1321 1755 148 1411 178
! 1000 5590 12.87   1606 1190 1295 1720 146 1383 175
! 1200 5527 12.77   1559 1155 1257 1670 141 1343 169
!

  c11=rho*vp*vp*1.d9/scale_fac
  c66=rho*vs*vs*1.d9/scale_fac

  A0=rho0*vp0*vp0*1.d9/scale_fac
  c33=c11+0.0349d0*A0
  c44=c66+0.00988d0*A0
  c12=c11-2.0d0*c66
  c13=c12-0.00881d0*A0

  end subroutine aniso_inner_core_model

