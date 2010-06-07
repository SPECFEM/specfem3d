!=====================================================================
!
!               S p e c f e m 3 D  V e r s i o n  1 . 4
!               ---------------------------------------
!
!                 Dimitri Komatitsch and Jeroen Tromp
!    Seismological Laboratory - California Institute of Technology
!         (c) California Institute of Technology September 2006
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

  subroutine hauksson_model(vp,vs,utm_x_eval,utm_y_eval,z_eval,vp_final,vs_final,MOHO_MAP_LUPEI)

  implicit none

  include "constants.h"
  include "constants_gocad.h"

!! DK DK UGLY one day, we should clarify the issue of merging Hauksson's Moho
!! DK DK UGLY with our Lupei Moho. Should not be a big issue because in
!! DK DK UGLY principle Hauksson used Lupei's map to build his model

  double precision utm_x_eval,utm_y_eval,z_eval
  double precision vp_final,vs_final
  logical MOHO_MAP_LUPEI

  double precision, dimension(NLAYERS_HAUKSSON,NGRID_NEW_HAUKSSON,NGRID_NEW_HAUKSSON) :: vp,vs
  double precision, dimension(NLAYERS_HAUKSSON) :: vp_interp,vs_interp

  integer ilayer
  integer icell_interp_x,icell_interp_y
  double precision spacing_x,spacing_y
  double precision utm_x_eval_copy,utm_y_eval_copy
  double precision gamma_interp_x,gamma_interp_y,gamma_interp_z
  double precision v1,v2,v3,v4
  double precision vp_upper,vs_upper,vp_lower,vs_lower,z_upper,z_lower

! copy input values
  utm_x_eval_copy = utm_x_eval
  utm_y_eval_copy = utm_y_eval

! make sure we stay inside Hauksson's grid
  if(utm_x_eval_copy < UTM_X_ORIG_HAUKSSON) utm_x_eval_copy = UTM_X_ORIG_HAUKSSON
  if(utm_y_eval_copy < UTM_Y_ORIG_HAUKSSON) utm_y_eval_copy = UTM_Y_ORIG_HAUKSSON

! determine spacing and cell for linear interpolation
  spacing_x = (utm_x_eval_copy - UTM_X_ORIG_HAUKSSON) / SPACING_UTM_X_HAUKSSON
  spacing_y = (utm_y_eval_copy - UTM_Y_ORIG_HAUKSSON) / SPACING_UTM_Y_HAUKSSON

  icell_interp_x = int(spacing_x) + 1
  icell_interp_y = int(spacing_y) + 1

  gamma_interp_x = spacing_x - int(spacing_x)
  gamma_interp_y = spacing_y - int(spacing_y)

! suppress edge effects for points outside of Hauksson's model
  if(icell_interp_x < 1) then
    icell_interp_x = 1
    gamma_interp_x = 0.d0
  endif
  if(icell_interp_x > NGRID_NEW_HAUKSSON-1) then
    icell_interp_x = NGRID_NEW_HAUKSSON-1
    gamma_interp_x = 1.d0
  endif

  if(icell_interp_y < 1) then
    icell_interp_y = 1
    gamma_interp_y = 0.d0
  endif
  if(icell_interp_y > NGRID_NEW_HAUKSSON-1) then
    icell_interp_y = NGRID_NEW_HAUKSSON-1
    gamma_interp_y = 1.d0
  endif

! make sure interpolation makes sense
  if(gamma_interp_x < -0.001d0 .or. gamma_interp_x > 1.001d0) &
        stop 'interpolation in x is incorrect in Hauksson'
  if(gamma_interp_y < -0.001d0 .or. gamma_interp_y > 1.001d0) &
        stop 'interpolation in y is incorrect in Hauksson'

! interpolate Hauksson's model at right location using bilinear interpolation
  do ilayer = 1,NLAYERS_HAUKSSON

! for Vp
  v1 = vp(ilayer,icell_interp_x,icell_interp_y)
  v2 = vp(ilayer,icell_interp_x+1,icell_interp_y)
  v3 = vp(ilayer,icell_interp_x+1,icell_interp_y+1)
  v4 = vp(ilayer,icell_interp_x,icell_interp_y+1)

  vp_interp(ilayer) = v1*(1.-gamma_interp_x)*(1.-gamma_interp_y) + &
                      v2*gamma_interp_x*(1.-gamma_interp_y) + &
                      v3*gamma_interp_x*gamma_interp_y + &
                      v4*(1.-gamma_interp_x)*gamma_interp_y

! for Vs
  v1 = vs(ilayer,icell_interp_x,icell_interp_y)
  v2 = vs(ilayer,icell_interp_x+1,icell_interp_y)
  v3 = vs(ilayer,icell_interp_x+1,icell_interp_y+1)
  v4 = vs(ilayer,icell_interp_x,icell_interp_y+1)

  vs_interp(ilayer) = v1*(1.-gamma_interp_x)*(1.-gamma_interp_y) + &
                      v2*gamma_interp_x*(1.-gamma_interp_y) + &
                      v3*gamma_interp_x*gamma_interp_y + &
                      v4*(1.-gamma_interp_x)*gamma_interp_y

  enddo

! choose right values depending on depth of target point
  if(z_eval >= Z_HAUKSSON_LAYER_1) then
    vp_final = vp_interp(1)
    vs_final = vs_interp(1)
    return

  else if(z_eval <= Z_HAUKSSON_LAYER_9) then
    vp_final = vp_interp(9)
    vs_final = vs_interp(9)
    return

  else if(z_eval >= Z_HAUKSSON_LAYER_2) then
    vp_upper = vp_interp(1)
    vs_upper = vs_interp(1)
    z_upper = Z_HAUKSSON_LAYER_1

    vp_lower = vp_interp(2)
    vs_lower = vs_interp(2)
    z_lower = Z_HAUKSSON_LAYER_2

  else if(z_eval >= Z_HAUKSSON_LAYER_3) then
    vp_upper = vp_interp(2)
    vs_upper = vs_interp(2)
    z_upper = Z_HAUKSSON_LAYER_2

    vp_lower = vp_interp(3)
    vs_lower = vs_interp(3)
    z_lower = Z_HAUKSSON_LAYER_3

  else if(z_eval >= Z_HAUKSSON_LAYER_4) then
    vp_upper = vp_interp(3)
    vs_upper = vs_interp(3)
    z_upper = Z_HAUKSSON_LAYER_3

    vp_lower = vp_interp(4)
    vs_lower = vs_interp(4)
    z_lower = Z_HAUKSSON_LAYER_4

  else if(z_eval >= Z_HAUKSSON_LAYER_5) then
    vp_upper = vp_interp(4)
    vs_upper = vs_interp(4)
    z_upper = Z_HAUKSSON_LAYER_4

    vp_lower = vp_interp(5)
    vs_lower = vs_interp(5)
    z_lower = Z_HAUKSSON_LAYER_5

 else if(z_eval >= Z_HAUKSSON_LAYER_6) then
    vp_upper = vp_interp(5)
    vs_upper = vs_interp(5)
    z_upper = Z_HAUKSSON_LAYER_5

    vp_lower = vp_interp(6)
    vs_lower = vs_interp(6)
    z_lower = Z_HAUKSSON_LAYER_6

  else if(z_eval >= Z_HAUKSSON_LAYER_7) then
    vp_upper = vp_interp(6)
    vs_upper = vs_interp(6)
    z_upper = Z_HAUKSSON_LAYER_6

    vp_lower = vp_interp(7)
    vs_lower = vs_interp(7)
    z_lower = Z_HAUKSSON_LAYER_7

  else if(z_eval >= Z_HAUKSSON_LAYER_8) then
    vp_upper = vp_interp(7)
    vs_upper = vs_interp(7)
    z_upper = Z_HAUKSSON_LAYER_7

    vp_lower = vp_interp(8)
    vs_lower = vs_interp(8)
    z_lower = Z_HAUKSSON_LAYER_8

  else
    if(.not. MOHO_MAP_LUPEI) then
      vp_upper = vp_interp(8)
      vs_upper = vs_interp(8)
      z_upper = Z_HAUKSSON_LAYER_8

      vp_lower = vp_interp(9)
      vs_lower = vs_interp(9)
      z_lower = Z_HAUKSSON_LAYER_9
   !!! waiting for better interpolation of Moho maps.
    else
      vp_upper = vp_interp(8)
      vs_upper = vs_interp(8)
      z_upper = Z_HAUKSSON_LAYER_8

      vp_lower = vp_interp(9)
      vs_lower = vs_interp(9)
      z_lower = Z_HAUKSSON_LAYER_9
    endif

  endif

    gamma_interp_z = (z_eval - z_lower) / (z_upper - z_lower)

    if(gamma_interp_z < -0.001d0 .or. gamma_interp_z > 1.001d0) &
        stop 'interpolation in z is incorrect in Hauksson'

    vp_final = vp_upper * gamma_interp_z + vp_lower * (1.-gamma_interp_z)
    vs_final = vs_upper * gamma_interp_z + vs_lower * (1.-gamma_interp_z)

  end subroutine hauksson_model

