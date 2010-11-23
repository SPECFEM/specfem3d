!=====================================================================
!
!          S p e c f e m 3 D  B a s i n  V e r s i o n  1 . 4
!          --------------------------------------------------
!
!          Main authors: Dimitri Komatitsch and Jeroen Tromp
!    Princeton University, USA and University of Pau / CNRS / INRIA
!         (c) California Institute of Technology July 2005
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

! constants.h.  Generated from constants.h.in by configure.

!
!--- larger Hauksson model for entire So-Cal, 15 km resolution
!

!! Hauksson (2000)
!! number of non-constant layers
!  integer, parameter :: NLAYERS_HAUKSSON = 9
!! depth of layers
!  double precision, parameter :: Z_HAUKSSON_LAYER_1 =  -1000.d0
!  double precision, parameter :: Z_HAUKSSON_LAYER_2 =  -4000.d0
!  double precision, parameter :: Z_HAUKSSON_LAYER_3 =  -6000.d0
!  double precision, parameter :: Z_HAUKSSON_LAYER_4 = -10000.d0
!  double precision, parameter :: Z_HAUKSSON_LAYER_5 = -15000.d0
!  double precision, parameter :: Z_HAUKSSON_LAYER_6 = -17000.d0
!  double precision, parameter :: Z_HAUKSSON_LAYER_7 = -22000.d0
!  double precision, parameter :: Z_HAUKSSON_LAYER_8 = -31000.d0
!  double precision, parameter :: Z_HAUKSSON_LAYER_9 = -33000.d0
!
!  integer, parameter :: NGRID_NEW_HAUKSSON = 201
!
!! corners of new Hauksson's interpolated grid
!  double precision, parameter :: UTM_X_ORIG_HAUKSSON = 122035.012d0
!  double precision, parameter :: UTM_X_END_HAUKSSON  = 766968.628d0
!  double precision, parameter :: UTM_Y_ORIG_HAUKSSON = 3547232.986d0
!  double precision, parameter :: UTM_Y_END_HAUKSSON  = 4098868.501d0

! Lin-Shearer-Hauksson-Thurber (2007)
! number of non-constant layers
   integer, parameter :: NLAYERS_HAUKSSON = 8
! depth of layers
  double precision, parameter :: Z_HAUKSSON_LAYER_1 =  0.d0
  double precision, parameter :: Z_HAUKSSON_LAYER_2 =  -3000.d0
  double precision, parameter :: Z_HAUKSSON_LAYER_3 =  -6000.d0
  double precision, parameter :: Z_HAUKSSON_LAYER_4 = -10000.d0
  double precision, parameter :: Z_HAUKSSON_LAYER_5 = -15000.d0
  double precision, parameter :: Z_HAUKSSON_LAYER_6 = -17000.d0
  double precision, parameter :: Z_HAUKSSON_LAYER_7 = -22000.d0
  double precision, parameter :: Z_HAUKSSON_LAYER_8 = -31000.d0

  integer, parameter :: NGRID_NEW_HAUKSSON = 241

! corners of new Hauksson's interpolated grid
  double precision, parameter :: UTM_X_ORIG_HAUKSSON = 88021.568d0
  double precision, parameter :: UTM_X_END_HAUKSSON  = 861517.886d0
  double precision, parameter :: UTM_Y_ORIG_HAUKSSON = 3404059.875d0
  double precision, parameter :: UTM_Y_END_HAUKSSON  = 4180234.582d0
