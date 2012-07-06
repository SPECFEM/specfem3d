!=====================================================================
!
!               S p e c f e m 3 D  V e r s i o n  2 . 1
!               ---------------------------------------
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

! depth at which we start to honor the basement interface
  double precision, parameter :: Z_THRESHOLD_HONOR_BASEMENT = -4700.d0

! min and max density in the model
  double precision, parameter :: DENSITY_MAX = 3000.d0
  double precision, parameter :: DENSITY_MIN = 2000.d0

! extend model below threshold and above topography to make sure
! there is no small gap between interpolated maps and sediments
  logical, parameter :: EXTEND_VOXET_BELOW_BASEMENT = .true.
  logical, parameter :: EXTEND_VOXET_ABOVE_TOPO = .true.
  double precision, parameter :: DISTMAX_ASSUME_SEDIMENTS = 210.d0
  integer, parameter :: NCELLS_EXTEND = 8

! define flag for elements - socal model
  integer, parameter :: IFLAG_ONE_LAYER_TOPOGRAPHY = 1
  integer, parameter :: IFLAG_BASEMENT_TOPO = 2
  integer, parameter :: IFLAG_16km_BASEMENT = 3
  integer, parameter :: IFLAG_MOHO_16km = 4
  integer, parameter :: IFLAG_HALFSPACE_MOHO = 5

!
!--- larger Hauksson model for entire So-Cal, 15 km resolution
!

! number of non-constant layers
integer, parameter :: NLAYERS_HAUKSSON = 9

! depth of layers
double precision, parameter :: Z_HAUKSSON_LAYER_1 =  -1000.d0
double precision, parameter :: Z_HAUKSSON_LAYER_2 =  -4000.d0
double precision, parameter :: Z_HAUKSSON_LAYER_3 =  -6000.d0
double precision, parameter :: Z_HAUKSSON_LAYER_4 = -10000.d0
double precision, parameter :: Z_HAUKSSON_LAYER_5 = -15000.d0
double precision, parameter :: Z_HAUKSSON_LAYER_6 = -17000.d0
double precision, parameter :: Z_HAUKSSON_LAYER_7 = -22000.d0
double precision, parameter :: Z_HAUKSSON_LAYER_8 = -31000.d0
double precision, parameter :: Z_HAUKSSON_LAYER_9 = -33000.d0

integer, parameter :: NGRID_NEW_HAUKSSON = 201

! corners of new Hauksson's interpolated grid
double precision, parameter :: UTM_X_ORIG_HAUKSSON = 122035.012d0
double precision, parameter :: UTM_X_END_HAUKSSON  = 766968.628d0
double precision, parameter :: UTM_Y_ORIG_HAUKSSON = 3547232.986d0
double precision, parameter :: UTM_Y_END_HAUKSSON  = 4098868.501d0

double precision, parameter :: SPACING_UTM_X_HAUKSSON = (UTM_X_END_HAUKSSON - UTM_X_ORIG_HAUKSSON) / (NGRID_NEW_HAUKSSON-1.d0)
double precision, parameter :: SPACING_UTM_Y_HAUKSSON = (UTM_Y_END_HAUKSSON - UTM_Y_ORIG_HAUKSSON) / (NGRID_NEW_HAUKSSON-1.d0)

!
! new Gocad Voxets Peter July 29, 2002 - high-res and medium-res blocks
!

! size of the medium-resolution Gocad voxet
integer, parameter :: NX_GOCAD_MR = 194, NY_GOCAD_MR = 196, NZ_GOCAD_MR = 100

double precision, parameter :: ORIG_X_GOCAD_MR = 283000.
double precision, parameter :: ORIG_Y_GOCAD_MR = 3655000.
double precision, parameter :: ORIG_Z_GOCAD_MR = -15000.

double precision, parameter :: SPACING_X_GOCAD_MR = 1000.
double precision, parameter :: SPACING_Y_GOCAD_MR = 1000.
double precision, parameter :: SPACING_Z_GOCAD_MR = 200.

! maximum size of model for tapering of transition between Hauksson and MR
double precision, parameter :: END_X_GOCAD_MR = ORIG_X_GOCAD_MR + SPACING_X_GOCAD_MR * (NX_GOCAD_MR - 1)
double precision, parameter :: END_Y_GOCAD_MR = ORIG_Y_GOCAD_MR + SPACING_Y_GOCAD_MR * (NY_GOCAD_MR - 1)

! size of the high-resolution Gocad voxet
integer, parameter :: NX_GOCAD_HR = 185, NY_GOCAD_HR = 196, NZ_GOCAD_HR = 100

double precision, parameter :: ORIG_X_GOCAD_HR = 371052.25
double precision, parameter :: ORIG_Y_GOCAD_HR = 3725250.
double precision, parameter :: ORIG_Z_GOCAD_HR = -9500.

double precision, parameter :: SPACING_X_GOCAD_HR = 250.
double precision, parameter :: SPACING_Y_GOCAD_HR = 250.
double precision, parameter :: SPACING_Z_GOCAD_HR = 100.

! maximum size of model for tapering of transition between HR and MR
double precision, parameter :: END_X_GOCAD_HR = ORIG_X_GOCAD_HR + SPACING_X_GOCAD_HR * (NX_GOCAD_HR - 1)
double precision, parameter :: END_Y_GOCAD_HR = ORIG_Y_GOCAD_HR + SPACING_Y_GOCAD_HR * (NY_GOCAD_HR - 1)

! implement smooth transition between Hauksson, HR and MR Gocad blocks
logical, parameter :: TAPER_GOCAD_TRANSITIONS = .true.

!  Salton Sea Gocad voxet
integer, parameter :: GOCAD_ST_NU = 638, GOCAD_ST_NV = 219, GOCAD_ST_NW = 76
double precision, parameter :: GOCAD_ST_O_X = 720844.0, GOCAD_ST_O_Y = 3401799.250, &
GOCAD_ST_O_Z =      -6354.334
double precision, parameter :: GOCAD_ST_U_X = -209197.89, GOCAD_ST_U_Y =  320741.71
double precision, parameter :: GOCAD_ST_V_X = 109670.74, GOCAD_ST_V_Y = 71530.72
double precision, parameter :: GOCAD_ST_W_Z =  7666.334
double precision, parameter :: GOCAD_ST_NO_DATA_VALUE = -99999

! unused parameters
!--------------------------
! size of Lupei Zhu's Moho map file for Southern California -- unused
!  integer, parameter :: NX_MOHO = 71,NY_MOHO = 51
!  double precision, parameter :: ORIG_LAT_MOHO = 32.d0
!  double precision, parameter :: ORIG_LONG_MOHO = -121.d0
!  double precision, parameter :: DEGREES_PER_CELL_MOHO = 0.1d0
!
! layers in the So-Cal regional model -- unused
! DEPTH_MOHO_SOCAL = -35 km was based on Dreger and Helmberger (1990)
! and is (July 2007) the preferred Moho depth for Dreger.
! The depth of 32 km is used in the standard processing (Wald et al., 1995)
! of SoCal events and is the value in the original Kanamori-Hadley (1975) model.
! double precision, parameter :: DEPTH_5p5km_SOCAL = -5500.d0
! double precision, parameter :: DEPTH_16km_SOCAL = -16000.d0
! double precision, parameter :: DEPTH_MOHO_SOCAL = -32000.d0
!
! magic ratio for heuristic rule -- unused
! this gives 120 degree angles in doubling
! standard value 0.5 gives 135-135-90, which is not optimal
!  double precision, parameter :: MAGIC_RATIO = 0.6056d0
