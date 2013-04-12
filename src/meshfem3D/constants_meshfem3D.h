!=====================================================================
!
!               S p e c f e m 3 D  V e r s i o n  2 . 1
!               ---------------------------------------
!
!          Main authors: Dimitri Komatitsch and Jeroen Tromp
!    Princeton University, USA and CNRS / INRIA / University of Pau
! (c) Princeton University / California Institute of Technology and CNRS / INRIA / University of Pau
!                             July 2012
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

!------------------------------------------------------
!----------- do not modify anything below -------------
!------------------------------------------------------


! number of GLL points not set in the mesher, do not modify this value
  integer, parameter :: NGLLX_M = 2
  integer, parameter :: NGLLY_M = NGLLX_M
  integer, parameter :: NGLLZ_M = NGLLX_M

! number of points per spectral element
  integer, parameter :: NGLLCUBE_M = NGLLX_M * NGLLY_M * NGLLZ_M

! define flag for elements
  integer, parameter :: IFLAG_ONE_LAYER_TOPOGRAPHY = 1
  integer, parameter :: IFLAG_BASEMENT_TOPO = 2

! flag for the four edges of each slice and for the bottom edge
  integer, parameter :: XI_MIN = 1
  integer, parameter :: XI_MAX = 2
  integer, parameter :: ETA_MIN = 3
  integer, parameter :: ETA_MAX = 4
  integer, parameter :: BOTTOM = 5

! size of topography and bathymetry file for Southern California
  integer, parameter :: NX_TOPO_SOCAL = 1401,NY_TOPO_SOCAL = 1001
  double precision, parameter :: ORIG_LAT_TOPO_SOCAL = 32.d0
  double precision, parameter :: ORIG_LONG_TOPO_SOCAL = -121.d0
  double precision, parameter :: DEGREES_PER_CELL_TOPO_SOCAL = 5.d0 / 1000.d0
  character(len=100), parameter :: TOPO_FILE_SOCAL = 'DATA/la_topography/topo_bathy_final.dat'

! size of Lupei Zhu's Moho map file for Southern California
  integer, parameter :: NX_MOHO = 71,NY_MOHO = 51
  double precision, parameter :: ORIG_LAT_MOHO = 32.d0
  double precision, parameter :: ORIG_LONG_MOHO = -121.d0
  double precision, parameter :: DEGREES_PER_CELL_MOHO = 0.1d0

! size of basement map file
  integer, parameter :: NX_BASEMENT = 161,NY_BASEMENT = 144
  double precision, parameter :: ORIG_X_BASEMENT = 316000.
  double precision, parameter :: ORIG_Y_BASEMENT = 3655000.
  double precision, parameter :: SPACING_X_BASEMENT = 1000.
  double precision, parameter :: SPACING_Y_BASEMENT = 1000.

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

! layers in the So-Cal regional model
! DEPTH_MOHO_SOCAL = -35 km was based on Dreger and Helmberger (1990)
! and is (July 2007) the preferred Moho depth for Dreger.
! The depth of 32 km is used in the standard processing (Wald et al., 1995)
! of SoCal events and is the value in the original Kanamori-Hadley (1975) model.
  double precision, parameter :: DEPTH_5p5km_SOCAL = -5500.d0
  double precision, parameter :: DEPTH_16km_SOCAL = -16000.d0
  double precision, parameter :: DEPTH_MOHO_SOCAL = -32000.d0

! reference surface of the model before adding topography
  double precision, parameter :: Z_SURFACE = 0.d0

! magic ratio for heuristic rule
! this gives 120 degree angles in doubling
! standard value 0.5 gives 135-135-90, which is not optimal
!  double precision, parameter :: MAGIC_RATIO = 0.6056d0

! type of elements for heuristic rule
  integer, parameter :: ITYPE_UNUSUAL_1  = 1
  integer, parameter :: ITYPE_UNUSUAL_1p = 2
  integer, parameter :: ITYPE_UNUSUAL_4  = 3
  integer, parameter :: ITYPE_UNUSUAL_4p = 4

! define number of spectral elements and points in basic symmetric mesh doubling superbrick
  integer, parameter :: NSPEC_DOUBLING_SUPERBRICK = 32
  integer, parameter :: NGLOB_DOUBLING_SUPERBRICK = 67
  integer, parameter :: NSPEC_SUPERBRICK_1L = 28
  integer, parameter :: NGLOB_SUPERBRICK_1L = 58

