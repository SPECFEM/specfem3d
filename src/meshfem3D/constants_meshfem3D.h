!=====================================================================
!
!               S p e c f e m 3 D  V e r s i o n  3 . 0
!               ---------------------------------------
!
!     Main historical authors: Dimitri Komatitsch and Jeroen Tromp
!                        Princeton University, USA
!                and CNRS / University of Marseille, France
!                 (there are currently many more authors!)
! (c) Princeton University and CNRS / University of Marseille, July 2012
!
! This program is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 3 of the License, or
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

! number of material properties + element ID in input file Mesh_Par_file
  integer, parameter :: NUMBER_OF_MATERIAL_PROPERTIES = 18

!------------------------------------------------------
!----------- do not modify anything below -------------
!------------------------------------------------------

! number of GLL points not set in the mesher, do not modify this value
! we can use 2 for NGNOD == 8 only and faster meshing; or 3 to allow also for NGNOD == 27 meshing
  integer, parameter :: NGLLX_M = 3
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

