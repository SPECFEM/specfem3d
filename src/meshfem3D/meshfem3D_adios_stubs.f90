!=====================================================================
!
!               S p e c f e m 3 D  V e r s i o n  3 . 0
!               ---------------------------------------
!
!     Main historical authors: Dimitri Komatitsch and Jeroen Tromp
!                              CNRS, France
!                       and Princeton University, USA
!                 (there are currently many more authors!)
!                           (c) October 2017
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

!==============================================================================
!> \file meshfem3D_adios_stubs.f90
!!
!!  Stubs for ADIOS functions. Avoid link error when not configured with
!!  ADIOS.
!!
!! \author MPBL
!==============================================================================


  subroutine save_databases_adios(LOCAL_PATH,sizeprocs, &
                                  nspec,nglob, &
                                  iMPIcut_xi,iMPIcut_eta, &
                                  nodes_coords,ispec_material_id, &
                                  nspec2D_xmin,nspec2D_xmax,nspec2D_ymin,nspec2D_ymax, &
                                  ibelm_xmin,ibelm_xmax,ibelm_ymin,ibelm_ymax,ibelm_bottom,ibelm_top)

  use constants

  use meshfem3D_par, only: NSPEC2DMAX_XMIN_XMAX,NSPEC2DMAX_YMIN_YMAX,NSPEC2D_BOTTOM,NSPEC2D_TOP

  use adios_manager_mod

  implicit none

  character(len=MAX_STRING_LEN) :: LOCAL_PATH
  integer :: sizeprocs

  integer :: nspec
  integer :: nglob

  logical, dimension(2,nspec):: iMPIcut_xi,iMPIcut_eta
  double precision, dimension(nglob,NDIM) :: nodes_coords
  integer, dimension(nspec) :: ispec_material_id

  integer :: nspec2D_xmin,nspec2D_xmax,nspec2D_ymin,nspec2D_ymax

  integer, dimension(NSPEC2DMAX_XMIN_XMAX) :: ibelm_xmin,ibelm_xmax
  integer, dimension(NSPEC2DMAX_YMIN_YMAX) :: ibelm_ymin,ibelm_ymax
  integer, dimension(NSPEC2D_BOTTOM) :: ibelm_bottom
  integer, dimension(NSPEC2D_TOP) :: ibelm_top

  ! local parameters
  integer(kind=4) :: unused_i4
  logical :: unused_l
  double precision :: unused_dp

  unused_i4 = len_trim(LOCAL_PATH)
  unused_i4 = sizeprocs

  unused_l  = iMPIcut_xi(1,1)
  unused_l  = iMPIcut_eta(1,1)

  unused_dp = nodes_coords(1,1)
  unused_i4 = ispec_material_id(1)

  unused_i4 = nspec2D_xmin
  unused_i4 = nspec2D_xmax
  unused_i4 = nspec2D_ymin
  unused_i4 = nspec2D_ymax
  unused_i4 = ibelm_xmin(1)
  unused_i4 = ibelm_xmax(1)
  unused_i4 = ibelm_ymin(1)
  unused_i4 = ibelm_ymax(1)
  unused_i4 = ibelm_bottom(1)
  unused_i4 = ibelm_top(1)

  call no_adios_err()

end subroutine save_databases_adios
