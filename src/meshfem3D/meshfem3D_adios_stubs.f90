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

subroutine save_databases_adios(LOCAL_PATH, myrank, sizeprocs, &
                                nspec,nglob,iproc_xi,iproc_eta, &
                                NPROC_XI,NPROC_ETA,addressing,iMPIcut_xi,iMPIcut_eta, &
                                ibool,nodes_coords,ispec_material_id, &
                                nspec2D_xmin,nspec2D_xmax,nspec2D_ymin,nspec2D_ymax, &
                                NSPEC2D_BOTTOM,NSPEC2D_TOP, NSPEC2DMAX_XMIN_XMAX,NSPEC2DMAX_YMIN_YMAX, &
                                ibelm_xmin,ibelm_xmax,ibelm_ymin,ibelm_ymax,ibelm_bottom,ibelm_top, &
                                NMATERIALS,material_properties, &
                                nspec_CPML,CPML_to_spec,CPML_regions,is_CPML)

  use constants
  use adios_manager_mod

  implicit none
  include "constants_meshfem3D.h"
  character(len=MAX_STRING_LEN) :: LOCAL_PATH
  integer :: myrank, sizeprocs
  integer :: nspec
  integer :: nglob
  integer :: iproc_xi,iproc_eta
  integer :: NPROC_XI,NPROC_ETA
  integer, dimension(0:NPROC_XI-1,0:NPROC_ETA-1) :: addressing
  logical, dimension(2,nspec):: iMPIcut_xi,iMPIcut_eta
  integer, dimension(NGLLX_M,NGLLY_M,NGLLZ_M,nspec) :: ibool
  double precision, dimension(nglob,3) :: nodes_coords
  integer, dimension(nspec) :: ispec_material_id
  integer :: nspec2D_xmin,nspec2D_xmax,nspec2D_ymin,nspec2D_ymax
  integer :: NSPEC2D_BOTTOM,NSPEC2D_TOP,NSPEC2DMAX_XMIN_XMAX,NSPEC2DMAX_YMIN_YMAX
  integer, dimension(NSPEC2DMAX_XMIN_XMAX) :: ibelm_xmin,ibelm_xmax
  integer, dimension(NSPEC2DMAX_YMIN_YMAX) :: ibelm_ymin,ibelm_ymax
  integer, dimension(NSPEC2D_BOTTOM) :: ibelm_bottom
  integer, dimension(NSPEC2D_TOP) :: ibelm_top
  integer :: NMATERIALS
  double precision, dimension(NMATERIALS,7) ::  material_properties

 ! CPML
  integer, intent(in) :: nspec_CPML
  integer, dimension(nspec_CPML), intent(in) :: CPML_to_spec,CPML_regions
  logical, dimension(nspec), intent(in) :: is_CPML

  integer(kind=4) :: unused_i4
  logical :: unused_l
  double precision :: unused_dp

  unused_i4 = len_trim(LOCAL_PATH)
  unused_i4 = myrank
  unused_i4 = sizeprocs
  unused_i4 = iproc_xi
  unused_i4 = iproc_eta
  unused_i4 = addressing(0,0)
  unused_l  = iMPIcut_xi(1,1)
  unused_l  = iMPIcut_eta(1,1)
  unused_i4 = ibool(1,1,1,1)
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
  unused_dp = material_properties(1,1)
  unused_i4 = nspec_CPML
  unused_i4 = CPML_to_spec(1)
  unused_i4 = CPML_regions(1)
  unused_l  = is_CPML(1)

  call no_adios_err()

end subroutine save_databases_adios
