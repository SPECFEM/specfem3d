!=====================================================================
!
!               S p e c f e m 3 D  V e r s i o n  2 . 1
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
   NPROC_XI,NPROC_ETA,addressing,iMPIcut_xi,iMPIcut_eta,&
   ibool,nodes_coords,true_material_num, &
   nspec2D_xmin,nspec2D_xmax,nspec2D_ymin,nspec2D_ymax, &
   NSPEC2D_BOTTOM,NSPEC2D_TOP, NSPEC2DMAX_XMIN_XMAX,NSPEC2DMAX_YMIN_YMAX, &
   ibelm_xmin,ibelm_xmax,ibelm_ymin,ibelm_ymax,ibelm_bottom,ibelm_top,&
   NMATERIALS,material_properties)

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
  integer, dimension(nspec) :: true_material_num
  integer :: nspec2D_xmin,nspec2D_xmax,nspec2D_ymin,nspec2D_ymax
  integer :: NSPEC2D_BOTTOM,NSPEC2D_TOP,NSPEC2DMAX_XMIN_XMAX,NSPEC2DMAX_YMIN_YMAX
  integer, dimension(NSPEC2DMAX_XMIN_XMAX) :: ibelm_xmin,ibelm_xmax
  integer, dimension(NSPEC2DMAX_YMIN_YMAX) :: ibelm_ymin,ibelm_ymax
  integer, dimension(NSPEC2D_BOTTOM) :: ibelm_bottom
  integer, dimension(NSPEC2D_TOP) :: ibelm_top
  integer :: NMATERIALS
  double precision, dimension(NMATERIALS,6) ::  material_properties

  call no_adios_err()
end subroutine save_databases_adios
