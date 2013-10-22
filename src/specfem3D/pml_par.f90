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
!
! United States and French Government Sponsorship Acknowledged.

module pml_par

  ! main parameter module for C-PML simulations

  use constants, only: CUSTOM_REAL

  implicit none

  ! number of C-PML spectral elements
  integer :: NSPEC_CPML

  ! C-PML spectral elements global indexing
  integer, dimension(:), allocatable :: CPML_to_spec

  ! C-PML spectral elements local indexing
  integer, dimension(:), allocatable :: spec_to_CPML

  ! C-PML regions
  integer, dimension(:), allocatable :: CPML_regions

  ! C-PML element type array: 1 = face, 2 = edge, 3 = corner
  integer, dimension(:), allocatable :: CPML_type

  ! mask of C-PML elements for the global mesh
  logical, dimension(:), allocatable :: is_CPML

  ! thickness of C-PML layers
  real(CUSTOM_REAL) :: CPML_width_x,CPML_width_y,CPML_width_z

  ! C-PML damping profile arrays
  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: d_store_x, d_store_y, d_store_z

  ! auxiliary parameters arrays
  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: K_store_x, K_store_y, K_store_z
  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: alpha_store_x,alpha_store_y,alpha_store_z

  !store the field of displ + 2 * deltat**2 * accel at time step n-1 for second order convolution scheme
  real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: displ_old

  ! derivatives of ux, uy and uz with respect to x, y and z
  real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable :: PML_dux_dxl,PML_dux_dyl,PML_dux_dzl
  real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable :: PML_duy_dxl,PML_duy_dyl,PML_duy_dzl
  real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable :: PML_duz_dxl,PML_duz_dyl,PML_duz_dzl

  real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable :: PML_dux_dxl_old,PML_dux_dyl_old,PML_dux_dzl_old
  real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable :: PML_duy_dxl_old,PML_duy_dyl_old,PML_duy_dzl_old
  real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable :: PML_duz_dxl_old,PML_duz_dyl_old,PML_duz_dzl_old

  ! derivatives of potential with respect to x, y and z
  real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable :: PML_dpotential_dxl,PML_dpotential_dyl,PML_dpotential_dzl

  real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable :: PML_dpotential_dxl_old,PML_dpotential_dyl_old,PML_dpotential_dzl_old

  ! C-PML memory variables
  real(kind=CUSTOM_REAL), dimension(:,:,:,:,:), allocatable :: rmemory_dux_dxl_x,rmemory_dux_dyl_x,rmemory_dux_dzl_x
  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: rmemory_duy_dxl_x,rmemory_duy_dyl_x
  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: rmemory_duz_dxl_x,rmemory_duz_dzl_x

  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: rmemory_dux_dxl_y,rmemory_dux_dyl_y
  real(kind=CUSTOM_REAL), dimension(:,:,:,:,:), allocatable :: rmemory_duy_dxl_y,rmemory_duy_dyl_y,rmemory_duy_dzl_y
  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: rmemory_duz_dzl_y,rmemory_duz_dyl_y

  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: rmemory_dux_dxl_z,rmemory_dux_dzl_z
  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: rmemory_duy_dyl_z,rmemory_duy_dzl_z
  real(kind=CUSTOM_REAL), dimension(:,:,:,:,:), allocatable :: rmemory_duz_dxl_z,rmemory_duz_dyl_z,rmemory_duz_dzl_z

  !store the field of displ + 2 * deltat**2 * accel at time step n-1 for second order convolution scheme
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: potential_acoustic_old

  !store the field accel at time step n-1 for second order convolution scheme
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: potential_dot_dot_acoustic_old

  real(kind=CUSTOM_REAL), dimension(:,:,:,:,:), allocatable :: rmemory_dpotential_dxl
  real(kind=CUSTOM_REAL), dimension(:,:,:,:,:), allocatable :: rmemory_dpotential_dyl
  real(kind=CUSTOM_REAL), dimension(:,:,:,:,:), allocatable :: rmemory_dpotential_dzl

  ! C-PML memory variable needed for displacement
  real(kind=CUSTOM_REAL), dimension(:,:,:,:,:,:), allocatable :: rmemory_displ_elastic

  ! C-PML memory variable needed for potential
  real(kind=CUSTOM_REAL), dimension(:,:,:,:,:), allocatable :: rmemory_potential_acoustic

  ! C-PML contribution to update acceleration to the global mesh
  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: accel_elastic_CPML

  ! C-PML contribution to update the second derivative of the potential to the global mesh
  real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable :: potential_dot_dot_acoustic_CPML

  ! C-PML contribution to update displacement on elastic/acoustic interface
  real(kind=CUSTOM_REAL), dimension(:,:,:,:,:,:), allocatable :: rmemory_coupling_ac_el_displ

  ! C-PML contribution to update displacement on elastic/acoustic interface
  real(kind=CUSTOM_REAL), dimension(:,:,:,:,:), allocatable :: rmemory_coupling_el_ac_potential
  real(kind=CUSTOM_REAL), dimension(:,:,:,:,:), allocatable :: rmemory_coupling_el_ac_potential_dot_dot

  ! --------------------------------------------------------------------------------------------
  ! for adjoint tomography
  ! need the array stored the points on interface between PML and interior computational domain
  ! --------------------------------------------------------------------------------------------
  integer :: nglob_interface_PML_acoustic,nglob_interface_PML_elastic
  integer, dimension(:), allocatable :: points_interface_PML_acoustic, points_interface_PML_elastic

  integer :: b_reclen_PML_field,b_reclen_PML_potential
  real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: b_PML_field,b_PML_potential

end module pml_par
