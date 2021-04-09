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
  real(kind=CUSTOM_REAL) :: CPML_width_x,CPML_width_y,CPML_width_z

  ! C-PML damping profile arrays
  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: d_store_x, d_store_y, d_store_z

  ! auxiliary parameters arrays
  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: K_store_x, K_store_y, K_store_z
  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: alpha_store_x,alpha_store_y,alpha_store_z

  real(kind=CUSTOM_REAL), dimension(:,:,:,:,:), allocatable :: convolution_coef_acoustic_alpha,convolution_coef_acoustic_beta

  ! minimum distance between parameters of CPML to avoid the singularities
  real(kind=CUSTOM_REAL) :: min_distance_between_CPML_parameter

  ! store the field of displ + (1-2 * \theta)/2*deltat * veloc + (1-\theta)/2*deltat**2 * accel
  ! for second order convolution scheme
  ! where displ, veloc, accel are defined at n-1 time step
  real(kind=CUSTOM_REAL), dimension(:,:,:,:,:), allocatable :: PML_displ_old
  real(kind=CUSTOM_REAL), parameter :: theta = 1.0_CUSTOM_REAL/8.0_CUSTOM_REAL

  ! store the field of displ + (1-2 * \theta)/2*deltat * veloc for second order convolution scheme
  ! where displ is defined at n time step, while veloc is predicted veloc at "n" time step
  real(kind=CUSTOM_REAL), dimension(:,:,:,:,:), allocatable :: PML_displ_new

  !store the field of
  !potential_acoustic + (1-2 * \theta)/2*deltat * potential_dot_acoustic + (1-\theta)/2*deltat**2 * potential_dot_dot_acoustic
  !where potential_acoustic, potential_dot_acoustic, potential_dot_dot_acoustic are defined at "n-1" time step
  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: PML_potential_acoustic_old

  !store the field of
  !potential_acoustic + (1-2 * \theta)/2*deltat * potential_dot_acoustic
  !potential_acoustic is defined at "n" time step,
  !and potential_dot_acoustic is predicted potential_dot_acoustic at "n" time step
  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: PML_potential_acoustic_new

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

  !store the field accel at time step n-1 for second order convolution scheme
!  real(kind=CUSTOM_REAL), dimension(:), allocatable :: potential_dot_dot_acoustic_old

  real(kind=CUSTOM_REAL), dimension(:,:,:,:,:), allocatable :: rmemory_dpotential_dxl
  real(kind=CUSTOM_REAL), dimension(:,:,:,:,:), allocatable :: rmemory_dpotential_dyl
  real(kind=CUSTOM_REAL), dimension(:,:,:,:,:), allocatable :: rmemory_dpotential_dzl

  ! C-PML memory variable needed for displacement
  real(kind=CUSTOM_REAL), dimension(:,:,:,:,:,:), allocatable :: rmemory_displ_elastic

  ! C-PML memory variable needed for potential
  real(kind=CUSTOM_REAL), dimension(:,:,:,:,:), allocatable :: rmemory_potential_acoustic

  ! C-PML contribution to update displacement on elastic/acoustic interface
  real(kind=CUSTOM_REAL), dimension(:,:,:,:,:,:), allocatable :: rmemory_coupling_ac_el_displ

  ! C-PML contribution to update displacement on elastic/acoustic interface
  real(kind=CUSTOM_REAL), dimension(:,:,:,:,:,:), allocatable :: rmemory_coupling_el_ac_potential
  real(kind=CUSTOM_REAL), dimension(:,:,:,:,:,:), allocatable :: rmemory_coupling_el_ac_potential_dot_dot

  ! --------------------------------------------------------------------------------------------
  ! for adjoint tomography
  ! need the array stored the points on interface between PML and interior computational domain
  ! --------------------------------------------------------------------------------------------
  integer :: nglob_interface_PML_acoustic,nglob_interface_PML_elastic
  integer, dimension(:), allocatable :: points_interface_PML_acoustic, points_interface_PML_elastic

  integer :: b_reclen_PML_field,b_reclen_PML_potential
  real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: b_PML_field,b_PML_potential

end module pml_par
