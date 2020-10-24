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

! for external mesh

  subroutine save_arrays_solver_ext_mesh_h5(nspec,nglob,APPROXIMATE_OCEAN_LOAD,ibool, &
                    num_interfaces_ext_mesh,my_neighbors_ext_mesh,nibool_interfaces_ext_mesh, &
                    max_interface_size_ext_mesh,ibool_interfaces_ext_mesh, &
                    SAVE_MESH_FILES,ANISOTROPY)
  use generate_databases_par, only: NGLLX,NGLLY,NGLLZ,IOUT, &
    USE_MESH_COLORING_GPU
  use create_regions_mesh_ext_par
  use constants, only: CUSTOM_REAL
  implicit none
  integer,intent(in) :: nspec,nglob
  ! ocean load
  logical,intent(in) :: APPROXIMATE_OCEAN_LOAD
  ! mesh coordinates
  integer, dimension(NGLLX,NGLLY,NGLLZ,nspec),intent(in) :: ibool
  ! MPI interfaces
  integer,intent(in) :: num_interfaces_ext_mesh
  integer, dimension(num_interfaces_ext_mesh),intent(in) :: my_neighbors_ext_mesh
  integer, dimension(num_interfaces_ext_mesh),intent(in) :: nibool_interfaces_ext_mesh
  integer,intent(in) :: max_interface_size_ext_mesh
  integer, dimension(NGLLX*NGLLX*max_interface_size_ext_mesh,num_interfaces_ext_mesh),intent(in) :: ibool_interfaces_ext_mesh
  logical,intent(in) :: SAVE_MESH_FILES
  logical,intent(in) :: ANISOTROPY
  end subroutine save_arrays_solver_ext_mesh_h5
