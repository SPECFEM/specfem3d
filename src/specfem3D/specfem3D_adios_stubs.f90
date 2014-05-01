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
!> \file specfem3D_adios_stubs.f90
!!
!!  Stubs for ADIOS functions. Avoid link error when not configured with
!!  ADIOS.
!!
!! \author MPBL
!==============================================================================

!------------------------------------------------.
! Subroutines from read_mesh_databases_adios.F90 |
!------------------------------------------------'

subroutine read_mesh_for_init()
  use adios_manager_mod

  call no_adios_err()
end subroutine read_mesh_for_init

subroutine read_mesh_databases_adios()
  use adios_manager_mod

  call no_adios_err()
end subroutine read_mesh_databases_adios

subroutine read_moho_mesh_adjoint_adios()
  use adios_manager_mod

  call no_adios_err()
end subroutine read_moho_mesh_adjoint_adios

!-----------------------------------------.
! Subroutines from save_kernels_adios.F90 |
!-----------------------------------------'

subroutine define_kernel_adios_variables()
  use adios_manager_mod

  call no_adios_err()
end subroutine define_kernel_adios_variables

subroutine perform_write_adios_kernels()
  use adios_manager_mod

  call no_adios_err()
end subroutine

subroutine save_kernels_acoustic_adios()
  use adios_manager_mod

  call no_adios_err()
end subroutine

subroutine save_kernels_elastic_adios()
  use adios_manager_mod

  call no_adios_err()
end subroutine

subroutine save_kernels_poroelastic_adios()
  use adios_manager_mod

  call no_adios_err()
end subroutine save_kernels_poroelastic_adios

subroutine save_kernels_hessian_adios()
  use adios_manager_mod

  call no_adios_err()
end subroutine save_kernels_hessian_adios

!------------------------------------------------.
! Subroutines from save_forward_arrays_adios.F90 |
!------------------------------------------------'

subroutine save_forward_arrays_adios()
  use adios_manager_mod

  call no_adios_err()
end subroutine save_forward_arrays_adios

!------------------------------------------------.
! Subroutines from read_forward_arrays_adios.F90 |
!------------------------------------------------'

subroutine read_forward_arrays_adios()
  use adios_manager_mod

  call no_adios_err()
end subroutine read_forward_arrays_adios

