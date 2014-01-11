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
!> \file generate_databases_adios_stubs.f90
!!
!!  Stubs for ADIOS functions. Avoid link error when not configured with
!!  ADIOS.
!!
!! \author MPBL
!==============================================================================

!--------------------------------------.
! Subroutines from model_gll_adios.F90 |
!--------------------------------------'

subroutine model_gll_adios(myrank,nspec,LOCAL_PATH)

  use adios_manager_mod

  implicit none

  integer, intent(in) :: myrank,nspec
  character(len=256) :: LOCAL_PATH

  ! dummy to avoid compiler warning about unused variables
  integer :: dummy
  character(len=256) :: dummy_line

  dummy = myrank
  dummy = nspec
  dummy_line = LOCAL_PATH

  call no_adios_err()
end subroutine

!----------------------------------------.
! Subroutines from model_ipati_adios.F90 |
!----------------------------------------'

module model_ipati_adios_mod
contains

subroutine model_ipati_adios(myrank,nspec,LOCAL_PATH)
  use adios_manager_mod

  implicit none

  integer, intent(in) :: myrank,nspec
  character(len=256), intent(in) :: LOCAL_PATH

  ! dummy to avoid compiler warning about unused variables
  integer :: dummy
  character(len=256) :: dummy_line

  dummy = myrank
  dummy = nspec
  dummy_line = LOCAL_PATH

  call no_adios_err()
end subroutine model_ipati_adios

subroutine model_ipati_water_adios(myrank,nspec,LOCAL_PATH)
  use adios_manager_mod

  implicit none

  integer, intent(in) :: myrank,nspec
  character(len=256), intent(in) :: LOCAL_PATH

  ! dummy to avoid compiler warning about unused variables
  integer :: dummy
  character(len=256) :: dummy_line

  dummy = myrank
  dummy = nspec
  dummy_line = LOCAL_PATH

  call no_adios_err()
end subroutine model_ipati_water_adios

subroutine read_model_vp_rho_adios(myrank, nspec, LOCAL_PATH, &
                                   rho_read, vp_read)
  use adios_manager_mod

  implicit none

  integer, intent(in) :: myrank,nspec
  character(len=256), intent(in) :: LOCAL_PATH
  real, dimension(:,:,:,:), intent(inout) :: vp_read,rho_read

  ! dummy to avoid compiler warning about unused variables
  integer :: dummy
  real :: dummy_r
  character(len=256) :: dummy_line

  dummy = myrank
  dummy = nspec
  dummy_line = LOCAL_PATH
  dummy_r = vp_read(1,1,1,1)
  dummy_r = rho_read(1,1,1,1)

  call no_adios_err()
end subroutine read_model_vp_rho_adios

end module model_ipati_adios_mod

!-------------------------------------------------.
! Subroutines from read_partition_files_adios.F90 |
!-------------------------------------------------'

subroutine read_partition_files_adios()
  use adios_manager_mod

  call no_adios_err()
end subroutine read_partition_files_adios

!-----------------------------------------------.
! Subroutines from save_arrays_solver_adios.F90 |
!-----------------------------------------------'

subroutine save_arrays_solver_ext_mesh_adios()
  use adios_manager_mod

  call no_adios_err()
end subroutine

subroutine save_arrays_solver_files_adios()
  use adios_manager_mod

  call no_adios_err()
end subroutine save_arrays_solver_files_adios

!--------------------------------------.
! Subroutines from save_moho_adios.F90 |
!--------------------------------------'

subroutine crm_save_moho_adios()
  use adios_manager_mod

  call no_adios_err()
end subroutine crm_save_moho_adios
