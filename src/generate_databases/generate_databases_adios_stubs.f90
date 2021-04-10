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
!> \file generate_databases_adios_stubs.f90
!!
!!  Stubs for ADIOS functions. Avoid link error when not configured with
!!  ADIOS.
!!
!! \author MPBL
!==============================================================================

!--------------------------------------.
! subroutines from model_gll_adios.F90 |
!--------------------------------------'

  subroutine model_gll_adios(myrank,nspec,LOCAL_PATH)

  use constants, only: MAX_STRING_LEN
  use adios_manager_mod

  implicit none

  integer, intent(in) :: myrank,nspec
  character(len=MAX_STRING_LEN) :: LOCAL_PATH

  integer(kind=4) :: unused_i4

  unused_i4 = myrank
  unused_i4 = nspec
  unused_i4 = len_trim(LOCAL_PATH)

  call no_adios_err()

  end subroutine

!----------------------------------------.
! subroutines from model_ipati_adios.F90 |
!----------------------------------------'

module model_ipati_adios_mod

  implicit none

contains

  subroutine model_ipati_adios(myrank,nspec,LOCAL_PATH)

  use constants, only: MAX_STRING_LEN
  use adios_manager_mod

  implicit none

  integer, intent(in) :: myrank,nspec
  character(len=MAX_STRING_LEN), intent(in) :: LOCAL_PATH

  integer(kind=4) :: unused_i4

  unused_i4 = myrank
  unused_i4 = nspec
  unused_i4 = len_trim(LOCAL_PATH)

  call no_adios_err()

  end subroutine model_ipati_adios

!
!-----------
!

  subroutine model_ipati_water_adios(myrank,nspec,LOCAL_PATH)

  use constants, only: MAX_STRING_LEN
  use adios_manager_mod

  implicit none

  integer, intent(in) :: myrank,nspec
  character(len=MAX_STRING_LEN), intent(in) :: LOCAL_PATH

  integer(kind=4) :: unused_i4

  unused_i4 = myrank
  unused_i4 = nspec
  unused_i4 = len_trim(LOCAL_PATH)

  call no_adios_err()

  end subroutine model_ipati_water_adios

!
!-----------
!

  subroutine read_model_vp_rho_adios(myrank, nspec, LOCAL_PATH, rho_read, vp_read)

  use constants, only: MAX_STRING_LEN
  use adios_manager_mod

  implicit none

  integer, intent(in) :: myrank,nspec
  character(len=MAX_STRING_LEN), intent(in) :: LOCAL_PATH
  real, dimension(:,:,:,:), intent(inout) :: vp_read,rho_read

  integer(kind=4) :: unused_i4
  real :: unused_r

  unused_i4 = myrank
  unused_i4 = nspec
  unused_i4 = len_trim(LOCAL_PATH)
  unused_r = vp_read(1,1,1,1)
  unused_r = rho_read(1,1,1,1)

  call no_adios_err()

  end subroutine read_model_vp_rho_adios

end module model_ipati_adios_mod

!-------------------------------------------------.
! subroutines from read_partition_files_adios.F90 |
!-------------------------------------------------'

  subroutine read_partition_files_adios()

  use adios_manager_mod

  implicit none

  call no_adios_err()

  end subroutine read_partition_files_adios

!-----------------------------------------------.
! subroutines from save_arrays_solver_adios.F90 |
!-----------------------------------------------'

  subroutine save_arrays_solver_ext_mesh_adios(nspec, ibool)

  use adios_manager_mod

  implicit none
  integer :: nspec
  integer, dimension(1,1,1,nspec) :: ibool

  integer(kind=4) :: unused_i4

  unused_i4 = ibool(1,1,1,1)

  call no_adios_err()

  end subroutine

!--------------------------------------.
! subroutines from save_moho_adios.F90 |
!--------------------------------------'

  subroutine crm_save_moho_adios()

  use adios_manager_mod

  implicit none

  call no_adios_err()

  end subroutine crm_save_moho_adios
