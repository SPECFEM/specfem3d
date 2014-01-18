
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

subroutine model_gll_adios()
  use adios_manager_mod

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

  call no_adios_err()
end subroutine model_ipati_adios

subroutine model_ipati_water_adios(myrank,nspec,LOCAL_PATH)
  use adios_manager_mod

  implicit none

  integer, intent(in) :: myrank,nspec
  character(len=256), intent(in) :: LOCAL_PATH

  call no_adios_err()
end subroutine model_ipati_water_adios

subroutine read_model_vp_rho_adios(myrank, nspec, LOCAL_PATH, &
                                   rho_read, vp_read)
  use adios_manager_mod

  implicit none

  integer, intent(in) :: myrank,nspec
  character(len=256), intent(in) :: LOCAL_PATH
  real, dimension(:,:,:,:), intent(inout) :: vp_read,rho_read

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
