
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

subroutine model_ipati_adios()
  use adios_manager_mod

  call no_adios_err()
end subroutine model_ipati_adios

subroutine model_ipati_water_adios()
  use adios_manager_mod

  call no_adios_err()
end subroutine model_ipati_water_adios

subroutine read_model_vp_rho_adios()
  use adios_manager_mod

  call no_adios_err()
end subroutine read_model_vp_rho_adios

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
