
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

