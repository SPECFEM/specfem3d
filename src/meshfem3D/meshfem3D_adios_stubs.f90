
!==============================================================================
!> \file meshfem3D_adios_stubs.f90
!!
!!  Stubs for ADIOS functions. Avoid link error when not configured with
!!  ADIOS.
!!
!! \author MPBL
!==============================================================================

subroutine save_databases_adios()
  use adios_manager_mod

  call no_adios_err()
end subroutine save_databases_adios
