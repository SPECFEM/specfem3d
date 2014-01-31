!==============================================================================
!> Tools to setup and cleanup ADIOS
!------------------------------------------------------------------------------
module adios_manager_mod

contains

!==============================================================================
!> Initialize ADIOS and setup the xml output file
subroutine adios_setup()

  use adios_write_mod, only: adios_init

  implicit none

  include 'constants.h'

  integer :: adios_err
  integer :: comm

  call world_get_comm(comm)

  call adios_init_noxml (comm, adios_err);
  call adios_allocate_buffer (ADIOS_BUFFER_SIZE_IN_MB, adios_err)

end subroutine adios_setup

!==============================================================================
!> Finalize ADIOS. Must be called once everything is written down.
subroutine adios_cleanup()

  use adios_write_mod, only: adios_finalize

  implicit none
  integer :: myrank
  integer :: adios_err

  call world_rank(myrank)
  call sync_all()

  call adios_finalize (myrank, adios_err)

end subroutine adios_cleanup

end module adios_manager_mod
