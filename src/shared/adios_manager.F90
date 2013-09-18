!==============================================================================
!> Tools to setup and cleanup ADIOS
!------------------------------------------------------------------------------
module adios_manager_mod

contains

!==============================================================================
!> Initialize ADIOS and setup the xml output file
subroutine adios_setup()
  use mpi
  use adios_write_mod, only: adios_init

  implicit none

  include 'constants.h'

  integer :: adios_err

  call adios_init_noxml (MPI_COMM_WORLD, adios_err);
  call adios_allocate_buffer (ADIOS_BUFFER_SIZE_IN_MB, adios_err)
end subroutine adios_setup

!==============================================================================
!> Finalize ADIOS. Must be called once everything is written down.
subroutine adios_cleanup()
  use mpi
  use adios_write_mod, only: adios_finalize

  implicit none
  integer :: myrank
  integer :: adios_err, ierr

  call MPI_Comm_rank(MPI_COMM_WORLD, myrank, ierr)
  call adios_finalize (myrank, adios_err)
end subroutine adios_cleanup

end module adios_manager_mod
