!==============================================================================
!> Tools to setup and cleanup ADIOS
!------------------------------------------------------------------------------
module adios_manager_mod

contains

subroutine no_adios_err()
  use mpi

  implicit none

  integer :: myrank, code, ier

  call MPI_Comm_rank(MPI_COMM_WORLD, myrank, ier)
  if (myrank == 0) then
    print *, "----------------------------------------------------"
    print *, "Not configure to be compile with ADIOS."
    print *, "Check your par_file and set ADIOS_ENABLED to .false."
    print *, "or reconfigure using --with-adios."
    print *, "----------------------------------------------------"
  endif
  call MPI_Abort(MPI_COMM_WORLD, code, ier)
end subroutine

!==============================================================================
!> Initialize ADIOS and setup the xml output file
subroutine adios_setup()
  call no_adios_err()
end subroutine adios_setup

!==============================================================================
!> Finalize ADIOS. Must be called once everything is written down.
subroutine adios_cleanup()
  call no_adios_err()
end subroutine adios_cleanup

end module adios_manager_mod
