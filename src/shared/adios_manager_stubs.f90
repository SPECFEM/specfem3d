!==============================================================================
!> Tools to setup and cleanup ADIOS
!------------------------------------------------------------------------------
module adios_manager_mod

contains

subroutine no_adios_err()
  implicit none

  integer :: myrank

  call world_rank(myrank)
  if (myrank == 0) then
    print *, "----------------------------------------------------"
    print *, "Not configure to be compile with ADIOS."
    print *, "Check your par_file and set ADIOS_ENABLED to .false."
    print *, "or reconfigure using --with-adios."
    print *, "----------------------------------------------------"
  endif
  call stop_all()
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
