
  program xmeshfem3D

  implicit none

#ifdef USE_MPI
! standard include of the MPI library
  include 'mpif.h'
#endif

#ifdef USE_MPI
  integer ier
#endif

#ifdef USE_MPI
! initialize the MPI communicator and start the NPROCTOT MPI processes.
  call MPI_INIT(ier)
#endif

! run the main program
  call meshfem3D

#ifdef USE_MPI
! stop all the MPI processes, and exit
  call MPI_FINALIZE(ier)
#endif

  end program xmeshfem3D
