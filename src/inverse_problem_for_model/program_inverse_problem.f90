 program inverse_problem

! MPI initialization
  call init_mpi()

! run the main program
  call inverse_problem_main()

! MPI finish
  call finalize_mpi()

end program inverse_problem
