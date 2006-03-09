
  program trampoline

! Perform Fortran mojo, and then run the Python script.

! With ifort v9 in particular, this function (i.e., MAIN__) will call
! the undocumented function __intel_new_proc_init or
! __intel_new_proc_init_P.  Without this, SPECFEM runs several
! times slower (!).

  call run_python_script

  end program trampoline
