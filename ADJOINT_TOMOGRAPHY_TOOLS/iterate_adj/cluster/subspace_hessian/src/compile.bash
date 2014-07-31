mpif90 -o subspace_hessian -O3 subspace_hessian.f90 exit_mpi.f90 gll_library.f90
rm *.o

#mpif90 -o sem_fun_resample -O3 sem_fun_resample.f90 exit_mpi.f90 rthetaphi_xyz.f90

#ifort -o sem_fun_resample -CB sem_fun_resample.f90 exit_mpi.f90 rthetaphi_xyz.f90
