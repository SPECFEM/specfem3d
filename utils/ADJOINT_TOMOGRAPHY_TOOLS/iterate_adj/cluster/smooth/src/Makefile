FC=mpif90
FFLAGS=-O3 -assume byterecl #-g -traceback

all: smooth_sem_fun

smooth_sem_fun: smooth_sem_fun.f90
	$(FC) -o smooth_sem_fun $(FFLAGS) smooth_sem_fun.f90 smooth_sub.f90 exit_mpi.f90 gll_library.f90


clean:
	rm -f smooth_sem_fun *.o

