FFLAGS=-O3 -assume byterecl

all: add_model

add_model: add_model.f90
	mpif90 -o add_model $(FFLAGS) add_model.f90 gll_library.f90 exit_mpi.f90


clean:
	rm -f add_model *.o *.mod


