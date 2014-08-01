FC=mpif90
FFLAGS=-O3 -assume byterecl

OBJS =  sum_kernels.o exit_mpi.o

all: sum_kernels

sum_kernels:  $(OBJS)
	$(FC) $(FFLAGS) -o $@ $(OBJS)

.SUFFIXES: .o .f90
.f90.o:
	$(FC) $(FFLAGS) -c $<


clean:
	rm -rf *.o *~ sum_kernels
