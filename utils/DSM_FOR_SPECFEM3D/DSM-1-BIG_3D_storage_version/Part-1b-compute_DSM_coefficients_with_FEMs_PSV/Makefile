
BIN=../bin

include ../config.h

program = $(BIN)/xTraPSV_write_ceof_mpi_PSV

objs =   trial.o rk3.o glu2.o dcsymbdl.o others.o calmat.o TraPSV-0.1.3.o

# add another directory to the search path for files
vpath %.f90 ../shared

.SUFFIXES: .o .f90

$(program): $(objs)
	$(FC)  -o $@ $(objs) $(FFLAGS)

.f90.o:
	$(FC) -o $@ -c $< $(FFLAGS)

.PHONY: clean

clean:
	rm -f $(program) $(objs) *.lst *.o

all: $(program)

