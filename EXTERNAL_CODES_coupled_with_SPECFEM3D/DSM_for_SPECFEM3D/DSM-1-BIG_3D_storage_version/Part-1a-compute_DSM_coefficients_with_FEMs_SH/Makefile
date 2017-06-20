
BIN=../bin

include ../config.h

program = $(BIN)/xTraPSV_write_ceof_mpi_SH

objs =   trial.o  glu2.o dclib.o dcsymbdl.o dcsymbdl_vector.o others.o calmat.o TraSH-0.1.0.o

# add another directory to the search path for files
#vpath %.f90 ../shared

.SUFFIXES: .o .f90

$(program): $(objs)
	$(FC)  -o $@ $(objs) $(FFLAGS)

.f90.o:
	$(FC) -o $@ -c $< $(FFLAGS)

.PHONY: clean

clean:
	rm -f $(program) $(objs) *.lst *.o

all: $(program)

