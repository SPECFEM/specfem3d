
BIN=../../bin

include ../../config.h

program = $(BIN)/xTraPSV_MPI_read_ceof_vert_PSV

objs =   trial.o  others.o TraPSV-0.1.3.o

# add another directory to the search path for files
vpath %.f90 ../../shared

.SUFFIXES: .o .f90

$(program): $(objs)
	$(FC)  -o $@ $(objs) $(FFLAGS)

.f90.o:
	$(FC) -o $@ -c $< $(FFLAGS)

.PHONY: clean

clean:
	rm -f $(program) $(objs) *.lst *.o

all: $(program)

