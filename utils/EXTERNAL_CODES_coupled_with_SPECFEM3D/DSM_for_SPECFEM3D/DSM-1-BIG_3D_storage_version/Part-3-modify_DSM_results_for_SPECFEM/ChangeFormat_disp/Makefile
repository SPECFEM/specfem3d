BIN=../../bin

include ../../config.h

program = $(BIN)/ChangeFormat_disp

objs = ChangeFormat.o  filters.o read_input.o

# add another directory to the search path for files
vpath %.f90 ../../shared

.SUFFIXES: .o .f90

$(program): $(objs)
	$(FC)  -o $@ $(objs) $(FFLAGS)
.f90.o:
	$(FC) -o $@ -c $< $(FFLAGS)

.PHONY: clean
clean:
	rm $(program) $(objs) *.lst *.o

all: $(program)
