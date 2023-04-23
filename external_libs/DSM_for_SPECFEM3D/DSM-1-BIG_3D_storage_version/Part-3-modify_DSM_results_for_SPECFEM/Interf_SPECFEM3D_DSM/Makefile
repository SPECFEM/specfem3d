BIN=../../bin

include ../../config.h
program = $(BIN)/xread_absorbing_interfaces
objs = read_absorbing_interfaces.o

.SUFFIXES: .o .f90

$(program): $(objs)
	$(FC)  -o $@ $(objs) $(FFLAGS)
.f90.o:
	$(FC) -o $@ -c $< $(FFLAGS)

.PHONY: clean
clean:
	rm $(program) $(objs) *.lst *.o

all: $(program)
