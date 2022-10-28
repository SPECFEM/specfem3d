
BIN=../bin

include ../config.h

program = $(BIN)/xcreate_inputs_files
program1 = $(BIN)/xcreate_all_inputs

objs =   param_reader.o create_inputs_files.o
objs1 =   param_reader.o create_all_input_files_for_Hybrid.o


# add another directory to the search path for files
vpath %.f90 ../shared

.SUFFIXES: .o .c .f90

$(program): $(objs)
	$(FC)  -o $@ $(objs) $(FFLAGS)

$(program1): $(objs1)
	$(FC)  -o $@ $(objs1) $(FFLAGS)

.f90.o:
	$(FC) -o $@ -c $< $(FFLAGS)

.c.o:
	$(CC) -o $@ -c $< $(CCLAGS)

.PHONY: clean

clean:
	rm -f $(program) $(program1) $(objs) *.lst *.o

all: $(program) $(program1)

