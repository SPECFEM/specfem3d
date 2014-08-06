
BIN=../bin

include ../config.h

program = $(BIN)/xmesh_chunk_vm

objs = ReadIasp91.o mesh_chunk.o gll_library.o

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

