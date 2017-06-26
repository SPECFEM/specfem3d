
BIN=.

program = $(BIN)/xexpand_2D_3D
program1= $(BIN)/xreformat
program2= $(BIN)/create_stations.x

include ./config.h

objs =  global_parameters_mod.o rotation_matrix_mod.o mpi_mod.o reading_inputs_mod.o  interpol_mod.o post_process_mod.o  writing_mod.o reading_fields_mod.o gll_library.o  program_expand_2D_to_3D.o

objs1 = reformat_output_files.o

objs2 = create_stations.o

# add another directory to the search path for files
#vpath %.f90 ../shared

.SUFFIXES: .o .f90

$(program): $(objs)
	$(FC)  -o $@ $(objs) $(FFLAGS)


$(program1) :  $(objs1)
	$(FC)  -o  $@ $(objs1) $(FFLAGS)


$(program2) :  $(objs2)
	$(FC)  -o  $@ $(objs2) $(FFLAGS)

.f90.o:

	$(FC) -o $@ -c $< $(FFLAGS)

.PHONY: clean

clean:

	rm -f $(program) $(objs) $(program1) $(objs1) $(program2) $(objs2) *.lst *.o *.mod

all: $(program) $(program1) $(program2)
