# includes default Makefile from previous configuration
include Makefile

# test target
default: test_valence

## compilation directories
O := ./obj

OBJECTS = \
	$O/decompose_mesh.dec.o \
	$O/fault_scotch.dec.o \
	$O/part_decompose_mesh.dec.o \
	$O/param_reader.cc.o \
	$O/shared_par.shared_module.o \
	$O/read_parameter_file.shared.o \
	$O/read_value_parameters.shared.o \
	$O/sort_array_coordinates.shared.o \
	$(EMPTY_MACRO)

test_valence:
	${FCCOMPILE_CHECK} ${FCFLAGS_f90} -o ./bin/test_valence test_valence.f90 -I./obj $(OBJECTS) $(SCOTCH_LIBS)

