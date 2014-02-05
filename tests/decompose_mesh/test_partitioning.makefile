# includes default Makefile from previous configuration
include Makefile

# test target
default: test_partitioning

OBJECTS = \
	./obj/decompose_mesh.dec.o \
	./obj/fault_scotch.dec.o \
	./obj/part_decompose_mesh.dec.o \
	./obj/get_value_parameters.shared.o \
	./obj/param_reader.cc.o \
	./obj/read_parameter_file.shared.o \
	./obj/read_value_parameters.shared.o \
	$(EMPTY_MACRO)

test_partitioning:
	${FCCOMPILE_CHECK} ${FCFLAGS_f90} -o ./bin/test_partitioning test_partitioning.f90 -I./obj $(OBJECTS) $(SCOTCH_LIBS)

