icc -c param_reader.c
ifort -c create_inputs_files_for_benchmark.f90
ifort -check all param_reader.o create_inputs_files_for_benchmark.o -o ../bin/xcreate_inputs_files

