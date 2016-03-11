
# Fortran90 compiler
F90 = gfortran
#F90 = ifort

FFLAGS = -O3

.NOTPARALLEL:

default: xconvert_external_layers_of_a_given_mesh_to_CPML_layers xadd_CPML_layers_to_an_existing_mesh xconvert_CUBIT_SPECFEM_mesh_to_OpenDX xshuffle_existing_mesh_for_tests xconvert_mesh_files_from_ASCII_to_binary xconvert_mesh_files_from_binary_to_ASCII cleanobj

all: clean default

xconvert_external_layers_of_a_given_mesh_to_CPML_layers: convert_external_layers_of_a_given_mesh_to_CPML_layers.o
	 ${F90} $(FFLAGS) -o xconvert_external_layers_of_a_given_mesh_to_CPML_layers convert_external_layers_of_a_given_mesh_to_CPML_layers.o

xadd_CPML_layers_to_an_existing_mesh: add_CPML_layers_to_an_existing_mesh.o
	 ${F90} $(FFLAGS) -o xadd_CPML_layers_to_an_existing_mesh add_CPML_layers_to_an_existing_mesh.o

xconvert_CUBIT_SPECFEM_mesh_to_OpenDX: convert_CUBIT_SPECFEM_mesh_to_OpenDX.o
	 ${F90} $(FFLAGS) -o xconvert_CUBIT_SPECFEM_mesh_to_OpenDX convert_CUBIT_SPECFEM_mesh_to_OpenDX.o

xshuffle_existing_mesh_for_tests: shuffle_existing_mesh_for_tests.o
	 ${F90} $(FFLAGS) -o xshuffle_existing_mesh_for_tests shuffle_existing_mesh_for_tests.o

xconvert_mesh_files_from_ASCII_to_binary: convert_mesh_files_from_ASCII_to_binary.o
	 ${F90} $(FFLAGS) -o xconvert_mesh_files_from_ASCII_to_binary convert_mesh_files_from_ASCII_to_binary.o

xconvert_mesh_files_from_binary_to_ASCII: convert_mesh_files_from_binary_to_ASCII.o
	 ${F90} $(FFLAGS) -o xconvert_mesh_files_from_binary_to_ASCII convert_mesh_files_from_binary_to_ASCII.o

clean:
	rm -rf *.o xconvert_external_layers_of_a_given_mesh_to_CPML_layers xadd_CPML_layers_to_an_existing_mesh xconvert_CUBIT_SPECFEM_mesh_to_OpenDX xshuffle_existing_mesh_for_tests xconvert_mesh_files_from_ASCII_to_binary xconvert_mesh_files_from_binary_to_ASCII

cleanobj:
	rm -rf *.o

# OBJECTS:  

convert_external_layers_of_a_given_mesh_to_CPML_layers.o: convert_external_layers_of_a_given_mesh_to_CPML_layers.f90 
	${F90} $(FFLAGS) -c -o convert_external_layers_of_a_given_mesh_to_CPML_layers.o convert_external_layers_of_a_given_mesh_to_CPML_layers.f90

add_CPML_layers_to_an_existing_mesh.o: add_CPML_layers_to_an_existing_mesh.f90 
	${F90} $(FFLAGS) -c -o add_CPML_layers_to_an_existing_mesh.o add_CPML_layers_to_an_existing_mesh.f90

convert_CUBIT_SPECFEM_mesh_to_OpenDX.o: convert_CUBIT_SPECFEM_mesh_to_OpenDX.f90 
	${F90} $(FFLAGS) -c -o convert_CUBIT_SPECFEM_mesh_to_OpenDX.o convert_CUBIT_SPECFEM_mesh_to_OpenDX.f90

shuffle_existing_mesh_for_tests.o: shuffle_existing_mesh_for_tests.f90 
	${F90} $(FFLAGS) -c -o shuffle_existing_mesh_for_tests.o shuffle_existing_mesh_for_tests.f90

convert_mesh_files_from_ASCII_to_binary.o: convert_mesh_files_from_ASCII_to_binary.f90 
	${F90} $(FFLAGS) -c -o convert_mesh_files_from_ASCII_to_binary.o convert_mesh_files_from_ASCII_to_binary.f90

convert_mesh_files_from_binary_to_ASCII.o: convert_mesh_files_from_binary_to_ASCII.f90 
	${F90} $(FFLAGS) -c -o convert_mesh_files_from_binary_to_ASCII.o convert_mesh_files_from_binary_to_ASCII.f90

