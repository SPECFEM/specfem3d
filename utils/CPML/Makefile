
# Fortran90 compiler
F90 = gfortran
#F90 = ifort

FFLAGS = -O3

default: xconvert_external_layers_to_CPML cleanobj

all: clean default

xconvert_external_layers_to_CPML: convert_external_layers_of_a_given_mesh_to_CPML_layers.o
	 ${F90} $(FFLAGS) -o xconvert_external_layers_to_CPML convert_external_layers_of_a_given_mesh_to_CPML_layers.o

clean:
	rm -rf *.o xconvert_external_layers_to_CPML

cleanobj:
	rm -rf *.o

# OBJECTS:  

convert_external_layers_of_a_given_mesh_to_CPML_layers.o: convert_external_layers_of_a_given_mesh_to_CPML_layers.f90 
	${F90} $(FFLAGS) -c -o convert_external_layers_of_a_given_mesh_to_CPML_layers.o convert_external_layers_of_a_given_mesh_to_CPML_layers.f90

