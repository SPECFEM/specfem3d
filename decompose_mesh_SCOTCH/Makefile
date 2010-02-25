# Makefile

#############################################################
## modify to match your compiler defaults 
## (which were used to compile SCOTCH libraries from below as well)
F90 = ifort      # use -warn
#F90 = gfortran    # use -Wall

## modify to match your library paths
SCOTCH_LIBS = -L/opt/seismo-util/source/scotch/lib -lscotch -lscotcherr
#SCOTCH_LIBS = -L/scratch/network/SCOTCH/lib -lscotch -lscotcherr							
#SCOTCH_LIBS = -L/n/sw/scotch-5.1/lib -lscotch -lscotcherr

#############################################################

LIBS = part_decompose_mesh_SCOTCH.o \
				decompose_mesh_SCOTCH.o \
				program_decompose_mesh_SCOTCH.o

# targets
all: xdecompose_mesh_SCOTCH

xdecompose_mesh_SCOTCH: $(LIBS)
	${F90} -o xdecompose_mesh_SCOTCH $(LIBS) $(SCOTCH_LIBS)


part_decompose_mesh_SCOTCH.o: part_decompose_mesh_SCOTCH.f90
	${F90} -c part_decompose_mesh_SCOTCH.f90

decompose_mesh_SCOTCH.o: decompose_mesh_SCOTCH.f90 part_decompose_mesh_SCOTCH.f90
	${F90} -c decompose_mesh_SCOTCH.f90

program_decompose_mesh_SCOTCH.o: program_decompose_mesh_SCOTCH.f90
	${F90} -c program_decompose_mesh_SCOTCH.f90

clean:
	rm -f *.o *.mod a.out xdecompose_mesh_SCOTCH
