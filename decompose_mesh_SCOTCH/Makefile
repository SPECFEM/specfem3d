# Makefile

F90 = gfortran

SCOTCH_LIBS = /data2/tarje/SCOTCH/lib/libscotch.a /data2/tarje/SCOTCH/lib/libscotcherr.a
#SCOTCH_LIBS = /scratch/network/SCOTCH/lib/libscotch.a /scratch/network/SCOTCH/lib/libscotcherr.a

LIBS = part_decompose_mesh_SCOTCH.o \
				decompose_mesh_SCOTCH.o \
				program_decompose_mesh_SCOTCH.o

# targets
all: xdecompose_mesh_SCOTCH

xdecompose_mesh_SCOTCH: $(LIBS)
	${F90} -Wall -o xdecompose_mesh_SCOTCH $(LIBS) $(SCOTCH_LIBS)


part_decompose_mesh_SCOTCH.o: part_decompose_mesh_SCOTCH.f90
	${F90} -Wall -c part_decompose_mesh_SCOTCH.f90

decompose_mesh_SCOTCH.o: decompose_mesh_SCOTCH.f90 part_decompose_mesh_SCOTCH.f90
	${F90} -Wall -c decompose_mesh_SCOTCH.f90

program_decompose_mesh_SCOTCH.o: program_decompose_mesh_SCOTCH.f90
	${F90} -Wall -c program_decompose_mesh_SCOTCH.f90

clean:
	rm -f *.o *.mod a.out xdecompose_mesh_SCOTCH
