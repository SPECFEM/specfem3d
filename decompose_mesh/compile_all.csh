#!/bin/sh

#. /opt/intel/fce/10.0.026/bin/ifortvars.sh
# export FCFLAGS="-g -traceback -implicitnone -warn stderrors -warn truncated_source -warn argument_checking -warn unused -warn declarations -std95 -assume byterecl -C"

rm *.o *.mod ./a.out

gfortran -c part_decompose_mesh.f90
gfortran -c decompose_mesh.f90
#gfortran decompose_mesh.o part_decompose_mesh.o ~/utils/metis-4.0/libmetis.a
#gfortran decompose_mesh.o part_decompose_mesh.o ~/utils/scotch_5.1/lib/libscotchmetis.a ~/utils/scotch_5.1/lib/libscotch.a ~/utils/scotch_5.1/lib/libscotcherr.a
gfortran decompose_mesh.o part_decompose_mesh.o ~/utils/scotch_5.1/lib/libscotch.a ~/utils/scotch_5.1/lib/libscotcherr.a

