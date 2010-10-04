#=====================================================================
#
#               S p e c f e m 3 D  V e r s i o n  1 . 4
#               ---------------------------------------
#
#                 Dimitri Komatitsch and Jeroen Tromp
#    Seismological Laboratory - California Institute of Technology
#         (c) California Institute of Technology July 2005
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
#
#=====================================================================
#
# United States Government Sponsorship Acknowledged.
#

# decompose_mesh_SCOTCH/Makefile.  Generated from Makefile.in by configure.

FC = gfortran
FCFLAGS = #-g -O2
MPIFC = mpif90
MPILIBS = 
FLAGS_CHECK = $(FLAGS_NO_CHECK)
FLAGS_NO_CHECK = -std=gnu -fimplicit-none -frange-check -O3 -pedantic -pedantic-errors -Waliasing -Wampersand -Wline-truncation -Wsurprising -Wunderflow -fno-trapping-math
FCFLAGS_f90 = 
SCOTCH_LIBS = -L/opt/scotch/gcc64/lib -lscotch -lscotcherr

FCCOMPILE_CHECK = ${FC} ${FCFLAGS} $(FLAGS_CHECK)
FCCOMPILE_NO_CHECK = ${FC} ${FCFLAGS} $(FLAGS_NO_CHECK)
MPIFCCOMPILE_CHECK = ${MPIFC} ${FCFLAGS} $(FLAGS_CHECK)
MPIFCCOMPILE_NO_CHECK = ${MPIFC} ${FCFLAGS} $(FLAGS_NO_CHECK)
#FCLINK = $(MPIFCCOMPILE_NO_CHECK)
##FCLINK = $(FCCOMPILE_NO_CHECK)
FCLINK = $(FCCOMPILE_NO_CHECK)

CC = gcc
CFLAGS = -g -O2 $(CPPFLAGS)
CPPFLAGS =  $(COND_MPI_CPPFLAGS)
COND_MPI_CPPFLAGS = -DWITH_MPI
#COND_MPI_CPPFLAGS =

AR = ar
ARFLAGS = cru
RANLIB = ranlib

O = ../obj

# objects toggled between the parallel and serial version
#COND_MPI_OBJECTS = $O/parallel.o
#COND_MPI_OBJECTS = $O/serial.o

####
#### targets
####

# default targets for the Pyrized version
DEFAULT = xdecompose_mesh_SCOTCH


default: $(DEFAULT)

all: clean default

decompose_mesh_SCOTCH: xdecompose_mesh_SCOTCH


####
#### rules for executables
####


#DEC_OBJETS = part_decompose_mesh_SCOTCH.o \
                                decompose_mesh_SCOTCH.o \
                                program_decompose_mesh_SCOTCH.o

DEC_OBJECTS = part_decompose_mesh_SCOTCH.o \
                decompose_mesh_SCOTCH.o \
                program_decompose_mesh_SCOTCH.o


# rules for the pure Fortran version
xdecompose_mesh_SCOTCH:: $(DEC_OBJECTS) 
	${FCLINK} -o xdecompose_mesh_SCOTCH $(DEC_OBJECTS) $(SCOTCH_LIBS) 

# solver also depends on values from mesher
#xspecfem3D: $(XSPECFEM_OBJECTS) $(COND_MPI_OBJECTS)
#	${FCLINK} -o xspecfem3D $(XSPECFEM_OBJECTS) $(COND_MPI_OBJECTS) $(MPILIBS)
#



clean:
	rm -f *.o *.mod xdecompose_mesh_SCOTCH 
####
#### rule to build each .o file below
####

part_decompose_mesh_SCOTCH.o: part_decompose_mesh_SCOTCH.f90 
	${FCCOMPILE_CHECK} -c part_decompose_mesh_SCOTCH.f90

decompose_mesh_SCOTCH.o: decompose_mesh_SCOTCH.f90 part_decompose_mesh_SCOTCH.f90 part_decompose_mesh_SCOTCH.o
	${FCCOMPILE_CHECK} -c decompose_mesh_SCOTCH.f90

program_decompose_mesh_SCOTCH.o: program_decompose_mesh_SCOTCH.f90 part_decompose_mesh_SCOTCH.o decompose_mesh_SCOTCH.o
	${FCCOMPILE_CHECK} -c program_decompose_mesh_SCOTCH.f90


