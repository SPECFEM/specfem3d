#=====================================================================
#
#                         S p e c f e m 3 D
#                         -----------------
#
#     Main historical authors: Dimitri Komatitsch and Jeroen Tromp
#                              CNRS, France
#                       and Princeton University, USA
#                 (there are currently many more authors!)
#                           (c) October 2017
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 3 of the License, or
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

FC = gfortran
FFLAGS = -g -std=gnu -fimplicit-none -frange-check -O3 -pedantic -pedantic-errors -Waliasing -Wampersand -Wline-truncation -Wsurprising -Wunderflow -Wunused -fno-trapping-math
FC_COMPILE = ${FC} ${FFLAGS}


CC = gcc
CFLAGS = -g -O3 -Werror -Wall -pedantic-errors
LDFLAGS = -lm
CC_COMPILE = ${CC} ${CFLAGS}
BIN = ../bin

# default targets for the pure Fortran version
DEFAULT = visualize

default: $(DEFAULT)

all: clean default


visualize: ${BIN}/visualize_par.o ${BIN}/string_process.o ${BIN}/visualize.o ${BIN}/read_input.o ${BIN}/visualize_collect.o ${BIN}/write_ensight.o ${BIN}/write_vtu.o ${BIN}/cfunc4fortran.o
	${FC_COMPILE} -g -o ${BIN}/visualize  ${BIN}/visualize_par.o ${BIN}/string_process.o ${BIN}/visualize.o ${BIN}/read_input.o ${BIN}/visualize_collect.o ${BIN}/write_ensight.o ${BIN}/write_vtu.o ${BIN}/cfunc4fortran.o

clean:
	rm -f ${BIN}/*.o ${BIN}/*.exe ${BIN}/visualize

${BIN}/visualize_par.o: visualize_par.f90
	${FC_COMPILE} -c -o ${BIN}/visualize_par.o visualize_par.f90

${BIN}/string_process.o: string_process.f90
	${FC_COMPILE} -c -o ${BIN}/string_process.o string_process.f90

${BIN}/visualize.o: visualize.f90
	${FC_COMPILE} -c -o ${BIN}/visualize.o visualize.f90

${BIN}/read_input.o: read_input.f90
	${FC_COMPILE} -c -o ${BIN}/read_input.o read_input.f90


${BIN}/visualize_collect.o: visualize_collect.f90
	${FC_COMPILE} -c -o ${BIN}/visualize_collect.o visualize_collect.f90

${BIN}/write_ensight.o: write_ensight.f90
	${FC_COMPILE} -c -o ${BIN}/write_ensight.o write_ensight.f90

${BIN}/write_vtu.o: write_vtu.f90
	${FC_COMPILE} -c -o ${BIN}/write_vtu.o write_vtu.f90

${BIN}/cfunc4fortran.o: cfunc4fortran.c
	${CC_COMPILE} -c -o ${BIN}/cfunc4fortran.o cfunc4fortran.c

