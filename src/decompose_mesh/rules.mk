#=====================================================================
#
#               S p e c f e m 3 D  V e r s i o n  3 . 0
#               ---------------------------------------
#
#     Main historical authors: Dimitri Komatitsch and Jeroen Tromp
#                        Princeton University, USA
#                and CNRS / University of Marseille, France
#                 (there are currently many more authors!)
# (c) Princeton University and CNRS / University of Marseille, July 2012
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

## compilation directories
S := ${S_TOP}/src/decompose_mesh
$(decompose_mesh_OBJECTS): S = ${S_TOP}/src/decompose_mesh

#######################################

####
#### targets
####

# default targets for the Pyrized version
decompose_mesh_TARGETS = \
	$E/xscotch \
	$E/xdecompose_mesh \
	$(EMPTY_MACRO)

decompose_mesh_OBJECTS = \
	$O/part_decompose_mesh.dec.o \
	$O/fault_scotch.dec.o \
	$O/decompose_mesh.dec.o \
	$O/program_decompose_mesh.dec.o \
	$(EMPTY_MACRO)

decompose_mesh_MODULES = \
	$(FC_MODDIR)/decompose_mesh.$(FC_MODEXT) \
	$(FC_MODDIR)/fault_scotch.$(FC_MODEXT) \
	$(FC_MODDIR)/part_decompose_mesh.$(FC_MODEXT) \
	$(EMPTY_MACRO)

decompose_mesh_SHARED_OBJECTS = \
	$O/shared_par.shared_module.o \
	$O/param_reader.cc.o \
	$O/read_parameter_file.shared.o \
	$O/read_value_parameters.shared.o \
	$O/sort_array_coordinates.shared.o \
	$(EMPTY_MACRO)


#######################################

####
#### rules for executables
####

dec: $(decompose_mesh_TARGETS)

decompose_mesh: xdecompose_mesh
xdecompose_mesh: $E/xdecompose_mesh

scotch: xscotch
xscotch: $E/xscotch

${SCOTCH_DIR}/include/scotchf.h: xscotch

# rules for the pure Fortran version
$E/xdecompose_mesh: ${SCOTCH_DIR}/include/scotchf.h $(decompose_mesh_SHARED_OBJECTS) $(decompose_mesh_OBJECTS) $(COND_MPI_OBJECTS)
	${FCLINK} -o  $E/xdecompose_mesh $(decompose_mesh_SHARED_OBJECTS) $(decompose_mesh_OBJECTS) $(SCOTCH_LIBS) $(COND_MPI_OBJECTS)

# scotch
$E/xscotch:
ifeq (${SCOTCH_BUNDLED},1)
	@echo "Using bundled Scotch in directory: ${SCOTCH_DIR}/src"
	$(MAKE) -C ${SCOTCH_DIR}/src
else
	@echo "Not using bundled Scotch"
endif



#######################################

###
### Module dependencies
###

$O/decompose_mesh.dec.o: $O/part_decompose_mesh.dec.o $O/fault_scotch.dec.o ${SCOTCH_DIR}/include/scotchf.h $O/shared_par.shared_module.o $(COND_MPI_OBJECTS)
$O/program_decompose_mesh.dec.o: $O/decompose_mesh.dec.o $O/shared_par.shared_module.o $(COND_MPI_OBJECTS)


#######################################

####
#### rule to build each .o file below
####

$O/%.dec.o: $S/%.f90 $O/shared_par.shared_module.o
	${FCCOMPILE_CHECK} ${FCFLAGS_f90} $(SCOTCH_INC) -c -o $@ $<

$O/%.dec.o: $S/%.F90 $O/shared_par.shared_module.o
	${FCCOMPILE_CHECK} ${FCFLAGS_f90} $(SCOTCH_INC) -c -o $@ $<



