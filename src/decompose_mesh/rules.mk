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
	$(xdecompose_mesh_OBJECTS) \
	$(EMPTY_MACRO)

decompose_mesh_MODULES = \
	$(FC_MODDIR)/decompose_mesh_par.$(FC_MODEXT) \
	$(FC_MODDIR)/fault_scotch.$(FC_MODEXT) \
	$(FC_MODDIR)/part_decompose_mesh.$(FC_MODEXT) \
	$(FC_MODDIR)/module_qsort.$(FC_MODEXT) \
	$(FC_MODDIR)/module_database.$(FC_MODEXT) \
	$(FC_MODDIR)/module_mesh.$(FC_MODEXT) \
	$(FC_MODDIR)/module_partition.$(FC_MODEXT) \
	$(FC_MODDIR)/scotch_par.$(FC_MODEXT) \
	$(EMPTY_MACRO)

decompose_mesh_SHARED_OBJECTS = \
	$O/shared_par.shared_module.o \
	$O/param_reader.cc.o \
	$O/exit_mpi.shared.o \
	$O/read_parameter_file.shared.o \
	$O/read_value_parameters.shared.o \
	$O/sort_array_coordinates.shared.o \
	$O/write_VTK_data.shared.o \
	$(EMPTY_MACRO)


#######################################

####
#### rules for executables
####

dec: $(decompose_mesh_TARGETS)

decompose_mesh: xdecompose_mesh
xdecompose_mesh: $E/xdecompose_mesh

decompose_mesh_mpi: xdecompose_mesh_mpi
xdecompose_mesh_mpi: $E/xdecompose_mesh_mpi

scotch: xscotch
xscotch: $E/xscotch

${SCOTCH_DIR}/include/scotchf.h: xscotch

#######################################

####
#### rules for each program follow
####

#######################################

##
## scotch
##
$E/xscotch:
ifeq (${USE_BUNDLED_SCOTCH},1)
	@echo "Using bundled Scotch in directory: ${SCOTCH_DIR}/src"
	$(MAKE) -C ${SCOTCH_DIR}/src
else
	@echo "Not using bundled Scotch"
endif


##
## xdecompose_mesh
##
xdecompose_mesh_OBJECTS = \
	$O/check_valence.dec.o \
	$O/decompose_mesh.dec.o \
	$O/decompose_mesh_par.dec_module.o \
	$O/fault_scotch.dec_module.o \
	$O/lts_helper.dec.o \
	$O/lts_setup_elements.dec.o \
	$O/part_decompose_mesh.dec_module.o \
	$O/partition_scotch.dec.o \
	$O/partition_metis.dec.o \
	$O/partition_patoh.dec.o \
	$O/partition_rows.dec.o \
	$O/program_decompose_mesh.dec.o \
	$O/read_mesh_files.dec.o \
	$O/wrap_patoh.o \
	$O/wrap_metis.o \
	$O/write_mesh_databases.dec.o \
	$(EMPTY_MACRO)

# rules for the pure Fortran version
$E/xdecompose_mesh: ${SCOTCH_DIR}/include/scotchf.h $(decompose_mesh_SHARED_OBJECTS) $(xdecompose_mesh_OBJECTS) $(COND_MPI_OBJECTS)
	@echo ""
	@echo "building xdecompose_mesh"
	@echo ""
	${FCLINK} -o  $E/xdecompose_mesh $(decompose_mesh_SHARED_OBJECTS) $(xdecompose_mesh_OBJECTS) $(PART_LIBS) $(COND_MPI_OBJECTS) $(LDFLAGS)
	@echo ""


##
## xdecompose_mesh_mpi
##
xdecompose_mesh_mpi_OBJECTS = \
	$O/fault_scotch.dec_module.o \
	$O/module_qsort.dec.o \
	$O/module_database.dec.o \
	$O/module_mesh.dec.o \
	$O/module_partition.dec.o \
	$O/program_decompose_mesh_mpi.mpidec.o \
	$(EMPTY_MACRO)

# conditional target, needs mpi
ifeq ($(MPI),yes)
decompose_mesh_TARGETS += \
	$E/xdecompose_mesh_mpi \
	$(EMPTY_MACRO)

decompose_mesh_OBJECTS += \
	$(xdecompose_mesh_mpi_OBJECTS) \
	$(EMPTY_MACRO)
endif

# parallel version of decompose_mesh
$E/xdecompose_mesh_mpi: $(decompose_mesh_SHARED_OBJECTS) $(xdecompose_mesh_mpi_OBJECTS) $(COND_MPI_OBJECTS)
	@echo ""
	@echo "building xdecompose_mesh_mpi"
	@echo ""
	${FCLINK} -o  $E/xdecompose_mesh_mpi $(decompose_mesh_SHARED_OBJECTS) $(xdecompose_mesh_mpi_OBJECTS) $(COND_MPI_OBJECTS) $(LDFLAGS)
	@echo ""


#######################################

###
### Module dependencies
###

# serial version
$O/program_decompose_mesh.dec.o: $(COND_MPI_OBJECTS)

$O/decompose_mesh.dec.o: $O/shared_par.shared_module.o $(COND_MPI_OBJECTS)
$O/decompose_mesh_par.dec_module.o: $O/part_decompose_mesh.dec_module.o $O/fault_scotch.dec_module.o

# mpi version
$O/program_decompose_mesh_mpi.mpidec.o: $O/shared_par.shared_module.o $O/module_mesh.dec.o $O/module_database.dec.o $O/module_partition.dec.o $(COND_MPI_OBJECTS)

$O/module_database.dec.o : $O/shared_par.shared_module.o
$O/module_partition.dec.o : $O/shared_par.shared_module.o $O/fault_scotch.dec_module.o $O/module_qsort.dec.o
$O/module_mesh.dec.o : $O/shared_par.shared_module.o $O/fault_scotch.dec_module.o

#######################################

####
#### rule to build each .o file below
####

$O/%.dec_module.o: $S/%.F90 $O/shared_par.shared_module.o
	${FCCOMPILE_CHECK} ${FCFLAGS_f90} $(PART_FLAGS) -c -o $@ $<

$O/%.dec_module.o: $S/%.f90 $O/shared_par.shared_module.o
	${FCCOMPILE_CHECK} ${FCFLAGS_f90} $(PART_FLAGS) -c -o $@ $<

$O/%.dec.o: $S/%.f90 $O/shared_par.shared_module.o $O/part_decompose_mesh.dec_module.o $O/decompose_mesh_par.dec_module.o
	${FCCOMPILE_CHECK} ${FCFLAGS_f90} $(PART_FLAGS) -c -o $@ $<

$O/%.dec.o: $S/%.F90 $O/shared_par.shared_module.o $O/part_decompose_mesh.dec_module.o  $O/decompose_mesh_par.dec_module.o
	${FCCOMPILE_CHECK} ${FCFLAGS_f90} $(PART_FLAGS) -c -o $@ $<

$O/%.mpidec.o: $S/%.f90 $O/shared_par.shared_module.o
	${MPIFCCOMPILE_CHECK} ${FCFLAGS_f90} -c -o $@ $<

$O/wrap_patoh.o: $S/wrap_patoh.cpp
	${CXX} -c $(CFLAGS) -o $@ $< $(PART_FLAGS)

$O/wrap_metis.o: $S/wrap_metis.c
	${CC} -c $(CFLAGS) -o $@ $< $(PART_FLAGS)

