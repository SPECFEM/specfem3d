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

## compilation directories
S := ${S_TOP}/src/shared
$(shared_OBJECTS): S = ${S_TOP}/src/shared

#######################################

shared_TARGETS = \
	$(shared_OBJECTS) \
	$(EMPTY_MACRO)

shared_OBJECTS = \
	$O/shared_par.shared_module.o \
	$O/adios_manager.shared_adios_module.o \
	$O/assemble_MPI_scalar.shared.o \
	$O/check_mesh_resolution.shared.o \
	$O/count_number_of_sources.shared.o \
	$O/create_name_database.shared.o \
	$O/define_derivation_matrices.shared.o \
	$O/detect_surface.shared.o \
	$O/exit_mpi.shared.o \
	$O/force_ftz.cc.o \
	$O/get_attenuation_model.shared.o \
	$O/get_element_face.shared.o \
	$O/get_global.shared.o \
	$O/get_jacobian_boundaries.shared.o \
	$O/get_shape2D.shared.o \
	$O/get_shape3D.shared.o \
	$O/gll_library.shared.o \
	$O/hex_nodes.shared.o \
	$O/heap_sort.shared.o \
	$O/init_openmp.shared.o \
	$O/lagrange_poly.shared.o \
	$O/merge_sort.shared.o \
	$O/netlib_specfun_erf.shared.o \
	$O/param_reader.cc.o \
	$O/prepare_assemble_MPI.shared.o \
	$O/read_parameter_file.shared.o \
	$O/read_topo_bathy_file.shared.o \
	$O/read_value_parameters.shared.o \
	$O/recompute_jacobian.shared.o \
	$O/safe_alloc_mod.shared.o \
	$O/save_header_file.shared.o \
	$O/search_kdtree.shared.o \
	$O/sort_array_coordinates.shared.o \
	$O/utm_geo.shared.o \
	$O/write_c_binary.cc.o \
	$O/write_VTK_data.shared.o \
	$(EMPTY_MACRO)


shared_MODULES = \
	$(FC_MODDIR)/constants.$(FC_MODEXT) \
	$(FC_MODDIR)/manager_adios.$(FC_MODEXT) \
	$(FC_MODDIR)/kdtree_search.$(FC_MODEXT) \
	$(FC_MODDIR)/safe_alloc_mod.$(FC_MODEXT) \
	$(FC_MODDIR)/shared_input_parameters.$(FC_MODEXT) \
	$(FC_MODDIR)/shared_compute_parameters.$(FC_MODEXT) \
	$(FC_MODDIR)/shared_parameters.$(FC_MODEXT) \
	$(FC_MODDIR)/attenuation_model.$(FC_MODEXT) \
	$(EMPTY_MACRO)


###
### MPI
###
shared_OBJECTS += $(COND_MPI_OBJECTS)
ifeq ($(MPI),yes)
shared_MODULES += $(FC_MODDIR)/my_mpi.$(FC_MODEXT)
endif

###
### ADIOS
###

adios_shared_OBJECTS = \
	$O/adios_helpers_addons.shared_adios_cc.o \
	$O/adios_helpers_definitions.shared_adios.o \
	$O/adios_helpers_readers.shared_adios.o \
	$O/adios_helpers_writers.shared_adios.o \
	$O/adios_helpers.shared_adios.o \
	$(EMPTY_MACRO)

adios_shared_MODULES = \
	$(FC_MODDIR)/adios_helpers_definitions_mod.$(FC_MODEXT) \
	$(FC_MODDIR)/adios_helpers_mod.$(FC_MODEXT) \
	$(FC_MODDIR)/adios_helpers_readers_mod.$(FC_MODEXT) \
	$(FC_MODDIR)/adios_helpers_writers_mod.$(FC_MODEXT) \
	$(EMPTY_MACRO)

adios_shared_STUBS = \
	$O/adios_method_stubs.cc.o \
	$(EMPTY_MACRO)

ifeq ($(ADIOS),yes)
shared_OBJECTS += $(adios_shared_OBJECTS)
shared_MODULES += $(adios_shared_MODULES)
else ifeq ($(ADIOS2),yes)
shared_OBJECTS += $(adios_shared_OBJECTS)
shared_MODULES += $(adios_shared_MODULES)
else
shared_OBJECTS += $(adios_shared_STUBS)
endif

## dependencies
ifeq ($(ADIOS),yes)
$O/adios_helpers.shared_adios.o: $O/adios_helpers_readers.shared_adios.o $O/adios_helpers_writers.shared_adios.o $O/adios_helpers_definitions.shared_adios.o
else ifeq ($(ADIOS2),yes)
$O/adios_helpers.shared_adios.o: $O/adios_helpers_readers.shared_adios.o $O/adios_helpers_writers.shared_adios.o $O/adios_helpers_definitions.shared_adios.o
endif

###
### ASDF
###

asdf_shared_OBJECTS = \
	$O/asdf_manager.shared_asdf.o \
	$(EMPTY_MACRO)

asdf_shared_STUBS = \
	$O/asdf_method_stubs.cc.o \
	$O/asdf_manager_stubs.shared_asdf.o \
	$(EMPTY_MACRO)

ifeq ($(ASDF),yes)
shared_OBJECTS += $(asdf_shared_OBJECTS)
else
shared_OBJECTS += $(asdf_shared_STUBS)
endif

#######################################

####
#### rule for each .o file below
####

##
## shared
##

$O/%.shared_module.o: $S/%.f90 ${SETUP}/constants.h
	${FCCOMPILE_CHECK} ${FCFLAGS_f90} -c -o $@ $<

$O/%.shared_module.o: $S/%.F90 ${SETUP}/constants.h
	${FCCOMPILE_CHECK} ${FCFLAGS_f90} -c -o $@ $<

$O/%.shared.o: $S/%.f90 $O/shared_par.shared_module.o
	${FCCOMPILE_CHECK} ${FCFLAGS_f90} -c -o $@ $<

$O/%.shared.o: $S/%.F90 $O/shared_par.shared_module.o
	${FCCOMPILE_CHECK} ${FCFLAGS_f90} -c -o $@ $<

$O/%.sharedmpi.o: $S/%.f90 ${SETUP}/constants.h $O/shared_par.shared_module.o
	${MPIFCCOMPILE_CHECK} ${FCFLAGS_f90} -c -o $@ $<

$O/%.sharedmpi.o: $S/%.F90 ${SETUP}/constants.h $O/shared_par.shared_module.o
	${MPIFCCOMPILE_CHECK} ${FCFLAGS_f90} -c -o $@ $<


##
## adios
##

$O/%.shared_adios_module.o: $S/%.f90 ${SETUP}/constants.h $O/shared_par.shared_module.o
	${FCCOMPILE_CHECK} ${FCFLAGS_f90} -c -o $@ $<

$O/%.shared_adios_module.o: $S/%.F90 ${SETUP}/constants.h $O/shared_par.shared_module.o
	${FCCOMPILE_CHECK} ${FCFLAGS_f90} -c -o $@ $<

$O/%.shared_adios.o: $S/%.f90 $O/adios_manager.shared_adios_module.o $S/adios_helpers_definitions_1d_generic.inc $S/adios_helpers_readers_1d_generic.inc $S/adios_helpers_writers_1d_generic.inc $S/adios_helpers_writers_1d_generic_offset.inc
	${FCCOMPILE_CHECK} ${FCFLAGS_f90} -c -o $@ $<

$O/%.shared_adios.o: $S/%.F90 $O/adios_manager.shared_adios_module.o $S/adios_helpers_definitions_1d_generic.inc $S/adios_helpers_readers_1d_generic.inc $S/adios_helpers_writers_1d_generic.inc $S/adios_helpers_writers_1d_generic_offset.inc
	${FCCOMPILE_CHECK} ${FCFLAGS_f90} -c -o $@ $<

$O/%.shared_adios_cc.o: $S/%.c ${SETUP}/config.h
	${MPICC} -c $(CPPFLAGS) $(CFLAGS) $(ADIOS2_DEF) $(ADIOS_DEF) -o $@ $<

$O/%.shared_noadios.o: $S/%.f90
	${FCCOMPILE_CHECK} ${FCFLAGS_f90} -c -o $@ $<

## asdf

$O/%.shared_asdf.o: $S/%.f90
	${FCCOMPILE_CHECK} ${FCFLAGS_f90} -c -o $@ $<

##
## C compilation
##

$O/%.cc.o: $S/%.c ${SETUP}/config.h
	${CC} -c $(CPPFLAGS) $(CFLAGS) -o $@ $<

