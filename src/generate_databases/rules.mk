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
S := ${S_TOP}/src/generate_databases
$(generate_databases_OBJECTS): S = ${S_TOP}/src/generate_databases

#######################################

####
#### targets
####

generate_databases_TARGETS = \
	$E/xgenerate_databases \
	$(EMPTY_MACRO)


generate_databases_OBJECTS = \
	$O/generate_databases_par.gen_mod.o \
	$O/calc_jacobian.gen.o \
	$O/fault_generate_databases.gen.o \
	$O/create_mass_matrices.gen.o \
	$O/create_regions_mesh.gen.o \
	$O/finalize_databases.gen.o \
	$O/generate_databases.gen.o \
	$O/get_absorbing_boundary.gen.o \
	$O/get_coupling_surfaces.gen.o \
	$O/get_model.gen.o \
	$O/get_MPI.gen.o \
	$O/get_perm_color.gen.o \
	$O/model_1d_cascadia.gen.o \
	$O/model_1d_prem.gen.o \
	$O/model_1d_socal.gen.o \
	$O/model_aniso.gen.o \
	$O/model_coupled.gen.o \
	$O/model_default.gen.o \
	$O/model_external_values.gen.o \
	$O/model_ipati.gen.o \
	$O/parse_sep.genc.o \
	$O/model_gll.gen.o \
	$O/model_salton_trough.gen.o \
	$O/model_tomography.gen.o \
	$O/pml_set_local_dampingcoeff.gen.o \
	$O/read_parameters.gen.o \
	$O/read_partition_files.gen.o \
	$O/save_arrays_solver.gen.o \
	$O/setup_color_perm.gen.o \
	$O/setup_mesh.gen.o \
	$O/memory_eval.gen.o \
	$(EMPTY_MACRO)


generate_databases_MODULES = \
	$(FC_MODDIR)/create_regions_mesh_ext_par.$(FC_MODEXT) \
	$(FC_MODDIR)/external_model.$(FC_MODEXT) \
	$(FC_MODDIR)/fault_generate_databases.$(FC_MODEXT) \
	$(FC_MODDIR)/generate_databases_par.$(FC_MODEXT) \
	$(FC_MODDIR)/manager_adios.$(FC_MODEXT) \
	$(FC_MODDIR)/model_coupled_par.$(FC_MODEXT) \
	$(FC_MODDIR)/model_sep_mod.$(FC_MODEXT) \
	$(FC_MODDIR)/model_tomography_par.$(FC_MODEXT) \
	$(FC_MODDIR)/salton_trough_par.$(FC_MODEXT) \
	$(EMPTY_MACRO)


generate_databases_SHARED_OBJECTS = \
	$O/assemble_MPI_scalar.shared.o \
	$O/shared_par.shared_module.o \
	$O/adios_manager.shared_adios_module.o \
	$O/check_mesh_resolution.shared.o \
	$O/count_number_of_sources.shared.o \
	$O/create_name_database.shared.o \
	$O/define_derivation_matrices.shared.o \
	$O/detect_surface.shared.o \
	$O/exit_mpi.shared.o \
	$O/get_attenuation_model.shared.o \
	$O/get_element_face.shared.o \
	$O/get_global.shared.o \
	$O/get_jacobian_boundaries.shared.o \
	$O/get_shape2D.shared.o \
	$O/get_shape3D.shared.o \
	$O/gll_library.shared.o \
	$O/hex_nodes.shared.o \
	$O/lagrange_poly.shared.o \
	$O/netlib_specfun_erf.shared.o \
	$O/param_reader.cc.o \
	$O/prepare_assemble_MPI.shared.o \
	$O/read_topo_bathy_file.shared.o \
	$O/read_parameter_file.shared.o \
	$O/read_value_parameters.shared.o \
	$O/recompute_jacobian.shared.o \
	$O/save_header_file.shared.o \
	$O/sort_array_coordinates.shared.o \
	$O/utm_geo.shared.o \
	$O/write_VTK_data.shared.o \
	$(EMPTY_MACRO)

# MPI stuffs
ifeq ($(MPI),no)
mpi_generate_databases_OBJECTS = $O/model_sep_nompi.gen.o
else
mpi_generate_databases_OBJECTS = $O/model_sep.mpi_gen.o
endif
generate_databases_OBJECTS += $(mpi_generate_databases_OBJECTS)

# using ADIOS files

adios_generate_databases_PREOBJECTS = \
	$O/adios_helpers_addons.shared_adios_cc.o \
	$O/adios_helpers_definitions.shared_adios.o \
	$O/adios_helpers_readers.shared_adios.o \
	$O/adios_helpers_writers.shared_adios.o \
	$O/adios_helpers.shared_adios.o \
	$(EMPTY_MACRO)

adios_generate_databases_OBJECTS = \
	$O/read_partition_files_adios.gen_adios.o \
	$O/save_arrays_solver_adios.gen_adios.o \
	$O/save_moho_adios.gen_adios.o \
	$O/model_gll_adios.gen_adios.o \
	$O/model_ipati_adios.gen_adios.o \
	$(EMPTY_MACRO)

adios_generate_databases_STUBS = \
	$O/adios_method_stubs.cc.o \
	$(EMPTY_MACRO)

# conditional adios linking
ifeq ($(ADIOS),yes)
generate_databases_OBJECTS += $(adios_generate_databases_OBJECTS)
generate_databases_SHARED_OBJECTS += $(adios_generate_databases_PREOBJECTS)
else ifeq ($(ADIOS2),yes)
generate_databases_OBJECTS += $(adios_generate_databases_OBJECTS)
generate_databases_SHARED_OBJECTS += $(adios_generate_databases_PREOBJECTS)
else
generate_databases_SHARED_OBJECTS += $(adios_generate_databases_STUBS)
endif

# objects for the pure Fortran version
XGENERATE_DATABASES_OBJECTS = \
	$(generate_databases_OBJECTS) $(generate_databases_SHARED_OBJECTS) \
	$(COND_MPI_OBJECTS)


#######################################


####
#### rules for executables
####

gen: $(generate_databases_TARGETS)

generate_databases: xgenerate_databases
xgenerate_databases: $E/xgenerate_databases

# rules for the pure Fortran version
$E/xgenerate_databases: $(XGENERATE_DATABASES_OBJECTS)
	@echo ""
	@echo "building xgenerate_databases"
	@echo ""
	${FCLINK} -o ${E}/xgenerate_databases $(XGENERATE_DATABASES_OBJECTS) $(MPILIBS)
	@echo ""



#######################################

###
### Module dependencies
###

# Version file
$O/generate_databases.gen.o: ${SETUP}/version.fh

ifeq ($(MPI),no)
$O/get_model.gen.o: $O/model_sep_nompi.gen.o
else
$O/get_model.gen.o: $O/model_sep.mpi_gen.o
endif

ifeq ($(MPI),yes)
$O/model_sep.mpi_gen.o: $O/parallel.sharedmpi.o
endif

$O/create_regions_mesh.gen.o: $O/fault_generate_databases.gen.o

## adios
$O/generate_databases.gen.o: $O/adios_manager.shared_adios_module.o

#######################################

####
#### rule to build each .o file below
####

$O/%.gen_mod.o: $S/%.F90 $O/shared_par.shared_module.o
	${FCCOMPILE_CHECK} ${FCFLAGS_f90} -c -o $@ $<

$O/%.gen.o: $S/%.f90 $O/shared_par.shared_module.o $O/generate_databases_par.gen_mod.o
	${FCCOMPILE_CHECK} ${FCFLAGS_f90} -c -o $@ $<

$O/%.gen.o: $S/%.F90 $O/shared_par.shared_module.o $O/generate_databases_par.gen_mod.o
	${FCCOMPILE_CHECK} ${FCFLAGS_f90} -c -o $@ $<

$O/%.mpi_gen.o: $S/%.f90 $O/shared_par.shared_module.o $O/generate_databases_par.gen_mod.o
	${MPIFCCOMPILE_CHECK} ${FCFLAGS_f90} -c -o $@ $<

$O/%.mpi_gen.o: $S/%.F90 $O/shared_par.shared_module.o $O/generate_databases_par.gen_mod.o
	${MPIFCCOMPILE_CHECK} ${FCFLAGS_f90} -c -o $@ $<

$O/%.genc.o: $S/%.c
	${CC} ${CFLAGS} -c -o $@ $<

###
### ADIOS compilation
###

$O/%.gen_adios.o: $S/%.F90 $O/shared_par.shared_module.o $O/generate_databases_par.gen_mod.o $O/adios_helpers.shared_adios.o
	${FCCOMPILE_CHECK} ${FCFLAGS_f90} -c -o $@ $<

$O/%.gen_adios.o: $S/%.f90 $O/shared_par.shared_module.o $O/generate_databases_par.gen_mod.o $O/adios_helpers.shared_adios.o
	${FCCOMPILE_CHECK} ${FCFLAGS_f90} -c -o $@ $<

$O/%.gen_noadios.o: $S/%.F90
	${FCCOMPILE_CHECK} ${FCFLAGS_f90} -c -o $@ $<

$O/%.gen_noadios.o: $S/%.f90
	${FCCOMPILE_CHECK} ${FCFLAGS_f90} -c -o $@ $<

