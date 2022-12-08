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
S := ${S_TOP}/src/specfem3D
$(specfem3D_OBJECTS): S = ${S_TOP}/src/specfem3D

#######################################
# solver objects - no statically allocated arrays anymore

####
#### targets
####

# default targets for the pure Fortran version
specfem3D_TARGETS = \
	$E/xspecfem3D \
	$(EMPTY_MACRO)


specfem3D_OBJECTS = \
	$O/specfem3D_par.spec_module.o \
	$O/asdf_data.spec_module.o \
	$O/assemble_MPI_vector.spec.o \
	$O/check_stability.spec.o \
	$O/comp_source_time_function.spec.o \
	$O/compute_add_sources_acoustic.spec.o \
	$O/compute_add_sources_viscoelastic.spec.o \
	$O/compute_add_sources_poroelastic.spec.o \
	$O/compute_adj_source_frechet.spec.o \
	$O/compute_arrays_source.spec.o \
	$O/compute_boundary_kernel.spec.o \
	$O/compute_coupling_acoustic_el.spec.o \
	$O/compute_coupling_acoustic_po.spec.o \
	$O/compute_coupling_viscoelastic_ac.spec.o \
	$O/compute_coupling_viscoelastic_po.spec.o \
	$O/compute_coupling_poroelastic_ac.spec.o \
	$O/compute_coupling_poroelastic_el.spec.o \
	$O/compute_forces_acoustic_calling_routine.spec.o \
	$O/compute_forces_acoustic.spec.o \
	$O/compute_forces_viscoelastic_calling_routine.spec.o \
	$O/compute_forces_viscoelastic.spec.o \
	$O/compute_element_att_memory.spec.o \
	$O/compute_element_strain.spec.o \
	$O/compute_forces_poro_fluid_part.spec.o \
	$O/compute_forces_poroelastic_calling_routine.spec.o \
	$O/compute_forces_poro_solid_part.spec.o \
	$O/compute_gradient_in_acoustic.spec.o \
	$O/compute_interpolated_dva.spec.o \
	$O/compute_kernels.spec.o \
	$O/compute_seismograms.spec.o \
	$O/compute_stacey_acoustic.spec.o \
	$O/compute_stacey_viscoelastic.spec.o \
	$O/compute_stacey_poroelastic.spec.o \
	$O/compute_energy.spec.o \
	$O/convert_time.spec.o \
	$O/couple_with_injection.spec.o \
	$O/calendar.spec.o \
	$O/create_color_image.spec.o \
	$O/detect_mesh_surfaces.spec.o \
	$O/fault_solver_common.spec.o \
	$O/fault_solver_dynamic.spec.o \
	$O/fault_solver_kinematic.spec.o \
	$O/finalize_simulation.spec.o \
	$O/get_cmt.spec.o \
	$O/get_elevation.spec.o \
	$O/get_force.spec.o \
	$O/gravity_perturbation.spec.o \
	$O/initialize_simulation.spec.o \
	$O/iterate_time.spec.o \
	$O/iterate_time_undoatt.spec.o \
	$O/locate_MPI_slice.spec.o \
	$O/locate_point.spec.o \
	$O/locate_receivers.spec.o \
	$O/locate_source.spec.o \
	$O/make_gravity.spec.o \
	$O/noise_tomography.spec.o \
	$O/pml_allocate_arrays.spec.o \
	$O/pml_output_VTKs.spec.o \
	$O/pml_compute_accel_contribution.spec.o \
	$O/pml_compute_memory_variables.spec.o \
	$O/pml_par.spec_module.o \
	$O/prepare_attenuation.spec.o \
	$O/prepare_gpu.spec.o \
	$O/prepare_gravity.spec.o \
	$O/prepare_noise.spec.o \
	$O/prepare_optimized_arrays.spec.o \
	$O/prepare_timerun.spec.o \
	$O/prepare_wavefields.spec.o \
	$O/print_stf_file.spec.o \
	$O/read_external_stf.spec.o \
	$O/read_forward_arrays.spec.o \
	$O/read_mesh_databases.spec.o \
	$O/read_stations.spec.o \
	$O/save_adjoint_kernels.spec.o \
	$O/save_forward_arrays.spec.o \
	$O/setup_GLL_points.spec.o \
	$O/setup_movie_meshes.spec.o \
	$O/setup_sources_receivers.spec.o \
	$O/specfem3D.spec.o \
	$O/station_filter.spec.o \
	$O/surface_or_volume_integral.spec.o \
	$O/update_displacement_scheme.spec.o \
	$O/update_displacement_LDDRK.spec.o \
	$O/write_movie_output.spec.o \
	$O/write_output_ASCII_or_binary.spec.o \
	$O/write_output_SU.spec.o \
	$O/write_seismograms.spec.o \
	$(EMPTY_MACRO)

specfem3D_SHARED_OBJECTS = \
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
	$O/get_jacobian_boundaries.shared.o \
	$O/get_shape3D.shared.o \
	$O/gll_library.shared.o \
	$O/heap_sort.shared.o \
	$O/hex_nodes.shared.o \
	$O/init_openmp.shared.o \
	$O/lagrange_poly.shared.o \
	$O/netlib_specfun_erf.shared.o \
	$O/param_reader.cc.o \
	$O/prepare_assemble_MPI.shared.o \
	$O/read_topo_bathy_file.shared.o \
	$O/read_parameter_file.shared.o \
	$O/read_value_parameters.shared.o \
	$O/recompute_jacobian.shared.o \
	$O/save_header_file.shared.o \
	$O/search_kdtree.shared.o \
	$O/sort_array_coordinates.shared.o \
	$O/utm_geo.shared.o \
	$O/write_VTK_data.shared.o \
	$O/write_c_binary.cc.o \
	$(EMPTY_MACRO)


specfem3D_MODULES = \
	$(FC_MODDIR)/asdf_data.$(FC_MODEXT) \
	$(FC_MODDIR)/asdf_manager_mod.$(FC_MODEXT) \
	$(FC_MODDIR)/fault_solver_common.$(FC_MODEXT) \
	$(FC_MODDIR)/fault_solver_dynamic.$(FC_MODEXT) \
	$(FC_MODDIR)/fault_solver_kinematic.$(FC_MODEXT) \
	$(FC_MODDIR)/gravity_perturbation.$(FC_MODEXT) \
	$(FC_MODDIR)/image_pnm_par.$(FC_MODEXT) \
	$(FC_MODDIR)/manager_adios.$(FC_MODEXT) \
	$(FC_MODDIR)/pml_par.$(FC_MODEXT) \
	$(FC_MODDIR)/specfem_par.$(FC_MODEXT) \
	$(FC_MODDIR)/specfem_par_acoustic.$(FC_MODEXT) \
	$(FC_MODDIR)/specfem_par_elastic.$(FC_MODEXT) \
	$(FC_MODDIR)/specfem_par_poroelastic.$(FC_MODEXT) \
	$(FC_MODDIR)/specfem_par_movie.$(FC_MODEXT) \
	$(FC_MODDIR)/specfem_par_coupling.$(FC_MODEXT) \
	$(FC_MODDIR)/specfem_par_noise.$(FC_MODEXT) \
	$(FC_MODDIR)/user_noise_distribution.$(FC_MODEXT) \
	$(EMPTY_MACRO)


###
### MPI
###
specfem3D_SHARED_OBJECTS += $(COND_MPI_OBJECTS)

###
### GPU
###

specfem3D_SHARED_OBJECTS += $(gpu_OBJECTS)

###
### ADIOS
###

# using ADIOS files

adios_specfem3D_OBJECTS= \
	$O/read_mesh_databases_adios.spec_adios.o \
	$O/save_forward_arrays_adios.spec_adios.o \
	$O/read_forward_arrays_adios.spec_adios.o \
	$O/save_kernels_adios.spec_adios.o \
	$(EMPTY_MACRO)

adios_specfem3D_PREOBJECTS = \
	$O/adios_helpers_addons.shared_adios_cc.o \
	$O/adios_helpers_definitions.shared_adios.o \
	$O/adios_helpers_readers.shared_adios.o \
	$O/adios_helpers_writers.shared_adios.o \
	$O/adios_helpers.shared_adios.o \
	$(EMPTY_MACRO)

adios_specfem3D_STUBS = \
	$O/adios_method_stubs.cc.o \
	$(EMPTY_MACRO)

# conditional adios linking
ifeq ($(ADIOS),yes)
specfem3D_OBJECTS += $(adios_specfem3D_OBJECTS)
specfem3D_SHARED_OBJECTS += $(adios_specfem3D_PREOBJECTS)
else ifeq ($(ADIOS2),yes)
specfem3D_OBJECTS += $(adios_specfem3D_OBJECTS)
specfem3D_SHARED_OBJECTS += $(adios_specfem3D_PREOBJECTS)
else
specfem3D_SHARED_OBJECTS += $(adios_specfem3D_STUBS)
endif

###
### ASDF
###

asdf_specfem3D_OBJECTS = \
	$O/write_output_ASDF.spec.o \
	$O/read_adjoint_sources_ASDF.spec.o \
	$(EMPTY_MACRO)

asdf_specfem3D_SHARED_OBJECTS = $(asdf_shared_OBJECTS)
asdf_specfem3D_SHARED_STUBS = $(asdf_shared_STUBS)

# conditional asdf linking
ifeq ($(ASDF),yes)
SPECFEM_LINK_FLAGS += $(ASDF_LIBS)
specfem3D_OBJECTS += $(asdf_specfem3D_OBJECTS)
specfem3D_SHARED_OBJECTS += $(asdf_specfem3D_SHARED_OBJECTS)
else
specfem3D_OBJECTS += $(asdf_specfem3D_STUBS)
specfem3D_SHARED_OBJECTS += $(asdf_specfem3D_SHARED_STUBS)
endif
#

###
### VTK
###

ifeq ($(VTK),yes)
specfem3D_OBJECTS += \
	$O/vtk_window.spec.o \
	$O/vtk_helper.visualcc.o \
	$(EMPTY_MACRO)
specfem3D_MODULES += \
	$(FC_MODDIR)/vtk_window_par.$(FC_MODEXT) \
	$(EMPTY_MACRO)
endif

#######################################

####
#### rules for executables
####

ifeq ($(HAS_GPU),yes)
INFO_GPU_SPECFEM="building xspecfem3D $(BUILD_VERSION_TXT)"

${E}/xspecfem3D: $(specfem3D_OBJECTS) $(specfem3D_SHARED_OBJECTS)
	@echo ""
	@echo $(INFO_GPU_SPECFEM)
	@echo ""
	${FCLINK} -o ${E}/xspecfem3D $(specfem3D_OBJECTS) $(specfem3D_SHARED_OBJECTS) $(MPILIBS) $(GPU_LINK) $(VTKLIBS) $(SPECFEM_LINK_FLAGS)
	@echo ""

else

## non-cuda version
${E}/xspecfem3D: $(specfem3D_OBJECTS) $(specfem3D_SHARED_OBJECTS)
	@echo ""
	@echo "building xspecfem3D"
	@echo ""
	${FCLINK} -o ${E}/xspecfem3D $(specfem3D_OBJECTS) $(specfem3D_SHARED_OBJECTS) $(MPILIBS) $(VTKLIBS) $(SPECFEM_LINK_FLAGS)
	@echo ""

endif



#######################################

###
### Module dependencies
###

# Version file
$O/initialize_simulation.spec.o: ${SETUP}/version.fh

## fault
$O/fault_solver_dynamic.spec.o: $O/fault_solver_common.spec.o
$O/fault_solver_kinematic.spec.o: $O/fault_solver_common.spec.o
$O/compute_forces_viscoelastic.spec.o: $O/fault_solver_dynamic.spec.o
$O/compute_forces_viscoelastic_calling_routine.spec.o: $O/fault_solver_dynamic.spec.o $O/fault_solver_kinematic.spec.o

$O/prepare_timerun.spec.o: $O/fault_solver_dynamic.spec.o $O/fault_solver_kinematic.spec.o
$O/prepare_gpu.spec.o: $O/fault_solver_dynamic.spec.o $O/fault_solver_kinematic.spec.o

## gravity
$O/finalize_simulation.spec.o: $O/gravity_perturbation.spec.o
$O/iterate_time.spec.o: $O/gravity_perturbation.spec.o
$O/iterate_time_undoatt.spec.o: $O/gravity_perturbation.spec.o
$O/prepare_gravity.spec.o: $O/gravity_perturbation.spec.o

## ADIOS
$O/finalize_simulation.spec.o: $O/adios_manager.shared_adios_module.o
$O/initialize_simulation.spec.o: $O/adios_manager.shared_adios_module.o

## ASDF compilation
$O/write_output_ASDF.spec.o: $O/asdf_data.spec_module.o

## kdtree
$O/locate_point.spec.o: $O/search_kdtree.shared.o
$O/setup_sources_receivers.spec.o: $O/search_kdtree.shared.o

####
#### rule to build each .o file below
####

## modules
$O/%.spec_module.o: $S/%.f90 $O/shared_par.shared_module.o
	${FCCOMPILE_CHECK} ${FCFLAGS_f90} -c -o $@ $<

$O/%.spec_module.o: $S/%.F90 $O/shared_par.shared_module.o
	${FCCOMPILE_CHECK} ${FCFLAGS_f90} -c -o $@ $<


$O/%.spec.o: $S/%.f90 $O/specfem3D_par.spec_module.o $O/pml_par.spec_module.o
	${FCCOMPILE_CHECK} ${FCFLAGS_f90} -c -o $@ $<

$O/%.spec.o: $S/%.F90 $O/specfem3D_par.spec_module.o $O/pml_par.spec_module.o
	${FCCOMPILE_CHECK} ${FCFLAGS_f90} -c -o $@ $<

###
### ADIOS compilation
###

$O/%.spec_adios.o: $S/%.F90 $O/specfem3D_par.spec_module.o $O/pml_par.spec_module.o $O/adios_helpers.shared_adios.o
	${FCCOMPILE_CHECK} ${FCFLAGS_f90} -c -o $@ $<

$O/%.spec_adios.o: $S/%.f90 $O/specfem3D_par.spec_module.o $O/pml_par.spec_module.o $O/adios_helpers.shared_adios.o
	${FCCOMPILE_CHECK} ${FCFLAGS_f90} -c -o $@ $<

$O/%.spec_noadios.o: $S/%.F90 $O/specfem3D_par.spec_module.o $O/pml_par.spec_module.o
	${FCCOMPILE_CHECK} ${FCFLAGS_f90} -c -o $@ $<

$O/%.spec_noadios.o: $S/%.f90 $O/specfem3D_par.spec_module.o $O/pml_par.spec_module.o
	${FCCOMPILE_CHECK} ${FCFLAGS_f90} -c -o $@ $<


###
### OpenMP compilation
###

$O/%.openmp.o: $S/%.f90 ${SETUP}/constants.h
	${FCCOMPILE_CHECK} ${FCFLAGS_f90} -c -o $@ $<


####
#### VTK file
####

$O/%.visualcc.o: $S/%.cpp ${SETUP}/config.h
	${CC} -c $(CPPFLAGS) $(CFLAGS) $(MPI_INCLUDES) -o $@ $<

$O/%.visualcc.o: $S/%.c ${SETUP}/config.h
	${CC} -c $(CPPFLAGS) $(CFLAGS) $(MPI_INCLUDES) -o $@ $<

