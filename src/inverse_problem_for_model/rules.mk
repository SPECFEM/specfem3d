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
S := ${S_TOP}/src/inverse_problem_for_model
$(inverse_problem_for_model_OBJECTS): S = ${S_TOP}/src/inverse_problem_for_model

#######################################
# solver objects - no statically allocated arrays anymore

####
#### targets
####

# default targets for the pure Fortran version
inverse_problem_for_model_TARGETS = \
	$E/xinverse_problem_for_model \
	$(EMPTY_MACRO)


## source folder objects
inv_top_OBJECTS = \
	$O/anisotropic_parametrisation_mod.inv.o \
	$O/elastic_tensor_tools_mod.inv.o \
	$O/inverse_problem_main.inv.o \
	$O/inverse_problem_par.inv_par.o \
	$O/parallel_for_inverse_problem.invmpi.o \
	$O/program_inverse_problem.inv.o \
	$(EMPTY_MACRO)

inv_adjoint_source_OBJECTS = \
	$O/adjoint_source_mod.inv_adjoint_source.o \
	$O/rotations_mod.inv_adjoint_source.o \
	$O/signal_processing_mod.inv_adjoint_source.o \
	$(EMPTY_MACRO)

inv_input_output_OBJECTS = \
	$O/input_output_mod.inv_input.o \
	$O/interpolation_mod.inv_input.o \
	$O/IO_model_mod.inv_input.o \
	$O/mesh_tools_mod.inv_input.o \
	$O/passive_imaging_format_mod.inv_input.o \
	$O/Teleseismic_IO.inv_inputmpi.o \
	$(EMPTY_MACRO)

inv_inversion_scheme_OBJECTS = \
	$O/family_parameter_mod.inv_inversion.o \
	$O/fwi_iteration_mod.inv_inversion.o \
	$O/inversion_scheme_mod.inv_inversion.o \
	$O/iso_parameters.inv_inversion.o \
	$O/PrecondFWI_mod.inv_inversion.o \
	$O/vti_parameters.inv_inversion.o \
	$(EMPTY_MACRO)

inv_parameterisation_OBJECTS = \
	$O/elastic_isotropic_mod.inv_parameterisation.o \
	$O/parameters_subfamily_mod.inv_parameterisation.o \
	$(EMPTY_MACRO)

inv_projection_on_FD_grid_OBJECTS = \
	$O/projection_on_FD_grid_mod.inv_projection.o \
	$(EMPTY_MACRO)

inv_regularization_OBJECTS = \
	$O/regularization_FD_mod.inv_regularization.o \
	$O/regularization_interface.inv_regularization.o \
	$O/regularization_SEM_mod.inv_regularization.o \
	$(EMPTY_MACRO)

inv_specfem_interface_OBJECTS = \
	$O/specfem_interface_mod.inv_specfem_interface.o \
	$(EMPTY_MACRO)

inverse_problem_for_model_OBJECTS = \
	$(inv_top_OBJECTS) \
	$(inv_adjoint_source_OBJECTS) \
	$(inv_input_output_OBJECTS) \
	$(inv_inversion_scheme_OBJECTS) \
	$(inv_parameterization_OBJECTS) \
	$(inv_projection_on_FD_grid_OBJECTS) \
	$(inv_regularization_OBJECTS) \
	$(inv_specfem_interface_OBJECTS) \
	$(EMPTY_MACRO)


## objects from other source directories
inverse_problem_for_model_OBJECTS += \
	$O/specfem3D_par.spec_module.o \
	$O/asdf_data.spec_module.o \
	$O/assemble_MPI_vector.spec.o \
	$O/calendar.spec.o \
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
	$O/station_filter.spec.o \
	$O/surface_or_volume_integral.spec.o \
	$O/update_displacement_scheme.spec.o \
	$O/update_displacement_LDDRK.spec.o \
	$O/write_movie_output.spec.o \
	$O/write_output_ASCII_or_binary.spec.o \
	$O/write_output_SU.spec.o \
	$O/write_seismograms.spec.o \
	$(EMPTY_MACRO)


inverse_problem_for_model_SHARED_OBJECTS = \
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


inverse_problem_for_model_MODULES = \
	$(FC_MODDIR)/user_noise_distribution.$(FC_MODEXT) \
	$(FC_MODDIR)/adjoint_source.$(FC_MODEXT) \
	$(FC_MODDIR)/input_output.$(FC_MODEXT) \
	$(FC_MODDIR)/mesh_tools.$(FC_MODEXT) \
	$(FC_MODDIR)/regularization_interface.$(FC_MODEXT) \
	$(FC_MODDIR)/anisotropic_parametrisation_mod.$(FC_MODEXT) \
	$(FC_MODDIR)/interpolation_mod.$(FC_MODEXT) \
	$(FC_MODDIR)/passive_imaging_format_mod.$(FC_MODEXT) \
	$(FC_MODDIR)/rotations_mod.$(FC_MODEXT) \
	$(FC_MODDIR)/inverse_problem_par.$(FC_MODEXT) \
	$(FC_MODDIR)/precond_mod.$(FC_MODEXT) \
	$(FC_MODDIR)/signal_processing.$(FC_MODEXT) \
	$(FC_MODDIR)/elastic_tensor_tools_mod.$(FC_MODEXT) \
	$(FC_MODDIR)/inversion_scheme.$(FC_MODEXT) \
	$(FC_MODDIR)/projection_on_fd_grid.$(FC_MODEXT) \
	$(FC_MODDIR)/specfem_interface.$(FC_MODEXT) \
	$(FC_MODDIR)/family_parameter.$(FC_MODEXT) \
	$(FC_MODDIR)/io_model.$(FC_MODEXT) \
	$(FC_MODDIR)/regularization.$(FC_MODEXT) \
	$(FC_MODDIR)/teleseismic_io_mod.$(FC_MODEXT) \
	$(FC_MODDIR)/fwi_iteration.$(FC_MODEXT) \
	$(FC_MODDIR)/iso_parameter_mod.$(FC_MODEXT) \
	$(FC_MODDIR)/regularization_fd_mod.$(FC_MODEXT) \
	$(FC_MODDIR)/vti_parameters_mod.$(FC_MODEXT) \
	$(EMPTY_MACRO)


###
### MPI
###
inverse_problem_for_model_SHARED_OBJECTS += $(COND_MPI_OBJECTS)

###
### OPENMP
###
inverse_problem_for_model_SHARED_OBJECTS += $(COND_OPENMP_OBJECTS)

###
### GPU
###

inverse_problem_for_model_SHARED_OBJECTS += $(gpu_OBJECTS)

###
### ADIOS
###

# using ADIOS files

adios_inverse_problem_for_model_OBJECTS= \
	$O/read_mesh_databases_adios.spec_adios.o \
	$O/save_forward_arrays_adios.spec_adios.o \
	$O/read_forward_arrays_adios.spec_adios.o \
	$O/save_kernels_adios.spec_adios.o \
	$(EMPTY_MACRO)

adios_inverse_problem_for_model_PREOBJECTS = \
	$O/adios_helpers_addons.shared_adios_cc.o \
	$O/adios_helpers_definitions.shared_adios.o \
	$O/adios_helpers_readers.shared_adios.o \
	$O/adios_helpers_writers.shared_adios.o \
	$O/adios_helpers.shared_adios.o \
	$(EMPTY_MACRO)

adios_inverse_problem_for_model_STUBS = \
	$O/adios_method_stubs.cc.o \
	$(EMPTY_MACRO)

# conditional adios linking
ifeq ($(ADIOS),yes)
inverse_problem_for_model_OBJECTS += $(adios_inverse_problem_for_model_OBJECTS)
inverse_problem_for_model_SHARED_OBJECTS += $(adios_inverse_problem_for_model_PREOBJECTS)
else ifeq ($(ADIOS2),yes)
inverse_problem_for_model_OBJECTS += $(adios_inverse_problem_for_model_OBJECTS)
inverse_problem_for_model_SHARED_OBJECTS += $(adios_inverse_problem_for_model_PREOBJECTS)
else
inverse_problem_for_model_OBJECTS += $(adios_inverse_problem_for_model_STUBS)
endif


# conditional asdf linking
ifeq ($(ASDF),yes)
INVERSE_LINK_FLAGS += $(ASDF_LIBS)
inverse_problem_for_model_OBJECTS += $(asdf_specfem3D_OBJECTS)
inverse_problem_for_model_SHARED_OBJECTS += $(asdf_specfem3D_SHARED_OBJECTS)
else
inverse_problem_for_model_OBJECTS += $(asdf_specfem3D_STUBS)
inverse_problem_for_model_SHARED_OBJECTS += $(asdf_specfem3D_SHARED_STUBS)
endif

# vtk
ifeq ($(VTK),yes)
inverse_problem_for_model_OBJECTS += \
	$O/vtk_window_stubs.visualcc.o \
	$(EMPTY_MACRO)
endif

#######################################

####
#### rules for executables
####

ifeq ($(HAS_GPU),yes)
INFO_GPU_INVERSE_PROBLEM="building xinverse_problem_for_model $(BUILD_VERSION_TXT)"

${E}/xinverse_problem_for_model: $(inverse_problem_for_model_OBJECTS) $(inverse_problem_for_model_SHARED_OBJECTS)
	@echo ""
	@echo $(INFO_GPU_INVERSE_PROBLEM)
	@echo ""
	${FCLINK} -o ${E}/xinverse_problem_for_model $(inverse_problem_for_model_OBJECTS) $(inverse_problem_for_model_SHARED_OBJECTS) $(MPILIBS) $(GPU_LINK) $(INVERSE_LINK_FLAGS)
	@echo ""

else

## non-cuda version
${E}/xinverse_problem_for_model: $(inverse_problem_for_model_OBJECTS) $(inverse_problem_for_model_SHARED_OBJECTS)
	@echo ""
	@echo "building xinverse_problem_for_model"
	@echo ""
	${FCLINK} -o ${E}/xinverse_problem_for_model $(inverse_problem_for_model_OBJECTS) $(inverse_problem_for_model_SHARED_OBJECTS) $(MPILIBS) $(INVERSE_LINK_FLAGS)
	@echo ""

endif



#######################################

###
### Module dependencies
###

$O/interpolation_mod.inv_input.o: $O/shared_par.shared_module.o
$O/signal_processing_mod.inv_adjoint_source.o: $O/specfem3D_par.spec_module.o

$O/mesh_tools_mod.inv_input.o: $O/inverse_problem_par.inv_par.o
$O/IO_model_mod.inv_input.o: \
	$O/mesh_tools_mod.inv_input.o \
	$O/shared_par.shared_module.o \
	$O/specfem3D_par.spec_module.o \
	$O/create_color_image.spec.o

$O/rotations_mod.inv_adjoint_source.o: $O/interpolation_mod.inv_input.o
$O/adjoint_source_mod.inv_adjoint_source.o: $O/signal_processing_mod.inv_adjoint_source.o $O/rotations_mod.inv_adjoint_source.o

$O/elastic_tensor_tools_mod.inv.o: $O/interpolation_mod.inv_input.o
$O/anisotropic_parametrisation_mod.inv.o: $O/elastic_tensor_tools_mod.inv.o

$O/passive_imaging_format_mod.inv_input.o: $O/rotations_mod.inv_adjoint_source.o
$O/Teleseismic_IO.inv_inputmpi.o: $O/mesh_tools_mod.inv_input.o $O/passive_imaging_format_mod.inv_input.o

$O/input_output_mod.inv_input.o: \
	$O/mesh_tools_mod.inv_input.o \
	$O/IO_model_mod.inv_input.o \
	$O/Teleseismic_IO.inv_inputmpi.o

$O/specfem_interface_mod.inv_specfem_interface.o: \
	$O/adjoint_source_mod.inv_adjoint_source.o \
	$O/input_output_mod.inv_input.o \
	$O/signal_processing_mod.inv_adjoint_source.o \
	$O/shared_par.shared_module.o \
	$O/specfem3D_par.spec_module.o \
	$O/create_color_image.spec.o

$O/regularization_interface.inv_regularization.o: $O/regularization_SEM_mod.inv_regularization.o $O/regularization_FD_mod.inv_regularization.o

$O/family_parameter_mod.inv_inversion.o: $O/input_output_mod.inv_input.o $O/iso_parameters.inv_inversion.o $O/vti_parameters.inv_inversion.o

$O/fwi_iteration_mod.inv_inversion.o: \
	$O/specfem_interface_mod.inv_specfem_interface.o \
	$O/family_parameter_mod.inv_inversion.o \
	$O/PrecondFWI_mod.inv_inversion.o \
	$O/regularization_interface.inv_regularization.o

$O/regularization_FD_mod.inv_regularization.o: $O/projection_on_FD_grid_mod.inv_projection.o

$O/inverse_problem_main.inv.o: \
	$O/elastic_tensor_tools_mod.inv.o \
	$O/adjoint_source_mod.inv_adjoint_source.o \
	$O/input_output_mod.inv_input.o \
	$O/regularization_SEM_mod.inv_regularization.o \
	$O/inversion_scheme_mod.inv_inversion.o \
	$O/projection_on_FD_grid_mod.inv_projection.o \
	$O/specfem_interface_mod.inv_specfem_interface.o \
	$O/fwi_iteration_mod.inv_inversion.o

$O/parallel_for_inverse_problem.invmpi.o: $(COND_MPI_OBJECTS)

####
#### rule to build each .o file
####

## main module
$O/%.inv_par.o: $S/%.f90 ${SETUP}/constants.h $O/shared_par.shared_module.o $O/specfem3D_par.spec_module.o
	${FCCOMPILE_CHECK} ${FCFLAGS_f90} -c -o $@ $<

## file object rules
$O/%.inv.o: $S/%.f90 ${SETUP}/constants.h $O/inverse_problem_par.inv_par.o
	${FCCOMPILE_CHECK} ${FCFLAGS_f90} -c -o $@ $<

$O/%.invmpi.o: $S/%.f90 ${SETUP}/constants.h $O/inverse_problem_par.inv_par.o
	${MPIFCCOMPILE_CHECK} ${FCFLAGS_f90} -c -o $@ $<

$O/%.inv_adjoint_source.o: $S/adjoint_source/%.f90 ${SETUP}/constants.h $O/inverse_problem_par.inv_par.o
	${FCCOMPILE_CHECK} ${FCFLAGS_f90} -c -o $@ $<

$O/%.inv_input.o: $S/input_output/%.f90 ${SETUP}/constants.h $O/inverse_problem_par.inv_par.o
	${MPIFCCOMPILE_CHECK} ${FCFLAGS_f90} -c -o $@ $<

$O/%.inv_inputmpi.o: $S/input_output/%.f90 ${SETUP}/constants.h $O/inverse_problem_par.inv_par.o
	${MPIFCCOMPILE_CHECK} ${FCFLAGS_f90} -c -o $@ $<

$O/%.inv_inversion.o: $S/inversion_scheme/%.f90 ${SETUP}/constants.h $O/inverse_problem_par.inv_par.o
	${FCCOMPILE_CHECK} ${FCFLAGS_f90} -c -o $@ $<

$O/%.inv_parameterisation.o: $S/parameterisation/%.f90 ${SETUP}/constants.h $O/inverse_problem_par.inv_par.o
	${FCCOMPILE_CHECK} ${FCFLAGS_f90} -c -o $@ $<

$O/%.inv_projection.o: $S/projection_on_FD_grid/%.f90 ${SETUP}/constants.h $O/inverse_problem_par.inv_par.o
	${FCCOMPILE_CHECK} ${FCFLAGS_f90} -c -o $@ $<

$O/%.inv_regularization.o: $S/regularization/%.f90 ${SETUP}/constants.h $O/inverse_problem_par.inv_par.o
	${FCCOMPILE_CHECK} ${FCFLAGS_f90} -c -o $@ $<

$O/%.inv_specfem_interface.o: $S/specfem_interface/%.F90 ${SETUP}/constants.h $O/inverse_problem_par.inv_par.o
	${FCCOMPILE_CHECK} ${FCFLAGS_f90} -c -o $@ $<

