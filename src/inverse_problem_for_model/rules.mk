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

# I make this Makefile serial for now because I do not know how to write a clean rules.mk here with all the right dependencies
# If someone knows how to do that then please do it (in this rules.mk file only; all the others in other directories are already OK)
.NOTPARALLEL:

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


inverse_problem_for_model_OBJECTS = \
	$O/specfem3D_par.spec_module.o \
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
	$O/compute_forces_acoustic_NGLL5_fast.spec.o \
	$O/compute_forces_acoustic_NGLLnot5_generic_slow.spec.o \
	$O/compute_forces_viscoelastic_calling_routine.spec.o \
	$O/compute_forces_viscoelastic.spec.o \
	$O/compute_element_att_memory.spec.o \
	$O/compute_forces_poro_fluid_part.spec.o \
	$O/compute_forces_poroelastic_calling_routine.spec.o \
	$O/compute_forces_poro_solid_part.spec.o \
	$O/compute_gradient_in_acoustic.spec.o \
	$O/compute_interpolated_dva.spec.o \
	$O/compute_kernels.spec.o \
	$O/compute_stacey_acoustic.spec.o \
	$O/compute_stacey_viscoelastic.spec.o \
	$O/compute_stacey_poroelastic.spec.o \
	$O/compute_energy.spec.o \
	$O/convert_time.spec.o \
	$O/calendar.spec.o \
	$O/create_color_image.spec.o \
	$O/detect_mesh_surfaces.spec.o \
	$O/fault_solver_common.spec.o \
	$O/fault_solver_dynamic.spec.o \
	$O/fault_solver_kinematic.spec.o \
	$O/finalize_simulation.spec.o \
	$O/get_cmt.spec.o \
	$O/get_force.spec.o \
	$O/gravity_perturbation.spec.o \
	$O/initialize_simulation.spec.o \
	$O/iterate_time.spec.o \
	$O/locate_receivers.spec.o \
	$O/locate_source.spec.o \
	$O/make_gravity.spec.o \
	$O/noise_tomography.spec.o \
	$O/pml_allocate_arrays.spec.o \
	$O/pml_output_VTKs.spec.o \
	$O/pml_compute_accel_contribution.spec.o \
	$O/pml_compute_memory_variables.spec.o \
	$O/pml_par.spec.o \
	$O/prepare_timerun.spec.o \
	$O/read_external_stf.spec.o \
	$O/read_mesh_databases.spec.o \
	$O/save_adjoint_kernels.spec.o \
	$O/setup_GLL_points.spec.o \
	$O/setup_movie_meshes.spec.o \
	$O/setup_sources_receivers.spec.o \
	$O/specfem3D.spec.o \
	$O/update_displacement_scheme.spec.o \
	$O/update_displacement_LDDRK.spec.o \
	$O/write_movie_output.spec.o \
	$O/write_output_ASCII_or_binary.spec.o \
	$O/write_output_SU.spec.o \
	$O/write_seismograms.spec.o \
	$O/inverse_problem_par.o \
	$O/signal_processing_mod.o \
	$O/parallel_for_inverse_problem.o \
	$O/adjoint_source_mod.o \
	$O/mesh_tools_mod.o \
	$O/interpolation_mod.o \
	$O/rotations_mod.o \
	$O/elastic_tensor_tools_mod.o \
	$O/anisotropic_parametrisation_mod.o \
	$O/passive_imaging_format_mod.o \
	$O/Teleseismic_IO.o \
	$O/IO_model_mod.o \
	$O/input_output_mod.o \
	$O/regularization_SEM_mod.o \
	$O/inversion_scheme_mod.o \
	$O/family_parameter_mod.o \
	$O/projection_on_FD_grid_mod.o \
	$O/regularization_FD_mod.o \
	$O/regularization_interface.o \
	$O/PrecondFWI_mod.o \
	$O/specfem_interface_mod.o \
	$O/fwi_iteration_mod.o \
	$O/inverse_problem_main.o \
	$O/program_inverse_problem.o \
	$(EMPTY_MACRO)

inverse_problem_for_model_SHARED_OBJECTS = \
	$O/shared_par.shared_module.o \
	$O/assemble_MPI_scalar.shared.o \
	$O/check_mesh_resolution.shared.o \
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
	$O/write_c_binary.cc.o \
	$(EMPTY_MACRO)


inverse_problem_for_model_MODULES = \
	$(FC_MODDIR)/fault_solver_common.$(FC_MODEXT) \
	$(FC_MODDIR)/fault_solver_dynamic.$(FC_MODEXT) \
	$(FC_MODDIR)/fault_solver_kinematic.$(FC_MODEXT) \
	$(FC_MODDIR)/gravity_perturbation.$(FC_MODEXT) \
	$(FC_MODDIR)/image_pnm_par.$(FC_MODEXT) \
	$(FC_MODDIR)/pml_par.$(FC_MODEXT) \
	$(FC_MODDIR)/specfem_par.$(FC_MODEXT) \
	$(FC_MODDIR)/specfem_par_acoustic.$(FC_MODEXT) \
	$(FC_MODDIR)/specfem_par_elastic.$(FC_MODEXT) \
	$(FC_MODDIR)/specfem_par_poroelastic.$(FC_MODEXT) \
	$(FC_MODDIR)/specfem_par_movie.$(FC_MODEXT) \
	$(FC_MODDIR)/user_noise_distribution.$(FC_MODEXT) \
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
### CUDA
###
cuda_inverse_problem_for_model_OBJECTS = \
	$O/assemble_MPI_scalar_cuda.cuda.o \
	$O/assemble_MPI_vector_cuda.cuda.o \
	$O/check_fields_cuda.cuda.o \
	$O/compute_add_sources_acoustic_cuda.cuda.o \
	$O/compute_add_sources_viscoelastic_cuda.cuda.o \
	$O/compute_coupling_cuda.cuda.o \
	$O/compute_forces_acoustic_cuda.cuda.o \
	$O/compute_forces_viscoelastic_cuda.cuda.o \
	$O/compute_kernels_cuda.cuda.o \
	$O/compute_stacey_acoustic_cuda.cuda.o \
	$O/compute_stacey_viscoelastic_cuda.cuda.o \
	$O/initialize_cuda.cuda.o \
	$O/noise_tomography_cuda.cuda.o \
	$O/prepare_mesh_constants_cuda.cuda.o \
	$O/save_and_compare_cpu_vs_gpu.cudacc.o \
	$O/transfer_fields_cuda.cuda.o \
	$O/update_displacement_cuda.cuda.o \
	$O/write_seismograms_cuda.cuda.o \
	$O/fault_solver_dynamics.cuda.o \
	$(EMPTY_MACRO)

cuda_inverse_problem_for_model_STUBS = \
	$O/specfem3D_gpu_cuda_method_stubs.cudacc.o \
	$(EMPTY_MACRO)

cuda_inverse_problem_for_model_DEVICE_OBJ = \
	$O/cuda_device_obj.o \
	$(EMPTY_MACRO)

ifeq ($(CUDA),yes)
inverse_problem_for_model_OBJECTS += $(cuda_inverse_problem_for_model_OBJECTS)
ifeq ($(CUDA_PLUS),yes)
inverse_problem_for_model_OBJECTS += $(cuda_inverse_problem_for_model_DEVICE_OBJ)
endif
else
inverse_problem_for_model_OBJECTS += $(cuda_inverse_problem_for_model_STUBS)
endif

###
### ADIOS
###

# using ADIOS files

adios_inverse_problem_for_model_OBJECTS= \
	$O/read_mesh_databases_adios.spec_adios.o \
	$O/save_forward_arrays_adios.spec_adios.o \
	$O/read_forward_arrays_adios.spec_adios.o \
	$O/save_kernels_adios.spec_adios.o

adios_inverse_problem_for_model_PREOBJECTS = \
	$O/adios_manager.shared_adios.o \
	$O/adios_helpers_definitions.shared_adios_module.o \
	$O/adios_helpers_writers.shared_adios_module.o \
	$O/adios_helpers.shared_adios.o

adios_inverse_problem_for_model_STUBS = \
	$O/specfem3D_adios_stubs.spec_noadios.o

adios_inverse_problem_for_model_PRESTUBS = \
	$O/adios_manager_stubs.shared_noadios.o

# conditional adios linking
ifeq ($(ADIOS),no)
adios_inverse_problem_for_model_OBJECTS = $(adios_inverse_problem_for_model_STUBS)
adios_inverse_problem_for_model_PREOBJECTS = $(adios_inverse_problem_for_model_PRESTUBS)
endif
inverse_problem_for_model_OBJECTS += $(adios_inverse_problem_for_model_OBJECTS)
inverse_problem_for_model_SHARED_OBJECTS += $(adios_inverse_problem_for_model_PREOBJECTS)



#######################################

####
#### rules for executables
####

ifeq ($(CUDA),yes)
## cuda version
ifeq ($(CUDA_PLUS),yes)
## cuda 5x & 6x version
INFO_CUDA_INVERSE_PROBLEM="building xinverse_problem_for_model with CUDA support"
else
## cuda 4 version
INFO_CUDA_INVERSE_PROBLEM="building xinverse_problem_for_model with CUDA 4 support"
endif

${E}/xinverse_problem_for_model: $(inverse_problem_for_model_OBJECTS) $(inverse_problem_for_model_SHARED_OBJECTS)
	@echo ""
	@echo $(INFO_CUDA_INVERSE_PROBLEM)
	@echo ""
	${FCLINK} -o ${E}/xinverse_problem_for_model $(inverse_problem_for_model_OBJECTS) $(inverse_problem_for_model_SHARED_OBJECTS) $(MPILIBS) $(CUDA_LINK)
	@echo ""

else

## non-cuda version
${E}/xinverse_problem_for_model: $(inverse_problem_for_model_OBJECTS) $(inverse_problem_for_model_SHARED_OBJECTS)
	@echo ""
	@echo "building xinverse_problem_for_model"
	@echo ""
	${FCLINK} -o ${E}/xinverse_problem_for_model $(inverse_problem_for_model_OBJECTS) $(inverse_problem_for_model_SHARED_OBJECTS) $(MPILIBS)
	@echo ""

endif



#######################################

###
### Module dependencies
###

# Version file
$O/initialize_simulation.spec.o: ${SETUP}/version.fh

## pml
$O/compute_coupling_acoustic_el.spec.o: $O/pml_par.spec.o
$O/compute_coupling_viscoelastic_ac.spec.o: $O/pml_par.spec.o
$O/compute_forces_acoustic_calling_routine.spec.o: $O/pml_par.spec.o
$O/compute_forces_acoustic_NGLL5_fast.spec.o: $O/pml_par.spec.o
$O/compute_forces_acoustic_NGLLnot5_generic_slow.spec.o: $O/pml_par.spec.o
$O/compute_energy.spec.o: $O/pml_par.spec.o
$O/pml_allocate_arrays.spec.o: $O/pml_par.spec.o
$O/pml_compute_accel_contribution.spec.o: $O/pml_par.spec.o
$O/pml_compute_memory_variables.spec.o: $O/pml_par.spec.o
$O/pml_output_VTKs.spec.o: $O/pml_par.spec.o
$O/read_mesh_databases.spec.o: $O/pml_par.spec.o
$O/update_displacement_scheme.spec.o: $O/pml_par.spec.o

## fault
$O/fault_solver_dynamic.spec.o: $O/fault_solver_common.spec.o
$O/fault_solver_kinematic.spec.o: $O/fault_solver_common.spec.o
$O/compute_forces_viscoelastic.spec.o: $O/pml_par.spec.o $O/fault_solver_dynamic.spec.o
$O/compute_forces_viscoelastic_calling_routine.spec.o: $O/pml_par.spec.o $O/fault_solver_dynamic.spec.o $O/fault_solver_kinematic.spec.o

## gravity
$O/iterate_time.spec.o: $O/gravity_perturbation.spec.o
$O/prepare_timerun.spec.o: $O/pml_par.spec.o $O/fault_solver_dynamic.spec.o $O/fault_solver_kinematic.spec.o $O/gravity_perturbation.spec.o

## adios
$O/read_forward_arrays_adios.spec_adios.o: $O/pml_par.spec.o
$O/read_mesh_databases_adios.spec_adios.o: $O/pml_par.spec.o
$O/initialize_simulation.spec.o: $(adios_inverse_problem_for_model_PREOBJECTS)
$O/save_kernels_adios.spec_adios.o: $(adios_inverse_problem_for_model_PREOBJECTS)
$O/save_forward_arrays_adios.spec_adios.o: $O/pml_par.spec.o $(adios_inverse_problem_for_model_PREOBJECTS)
$O/finalize_simulation.spec.o: $O/pml_par.spec.o $O/gravity_perturbation.spec.o $(adios_inverse_problem_for_model_PREOBJECTS)
$O/specfem3D_adios_stubs.spec_noadios.o: $O/adios_manager_stubs.shared_noadios.o
$O/adios_helpers.shared_adios.o: \
	$O/adios_helpers_definitions.shared_adios_module.o \
	$O/adios_helpers_writers.shared_adios_module.o


####
#### rule to build each .o file below
####

## modules
$O/%.spec_module.o: $S/%.f90 $O/shared_par.shared_module.o
	${FCCOMPILE_CHECK} ${FCFLAGS_f90} -c -o $@ $<

$O/%.spec_module.o: $S/%.F90 $O/shared_par.shared_module.o
	${FCCOMPILE_CHECK} ${FCFLAGS_f90} -c -o $@ $<


$O/%.spec.o: $S/%.f90 $O/specfem3D_par.spec_module.o
	${FCCOMPILE_CHECK} ${FCFLAGS_f90} -c -o $@ $<

$O/%.spec.o: $S/%.F90 $O/specfem3D_par.spec_module.o
	${FCCOMPILE_CHECK} ${FCFLAGS_f90} -c -o $@ $<

###
### ADIOS compilation
###

$O/%.spec_adios.o: $S/%.F90 $O/specfem3D_par.spec_module.o
	${FCCOMPILE_CHECK} ${FCFLAGS_f90} -c -o $@ $<

$O/%.spec_adios.o: $S/%.f90 $O/specfem3D_par.spec_module.o
	${FCCOMPILE_CHECK} ${FCFLAGS_f90} -c -o $@ $<

$O/%.spec_noadios.o: $S/%.F90 $O/specfem3D_par.spec_module.o
	${FCCOMPILE_CHECK} ${FCFLAGS_f90} -c -o $@ $<

$O/%.spec_noadios.o: $S/%.f90 $O/specfem3D_par.spec_module.o
	${FCCOMPILE_CHECK} ${FCFLAGS_f90} -c -o $@ $<


###
### OpenMP compilation
###

$O/%.openmp.o: $S/%.f90 ${SETUP}/constants.h
	${FCCOMPILE_CHECK} ${FCFLAGS_f90} -c -o $@ $<

############################################

###
### Main program compilation (part added by DK and VM)
###

$O/program_inverse_problem.o: $S/program_inverse_problem.f90 $(program_inverse_problem_SHARED_OBJECTS) ${SETUP}/constants.h $O/shared_par.shared_module.o $(inverse_problem_for_model_SHARED_OBJECTS)
	${FCCOMPILE_CHECK} $(FCFLAGS_f90) $(program_inverse_problem_SHARED_OBJECTS) $S/program_inverse_problem.f90 -c -o $O/program_inverse_problem.o

$O/inverse_problem_main.o: $S/inverse_problem_main.f90 $(program_inverse_problem_SHARED_OBJECTS) ${SETUP}/constants.h $O/shared_par.shared_module.o $(inverse_problem_for_model_SHARED_OBJECTS)
	${FCCOMPILE_CHECK} $(FCFLAGS_f90) $(program_inverse_problem_SHARED_OBJECTS) $S/inverse_problem_main.f90 -c -o $O/inverse_problem_main.o

# add here all complilation for local modules
$O/inverse_problem_par.o: $S/inverse_problem_par.f90 $(program_inverse_problem_SHARED_OBJECTS) ${SETUP}/constants.h $O/shared_par.shared_module.o $(inverse_problem_for_model_SHARED_OBJECTS)
	${FCCOMPILE_CHECK} $(FCFLAGS_f90) $(program_inverse_problem_SHARED_OBJECTS) $S/inverse_problem_par.f90 -c -o $O/inverse_problem_par.o

$O/parallel_for_inverse_problem.o: $S/parallel_for_inverse_problem.f90 ${SETUP}/constants.h $O/shared_par.shared_module.o $(inverse_problem_for_model_SHARED_OBJECTS)
	${MPIFCCOMPILE_CHECK} $(FCFLAGS_f90) $(program_inverse_problem_SHARED_OBJECTS) $S/parallel_for_inverse_problem.f90 -c -o $O/parallel_for_inverse_problem.o

# additionals tools for mesh
$O/mesh_tools_mod.o : $S/input_output/mesh_tools_mod.f90 ${SETUP}/constants.h $O/shared_par.shared_module.o $(inverse_problem_for_model_SHARED_OBJECTS)
	${FCCOMPILE_CHECK} $(FCFLAGS_f90) $S/input_output/mesh_tools_mod.f90 -c -o $O/mesh_tools_mod.o

$O/interpolation_mod.o : $S/input_output/interpolation_mod.f90 ${SETUP}/constants.h $O/shared_par.shared_module.o $(inverse_problem_for_model_SHARED_OBJECTS)
	${FCCOMPILE_CHECK} $(FCFLAGS_f90) $S/input_output/interpolation_mod.f90 -c -o $O/interpolation_mod.o

# IO for teleseismic case
$O/rotations_mod.o : $S/adjoint_source/rotations_mod.f90 ${SETUP}/constants.h $O/shared_par.shared_module.o $(inverse_problem_for_model_SHARED_OBJECTS)
	${FCCOMPILE_CHECK} $(FCFLAGS_f90) $S/adjoint_source/rotations_mod.f90 -c -o $O/rotations_mod.o

$O/passive_imaging_format_mod.o : $S/input_output/passive_imaging_format_mod.f90 ${SETUP}/constants.h $O/shared_par.shared_module.o $(inverse_problem_for_model_SHARED_OBJECTS)
	${FCCOMPILE_CHECK} $(FCFLAGS_f90) $S/input_output/passive_imaging_format_mod.f90 -c -o $O/passive_imaging_format_mod.o

$O/Teleseismic_IO.o : $S/input_output/Teleseismic_IO.f90 ${SETUP}/constants.h $O/shared_par.shared_module.o $(inverse_problem_for_model_SHARED_OBJECTS)
	${MPIFCCOMPILE_CHECK} $(FCFLAGS_f90) $S/input_output/Teleseismic_IO.f90 -c -o $O/Teleseismic_IO.o

# modules for full anisotropy and tensor related operations
$O/elastic_tensor_tools_mod.o : $S/elastic_tensor_tools_mod.f90 ${SETUP}/constants.h $O/shared_par.shared_module.o $(inverse_problem_for_model_SHARED_OBJECTS)
	${MPIFCCOMPILE_CHECK} $(FCFLAGS_f90) $S/elastic_tensor_tools_mod.f90 -c -o $O/elastic_tensor_tools_mod.o

$O/anisotropic_parametrisation_mod.o : $S/anisotropic_parametrisation_mod.f90 ${SETUP}/constants.h $O/shared_par.shared_module.o $(inverse_problem_for_model_SHARED_OBJECTS)
	${MPIFCCOMPILE_CHECK} $(FCFLAGS_f90) $S/anisotropic_parametrisation_mod.f90 -c -o $O/anisotropic_parametrisation_mod.o

# import or export model
$O/IO_model_mod.o : $S/input_output/IO_model_mod.f90 ${SETUP}/constants.h $O/shared_par.shared_module.o $(inverse_problem_for_model_SHARED_OBJECTS)
	${FCCOMPILE_CHECK} $(FCFLAGS_f90) $S/input_output/IO_model_mod.f90 -c -o $O/IO_model_mod.o

# input_output module
$O/input_output_mod.o: $S/input_output/input_output_mod.f90 ${SETUP}/constants.h $O/shared_par.shared_module.o $(inverse_problem_for_model_SHARED_OBJECTS)
	${MPIFCCOMPILE_CHECK} $(FCFLAGS_f90) $S/input_output/input_output_mod.f90 -c -o $O/input_output_mod.o

# some signal processing tools need for define adjoint sources
$O/signal_processing_mod.o: $S/adjoint_source/signal_processing_mod.f90 ${SETUP}/constants.h $O/shared_par.shared_module.o $(inverse_problem_for_model_SHARED_OBJECTS)
	${FCCOMPILE_CHECK} $(FCFLAGS_f90) $S/adjoint_source/signal_processing_mod.f90 -c -o $O/signal_processing_mod.o

# adjoint source module
$O/adjoint_source_mod.o: $S/adjoint_source/adjoint_source_mod.f90 ${SETUP}/constants.h $O/shared_par.shared_module.o $(inverse_problem_for_model_SHARED_OBJECTS)
	${FCCOMPILE_CHECK} $(FCFLAGS_f90) $S/adjoint_source/adjoint_source_mod.f90 -c -o $O/adjoint_source_mod.o

# regularization module
$O/regularization_SEM_mod.o: $S/regularization/regularization_SEM_mod.f90 $O/read_partition_files.gen.o  ${SETUP}/constants.h $O/shared_par.shared_module.o $(inverse_problem_for_model_SHARED_OBJECTS)
	${MPIFCCOMPILE_CHECK} $(FCFLAGS_f90) $S/regularization/regularization_SEM_mod.f90 -c -o $O/regularization_SEM_mod.o

# inversion schemes modules
$O/inversion_scheme_mod.o: $S/inversion_scheme/inversion_scheme_mod.f90 ${SETUP}/constants.h $O/shared_par.shared_module.o $(inverse_problem_for_model_SHARED_OBJECTS)
	${FCCOMPILE_CHECK} $(FCFLAGS_f90) $S/inversion_scheme/inversion_scheme_mod.f90 -c -o $O/inversion_scheme_mod.o

# generic inversion modules
$O/family_parameter_mod.o : $S/inversion_scheme/family_parameter_mod.f90 ${SETUP}/constants.h $O/shared_par.shared_module.o $(inverse_problem_for_model_SHARED_OBJECTS)
	${FCCOMPILE_CHECK} $(FCFLAGS_f90) $S/inversion_scheme/family_parameter_mod.f90 -c -o  $O/family_parameter_mod.o

$O/fwi_iteration_mod.o : $S/inversion_scheme/fwi_iteration_mod.f90 ${SETUP}/constants.h $O/shared_par.shared_module.o $(inverse_problem_for_model_SHARED_OBJECTS)
	${FCCOMPILE_CHECK} $(FCFLAGS_f90) $S/inversion_scheme/fwi_iteration_mod.f90 -c -o $O/fwi_iteration_mod.o

$O/PrecondFWI_mod.o : $S/inversion_scheme/PrecondFWI_mod.f90 ${SETUP}/constants.h $O/shared_par.shared_module.o $(inverse_problem_for_model_SHARED_OBJECTS)
	${FCCOMPILE_CHECK} $(FCFLAGS_f90) $S/inversion_scheme/PrecondFWI_mod.f90 -c -o $O/PrecondFWI_mod.o

# projection on finite-difference (FD) grid
$O/projection_on_FD_grid_mod.o: $S/projection_on_FD_grid/projection_on_FD_grid_mod.f90 ${SETUP}/constants.h $O/shared_par.shared_module.o $(inverse_problem_for_model_SHARED_OBJECTS)
	${FCCOMPILE_CHECK} $(FCFLAGS_f90) $S/projection_on_FD_grid/projection_on_FD_grid_mod.f90 -c -o $O/projection_on_FD_grid_mod.o

$O/regularization_FD_mod.o:  $S/regularization/regularization_FD_mod.f90  ${SETUP}/constants.h $O/shared_par.shared_module.o $(inverse_problem_for_model_SHARED_OBJECTS)
	${FCCOMPILE_CHECK} $(FCFLAGS_f90) $S/regularization/regularization_FD_mod.f90 -c -o $O/regularization_FD_mod.o

$O/regularization_interface.o:  $S/regularization/regularization_interface.f90 ${SETUP}/constants.h $O/shared_par.shared_module.o $(inverse_problem_for_model_SHARED_OBJECTS)
	${FCCOMPILE_CHECK} $(FCFLAGS_f90) $S/regularization/regularization_interface.f90 -c -o $O/regularization_interface.o

# specfem interface to compile at the end
$O/specfem_interface_mod.o: $S/specfem_interface/specfem_interface_mod.F90 ${SETUP}/constants.h $O/shared_par.shared_module.o $(inverse_problem_for_model_SHARED_OBJECTS)
	${FCCOMPILE_CHECK} $(FCFLAGS_f90) $S/specfem_interface/specfem_interface_mod.F90 -c -o $O/specfem_interface_mod.o

