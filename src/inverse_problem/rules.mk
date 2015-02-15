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
S := ${S_TOP}/src/inverse_problem
$(inverse_problem_OBJECTS): S := ${S_TOP}/src/inverse_problem

#######################################

inverse_problem_TARGETS = \
	$E/xprogram01_add_model_iso \
	$E/xprogram02_model_update \
	$E/xprogram03_smooth_sem \
	$E/xprogram04_sum_kernels \
	$E/xprogram05_sum_preconditioned_kernels \
	$(EMPTY_MACRO)

inverse_problem_OBJECTS = \
	$(xprogram01_add_model_iso_OBJECTS) \
	$(xprogram02_model_update_OBJECTS) \
	$(xprogram03_smooth_sem_OBJECTS) \
	$(xprogram04_sum_kernels_OBJECTS) \
	$(xprogram05_sum_preconditioned_kernels_OBJECTS) \
	$(EMPTY_MACRO)

# These files come from the shared directory
inverse_problem_SHARED_OBJECTS = \
	$(xprogram01_add_model_SHARED_OBJECTS) \
	$(xprogram02_model_update_SHARED_OBJECTS) \
	$(xprogram03_smooth_sem_SHARED_OBJECTS) \
	$(xprogram04_sum_kernels_SHARED_OBJECTS) \
	$(xprogram05_sum_preconditioned_kernels_SHARED_OBJECTS) \
	$(EMPTY_MACRO)

inverse_problem_MODULES = \
	$(FC_MODDIR)/inverse_problem_par.$(FC_MODEXT) \
	$(FC_MODDIR)/inverse_problem_kernels_iso.$(FC_MODEXT) \
	$(FC_MODDIR)/inverse_problem_kernels_tiso.$(FC_MODEXT) \
	$(FC_MODDIR)/inverse_problem_kernels_tiso_cg.$(FC_MODEXT) \
	$(FC_MODDIR)/inverse_problem_model_tiso.$(FC_MODEXT) \
	$(FC_MODDIR)/inverse_problem_model_iso.$(FC_MODEXT) \
	$(EMPTY_MACRO)

####
#### rules for executables
####

.PHONY: all_invprob invprob inverse_problem

all_invprob: $(inverse_problem_TARGETS)

invprob: $(inverse_problem_TARGETS)

inverse_problem: $(inverse_problem_TARGETS)

## single targets
program01_add_model_iso: xprogram01_add_model_iso
xprogram01_add_model_iso: $E/xprogram01_add_model_iso

program02_model_update: xprogram02_model_update
xprogram02_model_update: $E/xprogram02_model_update

program03_smooth_sem: xprogram03_smooth_sem
xprogram03_smooth_sem: $E/xprogram03_smooth_sem

program04_sum_kernels: xprogram04_sum_kernels
xprogram04_sum_kernels: $E/xprogram04_sum_kernels

program05_sum_preconditioned_kernels: xprogram05_sum_preconditioned_kernels
xprogram05_sum_preconditioned_kernels: $E/xprogram05_sum_preconditioned_kernels


#######################################

####
#### rules for each program follow
####

#######################################

##
## add_model
##
xprogram01_add_model_OBJECTS = \
	$O/inverse_problem_par.invprob_module.o \
	$O/subroutine01_compute_kernel_integral.invprob.o \
	$O/subroutine02_get_gradient_cg.invprob.o \
	$O/subroutine03_get_gradient_steepest.invprob.o \
	$O/subroutine04_read_kernels_cg.invprob.o \
	$O/subroutine05_read_kernels.invprob.o \
	$O/subroutine06_read_model.invprob.o \
	$O/subroutine07_read_parameters_invprob.invprob.o \
	$O/subroutine09_write_gradients.invprob.o \
	$O/subroutine10_write_new_model.invprob.o \
	$O/subroutine11_write_new_model_perturbations.invprob.o \
	$(EMPTY_MACRO)

xprogram01_add_model_SHARED_OBJECTS = \
	$O/specfem3D_par.spec.o \
	$O/pml_par.spec.o \
	$O/read_mesh_databases.spec.o \
	$O/shared_par.shared_module.o \
	$O/create_name_database.shared.o \
	$O/exit_mpi.shared.o \
	$O/gll_library.shared.o \
	$O/param_reader.cc.o \
	$O/read_parameter_file.shared.o \
	$O/read_value_parameters.shared.o \
	$(EMPTY_MACRO)

##
## xprogram01_add_model_iso
##
xprogram01_add_model_iso_OBJECTS = \
	$O/program01_add_model_iso.invprob.o \
	$(xprogram01_add_model_OBJECTS) \
	$(EMPTY_MACRO)

# extra dependencies
$O/program01_add_model_iso.invprob.o: $O/specfem3D_par.spec.o $O/inverse_problem_par.invprob_module.o

${E}/xprogram01_add_model_iso: $(xprogram01_add_model_iso_OBJECTS) $(xprogram01_add_model_SHARED_OBJECTS) $(COND_MPI_OBJECTS)
	${FCLINK} -o $@ $+ $(MPILIBS)


##
## xprogram02_model_update
##
xprogram02_model_update_OBJECTS = \
	$O/inverse_problem_par.invprob_module.o \
	$O/subroutine03_get_gradient_steepest.invprob.o \
	$O/program02_model_update.invprob.o \
	$O/subroutine05_read_kernels.invprob.o \
	$O/subroutine07_read_parameters_invprob.invprob.o \
	$O/subroutine08_save_external_bin_m_up.invprob.o \
	$O/subroutine09_write_gradients.invprob.o \
	$O/subroutine10_write_new_model.invprob.o \
	$O/subroutine11_write_new_model_perturbations.invprob.o \
	$(EMPTY_MACRO)

xprogram02_model_update_SHARED_OBJECTS = \
	$O/specfem3D_par.spec.o \
	$O/pml_par.spec.o \
	$O/initialize_simulation.spec.o \
	$O/read_mesh_databases.spec.o \
	$O/shared_par.shared_module.o \
	$O/check_mesh_resolution.shared.o \
	$O/create_name_database.shared.o \
	$O/exit_mpi.shared.o \
	$O/get_attenuation_model.shared.o \
	$O/gll_library.shared.o \
	$O/param_reader.cc.o \
	$O/read_parameter_file.shared.o \
	$O/read_value_parameters.shared.o \
	$O/write_VTK_data.shared.o \
	$(EMPTY_MACRO)

# cuda stubs
xprogram02_model_update_OBJECTS += $O/specfem3D_gpu_cuda_method_stubs.cudacc.o


# using ADIOS files
adios_program02_model_update_OBJECTS= \
	$O/read_mesh_databases_adios.spec_adios.o \
	$O/read_forward_arrays_adios.spec_adios.o

adios_program02_model_update_SHARED_OBJECTS = \
	$O/adios_manager.shared_adios.o \
	$O/adios_helpers_definitions.shared_adios_module.o \
	$O/adios_helpers_writers.shared_adios_module.o \
	$O/adios_helpers.shared_adios.o

adios_program02_model_update_STUBS = \
	$O/specfem3D_adios_stubs.spec_noadios.o

adios_program02_model_update_SHARED_STUBS = \
	$O/adios_manager_stubs.shared_noadios.o

# conditional adios linking
ifeq ($(ADIOS),yes)
xprogram02_model_update_OBJECTS += $(adios_program02_model_update_OBJECTS)
xprogram02_model_update_SHARED_OBJECTS += $(adios_program02_model_update_SHARED_OBJECTS)
else
xprogram02_model_update_OBJECTS += $(adios_program02_model_update_STUBS)
xprogram02_model_update_SHARED_OBJECTS += $(adios_program02_model_update_SHARED_STUBS)
endif

# extra dependencies
$O/program02_model_update.invprob.o: $O/specfem3D_par.spec.o $O/inverse_problem_par.invprob_module.o
$O/subroutine08_save_external_bin_m_up.invprob.o: $O/specfem3D_par.spec.o

${E}/xprogram02_model_update: $(xprogram02_model_update_OBJECTS) $(xprogram02_model_update_SHARED_OBJECTS) $(COND_MPI_OBJECTS)
	${FCLINK} -o $@ $+ $(MPILIBS)


##
## xprogram03_smooth_sem
##
xprogram03_smooth_sem_OBJECTS = \
	$O/inverse_problem_par.invprob_module.o \
	$O/program03_smooth_sem.invprob.o \
	$(EMPTY_MACRO)

xprogram03_smooth_sem_SHARED_OBJECTS = \
	$O/specfem3D_par.spec.o \
	$O/pml_par.spec.o \
	$O/read_mesh_databases.spec.o \
	$O/shared_par.shared_module.o \
	$O/check_mesh_resolution.shared.o \
	$O/create_name_database.shared.o \
	$O/exit_mpi.shared.o \
	$O/gll_library.shared.o \
	$O/param_reader.cc.o \
	$O/read_parameter_file.shared.o \
	$O/read_value_parameters.shared.o \
	$O/write_VTK_data.shared.o \
	$(EMPTY_MACRO)

# extra dependencies
$O/program03_smooth_sem.invprob.o: $O/specfem3D_par.spec.o $O/inverse_problem_par.invprob_module.o

${E}/xprogram03_smooth_sem: $(xprogram03_smooth_sem_OBJECTS) $(xprogram03_smooth_sem_SHARED_OBJECTS) $(COND_MPI_OBJECTS)
	${FCLINK} -o $@ $+ $(MPILIBS)


##
## xprogram04_sum_kernels
##
xprogram04_sum_kernels_OBJECTS = \
	$O/inverse_problem_par.invprob_module.o \
	$O/program04_sum_kernels.invprob.o \
	$(EMPTY_MACRO)

xprogram04_sum_kernels_SHARED_OBJECTS = \
	$O/shared_par.shared_module.o \
	$O/param_reader.cc.o \
	$O/read_parameter_file.shared.o \
	$O/read_value_parameters.shared.o \
	$(EMPTY_MACRO)

${E}/xprogram04_sum_kernels: $(xprogram04_sum_kernels_OBJECTS) $(xprogram04_sum_kernels_SHARED_OBJECTS) $(COND_MPI_OBJECTS)
	${FCLINK} -o $@ $+ $(MPILIBS)


##
## xprogram05_sum_preconditioned_kernels
##
xprogram05_sum_preconditioned_kernels_OBJECTS = \
	$O/inverse_problem_par.invprob_module.o \
	$O/program05_sum_preconditioned_kernels.invprob.o \
	$(EMPTY_MACRO)

xprogram05_sum_preconditioned_kernels_SHARED_OBJECTS = $(xprogram04_sum_kernels_SHARED_OBJECTS)

${E}/xprogram05_sum_preconditioned_kernels: $(xprogram05_sum_preconditioned_kernels_OBJECTS) $(xprogram05_sum_preconditioned_kernels_SHARED_OBJECTS) $(COND_MPI_OBJECTS)
	${FCLINK} -o $@ $+ $(MPILIBS)


#######################################

###
### Module dependencies
###
$O/inverse_problem_par.invprob_module.o: $O/shared_par.shared_module.o

####
#### rule for each .o file below
####

$O/%.invprob_module.o: $S/%.f90 ${SETUP}/constants_inverse_problem.h $O/shared_par.shared_module.o
	${FCCOMPILE_CHECK} ${FCFLAGS_f90} -c -o $@ $<

$O/%.invprob.o: $S/%.f90 ${SETUP}/constants_inverse_problem.h $O/inverse_problem_par.invprob_module.o
	${FCCOMPILE_CHECK} ${FCFLAGS_f90} -c -o $@ $<

