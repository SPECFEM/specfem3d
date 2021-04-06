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
S := ${S_TOP}/src/tomography
$(tomography_OBJECTS): S := ${S_TOP}/src/tomography

#######################################

tomography_TARGETS = \
	$E/xadd_model_iso \
	$E/xmodel_update \
	$E/xsum_kernels \
	$E/xsum_preconditioned_kernels \
	$(EMPTY_MACRO)

tomography_OBJECTS = \
	$(xadd_model_iso_OBJECTS) \
	$(xmodel_update_OBJECTS) \
	$(xsum_kernels_OBJECTS) \
	$(xsum_preconditioned_kernels_OBJECTS) \
	$(EMPTY_MACRO)

# These files come from the shared directory
tomography_SHARED_OBJECTS = \
	$(xadd_model_SHARED_OBJECTS) \
	$(xmodel_update_SHARED_OBJECTS) \
	$(xsum_kernels_SHARED_OBJECTS) \
	$(xsum_preconditioned_kernels_SHARED_OBJECTS) \
	$(EMPTY_MACRO)

tomography_MODULES = \
	$(FC_MODDIR)/tomography_par.$(FC_MODEXT) \
	$(FC_MODDIR)/tomography_kernels_iso.$(FC_MODEXT) \
	$(FC_MODDIR)/tomography_kernels_tiso.$(FC_MODEXT) \
	$(FC_MODDIR)/tomography_kernels_tiso_cg.$(FC_MODEXT) \
	$(FC_MODDIR)/tomography_model_tiso.$(FC_MODEXT) \
	$(FC_MODDIR)/tomography_model_iso.$(FC_MODEXT) \
	$(EMPTY_MACRO)

####
#### rules for executables
####

.PHONY: all_tomo tomo tomography

all_tomo: $(tomography_TARGETS)

tomo: $(tomography_TARGETS)

tomography: $(tomography_TARGETS)


### single targets
add_model_iso: xadd_model_iso
xadd_model_iso: $E/xadd_model_iso

model_update: xmodel_update
xmodel_update: $E/xmodel_update

sum_kernels: xsum_kernels
xsum_kernels: $E/xsum_kernels

sum_preconditioned_kernels: xsum_preconditioned_kernels
xsum_preconditioned_kernels: $E/xsum_preconditioned_kernels



#######################################

####
#### rules for each program follow
####

#######################################

##
## add_model
##
xadd_model_OBJECTS = \
	$O/tomography_par.tomo_module.o \
	$O/compute_kernel_integral.tomo.o \
	$O/get_cg_direction.tomo.o \
	$O/get_sd_direction.tomo.o \
	$O/read_kernels.tomo.o \
	$O/read_kernels_cg.tomo.o \
	$O/read_model.tomo.o \
	$O/read_parameters_tomo.tomo.o \
	$O/write_gradient.tomo.o \
	$O/write_new_model.tomo.o \
	$O/write_new_model_perturbations.tomo.o \
	$(EMPTY_MACRO)

xadd_model_SHARED_OBJECTS = \
	$O/specfem3D_par.spec_module.o \
	$O/pml_par.spec_module.o \
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
## xadd_model_iso
##
xadd_model_iso_OBJECTS = \
	$O/add_model_iso.tomo.o \
	$(xadd_model_OBJECTS) \
	$(EMPTY_MACRO)

# extra dependencies
$O/add_model_iso.tomo.o: $O/specfem3D_par.spec_module.o $O/tomography_par.tomo_module.o

${E}/xadd_model_iso: $(xadd_model_iso_OBJECTS) $(xadd_model_SHARED_OBJECTS) $(COND_MPI_OBJECTS)
	@echo ""
	@echo "building xadd_model_iso"
	@echo ""
	${FCLINK} -o $@ $+ $(MPILIBS)
	@echo ""


##
## xmodel_update
##
xmodel_update_OBJECTS = \
	$O/tomography_par.tomo_module.o \
	$O/get_sd_direction.tomo.o \
	$O/model_update.tomo.o \
	$O/read_kernels.tomo.o \
	$O/read_parameters_tomo.tomo.o \
	$O/save_external_bin_m_up.tomo.o \
	$O/write_gradient.tomo.o \
	$O/write_new_model.tomo.o \
	$O/write_new_model_perturbations.tomo.o \
	$(EMPTY_MACRO)

xmodel_update_SHARED_OBJECTS = \
	$O/specfem3D_par.spec_module.o \
	$O/pml_par.spec_module.o \
	$O/initialize_simulation.spec.o \
	$O/read_mesh_databases.spec.o \
	$O/shared_par.shared_module.o \
	$O/check_mesh_resolution.shared.o \
	$O/create_name_database.shared.o \
	$O/exit_mpi.shared.o \
	$O/get_attenuation_model.shared.o \
	$O/gll_library.shared.o \
	$O/init_openmp.shared.o \
	$O/param_reader.cc.o \
	$O/read_parameter_file.shared.o \
	$O/read_value_parameters.shared.o \
	$O/write_VTK_data.shared.o \
	$(EMPTY_MACRO)

# cuda stubs
xmodel_update_OBJECTS += $(gpu_specfem3D_STUBS)


# using ADIOS files
adios_model_update_OBJECTS= \
	$O/read_mesh_databases_adios.spec_adios.o \
	$O/read_forward_arrays_adios.spec_adios.o

adios_model_update_SHARED_OBJECTS = \
	$O/adios_manager.shared_adios.o \
	$O/adios_helpers_definitions.shared_adios_module.o \
	$O/adios_helpers_writers.shared_adios_module.o \
	$O/adios_helpers.shared_adios.o

adios_model_update_STUBS = \
	$O/specfem3D_adios_stubs.spec_noadios.o

adios_model_update_SHARED_STUBS = \
	$O/adios_manager_stubs.shared_noadios.o

# conditional adios linking
ifeq ($(ADIOS),yes)
xmodel_update_OBJECTS += $(adios_model_update_OBJECTS)
xmodel_update_SHARED_OBJECTS += $(adios_model_update_SHARED_OBJECTS)
else
xmodel_update_OBJECTS += $(adios_model_update_STUBS)
xmodel_update_SHARED_OBJECTS += $(adios_model_update_SHARED_STUBS)
endif

###
### ASDF
###

asdf_specfem3D_SHARED_STUBS = $(asdf_shared_STUBS)

# conditional asdf linking
xmodel_update_SHARED_OBJECTS += $(asdf_specfem3D_SHARED_STUBS)

# extra dependencies
$O/model_update.tomo.o: $O/specfem3D_par.spec_module.o $O/tomography_par.tomo_module.o
$O/save_external_bin_m_up.tomo.o: $O/specfem3D_par.spec_module.o

${E}/xmodel_update: $(xmodel_update_OBJECTS) $(xmodel_update_SHARED_OBJECTS) $(COND_MPI_OBJECTS)
	@echo ""
	@echo "building xmodel_update"
	@echo ""
	${FCLINK} -o $@ $+ $(MPILIBS)
	@echo ""


##
## xsum_kernels
##
xsum_kernels_OBJECTS = \
	$O/tomography_par.tomo_module.o \
	$O/sum_kernels.tomo.o \
	$(EMPTY_MACRO)

xsum_kernels_SHARED_OBJECTS = \
	$O/shared_par.shared_module.o \
	$O/param_reader.cc.o \
	$O/read_parameter_file.shared.o \
	$O/read_value_parameters.shared.o \
	$(EMPTY_MACRO)

${E}/xsum_kernels: $(xsum_kernels_OBJECTS) $(xsum_kernels_SHARED_OBJECTS) $(COND_MPI_OBJECTS)
	@echo ""
	@echo "building xsum_kernels"
	@echo ""
	${FCLINK} -o $@ $+ $(MPILIBS)
	@echo ""


##
## xsum_preconditioned_kernels
##
xsum_preconditioned_kernels_OBJECTS = \
	$O/tomography_par.tomo_module.o \
	$O/sum_preconditioned_kernels.tomo.o \
	$(EMPTY_MACRO)

xsum_preconditioned_kernels_SHARED_OBJECTS = $(xsum_kernels_SHARED_OBJECTS)

${E}/xsum_preconditioned_kernels: $(xsum_preconditioned_kernels_OBJECTS) $(xsum_preconditioned_kernels_SHARED_OBJECTS) $(COND_MPI_OBJECTS)
	@echo ""
	@echo "building xsum_preconditioned_kernels"
	@echo ""
	${FCLINK} -o $@ $+ $(MPILIBS)
	@echo ""


#######################################

###
### Module dependencies
###
$O/tomography_par.tomo_module.o: $O/shared_par.shared_module.o

####
#### rule for each .o file below
####

$O/%.tomo_module.o: $S/%.f90 ${SETUP}/constants_tomography.h $O/shared_par.shared_module.o
	${FCCOMPILE_CHECK} ${FCFLAGS_f90} -c -o $@ $<

$O/%.tomo.o: $S/%.f90 ${SETUP}/constants_tomography.h $O/tomography_par.tomo_module.o
	${FCCOMPILE_CHECK} ${FCFLAGS_f90} -c -o $@ $<

