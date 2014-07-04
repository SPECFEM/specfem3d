#=====================================================================
#
#               S p e c f e m 3 D  V e r s i o n  2 . 1
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
S := ${S_TOP}/src/auxiliaries
$(auxiliaries_OBJECTS): S := ${S_TOP}/src/auxiliaries

#######################################

auxiliaries_TARGETS = \
	$E/xcombine_surf_data \
	$E/xcombine_vol_data \
	$E/xconvolve_source_timefunction \
	$E/xcreate_movie_shakemap_AVS_DX_GMT \
	$E/xsmooth_vol_data \
	$E/xsum_kernels \
	$(EMPTY_MACRO)
##### DK DK July 2014 commented this out to temporarily avoid a dependency problem in the build system
###########	$E/xmodel_update \

auxiliaries_OBJECTS = \
	$O/combine_surf_data.aux.o \
	$O/combine_vol_data.aux.o \
	$O/convolve_source_timefunction.aux.o \
	$O/create_movie_shakemap_AVS_DX_GMT.aux.o \
	$O/smooth_vol_data.aux.o \
	$O/sum_kernels.aux.o \
	$(EMPTY_MACRO)
##### DK DK July 2014 commented this out to temporarily avoid a dependency problem in the build system
#########	$O/save_external_bin_m_up.aux.o \
#########	$O/model_update.aux.o \

# These files come from the shared directory
auxiliaries_SHARED_OBJECTS = \
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

auxiliaries_SHARED_OBJECTS += $(COND_MPI_OBJECTS)

auxiliaries_MODULES = \
	$(FC_MODDIR)/combine_vol_data_mod.$(FC_MODEXT) \
	$(FC_MODDIR)/combine_vtk.$(FC_MODEXT) \
	$(EMPTY_MACRO)

##
## model_update
##
model_upd_auxiliaries_OBJECTS = \
	$O/specfem3D_par.spec.o \
	$O/pml_par.spec.o \
	$O/initialize_simulation.spec.o \
	$O/read_mesh_databases.spec.o \
	$(EMPTY_MACRO)
##### DK DK July 2014 commented this out to temporarily avoid a dependency problem in the build system
#########	$O/model_update.aux.o \
########	$O/save_external_bin_m_up.aux.o \

model_upd_auxiliaries_SHARED_OBJECTS = \
	$O/check_mesh_resolution.shared.o \
	$O/create_name_database.shared.o \
	$O/exit_mpi.shared.o \
	$O/get_attenuation_model.shared.o \
	$O/gll_library.shared.o \
	$O/param_reader.cc.o \
	$O/read_parameter_file.shared.o \
	$O/read_value_parameters.shared.o \
	$O/unused_mod.shared_module.o \
	$O/write_VTK_data.shared.o \
	$(EMPTY_MACRO)

# cuda stubs
model_upd_auxiliaries_OBJECTS += $O/specfem3D_gpu_cuda_method_stubs.cudacc.o

# using ADIOS files
adios_model_upd_auxiliaries_OBJECTS= \
	$O/read_mesh_databases_adios.spec_adios.o \
	$O/read_forward_arrays_adios.spec_adios.o

adios_model_upd_auxiliaries_SHARED_OBJECTS = \
	$O/adios_manager.shared_adios.o \
	$O/adios_helpers_definitions.shared_adios_module.o \
	$O/adios_helpers_writers.shared_adios_module.o \
	$O/adios_helpers.shared_adios.o

adios_model_upd_auxiliaries_STUBS = \
	$O/specfem3D_adios_stubs.spec_noadios.o

adios_model_upd_auxiliaries_SHARED_STUBS = \
	$O/adios_manager_stubs.shared_noadios.o

# conditional adios linking
ifeq ($(ADIOS),yes)
model_upd_auxiliaries_OBJECTS += $(adios_model_upd_auxiliaries_OBJECTS)
model_upd_auxiliaries_SHARED_OBJECTS += $(adios_model_upd_auxiliaries_SHARED_OBJECTS)
else
model_upd_auxiliaries_OBJECTS += $(adios_model_upd_auxiliaries_STUBS)
model_upd_auxiliaries_SHARED_OBJECTS += $(adios_model_upd_auxiliaries_SHARED_STUBS)
endif

auxiliaries_OBJECTS += $(model_upd_auxiliaries_OBJECTS)
auxiliaries_SHARED_OBJECTS += $(model_upd_auxiliaries_SHARED_OBJECTS)


##
## sum_kernels
##
sum_kernels_auxiliaries_OBJECTS = \
	$O/sum_kernels.aux.o \
	$(EMPTY_MACRO)

sum_kernels_auxiliaries_SHARED_OBJECTS = \
	$O/exit_mpi.shared.o \
	$O/param_reader.cc.o \
	$O/read_parameter_file.shared.o \
	$O/read_value_parameters.shared.o \
	$(EMPTY_MACRO)

auxiliaries_OBJECTS += $(sum_kernels_auxiliaries_OBJECTS)
auxiliaries_SHARED_OBJECTS += $(sum_kernels_auxiliaries_SHARED_OBJECTS)
auxiliaries_MODULES += $(FC_MODDIR)/sum_par.$(FC_MODEXT)


##
## smooth_vol_data
##
smooth_vol_data_auxiliaries_OBJECTS = \
	$O/smooth_vol_data.aux.o \
	$(EMPTY_MACRO)

smooth_vol_data_auxiliaries_SHARED_OBJECTS = \
	$O/exit_mpi.shared.o \
	$O/gll_library.shared.o \
	$O/param_reader.cc.o \
	$O/read_parameter_file.shared.o \
	$O/read_value_parameters.shared.o \
	$(EMPTY_MACRO)

auxiliaries_OBJECTS += $(smooth_vol_data_auxiliaries_OBJECTS)
auxiliaries_SHARED_OBJECTS += $(smooth_vol_data_auxiliaries_SHARED_OBJECTS)


##
## combine_surf_data
##
combine_surf_data_auxiliaries_OBJECTS = \
	$O/combine_surf_data.aux.o \
	$(EMPTY_MACRO)

combine_surf_data_auxiliaries_SHARED_OBJECTS = \
	$O/param_reader.cc.o \
	$O/write_c_binary.cc.o \
	$(EMPTY_MACRO)

auxiliaries_OBJECTS += $(combine_surf_data_auxiliaries_OBJECTS)
auxiliaries_SHARED_OBJECTS += $(combine_surf_data_auxiliaries_SHARED_OBJECTS)

##
## combine_vol_data
##
combine_vol_data_auxiliaries_OBJECTS = \
	$O/combine_vol_data.aux.o \
	$O/combine_vol_data_impl.aux.o \
	$(EMPTY_MACRO)

combine_vol_data_auxiliaries_SHARED_OBJECTS = \
	$O/read_parameter_file.shared.o \
	$O/read_value_parameters.shared.o \
	$O/param_reader.cc.o \
	$O/unused_mod.shared_module.o \
	$O/write_c_binary.cc.o \
	$(EMPTY_MACRO)

## ADIOS
# conditional adios linking
ifeq ($(ADIOS),yes)
combine_vol_data_auxiliaries_OBJECTS += \
	$O/combine_vol_data_adios_impl.aux_adios.o
else
combine_vol_data_auxiliaries_OBJECTS += \
	$O/combine_vol_data_adios_stubs.aux_noadios.o
combine_vol_data_auxiliaries_SHARED_OBJECTS += \
	$O/adios_manager_stubs.shared_noadios.o
endif

auxiliaries_OBJECTS += $(combine_vol_data_auxiliaries_OBJECTS)
auxiliaries_SHARED_OBJECTS += $(combine_vol_data_auxiliaries_SHARED_OBJECTS)
auxiliaries_MODULES += $(FC_MODDIR)/combine_vol_data_adios_mod.$(FC_MODEXT)

##
## create_movie_shakemap_AVS_DX_GMT
##
create_movie_shakemap_AVS_DX_GMT_auxiliaries_OBJECTS = \
	$O/create_movie_shakemap_AVS_DX_GMT.aux.o \
	$(EMPTY_MACRO)

create_movie_shakemap_AVS_DX_GMT_auxiliaries_SHARED_OBJECTS = \
	$O/get_global.shared.o \
	$O/param_reader.cc.o \
	$O/read_parameter_file.shared.o \
	$O/read_value_parameters.shared.o \
	$O/sort_array_coordinates.shared.o \
	$O/utm_geo.shared.o \
	$(EMPTY_MACRO)

auxiliaries_OBJECTS += $(create_movie_shakemap_AVS_DX_GMT_auxiliaries_OBJECTS)
auxiliaries_SHARED_OBJECTS += $(create_movie_shakemap_AVS_DX_GMT_auxiliaries_SHARED_OBJECTS)

#######################################

####
#### rules for executables
####

aux: $(auxiliaries_TARGETS)

convolve_source_timefunction: xconvolve_source_timefunction
xconvolve_source_timefunction: $E/xconvolve_source_timefunction

combine_surf_data: xcombine_surf_data
xcombine_surf_data: $E/xcombine_surf_data

combine_vol_data: xcombine_vol_data
xcombine_vol_data: $E/xcombine_vol_data

create_movie_shakemap_AVS_DX_GMT: xcreate_movie_shakemap_AVS_DX_GMT
xcreate_movie_shakemap_AVS_DX_GMT: $E/xcreate_movie_shakemap_AVS_DX_GMT

model_update: xmodel_update
xmodel_update: $E/xmodel_update

smooth_vol_data: xsmooth_vol_data
xsmooth_vol_data: $E/xsmooth_vol_data

sum_kernels: xsum_kernels
xsum_kernels: $E/xsum_kernels



$E/xconvolve_source_timefunction: $O/convolve_source_timefunction.aux.o
	${FCCOMPILE_CHECK} -o  ${E}/xconvolve_source_timefunction $O/convolve_source_timefunction.aux.o

$E/xcombine_surf_data: $(combine_surf_data_auxiliaries_OBJECTS) $(combine_surf_data_auxiliaries_SHARED_OBJECTS)
	${FCLINK} -o $@ $+

$E/xcombine_vol_data: $(combine_vol_data_auxiliaries_OBJECTS) $(combine_vol_data_auxiliaries_SHARED_OBJECTS) $(COND_MPI_OBJECTS)
	${FCLINK} -o $@ $+ $(MPILIBS)

$E/xcreate_movie_shakemap_AVS_DX_GMT: $(create_movie_shakemap_AVS_DX_GMT_auxiliaries_OBJECTS) $(create_movie_shakemap_AVS_DX_GMT_auxiliaries_SHARED_OBJECTS)
	${FCLINK} -o $@ $+

$E/xmodel_update: $(model_upd_auxiliaries_OBJECTS) $(model_upd_auxiliaries_SHARED_OBJECTS) $(COND_MPI_OBJECTS)
	${FCLINK} -o $@ $+ $(MPILIBS)

$E/xsmooth_vol_data: $(smooth_vol_data_auxiliaries_OBJECTS) $(smooth_vol_data_auxiliaries_SHARED_OBJECTS) $(COND_MPI_OBJECTS)
	${FCLINK} -o $@ $+ $(MPILIBS)

$E/xsum_kernels: $(sum_kernels_auxiliaries_OBJECTS) $(sum_kernels_auxiliaries_SHARED_OBJECTS) $(COND_MPI_OBJECTS)
	${FCLINK} -o $@ $+ $(MPILIBS)


#######################################

###
### Module dependencies
###

$O/combine_vol_data_adios_stubs.aux_noadios.o: $O/unused_mod.shared_module.o $O/adios_manager_stubs.shared_noadios.o
ifeq ($(ADIOS),yes)
$O/combine_vol_data.aux.o: $O/combine_vol_data_impl.aux.o $O/combine_vol_data_adios_impl.aux_adios.o
$O/adios_helpers.shared_adios.o: \
	$O/adios_helpers_definitions.shared_adios_module.o \
	$O/adios_helpers_writers.shared_adios_module.o
else
$O/combine_vol_data.aux.o: $O/combine_vol_data_impl.aux.o $O/combine_vol_data_adios_stubs.aux_noadios.o
endif


#######################################

####
#### rule for each .o file below
####

##
## auxiliaries
##

$O/%.aux.o: $S/%.f90 $O/constants_mod.shared_module.o
	${FCCOMPILE_CHECK} ${FCFLAGS_f90} -c -o $@ $<


###
### ADIOS compilation
###

$O/%.aux_adios.o: $S/%.f90 $O/constants_mod.shared_module.o
	${FCCOMPILE_CHECK} ${FCFLAGS_f90} -c -o $@ $<

$O/%.aux_noadios.o: $S/%.f90
	${FCCOMPILE_CHECK} ${FCFLAGS_f90} -c -o $@ $<

