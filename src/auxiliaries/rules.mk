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
S := ${S_TOP}/src/auxiliaries
$(auxiliaries_OBJECTS): S := ${S_TOP}/src/auxiliaries

#######################################

auxiliaries_TARGETS = \
	$E/xcombine_surf_data \
	$E/xcombine_vol_data \
	$E/xconvolve_source_timefunction \
	$E/xdetect_duplicates_stations_file \
	$E/xcreate_movie_shakemap_AVS_DX_GMT \
	$(EMPTY_MACRO)

auxiliaries_OBJECTS = \
	$O/combine_surf_data.aux.o \
	$O/combine_vol_data.aux.o \
	$O/convolve_source_timefunction.aux.o \
	$O/detect_duplicates_stations_file.aux.o \
	$O/create_movie_shakemap_AVS_DX_GMT.aux.o \
	$(EMPTY_MACRO)

# These files come from the shared directory
auxiliaries_SHARED_OBJECTS = \
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

auxiliaries_SHARED_OBJECTS += $(COND_MPI_OBJECTS)

auxiliaries_MODULES = \
	$(FC_MODDIR)/combine_vol_data_mod.$(FC_MODEXT) \
	$(FC_MODDIR)/combine_vtk.$(FC_MODEXT) \
	$(EMPTY_MACRO)



##
## combine_surf_data
##
combine_surf_data_auxiliaries_OBJECTS = \
	$O/combine_surf_data.aux.o \
	$(EMPTY_MACRO)

combine_surf_data_auxiliaries_SHARED_OBJECTS = \
	$O/shared_par.shared_module.o \
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
	$O/shared_par.shared_module.o \
	$O/read_parameter_file.shared.o \
	$O/read_value_parameters.shared.o \
	$O/param_reader.cc.o \
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
	$O/shared_par.shared_module.o \
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

detect_duplicates_stations_file: xdetect_duplicates_stations_file
xdetect_duplicates_stations_file: $E/xdetect_duplicates_stations_file

combine_surf_data: xcombine_surf_data
xcombine_surf_data: $E/xcombine_surf_data

combine_vol_data: xcombine_vol_data
xcombine_vol_data: $E/xcombine_vol_data

create_movie_shakemap_AVS_DX_GMT: xcreate_movie_shakemap_AVS_DX_GMT
xcreate_movie_shakemap_AVS_DX_GMT: $E/xcreate_movie_shakemap_AVS_DX_GMT


$E/xconvolve_source_timefunction: $O/convolve_source_timefunction.aux.o $O/shared_par.shared_module.o
	${FCCOMPILE_CHECK} -o  ${E}/xconvolve_source_timefunction $O/convolve_source_timefunction.aux.o $O/shared_par.shared_module.o

$E/xdetect_duplicates_stations_file: $O/detect_duplicates_stations_file.aux.o $O/shared_par.shared_module.o
	${FCCOMPILE_CHECK} -o  ${E}/xdetect_duplicates_stations_file $O/detect_duplicates_stations_file.aux.o $O/shared_par.shared_module.o

$E/xcombine_surf_data: $(combine_surf_data_auxiliaries_OBJECTS) $(combine_surf_data_auxiliaries_SHARED_OBJECTS)
	${FCLINK} -o $@ $+

$E/xcombine_vol_data: $(combine_vol_data_auxiliaries_OBJECTS) $(combine_vol_data_auxiliaries_SHARED_OBJECTS) $(COND_MPI_OBJECTS)
	${FCLINK} -o $@ $+ $(MPILIBS)

$E/xcreate_movie_shakemap_AVS_DX_GMT: $(create_movie_shakemap_AVS_DX_GMT_auxiliaries_OBJECTS) $(create_movie_shakemap_AVS_DX_GMT_auxiliaries_SHARED_OBJECTS)
	${FCLINK} -o $@ $+




#######################################

###
### Module dependencies
###

# combine_vol_data
$O/combine_vol_data_adios_stubs.aux_noadios.o: $O/adios_manager_stubs.shared_noadios.o
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

$O/%.aux.o: $S/%.f90 $O/shared_par.shared_module.o
	${FCCOMPILE_CHECK} ${FCFLAGS_f90} -c -o $@ $<


###
### ADIOS compilation
###

$O/%.aux_adios.o: $S/%.f90 $O/shared_par.shared_module.o
	${FCCOMPILE_CHECK} ${FCFLAGS_f90} -c -o $@ $<

$O/%.aux_noadios.o: $S/%.f90 $O/shared_par.shared_module.o
	${FCCOMPILE_CHECK} ${FCFLAGS_f90} -c -o $@ $<

