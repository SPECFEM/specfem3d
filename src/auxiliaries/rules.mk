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
S := ${S_TOP}/src/auxiliaries
$(auxiliaries_OBJECTS): S := ${S_TOP}/src/auxiliaries

#######################################

auxiliaries_TARGETS = \
	$E/xcombine_surf_data \
	$E/xcombine_vol_data \
	$E/xcombine_vol_data_vtk \
	$E/xcombine_vol_data_vtu \
	$E/xcombine_vol_data_vtk_bin \
	$E/xconvolve_source_timefunction \
	$E/xdetect_duplicates_stations_file \
	$E/xcreate_movie_shakemap_AVS_DX_GMT \
	$(EMPTY_MACRO)

auxiliaries_OBJECTS = \
	$(xcombine_surf_data_OBJECTS) \
	$(xcombine_vol_data_OBJECTS) \
	$(xcombine_vol_data_vtk_OBJECTS) \
	$(xcombine_vol_data_vtu_OBJECTS) \
	$(xcombine_vol_data_vtk_bin_OBJECTS) \
	$(xcreate_movie_shakemap_AVS_DX_GMT_OBJECTS) \
	$(xconvolve_source_timefunction_OBJECTS) \
	$(xdetect_duplicates_stations_file_OBJECTS) \
	$(EMPTY_MACRO)

# These files come from the shared directory
auxiliaries_SHARED_OBJECTS = \
	$(xcombine_surf_data_SHARED_OBJECTS) \
	$(xcombine_vol_data_SHARED_OBJECTS) \
	$(xcombine_vol_data_vtk_SHARED_OBJECTS) \
	$(xcombine_vol_data_vtu_SHARED_OBJECTS) \
	$(xcombine_vol_data_vtk_bin_SHARED_OBJECTS) \
	$(xcreate_movie_shakemap_AVS_DX_GMT_SHARED_OBJECTS) \
	$(xconvolve_source_timefunction_SHARED_OBJECTS) \
	$(xdetect_duplicates_stations_file_SHARED_OBJECTS) \
	$(EMPTY_MACRO)

auxiliaries_SHARED_OBJECTS += $(COND_MPI_OBJECTS)

auxiliaries_MODULES = \
	$(FC_MODDIR)/combine_vol_data_mod.$(FC_MODEXT) \
	$(FC_MODDIR)/combine_vtk_par.$(FC_MODEXT) \
	$(EMPTY_MACRO)

## targets requiring MPI
ifeq ($(MPI),yes)
auxiliaries_TARGETS += $E/xproject_and_combine_vol_data_on_regular_grid
auxiliaries_OBJECTS += $(xproject_and_combine_vol_data_on_regular_grid_OBJECTS)
auxiliaries_SHARED_OBJECTS += $(xproject_and_combine_vol_data_on_regular_grid_SHARED_OBJECTS)
endif

## targets requiring ADIOS
adios_auxiliaries_TARGETS = \
	$E/xcombine_vol_data_adios \
	$E/xcombine_vol_data_vtk_adios \
	$E/xcombine_vol_data_vtu_adios \
	$(EMPTY_MACRO)

adios_auxiliaries_OBJECTS = \
	$(xcombine_vol_data_adios_OBJECTS) \
	$(xcombine_vol_data_vtk_adios_OBJECTS) \
	$(xcombine_vol_data_vtu_adios_OBJECTS) \
	$(EMPTY_MACRO)

adios_auxiliaries_MODULES = \
	$(FC_MODDIR)/combine_vol_data_adios_mod.$(FC_MODEXT) \
	$(EMPTY_MACRO)

ifeq ($(ADIOS),yes)
auxiliaries_TARGETS += $(adios_auxiliaries_TARGETS)
auxiliaries_OBJECTS += $(adios_auxiliaries_OBJECTS)
auxiliaries_MODULES += $(adios_auxiliaries_MODULES)
else ifeq ($(ADIOS2),yes)
auxiliaries_TARGETS += $(adios_auxiliaries_TARGETS)
auxiliaries_OBJECTS += $(adios_auxiliaries_OBJECTS)
auxiliaries_MODULES += $(adios_auxiliaries_MODULES)
endif


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

combine_vol_data_vtk: xcombine_vol_data_vtk
xcombine_vol_data_vtk: $E/xcombine_vol_data_vtk

combine_vol_data_vtu: xcombine_vol_data_vtu
xcombine_vol_data_vtu: $E/xcombine_vol_data_vtu

combine_vol_data_vtk_bin: xcombine_vol_data_vtk_bin
xcombine_vol_data_vtk_bin: $E/xcombine_vol_data_vtk_bin

combine_vol_data_adios: xcombine_vol_data_adios
xcombine_vol_data_adios: $E/xcombine_vol_data_adios

combine_vol_data_vtk_adios: xcombine_vol_data_vtk_adios
xcombine_vol_data_vtk_adios: $E/xcombine_vol_data_vtk_adios

combine_vol_data_vtu_adios: xcombine_vol_data_vtu_adios
xcombine_vol_data_vtu_adios: $E/xcombine_vol_data_vtu_adios

project_and_combine_vol_data_on_regular_grid: xproject_and_combine_vol_data_on_regular_grid
xproject_and_combine_vol_data_on_regular_grid: $E/xproject_and_combine_vol_data_on_regular_grid

create_movie_shakemap_AVS_DX_GMT: xcreate_movie_shakemap_AVS_DX_GMT
xcreate_movie_shakemap_AVS_DX_GMT: $E/xcreate_movie_shakemap_AVS_DX_GMT

#######################################

####
#### rules for each program follow
####

#######################################

##
## xcombine_surf_data
##
xcombine_surf_data_OBJECTS = \
	$O/combine_surf_data.aux.o \
	$(EMPTY_MACRO)

xcombine_surf_data_SHARED_OBJECTS = \
	$O/shared_par.shared_module.o \
	$O/count_number_of_sources.shared.o \
	$O/exit_mpi.shared.o \
	$O/read_parameter_file.shared.o \
	$O/read_value_parameters.shared.o \
	$O/param_reader.cc.o \
	$O/write_c_binary.cc.o \
	$(EMPTY_MACRO)

$E/xcombine_surf_data: $(xcombine_surf_data_OBJECTS) $(xcombine_surf_data_SHARED_OBJECTS) $(COND_MPI_OBJECTS)
	@echo ""
	@echo "building xcombine_surf_data"
	@echo ""
	${FCLINK} -o $@ $+ $(MPILIBS)
	@echo ""

#######################################

##
## xcombine_vol_data
##
xcombine_vol_data_OBJECTS = \
	$O/combine_vol_data.aux.o \
	$O/combine_vol_data_impl.aux.o \
	$(EMPTY_MACRO)

xcombine_vol_data_SHARED_OBJECTS = \
	$O/shared_par.shared_module.o \
	$O/count_number_of_sources.shared.o \
	$O/exit_mpi.shared.o \
	$O/read_parameter_file.shared.o \
	$O/read_value_parameters.shared.o \
	$O/param_reader.cc.o \
	$O/write_c_binary.cc.o \
	$O/write_VTK_data.shared.o \
	$(EMPTY_MACRO)

## MPI
xcombine_vol_data_SHARED_OBJECTS += $(COND_MPI_OBJECTS)

$E/xcombine_vol_data: $(xcombine_vol_data_OBJECTS) $(xcombine_vol_data_SHARED_OBJECTS)
	@echo ""
	@echo "building xcombine_vol_data"
	@echo ""
	${FCLINK} -o $@ $+ $(MPILIBS)
	@echo ""

### additional dependencies
$O/combine_vol_data.aux.o: $O/combine_vol_data_impl.aux.o


#######################################

## ADIOS
xcombine_vol_data_adios_OBJECTS = \
	$O/combine_vol_data.auxadios.o \
	$O/combine_vol_data_impl.aux.o \
	$O/combine_vol_data_adios_impl.auxadios.o \
	$(EMPTY_MACRO)

xcombine_vol_data_adios_SHARED_OBJECTS = $(xcombine_vol_data_SHARED_OBJECTS) \
	$O/adios_helpers_addons.shared_adios_cc.o \
	$O/adios_helpers_definitions.shared_adios.o \
	$O/adios_helpers_readers.shared_adios.o \
	$O/adios_helpers_writers.shared_adios.o \
	$O/adios_helpers.shared_adios.o \
	$O/adios_manager.shared_adios_module.o \
	$(EMPTY_MACRO)


${E}/xcombine_vol_data_adios: $(xcombine_vol_data_adios_OBJECTS) $(xcombine_vol_data_adios_SHARED_OBJECTS)
	@echo ""
	@echo "building xcombine_vol_data_adios"
	@echo ""
	${FCLINK} -o $@ $+ $(MPILIBS)
	@echo ""

### additional dependencies
$O/combine_vol_data.auxadios.o: $O/combine_vol_data_impl.aux.o $O/combine_vol_data_adios_impl.auxadios.o
$O/combine_vol_data_adios_impl.auxadios.o: $O/adios_manager.shared_adios_module.o $O/adios_helpers.shared_adios.o

#######################################

##
## xcombine_vol_data_vtk
##
xcombine_vol_data_vtk_OBJECTS = \
	$O/combine_vol_data.aux_vtk.o \
	$O/combine_vol_data_impl.aux.o \
	$(EMPTY_MACRO)

$E/xcombine_vol_data_vtk: $(xcombine_vol_data_vtk_OBJECTS) $(xcombine_vol_data_SHARED_OBJECTS)
	@echo ""
	@echo "building xcombine_vol_data_vtk"
	@echo ""
	${FCLINK} -o $@ $+ $(MPILIBS)
	@echo ""

### additional dependencies
$O/combine_vol_data.aux_vtk.o: $O/combine_vol_data_impl.aux.o

#######################################

## ADIOS
xcombine_vol_data_vtk_adios_OBJECTS = \
	$O/combine_vol_data.auxadios_vtk.o \
	$O/combine_vol_data_impl.aux.o \
	$O/combine_vol_data_adios_impl.auxadios.o \
	$(EMPTY_MACRO)

xcombine_vol_data_vtk_adios_SHARED_OBJECTS = $(xcombine_vol_data_adios_SHARED_OBJECTS)

$E/xcombine_vol_data_vtk_adios: $(xcombine_vol_data_vtk_adios_OBJECTS) $(xcombine_vol_data_vtk_adios_SHARED_OBJECTS)
	@echo ""
	@echo "building xcombine_vol_data_vtk_adios"
	@echo ""
	${FCLINK} -o $@ $+ $(MPILIBS)
	@echo ""

### additional dependencies
$O/combine_vol_data.auxadios_vtk.o: $O/combine_vol_data_impl.aux.o
$O/combine_vol_data.auxadios_vtk.o: $O/combine_vol_data_adios_impl.auxadios.o

#######################################

##
## xcombine_vol_data_vtu
##
xcombine_vol_data_vtu_OBJECTS = \
	$O/combine_vol_data.aux_vtu.o \
	$O/combine_vol_data_impl.aux.o \
	$(EMPTY_MACRO)

$E/xcombine_vol_data_vtu: $(xcombine_vol_data_vtu_OBJECTS) $(xcombine_vol_data_SHARED_OBJECTS)
	@echo ""
	@echo "building xcombine_vol_data_vtu"
	@echo ""
	${FCLINK} -o $@ $+ $(MPILIBS)
	@echo ""

### additional dependencies
$O/combine_vol_data.aux_vtu.o: $O/combine_vol_data_impl.aux.o

#######################################

## ADIOS
xcombine_vol_data_vtu_adios_OBJECTS = \
	$O/combine_vol_data.auxadios_vtu.o \
	$O/combine_vol_data_impl.aux.o \
	$O/combine_vol_data_adios_impl.auxadios.o \
	$(EMPTY_MACRO)

xcombine_vol_data_vtu_adios_SHARED_OBJECTS = $(xcombine_vol_data_adios_SHARED_OBJECTS)

$E/xcombine_vol_data_vtu_adios: $(xcombine_vol_data_vtu_adios_OBJECTS) $(xcombine_vol_data_vtu_adios_SHARED_OBJECTS)
	@echo ""
	@echo "building xcombine_vol_data_vtu_adios"
	@echo ""
	${FCLINK} -o $@ $+ $(MPILIBS)
	@echo ""

### additional dependencies
$O/combine_vol_data.auxadios_vtu.o: $O/combine_vol_data_impl.aux.o
$O/combine_vol_data.auxadios_vtu.o: $O/combine_vol_data_adios_impl.auxadios.o

#######################################

##
## xcombine_vol_data_vtk_bin
##
xcombine_vol_data_vtk_bin_OBJECTS = \
	$O/combine_vol_data_vtk_binary.aux.o \
	$O/combine_vol_data_impl.aux.o \
	$O/vtk_writer.aux.o \
	$(EMPTY_MACRO)

xcombine_vol_data_vtk_bin_SHARED_OBJECTS = $(xcombine_vol_data_SHARED_OBJECTS)

$E/xcombine_vol_data_vtk_bin: $(xcombine_vol_data_vtk_bin_OBJECTS) $(xcombine_vol_data_vtk_bin_SHARED_OBJECTS)
	@echo ""
	@echo "building xcombine_vol_data_vtk_bin"
	@echo ""
	${FCLINK} -o $@ $+ $(MPILIBS)
	@echo ""

### additional dependencies
$O/combine_vol_data_vtk_binary.aux.o:  $O/combine_vol_data_impl.aux.o


#######################################

##
## xproject_and_combine_vol_data_on_regular_grid
##

## requires compilation w/ MPI

xproject_and_combine_vol_data_on_regular_grid_OBJECTS = \
	$O/project_and_combine_vol_data_on_regular_grid.aux.o \
	$O/vtk_writer.aux.o \
	$(EMPTY_MACRO)

xproject_and_combine_vol_data_on_regular_grid_SHARED_OBJECTS = \
	$O/inverse_problem_par.inv_par.o \
	$O/parallel_for_inverse_problem.invmpi.o \
	$O/specfem3D_par.spec_module.o \
	$O/projection_on_FD_grid_mod.inv_projection.o \
	$O/assemble_MPI_scalar.shared.o \
	$O/check_mesh_resolution.shared.o \
	$O/create_name_database.shared.o \
	$O/define_derivation_matrices.shared.o \
	$O/detect_surface.shared.o \
	$O/force_ftz.cc.o \
	$O/get_attenuation_model.shared.o \
	$O/get_element_face.shared.o \
	$O/get_jacobian_boundaries.shared.o \
	$O/get_shape3D.shared.o \
	$O/gll_library.shared.o \
	$O/hex_nodes.shared.o \
	$O/lagrange_poly.shared.o \
	$O/netlib_specfun_erf.shared.o \
	$O/prepare_assemble_MPI.shared.o \
	$O/read_topo_bathy_file.shared.o \
	$O/recompute_jacobian.shared.o \
	$O/save_header_file.shared.o \
	$O/sort_array_coordinates.shared.o \
	$O/utm_geo.shared.o \
	$(EMPTY_MACRO)

xproject_and_combine_vol_data_on_regular_grid_SHARED_OBJECTS += $(xcombine_vol_data_SHARED_OBJECTS)

$E/xproject_and_combine_vol_data_on_regular_grid: $(xproject_and_combine_vol_data_on_regular_grid_OBJECTS) $(xproject_and_combine_vol_data_on_regular_grid_SHARED_OBJECTS)
	@echo ""
	@echo "building xproject_and_combine_vol_data_on_regular_grid"
	@echo ""
	${FCLINK} -o $@ $+ $(MPILIBS)
	@echo ""

### additional dependencies
$O/project_and_combine_vol_data_on_regular_grid.aux.o: \
	$O/specfem3D_par.spec_module.o \
	$O/inverse_problem_par.inv_par.o \
	$O/projection_on_FD_grid_mod.inv_projection.o


#######################################

##
## xcreate_movie_shakemap_AVS_DX_GMT
##
xcreate_movie_shakemap_AVS_DX_GMT_OBJECTS = \
	$O/create_movie_shakemap_AVS_DX_GMT.aux.o \
	$(EMPTY_MACRO)

xcreate_movie_shakemap_AVS_DX_GMT_SHARED_OBJECTS = \
	$O/shared_par.shared_module.o \
	$O/count_number_of_sources.shared.o \
	$O/get_global.shared.o \
	$O/param_reader.cc.o \
	$O/read_parameter_file.shared.o \
	$O/read_value_parameters.shared.o \
	$O/sort_array_coordinates.shared.o \
	$O/utm_geo.shared.o \
	$(EMPTY_MACRO)

$E/xcreate_movie_shakemap_AVS_DX_GMT: $(xcreate_movie_shakemap_AVS_DX_GMT_OBJECTS)  $(xcreate_movie_shakemap_AVS_DX_GMT_SHARED_OBJECTS) $(COND_MPI_OBJECTS)
	@echo ""
	@echo "building xcreate_movie_shakemap_AVS_DX_GMT"
	@echo ""
	${FCLINK} -o $@ $+ $(MPILIBS)
	@echo ""

#######################################

##
## xconvolve_source_timefunction
##
xconvolve_source_timefunction_OBJECTS = \
	$O/convolve_source_timefunction.aux.o \
	$(EMPTY_MACRO)

xconvolve_source_timefunction_SHARED_OBJECTS = \
	$O/shared_par.shared_module.o \
	$O/count_number_of_sources.shared.o \
	$O/exit_mpi.shared.o \
	$O/param_reader.cc.o \
	$O/read_parameter_file.shared.o \
	$O/read_value_parameters.shared.o \
	$(EMPTY_MACRO)

$E/xconvolve_source_timefunction: $(xconvolve_source_timefunction_OBJECTS) $(xconvolve_source_timefunction_SHARED_OBJECTS) $(COND_MPI_OBJECTS)
	@echo ""
	@echo "building xconvolve_source_timefunction"
	@echo ""
	${FCLINK} -o $@ $+ $(MPILIBS)
	@echo ""

#######################################

##
## xdetect_duplicates_stations_file
##

xdetect_duplicates_stations_file_OBJECTS = \
	$O/detect_duplicates_stations_file.aux.o \
	$(EMPTY_MACRO)

xdetect_duplicates_stations_file_SHARED_OBJECTS = \
	$O/shared_par.shared_module.o \
	$O/count_number_of_sources.shared.o \
	$O/exit_mpi.shared.o \
	$O/param_reader.cc.o \
	$O/read_parameter_file.shared.o \
	$O/read_value_parameters.shared.o \
	$(EMPTY_MACRO)

$E/xdetect_duplicates_stations_file: $(xdetect_duplicates_stations_file_OBJECTS) $(xdetect_duplicates_stations_file_SHARED_OBJECTS) $(COND_MPI_OBJECTS)
	@echo ""
	@echo "building xdetect_duplicates_stations_file"
	@echo ""
	${FCLINK} -o $@ $+ $(MPILIBS)
	@echo ""

#######################################

####
#### rule for each .o file below
####

##
## auxiliaries
##

$O/%.aux.o: $S/%.F90 $O/shared_par.shared_module.o
	${FCCOMPILE_CHECK} ${FCFLAGS_f90} -c -o $@ $<

$O/%.aux.o: $S/%.f90 $O/shared_par.shared_module.o
	${FCCOMPILE_CHECK} ${FCFLAGS_f90} -c -o $@ $<

$O/%.aux_vtk.o: $S/%.F90 $O/shared_par.shared_module.o
	${FCCOMPILE_CHECK} ${FCFLAGS_f90} -c -o $@ $< $(FC_DEFINE)USE_VTK_INSTEAD_OF_MESH

$O/%.aux_vtu.o: $S/%.F90 $O/shared_par.shared_module.o
	${FCCOMPILE_CHECK} ${FCFLAGS_f90} -c -o $@ $< $(FC_DEFINE)USE_VTU_INSTEAD_OF_MESH

$O/%.aux.o: $S/%.c ${SETUP}/config.h
	${CC} -c $(CPPFLAGS) $(CFLAGS) $(MPI_INCLUDES) -o $@ $<

###
### ADIOS compilation
###

$O/%.auxadios.o: $S/%.f90 $O/shared_par.shared_module.o
	${FCCOMPILE_CHECK} ${FCFLAGS_f90} -c -o $@ $< $(FC_DEFINE)USE_ADIOS_INSTEAD_OF_MESH

$O/%.auxadios.o: $S/%.F90 $O/shared_par.shared_module.o
	${FCCOMPILE_CHECK} ${FCFLAGS_f90} -c -o $@ $< $(FC_DEFINE)USE_ADIOS_INSTEAD_OF_MESH

$O/%.auxadios_vtk.o: $S/%.F90 $O/shared_par.shared_module.o
	${FCCOMPILE_CHECK} ${FCFLAGS_f90} -c -o $@ $< $(FC_DEFINE)USE_ADIOS_INSTEAD_OF_MESH $(FC_DEFINE)USE_VTK_INSTEAD_OF_MESH

$O/%.auxadios_vtu.o: $S/%.F90 $O/shared_par.shared_module.o
	${FCCOMPILE_CHECK} ${FCFLAGS_f90} -c -o $@ $< $(FC_DEFINE)USE_ADIOS_INSTEAD_OF_MESH $(FC_DEFINE)USE_VTU_INSTEAD_OF_MESH

$O/%.aux_noadios.o: $S/%.f90 $O/shared_par.shared_module.o
	${FCCOMPILE_CHECK} ${FCFLAGS_f90} -c -o $@ $<

