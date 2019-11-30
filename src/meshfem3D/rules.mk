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
S := ${S_TOP}/src/meshfem3D
$(meshfem3D_OBJECTS): S = ${S_TOP}/src/meshfem3D

#######################################

####
#### targets
####

meshfem3D_TARGETS = \
	$E/xmeshfem3D \
	$(EMPTY_MACRO)

meshfem3D_OBJECTS = \
	$O/calc_gll_points.mesh.o \
	$O/check_mesh_quality.mesh.o \
	$O/chunk_earth_mesh_mod.mesh.o \
	$O/compute_parameters.mesh.o \
	$O/create_meshfem_mesh.mesh.o \
	$O/create_CPML_regions.mesh.o \
	$O/create_interfaces_mesh.mesh.o \
	$O/create_visual_files.mesh.o \
	$O/define_subregions.mesh.o \
	$O/define_subregions_heuristic.mesh.o \
	$O/define_superbrick.mesh.o \
	$O/determine_cavity.mesh.o \
	$O/earth_chunk.mesh.o \
	$O/get_flags_boundaries.mesh.o \
	$O/get_MPI_cutplanes_eta.mesh.o \
	$O/get_MPI_cutplanes_xi.mesh.o \
	$O/meshfem3D.mesh.o \
	$O/meshfem3D_par.mesh_module.o \
	$O/read_mesh_parameter_file.mesh.o \
	$O/read_value_mesh_parameters.mesh.o \
	$O/save_databases.mesh.o \
	$O/store_boundaries.mesh.o \
	$O/store_coords.mesh.o \
	$(EMPTY_MACRO)

meshfem3D_MODULES = \
	$(FC_MODDIR)/constants_meshfem3d.$(FC_MODEXT) \
	$(FC_MODDIR)/meshfem3d_par.$(FC_MODEXT) \
	$(FC_MODDIR)/chunk_earth_mod.$(FC_MODEXT) \
	$(FC_MODDIR)/create_meshfem_par.$(FC_MODEXT) \
	$(EMPTY_MACRO)

meshfem3D_SHARED_OBJECTS = \
	$O/create_name_database.shared.o \
	$O/shared_par.shared_module.o \
	$O/exit_mpi.shared.o \
	$O/get_global.shared.o \
	$O/get_shape3D.shared.o \
	$O/gll_library.shared.o \
	$O/hex_nodes.shared.o \
	$O/param_reader.cc.o \
	$O/read_parameter_file.shared.o \
	$O/read_topo_bathy_file.shared.o \
	$O/read_value_parameters.shared.o \
	$O/safe_alloc_mod.shared.o \
	$O/sort_array_coordinates.shared.o \
	$O/utm_geo.shared.o \
	$O/write_VTK_data.shared.o \
	$(EMPTY_MACRO)


# using ADIOS files
adios_meshfem3D_PREOBJECTS= \
	$O/adios_manager.shared_adios.o \
	$O/adios_helpers_definitions.shared_adios_module.o \
	$O/adios_helpers_writers.shared_adios_module.o \
	$O/adios_helpers.shared_adios.o

adios_meshfem3D_OBJECTS= \
	$O/save_databases_adios.mesh_adios.o

adios_meshfem3D_PRESTUBS = \
	$O/adios_manager_stubs.shared_noadios.o

adios_meshfem3D_STUBS = \
	$O/meshfem3D_adios_stubs.mesh_noadios.o

# conditional adios linking
ifeq ($(ADIOS),no)
adios_meshfem3D_OBJECTS = $(adios_meshfem3D_STUBS)
adios_meshfem3D_PREOBJECTS = $(adios_meshfem3D_PRESTUBS)
endif
meshfem3D_OBJECTS += $(adios_meshfem3D_OBJECTS)
meshfem3D_SHARED_OBJECTS += $(adios_meshfem3D_PREOBJECTS)

# objects for the pure Fortran version
XMESHFEM_OBJECTS = \
	$(meshfem3D_OBJECTS) $(meshfem3D_SHARED_OBJECTS) \
	$(COND_MPI_OBJECTS)


#######################################


####
#### rules for executables
####

mesh: $(meshfem3D_TARGETS)

meshfem3D: xmeshfem3D
xmeshfem3D: $E/xmeshfem3D

# rules for the pure Fortran version
$E/xmeshfem3D: $(XMESHFEM_OBJECTS)
	@echo ""
	@echo "building xmeshfem3D"
	@echo ""
	${FCLINK} -o ${E}/xmeshfem3D $(XMESHFEM_OBJECTS) $(MPILIBS)
	@echo ""



#######################################


###
### Module dependencies
###

$O/meshfem3D.mesh.o: $O/chunk_earth_mesh_mod.mesh.o
$O/determine_cavity.mesh.o: $O/create_meshfem_mesh.mesh.o

## adios
$O/meshfem3D_adios_stubs.mesh_noadios.o: $O/shared_par.shared_module.o $O/adios_manager_stubs.shared_noadios.o

$O/save_databases_adios.mesh_adios.o: $O/safe_alloc_mod.shared.o $(adios_meshfem3D_PREOBJECTS)
$O/create_meshfem_mesh.mesh.o: $(adios_meshfem3D_PREOBJECTS)

$O/adios_helpers.shared_adios.o: \
	$O/adios_helpers_definitions.shared_adios_module.o \
	$O/adios_helpers_writers.shared_adios_module.o


####
#### rule to build each .o file below
####

$O/%.mesh_module.o: $S/%.f90 $O/shared_par.shared_module.o $S/constants_meshfem3D.h
	${FCCOMPILE_CHECK} ${FCFLAGS_f90} -c -o $@ $<

$O/%.mesh.o: $S/%.f90 $O/shared_par.shared_module.o $O/meshfem3D_par.mesh_module.o
	${FCCOMPILE_CHECK} ${FCFLAGS_f90} -c -o $@ $<

$O/%.mesh.o: $S/%.F90 $O/shared_par.shared_module.o $O/meshfem3D_par.mesh_module.o
	${FCCOMPILE_CHECK} ${FCFLAGS_f90} -c -o $@ $<


###
### ADIOS compilation
###

$O/%.mesh_adios.o: $S/%.F90 $O/shared_par.shared_module.o $O/meshfem3D_par.mesh_module.o
	${FCCOMPILE_CHECK} ${FCFLAGS_f90} -c -o $@ $<

$O/%.mesh_adios.o: $S/%.f90 $O/shared_par.shared_module.o $O/meshfem3D_par.mesh_module.o
	${FCCOMPILE_CHECK} ${FCFLAGS_f90}  -c -o $@ $<

$O/%.mesh_noadios.o: $S/%.F90 $O/meshfem3D_par.mesh_module.o
	${FCCOMPILE_CHECK} ${FCFLAGS_f90} -c -o $@ $<

$O/%.mesh_noadios.o: $S/%.f90 $O/meshfem3D_par.mesh_module.o
	${FCCOMPILE_CHECK} ${FCFLAGS_f90} -c -o $@ $<

