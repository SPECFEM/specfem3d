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
S := ${S_TOP}/src/tomography/postprocess_sensitivity_kernels
$(tomography/postprocess_sensitivity_kernels_OBJECTS): S := ${S_TOP}/src/tomography/postprocess_sensitivity_kernels

#######################################

tomography/postprocess_sensitivity_kernels_TARGETS = \
	$E/xclip_sem \
	$E/xcombine_sem \
	$E/xsmooth_sem \
	$E/xsmooth_sem_pde \
	$(EMPTY_MACRO)

tomography/postprocess_sensitivity_kernels_OBJECTS = \
	$(xclip_sem_OBJECTS) \
	$(xcombine_sem_OBJECTS) \
	$(xsmooth_sem_OBJECTS) \
	$(xsmooth_sem_pde_OBJECTS) \
	$(EMPTY_MACRO)

tomography/postprocess_sensitivity_kernels_SHARED_OBJECTS = \
	$(xclip_sem_SHARED_OBJECTS) \
	$(xcombine_sem_SHARED_OBJECTS) \
	$(xsmooth_sem_SHARED_OBJECTS) \
	$(xsmooth_sem_pde_SHARED_OBJECTS) \
	$(EMPTY_MACRO)

tomography/postprocess_sensitivity_kernels_MODULES = \
	$(FC_MODDIR)/postprocess_par.$(FC_MODEXT) \
	$(EMPTY_MACRO)

####
#### rules for executables
####

.PHONY: postprocess


postprocess: $(tomography/postprocess_sensitivity_kernels_TARGETS)

postprocess_sensitivity_kernels: postprocess

tomography/postprocess_sensitivity_kernels: postprocess

### single targets

clip_sem: xclip_sem
xclip_sem: $E/xclip_sem

combine_sem: xcombine_sem
xcombine_sem: $E/xcombine_sem

smooth_sem: xsmooth_sem
xsmooth_sem: $E/xsmooth_sem

smooth_sem_pde: xsmooth_sem_pde
xsmooth_sem_pde: $E/xsmooth_sem_pde

#######################################

####
#### rules for each program follow
####

#######################################


##
## xclip_sem
##
xclip_sem_OBJECTS = \
	$O/postprocess_par.postprocess_module.o \
	$O/parse_kernel_names.postprocess.o \
	$O/clip_sem.postprocess.o \
	$(EMPTY_MACRO)

xclip_sem_SHARED_OBJECTS = \
	$O/shared_par.shared_module.o \
	$O/count_number_of_sources.shared.o \
	$O/exit_mpi.shared.o \
	$O/param_reader.cc.o \
	$O/read_parameter_file.shared.o \
	$O/read_value_parameters.shared.o \
	$(EMPTY_MACRO)

${E}/xclip_sem: $(xclip_sem_OBJECTS) $(xclip_sem_SHARED_OBJECTS) $(COND_MPI_OBJECTS)
	@echo ""
	@echo "building xclip_sem"
	@echo ""
	${FCLINK} -o $@ $+ $(MPILIBS)
	@echo ""


##
## xcombine_sem
##
xcombine_sem_OBJECTS = \
	$O/postprocess_par.postprocess_module.o \
	$O/parse_kernel_names.postprocess.o \
	$O/combine_sem.postprocess.o \
	$(EMPTY_MACRO)

xcombine_sem_SHARED_OBJECTS = \
	$O/shared_par.shared_module.o \
	$O/count_number_of_sources.shared.o \
	$O/exit_mpi.shared.o \
	$O/param_reader.cc.o \
	$O/read_parameter_file.shared.o \
	$O/read_value_parameters.shared.o \
	$(EMPTY_MACRO)

${E}/xcombine_sem: $(xcombine_sem_OBJECTS) $(xcombine_sem_SHARED_OBJECTS) $(COND_MPI_OBJECTS)
	@echo ""
	@echo "building xcombine_sem"
	@echo ""
	${FCLINK} -o $@ $+ $(MPILIBS)
	@echo ""


##
## xsmooth_sem
##
xsmooth_sem_OBJECTS = \
	$O/postprocess_par.postprocess_module.o \
	$O/parse_kernel_names.postprocess.o \
	$O/smooth_sem.postprocess.o \
	$(EMPTY_MACRO)

xsmooth_sem_SHARED_OBJECTS = \
	$O/specfem3D_par.spec_module.o \
	$O/pml_par.spec_module.o \
	$O/read_mesh_databases.spec.o \
	$O/read_mesh_databases_hdf5.spec_hdf5.o \
	$O/shared_par.shared_module.o \
	$O/adios_manager.shared_adios_module.o \
	$O/check_mesh_resolution.shared.o \
	$O/count_number_of_sources.shared.o \
	$O/create_name_database.shared.o \
	$O/exit_mpi.shared.o \
	$O/gll_library.shared.o \
	$O/hdf5_manager.shared_hdf5_module.o \
	$O/heap_sort.shared.o \
	$O/param_reader.cc.o \
	$O/read_parameter_file.shared.o \
	$O/read_value_parameters.shared.o \
	$O/search_kdtree.shared.o \
	$O/write_VTK_data.shared.o \
	$(EMPTY_MACRO)


xsmooth_sem_pde_OBJECTS = \
	$O/postprocess_par.postprocess_module.o \
	$O/parse_kernel_names.postprocess.o \
	$O/smooth_sem_pde.postprocess.o \
	$(EMPTY_MACRO)


xsmooth_sem_pde_SHARED_OBJECTS = \
	$O/specfem3D_par.spec_module.o \
	$O/asdf_data.spec_module.o \
	$O/pml_par.spec_module.o \
	$O/calendar.spec.o \
	$O/comp_source_time_function.spec.o \
	$O/compute_add_sources_viscoelastic.spec.o \
	$O/compute_adj_source_frechet.spec.o \
	$O/compute_arrays_source.spec.o \
	$O/compute_element.spec.o \
	$O/compute_element_strain.spec.o \
	$O/compute_gradient_in_acoustic.spec.o \
	$O/compute_interpolated_dva.spec.o \
	$O/compute_seismograms.spec.o \
	$O/get_cmt.spec.o \
	$O/hdf5_io_server.spec_hdf5.o \
	$O/initialize_simulation.spec.o \
	$O/noise_tomography.spec.o \
	$O/read_external_stf.spec.o \
	$O/read_mesh_databases.spec.o \
	$O/read_mesh_databases_hdf5.spec_hdf5.o \
	$O/write_movie_output_HDF5.spec_hdf5.o \
	$O/write_output_HDF5.spec_hdf5.o \
	$O/write_seismograms.spec.o \
	$O/write_output_ASCII_or_binary.spec.o \
	$O/write_output_SU.spec.o \
	$(EMPTY_MACRO)

xsmooth_sem_pde_SHARED_OBJECTS += \
	$O/shared_par.shared_module.o \
	$O/adios_manager.shared_adios_module.o \
	$O/init_openmp.shared.o \
	$O/lagrange_poly.shared.o \
	$O/check_mesh_resolution.shared.o \
	$O/count_number_of_sources.shared.o \
	$O/create_name_database.shared.o \
	$O/exit_mpi.shared.o \
	$O/gll_library.shared.o \
	$O/hdf5_manager.shared_hdf5_module.o \
	$O/heap_sort.shared.o \
	$O/param_reader.cc.o \
	$O/read_parameter_file.shared.o \
	$O/read_value_parameters.shared.o \
	$O/search_kdtree.shared.o \
	$O/write_VTK_data.shared.o \
	$O/define_derivation_matrices.shared.o \
	$O/assemble_MPI_scalar.shared.o \
	$O/netlib_specfun_erf.shared.o \
	$O/write_c_binary.cc.o \
	$(EMPTY_MACRO)

###
### ADIOS
###

# conditional adios linking
ifeq ($(ADIOS),yes)
xsmooth_sem_OBJECTS += $(adios_specfem3D_OBJECTS)
xsmooth_sem_SHARED_OBJECTS += $(adios_specfem3D_PREOBJECTS)
xsmooth_sem_pde_OBJECTS += $(adios_specfem3D_OBJECTS)
xsmooth_sem_pde_SHARED_OBJECTS += $(adios_specfem3D_PREOBJECTS)
else ifeq ($(ADIOS2),yes)
xsmooth_sem_OBJECTS += $(adios_specfem3D_OBJECTS)
xsmooth_sem_SHARED_OBJECTS += $(adios_specfem3D_PREOBJECTS)
xsmooth_sem_pde_OBJECTS += $(adios_specfem3D_OBJECTS)
xsmooth_sem_pde_SHARED_OBJECTS += $(adios_specfem3D_PREOBJECTS)
else
xsmooth_sem_OBJECTS += $(adios_specfem3D_STUBS)
xsmooth_sem_pde_OBJECTS += $(adios_specfem3D_STUBS)
endif

###
### HDF 5
###

ifeq ($(HDF5),yes)
tomography/postprocess_sensitivity_kernels_MODULES += \
	$(FC_MODDIR)/specfem_par_movie_hdf5.$(FC_MODEXT) \
	$(EMPTY_MACRO)
endif

###
### ASDF
###

ifeq ($(ASDF),yes)
XSMOOTH_SEM_PDE_LINK_FLAGS += $(ASDF_LIBS)
xsmooth_sem_pde_OBJECTS += $(asdf_specfem3D_OBJECTS)
xsmooth_sem_pde_SHARED_OBJECTS += $(asdf_specfem3D_SHARED_OBJECTS)
else
xsmooth_sem_pde_OBJECTS += $(asdf_specfem3D_STUBS)
xsmooth_sem_pde_SHARED_OBJECTS += $(asdf_specfem3D_SHARED_STUBS)
endif

###
### GPU
###

xsmooth_sem_LIBS = $(MPILIBS)
xsmooth_sem_OBJECTS += $(gpu_OBJECTS)
xsmooth_sem_pde_LIBS = $(MPILIBS)
xsmooth_sem_pde_OBJECTS += $(gpu_OBJECTS)

## cuda
ifeq ($(HAS_GPU),yes)
xsmooth_sem_LIBS += $(GPU_LINK)
xsmooth_sem_pde_LIBS += $(GPU_LINK)
endif
INFO_SMOOTH="building xsmooth_sem $(BUILD_VERSION_TXT)"
INFO_SMOOTH_PDE="building xsmooth_sem_pde $(BUILD_VERSION_TXT)"

# extra dependencies
$O/smooth_sem.postprocess.o: $O/specfem3D_par.spec_module.o $O/postprocess_par.postprocess_module.o
$O/smooth_sem.postprocess.o: $O/search_kdtree.shared.o

$O/smooth_sem_pde.postprocess.o: $O/specfem3D_par.spec_module.o $O/postprocess_par.postprocess_module.o

${E}/xsmooth_sem: $(xsmooth_sem_OBJECTS) $(xsmooth_sem_SHARED_OBJECTS) $(COND_MPI_OBJECTS)
	@echo ""
	@echo $(INFO_SMOOTH)
	@echo ""
	${FCLINK} -o $@ $+ $(xsmooth_sem_LIBS)
	@echo ""

${E}/xsmooth_sem_pde: $(xsmooth_sem_pde_OBJECTS) $(xsmooth_sem_pde_SHARED_OBJECTS) $(COND_MPI_OBJECTS)
	@echo ""
	@echo $(INFO_SMOOTH_PDE)
	@echo ""
	${FCLINK} -o $@ $+ $(xsmooth_sem_pde_LIBS) $(XSMOOTH_SEM_PDE_LINK_FLAGS)
	@echo

#######################################

###
### Module dependencies
###
$O/postprocess_par.postprocess_module.o: $O/shared_par.shared_module.o

####
#### rule for each .o file below
####

$O/%.postprocess_module.o: $S/%.f90 ${SETUP}/constants_tomography.h $O/shared_par.shared_module.o
	${FCCOMPILE_CHECK} ${FCFLAGS_f90} -c -o $@ $<

$O/%.postprocess.o: $S/%.f90 ${SETUP}/constants_tomography.h $O/postprocess_par.postprocess_module.o
	${FCCOMPILE_CHECK} ${FCFLAGS_f90} -c -o $@ $<

$O/%.postprocess.o: $S/%.F90 ${SETUP}/constants_tomography.h $O/postprocess_par.postprocess_module.o
	${FCCOMPILE_CHECK} ${FCFLAGS_f90} -c -o $@ $<

$O/%.postprocess.o: $S/%.c ${SETUP}/config.h
	${CC} -c $(CPPFLAGS) $(CFLAGS) $(MPI_INCLUDES) -o $@ $<
