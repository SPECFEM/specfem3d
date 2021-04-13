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
S := ${S_TOP}/src/tomography/postprocess_sensitivity_kernels
$(tomography/postprocess_sensitivity_kernels_OBJECTS): S := ${S_TOP}/src/tomography/postprocess_sensitivity_kernels

#######################################

tomography/postprocess_sensitivity_kernels_TARGETS = \
	$E/xclip_sem \
	$E/xcombine_sem \
	$E/xsmooth_sem \
	$(EMPTY_MACRO)

tomography/postprocess_sensitivity_kernels_OBJECTS = \
	$(xclip_sem_OBJECTS) \
	$(xcombine_sem_OBJECTS) \
	$(xsmooth_sem_OBJECTS) \
	$(EMPTY_MACRO)

tomography/postprocess_sensitivity_kernels_SHARED_OBJECTS = \
	$(xclip_sem_SHARED_OBJECTS) \
	$(xcombine_sem_SHARED_OBJECTS) \
	$(xsmooth_sem_SHARED_OBJECTS) \
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

###
### GPU
###

xsmooth_sem_LIBS = $(MPILIBS)
xsmooth_sem_OBJECTS += $(gpu_OBJECTS)

## cuda
ifeq ($(HAS_GPU),yes)
xsmooth_sem_LIBS += $(GPU_LINK)
endif
INFO_SMOOTH="building xsmooth_sem $(BUILD_VERSION_TXT)"

# extra dependencies
$O/smooth_sem.postprocess.o: $O/specfem3D_par.spec_module.o $O/postprocess_par.postprocess_module.o

${E}/xsmooth_sem: $(xsmooth_sem_OBJECTS) $(xsmooth_sem_SHARED_OBJECTS) $(COND_MPI_OBJECTS)
	@echo ""
	@echo $(INFO_SMOOTH)
	@echo ""
	${FCLINK} -o $@ $+ $(xsmooth_sem_LIBS)
	@echo ""

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
