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
S := ${S_TOP}/src/gpu

KERNEL_DIR_NAME := kernels
KERNEL_DIR := ${S}/${KERNEL_DIR_NAME}

#######################################

gpu_specfem3D_TARGETS = \
	$(gpu_specfem3D_OBJECTS) \
	$(EMPTY_MACRO)

gpu_specfem3D_OBJECTS = \
	$O/assemble_MPI_scalar_cuda.o \
	$O/assemble_MPI_vector_cuda.o \
	$O/check_fields_cuda.o \
	$O/compute_add_sources_acoustic_cuda.o \
	$O/compute_add_sources_viscoelastic_cuda.o \
	$O/compute_coupling_cuda.o \
	$O/compute_forces_acoustic_cuda.o \
	$O/compute_forces_viscoelastic_cuda.o \
	$O/compute_kernels_cuda.o \
	$O/compute_stacey_acoustic_cuda.o \
	$O/compute_stacey_viscoelastic_cuda.o \
	$O/fault_solver_dynamics.o \
	$O/helper_functions.o \
	$O/initialize_cuda.o \
	$O/noise_tomography_cuda.o \
	$O/prepare_mesh_constants_cuda.o \
	$O/save_and_compare_cpu_vs_gpu.o \
	$O/smooth_cuda.o \
	$O/transfer_fields_cuda.o \
	$O/update_displacement_cuda.o \
	$O/write_seismograms_cuda.o \
	$(EMPTY_MACRO)

gpu_specfem3D_STUBS = \
	$O/specfem3D_gpu_cuda_method_stubs.gpu_cc.o \
	$(EMPTY_MACRO)


# CUDA
ifeq ($(CUDA),yes)
  cuda_specfem3D_DEVICE_OBJ =  $O/cuda_device_obj.o

  # defines $(cuda_kernels_OBJS)
  include $(KERNEL_DIR)/kernel_cuda.mk

	# replaces .o endings with .cuda.o for Cuda object files
	gpu_specfem3D_OBJECTS:=$(subst .o,.cuda.o,${gpu_specfem3D_OBJECTS})
endif

gpu_specfem3D_OBJECTS += $(cuda_specfem3D_DEVICE_OBJ) $(cuda_kernels_OBJS)

###
### variables
###

NVCC_CFLAGS := ${NVCC_FLAGS} -x cu

BUILD_VERSION_TXT := with
SELECTOR_CFLAG :=

ifeq ($(CUDA),yes)
  BUILD_VERSION_TXT += Cuda
  SELECTOR_CFLAG += $(FC_DEFINE)USE_CUDA

  ifeq ($(CUDA5),yes)
    BUILD_VERSION_TXT += (v5)
  endif
  ifeq ($(CUDA6),yes)
    BUILD_VERSION_TXT += (v6)
  endif
  ifeq ($(CUDA7),yes)
    BUILD_VERSION_TXT += (v7)
  endif
  ifeq ($(CUDA8),yes)
    BUILD_VERSION_TXT += (v8)
  endif
  ifeq ($(CUDA9),yes)
    BUILD_VERSION_TXT += (v9)
  endif
  ifeq ($(CUDA10),yes)
    BUILD_VERSION_TXT += (v10)
  endif
  ifeq ($(CUDA11),yes)
    BUILD_VERSION_TXT += (v11)
  endif
endif

BUILD_VERSION_TXT += support

#######################################

####
#### rule for each .o file below
####

###
### CUDA compilation
###

ifeq ($(CUDA),yes)
$O/%.cuda-kernel.o: $(KERNEL_DIR)/%.cu $S/mesh_constants_cuda.h $(KERNEL_DIR)/kernel_proto.cu.h #$S/mesh_constants_gpu.h 
	$(NVCC) -c $< -o $@ $(NVCC_CFLAGS) -I${SETUP} -I$(KERNEL_DIR) $(SELECTOR_CFLAG) -include $(word 2,$^)

$(cuda_specfem3D_DEVICE_OBJ): $(subst $(cuda_specfem3D_DEVICE_OBJ), ,$(gpu_specfem3D_OBJECTS)) $(cuda_kernels_OBJS)
	${NVCCLINK} -o $@ $^
endif

#print-%:
#	@echo '$*=$($*)'

$O/%.cuda.o: $S/%.cu ${SETUP}/config.h $S/mesh_constants_cuda.h $S/prepare_constants_cuda.h
	${NVCC} -c $< -o $@ $(NVCC_FLAGS) -I${SETUP} -I$(KERNEL_DIR) $(SELECTOR_CFLAG)

$O/%.cuda.o: $S/%.c ${SETUP}/config.h $S/mesh_constants_cuda.h
	$(NVCC) -c $< -o $@ $(NVCC_CFLAGS) -I${SETUP} -I$(KERNEL_DIR) $(SELECTOR_CFLAG)

$O/%.gpu_cc.o: $S/%.c ${SETUP}/config.h
	${CC} -c $(CPPFLAGS) $(CFLAGS) $(MPI_INCLUDES) -o $@ $<
