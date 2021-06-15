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
	$O/initialize_gpu.o \
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
endif

ifeq ($(HIP),yes)
  # defines $(cuda_kernels_OBJS)
  include $(KERNEL_DIR)/kernel_cuda.mk
  # renames endings
  cuda_kernels_OBJS:=$(subst .cuda-kernel.o,.hip-kernel.o,${cuda_kernels_OBJS})
endif


ifdef NO_GPU
gpu_OBJECTS = $(gpu_specfem3D_STUBS)
else
gpu_OBJECTS = $(gpu_specfem3D_OBJECTS)
endif


#######################################

# substitutes object endings to assign corresponding compilation rule
ifeq ($(HAS_GPU),yes)
	# CUDA kernels only
	ifeq ($(CUDA),yes)
		gpu_specfem3D_OBJECTS:=$(subst .o,.cuda.o,${gpu_specfem3D_OBJECTS})
	endif

	# HIP kernels only
	ifeq ($(HIP), yes)
		gpu_specfem3D_OBJECTS:=$(subst .o,.hip.o,${gpu_specfem3D_OBJECTS})
	endif
endif

gpu_specfem3D_OBJECTS += $(cuda_specfem3D_DEVICE_OBJ) $(cuda_kernels_OBJS)

###
### variables
###

BUILD_VERSION_TXT := with
SELECTOR_CFLAG :=

## CUDA compilation
NVCC_CFLAGS := ${NVCC_FLAGS} -x cu
ifeq ($(CUDA),yes)
  BUILD_VERSION_TXT += Cuda
  SELECTOR_CFLAG += $(FC_DEFINE)USE_CUDA
  GPU_LINK = $(CUDA_LINK)

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

## HIP compilation
HIPCC_CFLAGS := ${HIP_CFLAGS}      # ${HIP_CFLAG_ENDING} adds -x hip depending on platform
ifeq ($(HIP), yes)
  BUILD_VERSION_TXT += HIP
  SELECTOR_CFLAG += $(FC_DEFINE)USE_HIP
  ifneq ($(strip $(HIP_GPU_FLAGS)),)
    SELECTOR_CFLAG += -DHIP_GPU_CFLAGS="$(HIP_GPU_FLAGS)"
  endif
  GPU_LINK = $(HIP_LINK)

  # todo: compile hip with nvcc
  #ifeq ($(CUDA),yes)
  #  GPU_LINK += $(HIP_LINK)
  #  NVCC_CFLAGS += $(HIP_CPU_FLAGS)
  #endif
endif

BUILD_VERSION_TXT += support

#######################################

####
#### rule for each .o file below
####

###
### CUDA compilation
###

# source kernel files
ifeq ($(CUDA),yes)
$O/%.cuda-kernel.o: $(KERNEL_DIR)/%.cu $S/mesh_constants_gpu.h $S/mesh_constants_cuda.h $(KERNEL_DIR)/kernel_proto.cu.h
	$(NVCC) -c $< -o $@ $(NVCC_CFLAGS) -I${SETUP} -I$(KERNEL_DIR) $(SELECTOR_CFLAG) -include $(word 2,$^)

$(cuda_specfem3D_DEVICE_OBJ): $(subst $(cuda_specfem3D_DEVICE_OBJ), ,$(gpu_specfem3D_OBJECTS)) $(cuda_kernels_OBJS)
	${NVCCLINK} -o $@ $^
endif

ifeq ($(HIP),yes)
$O/%.hip-kernel.o: $(KERNEL_DIR)/%.cpp $S/mesh_constants_gpu.h $S/mesh_constants_hip.h $(KERNEL_DIR)/kernel_proto.cu.h
	$(HIPCC) -c $< -o $@ $(HIP_CFLAGS) -I${SETUP} -I$(KERNEL_DIR) $(SELECTOR_CFLAG) -include $(word 2,$^)
endif


# source files in src/gpu/
$O/%.cuda.o: $S/%.cu ${SETUP}/config.h $S/mesh_constants_gpu.h $S/mesh_constants_cuda.h $S/prepare_constants_cuda.h
	${NVCC} -c $< -o $@ $(NVCC_FLAGS) -I${SETUP} -I$(KERNEL_DIR) $(SELECTOR_CFLAG)

$O/%.cuda.o: $S/%.c ${SETUP}/config.h $S/mesh_constants_gpu.h $S/mesh_constants_cuda.h
	$(NVCC) -c $< -o $@ $(NVCC_CFLAGS) -I${SETUP} -I$(KERNEL_DIR) $(SELECTOR_CFLAG)


$O/%.hip.o: $S/%.cu ${SETUP}/config.h $S/mesh_constants_gpu.h  $S/mesh_constants_hip.h $S/prepare_constants_cuda.h
	${HIPCC} ${HIP_CFLAG_ENDING} -c $< -o $@ $(HIPCC_CFLAGS) -I${SETUP} -I$(KERNEL_DIR) $(SELECTOR_CFLAG)

$O/%.hip.o: $S/%.c ${SETUP}/config.h $S/mesh_constants_gpu.h  $S/mesh_constants_hip.h
	${HIPCC} ${HIP_CFLAG_ENDING} -c $< -o $@ $(HIPCC_CFLAGS) -I${SETUP} -I$(KERNEL_DIR) $(SELECTOR_CFLAG)


# C version
$O/%.gpu_cc.o: $S/%.c ${SETUP}/config.h
	${CC} -c $(CPPFLAGS) $(CFLAGS) $(MPI_INCLUDES) -o $@ $<

print-%:
	@echo '$*=$($*)'
