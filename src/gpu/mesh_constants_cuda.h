/*
 !=====================================================================
 !
 !               S p e c f e m 3 D  V e r s i o n  3 . 0
 !               ---------------------------------------
 !
 !     Main historical authors: Dimitri Komatitsch and Jeroen Tromp
 !                        Princeton University, USA
 !                and CNRS / University of Marseille, France
 !                 (there are currently many more authors!)
 ! (c) Princeton University and CNRS / University of Marseille, July 2012
 !
 ! This program is free software; you can redistribute it and/or modify
 ! it under the terms of the GNU General Public License as published by
 ! the Free Software Foundation; either version 3 of the License, or
 ! (at your option) any later version.
 !
 ! This program is distributed in the hope that it will be useful,
 ! but WITHOUT ANY WARRANTY; without even the implied warranty of
 ! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 ! GNU General Public License for more details.
 !
 ! You should have received a copy of the GNU General Public License along
 ! with this program; if not, write to the Free Software Foundation, Inc.,
 ! 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 !
 !=====================================================================
 */

#ifndef MESH_CONSTANTS_CUDA_H
#define MESH_CONSTANTS_CUDA_H


// CUDA specifics

// (optional) use launch_bounds specification to increase compiler optimization
// (depending on GPU type, register spilling might slow down the performance)
// elastic kernel
#ifdef GPU_DEVICE_K20
// specifics see: https://docs.nvidia.com/cuda/kepler-tuning-guide/index.html
// maximum shared memory per thread block 48KB
//
// note: main kernel is Kernel_2_***_impl() which is limited by shared memory usage to 8 active blocks
//       while register usage might use up to 9 blocks
//
//
// performance statistics: kernel Kernel_2_noatt_impl():
//       shared memory per block = 1700    for Kepler: total = 49152 -> limits active blocks to 16
//       registers per thread    = 48
//       registers per block     = 6144                total = 65536 -> limits active blocks to 10
//
// performance statistics: kernel Kernel_2_att_impl():
//       shared memory per block = 6100    for Kepler: total = 49152 -> limits active blocks to 8
//       registers per thread    = 59
//       registers per block     = 8192                total = 65536 -> limits active blocks to 8
//
// (uncomment if desired)
//#define USE_LAUNCH_BOUNDS
#define LAUNCH_MIN_BLOCKS 10
//#pragma message ("\nCompiling with: USE_LAUNCH_BOUNDS enabled for K20\n")
#endif

// add more card specific values
#ifdef GPU_DEVICE_Maxwell
// specifics see: https://docs.nvidia.com/cuda/maxwell-tuning-guide/index.html
// register file size 64k 32-bit registers per SM
// shared memory 64KB for GM107 and 96KB for GM204
#undef USE_LAUNCH_BOUNDS
#endif

#ifdef GPU_DEVICE_Pascal
// specifics see: https://docs.nvidia.com/cuda/pascal-tuning-guide/index.html
// register file size 64k 32-bit registers per SM
// shared memory 64KB for GP100 and 96KB for GP104
#undef USE_LAUNCH_BOUNDS
#endif

#ifdef GPU_DEVICE_Volta
// specifics see: https://docs.nvidia.com/cuda/volta-tuning-guide/index.html
// register file size 64k 32-bit registers per SM
// shared memory size 96KB per SM (maximum shared memory per thread block)
// maximum registers 255 per thread
#undef USE_LAUNCH_BOUNDS
#endif

#ifdef GPU_DEVICE_Turing
// specifics see: https://docs.nvidia.com/cuda/turing-tuning-guide/index.html
// register file size 64k 32-bit registers per SM
// shared memory size 64KB per SM (maximum shared memory per thread block)
// maximum registers 255 per thread
#undef USE_LAUNCH_BOUNDS
#endif

#ifdef GPU_DEVICE_Ampere
// specifics see: https://docs.nvidia.com/cuda/ampere-tuning-guide/index.html
// register file size 64k 32-bit registers per SM
// shared memory size 164KB per SM (maximum shared memory, 163KB per thread block)
// maximum registers 255 per thread
#undef USE_LAUNCH_BOUNDS
#endif


/* ----------------------------------------------------------------------------------------------- */

// CUDA specifics

#ifdef USE_CUDA
// definitions
typedef cudaEvent_t gpu_event;
typedef cudaStream_t gpu_stream;

// cuda header files
#include "kernel_proto.cu.h"

static inline void print_CUDA_error_if_any(cudaError_t err, int num) {
  if (cudaSuccess != err)
  {
    printf("\nCUDA error !!!!! <%s> !!!!! \nat CUDA call error code: # %d\n",cudaGetErrorString(err),num);
    fflush(stdout);

    // outputs error file
    FILE* fp;
    int myrank;
    char filename[BUFSIZ];
#ifdef WITH_MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
#else
    myrank = 0;
#endif
    sprintf(filename,OUTPUT_FILES"/error_message_%06d.txt",myrank);
    fp = fopen(filename,"a+");
    if (fp != NULL){
      fprintf(fp,"\nCUDA error !!!!! <%s> !!!!! \nat CUDA call error code: # %d\n",cudaGetErrorString(err),num);
      fclose(fp);
    }

    // stops program
#ifdef WITH_MPI
    MPI_Abort(MPI_COMM_WORLD,1);
#endif
    exit(EXIT_FAILURE);
  }
}

#endif  // USE_CUDA

#endif  // MESH_CONSTANTS_CUDA_H
