/*
 !=====================================================================
 !
 !               S p e c f e m 3 D  V e r s i o n  2 . 0
 !               ---------------------------------------
 !
 !          Main authors: Dimitri Komatitsch and Jeroen Tromp
 !    Princeton University, USA and University of Pau / CNRS / INRIA
 ! (c) Princeton University / California Institute of Technology and University of Pau / CNRS / INRIA
 !                            April 2011
 !
 ! This program is free software; you can redistribute it and/or modify
 ! it under the terms of the GNU General Public License as published by
 ! the Free Software Foundation; either version 2 of the License, or
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

#include <stdio.h>
#include <cuda.h>
#include <cublas.h>

#include "config.h"
#include "mesh_constants_cuda.h"


#define CUBLAS_ERROR(s,n)  if (s != CUBLAS_STATUS_SUCCESS) {  \
fprintf (stderr, "CUBLAS Memory Write Error @ %d\n",n); \
exit(EXIT_FAILURE); }

/* ----------------------------------------------------------------------------------------------- */

// elastic wavefield

/* ----------------------------------------------------------------------------------------------- */


__global__ void UpdateDispVeloc_kernel(realw* displ,
                                       realw* veloc,
                                       realw* accel,
                                       int size,
                                       realw deltat,
                                       realw deltatsqover2,
                                       realw deltatover2) {
  int id = threadIdx.x + blockIdx.x*blockDim.x + blockIdx.y*gridDim.x*blockDim.x;

  /* because of block and grid sizing problems, there is a small */
  /* amount of buffer at the end of the calculation */
  if(id < size) {
    displ[id] = displ[id] + deltat*veloc[id] + deltatsqover2*accel[id];
    veloc[id] = veloc[id] + deltatover2*accel[id];
    accel[id] = 0; // can do this using memset...not sure if faster
  }
}

/* ----------------------------------------------------------------------------------------------- */

extern "C"
void FC_FUNC_(it_update_displacement_cuda,
              IT_UPDATE_DISPLACMENT_CUDA)(long* Mesh_pointer_f,
                                                 int* size_F,
                                                 realw* deltat_F,
                                                 realw* deltatsqover2_F,
                                                 realw* deltatover2_F,
                                                 int* SIMULATION_TYPE,
                                                 realw* b_deltat_F,
                                                 realw* b_deltatsqover2_F,
                                                 realw* b_deltatover2_F) {

TRACE("it_update_displacement_cuda");

  Mesh* mp = (Mesh*)(*Mesh_pointer_f); // get Mesh from fortran integer wrapper

  //int i,device;

  int size = *size_F;
  realw deltat = *deltat_F;
  realw deltatsqover2 = *deltatsqover2_F;
  realw deltatover2 = *deltatover2_F;
  realw b_deltat = *b_deltat_F;
  realw b_deltatsqover2 = *b_deltatsqover2_F;
  realw b_deltatover2 = *b_deltatover2_F;
  //cublasStatus status;

  int blocksize = BLOCKSIZE_KERNEL1;
  int size_padded = ((int)ceil(((double)size)/((double)blocksize)))*blocksize;

  int num_blocks_x = size_padded/blocksize;
  int num_blocks_y = 1;
  while(num_blocks_x > 65535) {
    num_blocks_x = (int) ceil(num_blocks_x*0.5f);
    num_blocks_y = num_blocks_y*2;
  }

  dim3 grid(num_blocks_x,num_blocks_y);
  dim3 threads(blocksize,1,1);


//#ifdef ENABLE_VERY_SLOW_ERROR_CHECKING
//  exit_on_cuda_error("Before UpdateDispVeloc_kernel");
//#endif

  //launch kernel
  UpdateDispVeloc_kernel<<<grid,threads,0,mp->compute_stream>>>(mp->d_displ,mp->d_veloc,mp->d_accel,
                                           size,deltat,deltatsqover2,deltatover2);

  //cudaThreadSynchronize();
//#ifdef ENABLE_VERY_SLOW_ERROR_CHECKING
//  //printf("checking updatedispl_kernel launch...with %dx%d blocks\n",num_blocks_x,num_blocks_y);
//  // sync and check to catch errors from previous async operations
//  exit_on_cuda_error("UpdateDispVeloc_kernel");
//#endif

  // kernel for backward fields
  if(*SIMULATION_TYPE == 3) {

    UpdateDispVeloc_kernel<<<grid,threads,0,mp->compute_stream>>>(mp->d_b_displ,mp->d_b_veloc,mp->d_b_accel,
                                             size,b_deltat,b_deltatsqover2,b_deltatover2);

//#ifdef ENABLE_VERY_SLOW_ERROR_CHECKING
//    //printf("checking updatedispl_kernel launch...with %dx%d blocks\n",num_blocks_x,num_blocks_y);
//    exit_on_cuda_error("after SIM_TYPE==3 UpdateDispVeloc_kernel");
//#endif
  }

#ifdef ENABLE_VERY_SLOW_ERROR_CHECKING
  exit_on_cuda_error("it_update_displacement_cuda");
#endif
}

/* ----------------------------------------------------------------------------------------------- */

// acoustic wavefield

// KERNEL 1
/* ----------------------------------------------------------------------------------------------- */

__global__ void UpdatePotential_kernel(realw* potential_acoustic,
                                       realw* potential_dot_acoustic,
                                       realw* potential_dot_dot_acoustic,
                                       int size,
                                       realw deltat,
                                       realw deltatsqover2,
                                       realw deltatover2) {
  int id = threadIdx.x + blockIdx.x*blockDim.x + blockIdx.y*gridDim.x*blockDim.x;

  /* because of block and grid sizing problems, there is a small */
  /* amount of buffer at the end of the calculation */
  if(id < size) {
    potential_acoustic[id] = potential_acoustic[id]
                            + deltat*potential_dot_acoustic[id]
                            + deltatsqover2*potential_dot_dot_acoustic[id];

    potential_dot_acoustic[id] = potential_dot_acoustic[id]
                                + deltatover2*potential_dot_dot_acoustic[id];

    potential_dot_dot_acoustic[id] = 0;
  }
}

/* ----------------------------------------------------------------------------------------------- */

extern "C"
void FC_FUNC_(it_update_displacement_ac_cuda,
              it_update_displacement_ac_cuda)(long* Mesh_pointer_f,
                                               int* size_F,
                                               realw* deltat_F,
                                               realw* deltatsqover2_F,
                                               realw* deltatover2_F,
                                               int* SIMULATION_TYPE,
                                               realw* b_deltat_F,
                                               realw* b_deltatsqover2_F,
                                               realw* b_deltatover2_F) {
TRACE("it_update_displacement_ac_cuda");
  Mesh* mp = (Mesh*)(*Mesh_pointer_f); // get Mesh from fortran integer wrapper

  //int i,device;
  int size = *size_F;
  realw deltat = *deltat_F;
  realw deltatsqover2 = *deltatsqover2_F;
  realw deltatover2 = *deltatover2_F;
  realw b_deltat = *b_deltat_F;
  realw b_deltatsqover2 = *b_deltatsqover2_F;
  realw b_deltatover2 = *b_deltatover2_F;
  //cublasStatus status;

  int blocksize = BLOCKSIZE_KERNEL1;
  int size_padded = ((int)ceil(((double)size)/((double)blocksize)))*blocksize;

  int num_blocks_x = size_padded/blocksize;
  int num_blocks_y = 1;
  while(num_blocks_x > 65535) {
    num_blocks_x = (int) ceil(num_blocks_x*0.5f);
    num_blocks_y = num_blocks_y*2;
  }

  dim3 grid(num_blocks_x,num_blocks_y);
  dim3 threads(blocksize,1,1);

  //launch kernel
  UpdatePotential_kernel<<<grid,threads,0,mp->compute_stream>>>(mp->d_potential_acoustic,
                                           mp->d_potential_dot_acoustic,
                                           mp->d_potential_dot_dot_acoustic,
                                           size,deltat,deltatsqover2,deltatover2);

  if(*SIMULATION_TYPE == 3) {
    UpdatePotential_kernel<<<grid,threads,0,mp->compute_stream>>>(mp->d_b_potential_acoustic,
                                             mp->d_b_potential_dot_acoustic,
                                             mp->d_b_potential_dot_dot_acoustic,
                                             size,b_deltat,b_deltatsqover2,b_deltatover2);
  }

  //cudaThreadSynchronize();
#ifdef ENABLE_VERY_SLOW_ERROR_CHECKING
  //printf("checking updatedispl_kernel launch...with %dx%d blocks\n",num_blocks_x,num_blocks_y);
  exit_on_cuda_error("it_update_displacement_ac_cuda");
#endif
}
