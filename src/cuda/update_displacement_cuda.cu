/*
 !=====================================================================
 !
 !               S p e c f e m 3 D  V e r s i o n  3 . 0
 !               ---------------------------------------
 !
 !     Main historical authors: Dimitri Komatitsch and Jeroen Tromp
 !                              CNRS, France
 !                       and Princeton University, USA
 !                 (there are currently many more authors!)
 !                           (c) October 2017
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

#include "mesh_constants_cuda.h"

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

  // two dimensional array of blocks on grid where each block has one dimensional array of threads
  int id = threadIdx.x + blockIdx.x*blockDim.x + blockIdx.y*gridDim.x*blockDim.x;
  realw acc = accel[id];
  // because of block and grid sizing problems, there is a small
  // amount of buffer at the end of the calculation
  if (id < size) {
    displ[id] = displ[id] + deltat*veloc[id] + deltatsqover2*acc;
    veloc[id] = veloc[id] + deltatover2*acc;
    accel[id] = 0.0f; // can do this using memset...not sure if faster,probably not
  }

// -----------------
// total of: 6 FLOP per thread (without int id calculation at beginning)
//
//           8 * 4 BYTE = 32 DRAM accesses per thread
//
// arithmetic intensity: 6 FLOP / 32 BYTES ~ 0.19 FLOP/BYTE
// -----------------
// nvprof: 24599250 flops for 4099875 threads -> 6 FLOP per thread
}

/* ----------------------------------------------------------------------------------------------- */

extern "C"
void FC_FUNC_(update_displacement_cuda,
              UPDATE_DISPLACMENT_CUDA)(long* Mesh_pointer,
                                          realw* deltat_F,
                                          realw* deltatsqover2_F,
                                          realw* deltatover2_F,
                                          realw* b_deltat_F,
                                          realw* b_deltatsqover2_F,
                                          realw* b_deltatover2_F) {

  TRACE("\tupdate_displacement_cuda");

  Mesh* mp = (Mesh*)(*Mesh_pointer); // get Mesh from fortran integer wrapper

  realw deltat = *deltat_F;
  realw deltatsqover2 = *deltatsqover2_F;
  realw deltatover2 = *deltatover2_F;

  int size = NDIM * mp->NGLOB_AB;

  int blocksize = BLOCKSIZE_KERNEL1;
  int size_padded = ((int)ceil(((double)size)/((double)blocksize)))*blocksize;

  int num_blocks_x, num_blocks_y;
  get_blocks_xy(size_padded/blocksize,&num_blocks_x,&num_blocks_y);

  dim3 grid(num_blocks_x,num_blocks_y);
  dim3 threads(blocksize,1,1);

  // Cuda timing
  cudaEvent_t start,stop;
  if (CUDA_TIMING_UPDATE ){
    start_timing_cuda(&start,&stop);
  }

  // debug
  //realw max_d,max_v,max_a;
  //max_d = get_device_array_maximum_value(mp->d_displ, size);
  //max_v = get_device_array_maximum_value(mp->d_veloc, size);
  //max_a = get_device_array_maximum_value(mp->d_accel, size);
  //printf("rank %d - max displ: %f veloc: %f accel: %f\n",mp->myrank,max_d,max_v,max_a);

  //launch kernel
  UpdateDispVeloc_kernel<<<grid,threads,0,mp->compute_stream>>>(mp->d_displ,mp->d_veloc,mp->d_accel,
                                                                size,deltat,deltatsqover2,deltatover2);

  // kernel for backward fields
  if (mp->simulation_type == 3) {
    realw b_deltat = *b_deltat_F;
    realw b_deltatsqover2 = *b_deltatsqover2_F;
    realw b_deltatover2 = *b_deltatover2_F;

    UpdateDispVeloc_kernel<<<grid,threads,0,mp->compute_stream>>>(mp->d_b_displ,mp->d_b_veloc,mp->d_b_accel,
                                                                  size,b_deltat,b_deltatsqover2,b_deltatover2);
  }

  // Cuda timing
  if (CUDA_TIMING_UPDATE ){
    realw flops,time;
    stop_timing_cuda(&start,&stop,"UpdateDispVeloc_kernel",&time);
    // time in seconds
    time = time / 1000.;
    // performance: 6 FLOPS per thread
    flops = 6.0 * size;
    //printf("  performance: %f GFlop/s num_blocks x/y: %d %d threads: %d\n", flops/time * 1.e-9,num_blocks_x,num_blocks_y,size);
    printf("  performance: %f GFlop/s\n", flops/time * 1.e-9);
  }

#ifdef ENABLE_VERY_SLOW_ERROR_CHECKING
  exit_on_cuda_error("update_displacement_cuda");
#endif
}

/* ----------------------------------------------------------------------------------------------- */

// acoustic wavefield

// KERNEL 1
/* ----------------------------------------------------------------------------------------------- */

__global__ void UpdatePotential_kernel(field* potential_acoustic,
                                       field* potential_dot_acoustic,
                                       field* potential_dot_dot_acoustic,
                                       int size,
                                       realw deltat,
                                       realw deltatsqover2,
                                       realw deltatover2) {

  int id = threadIdx.x + blockIdx.x*blockDim.x + blockIdx.y*gridDim.x*blockDim.x;

  // because of block and grid sizing problems, there is a small
  // amount of buffer at the end of the calculation
  if (id < size) {
    field p_dot_dot = potential_dot_dot_acoustic[id];
    field p_dot = potential_dot_acoustic[id];
    potential_acoustic[id] +=   deltat*p_dot
                              + deltatsqover2*p_dot_dot;

    potential_dot_acoustic[id] = p_dot + deltatover2*p_dot_dot;

    potential_dot_dot_acoustic[id] = Make_field(0.f);
  }

// -----------------
// total of: 6 FLOP per thread (without id calculation)
//
//           8 * 4 BYTE = 32 DRAM accesses per thread
//
// arithmetic intensity: 6 FLOP / 32 BYTES ~ 0.19 FLOP/BYTE
// -----------------
//
// nvprof: nvprof --metrics flops_sp ./xspecfem3D
//          -> 8199750 FLOPS (Single) floating-point operations for 1366625 threads
//                                    1366625 (NGLOB) -> 10677 * 128 active threads- 31 ghost threads
//          -> 6 FLOP per thread


}

/* ----------------------------------------------------------------------------------------------- */

extern "C"
void FC_FUNC_(it_update_displacement_ac_cuda,
              it_update_displacement_ac_cuda)(long* Mesh_pointer,
                                               realw* deltat_F,
                                               realw* deltatsqover2_F,
                                               realw* deltatover2_F,
                                               realw* b_deltat_F,
                                               realw* b_deltatsqover2_F,
                                               realw* b_deltatover2_F) {
  TRACE("\tit_update_displacement_ac_cuda");
  Mesh* mp = (Mesh*)(*Mesh_pointer); // get Mesh from fortran integer wrapper

  int size = mp->NGLOB_AB;

  int blocksize = BLOCKSIZE_KERNEL1;
  int size_padded = ((int)ceil(((double)size)/((double)blocksize)))*blocksize;

  int num_blocks_x, num_blocks_y;
  get_blocks_xy(size_padded/blocksize,&num_blocks_x,&num_blocks_y);

  dim3 grid(num_blocks_x,num_blocks_y);
  dim3 threads(blocksize,1,1);

  //launch kernel
  // forward wavefields
  realw deltat = *deltat_F;
  realw deltatsqover2 = *deltatsqover2_F;
  realw deltatover2 = *deltatover2_F;

  // Cuda timing
  cudaEvent_t start,stop;
  if (CUDA_TIMING_UPDATE ){
    start_timing_cuda(&start,&stop);
  }

  UpdatePotential_kernel<<<grid,threads,0,mp->compute_stream>>>(mp->d_potential_acoustic,
                                                                 mp->d_potential_dot_acoustic,
                                                                 mp->d_potential_dot_dot_acoustic,
                                                                 size,deltat,deltatsqover2,deltatover2);

  // backward/reconstructed wavefields
  if (mp->simulation_type == 3) {
    realw b_deltat = *b_deltat_F;
    realw b_deltatsqover2 = *b_deltatsqover2_F;
    realw b_deltatover2 = *b_deltatover2_F;

    UpdatePotential_kernel<<<grid,threads,0,mp->compute_stream>>>(mp->d_b_potential_acoustic,
                                                                  mp->d_b_potential_dot_acoustic,
                                                                  mp->d_b_potential_dot_dot_acoustic,
                                                                  size,b_deltat,b_deltatsqover2,b_deltatover2);
  }

  // Cuda timing
  if (CUDA_TIMING_UPDATE ){
    realw flops,time;
    stop_timing_cuda(&start,&stop,"UpdatePotential_kernel",&time);
    // time in seconds
    time = time / 1000.;
    // performance
    // see with: nvprof --metrics flops_sp ./xspecfem3D
    //           -> using 8199750 FLOPS (Single) floating-point operations for 1366625 threads
    //              = 6 FLOPS per thread
    flops = 6.0 * size;
    //printf("  performance: %f GFlop/s num_blocks x/y: %d %d threads: %d\n", flops/time * 1.e-9,num_blocks_x,num_blocks_y,size);
    printf("  performance: %f GFlop/s\n", flops/time * 1.e-9);
  }

#ifdef ENABLE_VERY_SLOW_ERROR_CHECKING
  //printf("checking updatedispl_kernel launch...with %dx%d blocks\n",num_blocks_x,num_blocks_y);
  exit_on_cuda_error("it_update_displacement_ac_cuda");
#endif
}


/* ----------------------------------------------------------------------------------------------- */

// elastic domains

// KERNEL 3

/* ----------------------------------------------------------------------------------------------- */

__global__ void kernel_3_cuda_device(realw* veloc,
                                     realw* accel,
                                     realw* b_veloc,
                                     realw* b_accel,
                                     int size,
                                     int simulation_type,
                                     realw deltatover2,
                                     realw b_deltatover2,
                                     realw* rmassx,
                                     realw* rmassy,
                                     realw* rmassz) {

  int id = threadIdx.x + blockIdx.x*blockDim.x + blockIdx.y*gridDim.x*blockDim.x;
  realw rx,ry,rz;
  realw ax,ay,az;
  // because of block and grid sizing problems, there is a small
  // amount of buffer at the end of the calculation
  if (id < size) {
    rx = rmassx[id];
    ry = rmassy[id];
    rz = rmassz[id];
    ax = accel[3*id  ]*rx;
    ay = accel[3*id+1]*ry;
    az = accel[3*id+2]*rz;

    accel[3*id]   = ax;
    accel[3*id+1] = ay;
    accel[3*id+2] = az;

    veloc[3*id]   += deltatover2*ax;
    veloc[3*id+1] += deltatover2*ay;
    veloc[3*id+2] += deltatover2*az;

    if (simulation_type==3){
      ax = b_accel[3*id  ]*rx;
      ay = b_accel[3*id+1]*ry;
      az = b_accel[3*id+2]*rz;

      b_accel[3*id]   = ax;
      b_accel[3*id+1] = ay;
      b_accel[3*id+2] = az;

      b_veloc[3*id]   += b_deltatover2*ax;
      b_veloc[3*id+1] += b_deltatover2*ay;
      b_veloc[3*id+2] += b_deltatover2*az;

    }

  }

}

/* ----------------------------------------------------------------------------------------------- */

__global__ void kernel_3_accel_cuda_device(realw* accel,
                                           realw* b_accel,
                                           int size,
                                           int simulation_type,
                                           realw* rmassx,
                                           realw* rmassy,
                                           realw* rmassz) {
  int id = threadIdx.x + blockIdx.x*blockDim.x + blockIdx.y*gridDim.x*blockDim.x;

  realw rx,ry,rz;
  realw ax,ay,az;
  // because of block and grid sizing problems, there is a small
  // amount of buffer at the end of the calculation
  if (id < size) {
    rx = rmassx[id];
    ry = rmassy[id];
    rz = rmassz[id];
    ax = accel[3*id  ]*rx;
    ay = accel[3*id+1]*ry;
    az = accel[3*id+2]*rz;
    accel[3*id  ] = ax;
    accel[3*id+1] = ay;
    accel[3*id+2] = az;

    if (simulation_type==3){
      ax = b_accel[3*id  ]*rx;
      ay = b_accel[3*id+1]*ry;
      az = b_accel[3*id+2]*rz;

      b_accel[3*id]   = ax;
      b_accel[3*id+1] = ay;
      b_accel[3*id+2] = az;

    }

  }
}

/* ----------------------------------------------------------------------------------------------- */

__global__ void kernel_3_veloc_cuda_device(realw* veloc,
                                           realw* accel,
                                           int size,
                                           realw deltatover2) {

  int id = threadIdx.x + blockIdx.x*blockDim.x + blockIdx.y*gridDim.x*blockDim.x;

  // because of block and grid sizing problems, there is a small
  // amount of buffer at the end of the calculation
  if (id < size) {
    veloc[3*id] = veloc[3*id] + deltatover2*accel[3*id];
    veloc[3*id+1] = veloc[3*id+1] + deltatover2*accel[3*id+1];
    veloc[3*id+2] = veloc[3*id+2] + deltatover2*accel[3*id+2];
  }
}

/* ----------------------------------------------------------------------------------------------- */

extern "C"
void FC_FUNC_(kernel_3_a_cuda,
              KERNEL_3_A_CUDA)(long* Mesh_pointer,
                               realw* deltatover2_F,
                               realw* b_deltatover2_F,
                               int* APPROXIMATE_OCEAN_LOAD) {

  TRACE("\tkernel_3_a_cuda");

  Mesh* mp = (Mesh*)(*Mesh_pointer); // get Mesh from fortran integer wrapper

  int size = mp->NGLOB_AB;

  int blocksize = BLOCKSIZE_KERNEL3;
  int size_padded = ((int)ceil(((double)size)/((double)blocksize)))*blocksize;

  int num_blocks_x, num_blocks_y;
  get_blocks_xy(size_padded/blocksize,&num_blocks_x,&num_blocks_y);

  dim3 grid(num_blocks_x,num_blocks_y);
  dim3 threads(blocksize,1,1);

  // check whether we can update accel and veloc, or only accel at this point
  if (*APPROXIMATE_OCEAN_LOAD == 0){
   realw deltatover2 = *deltatover2_F;
   realw b_deltatover2 = *b_deltatover2_F;
   // updates both, accel and veloc
   kernel_3_cuda_device<<< grid, threads,0,mp->compute_stream>>>(mp->d_veloc,
                                                                 mp->d_accel,
                                                                 mp->d_b_veloc,
                                                                 mp->d_b_accel,
                                                                 size,mp->simulation_type,deltatover2,b_deltatover2,
                                                                 mp->d_rmassx,mp->d_rmassy,mp->d_rmassz);
  }else{
   // updates only accel
   kernel_3_accel_cuda_device<<< grid, threads,0,mp->compute_stream>>>(mp->d_accel,
                                                                       mp->d_b_accel,
                                                                       size,
                                                                       mp->simulation_type,
                                                                       mp->d_rmassx,
                                                                       mp->d_rmassy,
                                                                       mp->d_rmassz);
  }

#ifdef ENABLE_VERY_SLOW_ERROR_CHECKING
  //printf("checking updatedispl_kernel launch...with %dx%d blocks\n",num_blocks_x,num_blocks_y);
  exit_on_cuda_error("after kernel 3 a");
#endif
}

/* ----------------------------------------------------------------------------------------------- */

extern "C"
void FC_FUNC_(kernel_3_b_cuda,
              KERNEL_3_B_CUDA)(long* Mesh_pointer,
                               realw* deltatover2_F,
                               realw* b_deltatover2_F) {
  TRACE("\tkernel_3_b_cuda");

  Mesh* mp = (Mesh*)(*Mesh_pointer); // get Mesh from fortran integer wrapper

  int size = mp->NGLOB_AB;

  int blocksize = BLOCKSIZE_KERNEL3;
  int size_padded = ((int)ceil(((double)size)/((double)blocksize)))*blocksize;

  int num_blocks_x, num_blocks_y;
  get_blocks_xy(size_padded/blocksize,&num_blocks_x,&num_blocks_y);

  dim3 grid(num_blocks_x,num_blocks_y);
  dim3 threads(blocksize,1,1);

  realw deltatover2 = *deltatover2_F;
  // updates only veloc at this point
  kernel_3_veloc_cuda_device<<< grid, threads,0,mp->compute_stream>>>(mp->d_veloc,
                                                                      mp->d_accel,
                                                                      size,deltatover2);

  if (mp->simulation_type == 3) {
    realw b_deltatover2 = *b_deltatover2_F;
    kernel_3_veloc_cuda_device<<< grid, threads,0,mp->compute_stream>>>(mp->d_b_veloc,
                                                                        mp->d_b_accel,
                                                                        size,b_deltatover2);
  }

#ifdef ENABLE_VERY_SLOW_ERROR_CHECKING
  //printf("checking updatedispl_kernel launch...with %dx%d blocks\n",num_blocks_x,num_blocks_y);
  exit_on_cuda_error("after kernel 3 b");
#endif
}


/* ----------------------------------------------------------------------------------------------- */

// acoustic domains

// KERNEL 3

/* ----------------------------------------------------------------------------------------------- */


__global__ void kernel_3_acoustic_cuda_device(field* potential_dot_acoustic,
                                                field* potential_dot_dot_acoustic,
                                                field* b_potential_dot_acoustic,
                                                field* b_potential_dot_dot_acoustic,
                                                int simulation_type,
                                                int size,
                                                realw deltatover2,
                                                realw b_deltatover2,
                                                realw* rmass_acoustic) {

  int id = threadIdx.x + blockIdx.x*blockDim.x + blockIdx.y*gridDim.x*blockDim.x;
  realw rmass;
  field p_dot_dot;
  // because of block and grid sizing problems, there is a small
  // amount of buffer at the end of the calculation
  if (id < size) {
    rmass = rmass_acoustic[id];
    // multiplies pressure with the inverse of the mass matrix
    p_dot_dot = rmass*potential_dot_dot_acoustic[id];
    potential_dot_dot_acoustic[id] = p_dot_dot;
    potential_dot_acoustic[id] += deltatover2*p_dot_dot;
    if (simulation_type==3) {
      p_dot_dot = rmass*b_potential_dot_dot_acoustic[id];
      b_potential_dot_dot_acoustic[id] = p_dot_dot;
      b_potential_dot_acoustic[id] += b_deltatover2*p_dot_dot;
    }
  }
}

/* ----------------------------------------------------------------------------------------------- */

extern "C"
void FC_FUNC_(kernel_3_acoustic_cuda,
              KERNEL_3_ACOUSTIC_CUDA)(long* Mesh_pointer,
                                      realw* deltatover2_F,
                                      realw* b_deltatover2_F) {

TRACE("kernel_3_acoustic_cuda");

  Mesh* mp = (Mesh*)(*Mesh_pointer); // get Mesh from fortran integer wrapper

  int size = mp->NGLOB_AB;

  int blocksize = BLOCKSIZE_KERNEL3;
  int size_padded = ((int)ceil(((double)size)/((double)blocksize)))*blocksize;

  int num_blocks_x, num_blocks_y;
  get_blocks_xy(size_padded/blocksize,&num_blocks_x,&num_blocks_y);

  dim3 grid(num_blocks_x,num_blocks_y);
  dim3 threads(blocksize,1,1);

  realw deltaover2 = *deltatover2_F;
  realw b_deltaover2 = *b_deltatover2_F;

  kernel_3_acoustic_cuda_device<<< grid, threads>>>(mp->d_potential_dot_acoustic,
                                                    mp->d_potential_dot_dot_acoustic,
                                                    mp->d_b_potential_dot_acoustic,
                                                    mp->d_b_potential_dot_dot_acoustic,
                                                    mp->simulation_type,
                                                    size,
                                                    deltaover2,
                                                    b_deltaover2,
                                                    mp->d_rmass_acoustic);

#ifdef ENABLE_VERY_SLOW_ERROR_CHECKING
  //printf("checking updatedispl_kernel launch...with %dx%d blocks\n",num_blocks_x,num_blocks_y);
  exit_on_cuda_error("after kernel 3 ");
#endif
}

