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

#include <sys/time.h>
#include <sys/resource.h>

#include "config.h"
#include "mesh_constants_cuda.h"
// #include "epik_user.h"

//  cuda constant arrays
__device__ realw d_hprime_xx[NGLL2];

// only needed if NGLLX != NGLLY != NGLLZ
// __device__ realw d_hprime_yy[NGLL2];
// __device__ realw d_hprime_zz[NGLL2];
__device__ realw d_hprimewgll_xx[NGLL2];
__device__ realw d_hprimewgll_yy[NGLL2];
__device__ realw d_hprimewgll_zz[NGLL2];
__device__ realw d_wgllwgll_xy[NGLL2];
__device__ realw d_wgllwgll_xz[NGLL2];
__device__ realw d_wgllwgll_yz[NGLL2];

__constant__ realw d_wgll_cube[NGLL3]; // needed only for gravity case

//daniel: todo - check if necessary...
// prototype for the fortran function to do non-blocking mpi send
//extern "C"
//void assemble_mpi_vector_send_cuda_(void*,void*,void*,void*,void*,void*,void*,void*,void*); // {};


/* ----------------------------------------------------------------------------------------------- */

// prepares a device array with with all inter-element edge-nodes -- this
// is followed by a memcpy and MPI operations

__global__ void prepare_boundary_accel_on_device(realw* d_accel, realw* d_send_accel_buffer,
                                                 int num_interfaces_ext_mesh,
                                                 int max_nibool_interfaces_ext_mesh,
                                                 int* d_nibool_interfaces_ext_mesh,
                                                 int* d_ibool_interfaces_ext_mesh) {

  int id = threadIdx.x + blockIdx.x*blockDim.x + blockIdx.y*gridDim.x*blockDim.x;
  //int iinterface=0;

  for( int iinterface=0; iinterface < num_interfaces_ext_mesh; iinterface++) {
    if(id<d_nibool_interfaces_ext_mesh[iinterface]) {
      d_send_accel_buffer[3*(id + max_nibool_interfaces_ext_mesh*iinterface)] =
        d_accel[3*(d_ibool_interfaces_ext_mesh[id+max_nibool_interfaces_ext_mesh*iinterface]-1)];
      d_send_accel_buffer[3*(id + max_nibool_interfaces_ext_mesh*iinterface)+1] =
        d_accel[3*(d_ibool_interfaces_ext_mesh[id+max_nibool_interfaces_ext_mesh*iinterface]-1)+1];
      d_send_accel_buffer[3*(id + max_nibool_interfaces_ext_mesh*iinterface)+2] =
        d_accel[3*(d_ibool_interfaces_ext_mesh[id+max_nibool_interfaces_ext_mesh*iinterface]-1)+2];
    }
  }

}

/* ----------------------------------------------------------------------------------------------- */

// prepares and transfers the inter-element edge-nodes to the host to be MPI'd
// (elements on boundary)
extern "C"
void FC_FUNC_(transfer_boun_accel_from_device,
              TRANSFER_BOUN_ACCEL_FROM_DEVICE)(int* size, long* Mesh_pointer_f, realw* accel,
                                                    realw* send_accel_buffer,
                                                    int* num_interfaces_ext_mesh,
                                                    int* max_nibool_interfaces_ext_mesh,
                                                    int* nibool_interfaces_ext_mesh,
                                                    int* ibool_interfaces_ext_mesh,
                                                    int* FORWARD_OR_ADJOINT){
TRACE("transfer_boun_accel_from_device");

  Mesh* mp = (Mesh*)(*Mesh_pointer_f); //get mesh pointer out of fortran integer container

  if( *num_interfaces_ext_mesh == 0 ) return;

  int blocksize = BLOCKSIZE_TRANSFER;
  int size_padded = ((int)ceil(((double)mp->max_nibool_interfaces_ext_mesh)/((double)blocksize)))*blocksize;
  int num_blocks_x = size_padded/blocksize;
  int num_blocks_y = 1;
  while(num_blocks_x > 65535) {
    num_blocks_x = (int) ceil(num_blocks_x*0.5f);
    num_blocks_y = num_blocks_y*2;
  }

  dim3 grid(num_blocks_x,num_blocks_y);
  dim3 threads(blocksize,1,1);

  //timing for memory xfer
  // cudaEvent_t start, stop;
  // realw time;
  // cudaEventCreate(&start);
  // cudaEventCreate(&stop);
  // cudaEventRecord( start, 0 );
  if(*FORWARD_OR_ADJOINT == 1) {
    prepare_boundary_accel_on_device<<<grid,threads,0,mp->compute_stream>>>(mp->d_accel,mp->d_send_accel_buffer,
                                                       mp->num_interfaces_ext_mesh,
                                                       mp->max_nibool_interfaces_ext_mesh,
                                                       mp->d_nibool_interfaces_ext_mesh,
                                                       mp->d_ibool_interfaces_ext_mesh);
  }
  else if(*FORWARD_OR_ADJOINT == 3) {
    prepare_boundary_accel_on_device<<<grid,threads,0,mp->compute_stream>>>(mp->d_b_accel,mp->d_send_accel_buffer,
                                                       mp->num_interfaces_ext_mesh,
                                                       mp->max_nibool_interfaces_ext_mesh,
                                                       mp->d_nibool_interfaces_ext_mesh,
                                                       mp->d_ibool_interfaces_ext_mesh);
  }


  cudaMemcpy(send_accel_buffer,mp->d_send_accel_buffer,
             3*mp->max_nibool_interfaces_ext_mesh*mp->num_interfaces_ext_mesh*sizeof(realw),
             cudaMemcpyDeviceToHost);

  // finish timing of kernel+memcpy
  // cudaEventRecord( stop, 0 );
  // cudaEventSynchronize( stop );
  // cudaEventElapsedTime( &time, start, stop );
  // cudaEventDestroy( start );
  // cudaEventDestroy( stop );
  // printf("boundary xfer d->h Time: %f ms\n",time);
#ifdef ENABLE_VERY_SLOW_ERROR_CHECKING
  exit_on_cuda_error("transfer_boun_accel_from_device");
#endif
}

/* ----------------------------------------------------------------------------------------------- */

extern "C"
void FC_FUNC_(transfer_boundary_from_device_a,
              TRANSFER_BOUNDARY_FROM_DEVICE_A)(long* Mesh_pointer,
                                               int* nspec_outer_elastic) {

// asynchronous transfer from device to host

  TRACE("transfer_boundary_from_device_a");

  Mesh* mp = (Mesh*)(*Mesh_pointer); // get Mesh from fortran integer wrapper

//daniel: todo - check below with this...
 int blocksize = BLOCKSIZE_TRANSFER;
 int size_padded = ((int)ceil(((double)mp->max_nibool_interfaces_ext_mesh)/((double)blocksize)))*blocksize;
 int num_blocks_x = size_padded/blocksize;
 int num_blocks_y = 1;
 while(num_blocks_x > 65535) {
 num_blocks_x = (int) ceil(num_blocks_x*0.5f);
 num_blocks_y = num_blocks_y*2;
 }
 dim3 grid(num_blocks_x,num_blocks_y);
 dim3 threads(blocksize,1,1);

/*
//daniel: todo - check originally used...
  int num_blocks_x = *nspec_outer_elastic;
  int num_blocks_y = 1;
  while(num_blocks_x > 65535) {
    num_blocks_x = (int) ceil(num_blocks_x*0.5f);
    num_blocks_y = num_blocks_y*2;
  }

  int blocksize = NGLL3_PADDED;
  dim3 grid(num_blocks_x,num_blocks_y);
  dim3 threads(blocksize,1,1);
*/

  prepare_boundary_accel_on_device<<<grid,threads,0,mp->compute_stream>>>(mp->d_accel,mp->d_send_accel_buffer,
                                                                          mp->num_interfaces_ext_mesh,
                                                                          mp->max_nibool_interfaces_ext_mesh,
                                                                          mp->d_nibool_interfaces_ext_mesh,
                                                                          mp->d_ibool_interfaces_ext_mesh);
  // wait until kernel is finished before starting async memcpy
#if CUDA_VERSION >= 4000
  cudaDeviceSynchronize();
#else
  cudaThreadSynchronize();
#endif

  cudaMemcpyAsync(mp->h_send_accel_buffer,mp->d_send_accel_buffer,
                  3* mp->max_nibool_interfaces_ext_mesh* mp->num_interfaces_ext_mesh*sizeof(realw),
                  cudaMemcpyDeviceToHost,mp->copy_stream);
  // cudaMemcpyAsync(mp->h_send_accel_buffer,mp->d_send_accel_buffer,
  // 3* mp->max_nibool_interfaces_ext_mesh* mp->num_interfaces_ext_mesh*sizeof(realw),
  // cudaMemcpyDeviceToHost,mp->compute_stream);

}

/* ----------------------------------------------------------------------------------------------- */

__global__ void assemble_boundary_accel_on_device(realw* d_accel, realw* d_send_accel_buffer,
                                                  int num_interfaces_ext_mesh,
                                                  int max_nibool_interfaces_ext_mesh,
                                                  int* d_nibool_interfaces_ext_mesh,
                                                  int* d_ibool_interfaces_ext_mesh) {

  int id = threadIdx.x + blockIdx.x*blockDim.x + blockIdx.y*gridDim.x*blockDim.x;
  //int bx = blockIdx.y*gridDim.x+blockIdx.x;
  //int tx = threadIdx.x;
  //int iinterface=0;

  for( int iinterface=0; iinterface < num_interfaces_ext_mesh; iinterface++) {
    if(id < d_nibool_interfaces_ext_mesh[iinterface]) {

      // for testing atomic operations against not atomic operations (0.1ms vs. 0.04 ms)
      // d_accel[3*(d_ibool_interfaces_ext_mesh[id+max_nibool_interfaces_ext_mesh*iinterface]-1)] +=
      // d_send_accel_buffer[3*(id + max_nibool_interfaces_ext_mesh*iinterface)];
      // d_accel[3*(d_ibool_interfaces_ext_mesh[id+max_nibool_interfaces_ext_mesh*iinterface]-1)+1] +=
      // d_send_accel_buffer[3*(id + max_nibool_interfaces_ext_mesh*iinterface)+1];
      // d_accel[3*(d_ibool_interfaces_ext_mesh[id+max_nibool_interfaces_ext_mesh*iinterface]-1)+2] +=
      // d_send_accel_buffer[3*(id + max_nibool_interfaces_ext_mesh*iinterface)+2];


      atomicAdd(&d_accel[3*(d_ibool_interfaces_ext_mesh[id+max_nibool_interfaces_ext_mesh*iinterface]-1)],
                d_send_accel_buffer[3*(id + max_nibool_interfaces_ext_mesh*iinterface)]);
      atomicAdd(&d_accel[3*(d_ibool_interfaces_ext_mesh[id+max_nibool_interfaces_ext_mesh*iinterface]-1)+1],
                d_send_accel_buffer[3*(id + max_nibool_interfaces_ext_mesh*iinterface)+1]);
      atomicAdd(&d_accel[3*(d_ibool_interfaces_ext_mesh[id+max_nibool_interfaces_ext_mesh*iinterface]-1)+2],
                d_send_accel_buffer[3*(id + max_nibool_interfaces_ext_mesh*iinterface)+2]);
    }
  }
  // ! This step is done via previous function transfer_and_assemble...
  // ! do iinterface = 1, num_interfaces_ext_mesh
  // !   do ipoin = 1, nibool_interfaces_ext_mesh(iinterface)
  // !     array_val(:,ibool_interfaces_ext_mesh(ipoin,iinterface)) = &
  // !          array_val(:,ibool_interfaces_ext_mesh(ipoin,iinterface)) + buffer_recv_vector_ext_mesh(:,ipoin,iinterface)
  // !   enddo
  // ! enddo
}

/* ----------------------------------------------------------------------------------------------- */

extern "C"
void FC_FUNC_(transfer_boundary_to_device_a,
              TRANSFER_BOUNDARY_TO_DEVICE_A)(long* Mesh_pointer,
                                             realw* buffer_recv_vector_ext_mesh,
                                             int* num_interfaces_ext_mesh,
                                             int* max_nibool_interfaces_ext_mesh) {

// asynchronous transfer from host to device

  TRACE("transfer_boundary_to_device_a");

  Mesh* mp = (Mesh*)(*Mesh_pointer); //get mesh pointer out of fortran integer container

  memcpy(mp->h_recv_accel_buffer,buffer_recv_vector_ext_mesh,mp->size_mpi_recv_buffer*sizeof(realw));

  // cudaMemcpyAsync(mp->d_send_accel_buffer, buffer_recv_vector_ext_mesh,
  // 3*(mp->max_nibool_interfaces_ext_mesh)*(mp->num_interfaces_ext_mesh)*sizeof(realw),
  // cudaMemcpyHostToDevice,mp->compute_stream);
  //printf("xfer to device\n");
  cudaMemcpyAsync(mp->d_send_accel_buffer, buffer_recv_vector_ext_mesh,
                  3*(mp->max_nibool_interfaces_ext_mesh)*(mp->num_interfaces_ext_mesh)*sizeof(realw),
                  cudaMemcpyHostToDevice,mp->copy_stream);




}

/* ----------------------------------------------------------------------------------------------- */

//daniel: not used ...
//
//extern "C"
//void FC_FUNC_(assemble_accel_on_device,
//              ASSEMBLE_ACCEL_on_DEVICE)(long* Mesh_pointer, realw* accel,
//                                              realw* buffer_recv_vector_ext_mesh,
//                                              int* num_interfaces_ext_mesh,
//                                              int* max_nibool_interfaces_ext_mesh,
//                                              int* nibool_interfaces_ext_mesh,
//                                              int* ibool_interfaces_ext_mesh,
//                                              int* FORWARD_OR_ADJOINT) {
//  TRACE("assemble_accel_on_device");
//
//  Mesh* mp = (Mesh*)(*Mesh_pointer); //get mesh pointer out of fortran integer container
//
//  int blocksize = BLOCKSIZE_TRANSFER;
//  int size_padded = ((int)ceil(((double)mp->max_nibool_interfaces_ext_mesh)/((double)blocksize)))*blocksize;
//  int num_blocks_x = size_padded/blocksize;
//  int num_blocks_y = 1;
//  while(num_blocks_x > 65535) {
//    num_blocks_x = (int) ceil(num_blocks_x*0.5f);
//    num_blocks_y = num_blocks_y*2;
//  }
//
//  //double start_time = get_time();
//  dim3 grid(num_blocks_x,num_blocks_y);
//  dim3 threads(blocksize,1,1);
//  // cudaEvent_t start, stop;
//  // realw time;
//  // cudaEventCreate(&start);
//  // cudaEventCreate(&stop);
//  // cudaEventRecord( start, 0 );
//
//
//  // ***************************************************************************
//  // Wait until previous copy stream finishes. We assemble while other compute kernels execute.
//  cudaStreamSynchronize(mp->copy_stream);
//
//  // Assembling on the copy_stream breaks the solution and it "blows up"
//  if(*FORWARD_OR_ADJOINT == 1) { //assemble forward accel
//    assemble_boundary_accel_on_device<<<grid,threads,0,mp->compute_stream>>>(mp->d_accel, mp->d_send_accel_buffer,
//                                                                             mp->num_interfaces_ext_mesh,
//                                                                             mp->max_nibool_interfaces_ext_mesh,
//                                                                             mp->d_nibool_interfaces_ext_mesh,
//                                                                             mp->d_ibool_interfaces_ext_mesh);
//  }
//  else if(*FORWARD_OR_ADJOINT == 3) { //assemble adjoint accel
//    assemble_boundary_accel_on_device<<<grid,threads,0,mp->copy_stream>>>(mp->d_b_accel, mp->d_send_accel_buffer,
//                                                        mp->num_interfaces_ext_mesh,
//                                                        mp->max_nibool_interfaces_ext_mesh,
//                                                        mp->d_nibool_interfaces_ext_mesh,
//                                                        mp->d_ibool_interfaces_ext_mesh);
//  }
//
//  // cudaEventRecord( stop, 0 );
//  // cudaEventSynchronize( stop );
//  // cudaEventElapsedTime( &time, start, stop );
//  // cudaEventDestroy( start );
//  // cudaEventDestroy( stop );
//  // printf("Boundary Assemble Kernel Execution Time: %f ms\n",time);
//#ifdef ENABLE_VERY_SLOW_ERROR_CHECKING
//  //double end_time = get_time();
//  //printf("Elapsed time: %e\n",end_time-start_time);
//  exit_on_cuda_error("transfer_asmbl_accel_to_device");
//#endif
//}

/* ----------------------------------------------------------------------------------------------- */

// FORWARD_OR_ADJOINT == 1 for accel, and == 3 for b_accel
extern "C"
void FC_FUNC_(transfer_asmbl_accel_to_device,
              TRANSFER_ASMBL_ACCEL_TO_DEVICE)(long* Mesh_pointer, realw* accel,
                                                    realw* buffer_recv_vector_ext_mesh,
                                                    int* num_interfaces_ext_mesh,
                                                    int* max_nibool_interfaces_ext_mesh,
                                                    int* nibool_interfaces_ext_mesh,
                                                    int* ibool_interfaces_ext_mesh,
                                                    int* FORWARD_OR_ADJOINT) {
TRACE("transfer_asmbl_accel_to_device");

  Mesh* mp = (Mesh*)(*Mesh_pointer); //get mesh pointer out of fortran integer container

  //daniel: todo - check if this copy is only needed for adjoint simulation, otherwise it is called asynchronously?
  if(*FORWARD_OR_ADJOINT == 1 ){
    // Wait until previous copy stream finishes. We assemble while other compute kernels execute.
    cudaStreamSynchronize(mp->copy_stream);
  }
  else if(*FORWARD_OR_ADJOINT == 3 ){
    cudaMemcpy(mp->d_send_accel_buffer, buffer_recv_vector_ext_mesh,
             3*(mp->max_nibool_interfaces_ext_mesh)*(mp->num_interfaces_ext_mesh)*sizeof(realw),
             cudaMemcpyHostToDevice);
  }

  int blocksize = BLOCKSIZE_TRANSFER;
  int size_padded = ((int)ceil(((double)mp->max_nibool_interfaces_ext_mesh)/((double)blocksize)))*blocksize;
  int num_blocks_x = size_padded/blocksize;
  int num_blocks_y = 1;
  while(num_blocks_x > 65535) {
    num_blocks_x = (int) ceil(num_blocks_x*0.5f);
    num_blocks_y = num_blocks_y*2;
  }

  //double start_time = get_time();
  dim3 grid(num_blocks_x,num_blocks_y);
  dim3 threads(blocksize,1,1);
  // cudaEvent_t start, stop;
  // realw time;
  // cudaEventCreate(&start);
  // cudaEventCreate(&stop);
  // cudaEventRecord( start, 0 );
  if(*FORWARD_OR_ADJOINT == 1) { //assemble forward accel
    assemble_boundary_accel_on_device<<<grid,threads,0,mp->compute_stream>>>(mp->d_accel, mp->d_send_accel_buffer,
                                                        mp->num_interfaces_ext_mesh,
                                                        mp->max_nibool_interfaces_ext_mesh,
                                                        mp->d_nibool_interfaces_ext_mesh,
                                                        mp->d_ibool_interfaces_ext_mesh);
  }
  else if(*FORWARD_OR_ADJOINT == 3) { //assemble adjoint accel
    assemble_boundary_accel_on_device<<<grid,threads,0,mp->compute_stream>>>(mp->d_b_accel, mp->d_send_accel_buffer,
                                                        mp->num_interfaces_ext_mesh,
                                                        mp->max_nibool_interfaces_ext_mesh,
                                                        mp->d_nibool_interfaces_ext_mesh,
                                                        mp->d_ibool_interfaces_ext_mesh);
  }

  // cudaEventRecord( stop, 0 );
  // cudaEventSynchronize( stop );
  // cudaEventElapsedTime( &time, start, stop );
  // cudaEventDestroy( start );
  // cudaEventDestroy( stop );
  // printf("Boundary Assemble Kernel Execution Time: %f ms\n",time);
#ifdef ENABLE_VERY_SLOW_ERROR_CHECKING
  //double end_time = get_time();
  //printf("Elapsed time: %e\n",end_time-start_time);
  exit_on_cuda_error("transfer_asmbl_accel_to_device");
#endif
}


/* ----------------------------------------------------------------------------------------------- */

// KERNEL 2

/* ----------------------------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------------------------------- */

//__global__ void Kernel_test(realw* d_debug_output,int* d_phase_ispec_inner_elastic,
//                            int num_phase_ispec_elastic, int d_iphase, int* d_ibool) {
//  int bx = blockIdx.x;
//  int tx = threadIdx.x;
//  int working_element;
//  //int ispec;
//  //int NGLL3_ALIGN = 128;
//  if(tx==0 && bx==0) {
//
//    d_debug_output[0] = 420.0;
//
//    d_debug_output[2] = num_phase_ispec_elastic;
//    d_debug_output[3] = d_iphase;
//    working_element = d_phase_ispec_inner_elastic[bx + num_phase_ispec_elastic*(d_iphase-1)]-1;
//    d_debug_output[4] = working_element;
//    d_debug_output[5] = d_phase_ispec_inner_elastic[0];
//    /* d_debug_output[1] = d_ibool[working_element*NGLL3_ALIGN + tx]-1; */
//  }
//  /* d_debug_output[1+tx+128*bx] = 69.0; */
//
//}

/* ----------------------------------------------------------------------------------------------- */

// updates stress

__device__ void compute_element_att_stress(int tx,int working_element,int NSPEC,
                                           realw* R_xx,
                                           realw* R_yy,
                                           realw* R_xy,
                                           realw* R_xz,
                                           realw* R_yz,
                                           reald* sigma_xx,
                                           reald* sigma_yy,
                                           reald* sigma_zz,
                                           reald* sigma_xy,
                                           reald* sigma_xz,
                                           reald* sigma_yz) {

  int i_sls,offset_sls;
  reald R_xx_val,R_yy_val;

  for(i_sls = 0; i_sls < N_SLS; i_sls++){
    // index
    offset_sls = tx + NGLL3*(working_element + NSPEC*i_sls);

    R_xx_val = R_xx[offset_sls]; //(i,j,k,ispec,i_sls)
    R_yy_val = R_yy[offset_sls];

    *sigma_xx = *sigma_xx - R_xx_val;
    *sigma_yy = *sigma_yy - R_yy_val;
    *sigma_zz = *sigma_zz + R_xx_val + R_yy_val;
    *sigma_xy = *sigma_xy - R_xy[offset_sls];
    *sigma_xz = *sigma_xz - R_xz[offset_sls];
    *sigma_yz = *sigma_yz - R_yz[offset_sls];
  }
  return;
}

/* ----------------------------------------------------------------------------------------------- */

// updates R_memory

__device__ void compute_element_att_memory(int tx,int working_element,int NSPEC,
                                          realw* d_muv,
                                          realw* factor_common,
                                          realw* alphaval,realw* betaval,realw* gammaval,
                                          realw* R_xx,realw* R_yy,realw* R_xy,realw* R_xz,realw* R_yz,
                                          realw* epsilondev_xx,realw* epsilondev_yy,realw* epsilondev_xy,
                                          realw* epsilondev_xz,realw* epsilondev_yz,
                                          reald epsilondev_xx_loc,reald epsilondev_yy_loc,reald epsilondev_xy_loc,
                                          reald epsilondev_xz_loc,reald epsilondev_yz_loc
                                          ){

  int i_sls;
  int ijk_ispec;
  int offset_sls,offset_align,offset_common;
  reald mul;
  reald alphaval_loc,betaval_loc,gammaval_loc;
  reald factor_loc,Sn,Snp1;

  // indices
  offset_align = tx + NGLL3_PADDED * working_element;
  ijk_ispec = tx + NGLL3 * working_element;

  mul = d_muv[offset_align];

  // use Runge-Kutta scheme to march in time
  for(i_sls = 0; i_sls < N_SLS; i_sls++){

    // indices
    offset_common = i_sls + N_SLS*(tx + NGLL3*working_element); // (i_sls,i,j,k,ispec)
    offset_sls = tx + NGLL3*(working_element + NSPEC*i_sls);   // (i,j,k,ispec,i_sls)

    factor_loc = mul * factor_common[offset_common]; //mustore(i,j,k,ispec) * factor_common(i_sls,i,j,k,ispec)

    alphaval_loc = alphaval[i_sls]; // (i_sls)
    betaval_loc = betaval[i_sls];
    gammaval_loc = gammaval[i_sls];

    // term in xx
    Sn   = factor_loc * epsilondev_xx[ijk_ispec]; //(i,j,k,ispec)
    Snp1   = factor_loc * epsilondev_xx_loc; //(i,j,k)

    //R_xx(i,j,k,ispec,i_sls) = alphaval_loc * R_xx(i,j,k,ispec,i_sls) +
    //                          betaval_loc * Sn + gammaval_loc * Snp1;

    R_xx[offset_sls] = alphaval_loc * R_xx[offset_sls] +
                       betaval_loc * Sn + gammaval_loc * Snp1;

    // term in yy
    Sn   = factor_loc * epsilondev_yy[ijk_ispec];
    Snp1   = factor_loc * epsilondev_yy_loc;
    R_yy[offset_sls] = alphaval_loc * R_yy[offset_sls] +
                        betaval_loc * Sn + gammaval_loc * Snp1;
    // term in zz not computed since zero trace
    // term in xy
    Sn   = factor_loc * epsilondev_xy[ijk_ispec];
    Snp1   = factor_loc * epsilondev_xy_loc;
    R_xy[offset_sls] = alphaval_loc * R_xy[offset_sls] +
                        betaval_loc * Sn + gammaval_loc * Snp1;
    // term in xz
    Sn   = factor_loc * epsilondev_xz[ijk_ispec];
    Snp1   = factor_loc * epsilondev_xz_loc;
    R_xz[offset_sls] = alphaval_loc * R_xz[offset_sls] +
                        betaval_loc * Sn + gammaval_loc * Snp1;
    // term in yz
    Sn   = factor_loc * epsilondev_yz[ijk_ispec];
    Snp1   = factor_loc * epsilondev_yz_loc;
    R_yz[offset_sls] = alphaval_loc * R_yz[offset_sls] +
                        betaval_loc * Sn + gammaval_loc * Snp1;
  }
  return;
}

/* ----------------------------------------------------------------------------------------------- */

// pre-computes gravity term

__device__ void compute_element_gravity(int tx,int working_element,
                                        int* d_ibool,
                                        realw* d_minus_g,
                                        realw* d_minus_deriv_gravity,
                                        realw* d_rhostore,
                                        realw* wgll_cube,
                                        reald jacobianl,
                                        reald* s_dummyx_loc,
                                        reald* s_dummyy_loc,
                                        reald* s_dummyz_loc,
                                        reald* sigma_xx,
                                        reald* sigma_yy,
                                        reald* sigma_xz,
                                        reald* sigma_yz,
                                        reald* rho_s_H1,
                                        reald* rho_s_H2,
                                        reald* rho_s_H3){

  int iglob;
  reald minus_g,minus_dg;
  reald rhol;
  reald gzl; // gxl,gyl,
  reald sx_l,sy_l,sz_l;
  reald Hxxl,Hyyl,Hzzl; //,Hxyl,Hxzl,Hyzl;
  reald factor;

  // compute non-symmetric terms for gravity

  // get g, rho and dg/dr=dg
  iglob = d_ibool[working_element*NGLL3 + tx]-1;

  minus_g = d_minus_g[iglob];
  minus_dg = d_minus_deriv_gravity[iglob];

  // Cartesian components of the gravitational acceleration
  //gxl = 0.f;
  //gyl = 0.f;
  gzl = minus_g;

  // Cartesian components of gradient of gravitational acceleration
  // H = grad g
  // assumes g only acts in negative z-direction
  Hxxl = 0.f;
  Hyyl = 0.f;
  Hzzl = minus_dg;
  //Hxyl = 0.f;
  //Hxzl = 0.f;
  //Hyzl = 0.f;

  rhol = d_rhostore[working_element*NGLL3_PADDED + tx];

  // get displacement and multiply by density to compute G tensor
  // G = rho [ sg - (s * g) I  ]
  sx_l = rhol * s_dummyx_loc[tx]; // d_displ[iglob*3];
  sy_l = rhol * s_dummyy_loc[tx]; // d_displ[iglob*3 + 1];
  sz_l = rhol * s_dummyz_loc[tx]; // d_displ[iglob*3 + 2];

  // compute G tensor from s . g and add to sigma (not symmetric)
  //sigma_xx += sy_l*gyl + sz_l*gzl;
  *sigma_xx += sz_l*gzl;
  //sigma_yy += sx_l*gxl + sz_l*gzl;
  *sigma_yy += sz_l*gzl;
  //sigma_zz += sx_l*gxl + sy_l*gyl;

  //sigma_xy -= sx_l*gyl;
  //sigma_yx -= sy_l*gxl;

  *sigma_xz -= sx_l*gzl;
  //sigma_zx -= sz_l*gxl;

  *sigma_yz -= sy_l*gzl;
  //sigma_zy -= sz_l*gyl;

  // precompute vector
  factor = jacobianl * wgll_cube[tx];

  //rho_s_H1 = fac1 * (sx_l * Hxxl + sy_l * Hxyl + sz_l * Hxzl);
  //rho_s_H2 = fac1 * (sx_l * Hxyl + sy_l * Hyyl + sz_l * Hyzl);
  //rho_s_H3 = fac1 * (sx_l * Hxzl + sy_l * Hyzl + sz_l * Hzzl);

  // only non-zero z-direction
  *rho_s_H1 = factor * sx_l * Hxxl ; // 0.f;
  *rho_s_H2 = factor * sy_l * Hyyl ; // 0.f;
  *rho_s_H3 = factor * sz_l * Hzzl ;

  // debug
  //*rho_s_H1 = 0.f;
  //*rho_s_H2 = 0.f;
  //*rho_s_H3 = 0.f ;

}

/* ----------------------------------------------------------------------------------------------- */

// double precision temporary variables leads to 10% performance
// decrease in Kernel_2_impl (not very much..)
//typedef realw reald;
#ifdef USE_TEXTURES_FIELDS
texture<realw, cudaTextureType1D, cudaReadModeElementType> d_displ_tex;
texture<realw, cudaTextureType1D, cudaReadModeElementType> d_accel_tex;
#endif

#ifdef USE_TEXTURES_CONSTANTS
texture<realw, cudaTextureType1D, cudaReadModeElementType> d_hprime_xx_tex;
#endif

__global__ void Kernel_2_impl(int nb_blocks_to_compute,
                              int NGLOB,
                              int* d_ibool,
                              int* d_phase_ispec_inner_elastic, int num_phase_ispec_elastic,
                              int d_iphase,
                              int use_mesh_coloring_gpu,
                              realw* d_displ, realw* d_accel,
                              realw* d_xix, realw* d_xiy, realw* d_xiz,
                              realw* d_etax, realw* d_etay, realw* d_etaz,
                              realw* d_gammax, realw* d_gammay, realw* d_gammaz,
                              realw* d_kappav, realw* d_muv,
                              int COMPUTE_AND_STORE_STRAIN,
                              realw* epsilondev_xx,realw* epsilondev_yy,realw* epsilondev_xy,
                              realw* epsilondev_xz,realw* epsilondev_yz,
                              realw* epsilon_trace_over_3,
                              int SIMULATION_TYPE,
                              int ATTENUATION,
                              int NSPEC,
                              realw* one_minus_sum_beta,realw* factor_common,
                              realw* R_xx, realw* R_yy, realw* R_xy, realw* R_xz, realw* R_yz,
                              realw* alphaval,realw* betaval,realw* gammaval,
                              int ANISOTROPY,
                              realw* d_c11store,
                              realw* d_c12store,
                              realw* d_c13store,
                              realw* d_c14store,
                              realw* d_c15store,
                              realw* d_c16store,
                              realw* d_c22store,
                              realw* d_c23store,
                              realw* d_c24store,
                              realw* d_c25store,
                              realw* d_c26store,
                              realw* d_c33store,
                              realw* d_c34store,
                              realw* d_c35store,
                              realw* d_c36store,
                              realw* d_c44store,
                              realw* d_c45store,
                              realw* d_c46store,
                              realw* d_c55store,
                              realw* d_c56store,
                              realw* d_c66store,
                              int gravity,
                              realw* d_minus_g,
                              realw* d_minus_deriv_gravity,
                              realw* d_rhostore,
                              realw* wgll_cube){

  /* int bx = blockIdx.y*blockDim.x+blockIdx.x; //possible bug in original code*/
  int bx = blockIdx.y*gridDim.x+blockIdx.x;
  /* int bx = blockIdx.x; */
  int tx = threadIdx.x;

  //const int NGLLX = 5;
  // const int NGLL2 = 25;
  //const int NGLL3 = NGLL3;
  const int NGLL3_ALIGN = NGLL3_PADDED;

  int K = (tx/NGLL2);
  int J = ((tx-K*NGLL2)/NGLLX);
  int I = (tx-K*NGLL2-J*NGLLX);

  int active,offset;
  int iglob = 0;
  int working_element;

  reald tempx1l,tempx2l,tempx3l,tempy1l,tempy2l,tempy3l,tempz1l,tempz2l,tempz3l;
  reald xixl,xiyl,xizl,etaxl,etayl,etazl,gammaxl,gammayl,gammazl,jacobianl;
  reald duxdxl,duxdyl,duxdzl,duydxl,duydyl,duydzl,duzdxl,duzdyl,duzdzl;
  reald duxdxl_plus_duydyl,duxdxl_plus_duzdzl,duydyl_plus_duzdzl;
  reald duxdyl_plus_duydxl,duzdxl_plus_duxdzl,duzdyl_plus_duydzl;
  reald fac1,fac2,fac3,lambdal,mul,lambdalplus2mul,kappal;
  reald sigma_xx,sigma_yy,sigma_zz,sigma_xy,sigma_xz,sigma_yz;
  reald epsilondev_xx_loc,epsilondev_yy_loc,epsilondev_xy_loc,epsilondev_xz_loc,epsilondev_yz_loc;
  reald c11,c12,c13,c14,c15,c16,c22,c23,c24,c25,c26,c33,c34,c35,c36,c44,c45,c46,c55,c56,c66;
  reald sum_terms1,sum_terms2,sum_terms3;

  // gravity variables
  reald sigma_yx,sigma_zx,sigma_zy;
  reald rho_s_H1,rho_s_H2,rho_s_H3;

#ifndef MANUALLY_UNROLLED_LOOPS
    int l;
    realw hp1,hp2,hp3;
#endif

    __shared__ reald s_dummyx_loc[NGLL3];
    __shared__ reald s_dummyy_loc[NGLL3];
    __shared__ reald s_dummyz_loc[NGLL3];

    __shared__ reald s_tempx1[NGLL3];
    __shared__ reald s_tempx2[NGLL3];
    __shared__ reald s_tempx3[NGLL3];
    __shared__ reald s_tempy1[NGLL3];
    __shared__ reald s_tempy2[NGLL3];
    __shared__ reald s_tempy3[NGLL3];
    __shared__ reald s_tempz1[NGLL3];
    __shared__ reald s_tempz2[NGLL3];
    __shared__ reald s_tempz3[NGLL3];

    __shared__ reald sh_hprime_xx[NGLL2];

// use only NGLL^3 = 125 active threads, plus 3 inactive/ghost threads,
// because we used memory padding from NGLL^3 = 125 to 128 to get coalescent memory accesses
    active = (tx < NGLL3 && bx < nb_blocks_to_compute) ? 1:0;

// copy from global memory to shared memory
// each thread writes one of the NGLL^3 = 125 data points
    if (active) {

#ifdef USE_MESH_COLORING_GPU
      working_element = bx;
#else
      //mesh coloring
      if( use_mesh_coloring_gpu ){
        working_element = bx;
      }else{
        // iphase-1 and working_element-1 for Fortran->C array conventions
        working_element = d_phase_ispec_inner_elastic[bx + num_phase_ispec_elastic*(d_iphase-1)]-1;
      }
#endif

      // iglob = d_ibool[working_element*NGLL3_ALIGN + tx]-1;
      iglob = d_ibool[working_element*NGLL3 + tx]-1;

#ifdef USE_TEXTURES_FIELDS
      s_dummyx_loc[tx] = tex1Dfetch(d_displ_tex, iglob*3);
      s_dummyy_loc[tx] = tex1Dfetch(d_displ_tex, iglob*3 + 1);
      s_dummyz_loc[tx] = tex1Dfetch(d_displ_tex, iglob*3 + 2);
#else
      // changing iglob indexing to match fortran row changes fast style
      s_dummyx_loc[tx] = d_displ[iglob*3];
      s_dummyy_loc[tx] = d_displ[iglob*3 + 1];
      s_dummyz_loc[tx] = d_displ[iglob*3 + 2];
#endif
    }

    if (tx < NGLL2) {
      #ifdef USE_TEXTURES_CONSTANTS
      sh_hprime_xx[tx] = tex1Dfetch(d_hprime_xx_tex,tx);
      #else
      sh_hprime_xx[tx] = d_hprime_xx[tx];
      #endif
    }

// synchronize all the threads (one thread for each of the NGLL grid points of the
// current spectral element) because we need the whole element to be ready in order
// to be able to compute the matrix products along cut planes of the 3D element below
    __syncthreads();

    if (active) {

#ifndef MANUALLY_UNROLLED_LOOPS

      tempx1l = 0.f;
      tempx2l = 0.f;
      tempx3l = 0.f;

      tempy1l = 0.f;
      tempy2l = 0.f;
      tempy3l = 0.f;

      tempz1l = 0.f;
      tempz2l = 0.f;
      tempz3l = 0.f;

      for (l=0;l<NGLLX;l++) {
          hp1 = sh_hprime_xx[l*NGLLX+I];
          offset = K*NGLL2+J*NGLLX+l;
          tempx1l += s_dummyx_loc[offset]*hp1;
          tempy1l += s_dummyy_loc[offset]*hp1;
          tempz1l += s_dummyz_loc[offset]*hp1;

          hp2 = sh_hprime_xx[l*NGLLX+J];
          offset = K*NGLL2+l*NGLLX+I;
          tempx2l += s_dummyx_loc[offset]*hp2;
          tempy2l += s_dummyy_loc[offset]*hp2;
          tempz2l += s_dummyz_loc[offset]*hp2;

          hp3 = sh_hprime_xx[l*NGLLX+K];
          offset = l*NGLL2+J*NGLLX+I;
          tempx3l += s_dummyx_loc[offset]*hp3;
          tempy3l += s_dummyy_loc[offset]*hp3;
          tempz3l += s_dummyz_loc[offset]*hp3;

      }
#else

      tempx1l = s_dummyx_loc[K*NGLL2+J*NGLLX]*d_hprime_xx[I]
              + s_dummyx_loc[K*NGLL2+J*NGLLX+1]*d_hprime_xx[NGLLX+I]
              + s_dummyx_loc[K*NGLL2+J*NGLLX+2]*d_hprime_xx[2*NGLLX+I]
              + s_dummyx_loc[K*NGLL2+J*NGLLX+3]*d_hprime_xx[3*NGLLX+I]
              + s_dummyx_loc[K*NGLL2+J*NGLLX+4]*d_hprime_xx[4*NGLLX+I];

      tempy1l = s_dummyy_loc[K*NGLL2+J*NGLLX]*d_hprime_xx[I]
              + s_dummyy_loc[K*NGLL2+J*NGLLX+1]*d_hprime_xx[NGLLX+I]
              + s_dummyy_loc[K*NGLL2+J*NGLLX+2]*d_hprime_xx[2*NGLLX+I]
              + s_dummyy_loc[K*NGLL2+J*NGLLX+3]*d_hprime_xx[3*NGLLX+I]
              + s_dummyy_loc[K*NGLL2+J*NGLLX+4]*d_hprime_xx[4*NGLLX+I];

      tempz1l = s_dummyz_loc[K*NGLL2+J*NGLLX]*d_hprime_xx[I]
              + s_dummyz_loc[K*NGLL2+J*NGLLX+1]*d_hprime_xx[NGLLX+I]
              + s_dummyz_loc[K*NGLL2+J*NGLLX+2]*d_hprime_xx[2*NGLLX+I]
              + s_dummyz_loc[K*NGLL2+J*NGLLX+3]*d_hprime_xx[3*NGLLX+I]
              + s_dummyz_loc[K*NGLL2+J*NGLLX+4]*d_hprime_xx[4*NGLLX+I];

      tempx2l = s_dummyx_loc[K*NGLL2+I]*d_hprime_xx[J]
              + s_dummyx_loc[K*NGLL2+NGLLX+I]*d_hprime_xx[NGLLX+J]
              + s_dummyx_loc[K*NGLL2+2*NGLLX+I]*d_hprime_xx[2*NGLLX+J]
              + s_dummyx_loc[K*NGLL2+3*NGLLX+I]*d_hprime_xx[3*NGLLX+J]
              + s_dummyx_loc[K*NGLL2+4*NGLLX+I]*d_hprime_xx[4*NGLLX+J];

      tempy2l = s_dummyy_loc[K*NGLL2+I]*d_hprime_xx[J]
              + s_dummyy_loc[K*NGLL2+NGLLX+I]*d_hprime_xx[NGLLX+J]
              + s_dummyy_loc[K*NGLL2+2*NGLLX+I]*d_hprime_xx[2*NGLLX+J]
              + s_dummyy_loc[K*NGLL2+3*NGLLX+I]*d_hprime_xx[3*NGLLX+J]
              + s_dummyy_loc[K*NGLL2+4*NGLLX+I]*d_hprime_xx[4*NGLLX+J];

      tempz2l = s_dummyz_loc[K*NGLL2+I]*d_hprime_xx[J]
              + s_dummyz_loc[K*NGLL2+NGLLX+I]*d_hprime_xx[NGLLX+J]
              + s_dummyz_loc[K*NGLL2+2*NGLLX+I]*d_hprime_xx[2*NGLLX+J]
              + s_dummyz_loc[K*NGLL2+3*NGLLX+I]*d_hprime_xx[3*NGLLX+J]
              + s_dummyz_loc[K*NGLL2+4*NGLLX+I]*d_hprime_xx[4*NGLLX+J];

      tempx3l = s_dummyx_loc[J*NGLLX+I]*d_hprime_xx[K]
              + s_dummyx_loc[NGLL2+J*NGLLX+I]*d_hprime_xx[NGLLX+K]
              + s_dummyx_loc[2*NGLL2+J*NGLLX+I]*d_hprime_xx[2*NGLLX+K]
              + s_dummyx_loc[3*NGLL2+J*NGLLX+I]*d_hprime_xx[3*NGLLX+K]
              + s_dummyx_loc[4*NGLL2+J*NGLLX+I]*d_hprime_xx[4*NGLLX+K];

      tempy3l = s_dummyy_loc[J*NGLLX+I]*d_hprime_xx[K]
              + s_dummyy_loc[NGLL2+J*NGLLX+I]*d_hprime_xx[NGLLX+K]
              + s_dummyy_loc[2*NGLL2+J*NGLLX+I]*d_hprime_xx[2*NGLLX+K]
              + s_dummyy_loc[3*NGLL2+J*NGLLX+I]*d_hprime_xx[3*NGLLX+K]
              + s_dummyy_loc[4*NGLL2+J*NGLLX+I]*d_hprime_xx[4*NGLLX+K];

      tempz3l = s_dummyz_loc[J*NGLLX+I]*d_hprime_xx[K]
              + s_dummyz_loc[NGLL2+J*NGLLX+I]*d_hprime_xx[NGLLX+K]
              + s_dummyz_loc[2*NGLL2+J*NGLLX+I]*d_hprime_xx[2*NGLLX+K]
              + s_dummyz_loc[3*NGLL2+J*NGLLX+I]*d_hprime_xx[3*NGLLX+K]
              + s_dummyz_loc[4*NGLL2+J*NGLLX+I]*d_hprime_xx[4*NGLLX+K];

#endif

// compute derivatives of ux, uy and uz with respect to x, y and z
      offset = working_element*NGLL3_ALIGN + tx;

      xixl = d_xix[offset];
      xiyl = d_xiy[offset];
      xizl = d_xiz[offset];
      etaxl = d_etax[offset];
      etayl = d_etay[offset];
      etazl = d_etaz[offset];
      gammaxl = d_gammax[offset];
      gammayl = d_gammay[offset];
      gammazl = d_gammaz[offset];

      duxdxl = xixl*tempx1l + etaxl*tempx2l + gammaxl*tempx3l;
      duxdyl = xiyl*tempx1l + etayl*tempx2l + gammayl*tempx3l;
      duxdzl = xizl*tempx1l + etazl*tempx2l + gammazl*tempx3l;

      duydxl = xixl*tempy1l + etaxl*tempy2l + gammaxl*tempy3l;
      duydyl = xiyl*tempy1l + etayl*tempy2l + gammayl*tempy3l;
      duydzl = xizl*tempy1l + etazl*tempy2l + gammazl*tempy3l;

      duzdxl = xixl*tempz1l + etaxl*tempz2l + gammaxl*tempz3l;
      duzdyl = xiyl*tempz1l + etayl*tempz2l + gammayl*tempz3l;
      duzdzl = xizl*tempz1l + etazl*tempz2l + gammazl*tempz3l;

      // precompute some sums to save CPU time
      duxdxl_plus_duydyl = duxdxl + duydyl;
      duxdxl_plus_duzdzl = duxdxl + duzdzl;
      duydyl_plus_duzdzl = duydyl + duzdzl;
      duxdyl_plus_duydxl = duxdyl + duydxl;
      duzdxl_plus_duxdzl = duzdxl + duxdzl;
      duzdyl_plus_duydzl = duzdyl + duydzl;

      // computes deviatoric strain attenuation and/or for kernel calculations
      if(COMPUTE_AND_STORE_STRAIN) {
        realw templ = 0.33333333333333333333f * (duxdxl + duydyl + duzdzl); // 1./3. = 0.33333
        /*
        epsilondev_xx[offset] = duxdxl - templ;
        epsilondev_yy[offset] = duydyl - templ;
        epsilondev_xy[offset] = 0.5f * duxdyl_plus_duydxl;
        epsilondev_xz[offset] = 0.5f * duzdxl_plus_duxdzl;
        epsilondev_yz[offset] = 0.5f * duzdyl_plus_duydzl;
              */
        // local storage: stresses at this current time step
        epsilondev_xx_loc = duxdxl - templ;
        epsilondev_yy_loc = duydyl - templ;
        epsilondev_xy_loc = 0.5f * duxdyl_plus_duydxl;
        epsilondev_xz_loc = 0.5f * duzdxl_plus_duxdzl;
        epsilondev_yz_loc = 0.5f * duzdyl_plus_duydzl;

        if(SIMULATION_TYPE == 3) {
          epsilon_trace_over_3[tx + working_element*NGLL3] = templ;
        }
      }

      // compute elements with an elastic isotropic rheology
      kappal = d_kappav[offset];
      mul = d_muv[offset];

      // attenuation
      if(ATTENUATION){
        // use unrelaxed parameters if attenuation
        mul  = mul * one_minus_sum_beta[tx+working_element*NGLL3]; // (i,j,k,ispec)
      }

      // full anisotropic case, stress calculations
      if(ANISOTROPY){

        c11 = d_c11store[offset];
        c12 = d_c12store[offset];
        c13 = d_c13store[offset];
        c14 = d_c14store[offset];
        c15 = d_c15store[offset];
        c16 = d_c16store[offset];
        c22 = d_c22store[offset];
        c23 = d_c23store[offset];
        c24 = d_c24store[offset];
        c25 = d_c25store[offset];
        c26 = d_c26store[offset];
        c33 = d_c33store[offset];
        c34 = d_c34store[offset];
        c35 = d_c35store[offset];
        c36 = d_c36store[offset];
        c44 = d_c44store[offset];
        c45 = d_c45store[offset];
        c46 = d_c46store[offset];
        c55 = d_c55store[offset];
        c56 = d_c56store[offset];
        c66 = d_c66store[offset];

        sigma_xx = c11*duxdxl + c16*duxdyl_plus_duydxl + c12*duydyl +
                   c15*duzdxl_plus_duxdzl + c14*duzdyl_plus_duydzl + c13*duzdzl;
        sigma_yy = c12*duxdxl + c26*duxdyl_plus_duydxl + c22*duydyl +
                   c25*duzdxl_plus_duxdzl + c24*duzdyl_plus_duydzl + c23*duzdzl;
        sigma_zz = c13*duxdxl + c36*duxdyl_plus_duydxl + c23*duydyl +
                   c35*duzdxl_plus_duxdzl + c34*duzdyl_plus_duydzl + c33*duzdzl;
        sigma_xy = c16*duxdxl + c66*duxdyl_plus_duydxl + c26*duydyl +
                   c56*duzdxl_plus_duxdzl + c46*duzdyl_plus_duydzl + c36*duzdzl;
        sigma_xz = c15*duxdxl + c56*duxdyl_plus_duydxl + c25*duydyl +
                   c55*duzdxl_plus_duxdzl + c45*duzdyl_plus_duydzl + c35*duzdzl;
        sigma_yz = c14*duxdxl + c46*duxdyl_plus_duydxl + c24*duydyl +
                   c45*duzdxl_plus_duxdzl + c44*duzdyl_plus_duydzl + c34*duzdzl;

      }else{

        // isotropic case

        lambdalplus2mul = kappal + 1.33333333333333333333f * mul;  // 4./3. = 1.3333333
        lambdal = lambdalplus2mul - 2.0f * mul;

        // compute the six components of the stress tensor sigma
        sigma_xx = lambdalplus2mul*duxdxl + lambdal*duydyl_plus_duzdzl;
        sigma_yy = lambdalplus2mul*duydyl + lambdal*duxdxl_plus_duzdzl;
        sigma_zz = lambdalplus2mul*duzdzl + lambdal*duxdxl_plus_duydyl;

        sigma_xy = mul*duxdyl_plus_duydxl;
        sigma_xz = mul*duzdxl_plus_duxdzl;
        sigma_yz = mul*duzdyl_plus_duydzl;
      }

      if(ATTENUATION){
        // subtracts memory variables if attenuation
        compute_element_att_stress(tx,working_element,NSPEC,
                                   R_xx,R_yy,R_xy,R_xz,R_yz,
                                   &sigma_xx,&sigma_yy,&sigma_zz,&sigma_xy,&sigma_xz,&sigma_yz);
      }

      jacobianl = 1.0f / (xixl*(etayl*gammazl-etazl*gammayl)-xiyl*(etaxl*gammazl-etazl*gammaxl)+xizl*(etaxl*gammayl-etayl*gammaxl));

      // define symmetric components (needed for non-symmetric dot product and sigma for gravity)
      sigma_yx = sigma_xy;
      sigma_zx = sigma_xz;
      sigma_zy = sigma_yz;

      if( gravity ){
        //  computes non-symmetric terms for gravity
        compute_element_gravity(tx,working_element,d_ibool,d_minus_g,d_minus_deriv_gravity,
                                d_rhostore,wgll_cube,jacobianl,
                                s_dummyx_loc,s_dummyy_loc,s_dummyz_loc,
                                &sigma_xx,&sigma_yy,&sigma_xz,&sigma_yz,
                                &rho_s_H1,&rho_s_H2,&rho_s_H3);
      }

      // form dot product with test vector, non-symmetric form
      s_tempx1[tx] = jacobianl * (sigma_xx*xixl + sigma_yx*xiyl + sigma_zx*xizl);
      s_tempy1[tx] = jacobianl * (sigma_xy*xixl + sigma_yy*xiyl + sigma_zy*xizl);
      s_tempz1[tx] = jacobianl * (sigma_xz*xixl + sigma_yz*xiyl + sigma_zz*xizl);

      s_tempx2[tx] = jacobianl * (sigma_xx*etaxl + sigma_yx*etayl + sigma_zx*etazl);
      s_tempy2[tx] = jacobianl * (sigma_xy*etaxl + sigma_yy*etayl + sigma_zy*etazl);
      s_tempz2[tx] = jacobianl * (sigma_xz*etaxl + sigma_yz*etayl + sigma_zz*etazl);

      s_tempx3[tx] = jacobianl * (sigma_xx*gammaxl + sigma_yx*gammayl + sigma_zx*gammazl);
      s_tempy3[tx] = jacobianl * (sigma_xy*gammaxl + sigma_yy*gammayl + sigma_zy*gammazl);
      s_tempz3[tx] = jacobianl * (sigma_xz*gammaxl + sigma_yz*gammayl + sigma_zz*gammazl);

    }

// synchronize all the threads (one thread for each of the NGLL grid points of the
// current spectral element) because we need the whole element to be ready in order
// to be able to compute the matrix products along cut planes of the 3D element below
    __syncthreads();

    if (active) {

#ifndef MANUALLY_UNROLLED_LOOPS

      tempx1l = 0.f;
      tempy1l = 0.f;
      tempz1l = 0.f;

      tempx2l = 0.f;
      tempy2l = 0.f;
      tempz2l = 0.f;

      tempx3l = 0.f;
      tempy3l = 0.f;
      tempz3l = 0.f;

      for (l=0;l<NGLLX;l++) {

        fac1 = d_hprimewgll_xx[I*NGLLX+l];
        offset = K*NGLL2+J*NGLLX+l;
        tempx1l += s_tempx1[offset]*fac1;
        tempy1l += s_tempy1[offset]*fac1;
        tempz1l += s_tempz1[offset]*fac1;

        fac2 = d_hprimewgll_yy[J*NGLLX+l];
        offset = K*NGLL2+l*NGLLX+I;
        tempx2l += s_tempx2[offset]*fac2;
        tempy2l += s_tempy2[offset]*fac2;
        tempz2l += s_tempz2[offset]*fac2;

        fac3 = d_hprimewgll_zz[K*NGLLX+l];
        offset = l*NGLL2+J*NGLLX+I;
        tempx3l += s_tempx3[offset]*fac3;
        tempy3l += s_tempy3[offset]*fac3;
        tempz3l += s_tempz3[offset]*fac3;

      }
#else

      tempx1l = s_tempx1[K*NGLL2+J*NGLLX]*d_hprimewgll_xx[I*NGLLX]
              + s_tempx1[K*NGLL2+J*NGLLX+1]*d_hprimewgll_xx[I*NGLLX+1]
              + s_tempx1[K*NGLL2+J*NGLLX+2]*d_hprimewgll_xx[I*NGLLX+2]
              + s_tempx1[K*NGLL2+J*NGLLX+3]*d_hprimewgll_xx[I*NGLLX+3]
              + s_tempx1[K*NGLL2+J*NGLLX+4]*d_hprimewgll_xx[I*NGLLX+4];

      tempy1l = s_tempy1[K*NGLL2+J*NGLLX]*d_hprimewgll_xx[I*NGLLX]
              + s_tempy1[K*NGLL2+J*NGLLX+1]*d_hprimewgll_xx[I*NGLLX+1]
              + s_tempy1[K*NGLL2+J*NGLLX+2]*d_hprimewgll_xx[I*NGLLX+2]
              + s_tempy1[K*NGLL2+J*NGLLX+3]*d_hprimewgll_xx[I*NGLLX+3]
              + s_tempy1[K*NGLL2+J*NGLLX+4]*d_hprimewgll_xx[I*NGLLX+4];

      tempz1l = s_tempz1[K*NGLL2+J*NGLLX]*d_hprimewgll_xx[I*NGLLX]
              + s_tempz1[K*NGLL2+J*NGLLX+1]*d_hprimewgll_xx[I*NGLLX+1]
              + s_tempz1[K*NGLL2+J*NGLLX+2]*d_hprimewgll_xx[I*NGLLX+2]
              + s_tempz1[K*NGLL2+J*NGLLX+3]*d_hprimewgll_xx[I*NGLLX+3]
              + s_tempz1[K*NGLL2+J*NGLLX+4]*d_hprimewgll_xx[I*NGLLX+4];

      tempx2l = s_tempx2[K*NGLL2+I]*d_hprimewgll_yy[J*NGLLX]
              + s_tempx2[K*NGLL2+NGLLX+I]*d_hprimewgll_yy[J*NGLLX+1]
              + s_tempx2[K*NGLL2+2*NGLLX+I]*d_hprimewgll_yy[J*NGLLX+2]
              + s_tempx2[K*NGLL2+3*NGLLX+I]*d_hprimewgll_yy[J*NGLLX+3]
              + s_tempx2[K*NGLL2+4*NGLLX+I]*d_hprimewgll_yy[J*NGLLX+4];

      tempy2l = s_tempy2[K*NGLL2+I]*d_hprimewgll_yy[J*NGLLX]
              + s_tempy2[K*NGLL2+NGLLX+I]*d_hprimewgll_yy[J*NGLLX+1]
              + s_tempy2[K*NGLL2+2*NGLLX+I]*d_hprimewgll_yy[J*NGLLX+2]
              + s_tempy2[K*NGLL2+3*NGLLX+I]*d_hprimewgll_yy[J*NGLLX+3]
              + s_tempy2[K*NGLL2+4*NGLLX+I]*d_hprimewgll_yy[J*NGLLX+4];

      tempz2l = s_tempz2[K*NGLL2+I]*d_hprimewgll_yy[J*NGLLX]
              + s_tempz2[K*NGLL2+NGLLX+I]*d_hprimewgll_yy[J*NGLLX+1]
              + s_tempz2[K*NGLL2+2*NGLLX+I]*d_hprimewgll_yy[J*NGLLX+2]
              + s_tempz2[K*NGLL2+3*NGLLX+I]*d_hprimewgll_yy[J*NGLLX+3]
              + s_tempz2[K*NGLL2+4*NGLLX+I]*d_hprimewgll_yy[J*NGLLX+4];

      tempx3l = s_tempx3[J*NGLLX+I]*d_hprimewgll_zz[K*NGLLX]
              + s_tempx3[NGLL2+J*NGLLX+I]*d_hprimewgll_zz[K*NGLLX+1]
              + s_tempx3[2*NGLL2+J*NGLLX+I]*d_hprimewgll_zz[K*NGLLX+2]
              + s_tempx3[3*NGLL2+J*NGLLX+I]*d_hprimewgll_zz[K*NGLLX+3]
              + s_tempx3[4*NGLL2+J*NGLLX+I]*d_hprimewgll_zz[K*NGLLX+4];

      tempy3l = s_tempy3[J*NGLLX+I]*d_hprimewgll_zz[K*NGLLX]
              + s_tempy3[NGLL2+J*NGLLX+I]*d_hprimewgll_zz[K*NGLLX+1]
              + s_tempy3[2*NGLL2+J*NGLLX+I]*d_hprimewgll_zz[K*NGLLX+2]
              + s_tempy3[3*NGLL2+J*NGLLX+I]*d_hprimewgll_zz[K*NGLLX+3]
              + s_tempy3[4*NGLL2+J*NGLLX+I]*d_hprimewgll_zz[K*NGLLX+4];

      tempz3l = s_tempz3[J*NGLLX+I]*d_hprimewgll_zz[K*NGLLX]
              + s_tempz3[NGLL2+J*NGLLX+I]*d_hprimewgll_zz[K*NGLLX+1]
              + s_tempz3[2*NGLL2+J*NGLLX+I]*d_hprimewgll_zz[K*NGLLX+2]
              + s_tempz3[3*NGLL2+J*NGLLX+I]*d_hprimewgll_zz[K*NGLLX+3]
              + s_tempz3[4*NGLL2+J*NGLLX+I]*d_hprimewgll_zz[K*NGLLX+4];

#endif

      fac1 = d_wgllwgll_yz[K*NGLLX+J];
      fac2 = d_wgllwgll_xz[K*NGLLX+I];
      fac3 = d_wgllwgll_xy[J*NGLLX+I];

      sum_terms1 = - (fac1*tempx1l + fac2*tempx2l + fac3*tempx3l);
      sum_terms2 = - (fac1*tempy1l + fac2*tempy2l + fac3*tempy3l);
      sum_terms3 = - (fac1*tempz1l + fac2*tempz2l + fac3*tempz3l);

      // adds gravity term
      if( gravity ){
        sum_terms1 += rho_s_H1;
        sum_terms2 += rho_s_H2;
        sum_terms3 += rho_s_H3;
      }


#ifdef USE_MESH_COLORING_GPU
      // no atomic operation needed, colors don't share global points between elements

#ifdef USE_TEXTURES_FIELDS
      d_accel[iglob*3]     = tex1Dfetch(d_accel_tex, iglob*3) + sum_terms1;
      d_accel[iglob*3 + 1] = tex1Dfetch(d_accel_tex, iglob*3 + 1) + sum_terms2;
      d_accel[iglob*3 + 2] = tex1Dfetch(d_accel_tex, iglob*3 + 2) + sum_terms3;
#else
      d_accel[iglob*3]     += sum_terms1;
      d_accel[iglob*3 + 1] += sum_terms2;
      d_accel[iglob*3 + 2] += sum_terms3;
#endif // USE_TEXTURES_FIELDS

#else // MESH_COLORING
      //mesh coloring
      if( use_mesh_coloring_gpu ){

        // no atomic operation needed, colors don't share global points between elements
        // d_accel[iglob*3]     += sum_terms1;
        // d_accel[iglob*3 + 1] += sum_terms2;
        // d_accel[iglob*3 + 2] += sum_terms3;
#ifdef USE_TEXTURES_FIELDS
        d_accel[iglob*3]     = tex1Dfetch(d_accel_tex, iglob*3) + sum_terms1;
        d_accel[iglob*3 + 1] = tex1Dfetch(d_accel_tex, iglob*3 + 1) + sum_terms2;
        d_accel[iglob*3 + 2] = tex1Dfetch(d_accel_tex, iglob*3 + 2) + sum_terms3;
#else
        d_accel[iglob*3]     += sum_terms1;
        d_accel[iglob*3 + 1] += sum_terms2;
        d_accel[iglob*3 + 2] += sum_terms3;
#endif // USE_TEXTURES_FIELDS

      }
      else {

        // for testing purposes only: w/out atomic updates
        //d_accel[iglob*3] -= (0.00000001f*tempx1l + 0.00000001f*tempx2l + 0.00000001f*tempx3l);
        //d_accel[iglob*3 + 1] -= (0.00000001f*tempy1l + 0.00000001f*tempy2l + 0.00000001f*tempy3l);
        //d_accel[iglob*3 + 2] -= (0.00000001f*tempz1l + 0.00000001f*tempz2l + 0.00000001f*tempz3l);

        atomicAdd(&d_accel[iglob*3], sum_terms1);
        atomicAdd(&d_accel[iglob*3+1], sum_terms2);
        atomicAdd(&d_accel[iglob*3+2], sum_terms3);

      } // if(use_mesh_coloring_gpu)

#endif // MESH_COLORING


      // update memory variables based upon the Runge-Kutta scheme
      if( ATTENUATION ){
        compute_element_att_memory(tx,working_element,NSPEC,
                                  d_muv,
                                  factor_common,alphaval,betaval,gammaval,
                                  R_xx,R_yy,R_xy,R_xz,R_yz,
                                  epsilondev_xx,epsilondev_yy,epsilondev_xy,epsilondev_xz,epsilondev_yz,
                                  epsilondev_xx_loc,epsilondev_yy_loc,epsilondev_xy_loc,epsilondev_xz_loc,epsilondev_yz_loc);
      }

      // save deviatoric strain for Runge-Kutta scheme
      if( COMPUTE_AND_STORE_STRAIN ){
        int ijk_ispec = tx + working_element*NGLL3;

        // fortran: epsilondev_xx(:,:,:,ispec) = epsilondev_xx_loc(:,:,:)
        epsilondev_xx[ijk_ispec] = epsilondev_xx_loc;
        epsilondev_yy[ijk_ispec] = epsilondev_yy_loc;
        epsilondev_xy[ijk_ispec] = epsilondev_xy_loc;
        epsilondev_xz[ijk_ispec] = epsilondev_xz_loc;
        epsilondev_yz[ijk_ispec] = epsilondev_yz_loc;
      }

    } // if(active)

} // kernel_2_impl()

/* ----------------------------------------------------------------------------------------------- */

void Kernel_2(int nb_blocks_to_compute,Mesh* mp,int d_iphase,
              int COMPUTE_AND_STORE_STRAIN,int SIMULATION_TYPE,
              int ATTENUATION,int ANISOTROPY,
              int* d_ibool,
              realw* d_xix,
              realw* d_xiy,
              realw* d_xiz,
              realw* d_etax,
              realw* d_etay,
              realw* d_etaz,
              realw* d_gammax,
              realw* d_gammay,
              realw* d_gammaz,
              realw* d_kappav,
              realw* d_muv,
              realw* d_epsilondev_xx,
              realw* d_epsilondev_yy,
              realw* d_epsilondev_xy,
              realw* d_epsilondev_xz,
              realw* d_epsilondev_yz,
              realw* d_epsilon_trace_over_3,
              realw* d_one_minus_sum_beta,
              realw* d_factor_common,
              realw* d_R_xx,
              realw* d_R_yy,
              realw* d_R_xy,
              realw* d_R_xz,
              realw* d_R_yz,
              realw* d_b_epsilondev_xx,
              realw* d_b_epsilondev_yy,
              realw* d_b_epsilondev_xy,
              realw* d_b_epsilondev_xz,
              realw* d_b_epsilondev_yz,
              realw* d_b_epsilon_trace_over_3,
              realw* d_b_R_xx,
              realw* d_b_R_yy,
              realw* d_b_R_xy,
              realw* d_b_R_xz,
              realw* d_b_R_yz,
              realw* d_c11store,
              realw* d_c12store,
              realw* d_c13store,
              realw* d_c14store,
              realw* d_c15store,
              realw* d_c16store,
              realw* d_c22store,
              realw* d_c23store,
              realw* d_c24store,
              realw* d_c25store,
              realw* d_c26store,
              realw* d_c33store,
              realw* d_c34store,
              realw* d_c35store,
              realw* d_c36store,
              realw* d_c44store,
              realw* d_c45store,
              realw* d_c46store,
              realw* d_c55store,
              realw* d_c56store,
              realw* d_c66store,
              realw* d_rhostore){

#ifdef ENABLE_VERY_SLOW_ERROR_CHECKING
  exit_on_cuda_error("before kernel Kernel 2");
#endif

  /* if the grid can handle the number of blocks, we let it be 1D */
  /* grid_2_x = nb_elem_color; */
  /* nb_elem_color is just how many blocks we are computing now */

  int num_blocks_x = nb_blocks_to_compute;
  int num_blocks_y = 1;
  while(num_blocks_x > 65535) {
    num_blocks_x = (int) ceil(num_blocks_x*0.5f);
    num_blocks_y = num_blocks_y*2;
  }

  int blocksize = NGLL3_PADDED;
  dim3 grid(num_blocks_x,num_blocks_y);
  dim3 threads(blocksize,1,1);

  // Cuda timing
  // cudaEvent_t start, stop;
  // realw time;
  // cudaEventCreate(&start);
  // cudaEventCreate(&stop);
  // cudaEventRecord( start, 0 );

  Kernel_2_impl<<<grid,threads,0,mp->compute_stream>>>(nb_blocks_to_compute,
                                  mp->NGLOB_AB,
                                  d_ibool,
                                  mp->d_phase_ispec_inner_elastic,
                                  mp->num_phase_ispec_elastic,
                                  d_iphase,
                                  mp->use_mesh_coloring_gpu,
                                  mp->d_displ, mp->d_accel,
                                  d_xix, d_xiy, d_xiz,
                                  d_etax, d_etay, d_etaz,
                                  d_gammax, d_gammay, d_gammaz,
                                  d_kappav, d_muv,
                                  COMPUTE_AND_STORE_STRAIN,
                                  d_epsilondev_xx,
                                  d_epsilondev_yy,
                                  d_epsilondev_xy,
                                  d_epsilondev_xz,
                                  d_epsilondev_yz,
                                  d_epsilon_trace_over_3,
                                  SIMULATION_TYPE,
                                  ATTENUATION,mp->NSPEC_AB,
                                  d_one_minus_sum_beta,
                                  d_factor_common,
                                  d_R_xx,d_R_yy,d_R_xy,d_R_xz,d_R_yz,
                                  mp->d_alphaval,mp->d_betaval,mp->d_gammaval,
                                  ANISOTROPY,
                                  d_c11store,
                                  d_c12store,
                                  d_c13store,
                                  d_c14store,
                                  d_c15store,
                                  d_c16store,
                                  d_c22store,
                                  d_c23store,
                                  d_c24store,
                                  d_c25store,
                                  d_c26store,
                                  d_c33store,
                                  d_c34store,
                                  d_c35store,
                                  d_c36store,
                                  d_c44store,
                                  d_c45store,
                                  d_c46store,
                                  d_c55store,
                                  d_c56store,
                                  d_c66store,
                                  mp->gravity,
                                  mp->d_minus_g,
                                  mp->d_minus_deriv_gravity,
                                  d_rhostore,
                                  mp->d_wgll_cube);


  if(SIMULATION_TYPE == 3) {
    Kernel_2_impl<<< grid,threads,0,mp->compute_stream>>>(nb_blocks_to_compute,
                                     mp->NGLOB_AB,
                                     d_ibool,
                                     mp->d_phase_ispec_inner_elastic,
                                     mp->num_phase_ispec_elastic,
                                     d_iphase,
                                     mp->use_mesh_coloring_gpu,
                                     mp->d_b_displ, mp->d_b_accel,
                                     d_xix, d_xiy, d_xiz,
                                     d_etax, d_etay, d_etaz,
                                     d_gammax, d_gammay, d_gammaz,
                                     d_kappav, d_muv,
                                     COMPUTE_AND_STORE_STRAIN,
                                     d_b_epsilondev_xx,
                                     d_b_epsilondev_yy,
                                     d_b_epsilondev_xy,
                                     d_b_epsilondev_xz,
                                     d_b_epsilondev_yz,
                                     d_b_epsilon_trace_over_3,
                                     SIMULATION_TYPE,
                                     ATTENUATION,mp->NSPEC_AB,
                                     d_one_minus_sum_beta,
                                     d_factor_common,
                                     d_b_R_xx,d_b_R_yy,d_b_R_xy,d_b_R_xz,d_b_R_yz,
                                     mp->d_b_alphaval,mp->d_b_betaval,mp->d_b_gammaval,
                                     ANISOTROPY,
                                     d_c11store,
                                     d_c12store,
                                     d_c13store,
                                     d_c14store,
                                     d_c15store,
                                     d_c16store,
                                     d_c22store,
                                     d_c23store,
                                     d_c24store,
                                     d_c25store,
                                     d_c26store,
                                     d_c33store,
                                     d_c34store,
                                     d_c35store,
                                     d_c36store,
                                     d_c44store,
                                     d_c45store,
                                     d_c46store,
                                     d_c55store,
                                     d_c56store,
                                     d_c66store,
                                     mp->gravity,
                                     mp->d_minus_g,
                                     mp->d_minus_deriv_gravity,
                                     d_rhostore,
                                     mp->d_wgll_cube);
  }

  // cudaEventRecord( stop, 0 );
  // cudaEventSynchronize( stop );
  // cudaEventElapsedTime( &time, start, stop );
  // cudaEventDestroy( start );
  // cudaEventDestroy( stop );
  // printf("Kernel2 Execution Time: %f ms\n",time);

  /* cudaThreadSynchronize(); */
  /* LOG("Kernel 2 finished"); */
#ifdef ENABLE_VERY_SLOW_ERROR_CHECKING
  exit_on_cuda_error("Kernel_2_impl ");
#endif
}

/* ----------------------------------------------------------------------------------------------- */


extern "C"
void FC_FUNC_(compute_forces_elastic_cuda,
              COMPUTE_FORCES_ELASTIC_CUDA)(long* Mesh_pointer_f,
                                           int* iphase,
                                           int* nspec_outer_elastic,
                                           int* nspec_inner_elastic,
                                           int* SIMULATION_TYPE,
                                           int* COMPUTE_AND_STORE_STRAIN,
                                           int* ATTENUATION,
                                           int* ANISOTROPY) {

  TRACE("compute_forces_elastic_cuda");
  // EPIK_TRACER("compute_forces_elastic_cuda");
  //printf("Running compute_forces\n");
  //double start_time = get_time();

  Mesh* mp = (Mesh*)(*Mesh_pointer_f); // get Mesh from fortran integer wrapper

  int num_elements;

  if( *iphase == 1 )
    num_elements = *nspec_outer_elastic;
  else
    num_elements = *nspec_inner_elastic;

  // checks if anything to do
  if( num_elements == 0 ) return;

  // mesh coloring
  if( mp->use_mesh_coloring_gpu ){

    // note: array offsets require sorted arrays, such that e.g. ibool starts with elastic elements
    //         and followed by acoustic ones.
    //         elastic elements also start with outer than inner element ordering

    int nb_colors,nb_blocks_to_compute;
    int istart;
    int color_offset,color_offset_nonpadded,color_offset_nonpadded_att2;

    // sets up color loop
    if( *iphase == 1 ){
      // outer elements
      nb_colors = mp->num_colors_outer_elastic;
      istart = 0;

      // array offsets
      color_offset = 0;
      color_offset_nonpadded = 0;
      color_offset_nonpadded_att2 = 0;
    }else{
      // inner elements (start after outer elements)
      nb_colors = mp->num_colors_outer_elastic + mp->num_colors_inner_elastic;
      istart = mp->num_colors_outer_elastic;

      // array offsets
      color_offset = (*nspec_outer_elastic) * NGLL3_PADDED;
      color_offset_nonpadded = (*nspec_outer_elastic) * NGLL3;
      color_offset_nonpadded_att2 = (*nspec_outer_elastic) * NGLL3 * N_SLS;
    }

    // loops over colors
    for(int icolor = istart; icolor < nb_colors; icolor++){

      nb_blocks_to_compute = mp->h_num_elem_colors_elastic[icolor];

      // checks
      //if( nb_blocks_to_compute <= 0 ){
      //  printf("error number of elastic color blocks: %d -- color = %d \n",nb_blocks_to_compute,icolor);
      //  exit(EXIT_FAILURE);
      //}

      Kernel_2(nb_blocks_to_compute,mp,*iphase,
               *COMPUTE_AND_STORE_STRAIN,*SIMULATION_TYPE,
               *ATTENUATION,*ANISOTROPY,
               mp->d_ibool + color_offset_nonpadded,
               mp->d_xix + color_offset,
               mp->d_xiy + color_offset,
               mp->d_xiz + color_offset,
               mp->d_etax + color_offset,
               mp->d_etay + color_offset,
               mp->d_etaz + color_offset,
               mp->d_gammax + color_offset,
               mp->d_gammay + color_offset,
               mp->d_gammaz + color_offset,
               mp->d_kappav + color_offset,
               mp->d_muv + color_offset,
               mp->d_epsilondev_xx + color_offset_nonpadded,
               mp->d_epsilondev_yy + color_offset_nonpadded,
               mp->d_epsilondev_xy + color_offset_nonpadded,
               mp->d_epsilondev_xz + color_offset_nonpadded,
               mp->d_epsilondev_yz + color_offset_nonpadded,
               mp->d_epsilon_trace_over_3 + color_offset_nonpadded,
               mp->d_one_minus_sum_beta + color_offset_nonpadded,
               mp->d_factor_common + color_offset_nonpadded_att2,
               mp->d_R_xx + color_offset_nonpadded,
               mp->d_R_yy + color_offset_nonpadded,
               mp->d_R_xy + color_offset_nonpadded,
               mp->d_R_xz + color_offset_nonpadded,
               mp->d_R_yz + color_offset_nonpadded,
               mp->d_b_epsilondev_xx + color_offset_nonpadded,
               mp->d_b_epsilondev_yy + color_offset_nonpadded,
               mp->d_b_epsilondev_xy + color_offset_nonpadded,
               mp->d_b_epsilondev_xz + color_offset_nonpadded,
               mp->d_b_epsilondev_yz + color_offset_nonpadded,
               mp->d_b_epsilon_trace_over_3 + color_offset_nonpadded,
               mp->d_b_R_xx + color_offset_nonpadded,
               mp->d_b_R_yy + color_offset_nonpadded,
               mp->d_b_R_xy + color_offset_nonpadded,
               mp->d_b_R_xz + color_offset_nonpadded,
               mp->d_b_R_yz + color_offset_nonpadded,
               mp->d_c11store + color_offset,
               mp->d_c12store + color_offset,
               mp->d_c13store + color_offset,
               mp->d_c14store + color_offset,
               mp->d_c15store + color_offset,
               mp->d_c16store + color_offset,
               mp->d_c22store + color_offset,
               mp->d_c23store + color_offset,
               mp->d_c24store + color_offset,
               mp->d_c25store + color_offset,
               mp->d_c26store + color_offset,
               mp->d_c33store + color_offset,
               mp->d_c34store + color_offset,
               mp->d_c35store + color_offset,
               mp->d_c36store + color_offset,
               mp->d_c44store + color_offset,
               mp->d_c45store + color_offset,
               mp->d_c46store + color_offset,
               mp->d_c55store + color_offset,
               mp->d_c56store + color_offset,
               mp->d_c66store + color_offset,
               mp->d_rhostore + color_offset);

      // for padded and aligned arrays
      color_offset += nb_blocks_to_compute * NGLL3_PADDED;
      // for no-aligned arrays
      color_offset_nonpadded += nb_blocks_to_compute * NGLL3;
      // for factor_common array
      color_offset_nonpadded_att2 += nb_blocks_to_compute * NGLL3 * N_SLS;
    }

  }else{

    // no mesh coloring: uses atomic updates

    Kernel_2(num_elements,mp,*iphase,
             *COMPUTE_AND_STORE_STRAIN,*SIMULATION_TYPE,
             *ATTENUATION,*ANISOTROPY,
             mp->d_ibool,
             mp->d_xix,
             mp->d_xiy,
             mp->d_xiz,
             mp->d_etax,
             mp->d_etay,
             mp->d_etaz,
             mp->d_gammax,
             mp->d_gammay,
             mp->d_gammaz,
             mp->d_kappav,
             mp->d_muv,
             mp->d_epsilondev_xx,
             mp->d_epsilondev_yy,
             mp->d_epsilondev_xy,
             mp->d_epsilondev_xz,
             mp->d_epsilondev_yz,
             mp->d_epsilon_trace_over_3,
             mp->d_one_minus_sum_beta,
             mp->d_factor_common,
             mp->d_R_xx,
             mp->d_R_yy,
             mp->d_R_xy,
             mp->d_R_xz,
             mp->d_R_yz,
             mp->d_b_epsilondev_xx,
             mp->d_b_epsilondev_yy,
             mp->d_b_epsilondev_xy,
             mp->d_b_epsilondev_xz,
             mp->d_b_epsilondev_yz,
             mp->d_b_epsilon_trace_over_3,
             mp->d_b_R_xx,
             mp->d_b_R_yy,
             mp->d_b_R_xy,
             mp->d_b_R_xz,
             mp->d_b_R_yz,
             mp->d_c11store,
             mp->d_c12store,
             mp->d_c13store,
             mp->d_c14store,
             mp->d_c15store,
             mp->d_c16store,
             mp->d_c22store,
             mp->d_c23store,
             mp->d_c24store,
             mp->d_c25store,
             mp->d_c26store,
             mp->d_c33store,
             mp->d_c34store,
             mp->d_c35store,
             mp->d_c36store,
             mp->d_c44store,
             mp->d_c45store,
             mp->d_c46store,
             mp->d_c55store,
             mp->d_c56store,
             mp->d_c66store,
             mp->d_rhostore);
  }

  //daniel: todo - check with routine sync_copy_from_device below...
//  // Wait until async-memcpy of outer elements is finished and start MPI.
//  if(*iphase==2) {
//    cudaStreamSynchronize(mp->copy_stream);
//
//    // There have been problems using the pinned-memory with MPI, so
//    // we copy the buffer into a non-pinned region.
//    memcpy(mp->send_buffer,mp->h_send_accel_buffer,
//           mp->size_mpi_send_buffer*sizeof(float));
//
//    // memory copy is now finished, so non-blocking MPI send can proceed
//    // MPI based halo exchange
//
//    assemble_mpi_vector_send_cuda_(&(mp->NPROCS),
//                                   mp->send_buffer, /* "regular" memory */
//                                   // mp->h_send_accel_buffer, /* pinned memory **CRASH** */
//                                   mp->buffer_recv_vector_ext_mesh,
//                                   &mp->num_interfaces_ext_mesh,
//                                   &mp->max_nibool_interfaces_ext_mesh,
//                                   mp->nibool_interfaces_ext_mesh,
//                                   mp->my_neighbours_ext_mesh,
//                                   mp->request_send_vector_ext_mesh,
//                                   mp->request_recv_vector_ext_mesh);
//
//    // Decided to keep launching kernels and to wait for MPI & do memcpy while other kernels launch.
//    // cudaDeviceSynchronize();
//  }

}

/* ----------------------------------------------------------------------------------------------- */

//daniel: todo - use this instead above call to fortran routine to avoid compilation problems
extern "C"
void FC_FUNC_(sync_copy_from_device,
              SYNC_copy_FROM_DEVICE)(long* Mesh_pointer_f,
                                     int* iphase,
                                     realw* send_buffer) {

  TRACE("sync_copy_from_device");

  Mesh* mp = (Mesh*)(*Mesh_pointer_f); // get Mesh from fortran integer wrapper

  // Wait until async-memcpy of outer elements is finished and start MPI.
  if( *iphase != 2 ){ exit_on_cuda_error("sync_copy_from_device must be called for iphase == 2"); }

  //if(*iphase==2) {

  // waits for asynchronous copy to finish
  cudaStreamSynchronize(mp->copy_stream);

  // There have been problems using the pinned-memory with MPI, so
  // we copy the buffer into a non-pinned region.
  memcpy(send_buffer,mp->h_send_accel_buffer,
         mp->size_mpi_send_buffer*sizeof(float));

  // memory copy is now finished, so non-blocking MPI send can proceed

  //}
}

/* ----------------------------------------------------------------------------------------------- */

// KERNEL 3

/* ----------------------------------------------------------------------------------------------- */

__global__ void kernel_3_cuda_device(realw* veloc,
                                     realw* accel, int size,
                                     realw deltatover2,
                                     realw* rmass) {
  int id = threadIdx.x + blockIdx.x*blockDim.x + blockIdx.y*gridDim.x*blockDim.x;

  /* because of block and grid sizing problems, there is a small */
  /* amount of buffer at the end of the calculation */
  if(id < size) {
    realw new_accel = accel[id] * rmass[id / 3];
    veloc[id] += deltatover2 * new_accel;
    accel[id] = new_accel;
/*
    accel[3*id] = accel[3*id]*rmass[id];
    accel[3*id+1] = accel[3*id+1]*rmass[id];
    accel[3*id+2] = accel[3*id+2]*rmass[id];

    veloc[3*id] = veloc[3*id] + deltatover2*accel[3*id];
    veloc[3*id+1] = veloc[3*id+1] + deltatover2*accel[3*id+1];
    veloc[3*id+2] = veloc[3*id+2] + deltatover2*accel[3*id+2];
*/
  }
}

/* ----------------------------------------------------------------------------------------------- */

__global__ void kernel_3_accel_cuda_device(realw* accel,
                                           int size,
                                           realw* rmass) {
  int id = threadIdx.x + blockIdx.x*blockDim.x + blockIdx.y*gridDim.x*blockDim.x;

  /* because of block and grid sizing problems, there is a small */
  /* amount of buffer at the end of the calculation */
  if(id < size) {
    accel[id] *= rmass[id / 3];
/*
    accel[3*id] = accel[3*id]*rmass[id];
    accel[3*id+1] = accel[3*id+1]*rmass[id];
    accel[3*id+2] = accel[3*id+2]*rmass[id];
*/
  }
}

/* ----------------------------------------------------------------------------------------------- */

__global__ void kernel_3_veloc_cuda_device(realw* veloc,
                                           realw* accel,
                                           int size,
                                           realw deltatover2) {
  int id = threadIdx.x + blockIdx.x*blockDim.x + blockIdx.y*gridDim.x*blockDim.x;

  /* because of block and grid sizing problems, there is a small */
  /* amount of buffer at the end of the calculation */
  if(id < size) {
    veloc[3*id] = veloc[3*id] + deltatover2*accel[3*id];
    veloc[3*id+1] = veloc[3*id+1] + deltatover2*accel[3*id+1];
    veloc[3*id+2] = veloc[3*id+2] + deltatover2*accel[3*id+2];
  }
}

/* ----------------------------------------------------------------------------------------------- */

extern "C"
void FC_FUNC_(kernel_3_a_cuda,
              KERNEL_3_A_CUDA)(long* Mesh_pointer,
                               int* size_F,
                               realw* deltatover2_F,
                               int* SIMULATION_TYPE_f,
                               realw* b_deltatover2_F,
                               int* OCEANS) {
TRACE("kernel_3_a_cuda");

   Mesh* mp = (Mesh*)(*Mesh_pointer); // get Mesh from fortran integer wrapper
   //int size = *size_F;
   int size = *size_F * 3;
   int SIMULATION_TYPE = *SIMULATION_TYPE_f;
   realw deltatover2 = *deltatover2_F;
   realw b_deltatover2 = *b_deltatover2_F;

   int blocksize = BLOCKSIZE_KERNEL3;
   int size_padded = ((int)ceil(((double)size)/((double)blocksize)))*blocksize;

   int num_blocks_x = size_padded/blocksize;
   int num_blocks_y = 1;
   while(num_blocks_x > 65535) {
     num_blocks_x = (int) ceil(num_blocks_x*0.5f);
     num_blocks_y = num_blocks_y*2;
   }

   dim3 grid(num_blocks_x,num_blocks_y);
   dim3 threads(blocksize,1,1);

   // check whether we can update accel and veloc, or only accel at this point
   if( *OCEANS == 0 ){
     // updates both, accel and veloc
     kernel_3_cuda_device<<< grid, threads,0,mp->compute_stream>>>(mp->d_veloc, mp->d_accel, size, deltatover2, mp->d_rmass);

     if(SIMULATION_TYPE == 3) {
       kernel_3_cuda_device<<< grid, threads,0,mp->compute_stream>>>(mp->d_b_veloc, mp->d_b_accel, size, b_deltatover2,mp->d_rmass);
     }
   }else{
     // updates only accel
     kernel_3_accel_cuda_device<<< grid, threads,0,mp->compute_stream>>>(mp->d_accel, size, mp->d_rmass);

     if(SIMULATION_TYPE == 3) {
       kernel_3_accel_cuda_device<<< grid, threads,0,mp->compute_stream>>>(mp->d_b_accel, size, mp->d_rmass);
     }
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
                             int* size_F,
                             realw* deltatover2_F,
                             int* SIMULATION_TYPE_f,
                             realw* b_deltatover2_F) {
  TRACE("kernel_3_b_cuda");

  Mesh* mp = (Mesh*)(*Mesh_pointer); // get Mesh from fortran integer wrapper
  int size = *size_F;
  int SIMULATION_TYPE = *SIMULATION_TYPE_f;
  realw deltatover2 = *deltatover2_F;
  realw b_deltatover2 = *b_deltatover2_F;

  int blocksize = BLOCKSIZE_KERNEL3;
  int size_padded = ((int)ceil(((double)size)/((double)blocksize)))*blocksize;

  int num_blocks_x = size_padded/blocksize;
  int num_blocks_y = 1;
  while(num_blocks_x > 65535) {
    num_blocks_x = (int) ceil(num_blocks_x*0.5f);
    num_blocks_y = num_blocks_y*2;
  }

  dim3 grid(num_blocks_x,num_blocks_y);
  dim3 threads(blocksize,1,1);

  // updates only veloc at this point
  kernel_3_veloc_cuda_device<<< grid, threads,0,mp->compute_stream>>>(mp->d_veloc,mp->d_accel,size,deltatover2);

  if(SIMULATION_TYPE == 3) {
    kernel_3_veloc_cuda_device<<< grid, threads,0,mp->compute_stream>>>(mp->d_b_veloc,mp->d_b_accel,size,b_deltatover2);
  }

#ifdef ENABLE_VERY_SLOW_ERROR_CHECKING
  //printf("checking updatedispl_kernel launch...with %dx%d blocks\n",num_blocks_x,num_blocks_y);
  exit_on_cuda_error("after kernel 3 b");
#endif
}


/* ----------------------------------------------------------------------------------------------- */

/* OCEANS load on free surface */

/* ----------------------------------------------------------------------------------------------- */


__global__ void elastic_ocean_load_cuda_kernel(realw* accel,
                                               realw* rmass,
                                               realw* rmass_ocean_load,
                                               int num_free_surface_faces,
                                               int* free_surface_ispec,
                                               int* free_surface_ijk,
                                               realw* free_surface_normal,
                                               int* ibool,
                                               int* updated_dof_ocean_load) {
  // gets spectral element face id
  int igll = threadIdx.x ;  //  threadIdx.y*blockDim.x will be always = 0 for thread block (25,1,1)
  int iface = blockIdx.x + gridDim.x*blockIdx.y;
  realw nx,ny,nz;
  realw force_normal_comp,additional_term;

  // for all faces on free surface
  if( iface < num_free_surface_faces ){

    int ispec = free_surface_ispec[iface]-1;

    // gets global point index
    int i = free_surface_ijk[INDEX3(NDIM,NGLL2,0,igll,iface)] - 1; // (1,igll,iface)
    int j = free_surface_ijk[INDEX3(NDIM,NGLL2,1,igll,iface)] - 1;
    int k = free_surface_ijk[INDEX3(NDIM,NGLL2,2,igll,iface)] - 1;

    int iglob = ibool[INDEX4(5,5,5,i,j,k,ispec)] - 1;

    //if(igll == 0 ) printf("igll %d %d %d %d\n",igll,i,j,k,iglob);

    // only update this global point once

    // daniel: TODO - there might be better ways to implement a mutex like below,
    //            and find a workaround to not use the temporary update array.
    //            atomicExch: returns the old value, i.e. 0 indicates that we still have to do this point

    if( atomicExch(&updated_dof_ocean_load[iglob],1) == 0){

      // get normal
      nx = free_surface_normal[INDEX3(NDIM,NGLL2,0,igll,iface)]; //(1,igll,iface)
      ny = free_surface_normal[INDEX3(NDIM,NGLL2,1,igll,iface)];
      nz = free_surface_normal[INDEX3(NDIM,NGLL2,2,igll,iface)];

      // make updated component of right-hand side
      // we divide by rmass() which is 1 / M
      // we use the total force which includes the Coriolis term above
      force_normal_comp = ( accel[iglob*3]*nx + accel[iglob*3+1]*ny + accel[iglob*3+2]*nz ) / rmass[iglob];

      additional_term = (rmass_ocean_load[iglob] - rmass[iglob]) * force_normal_comp;

      // probably wouldn't need atomicAdd anymore, but just to be sure...
      atomicAdd(&accel[iglob*3], + additional_term * nx);
      atomicAdd(&accel[iglob*3+1], + additional_term * ny);
      atomicAdd(&accel[iglob*3+2], + additional_term * nz);
    }
  }
}

/* ----------------------------------------------------------------------------------------------- */

extern "C"
void FC_FUNC_(elastic_ocean_load_cuda,
              ELASTIC_OCEAN_LOAD_CUDA)(long* Mesh_pointer_f,
                                       int* SIMULATION_TYPE) {

  TRACE("elastic_ocean_load_cuda");

  Mesh* mp = (Mesh*)(*Mesh_pointer_f); //get mesh pointer out of fortran integer container

  // checks if anything to do
  if( mp->num_free_surface_faces == 0 ) return;

  // block sizes: exact blocksize to match NGLLSQUARE
  int blocksize = NGLL2;

  int num_blocks_x = mp->num_free_surface_faces;
  int num_blocks_y = 1;
  while(num_blocks_x > 65535) {
    num_blocks_x = (int) ceil(num_blocks_x*0.5f);
    num_blocks_y = num_blocks_y*2;
  }

  dim3 grid(num_blocks_x,num_blocks_y);
  dim3 threads(blocksize,1,1);


  // initializes temporary array to zero
  print_CUDA_error_if_any(cudaMemset(mp->d_updated_dof_ocean_load,0,
                                     sizeof(int)*mp->NGLOB_AB),88501);

#ifdef ENABLE_VERY_SLOW_ERROR_CHECKING
  exit_on_cuda_error("before kernel elastic_ocean_load_cuda");
#endif

  elastic_ocean_load_cuda_kernel<<<grid,threads,0,mp->compute_stream>>>(mp->d_accel,
                                                   mp->d_rmass,
                                                   mp->d_rmass_ocean_load,
                                                   mp->num_free_surface_faces,
                                                   mp->d_free_surface_ispec,
                                                   mp->d_free_surface_ijk,
                                                   mp->d_free_surface_normal,
                                                   mp->d_ibool,
                                                   mp->d_updated_dof_ocean_load);
  // for backward/reconstructed potentials
  if(*SIMULATION_TYPE == 3) {
    // re-initializes array
    print_CUDA_error_if_any(cudaMemset(mp->d_updated_dof_ocean_load,0,
                                       sizeof(int)*mp->NGLOB_AB),88502);

    elastic_ocean_load_cuda_kernel<<<grid,threads,0,mp->compute_stream>>>(mp->d_b_accel,
                                                     mp->d_rmass,
                                                     mp->d_rmass_ocean_load,
                                                     mp->num_free_surface_faces,
                                                     mp->d_free_surface_ispec,
                                                     mp->d_free_surface_ijk,
                                                     mp->d_free_surface_normal,
                                                     mp->d_ibool,
                                                     mp->d_updated_dof_ocean_load);

  }


#ifdef ENABLE_VERY_SLOW_ERROR_CHECKING
  exit_on_cuda_error("elastic_ocean_load_cuda");
#endif
}


/* ----------------------------------------------------------------------------------------------- */

/* note:
 constant arrays when used in compute_forces_acoustic_cuda.cu routines stay zero,
 constant declaration and cudaMemcpyToSymbol would have to be in the same file...

 extern keyword doesn't work for __constant__ declarations.

 also:
 cudaMemcpyToSymbol("deviceCaseParams", caseParams, sizeof(CaseParams));
 ..
 and compile with -arch=sm_20

 see also: http://stackoverflow.com/questions/4008031/how-to-use-cuda-constant-memory-in-a-programmer-pleasant-way
 doesn't seem to work.

 we could keep arrays separated for acoustic and elastic routines...

 for now, we store pointers with cudaGetSymbolAddress() function calls.

 */


// constant arrays

void setConst_hprime_xx(realw* array,Mesh* mp)
{

  cudaError_t err = cudaMemcpyToSymbol(d_hprime_xx, array, NGLL2*sizeof(realw));
  if (err != cudaSuccess)
  {
    fprintf(stderr, "Error in setConst_hprime_xx: %s\n", cudaGetErrorString(err));
    fprintf(stderr, "The problem is maybe -arch sm_13 instead of -arch sm_11 in the Makefile, please doublecheck\n");
    exit(1);
  }

  err = cudaGetSymbolAddress((void**)&(mp->d_hprime_xx),"d_hprime_xx");
  if(err != cudaSuccess) {
    fprintf(stderr, "Error with d_hprime_xx: %s\n", cudaGetErrorString(err));
    exit(1);
  }
}

// void setConst_hprime_yy(realw* array,Mesh* mp)
// {

//   cudaError_t err = cudaMemcpyToSymbol(d_hprime_yy, array, NGLL2*sizeof(realw));
//   if (err != cudaSuccess)
//   {
//     fprintf(stderr, "Error in setConst_hprime_yy: %s\n", cudaGetErrorString(err));
//     fprintf(stderr, "The problem is maybe -arch sm_13 instead of -arch sm_11 in the Makefile, please doublecheck\n");
//     exit(1);
//   }

//   err = cudaGetSymbolAddress((void**)&(mp->d_hprime_yy),"d_hprime_yy");
//   if(err != cudaSuccess) {
//     fprintf(stderr, "Error with d_hprime_yy: %s\n", cudaGetErrorString(err));
//     exit(1);
//   }
// }

// void setConst_hprime_zz(realw* array,Mesh* mp)
// {

//   cudaError_t err = cudaMemcpyToSymbol(d_hprime_zz, array, NGLL2*sizeof(realw));
//   if (err != cudaSuccess)
//   {
//     fprintf(stderr, "Error in setConst_hprime_zz: %s\n", cudaGetErrorString(err));
//     fprintf(stderr, "The problem is maybe -arch sm_13 instead of -arch sm_11 in the Makefile, please doublecheck\n");
//     exit(1);
//   }

//   err = cudaGetSymbolAddress((void**)&(mp->d_hprime_zz),"d_hprime_zz");
//   if(err != cudaSuccess) {
//     fprintf(stderr, "Error with d_hprime_zz: %s\n", cudaGetErrorString(err));
//     exit(1);
//   }
// }


void setConst_hprimewgll_xx(realw* array,Mesh* mp)
{
  cudaError_t err = cudaMemcpyToSymbol(d_hprimewgll_xx, array, NGLL2*sizeof(realw));
  if (err != cudaSuccess)
  {
    fprintf(stderr, "Error in setConst_hprimewgll_xx: %s\n", cudaGetErrorString(err));
    exit(1);
  }

  err = cudaGetSymbolAddress((void**)&(mp->d_hprimewgll_xx),"d_hprimewgll_xx");
  if(err != cudaSuccess) {
    fprintf(stderr, "Error with d_hprimewgll_xx: %s\n", cudaGetErrorString(err));
    exit(1);
  }
}

void setConst_hprimewgll_yy(realw* array,Mesh* mp)
{
  cudaError_t err = cudaMemcpyToSymbol(d_hprimewgll_yy, array, NGLL2*sizeof(realw));
  if (err != cudaSuccess)
  {
    fprintf(stderr, "Error in setConst_hprimewgll_yy: %s\n", cudaGetErrorString(err));
    exit(1);
  }

  err = cudaGetSymbolAddress((void**)&(mp->d_hprimewgll_yy),"d_hprimewgll_yy");
  if(err != cudaSuccess) {
    fprintf(stderr, "Error with d_hprimewgll_yy: %s\n", cudaGetErrorString(err));
    exit(1);
  }
}

void setConst_hprimewgll_zz(realw* array,Mesh* mp)
{
  cudaError_t err = cudaMemcpyToSymbol(d_hprimewgll_zz, array, NGLL2*sizeof(realw));
  if (err != cudaSuccess)
  {
    fprintf(stderr, "Error in setConst_hprimewgll_zz: %s\n", cudaGetErrorString(err));
    exit(1);
  }

  err = cudaGetSymbolAddress((void**)&(mp->d_hprimewgll_zz),"d_hprimewgll_zz");
  if(err != cudaSuccess) {
    fprintf(stderr, "Error with d_hprimewgll_zz: %s\n", cudaGetErrorString(err));
    exit(1);
  }
}

void setConst_wgllwgll_xy(realw* array,Mesh* mp)
{
  cudaError_t err = cudaMemcpyToSymbol(d_wgllwgll_xy, array, NGLL2*sizeof(realw));
  if (err != cudaSuccess)
  {
    fprintf(stderr, "Error in setConst_wgllwgll_xy: %s\n", cudaGetErrorString(err));
    exit(1);
  }
  //mp->d_wgllwgll_xy = d_wgllwgll_xy;
  err = cudaGetSymbolAddress((void**)&(mp->d_wgllwgll_xy),"d_wgllwgll_xy");
  if(err != cudaSuccess) {
    fprintf(stderr, "Error with d_wgllwgll_xy: %s\n", cudaGetErrorString(err));
    exit(1);
  }

}

void setConst_wgllwgll_xz(realw* array,Mesh* mp)
{
  cudaError_t err = cudaMemcpyToSymbol(d_wgllwgll_xz, array, NGLL2*sizeof(realw));
  if (err != cudaSuccess)
  {
    fprintf(stderr, "Error in  setConst_wgllwgll_xz: %s\n", cudaGetErrorString(err));
    exit(1);
  }
  //mp->d_wgllwgll_xz = d_wgllwgll_xz;
  err = cudaGetSymbolAddress((void**)&(mp->d_wgllwgll_xz),"d_wgllwgll_xz");
  if(err != cudaSuccess) {
    fprintf(stderr, "Error with d_wgllwgll_xz: %s\n", cudaGetErrorString(err));
    exit(1);
  }

}

void setConst_wgllwgll_yz(realw* array,Mesh* mp)
{
  cudaError_t err = cudaMemcpyToSymbol(d_wgllwgll_yz, array, NGLL2*sizeof(realw));
  if (err != cudaSuccess)
  {
    fprintf(stderr, "Error in setConst_wgllwgll_yz: %s\n", cudaGetErrorString(err));
    exit(1);
  }
  //mp->d_wgllwgll_yz = d_wgllwgll_yz;
  err = cudaGetSymbolAddress((void**)&(mp->d_wgllwgll_yz),"d_wgllwgll_yz");
  if(err != cudaSuccess) {
    fprintf(stderr, "Error with d_wgllwgll_yz: %s\n", cudaGetErrorString(err));
    exit(1);
  }

}

void setConst_wgll_cube(realw* array,Mesh* mp)
{
  cudaError_t err = cudaMemcpyToSymbol(d_wgll_cube, array, NGLL3*sizeof(realw));
  if (err != cudaSuccess)
  {
    fprintf(stderr, "Error in setConst_wgll_cube: %s\n", cudaGetErrorString(err));
    exit(1);
  }
  //mp->d_wgll_cube = d_wgll_cube;
  err = cudaGetSymbolAddress((void**)&(mp->d_wgll_cube),"d_wgll_cube");
  if(err != cudaSuccess) {
    fprintf(stderr, "Error with d_wgll_cube: %s\n", cudaGetErrorString(err));
    exit(1);
  }

}
