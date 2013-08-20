/*
 !=====================================================================
 !
 !               S p e c f e m 3 D  V e r s i o n  2 . 1
 !               ---------------------------------------
 !
 !          Main authors: Dimitri Komatitsch and Jeroen Tromp
 !    Princeton University, USA and CNRS / INRIA / University of Pau
 ! (c) Princeton University / California Institute of Technology and CNRS / INRIA / University of Pau
 !                             July 2012
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


#ifdef USE_TEXTURES_FIELDS
realw_texture d_displ_tex;
realw_texture d_veloc_tex;
realw_texture d_accel_tex;
//backward/reconstructed
realw_texture d_b_displ_tex;
realw_texture d_b_veloc_tex;
realw_texture d_b_accel_tex;

//note: texture variables are implicitly static, and cannot be passed as arguments to cuda kernels;
//      thus, 1) we thus use if-statements (FORWARD_OR_ADJOINT) to determine from which texture to fetch from
//            2) we use templates
//      since if-statements are a bit slower as the variable is only known at runtime, we use option 2)

// templates definitions
template<int FORWARD_OR_ADJOINT> __device__ float texfetch_displ(int x);
template<int FORWARD_OR_ADJOINT> __device__ float texfetch_veloc(int x);
template<int FORWARD_OR_ADJOINT> __device__ float texfetch_accel(int x);

// templates for texture fetching
// FORWARD_OR_ADJOINT == 1 <- forward arrays
template<> __device__ float texfetch_displ<1>(int x) { return tex1Dfetch(d_displ_tex, x); }
template<> __device__ float texfetch_veloc<1>(int x) { return tex1Dfetch(d_veloc_tex, x); }
template<> __device__ float texfetch_accel<1>(int x) { return tex1Dfetch(d_accel_tex, x); }
// FORWARD_OR_ADJOINT == 3 <- backward/reconstructed arrays
template<> __device__ float texfetch_displ<3>(int x) { return tex1Dfetch(d_b_displ_tex, x); }
template<> __device__ float texfetch_veloc<3>(int x) { return tex1Dfetch(d_b_veloc_tex, x); }
template<> __device__ float texfetch_accel<3>(int x) { return tex1Dfetch(d_b_accel_tex, x); }

#endif

#ifdef USE_TEXTURES_CONSTANTS
realw_texture d_hprime_xx_tex;
#endif


/* ----------------------------------------------------------------------------------------------- */

// prepares a device array with with all inter-element edge-nodes -- this
// is followed by a memcpy and MPI operations

__global__ void prepare_boundary_accel_on_device(realw* d_accel, realw* d_send_accel_buffer,
                                                 int num_interfaces_ext_mesh,
                                                 int max_nibool_interfaces_ext_mesh,
                                                 int* d_nibool_interfaces_ext_mesh,
                                                 int* d_ibool_interfaces_ext_mesh) {

  int id = threadIdx.x + blockIdx.x*blockDim.x + blockIdx.y*gridDim.x*blockDim.x;
  int ientry,iglob;

  for( int iinterface=0; iinterface < num_interfaces_ext_mesh; iinterface++) {
    if( id < d_nibool_interfaces_ext_mesh[iinterface] ) {

      // entry in interface array
      ientry = id + max_nibool_interfaces_ext_mesh*iinterface;
      // global index in wavefield
      iglob = d_ibool_interfaces_ext_mesh[ientry] - 1;

      d_send_accel_buffer[3*ientry] = d_accel[3*iglob];
      d_send_accel_buffer[3*ientry + 1 ] = d_accel[3*iglob + 1];
      d_send_accel_buffer[3*ientry + 2 ] = d_accel[3*iglob + 2];
    }
  }

}

/* ----------------------------------------------------------------------------------------------- */

// prepares and transfers the inter-element edge-nodes to the host to be MPI'd
// (elements on boundary)
extern "C"
void FC_FUNC_(transfer_boun_accel_from_device,
              TRANSFER_BOUN_ACCEL_FROM_DEVICE)(long* Mesh_pointer,
                                               realw* accel,
                                               realw* send_accel_buffer,
                                               int* FORWARD_OR_ADJOINT){
TRACE("\ttransfer_boun_accel_from_device");

  Mesh* mp = (Mesh*)(*Mesh_pointer); //get mesh pointer out of fortran integer container

  // checks if anything to do
  if( mp->size_mpi_buffer > 0 ){

    int blocksize = BLOCKSIZE_TRANSFER;
    int size_padded = ((int)ceil(((double)mp->max_nibool_interfaces_ext_mesh)/((double)blocksize)))*blocksize;

    int num_blocks_x, num_blocks_y;
    get_blocks_xy(size_padded/blocksize,&num_blocks_x,&num_blocks_y);

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
      // synchronizes
      //synchronize_cuda();
      // explicitly waits until previous compute stream finishes
      // (cudaMemcpy implicitly synchronizes all other cuda operations)
      cudaStreamSynchronize(mp->compute_stream);

      // copies buffer from GPU to CPU host
      print_CUDA_error_if_any(cudaMemcpy(send_accel_buffer,mp->d_send_accel_buffer,
                              mp->size_mpi_buffer*sizeof(realw),cudaMemcpyDeviceToHost),97001);

    }
    else if(*FORWARD_OR_ADJOINT == 3) {
      prepare_boundary_accel_on_device<<<grid,threads,0,mp->compute_stream>>>(mp->d_b_accel,mp->d_b_send_accel_buffer,
                                                                              mp->num_interfaces_ext_mesh,
                                                                              mp->max_nibool_interfaces_ext_mesh,
                                                                              mp->d_nibool_interfaces_ext_mesh,
                                                                              mp->d_ibool_interfaces_ext_mesh);
      // synchronizes
      //synchronize_cuda();
      // explicitly waits until previous compute stream finishes
      // (cudaMemcpy implicitly synchronizes all other cuda operations)
      cudaStreamSynchronize(mp->compute_stream);

      // copies buffer from GPU to CPU host
      print_CUDA_error_if_any(cudaMemcpy(send_accel_buffer,mp->d_b_send_accel_buffer,
                              mp->size_mpi_buffer*sizeof(realw),cudaMemcpyDeviceToHost),97002);
    }

    // finish timing of kernel+memcpy
    // cudaEventRecord( stop, 0 );
    // cudaEventSynchronize( stop );
    // cudaEventElapsedTime( &time, start, stop );
    // cudaEventDestroy( start );
    // cudaEventDestroy( stop );
    // printf("boundary xfer d->h Time: %f ms\n",time);
  }

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

  TRACE("\ttransfer_boundary_from_device_a");

  Mesh* mp = (Mesh*)(*Mesh_pointer); // get Mesh from fortran integer wrapper

  if( mp->size_mpi_buffer > 0 ){

    int blocksize = BLOCKSIZE_TRANSFER;
    int size_padded = ((int)ceil(((double)mp->max_nibool_interfaces_ext_mesh)/((double)blocksize)))*blocksize;

    int num_blocks_x, num_blocks_y;
    get_blocks_xy(size_padded/blocksize,&num_blocks_x,&num_blocks_y);

    dim3 grid(num_blocks_x,num_blocks_y);
    dim3 threads(blocksize,1,1);

    prepare_boundary_accel_on_device<<<grid,threads,0,mp->compute_stream>>>(mp->d_accel,mp->d_send_accel_buffer,
                                                                            mp->num_interfaces_ext_mesh,
                                                                            mp->max_nibool_interfaces_ext_mesh,
                                                                            mp->d_nibool_interfaces_ext_mesh,
                                                                            mp->d_ibool_interfaces_ext_mesh);
    // waits until kernel is finished before starting async memcpy
    //synchronize_cuda();
    // waits until previous compute stream finishes
    cudaStreamSynchronize(mp->compute_stream);

    cudaMemcpyAsync(mp->h_send_accel_buffer,mp->d_send_accel_buffer,
                    mp->size_mpi_buffer*sizeof(realw),cudaMemcpyDeviceToHost,mp->copy_stream);
  }
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

  if( mp->size_mpi_buffer > 0 ){
    // copy on host memory
    memcpy(mp->h_recv_accel_buffer,buffer_recv_vector_ext_mesh,mp->size_mpi_buffer*sizeof(realw));

    // asynchronous copy to GPU using copy_stream
    cudaMemcpyAsync(mp->d_send_accel_buffer, buffer_recv_vector_ext_mesh,
                    mp->size_mpi_buffer*sizeof(realw),cudaMemcpyHostToDevice,mp->copy_stream);
  }
}


/* ----------------------------------------------------------------------------------------------- */

// Assembly

/* ----------------------------------------------------------------------------------------------- */

__global__ void assemble_boundary_accel_on_device(realw* d_accel, realw* d_send_accel_buffer,
                                                  int num_interfaces_ext_mesh,
                                                  int max_nibool_interfaces_ext_mesh,
                                                  int* d_nibool_interfaces_ext_mesh,
                                                  int* d_ibool_interfaces_ext_mesh) {

  //int bx = blockIdx.y*gridDim.x+blockIdx.x;
  //int tx = threadIdx.x;
  int id = threadIdx.x + blockIdx.x*blockDim.x + blockIdx.y*gridDim.x*blockDim.x;

  int ientry,iglob;

  for( int iinterface=0; iinterface < num_interfaces_ext_mesh; iinterface++) {
    if( id < d_nibool_interfaces_ext_mesh[iinterface] ) {

      // entry in interface array
      ientry = id + max_nibool_interfaces_ext_mesh*iinterface;
      // global index in wavefield
      iglob = d_ibool_interfaces_ext_mesh[ientry] - 1;

      // for testing atomic operations against not atomic operations (0.1ms vs. 0.04 ms)
      // d_accel[3*(iglob)] += d_send_accel_buffer[3*(ientry)];
      // d_accel[3*(iglob)+1] += d_send_accel_buffer[3*(ientry)+1];
      // d_accel[3*(iglob)+2] += d_send_accel_buffer[3*(ientry)+2];

      atomicAdd(&d_accel[3*iglob],d_send_accel_buffer[3*ientry]);
      atomicAdd(&d_accel[3*iglob + 1],d_send_accel_buffer[3*ientry + 1]);
      atomicAdd(&d_accel[3*iglob + 2],d_send_accel_buffer[3*ientry + 2]);
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
TRACE("\ttransfer_asmbl_accel_to_device");

  Mesh* mp = (Mesh*)(*Mesh_pointer); //get mesh pointer out of fortran integer container

  if( mp->size_mpi_buffer > 0 ){

    //daniel: todo - check if this copy is only needed for adjoint simulation, otherwise it is called asynchronously?
    if(*FORWARD_OR_ADJOINT == 1 ){
      // Wait until previous copy stream finishes. We assemble while other compute kernels execute.
      cudaStreamSynchronize(mp->copy_stream);
    }
    else if(*FORWARD_OR_ADJOINT == 3 ){
      // explicitly synchronizes
      // (cudaMemcpy implicitly synchronizes all other cuda operations)
      synchronize_cuda();

      print_CUDA_error_if_any(cudaMemcpy(mp->d_b_send_accel_buffer, buffer_recv_vector_ext_mesh,
                              mp->size_mpi_buffer*sizeof(realw),cudaMemcpyHostToDevice),97001);
    }

    int blocksize = BLOCKSIZE_TRANSFER;
    int size_padded = ((int)ceil(((double)mp->max_nibool_interfaces_ext_mesh)/((double)blocksize)))*blocksize;

    int num_blocks_x, num_blocks_y;
    get_blocks_xy(size_padded/blocksize,&num_blocks_x,&num_blocks_y);

    dim3 grid(num_blocks_x,num_blocks_y);
    dim3 threads(blocksize,1,1);

    //double start_time = get_time();
    // cudaEvent_t start, stop;
    // realw time;
    // cudaEventCreate(&start);
    // cudaEventCreate(&stop);
    // cudaEventRecord( start, 0 );

    if(*FORWARD_OR_ADJOINT == 1) {
      //assemble forward accel
      assemble_boundary_accel_on_device<<<grid,threads,0,mp->compute_stream>>>(mp->d_accel, mp->d_send_accel_buffer,
                                                                               mp->num_interfaces_ext_mesh,
                                                                               mp->max_nibool_interfaces_ext_mesh,
                                                                               mp->d_nibool_interfaces_ext_mesh,
                                                                               mp->d_ibool_interfaces_ext_mesh);
    }
    else if(*FORWARD_OR_ADJOINT == 3) {
      //assemble adjoint accel
      assemble_boundary_accel_on_device<<<grid,threads,0,mp->compute_stream>>>(mp->d_b_accel, mp->d_b_send_accel_buffer,
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
  }

#ifdef ENABLE_VERY_SLOW_ERROR_CHECKING
  //double end_time = get_time();
  //printf("Elapsed time: %e\n",end_time-start_time);
  exit_on_cuda_error("transfer_asmbl_accel_to_device");
#endif
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
//
//  int num_blocks_x, num_blocks_y;
//  get_blocks_xy(size_padded/blocksize,&num_blocks_x,&num_blocks_y);
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
//  exit_on_cuda_error("assemble_accel_on_device");
//#endif
//}



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
                                           realw* R_xx,realw* R_yy,realw* R_xy,
                                           realw* R_xz,realw* R_yz,
                                           realw* sigma_xx,realw* sigma_yy,realw* sigma_zz,
                                           realw* sigma_xy,realw* sigma_xz,realw* sigma_yz) {

  int i_sls,offset_sls;
  realw R_xx_val,R_yy_val;

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
                                          realw epsilondev_xx_loc,realw epsilondev_yy_loc,realw epsilondev_xy_loc,
                                          realw epsilondev_xz_loc,realw epsilondev_yz_loc
                                          ){

  int i_sls;
  int ijk_ispec;
  int offset_sls,offset_align,offset_common;
  realw mul;
  realw alphaval_loc,betaval_loc,gammaval_loc;
  realw factor_loc,Sn,Snp1;

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
                                        realw jacobianl,
                                        realw* s_dummyx_loc,
                                        realw* s_dummyy_loc,
                                        realw* s_dummyz_loc,
                                        realw* sigma_xx,
                                        realw* sigma_yy,
                                        realw* sigma_xz,
                                        realw* sigma_yz,
                                        realw* rho_s_H1,
                                        realw* rho_s_H2,
                                        realw* rho_s_H3){

  int iglob;
  realw minus_g,minus_dg;
  realw rhol;
  realw gzl; // gxl,gyl,
  realw sx_l,sy_l,sz_l;
  realw Hxxl,Hyyl,Hzzl; //,Hxyl,Hxzl,Hyzl;
  realw factor;

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

// KERNEL 2
//
// for elastic domains

/* ----------------------------------------------------------------------------------------------- */

/*

// unused
// original elastic kernel, please leave this code here for reference...

__global__ void Kernel_2_impl(int nb_blocks_to_compute,
                              int NGLOB,
                              int* d_ibool,
                              int* d_phase_ispec_inner_elastic, int num_phase_ispec_elastic,
                              int d_iphase,
                              int use_mesh_coloring_gpu,
                              realw d_deltat,
                              realw* d_displ,realw* d_veloc,realw* d_accel,
                              realw* d_xix, realw* d_xiy, realw* d_xiz,
                              realw* d_etax, realw* d_etay, realw* d_etaz,
                              realw* d_gammax, realw* d_gammay, realw* d_gammaz,
                              realw* d_hprime_xx,
                              realw* d_hprimewgll_xx,
                              realw* d_wgllwgll_xy,realw* d_wgllwgll_xz,realw* d_wgllwgll_yz,
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
                              realw* d_c11store,realw* d_c12store,realw* d_c13store,
                              realw* d_c14store,realw* d_c15store,realw* d_c16store,
                              realw* d_c22store,realw* d_c23store,realw* d_c24store,
                              realw* d_c25store,realw* d_c26store,realw* d_c33store,
                              realw* d_c34store,realw* d_c35store,realw* d_c36store,
                              realw* d_c44store,realw* d_c45store,realw* d_c46store,
                              realw* d_c55store,realw* d_c56store,realw* d_c66store,
                              int gravity,
                              realw* d_minus_g,
                              realw* d_minus_deriv_gravity,
                              realw* d_rhostore,
                              realw* wgll_cube){

  int bx = blockIdx.y*gridDim.x + blockIdx.x;
  int tx = threadIdx.x;

  const int NGLL3_ALIGN = NGLL3_PADDED;

  int K = (tx/NGLL2);
  int J = ((tx-K*NGLL2)/NGLLX);
  int I = (tx-K*NGLL2-J*NGLLX);

  int active,offset;
  int iglob = 0;
  int working_element;

  realw tempx1l,tempx2l,tempx3l,tempy1l,tempy2l,tempy3l,tempz1l,tempz2l,tempz3l;
  realw xixl,xiyl,xizl,etaxl,etayl,etazl,gammaxl,gammayl,gammazl,jacobianl;
  realw duxdxl,duxdyl,duxdzl,duydxl,duydyl,duydzl,duzdxl,duzdyl,duzdzl;
  realw duxdxl_plus_duydyl,duxdxl_plus_duzdzl,duydyl_plus_duzdzl;
  realw duxdyl_plus_duydxl,duzdxl_plus_duxdzl,duzdyl_plus_duydzl;

  realw tempx1l_att,tempx2l_att,tempx3l_att,tempy1l_att,tempy2l_att,tempy3l_att,tempz1l_att,tempz2l_att,tempz3l_att;
  realw duxdxl_att,duxdyl_att,duxdzl_att,duydxl_att,duydyl_att,duydzl_att,duzdxl_att,duzdyl_att,duzdzl_att;
  realw duxdyl_plus_duydxl_att,duzdxl_plus_duxdzl_att,duzdyl_plus_duydzl_att;

  realw fac1,fac2,fac3,lambdal,mul,lambdalplus2mul,kappal;
  realw sigma_xx,sigma_yy,sigma_zz,sigma_xy,sigma_xz,sigma_yz;
  realw epsilondev_xx_loc,epsilondev_yy_loc,epsilondev_xy_loc,epsilondev_xz_loc,epsilondev_yz_loc;

  realw c11,c12,c13,c14,c15,c16,c22,c23,c24,c25,c26,c33,c34,c35,c36,c44,c45,c46,c55,c56,c66;
  realw sum_terms1,sum_terms2,sum_terms3;

  // gravity variables
  realw sigma_yx,sigma_zx,sigma_zy;
  realw rho_s_H1,rho_s_H2,rho_s_H3;

#ifndef MANUALLY_UNROLLED_LOOPS
  int l;
  realw hp1,hp2,hp3;
#endif

  __shared__ realw s_dummyx_loc[NGLL3];
  __shared__ realw s_dummyy_loc[NGLL3];
  __shared__ realw s_dummyz_loc[NGLL3];

  __shared__ realw s_dummyx_loc_att[NGLL3];
  __shared__ realw s_dummyy_loc_att[NGLL3];
  __shared__ realw s_dummyz_loc_att[NGLL3];

  __shared__ realw s_tempx1[NGLL3];
  __shared__ realw s_tempx2[NGLL3];
  __shared__ realw s_tempx3[NGLL3];
  __shared__ realw s_tempy1[NGLL3];
  __shared__ realw s_tempy2[NGLL3];
  __shared__ realw s_tempy3[NGLL3];
  __shared__ realw s_tempz1[NGLL3];
  __shared__ realw s_tempz2[NGLL3];
  __shared__ realw s_tempz3[NGLL3];

  __shared__ realw sh_hprime_xx[NGLL2];

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

  // JC JC here we will need to add GPU support for the new C-PML routines
  if(ATTENUATION){
    // use first order Taylor expansion of displacement for local storage of stresses
    // at this current time step, to fix attenuation in a consistent way
#ifdef USE_TEXTURES_FIELDS
    s_dummyx_loc_att[tx] = s_dummyx_loc[tx] + d_deltat * tex1Dfetch(d_veloc_tex, iglob);
    s_dummyy_loc_att[tx] = s_dummyy_loc[tx] + d_deltat * tex1Dfetch(d_veloc_tex, iglob + NGLOB);
    s_dummyz_loc_att[tx] = s_dummyz_loc[tx] + d_deltat * tex1Dfetch(d_veloc_tex, iglob + 2*NGLOB);
#else
    s_dummyx_loc_att[tx] = s_dummyx_loc[tx] + d_deltat * d_veloc[iglob*3];
    s_dummyy_loc_att[tx] = s_dummyy_loc[tx] + d_deltat * d_veloc[iglob*3 + 1];
    s_dummyz_loc_att[tx] = s_dummyz_loc[tx] + d_deltat * d_veloc[iglob*3 + 2];
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

        //assumes that hprime_xx = hprime_yy = hprime_zz
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

    // JC JC here we will need to add GPU support for the new C-PML routines
    if( ATTENUATION){
      // temporary variables used for fixing attenuation in a consistent way

      tempx1l_att = 0.f;
      tempx2l_att = 0.f;
      tempx3l_att = 0.f;

      tempy1l_att = 0.f;
      tempy2l_att = 0.f;
      tempy3l_att = 0.f;

      tempz1l_att = 0.f;
      tempz2l_att = 0.f;
      tempz3l_att = 0.f;

      for (l=0;l<NGLLX;l++) {
  hp1 = sh_hprime_xx[l*NGLLX+I];
  offset = K*NGLL2+J*NGLLX+l;
  tempx1l_att += s_dummyx_loc_att[offset]*hp1;
  tempy1l_att += s_dummyy_loc_att[offset]*hp1;
  tempz1l_att += s_dummyz_loc_att[offset]*hp1;

  hp2 = sh_hprime_xx[l*NGLLX+J];
  offset = K*NGLL2+l*NGLLX+I;
  tempx2l_att += s_dummyx_loc_att[offset]*hp2;
  tempy2l_att += s_dummyy_loc_att[offset]*hp2;
  tempz2l_att += s_dummyz_loc_att[offset]*hp2;

  hp3 = sh_hprime_xx[l*NGLLX+K];
  offset = l*NGLL2+J*NGLLX+I;
  tempx3l_att += s_dummyx_loc_att[offset]*hp3;
  tempy3l_att += s_dummyy_loc_att[offset]*hp3;
  tempz3l_att += s_dummyz_loc_att[offset]*hp3;

      }
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

    // JC JC here we will need to add GPU support for the new C-PML routines
    if( ATTENUATION){
      // temporary variables used for fixing attenuation in a consistent way

      tempx1l_att = s_dummyx_loc_att[K*NGLL2+J*NGLLX]*d_hprime_xx[I]
  + s_dummyx_loc_att[K*NGLL2+J*NGLLX+1]*d_hprime_xx[NGLLX+I]
  + s_dummyx_loc_att[K*NGLL2+J*NGLLX+2]*d_hprime_xx[2*NGLLX+I]
  + s_dummyx_loc_att[K*NGLL2+J*NGLLX+3]*d_hprime_xx[3*NGLLX+I]
  + s_dummyx_loc_att[K*NGLL2+J*NGLLX+4]*d_hprime_xx[4*NGLLX+I];

      tempy1l_att = s_dummyy_loc_att[K*NGLL2+J*NGLLX]*d_hprime_xx[I]
  + s_dummyy_loc_att[K*NGLL2+J*NGLLX+1]*d_hprime_xx[NGLLX+I]
  + s_dummyy_loc_att[K*NGLL2+J*NGLLX+2]*d_hprime_xx[2*NGLLX+I]
  + s_dummyy_loc_att[K*NGLL2+J*NGLLX+3]*d_hprime_xx[3*NGLLX+I]
  + s_dummyy_loc_att[K*NGLL2+J*NGLLX+4]*d_hprime_xx[4*NGLLX+I];

      tempz1l_att = s_dummyz_loc_att[K*NGLL2+J*NGLLX]*d_hprime_xx[I]
  + s_dummyz_loc_att[K*NGLL2+J*NGLLX+1]*d_hprime_xx[NGLLX+I]
  + s_dummyz_loc_att[K*NGLL2+J*NGLLX+2]*d_hprime_xx[2*NGLLX+I]
  + s_dummyz_loc_att[K*NGLL2+J*NGLLX+3]*d_hprime_xx[3*NGLLX+I]
  + s_dummyz_loc_att[K*NGLL2+J*NGLLX+4]*d_hprime_xx[4*NGLLX+I];

      tempx2l_att = s_dummyx_loc_att[K*NGLL2+I]*d_hprime_xx[J]
  + s_dummyx_loc_att[K*NGLL2+NGLLX+I]*d_hprime_xx[NGLLX+J]
  + s_dummyx_loc_att[K*NGLL2+2*NGLLX+I]*d_hprime_xx[2*NGLLX+J]
  + s_dummyx_loc_att[K*NGLL2+3*NGLLX+I]*d_hprime_xx[3*NGLLX+J]
  + s_dummyx_loc_att[K*NGLL2+4*NGLLX+I]*d_hprime_xx[4*NGLLX+J];

      tempy2l_att = s_dummyy_loc_att[K*NGLL2+I]*d_hprime_xx[J]
  + s_dummyy_loc_att[K*NGLL2+NGLLX+I]*d_hprime_xx[NGLLX+J]
  + s_dummyy_loc_att[K*NGLL2+2*NGLLX+I]*d_hprime_xx[2*NGLLX+J]
  + s_dummyy_loc_att[K*NGLL2+3*NGLLX+I]*d_hprime_xx[3*NGLLX+J]
  + s_dummyy_loc_att[K*NGLL2+4*NGLLX+I]*d_hprime_xx[4*NGLLX+J];

      tempz2l_att = s_dummyz_loc_att[K*NGLL2+I]*d_hprime_xx[J]
  + s_dummyz_loc_att[K*NGLL2+NGLLX+I]*d_hprime_xx[NGLLX+J]
  + s_dummyz_loc_att[K*NGLL2+2*NGLLX+I]*d_hprime_xx[2*NGLLX+J]
  + s_dummyz_loc_att[K*NGLL2+3*NGLLX+I]*d_hprime_xx[3*NGLLX+J]
  + s_dummyz_loc_att[K*NGLL2+4*NGLLX+I]*d_hprime_xx[4*NGLLX+J];

      tempx3l_att = s_dummyx_loc_att[J*NGLLX+I]*d_hprime_xx[K]
  + s_dummyx_loc_att[NGLL2+J*NGLLX+I]*d_hprime_xx[NGLLX+K]
  + s_dummyx_loc_att[2*NGLL2+J*NGLLX+I]*d_hprime_xx[2*NGLLX+K]
  + s_dummyx_loc_att[3*NGLL2+J*NGLLX+I]*d_hprime_xx[3*NGLLX+K]
  + s_dummyx_loc_att[4*NGLL2+J*NGLLX+I]*d_hprime_xx[4*NGLLX+K];

      tempy3l_att = s_dummyy_loc_att[J*NGLLX+I]*d_hprime_xx[K]
  + s_dummyy_loc_att[NGLL2+J*NGLLX+I]*d_hprime_xx[NGLLX+K]
  + s_dummyy_loc_att[2*NGLL2+J*NGLLX+I]*d_hprime_xx[2*NGLLX+K]
  + s_dummyy_loc_att[3*NGLL2+J*NGLLX+I]*d_hprime_xx[3*NGLLX+K]
  + s_dummyy_loc_att[4*NGLL2+J*NGLLX+I]*d_hprime_xx[4*NGLLX+K];

      tempz3l_att = s_dummyz_loc_att[J*NGLLX+I]*d_hprime_xx[K]
  + s_dummyz_loc_att[NGLL2+J*NGLLX+I]*d_hprime_xx[NGLLX+K]
  + s_dummyz_loc_att[2*NGLL2+J*NGLLX+I]*d_hprime_xx[2*NGLLX+K]
  + s_dummyz_loc_att[3*NGLL2+J*NGLLX+I]*d_hprime_xx[3*NGLLX+K]
  + s_dummyz_loc_att[4*NGLL2+J*NGLLX+I]*d_hprime_xx[4*NGLLX+K];
    }

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

    // JC JC here we will need to add GPU support for the new C-PML routines

    // precompute some sums to save CPU time
    duxdxl_plus_duydyl = duxdxl + duydyl;
    duxdxl_plus_duzdzl = duxdxl + duzdzl;
    duydyl_plus_duzdzl = duydyl + duzdzl;
    duxdyl_plus_duydxl = duxdyl + duydxl;
    duzdxl_plus_duxdzl = duzdxl + duxdzl;
    duzdyl_plus_duydzl = duzdyl + duydzl;

    // JC JC here we will need to add GPU support for the new C-PML routines
    if( ATTENUATION){
      // temporary variables used for fixing attenuation in a consistent way

      duxdxl_att = xixl*tempx1l_att + etaxl*tempx2l_att + gammaxl*tempx3l_att;
      duxdyl_att = xiyl*tempx1l_att + etayl*tempx2l_att + gammayl*tempx3l_att;
      duxdzl_att = xizl*tempx1l_att + etazl*tempx2l_att + gammazl*tempx3l_att;

      duydxl_att = xixl*tempy1l_att + etaxl*tempy2l_att + gammaxl*tempy3l_att;
      duydyl_att = xiyl*tempy1l_att + etayl*tempy2l_att + gammayl*tempy3l_att;
      duydzl_att = xizl*tempy1l_att + etazl*tempy2l_att + gammazl*tempy3l_att;

      duzdxl_att = xixl*tempz1l_att + etaxl*tempz2l_att + gammaxl*tempz3l_att;
      duzdyl_att = xiyl*tempz1l_att + etayl*tempz2l_att + gammayl*tempz3l_att;
      duzdzl_att = xizl*tempz1l_att + etazl*tempz2l_att + gammazl*tempz3l_att;

      // precompute some sums to save CPU time
      duxdyl_plus_duydxl_att = duxdyl_att + duydxl_att;
      duzdxl_plus_duxdzl_att = duzdxl_att + duxdzl_att;
      duzdyl_plus_duydzl_att = duzdyl_att + duydzl_att;

      // computes deviatoric strain attenuation and/or for kernel calculations
      if(COMPUTE_AND_STORE_STRAIN) {
  realw templ = 0.33333333333333333333f * (duxdxl_att + duydyl_att + duzdzl_att); // 1./3. = 0.33333

  // local storage: stresses at this current time step
  epsilondev_xx_loc = duxdxl_att - templ;
  epsilondev_yy_loc = duydyl_att - templ;
  epsilondev_xy_loc = 0.5f * duxdyl_plus_duydxl_att;
  epsilondev_xz_loc = 0.5f * duzdxl_plus_duxdzl_att;
  epsilondev_yz_loc = 0.5f * duzdyl_plus_duydzl_att;

  if(SIMULATION_TYPE == 3) {
    epsilon_trace_over_3[tx + working_element*NGLL3] = templ;
  }

  // JC JC here we will need to add GPU support for the new C-PML routines
      }
    }else{
      // computes deviatoric strain attenuation and/or for kernel calculations
      if(COMPUTE_AND_STORE_STRAIN) {
  realw templ = 0.33333333333333333333f * (duxdxl + duydyl + duzdzl); // 1./3. = 0.33333

  //  epsilondev_xx[offset] = duxdxl - templ;
  //  epsilondev_yy[offset] = duydyl - templ;
  //  epsilondev_xy[offset] = 0.5f * duxdyl_plus_duydxl;
  //  epsilondev_xz[offset] = 0.5f * duzdxl_plus_duxdzl;
  //  epsilondev_yz[offset] = 0.5f * duzdyl_plus_duydzl;

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

  // JC JC here we will need to add GPU support for the new C-PML routines

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

      // assumes hprimewgll_xx == hprimewgll_yy == hprimewgll_zz
      fac2 = d_hprimewgll_xx[J*NGLLX+l];
      offset = K*NGLL2+l*NGLLX+I;
      tempx2l += s_tempx2[offset]*fac2;
      tempy2l += s_tempy2[offset]*fac2;
      tempz2l += s_tempz2[offset]*fac2;

      fac3 = d_hprimewgll_xx[K*NGLLX+l];
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

    tempx2l = s_tempx2[K*NGLL2+I]*d_hprimewgll_xx[J*NGLLX]
            + s_tempx2[K*NGLL2+NGLLX+I]*d_hprimewgll_xx[J*NGLLX+1]
            + s_tempx2[K*NGLL2+2*NGLLX+I]*d_hprimewgll_xx[J*NGLLX+2]
            + s_tempx2[K*NGLL2+3*NGLLX+I]*d_hprimewgll_xx[J*NGLLX+3]
            + s_tempx2[K*NGLL2+4*NGLLX+I]*d_hprimewgll_xx[J*NGLLX+4];

    tempy2l = s_tempy2[K*NGLL2+I]*d_hprimewgll_xx[J*NGLLX]
            + s_tempy2[K*NGLL2+NGLLX+I]*d_hprimewgll_xx[J*NGLLX+1]
            + s_tempy2[K*NGLL2+2*NGLLX+I]*d_hprimewgll_xx[J*NGLLX+2]
            + s_tempy2[K*NGLL2+3*NGLLX+I]*d_hprimewgll_xx[J*NGLLX+3]
            + s_tempy2[K*NGLL2+4*NGLLX+I]*d_hprimewgll_xx[J*NGLLX+4];

    tempz2l = s_tempz2[K*NGLL2+I]*d_hprimewgll_xx[J*NGLLX]
            + s_tempz2[K*NGLL2+NGLLX+I]*d_hprimewgll_xx[J*NGLLX+1]
            + s_tempz2[K*NGLL2+2*NGLLX+I]*d_hprimewgll_xx[J*NGLLX+2]
            + s_tempz2[K*NGLL2+3*NGLLX+I]*d_hprimewgll_xx[J*NGLLX+3]
            + s_tempz2[K*NGLL2+4*NGLLX+I]*d_hprimewgll_xx[J*NGLLX+4];

    tempx3l = s_tempx3[J*NGLLX+I]*d_hprimewgll_xx[K*NGLLX]
            + s_tempx3[NGLL2+J*NGLLX+I]*d_hprimewgll_xx[K*NGLLX+1]
            + s_tempx3[2*NGLL2+J*NGLLX+I]*d_hprimewgll_xx[K*NGLLX+2]
            + s_tempx3[3*NGLL2+J*NGLLX+I]*d_hprimewgll_xx[K*NGLLX+3]
            + s_tempx3[4*NGLL2+J*NGLLX+I]*d_hprimewgll_xx[K*NGLLX+4];

    tempy3l = s_tempy3[J*NGLLX+I]*d_hprimewgll_xx[K*NGLLX]
            + s_tempy3[NGLL2+J*NGLLX+I]*d_hprimewgll_xx[K*NGLLX+1]
            + s_tempy3[2*NGLL2+J*NGLLX+I]*d_hprimewgll_xx[K*NGLLX+2]
            + s_tempy3[3*NGLL2+J*NGLLX+I]*d_hprimewgll_xx[K*NGLLX+3]
            + s_tempy3[4*NGLL2+J*NGLLX+I]*d_hprimewgll_xx[K*NGLLX+4];

    tempz3l = s_tempz3[J*NGLLX+I]*d_hprimewgll_xx[K*NGLLX]
            + s_tempz3[NGLL2+J*NGLLX+I]*d_hprimewgll_xx[K*NGLLX+1]
            + s_tempz3[2*NGLL2+J*NGLLX+I]*d_hprimewgll_xx[K*NGLLX+2]
            + s_tempz3[3*NGLL2+J*NGLLX+I]*d_hprimewgll_xx[K*NGLLX+3]
            + s_tempz3[4*NGLL2+J*NGLLX+I]*d_hprimewgll_xx[K*NGLLX+4];

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

    // JC JC here we will need to add GPU support for the new C-PML routines

#else // MESH_COLORING

    //mesh coloring
    if( use_mesh_coloring_gpu ){

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

  // JC JC here we will need to add GPU support for the new C-PML routines

} // kernel_2_impl()

*/

/* ----------------------------------------------------------------------------------------------- */

// note: kernel_2 is split into two kernels:
//       - a kernel without attenuation Kernel_2_noatt_impl() and
//       - a kernel including attenuation Kernel_2_att_impl()
//       this separation should help with performance


// kernel without attenuation
//
// we use templates to distinguish between calls with forward or adjoint texture fields

template<int FORWARD_OR_ADJOINT> __global__ void Kernel_2_noatt_impl(int nb_blocks_to_compute,
                                                                      int NGLOB,
                                                                      int* d_ibool,
                                                                      int* d_phase_ispec_inner_elastic, int num_phase_ispec_elastic,
                                                                      int d_iphase,
                                                                      int use_mesh_coloring_gpu,
                                                                      realw* d_displ,realw* d_veloc,realw* d_accel,
                                                                      realw* d_xix, realw* d_xiy, realw* d_xiz,
                                                                      realw* d_etax, realw* d_etay, realw* d_etaz,
                                                                      realw* d_gammax, realw* d_gammay, realw* d_gammaz,
                                                                      realw* d_hprime_xx,
                                                                      realw* d_hprimewgll_xx,
                                                                      realw* d_wgllwgll_xy,realw* d_wgllwgll_xz,realw* d_wgllwgll_yz,
                                                                      realw* d_kappav, realw* d_muv,
                                                                      int COMPUTE_AND_STORE_STRAIN,
                                                                      realw* epsilondev_xx,realw* epsilondev_yy,realw* epsilondev_xy,
                                                                      realw* epsilondev_xz,realw* epsilondev_yz,
                                                                      realw* epsilon_trace_over_3,
                                                                      int SIMULATION_TYPE,
                                                                      int NSPEC,
                                                                      realw* one_minus_sum_beta,realw* factor_common,
                                                                      realw* R_xx, realw* R_yy, realw* R_xy, realw* R_xz, realw* R_yz,
                                                                      realw* alphaval,realw* betaval,realw* gammaval,
                                                                      int ANISOTROPY,
                                                                      realw* d_c11store,realw* d_c12store,realw* d_c13store,
                                                                      realw* d_c14store,realw* d_c15store,realw* d_c16store,
                                                                      realw* d_c22store,realw* d_c23store,realw* d_c24store,
                                                                      realw* d_c25store,realw* d_c26store,realw* d_c33store,
                                                                      realw* d_c34store,realw* d_c35store,realw* d_c36store,
                                                                      realw* d_c44store,realw* d_c45store,realw* d_c46store,
                                                                      realw* d_c55store,realw* d_c56store,realw* d_c66store,
                                                                      int gravity,
                                                                      realw* d_minus_g,
                                                                      realw* d_minus_deriv_gravity,
                                                                      realw* d_rhostore,
                                                                      realw* wgll_cube ){

// elastic compute kernel without attenuation
// holds for: ATTENUATION = .false.
//            COMPUTE_AND_STORE_STRAIN = .true. or .false. (true for kernel simulations)

  int bx = blockIdx.y*gridDim.x+blockIdx.x;
  int tx = threadIdx.x;

  const int NGLL3_ALIGN = NGLL3_PADDED;

  int K = (tx/NGLL2);
  int J = ((tx-K*NGLL2)/NGLLX);
  int I = (tx-K*NGLL2-J*NGLLX);

  int active,offset;
  int iglob = 0;
  int working_element;

  realw tempx1l,tempx2l,tempx3l,tempy1l,tempy2l,tempy3l,tempz1l,tempz2l,tempz3l;
  realw xixl,xiyl,xizl,etaxl,etayl,etazl,gammaxl,gammayl,gammazl,jacobianl;
  realw duxdxl,duxdyl,duxdzl,duydxl,duydyl,duydzl,duzdxl,duzdyl,duzdzl;
  realw duxdxl_plus_duydyl,duxdxl_plus_duzdzl,duydyl_plus_duzdzl;
  realw duxdyl_plus_duydxl,duzdxl_plus_duxdzl,duzdyl_plus_duydzl;

  realw fac1,fac2,fac3,lambdal,mul,lambdalplus2mul,kappal;
  realw sigma_xx,sigma_yy,sigma_zz,sigma_xy,sigma_xz,sigma_yz;
  realw epsilondev_xx_loc,epsilondev_yy_loc,epsilondev_xy_loc,epsilondev_xz_loc,epsilondev_yz_loc;

  realw c11,c12,c13,c14,c15,c16,c22,c23,c24,c25,c26,c33,c34,c35,c36,c44,c45,c46,c55,c56,c66;
  realw sum_terms1,sum_terms2,sum_terms3;

  // gravity variables
  realw sigma_yx,sigma_zx,sigma_zy;
  realw rho_s_H1,rho_s_H2,rho_s_H3;

#ifndef MANUALLY_UNROLLED_LOOPS
  int l;
  realw hp1,hp2,hp3;
#endif

  __shared__ realw s_dummyx_loc[NGLL3];
  __shared__ realw s_dummyy_loc[NGLL3];
  __shared__ realw s_dummyz_loc[NGLL3];

  __shared__ realw s_tempx1[NGLL3];
  __shared__ realw s_tempx2[NGLL3];
  __shared__ realw s_tempx3[NGLL3];

  __shared__ realw s_tempy1[NGLL3];
  __shared__ realw s_tempy2[NGLL3];
  __shared__ realw s_tempy3[NGLL3];

  __shared__ realw s_tempz1[NGLL3];
  __shared__ realw s_tempz2[NGLL3];
  __shared__ realw s_tempz3[NGLL3];

  __shared__ realw sh_hprime_xx[NGLL2];

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

    iglob = d_ibool[working_element*NGLL3 + tx]-1;
    // debug
    //if( iglob < 0 || iglob >= NGLOB ){ printf("wrong iglob %d\n",iglob);  }

#ifdef USE_TEXTURES_FIELDS
    s_dummyx_loc[tx] = texfetch_displ<FORWARD_OR_ADJOINT>(iglob*3);
    s_dummyy_loc[tx] = texfetch_displ<FORWARD_OR_ADJOINT>(iglob*3 + 1);
    s_dummyz_loc[tx] = texfetch_displ<FORWARD_OR_ADJOINT>(iglob*3 + 2);
#else
    // changing iglob indexing to match fortran row changes fast style
    s_dummyx_loc[tx] = d_displ[iglob*3];
    s_dummyy_loc[tx] = d_displ[iglob*3 + 1];
    s_dummyz_loc[tx] = d_displ[iglob*3 + 2];
#endif
  }

  // JC JC here we will need to add GPU support for the new C-PML routines

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

      //assumes that hprime_xx = hprime_yy = hprime_zz
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

    // JC JC here we will need to add GPU support for the new C-PML routines

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

    // JC JC here we will need to add GPU support for the new C-PML routines

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

    // JC JC here we will need to add GPU support for the new C-PML routines

    // precompute some sums to save CPU time
    duxdxl_plus_duydyl = duxdxl + duydyl;
    duxdxl_plus_duzdzl = duxdxl + duzdzl;
    duydyl_plus_duzdzl = duydyl + duzdzl;
    duxdyl_plus_duydxl = duxdyl + duydxl;
    duzdxl_plus_duxdzl = duzdxl + duxdzl;
    duzdyl_plus_duydzl = duzdyl + duydzl;

    // JC JC here we will need to add GPU support for the new C-PML routines

    // computes deviatoric strain for kernel calculations
    if(COMPUTE_AND_STORE_STRAIN) {
      realw templ = 0.33333333333333333333f * (duxdxl + duydyl + duzdzl); // 1./3. = 0.33333
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

  // JC JC here we will need to add GPU support for the new C-PML routines

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

      // assumes hprimewgll_xx == hprimewgll_yy == hprimewgll_zz
      fac2 = d_hprimewgll_xx[J*NGLLX+l];
      offset = K*NGLL2+l*NGLLX+I;
      tempx2l += s_tempx2[offset]*fac2;
      tempy2l += s_tempy2[offset]*fac2;
      tempz2l += s_tempz2[offset]*fac2;

      fac3 = d_hprimewgll_xx[K*NGLLX+l];
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

    tempx2l = s_tempx2[K*NGLL2+I]*d_hprimewgll_xx[J*NGLLX]
            + s_tempx2[K*NGLL2+NGLLX+I]*d_hprimewgll_xx[J*NGLLX+1]
            + s_tempx2[K*NGLL2+2*NGLLX+I]*d_hprimewgll_xx[J*NGLLX+2]
            + s_tempx2[K*NGLL2+3*NGLLX+I]*d_hprimewgll_xx[J*NGLLX+3]
            + s_tempx2[K*NGLL2+4*NGLLX+I]*d_hprimewgll_xx[J*NGLLX+4];

    tempy2l = s_tempy2[K*NGLL2+I]*d_hprimewgll_xx[J*NGLLX]
            + s_tempy2[K*NGLL2+NGLLX+I]*d_hprimewgll_xx[J*NGLLX+1]
            + s_tempy2[K*NGLL2+2*NGLLX+I]*d_hprimewgll_xx[J*NGLLX+2]
            + s_tempy2[K*NGLL2+3*NGLLX+I]*d_hprimewgll_xx[J*NGLLX+3]
            + s_tempy2[K*NGLL2+4*NGLLX+I]*d_hprimewgll_xx[J*NGLLX+4];

    tempz2l = s_tempz2[K*NGLL2+I]*d_hprimewgll_xx[J*NGLLX]
            + s_tempz2[K*NGLL2+NGLLX+I]*d_hprimewgll_xx[J*NGLLX+1]
            + s_tempz2[K*NGLL2+2*NGLLX+I]*d_hprimewgll_xx[J*NGLLX+2]
            + s_tempz2[K*NGLL2+3*NGLLX+I]*d_hprimewgll_xx[J*NGLLX+3]
            + s_tempz2[K*NGLL2+4*NGLLX+I]*d_hprimewgll_xx[J*NGLLX+4];

    tempx3l = s_tempx3[J*NGLLX+I]*d_hprimewgll_xx[K*NGLLX]
            + s_tempx3[NGLL2+J*NGLLX+I]*d_hprimewgll_xx[K*NGLLX+1]
            + s_tempx3[2*NGLL2+J*NGLLX+I]*d_hprimewgll_xx[K*NGLLX+2]
            + s_tempx3[3*NGLL2+J*NGLLX+I]*d_hprimewgll_xx[K*NGLLX+3]
            + s_tempx3[4*NGLL2+J*NGLLX+I]*d_hprimewgll_xx[K*NGLLX+4];

    tempy3l = s_tempy3[J*NGLLX+I]*d_hprimewgll_xx[K*NGLLX]
            + s_tempy3[NGLL2+J*NGLLX+I]*d_hprimewgll_xx[K*NGLLX+1]
            + s_tempy3[2*NGLL2+J*NGLLX+I]*d_hprimewgll_xx[K*NGLLX+2]
            + s_tempy3[3*NGLL2+J*NGLLX+I]*d_hprimewgll_xx[K*NGLLX+3]
            + s_tempy3[4*NGLL2+J*NGLLX+I]*d_hprimewgll_xx[K*NGLLX+4];

    tempz3l = s_tempz3[J*NGLLX+I]*d_hprimewgll_xx[K*NGLLX]
            + s_tempz3[NGLL2+J*NGLLX+I]*d_hprimewgll_xx[K*NGLLX+1]
            + s_tempz3[2*NGLL2+J*NGLLX+I]*d_hprimewgll_xx[K*NGLLX+2]
            + s_tempz3[3*NGLL2+J*NGLLX+I]*d_hprimewgll_xx[K*NGLLX+3]
            + s_tempz3[4*NGLL2+J*NGLLX+I]*d_hprimewgll_xx[K*NGLLX+4];

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
    d_accel[iglob*3]     = texfetch_accel<FORWARD_OR_ADJOINT>(iglob*3) + sum_terms1;
    d_accel[iglob*3 + 1] = texfetch_accel<FORWARD_OR_ADJOINT>(iglob*3 + 1) + sum_terms2;
    d_accel[iglob*3 + 2] = texfetch_accel<FORWARD_OR_ADJOINT>(iglob*3 + 2) + sum_terms3;
#else
    d_accel[iglob*3]     += sum_terms1;
    d_accel[iglob*3 + 1] += sum_terms2;
    d_accel[iglob*3 + 2] += sum_terms3;
#endif // USE_TEXTURES_FIELDS

    // JC JC here we will need to add GPU support for the new C-PML routines

#else // MESH_COLORING

    //mesh coloring
    if( use_mesh_coloring_gpu ){

      // no atomic operation needed, colors don't share global points between elements
#ifdef USE_TEXTURES_FIELDS
      d_accel[iglob*3]     = texfetch_accel<FORWARD_OR_ADJOINT>(iglob*3) + sum_terms1;
      d_accel[iglob*3 + 1] = texfetch_accel<FORWARD_OR_ADJOINT>(iglob*3 + 1) + sum_terms2;
      d_accel[iglob*3 + 2] = texfetch_accel<FORWARD_OR_ADJOINT>(iglob*3 + 2) + sum_terms3;
#else
      d_accel[iglob*3]     += sum_terms1;
      d_accel[iglob*3 + 1] += sum_terms2;
      d_accel[iglob*3 + 2] += sum_terms3;
#endif // USE_TEXTURES_FIELDS

    }else {

      // for testing purposes only: w/out atomic updates
      //d_accel[iglob*3] -= (0.00000001f*tempx1l + 0.00000001f*tempx2l + 0.00000001f*tempx3l);
      //d_accel[iglob*3 + 1] -= (0.00000001f*tempy1l + 0.00000001f*tempy2l + 0.00000001f*tempy3l);
      //d_accel[iglob*3 + 2] -= (0.00000001f*tempz1l + 0.00000001f*tempz2l + 0.00000001f*tempz3l);
      // w/out atomic update
      //d_accel[iglob*3]     += sum_terms1;
      //d_accel[iglob*3 + 1] += sum_terms2;
      //d_accel[iglob*3 + 2] += sum_terms3;

      atomicAdd(&d_accel[iglob*3], sum_terms1);
      atomicAdd(&d_accel[iglob*3+1], sum_terms2);
      atomicAdd(&d_accel[iglob*3+2], sum_terms3);

    } // if(use_mesh_coloring_gpu)

#endif // MESH_COLORING

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

  // JC JC here we will need to add GPU support for the new C-PML routines

} // kernel_2_noatt_impl()


/* ----------------------------------------------------------------------------------------------- */

// kernel with attenuation
//
// we use templates to distinguish between calls with forward or adjoint texture fields

template<int FORWARD_OR_ADJOINT> __global__ void Kernel_2_att_impl(int nb_blocks_to_compute,
                                                                  int NGLOB,
                                                                  int* d_ibool,
                                                                  int* d_phase_ispec_inner_elastic, int num_phase_ispec_elastic,
                                                                  int d_iphase,
                                                                  int use_mesh_coloring_gpu,
                                                                  realw d_deltat,
                                                                  realw* d_displ,realw* d_veloc,realw* d_accel,
                                                                  realw* d_xix, realw* d_xiy, realw* d_xiz,
                                                                  realw* d_etax, realw* d_etay, realw* d_etaz,
                                                                  realw* d_gammax, realw* d_gammay, realw* d_gammaz,
                                                                  realw* d_hprime_xx,
                                                                  realw* d_hprimewgll_xx,
                                                                  realw* d_wgllwgll_xy,realw* d_wgllwgll_xz,realw* d_wgllwgll_yz,
                                                                  realw* d_kappav, realw* d_muv,
                                                                  realw* epsilondev_xx,realw* epsilondev_yy,realw* epsilondev_xy,
                                                                  realw* epsilondev_xz,realw* epsilondev_yz,
                                                                  realw* epsilon_trace_over_3,
                                                                  int SIMULATION_TYPE,
                                                                  int NSPEC,
                                                                  realw* one_minus_sum_beta,realw* factor_common,
                                                                  realw* R_xx, realw* R_yy, realw* R_xy, realw* R_xz, realw* R_yz,
                                                                  realw* alphaval,realw* betaval,realw* gammaval,
                                                                  int ANISOTROPY,
                                                                  realw* d_c11store,realw* d_c12store,realw* d_c13store,
                                                                  realw* d_c14store,realw* d_c15store,realw* d_c16store,
                                                                  realw* d_c22store,realw* d_c23store,realw* d_c24store,
                                                                  realw* d_c25store,realw* d_c26store,realw* d_c33store,
                                                                  realw* d_c34store,realw* d_c35store,realw* d_c36store,
                                                                  realw* d_c44store,realw* d_c45store,realw* d_c46store,
                                                                  realw* d_c55store,realw* d_c56store,realw* d_c66store,
                                                                  int gravity,
                                                                  realw* d_minus_g,
                                                                  realw* d_minus_deriv_gravity,
                                                                  realw* d_rhostore,
                                                                  realw* wgll_cube){


// elastic compute kernel with attenuation
// holds for: ATTENUATION = .true.
//            COMPUTE_AND_STORE_STRAIN = .true. (always true for attenuation)

  int bx = blockIdx.y*gridDim.x+blockIdx.x;
  int tx = threadIdx.x;

  const int NGLL3_ALIGN = NGLL3_PADDED;

  int K = (tx/NGLL2);
  int J = ((tx-K*NGLL2)/NGLLX);
  int I = (tx-K*NGLL2-J*NGLLX);

  int active,offset;
  int iglob = 0;
  int working_element;

  realw tempx1l,tempx2l,tempx3l,tempy1l,tempy2l,tempy3l,tempz1l,tempz2l,tempz3l;
  realw xixl,xiyl,xizl,etaxl,etayl,etazl,gammaxl,gammayl,gammazl,jacobianl;
  realw duxdxl,duxdyl,duxdzl,duydxl,duydyl,duydzl,duzdxl,duzdyl,duzdzl;
  realw duxdxl_plus_duydyl,duxdxl_plus_duzdzl,duydyl_plus_duzdzl;
  realw duxdyl_plus_duydxl,duzdxl_plus_duxdzl,duzdyl_plus_duydzl;

  realw tempx1l_att,tempx2l_att,tempx3l_att,tempy1l_att,tempy2l_att,tempy3l_att,tempz1l_att,tempz2l_att,tempz3l_att;
  realw duxdxl_att,duxdyl_att,duxdzl_att,duydxl_att,duydyl_att,duydzl_att,duzdxl_att,duzdyl_att,duzdzl_att;
  realw duxdyl_plus_duydxl_att,duzdxl_plus_duxdzl_att,duzdyl_plus_duydzl_att;

  realw fac1,fac2,fac3,lambdal,mul,lambdalplus2mul,kappal;
  realw sigma_xx,sigma_yy,sigma_zz,sigma_xy,sigma_xz,sigma_yz;
  realw epsilondev_xx_loc,epsilondev_yy_loc,epsilondev_xy_loc,epsilondev_xz_loc,epsilondev_yz_loc;

  realw c11,c12,c13,c14,c15,c16,c22,c23,c24,c25,c26,c33,c34,c35,c36,c44,c45,c46,c55,c56,c66;
  realw sum_terms1,sum_terms2,sum_terms3;

  // gravity variables
  realw sigma_yx,sigma_zx,sigma_zy;
  realw rho_s_H1,rho_s_H2,rho_s_H3;

#ifndef MANUALLY_UNROLLED_LOOPS
  int l;
  realw hp1,hp2,hp3;
#endif

  __shared__ realw s_dummyx_loc[NGLL3];
  __shared__ realw s_dummyy_loc[NGLL3];
  __shared__ realw s_dummyz_loc[NGLL3];

  __shared__ realw s_dummyx_loc_att[NGLL3];
  __shared__ realw s_dummyy_loc_att[NGLL3];
  __shared__ realw s_dummyz_loc_att[NGLL3];

  __shared__ realw s_tempx1[NGLL3];
  __shared__ realw s_tempx2[NGLL3];
  __shared__ realw s_tempx3[NGLL3];

  __shared__ realw s_tempy1[NGLL3];
  __shared__ realw s_tempy2[NGLL3];
  __shared__ realw s_tempy3[NGLL3];

  __shared__ realw s_tempz1[NGLL3];
  __shared__ realw s_tempz2[NGLL3];
  __shared__ realw s_tempz3[NGLL3];

  __shared__ realw sh_hprime_xx[NGLL2];

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

    iglob = d_ibool[working_element*NGLL3 + tx]-1;

#ifdef USE_TEXTURES_FIELDS
    s_dummyx_loc[tx] = texfetch_displ<FORWARD_OR_ADJOINT>(iglob*3);
    s_dummyy_loc[tx] = texfetch_displ<FORWARD_OR_ADJOINT>(iglob*3 + 1);
    s_dummyz_loc[tx] = texfetch_displ<FORWARD_OR_ADJOINT>(iglob*3 + 2);
#else
    // changing iglob indexing to match fortran row changes fast style
    s_dummyx_loc[tx] = d_displ[iglob*3];
    s_dummyy_loc[tx] = d_displ[iglob*3 + 1];
    s_dummyz_loc[tx] = d_displ[iglob*3 + 2];
#endif

  // JC JC here we will need to add GPU support for the new C-PML routines

  // attenuation
  // use first order Taylor expansion of displacement for local storage of stresses
  // at this current time step, to fix attenuation in a consistent way
#ifdef USE_TEXTURES_FIELDS
  s_dummyx_loc_att[tx] = s_dummyx_loc[tx] + d_deltat * texfetch_veloc<FORWARD_OR_ADJOINT>(iglob*3);
  s_dummyy_loc_att[tx] = s_dummyy_loc[tx] + d_deltat * texfetch_veloc<FORWARD_OR_ADJOINT>(iglob*3 + 1);
  s_dummyz_loc_att[tx] = s_dummyz_loc[tx] + d_deltat * texfetch_veloc<FORWARD_OR_ADJOINT>(iglob*3 + 2);
#else
  s_dummyx_loc_att[tx] = s_dummyx_loc[tx] + d_deltat * d_veloc[iglob*3];
  s_dummyy_loc_att[tx] = s_dummyy_loc[tx] + d_deltat * d_veloc[iglob*3 + 1];
  s_dummyz_loc_att[tx] = s_dummyz_loc[tx] + d_deltat * d_veloc[iglob*3 + 2];
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

      //assumes that hprime_xx = hprime_yy = hprime_zz
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

    // JC JC here we will need to add GPU support for the new C-PML routines

    // attenuation
    // temporary variables used for fixing attenuation in a consistent way
    tempx1l_att = 0.f;
    tempx2l_att = 0.f;
    tempx3l_att = 0.f;

    tempy1l_att = 0.f;
    tempy2l_att = 0.f;
    tempy3l_att = 0.f;

    tempz1l_att = 0.f;
    tempz2l_att = 0.f;
    tempz3l_att = 0.f;

    for (l=0;l<NGLLX;l++) {
      hp1 = sh_hprime_xx[l*NGLLX+I];
      offset = K*NGLL2+J*NGLLX+l;
      tempx1l_att += s_dummyx_loc_att[offset]*hp1;
      tempy1l_att += s_dummyy_loc_att[offset]*hp1;
      tempz1l_att += s_dummyz_loc_att[offset]*hp1;

      hp2 = sh_hprime_xx[l*NGLLX+J];
      offset = K*NGLL2+l*NGLLX+I;
      tempx2l_att += s_dummyx_loc_att[offset]*hp2;
      tempy2l_att += s_dummyy_loc_att[offset]*hp2;
      tempz2l_att += s_dummyz_loc_att[offset]*hp2;

      hp3 = sh_hprime_xx[l*NGLLX+K];
      offset = l*NGLL2+J*NGLLX+I;
      tempx3l_att += s_dummyx_loc_att[offset]*hp3;
      tempy3l_att += s_dummyy_loc_att[offset]*hp3;
      tempz3l_att += s_dummyz_loc_att[offset]*hp3;
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

    // JC JC here we will need to add GPU support for the new C-PML routines

    // attenuation
    // temporary variables used for fixing attenuation in a consistent way
    tempx1l_att = s_dummyx_loc_att[K*NGLL2+J*NGLLX]*d_hprime_xx[I]
                  + s_dummyx_loc_att[K*NGLL2+J*NGLLX+1]*d_hprime_xx[NGLLX+I]
                  + s_dummyx_loc_att[K*NGLL2+J*NGLLX+2]*d_hprime_xx[2*NGLLX+I]
                  + s_dummyx_loc_att[K*NGLL2+J*NGLLX+3]*d_hprime_xx[3*NGLLX+I]
                  + s_dummyx_loc_att[K*NGLL2+J*NGLLX+4]*d_hprime_xx[4*NGLLX+I];

    tempy1l_att = s_dummyy_loc_att[K*NGLL2+J*NGLLX]*d_hprime_xx[I]
                  + s_dummyy_loc_att[K*NGLL2+J*NGLLX+1]*d_hprime_xx[NGLLX+I]
                  + s_dummyy_loc_att[K*NGLL2+J*NGLLX+2]*d_hprime_xx[2*NGLLX+I]
                  + s_dummyy_loc_att[K*NGLL2+J*NGLLX+3]*d_hprime_xx[3*NGLLX+I]
                  + s_dummyy_loc_att[K*NGLL2+J*NGLLX+4]*d_hprime_xx[4*NGLLX+I];

    tempz1l_att = s_dummyz_loc_att[K*NGLL2+J*NGLLX]*d_hprime_xx[I]
                  + s_dummyz_loc_att[K*NGLL2+J*NGLLX+1]*d_hprime_xx[NGLLX+I]
                  + s_dummyz_loc_att[K*NGLL2+J*NGLLX+2]*d_hprime_xx[2*NGLLX+I]
                  + s_dummyz_loc_att[K*NGLL2+J*NGLLX+3]*d_hprime_xx[3*NGLLX+I]
                  + s_dummyz_loc_att[K*NGLL2+J*NGLLX+4]*d_hprime_xx[4*NGLLX+I];

    tempx2l_att = s_dummyx_loc_att[K*NGLL2+I]*d_hprime_xx[J]
                  + s_dummyx_loc_att[K*NGLL2+NGLLX+I]*d_hprime_xx[NGLLX+J]
                  + s_dummyx_loc_att[K*NGLL2+2*NGLLX+I]*d_hprime_xx[2*NGLLX+J]
                  + s_dummyx_loc_att[K*NGLL2+3*NGLLX+I]*d_hprime_xx[3*NGLLX+J]
                  + s_dummyx_loc_att[K*NGLL2+4*NGLLX+I]*d_hprime_xx[4*NGLLX+J];

    tempy2l_att = s_dummyy_loc_att[K*NGLL2+I]*d_hprime_xx[J]
                  + s_dummyy_loc_att[K*NGLL2+NGLLX+I]*d_hprime_xx[NGLLX+J]
                  + s_dummyy_loc_att[K*NGLL2+2*NGLLX+I]*d_hprime_xx[2*NGLLX+J]
                  + s_dummyy_loc_att[K*NGLL2+3*NGLLX+I]*d_hprime_xx[3*NGLLX+J]
                  + s_dummyy_loc_att[K*NGLL2+4*NGLLX+I]*d_hprime_xx[4*NGLLX+J];

    tempz2l_att = s_dummyz_loc_att[K*NGLL2+I]*d_hprime_xx[J]
                  + s_dummyz_loc_att[K*NGLL2+NGLLX+I]*d_hprime_xx[NGLLX+J]
                  + s_dummyz_loc_att[K*NGLL2+2*NGLLX+I]*d_hprime_xx[2*NGLLX+J]
                  + s_dummyz_loc_att[K*NGLL2+3*NGLLX+I]*d_hprime_xx[3*NGLLX+J]
                  + s_dummyz_loc_att[K*NGLL2+4*NGLLX+I]*d_hprime_xx[4*NGLLX+J];

    tempx3l_att = s_dummyx_loc_att[J*NGLLX+I]*d_hprime_xx[K]
                  + s_dummyx_loc_att[NGLL2+J*NGLLX+I]*d_hprime_xx[NGLLX+K]
                  + s_dummyx_loc_att[2*NGLL2+J*NGLLX+I]*d_hprime_xx[2*NGLLX+K]
                  + s_dummyx_loc_att[3*NGLL2+J*NGLLX+I]*d_hprime_xx[3*NGLLX+K]
                  + s_dummyx_loc_att[4*NGLL2+J*NGLLX+I]*d_hprime_xx[4*NGLLX+K];

    tempy3l_att = s_dummyy_loc_att[J*NGLLX+I]*d_hprime_xx[K]
                  + s_dummyy_loc_att[NGLL2+J*NGLLX+I]*d_hprime_xx[NGLLX+K]
                  + s_dummyy_loc_att[2*NGLL2+J*NGLLX+I]*d_hprime_xx[2*NGLLX+K]
                  + s_dummyy_loc_att[3*NGLL2+J*NGLLX+I]*d_hprime_xx[3*NGLLX+K]
                  + s_dummyy_loc_att[4*NGLL2+J*NGLLX+I]*d_hprime_xx[4*NGLLX+K];

    tempz3l_att = s_dummyz_loc_att[J*NGLLX+I]*d_hprime_xx[K]
                  + s_dummyz_loc_att[NGLL2+J*NGLLX+I]*d_hprime_xx[NGLLX+K]
                  + s_dummyz_loc_att[2*NGLL2+J*NGLLX+I]*d_hprime_xx[2*NGLLX+K]
                  + s_dummyz_loc_att[3*NGLL2+J*NGLLX+I]*d_hprime_xx[3*NGLLX+K]
                  + s_dummyz_loc_att[4*NGLL2+J*NGLLX+I]*d_hprime_xx[4*NGLLX+K];

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

    // JC JC here we will need to add GPU support for the new C-PML routines

    // precompute some sums to save CPU time
    duxdxl_plus_duydyl = duxdxl + duydyl;
    duxdxl_plus_duzdzl = duxdxl + duzdzl;
    duydyl_plus_duzdzl = duydyl + duzdzl;
    duxdyl_plus_duydxl = duxdyl + duydxl;
    duzdxl_plus_duxdzl = duzdxl + duxdzl;
    duzdyl_plus_duydzl = duzdyl + duydzl;

    // JC JC here we will need to add GPU support for the new C-PML routines

    // attenuation
    // temporary variables used for fixing attenuation in a consistent way
    duxdxl_att = xixl*tempx1l_att + etaxl*tempx2l_att + gammaxl*tempx3l_att;
    duxdyl_att = xiyl*tempx1l_att + etayl*tempx2l_att + gammayl*tempx3l_att;
    duxdzl_att = xizl*tempx1l_att + etazl*tempx2l_att + gammazl*tempx3l_att;

    duydxl_att = xixl*tempy1l_att + etaxl*tempy2l_att + gammaxl*tempy3l_att;
    duydyl_att = xiyl*tempy1l_att + etayl*tempy2l_att + gammayl*tempy3l_att;
    duydzl_att = xizl*tempy1l_att + etazl*tempy2l_att + gammazl*tempy3l_att;

    duzdxl_att = xixl*tempz1l_att + etaxl*tempz2l_att + gammaxl*tempz3l_att;
    duzdyl_att = xiyl*tempz1l_att + etayl*tempz2l_att + gammayl*tempz3l_att;
    duzdzl_att = xizl*tempz1l_att + etazl*tempz2l_att + gammazl*tempz3l_att;

    // precompute some sums to save CPU time
    duxdyl_plus_duydxl_att = duxdyl_att + duydxl_att;
    duzdxl_plus_duxdzl_att = duzdxl_att + duxdzl_att;
    duzdyl_plus_duydzl_att = duzdyl_att + duydzl_att;

    // attenuation
    // computes deviatoric strain attenuation and/or for kernel calculations
    realw templ = 0.33333333333333333333f * (duxdxl_att + duydyl_att + duzdzl_att); // 1./3. = 0.33333
    // local storage: stresses at this current time step
    epsilondev_xx_loc = duxdxl_att - templ;
    epsilondev_yy_loc = duydyl_att - templ;
    epsilondev_xy_loc = 0.5f * duxdyl_plus_duydxl_att;
    epsilondev_xz_loc = 0.5f * duzdxl_plus_duxdzl_att;
    epsilondev_yz_loc = 0.5f * duzdyl_plus_duydzl_att;

    if(SIMULATION_TYPE == 3) {
      epsilon_trace_over_3[tx + working_element*NGLL3] = templ;
    }

    // compute elements with an elastic isotropic rheology
    kappal = d_kappav[offset];
    mul = d_muv[offset];

    // attenuation
    // use unrelaxed parameters if attenuation
    mul  = mul * one_minus_sum_beta[tx+working_element*NGLL3]; // (i,j,k,ispec)

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

    // attenuation
    // subtracts memory variables if attenuation
    compute_element_att_stress(tx,working_element,NSPEC,
                               R_xx,R_yy,R_xy,R_xz,R_yz,
                               &sigma_xx,&sigma_yy,&sigma_zz,&sigma_xy,&sigma_xz,&sigma_yz);

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

  // JC JC here we will need to add GPU support for the new C-PML routines

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

      // assumes hprimewgll_xx == hprimewgll_yy == hprimewgll_zz
      fac2 = d_hprimewgll_xx[J*NGLLX+l];
      offset = K*NGLL2+l*NGLLX+I;
      tempx2l += s_tempx2[offset]*fac2;
      tempy2l += s_tempy2[offset]*fac2;
      tempz2l += s_tempz2[offset]*fac2;

      fac3 = d_hprimewgll_xx[K*NGLLX+l];
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

    tempx2l = s_tempx2[K*NGLL2+I]*d_hprimewgll_xx[J*NGLLX]
            + s_tempx2[K*NGLL2+NGLLX+I]*d_hprimewgll_xx[J*NGLLX+1]
            + s_tempx2[K*NGLL2+2*NGLLX+I]*d_hprimewgll_xx[J*NGLLX+2]
            + s_tempx2[K*NGLL2+3*NGLLX+I]*d_hprimewgll_xx[J*NGLLX+3]
            + s_tempx2[K*NGLL2+4*NGLLX+I]*d_hprimewgll_xx[J*NGLLX+4];

    tempy2l = s_tempy2[K*NGLL2+I]*d_hprimewgll_xx[J*NGLLX]
            + s_tempy2[K*NGLL2+NGLLX+I]*d_hprimewgll_xx[J*NGLLX+1]
            + s_tempy2[K*NGLL2+2*NGLLX+I]*d_hprimewgll_xx[J*NGLLX+2]
            + s_tempy2[K*NGLL2+3*NGLLX+I]*d_hprimewgll_xx[J*NGLLX+3]
            + s_tempy2[K*NGLL2+4*NGLLX+I]*d_hprimewgll_xx[J*NGLLX+4];

    tempz2l = s_tempz2[K*NGLL2+I]*d_hprimewgll_xx[J*NGLLX]
            + s_tempz2[K*NGLL2+NGLLX+I]*d_hprimewgll_xx[J*NGLLX+1]
            + s_tempz2[K*NGLL2+2*NGLLX+I]*d_hprimewgll_xx[J*NGLLX+2]
            + s_tempz2[K*NGLL2+3*NGLLX+I]*d_hprimewgll_xx[J*NGLLX+3]
            + s_tempz2[K*NGLL2+4*NGLLX+I]*d_hprimewgll_xx[J*NGLLX+4];

    tempx3l = s_tempx3[J*NGLLX+I]*d_hprimewgll_xx[K*NGLLX]
            + s_tempx3[NGLL2+J*NGLLX+I]*d_hprimewgll_xx[K*NGLLX+1]
            + s_tempx3[2*NGLL2+J*NGLLX+I]*d_hprimewgll_xx[K*NGLLX+2]
            + s_tempx3[3*NGLL2+J*NGLLX+I]*d_hprimewgll_xx[K*NGLLX+3]
            + s_tempx3[4*NGLL2+J*NGLLX+I]*d_hprimewgll_xx[K*NGLLX+4];

    tempy3l = s_tempy3[J*NGLLX+I]*d_hprimewgll_xx[K*NGLLX]
            + s_tempy3[NGLL2+J*NGLLX+I]*d_hprimewgll_xx[K*NGLLX+1]
            + s_tempy3[2*NGLL2+J*NGLLX+I]*d_hprimewgll_xx[K*NGLLX+2]
            + s_tempy3[3*NGLL2+J*NGLLX+I]*d_hprimewgll_xx[K*NGLLX+3]
            + s_tempy3[4*NGLL2+J*NGLLX+I]*d_hprimewgll_xx[K*NGLLX+4];

    tempz3l = s_tempz3[J*NGLLX+I]*d_hprimewgll_xx[K*NGLLX]
            + s_tempz3[NGLL2+J*NGLLX+I]*d_hprimewgll_xx[K*NGLLX+1]
            + s_tempz3[2*NGLL2+J*NGLLX+I]*d_hprimewgll_xx[K*NGLLX+2]
            + s_tempz3[3*NGLL2+J*NGLLX+I]*d_hprimewgll_xx[K*NGLLX+3]
            + s_tempz3[4*NGLL2+J*NGLLX+I]*d_hprimewgll_xx[K*NGLLX+4];

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
    d_accel[iglob*3]     = texfetch_accel<FORWARD_OR_ADJOINT>(iglob*3) + sum_terms1;
    d_accel[iglob*3 + 1] = texfetch_accel<FORWARD_OR_ADJOINT>(iglob*3 + 1) + sum_terms2;
    d_accel[iglob*3 + 2] = texfetch_accel<FORWARD_OR_ADJOINT>(iglob*3 + 2) + sum_terms3;
#else
    d_accel[iglob*3]     += sum_terms1;
    d_accel[iglob*3 + 1] += sum_terms2;
    d_accel[iglob*3 + 2] += sum_terms3;
#endif // USE_TEXTURES_FIELDS

    // JC JC here we will need to add GPU support for the new C-PML routines

#else // MESH_COLORING

    //mesh coloring
    if( use_mesh_coloring_gpu ){

      // no atomic operation needed, colors don't share global points between elements
#ifdef USE_TEXTURES_FIELDS
      d_accel[iglob*3]     = texfetch_accel<FORWARD_OR_ADJOINT>(iglob*3) + sum_terms1;
      d_accel[iglob*3 + 1] = texfetch_accel<FORWARD_OR_ADJOINT>(iglob*3 + 1) + sum_terms2;
      d_accel[iglob*3 + 2] = texfetch_accel<FORWARD_OR_ADJOINT>(iglob*3 + 2) + sum_terms3;
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
      // w/out atomic update
      //d_accel[iglob*3]     += sum_terms1;
      //d_accel[iglob*3 + 1] += sum_terms2;
      //d_accel[iglob*3 + 2] += sum_terms3;

      atomicAdd(&d_accel[iglob*3], sum_terms1);
      atomicAdd(&d_accel[iglob*3+1], sum_terms2);
      atomicAdd(&d_accel[iglob*3+2], sum_terms3);

    } // if(use_mesh_coloring_gpu)

#endif // MESH_COLORING

    // attenuation
    // update memory variables based upon the Runge-Kutta scheme
    compute_element_att_memory(tx,working_element,NSPEC,
                              d_muv,
                              factor_common,alphaval,betaval,gammaval,
                              R_xx,R_yy,R_xy,R_xz,R_yz,
                              epsilondev_xx,epsilondev_yy,epsilondev_xy,epsilondev_xz,epsilondev_yz,
                              epsilondev_xx_loc,epsilondev_yy_loc,epsilondev_xy_loc,epsilondev_xz_loc,epsilondev_yz_loc);

    // save deviatoric strain for Runge-Kutta scheme
    int ijk_ispec = tx + working_element*NGLL3;

    // fortran: epsilondev_xx(:,:,:,ispec) = epsilondev_xx_loc(:,:,:)
    epsilondev_xx[ijk_ispec] = epsilondev_xx_loc;
    epsilondev_yy[ijk_ispec] = epsilondev_yy_loc;
    epsilondev_xy[ijk_ispec] = epsilondev_xy_loc;
    epsilondev_xz[ijk_ispec] = epsilondev_xz_loc;
    epsilondev_yz[ijk_ispec] = epsilondev_yz_loc;

  } // if(active)

  // JC JC here we will need to add GPU support for the new C-PML routines

} // kernel_2_att_impl()


/* ----------------------------------------------------------------------------------------------- */

void Kernel_2(int nb_blocks_to_compute,Mesh* mp,int d_iphase,realw d_deltat,
              int COMPUTE_AND_STORE_STRAIN,
              int ATTENUATION,int ANISOTROPY,
              int* d_ibool,
              realw* d_xix,realw* d_xiy,realw* d_xiz,
              realw* d_etax,realw* d_etay,realw* d_etaz,
              realw* d_gammax,realw* d_gammay,realw* d_gammaz,
              realw* d_kappav,
              realw* d_muv,
              realw* d_epsilondev_xx,realw* d_epsilondev_yy,realw* d_epsilondev_xy,
              realw* d_epsilondev_xz,realw* d_epsilondev_yz,
              realw* d_epsilon_trace_over_3,
              realw* d_one_minus_sum_beta,
              realw* d_factor_common,
              realw* d_R_xx,realw* d_R_yy,realw* d_R_xy,
              realw* d_R_xz,realw* d_R_yz,
              realw* d_b_epsilondev_xx,realw* d_b_epsilondev_yy,realw* d_b_epsilondev_xy,
              realw* d_b_epsilondev_xz,realw* d_b_epsilondev_yz,
              realw* d_b_epsilon_trace_over_3,
              realw* d_b_R_xx,realw* d_b_R_yy,realw* d_b_R_xy,
              realw* d_b_R_xz,realw* d_b_R_yz,
              realw* d_c11store,realw* d_c12store,realw* d_c13store,
              realw* d_c14store,realw* d_c15store,realw* d_c16store,
              realw* d_c22store,realw* d_c23store,realw* d_c24store,
              realw* d_c25store,realw* d_c26store,realw* d_c33store,
              realw* d_c34store,realw* d_c35store,realw* d_c36store,
              realw* d_c44store,realw* d_c45store,realw* d_c46store,
              realw* d_c55store,realw* d_c56store,realw* d_c66store,
              realw* d_rhostore){

  TRACE("\tKernel_2");

#ifdef ENABLE_VERY_SLOW_ERROR_CHECKING
  exit_on_cuda_error("before kernel Kernel 2");
#endif

  /* if the grid can handle the number of blocks, we let it be 1D */
  /* grid_2_x = nb_elem_color; */
  /* nb_elem_color is just how many blocks we are computing now */

  int blocksize = NGLL3_PADDED;

  int num_blocks_x, num_blocks_y;
  get_blocks_xy(nb_blocks_to_compute,&num_blocks_x,&num_blocks_y);

  dim3 grid(num_blocks_x,num_blocks_y);
  dim3 threads(blocksize,1,1);

  // Cuda timing
  // cudaEvent_t start, stop;
  // realw time;
  // cudaEventCreate(&start);
  // cudaEventCreate(&stop);
  // cudaEventRecord( start, 0 );

  if( ATTENUATION ){
    // debug
    //printf("Running Kernel_2 with attenuation\n");

    // compute kernels with attenuation
    // forward wavefields -> FORWARD_OR_ADJOINT == 1
    Kernel_2_att_impl<1><<<grid,threads,0,mp->compute_stream>>>(nb_blocks_to_compute,
                                                                mp->NGLOB_AB,
                                                                d_ibool,
                                                                mp->d_phase_ispec_inner_elastic,
                                                                mp->num_phase_ispec_elastic,
                                                                d_iphase,
                                                                mp->use_mesh_coloring_gpu,
                                                                d_deltat,
                                                                mp->d_displ,mp->d_veloc,mp->d_accel,
                                                                d_xix, d_xiy, d_xiz,
                                                                d_etax, d_etay, d_etaz,
                                                                d_gammax, d_gammay, d_gammaz,
                                                                mp->d_hprime_xx,
                                                                mp->d_hprimewgll_xx,
                                                                mp->d_wgllwgll_xy, mp->d_wgllwgll_xz, mp->d_wgllwgll_yz,
                                                                d_kappav, d_muv,
                                                                d_epsilondev_xx,d_epsilondev_yy,d_epsilondev_xy,
                                                                d_epsilondev_xz,d_epsilondev_yz,
                                                                d_epsilon_trace_over_3,
                                                                mp->simulation_type,
                                                                mp->NSPEC_AB,
                                                                d_one_minus_sum_beta,
                                                                d_factor_common,
                                                                d_R_xx,d_R_yy,d_R_xy,d_R_xz,d_R_yz,
                                                                mp->d_alphaval,mp->d_betaval,mp->d_gammaval,
                                                                ANISOTROPY,
                                                                d_c11store,d_c12store,d_c13store,
                                                                d_c14store,d_c15store,d_c16store,
                                                                d_c22store,d_c23store,d_c24store,
                                                                d_c25store,d_c26store,d_c33store,
                                                                d_c34store,d_c35store,d_c36store,
                                                                d_c44store,d_c45store,d_c46store,
                                                                d_c55store,d_c56store,d_c66store,
                                                                mp->gravity,
                                                                mp->d_minus_g,
                                                                mp->d_minus_deriv_gravity,
                                                                d_rhostore,
                                                                mp->d_wgll_cube);

    if(mp->simulation_type == 3) {
      // backward/reconstructed wavefields -> FORWARD_OR_ADJOINT == 3
      Kernel_2_att_impl<3><<< grid,threads,0,mp->compute_stream>>>(nb_blocks_to_compute,
                                                                   mp->NGLOB_AB,
                                                                   d_ibool,
                                                                   mp->d_phase_ispec_inner_elastic,
                                                                   mp->num_phase_ispec_elastic,
                                                                   d_iphase,
                                                                   mp->use_mesh_coloring_gpu,
                                                                   d_deltat,
                                                                   mp->d_b_displ,mp->d_b_veloc,mp->d_b_accel,
                                                                   d_xix, d_xiy, d_xiz,
                                                                   d_etax, d_etay, d_etaz,
                                                                   d_gammax, d_gammay, d_gammaz,
                                                                   mp->d_hprime_xx,
                                                                   mp->d_hprimewgll_xx,
                                                                   mp->d_wgllwgll_xy, mp->d_wgllwgll_xz, mp->d_wgllwgll_yz,
                                                                   d_kappav, d_muv,
                                                                   d_b_epsilondev_xx,d_b_epsilondev_yy,d_b_epsilondev_xy,
                                                                   d_b_epsilondev_xz,d_b_epsilondev_yz,
                                                                   d_b_epsilon_trace_over_3,
                                                                   mp->simulation_type,
                                                                   mp->NSPEC_AB,
                                                                   d_one_minus_sum_beta,
                                                                   d_factor_common,
                                                                   d_b_R_xx,d_b_R_yy,d_b_R_xy,d_b_R_xz,d_b_R_yz,
                                                                   mp->d_b_alphaval,mp->d_b_betaval,mp->d_b_gammaval,
                                                                   ANISOTROPY,
                                                                   d_c11store,d_c12store,d_c13store,
                                                                   d_c14store,d_c15store,d_c16store,
                                                                   d_c22store,d_c23store,d_c24store,
                                                                   d_c25store,d_c26store,d_c33store,
                                                                   d_c34store,d_c35store,d_c36store,
                                                                   d_c44store,d_c45store,d_c46store,
                                                                   d_c55store,d_c56store,d_c66store,
                                                                   mp->gravity,
                                                                   mp->d_minus_g,
                                                                   mp->d_minus_deriv_gravity,
                                                                   d_rhostore,
                                                                   mp->d_wgll_cube);
    }
  }else{
    // debug
    //printf("Running Kernel_2 without attenuation\n");

    // compute kernels without attenuation
    // forward wavefields -> FORWARD_OR_ADJOINT == 1
    Kernel_2_noatt_impl<1><<<grid,threads,0,mp->compute_stream>>>(nb_blocks_to_compute,
                                                                  mp->NGLOB_AB,
                                                                  d_ibool,
                                                                  mp->d_phase_ispec_inner_elastic,mp->num_phase_ispec_elastic,
                                                                  d_iphase,
                                                                  mp->use_mesh_coloring_gpu,
                                                                  mp->d_displ,mp->d_veloc,mp->d_accel,
                                                                  d_xix, d_xiy, d_xiz,
                                                                  d_etax, d_etay, d_etaz,
                                                                  d_gammax, d_gammay, d_gammaz,
                                                                  mp->d_hprime_xx,
                                                                  mp->d_hprimewgll_xx,
                                                                  mp->d_wgllwgll_xy, mp->d_wgllwgll_xz, mp->d_wgllwgll_yz,
                                                                  d_kappav, d_muv,
                                                                  COMPUTE_AND_STORE_STRAIN,
                                                                  d_epsilondev_xx,d_epsilondev_yy,d_epsilondev_xy,
                                                                  d_epsilondev_xz,d_epsilondev_yz,
                                                                  d_epsilon_trace_over_3,
                                                                  mp->simulation_type,
                                                                  mp->NSPEC_AB,
                                                                  d_one_minus_sum_beta,d_factor_common,
                                                                  d_R_xx,d_R_yy,d_R_xy,d_R_xz,d_R_yz,
                                                                  mp->d_alphaval,mp->d_betaval,mp->d_gammaval,
                                                                  ANISOTROPY,
                                                                  d_c11store,d_c12store,d_c13store,
                                                                  d_c14store,d_c15store,d_c16store,
                                                                  d_c22store,d_c23store,d_c24store,
                                                                  d_c25store,d_c26store,d_c33store,
                                                                  d_c34store,d_c35store,d_c36store,
                                                                  d_c44store,d_c45store,d_c46store,
                                                                  d_c55store,d_c56store,d_c66store,
                                                                  mp->gravity,
                                                                  mp->d_minus_g,
                                                                  mp->d_minus_deriv_gravity,
                                                                  d_rhostore,
                                                                  mp->d_wgll_cube );

    // backward/reconstructed wavefield
    if(mp->simulation_type == 3) {
      // backward/reconstructed wavefields -> FORWARD_OR_ADJOINT == 3
      Kernel_2_noatt_impl<3><<< grid,threads,0,mp->compute_stream>>>(nb_blocks_to_compute,
                                                                     mp->NGLOB_AB,
                                                                     d_ibool,
                                                                     mp->d_phase_ispec_inner_elastic,mp->num_phase_ispec_elastic,
                                                                     d_iphase,
                                                                     mp->use_mesh_coloring_gpu,
                                                                     mp->d_b_displ,mp->d_b_veloc,mp->d_b_accel,
                                                                     d_xix, d_xiy, d_xiz,
                                                                     d_etax, d_etay, d_etaz,
                                                                     d_gammax, d_gammay, d_gammaz,
                                                                     mp->d_hprime_xx,
                                                                     mp->d_hprimewgll_xx,
                                                                     mp->d_wgllwgll_xy, mp->d_wgllwgll_xz, mp->d_wgllwgll_yz,
                                                                     d_kappav, d_muv,
                                                                     COMPUTE_AND_STORE_STRAIN,
                                                                     d_b_epsilondev_xx,d_b_epsilondev_yy,d_b_epsilondev_xy,
                                                                     d_b_epsilondev_xz,d_b_epsilondev_yz,
                                                                     d_b_epsilon_trace_over_3,
                                                                     mp->simulation_type,
                                                                     mp->NSPEC_AB,
                                                                     d_one_minus_sum_beta,d_factor_common,
                                                                     d_b_R_xx,d_b_R_yy,d_b_R_xy,d_b_R_xz,d_b_R_yz,
                                                                     mp->d_b_alphaval,mp->d_b_betaval,mp->d_b_gammaval,
                                                                     ANISOTROPY,
                                                                     d_c11store,d_c12store,d_c13store,
                                                                     d_c14store,d_c15store,d_c16store,
                                                                     d_c22store,d_c23store,d_c24store,
                                                                     d_c25store,d_c26store,d_c33store,
                                                                     d_c34store,d_c35store,d_c36store,
                                                                     d_c44store,d_c45store,d_c46store,
                                                                     d_c55store,d_c56store,d_c66store,
                                                                     mp->gravity,
                                                                     mp->d_minus_g,
                                                                     mp->d_minus_deriv_gravity,
                                                                     d_rhostore,
                                                                     mp->d_wgll_cube );
    }
  }
  // cudaEventRecord( stop, 0 );
  // cudaEventSynchronize( stop );
  // cudaEventElapsedTime( &time, start, stop );
  // cudaEventDestroy( start );
  // cudaEventDestroy( stop );
  // printf("Kernel2 Execution Time: %f ms\n",time);

  // cudaThreadSynchronize(); //
  // LOG("Kernel 2 finished"); //
#ifdef ENABLE_VERY_SLOW_ERROR_CHECKING
  exit_on_cuda_error("Kernel_2_impl");
#endif
}

/* ----------------------------------------------------------------------------------------------- */


extern "C"
void FC_FUNC_(compute_forces_viscoelastic_cuda,
              COMPUTE_FORCES_VISCOELASTIC_CUDA)(long* Mesh_pointer,
                                                int* iphase,
                                                realw* deltat,
                                                int* nspec_outer_elastic,
                                                int* nspec_inner_elastic,
                                                int* COMPUTE_AND_STORE_STRAIN,
                                                int* ATTENUATION,
                                                int* ANISOTROPY) {

  TRACE("\tcompute_forces_viscoelastic_cuda");
  // EPIK_TRACER("compute_forces_viscoelastic_cuda");
  //printf("Running compute_forces\n");
  //double start_time = get_time();

  Mesh* mp = (Mesh*)(*Mesh_pointer); // get Mesh from fortran integer wrapper

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
    int offset,offset_nonpadded,offset_nonpadded_att2;

    // sets up color loop
    if( *iphase == 1 ){
      // outer elements
      nb_colors = mp->num_colors_outer_elastic;
      istart = 0;

      // array offsets
      offset = 0;
      offset_nonpadded = 0;
      offset_nonpadded_att2 = 0;
    }else{
      // inner elements (start after outer elements)
      nb_colors = mp->num_colors_outer_elastic + mp->num_colors_inner_elastic;
      istart = mp->num_colors_outer_elastic;

      // array offsets
      offset = (*nspec_outer_elastic) * NGLL3_PADDED;
      offset_nonpadded = (*nspec_outer_elastic) * NGLL3;
      offset_nonpadded_att2 = (*nspec_outer_elastic) * NGLL3 * N_SLS;
    }

    // loops over colors
    for(int icolor = istart; icolor < nb_colors; icolor++){

      nb_blocks_to_compute = mp->h_num_elem_colors_elastic[icolor];

      // checks
      //if( nb_blocks_to_compute <= 0 ){
      //  printf("error number of elastic color blocks: %d -- color = %d \n",nb_blocks_to_compute,icolor);
      //  exit(EXIT_FAILURE);
      //}

      Kernel_2(nb_blocks_to_compute,mp,*iphase,*deltat,
               *COMPUTE_AND_STORE_STRAIN,
               *ATTENUATION,*ANISOTROPY,
               mp->d_ibool + offset_nonpadded,
               mp->d_xix + offset,mp->d_xiy + offset,mp->d_xiz + offset,
               mp->d_etax + offset,mp->d_etay + offset,mp->d_etaz + offset,
               mp->d_gammax + offset,mp->d_gammay + offset,mp->d_gammaz + offset,
               mp->d_kappav + offset,
               mp->d_muv + offset,
               mp->d_epsilondev_xx + offset_nonpadded,mp->d_epsilondev_yy + offset_nonpadded,mp->d_epsilondev_xy + offset_nonpadded,
               mp->d_epsilondev_xz + offset_nonpadded,mp->d_epsilondev_yz + offset_nonpadded,
               mp->d_epsilon_trace_over_3 + offset_nonpadded,
               mp->d_one_minus_sum_beta + offset_nonpadded,
               mp->d_factor_common + offset_nonpadded_att2,
               mp->d_R_xx + offset_nonpadded,mp->d_R_yy + offset_nonpadded,mp->d_R_xy + offset_nonpadded,
               mp->d_R_xz + offset_nonpadded,mp->d_R_yz + offset_nonpadded,
               mp->d_b_epsilondev_xx + offset_nonpadded,mp->d_b_epsilondev_yy + offset_nonpadded,mp->d_b_epsilondev_xy + offset_nonpadded,
               mp->d_b_epsilondev_xz + offset_nonpadded,mp->d_b_epsilondev_yz + offset_nonpadded,
               mp->d_b_epsilon_trace_over_3 + offset_nonpadded,
               mp->d_b_R_xx + offset_nonpadded,mp->d_b_R_yy + offset_nonpadded,mp->d_b_R_xy + offset_nonpadded,
               mp->d_b_R_xz + offset_nonpadded,mp->d_b_R_yz + offset_nonpadded,
               mp->d_c11store + offset,mp->d_c12store + offset,mp->d_c13store + offset,
               mp->d_c14store + offset,mp->d_c15store + offset,mp->d_c16store + offset,
               mp->d_c22store + offset,mp->d_c23store + offset,mp->d_c24store + offset,
               mp->d_c25store + offset,mp->d_c26store + offset,mp->d_c33store + offset,
               mp->d_c34store + offset,mp->d_c35store + offset,mp->d_c36store + offset,
               mp->d_c44store + offset,mp->d_c45store + offset,mp->d_c46store + offset,
               mp->d_c55store + offset,mp->d_c56store + offset,mp->d_c66store + offset,
               mp->d_rhostore + offset);

      // for padded and aligned arrays
      offset += nb_blocks_to_compute * NGLL3_PADDED;
      // for no-aligned arrays
      offset_nonpadded += nb_blocks_to_compute * NGLL3;
      // for factor_common array
      offset_nonpadded_att2 += nb_blocks_to_compute * NGLL3 * N_SLS;

      //note: we use the same stream, so kernels are executed one after the other
      //      thus, there should be no need to synchronize in case we run on only 1 process to avoid race-conditions

    }

  }else{
    // no mesh coloring: uses atomic updates
    Kernel_2(num_elements,mp,*iphase,*deltat,
             *COMPUTE_AND_STORE_STRAIN,
             *ATTENUATION,*ANISOTROPY,
             mp->d_ibool,
             mp->d_xix,mp->d_xiy,mp->d_xiz,
             mp->d_etax,mp->d_etay,mp->d_etaz,
             mp->d_gammax,mp->d_gammay,mp->d_gammaz,
             mp->d_kappav,
             mp->d_muv,
             mp->d_epsilondev_xx,mp->d_epsilondev_yy,mp->d_epsilondev_xy,
             mp->d_epsilondev_xz,mp->d_epsilondev_yz,
             mp->d_epsilon_trace_over_3,
             mp->d_one_minus_sum_beta,
             mp->d_factor_common,
             mp->d_R_xx,mp->d_R_yy,mp->d_R_xy,
             mp->d_R_xz,mp->d_R_yz,
             mp->d_b_epsilondev_xx,mp->d_b_epsilondev_yy,mp->d_b_epsilondev_xy,
             mp->d_b_epsilondev_xz,mp->d_b_epsilondev_yz,
             mp->d_b_epsilon_trace_over_3,
             mp->d_b_R_xx,mp->d_b_R_yy,mp->d_b_R_xy,
             mp->d_b_R_xz,mp->d_b_R_yz,
             mp->d_c11store,mp->d_c12store,mp->d_c13store,
             mp->d_c14store,mp->d_c15store,mp->d_c16store,
             mp->d_c22store,mp->d_c23store,mp->d_c24store,
             mp->d_c25store,mp->d_c26store,mp->d_c33store,
             mp->d_c34store,mp->d_c35store,mp->d_c36store,
             mp->d_c44store,mp->d_c45store,mp->d_c46store,
             mp->d_c55store,mp->d_c56store,mp->d_c66store,
             mp->d_rhostore);
  }
}

/* ----------------------------------------------------------------------------------------------- */

extern "C"
void FC_FUNC_(sync_copy_from_device,
              SYNC_copy_FROM_DEVICE)(long* Mesh_pointer,
                                     int* iphase,
                                     realw* send_buffer) {

  TRACE("sync_copy_from_device");

  Mesh* mp = (Mesh*)(*Mesh_pointer); // get Mesh from fortran integer wrapper

  // Wait until async-memcpy of outer elements is finished and start MPI.
  if( *iphase != 2 ){ exit_on_cuda_error("sync_copy_from_device must be called for iphase == 2"); }

  if( mp->size_mpi_buffer > 0 ){
    // waits for asynchronous copy to finish
    cudaStreamSynchronize(mp->copy_stream);

    // There have been problems using the pinned-memory with MPI, so
    // we copy the buffer into a non-pinned region.
    memcpy(send_buffer,mp->h_send_accel_buffer,mp->size_mpi_buffer*sizeof(float));
  }
  // memory copy is now finished, so non-blocking MPI send can proceed
}

