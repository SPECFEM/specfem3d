/*
 !=====================================================================
 !
 !                         S p e c f e m 3 D
 !                         -----------------
 !
 !    Main historical authors: Dimitri Komatitsch and Jeroen Tromp
 !                             CNRS, France
 !                      and Princeton University, USA
 !                (there are currently many more authors!)
 !                          (c) October 2017
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

#include "mesh_constants_gpu.h"

/* ----------------------------------------------------------------------------------------------- */

// ASSEMBLY - mpi data transfer between CPU-GPU

/* ----------------------------------------------------------------------------------------------- */

// prepares and transfers the inter-element edge-nodes to the host to be MPI'd
// (elements on boundary)

extern EXTERN_LANG
void FC_FUNC_(transfer_boun_accel_from_device,
              TRANSFER_BOUN_ACCEL_FROM_DEVICE)(long* Mesh_pointer,
                                               realw* send_accel_buffer,
                                               const int* FORWARD_OR_ADJOINT){
TRACE("transfer_boun_accel_from_device");

  Mesh* mp = (Mesh*)(*Mesh_pointer); //get mesh pointer out of fortran integer container

  // checks if anything to do
  if (mp->size_mpi_buffer > 0){

    int blocksize = BLOCKSIZE_TRANSFER;
    int size_padded = ((int)ceil(((double)mp->max_nibool_interfaces_ext_mesh)/((double)blocksize)))*blocksize;

    int num_blocks_x, num_blocks_y;
    get_blocks_xy(size_padded/blocksize,&num_blocks_x,&num_blocks_y);

    dim3 grid(num_blocks_x,num_blocks_y);
    dim3 threads(blocksize,1,1);

    // selects arrays
    realw* d_accel = NULL;
    realw* d_send_buffer = NULL;
    if (*FORWARD_OR_ADJOINT == 1) {
      // forward wavefield
      d_accel = mp->d_accel;
      d_send_buffer = mp->d_send_accel_buffer;
    } else if (*FORWARD_OR_ADJOINT == 3) {
      // backward/reconstructed wavefield
      d_accel = mp->d_b_accel;
      d_send_buffer = mp->d_b_send_accel_buffer;
    }

    // Cuda timing
    //gpu_event start, stop;
    //start_timing_gpu(&start,&stop);

    // fills mpi boundary buffer
#ifdef USE_CUDA
    if (run_cuda){
      prepare_boundary_accel_on_device<<<grid,threads,0,mp->compute_stream>>>(d_accel,d_send_buffer,
                                                                              mp->num_interfaces_ext_mesh,
                                                                              mp->max_nibool_interfaces_ext_mesh,
                                                                              mp->d_nibool_interfaces_ext_mesh,
                                                                              mp->d_ibool_interfaces_ext_mesh);
    }
#endif
#ifdef USE_HIP
    if (run_hip){
      hipLaunchKernelGGL(prepare_boundary_accel_on_device, dim3(grid), dim3(threads), 0, mp->compute_stream,
                                                           d_accel,d_send_buffer,
                                                           mp->num_interfaces_ext_mesh,
                                                           mp->max_nibool_interfaces_ext_mesh,
                                                           mp->d_nibool_interfaces_ext_mesh,
                                                           mp->d_ibool_interfaces_ext_mesh);
    }
#endif

    // synchronizes
    //gpuSynchronize();
    // explicitly waits until previous compute stream finishes
    // (cudaMemcpy implicitly synchronizes all other cuda operations)
    gpuStreamSynchronize(mp->compute_stream);

    // copies buffer from GPU to CPU host
    gpuMemcpy_tohost_realw(send_accel_buffer,d_send_buffer,mp->size_mpi_buffer);

    // kernel timing
    // finish timing of kernel+memcpy
    //stop_timing_gpu(&start,&stop,"prepare_boundary_accel_on_device");
  }

  GPU_ERROR_CHECKING("transfer_boun_accel_from_device");
}

/* ----------------------------------------------------------------------------------------------- */

extern EXTERN_LANG
void FC_FUNC_(transfer_boundary_from_device_a,
              TRANSFER_BOUNDARY_FROM_DEVICE_A)(long* Mesh_pointer,
                                               const int* nspec_outer_elastic) {

// asynchronous transfer from device to host

  TRACE("transfer_boundary_from_device_a");

  Mesh* mp = (Mesh*)(*Mesh_pointer); // get Mesh from fortran integer wrapper

  if (mp->size_mpi_buffer > 0){

    int blocksize = BLOCKSIZE_TRANSFER;
    int size_padded = ((int)ceil(((double)mp->max_nibool_interfaces_ext_mesh)/((double)blocksize)))*blocksize;

    int num_blocks_x, num_blocks_y;
    get_blocks_xy(size_padded/blocksize,&num_blocks_x,&num_blocks_y);

    dim3 grid(num_blocks_x,num_blocks_y);
    dim3 threads(blocksize,1,1);

    // prepares boundary buffer
#ifdef USE_CUDA
    if (run_cuda){
      prepare_boundary_accel_on_device<<<grid,threads,0,mp->compute_stream>>>(mp->d_accel,mp->d_send_accel_buffer,
                                                                              mp->num_interfaces_ext_mesh,
                                                                              mp->max_nibool_interfaces_ext_mesh,
                                                                              mp->d_nibool_interfaces_ext_mesh,
                                                                              mp->d_ibool_interfaces_ext_mesh);
    }
#endif
#ifdef USE_HIP
    if (run_hip){
      hipLaunchKernelGGL(prepare_boundary_accel_on_device, dim3(grid), dim3(threads), 0, mp->compute_stream,
                                                           mp->d_accel,mp->d_send_accel_buffer,
                                                           mp->num_interfaces_ext_mesh,
                                                           mp->max_nibool_interfaces_ext_mesh,
                                                           mp->d_nibool_interfaces_ext_mesh,
                                                           mp->d_ibool_interfaces_ext_mesh);
    }
#endif

    // waits until kernel is finished before starting async memcpy
    //gpuSynchronize();
    // waits until previous compute stream finishes
    gpuStreamSynchronize(mp->compute_stream);

    gpuMemcpyAsync_tohost_realw(mp->h_send_accel_buffer,mp->d_send_accel_buffer,mp->size_mpi_buffer,mp->copy_stream);
  }
}

/* ----------------------------------------------------------------------------------------------- */

extern EXTERN_LANG
void FC_FUNC_(transfer_boundary_to_device_a,
              TRANSFER_BOUNDARY_TO_DEVICE_A)(long* Mesh_pointer,
                                             realw* buffer_recv_vector_ext_mesh,
                                             const int* num_interfaces_ext_mesh,
                                             const int* max_nibool_interfaces_ext_mesh) {

// asynchronous transfer from host to device

  TRACE("transfer_boundary_to_device_a");

  Mesh* mp = (Mesh*)(*Mesh_pointer); //get mesh pointer out of fortran integer container

  if (mp->size_mpi_buffer > 0){
    // copy on host memory
    memcpy(mp->h_recv_accel_buffer,buffer_recv_vector_ext_mesh,mp->size_mpi_buffer*sizeof(realw));

    // asynchronous copy to GPU using copy_stream
    gpuMemcpyAsync_todevice_realw(mp->d_send_accel_buffer, mp->h_recv_accel_buffer,mp->size_mpi_buffer,mp->copy_stream);
  }
}


/* ----------------------------------------------------------------------------------------------- */

// Assembly

/* ----------------------------------------------------------------------------------------------- */

// FORWARD_OR_ADJOINT == 1 for accel, and == 3 for b_accel

extern EXTERN_LANG
void FC_FUNC_(transfer_asmbl_accel_to_device,
              TRANSFER_ASMBL_ACCEL_TO_DEVICE)(long* Mesh_pointer,
                                              realw* buffer_recv_vector_ext_mesh,
                                              const int* num_interfaces_ext_mesh,
                                              const int* max_nibool_interfaces_ext_mesh,
                                              const int* nibool_interfaces_ext_mesh,
                                              const int* ibool_interfaces_ext_mesh,
                                              const int* FORWARD_OR_ADJOINT) {
TRACE("transfer_asmbl_accel_to_device");

  Mesh* mp = (Mesh*)(*Mesh_pointer); //get mesh pointer out of fortran integer container

  if (mp->size_mpi_buffer > 0){

    //daniel: todo - check if this copy is only needed for adjoint simulation, otherwise it is called asynchronously?
    if (*FORWARD_OR_ADJOINT == 1){
      // Wait until previous copy stream finishes. We assemble while other compute kernels execute.
      gpuStreamSynchronize(mp->copy_stream);
    }
    else if (*FORWARD_OR_ADJOINT == 3){
      // explicitly synchronizes
      // (cudaMemcpy implicitly synchronizes all other cuda operations)
      gpuSynchronize();

      gpuMemcpy_todevice_realw(mp->d_b_send_accel_buffer, buffer_recv_vector_ext_mesh,mp->size_mpi_buffer);
    }

    int blocksize = BLOCKSIZE_TRANSFER;
    int size_padded = ((int)ceil(((double)mp->max_nibool_interfaces_ext_mesh)/((double)blocksize)))*blocksize;

    int num_blocks_x, num_blocks_y;
    get_blocks_xy(size_padded/blocksize,&num_blocks_x,&num_blocks_y);

    dim3 grid(num_blocks_x,num_blocks_y);
    dim3 threads(blocksize,1,1);

    // selects arrays
    realw* d_accel = NULL;
    realw* d_send_buffer = NULL;
    if (*FORWARD_OR_ADJOINT == 1) {
      // forward wavefield
      d_accel = mp->d_accel;
      d_send_buffer = mp->d_send_accel_buffer;
    } else if (*FORWARD_OR_ADJOINT == 3) {
      // backward/reconstructed wavefield
      d_accel = mp->d_b_accel;
      d_send_buffer = mp->d_b_send_accel_buffer;
    }

    //double start_time = get_time_val();
    // gpu_event start, stop;
    // realw time;
    // cudaEventCreate(&start);
    // cudaEventCreate(&stop);
    // cudaEventRecord( start, 0);

    // assembles accel
#ifdef USE_CUDA
    if (run_cuda){
      assemble_boundary_accel_on_device<<<grid,threads,0,mp->compute_stream>>>(d_accel, d_send_buffer,
                                                                               mp->num_interfaces_ext_mesh,
                                                                               mp->max_nibool_interfaces_ext_mesh,
                                                                               mp->d_nibool_interfaces_ext_mesh,
                                                                               mp->d_ibool_interfaces_ext_mesh);
    }
#endif
#ifdef USE_HIP
    if (run_hip){
      hipLaunchKernelGGL(assemble_boundary_accel_on_device, dim3(grid), dim3(threads), 0, mp->compute_stream,
                                                            d_accel, d_send_buffer,
                                                            mp->num_interfaces_ext_mesh,
                                                            mp->max_nibool_interfaces_ext_mesh,
                                                            mp->d_nibool_interfaces_ext_mesh,
                                                            mp->d_ibool_interfaces_ext_mesh);
    }
#endif

    // cudaEventRecord( stop, 0);
    // cudaEventSynchronize( stop );
    // cudaEventElapsedTime( &time, start, stop );
    // cudaEventDestroy( start );
    // cudaEventDestroy( stop );
    // printf("Boundary Assemble Kernel Execution Time: %f ms\n",time);
  }

  //double end_time = get_time_val();
  //printf("Elapsed time: %e\n",end_time-start_time);

  GPU_ERROR_CHECKING("transfer_asmbl_accel_to_device");
}

/* ----------------------------------------------------------------------------------------------- */



/* ----------------------------------------------------------------------------------------------- */

// FORWARD_OR_ADJOINT == 1 for accel, and == 3 for b_accel
// This sync function is for FAULT_SOLVER
extern EXTERN_LANG
void FC_FUNC_(transfer_sync_accel_to_device,
              TRANSFER_ASMBL_ACCEL_TO_DEVICE)(long* Mesh_pointer,
                                              realw* buffer_recv_vector_ext_mesh,
                                              const int* num_interfaces_ext_mesh,
                                              const int* max_nibool_interfaces_ext_mesh,
                                              const int* nibool_interfaces_ext_mesh,
                                              const int* ibool_interfaces_ext_mesh,
                                              const int* FORWARD_OR_ADJOINT) {
TRACE("transfer_sync_accel_to_device");

  Mesh* mp = (Mesh*)(*Mesh_pointer); //get mesh pointer out of fortran integer container

  if (mp->size_mpi_buffer > 0){

    //daniel: todo - check if this copy is only needed for adjoint simulation, otherwise it is called asynchronously?
    if (*FORWARD_OR_ADJOINT == 1){
      // Wait until previous copy stream finishes. We assemble while other compute kernels execute.
      gpuStreamSynchronize(mp->copy_stream);
    }
    else if (*FORWARD_OR_ADJOINT == 3){
      // explicitly synchronizes
      // (cudaMemcpy implicitly synchronizes all other cuda operations)
      gpuSynchronize();

      gpuMemcpy_todevice_realw(mp->d_b_send_accel_buffer, buffer_recv_vector_ext_mesh,mp->size_mpi_buffer);
    }

    int blocksize = BLOCKSIZE_TRANSFER;
    int size_padded = ((int)ceil(((double)mp->max_nibool_interfaces_ext_mesh)/((double)blocksize)))*blocksize;

    int num_blocks_x, num_blocks_y;
    get_blocks_xy(size_padded/blocksize,&num_blocks_x,&num_blocks_y);

    dim3 grid(num_blocks_x,num_blocks_y);
    dim3 threads(blocksize,1,1);

    //double start_time = get_time_val();
    // cudaEvent_t start, stop;
    // realw time;
    // cudaEventCreate(&start);
    // cudaEventCreate(&stop);
    // cudaEventRecord( start, 0);

    // selects arrays
    realw* d_accel = NULL;
    realw* d_send_buffer = NULL;
    if (*FORWARD_OR_ADJOINT == 1) {
      // forward wavefield
      d_accel = mp->d_accel;
      d_send_buffer = mp->d_send_accel_buffer;
    } else if (*FORWARD_OR_ADJOINT == 3) {
      // backward/reconstructed wavefield
      d_accel = mp->d_b_accel;
      d_send_buffer = mp->d_b_send_accel_buffer;
    }

    //assembles accel
#ifdef USE_CUDA
    if (run_cuda){
      synchronize_boundary_accel_on_device<<<grid,threads,0,mp->compute_stream>>>(d_accel, d_send_buffer,
                                                                                  mp->num_interfaces_ext_mesh,
                                                                                  mp->max_nibool_interfaces_ext_mesh,
                                                                                  mp->d_nibool_interfaces_ext_mesh,
                                                                                  mp->d_ibool_interfaces_ext_mesh);
    }
#endif
#ifdef USE_HIP
    if (run_hip){
      hipLaunchKernelGGL(synchronize_boundary_accel_on_device, dim3(grid), dim3(threads), 0, mp->compute_stream,
                                                               d_accel, d_send_buffer,
                                                               mp->num_interfaces_ext_mesh,
                                                               mp->max_nibool_interfaces_ext_mesh,
                                                               mp->d_nibool_interfaces_ext_mesh,
                                                               mp->d_ibool_interfaces_ext_mesh);
    }
#endif

    // cudaEventRecord( stop, 0);
    // cudaEventSynchronize( stop );
    // cudaEventElapsedTime( &time, start, stop );
    // cudaEventDestroy( start );
    // cudaEventDestroy( stop );
    // printf("Boundary Assemble Kernel Execution Time: %f ms\n",time);
  }

  //double end_time = get_time_val();
  //printf("Elapsed time: %e\n",end_time-start_time);

  GPU_ERROR_CHECKING("transfer_sync_accel_to_device");
}

/* ----------------------------------------------------------------------------------------------- */

//daniel: not used ...
//
//extern EXTERN_LANG
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
//  //double start_time = get_time_val();
//  dim3 grid(num_blocks_x,num_blocks_y);
//  dim3 threads(blocksize,1,1);
//
//  // ***************************************************************************
//  // Wait until previous copy stream finishes. We assemble while other compute kernels execute.
//  cudaStreamSynchronize(mp->copy_stream);
//
//  // Assembling on the copy_stream breaks the solution and it "blows up"
//  if (*FORWARD_OR_ADJOINT == 1) { //assemble forward accel
//    assemble_boundary_accel_on_device<<<grid,threads,0,mp->compute_stream>>>(mp->d_accel, mp->d_send_accel_buffer,
//                                                                             mp->num_interfaces_ext_mesh,
//                                                                             mp->max_nibool_interfaces_ext_mesh,
//                                                                             mp->d_nibool_interfaces_ext_mesh,
//                                                                             mp->d_ibool_interfaces_ext_mesh);
//  }
//  else if (*FORWARD_OR_ADJOINT == 3) { //assemble adjoint accel
//    assemble_boundary_accel_on_device<<<grid,threads,0,mp->copy_stream>>>(mp->d_b_accel, mp->d_send_accel_buffer,
//                                                        mp->num_interfaces_ext_mesh,
//                                                        mp->max_nibool_interfaces_ext_mesh,
//                                                        mp->d_nibool_interfaces_ext_mesh,
//                                                        mp->d_ibool_interfaces_ext_mesh);
//  }
//
//GPU_ERROR_CHECKING("assemble_accel_on_device");
//}

/* ----------------------------------------------------------------------------------------------- */

extern EXTERN_LANG
void FC_FUNC_(sync_copy_from_device,
              SYNC_copy_FROM_DEVICE)(long* Mesh_pointer,
                                     int* iphase,
                                     realw* send_buffer) {

  TRACE("sync_copy_from_device");

  Mesh* mp = (Mesh*)(*Mesh_pointer); // get Mesh from fortran integer wrapper

  // Wait until async-memcpy of outer elements is finished and start MPI.
  if (*iphase != 2){ exit_on_gpu_error("sync_copy_from_device must be called for iphase == 2"); }

  if (mp->size_mpi_buffer > 0){
    // waits for asynchronous copy to finish
    gpuStreamSynchronize(mp->copy_stream);

    // There have been problems using the pinned-memory with MPI, so
    // we copy the buffer into a non-pinned region.
    memcpy(send_buffer,mp->h_send_accel_buffer,mp->size_mpi_buffer*sizeof(float));
  }
  // memory copy is now finished, so non-blocking MPI send can proceed
}

