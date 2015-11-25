/*
 !=====================================================================
 !
 !               S p e c f e m 3 D  V e r s i o n  3 . 0
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

#include "mesh_constants_cuda.h"

/* ----------------------------------------------------------------------------------------------- */

// ASSEMBLY - mpi data transfer between CPU-GPU

/* ----------------------------------------------------------------------------------------------- */

// prepares a device array with with all inter-element edge-nodes -- this
// is followed by a memcpy and MPI operations

__global__ void prepare_boundary_accel_on_device(realw* d_accel, realw* d_send_accel_buffer,
                                                 const int num_interfaces_ext_mesh,
                                                 const int max_nibool_interfaces_ext_mesh,
                                                 const int* d_nibool_interfaces_ext_mesh,
                                                 const int* d_ibool_interfaces_ext_mesh) {

  int id = threadIdx.x + blockIdx.x*blockDim.x + blockIdx.y*gridDim.x*blockDim.x;
  int ientry,iglob;

  for( int iinterface=0; iinterface < num_interfaces_ext_mesh; iinterface++) {
    if (id < d_nibool_interfaces_ext_mesh[iinterface]) {

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
                                               const int* FORWARD_OR_ADJOINT){
TRACE("\ttransfer_boun_accel_from_device");

  Mesh* mp = (Mesh*)(*Mesh_pointer); //get mesh pointer out of fortran integer container

  // checks if anything to do
  if (mp->size_mpi_buffer > 0){

    int blocksize = BLOCKSIZE_TRANSFER;
    int size_padded = ((int)ceil(((double)mp->max_nibool_interfaces_ext_mesh)/((double)blocksize)))*blocksize;

    int num_blocks_x, num_blocks_y;
    get_blocks_xy(size_padded/blocksize,&num_blocks_x,&num_blocks_y);

    dim3 grid(num_blocks_x,num_blocks_y);
    dim3 threads(blocksize,1,1);

    // Cuda timing
    //cudaEvent_t start, stop;
    //start_timing_cuda(&start,&stop);

    if (*FORWARD_OR_ADJOINT == 1) {
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
    else if (*FORWARD_OR_ADJOINT == 3) {
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

    // Cuda timing
    // finish timing of kernel+memcpy
    //stop_timing_cuda(&start,&stop,"prepare_boundary_accel_on_device");
  }

#ifdef ENABLE_VERY_SLOW_ERROR_CHECKING
  exit_on_cuda_error("transfer_boun_accel_from_device");
#endif
}

/* ----------------------------------------------------------------------------------------------- */

extern "C"
void FC_FUNC_(transfer_boundary_from_device_a,
              TRANSFER_BOUNDARY_FROM_DEVICE_A)(long* Mesh_pointer,
                                               const int* nspec_outer_elastic) {

// asynchronous transfer from device to host

  TRACE("\ttransfer_boundary_from_device_a");

  Mesh* mp = (Mesh*)(*Mesh_pointer); // get Mesh from fortran integer wrapper

  if (mp->size_mpi_buffer > 0){

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
                                             const int* num_interfaces_ext_mesh,
                                             const int* max_nibool_interfaces_ext_mesh) {

// asynchronous transfer from host to device

  TRACE("transfer_boundary_to_device_a");

  Mesh* mp = (Mesh*)(*Mesh_pointer); //get mesh pointer out of fortran integer container

  if (mp->size_mpi_buffer > 0){
    // copy on host memory
    memcpy(mp->h_recv_accel_buffer,buffer_recv_vector_ext_mesh,mp->size_mpi_buffer*sizeof(realw));

    // asynchronous copy to GPU using copy_stream
    cudaMemcpyAsync(mp->d_send_accel_buffer, mp->h_recv_accel_buffer,
                    mp->size_mpi_buffer*sizeof(realw),cudaMemcpyHostToDevice,mp->copy_stream);
  }
}


/* ----------------------------------------------------------------------------------------------- */

// Assembly

/* ----------------------------------------------------------------------------------------------- */

__global__ void assemble_boundary_accel_on_device(realw* d_accel, realw* d_send_accel_buffer,
                                                  const int num_interfaces_ext_mesh,
                                                  const int max_nibool_interfaces_ext_mesh,
                                                  const int* d_nibool_interfaces_ext_mesh,
                                                  const int* d_ibool_interfaces_ext_mesh) {

  //int bx = blockIdx.y*gridDim.x+blockIdx.x;
  //int tx = threadIdx.x;
  int id = threadIdx.x + blockIdx.x*blockDim.x + blockIdx.y*gridDim.x*blockDim.x;

  int ientry,iglob;

  for( int iinterface=0; iinterface < num_interfaces_ext_mesh; iinterface++) {
    if (id < d_nibool_interfaces_ext_mesh[iinterface]) {

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



__global__ void synchronize_boundary_accel_on_device(realw* d_accel, realw* d_send_accel_buffer,
                                                  const int num_interfaces_ext_mesh,
                                                  const int max_nibool_interfaces_ext_mesh,
                                                  const int* d_nibool_interfaces_ext_mesh,
                                                  const int* d_ibool_interfaces_ext_mesh) {

  //int bx = blockIdx.y*gridDim.x+blockIdx.x;
  //int tx = threadIdx.x;
  int id = threadIdx.x + blockIdx.x*blockDim.x + blockIdx.y*gridDim.x*blockDim.x;
  // printf("inside the synchronization!\n");

  int ientry,iglob;

  for( int iinterface=0; iinterface < num_interfaces_ext_mesh; iinterface++) {
    if (id < d_nibool_interfaces_ext_mesh[iinterface]) {

      // entry in interface array
      ientry = id + max_nibool_interfaces_ext_mesh*iinterface;
      // global index in wavefield
      iglob = d_ibool_interfaces_ext_mesh[ientry] - 1;

      // for testing atomic operations against not atomic operations (0.1ms vs. 0.04 ms)
      // d_accel[3*(iglob)] += d_send_accel_buffer[3*(ientry)];
      // d_accel[3*(iglob)+1] += d_send_accel_buffer[3*(ientry)+1];
      // d_accel[3*(iglob)+2] += d_send_accel_buffer[3*(ientry)+2];
        //    printf("\n1:%f,2:%f,max:%f\n",d_accel[3*iglob],d_send_accel_buffer[3*ientry],fmax(d_accel[3*iglob],d_send_accel_buffer[3*ientry]));


      d_accel[3*iglob]  = fmax(d_accel[3*iglob],d_send_accel_buffer[3*ientry]);
      d_accel[3*iglob+1]= fmax(d_accel[3*iglob + 1],d_send_accel_buffer[3*ientry + 1]);
      d_accel[3*iglob+2]= fmax(d_accel[3*iglob + 2],d_send_accel_buffer[3*ientry + 2]);
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
                                              const int* num_interfaces_ext_mesh,
                                              const int* max_nibool_interfaces_ext_mesh,
                                              const int* nibool_interfaces_ext_mesh,
                                              const int* ibool_interfaces_ext_mesh,
                                              const int* FORWARD_OR_ADJOINT) {
TRACE("\ttransfer_asmbl_accel_to_device");

  Mesh* mp = (Mesh*)(*Mesh_pointer); //get mesh pointer out of fortran integer container

  if (mp->size_mpi_buffer > 0){

    //daniel: todo - check if this copy is only needed for adjoint simulation, otherwise it is called asynchronously?
    if (*FORWARD_OR_ADJOINT == 1){
      // Wait until previous copy stream finishes. We assemble while other compute kernels execute.
      cudaStreamSynchronize(mp->copy_stream);
    }
    else if (*FORWARD_OR_ADJOINT == 3){
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
    // cudaEventRecord( start, 0);

    if (*FORWARD_OR_ADJOINT == 1) {
      //assemble forward accel
      assemble_boundary_accel_on_device<<<grid,threads,0,mp->compute_stream>>>(mp->d_accel, mp->d_send_accel_buffer,
                                                                               mp->num_interfaces_ext_mesh,
                                                                               mp->max_nibool_interfaces_ext_mesh,
                                                                               mp->d_nibool_interfaces_ext_mesh,
                                                                               mp->d_ibool_interfaces_ext_mesh);
    }
    else if (*FORWARD_OR_ADJOINT == 3) {
      //assemble adjoint accel
      assemble_boundary_accel_on_device<<<grid,threads,0,mp->compute_stream>>>(mp->d_b_accel, mp->d_b_send_accel_buffer,
                                                                               mp->num_interfaces_ext_mesh,
                                                                               mp->max_nibool_interfaces_ext_mesh,
                                                                               mp->d_nibool_interfaces_ext_mesh,
                                                                               mp->d_ibool_interfaces_ext_mesh);
    }

    // cudaEventRecord( stop, 0);
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



/* ----------------------------------------------------------------------------------------------- */

// FORWARD_OR_ADJOINT == 1 for accel, and == 3 for b_accel
// This sync function is for FAULT_SOLVER
extern "C"
void FC_FUNC_(transfer_sync_accel_to_device,
              TRANSFER_ASMBL_ACCEL_TO_DEVICE)(long* Mesh_pointer, realw* accel,
                                              realw* buffer_recv_vector_ext_mesh,
                                              const int* num_interfaces_ext_mesh,
                                              const int* max_nibool_interfaces_ext_mesh,
                                              const int* nibool_interfaces_ext_mesh,
                                              const int* ibool_interfaces_ext_mesh,
                                              const int* FORWARD_OR_ADJOINT) {
TRACE("\ttransfer_sync_accel_to_device");

  Mesh* mp = (Mesh*)(*Mesh_pointer); //get mesh pointer out of fortran integer container

  if (mp->size_mpi_buffer > 0){

    //daniel: todo - check if this copy is only needed for adjoint simulation, otherwise it is called asynchronously?
    if (*FORWARD_OR_ADJOINT == 1){
      // Wait until previous copy stream finishes. We assemble while other compute kernels execute.
      cudaStreamSynchronize(mp->copy_stream);
    }
    else if (*FORWARD_OR_ADJOINT == 3){
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
    // cudaEventRecord( start, 0);

    if (*FORWARD_OR_ADJOINT == 1) {
      //assemble forward accel
      synchronize_boundary_accel_on_device<<<grid,threads,0,mp->compute_stream>>>(mp->d_accel, mp->d_send_accel_buffer,
                                                                               mp->num_interfaces_ext_mesh,
                                                                               mp->max_nibool_interfaces_ext_mesh,
                                                                               mp->d_nibool_interfaces_ext_mesh,
                                                                               mp->d_ibool_interfaces_ext_mesh);
    }
    else if (*FORWARD_OR_ADJOINT == 3) {
      //assemble adjoint accel
      synchronize_boundary_accel_on_device<<<grid,threads,0,mp->compute_stream>>>(mp->d_b_accel, mp->d_b_send_accel_buffer,
                                                                               mp->num_interfaces_ext_mesh,
                                                                               mp->max_nibool_interfaces_ext_mesh,
                                                                               mp->d_nibool_interfaces_ext_mesh,
                                                                               mp->d_ibool_interfaces_ext_mesh);
    }

    // cudaEventRecord( stop, 0);
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
//#ifdef ENABLE_VERY_SLOW_ERROR_CHECKING
//  exit_on_cuda_error("assemble_accel_on_device");
//#endif
//}

/* ----------------------------------------------------------------------------------------------- */

extern "C"
void FC_FUNC_(sync_copy_from_device,
              SYNC_copy_FROM_DEVICE)(long* Mesh_pointer,
                                     int* iphase,
                                     realw* send_buffer) {

  TRACE("sync_copy_from_device");

  Mesh* mp = (Mesh*)(*Mesh_pointer); // get Mesh from fortran integer wrapper

  // Wait until async-memcpy of outer elements is finished and start MPI.
  if (*iphase != 2){ exit_on_cuda_error("sync_copy_from_device must be called for iphase == 2"); }

  if (mp->size_mpi_buffer > 0){
    // waits for asynchronous copy to finish
    cudaStreamSynchronize(mp->copy_stream);

    // There have been problems using the pinned-memory with MPI, so
    // we copy the buffer into a non-pinned region.
    memcpy(send_buffer,mp->h_send_accel_buffer,mp->size_mpi_buffer*sizeof(float));
  }
  // memory copy is now finished, so non-blocking MPI send can proceed
}

