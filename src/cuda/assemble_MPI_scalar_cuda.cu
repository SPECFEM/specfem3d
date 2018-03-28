/*
 !=====================================================================
 !
 !               S p e c f e m 3 D  V e r s i o n  3 . 0
 !               ---------------------------------------
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

#include "mesh_constants_cuda.h"

/* ----------------------------------------------------------------------------------------------- */

// ASSEMBLY - mpi data transfer between CPU-GPU

/* ----------------------------------------------------------------------------------------------- */

// prepares a device array with with all inter-element edge-nodes -- this
// is followed by a memcpy and MPI operations
__global__ void prepare_boundary_potential_on_device(field* d_potential_dot_dot_acoustic,
                                                     field* d_send_potential_dot_dot_buffer,
                                                     const int num_interfaces_ext_mesh,
                                                     const int max_nibool_interfaces_ext_mesh,
                                                     const int* d_nibool_interfaces_ext_mesh,
                                                     const int* d_ibool_interfaces_ext_mesh) {

  int id = threadIdx.x + blockIdx.x*blockDim.x + blockIdx.y*gridDim.x*blockDim.x;
  int ientry,iglob;

  for(int iinterface=0; iinterface < num_interfaces_ext_mesh; iinterface++) {
    if (id<d_nibool_interfaces_ext_mesh[iinterface]) {

      // entry in interface array
      ientry = id + max_nibool_interfaces_ext_mesh*iinterface;
      // global index in wavefield
      iglob = d_ibool_interfaces_ext_mesh[ientry] - 1;

      d_send_potential_dot_dot_buffer[ientry] = d_potential_dot_dot_acoustic[iglob];
    }
  }

}


/* ----------------------------------------------------------------------------------------------- */

// prepares and transfers the inter-element edge-nodes to the host to be MPI'd
extern "C"
void FC_FUNC_(transfer_boun_pot_from_device,
              TRANSFER_BOUN_POT_FROM_DEVICE)(long* Mesh_pointer,
                                             field* potential_dot_dot_acoustic,
                                             field* send_potential_dot_dot_buffer,
                                             const int* FORWARD_OR_ADJOINT){

TRACE("transfer_boun_pot_from_device");

  Mesh* mp = (Mesh*)(*Mesh_pointer); //get mesh pointer out of fortran integer container

  // checks if anything to do
  if (mp->size_mpi_buffer_potential > 0){

    int blocksize = BLOCKSIZE_TRANSFER;
    int size_padded = ((int)ceil(((double)(mp->max_nibool_interfaces_ext_mesh))/((double)blocksize)))*blocksize;

    int num_blocks_x, num_blocks_y;
    get_blocks_xy(size_padded/blocksize,&num_blocks_x,&num_blocks_y);

    dim3 grid(num_blocks_x,num_blocks_y);
    dim3 threads(blocksize,1,1);

    if (*FORWARD_OR_ADJOINT == 1) {
      prepare_boundary_potential_on_device<<<grid,threads,0,mp->compute_stream>>>(mp->d_potential_dot_dot_acoustic,
                                                                                   mp->d_send_potential_dot_dot_buffer,
                                                                                   mp->num_interfaces_ext_mesh,
                                                                                   mp->max_nibool_interfaces_ext_mesh,
                                                                                   mp->d_nibool_interfaces_ext_mesh,
                                                                                   mp->d_ibool_interfaces_ext_mesh);

      // synchronizes
      //synchronize_cuda();
      // explicitly waits until previous compute stream finishes
      // (cudaMemcpy implicitly synchronizes all other cuda operations)
      cudaStreamSynchronize(mp->compute_stream);

      print_CUDA_error_if_any(cudaMemcpy(send_potential_dot_dot_buffer,mp->d_send_potential_dot_dot_buffer,
                                         mp->size_mpi_buffer_potential*sizeof(field),cudaMemcpyDeviceToHost),98000);
    }
    else if (*FORWARD_OR_ADJOINT == 3) {
      // backward/reconstructed wavefield buffer
      prepare_boundary_potential_on_device<<<grid,threads,0,mp->compute_stream>>>(mp->d_b_potential_dot_dot_acoustic,
                                                                                   mp->d_b_send_potential_dot_dot_buffer,
                                                                                   mp->num_interfaces_ext_mesh,
                                                                                   mp->max_nibool_interfaces_ext_mesh,
                                                                                   mp->d_nibool_interfaces_ext_mesh,
                                                                                   mp->d_ibool_interfaces_ext_mesh);

      // synchronizes
      //synchronize_cuda();
      // explicitly waits until previous compute stream finishes
      // (cudaMemcpy implicitly synchronizes all other cuda operations)
      cudaStreamSynchronize(mp->compute_stream);

      print_CUDA_error_if_any(cudaMemcpy(send_potential_dot_dot_buffer,mp->d_b_send_potential_dot_dot_buffer,
                                         mp->size_mpi_buffer_potential*sizeof(field),cudaMemcpyDeviceToHost),98000);
    }
  }

#ifdef ENABLE_VERY_SLOW_ERROR_CHECKING
  exit_on_cuda_error("after prepare_boundary_potential_on_device");
#endif


  // finish timing of kernel+memcpy
  // cudaEventRecord( stop, 0);
  // cudaEventSynchronize( stop );
  // cudaEventElapsedTime( &time, start, stop );
  // cudaEventDestroy( start );
  // cudaEventDestroy( stop );
  // printf("boundary xfer d->h Time: %f ms\n",time);
#ifdef ENABLE_VERY_SLOW_ERROR_CHECKING
  exit_on_cuda_error("transfer_boun_pot_from_device");
#endif
}

/* ----------------------------------------------------------------------------------------------- */

// Assembly

/* ----------------------------------------------------------------------------------------------- */


__global__ void assemble_boundary_potential_on_device(field* d_potential_dot_dot_acoustic,
                                                      field* d_send_potential_dot_dot_buffer,
                                                      const int num_interfaces_ext_mesh,
                                                      const int max_nibool_interfaces_ext_mesh,
                                                      const int* d_nibool_interfaces_ext_mesh,
                                                      const int* d_ibool_interfaces_ext_mesh) {

  int id = threadIdx.x + blockIdx.x*blockDim.x + blockIdx.y*gridDim.x*blockDim.x;
  int ientry,iglob;

  for( int iinterface=0; iinterface < num_interfaces_ext_mesh; iinterface++) {
    if (id<d_nibool_interfaces_ext_mesh[iinterface]) {

      // entry in interface array
      ientry = id + max_nibool_interfaces_ext_mesh*iinterface;
      // global index in wavefield
      iglob = d_ibool_interfaces_ext_mesh[ientry] - 1;

      // for testing atomic operations against not atomic operations (0.1ms vs. 0.04 ms)
      // d_potential_dot_dot_acoustic[3*(d_ibool_interfaces_ext_mesh[id+max_nibool_interfaces_ext_mesh*iinterface]-1)] +=
      // d_send_potential_dot_dot_buffer[3*(id + max_nibool_interfaces_ext_mesh*iinterface)];

      atomicAdd(&d_potential_dot_dot_acoustic[iglob],d_send_potential_dot_dot_buffer[ientry]);
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
void FC_FUNC_(transfer_asmbl_pot_to_device,
              TRANSFER_ASMBL_POT_TO_DEVICE)(long* Mesh_pointer,
                                            field* potential_dot_dot_acoustic,
                                            field* buffer_recv_scalar_ext_mesh,
                                            const int* FORWARD_OR_ADJOINT) {

TRACE("transfer_asmbl_pot_to_device");

  Mesh* mp = (Mesh*)(*Mesh_pointer); //get mesh pointer out of fortran integer container

  // Cuda timing
  //cudaEvent_t start, stop;
  //start_timing_cuda(&start,&stop);

  // checks if anything to do
  if (mp->size_mpi_buffer_potential > 0){

    // assembles on GPU
    int blocksize = BLOCKSIZE_TRANSFER;
    int size_padded = ((int)ceil(((double)mp->max_nibool_interfaces_ext_mesh)/((double)blocksize)))*blocksize;

    int num_blocks_x, num_blocks_y;
    get_blocks_xy(size_padded/blocksize,&num_blocks_x,&num_blocks_y);

    dim3 grid(num_blocks_x,num_blocks_y);
    dim3 threads(blocksize,1,1);

    // synchronizes
    synchronize_cuda();

    if (*FORWARD_OR_ADJOINT == 1) {
      // copies buffer onto GPU
      print_CUDA_error_if_any(cudaMemcpy(mp->d_send_potential_dot_dot_buffer, buffer_recv_scalar_ext_mesh,
                                         mp->size_mpi_buffer_potential*sizeof(field), cudaMemcpyHostToDevice),98010);

      //assemble forward field
      assemble_boundary_potential_on_device<<<grid,threads,0,mp->compute_stream>>>(mp->d_potential_dot_dot_acoustic,
                                                                                    mp->d_send_potential_dot_dot_buffer,
                                                                                    mp->num_interfaces_ext_mesh,
                                                                                    mp->max_nibool_interfaces_ext_mesh,
                                                                                    mp->d_nibool_interfaces_ext_mesh,
                                                                                    mp->d_ibool_interfaces_ext_mesh);
    }
    else if (*FORWARD_OR_ADJOINT == 3) {
      // copies buffer onto GPU
      print_CUDA_error_if_any(cudaMemcpy(mp->d_b_send_potential_dot_dot_buffer, buffer_recv_scalar_ext_mesh,
                                         mp->size_mpi_buffer_potential*sizeof(field), cudaMemcpyHostToDevice),98011);

      //assemble reconstructed/backward field
      assemble_boundary_potential_on_device<<<grid,threads,0,mp->compute_stream>>>(mp->d_b_potential_dot_dot_acoustic,
                                                                                    mp->d_b_send_potential_dot_dot_buffer,
                                                                                    mp->num_interfaces_ext_mesh,
                                                                                    mp->max_nibool_interfaces_ext_mesh,
                                                                                    mp->d_nibool_interfaces_ext_mesh,
                                                                                    mp->d_ibool_interfaces_ext_mesh);
    }
  }

  // Cuda timing
  //stop_timing_cuda(&start,&stop,"assemble_boundary_potential_on_device");

#ifdef ENABLE_VERY_SLOW_ERROR_CHECKING
  exit_on_cuda_error("transfer_asmbl_pot_to_device");
#endif
}

