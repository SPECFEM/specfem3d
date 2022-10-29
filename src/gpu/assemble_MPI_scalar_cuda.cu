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

extern EXTERN_LANG
void FC_FUNC_(transfer_boun_pot_from_device,
              TRANSFER_BOUN_POT_FROM_DEVICE)(long* Mesh_pointer,
                                             field* send_potential_dot_dot_buffer,
                                             const int* FORWARD_OR_ADJOINT){

  TRACE("transfer_boun_pot_from_device");
  //get mesh pointer out of fortran integer container
  Mesh* mp = (Mesh*)(*Mesh_pointer);

  // checks if anything to do
  if (mp->size_mpi_buffer_potential > 0){

    int blocksize = BLOCKSIZE_TRANSFER;
    int size_padded = ((int)ceil(((double)(mp->max_nibool_interfaces_ext_mesh))/((double)blocksize)))*blocksize;

    int num_blocks_x, num_blocks_y;
    get_blocks_xy(size_padded/blocksize,&num_blocks_x,&num_blocks_y);

    dim3 grid(num_blocks_x,num_blocks_y);
    dim3 threads(blocksize,1,1);

    // selects arrays
    field* d_potential_dot_dot = NULL;
    field* d_send_buffer = NULL;
    if (*FORWARD_OR_ADJOINT == 1) {
      // forward wavefield
      d_potential_dot_dot = mp->d_potential_dot_dot_acoustic;
      d_send_buffer = mp->d_send_potential_dot_dot_buffer;
    } else if (*FORWARD_OR_ADJOINT == 3) {
      // backward/reconstructed wavefield
      d_potential_dot_dot = mp->d_b_potential_dot_dot_acoustic;
      d_send_buffer = mp->d_b_send_potential_dot_dot_buffer;
    }

#ifdef USE_CUDA
    if (run_cuda){
      // fills mpi boundary buffer
      prepare_boundary_potential_on_device<<<grid,threads,0,mp->compute_stream>>>(d_potential_dot_dot,
                                                                                  d_send_buffer,
                                                                                  mp->num_interfaces_ext_mesh,
                                                                                  mp->max_nibool_interfaces_ext_mesh,
                                                                                  mp->d_nibool_interfaces_ext_mesh,
                                                                                  mp->d_ibool_interfaces_ext_mesh);
    }
#endif
#ifdef USE_HIP
    if (run_hip){
      // fills mpi boundary buffer
      hipLaunchKernelGGL(prepare_boundary_potential_on_device, dim3(grid), dim3(threads), 0, mp->compute_stream,
                                                               d_potential_dot_dot,
                                                               d_send_buffer,
                                                               mp->num_interfaces_ext_mesh,
                                                               mp->max_nibool_interfaces_ext_mesh,
                                                               mp->d_nibool_interfaces_ext_mesh,
                                                               mp->d_ibool_interfaces_ext_mesh);
    }
#endif

    //GPU_ERROR_CHECKING("after prepare_boundary_potential_on_device");

    // synchronizes
    //gpuSynchronize();
    // explicitly waits until previous compute stream finishes
    // (cudaMemcpy implicitly synchronizes all other cuda operations)
    gpuStreamSynchronize(mp->compute_stream);

    // copies buffer to CPU
    gpuMemcpy_tohost_field(send_potential_dot_dot_buffer,d_send_buffer,mp->size_mpi_buffer_potential);
  }

  // finish timing of kernel+memcpy
  // cudaEventRecord( stop, 0);
  // cudaEventSynchronize( stop );
  // cudaEventElapsedTime( &time, start, stop );
  // cudaEventDestroy( start );
  // cudaEventDestroy( stop );
  // printf("boundary xfer d->h Time: %f ms\n",time);
  GPU_ERROR_CHECKING("transfer_boun_pot_from_device");
}

/* ----------------------------------------------------------------------------------------------- */

// Assembly

/* ----------------------------------------------------------------------------------------------- */

extern EXTERN_LANG
void FC_FUNC_(transfer_asmbl_pot_to_device,
              TRANSFER_ASMBL_POT_TO_DEVICE)(long* Mesh_pointer,
                                            field* buffer_recv_scalar_ext_mesh,
                                            const int* FORWARD_OR_ADJOINT) {

  TRACE("transfer_asmbl_pot_to_device");
  //get mesh pointer out of fortran integer container
  Mesh* mp = (Mesh*)(*Mesh_pointer);

  // Cuda timing
  //gpu_event start, stop;
  //start_timing_gpu(&start,&stop);

  // checks if anything to do
  if (mp->size_mpi_buffer_potential > 0){

    // assembles on GPU
    int blocksize = BLOCKSIZE_TRANSFER;
    int size_padded = ((int)ceil(((double)mp->max_nibool_interfaces_ext_mesh)/((double)blocksize)))*blocksize;

    int num_blocks_x, num_blocks_y;
    get_blocks_xy(size_padded/blocksize,&num_blocks_x,&num_blocks_y);

    dim3 grid(num_blocks_x,num_blocks_y);
    dim3 threads(blocksize,1,1);

    // selects arrays
    field* d_potential_dot_dot = NULL;
    field* d_send_buffer = NULL;
    if (*FORWARD_OR_ADJOINT == 1) {
      // forward wavefield
      d_potential_dot_dot = mp->d_potential_dot_dot_acoustic;
      d_send_buffer = mp->d_send_potential_dot_dot_buffer;
    } else if (*FORWARD_OR_ADJOINT == 3) {
      // backward/reconstructed wavefield
      d_potential_dot_dot = mp->d_b_potential_dot_dot_acoustic;
      d_send_buffer = mp->d_b_send_potential_dot_dot_buffer;
    }

    // synchronizes
    gpuSynchronize();

    // copies buffer onto GPU
    gpuMemcpy_todevice_field(d_send_buffer, buffer_recv_scalar_ext_mesh,mp->size_mpi_buffer_potential);

    // assembles field
#ifdef USE_CUDA
    if (run_cuda){
      assemble_boundary_potential_on_device<<<grid,threads,0,mp->compute_stream>>>(d_potential_dot_dot,
                                                                                   d_send_buffer,
                                                                                   mp->num_interfaces_ext_mesh,
                                                                                   mp->max_nibool_interfaces_ext_mesh,
                                                                                   mp->d_nibool_interfaces_ext_mesh,
                                                                                   mp->d_ibool_interfaces_ext_mesh);
    }
#endif
#ifdef USE_HIP
    if (run_hip){
      hipLaunchKernelGGL(assemble_boundary_potential_on_device, dim3(grid), dim3(threads), 0, mp->compute_stream,
                                                                d_potential_dot_dot,
                                                                d_send_buffer,
                                                                mp->num_interfaces_ext_mesh,
                                                                mp->max_nibool_interfaces_ext_mesh,
                                                                mp->d_nibool_interfaces_ext_mesh,
                                                                mp->d_ibool_interfaces_ext_mesh);
    }
#endif
  }
  // kernel timing
  //stop_timing_gpu(&start,&stop,"assemble_boundary_potential_on_device");

  GPU_ERROR_CHECKING("transfer_asmbl_pot_to_device");
}

