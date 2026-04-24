/*
 !=====================================================================
 !
 !                         S p e c f e m 3 D
 !                         -----------------
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

#include "mesh_constants_gpu.h"

/* ----------------------------------------------------------------------------------------------- */

extern EXTERN_LANG
void FC_FUNC_(sync_copy_reduced_from_device,
              SYNC_copy_reduced_FROM_DEVICE)(long* Mesh_pointer,
                                             int* iphase,
                                             realw* send_buffer,
                                             int* num_interface_p_refine_boundary_f) {

  TRACE("sync_copy_reduced_from_device");

  Mesh* mp = (Mesh*)(*Mesh_pointer); // get Mesh from fortran integer wrapper
  int num_interface_p_refine_boundary = *num_interface_p_refine_boundary_f;

  // checks if anything to do
  if (num_interface_p_refine_boundary == 0) return;

  // Wait until async-memcpy of outer elements is finished and start MPI.
  if (*iphase != 2){ exit_on_error("sync_copy_reduced_from_device must be called for iphase == 2"); }

  if (mp->size_mpi_buffer > 0){
    // waits for asynchronous copy to finish
    gpuStreamSynchronize(mp->copy_stream);

    // There have been problems using the pinned-memory with MPI, so
    // we copy the buffer into a non-pinned region.
    memcpy(send_buffer, mp->h_send_accel_buffer, NDIM * num_interface_p_refine_boundary * sizeof(realw));
  }
  // memory copy is now finished, so non-blocking MPI send can proceed
}


/* ----------------------------------------------------------------------------------------------- */
// transfer routines
/* ----------------------------------------------------------------------------------------------- */

extern EXTERN_LANG
void FC_FUNC_(test_boundary_transfer_lts,
              TEST_BOUNDARY_TRANSFER_LTS)(long* Mesh_pointer,
                                          int* ilevel_f,
                                          int* max_num_interface_p_refine_ibool_f,
                                          realw* copy_buffer) {

  TRACE("test_boundary_transfer_lts");

  Mesh* mp = (Mesh*)(*Mesh_pointer); // get Mesh from fortran integer wrapper
  int ilevel = *ilevel_f;
  int max_num_interface_p_refine_ibool = *max_num_interface_p_refine_ibool_f;

  if (mp->size_mpi_buffer > 0){
    int blocksize = BLOCKSIZE_TRANSFER;

    int size_padded = ((int)ceil(((double)max_num_interface_p_refine_ibool)/((double)blocksize)))*blocksize;

    int num_blocks_x, num_blocks_y;
    get_blocks_xy(size_padded/blocksize,&num_blocks_x,&num_blocks_y);

    dim3 grid(num_blocks_x,num_blocks_y);
    dim3 threads(blocksize,1,1);

    gpuMemset_realw(mp->d_send_accel_buffer, NDIM * mp->max_nibool_interfaces_ext_mesh * mp->num_interfaces_ext_mesh, 0);

#ifdef USE_CUDA
    if (run_cuda){
      prepare_boundary_lts_accel_on_device<<<grid,threads,0,mp->compute_stream>>>(mp->d_accel,
                                                                                  mp->d_send_accel_buffer,
                                                                                  mp->num_interfaces_ext_mesh,
                                                                                  mp->d_lts_num_interface_p_refine_ibool,
                                                                                  mp->d_lts_interface_p_refine_ibool,
                                                                                  mp->max_nibool_interfaces_ext_mesh,
                                                                                  ilevel);
    }
#endif
#ifdef USE_HIP
    if (run_hip){
      hipLaunchKernelGGL(prepare_boundary_lts_accel_on_device, dim3(grid), dim3(threads), 0, mp->compute_stream,
                                                                                  mp->d_accel,
                                                                                  mp->d_send_accel_buffer,
                                                                                  mp->num_interfaces_ext_mesh,
                                                                                  mp->d_lts_num_interface_p_refine_ibool,
                                                                                  mp->d_lts_interface_p_refine_ibool,
                                                                                  mp->max_nibool_interfaces_ext_mesh,
                                                                                  ilevel);
    }
#endif
    GPU_ERROR_CHECKING("prepare_boundary_lts_accel_on_device");

    // wait until kernel is finished before starting async memcpy
    gpuSynchronize();

    // this vector is bigger than necessary, however building it to fit properly is quite difficult
    gpuMemcpy_tohost_realw(copy_buffer, mp->d_send_accel_buffer, NDIM * mp->max_nibool_interfaces_ext_mesh * mp->num_interfaces_ext_mesh);
  }
}

/* ----------------------------------------------------------------------------------------------- */

extern EXTERN_LANG
void FC_FUNC_(transfer_reduced_boundary_from_device_async_lts,
              TRANSFER_REDUCED_BOUNDARY_FROM_DEVICE_ASYNC_LTS)(long* Mesh_pointer,
                                                               int* ilevel_f,
                                                               int* num_interface_p_refine_boundary_f) {

// asynchronous transfer from device to host

  TRACE("transfer_reduced_boundary_from_device_async_lts");

  Mesh* mp = (Mesh*)(*Mesh_pointer); // get Mesh from fortran integer wrapper
  int ilevel = *ilevel_f;
  int num_interface_p_refine_boundary = *num_interface_p_refine_boundary_f;

  // checks if anything to do
  if (num_interface_p_refine_boundary == 0) return;

  if (mp->size_mpi_buffer > 0){
    int blocksize = BLOCKSIZE_TRANSFER;

    int size_padded = ((int)ceil(((double)num_interface_p_refine_boundary)/((double)blocksize)))*blocksize;

    int num_blocks_x, num_blocks_y;
    get_blocks_xy(size_padded/blocksize,&num_blocks_x,&num_blocks_y);

    dim3 grid(num_blocks_x,num_blocks_y);
    dim3 threads(blocksize,1,1);

#ifdef USE_CUDA
    if (run_cuda){
      prepare_reduced_boundary_lts_accel_on_device<<<grid,threads,0,mp->compute_stream>>>(mp->d_accel,
                                                                                          mp->d_send_accel_buffer,
                                                                                          num_interface_p_refine_boundary,
                                                                                          mp->d_lts_interface_p_refine_boundary,
                                                                                          mp->lts_max_nibool_interfaces_boundary,
                                                                                          ilevel);
    }
#endif
#ifdef USE_HIP
    if (run_hip){
      hipLaunchKernelGGL(prepare_reduced_boundary_lts_accel_on_device, dim3(grid), dim3(threads), 0, mp->compute_stream,
                                                                                          mp->d_accel,
                                                                                          mp->d_send_accel_buffer,
                                                                                          num_interface_p_refine_boundary,
                                                                                          mp->d_lts_interface_p_refine_boundary,
                                                                                          mp->lts_max_nibool_interfaces_boundary,
                                                                                          ilevel);
    }
#endif

    // wait until kernel is finished before starting async memcpy
    gpuSynchronize();

    // this vector is only as big as mpi-boundary for current p-level
    gpuMemcpyAsync_tohost_realw(mp->h_send_accel_buffer, mp->d_send_accel_buffer,
                                NDIM * num_interface_p_refine_boundary, mp->copy_stream);
  }
}

/* ----------------------------------------------------------------------------------------------- */

extern EXTERN_LANG
void FC_FUNC_(transfer_boundary_from_device_async_lts,
              TRANSFER_BOUNDARY_FROM_DEVICE_ASYNC_LTS)(long* Mesh_pointer,
                                                       int* ilevel_f,
                                                       int* max_num_interface_p_refine_ibool_f) {

// asynchronous transfer from device to host

  TRACE("transfer_boundary_from_device_async_lts");

  Mesh* mp = (Mesh*)(*Mesh_pointer); // get Mesh from fortran integer wrapper
  int ilevel = *ilevel_f;
  int max_num_interface_p_refine_ibool = *max_num_interface_p_refine_ibool_f;

  if (mp->size_mpi_buffer > 0){
    int blocksize = BLOCKSIZE_TRANSFER;

    int size_padded = ((int)ceil(((double)max_num_interface_p_refine_ibool)/((double)blocksize)))*blocksize;

    int num_blocks_x, num_blocks_y;
    get_blocks_xy(size_padded/blocksize,&num_blocks_x,&num_blocks_y);

    dim3 grid(num_blocks_x,num_blocks_y);
    dim3 threads(blocksize,1,1);

#ifdef USE_CUDA
    if (run_cuda){
      prepare_boundary_lts_accel_on_device<<<grid,threads,0,mp->compute_stream>>>(mp->d_accel,
                                                                                  mp->d_send_accel_buffer,
                                                                                  mp->num_interfaces_ext_mesh,
                                                                                  mp->d_lts_num_interface_p_refine_ibool,
                                                                                  mp->d_lts_interface_p_refine_ibool,
                                                                                  mp->max_nibool_interfaces_ext_mesh,
                                                                                  ilevel);
    }
#endif
#ifdef USE_HIP
    if (run_hip){
      hipLaunchKernelGGL(prepare_boundary_lts_accel_on_device, dim3(grid), dim3(threads), 0, mp->compute_stream,
                                                                                  mp->d_accel,
                                                                                  mp->d_send_accel_buffer,
                                                                                  mp->num_interfaces_ext_mesh,
                                                                                  mp->d_lts_num_interface_p_refine_ibool,
                                                                                  mp->d_lts_interface_p_refine_ibool,
                                                                                  mp->max_nibool_interfaces_ext_mesh,
                                                                                  ilevel);
    }
#endif

    // wait until kernel is finished before starting async memcpy
    gpuSynchronize();

    // this vector is bigger than necessary, however building it to fit properly is quite difficult
    gpuMemcpyAsync_tohost_realw(mp->h_send_accel_buffer, mp->d_send_accel_buffer,
                                NDIM * mp->max_nibool_interfaces_ext_mesh * mp->num_interfaces_ext_mesh, mp->copy_stream);
  }
}

/* ----------------------------------------------------------------------------------------------- */

extern EXTERN_LANG
void FC_FUNC_(transfer_reduced_boundary_to_device_async_lts,
              TRANSFER_REDUCED_BOUNDARY_TO_DEVICE_ASYNC_LTS)(long* Mesh_pointer,
                                                             realw* buffer_reduced_recv_vector,
                                                             int* num_interface_p_refine_boundary_f) {

// asynchronous transfer from host to device

  TRACE("transfer_reduced_boundary_to_device_async_lts");

  Mesh* mp = (Mesh*)(*Mesh_pointer); //get mesh pointer out of fortran integer container
  int num_interface_p_refine_boundary = *num_interface_p_refine_boundary_f;

  // checks if anything to do
  if (num_interface_p_refine_boundary == 0) return;

  if (mp->size_mpi_buffer > 0){
    // copy on host memory
    memcpy(mp->h_recv_accel_buffer,buffer_reduced_recv_vector,3*num_interface_p_refine_boundary*sizeof(realw));

    // asynchronous copy to GPU using copy_stream
    gpuMemcpyAsync_todevice_realw(mp->d_send_accel_buffer, mp->h_recv_accel_buffer, NDIM * num_interface_p_refine_boundary, mp->copy_stream);
  }
}

/* ----------------------------------------------------------------------------------------------- */
// assemble routines
/* ----------------------------------------------------------------------------------------------- */

extern EXTERN_LANG
void FC_FUNC_(assemble_mpi_device_lts,
              ASSEMBLE_MPI_DEVICE_LTS)(long* Mesh_pointer,
                                       int* ilevel_f,
                                       int* max_num_interface_p_refine_ibool_f) {

// asynchronous transfer from device to host

  TRACE("assemble_mpi_device_lts");

  Mesh* mp = (Mesh*)(*Mesh_pointer); // get Mesh from fortran integer wrapper
  int ilevel = *ilevel_f;
  int max_num_interface_p_refine_ibool = *max_num_interface_p_refine_ibool_f;

  if (mp->size_mpi_buffer > 0){

    gpuStreamSynchronize(mp->copy_stream);

    int blocksize = BLOCKSIZE_TRANSFER;

    int size_padded = ((int)ceil(((double)max_num_interface_p_refine_ibool)/((double)blocksize)))*blocksize;

    int num_blocks_x, num_blocks_y;
    get_blocks_xy(size_padded/blocksize,&num_blocks_x,&num_blocks_y);

    dim3 grid(num_blocks_x,num_blocks_y);
    dim3 threads(blocksize,1,1);

#ifdef USE_CUDA
    if (run_cuda){
      assemble_boundary_lts_accel_on_device<<<grid,threads,0,mp->compute_stream>>>(mp->d_accel,
                                                                                   mp->d_send_accel_buffer,
                                                                                   mp->num_interfaces_ext_mesh,
                                                                                   mp->d_lts_num_interface_p_refine_ibool,
                                                                                   mp->d_lts_interface_p_refine_ibool,
                                                                                   mp->max_nibool_interfaces_ext_mesh,
                                                                                   ilevel);
    }
#endif
#ifdef USE_HIP
    if (run_hip){
      hipLaunchKernelGGL(assemble_boundary_lts_accel_on_device, dim3(grid), dim3(threads), 0, mp->compute_stream,
                                                                                   mp->d_accel,
                                                                                   mp->d_send_accel_buffer,
                                                                                   mp->num_interfaces_ext_mesh,
                                                                                   mp->d_lts_num_interface_p_refine_ibool,
                                                                                   mp->d_lts_interface_p_refine_ibool,
                                                                                   mp->max_nibool_interfaces_ext_mesh,
                                                                                   ilevel);
    }
#endif
    GPU_ERROR_CHECKING("assemble_boundary_lts_accel_on_device");

    // wait until kernel is finished before starting async memcpy
    gpuSynchronize();
  }
}

/* ----------------------------------------------------------------------------------------------- */

extern EXTERN_LANG
void FC_FUNC_(assemble_reduced_mpi_device_lts,
              ASSEMBLE_REDUCED_MPI_DEVICE_LTS)(long* Mesh_pointer,
                                               int* ilevel_f,
                                               int* num_interface_p_refine_boundary_f) {

  // asynchronous transfer from device to host

  TRACE("assemble_reduced_mpi_device_lts");

  Mesh* mp = (Mesh*)(*Mesh_pointer); // get Mesh from fortran integer wrapper
  int ilevel = *ilevel_f;
  int num_interface_p_refine_boundary = *num_interface_p_refine_boundary_f;

  // checks if anything to do
  if (num_interface_p_refine_boundary == 0) return;

  if (mp->size_mpi_buffer > 0){

    gpuStreamSynchronize(mp->copy_stream);

    int blocksize = BLOCKSIZE_TRANSFER;

    int size_padded = ((int)ceil(((double)num_interface_p_refine_boundary)/((double)blocksize)))*blocksize;

    int num_blocks_x, num_blocks_y;
    get_blocks_xy(size_padded/blocksize,&num_blocks_x,&num_blocks_y);

    dim3 grid(num_blocks_x,num_blocks_y);
    dim3 threads(blocksize,1,1);

#ifdef USE_CUDA
    if (run_cuda){
      assemble_reduced_boundary_lts_accel_on_device<<<grid,threads,0,mp->compute_stream>>>(mp->d_accel,
                                                                                           mp->d_send_accel_buffer,
                                                                                           num_interface_p_refine_boundary,
                                                                                           mp->d_lts_interface_p_refine_boundary,
                                                                                           mp->lts_max_nibool_interfaces_boundary,
                                                                                           ilevel);
    }
#endif
#ifdef USE_HIP
    if (run_hip){
      hipLaunchKernelGGL(assemble_reduced_boundary_lts_accel_on_device, dim3(grid), dim3(threads), 0, mp->compute_stream,
                                                                                           mp->d_accel,
                                                                                           mp->d_send_accel_buffer,
                                                                                           num_interface_p_refine_boundary,
                                                                                           mp->d_lts_interface_p_refine_boundary,
                                                                                           mp->lts_max_nibool_interfaces_boundary,
                                                                                           ilevel);
    }
#endif
    GPU_ERROR_CHECKING("assemble_reduced_boundary_lts_accel_on_device");

    // wait until kernel is finished before starting async memcpy
    gpuSynchronize();
  }
}
