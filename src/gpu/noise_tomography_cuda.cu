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

// randomize displ for testing
extern EXTERN_LANG
void FC_FUNC_(make_displ_rand,MAKE_DISPL_RAND)(long* Mesh_pointer,realw* h_displ) {
TRACE("make_displ_rand");

  Mesh* mp = (Mesh*)(*Mesh_pointer); // get Mesh from fortran integer wrapper
  // realw* displ_rnd = (realw*)malloc(mp->NGLOB_AB*3*sizeof(realw));
  for(int i=0;i<mp->NGLOB_AB*3;i++) {
    h_displ[i] = rand();
  }
  gpuMemcpy_todevice_realw(mp->d_displ,h_displ,mp->NGLOB_AB*3);
}

/* ----------------------------------------------------------------------------------------------- */


extern EXTERN_LANG
void FC_FUNC_(transfer_surface_to_host,
              TRANSFER_SURFACE_TO_HOST)(long* Mesh_pointer,
                                        realw* h_noise_surface_movie) {
TRACE("transfer_surface_to_host");

  Mesh* mp = (Mesh*)(*Mesh_pointer); // get Mesh from fortran integer wrapper

  int num_blocks_x, num_blocks_y;
  get_blocks_xy(mp->num_free_surface_faces,&num_blocks_x,&num_blocks_y);

  dim3 grid(num_blocks_x,num_blocks_y,1);
  dim3 threads(NGLL2,1,1);

#ifdef USE_CUDA
  if (run_cuda){
    transfer_surface_to_host_kernel<<<grid,threads>>>(mp->d_free_surface_ispec,
                                                      mp->d_free_surface_ijk,
                                                      mp->num_free_surface_faces,
                                                      mp->d_ibool,
                                                      mp->d_displ,
                                                      mp->d_noise_surface_movie);
  }
#endif
#ifdef USE_HIP
  if (run_hip){
    hipLaunchKernelGGL(transfer_surface_to_host_kernel, dim3(grid), dim3(threads), 0, 0,
                                                        mp->d_free_surface_ispec,
                                                        mp->d_free_surface_ijk,
                                                        mp->num_free_surface_faces,
                                                        mp->d_ibool,
                                                        mp->d_displ,
                                                        mp->d_noise_surface_movie);
  }
#endif

  gpuMemcpy_tohost_realw(h_noise_surface_movie,mp->d_noise_surface_movie,3*NGLL2*(mp->num_free_surface_faces));

  GPU_ERROR_CHECKING("transfer_surface_to_host");
}

/* ----------------------------------------------------------------------------------------------- */


extern EXTERN_LANG
void FC_FUNC_(noise_read_add_surface_movie_cu,
              NOISE_READ_ADD_SURFACE_MOVIE_CU)(long* Mesh_pointer,
                                               realw* h_noise_surface_movie,
                                               int* NOISE_TOMOGRAPHYf) {
  TRACE("noise_read_add_surface_movie_cu");

  // EPIK_TRACER("noise_read_add_surface_movie_cu");

  Mesh* mp = (Mesh*)(*Mesh_pointer); //get mesh pointer out of fortran integer container
  int NOISE_TOMOGRAPHY = *NOISE_TOMOGRAPHYf;

  gpuMemcpy_todevice_realw(mp->d_noise_surface_movie,h_noise_surface_movie,3*NGLL2*(mp->num_free_surface_faces));

  int num_blocks_x, num_blocks_y;
  get_blocks_xy(mp->num_free_surface_faces,&num_blocks_x,&num_blocks_y);

  dim3 grid(num_blocks_x,num_blocks_y,1);
  dim3 threads(NGLL2,1,1);

  if (NOISE_TOMOGRAPHY == 2) {
    // add surface source to forward field
#ifdef USE_CUDA
    if (run_cuda){
      noise_read_add_surface_movie_cuda_kernel<<<grid,threads>>>(mp->d_accel,
                                                                 mp->d_ibool,
                                                                 mp->d_free_surface_ispec,
                                                                 mp->d_free_surface_ijk,
                                                                 mp->num_free_surface_faces,
                                                                 mp->d_noise_surface_movie,
                                                                 mp->d_normal_x_noise,
                                                                 mp->d_normal_y_noise,
                                                                 mp->d_normal_z_noise,
                                                                 mp->d_mask_noise,
                                                                 mp->d_free_surface_jacobian2Dw);
    }
#endif
#ifdef USE_HIP
    if (run_hip){
      hipLaunchKernelGGL(noise_read_add_surface_movie_cuda_kernel, dim3(grid), dim3(threads), 0, 0,
                                                                   mp->d_accel,
                                                                   mp->d_ibool,
                                                                   mp->d_free_surface_ispec,
                                                                   mp->d_free_surface_ijk,
                                                                   mp->num_free_surface_faces,
                                                                   mp->d_noise_surface_movie,
                                                                   mp->d_normal_x_noise,
                                                                   mp->d_normal_y_noise,
                                                                   mp->d_normal_z_noise,
                                                                   mp->d_mask_noise,
                                                                   mp->d_free_surface_jacobian2Dw);
    }
#endif

  }
  else if (NOISE_TOMOGRAPHY == 3) {
    // add surface source to adjoint (backward) field
#ifdef USE_CUDA
    if (run_cuda){
      noise_read_add_surface_movie_cuda_kernel<<<grid,threads>>>(mp->d_b_accel,
                                                                 mp->d_ibool,
                                                                 mp->d_free_surface_ispec,
                                                                 mp->d_free_surface_ijk,
                                                                 mp->num_free_surface_faces,
                                                                 mp->d_noise_surface_movie,
                                                                 mp->d_normal_x_noise,
                                                                 mp->d_normal_y_noise,
                                                                 mp->d_normal_z_noise,
                                                                 mp->d_mask_noise,
                                                                 mp->d_free_surface_jacobian2Dw);
    }
#endif
#ifdef USE_HIP
    if (run_hip){
      hipLaunchKernelGGL(noise_read_add_surface_movie_cuda_kernel, dim3(grid), dim3(threads), 0, 0,
                                                                   mp->d_b_accel,
                                                                   mp->d_ibool,
                                                                   mp->d_free_surface_ispec,
                                                                   mp->d_free_surface_ijk,
                                                                   mp->num_free_surface_faces,
                                                                   mp->d_noise_surface_movie,
                                                                   mp->d_normal_x_noise,
                                                                   mp->d_normal_y_noise,
                                                                   mp->d_normal_z_noise,
                                                                   mp->d_mask_noise,
                                                                   mp->d_free_surface_jacobian2Dw);
    }
#endif

  }

  GPU_ERROR_CHECKING("noise_read_add_surface_movie_cuda_kernel");
}
