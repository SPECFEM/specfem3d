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

// Stacey absorbing boundary - acoustic domains

/* ----------------------------------------------------------------------------------------------- */

extern EXTERN_LANG
void FC_FUNC_(compute_stacey_acoustic_cuda,
              COMPUTE_STACEY_ACOUSTIC_CUDA)(long* Mesh_pointer,
                                            int* iphasef,
                                            realw* h_b_absorb_potential,
                                            int* FORWARD_OR_ADJOINT_f) {

  TRACE("compute_stacey_acoustic_cuda");

  //double start_time = get_time_val();

  Mesh* mp = (Mesh*)(*Mesh_pointer); //get mesh pointer out of fortran integer container
  int FORWARD_OR_ADJOINT = *FORWARD_OR_ADJOINT_f;

  // checks if anything to do
  if (mp->d_num_abs_boundary_faces == 0) return;

  int iphase          = *iphasef;

  // only add these contributions in first pass
  if (iphase != 1) return;

  // way 1: Elapsed time: 4.385948e-03
  // > NGLLSQUARE==NGLL2==25, but we handle this inside kernel
  //  int blocksize = 32;

  // way 2: Elapsed time: 4.379034e-03
  // > NGLLSQUARE==NGLL2==25, no further check inside kernel
  int blocksize = NGLL2;

  int num_blocks_x, num_blocks_y;
  get_blocks_xy(mp->d_num_abs_boundary_faces,&num_blocks_x,&num_blocks_y);

  dim3 grid(num_blocks_x,num_blocks_y);
  dim3 threads(blocksize,1,1);

  // safety check
  if (FORWARD_OR_ADJOINT != 0 && FORWARD_OR_ADJOINT != 1 && FORWARD_OR_ADJOINT != 3) {
    exit_on_error("Error invalid FORWARD_OR_ADJOINT in compute_stacey_acoustic_cuda() routine");
  }

  //  adjoint simulations: reads in absorbing boundary
  if (mp->simulation_type == 3 && FORWARD_OR_ADJOINT != 1){
    // copies array to GPU
    gpuMemcpy_todevice_void((void*)mp->d_b_absorb_potential,(void*)h_b_absorb_potential,mp->d_b_reclen_potential);
  }

  if (FORWARD_OR_ADJOINT == 0){
    // combined forward/backward fields
#ifdef USE_CUDA
    if (run_cuda){
      compute_stacey_acoustic_kernel<<<grid,threads>>>(mp->d_potential_dot_acoustic,
                                                       mp->d_potential_dot_dot_acoustic,
                                                       mp->d_abs_boundary_ispec,
                                                       mp->d_abs_boundary_ijk,
                                                       mp->d_abs_boundary_jacobian2Dw,
                                                       mp->d_ibool,
                                                       mp->d_rhostore,
                                                       mp->d_kappastore,
                                                       mp->d_ispec_is_acoustic,
                                                       mp->simulation_type,
                                                       mp->save_forward,
                                                       mp->d_num_abs_boundary_faces,
                                                       mp->d_b_potential_dot_acoustic,
                                                       mp->d_b_potential_dot_dot_acoustic,
                                                       mp->d_b_absorb_potential,
                                                       mp->gravity);
    }
#endif
#ifdef USE_HIP
    if (run_hip){
      hipLaunchKernelGGL(compute_stacey_acoustic_kernel, dim3(grid), dim3(threads), 0, 0,
                                                         mp->d_potential_dot_acoustic,
                                                         mp->d_potential_dot_dot_acoustic,
                                                         mp->d_abs_boundary_ispec,
                                                         mp->d_abs_boundary_ijk,
                                                         mp->d_abs_boundary_jacobian2Dw,
                                                         mp->d_ibool,
                                                         mp->d_rhostore,
                                                         mp->d_kappastore,
                                                         mp->d_ispec_is_acoustic,
                                                         mp->simulation_type,
                                                         mp->save_forward,
                                                         mp->d_num_abs_boundary_faces,
                                                         mp->d_b_potential_dot_acoustic,
                                                         mp->d_b_potential_dot_dot_acoustic,
                                                         mp->d_b_absorb_potential,
                                                         mp->gravity);
    }
#endif

  }else{
    // sets gpu arrays
    field *potential_dot, *potential_dot_dot;
    if (FORWARD_OR_ADJOINT == 1) {
      potential_dot = mp->d_potential_dot_acoustic;
      potential_dot_dot = mp->d_potential_dot_dot_acoustic;
    } else {
      // for backward/reconstructed fields
      potential_dot = mp->d_b_potential_dot_acoustic;
      potential_dot_dot = mp->d_b_potential_dot_dot_acoustic;
    }
    // single forward or backward fields
#ifdef USE_CUDA
    if (run_cuda){
      compute_stacey_acoustic_single_kernel<<<grid,threads>>>(potential_dot,
                                                              potential_dot_dot,
                                                              mp->d_abs_boundary_ispec,
                                                              mp->d_abs_boundary_ijk,
                                                              mp->d_abs_boundary_jacobian2Dw,
                                                              mp->d_ibool,
                                                              mp->d_rhostore,
                                                              mp->d_kappastore,
                                                              mp->d_ispec_is_acoustic,
                                                              FORWARD_OR_ADJOINT,
                                                              mp->simulation_type,
                                                              mp->save_forward,
                                                              mp->d_num_abs_boundary_faces,
                                                              mp->d_b_absorb_potential,
                                                              mp->gravity);
    }
#endif
#ifdef USE_HIP
    if (run_hip){
      hipLaunchKernelGGL(compute_stacey_acoustic_single_kernel, dim3(grid), dim3(threads), 0, 0,
                                                                potential_dot,
                                                                potential_dot_dot,
                                                                mp->d_abs_boundary_ispec,
                                                                mp->d_abs_boundary_ijk,
                                                                mp->d_abs_boundary_jacobian2Dw,
                                                                mp->d_ibool,
                                                                mp->d_rhostore,
                                                                mp->d_kappastore,
                                                                mp->d_ispec_is_acoustic,
                                                                FORWARD_OR_ADJOINT,
                                                                mp->simulation_type,
                                                                mp->save_forward,
                                                                mp->d_num_abs_boundary_faces,
                                                                mp->d_b_absorb_potential,
                                                                mp->gravity);
    }
#endif
  }

  //  adjoint simulations: stores absorbed wavefield part
  if (mp->simulation_type == 1 && mp->save_forward){
    // (cudaMemcpy implicitly synchronizes all other cuda operations)
    // copies array to CPU
    gpuMemcpy_tohost_void((void*)h_b_absorb_potential,(void*)mp->d_b_absorb_potential,mp->d_b_reclen_potential);
  }

  GPU_ERROR_CHECKING("compute_stacey_acoustic_kernel");
}

/* ----------------------------------------------------------------------------------------------- */

extern EXTERN_LANG
void FC_FUNC_(compute_stacey_acoustic_undoatt_cuda,
              COMPUTE_STACEY_ACOUSTIC_UNDOATT_CUDA)(long* Mesh_pointer,
                                                     int* iphasef,
                                                     int* FORWARD_OR_ADJOINT_f) {

  TRACE("compute_stacey_acoustic_undoatt_cuda");

  //double start_time = get_time_val();

  Mesh* mp = (Mesh*)(*Mesh_pointer); //get mesh pointer out of fortran integer container
  int FORWARD_OR_ADJOINT = *FORWARD_OR_ADJOINT_f;

  // safety check
  if (FORWARD_OR_ADJOINT != 1 && FORWARD_OR_ADJOINT != 3) {
    exit_on_error("Error invalid FORWARD_OR_ADJOINT in compute_stacey_acoustic_undoatt_cuda() routine");
  }

  // checks if anything to do
  if (mp->d_num_abs_boundary_faces == 0) return;

  int iphase = *iphasef;

  // only add these contributions in first pass
  if (iphase != 1) return;

  // way 1: Elapsed time: 4.385948e-03
  // > NGLLSQUARE==NGLL2==25, but we handle this inside kernel
  //  int blocksize = 32;

  // way 2: Elapsed time: 4.379034e-03
  // > NGLLSQUARE==NGLL2==25, no further check inside kernel
  int blocksize = NGLL2;

  int num_blocks_x, num_blocks_y;
  get_blocks_xy(mp->d_num_abs_boundary_faces,&num_blocks_x,&num_blocks_y);

  dim3 grid(num_blocks_x,num_blocks_y);
  dim3 threads(blocksize,1,1);

  // no absorbing boundary need to be stored, only propagates forward in time
  // sets gpu arrays
  field *potential_dot, *potential_dot_dot;
  if (FORWARD_OR_ADJOINT == 1) {
    potential_dot = mp->d_potential_dot_acoustic;
    potential_dot_dot = mp->d_potential_dot_dot_acoustic;
  } else {
    // for backward/reconstructed fields
    potential_dot = mp->d_b_potential_dot_acoustic;
    potential_dot_dot = mp->d_b_potential_dot_dot_acoustic;
  }

  // undoatt: single forward or backward fields
#ifdef USE_CUDA
  if (run_cuda){
    compute_stacey_acoustic_undoatt_kernel<<<grid,threads>>>(potential_dot,
                                                             potential_dot_dot,
                                                             mp->d_abs_boundary_ispec,
                                                             mp->d_abs_boundary_ijk,
                                                             mp->d_abs_boundary_jacobian2Dw,
                                                             mp->d_ibool,
                                                             mp->d_rhostore,
                                                             mp->d_kappastore,
                                                             mp->d_ispec_is_acoustic,
                                                             mp->d_num_abs_boundary_faces,
                                                             mp->gravity);
  }
#endif
#ifdef USE_HIP
  if (run_hip){
    hipLaunchKernelGGL(compute_stacey_acoustic_undoatt_kernel, dim3(grid), dim3(threads), 0, 0,
                                                               potential_dot,
                                                               potential_dot_dot,
                                                               mp->d_abs_boundary_ispec,
                                                               mp->d_abs_boundary_ijk,
                                                               mp->d_abs_boundary_jacobian2Dw,
                                                               mp->d_ibool,
                                                               mp->d_rhostore,
                                                               mp->d_kappastore,
                                                               mp->d_ispec_is_acoustic,
                                                               mp->d_num_abs_boundary_faces,
                                                               mp->gravity);
  }
#endif

  //  no need to store absorbed wavefield part

  GPU_ERROR_CHECKING("compute_stacey_acoustic_undoatt_cuda");
}

