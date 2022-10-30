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

// Stacey absorbing boundary - elastic domains

/* ----------------------------------------------------------------------------------------------- */

extern EXTERN_LANG
void FC_FUNC_(compute_stacey_viscoelastic_cuda,
              COMPUTE_STACEY_VISCOELASTIC_CUDA)(long* Mesh_pointer,
                                                int* iphasef,
                                                realw* b_absorb_field,
                                                int* FORWARD_OR_ADJOINT_f) {

  TRACE("compute_stacey_viscoelastic_cuda");

  Mesh* mp = (Mesh*)(*Mesh_pointer); //get mesh pointer out of fortran integer container
  int FORWARD_OR_ADJOINT = *FORWARD_OR_ADJOINT_f;

  // checks if anything to do
  if (mp->d_num_abs_boundary_faces == 0) return;

  int iphase = *iphasef;

  // only add these contributions in first pass
  if (iphase != 1) return;

  // way 1
  // > NGLLSQUARE==NGLL2==25, but we handle this inside kernel
  //int blocksize = 32;

  // way 2: seems sligthly faster
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
    // reading is done in fortran routine
    gpuMemcpy_todevice_void((void*)mp->d_b_absorb_field,(void*)b_absorb_field,mp->d_b_reclen_field);
  }

  GPU_ERROR_CHECKING("between cudamemcpy and compute_stacey_elastic_kernel");

  if (FORWARD_OR_ADJOINT == 0){
    // combined forward/backward fields
#ifdef USE_CUDA
    if (run_cuda){
      compute_stacey_elastic_kernel<<<grid,threads,0,mp->compute_stream>>>(mp->d_veloc,
                                                                           mp->d_accel,
                                                                           mp->d_abs_boundary_ispec,
                                                                           mp->d_abs_boundary_ijk,
                                                                           mp->d_abs_boundary_normal,
                                                                           mp->d_abs_boundary_jacobian2Dw,
                                                                           mp->d_ibool,
                                                                           mp->d_rho_vp,
                                                                           mp->d_rho_vs,
                                                                           mp->d_ispec_is_elastic,
                                                                           mp->simulation_type,
                                                                           mp->save_forward,
                                                                           mp->d_num_abs_boundary_faces,
                                                                           mp->d_b_absorb_field);
    }
#endif
#ifdef USE_HIP
    if (run_hip){
      hipLaunchKernelGGL(compute_stacey_elastic_kernel, dim3(grid), dim3(threads), 0, mp->compute_stream,
                                                        mp->d_veloc,
                                                        mp->d_accel,
                                                        mp->d_abs_boundary_ispec,
                                                        mp->d_abs_boundary_ijk,
                                                        mp->d_abs_boundary_normal,
                                                        mp->d_abs_boundary_jacobian2Dw,
                                                        mp->d_ibool,
                                                        mp->d_rho_vp,
                                                        mp->d_rho_vs,
                                                        mp->d_ispec_is_elastic,
                                                        mp->simulation_type,
                                                        mp->save_forward,
                                                        mp->d_num_abs_boundary_faces,
                                                        mp->d_b_absorb_field);
    }
#endif

    // adjoint simulations
    if (mp->simulation_type == 3){
#ifdef USE_CUDA
      if (run_cuda){
        compute_stacey_elastic_sim3_kernel<<<grid,threads,0,mp->compute_stream>>>(mp->d_abs_boundary_ispec,
                                                                                  mp->d_abs_boundary_ijk,
                                                                                  mp->d_ibool,
                                                                                  mp->d_ispec_is_elastic,
                                                                                  mp->d_num_abs_boundary_faces,
                                                                                  mp->d_b_accel,
                                                                                  mp->d_b_absorb_field);
      }
#endif
#ifdef USE_HIP
      if (run_hip){
        hipLaunchKernelGGL(compute_stacey_elastic_sim3_kernel, dim3(grid), dim3(threads), 0, mp->compute_stream,
                                                               mp->d_abs_boundary_ispec,
                                                               mp->d_abs_boundary_ijk,
                                                               mp->d_ibool,
                                                               mp->d_ispec_is_elastic,
                                                               mp->d_num_abs_boundary_faces,
                                                               mp->d_b_accel,
                                                               mp->d_b_absorb_field);
      }
#endif

    }
  }else{
    // sets gpu arrays
    realw *veloc, *accel;
    if (FORWARD_OR_ADJOINT == 1) {
      veloc = mp->d_veloc;
      accel = mp->d_accel;
    } else {
      // for backward/reconstructed fields
      veloc = mp->d_b_veloc;
      accel = mp->d_b_accel;
    }
    // single forward or backward fields
#ifdef USE_CUDA
    if (run_cuda){
      compute_stacey_elastic_single_kernel<<<grid,threads,0,mp->compute_stream>>>(veloc,
                                                                                  accel,
                                                                                  mp->d_abs_boundary_ispec,
                                                                                  mp->d_abs_boundary_ijk,
                                                                                  mp->d_abs_boundary_normal,
                                                                                  mp->d_abs_boundary_jacobian2Dw,
                                                                                  mp->d_ibool,
                                                                                  mp->d_rho_vp,
                                                                                  mp->d_rho_vs,
                                                                                  mp->d_ispec_is_elastic,
                                                                                  FORWARD_OR_ADJOINT,
                                                                                  mp->simulation_type,
                                                                                  mp->save_forward,
                                                                                  mp->d_num_abs_boundary_faces,
                                                                                  mp->d_b_absorb_field);
    }
#endif
#ifdef USE_HIP
    if (run_hip){
      hipLaunchKernelGGL(compute_stacey_elastic_single_kernel, dim3(grid), dim3(threads), 0, mp->compute_stream,
                                                               veloc,
                                                               accel,
                                                               mp->d_abs_boundary_ispec,
                                                               mp->d_abs_boundary_ijk,
                                                               mp->d_abs_boundary_normal,
                                                               mp->d_abs_boundary_jacobian2Dw,
                                                               mp->d_ibool,
                                                               mp->d_rho_vp,
                                                               mp->d_rho_vs,
                                                               mp->d_ispec_is_elastic,
                                                               FORWARD_OR_ADJOINT,
                                                               mp->simulation_type,
                                                               mp->save_forward,
                                                               mp->d_num_abs_boundary_faces,
                                                               mp->d_b_absorb_field);
    }
#endif
  }

  GPU_ERROR_CHECKING("compute_stacey_elastic_kernel");

  if (mp->simulation_type == 1 && mp->save_forward) {
    // explicitly wait until compute stream is done
    // (cudaMemcpy implicitly synchronizes all other cuda operations)
    gpuStreamSynchronize(mp->compute_stream);

    // copies absorb_field values to CPU
    gpuMemcpy_tohost_void((void*)b_absorb_field,(void*)mp->d_b_absorb_field,mp->d_b_reclen_field);
    // writing is done in fortran routine
  }

  GPU_ERROR_CHECKING("after compute_stacey_elastic after cudamemcpy");
}

/* ----------------------------------------------------------------------------------------------- */

extern EXTERN_LANG
void FC_FUNC_(compute_stacey_viscoelastic_undoatt_cuda,
              COMPUTE_STACEY_VISCOELASTIC_UNDOATT_CUDA)(long* Mesh_pointer,
                                                        int* iphasef,
                                                        int* FORWARD_OR_ADJOINT_f) {

  TRACE("compute_stacey_viscoelastic_undoatt_cuda");

  Mesh* mp = (Mesh*)(*Mesh_pointer); //get mesh pointer out of fortran integer container
  int FORWARD_OR_ADJOINT = *FORWARD_OR_ADJOINT_f;
  // safety check
  if (FORWARD_OR_ADJOINT != 1 && FORWARD_OR_ADJOINT != 3) {
    exit_on_error("Error invalid FORWARD_OR_ADJOINT in compute_stacey_acoustic_undoatt_cuda() routine");
  }

  // checks if anything to do
  if (mp->d_num_abs_boundary_faces == 0) return;

  int iphase    = *iphasef;

  // only add these contributions in first pass
  if (iphase != 1) return;

  // way 1
  // > NGLLSQUARE==NGLL2==25, but we handle this inside kernel
  //int blocksize = 32;

  // way 2: seems sligthly faster
  // > NGLLSQUARE==NGLL2==25, no further check inside kernel
  int blocksize = NGLL2;

  int num_blocks_x, num_blocks_y;
  get_blocks_xy(mp->d_num_abs_boundary_faces,&num_blocks_x,&num_blocks_y);

  dim3 grid(num_blocks_x,num_blocks_y);
  dim3 threads(blocksize,1,1);

  realw *veloc, *accel;
  if (FORWARD_OR_ADJOINT == 1) {
    veloc = mp->d_veloc;
    accel = mp->d_accel;
  } else {
    // for backward/reconstructed fields
    veloc = mp->d_b_veloc;
    accel = mp->d_b_accel;
  }

  // undoatt: single forward or backward fields
#ifdef USE_CUDA
  if (run_cuda){
    compute_stacey_elastic_undoatt_kernel<<<grid,threads,0,mp->compute_stream>>>(veloc,
                                                                                 accel,
                                                                                 mp->d_abs_boundary_ispec,
                                                                                 mp->d_abs_boundary_ijk,
                                                                                 mp->d_abs_boundary_normal,
                                                                                 mp->d_abs_boundary_jacobian2Dw,
                                                                                 mp->d_ibool,
                                                                                 mp->d_rho_vp,
                                                                                 mp->d_rho_vs,
                                                                                 mp->d_ispec_is_elastic,
                                                                                 mp->d_num_abs_boundary_faces);
  }
#endif
#ifdef USE_HIP
  if (run_hip){
    hipLaunchKernelGGL(compute_stacey_elastic_undoatt_kernel, dim3(grid), dim3(threads), 0, mp->compute_stream,
                                                              veloc,
                                                              accel,
                                                              mp->d_abs_boundary_ispec,
                                                              mp->d_abs_boundary_ijk,
                                                              mp->d_abs_boundary_normal,
                                                              mp->d_abs_boundary_jacobian2Dw,
                                                              mp->d_ibool,
                                                              mp->d_rho_vp,
                                                              mp->d_rho_vs,
                                                              mp->d_ispec_is_elastic,
                                                              mp->d_num_abs_boundary_faces);
  }
#endif

  GPU_ERROR_CHECKING("after compute_stacey_elastic_undoatt_kernel");
}

