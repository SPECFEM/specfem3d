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

// ELASTIC SIMULATIONS

/* ----------------------------------------------------------------------------------------------- */

extern EXTERN_LANG
void FC_FUNC_(compute_kernels_elastic_cuda,
              COMPUTE_KERNELS_ELASTIC_CUDA)(long* Mesh_pointer,
                                            realw* deltat_f) {

  TRACE("compute_kernels_elastic_cuda");

  Mesh* mp = (Mesh*)(*Mesh_pointer); //get mesh pointer out of fortran integer container

  int blocksize = NGLL3; // NGLLX*NGLLY*NGLLZ
  realw deltat = *deltat_f;

  int num_blocks_x, num_blocks_y;
  get_blocks_xy(mp->NSPEC_AB,&num_blocks_x,&num_blocks_y);

  dim3 grid(num_blocks_x,num_blocks_y);
  dim3 threads(blocksize,1,1);

  // current strain field
  if (mp->undo_attenuation){
    // simulations with UNDO_ATTENUATION save as much memory as possible;
    // backward/reconstructed wavefield strain will be re-computed locally here
#ifdef USE_CUDA
    if (run_cuda){
      compute_element_strain_cudakernel<<<grid,threads>>>(mp->d_ispec_is_elastic,mp->d_ibool,
                                                          mp->d_b_displ,
                                                          mp->d_b_epsilondev_xx,
                                                          mp->d_b_epsilondev_yy,
                                                          mp->d_b_epsilondev_xy,
                                                          mp->d_b_epsilondev_xz,
                                                          mp->d_b_epsilondev_yz,
                                                          mp->d_b_epsilondev_trace,
                                                          mp->d_b_epsilon_trace_over_3,
                                                          mp->d_xix,mp->d_xiy,mp->d_xiz,
                                                          mp->d_etax,mp->d_etay,mp->d_etaz,
                                                          mp->d_gammax,mp->d_gammay,mp->d_gammaz,
                                                          mp->d_irregular_element_number,
                                                          mp->xix_regular,
                                                          mp->d_hprime_xx,
                                                          mp->NSPEC_AB);
    }
#endif
#ifdef USE_HIP
    if (run_hip){
      hipLaunchKernelGGL(compute_element_strain_cudakernel, dim3(grid), dim3(threads), 0, 0,
                                                            mp->d_ispec_is_elastic,mp->d_ibool,
                                                            mp->d_b_displ,
                                                            mp->d_b_epsilondev_xx,
                                                            mp->d_b_epsilondev_yy,
                                                            mp->d_b_epsilondev_xy,
                                                            mp->d_b_epsilondev_xz,
                                                            mp->d_b_epsilondev_yz,
                                                            mp->d_b_epsilondev_trace,
                                                            mp->d_b_epsilon_trace_over_3,
                                                            mp->d_xix,mp->d_xiy,mp->d_xiz,
                                                            mp->d_etax,mp->d_etay,mp->d_etaz,
                                                            mp->d_gammax,mp->d_gammay,mp->d_gammaz,
                                                            mp->d_irregular_element_number,
                                                            mp->xix_regular,
                                                            mp->d_hprime_xx,
                                                            mp->NSPEC_AB);
    }
#endif
  }

  // elastic kernels
  if (mp->anisotropic_kl ){
    // anisotropic kernel
#ifdef USE_CUDA
    if (run_cuda){
      compute_kernels_ani_cudakernel<<<grid,threads>>>(mp->d_ispec_is_elastic,mp->d_ibool,
                                                       mp->d_accel, mp->d_b_displ,
                                                       mp->d_epsilondev_xx,
                                                       mp->d_epsilondev_yy,
                                                       mp->d_epsilondev_xy,
                                                       mp->d_epsilondev_xz,
                                                       mp->d_epsilondev_yz,
                                                       mp->d_b_epsilondev_xx,
                                                       mp->d_b_epsilondev_yy,
                                                       mp->d_b_epsilondev_xy,
                                                       mp->d_b_epsilondev_xz,
                                                       mp->d_b_epsilondev_yz,
                                                       mp->d_rho_kl,
                                                       deltat,
                                                       mp->d_cijkl_kl,
                                                       mp->d_epsilon_trace_over_3,
                                                       mp->d_b_epsilon_trace_over_3,
                                                       mp->NSPEC_AB);
    }
#endif
#ifdef USE_HIP
    if (run_hip){
      hipLaunchKernelGGL(compute_kernels_ani_cudakernel, dim3(grid), dim3(threads), 0, 0,
                                                         mp->d_ispec_is_elastic,mp->d_ibool,
                                                         mp->d_accel, mp->d_b_displ,
                                                         mp->d_epsilondev_xx,
                                                         mp->d_epsilondev_yy,
                                                         mp->d_epsilondev_xy,
                                                         mp->d_epsilondev_xz,
                                                         mp->d_epsilondev_yz,
                                                         mp->d_b_epsilondev_xx,
                                                         mp->d_b_epsilondev_yy,
                                                         mp->d_b_epsilondev_xy,
                                                         mp->d_b_epsilondev_xz,
                                                         mp->d_b_epsilondev_yz,
                                                         mp->d_rho_kl,
                                                         deltat,
                                                         mp->d_cijkl_kl,
                                                         mp->d_epsilon_trace_over_3,
                                                         mp->d_b_epsilon_trace_over_3,
                                                         mp->NSPEC_AB);
    }
#endif

  }else{
    // isotropic kernel
#ifdef USE_CUDA
    if (run_cuda){
      compute_kernels_cudakernel<<<grid,threads>>>(mp->d_ispec_is_elastic,mp->d_ibool,
                                                   mp->d_accel, mp->d_b_displ,
                                                   mp->d_epsilondev_xx,
                                                   mp->d_epsilondev_yy,
                                                   mp->d_epsilondev_xy,
                                                   mp->d_epsilondev_xz,
                                                   mp->d_epsilondev_yz,
                                                   mp->d_b_epsilondev_xx,
                                                   mp->d_b_epsilondev_yy,
                                                   mp->d_b_epsilondev_xy,
                                                   mp->d_b_epsilondev_xz,
                                                   mp->d_b_epsilondev_yz,
                                                   mp->d_rho_kl,
                                                   deltat,
                                                   mp->d_mu_kl,
                                                   mp->d_kappa_kl,
                                                   mp->d_epsilon_trace_over_3,
                                                   mp->d_b_epsilon_trace_over_3,
                                                   mp->NSPEC_AB);
    }
#endif
#ifdef USE_HIP
    if (run_hip){
      hipLaunchKernelGGL(compute_kernels_cudakernel, dim3(grid), dim3(threads), 0, 0,
                                                     mp->d_ispec_is_elastic,mp->d_ibool,
                                                     mp->d_accel, mp->d_b_displ,
                                                     mp->d_epsilondev_xx,
                                                     mp->d_epsilondev_yy,
                                                     mp->d_epsilondev_xy,
                                                     mp->d_epsilondev_xz,
                                                     mp->d_epsilondev_yz,
                                                     mp->d_b_epsilondev_xx,
                                                     mp->d_b_epsilondev_yy,
                                                     mp->d_b_epsilondev_xy,
                                                     mp->d_b_epsilondev_xz,
                                                     mp->d_b_epsilondev_yz,
                                                     mp->d_rho_kl,
                                                     deltat,
                                                     mp->d_mu_kl,
                                                     mp->d_kappa_kl,
                                                     mp->d_epsilon_trace_over_3,
                                                     mp->d_b_epsilon_trace_over_3,
                                                     mp->NSPEC_AB);
    }
#endif

  }

  GPU_ERROR_CHECKING("compute_kernels_elastic_cuda");
}


/* ----------------------------------------------------------------------------------------------- */

// NOISE SIMULATIONS

/* ----------------------------------------------------------------------------------------------- */

extern EXTERN_LANG
void FC_FUNC_(compute_kernels_strgth_noise_cu,
              COMPUTE_KERNELS_STRGTH_NOISE_CU)(long* Mesh_pointer,
                                                    realw* h_noise_surface_movie,
                                                    realw* deltat) {

  TRACE("compute_kernels_strgth_noise_cu");

  Mesh* mp = (Mesh*)(*Mesh_pointer); //get mesh pointer out of fortran integer container

  // checks if anything to do
  if (mp->num_free_surface_faces == 0) return;

  int num_blocks_x, num_blocks_y;
  get_blocks_xy(mp->num_free_surface_faces,&num_blocks_x,&num_blocks_y);

  dim3 grid(num_blocks_x,num_blocks_y);
  dim3 threads(NGLL2,1,1);

  gpuMemcpy_todevice_realw(mp->d_noise_surface_movie,h_noise_surface_movie,NDIM*NGLL2*(mp->num_free_surface_faces));

#ifdef USE_CUDA
  if (run_cuda){
    compute_kernels_strength_noise_cuda_kernel<<<grid,threads>>>(mp->d_displ,
                                                                 mp->d_free_surface_ispec,
                                                                 mp->d_free_surface_ijk,
                                                                 mp->d_ibool,
                                                                 mp->d_noise_surface_movie,
                                                                 mp->d_normal_x_noise,
                                                                 mp->d_normal_y_noise,
                                                                 mp->d_normal_z_noise,
                                                                 mp->d_sigma_kl,*deltat,
                                                                 mp->num_free_surface_faces);
  }
#endif
#ifdef USE_HIP
  if (run_hip){
    hipLaunchKernelGGL(compute_kernels_strength_noise_cuda_kernel, dim3(grid), dim3(threads), 0, 0,
                                                                   mp->d_displ,
                                                                   mp->d_free_surface_ispec,
                                                                   mp->d_free_surface_ijk,
                                                                   mp->d_ibool,
                                                                   mp->d_noise_surface_movie,
                                                                   mp->d_normal_x_noise,
                                                                   mp->d_normal_y_noise,
                                                                   mp->d_normal_z_noise,
                                                                   mp->d_sigma_kl,*deltat,
                                                                   mp->num_free_surface_faces);
  }
#endif

  GPU_ERROR_CHECKING("compute_kernels_strength_noise_cuda_kernel");
}



/* ----------------------------------------------------------------------------------------------- */

// ACOUSTIC SIMULATIONS

/* ----------------------------------------------------------------------------------------------- */

extern EXTERN_LANG
void FC_FUNC_(compute_kernels_acoustic_cuda,
              COMPUTE_KERNELS_ACOUSTIC_CUDA)(long* Mesh_pointer,
                                             realw* deltat_f) {

  TRACE("compute_kernels_acoustic_cuda");

  Mesh* mp = (Mesh*)(*Mesh_pointer); //get mesh pointer out of fortran integer container

  int blocksize = NGLL3; // NGLLX*NGLLY*NGLLZ
  realw deltat = *deltat_f;

  int num_blocks_x, num_blocks_y;
  get_blocks_xy(mp->NSPEC_AB,&num_blocks_x,&num_blocks_y);

  dim3 grid(num_blocks_x,num_blocks_y);
  dim3 threads(blocksize,1,1);

#ifdef USE_CUDA
  if (run_cuda){
    compute_kernels_acoustic_kernel<<<grid,threads>>>(mp->d_ispec_is_acoustic,
                                                      mp->d_ibool,
                                                      mp->d_rhostore,
                                                      mp->d_hprime_xx,
                                                      mp->d_irregular_element_number,
                                                      mp->d_xix,mp->d_xiy,mp->d_xiz,
                                                      mp->d_etax,mp->d_etay,mp->d_etaz,
                                                      mp->d_gammax,mp->d_gammay,mp->d_gammaz,
                                                      mp->xix_regular,
                                                      mp->d_potential_acoustic,
                                                      mp->d_potential_dot_dot_acoustic,
                                                      mp->d_b_potential_acoustic,
                                                      mp->d_b_potential_dot_dot_acoustic,
                                                      mp->d_rho_ac_kl,
                                                      mp->d_kappa_ac_kl,
                                                      deltat,
                                                      mp->NSPEC_AB,
                                                      mp->gravity);
  }
#endif
#ifdef USE_HIP
  if (run_hip){
    hipLaunchKernelGGL(compute_kernels_acoustic_kernel, dim3(grid), dim3(threads), 0, 0,
                                                      mp->d_ispec_is_acoustic,
                                                      mp->d_ibool,
                                                      mp->d_rhostore,
                                                      mp->d_hprime_xx,
                                                      mp->d_irregular_element_number,
                                                      mp->d_xix,mp->d_xiy,mp->d_xiz,
                                                      mp->d_etax,mp->d_etay,mp->d_etaz,
                                                      mp->d_gammax,mp->d_gammay,mp->d_gammaz,
                                                      mp->xix_regular,
                                                      mp->d_potential_acoustic,
                                                      mp->d_potential_dot_dot_acoustic,
                                                      mp->d_b_potential_acoustic,
                                                      mp->d_b_potential_dot_dot_acoustic,
                                                      mp->d_rho_ac_kl,
                                                      mp->d_kappa_ac_kl,
                                                      deltat,
                                                      mp->NSPEC_AB,
                                                      mp->gravity);
  }
#endif

  GPU_ERROR_CHECKING("compute_kernels_acoustic_kernel");
}

/* ----------------------------------------------------------------------------------------------- */

// preconditioner (approximate Hessian kernel)

/* ----------------------------------------------------------------------------------------------- */

extern EXTERN_LANG
void FC_FUNC_(compute_kernels_hess_cuda,
              COMPUTE_KERNELS_HESS_CUDA)(long* Mesh_pointer,
                                         realw* deltat_f,
                                         int* ELASTIC_SIMULATION,
                                         int* ACOUSTIC_SIMULATION) {
  TRACE("compute_kernels_hess_cuda");

  Mesh* mp = (Mesh*)(*Mesh_pointer); //get mesh pointer out of fortran integer container

  int blocksize = NGLL3; // NGLLX*NGLLY*NGLLZ
  realw deltat = *deltat_f;

  int num_blocks_x, num_blocks_y;
  get_blocks_xy(mp->NSPEC_AB,&num_blocks_x,&num_blocks_y);

  dim3 grid(num_blocks_x,num_blocks_y);
  dim3 threads(blocksize,1,1);

  if (*ELASTIC_SIMULATION) {
#ifdef USE_CUDA
    if (run_cuda){
      compute_kernels_hess_el_cudakernel<<<grid,threads>>>(mp->d_ispec_is_elastic,
                                                           mp->d_ibool,
                                                           mp->d_accel,
                                                           mp->d_b_accel,
                                                           mp->d_b_veloc,
                                                           mp->d_b_epsilondev_xx,
                                                           mp->d_b_epsilondev_yy,
                                                           mp->d_b_epsilondev_xy,
                                                           mp->d_b_epsilondev_xz,
                                                           mp->d_b_epsilondev_yz,
                                                           mp->d_b_epsilon_trace_over_3,
                                                           mp->d_hess_el_kl,
                                                           mp->d_hess_rho_el_kl,
                                                           mp->d_hess_kappa_el_kl,
                                                           mp->d_hess_mu_el_kl,
                                                           deltat,
                                                           mp->NSPEC_AB);
    }
#endif
#ifdef USE_HIP
    if (run_hip){
      hipLaunchKernelGGL(compute_kernels_hess_el_cudakernel, dim3(grid), dim3(threads), 0, 0,
                                                             mp->d_ispec_is_elastic,
                                                             mp->d_ibool,
                                                             mp->d_accel,
                                                             mp->d_b_accel,
                                                             mp->d_b_veloc,
                                                             mp->d_b_epsilondev_xx,
                                                             mp->d_b_epsilondev_yy,
                                                             mp->d_b_epsilondev_xy,
                                                             mp->d_b_epsilondev_xz,
                                                             mp->d_b_epsilondev_yz,
                                                             mp->d_b_epsilon_trace_over_3,
                                                             mp->d_hess_el_kl,
                                                             mp->d_hess_rho_el_kl,
                                                             mp->d_hess_kappa_el_kl,
                                                             mp->d_hess_mu_el_kl,
                                                             deltat,
                                                             mp->NSPEC_AB);
    }
#endif
  }

  if (*ACOUSTIC_SIMULATION) {
#ifdef USE_CUDA
    if (run_cuda){
      compute_kernels_hess_ac_cudakernel<<<grid,threads>>>(mp->d_ispec_is_acoustic,
                                                           mp->d_ibool,
                                                           mp->d_potential_dot_dot_acoustic,
                                                           mp->d_b_potential_dot_dot_acoustic,
                                                           mp->d_b_potential_dot_acoustic,
                                                           mp->d_rhostore,
                                                           mp->d_kappastore,
                                                           mp->d_hprime_xx,
                                                           mp->d_irregular_element_number,
                                                           mp->d_xix,mp->d_xiy,mp->d_xiz,
                                                           mp->d_etax,mp->d_etay,mp->d_etaz,
                                                           mp->d_gammax,mp->d_gammay,mp->d_gammaz,
                                                           mp->xix_regular,
                                                           mp->d_hess_ac_kl,
                                                           mp->d_hess_rho_ac_kl,
                                                           mp->d_hess_kappa_ac_kl,
                                                           deltat,
                                                           mp->NSPEC_AB,
                                                           mp->gravity);
    }
#endif
#ifdef USE_HIP
    if (run_hip){
      hipLaunchKernelGGL(compute_kernels_hess_ac_cudakernel, dim3(grid), dim3(threads), 0, 0,
                                                             mp->d_ispec_is_acoustic,
                                                             mp->d_ibool,
                                                             mp->d_potential_dot_dot_acoustic,
                                                             mp->d_b_potential_dot_dot_acoustic,
                                                             mp->d_b_potential_dot_acoustic,
                                                             mp->d_rhostore,
                                                             mp->d_kappastore,
                                                             mp->d_hprime_xx,
                                                             mp->d_irregular_element_number,
                                                             mp->d_xix,mp->d_xiy,mp->d_xiz,
                                                             mp->d_etax,mp->d_etay,mp->d_etaz,
                                                             mp->d_gammax,mp->d_gammay,mp->d_gammaz,
                                                             mp->xix_regular,
                                                             mp->d_hess_ac_kl,
                                                             mp->d_hess_rho_ac_kl,
                                                             mp->d_hess_kappa_ac_kl,
                                                             deltat,
                                                             mp->NSPEC_AB,
                                                             mp->gravity);
    }
#endif
  }

  GPU_ERROR_CHECKING("compute_kernels_hess_cuda");
}

