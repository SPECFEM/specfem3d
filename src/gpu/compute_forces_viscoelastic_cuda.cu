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

// KERNEL 2

/* ----------------------------------------------------------------------------------------------- */

void Kernel_2(int nb_blocks_to_compute,Mesh* mp,int d_iphase,realw d_deltat,
              const int COMPUTE_AND_STORE_STRAIN,
              const int ATTENUATION,
              int* d_ibool,
              realw* d_xix,realw* d_xiy,realw* d_xiz,
              realw* d_etax,realw* d_etay,realw* d_etaz,
              realw* d_gammax,realw* d_gammay,realw* d_gammaz,
              realw* d_kappav,
              realw* d_muv,
              realw* d_epsilondev_xx,realw* d_epsilondev_yy,realw* d_epsilondev_xy,
              realw* d_epsilondev_xz,realw* d_epsilondev_yz,
              realw* d_epsilon_trace_over_3,
              realw* d_factor_common,
              realw* d_R_xx,realw* d_R_yy,realw* d_R_xy,
              realw* d_R_xz,realw* d_R_yz,
              realw* d_factor_common_kappa,
              realw* d_R_trace,realw* d_epsilondev_trace,
              realw* d_b_epsilondev_xx,realw* d_b_epsilondev_yy,realw* d_b_epsilondev_xy,
              realw* d_b_epsilondev_xz,realw* d_b_epsilondev_yz,
              realw* d_b_epsilon_trace_over_3,
              realw* d_b_R_xx,realw* d_b_R_yy,realw* d_b_R_xy,
              realw* d_b_R_xz,realw* d_b_R_yz,
              realw* d_b_R_trace,realw* d_b_epsilondev_trace,
              realw* d_c11store,realw* d_c12store,realw* d_c13store,
              realw* d_c14store,realw* d_c15store,realw* d_c16store,
              realw* d_c22store,realw* d_c23store,realw* d_c24store,
              realw* d_c25store,realw* d_c26store,realw* d_c33store,
              realw* d_c34store,realw* d_c35store,realw* d_c36store,
              realw* d_c44store,realw* d_c45store,realw* d_c46store,
              realw* d_c55store,realw* d_c56store,realw* d_c66store,
              realw* d_rhostore,
              const int FORWARD_OR_ADJOINT){

  TRACE("Kernel_2");

  GPU_ERROR_CHECKING("before kernel Kernel 2");

  // safety check
  if (FORWARD_OR_ADJOINT != 0 && FORWARD_OR_ADJOINT != 1 && FORWARD_OR_ADJOINT != 3) {
    exit_on_error("Error invalid FORWARD_OR_ADJOINT in Kernel_2_acoustic() routine");
  }

  // if the grid can handle the number of blocks, we let it be 1D
  // grid_2_x = nb_elem_color;
  // nb_elem_color is just how many blocks we are computing now

  int blocksize = NGLL3_PADDED;

  int num_blocks_x, num_blocks_y;
  get_blocks_xy(nb_blocks_to_compute,&num_blocks_x,&num_blocks_y);

  dim3 grid(num_blocks_x,num_blocks_y);
  dim3 threads(blocksize,1,1);

  // kernel timing
  gpu_event start,stop;
  if (CUDA_TIMING){ start_timing_gpu(&start,&stop); }

  // defines local parameters for forward/adjoint function calls
  realw *displ,*veloc,*accel;
  realw *epsilondev_xx,*epsilondev_yy,*epsilondev_xy,*epsilondev_xz,*epsilondev_yz;
  realw *epsilondev_trace,*epsilon_trace_over_3;
  realw *R_xx,*R_yy,*R_xy,*R_xz,*R_yz,*R_trace;
  realw *alphaval,*betaval,*gammaval;

  // sets gpu arrays
  if (FORWARD_OR_ADJOINT == 1 || FORWARD_OR_ADJOINT == 0) {
    displ = mp->d_displ;
    veloc = mp->d_veloc;
    accel = mp->d_accel;
    epsilondev_xx = d_epsilondev_xx;
    epsilondev_yy = d_epsilondev_yy;
    epsilondev_xy = d_epsilondev_xy;
    epsilondev_xz = d_epsilondev_xz;
    epsilondev_yz = d_epsilondev_yz;
    epsilondev_trace = d_epsilondev_trace;
    epsilon_trace_over_3 = d_epsilon_trace_over_3;
    R_xx = d_R_xx;
    R_yy = d_R_yy;
    R_xy = d_R_xy;
    R_xz = d_R_xz;
    R_yz = d_R_yz;
    R_trace = d_R_trace;
    alphaval = mp->d_alphaval;
    betaval = mp->d_betaval;
    gammaval = mp->d_gammaval;
  } else {
    // backward/reconstructed wavefields -> FORWARD_OR_ADJOINT == 3
    displ = mp->d_b_displ;
    veloc = mp->d_b_veloc;
    accel = mp->d_b_accel;
    epsilondev_xx = d_b_epsilondev_xx;
    epsilondev_yy = d_b_epsilondev_yy;
    epsilondev_xy = d_b_epsilondev_xy;
    epsilondev_xz = d_b_epsilondev_xz;
    epsilondev_yz = d_b_epsilondev_yz;
    epsilondev_trace = d_b_epsilondev_trace;
    epsilon_trace_over_3 = d_b_epsilon_trace_over_3;
    R_xx = d_b_R_xx;
    R_yy = d_b_R_yy;
    R_xy = d_b_R_xy;
    R_xz = d_b_R_xz;
    R_yz = d_b_R_yz;
    R_trace = d_b_R_trace;
    alphaval = mp->d_b_alphaval;
    betaval = mp->d_b_betaval;
    gammaval = mp->d_b_gammaval;
  }

  // cuda kernel call
  if (ATTENUATION ){
    TRACE("\tKernel_2: Kernel_2_att_impl");
    // compute kernels with attenuation
    // forward wavefields -> FORWARD_OR_ADJOINT == 1
#ifdef USE_CUDA
    if (run_cuda){
      Kernel_2_att_impl<<<grid,threads,0,mp->compute_stream>>>( nb_blocks_to_compute,
                                                                d_ibool,
                                                                mp->d_phase_ispec_inner_elastic,
                                                                mp->num_phase_ispec_elastic,
                                                                d_iphase,
                                                                mp->d_irregular_element_number,
                                                                mp->use_mesh_coloring_gpu,
                                                                d_deltat,
                                                                displ,veloc,accel,
                                                                d_xix, d_xiy, d_xiz,
                                                                d_etax, d_etay, d_etaz,
                                                                d_gammax, d_gammay, d_gammaz,
                                                                mp->xix_regular,mp->jacobian_regular,
                                                                mp->d_hprime_xx,
                                                                mp->d_hprimewgll_xx,
                                                                mp->d_wgllwgll_xy, mp->d_wgllwgll_xz, mp->d_wgllwgll_yz,
                                                                d_kappav, d_muv,
                                                                epsilondev_xx,epsilondev_yy,epsilondev_xy,
                                                                epsilondev_xz,epsilondev_yz,
                                                                epsilondev_trace,
                                                                epsilon_trace_over_3,
                                                                mp->simulation_type,
                                                                mp->NSPEC_AB,
                                                                d_factor_common,
                                                                R_xx,R_yy,R_xy,R_xz,R_yz,
                                                                R_trace,
                                                                d_factor_common_kappa,
                                                                alphaval,betaval,gammaval,
                                                                mp->ANISOTROPY,
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
                                                                mp->d_wgll_cube,
                                                                FORWARD_OR_ADJOINT);
    }
#endif
#ifdef USE_HIP
    if (run_hip){
      hipLaunchKernelGGL(Kernel_2_att_impl, dim3(grid), dim3(threads), 0, mp->compute_stream,
                                            nb_blocks_to_compute,
                                            d_ibool,
                                            mp->d_phase_ispec_inner_elastic,
                                            mp->num_phase_ispec_elastic,
                                            d_iphase,
                                            mp->d_irregular_element_number,
                                            mp->use_mesh_coloring_gpu,
                                            d_deltat,
                                            displ,veloc,accel,
                                            d_xix, d_xiy, d_xiz,
                                            d_etax, d_etay, d_etaz,
                                            d_gammax, d_gammay, d_gammaz,
                                            mp->xix_regular,mp->jacobian_regular,
                                            mp->d_hprime_xx,
                                            mp->d_hprimewgll_xx,
                                            mp->d_wgllwgll_xy, mp->d_wgllwgll_xz, mp->d_wgllwgll_yz,
                                            d_kappav, d_muv,
                                            epsilondev_xx,epsilondev_yy,epsilondev_xy,
                                            epsilondev_xz,epsilondev_yz,
                                            epsilondev_trace,
                                            epsilon_trace_over_3,
                                            mp->simulation_type,
                                            mp->NSPEC_AB,
                                            d_factor_common,
                                            R_xx,R_yy,R_xy,R_xz,R_yz,
                                            R_trace,
                                            d_factor_common_kappa,
                                            alphaval,betaval,gammaval,
                                            mp->ANISOTROPY,
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
                                            mp->d_wgll_cube,
                                            FORWARD_OR_ADJOINT);
    }
#endif

    if (mp->simulation_type == 3 && FORWARD_OR_ADJOINT == 0) {
      // backward/reconstructed wavefields -> FORWARD_OR_ADJOINT == 3
#ifdef USE_CUDA
      if (run_cuda){
        Kernel_2_att_impl<<<grid,threads,0,mp->compute_stream>>>(nb_blocks_to_compute,
                                                                 d_ibool,
                                                                 mp->d_phase_ispec_inner_elastic,
                                                                 mp->num_phase_ispec_elastic,
                                                                 d_iphase,
                                                                 mp->d_irregular_element_number,
                                                                 mp->use_mesh_coloring_gpu,
                                                                 d_deltat,
                                                                 mp->d_b_displ,mp->d_b_veloc,mp->d_b_accel,
                                                                 d_xix, d_xiy, d_xiz,
                                                                 d_etax, d_etay, d_etaz,
                                                                 d_gammax, d_gammay, d_gammaz,
                                                                 mp->xix_regular,mp->jacobian_regular,
                                                                 mp->d_hprime_xx,
                                                                 mp->d_hprimewgll_xx,
                                                                 mp->d_wgllwgll_xy, mp->d_wgllwgll_xz, mp->d_wgllwgll_yz,
                                                                 d_kappav, d_muv,
                                                                 d_b_epsilondev_xx,d_b_epsilondev_yy,d_b_epsilondev_xy,
                                                                 d_b_epsilondev_xz,d_b_epsilondev_yz,
                                                                 d_b_epsilondev_trace,
                                                                 d_b_epsilon_trace_over_3,
                                                                 mp->simulation_type,
                                                                 mp->NSPEC_AB,
                                                                 d_factor_common,
                                                                 d_b_R_xx,d_b_R_yy,d_b_R_xy,d_b_R_xz,d_b_R_yz,
                                                                 d_b_R_trace,
                                                                 d_factor_common_kappa,
                                                                 mp->d_b_alphaval,mp->d_b_betaval,mp->d_b_gammaval,
                                                                 mp->ANISOTROPY,
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
                                                                 mp->d_wgll_cube,
                                                                 3);  // 3 == backward
      }
#endif
#ifdef USE_HIP
      if (run_hip){
        hipLaunchKernelGGL(Kernel_2_att_impl, dim3(grid), dim3(threads), 0, mp->compute_stream,
                                              nb_blocks_to_compute,
                                              d_ibool,
                                              mp->d_phase_ispec_inner_elastic,
                                              mp->num_phase_ispec_elastic,
                                              d_iphase,
                                              mp->d_irregular_element_number,
                                              mp->use_mesh_coloring_gpu,
                                              d_deltat,
                                              mp->d_b_displ,mp->d_b_veloc,mp->d_b_accel,
                                              d_xix, d_xiy, d_xiz,
                                              d_etax, d_etay, d_etaz,
                                              d_gammax, d_gammay, d_gammaz,
                                              mp->xix_regular,mp->jacobian_regular,
                                              mp->d_hprime_xx,
                                              mp->d_hprimewgll_xx,
                                              mp->d_wgllwgll_xy, mp->d_wgllwgll_xz, mp->d_wgllwgll_yz,
                                              d_kappav, d_muv,
                                              d_b_epsilondev_xx,d_b_epsilondev_yy,d_b_epsilondev_xy,
                                              d_b_epsilondev_xz,d_b_epsilondev_yz,
                                              d_b_epsilondev_trace,
                                              d_b_epsilon_trace_over_3,
                                              mp->simulation_type,
                                              mp->NSPEC_AB,
                                              d_factor_common,
                                              d_b_R_xx,d_b_R_yy,d_b_R_xy,d_b_R_xz,d_b_R_yz,
                                              d_b_R_trace,
                                              d_factor_common_kappa,
                                              mp->d_b_alphaval,mp->d_b_betaval,mp->d_b_gammaval,
                                              mp->ANISOTROPY,
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
                                              mp->d_wgll_cube,
                                              3);  // 3 == backward
      }
#endif

    }
  }else{
    // compute kernels without attenuation
    if (mp->ANISOTROPY){
      TRACE("\tKernel_2: Kernel_2_noatt_ani_impl");
      // full anisotropy
      // forward wavefields -> FORWARD_OR_ADJOINT == 1
#ifdef USE_CUDA
      if (run_cuda){
        Kernel_2_noatt_ani_impl<<<grid,threads,0,mp->compute_stream>>>( nb_blocks_to_compute,
                                                                        d_ibool,
                                                                        mp->d_phase_ispec_inner_elastic,mp->num_phase_ispec_elastic,
                                                                        d_iphase,
                                                                        mp->d_irregular_element_number,
                                                                        mp->use_mesh_coloring_gpu,
                                                                        displ,
                                                                        accel,
                                                                        d_xix, d_xiy, d_xiz,
                                                                        d_etax, d_etay, d_etaz,
                                                                        d_gammax, d_gammay, d_gammaz,
                                                                        mp->xix_regular,mp->jacobian_regular,
                                                                        mp->d_hprime_xx,
                                                                        mp->d_hprimewgll_xx,
                                                                        mp->d_wgllwgll_xy, mp->d_wgllwgll_xz, mp->d_wgllwgll_yz,
                                                                        d_kappav, d_muv,
                                                                        COMPUTE_AND_STORE_STRAIN,
                                                                        epsilondev_xx,epsilondev_yy,epsilondev_xy,
                                                                        epsilondev_xz,epsilondev_yz,
                                                                        epsilon_trace_over_3,
                                                                        mp->simulation_type,
                                                                        mp->ANISOTROPY,
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
                                                                        mp->d_wgll_cube,
                                                                        FORWARD_OR_ADJOINT);
      }
#endif
#ifdef USE_HIP
      if (run_hip){
        hipLaunchKernelGGL(Kernel_2_noatt_ani_impl, dim3(grid), dim3(threads), 0, mp->compute_stream,
                                                    nb_blocks_to_compute,
                                                    d_ibool,
                                                    mp->d_phase_ispec_inner_elastic,mp->num_phase_ispec_elastic,
                                                    d_iphase,
                                                    mp->d_irregular_element_number,
                                                    mp->use_mesh_coloring_gpu,
                                                    displ,
                                                    accel,
                                                    d_xix, d_xiy, d_xiz,
                                                    d_etax, d_etay, d_etaz,
                                                    d_gammax, d_gammay, d_gammaz,
                                                    mp->xix_regular,mp->jacobian_regular,
                                                    mp->d_hprime_xx,
                                                    mp->d_hprimewgll_xx,
                                                    mp->d_wgllwgll_xy, mp->d_wgllwgll_xz, mp->d_wgllwgll_yz,
                                                    d_kappav, d_muv,
                                                    COMPUTE_AND_STORE_STRAIN,
                                                    epsilondev_xx,epsilondev_yy,epsilondev_xy,
                                                    epsilondev_xz,epsilondev_yz,
                                                    epsilon_trace_over_3,
                                                    mp->simulation_type,
                                                    mp->ANISOTROPY,
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
                                                    mp->d_wgll_cube,
                                                    FORWARD_OR_ADJOINT);
      }
#endif

      // backward/reconstructed wavefield
      if (mp->simulation_type == 3 && FORWARD_OR_ADJOINT == 0) {
        // backward/reconstructed wavefields -> FORWARD_OR_ADJOINT == 3
#ifdef USE_CUDA
        if (run_cuda){
          Kernel_2_noatt_ani_impl<<< grid,threads,0,mp->compute_stream>>>( nb_blocks_to_compute,
                                                                           d_ibool,
                                                                           mp->d_phase_ispec_inner_elastic,mp->num_phase_ispec_elastic,
                                                                           d_iphase,
                                                                           mp->d_irregular_element_number,
                                                                           mp->use_mesh_coloring_gpu,
                                                                           mp->d_b_displ,
                                                                           mp->d_b_accel,
                                                                           d_xix, d_xiy, d_xiz,
                                                                           d_etax, d_etay, d_etaz,
                                                                           d_gammax, d_gammay, d_gammaz,
                                                                           mp->xix_regular,mp->jacobian_regular,
                                                                           mp->d_hprime_xx,
                                                                           mp->d_hprimewgll_xx,
                                                                           mp->d_wgllwgll_xy, mp->d_wgllwgll_xz, mp->d_wgllwgll_yz,
                                                                           d_kappav, d_muv,
                                                                           COMPUTE_AND_STORE_STRAIN,
                                                                           d_b_epsilondev_xx,d_b_epsilondev_yy,d_b_epsilondev_xy,
                                                                           d_b_epsilondev_xz,d_b_epsilondev_yz,
                                                                           d_b_epsilon_trace_over_3,
                                                                           mp->simulation_type,
                                                                           mp->ANISOTROPY,
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
                                                                           mp->d_wgll_cube,
                                                                           3); // 3 == backward
        }
#endif
#ifdef USE_HIP
        if (run_hip){
          hipLaunchKernelGGL(Kernel_2_noatt_ani_impl, dim3(grid), dim3(threads), 0, mp->compute_stream,
                                                      nb_blocks_to_compute,
                                                      d_ibool,
                                                      mp->d_phase_ispec_inner_elastic,mp->num_phase_ispec_elastic,
                                                      d_iphase,
                                                      mp->d_irregular_element_number,
                                                      mp->use_mesh_coloring_gpu,
                                                      mp->d_b_displ,
                                                      mp->d_b_accel,
                                                      d_xix, d_xiy, d_xiz,
                                                      d_etax, d_etay, d_etaz,
                                                      d_gammax, d_gammay, d_gammaz,
                                                      mp->xix_regular,mp->jacobian_regular,
                                                      mp->d_hprime_xx,
                                                      mp->d_hprimewgll_xx,
                                                      mp->d_wgllwgll_xy, mp->d_wgllwgll_xz, mp->d_wgllwgll_yz,
                                                      d_kappav, d_muv,
                                                      COMPUTE_AND_STORE_STRAIN,
                                                      d_b_epsilondev_xx,d_b_epsilondev_yy,d_b_epsilondev_xy,
                                                      d_b_epsilondev_xz,d_b_epsilondev_yz,
                                                      d_b_epsilon_trace_over_3,
                                                      mp->simulation_type,
                                                      mp->ANISOTROPY,
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
                                                      mp->d_wgll_cube,
                                                      3); // 3 == backward
        }
#endif
      }
    }else{
      // isotropic case
      if (mp->gravity){
        TRACE("\tKernel_2: Kernel_2_noatt_iso_grav_impl");
        // with gravity
        // forward wavefields -> FORWARD_OR_ADJOINT == 1
#ifdef USE_CUDA
        if (run_cuda){
          Kernel_2_noatt_iso_grav_impl<<<grid,threads,0,mp->compute_stream>>>(nb_blocks_to_compute,
                                                                              d_ibool,
                                                                              mp->d_phase_ispec_inner_elastic,mp->num_phase_ispec_elastic,
                                                                              d_iphase,
                                                                              mp->d_irregular_element_number,
                                                                              mp->use_mesh_coloring_gpu,
                                                                              displ,
                                                                              accel,
                                                                              d_xix, d_xiy, d_xiz,
                                                                              d_etax, d_etay, d_etaz,
                                                                              d_gammax, d_gammay, d_gammaz,
                                                                              mp->xix_regular,mp->jacobian_regular,
                                                                              mp->d_hprime_xx,
                                                                              mp->d_hprimewgll_xx,
                                                                              mp->d_wgllwgll_xy, mp->d_wgllwgll_xz, mp->d_wgllwgll_yz,
                                                                              d_kappav, d_muv,
                                                                              COMPUTE_AND_STORE_STRAIN,
                                                                              epsilondev_xx,epsilondev_yy,epsilondev_xy,
                                                                              epsilondev_xz,epsilondev_yz,
                                                                              epsilon_trace_over_3,
                                                                              mp->simulation_type,
                                                                              mp->gravity,
                                                                              mp->d_minus_g,
                                                                              mp->d_minus_deriv_gravity,
                                                                              d_rhostore,
                                                                              mp->d_wgll_cube,
                                                                              FORWARD_OR_ADJOINT);
        }
#endif
#ifdef USE_HIP
        if (run_hip){
          hipLaunchKernelGGL(Kernel_2_noatt_iso_grav_impl, dim3(grid), dim3(threads), 0, mp->compute_stream,
                                                           nb_blocks_to_compute,
                                                           d_ibool,
                                                           mp->d_phase_ispec_inner_elastic,mp->num_phase_ispec_elastic,
                                                           d_iphase,
                                                           mp->d_irregular_element_number,
                                                           mp->use_mesh_coloring_gpu,
                                                           displ,
                                                           accel,
                                                           d_xix, d_xiy, d_xiz,
                                                           d_etax, d_etay, d_etaz,
                                                           d_gammax, d_gammay, d_gammaz,
                                                           mp->xix_regular,mp->jacobian_regular,
                                                           mp->d_hprime_xx,
                                                           mp->d_hprimewgll_xx,
                                                           mp->d_wgllwgll_xy, mp->d_wgllwgll_xz, mp->d_wgllwgll_yz,
                                                           d_kappav, d_muv,
                                                           COMPUTE_AND_STORE_STRAIN,
                                                           epsilondev_xx,epsilondev_yy,epsilondev_xy,
                                                           epsilondev_xz,epsilondev_yz,
                                                           epsilon_trace_over_3,
                                                           mp->simulation_type,
                                                           mp->gravity,
                                                           mp->d_minus_g,
                                                           mp->d_minus_deriv_gravity,
                                                           d_rhostore,
                                                           mp->d_wgll_cube,
                                                           FORWARD_OR_ADJOINT);
        }
#endif

        // backward/reconstructed wavefield
        if (mp->simulation_type == 3 && FORWARD_OR_ADJOINT == 0) {
          // backward/reconstructed wavefields -> FORWARD_OR_ADJOINT == 3
#ifdef USE_CUDA
          if (run_cuda){
            Kernel_2_noatt_iso_grav_impl<<< grid,threads,0,mp->compute_stream>>>(nb_blocks_to_compute,
                                                                                 d_ibool,
                                                                                 mp->d_phase_ispec_inner_elastic,mp->num_phase_ispec_elastic,
                                                                                 d_iphase,
                                                                                 mp->d_irregular_element_number,
                                                                                 mp->use_mesh_coloring_gpu,
                                                                                 mp->d_b_displ,
                                                                                 mp->d_b_accel,
                                                                                 d_xix, d_xiy, d_xiz,
                                                                                 d_etax, d_etay, d_etaz,
                                                                                 d_gammax, d_gammay, d_gammaz,
                                                                                 mp->xix_regular,mp->jacobian_regular,
                                                                                 mp->d_hprime_xx,
                                                                                 mp->d_hprimewgll_xx,
                                                                                 mp->d_wgllwgll_xy, mp->d_wgllwgll_xz, mp->d_wgllwgll_yz,
                                                                                 d_kappav, d_muv,
                                                                                 COMPUTE_AND_STORE_STRAIN,
                                                                                 d_b_epsilondev_xx,d_b_epsilondev_yy,d_b_epsilondev_xy,
                                                                                 d_b_epsilondev_xz,d_b_epsilondev_yz,
                                                                                 d_b_epsilon_trace_over_3,
                                                                                 mp->simulation_type,
                                                                                 mp->gravity,
                                                                                 mp->d_minus_g,
                                                                                 mp->d_minus_deriv_gravity,
                                                                                 d_rhostore,
                                                                                 mp->d_wgll_cube,
                                                                                 3); // 3 == backward
          }
#endif
#ifdef USE_HIP
          if (run_hip){
            hipLaunchKernelGGL(Kernel_2_noatt_iso_grav_impl, dim3(grid), dim3(threads), 0, mp->compute_stream,
                                                             nb_blocks_to_compute,
                                                             d_ibool,
                                                             mp->d_phase_ispec_inner_elastic,mp->num_phase_ispec_elastic,
                                                             d_iphase,
                                                             mp->d_irregular_element_number,
                                                             mp->use_mesh_coloring_gpu,
                                                             mp->d_b_displ,
                                                             mp->d_b_accel,
                                                             d_xix, d_xiy, d_xiz,
                                                             d_etax, d_etay, d_etaz,
                                                             d_gammax, d_gammay, d_gammaz,
                                                             mp->xix_regular,mp->jacobian_regular,
                                                             mp->d_hprime_xx,
                                                             mp->d_hprimewgll_xx,
                                                             mp->d_wgllwgll_xy, mp->d_wgllwgll_xz, mp->d_wgllwgll_yz,
                                                             d_kappav, d_muv,
                                                             COMPUTE_AND_STORE_STRAIN,
                                                             d_b_epsilondev_xx,d_b_epsilondev_yy,d_b_epsilondev_xy,
                                                             d_b_epsilondev_xz,d_b_epsilondev_yz,
                                                             d_b_epsilon_trace_over_3,
                                                             mp->simulation_type,
                                                             mp->gravity,
                                                             mp->d_minus_g,
                                                             mp->d_minus_deriv_gravity,
                                                             d_rhostore,
                                                             mp->d_wgll_cube,
                                                             3); // 3 == backward
          }
#endif

        }
      }else{
        // without gravity
        if (mp->use_mesh_coloring_gpu) {
          TRACE("\tKernel_2: Kernel_2_noatt_iso_col_impl");
          // with mesh coloring
          // forward wavefields -> FORWARD_OR_ADJOINT == 1
#ifdef USE_CUDA
          if (run_cuda){
            Kernel_2_noatt_iso_col_impl<<<grid,threads,0,mp->compute_stream>>>( nb_blocks_to_compute,
                                                                                d_ibool,
                                                                                mp->d_phase_ispec_inner_elastic,mp->num_phase_ispec_elastic,
                                                                                d_iphase,
                                                                                mp->use_mesh_coloring_gpu,
                                                                                displ,
                                                                                accel,
                                                                                d_xix, d_xiy, d_xiz,
                                                                                d_etax, d_etay, d_etaz,
                                                                                d_gammax, d_gammay, d_gammaz,
                                                                                mp->d_hprime_xx,
                                                                                mp->d_hprimewgll_xx,
                                                                                mp->d_wgllwgll_xy, mp->d_wgllwgll_xz, mp->d_wgllwgll_yz,
                                                                                d_kappav, d_muv,
                                                                                COMPUTE_AND_STORE_STRAIN,
                                                                                epsilondev_xx,epsilondev_yy,epsilondev_xy,
                                                                                epsilondev_xz,epsilondev_yz,
                                                                                epsilon_trace_over_3,
                                                                                mp->simulation_type,
                                                                                FORWARD_OR_ADJOINT);
          }
#endif
#ifdef USE_HIP
          if (run_hip){
            hipLaunchKernelGGL(Kernel_2_noatt_iso_col_impl, dim3(grid), dim3(threads), 0, mp->compute_stream,
                                                            nb_blocks_to_compute,
                                                            d_ibool,
                                                            mp->d_phase_ispec_inner_elastic,mp->num_phase_ispec_elastic,
                                                            d_iphase,
                                                            mp->use_mesh_coloring_gpu,
                                                            displ,
                                                            accel,
                                                            d_xix, d_xiy, d_xiz,
                                                            d_etax, d_etay, d_etaz,
                                                            d_gammax, d_gammay, d_gammaz,
                                                            mp->d_hprime_xx,
                                                            mp->d_hprimewgll_xx,
                                                            mp->d_wgllwgll_xy, mp->d_wgllwgll_xz, mp->d_wgllwgll_yz,
                                                            d_kappav, d_muv,
                                                            COMPUTE_AND_STORE_STRAIN,
                                                            epsilondev_xx,epsilondev_yy,epsilondev_xy,
                                                            epsilondev_xz,epsilondev_yz,
                                                            epsilon_trace_over_3,
                                                            mp->simulation_type,
                                                            FORWARD_OR_ADJOINT);
          }
#endif

          // backward/reconstructed wavefield
          if (mp->simulation_type == 3 && FORWARD_OR_ADJOINT == 0) {
            // backward/reconstructed wavefields -> FORWARD_OR_ADJOINT == 3
#ifdef USE_CUDA
            if (run_cuda){
              Kernel_2_noatt_iso_col_impl<<< grid,threads,0,mp->compute_stream>>>( nb_blocks_to_compute,
                                                                                   d_ibool,
                                                                                   mp->d_phase_ispec_inner_elastic,mp->num_phase_ispec_elastic,
                                                                                   d_iphase,
                                                                                   mp->use_mesh_coloring_gpu,
                                                                                   mp->d_b_displ,
                                                                                   mp->d_b_accel,
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
                                                                                   3); // 3 == backward
            }
#endif
#ifdef USE_HIP
            if (run_hip){
              hipLaunchKernelGGL(Kernel_2_noatt_iso_col_impl, dim3(grid), dim3(threads), 0, mp->compute_stream,
                                                              nb_blocks_to_compute,
                                                              d_ibool,
                                                              mp->d_phase_ispec_inner_elastic,mp->num_phase_ispec_elastic,
                                                              d_iphase,
                                                              mp->use_mesh_coloring_gpu,
                                                              mp->d_b_displ,
                                                              mp->d_b_accel,
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
                                                              3); // 3 == backward
            }
#endif
          }
        }else{
          // without mesh coloring
          if (COMPUTE_AND_STORE_STRAIN){
            TRACE("\tKernel_2: Kernel_2_noatt_iso_strain_impl");
            // stores strains
            // forward wavefields -> FORWARD_OR_ADJOINT == 1
#ifdef USE_CUDA
            if (run_cuda){
              Kernel_2_noatt_iso_strain_impl<<<grid,threads,0,mp->compute_stream>>>(nb_blocks_to_compute,
                                                                                    d_ibool,
                                                                                    mp->d_phase_ispec_inner_elastic,
                                                                                    mp->num_phase_ispec_elastic,
                                                                                    d_iphase,
                                                                                    mp->d_irregular_element_number,
                                                                                    displ,
                                                                                    accel,
                                                                                    d_xix, d_xiy, d_xiz,
                                                                                    d_etax, d_etay, d_etaz,
                                                                                    d_gammax, d_gammay, d_gammaz,
                                                                                    mp->xix_regular,mp->jacobian_regular,
                                                                                    mp->d_hprime_xx,
                                                                                    mp->d_hprimewgll_xx,
                                                                                    mp->d_wgllwgll_xy, mp->d_wgllwgll_xz, mp->d_wgllwgll_yz,
                                                                                    d_kappav, d_muv,
                                                                                    COMPUTE_AND_STORE_STRAIN,
                                                                                    epsilondev_xx,epsilondev_yy,epsilondev_xy,
                                                                                    epsilondev_xz,epsilondev_yz,
                                                                                    epsilon_trace_over_3,
                                                                                    mp->simulation_type,
                                                                                    FORWARD_OR_ADJOINT);
            }
#endif
#ifdef USE_HIP
            if (run_hip){
              hipLaunchKernelGGL(Kernel_2_noatt_iso_strain_impl, dim3(grid), dim3(threads), 0, mp->compute_stream,
                                                                 nb_blocks_to_compute,
                                                                 d_ibool,
                                                                 mp->d_phase_ispec_inner_elastic,
                                                                 mp->num_phase_ispec_elastic,
                                                                 d_iphase,
                                                                 mp->d_irregular_element_number,
                                                                 displ,
                                                                 accel,
                                                                 d_xix, d_xiy, d_xiz,
                                                                 d_etax, d_etay, d_etaz,
                                                                 d_gammax, d_gammay, d_gammaz,
                                                                 mp->xix_regular,mp->jacobian_regular,
                                                                 mp->d_hprime_xx,
                                                                 mp->d_hprimewgll_xx,
                                                                 mp->d_wgllwgll_xy, mp->d_wgllwgll_xz, mp->d_wgllwgll_yz,
                                                                 d_kappav, d_muv,
                                                                 COMPUTE_AND_STORE_STRAIN,
                                                                 epsilondev_xx,epsilondev_yy,epsilondev_xy,
                                                                 epsilondev_xz,epsilondev_yz,
                                                                 epsilon_trace_over_3,
                                                                 mp->simulation_type,
                                                                 FORWARD_OR_ADJOINT);
            }
#endif
            // backward/reconstructed wavefield
            if (mp->simulation_type == 3 && FORWARD_OR_ADJOINT == 0) {
              // backward/reconstructed wavefields -> FORWARD_OR_ADJOINT == 3
#ifdef USE_CUDA
              if (run_cuda){
                Kernel_2_noatt_iso_strain_impl<<< grid,threads,0,mp->compute_stream>>>(nb_blocks_to_compute,
                                                                                       d_ibool,
                                                                                       mp->d_phase_ispec_inner_elastic,
                                                                                       mp->num_phase_ispec_elastic,
                                                                                       d_iphase,
                                                                                       mp->d_irregular_element_number,
                                                                                       mp->d_b_displ,
                                                                                       mp->d_b_accel,
                                                                                       d_xix, d_xiy, d_xiz,
                                                                                       d_etax, d_etay, d_etaz,
                                                                                       d_gammax, d_gammay, d_gammaz,
                                                                                       mp->xix_regular,mp->jacobian_regular,
                                                                                       mp->d_hprime_xx,
                                                                                       mp->d_hprimewgll_xx,
                                                                                       mp->d_wgllwgll_xy, mp->d_wgllwgll_xz, mp->d_wgllwgll_yz,
                                                                                       d_kappav, d_muv,
                                                                                       COMPUTE_AND_STORE_STRAIN,
                                                                                       d_b_epsilondev_xx,d_b_epsilondev_yy,d_b_epsilondev_xy,
                                                                                       d_b_epsilondev_xz,d_b_epsilondev_yz,
                                                                                       d_b_epsilon_trace_over_3,
                                                                                       mp->simulation_type,
                                                                                       3); // 3 == backward
              }
#endif
#ifdef USE_HIP
              if (run_hip){
                hipLaunchKernelGGL(Kernel_2_noatt_iso_strain_impl, dim3(grid), dim3(threads), 0, mp->compute_stream,
                                                                   nb_blocks_to_compute,
                                                                   d_ibool,
                                                                   mp->d_phase_ispec_inner_elastic,
                                                                   mp->num_phase_ispec_elastic,
                                                                   d_iphase,
                                                                   mp->d_irregular_element_number,
                                                                   mp->d_b_displ,
                                                                   mp->d_b_accel,
                                                                   d_xix, d_xiy, d_xiz,
                                                                   d_etax, d_etay, d_etaz,
                                                                   d_gammax, d_gammay, d_gammaz,
                                                                   mp->xix_regular,mp->jacobian_regular,
                                                                   mp->d_hprime_xx,
                                                                   mp->d_hprimewgll_xx,
                                                                   mp->d_wgllwgll_xy, mp->d_wgllwgll_xz, mp->d_wgllwgll_yz,
                                                                   d_kappav, d_muv,
                                                                   COMPUTE_AND_STORE_STRAIN,
                                                                   d_b_epsilondev_xx,d_b_epsilondev_yy,d_b_epsilondev_xy,
                                                                   d_b_epsilondev_xz,d_b_epsilondev_yz,
                                                                   d_b_epsilon_trace_over_3,
                                                                   mp->simulation_type,
                                                                   3); // 3 == backward
              }
#endif

            }
          }else{
            // dynamic rupture simulation
            if (mp->use_Kelvin_Voigt_damping) {
              TRACE("\tKernel_2: Kernel_2_noatt_iso_kelvinvoigt_impl");
              // use_Kelvin_Voigt_damping == true means there is fault in this partition
#ifdef USE_CUDA
              if (run_cuda){
                Kernel_2_noatt_iso_kelvinvoigt_impl<<<grid,threads,0,mp->compute_stream>>>(nb_blocks_to_compute,
                                                                                           d_ibool,
                                                                                           mp->d_phase_ispec_inner_elastic,
                                                                                           mp->num_phase_ispec_elastic,
                                                                                           d_iphase,
                                                                                           mp->d_irregular_element_number,
                                                                                           mp->d_Kelvin_Voigt_eta,
                                                                                           displ,
                                                                                           veloc,
                                                                                           accel,
                                                                                           d_xix, d_xiy, d_xiz,
                                                                                           d_etax, d_etay, d_etaz,
                                                                                           d_gammax, d_gammay, d_gammaz,
                                                                                           mp->xix_regular,mp->jacobian_regular,
                                                                                           mp->d_hprime_xx,
                                                                                           mp->d_hprimewgll_xx,
                                                                                           mp->d_wgllwgll_xy, mp->d_wgllwgll_xz, mp->d_wgllwgll_yz,
                                                                                           d_kappav, d_muv,
                                                                                           FORWARD_OR_ADJOINT);
              }
#endif
#ifdef USE_HIP
              if (run_hip){
                hipLaunchKernelGGL(Kernel_2_noatt_iso_kelvinvoigt_impl, dim3(grid), dim3(threads), 0, mp->compute_stream,
                                                                        nb_blocks_to_compute,
                                                                        d_ibool,
                                                                        mp->d_phase_ispec_inner_elastic,
                                                                        mp->num_phase_ispec_elastic,
                                                                        d_iphase,
                                                                        mp->d_irregular_element_number,
                                                                        mp->d_Kelvin_Voigt_eta,
                                                                        displ,
                                                                        veloc,
                                                                        accel,
                                                                        d_xix, d_xiy, d_xiz,
                                                                        d_etax, d_etay, d_etaz,
                                                                        d_gammax, d_gammay, d_gammaz,
                                                                        mp->xix_regular,mp->jacobian_regular,
                                                                        mp->d_hprime_xx,
                                                                        mp->d_hprimewgll_xx,
                                                                        mp->d_wgllwgll_xy, mp->d_wgllwgll_xz, mp->d_wgllwgll_yz,
                                                                        d_kappav, d_muv,
                                                                        FORWARD_OR_ADJOINT);
              }
#endif
            }else{
              TRACE("\tKernel_2: Kernel_2_noatt_iso_impl");
              // without storing strains
              // forward wavefields -> FORWARD_OR_ADJOINT == 1
#ifdef USE_CUDA
              if (run_cuda){
                Kernel_2_noatt_iso_impl<<<grid,threads,0,mp->compute_stream>>>(nb_blocks_to_compute,
                                                                               d_ibool,
                                                                               mp->d_phase_ispec_inner_elastic,mp->num_phase_ispec_elastic,
                                                                               d_iphase,
                                                                               mp->d_irregular_element_number,
                                                                               displ,
                                                                               accel,
                                                                               d_xix, d_xiy, d_xiz,
                                                                               d_etax, d_etay, d_etaz,
                                                                               d_gammax, d_gammay, d_gammaz,
                                                                               mp->xix_regular,mp->jacobian_regular,
                                                                               mp->d_hprime_xx,
                                                                               mp->d_hprimewgll_xx,
                                                                               mp->d_wgllwgll_xy, mp->d_wgllwgll_xz, mp->d_wgllwgll_yz,
                                                                               d_kappav, d_muv,
                                                                               FORWARD_OR_ADJOINT);
              }
#endif
#ifdef USE_HIP
              if (run_hip){
                 hipLaunchKernelGGL(Kernel_2_noatt_iso_impl, dim3(grid), dim3(threads), 0, mp->compute_stream,
                                                             nb_blocks_to_compute,
                                                             d_ibool,
                                                             mp->d_phase_ispec_inner_elastic,mp->num_phase_ispec_elastic,
                                                             d_iphase,
                                                             mp->d_irregular_element_number,
                                                             displ,
                                                             accel,
                                                             d_xix, d_xiy, d_xiz,
                                                             d_etax, d_etay, d_etaz,
                                                             d_gammax, d_gammay, d_gammaz,
                                                             mp->xix_regular,mp->jacobian_regular,
                                                             mp->d_hprime_xx,
                                                             mp->d_hprimewgll_xx,
                                                             mp->d_wgllwgll_xy, mp->d_wgllwgll_xz, mp->d_wgllwgll_yz,
                                                             d_kappav, d_muv,
                                                             FORWARD_OR_ADJOINT);
              }
#endif
              // backward/reconstructed wavefield
              if (mp->simulation_type == 3 && FORWARD_OR_ADJOINT == 0) {
                // backward/reconstructed wavefields -> FORWARD_OR_ADJOINT == 3
#ifdef USE_CUDA
                if (run_cuda){
                  Kernel_2_noatt_iso_impl<<< grid,threads,0,mp->compute_stream>>>(nb_blocks_to_compute,
                                                                                  d_ibool,
                                                                                  mp->d_phase_ispec_inner_elastic,
                                                                                  mp->num_phase_ispec_elastic,
                                                                                  d_iphase,
                                                                                  mp->d_irregular_element_number,
                                                                                  mp->d_b_displ,
                                                                                  mp->d_b_accel,
                                                                                  d_xix, d_xiy, d_xiz,
                                                                                  d_etax, d_etay, d_etaz,
                                                                                  d_gammax, d_gammay, d_gammaz,
                                                                                  mp->xix_regular,mp->jacobian_regular,
                                                                                  mp->d_hprime_xx,
                                                                                  mp->d_hprimewgll_xx,
                                                                                  mp->d_wgllwgll_xy, mp->d_wgllwgll_xz, mp->d_wgllwgll_yz,
                                                                                  d_kappav, d_muv,
                                                                                  3); // 3 == backward
                }
#endif
#ifdef USE_HIP
                if (run_hip){
                  hipLaunchKernelGGL(Kernel_2_noatt_iso_impl, dim3(grid), dim3(threads), 0, mp->compute_stream,
                                                              nb_blocks_to_compute,
                                                              d_ibool,
                                                              mp->d_phase_ispec_inner_elastic,
                                                              mp->num_phase_ispec_elastic,
                                                              d_iphase,
                                                              mp->d_irregular_element_number,
                                                              mp->d_b_displ,
                                                              mp->d_b_accel,
                                                              d_xix, d_xiy, d_xiz,
                                                              d_etax, d_etay, d_etaz,
                                                              d_gammax, d_gammay, d_gammaz,
                                                              mp->xix_regular,mp->jacobian_regular,
                                                              mp->d_hprime_xx,
                                                              mp->d_hprimewgll_xx,
                                                              mp->d_wgllwgll_xy, mp->d_wgllwgll_xz, mp->d_wgllwgll_yz,
                                                              d_kappav, d_muv,
                                                              3); // 3 == backward
                }
#endif
              }
            }
          } // COMPUTE_AND_STORE_STRAIN
        } // use_mesh_coloring_gpu
      } // gravity
    } // ANISOTROPY
  } // ATTENUATION

  // Cuda timing
  if (CUDA_TIMING){
    if (ATTENUATION){
      stop_timing_gpu(&start,&stop,"Kernel_2_att_impl");
    }else{
      if (mp->ANISOTROPY){
        stop_timing_gpu(&start,&stop,"Kernel_2_noatt_ani_impl");
      }else{
        if (mp->gravity){
          stop_timing_gpu(&start,&stop,"Kernel_2_noatt_iso_grav_impl");
        }else{
          if (COMPUTE_AND_STORE_STRAIN){
            stop_timing_gpu(&start,&stop,"Kernel_2_noatt_iso_strain_impl");
          }else{
            realw time;
            stop_timing_gpu(&start,&stop,"Kernel_2_noatt_iso_impl",&time);
            // time in seconds
            time = time / 1000.;
            // performance
            // see with: nvprof --metrics flops_sp ./xspecfem3D -> using 883146240 FLOPS (Single) floating-point operations
            // hand-counts: 89344 * number-of-blocks
            realw flops = 89344 * nb_blocks_to_compute;
            printf("  performance: %f GFlops/s\n", flops/time *(1./1000./1000./1000.));
          }
        }
      }
    }
  }

  GPU_ERROR_CHECKING("Kernel_2_impl");
}

/* ----------------------------------------------------------------------------------------------- */


extern EXTERN_LANG
void FC_FUNC_(compute_forces_viscoelastic_cuda,
              COMPUTE_FORCES_VISCOELASTIC_CUDA)(long* Mesh_pointer,
                                                int* iphase,
                                                realw* deltat,
                                                int* nspec_outer_elastic,
                                                int* nspec_inner_elastic,
                                                int* COMPUTE_AND_STORE_STRAIN,
                                                int* ATTENUATION,
                                                int* FORWARD_OR_ADJOINT_f) {

  TRACE("compute_forces_viscoelastic_cuda");
  // EPIK_TRACER("compute_forces_viscoelastic_cuda");
  //printf("Running compute_forces\n");
  //double start_time = get_time_val();

  Mesh* mp = (Mesh*)(*Mesh_pointer); // get Mesh from fortran integer wrapper
  int FORWARD_OR_ADJOINT = *FORWARD_OR_ADJOINT_f;

  int num_elements;

  if (*iphase == 1)
    num_elements = *nspec_outer_elastic;
  else
    num_elements = *nspec_inner_elastic;

  // checks if anything to do
  if (num_elements == 0) return;

  // mesh coloring
  if (mp->use_mesh_coloring_gpu ){
    // note: array offsets require sorted arrays, such that e.g. ibool starts with elastic elements
    //         and followed by acoustic ones.
    //         elastic elements also start with outer than inner element ordering
    int nb_colors,nb_blocks_to_compute;
    int istart;
    int offset,offset_nonpadded,offset_nonpadded_att2;

    // sets up color loop
    if (*iphase == 1){
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
      //if (nb_blocks_to_compute <= 0){
      //  printf("error number of elastic color blocks: %d -- color = %d \n",nb_blocks_to_compute,icolor);
      //  exit(EXIT_FAILURE);
      //}

      Kernel_2(nb_blocks_to_compute,mp,*iphase,*deltat,
               *COMPUTE_AND_STORE_STRAIN,
               *ATTENUATION,
               mp->d_ibool + offset,
               mp->d_xix + offset,mp->d_xiy + offset,mp->d_xiz + offset,
               mp->d_etax + offset,mp->d_etay + offset,mp->d_etaz + offset,
               mp->d_gammax + offset,mp->d_gammay + offset,mp->d_gammaz + offset,
               mp->d_kappav + offset,
               mp->d_muv + offset,
               mp->d_epsilondev_xx + offset_nonpadded,mp->d_epsilondev_yy + offset_nonpadded,mp->d_epsilondev_xy + offset_nonpadded,
               mp->d_epsilondev_xz + offset_nonpadded,mp->d_epsilondev_yz + offset_nonpadded,
               mp->d_epsilon_trace_over_3 + offset_nonpadded,
               mp->d_factor_common + offset_nonpadded_att2,
               mp->d_R_xx + offset_nonpadded,mp->d_R_yy + offset_nonpadded,mp->d_R_xy + offset_nonpadded,
               mp->d_R_xz + offset_nonpadded,mp->d_R_yz + offset_nonpadded,
               mp->d_factor_common_kappa + offset_nonpadded_att2,
               mp->d_R_trace + offset_nonpadded,
               mp->d_epsilondev_trace + offset_nonpadded,
               mp->d_b_epsilondev_xx + offset_nonpadded,mp->d_b_epsilondev_yy + offset_nonpadded,mp->d_b_epsilondev_xy + offset_nonpadded,
               mp->d_b_epsilondev_xz + offset_nonpadded,mp->d_b_epsilondev_yz + offset_nonpadded,
               mp->d_b_epsilon_trace_over_3 + offset_nonpadded,
               mp->d_b_R_xx + offset_nonpadded,mp->d_b_R_yy + offset_nonpadded,mp->d_b_R_xy + offset_nonpadded,
               mp->d_b_R_xz + offset_nonpadded,mp->d_b_R_yz + offset_nonpadded,
               mp->d_b_R_trace + offset_nonpadded,
               mp->d_b_epsilondev_trace + offset_nonpadded,
               mp->d_c11store + offset,mp->d_c12store + offset,mp->d_c13store + offset,
               mp->d_c14store + offset,mp->d_c15store + offset,mp->d_c16store + offset,
               mp->d_c22store + offset,mp->d_c23store + offset,mp->d_c24store + offset,
               mp->d_c25store + offset,mp->d_c26store + offset,mp->d_c33store + offset,
               mp->d_c34store + offset,mp->d_c35store + offset,mp->d_c36store + offset,
               mp->d_c44store + offset,mp->d_c45store + offset,mp->d_c46store + offset,
               mp->d_c55store + offset,mp->d_c56store + offset,mp->d_c66store + offset,
               mp->d_rhostore + offset,
               FORWARD_OR_ADJOINT);

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
             *ATTENUATION,
             mp->d_ibool,
             mp->d_xix,mp->d_xiy,mp->d_xiz,
             mp->d_etax,mp->d_etay,mp->d_etaz,
             mp->d_gammax,mp->d_gammay,mp->d_gammaz,
             mp->d_kappav,
             mp->d_muv,
             mp->d_epsilondev_xx,mp->d_epsilondev_yy,mp->d_epsilondev_xy,
             mp->d_epsilondev_xz,mp->d_epsilondev_yz,
             mp->d_epsilon_trace_over_3,
             mp->d_factor_common,
             mp->d_R_xx,mp->d_R_yy,mp->d_R_xy,
             mp->d_R_xz,mp->d_R_yz,
             mp->d_factor_common_kappa,
             mp->d_R_trace,mp->d_epsilondev_trace,
             mp->d_b_epsilondev_xx,mp->d_b_epsilondev_yy,mp->d_b_epsilondev_xy,
             mp->d_b_epsilondev_xz,mp->d_b_epsilondev_yz,
             mp->d_b_epsilon_trace_over_3,
             mp->d_b_R_xx,mp->d_b_R_yy,mp->d_b_R_xy,
             mp->d_b_R_xz,mp->d_b_R_yz,
             mp->d_b_R_trace,mp->d_b_epsilondev_trace,
             mp->d_c11store,mp->d_c12store,mp->d_c13store,
             mp->d_c14store,mp->d_c15store,mp->d_c16store,
             mp->d_c22store,mp->d_c23store,mp->d_c24store,
             mp->d_c25store,mp->d_c26store,mp->d_c33store,
             mp->d_c34store,mp->d_c35store,mp->d_c36store,
             mp->d_c44store,mp->d_c45store,mp->d_c46store,
             mp->d_c55store,mp->d_c56store,mp->d_c66store,
             mp->d_rhostore,
             FORWARD_OR_ADJOINT);
  }
}

