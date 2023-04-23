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


extern EXTERN_LANG
void FC_FUNC_(compute_strain_cuda,
              COMPUTE_strain_CUDA)(long* Mesh_pointer,
                                   int* FORWARD_OR_ADJOINT) {

  TRACE("compute_strain_cuda");

  Mesh* mp = (Mesh*)(*Mesh_pointer); //get mesh pointer out of fortran integer container

  int size = mp->NSPEC_AB;
  int blocksize = NGLL3; // NGLLX*NGLLY*NGLLZ

  realw *displ;
  realw *epsilondev_xx,*epsilondev_yy,*epsilondev_xy,*epsilondev_xz,*epsilondev_yz;
  realw *epsilondev_trace,*epsilon_trace_over_3;

  // safety check
  if (*FORWARD_OR_ADJOINT != 1 && *FORWARD_OR_ADJOINT != 3) {
    exit_on_error("Error invalid FORWARD_OR_ADJOINT in compute_strain_cuda() routine");
  }
  // only for kernel runs (needs epsilon_trace_over_3 arrays allocated)
  if (mp->simulation_type != 3) exit_on_error("Error invalid SIMULATION_TYPE in compute_strain_cuda() routine");

  // checks flag
  if (! mp->undo_attenuation) exit_on_error("Error invalid UNDO_ATTENUATION flag in compute_strain_cuda() routine");

  int num_blocks_x, num_blocks_y;
  get_blocks_xy(size,&num_blocks_x,&num_blocks_y);

  dim3 grid(num_blocks_x,num_blocks_y);
  dim3 threads(blocksize,1,1);

  // sets gpu arrays
  if (*FORWARD_OR_ADJOINT == 1) {
    // forward fields
    displ = mp->d_displ;
    epsilondev_xx = mp->d_epsilondev_xx;
    epsilondev_yy = mp->d_epsilondev_yy;
    epsilondev_xy = mp->d_epsilondev_xy;
    epsilondev_xz = mp->d_epsilondev_xz;
    epsilondev_yz = mp->d_epsilondev_yz;
    epsilondev_trace = mp->d_epsilondev_trace;
    epsilon_trace_over_3 = mp->d_epsilon_trace_over_3;
  } else {
    // for backward/reconstructed fields
    displ = mp->d_b_displ;
    epsilondev_xx = mp->d_b_epsilondev_xx;
    epsilondev_yy = mp->d_b_epsilondev_yy;
    epsilondev_xy = mp->d_b_epsilondev_xy;
    epsilondev_xz = mp->d_b_epsilondev_xz;
    epsilondev_yz = mp->d_b_epsilondev_yz;
    epsilondev_trace = mp->d_b_epsilondev_trace;
    epsilon_trace_over_3 = mp->d_b_epsilon_trace_over_3;
  }

  // computes strain
#ifdef USE_CUDA
  if (run_cuda){
    compute_element_strain_cudakernel<<<grid,threads>>>(mp->d_ispec_is_elastic,mp->d_ibool,
                                                        displ,
                                                        epsilondev_xx,
                                                        epsilondev_yy,
                                                        epsilondev_xy,
                                                        epsilondev_xz,
                                                        epsilondev_yz,
                                                        epsilondev_trace,
                                                        epsilon_trace_over_3,
                                                        mp->d_xix,mp->d_xiy,mp->d_xiz,
                                                        mp->d_etax,mp->d_etay,mp->d_etaz,
                                                        mp->d_gammax,mp->d_gammay,mp->d_gammaz,
                                                        mp->d_irregular_element_number,
                                                        mp->xix_regular,
                                                        mp->d_hprime_xx,
                                                        size);
  }
#endif
#ifdef USE_HIP
  if (run_hip){
    hipLaunchKernelGGL(compute_element_strain_cudakernel, dim3(grid), dim3(threads), 0, 0,
                                                          mp->d_ispec_is_elastic,mp->d_ibool,
                                                          displ,
                                                          epsilondev_xx,
                                                          epsilondev_yy,
                                                          epsilondev_xy,
                                                          epsilondev_xz,
                                                          epsilondev_yz,
                                                          epsilondev_trace,
                                                          epsilon_trace_over_3,
                                                          mp->d_xix,mp->d_xiy,mp->d_xiz,
                                                          mp->d_etax,mp->d_etay,mp->d_etaz,
                                                          mp->d_gammax,mp->d_gammay,mp->d_gammaz,
                                                          mp->d_irregular_element_number,
                                                          mp->xix_regular,
                                                          mp->d_hprime_xx,
                                                          size);
  }
#endif

  GPU_ERROR_CHECKING("compute_strain_cuda");
}

