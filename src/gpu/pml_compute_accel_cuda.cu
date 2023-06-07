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

/*

extern EXTERN_LANG
void FC_FUNC_(compute_forces_viscoelastic_pml_cuda,
              COMPUTE_FORCES_VISCOELASTIC_PML_CUDA)(long* Mesh_pointer,
                                                    int* iphase,
                                                    realw* deltat,
                                                    int* nspec_outer_elastic,
                                                    int* nspec_inner_elastic,
                                                    int* COMPUTE_AND_STORE_STRAIN,
                                                    int* ATTENUATION,
                                                    int* FORWARD_OR_ADJOINT_f) {

  TRACE("compute_forces_viscoelastic_pml_cuda");

  Mesh* mp = (Mesh*)(*Mesh_pointer); // get Mesh from fortran integer wrapper

  int FORWARD_OR_ADJOINT = *FORWARD_OR_ADJOINT_f;
  // safety check
  if (FORWARD_OR_ADJOINT != 0 && FORWARD_OR_ADJOINT != 1 && FORWARD_OR_ADJOINT != 3) {
    exit_on_error("Error invalid FORWARD_OR_ADJOINT in compute_forces_viscoelastic_pml_cuda() routine");
  }
  // mesh coloring
  if (mp->use_mesh_coloring_gpu ){ exit_on_error("PML kernel w/ mesh coloring not supported yet"); }

  int num_elements;
  if (*iphase == 1)
    num_elements = *nspec_outer_elastic;
  else
    num_elements = *nspec_inner_elastic;

  // checks if anything to do
  if (num_elements == 0) return;

  int blocksize = NGLL3_PADDED;

  int num_blocks_x, num_blocks_y;
  get_blocks_xy(nb_blocks_to_compute,&num_blocks_x,&num_blocks_y);

  dim3 grid(num_blocks_x,num_blocks_y);
  dim3 threads(blocksize,1,1);

  // cuda kernel call
#ifdef USE_CUDA
  if (run_cuda){
    pml_compute_accel_cuda_kernel<<<grid,threads,0,mp->compute_stream>>>( nb_blocks_to_compute,
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
    hipLaunchKernelGGL(pml_compute_accel_cuda_kernel, dim3(grid), dim3(threads), 0, mp->compute_stream,
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

  GPU_ERROR_CHECKING("compute_forces_viscoelastic_pml_cuda");
}

*/

/* ----------------------------------------------------------------------------------------------- */


extern EXTERN_LANG
void FC_FUNC_(pml_impose_boundary_condition_elastic_cuda,
              PML_IMPOSE_BOUNDARY_CONDITION_ELASTIC_CUDA)(long* Mesh_pointer) {

  TRACE("pml_impose_boundary_condition_elastic_cuda");

  Mesh* mp = (Mesh*)(*Mesh_pointer); //get mesh pointer out of fortran integer container

  // checks if anything to do
  if (! mp->pml_conditions) return;

  if (mp->d_num_abs_boundary_faces == 0) return;
  if (mp->NSPEC_CPML == 0) return;

  int num_blocks_x, num_blocks_y;
  get_blocks_xy(mp->d_num_abs_boundary_faces,&num_blocks_x,&num_blocks_y);

  dim3 grid(num_blocks_x,num_blocks_y);
  dim3 threads(NGLL2,1,1);

  // impose Dirichlet conditions on the outer edges of the C-PML layers
#ifdef USE_CUDA
  if (run_cuda){
    pml_impose_boundary_condition_cuda_kernel<<<grid,threads,0,mp->compute_stream>>>(mp->d_accel,mp->d_veloc,mp->d_displ,
                                                                                     mp->d_PML_displ_old,mp->d_PML_displ_new,
                                                                                     mp->d_abs_boundary_ispec,
                                                                                     mp->d_abs_boundary_ijk,
                                                                                     mp->d_num_abs_boundary_faces,
                                                                                     mp->d_ibool,
                                                                                     mp->d_ispec_is_elastic,
                                                                                     mp->d_is_CPML,
                                                                                     mp->d_spec_to_CPML);
  }
#endif
#ifdef USE_HIP
  if (run_hip){
    hipLaunchKernelGGL(pml_impose_boundary_condition_cuda_kernel, dim3(grid), dim3(threads), 0, mp->compute_stream,
                                                                 mp->d_accel,mp->d_veloc,mp->d_displ,
                                                                 mp->d_PML_displ_old,mp->d_PML_displ_new,
                                                                 mp->d_abs_boundary_ispec,
                                                                 mp->d_abs_boundary_ijk,
                                                                 mp->d_num_abs_boundary_faces,
                                                                 mp->d_ibool,
                                                                 mp->d_ispec_is_elastic,
                                                                 mp->d_is_CPML,
                                                                 mp->d_spec_to_CPML);
  }
#endif

  GPU_ERROR_CHECKING("pml_impose_boundary_condition_elastic_cuda");
}
