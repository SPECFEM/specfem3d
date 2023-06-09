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

// PML compute forces contribution

/* ----------------------------------------------------------------------------------------------- */

void compute_forces_viscoelastic_pml_cuda(long* Mesh_pointer,
                                          int* iphase,
                                          int* nspec_outer_elastic,
                                          int* nspec_inner_elastic,
                                          int* COMPUTE_AND_STORE_STRAIN_f,
                                          int* FORWARD_OR_ADJOINT_f) {

  TRACE("compute_forces_viscoelastic_pml_cuda");

  Mesh* mp = (Mesh*)(*Mesh_pointer); // get Mesh from fortran integer wrapper

  // needs PML conditions set to .true.
  if (! mp->pml_conditions) return;
  if (mp->NSPEC_CPML == 0) return;

  int FORWARD_OR_ADJOINT = *FORWARD_OR_ADJOINT_f;
  int COMPUTE_AND_STORE_STRAIN = *COMPUTE_AND_STORE_STRAIN_f;

  // safety check
  if (FORWARD_OR_ADJOINT != 0 && FORWARD_OR_ADJOINT != 1 && FORWARD_OR_ADJOINT != 3) {
    exit_on_error("Invalid FORWARD_OR_ADJOINT in PML compute_forces_viscoelastic_pml_cuda() routine");
  }
  // only for forward wavefields
  if (FORWARD_OR_ADJOINT == 3) return;

  // mesh coloring
  if (mp->use_mesh_coloring_gpu){ exit_on_error("PML compute_forces_viscoelastic_pml_cuda() routine w/ mesh coloring not supported yet"); }

  int num_elements;
  if (*iphase == 1)
    num_elements = *nspec_outer_elastic;
  else
    num_elements = *nspec_inner_elastic;

  // checks if anything to do
  if (num_elements == 0) return;

  int blocksize = NGLL3_PADDED;

  int num_blocks_x, num_blocks_y;
  get_blocks_xy(num_elements,&num_blocks_x,&num_blocks_y);

  dim3 grid(num_blocks_x,num_blocks_y);
  dim3 threads(blocksize,1,1);

  // defines local parameters for forward/adjoint function calls
  realw *displ,*veloc,*accel;
  realw *epsilondev_xx,*epsilondev_yy,*epsilondev_xy,*epsilondev_xz,*epsilondev_yz;
  realw *epsilon_trace_over_3;

  int d_iphase = *iphase;

  // sets gpu arrays
  if (FORWARD_OR_ADJOINT == 1 || FORWARD_OR_ADJOINT == 0) {
    displ = mp->d_displ;
    veloc = mp->d_veloc;
    accel = mp->d_accel;
    epsilondev_xx = mp->d_epsilondev_xx;
    epsilondev_yy = mp->d_epsilondev_yy;
    epsilondev_xy = mp->d_epsilondev_xy;
    epsilondev_xz = mp->d_epsilondev_xz;
    epsilondev_yz = mp->d_epsilondev_yz;
    epsilon_trace_over_3 = mp->d_epsilon_trace_over_3;
  } else {
    // backward/reconstructed wavefields -> FORWARD_OR_ADJOINT == 3
    // safety stop
    exit_on_error("PML compute_forces_viscoelastic_pml_cuda() routine for backward wavefields not supported yet");
    // would need also b_PML_displ_new,.. arrays
    displ = mp->d_b_displ;
    veloc = mp->d_b_veloc;
    accel = mp->d_b_accel;
    epsilondev_xx = mp->d_b_epsilondev_xx;
    epsilondev_yy = mp->d_b_epsilondev_yy;
    epsilondev_xy = mp->d_b_epsilondev_xy;
    epsilondev_xz = mp->d_b_epsilondev_xz;
    epsilondev_yz = mp->d_b_epsilondev_yz;
    epsilon_trace_over_3 = mp->d_b_epsilon_trace_over_3;
  }

  // cuda kernel call
#ifdef USE_CUDA
  if (run_cuda){
    pml_kernel_2_impl<<<grid,threads,0,mp->compute_stream>>>(num_elements,
                                                             mp->d_ibool,
                                                             mp->d_phase_ispec_inner_elastic,
                                                             mp->num_phase_ispec_elastic,
                                                             d_iphase,
                                                             mp->d_irregular_element_number,
                                                             displ,veloc,accel,
                                                             mp->d_xix, mp->d_xiy, mp->d_xiz,
                                                             mp->d_etax, mp->d_etay, mp->d_etaz,
                                                             mp->d_gammax, mp->d_gammay, mp->d_gammaz,
                                                             mp->xix_regular,mp->jacobian_regular,
                                                             mp->d_hprime_xx,
                                                             mp->d_hprimewgll_xx,
                                                             mp->d_wgllwgll_xy, mp->d_wgllwgll_xz, mp->d_wgllwgll_yz,
                                                             mp->d_kappav, mp->d_muv,
                                                             COMPUTE_AND_STORE_STRAIN,
                                                             epsilondev_xx,epsilondev_yy,epsilondev_xy,
                                                             epsilondev_xz,epsilondev_yz,
                                                             epsilon_trace_over_3,
                                                             mp->simulation_type,
                                                             mp->d_rhostore,
                                                             mp->d_wgll_cube,
                                                             mp->NSPEC_CPML,
                                                             mp->d_is_CPML,
                                                             mp->d_spec_to_CPML,
                                                             mp->d_PML_displ_new,
                                                             mp->d_PML_displ_old,
                                                             mp->d_rmemory_displ_elastic,
                                                             mp->d_rmemory_dux_dxl_x,
                                                             mp->d_rmemory_duy_dxl_y,
                                                             mp->d_rmemory_duz_dxl_z,
                                                             mp->d_rmemory_dux_dyl_x,
                                                             mp->d_rmemory_duy_dyl_y,
                                                             mp->d_rmemory_duz_dyl_z,
                                                             mp->d_rmemory_dux_dzl_x,
                                                             mp->d_rmemory_duy_dzl_y,
                                                             mp->d_rmemory_duz_dzl_z,
                                                             mp->d_rmemory_dux_dxl_y,
                                                             mp->d_rmemory_dux_dxl_z,
                                                             mp->d_rmemory_duy_dxl_x,
                                                             mp->d_rmemory_duz_dxl_x,
                                                             mp->d_rmemory_dux_dyl_y,
                                                             mp->d_rmemory_duy_dyl_x,
                                                             mp->d_rmemory_duy_dyl_z,
                                                             mp->d_rmemory_duz_dyl_y,
                                                             mp->d_rmemory_dux_dzl_z,
                                                             mp->d_rmemory_duy_dzl_z,
                                                             mp->d_rmemory_duz_dzl_x,
                                                             mp->d_rmemory_duz_dzl_y,
                                                             mp->d_pml_convolution_coef_alpha,
                                                             mp->d_pml_convolution_coef_beta,
                                                             mp->d_pml_convolution_coef_strain,
                                                             mp->d_pml_convolution_coef_abar);
  }
#endif
#ifdef USE_HIP
  if (run_hip){
    hipLaunchKernelGGL(pml_kernel_2_impl, dim3(grid), dim3(threads), 0, mp->compute_stream,
                                                            num_elements,
                                                            mp->d_ibool,
                                                            mp->d_phase_ispec_inner_elastic,
                                                            mp->num_phase_ispec_elastic,
                                                            d_iphase,
                                                            mp->d_irregular_element_number,
                                                            displ,veloc,accel,
                                                            mp->d_xix, mp->d_xiy, mp->d_xiz,
                                                            mp->d_etax, mp->d_etay, mp->d_etaz,
                                                            mp->d_gammax, mp->d_gammay, mp->d_gammaz,
                                                            mp->xix_regular,mp->jacobian_regular,
                                                            mp->d_hprime_xx,
                                                            mp->d_hprimewgll_xx,
                                                            mp->d_wgllwgll_xy, mp->d_wgllwgll_xz, mp->d_wgllwgll_yz,
                                                            mp->d_kappav, mp->d_muv,
                                                            COMPUTE_AND_STORE_STRAIN,
                                                            epsilondev_xx,epsilondev_yy,epsilondev_xy,
                                                            epsilondev_xz,epsilondev_yz,
                                                            epsilon_trace_over_3,
                                                            mp->simulation_type,
                                                            mp->d_rhostore,
                                                            mp->d_wgll_cube,
                                                            mp->NSPEC_CPML,
                                                            mp->d_is_CPML,
                                                            mp->d_spec_to_CPML,
                                                            mp->d_PML_displ_new,
                                                            mp->d_PML_displ_old,
                                                            mp->d_rmemory_displ_elastic,
                                                            mp->d_rmemory_dux_dxl_x,
                                                            mp->d_rmemory_duy_dxl_y,
                                                            mp->d_rmemory_duz_dxl_z,
                                                            mp->d_rmemory_dux_dyl_x,
                                                            mp->d_rmemory_duy_dyl_y,
                                                            mp->d_rmemory_duz_dyl_z,
                                                            mp->d_rmemory_dux_dzl_x,
                                                            mp->d_rmemory_duy_dzl_y,
                                                            mp->d_rmemory_duz_dzl_z,
                                                            mp->d_rmemory_dux_dxl_y,
                                                            mp->d_rmemory_dux_dxl_z,
                                                            mp->d_rmemory_duy_dxl_x,
                                                            mp->d_rmemory_duz_dxl_x,
                                                            mp->d_rmemory_dux_dyl_y,
                                                            mp->d_rmemory_duy_dyl_x,
                                                            mp->d_rmemory_duy_dyl_z,
                                                            mp->d_rmemory_duz_dyl_y,
                                                            mp->d_rmemory_dux_dzl_z,
                                                            mp->d_rmemory_duy_dzl_z,
                                                            mp->d_rmemory_duz_dzl_x,
                                                            mp->d_rmemory_duz_dzl_y,
                                                            mp->d_pml_convolution_coef_alpha,
                                                            mp->d_pml_convolution_coef_beta,
                                                            mp->d_pml_convolution_coef_strain,
                                                            mp->d_pml_convolution_coef_abar);
  }
#endif

  GPU_ERROR_CHECKING("compute_forces_viscoelastic_pml_cuda");
}


/* ----------------------------------------------------------------------------------------------- */

// PML Dirichlet condition on boundary interface

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
