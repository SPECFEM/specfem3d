/*
 !=====================================================================
 !
 !               S p e c f e m 3 D  V e r s i o n  3 . 0
 !               ---------------------------------------
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

// KERNEL 2 - acoustic compute forces kernel

/* ----------------------------------------------------------------------------------------------- */

void Kernel_2_acoustic(int nb_blocks_to_compute, Mesh* mp, int d_iphase,
                       int* d_ibool,
                       realw* d_xix,realw* d_xiy,realw* d_xiz,
                       realw* d_etax,realw* d_etay,realw* d_etaz,
                       realw* d_gammax,realw* d_gammay,realw* d_gammaz,
                       realw* d_rhostore,
                       realw* d_kappastore,
                       const int FORWARD_OR_ADJOINT){

  GPU_ERROR_CHECKING("before acoustic kernel Kernel 2");

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

  // note: for computational efficienty, the FORWARD_OR_ADJOINT variable here can have a special case (== 0)
  //       to combine forward and backward wavefield in the same kernel call

  // kernel timing
  gpu_event start, stop;
  if (CUDA_TIMING ){
    start_timing_gpu(&start,&stop);
  }
  int nb_field = mp->simulation_type == 3 ? 2 : 1 ;

  // sets gpu arrays
  field* potential, *potential_dot_dot;
  if (FORWARD_OR_ADJOINT == 1 || FORWARD_OR_ADJOINT == 0) {
    // forward wavefields -> FORWARD_OR_ADJOINT == 1
    potential = mp->d_potential_acoustic;
    potential_dot_dot = mp->d_potential_dot_dot_acoustic;
  } else {
    // for backward/reconstructed fields
    // backward wavefields -> FORWARD_OR_ADJOINT == 3
    potential = mp->d_b_potential_acoustic;
    potential_dot_dot = mp->d_b_potential_dot_dot_acoustic;
  }

  // acoustic kernel
  if (FORWARD_OR_ADJOINT == 0){
    // This kernel treats both forward and adjoint wavefield within the same call, to increase performance
    // ( ~37% faster for pure acoustic simulations )
#ifdef USE_CUDA
    if (run_cuda){
      Kernel_2_acoustic_impl<1><<<grid,threads,0,mp->compute_stream>>>(nb_blocks_to_compute,
                                                                        d_ibool,
                                                                        mp->d_irregular_element_number,
                                                                        mp->d_phase_ispec_inner_acoustic,
                                                                        mp->num_phase_ispec_acoustic,
                                                                        d_iphase,
                                                                        mp->d_potential_acoustic, mp->d_potential_dot_dot_acoustic,
                                                                        mp->d_b_potential_acoustic, mp->d_b_potential_dot_dot_acoustic,
                                                                        nb_field,
                                                                        d_xix, d_xiy, d_xiz,
                                                                        d_etax, d_etay, d_etaz,
                                                                        d_gammax, d_gammay, d_gammaz,
                                                                        mp->xix_regular,mp->jacobian_regular,
                                                                        mp->d_hprime_xx,
                                                                        mp->d_hprimewgll_xx,
                                                                        mp->d_wgllwgll_xy, mp->d_wgllwgll_xz, mp->d_wgllwgll_yz,
                                                                        d_rhostore,
                                                                        mp->use_mesh_coloring_gpu,
                                                                        mp->gravity,
                                                                        mp->d_minus_g,
                                                                        d_kappastore,
                                                                        mp->d_wgll_cube);
    }
#endif
#ifdef USE_HIP
    if (run_hip){
      hipLaunchKernelGGL(HIP_KERNEL_NAME(Kernel_2_acoustic_impl<1>), dim3(grid), dim3(threads), 0, mp->compute_stream,
                                                                     nb_blocks_to_compute,
                                                                     d_ibool,
                                                                     mp->d_irregular_element_number,
                                                                     mp->d_phase_ispec_inner_acoustic,
                                                                     mp->num_phase_ispec_acoustic,
                                                                     d_iphase,
                                                                     mp->d_potential_acoustic, mp->d_potential_dot_dot_acoustic,
                                                                     mp->d_b_potential_acoustic, mp->d_b_potential_dot_dot_acoustic,
                                                                     nb_field,
                                                                     d_xix, d_xiy, d_xiz,
                                                                     d_etax, d_etay, d_etaz,
                                                                     d_gammax, d_gammay, d_gammaz,
                                                                     mp->xix_regular,mp->jacobian_regular,
                                                                     mp->d_hprime_xx,
                                                                     mp->d_hprimewgll_xx,
                                                                     mp->d_wgllwgll_xy, mp->d_wgllwgll_xz, mp->d_wgllwgll_yz,
                                                                     d_rhostore,
                                                                     mp->use_mesh_coloring_gpu,
                                                                     mp->gravity,
                                                                     mp->d_minus_g,
                                                                     d_kappastore,
                                                                     mp->d_wgll_cube);
    }
#endif

  } else {
    // solving a single wavefield
#ifdef USE_CUDA
    if (run_cuda){
      Kernel_2_acoustic_single_impl<<<grid,threads,0,mp->compute_stream>>>(nb_blocks_to_compute,
                                                                           d_ibool,
                                                                           mp->d_phase_ispec_inner_acoustic,
                                                                           mp->num_phase_ispec_acoustic,
                                                                           d_iphase,
                                                                           potential,
                                                                           potential_dot_dot,
                                                                           d_xix, d_xiy, d_xiz,
                                                                           d_etax, d_etay, d_etaz,
                                                                           d_gammax, d_gammay, d_gammaz,
                                                                           mp->d_irregular_element_number,
                                                                           mp->xix_regular,mp->jacobian_regular,
                                                                           mp->d_hprime_xx,
                                                                           mp->d_hprimewgll_xx,
                                                                           mp->d_wgllwgll_xy, mp->d_wgllwgll_xz, mp->d_wgllwgll_yz,
                                                                           d_rhostore,
                                                                           mp->use_mesh_coloring_gpu,
                                                                           mp->gravity,
                                                                           mp->d_minus_g,
                                                                           d_kappastore,
                                                                           mp->d_wgll_cube,
                                                                           FORWARD_OR_ADJOINT);
    }
#endif
#ifdef USE_HIP
    if (run_hip){
      hipLaunchKernelGGL(Kernel_2_acoustic_single_impl, dim3(grid), dim3(threads), 0, mp->compute_stream,
                                                        nb_blocks_to_compute,
                                                        d_ibool,
                                                        mp->d_phase_ispec_inner_acoustic,
                                                        mp->num_phase_ispec_acoustic,
                                                        d_iphase,
                                                        potential,
                                                        potential_dot_dot,
                                                        d_xix, d_xiy, d_xiz,
                                                        d_etax, d_etay, d_etaz,
                                                        d_gammax, d_gammay, d_gammaz,
                                                        mp->d_irregular_element_number,
                                                        mp->xix_regular,mp->jacobian_regular,
                                                        mp->d_hprime_xx,
                                                        mp->d_hprimewgll_xx,
                                                        mp->d_wgllwgll_xy, mp->d_wgllwgll_xz, mp->d_wgllwgll_yz,
                                                        d_rhostore,
                                                        mp->use_mesh_coloring_gpu,
                                                        mp->gravity,
                                                        mp->d_minus_g,
                                                        d_kappastore,
                                                        mp->d_wgll_cube,
                                                        FORWARD_OR_ADJOINT);
    }
#endif

  }

  // Cuda timing
  if (CUDA_TIMING ){
    realw flops,time;
    stop_timing_gpu(&start,&stop,"Kernel_2_acoustic_impl",&time);
    // time in seconds
    time = time / 1000.;
    // performance
    if (! mp->gravity) {
      if (! mp->use_mesh_coloring_gpu ){
        // see with: nvprof --metrics flops_sp ./xspecfem3D
        //           -> using 322631424 FLOPS (Single) floating-point operations for 20736 elements
        //              = 15559 FLOPS per block
        flops = 15559 * nb_blocks_to_compute;
      }else{
        // coloring
        flops = 15559 * nb_blocks_to_compute;
      }
    }else{
      // gravity
      flops = 15559 * nb_blocks_to_compute;
    }
    printf("  performance: %f GFlop/s\n", flops/time * 1.e-9);
  }

  GPU_ERROR_CHECKING("kernel Kernel_2");
}

/* ----------------------------------------------------------------------------------------------- */

// main compute_forces_acoustic CUDA routine

/* ----------------------------------------------------------------------------------------------- */

extern EXTERN_LANG
void FC_FUNC_(compute_forces_acoustic_cuda,
              COMPUTE_FORCES_ACOUSTIC_CUDA)(long* Mesh_pointer,
                                            int* iphase,
                                            int* nspec_outer_acoustic,
                                            int* nspec_inner_acoustic,
                                            int* FORWARD_OR_ADJOINT_f) {

  TRACE("compute_forces_acoustic_cuda");
  //double start_time = get_time_val();

  Mesh* mp = (Mesh*)(*Mesh_pointer); // get Mesh from fortran integer wrapper
  int FORWARD_OR_ADJOINT = *FORWARD_OR_ADJOINT_f;

  // determines number of elements to loop over (inner/outer elements)
  int num_elements;
  if (*iphase == 1)
    num_elements = *nspec_outer_acoustic;
  else
    num_elements = *nspec_inner_acoustic;

  if (num_elements == 0) return;

  // mesh coloring
  if (mp->use_mesh_coloring_gpu ){

    // note: array offsets require sorted arrays, such that e.g. ibool starts with elastic elements
    //         and followed by acoustic ones.
    //         acoustic elements also start with outer than inner element ordering

    int nb_colors,nb_blocks_to_compute;
    int istart;
    int offset,offset_nonpadded;

    // sets up color loop
    if (*iphase == 1){
      // outer elements
      nb_colors = mp->num_colors_outer_acoustic;
      istart = 0;

      // array offsets (acoustic elements start after elastic ones)
      offset = mp->nspec_elastic * NGLL3_PADDED;
      offset_nonpadded = mp->nspec_elastic * NGLL3;
    }else{
      // inner element colors (start after outer elements)
      nb_colors = mp->num_colors_outer_acoustic + mp->num_colors_inner_acoustic;
      istart = mp->num_colors_outer_acoustic;

      // array offsets (inner elements start after outer ones)
      offset = ( mp->nspec_elastic + (*nspec_outer_acoustic) ) * NGLL3_PADDED;
      offset_nonpadded = ( mp->nspec_elastic + (*nspec_outer_acoustic) ) * NGLL3;
    }

    // loops over colors
    for(int icolor = istart; icolor < nb_colors; icolor++){

      nb_blocks_to_compute = mp->h_num_elem_colors_acoustic[icolor];

      Kernel_2_acoustic(nb_blocks_to_compute,mp,*iphase,
                         mp->d_ibool + offset,
                         mp->d_xix + offset,mp->d_xiy + offset,mp->d_xiz + offset,
                         mp->d_etax + offset,mp->d_etay + offset,mp->d_etaz + offset,
                         mp->d_gammax + offset,mp->d_gammay + offset,mp->d_gammaz + offset,
                         mp->d_rhostore + offset,
                         mp->d_kappastore + offset_nonpadded,
                         FORWARD_OR_ADJOINT);

      // for padded and aligned arrays
      offset += nb_blocks_to_compute * NGLL3_PADDED;
      // for no-aligned arrays
      offset_nonpadded += nb_blocks_to_compute * NGLL3;
    }

  }else{

    // no mesh coloring: uses atomic updates
    Kernel_2_acoustic(num_elements, mp, *iphase,
                      mp->d_ibool,
                      mp->d_xix,mp->d_xiy,mp->d_xiz,
                      mp->d_etax,mp->d_etay,mp->d_etaz,
                      mp->d_gammax,mp->d_gammay,mp->d_gammaz,
                      mp->d_rhostore,
                      mp->d_kappastore,
                      FORWARD_OR_ADJOINT);

  }
}



/* ----------------------------------------------------------------------------------------------- */

/* KERNEL for enforce free surface */

/* ----------------------------------------------------------------------------------------------- */

extern EXTERN_LANG
void FC_FUNC_(acoustic_enforce_free_surf_cuda,
              ACOUSTIC_ENFORCE_FREE_SURF_CUDA)(long* Mesh_pointer,
                                               int* ABSORB_INSTEAD_OF_FREE_SURFACE,
                                               int* FORWARD_OR_ADJOINT) {

TRACE("acoustic_enforce_free_surf_cuda");

  Mesh* mp = (Mesh*)(*Mesh_pointer); //get mesh pointer out of fortran integer container

  // checks if anything to do
  if (*ABSORB_INSTEAD_OF_FREE_SURFACE == 0){

    // does not absorb free surface, thus we enforce the potential to be zero at surface

    // block sizes
    int num_blocks_x, num_blocks_y;
    get_blocks_xy(mp->num_free_surface_faces,&num_blocks_x,&num_blocks_y);

    dim3 grid(num_blocks_x,num_blocks_y,1);
    dim3 threads(NGLL2,1,1);

    // sets gpu arrays
    field* potential, *potential_dot, *potential_dot_dot;
    if (*FORWARD_OR_ADJOINT == 1) {
      potential = mp->d_potential_acoustic;
      potential_dot = mp->d_potential_dot_acoustic;
      potential_dot_dot = mp->d_potential_dot_dot_acoustic;
    } else {
      // for backward/reconstructed fields
      potential = mp->d_b_potential_acoustic;
      potential_dot = mp->d_b_potential_dot_acoustic;
      potential_dot_dot = mp->d_b_potential_dot_dot_acoustic;
    }

    // sets potentials to zero at free surface
#ifdef USE_CUDA
    if (run_cuda){
      enforce_free_surface_cuda_kernel<<<grid,threads,0,mp->compute_stream>>>(potential,
                                                                              potential_dot,
                                                                              potential_dot_dot,
                                                                              mp->num_free_surface_faces,
                                                                              mp->d_free_surface_ispec,
                                                                              mp->d_free_surface_ijk,
                                                                              mp->d_ibool,
                                                                              mp->d_ispec_is_acoustic);
    }
#endif
#ifdef USE_HIP
    if (run_hip){
      hipLaunchKernelGGL(enforce_free_surface_cuda_kernel, dim3(grid), dim3(threads), 0, mp->compute_stream,
                                                           potential,
                                                           potential_dot,
                                                           potential_dot_dot,
                                                           mp->num_free_surface_faces,
                                                           mp->d_free_surface_ispec,
                                                           mp->d_free_surface_ijk,
                                                           mp->d_ibool,
                                                           mp->d_ispec_is_acoustic);
    }
#endif
  }

  GPU_ERROR_CHECKING("enforce_free_surface_cuda");
}

