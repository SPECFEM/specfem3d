/*
 !=====================================================================
 !
 !               S p e c f e m 3 D  V e r s i o n  2 . 1
 !               ---------------------------------------
 !
 !          Main authors: Dimitri Komatitsch and Jeroen Tromp
 !    Princeton University, USA and CNRS / INRIA / University of Pau
 ! (c) Princeton University / California Institute of Technology and CNRS / INRIA / University of Pau
 !                             July 2012
 !
 ! This program is free software; you can redistribute it and/or modify
 ! it under the terms of the GNU General Public License as published by
 ! the Free Software Foundation; either version 2 of the License, or
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

#include <stdio.h>
#include <cuda.h>
#include <cublas.h>

#include <sys/time.h>
#include <sys/resource.h>

#include "config.h"
#include "mesh_constants_cuda.h"

#ifdef USE_TEXTURES_FIELDS
realw_texture d_potential_tex;
realw_texture d_potential_dot_dot_tex;
//backward/reconstructed
realw_texture d_b_potential_tex;
realw_texture d_b_potential_dot_dot_tex;

//note: texture variables are implicitly static, and cannot be passed as arguments to cuda kernels;
//      thus, 1) we thus use if-statements (FORWARD_OR_ADJOINT) to determine from which texture to fetch from
//            2) we use templates
//      since if-statements are a bit slower as the variable is only known at runtime, we use option 2)

// templates definitions
template<int FORWARD_OR_ADJOINT> __device__ float texfetch_potential(int x);
template<int FORWARD_OR_ADJOINT> __device__ float texfetch_potential_dot_dot(int x);

// templates for texture fetching
// FORWARD_OR_ADJOINT == 1 <- forward arrays
template<> __device__ float texfetch_potential<1>(int x) { return tex1Dfetch(d_potential_tex, x); }
template<> __device__ float texfetch_potential_dot_dot<1>(int x) { return tex1Dfetch(d_potential_dot_dot_tex, x); }
// FORWARD_OR_ADJOINT == 3 <- backward/reconstructed arrays
template<> __device__ float texfetch_potential<3>(int x) { return tex1Dfetch(d_b_potential_tex, x); }
template<> __device__ float texfetch_potential_dot_dot<3>(int x) { return tex1Dfetch(d_b_potential_dot_dot_tex, x); }

#endif

#ifdef USE_TEXTURES_CONSTANTS
extern realw_texture d_hprime_xx_tex;
#endif


/* ----------------------------------------------------------------------------------------------- */

// KERNEL 2 - acoustic compute forces kernel

/* ----------------------------------------------------------------------------------------------- */


// acoustic kernel without gravity and without mesh coloring

template<int FORWARD_OR_ADJOINT> __global__ void
#ifdef USE_LAUNCH_BOUNDS
// adds compiler specification
__launch_bounds__(NGLL3_PADDED,LAUNCH_MIN_BLOCKS_ACOUSTIC)
#endif
Kernel_2_acoustic_impl(int nb_blocks_to_compute,
                       const int NGLOB,
                       const int* d_ibool,
                       const int* d_phase_ispec_inner_acoustic,
                       const int num_phase_ispec_acoustic,
                       const int d_iphase,
                       realw_const_p d_potential_acoustic,
                       realw_p d_potential_dot_dot_acoustic,
                       realw* d_xix,realw* d_xiy,realw* d_xiz,
                       realw* d_etax,realw* d_etay,realw* d_etaz,
                       realw* d_gammax,realw* d_gammay,realw* d_gammaz,
                       realw_const_p d_hprime_xx,
                       realw_const_p hprimewgll_xx,
                       realw_const_p wgllwgll_xy,realw_const_p wgllwgll_xz,realw_const_p wgllwgll_yz,
                       realw* d_rhostore){

  // block-id == number of local element id in phase_ispec array
  int bx = blockIdx.y*gridDim.x+blockIdx.x;

  // checks if anything to do
  if( bx >= nb_blocks_to_compute ) return;

  // thread-id == GLL node id
  // note: use only NGLL^3 = 125 active threads, plus 3 inactive/ghost threads,
  //       because we used memory padding from NGLL^3 = 125 to 128 to get coalescent memory accesses;
  //       to avoid execution branching and the need of registers to store an active state variable,
  //       the thread ids are put in valid range
  int tx = threadIdx.x;
  if( tx >= NGLL3 ) tx = NGLL3-1;
  
  int K = (tx/NGLL2);
  int J = ((tx-K*NGLL2)/NGLLX);
  int I = (tx-K*NGLL2-J*NGLLX);

  int iglob,offset;
  int working_element;

  realw temp1l,temp2l,temp3l;
  realw xixl,xiyl,xizl,etaxl,etayl,etazl,gammaxl,gammayl,gammazl,jacobianl;

  realw dpotentialdxl,dpotentialdyl,dpotentialdzl;
  realw fac1,fac2,fac3;
  realw rho_invl;

  realw sum_terms;

  __shared__ realw s_dummy_loc[NGLL3];

  __shared__ realw s_temp1[NGLL3];
  __shared__ realw s_temp2[NGLL3];
  __shared__ realw s_temp3[NGLL3];

  __shared__ realw sh_hprime_xx[NGLL2];
  __shared__ realw sh_hprimewgll_xx[NGLL2];

  // loads hprime into shared memory
  if (tx < NGLL2) {
#ifdef USE_TEXTURES_CONSTANTS
    sh_hprime_xx[tx] = tex1Dfetch(d_hprime_xx_tex,tx);
#else
    sh_hprime_xx[tx] = d_hprime_xx[tx];
#endif
    // loads hprimewgll into shared memory
    sh_hprimewgll_xx[tx] = hprimewgll_xx[tx];
  }

  // spectral-element id
  // iphase-1 and working_element-1 for Fortran->C array conventions
  working_element = d_phase_ispec_inner_acoustic[bx + num_phase_ispec_acoustic*(d_iphase-1)] - 1;

  // local padded index
  offset = working_element*NGLL3_PADDED + tx;

  // loads potential values into shared memory
  if(threadIdx.x < NGLL3) {
    // global index
    iglob = d_ibool[offset] - 1;

    // loads potentials
#ifdef USE_TEXTURES_FIELDS
    s_dummy_loc[tx] = texfetch_potential<FORWARD_OR_ADJOINT>(iglob);
#else
    // changing iglob indexing to match fortran row changes fast style
    s_dummy_loc[tx] = d_potential_acoustic[iglob];
#endif
  }

  // loads mesh values here to give compiler possibility to overlap memory fetches with some computations
  // note: arguments defined as realw* instead of const realw* __restrict__ to avoid that the compiler
  //       loads all memory by texture loads
  //       we only use the first loads explicitly by texture loads, all subsequent without. this should lead/trick
  //       the compiler to use global memory loads for all the subsequent accesses.
  //
  // calculates laplacian
  xixl = get_global_cr( &d_xix[offset] ); // first array with texture load
  xiyl = d_xiy[offset]; // all subsequent without to avoid over-use of texture for coalescent access
  xizl = d_xiz[offset];
  etaxl = d_etax[offset];
  etayl = d_etay[offset];
  etazl = d_etaz[offset];
  gammaxl = d_gammax[offset];
  gammayl = d_gammay[offset];
  gammazl = d_gammaz[offset];

  jacobianl = 1.f / (xixl*(etayl*gammazl-etazl*gammayl)
                    -xiyl*(etaxl*gammazl-etazl*gammaxl)
                    +xizl*(etaxl*gammayl-etayl*gammaxl));

  // density (reciproc)
  rho_invl = 1.f / d_rhostore[offset];

  // synchronize all the threads (one thread for each of the NGLL grid points of the
  // current spectral element) because we need the whole element to be ready in order
  // to be able to compute the matrix products along cut planes of the 3D element below
  __syncthreads();

  // computes first matrix product
  temp1l = 0.f;
  temp2l = 0.f;
  temp3l = 0.f;
  for (int l=0;l<NGLLX;l++) {
    //assumes that hprime_xx = hprime_yy = hprime_zz
    fac1 = sh_hprime_xx[l*NGLLX+I];
    temp1l += s_dummy_loc[K*NGLL2+J*NGLLX+l]*fac1;

    fac2 = sh_hprime_xx[l*NGLLX+J];
    temp2l += s_dummy_loc[K*NGLL2+l*NGLLX+I]*fac2;

    fac3 = sh_hprime_xx[l*NGLLX+K];
    temp3l += s_dummy_loc[l*NGLL2+J*NGLLX+I]*fac3;
  }

  // compute derivatives of ux, uy and uz with respect to x, y and z
  // derivatives of potential
  dpotentialdxl = xixl*temp1l + etaxl*temp2l + gammaxl*temp3l;
  dpotentialdyl = xiyl*temp1l + etayl*temp2l + gammayl*temp3l;
  dpotentialdzl = xizl*temp1l + etazl*temp2l + gammazl*temp3l;

  // form the dot product with the test vector
  if( threadIdx.x < NGLL3 ) {
    s_temp1[tx] = jacobianl * rho_invl * (dpotentialdxl*xixl + dpotentialdyl*xiyl + dpotentialdzl*xizl);
    s_temp2[tx] = jacobianl * rho_invl * (dpotentialdxl*etaxl + dpotentialdyl*etayl + dpotentialdzl*etazl);
    s_temp3[tx] = jacobianl * rho_invl * (dpotentialdxl*gammaxl + dpotentialdyl*gammayl + dpotentialdzl*gammazl);
  }

  // synchronize all the threads (one thread for each of the NGLL grid points of the
  // current spectral element) because we need the whole element to be ready in order
  // to be able to compute the matrix products along cut planes of the 3D element below
  __syncthreads();

  // computes second matrix product
  temp1l = 0.f;
  temp2l = 0.f;
  temp3l = 0.f;
  for (int l=0;l<NGLLX;l++) {
    //assumes hprimewgll_xx = hprimewgll_yy = hprimewgll_zz
    fac1 = sh_hprimewgll_xx[I*NGLLX+l]; // hprimewgll_xx[I*NGLLX+l];
    temp1l += s_temp1[K*NGLL2+J*NGLLX+l]*fac1;

    fac2 = sh_hprimewgll_xx[J*NGLLX+l]; // hprimewgll_xx[J*NGLLX+l];
    temp2l += s_temp2[K*NGLL2+l*NGLLX+I]*fac2;

    fac3 = sh_hprimewgll_xx[K*NGLLX+l]; // hprimewgll_xx[K*NGLLX+l];
    temp3l += s_temp3[l*NGLL2+J*NGLLX+I]*fac3;
  }

  // summed terms with added gll weights
  sum_terms = -(wgllwgll_yz[K*NGLLX+J]*temp1l + wgllwgll_xz[K*NGLLX+I]*temp2l + wgllwgll_xy[J*NGLLX+I]*temp3l);

  // assembles potential array
  if(threadIdx.x < NGLL3) {
      atomicAdd(&d_potential_dot_dot_acoustic[iglob],sum_terms);
  }

}

/* ----------------------------------------------------------------------------------------------- */

// acoustic kernel without gravity but with coloring

template<int FORWARD_OR_ADJOINT> __global__ void
#ifdef USE_LAUNCH_BOUNDS
// adds compiler specification
__launch_bounds__(NGLL3_PADDED,LAUNCH_MIN_BLOCKS_ACOUSTIC)
#endif
Kernel_2_acoustic_col_impl(int nb_blocks_to_compute,
                           const int NGLOB,
                           const int* d_ibool,
                           const int* d_phase_ispec_inner_acoustic,
                           const int num_phase_ispec_acoustic,
                           const int d_iphase,
                           realw_const_p d_potential_acoustic,
                           realw_p d_potential_dot_dot_acoustic,
                           realw* d_xix,realw* d_xiy,realw* d_xiz,
                           realw* d_etax,realw* d_etay,realw* d_etaz,
                           realw* d_gammax,realw* d_gammay,realw* d_gammaz,
                           realw_const_p d_hprime_xx,
                           realw_const_p hprimewgll_xx,
                           realw_const_p wgllwgll_xy,realw_const_p wgllwgll_xz,realw_const_p wgllwgll_yz,
                           realw* d_rhostore,
                           const int use_mesh_coloring_gpu){

  // block-id == spectral-element id
  int bx = blockIdx.y*gridDim.x+blockIdx.x;

  // checks if anything to do
  if( bx >= nb_blocks_to_compute ) return;

  // thread-id == GLL node id
  // note: use only NGLL^3 = 125 active threads, plus 3 inactive/ghost threads,
  //       because we used memory padding from NGLL^3 = 125 to 128 to get coalescent memory accesses;
  //       to avoid execution branching and the need of registers to store an active state variable,
  //       the thread ids are put in valid range
  int tx = threadIdx.x;
  if( tx >= NGLL3 ) tx = NGLL3-1;
  
  int K = (tx/NGLL2);
  int J = ((tx-K*NGLL2)/NGLLX);
  int I = (tx-K*NGLL2-J*NGLLX);

  int iglob,offset;
  int working_element;

  realw temp1l,temp2l,temp3l;
  realw xixl,xiyl,xizl,etaxl,etayl,etazl,gammaxl,gammayl,gammazl,jacobianl;

  realw dpotentialdxl,dpotentialdyl,dpotentialdzl;
  realw fac1,fac2,fac3;
  realw rho_invl;

  realw sum_terms;

  __shared__ realw s_dummy_loc[NGLL3];

  __shared__ realw s_temp1[NGLL3];
  __shared__ realw s_temp2[NGLL3];
  __shared__ realw s_temp3[NGLL3];

  __shared__ realw sh_hprime_xx[NGLL2];
  __shared__ realw sh_hprimewgll_xx[NGLL2];

  // loads hprime into shared memory
  if (tx < NGLL2) {
#ifdef USE_TEXTURES_CONSTANTS
    sh_hprime_xx[tx] = tex1Dfetch(d_hprime_xx_tex,tx);
#else
    sh_hprime_xx[tx] = d_hprime_xx[tx];
#endif
    // loads hprimewgll into shared memory
    sh_hprimewgll_xx[tx] = hprimewgll_xx[tx];
  }

  // copy from global memory to shared memory
  // each thread writes one of the NGLL^3 = 125 data points
#ifdef USE_MESH_COLORING_GPU
  working_element = bx;
#else
  //mesh coloring
  if( use_mesh_coloring_gpu ){
    working_element = bx;
  }else{
    // iphase-1 and working_element-1 for Fortran->C array conventions
    working_element = d_phase_ispec_inner_acoustic[bx + num_phase_ispec_acoustic*(d_iphase-1)]-1;
  }
#endif
  // local padded index
  offset = working_element*NGLL3_PADDED + tx;

  // loads potential values into shared memory
  if(threadIdx.x < NGLL3) {
    // global index
    iglob = d_ibool[offset] - 1;

    // loads potentials
#ifdef USE_TEXTURES_FIELDS
    s_dummy_loc[tx] = texfetch_potential<FORWARD_OR_ADJOINT>(iglob);
#else
    // changing iglob indexing to match fortran row changes fast style
    s_dummy_loc[tx] = d_potential_acoustic[iglob];
#endif
  }

  // loads mesh values here to give compiler possibility to overlap memory fetches with some computations
  // note: arguments defined as realw* instead of const realw* __restrict__ to avoid that the compiler
  //       loads all memory by texture loads
  //       we only use the first loads explicitly by texture loads, all subsequent without. this should lead/trick
  //       the compiler to use global memory loads for all the subsequent accesses.
  //
  // calculates laplacian
  xixl = get_global_cr( &d_xix[offset] ); // first array with texture load
  xiyl = d_xiy[offset]; // all subsequent without to avoid over-use of texture for coalescent access
  xizl = d_xiz[offset];
  etaxl = d_etax[offset];
  etayl = d_etay[offset];
  etazl = d_etaz[offset];
  gammaxl = d_gammax[offset];
  gammayl = d_gammay[offset];
  gammazl = d_gammaz[offset];

  jacobianl = 1.f / (xixl*(etayl*gammazl-etazl*gammayl)
                    -xiyl*(etaxl*gammazl-etazl*gammaxl)
                    +xizl*(etaxl*gammayl-etayl*gammaxl));

  // density (reciproc)
  rho_invl = 1.f / d_rhostore[offset];

  // synchronize all the threads (one thread for each of the NGLL grid points of the
  // current spectral element) because we need the whole element to be ready in order
  // to be able to compute the matrix products along cut planes of the 3D element below
  __syncthreads();

  // computes first matrix product
  temp1l = 0.f;
  temp2l = 0.f;
  temp3l = 0.f;
  for (int l=0;l<NGLLX;l++) {
    //assumes that hprime_xx = hprime_yy = hprime_zz
    fac1 = sh_hprime_xx[l*NGLLX+I];
    temp1l += s_dummy_loc[K*NGLL2+J*NGLLX+l]*fac1;

    fac2 = sh_hprime_xx[l*NGLLX+J];
    temp2l += s_dummy_loc[K*NGLL2+l*NGLLX+I]*fac2;

    fac3 = sh_hprime_xx[l*NGLLX+K];
    temp3l += s_dummy_loc[l*NGLL2+J*NGLLX+I]*fac3;
  }

  // compute derivatives of ux, uy and uz with respect to x, y and z
  // derivatives of potential
  dpotentialdxl = xixl*temp1l + etaxl*temp2l + gammaxl*temp3l;
  dpotentialdyl = xiyl*temp1l + etayl*temp2l + gammayl*temp3l;
  dpotentialdzl = xizl*temp1l + etazl*temp2l + gammazl*temp3l;

  // form the dot product with the test vector
  if( threadIdx.x < NGLL3 ) {
    s_temp1[tx] = jacobianl * rho_invl * (dpotentialdxl*xixl + dpotentialdyl*xiyl + dpotentialdzl*xizl);
    s_temp2[tx] = jacobianl * rho_invl * (dpotentialdxl*etaxl + dpotentialdyl*etayl + dpotentialdzl*etazl);
    s_temp3[tx] = jacobianl * rho_invl * (dpotentialdxl*gammaxl + dpotentialdyl*gammayl + dpotentialdzl*gammazl);
  }

  // synchronize all the threads (one thread for each of the NGLL grid points of the
  // current spectral element) because we need the whole element to be ready in order
  // to be able to compute the matrix products along cut planes of the 3D element below
  __syncthreads();

  // computes second matrix product
  temp1l = 0.f;
  temp2l = 0.f;
  temp3l = 0.f;
  for (int l=0;l<NGLLX;l++) {
    //assumes hprimewgll_xx = hprimewgll_yy = hprimewgll_zz
    fac1 = sh_hprimewgll_xx[I*NGLLX+l]; // hprimewgll_xx[I*NGLLX+l];
    temp1l += s_temp1[K*NGLL2+J*NGLLX+l]*fac1;

    fac2 = sh_hprimewgll_xx[J*NGLLX+l]; // hprimewgll_xx[J*NGLLX+l];
    temp2l += s_temp2[K*NGLL2+l*NGLLX+I]*fac2;

    fac3 = sh_hprimewgll_xx[K*NGLLX+l]; // hprimewgll_xx[K*NGLLX+l];
    temp3l += s_temp3[l*NGLL2+J*NGLLX+I]*fac3;
  }

  // summed terms with added gll weights
  sum_terms = -(wgllwgll_yz[K*NGLLX+J]*temp1l + wgllwgll_xz[K*NGLLX+I]*temp2l + wgllwgll_xy[J*NGLLX+I]*temp3l);

  // assembles potential array
  if(threadIdx.x < NGLL3) {
#ifdef USE_MESH_COLORING_GPU
    // no atomic operation needed, colors don't share global points between elements
#ifdef USE_TEXTURES_FIELDS
    d_potential_dot_dot_acoustic[iglob] = texfetch_potential_dot_dot<FORWARD_OR_ADJOINT>(iglob) + sum_terms;
#else
    d_potential_dot_dot_acoustic[iglob] += sum_terms;
#endif // USE_TEXTURES_FIELDS

#else  // MESH_COLORING

    //mesh coloring
    if( use_mesh_coloring_gpu ){
      // no atomic operation needed, colors don't share global points between elements
#ifdef USE_TEXTURES_FIELDS
      d_potential_dot_dot_acoustic[iglob] = texfetch_potential_dot_dot<FORWARD_OR_ADJOINT>(iglob) + sum_terms;
#else
      d_potential_dot_dot_acoustic[iglob] += sum_terms;
#endif // USE_TEXTURES_FIELDS

    }else{
      atomicAdd(&d_potential_dot_dot_acoustic[iglob],sum_terms);
    }
#endif // MESH_COLORING
  } // threadIdx.x

}



/* ----------------------------------------------------------------------------------------------- */

// acoustic kernel with gravity

template<int FORWARD_OR_ADJOINT> __global__ void
#ifdef USE_LAUNCH_BOUNDS
// adds compiler specification
__launch_bounds__(NGLL3_PADDED,LAUNCH_MIN_BLOCKS_ACOUSTIC)
#endif
Kernel_2_acoustic_grav_impl(int nb_blocks_to_compute,
                            const int NGLOB,
                            const int* d_ibool,
                            const int* d_phase_ispec_inner_acoustic,
                            const int num_phase_ispec_acoustic,
                            const int d_iphase,
                            realw_const_p d_potential_acoustic,
                            realw_p d_potential_dot_dot_acoustic,
                            realw* d_xix,realw* d_xiy,realw* d_xiz,
                            realw* d_etax,realw* d_etay,realw* d_etaz,
                            realw* d_gammax,realw* d_gammay,realw* d_gammaz,
                            realw_const_p d_hprime_xx,
                            realw_const_p hprimewgll_xx,
                            realw_const_p wgllwgll_xy,realw_const_p wgllwgll_xz,realw_const_p wgllwgll_yz,
                            realw* d_rhostore,
                            const int use_mesh_coloring_gpu,
                            const int gravity,
                            realw_const_p d_minus_g,
                            realw* d_kappastore,
                            realw_const_p d_wgll_cube){

  // block-id == spectral-element id
  int bx = blockIdx.y*gridDim.x+blockIdx.x;

  // checks if anything to do
  if( bx >= nb_blocks_to_compute ) return;

  // thread-id == GLL node id
  // use only NGLL^3 = 125 active threads, plus 3 inactive/ghost threads,
  // because we used memory padding from NGLL^3 = 125 to 128 to get coalescent memory accesses
  //
  // to avoid execution branching, thread ids are put in valid range
  int tx = threadIdx.x;
  if( tx >= NGLL3 ) tx = NGLL3-1;
  
  int K = (tx/NGLL2);
  int J = ((tx-K*NGLL2)/NGLLX);
  int I = (tx-K*NGLL2-J*NGLLX);

  int iglob,offset;
  int working_element;

  realw temp1l,temp2l,temp3l;
  realw xixl,xiyl,xizl,etaxl,etayl,etazl,gammaxl,gammayl,gammazl,jacobianl;

  realw dpotentialdxl,dpotentialdyl,dpotentialdzl;
  realw fac1,fac2,fac3;
  realw rho_invl,kappa_invl;

  realw sum_terms;
  realw gravity_term;

  __shared__ realw s_dummy_loc[NGLL3];

  __shared__ realw s_temp1[NGLL3];
  __shared__ realw s_temp2[NGLL3];
  __shared__ realw s_temp3[NGLL3];

  __shared__ realw sh_hprime_xx[NGLL2];
  __shared__ realw sh_hprimewgll_xx[NGLL2];

  // loads hprime into shared memory
  if (tx < NGLL2) {
#ifdef USE_TEXTURES_CONSTANTS
    sh_hprime_xx[tx] = tex1Dfetch(d_hprime_xx_tex,tx);
#else
    sh_hprime_xx[tx] = d_hprime_xx[tx];
#endif
    // loads hprimewgll into shared memory
    sh_hprimewgll_xx[tx] = hprimewgll_xx[tx];
  }

  // copy from global memory to shared memory
  // each thread writes one of the NGLL^3 = 125 data points
#ifdef USE_MESH_COLORING_GPU
  working_element = bx;
#else
  //mesh coloring
  if( use_mesh_coloring_gpu ){
    working_element = bx;
  }else{
    // iphase-1 and working_element-1 for Fortran->C array conventions
    working_element = d_phase_ispec_inner_acoustic[bx + num_phase_ispec_acoustic*(d_iphase-1)]-1;
  }
#endif
  // local padded index
  offset = working_element*NGLL3_PADDED + tx;

  // loads potential values into shared memory
  if(threadIdx.x < NGLL3) {
    // global index
    iglob = d_ibool[offset] - 1;

#ifdef USE_TEXTURES_FIELDS
    s_dummy_loc[tx] = texfetch_potential<FORWARD_OR_ADJOINT>(iglob);
#else
    // changing iglob indexing to match fortran row changes fast style
    s_dummy_loc[tx] = d_potential_acoustic[iglob];
#endif
  }

  // loads mesh values here to give compiler possibility to overlap memory fetches with some computations
  // note: arguments defined as realw* instead of const realw* __restrict__ to avoid that the compiler
  //       loads all memory by texture loads
  //       we only use the first loads explicitly by texture loads, all subsequent without. this should lead/trick
  //       the compiler to use global memory loads for all the subsequent accesses.
  //
  // calculates laplacian
  xixl = get_global_cr( &d_xix[offset] ); // first array with texture load
  xiyl = d_xiy[offset]; // all subsequent without to avoid over-use of texture for coalescent access
  xizl = d_xiz[offset];
  etaxl = d_etax[offset];
  etayl = d_etay[offset];
  etazl = d_etaz[offset];
  gammaxl = d_gammax[offset];
  gammayl = d_gammay[offset];
  gammazl = d_gammaz[offset];

  jacobianl = 1.f / (xixl*(etayl*gammazl-etazl*gammayl)
                    -xiyl*(etaxl*gammazl-etazl*gammaxl)
                    +xizl*(etaxl*gammayl-etayl*gammaxl));

  // density (reciproc)
  rho_invl = 1.f / d_rhostore[offset];

  // gravity
  if( gravity ){
    // gravity term: 1/kappa grad(chi) * g
    // assumes that g only acts in (negative) z-direction
    kappa_invl = 1.f / d_kappastore[working_element*NGLL3 + tx];
    // daniel: TODO - check gravity
    //if( kappa_invl <= 0.0f ){
    //  printf("kappa error: %f %f\n",kappa_invl,d_kappastore[working_element*NGLL3 + tx]);
    //  printf("kappa error: thread %d %d \n",tx,working_element);
    //  asm("trap;");
    //}
    //if( iglob <= 0 ){
    //  printf("iglob error: %d %d %d \n",iglob,tx,working_element);
    //  asm("trap;");
    //}
  }
  
  // synchronize all the threads (one thread for each of the NGLL grid points of the
  // current spectral element) because we need the whole element to be ready in order
  // to be able to compute the matrix products along cut planes of the 3D element below
  __syncthreads();

  // computes first matrix product
  temp1l = 0.f;
  temp2l = 0.f;
  temp3l = 0.f;

  for (int l=0;l<NGLLX;l++) {
    //assumes that hprime_xx = hprime_yy = hprime_zz
    fac1 = sh_hprime_xx[l*NGLLX+I];
    temp1l += s_dummy_loc[K*NGLL2+J*NGLLX+l]*fac1;

    fac2 = sh_hprime_xx[l*NGLLX+J];
    temp2l += s_dummy_loc[K*NGLL2+l*NGLLX+I]*fac2;

    fac3 = sh_hprime_xx[l*NGLLX+K];
    temp3l += s_dummy_loc[l*NGLL2+J*NGLLX+I]*fac3;
  }

  // compute derivatives of ux, uy and uz with respect to x, y and z
  // derivatives of potential
  dpotentialdxl = xixl*temp1l + etaxl*temp2l + gammaxl*temp3l;
  dpotentialdyl = xiyl*temp1l + etayl*temp2l + gammayl*temp3l;
  dpotentialdzl = xizl*temp1l + etazl*temp2l + gammazl*temp3l;

  // pre-computes gravity sum term
  if( gravity ){
    // uses potential definition: s = grad(chi)

    // gravity term: 1/kappa grad(chi) * g
    // assumes that g only acts in (negative) z-direction
    gravity_term = d_minus_g[iglob] * kappa_invl * jacobianl * d_wgll_cube[tx] * dpotentialdzl;

    // daniel: TODO - check gravity
    //gravity_term = 0.f;
    //if( iglob == 5 ){
    //  printf("iglob infos: %f %f %f %f %f \n",d_minus_g[iglob],kappa_invl,jacobianl,d_wgll_cube[tx],dpotentialdzl);
    //}
  }

  // form the dot product with the test vector
  if( threadIdx.x < NGLL3 ) {
    s_temp1[tx] = jacobianl * rho_invl * (dpotentialdxl*xixl + dpotentialdyl*xiyl + dpotentialdzl*xizl);
    s_temp2[tx] = jacobianl * rho_invl * (dpotentialdxl*etaxl + dpotentialdyl*etayl + dpotentialdzl*etazl);
    s_temp3[tx] = jacobianl * rho_invl * (dpotentialdxl*gammaxl + dpotentialdyl*gammayl + dpotentialdzl*gammazl);
  }

  // synchronize all the threads (one thread for each of the NGLL grid points of the
  // current spectral element) because we need the whole element to be ready in order
  // to be able to compute the matrix products along cut planes of the 3D element below
  __syncthreads();

  // computes second matrix product
  temp1l = 0.f;
  temp2l = 0.f;
  temp3l = 0.f;

  for (int l=0;l<NGLLX;l++) {
    //assumes hprimewgll_xx = hprimewgll_yy = hprimewgll_zz
    fac1 = sh_hprimewgll_xx[I*NGLLX+l]; // hprimewgll_xx[I*NGLLX+l];
    temp1l += s_temp1[K*NGLL2+J*NGLLX+l]*fac1;

    fac2 = sh_hprimewgll_xx[J*NGLLX+l]; // hprimewgll_xx[J*NGLLX+l];
    temp2l += s_temp2[K*NGLL2+l*NGLLX+I]*fac2;

    fac3 = sh_hprimewgll_xx[K*NGLLX+l]; // hprimewgll_xx[K*NGLLX+l];
    temp3l += s_temp3[l*NGLL2+J*NGLLX+I]*fac3;
  }

  // summed terms with added gll weights
  sum_terms = -(wgllwgll_yz[K*NGLLX+J]*temp1l + wgllwgll_xz[K*NGLLX+I]*temp2l + wgllwgll_xy[J*NGLLX+I]*temp3l);

  // adds gravity contribution
  if( gravity ) sum_terms += gravity_term;

  // assembles potential array
  if(threadIdx.x < NGLL3) {
#ifdef USE_MESH_COLORING_GPU
    // no atomic operation needed, colors don't share global points between elements
#ifdef USE_TEXTURES_FIELDS
    d_potential_dot_dot_acoustic[iglob] = texfetch_potential_dot_dot<FORWARD_OR_ADJOINT>(iglob) + sum_terms;
#else
    d_potential_dot_dot_acoustic[iglob] += sum_terms;
#endif // USE_TEXTURES_FIELDS

#else  // MESH_COLORING

    //mesh coloring
    if( use_mesh_coloring_gpu ){
      // no atomic operation needed, colors don't share global points between elements
#ifdef USE_TEXTURES_FIELDS
      d_potential_dot_dot_acoustic[iglob] = texfetch_potential_dot_dot<FORWARD_OR_ADJOINT>(iglob) + sum_terms;
#else
      d_potential_dot_dot_acoustic[iglob] += sum_terms;
#endif // USE_TEXTURES_FIELDS

    }else{
      atomicAdd(&d_potential_dot_dot_acoustic[iglob],sum_terms);
    }
#endif // MESH_COLORING
  } // threadIdx.x

}

/* ----------------------------------------------------------------------------------------------- */

/*

// original kernel
// please leave it here for reference...

//template<int FORWARD_OR_ADJOINT> __global__ void
//Kernel_2_acoustic_org_impl(int nb_blocks_to_compute,
                       const int NGLOB,
                       const int* d_ibool,
                       const int* d_phase_ispec_inner_acoustic,
                       const int num_phase_ispec_acoustic,
                       const int d_iphase,
                       const int use_mesh_coloring_gpu,
                       realw_const_p d_potential_acoustic,
                       realw_p d_potential_dot_dot_acoustic,
                       realw_const_p d_xix,realw_const_p d_xiy,realw_const_p d_xiz,
                       realw_const_p d_etax,realw_const_p d_etay,realw_const_p d_etaz,
                       realw_const_p d_gammax,realw_const_p d_gammay,realw_const_p d_gammaz,
                       realw_const_p d_hprime_xx,
                       realw_const_p hprimewgll_xx,
                       realw_const_p wgllwgll_xy,realw_const_p wgllwgll_xz,realw_const_p wgllwgll_yz,
                       realw_const_p d_rhostore,
                       const int gravity,
                       realw_const_p minus_g,
                       realw_const_p d_kappastore,
                       realw_const_p wgll_cube){

  int bx = blockIdx.y*gridDim.x+blockIdx.x;
  int tx = threadIdx.x;

  int K = (tx/NGLL2);
  int J = ((tx-K*NGLL2)/NGLLX);
  int I = (tx-K*NGLL2-J*NGLLX);

  unsigned short int active;
  int iglob,offset;
  int working_element;

  realw temp1l,temp2l,temp3l;
  realw xixl,xiyl,xizl,etaxl,etayl,etazl,gammaxl,gammayl,gammazl,jacobianl;

  realw dpotentialdxl,dpotentialdyl,dpotentialdzl;
  realw fac1,fac2,fac3;
  realw rho_invl,kappa_invl;

  realw sum_terms;
  realw gravity_term;

#ifndef MANUALLY_UNROLLED_LOOPS
  int l;
#endif

  __shared__ realw s_dummy_loc[NGLL3];

  __shared__ realw s_temp1[NGLL3];
  __shared__ realw s_temp2[NGLL3];
  __shared__ realw s_temp3[NGLL3];

  __shared__ realw sh_hprime_xx[NGLL2];

// use only NGLL^3 = 125 active threads, plus 3 inactive/ghost threads,
// because we used memory padding from NGLL^3 = 125 to 128 to get coalescent memory accesses
  active = (tx < NGLL3 && bx < nb_blocks_to_compute) ? 1:0;

// copy from global memory to shared memory
// each thread writes one of the NGLL^3 = 125 data points
  if (active) {

#ifdef USE_MESH_COLORING_GPU
    working_element = bx;
#else
    //mesh coloring
    if( use_mesh_coloring_gpu ){
      working_element = bx;
    }else{
      // iphase-1 and working_element-1 for Fortran->C array conventions
      working_element = d_phase_ispec_inner_acoustic[bx + num_phase_ispec_acoustic*(d_iphase-1)]-1;
    }
#endif
    // local padded index
    offset = working_element*NGLL3_PADDED + tx;

    // global index
    iglob = d_ibool[offset] - 1;

#ifdef USE_TEXTURES_FIELDS
    s_dummy_loc[tx] = texfetch_potential<FORWARD_OR_ADJOINT>(iglob);
#else
    // changing iglob indexing to match fortran row changes fast style
    s_dummy_loc[tx] = d_potential_acoustic[iglob];
#endif
  }

  if (tx < NGLL2) {
#ifdef USE_TEXTURES_CONSTANTS
    sh_hprime_xx[tx] = tex1Dfetch(d_hprime_xx_tex,tx);
#else
    sh_hprime_xx[tx] = d_hprime_xx[tx];
#endif
  }

// synchronize all the threads (one thread for each of the NGLL grid points of the
// current spectral element) because we need the whole element to be ready in order
// to be able to compute the matrix products along cut planes of the 3D element below
  __syncthreads();

  if (active) {

#ifndef MANUALLY_UNROLLED_LOOPS

    temp1l = 0.f;
    temp2l = 0.f;
    temp3l = 0.f;

    for (l=0;l<NGLLX;l++) {
      //assumes that hprime_xx = hprime_yy = hprime_zz
      fac1 = sh_hprime_xx[l*NGLLX+I];
      temp1l += s_dummy_loc[K*NGLL2+J*NGLLX+l]*fac1;

      fac2 = sh_hprime_xx[l*NGLLX+J];
      temp2l += s_dummy_loc[K*NGLL2+l*NGLLX+I]*fac2;

      fac3 = sh_hprime_xx[l*NGLLX+K];
      temp3l += s_dummy_loc[l*NGLL2+J*NGLLX+I]*fac3;
    }
#else
    temp1l = s_dummy_loc[K*NGLL2+J*NGLLX]*d_hprime_xx[I]
            + s_dummy_loc[K*NGLL2+J*NGLLX+1]*d_hprime_xx[NGLLX+I]
            + s_dummy_loc[K*NGLL2+J*NGLLX+2]*d_hprime_xx[2*NGLLX+I]
            + s_dummy_loc[K*NGLL2+J*NGLLX+3]*d_hprime_xx[3*NGLLX+I]
            + s_dummy_loc[K*NGLL2+J*NGLLX+4]*d_hprime_xx[4*NGLLX+I];

    temp2l = s_dummy_loc[K*NGLL2+I]*d_hprime_xx[J]
            + s_dummy_loc[K*NGLL2+NGLLX+I]*d_hprime_xx[NGLLX+J]
            + s_dummy_loc[K*NGLL2+2*NGLLX+I]*d_hprime_xx[2*NGLLX+J]
            + s_dummy_loc[K*NGLL2+3*NGLLX+I]*d_hprime_xx[3*NGLLX+J]
            + s_dummy_loc[K*NGLL2+4*NGLLX+I]*d_hprime_xx[4*NGLLX+J];

    temp3l = s_dummy_loc[J*NGLLX+I]*d_hprime_xx[K]
            + s_dummy_loc[NGLL2+J*NGLLX+I]*d_hprime_xx[NGLLX+K]
            + s_dummy_loc[2*NGLL2+J*NGLLX+I]*d_hprime_xx[2*NGLLX+K]
            + s_dummy_loc[3*NGLL2+J*NGLLX+I]*d_hprime_xx[3*NGLLX+K]
            + s_dummy_loc[4*NGLL2+J*NGLLX+I]*d_hprime_xx[4*NGLLX+K];
#endif

    // compute derivatives of ux, uy and uz with respect to x, y and z
    xixl = d_xix[offset];
    xiyl = d_xiy[offset];
    xizl = d_xiz[offset];
    etaxl = d_etax[offset];
    etayl = d_etay[offset];
    etazl = d_etaz[offset];
    gammaxl = d_gammax[offset];
    gammayl = d_gammay[offset];
    gammazl = d_gammaz[offset];

    jacobianl = 1.f / (xixl*(etayl*gammazl-etazl*gammayl)
                      -xiyl*(etaxl*gammazl-etazl*gammaxl)
                      +xizl*(etaxl*gammayl-etayl*gammaxl));

    // derivatives of potential
    dpotentialdxl = xixl*temp1l + etaxl*temp2l + gammaxl*temp3l;
    dpotentialdyl = xiyl*temp1l + etayl*temp2l + gammayl*temp3l;
    dpotentialdzl = xizl*temp1l + etazl*temp2l + gammazl*temp3l;

    // pre-computes gravity sum term
    if( gravity ){
      // uses potential definition: s = grad(chi)

      // gravity term: 1/kappa grad(chi) * g
      // assumes that g only acts in (negative) z-direction
      kappa_invl = 1.f / d_kappastore[working_element*NGLL3 + tx];

      gravity_term = minus_g[iglob] * kappa_invl * jacobianl * wgll_cube[tx] * dpotentialdzl;
    }

    // density (reciproc)
    rho_invl = 1.f / d_rhostore[offset];

    // form the dot product with the test vector
    s_temp1[tx] = jacobianl * rho_invl * (dpotentialdxl*xixl + dpotentialdyl*xiyl + dpotentialdzl*xizl);
    s_temp2[tx] = jacobianl * rho_invl * (dpotentialdxl*etaxl + dpotentialdyl*etayl + dpotentialdzl*etazl);
    s_temp3[tx] = jacobianl * rho_invl * (dpotentialdxl*gammaxl + dpotentialdyl*gammayl + dpotentialdzl*gammazl);
  }

// synchronize all the threads (one thread for each of the NGLL grid points of the
// current spectral element) because we need the whole element to be ready in order
// to be able to compute the matrix products along cut planes of the 3D element below
  __syncthreads();

  if (active) {

#ifndef MANUALLY_UNROLLED_LOOPS

    temp1l = 0.f;
    temp2l = 0.f;
    temp3l = 0.f;

    for (l=0;l<NGLLX;l++) {
      //assumes hprimewgll_xx = hprimewgll_yy = hprimewgll_zz
      fac1 = hprimewgll_xx[I*NGLLX+l];
      temp1l += s_temp1[K*NGLL2+J*NGLLX+l]*fac1;

      fac2 = hprimewgll_xx[J*NGLLX+l];
      temp2l += s_temp2[K*NGLL2+l*NGLLX+I]*fac2;

      fac3 = hprimewgll_xx[K*NGLLX+l];
      temp3l += s_temp3[l*NGLL2+J*NGLLX+I]*fac3;
    }
#else

    temp1l = s_temp1[K*NGLL2+J*NGLLX]*hprimewgll_xx[I*NGLLX]
            + s_temp1[K*NGLL2+J*NGLLX+1]*hprimewgll_xx[I*NGLLX+1]
            + s_temp1[K*NGLL2+J*NGLLX+2]*hprimewgll_xx[I*NGLLX+2]
            + s_temp1[K*NGLL2+J*NGLLX+3]*hprimewgll_xx[I*NGLLX+3]
            + s_temp1[K*NGLL2+J*NGLLX+4]*hprimewgll_xx[I*NGLLX+4];


    temp2l = s_temp2[K*NGLL2+I]*hprimewgll_xx[J*NGLLX]
            + s_temp2[K*NGLL2+NGLLX+I]*hprimewgll_xx[J*NGLLX+1]
            + s_temp2[K*NGLL2+2*NGLLX+I]*hprimewgll_xx[J*NGLLX+2]
            + s_temp2[K*NGLL2+3*NGLLX+I]*hprimewgll_xx[J*NGLLX+3]
            + s_temp2[K*NGLL2+4*NGLLX+I]*hprimewgll_xx[J*NGLLX+4];


    temp3l = s_temp3[J*NGLLX+I]*hprimewgll_xx[K*NGLLX]
            + s_temp3[NGLL2+J*NGLLX+I]*hprimewgll_xx[K*NGLLX+1]
            + s_temp3[2*NGLL2+J*NGLLX+I]*hprimewgll_xx[K*NGLLX+2]
            + s_temp3[3*NGLL2+J*NGLLX+I]*hprimewgll_xx[K*NGLLX+3]
            + s_temp3[4*NGLL2+J*NGLLX+I]*hprimewgll_xx[K*NGLLX+4];


#endif

    fac1 = wgllwgll_yz[K*NGLLX+J];
    fac2 = wgllwgll_xz[K*NGLLX+I];
    fac3 = wgllwgll_xy[J*NGLLX+I];

    sum_terms = -(fac1*temp1l + fac2*temp2l + fac3*temp3l);
    if( gravity ) sum_terms += gravity_term;

#ifdef USE_MESH_COLORING_GPU
    // no atomic operation needed, colors don't share global points between elements

#ifdef USE_TEXTURES_FIELDS
    d_potential_dot_dot_acoustic[iglob] = texfetch_potential_dot_dot<FORWARD_OR_ADJOINT>(iglob) + sum_terms;
#else
    d_potential_dot_dot_acoustic[iglob] += sum_terms;
#endif // USE_TEXTURES_FIELDS


#else  // MESH_COLORING

    //mesh coloring
    if( use_mesh_coloring_gpu ){
      // no atomic operation needed, colors don't share global points between elements
#ifdef USE_TEXTURES_FIELDS
      d_potential_dot_dot_acoustic[iglob] = texfetch_potential_dot_dot<FORWARD_OR_ADJOINT>(iglob) + sum_terms;
#else
      d_potential_dot_dot_acoustic[iglob] += sum_terms;
#endif // USE_TEXTURES_FIELDS

    }else{
      atomicAdd(&d_potential_dot_dot_acoustic[iglob],sum_terms);
    }
#endif // MESH_COLORING

  }
}
*/



/* ----------------------------------------------------------------------------------------------- */

void Kernel_2_acoustic(int nb_blocks_to_compute, Mesh* mp, int d_iphase,
                       int* d_ibool,
                       realw* d_xix,realw* d_xiy,realw* d_xiz,
                       realw* d_etax,realw* d_etay,realw* d_etaz,
                       realw* d_gammax,realw* d_gammay,realw* d_gammaz,
                       realw* d_rhostore,
                       realw* d_kappastore){

#ifdef ENABLE_VERY_SLOW_ERROR_CHECKING
  exit_on_cuda_error("before acoustic kernel Kernel 2");
#endif

  // if the grid can handle the number of blocks, we let it be 1D
  // grid_2_x = nb_elem_color;
  // nb_elem_color is just how many blocks we are computing now

  int blocksize = NGLL3_PADDED;

  int num_blocks_x, num_blocks_y;
  get_blocks_xy(nb_blocks_to_compute,&num_blocks_x,&num_blocks_y);

  dim3 grid(num_blocks_x,num_blocks_y);
  dim3 threads(blocksize,1,1);

  // Cuda timing
  //cudaEvent_t start, stop;
  //start_timing_cuda(&start,&stop);

  if( ! mp->gravity ) {
    // without gravity

    if( ! mp->use_mesh_coloring_gpu ){
      // without mesh coloring
      // forward wavefields -> FORWARD_OR_ADJOINT == 1
      Kernel_2_acoustic_impl<1><<<grid,threads,0,mp->compute_stream>>>(nb_blocks_to_compute,
                                                                mp->NGLOB_AB,
                                                                d_ibool,
                                                                mp->d_phase_ispec_inner_acoustic,
                                                                mp->num_phase_ispec_acoustic,
                                                                d_iphase,
                                                                mp->d_potential_acoustic, mp->d_potential_dot_dot_acoustic,
                                                                d_xix, d_xiy, d_xiz,
                                                                d_etax, d_etay, d_etaz,
                                                                d_gammax, d_gammay, d_gammaz,
                                                                mp->d_hprime_xx,
                                                                mp->d_hprimewgll_xx,
                                                                mp->d_wgllwgll_xy, mp->d_wgllwgll_xz, mp->d_wgllwgll_yz,
                                                                d_rhostore);

      if(mp->simulation_type == 3) {
        // backward/reconstructed wavefields -> FORWARD_OR_ADJOINT == 3
        Kernel_2_acoustic_impl<3><<<grid,threads,0,mp->compute_stream>>>(nb_blocks_to_compute,
                                                              mp->NGLOB_AB,
                                                              d_ibool,
                                                              mp->d_phase_ispec_inner_acoustic,
                                                              mp->num_phase_ispec_acoustic,
                                                              d_iphase,
                                                              mp->d_b_potential_acoustic, mp->d_b_potential_dot_dot_acoustic,
                                                              d_xix, d_xiy, d_xiz,
                                                              d_etax, d_etay, d_etaz,
                                                              d_gammax, d_gammay, d_gammaz,
                                                              mp->d_hprime_xx,
                                                              mp->d_hprimewgll_xx,
                                                              mp->d_wgllwgll_xy, mp->d_wgllwgll_xz, mp->d_wgllwgll_yz,
                                                              d_rhostore);
      }
    }else{
      // with mesh coloring
      // forward wavefields -> FORWARD_OR_ADJOINT == 1
      Kernel_2_acoustic_col_impl<1><<<grid,threads,0,mp->compute_stream>>>(nb_blocks_to_compute,
                                                                mp->NGLOB_AB,
                                                                d_ibool,
                                                                mp->d_phase_ispec_inner_acoustic,
                                                                mp->num_phase_ispec_acoustic,
                                                                d_iphase,
                                                                mp->d_potential_acoustic, mp->d_potential_dot_dot_acoustic,
                                                                d_xix, d_xiy, d_xiz,
                                                                d_etax, d_etay, d_etaz,
                                                                d_gammax, d_gammay, d_gammaz,
                                                                mp->d_hprime_xx,
                                                                mp->d_hprimewgll_xx,
                                                                mp->d_wgllwgll_xy, mp->d_wgllwgll_xz, mp->d_wgllwgll_yz,
                                                                d_rhostore,
                                                                mp->use_mesh_coloring_gpu);

      if(mp->simulation_type == 3) {
        // backward/reconstructed wavefields -> FORWARD_OR_ADJOINT == 3
        Kernel_2_acoustic_col_impl<3><<<grid,threads,0,mp->compute_stream>>>(nb_blocks_to_compute,
                                                              mp->NGLOB_AB,
                                                              d_ibool,
                                                              mp->d_phase_ispec_inner_acoustic,
                                                              mp->num_phase_ispec_acoustic,
                                                              d_iphase,
                                                              mp->d_b_potential_acoustic, mp->d_b_potential_dot_dot_acoustic,
                                                              d_xix, d_xiy, d_xiz,
                                                              d_etax, d_etay, d_etaz,
                                                              d_gammax, d_gammay, d_gammaz,
                                                              mp->d_hprime_xx,
                                                              mp->d_hprimewgll_xx,
                                                              mp->d_wgllwgll_xy, mp->d_wgllwgll_xz, mp->d_wgllwgll_yz,
                                                              d_rhostore,
                                                              mp->use_mesh_coloring_gpu);
      }
    } // use_mesh_coloring_gpu
    
  }else{
    // with gravity
    
    // forward wavefields -> FORWARD_OR_ADJOINT == 1
    Kernel_2_acoustic_grav_impl<1><<<grid,threads,0,mp->compute_stream>>>(nb_blocks_to_compute,
                                                              mp->NGLOB_AB,
                                                              d_ibool,
                                                              mp->d_phase_ispec_inner_acoustic,
                                                              mp->num_phase_ispec_acoustic,
                                                              d_iphase,
                                                              mp->d_potential_acoustic, mp->d_potential_dot_dot_acoustic,
                                                              d_xix, d_xiy, d_xiz,
                                                              d_etax, d_etay, d_etaz,
                                                              d_gammax, d_gammay, d_gammaz,
                                                              mp->d_hprime_xx,
                                                              mp->d_hprimewgll_xx,
                                                              mp->d_wgllwgll_xy, mp->d_wgllwgll_xz, mp->d_wgllwgll_yz,
                                                              d_rhostore,
                                                              mp->use_mesh_coloring_gpu,
                                                              mp->gravity,
                                                              mp->d_minus_g,
                                                              d_kappastore,
                                                              mp->d_wgll_cube);

    if(mp->simulation_type == 3) {
      // backward/reconstructed wavefields -> FORWARD_OR_ADJOINT == 3
      Kernel_2_acoustic_grav_impl<3><<<grid,threads,0,mp->compute_stream>>>(nb_blocks_to_compute,
                                                            mp->NGLOB_AB,
                                                            d_ibool,
                                                            mp->d_phase_ispec_inner_acoustic,
                                                            mp->num_phase_ispec_acoustic,
                                                            d_iphase,
                                                            mp->d_b_potential_acoustic, mp->d_b_potential_dot_dot_acoustic,
                                                            d_xix, d_xiy, d_xiz,
                                                            d_etax, d_etay, d_etaz,
                                                            d_gammax, d_gammay, d_gammaz,
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
  } // gravity

  // Cuda timing
  //if( ! mp->gravity ) {
  //  if( ! mp->use_mesh_coloring_gpu ){
  //    stop_timing_cuda(&start,&stop,"Kernel_2_acoustic_impl");
  //  }else{
  //    stop_timing_cuda(&start,&stop,"Kernel_2_acoustic_col_impl");
  //  }
  //}else{
  //  stop_timing_cuda(&start,&stop,"Kernel_2_acoustic_grav_impl");
  //}
  
#ifdef ENABLE_VERY_SLOW_ERROR_CHECKING
  exit_on_cuda_error("kernel Kernel_2");
#endif
}

/* ----------------------------------------------------------------------------------------------- */

// main compute_forces_acoustic CUDA routine

/* ----------------------------------------------------------------------------------------------- */

extern "C"
void FC_FUNC_(compute_forces_acoustic_cuda,
              COMPUTE_FORCES_ACOUSTIC_CUDA)(long* Mesh_pointer,
                                            int* iphase,
                                            int* nspec_outer_acoustic,
                                            int* nspec_inner_acoustic) {

  TRACE("compute_forces_acoustic_cuda");
  //double start_time = get_time();

  Mesh* mp = (Mesh*)(*Mesh_pointer); // get Mesh from fortran integer wrapper

  int num_elements;

  if( *iphase == 1 )
    num_elements = *nspec_outer_acoustic;
  else
    num_elements = *nspec_inner_acoustic;

  if( num_elements == 0 ) return;

  // mesh coloring
  if( mp->use_mesh_coloring_gpu ){

    // note: array offsets require sorted arrays, such that e.g. ibool starts with elastic elements
    //         and followed by acoustic ones.
    //         acoustic elements also start with outer than inner element ordering

    int nb_colors,nb_blocks_to_compute;
    int istart;
    int offset,offset_nonpadded;

    // sets up color loop
    if( *iphase == 1 ){
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
                         mp->d_ibool + offset_nonpadded,
                         mp->d_xix + offset,mp->d_xiy + offset,mp->d_xiz + offset,
                         mp->d_etax + offset,mp->d_etay + offset,mp->d_etaz + offset,
                         mp->d_gammax + offset,mp->d_gammay + offset,mp->d_gammaz + offset,
                         mp->d_rhostore + offset,
                         mp->d_kappastore + offset_nonpadded);

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
                      mp->d_kappastore);

  }
}



/* ----------------------------------------------------------------------------------------------- */

/* KERNEL for enforce free surface */

/* ----------------------------------------------------------------------------------------------- */


__global__ void enforce_free_surface_cuda_kernel(
                                       realw_p potential_acoustic,
                                       realw_p potential_dot_acoustic,
                                       realw_p potential_dot_dot_acoustic,
                                       const int num_free_surface_faces,
                                       const int* free_surface_ispec,
                                       const int* free_surface_ijk,
                                       const int* d_ibool,
                                       const int* ispec_is_acoustic) {
  // gets spectral element face id
  int iface = blockIdx.x + gridDim.x*blockIdx.y;

  // for all faces on free surface
  if( iface < num_free_surface_faces ){

    int ispec = free_surface_ispec[iface]-1;

    // checks if element is in acoustic domain
    if( ispec_is_acoustic[ispec] ){

      // gets global point index
      int igll = threadIdx.x + threadIdx.y*blockDim.x;

      int i = free_surface_ijk[INDEX3(NDIM,NGLL2,0,igll,iface)] - 1; // (1,igll,iface)
      int j = free_surface_ijk[INDEX3(NDIM,NGLL2,1,igll,iface)] - 1;
      int k = free_surface_ijk[INDEX3(NDIM,NGLL2,2,igll,iface)] - 1;

      int iglob = d_ibool[INDEX4_PADDED(NGLLX,NGLLX,NGLLX,i,j,k,ispec)] - 1;

      // sets potentials to zero at free surface
      potential_acoustic[iglob] = 0;
      potential_dot_acoustic[iglob] = 0;
      potential_dot_dot_acoustic[iglob] = 0;
    }
  }
}

/* ----------------------------------------------------------------------------------------------- */

extern "C"
void FC_FUNC_(acoustic_enforce_free_surf_cuda,
              ACOUSTIC_ENFORCE_FREE_SURF_CUDA)(long* Mesh_pointer,
                                               int* ABSORB_INSTEAD_OF_FREE_SURFACE) {

TRACE("acoustic_enforce_free_surf_cuda");

  Mesh* mp = (Mesh*)(*Mesh_pointer); //get mesh pointer out of fortran integer container

  // checks if anything to do
  if( *ABSORB_INSTEAD_OF_FREE_SURFACE == 0 ){

    // does not absorb free surface, thus we enforce the potential to be zero at surface

    // block sizes
    int num_blocks_x, num_blocks_y;
    get_blocks_xy(mp->num_free_surface_faces,&num_blocks_x,&num_blocks_y);

    dim3 grid(num_blocks_x,num_blocks_y,1);
    dim3 threads(NGLL2,1,1);

    // sets potentials to zero at free surface
    enforce_free_surface_cuda_kernel<<<grid,threads,0,mp->compute_stream>>>(mp->d_potential_acoustic,
                                                                             mp->d_potential_dot_acoustic,
                                                                             mp->d_potential_dot_dot_acoustic,
                                                                             mp->num_free_surface_faces,
                                                                             mp->d_free_surface_ispec,
                                                                             mp->d_free_surface_ijk,
                                                                             mp->d_ibool,
                                                                             mp->d_ispec_is_acoustic);
    // for backward/reconstructed potentials
    if(mp->simulation_type == 3) {
      enforce_free_surface_cuda_kernel<<<grid,threads,0,mp->compute_stream>>>(mp->d_b_potential_acoustic,
                                                                               mp->d_b_potential_dot_acoustic,
                                                                               mp->d_b_potential_dot_dot_acoustic,
                                                                               mp->num_free_surface_faces,
                                                                               mp->d_free_surface_ispec,
                                                                               mp->d_free_surface_ijk,
                                                                               mp->d_ibool,
                                                                               mp->d_ispec_is_acoustic);
    }
  }

#ifdef ENABLE_VERY_SLOW_ERROR_CHECKING
  exit_on_cuda_error("enforce_free_surface_cuda");
#endif
}

