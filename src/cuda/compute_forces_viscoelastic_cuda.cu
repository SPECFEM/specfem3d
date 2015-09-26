/*
 !=====================================================================
 !
 !               S p e c f e m 3 D  V e r s i o n  3 . 0
 !               ---------------------------------------
 !
 !     Main historical authors: Dimitri Komatitsch and Jeroen Tromp
 !                        Princeton University, USA
 !                and CNRS / University of Marseille, France
 !                 (there are currently many more authors!)
 ! (c) Princeton University and CNRS / University of Marseille, July 2012
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

#include "mesh_constants_cuda.h"

/* ----------------------------------------------------------------------------------------------- */

#ifdef USE_TEXTURES_FIELDS
realw_texture d_displ_tex;
realw_texture d_veloc_tex;
realw_texture d_accel_tex;
//backward/reconstructed
realw_texture d_b_displ_tex;
realw_texture d_b_veloc_tex;
realw_texture d_b_accel_tex;

//note: texture variables are implicitly static, and cannot be passed as arguments to cuda kernels;
//      thus, 1) we thus use if-statements (FORWARD_OR_ADJOINT) to determine from which texture to fetch from
//            2) we use templates
//      since if-statements are a bit slower as the variable is only known at runtime, we use option 2)

// templates definitions
template<int FORWARD_OR_ADJOINT> __device__ float texfetch_displ(int x);
template<int FORWARD_OR_ADJOINT> __device__ float texfetch_veloc(int x);
template<int FORWARD_OR_ADJOINT> __device__ float texfetch_accel(int x);

// templates for texture fetching
// FORWARD_OR_ADJOINT == 1 <- forward arrays
template<> __device__ float texfetch_displ<1>(int x) { return tex1Dfetch(d_displ_tex, x); }
template<> __device__ float texfetch_veloc<1>(int x) { return tex1Dfetch(d_veloc_tex, x); }
template<> __device__ float texfetch_accel<1>(int x) { return tex1Dfetch(d_accel_tex, x); }
// FORWARD_OR_ADJOINT == 3 <- backward/reconstructed arrays
template<> __device__ float texfetch_displ<3>(int x) { return tex1Dfetch(d_b_displ_tex, x); }
template<> __device__ float texfetch_veloc<3>(int x) { return tex1Dfetch(d_b_veloc_tex, x); }
template<> __device__ float texfetch_accel<3>(int x) { return tex1Dfetch(d_b_accel_tex, x); }

#endif

#ifdef USE_TEXTURES_CONSTANTS
realw_texture d_hprime_xx_tex;
#endif


/* ----------------------------------------------------------------------------------------------- */

// KERNEL 2

/* ----------------------------------------------------------------------------------------------- */


// updates stress

__device__ __forceinline__ void compute_element_att_stress(int tx,int working_element,const int NSPEC,
                                           realw* R_xx,realw* R_yy,realw* R_xy,
                                           realw* R_xz,realw* R_yz,
                                           realw* sigma_xx,realw* sigma_yy,realw* sigma_zz,
                                           realw* sigma_xy,realw* sigma_xz,realw* sigma_yz) {

  int offset_sls;
  realw R_xx_val,R_yy_val;

  for(int i_sls = 0; i_sls < N_SLS; i_sls++){
    // index
    offset_sls = tx + NGLL3*(working_element + NSPEC*i_sls);

    R_xx_val = get_global_cr( &R_xx[offset_sls] ); //(i,j,k,ispec,i_sls)
    R_yy_val = get_global_cr( &R_yy[offset_sls] );

    *sigma_xx = *sigma_xx - R_xx_val;
    *sigma_yy = *sigma_yy - R_yy_val;
    *sigma_zz = *sigma_zz + R_xx_val + R_yy_val;
    *sigma_xy = *sigma_xy - get_global_cr( &R_xy[offset_sls] );
    *sigma_xz = *sigma_xz - get_global_cr( &R_xz[offset_sls] );
    *sigma_yz = *sigma_yz - get_global_cr( &R_yz[offset_sls] );
  }
  return;
}

/* ----------------------------------------------------------------------------------------------- */

// updates R_memory

__device__  __forceinline__ void compute_element_att_memory(int tx,int working_element,const int NSPEC,
                                          realw_const_p d_muv,
                                          realw_const_p factor_common,
                                          realw_const_p alphaval,realw_const_p betaval,realw_const_p gammaval,
                                          realw_p R_xx,realw_p R_yy,realw_p R_xy,realw_p R_xz,realw_p R_yz,
                                          realw_p epsilondev_xx,realw_p epsilondev_yy,realw_p epsilondev_xy,
                                          realw_p epsilondev_xz,realw_p epsilondev_yz,
                                          realw epsilondev_xx_loc,realw epsilondev_yy_loc,realw epsilondev_xy_loc,
                                          realw epsilondev_xz_loc,realw epsilondev_yz_loc
                                          ){

  int ijk_ispec;
  int offset_sls,offset_align,offset_common;
  realw mul;
  realw alphaval_loc,betaval_loc,gammaval_loc;
  realw factor_loc;

  realw rxxl,ryyl,rxyl,rxzl,ryzl;
  realw Sn_xx,Sn_yy,Sn_xy,Sn_xz,Sn_yz;
  realw Snp1_xx,Snp1_yy,Snp1_xy,Snp1_xz,Snp1_yz;

  // indices
  offset_align = tx + NGLL3_PADDED * working_element;
  ijk_ispec = tx + NGLL3 * working_element;

  mul = get_global_cr( &d_muv[offset_align] );

  // use Runge-Kutta scheme to march in time
  for(int i_sls = 0; i_sls < N_SLS; i_sls++){

    // indices
    offset_common = i_sls + N_SLS*(tx + NGLL3*working_element); // (i_sls,i,j,k,ispec)
    offset_sls = tx + NGLL3*(working_element + NSPEC*i_sls);   // (i,j,k,ispec,i_sls)

    factor_loc = mul * get_global_cr( &factor_common[offset_common] ); //mustore(i,j,k,ispec) * factor_common(i_sls,i,j,k,ispec)

    alphaval_loc = alphaval[i_sls]; // (i_sls)
    betaval_loc = betaval[i_sls];
    gammaval_loc = gammaval[i_sls];

    // term in xx
    Sn_xx   = factor_loc * get_global_cr( &epsilondev_xx[ijk_ispec] ); //(i,j,k,ispec)
    Snp1_xx   = factor_loc * epsilondev_xx_loc; //(i,j,k)

    rxxl = get_global_cr( &R_xx[offset_sls] );
    R_xx[offset_sls] = alphaval_loc * rxxl + betaval_loc * Sn_xx + gammaval_loc * Snp1_xx;

    // term in yy
    Sn_yy   = factor_loc * get_global_cr( &epsilondev_yy[ijk_ispec] );
    Snp1_yy   = factor_loc * epsilondev_yy_loc;

    ryyl = get_global_cr( &R_yy[offset_sls] );
    R_yy[offset_sls] = alphaval_loc * ryyl + betaval_loc * Sn_yy + gammaval_loc * Snp1_yy;

    // term in zz not computed since zero trace
    // term in xy
    Sn_xy   = factor_loc * get_global_cr( &epsilondev_xy[ijk_ispec] );
    Snp1_xy   = factor_loc * epsilondev_xy_loc;

    rxyl = get_global_cr( &R_xy[offset_sls] );
    R_xy[offset_sls] = alphaval_loc * rxyl + betaval_loc * Sn_xy + gammaval_loc * Snp1_xy;

    // term in xz
    Sn_xz   = factor_loc * get_global_cr( &epsilondev_xz[ijk_ispec] );
    Snp1_xz   = factor_loc * epsilondev_xz_loc;

    rxzl = get_global_cr( &R_xz[offset_sls] );
    R_xz[offset_sls] = alphaval_loc * rxzl + betaval_loc * Sn_xz + gammaval_loc * Snp1_xz;

    // term in yz
    Sn_yz   = factor_loc * get_global_cr( &epsilondev_yz[ijk_ispec] );
    Snp1_yz   = factor_loc * epsilondev_yz_loc;

    ryzl = get_global_cr( &R_yz[offset_sls] );
    R_yz[offset_sls] = alphaval_loc * ryzl + betaval_loc * Sn_yz + gammaval_loc * Snp1_yz;

  }
  return;
}

/* ----------------------------------------------------------------------------------------------- */

// pre-computes gravity term

__device__ __forceinline__ void compute_element_gravity(int tx,int working_element,
                                        const int* iglob,
                                        realw_const_p d_minus_g,
                                        realw_const_p d_minus_deriv_gravity,
                                        realw_const_p d_rhostore,
                                        realw_const_p wgll_cube,
                                        realw jacobianl,
                                        realw* sh_tempx,
                                        realw* sh_tempy,
                                        realw* sh_tempz,
                                        realw* sigma_xx,
                                        realw* sigma_yy,
                                        realw* sigma_xz,
                                        realw* sigma_yz,
                                        realw* rho_s_H1,
                                        realw* rho_s_H2,
                                        realw* rho_s_H3){

  realw minus_g,minus_dg;
  realw rhol;
  realw gzl; // gxl,gyl,
  realw sx_l,sy_l,sz_l;
  realw Hxxl,Hyyl,Hzzl; //,Hxyl,Hxzl,Hyzl;
  realw factor;

  // compute non-symmetric terms for gravity

  // get g, rho and dg/dr=dg
  minus_g = d_minus_g[*iglob];
  minus_dg = d_minus_deriv_gravity[*iglob];

  // Cartesian components of the gravitational acceleration
  //gxl = 0.f;
  //gyl = 0.f;
  gzl = minus_g;

  // Cartesian components of gradient of gravitational acceleration
  // H = grad g
  // assumes g only acts in negative z-direction
  Hxxl = 0.f;
  Hyyl = 0.f;
  Hzzl = minus_dg;
  //Hxyl = 0.f;
  //Hxzl = 0.f;
  //Hyzl = 0.f;

  rhol = get_global_cr( &d_rhostore[working_element*NGLL3_PADDED + tx] );

  // get displacement and multiply by density to compute G tensor
  // G = rho [ sg - (s * g) I  ]
  sx_l = rhol * sh_tempx[tx]; // d_displ[iglob*3];
  sy_l = rhol * sh_tempy[tx]; // d_displ[iglob*3 + 1];
  sz_l = rhol * sh_tempz[tx]; // d_displ[iglob*3 + 2];

  // compute G tensor from s . g and add to sigma (not symmetric)
  //sigma_xx += sy_l*gyl + sz_l*gzl;
  *sigma_xx += sz_l*gzl;
  //sigma_yy += sx_l*gxl + sz_l*gzl;
  *sigma_yy += sz_l*gzl;
  //sigma_zz += sx_l*gxl + sy_l*gyl;

  //sigma_xy -= sx_l*gyl;
  //sigma_yx -= sy_l*gxl;

  *sigma_xz -= sx_l*gzl;
  //sigma_zx -= sz_l*gxl;

  *sigma_yz -= sy_l*gzl;
  //sigma_zy -= sz_l*gyl;

  // precompute vector
  factor = jacobianl * wgll_cube[tx];

  //rho_s_H1 = fac1 * (sx_l * Hxxl + sy_l * Hxyl + sz_l * Hxzl);
  //rho_s_H2 = fac1 * (sx_l * Hxyl + sy_l * Hyyl + sz_l * Hyzl);
  //rho_s_H3 = fac1 * (sx_l * Hxzl + sy_l * Hyzl + sz_l * Hzzl);

  // only non-zero z-direction
  *rho_s_H1 = factor * sx_l * Hxxl ; // 0.f;
  *rho_s_H2 = factor * sy_l * Hyyl ; // 0.f;
  *rho_s_H3 = factor * sz_l * Hzzl ;

  // debug
  //*rho_s_H1 = 0.f;
  //*rho_s_H2 = 0.f;
  //*rho_s_H3 = 0.f ;

}

/* ----------------------------------------------------------------------------------------------- */

// loads displacement into shared memory for element

template<int FORWARD_OR_ADJOINT>
__device__  __forceinline__ void load_shared_memory_displ(const int* tx, const int* iglob,
                                                          realw_const_p d_displ,
                                                          realw* sh_tempx,
                                                          realw* sh_tempy,
                                                          realw* sh_tempz){

  // copy from global memory to shared memory
  // each thread writes one of the NGLL^3 = 125 data points
#ifdef USE_TEXTURES_FIELDS
  sh_tempx[(*tx)] = texfetch_displ<FORWARD_OR_ADJOINT>((*iglob)*3);
  sh_tempy[(*tx)] = texfetch_displ<FORWARD_OR_ADJOINT>((*iglob)*3 + 1);
  sh_tempz[(*tx)] = texfetch_displ<FORWARD_OR_ADJOINT>((*iglob)*3 + 2);
#else
  // changing iglob indexing to match fortran row changes fast style
  sh_tempx[(*tx)] = d_displ[(*iglob)*3];
  sh_tempy[(*tx)] = d_displ[(*iglob)*3 + 1];
  sh_tempz[(*tx)] = d_displ[(*iglob)*3 + 2];
#endif
}

/* ----------------------------------------------------------------------------------------------- */

// loads displacement + viscosity * velocity into shared memory for element

template<int FORWARD_OR_ADJOINT>
__device__  __forceinline__ void load_shared_memory_displ_visco(const int* tx, const int* iglob,
                                                          realw_const_p d_displ,
                                                          realw_const_p d_veloc,
                                                          realw visco,
                                                          realw* sh_tempx,
                                                          realw* sh_tempy,
                                                          realw* sh_tempz){

  // copy from global memory to shared memory
  // each thread writes one of the NGLL^3 = 125 data points
#ifdef USE_TEXTURES_FIELDS
  sh_tempx[(*tx)] = texfetch_displ<FORWARD_OR_ADJOINT>((*iglob)*3) + visco*texfetch_veloc<FORWARD_OR_ADJOINT>((*iglob)*3);
  sh_tempy[(*tx)] = texfetch_displ<FORWARD_OR_ADJOINT>((*iglob)*3 + 1) + visco*texfetch_veloc<FORWARD_OR_ADJOINT>((*iglob)*3 + 1);
  sh_tempz[(*tx)] = texfetch_displ<FORWARD_OR_ADJOINT>((*iglob)*3 + 2) + visco*texfetch_veloc<FORWARD_OR_ADJOINT>((*iglob)*3 + 2);
#else
  // changing iglob indexing to match fortran row changes fast style
  sh_tempx[(*tx)] = d_displ[(*iglob)*3] + visco * d_veloc[(*iglob)*3];
  sh_tempy[(*tx)] = d_displ[(*iglob)*3 + 1] + visco * d_veloc[(*iglob)*3 + 1];
  sh_tempz[(*tx)] = d_displ[(*iglob)*3 + 2] + visco * d_veloc[(*iglob)*3 + 2];
#endif
}

/* ----------------------------------------------------------------------------------------------- */

// loads hprime into shared memory for element

__device__  __forceinline__ void load_shared_memory_hprime(const int* tx,
                                                           realw_const_p d_hprime_xx,
                                                           realw* sh_hprime_xx){

  // each thread reads its corresponding value
  // (might be faster sometimes...)
#ifdef USE_TEXTURES_CONSTANTS
  // hprime
  sh_hprime_xx[(*tx)] = tex1Dfetch(d_hprime_xx_tex,tx + d_hprime_xx_tex_offset);
#else
  // hprime
  sh_hprime_xx[(*tx)] = d_hprime_xx[(*tx)];
#endif
}


/* ----------------------------------------------------------------------------------------------- */

// loads hprimewgll into shared memory for element

__device__  __forceinline__ void load_shared_memory_hprimewgll(const int* tx,
                                                               realw_const_p d_hprimewgll_xx,
                                                               realw* sh_hprimewgll_xx ){

  // each thread reads its corresponding value
  // weighted hprime
  sh_hprimewgll_xx[(*tx)] = d_hprimewgll_xx[(*tx)];
}

/* ----------------------------------------------------------------------------------------------- */

// computes a 3D matrix-vector product along a 2D cut-plane

__device__  __forceinline__ void sum_hprime_xi(int I, int J, int K,
                                              realw* tempxl,realw* tempyl,realw* tempzl,
                                              realw* sh_tempx,realw* sh_tempy,realw* sh_tempz, realw* sh_hprime ){

  realw fac;

  // initializes
  realw sumx = 0.f;
  realw sumy = 0.f;
  realw sumz = 0.f;

  // 1. cut-plane along xi-direction
  #pragma unroll
  for (int l=0;l<NGLLX;l++) {
    fac = sh_hprime[l*NGLLX+I];

    sumx += sh_tempx[K*NGLL2+J*NGLLX+l] * fac;
    sumy += sh_tempy[K*NGLL2+J*NGLLX+l] * fac;
    sumz += sh_tempz[K*NGLL2+J*NGLLX+l] * fac;
  }

// counts:
// + NGLLX * ( 2 + 3*6) FLOP = 100 FLOP
//
// + 0 BYTE

  *tempxl = sumx;
  *tempyl = sumy;
  *tempzl = sumz;
}

/* ----------------------------------------------------------------------------------------------- */

// computes a 3D matrix-vector product along a 2D cut-plane

__device__  __forceinline__ void sum_hprime_eta(int I, int J, int K,
                                               realw* tempxl,realw* tempyl,realw* tempzl,
                                               realw* sh_tempx,realw* sh_tempy,realw* sh_tempz, realw* sh_hprime ){

  realw fac;

  // initializes
  realw sumx = 0.f;
  realw sumy = 0.f;
  realw sumz = 0.f;

  // 2. cut-plane along eta-direction
  #pragma unroll
  for (int l=0;l<NGLLX;l++) {
    fac = sh_hprime[l*NGLLX+J];

    sumx += sh_tempx[K*NGLL2+l*NGLLX+I] * fac;
    sumy += sh_tempy[K*NGLL2+l*NGLLX+I] * fac;
    sumz += sh_tempz[K*NGLL2+l*NGLLX+I] * fac;
  }

  *tempxl = sumx;
  *tempyl = sumy;
  *tempzl = sumz;
}

/* ----------------------------------------------------------------------------------------------- */

// computes a 3D matrix-vector product along a 2D cut-plane

__device__  __forceinline__ void sum_hprime_gamma(int I, int J, int K,
                                                 realw* tempxl,realw* tempyl,realw* tempzl,
                                                 realw* sh_tempx,realw* sh_tempy,realw* sh_tempz, realw* sh_hprime ){

  realw fac;

  // initializes
  realw sumx = 0.f;
  realw sumy = 0.f;
  realw sumz = 0.f;

  // 3. cut-plane along gamma-direction
  #pragma unroll
  for (int l=0;l<NGLLX;l++) {
    fac = sh_hprime[l*NGLLX+K];

    sumx += sh_tempx[l*NGLL2+J*NGLLX+I] * fac;
    sumy += sh_tempy[l*NGLL2+J*NGLLX+I] * fac;
    sumz += sh_tempz[l*NGLL2+J*NGLLX+I] * fac;
  }

  *tempxl = sumx;
  *tempyl = sumy;
  *tempzl = sumz;
}

/* ----------------------------------------------------------------------------------------------- */

// computes a 3D matrix-vector product along a 2D cut-plane

__device__  __forceinline__ void sum_hprimewgll_xi(int I, int J, int K,
                                                   realw* tempxl,realw* tempyl,realw* tempzl,
                                                   realw* sh_tempx,realw* sh_tempy,realw* sh_tempz, realw* sh_hprimewgll ){

  realw fac;

  // initializes
  realw sumx = 0.f;
  realw sumy = 0.f;
  realw sumz = 0.f;

  // 1. cut-plane along xi-direction
  #pragma unroll
  for (int l=0;l<NGLLX;l++) {
    fac = sh_hprimewgll[I*NGLLX+l]; //  d_hprimewgll_xx[I*NGLLX+l];

    sumx += sh_tempx[K*NGLL2+J*NGLLX+l] * fac;
    sumy += sh_tempy[K*NGLL2+J*NGLLX+l] * fac;
    sumz += sh_tempz[K*NGLL2+J*NGLLX+l] * fac;
  }

  *tempxl = sumx;
  *tempyl = sumy;
  *tempzl = sumz;
}

/* ----------------------------------------------------------------------------------------------- */

// computes a 3D matrix-vector product along a 2D cut-plane

__device__  __forceinline__ void sum_hprimewgll_eta(int I, int J, int K,
                                               realw* tempxl,realw* tempyl,realw* tempzl,
                                               realw* sh_tempx,realw* sh_tempy,realw* sh_tempz, realw* sh_hprimewgll ){

  realw fac;

  // initializes
  realw sumx = 0.f;
  realw sumy = 0.f;
  realw sumz = 0.f;

  // 2. cut-plane along eta-direction
  #pragma unroll
  for (int l=0;l<NGLLX;l++) {
    fac = sh_hprimewgll[J*NGLLX+l]; // d_hprimewgll_xx[J*NGLLX+l];

    sumx += sh_tempx[K*NGLL2+l*NGLLX+I] * fac;
    sumy += sh_tempy[K*NGLL2+l*NGLLX+I] * fac;
    sumz += sh_tempz[K*NGLL2+l*NGLLX+I] * fac;
  }

  *tempxl = sumx;
  *tempyl = sumy;
  *tempzl = sumz;
}

/* ----------------------------------------------------------------------------------------------- */

// computes a 3D matrix-vector product along a 2D cut-plane

__device__  __forceinline__ void sum_hprimewgll_gamma(int I, int J, int K,
                                                 realw* tempxl,realw* tempyl,realw* tempzl,
                                                 realw* sh_tempx,realw* sh_tempy,realw* sh_tempz, realw* sh_hprimewgll ){

  realw fac;

  // initializes
  realw sumx = 0.f;
  realw sumy = 0.f;
  realw sumz = 0.f;

  // 3. cut-plane along gamma-direction
  #pragma unroll
  for (int l=0;l<NGLLX;l++) {
    fac = sh_hprimewgll[K*NGLLX+l]; // d_hprimewgll_xx[K*NGLLX+l];

    sumx += sh_tempx[l*NGLL2+J*NGLLX+I] * fac;
    sumy += sh_tempy[l*NGLL2+J*NGLLX+I] * fac;
    sumz += sh_tempz[l*NGLL2+J*NGLLX+I] * fac;
  }

  *tempxl = sumx;
  *tempyl = sumy;
  *tempzl = sumz;
}


/* ----------------------------------------------------------------------------------------------- */

// KERNEL 2
//
// for elastic domains

/* ----------------------------------------------------------------------------------------------- */

// note:
// kernel_2 is split into several kernels:
//  - a kernel without attenuation and for isotropic media: Kernel_2_noatt_iso_impl()
//  - a kernel without attenuation and for isotropic media with coloring: Kernel_2_noatt_iso_col_impl()
//  - a kernel without attenuation and for isotropic media with gravity: Kernel_2_noatt_iso_grav_impl()
//  - a kernel without attenuation and for anisotropic media: Kernel_2_noatt_ani_impl()
//  - a kernel including attenuation: Kernel_2_att_impl()
//
// this should help with performance:
// the high number of registers needed for our kernels limits the occupancy; separation tries to reduce this.


// kernel without attenuation
//
// we use templates to distinguish between calls with forward or adjoint texture fields

template<int FORWARD_OR_ADJOINT> __global__ void
#ifdef USE_LAUNCH_BOUNDS
// adds compiler specification
__launch_bounds__(NGLL3_PADDED,LAUNCH_MIN_BLOCKS)
#endif
// main kernel
Kernel_2_noatt_iso_impl(const int nb_blocks_to_compute,
                        const int* d_ibool,
                        const int* d_phase_ispec_inner_elastic,const int num_phase_ispec_elastic,
                        const int d_iphase,
                        realw_const_p d_displ,
                        realw_p d_accel,
                        realw* d_xix,realw* d_xiy,realw* d_xiz,
                        realw* d_etax,realw* d_etay,realw* d_etaz,
                        realw* d_gammax,realw* d_gammay,realw* d_gammaz,
                        realw_const_p d_hprime_xx,
                        realw_const_p d_hprimewgll_xx,
                        realw_const_p d_wgllwgll_xy,realw_const_p d_wgllwgll_xz,realw_const_p d_wgllwgll_yz,
                        realw* d_kappav,realw* d_muv){

// elastic compute kernel without attenuation for isotropic elements
//
// holds for:
//  ATTENUATION               = .false.
//  ANISOTROPY                = .false.
//  COMPUTE_AND_STORE_STRAIN  = .true. or .false. (true for kernel simulations)
//  gravity                   = .false.
//  use_mesh_coloring_gpu     = .false.
//  COMPUTE_AND_STORE_STRAIN  = .false.

  // block-id == number of local element id in phase_ispec array
  int bx = blockIdx.y*gridDim.x+blockIdx.x;

  // thread-id == GLL node id
  // note: use only NGLL^3 = 125 active threads, plus 3 inactive/ghost threads,
  //       because we used memory padding from NGLL^3 = 125 to 128 to get coalescent memory accesses;
  //       to avoid execution branching and the need of registers to store an active state variable,
  //       the thread ids are put in valid range
  int tx = threadIdx.x;

  int I,J,K;
  int iglob,offset;
  int working_element;

  realw tempx1l,tempx2l,tempx3l,tempy1l,tempy2l,tempy3l,tempz1l,tempz2l,tempz3l;
  realw xixl,xiyl,xizl,etaxl,etayl,etazl,gammaxl,gammayl,gammazl,jacobianl;
  realw duxdxl,duxdyl,duxdzl,duydxl,duydyl,duydzl,duzdxl,duzdyl,duzdzl;
  realw duxdxl_plus_duydyl,duxdxl_plus_duzdzl,duydyl_plus_duzdzl;
  realw duxdyl_plus_duydxl,duzdxl_plus_duxdzl,duzdyl_plus_duydzl;

  realw fac1,fac2,fac3;
  realw lambdal,mul,lambdalplus2mul,kappal;
  realw sigma_xx,sigma_yy,sigma_zz,sigma_xy,sigma_xz,sigma_yz;
  realw sum_terms1,sum_terms2,sum_terms3;

  // shared memory
  __shared__ realw sh_tempx[NGLL3];
  __shared__ realw sh_tempy[NGLL3];
  __shared__ realw sh_tempz[NGLL3];

  // note: using shared memory for hprime's improves performance
  //       (but could tradeoff with occupancy)
  __shared__ realw sh_hprime_xx[NGLL2];
  __shared__ realw sh_hprimewgll_xx[NGLL2];

// arithmetic intensity: ratio of number-of-arithmetic-operations / number-of-bytes-accessed-on-DRAM
//
// hand-counts on floating-point operations: counts addition/subtraction/multiplication/division
//                                           no counts for operations on indices in for-loops (compiler will likely unrool loops)
//
//                                           counts accesses to global memory, but no shared memory or register loads/stores
//                                           float has 4 bytes

// counts:
// 2 FLOP

  // checks if anything to do
  if (bx >= nb_blocks_to_compute) return;

  // limits thread ids to range [0,125-1]
  if (tx >= NGLL3) tx = NGLL3 - 1;

// counts:
// + 1 FLOP
//
// + 0 BYTE

  // loads hprime's into shared memory
  if (tx < NGLL2) {
    // copy hprime from global memory to shared memory
    load_shared_memory_hprime(&tx,d_hprime_xx,sh_hprime_xx);

    // copy hprimewgll from global memory to shared memory
    load_shared_memory_hprimewgll(&tx,d_hprimewgll_xx,sh_hprimewgll_xx);
  }

// counts:
// + 0 FLOP
//
// 2 * 1 float * 25 threads = 200 BYTE

  // spectral-element id
  // iphase-1 and working_element-1 for Fortran->C array conventions
  working_element = d_phase_ispec_inner_elastic[bx + num_phase_ispec_elastic*(d_iphase-1)] - 1;

  // local padded index
  offset = working_element*NGLL3_PADDED + tx;

  // global index
  iglob = d_ibool[offset] - 1 ;

// counts:
// + 7 FLOP
//
// + 2 float * 128 threads = 1024 BYTE

  // copy from global memory to shared memory
  // each thread writes one of the NGLL^3 = 125 data points
  if (threadIdx.x < NGLL3 ){
    // copy displacement from global memory to shared memory
    load_shared_memory_displ<FORWARD_OR_ADJOINT>(&tx,&iglob,d_displ,sh_tempx,sh_tempy,sh_tempz);
  }

// counts:
// + 5 FLOP
//
// + 3 float * 125 threads = 1500 BYTE

  kappal = d_kappav[offset];
  mul = d_muv[offset];

// counts:
// + 0 FLOP
//
// + 2 * 1 float * 128 threads = 1024 BYTE

  // loads mesh values here to give compiler possibility to overlap memory fetches with some computations
  // note: arguments defined as realw* instead of const realw* __restrict__ to avoid that the compiler
  //       loads all memory by texture loads
  //       we only use the first loads explicitly by texture loads, all subsequent without. this should lead/trick
  //       the compiler to use global memory loads for all the subsequent accesses.
  //
  // calculates laplacian
  xixl = get_global_cr( &d_xix[offset] ); // first array with texture load
  xiyl = get_global_cr( &d_xiy[offset] ); // first array with texture load
  xizl = get_global_cr( &d_xiz[offset] ); // first array with texture load

//  xixl = d_xix[offset]; // first array with texture load
//  xiyl = d_xiy[offset]; // all subsequent without to avoid over-use of texture for coalescent access
//  xizl = d_xiz[offset];

  etaxl = d_etax[offset];
  etayl = d_etay[offset];
  etazl = d_etaz[offset];

  gammaxl = d_gammax[offset];
  gammayl = d_gammay[offset];
  gammazl = d_gammaz[offset];

  jacobianl = 1.f / (xixl*(etayl*gammazl-etazl*gammayl)
                    -xiyl*(etaxl*gammazl-etazl*gammaxl)
                    +xizl*(etaxl*gammayl-etayl*gammaxl));

// counts:
// + 15 FLOP
//
// + 9 float * 128 threads = 4608 BYTE



  // local index
  K = (tx/NGLL2);
  J = ((tx-K*NGLL2)/NGLLX);
  I = (tx-K*NGLL2-J*NGLLX);

// counts:
// + 8 FLOP
//
// + 0 BYTE

  // synchronize all the threads (one thread for each of the NGLL grid points of the
  // current spectral element) because we need the whole element to be ready in order
  // to be able to compute the matrix products along cut planes of the 3D element below
  __syncthreads();

  // computes first matrix products
  // 1. cut-plane
  sum_hprime_xi(I,J,K,&tempx1l,&tempy1l,&tempz1l,sh_tempx,sh_tempy,sh_tempz,sh_hprime_xx);
  // 2. cut-plane
  sum_hprime_eta(I,J,K,&tempx2l,&tempy2l,&tempz2l,sh_tempx,sh_tempy,sh_tempz,sh_hprime_xx);
  // 3. cut-plane
  sum_hprime_gamma(I,J,K,&tempx3l,&tempy3l,&tempz3l,sh_tempx,sh_tempy,sh_tempz,sh_hprime_xx);

// counts:
// + 3 * 100 FLOP = 300 FLOP
//
// + 0 BYTE

  // compute derivatives of ux, uy and uz with respect to x, y and z
  duxdxl = xixl*tempx1l + etaxl*tempx2l + gammaxl*tempx3l;
  duxdyl = xiyl*tempx1l + etayl*tempx2l + gammayl*tempx3l;
  duxdzl = xizl*tempx1l + etazl*tempx2l + gammazl*tempx3l;

  duydxl = xixl*tempy1l + etaxl*tempy2l + gammaxl*tempy3l;
  duydyl = xiyl*tempy1l + etayl*tempy2l + gammayl*tempy3l;
  duydzl = xizl*tempy1l + etazl*tempy2l + gammazl*tempy3l;

  duzdxl = xixl*tempz1l + etaxl*tempz2l + gammaxl*tempz3l;
  duzdyl = xiyl*tempz1l + etayl*tempz2l + gammayl*tempz3l;
  duzdzl = xizl*tempz1l + etazl*tempz2l + gammazl*tempz3l;

// counts:
// + 9 * 5 FLOP = 45 FLOP
//
// + 0 BYTE

  // precompute some sums to save CPU time
  duxdxl_plus_duydyl = duxdxl + duydyl;
  duxdxl_plus_duzdzl = duxdxl + duzdzl;
  duydyl_plus_duzdzl = duydyl + duzdzl;
  duxdyl_plus_duydxl = duxdyl + duydxl;
  duzdxl_plus_duxdzl = duzdxl + duxdzl;
  duzdyl_plus_duydzl = duzdyl + duydzl;

  // stress calculations

  // isotropic case
  // compute elements with an elastic isotropic rheology

  lambdalplus2mul = kappal + 1.33333333333333333333f * mul;  // 4./3. = 1.3333333
  lambdal = lambdalplus2mul - 2.0f * mul;

  // compute the six components of the stress tensor sigma
  sigma_xx = lambdalplus2mul*duxdxl + lambdal*duydyl_plus_duzdzl;
  sigma_yy = lambdalplus2mul*duydyl + lambdal*duxdxl_plus_duzdzl;
  sigma_zz = lambdalplus2mul*duzdzl + lambdal*duxdxl_plus_duydyl;

  sigma_xy = mul*duxdyl_plus_duydxl;
  sigma_xz = mul*duzdxl_plus_duxdzl;
  sigma_yz = mul*duzdyl_plus_duydzl;

// counts:
// + 22 FLOP
//
// + 0 BYTE

  // form dot product with test vector, symmetric form
  // 1. cut-plane xi
  __syncthreads();
  // fills shared memory arrays
  if (threadIdx.x < NGLL3) {
    sh_tempx[tx] = jacobianl * (sigma_xx*xixl + sigma_xy*xiyl + sigma_xz*xizl); // sh_tempx1
    sh_tempy[tx] = jacobianl * (sigma_xy*xixl + sigma_yy*xiyl + sigma_yz*xizl); // sh_tempy1
    sh_tempz[tx] = jacobianl * (sigma_xz*xixl + sigma_yz*xiyl + sigma_zz*xizl); // sh_tempz1
  }
  __syncthreads();
  // 1. cut-plane xi
  sum_hprimewgll_xi(I,J,K,&tempx1l,&tempy1l,&tempz1l,sh_tempx,sh_tempy,sh_tempz,sh_hprimewgll_xx);

  // 2. cut-plane eta
  __syncthreads();
  // fills shared memory arrays
  if (threadIdx.x < NGLL3) {
    sh_tempx[tx] = jacobianl * (sigma_xx*etaxl + sigma_xy*etayl + sigma_xz*etazl); // sh_tempx2
    sh_tempy[tx] = jacobianl * (sigma_xy*etaxl + sigma_yy*etayl + sigma_yz*etazl); // sh_tempy2
    sh_tempz[tx] = jacobianl * (sigma_xz*etaxl + sigma_yz*etayl + sigma_zz*etazl); // sh_tempz2
  }
  __syncthreads();
  // 2. cut-plane eta
  sum_hprimewgll_eta(I,J,K,&tempx2l,&tempy2l,&tempz2l,sh_tempx,sh_tempy,sh_tempz,sh_hprimewgll_xx);

  // 3. cut-plane gamma
  __syncthreads();
  // fills shared memory arrays
  if (threadIdx.x < NGLL3) {
    sh_tempx[tx] = jacobianl * (sigma_xx*gammaxl + sigma_xy*gammayl + sigma_xz*gammazl); // sh_tempx3
    sh_tempy[tx] = jacobianl * (sigma_xy*gammaxl + sigma_yy*gammayl + sigma_yz*gammazl); // sh_tempy3
    sh_tempz[tx] = jacobianl * (sigma_xz*gammaxl + sigma_yz*gammayl + sigma_zz*gammazl); // sh_tempz3
  }
  __syncthreads();
  // 3. cut-plane gamma
  sum_hprimewgll_gamma(I,J,K,&tempx3l,&tempy3l,&tempz3l,sh_tempx,sh_tempy,sh_tempz,sh_hprimewgll_xx);

// counts:
// + 3 * 3 * 6 FLOP = 54 FLOP
// + 3 * 100 FLOP = 300 FLOP
//
// + 0 BYTE

  // gets double weights
  fac1 = d_wgllwgll_yz[K*NGLLX+J];
  fac2 = d_wgllwgll_xz[K*NGLLX+I];
  fac3 = d_wgllwgll_xy[J*NGLLX+I];

// counts:
// + 3 * 2 FLOP = 6 FLOP
//
// + 3 float * 128 threads = 1536 BYTE

  sum_terms1 = - (fac1*tempx1l + fac2*tempx2l + fac3*tempx3l);
  sum_terms2 = - (fac1*tempy1l + fac2*tempy2l + fac3*tempy3l);
  sum_terms3 = - (fac1*tempz1l + fac2*tempz2l + fac3*tempz3l);

// counts:
// + 3 * 6 FLOP = 18 FLOP
//
// + 0 BYTE

  // assembles acceleration array
  if (threadIdx.x < NGLL3) {
    atomicAdd(&d_accel[iglob*3], sum_terms1);
    atomicAdd(&d_accel[iglob*3+1], sum_terms2);
    atomicAdd(&d_accel[iglob*3+2], sum_terms3);
  }

// counts:
// + 8 FLOP
//
// + 3 float * 125 threads = 1500 BYTE


// counts:
// -----------------
// total of: 790 FLOP per thread
//           ~ 128 * 790 = 101120 FLOP per block
//
//           11392 BYTE DRAM accesses per block
//
// arithmetic intensity: 101120 FLOP / 11392 BYTES ~ 8.9 FLOP/BYTE
// -----------------
//
// nvprof: nvprof --metrics flops_sp ./xspecfem3D
//          -> 883146240 FLOPS (Single) floating-point operations for 20736 elements
//          -> 42590 FLOP per block
// arithmetic intensity: 42590 FLOP / 11392 BYTES ~ 3.74 FLOP/BYTE
//
// roofline model: Kepler K20x
// ---------------------------
//   for a Kepler K20x card, the peak single-precision performance is about 3.95 TFlop/s.
//   global memory access has a bandwidth of ~ 250 GB/s.
//
//   memory bandwidth: 250 GB/s
//   single-precision peak performance: 3.95 TFlop/s -> corner arithmetic intensity = 3950./250. ~ 15.8 flop/byte
//
//   elastic kernel has an arithmetic intensity of: hand-counts   ~ 8.9 flop/byte
//                                                  nvprof-counts ~ 42590./11392. flop/byte = 3.74 flop/byte
//
//   -> we can only achieve about: (hand-counts)   56% of the peak performance
//                                 (nvprof-counts) 24% of the peak performance -> 935.0 GFlop/s
//
// roofline model: Tesla K20c (Kepler architecture: http://www.nvidia.com/content/tesla/pdf/Tesla-KSeries-Overview-LR.pdf)
// ---------------------------
//   memory bandwidth: 208 GB/s
//   single-precision peak performance: 3.52 TFlop/s -> corner arithmetic intensity = 3520 / 208 ~ 16.9 flop/byte
//
//   we can only achieve about: (hand-counts)   52% of the peak performance
//                              (nvprof-counts) 22% of the peak performance -> 779.0 GFlop/s - measured: 647.3 GFlop/s


} // kernel_2_noatt_iso_impl()

/* ----------------------------------------------------------------------------------------------- */


template<int FORWARD_OR_ADJOINT> __global__ void
#ifdef USE_LAUNCH_BOUNDS
// adds compiler specification
__launch_bounds__(NGLL3_PADDED,LAUNCH_MIN_BLOCKS)
#endif
// main kernel
Kernel_2_noatt_iso_strain_impl(int nb_blocks_to_compute,
                              const int* d_ibool,
                              const int* d_phase_ispec_inner_elastic,const int num_phase_ispec_elastic,
                              const int d_iphase,
                              realw_const_p d_displ,
                              realw_p d_accel,
                              realw* d_xix,realw* d_xiy,realw* d_xiz,
                              realw* d_etax,realw* d_etay,realw* d_etaz,
                              realw* d_gammax,realw* d_gammay,realw* d_gammaz,
                              realw_const_p d_hprime_xx,
                              realw_const_p d_hprimewgll_xx,
                              realw_const_p d_wgllwgll_xy,realw_const_p d_wgllwgll_xz,realw_const_p d_wgllwgll_yz,
                              realw_const_p d_kappav,realw_const_p d_muv,
                              const int COMPUTE_AND_STORE_STRAIN,
                              realw_p epsilondev_xx,realw_p epsilondev_yy,realw_p epsilondev_xy,
                              realw_p epsilondev_xz,realw_p epsilondev_yz,
                              realw_p epsilon_trace_over_3,
                              const int SIMULATION_TYPE){

// elastic compute kernel without attenuation for isotropic elements
//
// holds for:
//  ATTENUATION               = .false.
//  ANISOTROPY                = .false.
//  COMPUTE_AND_STORE_STRAIN  = .true. or .false. (true for kernel simulations)
//  gravity                   = .false.
//  use_mesh_coloring_gpu     = .false.
//  COMPUTE_AND_STORE_STRAIN  = .true.

  // block-id == number of local element id in phase_ispec array
  int bx = blockIdx.y*gridDim.x+blockIdx.x;

  // checks if anything to do
  if (bx >= nb_blocks_to_compute) return;

  // thread-id == GLL node id
  // note: use only NGLL^3 = 125 active threads, plus 3 inactive/ghost threads,
  //       because we used memory padding from NGLL^3 = 125 to 128 to get coalescent memory accesses;
  //       to avoid execution branching and the need of registers to store an active state variable,
  //       the thread ids are put in valid range
  int tx = threadIdx.x;
  if (tx >= NGLL3) tx = NGLL3 - 1;

  int K = (tx/NGLL2);
  int J = ((tx-K*NGLL2)/NGLLX);
  int I = (tx-K*NGLL2-J*NGLLX);

  int iglob,offset;
  int working_element;

  realw tempx1l,tempx2l,tempx3l,tempy1l,tempy2l,tempy3l,tempz1l,tempz2l,tempz3l;
  realw xixl,xiyl,xizl,etaxl,etayl,etazl,gammaxl,gammayl,gammazl,jacobianl;
  realw duxdxl,duxdyl,duxdzl,duydxl,duydyl,duydzl,duzdxl,duzdyl,duzdzl;
  realw duxdxl_plus_duydyl,duxdxl_plus_duzdzl,duydyl_plus_duzdzl;
  realw duxdyl_plus_duydxl,duzdxl_plus_duxdzl,duzdyl_plus_duydzl;

  realw fac1,fac2,fac3;
  realw lambdal,mul,lambdalplus2mul,kappal;
  realw sigma_xx,sigma_yy,sigma_zz,sigma_xy,sigma_xz,sigma_yz;
  realw sum_terms1,sum_terms2,sum_terms3;

  // shared memory
  __shared__ realw sh_tempx[NGLL3];
  __shared__ realw sh_tempy[NGLL3];
  __shared__ realw sh_tempz[NGLL3];

  // note: using shared memory for hprime's improves performance
  //       (but could tradeoff with occupancy)
  __shared__ realw sh_hprime_xx[NGLL2];
  __shared__ realw sh_hprimewgll_xx[NGLL2];

  // loads hprime's into shared memory
  if (tx < NGLL2) {
    // copy hprime from global memory to shared memory
    load_shared_memory_hprime(&tx,d_hprime_xx,sh_hprime_xx);
    // copy hprime from global memory to shared memory
    load_shared_memory_hprimewgll(&tx,d_hprimewgll_xx,sh_hprimewgll_xx);
  }

  // spectral-element id
  // iphase-1 and working_element-1 for Fortran->C array conventions
  working_element = d_phase_ispec_inner_elastic[bx + num_phase_ispec_elastic*(d_iphase-1)] - 1;

  // local padded index
  offset = working_element*NGLL3_PADDED + tx;

  // global index
  iglob = d_ibool[offset] - 1 ;

  // copy from global memory to shared memory
  // each thread writes one of the NGLL^3 = 125 data points
  if (threadIdx.x < NGLL3 ){
    // copy displacement from global memory to shared memory
    load_shared_memory_displ<FORWARD_OR_ADJOINT>(&tx,&iglob,d_displ,sh_tempx,sh_tempy,sh_tempz);
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

  // synchronize all the threads (one thread for each of the NGLL grid points of the
  // current spectral element) because we need the whole element to be ready in order
  // to be able to compute the matrix products along cut planes of the 3D element below
  __syncthreads();

  // computes first matrix products
  // 1. cut-plane
  sum_hprime_xi(I,J,K,&tempx1l,&tempy1l,&tempz1l,sh_tempx,sh_tempy,sh_tempz,sh_hprime_xx);
  // 2. cut-plane
  sum_hprime_eta(I,J,K,&tempx2l,&tempy2l,&tempz2l,sh_tempx,sh_tempy,sh_tempz,sh_hprime_xx);
  // 3. cut-plane
  sum_hprime_gamma(I,J,K,&tempx3l,&tempy3l,&tempz3l,sh_tempx,sh_tempy,sh_tempz,sh_hprime_xx);

  // compute derivatives of ux, uy and uz with respect to x, y and z
  duxdxl = xixl*tempx1l + etaxl*tempx2l + gammaxl*tempx3l;
  duxdyl = xiyl*tempx1l + etayl*tempx2l + gammayl*tempx3l;
  duxdzl = xizl*tempx1l + etazl*tempx2l + gammazl*tempx3l;

  duydxl = xixl*tempy1l + etaxl*tempy2l + gammaxl*tempy3l;
  duydyl = xiyl*tempy1l + etayl*tempy2l + gammayl*tempy3l;
  duydzl = xizl*tempy1l + etazl*tempy2l + gammazl*tempy3l;

  duzdxl = xixl*tempz1l + etaxl*tempz2l + gammaxl*tempz3l;
  duzdyl = xiyl*tempz1l + etayl*tempz2l + gammayl*tempz3l;
  duzdzl = xizl*tempz1l + etazl*tempz2l + gammazl*tempz3l;

  // precompute some sums to save CPU time
  duxdxl_plus_duydyl = duxdxl + duydyl;
  duxdxl_plus_duzdzl = duxdxl + duzdzl;
  duydyl_plus_duzdzl = duydyl + duzdzl;
  duxdyl_plus_duydxl = duxdyl + duydxl;
  duzdxl_plus_duxdzl = duzdxl + duxdzl;
  duzdyl_plus_duydzl = duzdyl + duydzl;

  // computes deviatoric strain for kernel calculations
  if (COMPUTE_AND_STORE_STRAIN) {
    // save deviatoric strain for Runge-Kutta scheme
    if (threadIdx.x < NGLL3) {
      realw templ = 0.33333333333333333333f * (duxdxl + duydyl + duzdzl); // 1./3. = 0.33333
      // local storage: stresses at this current time step
      // fortran: epsilondev_xx(:,:,:,ispec) = epsilondev_xx_loc(:,:,:)
      epsilondev_xx[tx + working_element*NGLL3] = duxdxl - templ; // epsilondev_xx_loc;
      epsilondev_yy[tx + working_element*NGLL3] = duydyl - templ; // epsilondev_yy_loc;
      epsilondev_xy[tx + working_element*NGLL3] = 0.5f * duxdyl_plus_duydxl; // epsilondev_xy_loc;
      epsilondev_xz[tx + working_element*NGLL3] = 0.5f * duzdxl_plus_duxdzl; // epsilondev_xz_loc;
      epsilondev_yz[tx + working_element*NGLL3] = 0.5f * duzdyl_plus_duydzl; //epsilondev_yz_loc;
      // kernel simulations
      if (SIMULATION_TYPE == 3){
        epsilon_trace_over_3[tx + working_element*NGLL3] = templ;
      }
    } // threadIdx.x
  }

  // stress calculations

  // isotropic case
  // compute elements with an elastic isotropic rheology
  kappal = d_kappav[offset];
  mul = d_muv[offset];

  lambdalplus2mul = kappal + 1.33333333333333333333f * mul;  // 4./3. = 1.3333333
  lambdal = lambdalplus2mul - 2.0f * mul;

  // compute the six components of the stress tensor sigma
  sigma_xx = lambdalplus2mul*duxdxl + lambdal*duydyl_plus_duzdzl;
  sigma_yy = lambdalplus2mul*duydyl + lambdal*duxdxl_plus_duzdzl;
  sigma_zz = lambdalplus2mul*duzdzl + lambdal*duxdxl_plus_duydyl;

  sigma_xy = mul*duxdyl_plus_duydxl;
  sigma_xz = mul*duzdxl_plus_duxdzl;
  sigma_yz = mul*duzdyl_plus_duydzl;

  // form dot product with test vector, symmetric form
  // 1. cut-plane xi
  __syncthreads();
  // fills shared memory arrays
  if (threadIdx.x < NGLL3) {
    sh_tempx[tx] = jacobianl * (sigma_xx*xixl + sigma_xy*xiyl + sigma_xz*xizl); // sh_tempx1
    sh_tempy[tx] = jacobianl * (sigma_xy*xixl + sigma_yy*xiyl + sigma_yz*xizl); // sh_tempy1
    sh_tempz[tx] = jacobianl * (sigma_xz*xixl + sigma_yz*xiyl + sigma_zz*xizl); // sh_tempz1
  }
  __syncthreads();
  // 1. cut-plane xi
  sum_hprimewgll_xi(I,J,K,&tempx1l,&tempy1l,&tempz1l,sh_tempx,sh_tempy,sh_tempz,sh_hprimewgll_xx);

  // 2. cut-plane eta
  __syncthreads();
  // fills shared memory arrays
  if (threadIdx.x < NGLL3) {
    sh_tempx[tx] = jacobianl * (sigma_xx*etaxl + sigma_xy*etayl + sigma_xz*etazl); // sh_tempx2
    sh_tempy[tx] = jacobianl * (sigma_xy*etaxl + sigma_yy*etayl + sigma_yz*etazl); // sh_tempy2
    sh_tempz[tx] = jacobianl * (sigma_xz*etaxl + sigma_yz*etayl + sigma_zz*etazl); // sh_tempz2
  }
  __syncthreads();
  // 2. cut-plane eta
  sum_hprimewgll_eta(I,J,K,&tempx2l,&tempy2l,&tempz2l,sh_tempx,sh_tempy,sh_tempz,sh_hprimewgll_xx);

  // 3. cut-plane gamma
  __syncthreads();
  // fills shared memory arrays
  if (threadIdx.x < NGLL3) {
    sh_tempx[tx] = jacobianl * (sigma_xx*gammaxl + sigma_xy*gammayl + sigma_xz*gammazl); // sh_tempx3
    sh_tempy[tx] = jacobianl * (sigma_xy*gammaxl + sigma_yy*gammayl + sigma_yz*gammazl); // sh_tempy3
    sh_tempz[tx] = jacobianl * (sigma_xz*gammaxl + sigma_yz*gammayl + sigma_zz*gammazl); // sh_tempz3
  }
  __syncthreads();
  // 3. cut-plane gamma
  sum_hprimewgll_gamma(I,J,K,&tempx3l,&tempy3l,&tempz3l,sh_tempx,sh_tempy,sh_tempz,sh_hprimewgll_xx);

  // gets double weights
  fac1 = d_wgllwgll_yz[K*NGLLX+J];
  fac2 = d_wgllwgll_xz[K*NGLLX+I];
  fac3 = d_wgllwgll_xy[J*NGLLX+I];

  sum_terms1 = - (fac1*tempx1l + fac2*tempx2l + fac3*tempx3l);
  sum_terms2 = - (fac1*tempy1l + fac2*tempy2l + fac3*tempy3l);
  sum_terms3 = - (fac1*tempz1l + fac2*tempz2l + fac3*tempz3l);

  // assembles acceleration array
  if (threadIdx.x < NGLL3) {
    atomicAdd(&d_accel[iglob*3], sum_terms1);
    atomicAdd(&d_accel[iglob*3+1], sum_terms2);
    atomicAdd(&d_accel[iglob*3+2], sum_terms3);
  } // threadIdx.x

} // kernel_2_noatt_iso_strain_impl()

/* ----------------------------------------------------------------------------------------------- */

template<int FORWARD_OR_ADJOINT> __global__ void
#ifdef USE_LAUNCH_BOUNDS
// adds compiler specification
__launch_bounds__(NGLL3_PADDED,LAUNCH_MIN_BLOCKS)
#endif
// main kernel
Kernel_2_noatt_iso_col_impl(int nb_blocks_to_compute,
                        const int* d_ibool,
                        const int* d_phase_ispec_inner_elastic,const int num_phase_ispec_elastic,
                        const int d_iphase,
                        const int use_mesh_coloring_gpu,
                        realw_const_p d_displ,
                        realw_p d_accel,
                        realw* d_xix,realw* d_xiy,realw* d_xiz,
                        realw* d_etax,realw* d_etay,realw* d_etaz,
                        realw* d_gammax,realw* d_gammay,realw* d_gammaz,
                        realw_const_p d_hprime_xx,
                        realw_const_p d_hprimewgll_xx,
                        realw_const_p d_wgllwgll_xy,realw_const_p d_wgllwgll_xz,realw_const_p d_wgllwgll_yz,
                        realw_const_p d_kappav,realw_const_p d_muv,
                        const int COMPUTE_AND_STORE_STRAIN,
                        realw_p epsilondev_xx,realw_p epsilondev_yy,realw_p epsilondev_xy,
                        realw_p epsilondev_xz,realw_p epsilondev_yz,
                        realw_p epsilon_trace_over_3,
                        const int SIMULATION_TYPE){

// elastic compute kernel without attenuation for isotropic elements
//
// holds for:
//  ATTENUATION               = .false.
//  ANISOTROPY                = .false.
//  COMPUTE_AND_STORE_STRAIN  = .true. or .false. (true for kernel simulations)
//  gravity                   = .false.
//  use_mesh_coloring_gpu     = .true.

  // block-id == number of local element id in phase_ispec array
  int bx = blockIdx.y*gridDim.x+blockIdx.x;

  // checks if anything to do
  if (bx >= nb_blocks_to_compute) return;

  // thread-id == GLL node id
  // note: use only NGLL^3 = 125 active threads, plus 3 inactive/ghost threads,
  //       because we used memory padding from NGLL^3 = 125 to 128 to get coalescent memory accesses;
  //       to avoid execution branching and the need of registers to store an active state variable,
  //       the thread ids are put in valid range
  int tx = threadIdx.x;
  if (tx >= NGLL3) tx = NGLL3-1;

  int K = (tx/NGLL2);
  int J = ((tx-K*NGLL2)/NGLLX);
  int I = (tx-K*NGLL2-J*NGLLX);

  int iglob,offset;
  int working_element;

  realw tempx1l,tempx2l,tempx3l,tempy1l,tempy2l,tempy3l,tempz1l,tempz2l,tempz3l;
  realw xixl,xiyl,xizl,etaxl,etayl,etazl,gammaxl,gammayl,gammazl,jacobianl;
  realw duxdxl,duxdyl,duxdzl,duydxl,duydyl,duydzl,duzdxl,duzdyl,duzdzl;
  realw duxdxl_plus_duydyl,duxdxl_plus_duzdzl,duydyl_plus_duzdzl;
  realw duxdyl_plus_duydxl,duzdxl_plus_duxdzl,duzdyl_plus_duydzl;

  realw fac1,fac2,fac3;
  realw lambdal,mul,lambdalplus2mul,kappal;
  realw sigma_xx,sigma_yy,sigma_zz,sigma_xy,sigma_xz,sigma_yz;
  realw sum_terms1,sum_terms2,sum_terms3;

  // shared memory
  __shared__ realw sh_tempx[NGLL3];
  __shared__ realw sh_tempy[NGLL3];
  __shared__ realw sh_tempz[NGLL3];

  // note: using shared memory for hprime's improves performance
  //       (but could tradeoff with occupancy)
  __shared__ realw sh_hprime_xx[NGLL2];
  __shared__ realw sh_hprimewgll_xx[NGLL2];

  // loads hprime's into shared memory
  if (tx < NGLL2) {
    // copy hprime from global memory to shared memory
    load_shared_memory_hprime(&tx,d_hprime_xx,sh_hprime_xx);
    // copy hprime from global memory to shared memory
    load_shared_memory_hprimewgll(&tx,d_hprimewgll_xx,sh_hprimewgll_xx);
  }

  // spectral-element id
  // iphase-1 and working_element-1 for Fortran->C array conventions
#ifdef USE_MESH_COLORING_GPU
  working_element = bx;
#else
  //mesh coloring
  if (use_mesh_coloring_gpu ){
    working_element = bx;
  }else{
    // iphase-1 and working_element-1 for Fortran->C array conventions
    working_element = d_phase_ispec_inner_elastic[bx + num_phase_ispec_elastic*(d_iphase-1)] - 1;
  }
#endif
  // local padded index
  offset = working_element*NGLL3_PADDED + tx;

  // global index
  iglob = d_ibool[offset] - 1 ;

  // copy from global memory to shared memory
  // each thread writes one of the NGLL^3 = 125 data points
  if (threadIdx.x < NGLL3 ){
    // copy displacement from global memory to shared memory
    load_shared_memory_displ<FORWARD_OR_ADJOINT>(&tx,&iglob,d_displ,sh_tempx,sh_tempy,sh_tempz);
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

  // synchronize all the threads (one thread for each of the NGLL grid points of the
  // current spectral element) because we need the whole element to be ready in order
  // to be able to compute the matrix products along cut planes of the 3D element below
  __syncthreads();

  // computes first matrix products
  // 1. cut-plane
  sum_hprime_xi(I,J,K,&tempx1l,&tempy1l,&tempz1l,sh_tempx,sh_tempy,sh_tempz,sh_hprime_xx);
  // 2. cut-plane
  sum_hprime_eta(I,J,K,&tempx2l,&tempy2l,&tempz2l,sh_tempx,sh_tempy,sh_tempz,sh_hprime_xx);
  // 3. cut-plane
  sum_hprime_gamma(I,J,K,&tempx3l,&tempy3l,&tempz3l,sh_tempx,sh_tempy,sh_tempz,sh_hprime_xx);

  // compute derivatives of ux, uy and uz with respect to x, y and z
  duxdxl = xixl*tempx1l + etaxl*tempx2l + gammaxl*tempx3l;
  duxdyl = xiyl*tempx1l + etayl*tempx2l + gammayl*tempx3l;
  duxdzl = xizl*tempx1l + etazl*tempx2l + gammazl*tempx3l;

  duydxl = xixl*tempy1l + etaxl*tempy2l + gammaxl*tempy3l;
  duydyl = xiyl*tempy1l + etayl*tempy2l + gammayl*tempy3l;
  duydzl = xizl*tempy1l + etazl*tempy2l + gammazl*tempy3l;

  duzdxl = xixl*tempz1l + etaxl*tempz2l + gammaxl*tempz3l;
  duzdyl = xiyl*tempz1l + etayl*tempz2l + gammayl*tempz3l;
  duzdzl = xizl*tempz1l + etazl*tempz2l + gammazl*tempz3l;

  // precompute some sums to save CPU time
  duxdxl_plus_duydyl = duxdxl + duydyl;
  duxdxl_plus_duzdzl = duxdxl + duzdzl;
  duydyl_plus_duzdzl = duydyl + duzdzl;
  duxdyl_plus_duydxl = duxdyl + duydxl;
  duzdxl_plus_duxdzl = duzdxl + duxdzl;
  duzdyl_plus_duydzl = duzdyl + duydzl;

  // computes deviatoric strain for kernel calculations
  if (COMPUTE_AND_STORE_STRAIN) {
    // save deviatoric strain for Runge-Kutta scheme
    if (threadIdx.x < NGLL3) {
      realw templ = 0.33333333333333333333f * (duxdxl + duydyl + duzdzl); // 1./3. = 0.33333
      // local storage: stresses at this current time step
      // fortran: epsilondev_xx(:,:,:,ispec) = epsilondev_xx_loc(:,:,:)
      epsilondev_xx[tx + working_element*NGLL3] = duxdxl - templ; // epsilondev_xx_loc;
      epsilondev_yy[tx + working_element*NGLL3] = duydyl - templ; // epsilondev_yy_loc;
      epsilondev_xy[tx + working_element*NGLL3] = 0.5f * duxdyl_plus_duydxl; // epsilondev_xy_loc;
      epsilondev_xz[tx + working_element*NGLL3] = 0.5f * duzdxl_plus_duxdzl; // epsilondev_xz_loc;
      epsilondev_yz[tx + working_element*NGLL3] = 0.5f * duzdyl_plus_duydzl; //epsilondev_yz_loc;
      // kernel simulations
      if (SIMULATION_TYPE == 3){
        epsilon_trace_over_3[tx + working_element*NGLL3] = templ;
      }
    } // threadIdx.x
  }

  // stress calculations

  // isotropic case
  // compute elements with an elastic isotropic rheology
  kappal = d_kappav[offset];
  mul = d_muv[offset];

  lambdalplus2mul = kappal + 1.33333333333333333333f * mul;  // 4./3. = 1.3333333
  lambdal = lambdalplus2mul - 2.0f * mul;

  // compute the six components of the stress tensor sigma
  sigma_xx = lambdalplus2mul*duxdxl + lambdal*duydyl_plus_duzdzl;
  sigma_yy = lambdalplus2mul*duydyl + lambdal*duxdxl_plus_duzdzl;
  sigma_zz = lambdalplus2mul*duzdzl + lambdal*duxdxl_plus_duydyl;

  sigma_xy = mul*duxdyl_plus_duydxl;
  sigma_xz = mul*duzdxl_plus_duxdzl;
  sigma_yz = mul*duzdyl_plus_duydzl;

  // form dot product with test vector, symmetric form
  // 1. cut-plane xi
  __syncthreads();
  // fills shared memory arrays
  if (threadIdx.x < NGLL3) {
    sh_tempx[tx] = jacobianl * (sigma_xx*xixl + sigma_xy*xiyl + sigma_xz*xizl); // sh_tempx1
    sh_tempy[tx] = jacobianl * (sigma_xy*xixl + sigma_yy*xiyl + sigma_yz*xizl); // sh_tempy1
    sh_tempz[tx] = jacobianl * (sigma_xz*xixl + sigma_yz*xiyl + sigma_zz*xizl); // sh_tempz1
  }
  __syncthreads();
  // 1. cut-plane xi
  sum_hprimewgll_xi(I,J,K,&tempx1l,&tempy1l,&tempz1l,sh_tempx,sh_tempy,sh_tempz,sh_hprimewgll_xx);

  // 2. cut-plane eta
  __syncthreads();
  // fills shared memory arrays
  if (threadIdx.x < NGLL3) {
    sh_tempx[tx] = jacobianl * (sigma_xx*etaxl + sigma_xy*etayl + sigma_xz*etazl); // sh_tempx2
    sh_tempy[tx] = jacobianl * (sigma_xy*etaxl + sigma_yy*etayl + sigma_yz*etazl); // sh_tempy2
    sh_tempz[tx] = jacobianl * (sigma_xz*etaxl + sigma_yz*etayl + sigma_zz*etazl); // sh_tempz2
  }
  __syncthreads();
  // 2. cut-plane eta
  sum_hprimewgll_eta(I,J,K,&tempx2l,&tempy2l,&tempz2l,sh_tempx,sh_tempy,sh_tempz,sh_hprimewgll_xx);

  // 3. cut-plane gamma
  __syncthreads();
  // fills shared memory arrays
  if (threadIdx.x < NGLL3) {
    sh_tempx[tx] = jacobianl * (sigma_xx*gammaxl + sigma_xy*gammayl + sigma_xz*gammazl); // sh_tempx3
    sh_tempy[tx] = jacobianl * (sigma_xy*gammaxl + sigma_yy*gammayl + sigma_yz*gammazl); // sh_tempy3
    sh_tempz[tx] = jacobianl * (sigma_xz*gammaxl + sigma_yz*gammayl + sigma_zz*gammazl); // sh_tempz3
  }
  __syncthreads();
  // 3. cut-plane gamma
  sum_hprimewgll_gamma(I,J,K,&tempx3l,&tempy3l,&tempz3l,sh_tempx,sh_tempy,sh_tempz,sh_hprimewgll_xx);

  // gets double weights
  fac1 = d_wgllwgll_yz[K*NGLLX+J];
  fac2 = d_wgllwgll_xz[K*NGLLX+I];
  fac3 = d_wgllwgll_xy[J*NGLLX+I];

  sum_terms1 = - (fac1*tempx1l + fac2*tempx2l + fac3*tempx3l);
  sum_terms2 = - (fac1*tempy1l + fac2*tempy2l + fac3*tempy3l);
  sum_terms3 = - (fac1*tempz1l + fac2*tempz2l + fac3*tempz3l);

  // assembles acceleration array
  if (threadIdx.x < NGLL3) {

#ifdef USE_MESH_COLORING_GPU
    // no atomic operation needed, colors don't share global points between elements

#ifdef USE_TEXTURES_FIELDS
    d_accel[iglob*3]     = texfetch_accel<FORWARD_OR_ADJOINT>(iglob*3) + sum_terms1;
    d_accel[iglob*3 + 1] = texfetch_accel<FORWARD_OR_ADJOINT>(iglob*3 + 1) + sum_terms2;
    d_accel[iglob*3 + 2] = texfetch_accel<FORWARD_OR_ADJOINT>(iglob*3 + 2) + sum_terms3;
#else
    d_accel[iglob*3]     += sum_terms1;
    d_accel[iglob*3 + 1] += sum_terms2;
    d_accel[iglob*3 + 2] += sum_terms3;
#endif // USE_TEXTURES_FIELDS

#else // MESH_COLORING

    //mesh coloring
    if (use_mesh_coloring_gpu ){

      // no atomic operation needed, colors don't share global points between elements
#ifdef USE_TEXTURES_FIELDS
      d_accel[iglob*3]     = texfetch_accel<FORWARD_OR_ADJOINT>(iglob*3) + sum_terms1;
      d_accel[iglob*3 + 1] = texfetch_accel<FORWARD_OR_ADJOINT>(iglob*3 + 1) + sum_terms2;
      d_accel[iglob*3 + 2] = texfetch_accel<FORWARD_OR_ADJOINT>(iglob*3 + 2) + sum_terms3;
#else
      d_accel[iglob*3]     += sum_terms1;
      d_accel[iglob*3 + 1] += sum_terms2;
      d_accel[iglob*3 + 2] += sum_terms3;
#endif // USE_TEXTURES_FIELDS

    }else {
      atomicAdd(&d_accel[iglob*3], sum_terms1);
      atomicAdd(&d_accel[iglob*3+1], sum_terms2);
      atomicAdd(&d_accel[iglob*3+2], sum_terms3);
    } // if (use_mesh_coloring_gpu)

#endif // MESH_COLORING

  } // threadIdx.x

} // kernel_2_noatt_iso_col_impl()

/* ----------------------------------------------------------------------------------------------- */


template<int FORWARD_OR_ADJOINT> __global__ void
#ifdef USE_LAUNCH_BOUNDS
// adds compiler specification
__launch_bounds__(NGLL3_PADDED,LAUNCH_MIN_BLOCKS)
#endif
// main kernel
Kernel_2_noatt_iso_grav_impl(int nb_blocks_to_compute,
                        const int* d_ibool,
                        const int* d_phase_ispec_inner_elastic,const int num_phase_ispec_elastic,
                        const int d_iphase,
                        const int use_mesh_coloring_gpu,
                        realw_const_p d_displ,
                        realw_p d_accel,
                        realw* d_xix,realw* d_xiy,realw* d_xiz,
                        realw* d_etax,realw* d_etay,realw* d_etaz,
                        realw* d_gammax,realw* d_gammay,realw* d_gammaz,
                        realw_const_p d_hprime_xx,
                        realw_const_p d_hprimewgll_xx,
                        realw_const_p d_wgllwgll_xy,realw_const_p d_wgllwgll_xz,realw_const_p d_wgllwgll_yz,
                        realw_const_p d_kappav,realw_const_p d_muv,
                        const int COMPUTE_AND_STORE_STRAIN,
                        realw_p epsilondev_xx,realw_p epsilondev_yy,realw_p epsilondev_xy,
                        realw_p epsilondev_xz,realw_p epsilondev_yz,
                        realw_p epsilon_trace_over_3,
                        const int SIMULATION_TYPE,
                        const int gravity,
                        realw_const_p d_minus_g,
                        realw_const_p d_minus_deriv_gravity,
                        realw_const_p d_rhostore,
                        realw_const_p wgll_cube ){

// elastic compute kernel without attenuation for isotropic elements
//
// holds for:
//  ATTENUATION               = .false.
//  ANISOTROPY                = .false.
//  COMPUTE_AND_STORE_STRAIN  = .true. or .false. (true for kernel simulations)
//  gravity                   = .true.
//

  // block-id == number of local element id in phase_ispec array
  int bx = blockIdx.y*gridDim.x+blockIdx.x;

  // checks if anything to do
  if (bx >= nb_blocks_to_compute) return;

  // thread-id == GLL node id
  // note: use only NGLL^3 = 125 active threads, plus 3 inactive/ghost threads,
  //       because we used memory padding from NGLL^3 = 125 to 128 to get coalescent memory accesses;
  //       to avoid execution branching and the need of registers to store an active state variable,
  //       the thread ids are put in valid range
  int tx = threadIdx.x;
  if (tx >= NGLL3) tx = NGLL3-1;

  int K = (tx/NGLL2);
  int J = ((tx-K*NGLL2)/NGLLX);
  int I = (tx-K*NGLL2-J*NGLLX);

  int iglob,offset;
  int working_element;

  realw tempx1l,tempx2l,tempx3l,tempy1l,tempy2l,tempy3l,tempz1l,tempz2l,tempz3l;
  realw xixl,xiyl,xizl,etaxl,etayl,etazl,gammaxl,gammayl,gammazl,jacobianl;
  realw duxdxl,duxdyl,duxdzl,duydxl,duydyl,duydzl,duzdxl,duzdyl,duzdzl;
  realw duxdxl_plus_duydyl,duxdxl_plus_duzdzl,duydyl_plus_duzdzl;
  realw duxdyl_plus_duydxl,duzdxl_plus_duxdzl,duzdyl_plus_duydzl;

  realw fac1,fac2,fac3;
  realw lambdal,mul,lambdalplus2mul,kappal;
  realw sigma_xx,sigma_yy,sigma_zz,sigma_xy,sigma_xz,sigma_yz;
  realw sum_terms1,sum_terms2,sum_terms3;

  // gravity variables
  realw sigma_yx,sigma_zx,sigma_zy;
  realw rho_s_H1 = 0.f;
  realw rho_s_H2 = 0.f;
  realw rho_s_H3 = 0.f;

  // shared memory
  __shared__ realw sh_tempx[NGLL3];
  __shared__ realw sh_tempy[NGLL3];
  __shared__ realw sh_tempz[NGLL3];

  // note: using shared memory for hprime's improves performance
  //       (but could tradeoff with occupancy)
  __shared__ realw sh_hprime_xx[NGLL2];
  __shared__ realw sh_hprimewgll_xx[NGLL2];

  // loads hprime's into shared memory
  if (tx < NGLL2) {
    // copy hprime from global memory to shared memory
    load_shared_memory_hprime(&tx,d_hprime_xx,sh_hprime_xx);
    // copy hprime from global memory to shared memory
    load_shared_memory_hprimewgll(&tx,d_hprimewgll_xx,sh_hprimewgll_xx);
  }

  // spectral-element id
  // iphase-1 and working_element-1 for Fortran->C array conventions
#ifdef USE_MESH_COLORING_GPU
  working_element = bx;
#else
  //mesh coloring
  if (use_mesh_coloring_gpu ){
    working_element = bx;
  }else{
    // iphase-1 and working_element-1 for Fortran->C array conventions
    working_element = d_phase_ispec_inner_elastic[bx + num_phase_ispec_elastic*(d_iphase-1)] - 1;
  }
#endif
  // local padded index
  offset = working_element*NGLL3_PADDED + tx;

  // global index
  iglob = d_ibool[offset] - 1 ;

  // copy from global memory to shared memory
  // each thread writes one of the NGLL^3 = 125 data points
  if (threadIdx.x < NGLL3 ){
    // copy displacement from global memory to shared memory
    load_shared_memory_displ<FORWARD_OR_ADJOINT>(&tx,&iglob,d_displ,sh_tempx,sh_tempy,sh_tempz);
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

  // synchronize all the threads (one thread for each of the NGLL grid points of the
  // current spectral element) because we need the whole element to be ready in order
  // to be able to compute the matrix products along cut planes of the 3D element below
  __syncthreads();

  // computes first matrix products
  // 1. cut-plane
  sum_hprime_xi(I,J,K,&tempx1l,&tempy1l,&tempz1l,sh_tempx,sh_tempy,sh_tempz,sh_hprime_xx);
  // 2. cut-plane
  sum_hprime_eta(I,J,K,&tempx2l,&tempy2l,&tempz2l,sh_tempx,sh_tempy,sh_tempz,sh_hprime_xx);
  // 3. cut-plane
  sum_hprime_gamma(I,J,K,&tempx3l,&tempy3l,&tempz3l,sh_tempx,sh_tempy,sh_tempz,sh_hprime_xx);

  // compute derivatives of ux, uy and uz with respect to x, y and z
  duxdxl = xixl*tempx1l + etaxl*tempx2l + gammaxl*tempx3l;
  duxdyl = xiyl*tempx1l + etayl*tempx2l + gammayl*tempx3l;
  duxdzl = xizl*tempx1l + etazl*tempx2l + gammazl*tempx3l;

  duydxl = xixl*tempy1l + etaxl*tempy2l + gammaxl*tempy3l;
  duydyl = xiyl*tempy1l + etayl*tempy2l + gammayl*tempy3l;
  duydzl = xizl*tempy1l + etazl*tempy2l + gammazl*tempy3l;

  duzdxl = xixl*tempz1l + etaxl*tempz2l + gammaxl*tempz3l;
  duzdyl = xiyl*tempz1l + etayl*tempz2l + gammayl*tempz3l;
  duzdzl = xizl*tempz1l + etazl*tempz2l + gammazl*tempz3l;

  // precompute some sums to save CPU time
  duxdxl_plus_duydyl = duxdxl + duydyl;
  duxdxl_plus_duzdzl = duxdxl + duzdzl;
  duydyl_plus_duzdzl = duydyl + duzdzl;
  duxdyl_plus_duydxl = duxdyl + duydxl;
  duzdxl_plus_duxdzl = duzdxl + duxdzl;
  duzdyl_plus_duydzl = duzdyl + duydzl;

  // computes deviatoric strain for kernel calculations
  if (COMPUTE_AND_STORE_STRAIN) {
    // save deviatoric strain for Runge-Kutta scheme
    if (threadIdx.x < NGLL3) {
      realw templ = 0.33333333333333333333f * (duxdxl + duydyl + duzdzl); // 1./3. = 0.33333
      // local storage: stresses at this current time step
      // fortran: epsilondev_xx(:,:,:,ispec) = epsilondev_xx_loc(:,:,:)
      epsilondev_xx[tx + working_element*NGLL3] = duxdxl - templ; // epsilondev_xx_loc;
      epsilondev_yy[tx + working_element*NGLL3] = duydyl - templ; // epsilondev_yy_loc;
      epsilondev_xy[tx + working_element*NGLL3] = 0.5f * duxdyl_plus_duydxl; // epsilondev_xy_loc;
      epsilondev_xz[tx + working_element*NGLL3] = 0.5f * duzdxl_plus_duxdzl; // epsilondev_xz_loc;
      epsilondev_yz[tx + working_element*NGLL3] = 0.5f * duzdyl_plus_duydzl; //epsilondev_yz_loc;
      // kernel simulations
      if (SIMULATION_TYPE == 3){
        epsilon_trace_over_3[tx + working_element*NGLL3] = templ;
      }
    } // threadIdx.x
  }

  // stress calculations

  // isotropic case
  // compute elements with an elastic isotropic rheology
  kappal = d_kappav[offset];
  mul = d_muv[offset];

  lambdalplus2mul = kappal + 1.33333333333333333333f * mul;  // 4./3. = 1.3333333
  lambdal = lambdalplus2mul - 2.0f * mul;

  // compute the six components of the stress tensor sigma
  sigma_xx = lambdalplus2mul*duxdxl + lambdal*duydyl_plus_duzdzl;
  sigma_yy = lambdalplus2mul*duydyl + lambdal*duxdxl_plus_duzdzl;
  sigma_zz = lambdalplus2mul*duzdzl + lambdal*duxdxl_plus_duydyl;

  sigma_xy = mul*duxdyl_plus_duydxl;
  sigma_xz = mul*duzdxl_plus_duxdzl;
  sigma_yz = mul*duzdyl_plus_duydzl;

  // define symmetric components (needed for non-symmetric dot product and sigma for gravity)
  sigma_yx = sigma_xy;
  sigma_zx = sigma_xz;
  sigma_zy = sigma_yz;

  if (gravity ){
    //  computes non-symmetric terms for gravity
    compute_element_gravity(tx,working_element,&iglob,d_minus_g,d_minus_deriv_gravity,
                            d_rhostore,wgll_cube,jacobianl,
                            sh_tempx,sh_tempy,sh_tempz,
                            &sigma_xx,&sigma_yy,&sigma_xz,&sigma_yz,
                            &rho_s_H1,&rho_s_H2,&rho_s_H3);
  }

  // form dot product with test vector, non-symmetric form
  // 1. cut-plane xi
  __syncthreads();
  if (threadIdx.x < NGLL3) {
    sh_tempx[tx] = jacobianl * (sigma_xx*xixl + sigma_yx*xiyl + sigma_zx*xizl); // sh_tempx1
    sh_tempy[tx] = jacobianl * (sigma_xy*xixl + sigma_yy*xiyl + sigma_zy*xizl); // sh_tempy1
    sh_tempz[tx] = jacobianl * (sigma_xz*xixl + sigma_yz*xiyl + sigma_zz*xizl); // sh_tempz1
  }
  __syncthreads();
  // 1. cut-plane xi
  sum_hprimewgll_xi(I,J,K,&tempx1l,&tempy1l,&tempz1l,sh_tempx,sh_tempy,sh_tempz,sh_hprimewgll_xx);

  // 2. cut-plane eta
  __syncthreads();
  if (threadIdx.x < NGLL3) {
    sh_tempx[tx] = jacobianl * (sigma_xx*etaxl + sigma_yx*etayl + sigma_zx*etazl); // sh_tempx2
    sh_tempy[tx] = jacobianl * (sigma_xy*etaxl + sigma_yy*etayl + sigma_zy*etazl); // sh_tempy2
    sh_tempz[tx] = jacobianl * (sigma_xz*etaxl + sigma_yz*etayl + sigma_zz*etazl); // sh_tempz2
  }
  __syncthreads();
  // 2. cut-plane eta
  sum_hprimewgll_eta(I,J,K,&tempx2l,&tempy2l,&tempz2l,sh_tempx,sh_tempy,sh_tempz,sh_hprimewgll_xx);

  // 3. cut-plane gamma
  __syncthreads();
  if (threadIdx.x < NGLL3) {
    sh_tempx[tx] = jacobianl * (sigma_xx*gammaxl + sigma_yx*gammayl + sigma_zx*gammazl); // sh_tempx3
    sh_tempy[tx] = jacobianl * (sigma_xy*gammaxl + sigma_yy*gammayl + sigma_zy*gammazl); // sh_tempy3
    sh_tempz[tx] = jacobianl * (sigma_xz*gammaxl + sigma_yz*gammayl + sigma_zz*gammazl); // sh_tempz3
  }
  __syncthreads();
  // 3. cut-plane gamma
  sum_hprimewgll_gamma(I,J,K,&tempx3l,&tempy3l,&tempz3l,sh_tempx,sh_tempy,sh_tempz,sh_hprimewgll_xx);

  // gets double weights
  fac1 = d_wgllwgll_yz[K*NGLLX+J];
  fac2 = d_wgllwgll_xz[K*NGLLX+I];
  fac3 = d_wgllwgll_xy[J*NGLLX+I];

  sum_terms1 = - (fac1*tempx1l + fac2*tempx2l + fac3*tempx3l);
  sum_terms2 = - (fac1*tempy1l + fac2*tempy2l + fac3*tempy3l);
  sum_terms3 = - (fac1*tempz1l + fac2*tempz2l + fac3*tempz3l);

  // adds gravity term
  sum_terms1 += rho_s_H1;
  sum_terms2 += rho_s_H2;
  sum_terms3 += rho_s_H3;

  // assembles acceleration array
  if (threadIdx.x < NGLL3) {

#ifdef USE_MESH_COLORING_GPU
    // no atomic operation needed, colors don't share global points between elements

#ifdef USE_TEXTURES_FIELDS
    d_accel[iglob*3]     = texfetch_accel<FORWARD_OR_ADJOINT>(iglob*3) + sum_terms1;
    d_accel[iglob*3 + 1] = texfetch_accel<FORWARD_OR_ADJOINT>(iglob*3 + 1) + sum_terms2;
    d_accel[iglob*3 + 2] = texfetch_accel<FORWARD_OR_ADJOINT>(iglob*3 + 2) + sum_terms3;
#else
    d_accel[iglob*3]     += sum_terms1;
    d_accel[iglob*3 + 1] += sum_terms2;
    d_accel[iglob*3 + 2] += sum_terms3;
#endif // USE_TEXTURES_FIELDS

#else // MESH_COLORING

    //mesh coloring
    if (use_mesh_coloring_gpu ){

      // no atomic operation needed, colors don't share global points between elements
#ifdef USE_TEXTURES_FIELDS
      d_accel[iglob*3]     = texfetch_accel<FORWARD_OR_ADJOINT>(iglob*3) + sum_terms1;
      d_accel[iglob*3 + 1] = texfetch_accel<FORWARD_OR_ADJOINT>(iglob*3 + 1) + sum_terms2;
      d_accel[iglob*3 + 2] = texfetch_accel<FORWARD_OR_ADJOINT>(iglob*3 + 2) + sum_terms3;
#else
      d_accel[iglob*3]     += sum_terms1;
      d_accel[iglob*3 + 1] += sum_terms2;
      d_accel[iglob*3 + 2] += sum_terms3;
#endif // USE_TEXTURES_FIELDS

    }else {
      atomicAdd(&d_accel[iglob*3], sum_terms1);
      atomicAdd(&d_accel[iglob*3+1], sum_terms2);
      atomicAdd(&d_accel[iglob*3+2], sum_terms3);
    } // if (use_mesh_coloring_gpu)

#endif // MESH_COLORING

  } // threadIdx.x

} // kernel_2_noatt_iso_grav_impl()

/* ----------------------------------------------------------------------------------------------- */


template<int FORWARD_OR_ADJOINT> __global__ void
#ifdef USE_LAUNCH_BOUNDS
// adds compiler specification
__launch_bounds__(NGLL3_PADDED,LAUNCH_MIN_BLOCKS)
#endif
// main kernel
Kernel_2_noatt_ani_impl(int nb_blocks_to_compute,
                        const int* d_ibool,
                        const int* d_phase_ispec_inner_elastic,const int num_phase_ispec_elastic,
                        const int d_iphase,
                        const int use_mesh_coloring_gpu,
                        realw_const_p d_displ,
                        realw_p d_accel,
                        realw* d_xix,realw* d_xiy,realw* d_xiz,
                        realw* d_etax,realw* d_etay,realw* d_etaz,
                        realw* d_gammax,realw* d_gammay,realw* d_gammaz,
                        realw_const_p d_hprime_xx,
                        realw_const_p d_hprimewgll_xx,
                        realw_const_p d_wgllwgll_xy,realw_const_p d_wgllwgll_xz,realw_const_p d_wgllwgll_yz,
                        realw_const_p d_kappav,realw_const_p d_muv,
                        const int COMPUTE_AND_STORE_STRAIN,
                        realw_p epsilondev_xx,realw_p epsilondev_yy,realw_p epsilondev_xy,
                        realw_p epsilondev_xz,realw_p epsilondev_yz,
                        realw_p epsilon_trace_over_3,
                        const int SIMULATION_TYPE,
                        const int ANISOTROPY,
                        realw* d_c11store,realw* d_c12store,realw* d_c13store,
                        realw* d_c14store,realw* d_c15store,realw* d_c16store,
                        realw* d_c22store,realw* d_c23store,realw* d_c24store,
                        realw* d_c25store,realw* d_c26store,realw* d_c33store,
                        realw* d_c34store,realw* d_c35store,realw* d_c36store,
                        realw* d_c44store,realw* d_c45store,realw* d_c46store,
                        realw* d_c55store,realw* d_c56store,realw* d_c66store,
                        const int gravity,
                        realw_const_p d_minus_g,
                        realw_const_p d_minus_deriv_gravity,
                        realw_const_p d_rhostore,
                        realw_const_p wgll_cube ){

// elastic compute kernel without attenuation for anisotropic elements
//
// holds for:
//  ATTENUATION               = .false.
//  ANISOTROPY                = .true.
//  COMPUTE_AND_STORE_STRAIN  = .true. or .false. (true for kernel simulations)


  // block-id == number of local element id in phase_ispec array
  int bx = blockIdx.y*gridDim.x+blockIdx.x;

  // checks if anything to do
  if (bx >= nb_blocks_to_compute) return;

  // thread-id == GLL node id
  // note: use only NGLL^3 = 125 active threads, plus 3 inactive/ghost threads,
  //       because we used memory padding from NGLL^3 = 125 to 128 to get coalescent memory accesses;
  //       to avoid execution branching and the need of registers to store an active state variable,
  //       the thread ids are put in valid range
  int tx = threadIdx.x;
  if (tx >= NGLL3) tx = NGLL3-1;

  int K = (tx/NGLL2);
  int J = ((tx-K*NGLL2)/NGLLX);
  int I = (tx-K*NGLL2-J*NGLLX);

  int iglob,offset;
  int working_element;

  realw tempx1l,tempx2l,tempx3l,tempy1l,tempy2l,tempy3l,tempz1l,tempz2l,tempz3l;
  realw xixl,xiyl,xizl,etaxl,etayl,etazl,gammaxl,gammayl,gammazl,jacobianl;
  realw duxdxl,duxdyl,duxdzl,duydxl,duydyl,duydzl,duzdxl,duzdyl,duzdzl;
  realw duxdxl_plus_duydyl,duxdxl_plus_duzdzl,duydyl_plus_duzdzl;
  realw duxdyl_plus_duydxl,duzdxl_plus_duxdzl,duzdyl_plus_duydzl;

  realw fac1,fac2,fac3;
  realw lambdal,mul,lambdalplus2mul,kappal;

  realw sigma_xx,sigma_yy,sigma_zz,sigma_xy,sigma_xz,sigma_yz;
  realw epsilondev_xx_loc,epsilondev_yy_loc,epsilondev_xy_loc,epsilondev_xz_loc,epsilondev_yz_loc;

  realw c11,c12,c13,c14,c15,c16,c22,c23,c24,c25,c26,c33,c34,c35,c36,c44,c45,c46,c55,c56,c66;
  realw sum_terms1,sum_terms2,sum_terms3;

  // gravity variables
  realw sigma_yx,sigma_zx,sigma_zy;
  realw rho_s_H1,rho_s_H2,rho_s_H3;

  // shared memory
  __shared__ realw sh_tempx[NGLL3];
  __shared__ realw sh_tempy[NGLL3];
  __shared__ realw sh_tempz[NGLL3];

  // note: using shared memory for hprime's improves performance
  //       (but could tradeoff with occupancy)
  __shared__ realw sh_hprime_xx[NGLL2];
  __shared__ realw sh_hprimewgll_xx[NGLL2];

  // loads hprime's into shared memory
  if (tx < NGLL2) {
    // copy hprime from global memory to shared memory
    load_shared_memory_hprime(&tx,d_hprime_xx,sh_hprime_xx);
    // copy hprime from global memory to shared memory
    load_shared_memory_hprimewgll(&tx,d_hprimewgll_xx,sh_hprimewgll_xx);
  }

  // spectral-element id
  // iphase-1 and working_element-1 for Fortran->C array conventions
#ifdef USE_MESH_COLORING_GPU
  working_element = bx;
#else
  //mesh coloring
  if (use_mesh_coloring_gpu ){
    working_element = bx;
  }else{
    // iphase-1 and working_element-1 for Fortran->C array conventions
    working_element = d_phase_ispec_inner_elastic[bx + num_phase_ispec_elastic*(d_iphase-1)] - 1;
  }
#endif
  // local padded index
  offset = working_element*NGLL3_PADDED + tx;

  // global index
  iglob = d_ibool[offset] - 1 ;

  // copy from global memory to shared memory
  // each thread writes one of the NGLL^3 = 125 data points
  if (threadIdx.x < NGLL3 ){
    // copy displacement from global memory to shared memory
    load_shared_memory_displ<FORWARD_OR_ADJOINT>(&tx,&iglob,d_displ,sh_tempx,sh_tempy,sh_tempz);
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

  // synchronize all the threads (one thread for each of the NGLL grid points of the
  // current spectral element) because we need the whole element to be ready in order
  // to be able to compute the matrix products along cut planes of the 3D element below
  __syncthreads();

  // computes first matrix products
  // 1. cut-plane
  sum_hprime_xi(I,J,K,&tempx1l,&tempy1l,&tempz1l,sh_tempx,sh_tempy,sh_tempz,sh_hprime_xx);
  // 2. cut-plane
  sum_hprime_eta(I,J,K,&tempx2l,&tempy2l,&tempz2l,sh_tempx,sh_tempy,sh_tempz,sh_hprime_xx);
  // 3. cut-plane
  sum_hprime_gamma(I,J,K,&tempx3l,&tempy3l,&tempz3l,sh_tempx,sh_tempy,sh_tempz,sh_hprime_xx);

  // compute derivatives of ux, uy and uz with respect to x, y and z
  duxdxl = xixl*tempx1l + etaxl*tempx2l + gammaxl*tempx3l;
  duxdyl = xiyl*tempx1l + etayl*tempx2l + gammayl*tempx3l;
  duxdzl = xizl*tempx1l + etazl*tempx2l + gammazl*tempx3l;

  duydxl = xixl*tempy1l + etaxl*tempy2l + gammaxl*tempy3l;
  duydyl = xiyl*tempy1l + etayl*tempy2l + gammayl*tempy3l;
  duydzl = xizl*tempy1l + etazl*tempy2l + gammazl*tempy3l;

  duzdxl = xixl*tempz1l + etaxl*tempz2l + gammaxl*tempz3l;
  duzdyl = xiyl*tempz1l + etayl*tempz2l + gammayl*tempz3l;
  duzdzl = xizl*tempz1l + etazl*tempz2l + gammazl*tempz3l;

  // precompute some sums to save CPU time
  duxdxl_plus_duydyl = duxdxl + duydyl;
  duxdxl_plus_duzdzl = duxdxl + duzdzl;
  duydyl_plus_duzdzl = duydyl + duzdzl;
  duxdyl_plus_duydxl = duxdyl + duydxl;
  duzdxl_plus_duxdzl = duzdxl + duxdzl;
  duzdyl_plus_duydzl = duzdyl + duydzl;

  // computes deviatoric strain for kernel calculations
  if (COMPUTE_AND_STORE_STRAIN) {
    realw templ = 0.33333333333333333333f * (duxdxl + duydyl + duzdzl); // 1./3. = 0.33333
    // local storage: stresses at this current time step
    epsilondev_xx_loc = duxdxl - templ;
    epsilondev_yy_loc = duydyl - templ;
    epsilondev_xy_loc = 0.5f * duxdyl_plus_duydxl;
    epsilondev_xz_loc = 0.5f * duzdxl_plus_duxdzl;
    epsilondev_yz_loc = 0.5f * duzdyl_plus_duydzl;

    if (SIMULATION_TYPE == 3){
      if (threadIdx.x < NGLL3) {
        epsilon_trace_over_3[tx + working_element*NGLL3] = templ;
      }
    }
  }

  // full anisotropic case, stress calculations
  if (ANISOTROPY){
    c11 = d_c11store[offset];
    c12 = d_c12store[offset];
    c13 = d_c13store[offset];
    c14 = d_c14store[offset];
    c15 = d_c15store[offset];
    c16 = d_c16store[offset];
    c22 = d_c22store[offset];
    c23 = d_c23store[offset];
    c24 = d_c24store[offset];
    c25 = d_c25store[offset];
    c26 = d_c26store[offset];
    c33 = d_c33store[offset];
    c34 = d_c34store[offset];
    c35 = d_c35store[offset];
    c36 = d_c36store[offset];
    c44 = d_c44store[offset];
    c45 = d_c45store[offset];
    c46 = d_c46store[offset];
    c55 = d_c55store[offset];
    c56 = d_c56store[offset];
    c66 = d_c66store[offset];

    sigma_xx = c11*duxdxl + c16*duxdyl_plus_duydxl + c12*duydyl +
               c15*duzdxl_plus_duxdzl + c14*duzdyl_plus_duydzl + c13*duzdzl;
    sigma_yy = c12*duxdxl + c26*duxdyl_plus_duydxl + c22*duydyl +
               c25*duzdxl_plus_duxdzl + c24*duzdyl_plus_duydzl + c23*duzdzl;
    sigma_zz = c13*duxdxl + c36*duxdyl_plus_duydxl + c23*duydyl +
               c35*duzdxl_plus_duxdzl + c34*duzdyl_plus_duydzl + c33*duzdzl;
    sigma_xy = c16*duxdxl + c66*duxdyl_plus_duydxl + c26*duydyl +
               c56*duzdxl_plus_duxdzl + c46*duzdyl_plus_duydzl + c36*duzdzl;
    sigma_xz = c15*duxdxl + c56*duxdyl_plus_duydxl + c25*duydyl +
               c55*duzdxl_plus_duxdzl + c45*duzdyl_plus_duydzl + c35*duzdzl;
    sigma_yz = c14*duxdxl + c46*duxdyl_plus_duydxl + c24*duydyl +
               c45*duzdxl_plus_duxdzl + c44*duzdyl_plus_duydzl + c34*duzdzl;

  }else{

    // isotropic case

    // compute elements with an elastic isotropic rheology
    kappal = d_kappav[offset];
    mul = d_muv[offset];

    lambdalplus2mul = kappal + 1.33333333333333333333f * mul;  // 4./3. = 1.3333333
    lambdal = lambdalplus2mul - 2.0f * mul;

    // compute the six components of the stress tensor sigma
    sigma_xx = lambdalplus2mul*duxdxl + lambdal*duydyl_plus_duzdzl;
    sigma_yy = lambdalplus2mul*duydyl + lambdal*duxdxl_plus_duzdzl;
    sigma_zz = lambdalplus2mul*duzdzl + lambdal*duxdxl_plus_duydyl;

    sigma_xy = mul*duxdyl_plus_duydxl;
    sigma_xz = mul*duzdxl_plus_duxdzl;
    sigma_yz = mul*duzdyl_plus_duydzl;
  }

  // define symmetric components (needed for non-symmetric dot product and sigma for gravity)
  sigma_yx = sigma_xy;
  sigma_zx = sigma_xz;
  sigma_zy = sigma_yz;

  if (gravity ){
    //  computes non-symmetric terms for gravity
    compute_element_gravity(tx,working_element,&iglob,d_minus_g,d_minus_deriv_gravity,
                            d_rhostore,wgll_cube,jacobianl,
                            sh_tempx,sh_tempy,sh_tempz,
                            &sigma_xx,&sigma_yy,&sigma_xz,&sigma_yz,
                            &rho_s_H1,&rho_s_H2,&rho_s_H3);
  }

  // form dot product with test vector, non-symmetric form
  // 1. cut-plane xi
  __syncthreads();
  if (threadIdx.x < NGLL3) {
    sh_tempx[tx] = jacobianl * (sigma_xx*xixl + sigma_yx*xiyl + sigma_zx*xizl); // sh_tempx1
    sh_tempy[tx] = jacobianl * (sigma_xy*xixl + sigma_yy*xiyl + sigma_zy*xizl); // sh_tempy1
    sh_tempz[tx] = jacobianl * (sigma_xz*xixl + sigma_yz*xiyl + sigma_zz*xizl); // sh_tempz1
  }
  __syncthreads();
  // 1. cut-plane xi
  sum_hprimewgll_xi(I,J,K,&tempx1l,&tempy1l,&tempz1l,sh_tempx,sh_tempy,sh_tempz,sh_hprimewgll_xx);

  // 2. cut-plane eta
  __syncthreads();
  if (threadIdx.x < NGLL3) {
    sh_tempx[tx] = jacobianl * (sigma_xx*etaxl + sigma_yx*etayl + sigma_zx*etazl); // sh_tempx2
    sh_tempy[tx] = jacobianl * (sigma_xy*etaxl + sigma_yy*etayl + sigma_zy*etazl); // sh_tempy2
    sh_tempz[tx] = jacobianl * (sigma_xz*etaxl + sigma_yz*etayl + sigma_zz*etazl); // sh_tempz2
  }
  __syncthreads();
  // 2. cut-plane eta
  sum_hprimewgll_eta(I,J,K,&tempx2l,&tempy2l,&tempz2l,sh_tempx,sh_tempy,sh_tempz,sh_hprimewgll_xx);

  // 3. cut-plane gamma
  __syncthreads();
  if (threadIdx.x < NGLL3) {
    sh_tempx[tx] = jacobianl * (sigma_xx*gammaxl + sigma_yx*gammayl + sigma_zx*gammazl); // sh_tempx3
    sh_tempy[tx] = jacobianl * (sigma_xy*gammaxl + sigma_yy*gammayl + sigma_zy*gammazl); // sh_tempy3
    sh_tempz[tx] = jacobianl * (sigma_xz*gammaxl + sigma_yz*gammayl + sigma_zz*gammazl); // sh_tempz3
  }
  __syncthreads();
  // 3. cut-plane gamma
  sum_hprimewgll_gamma(I,J,K,&tempx3l,&tempy3l,&tempz3l,sh_tempx,sh_tempy,sh_tempz,sh_hprimewgll_xx);

  // gets double weights
  fac1 = d_wgllwgll_yz[K*NGLLX+J];
  fac2 = d_wgllwgll_xz[K*NGLLX+I];
  fac3 = d_wgllwgll_xy[J*NGLLX+I];

  sum_terms1 = - (fac1*tempx1l + fac2*tempx2l + fac3*tempx3l);
  sum_terms2 = - (fac1*tempy1l + fac2*tempy2l + fac3*tempy3l);
  sum_terms3 = - (fac1*tempz1l + fac2*tempz2l + fac3*tempz3l);

  // adds gravity term
  if (gravity ){
    sum_terms1 += rho_s_H1;
    sum_terms2 += rho_s_H2;
    sum_terms3 += rho_s_H3;
  }

  // assembles acceleration array
  if (threadIdx.x < NGLL3) {

#ifdef USE_MESH_COLORING_GPU
    // no atomic operation needed, colors don't share global points between elements

#ifdef USE_TEXTURES_FIELDS
    d_accel[iglob*3]     = texfetch_accel<FORWARD_OR_ADJOINT>(iglob*3) + sum_terms1;
    d_accel[iglob*3 + 1] = texfetch_accel<FORWARD_OR_ADJOINT>(iglob*3 + 1) + sum_terms2;
    d_accel[iglob*3 + 2] = texfetch_accel<FORWARD_OR_ADJOINT>(iglob*3 + 2) + sum_terms3;
#else
    d_accel[iglob*3]     += sum_terms1;
    d_accel[iglob*3 + 1] += sum_terms2;
    d_accel[iglob*3 + 2] += sum_terms3;
#endif // USE_TEXTURES_FIELDS

#else // MESH_COLORING

    //mesh coloring
    if (use_mesh_coloring_gpu ){

      // no atomic operation needed, colors don't share global points between elements
#ifdef USE_TEXTURES_FIELDS
      d_accel[iglob*3]     = texfetch_accel<FORWARD_OR_ADJOINT>(iglob*3) + sum_terms1;
      d_accel[iglob*3 + 1] = texfetch_accel<FORWARD_OR_ADJOINT>(iglob*3 + 1) + sum_terms2;
      d_accel[iglob*3 + 2] = texfetch_accel<FORWARD_OR_ADJOINT>(iglob*3 + 2) + sum_terms3;
#else
      d_accel[iglob*3]     += sum_terms1;
      d_accel[iglob*3 + 1] += sum_terms2;
      d_accel[iglob*3 + 2] += sum_terms3;
#endif // USE_TEXTURES_FIELDS

    }else {
      atomicAdd(&d_accel[iglob*3], sum_terms1);
      atomicAdd(&d_accel[iglob*3+1], sum_terms2);
      atomicAdd(&d_accel[iglob*3+2], sum_terms3);
    } // if (use_mesh_coloring_gpu)

#endif // MESH_COLORING

    // save deviatoric strain for Runge-Kutta scheme
    if (COMPUTE_AND_STORE_STRAIN ){
      // fortran: epsilondev_xx(:,:,:,ispec) = epsilondev_xx_loc(:,:,:)
      epsilondev_xx[tx + working_element*NGLL3] = epsilondev_xx_loc;
      epsilondev_yy[tx + working_element*NGLL3] = epsilondev_yy_loc;
      epsilondev_xy[tx + working_element*NGLL3] = epsilondev_xy_loc;
      epsilondev_xz[tx + working_element*NGLL3] = epsilondev_xz_loc;
      epsilondev_yz[tx + working_element*NGLL3] = epsilondev_yz_loc;
    }

  } // threadIdx.x

} // kernel_2_noatt_ani_impl()


/* ----------------------------------------------------------------------------------------------- */

// kernel with attenuation
//
// we use templates to distinguish between calls with forward or adjoint texture fields

template<int FORWARD_OR_ADJOINT> __global__ void
#ifdef USE_LAUNCH_BOUNDS
// adds compiler specification
__launch_bounds__(NGLL3_PADDED,LAUNCH_MIN_BLOCKS)
#endif
// main kernel
Kernel_2_att_impl(int nb_blocks_to_compute,
                  const int* d_ibool,
                  const int* d_phase_ispec_inner_elastic,const int num_phase_ispec_elastic,
                  const int d_iphase,
                  const int use_mesh_coloring_gpu,
                  realw d_deltat,
                  realw_const_p d_displ,
                  realw_const_p d_veloc,
                  realw_p d_accel,
                  realw_const_p d_xix,realw_const_p d_xiy,realw_const_p d_xiz,
                  realw_const_p d_etax,realw_const_p d_etay,realw_const_p d_etaz,
                  realw_const_p d_gammax,realw_const_p d_gammay,realw_const_p d_gammaz,
                  realw_const_p d_hprime_xx,
                  realw_const_p d_hprimewgll_xx,
                  realw_const_p d_wgllwgll_xy,realw_const_p d_wgllwgll_xz,realw_const_p d_wgllwgll_yz,
                  realw_const_p d_kappav,realw_const_p d_muv,
                  realw_p epsilondev_xx,realw_p epsilondev_yy,realw_p epsilondev_xy,
                  realw_p epsilondev_xz,realw_p epsilondev_yz,
                  realw_p epsilon_trace_over_3,
                  const int SIMULATION_TYPE,
                  const int NSPEC,
                  realw_const_p one_minus_sum_beta,realw_const_p factor_common,
                  realw_p R_xx,realw_p R_yy,realw_p R_xy,realw_p R_xz,realw_p R_yz,
                  realw_const_p alphaval,realw_const_p betaval,realw_const_p gammaval,
                  const int ANISOTROPY,
                  realw_const_p d_c11store,realw_const_p d_c12store,realw_const_p d_c13store,
                  realw_const_p d_c14store,realw_const_p d_c15store,realw_const_p d_c16store,
                  realw_const_p d_c22store,realw_const_p d_c23store,realw_const_p d_c24store,
                  realw_const_p d_c25store,realw_const_p d_c26store,realw_const_p d_c33store,
                  realw_const_p d_c34store,realw_const_p d_c35store,realw_const_p d_c36store,
                  realw_const_p d_c44store,realw_const_p d_c45store,realw_const_p d_c46store,
                  realw_const_p d_c55store,realw_const_p d_c56store,realw_const_p d_c66store,
                  const int gravity,
                  realw_const_p d_minus_g,
                  realw_const_p d_minus_deriv_gravity,
                  realw_const_p d_rhostore,
                  realw_const_p wgll_cube){


// elastic compute kernel with attenuation
// holds for: ATTENUATION = .true.
//            COMPUTE_AND_STORE_STRAIN = .true. (always true for attenuation)

  int bx = blockIdx.y*gridDim.x+blockIdx.x;
  int tx = threadIdx.x;

  int K = (tx/NGLL2);
  int J = ((tx-K*NGLL2)/NGLLX);
  int I = (tx-K*NGLL2-J*NGLLX);

  unsigned short int active;
  int iglob,offset;
  int working_element;

  realw tempx1l,tempx2l,tempx3l,tempy1l,tempy2l,tempy3l,tempz1l,tempz2l,tempz3l;
  realw xixl,xiyl,xizl,etaxl,etayl,etazl,gammaxl,gammayl,gammazl,jacobianl;
  realw duxdxl,duxdyl,duxdzl,duydxl,duydyl,duydzl,duzdxl,duzdyl,duzdzl;
  realw duxdxl_plus_duydyl,duxdxl_plus_duzdzl,duydyl_plus_duzdzl;
  realw duxdyl_plus_duydxl,duzdxl_plus_duxdzl,duzdyl_plus_duydzl;
  realw templ;

  realw tempx1l_att,tempx2l_att,tempx3l_att,tempy1l_att,tempy2l_att,tempy3l_att,tempz1l_att,tempz2l_att,tempz3l_att;
  realw duxdxl_att,duxdyl_att,duxdzl_att,duydxl_att,duydyl_att,duydzl_att,duzdxl_att,duzdyl_att,duzdzl_att;
  realw duxdyl_plus_duydxl_att,duzdxl_plus_duxdzl_att,duzdyl_plus_duydzl_att;

  realw fac1,fac2,fac3;
  realw lambdal,mul,lambdalplus2mul,kappal;

  realw sigma_xx,sigma_yy,sigma_zz,sigma_xy,sigma_xz,sigma_yz;
  realw epsilondev_xx_loc,epsilondev_yy_loc,epsilondev_xy_loc,epsilondev_xz_loc,epsilondev_yz_loc;

  realw c11,c12,c13,c14,c15,c16,c22,c23,c24,c25,c26,c33,c34,c35,c36,c44,c45,c46,c55,c56,c66;
  realw sum_terms1,sum_terms2,sum_terms3;

  // gravity variables
  realw sigma_yx,sigma_zx,sigma_zy;
  realw rho_s_H1,rho_s_H2,rho_s_H3;

#ifndef MANUALLY_UNROLLED_LOOPS
  int l;
#endif

  __shared__ realw sh_tempx1[NGLL3];
  __shared__ realw sh_tempx2[NGLL3];
  __shared__ realw sh_tempx3[NGLL3];

  __shared__ realw sh_tempy1[NGLL3];
  __shared__ realw sh_tempy2[NGLL3];
  __shared__ realw sh_tempy3[NGLL3];

  __shared__ realw sh_tempz1[NGLL3];
  __shared__ realw sh_tempz2[NGLL3];
  __shared__ realw sh_tempz3[NGLL3];

  __shared__ realw sh_hprime_xx[NGLL2];

  __shared__ realw s_dummyx_loc[NGLL3];
  __shared__ realw s_dummyy_loc[NGLL3];
  __shared__ realw s_dummyz_loc[NGLL3];

  // re-assigns shared array to decrease shared memory usage
  // note: this will re-use s_temp arrays from above to save shared memory
  realw* s_dummyx_loc_att = (realw*) sh_tempx1;
  realw* s_dummyy_loc_att = (realw*) sh_tempx2;
  realw* s_dummyz_loc_att = (realw*) sh_tempx3;

  // use only NGLL^3 = 125 active threads, plus 3 inactive/ghost threads,
  // because we used memory padding from NGLL^3 = 125 to 128 to get coalescent memory accesses
  active = (tx < NGLL3 && bx < nb_blocks_to_compute) ? 1:0;

  // copy from global memory to shared memory
  // each thread writes one of the NGLL^3 = 125 data points
  if (active ){

#ifdef USE_MESH_COLORING_GPU
    working_element = bx;
#else
    //mesh coloring
    if (use_mesh_coloring_gpu ){
      working_element = bx;
    }else{
      // iphase-1 and working_element-1 for Fortran->C array conventions
      working_element = d_phase_ispec_inner_elastic[bx + num_phase_ispec_elastic*(d_iphase-1)]-1;
    }
#endif
    // local padded index
    offset = working_element*NGLL3_PADDED + tx;

    // global index
    iglob = d_ibool[offset] - 1;

    // copy displacement from global memory to shared memory
    load_shared_memory_displ<FORWARD_OR_ADJOINT>(&tx,&iglob,d_displ,s_dummyx_loc,s_dummyy_loc,s_dummyz_loc);

// JC JC here we will need to add GPU support for the new C-PML routines

    // attenuation
    // use first order Taylor expansion of displacement for local storage of stresses
    // at this current time step, to fix attenuation in a consistent way
#ifdef USE_TEXTURES_FIELDS
    s_dummyx_loc_att[tx] = s_dummyx_loc[tx] + d_deltat * texfetch_veloc<FORWARD_OR_ADJOINT>(iglob*3);
    s_dummyy_loc_att[tx] = s_dummyy_loc[tx] + d_deltat * texfetch_veloc<FORWARD_OR_ADJOINT>(iglob*3 + 1);
    s_dummyz_loc_att[tx] = s_dummyz_loc[tx] + d_deltat * texfetch_veloc<FORWARD_OR_ADJOINT>(iglob*3 + 2);
#else
    s_dummyx_loc_att[tx] = s_dummyx_loc[tx] + d_deltat * d_veloc[iglob*3];
    s_dummyy_loc_att[tx] = s_dummyy_loc[tx] + d_deltat * d_veloc[iglob*3 + 1];
    s_dummyz_loc_att[tx] = s_dummyz_loc[tx] + d_deltat * d_veloc[iglob*3 + 2];
#endif
  }// active

  // loads hprime's into shared memory
  if (tx < NGLL2) {
    // copy hprime from global memory to shared memory
    load_shared_memory_hprime(&tx,d_hprime_xx,sh_hprime_xx);
  }

  // synchronize all the threads (one thread for each of the NGLL grid points of the
  // current spectral element) because we need the whole element to be ready in order
  // to be able to compute the matrix products along cut planes of the 3D element below
  __syncthreads();

  if (active ){

#ifndef MANUALLY_UNROLLED_LOOPS

    tempx1l = 0.f;
    tempx2l = 0.f;
    tempx3l = 0.f;

    tempy1l = 0.f;
    tempy2l = 0.f;
    tempy3l = 0.f;

    tempz1l = 0.f;
    tempz2l = 0.f;
    tempz3l = 0.f;

    for (l=0;l<NGLLX;l++) {
      //assumes that hprime_xx = hprime_yy = hprime_zz
      fac1 = sh_hprime_xx[l*NGLLX+I];
      tempx1l += s_dummyx_loc[K*NGLL2+J*NGLLX+l]*fac1;
      tempy1l += s_dummyy_loc[K*NGLL2+J*NGLLX+l]*fac1;
      tempz1l += s_dummyz_loc[K*NGLL2+J*NGLLX+l]*fac1;

      fac2 = sh_hprime_xx[l*NGLLX+J];
      tempx2l += s_dummyx_loc[K*NGLL2+l*NGLLX+I]*fac2;
      tempy2l += s_dummyy_loc[K*NGLL2+l*NGLLX+I]*fac2;
      tempz2l += s_dummyz_loc[K*NGLL2+l*NGLLX+I]*fac2;

      fac3 = sh_hprime_xx[l*NGLLX+K];
      tempx3l += s_dummyx_loc[l*NGLL2+J*NGLLX+I]*fac3;
      tempy3l += s_dummyy_loc[l*NGLL2+J*NGLLX+I]*fac3;
      tempz3l += s_dummyz_loc[l*NGLL2+J*NGLLX+I]*fac3;
    }

// JC JC here we will need to add GPU support for the new C-PML routines

    // attenuation
    // temporary variables used for fixing attenuation in a consistent way
    tempx1l_att = 0.f;
    tempx2l_att = 0.f;
    tempx3l_att = 0.f;

    tempy1l_att = 0.f;
    tempy2l_att = 0.f;
    tempy3l_att = 0.f;

    tempz1l_att = 0.f;
    tempz2l_att = 0.f;
    tempz3l_att = 0.f;

    for (l=0;l<NGLLX;l++) {
      fac1 = sh_hprime_xx[l*NGLLX+I];
      tempx1l_att += s_dummyx_loc_att[K*NGLL2+J*NGLLX+l]*fac1;
      tempy1l_att += s_dummyy_loc_att[K*NGLL2+J*NGLLX+l]*fac1;
      tempz1l_att += s_dummyz_loc_att[K*NGLL2+J*NGLLX+l]*fac1;

      fac2 = sh_hprime_xx[l*NGLLX+J];
      tempx2l_att += s_dummyx_loc_att[K*NGLL2+l*NGLLX+I]*fac2;
      tempy2l_att += s_dummyy_loc_att[K*NGLL2+l*NGLLX+I]*fac2;
      tempz2l_att += s_dummyz_loc_att[K*NGLL2+l*NGLLX+I]*fac2;

      fac3 = sh_hprime_xx[l*NGLLX+K];
      tempx3l_att += s_dummyx_loc_att[l*NGLL2+J*NGLLX+I]*fac3;
      tempy3l_att += s_dummyy_loc_att[l*NGLL2+J*NGLLX+I]*fac3;
      tempz3l_att += s_dummyz_loc_att[l*NGLL2+J*NGLLX+I]*fac3;
    }

#else
    tempx1l = s_dummyx_loc[K*NGLL2+J*NGLLX]*d_hprime_xx[I]
            + s_dummyx_loc[K*NGLL2+J*NGLLX+1]*d_hprime_xx[NGLLX+I]
            + s_dummyx_loc[K*NGLL2+J*NGLLX+2]*d_hprime_xx[2*NGLLX+I]
            + s_dummyx_loc[K*NGLL2+J*NGLLX+3]*d_hprime_xx[3*NGLLX+I]
            + s_dummyx_loc[K*NGLL2+J*NGLLX+4]*d_hprime_xx[4*NGLLX+I];

    tempy1l = s_dummyy_loc[K*NGLL2+J*NGLLX]*d_hprime_xx[I]
            + s_dummyy_loc[K*NGLL2+J*NGLLX+1]*d_hprime_xx[NGLLX+I]
            + s_dummyy_loc[K*NGLL2+J*NGLLX+2]*d_hprime_xx[2*NGLLX+I]
            + s_dummyy_loc[K*NGLL2+J*NGLLX+3]*d_hprime_xx[3*NGLLX+I]
            + s_dummyy_loc[K*NGLL2+J*NGLLX+4]*d_hprime_xx[4*NGLLX+I];

    tempz1l = s_dummyz_loc[K*NGLL2+J*NGLLX]*d_hprime_xx[I]
            + s_dummyz_loc[K*NGLL2+J*NGLLX+1]*d_hprime_xx[NGLLX+I]
            + s_dummyz_loc[K*NGLL2+J*NGLLX+2]*d_hprime_xx[2*NGLLX+I]
            + s_dummyz_loc[K*NGLL2+J*NGLLX+3]*d_hprime_xx[3*NGLLX+I]
            + s_dummyz_loc[K*NGLL2+J*NGLLX+4]*d_hprime_xx[4*NGLLX+I];

    tempx2l = s_dummyx_loc[K*NGLL2+I]*d_hprime_xx[J]
            + s_dummyx_loc[K*NGLL2+NGLLX+I]*d_hprime_xx[NGLLX+J]
            + s_dummyx_loc[K*NGLL2+2*NGLLX+I]*d_hprime_xx[2*NGLLX+J]
            + s_dummyx_loc[K*NGLL2+3*NGLLX+I]*d_hprime_xx[3*NGLLX+J]
            + s_dummyx_loc[K*NGLL2+4*NGLLX+I]*d_hprime_xx[4*NGLLX+J];

    tempy2l = s_dummyy_loc[K*NGLL2+I]*d_hprime_xx[J]
            + s_dummyy_loc[K*NGLL2+NGLLX+I]*d_hprime_xx[NGLLX+J]
            + s_dummyy_loc[K*NGLL2+2*NGLLX+I]*d_hprime_xx[2*NGLLX+J]
            + s_dummyy_loc[K*NGLL2+3*NGLLX+I]*d_hprime_xx[3*NGLLX+J]
            + s_dummyy_loc[K*NGLL2+4*NGLLX+I]*d_hprime_xx[4*NGLLX+J];

    tempz2l = s_dummyz_loc[K*NGLL2+I]*d_hprime_xx[J]
            + s_dummyz_loc[K*NGLL2+NGLLX+I]*d_hprime_xx[NGLLX+J]
            + s_dummyz_loc[K*NGLL2+2*NGLLX+I]*d_hprime_xx[2*NGLLX+J]
            + s_dummyz_loc[K*NGLL2+3*NGLLX+I]*d_hprime_xx[3*NGLLX+J]
            + s_dummyz_loc[K*NGLL2+4*NGLLX+I]*d_hprime_xx[4*NGLLX+J];

    tempx3l = s_dummyx_loc[J*NGLLX+I]*d_hprime_xx[K]
            + s_dummyx_loc[NGLL2+J*NGLLX+I]*d_hprime_xx[NGLLX+K]
            + s_dummyx_loc[2*NGLL2+J*NGLLX+I]*d_hprime_xx[2*NGLLX+K]
            + s_dummyx_loc[3*NGLL2+J*NGLLX+I]*d_hprime_xx[3*NGLLX+K]
            + s_dummyx_loc[4*NGLL2+J*NGLLX+I]*d_hprime_xx[4*NGLLX+K];

    tempy3l = s_dummyy_loc[J*NGLLX+I]*d_hprime_xx[K]
            + s_dummyy_loc[NGLL2+J*NGLLX+I]*d_hprime_xx[NGLLX+K]
            + s_dummyy_loc[2*NGLL2+J*NGLLX+I]*d_hprime_xx[2*NGLLX+K]
            + s_dummyy_loc[3*NGLL2+J*NGLLX+I]*d_hprime_xx[3*NGLLX+K]
            + s_dummyy_loc[4*NGLL2+J*NGLLX+I]*d_hprime_xx[4*NGLLX+K];

    tempz3l = s_dummyz_loc[J*NGLLX+I]*d_hprime_xx[K]
            + s_dummyz_loc[NGLL2+J*NGLLX+I]*d_hprime_xx[NGLLX+K]
            + s_dummyz_loc[2*NGLL2+J*NGLLX+I]*d_hprime_xx[2*NGLLX+K]
            + s_dummyz_loc[3*NGLL2+J*NGLLX+I]*d_hprime_xx[3*NGLLX+K]
            + s_dummyz_loc[4*NGLL2+J*NGLLX+I]*d_hprime_xx[4*NGLLX+K];

// JC JC here we will need to add GPU support for the new C-PML routines

    // attenuation
    // temporary variables used for fixing attenuation in a consistent way
    tempx1l_att = s_dummyx_loc_att[K*NGLL2+J*NGLLX]*d_hprime_xx[I]
                  + s_dummyx_loc_att[K*NGLL2+J*NGLLX+1]*d_hprime_xx[NGLLX+I]
                  + s_dummyx_loc_att[K*NGLL2+J*NGLLX+2]*d_hprime_xx[2*NGLLX+I]
                  + s_dummyx_loc_att[K*NGLL2+J*NGLLX+3]*d_hprime_xx[3*NGLLX+I]
                  + s_dummyx_loc_att[K*NGLL2+J*NGLLX+4]*d_hprime_xx[4*NGLLX+I];

    tempy1l_att = s_dummyy_loc_att[K*NGLL2+J*NGLLX]*d_hprime_xx[I]
                  + s_dummyy_loc_att[K*NGLL2+J*NGLLX+1]*d_hprime_xx[NGLLX+I]
                  + s_dummyy_loc_att[K*NGLL2+J*NGLLX+2]*d_hprime_xx[2*NGLLX+I]
                  + s_dummyy_loc_att[K*NGLL2+J*NGLLX+3]*d_hprime_xx[3*NGLLX+I]
                  + s_dummyy_loc_att[K*NGLL2+J*NGLLX+4]*d_hprime_xx[4*NGLLX+I];

    tempz1l_att = s_dummyz_loc_att[K*NGLL2+J*NGLLX]*d_hprime_xx[I]
                  + s_dummyz_loc_att[K*NGLL2+J*NGLLX+1]*d_hprime_xx[NGLLX+I]
                  + s_dummyz_loc_att[K*NGLL2+J*NGLLX+2]*d_hprime_xx[2*NGLLX+I]
                  + s_dummyz_loc_att[K*NGLL2+J*NGLLX+3]*d_hprime_xx[3*NGLLX+I]
                  + s_dummyz_loc_att[K*NGLL2+J*NGLLX+4]*d_hprime_xx[4*NGLLX+I];

    tempx2l_att = s_dummyx_loc_att[K*NGLL2+I]*d_hprime_xx[J]
                  + s_dummyx_loc_att[K*NGLL2+NGLLX+I]*d_hprime_xx[NGLLX+J]
                  + s_dummyx_loc_att[K*NGLL2+2*NGLLX+I]*d_hprime_xx[2*NGLLX+J]
                  + s_dummyx_loc_att[K*NGLL2+3*NGLLX+I]*d_hprime_xx[3*NGLLX+J]
                  + s_dummyx_loc_att[K*NGLL2+4*NGLLX+I]*d_hprime_xx[4*NGLLX+J];

    tempy2l_att = s_dummyy_loc_att[K*NGLL2+I]*d_hprime_xx[J]
                  + s_dummyy_loc_att[K*NGLL2+NGLLX+I]*d_hprime_xx[NGLLX+J]
                  + s_dummyy_loc_att[K*NGLL2+2*NGLLX+I]*d_hprime_xx[2*NGLLX+J]
                  + s_dummyy_loc_att[K*NGLL2+3*NGLLX+I]*d_hprime_xx[3*NGLLX+J]
                  + s_dummyy_loc_att[K*NGLL2+4*NGLLX+I]*d_hprime_xx[4*NGLLX+J];

    tempz2l_att = s_dummyz_loc_att[K*NGLL2+I]*d_hprime_xx[J]
                  + s_dummyz_loc_att[K*NGLL2+NGLLX+I]*d_hprime_xx[NGLLX+J]
                  + s_dummyz_loc_att[K*NGLL2+2*NGLLX+I]*d_hprime_xx[2*NGLLX+J]
                  + s_dummyz_loc_att[K*NGLL2+3*NGLLX+I]*d_hprime_xx[3*NGLLX+J]
                  + s_dummyz_loc_att[K*NGLL2+4*NGLLX+I]*d_hprime_xx[4*NGLLX+J];

    tempx3l_att = s_dummyx_loc_att[J*NGLLX+I]*d_hprime_xx[K]
                  + s_dummyx_loc_att[NGLL2+J*NGLLX+I]*d_hprime_xx[NGLLX+K]
                  + s_dummyx_loc_att[2*NGLL2+J*NGLLX+I]*d_hprime_xx[2*NGLLX+K]
                  + s_dummyx_loc_att[3*NGLL2+J*NGLLX+I]*d_hprime_xx[3*NGLLX+K]
                  + s_dummyx_loc_att[4*NGLL2+J*NGLLX+I]*d_hprime_xx[4*NGLLX+K];

    tempy3l_att = s_dummyy_loc_att[J*NGLLX+I]*d_hprime_xx[K]
                  + s_dummyy_loc_att[NGLL2+J*NGLLX+I]*d_hprime_xx[NGLLX+K]
                  + s_dummyy_loc_att[2*NGLL2+J*NGLLX+I]*d_hprime_xx[2*NGLLX+K]
                  + s_dummyy_loc_att[3*NGLL2+J*NGLLX+I]*d_hprime_xx[3*NGLLX+K]
                  + s_dummyy_loc_att[4*NGLL2+J*NGLLX+I]*d_hprime_xx[4*NGLLX+K];

    tempz3l_att = s_dummyz_loc_att[J*NGLLX+I]*d_hprime_xx[K]
                  + s_dummyz_loc_att[NGLL2+J*NGLLX+I]*d_hprime_xx[NGLLX+K]
                  + s_dummyz_loc_att[2*NGLL2+J*NGLLX+I]*d_hprime_xx[2*NGLLX+K]
                  + s_dummyz_loc_att[3*NGLL2+J*NGLLX+I]*d_hprime_xx[3*NGLLX+K]
                  + s_dummyz_loc_att[4*NGLL2+J*NGLLX+I]*d_hprime_xx[4*NGLLX+K];
#endif

    // compute derivatives of ux, uy and uz with respect to x, y and z
    xixl = get_global_cr( &d_xix[offset] );
    xiyl = get_global_cr( &d_xiy[offset] );
    xizl = get_global_cr( &d_xiz[offset] );
    etaxl = get_global_cr( &d_etax[offset] );
    etayl = get_global_cr( &d_etay[offset] );
    etazl = get_global_cr( &d_etaz[offset] );
    gammaxl = get_global_cr( &d_gammax[offset] );
    gammayl = get_global_cr( &d_gammay[offset] );
    gammazl = get_global_cr( &d_gammaz[offset] );

    jacobianl = 1.0f / (xixl*(etayl*gammazl-etazl*gammayl)
                       -xiyl*(etaxl*gammazl-etazl*gammaxl)
                       +xizl*(etaxl*gammayl-etayl*gammaxl));

    duxdxl = xixl*tempx1l + etaxl*tempx2l + gammaxl*tempx3l;
    duxdyl = xiyl*tempx1l + etayl*tempx2l + gammayl*tempx3l;
    duxdzl = xizl*tempx1l + etazl*tempx2l + gammazl*tempx3l;

    duydxl = xixl*tempy1l + etaxl*tempy2l + gammaxl*tempy3l;
    duydyl = xiyl*tempy1l + etayl*tempy2l + gammayl*tempy3l;
    duydzl = xizl*tempy1l + etazl*tempy2l + gammazl*tempy3l;

    duzdxl = xixl*tempz1l + etaxl*tempz2l + gammaxl*tempz3l;
    duzdyl = xiyl*tempz1l + etayl*tempz2l + gammayl*tempz3l;
    duzdzl = xizl*tempz1l + etazl*tempz2l + gammazl*tempz3l;

// JC JC here we will need to add GPU support for the new C-PML routines

    // precompute some sums to save CPU time
    duxdxl_plus_duydyl = duxdxl + duydyl;
    duxdxl_plus_duzdzl = duxdxl + duzdzl;
    duydyl_plus_duzdzl = duydyl + duzdzl;
    duxdyl_plus_duydxl = duxdyl + duydxl;
    duzdxl_plus_duxdzl = duzdxl + duxdzl;
    duzdyl_plus_duydzl = duzdyl + duydzl;

// JC JC here we will need to add GPU support for the new C-PML routines

    // attenuation
    // temporary variables used for fixing attenuation in a consistent way
    duxdxl_att = xixl*tempx1l_att + etaxl*tempx2l_att + gammaxl*tempx3l_att;
    duxdyl_att = xiyl*tempx1l_att + etayl*tempx2l_att + gammayl*tempx3l_att;
    duxdzl_att = xizl*tempx1l_att + etazl*tempx2l_att + gammazl*tempx3l_att;

    duydxl_att = xixl*tempy1l_att + etaxl*tempy2l_att + gammaxl*tempy3l_att;
    duydyl_att = xiyl*tempy1l_att + etayl*tempy2l_att + gammayl*tempy3l_att;
    duydzl_att = xizl*tempy1l_att + etazl*tempy2l_att + gammazl*tempy3l_att;

    duzdxl_att = xixl*tempz1l_att + etaxl*tempz2l_att + gammaxl*tempz3l_att;
    duzdyl_att = xiyl*tempz1l_att + etayl*tempz2l_att + gammayl*tempz3l_att;
    duzdzl_att = xizl*tempz1l_att + etazl*tempz2l_att + gammazl*tempz3l_att;

    // precompute some sums to save CPU time
    duxdyl_plus_duydxl_att = duxdyl_att + duydxl_att;
    duzdxl_plus_duxdzl_att = duzdxl_att + duxdzl_att;
    duzdyl_plus_duydzl_att = duzdyl_att + duydzl_att;

    // attenuation
    // computes deviatoric strain attenuation and/or for kernel calculations
    templ = 0.33333333333333333333f * (duxdxl_att + duydyl_att + duzdzl_att); // 1./3. = 0.33333
    // local storage: stresses at this current time step
    epsilondev_xx_loc = duxdxl_att - templ;
    epsilondev_yy_loc = duydyl_att - templ;
    epsilondev_xy_loc = 0.5f * duxdyl_plus_duydxl_att;
    epsilondev_xz_loc = 0.5f * duzdxl_plus_duxdzl_att;
    epsilondev_yz_loc = 0.5f * duzdyl_plus_duydzl_att;

    if (SIMULATION_TYPE == 3) {
      epsilon_trace_over_3[tx + working_element*NGLL3] = templ;
    }

    // full anisotropic case, stress calculations
    if (ANISOTROPY){
      c11 = d_c11store[offset];
      c12 = d_c12store[offset];
      c13 = d_c13store[offset];
      c14 = d_c14store[offset];
      c15 = d_c15store[offset];
      c16 = d_c16store[offset];
      c22 = d_c22store[offset];
      c23 = d_c23store[offset];
      c24 = d_c24store[offset];
      c25 = d_c25store[offset];
      c26 = d_c26store[offset];
      c33 = d_c33store[offset];
      c34 = d_c34store[offset];
      c35 = d_c35store[offset];
      c36 = d_c36store[offset];
      c44 = d_c44store[offset];
      c45 = d_c45store[offset];
      c46 = d_c46store[offset];
      c55 = d_c55store[offset];
      c56 = d_c56store[offset];
      c66 = d_c66store[offset];

      sigma_xx = c11*duxdxl + c16*duxdyl_plus_duydxl + c12*duydyl +
                 c15*duzdxl_plus_duxdzl + c14*duzdyl_plus_duydzl + c13*duzdzl;
      sigma_yy = c12*duxdxl + c26*duxdyl_plus_duydxl + c22*duydyl +
                 c25*duzdxl_plus_duxdzl + c24*duzdyl_plus_duydzl + c23*duzdzl;
      sigma_zz = c13*duxdxl + c36*duxdyl_plus_duydxl + c23*duydyl +
                 c35*duzdxl_plus_duxdzl + c34*duzdyl_plus_duydzl + c33*duzdzl;
      sigma_xy = c16*duxdxl + c66*duxdyl_plus_duydxl + c26*duydyl +
                 c56*duzdxl_plus_duxdzl + c46*duzdyl_plus_duydzl + c36*duzdzl;
      sigma_xz = c15*duxdxl + c56*duxdyl_plus_duydxl + c25*duydyl +
                 c55*duzdxl_plus_duxdzl + c45*duzdyl_plus_duydzl + c35*duzdzl;
      sigma_yz = c14*duxdxl + c46*duxdyl_plus_duydxl + c24*duydyl +
                 c45*duzdxl_plus_duxdzl + c44*duzdyl_plus_duydzl + c34*duzdzl;
    }else{

      // isotropic case

      // compute elements with an elastic isotropic rheology
      kappal = get_global_cr( &d_kappav[offset] );
      mul = get_global_cr( &d_muv[offset] );

      // attenuation
      // use unrelaxed parameters if attenuation
      mul  = mul * get_global_cr( &one_minus_sum_beta[tx+working_element*NGLL3] ); // (i,j,k,ispec)

      lambdalplus2mul = kappal + 1.33333333333333333333f * mul;  // 4./3. = 1.3333333
      lambdal = lambdalplus2mul - 2.0f * mul;

      // compute the six components of the stress tensor sigma
      sigma_xx = lambdalplus2mul*duxdxl + lambdal*duydyl_plus_duzdzl;
      sigma_yy = lambdalplus2mul*duydyl + lambdal*duxdxl_plus_duzdzl;
      sigma_zz = lambdalplus2mul*duzdzl + lambdal*duxdxl_plus_duydyl;

      sigma_xy = mul*duxdyl_plus_duydxl;
      sigma_xz = mul*duzdxl_plus_duxdzl;
      sigma_yz = mul*duzdyl_plus_duydzl;
    }

    // attenuation
    // subtracts memory variables if attenuation
    compute_element_att_stress(tx,working_element,NSPEC,
                               R_xx,R_yy,R_xy,R_xz,R_yz,
                               &sigma_xx,&sigma_yy,&sigma_zz,&sigma_xy,&sigma_xz,&sigma_yz);

    // define symmetric components (needed for non-symmetric dot product and sigma for gravity)
    sigma_yx = sigma_xy;
    sigma_zx = sigma_xz;
    sigma_zy = sigma_yz;

    if (gravity ){
      //  computes non-symmetric terms for gravity
      compute_element_gravity(tx,working_element,&iglob,d_minus_g,d_minus_deriv_gravity,
                              d_rhostore,wgll_cube,jacobianl,
                              s_dummyx_loc,s_dummyy_loc,s_dummyz_loc,
                              &sigma_xx,&sigma_yy,&sigma_xz,&sigma_yz,
                              &rho_s_H1,&rho_s_H2,&rho_s_H3);
    }
  } // active

  //note: due to re-assignement of s_dummyx_loc_att,..,we need to sync before updating sh_tempx1...
  __syncthreads();

  if (active ){
    // form dot product with test vector, non-symmetric form
    sh_tempx1[tx] = jacobianl * (sigma_xx*xixl + sigma_yx*xiyl + sigma_zx*xizl);
    sh_tempy1[tx] = jacobianl * (sigma_xy*xixl + sigma_yy*xiyl + sigma_zy*xizl);
    sh_tempz1[tx] = jacobianl * (sigma_xz*xixl + sigma_yz*xiyl + sigma_zz*xizl);

    sh_tempx2[tx] = jacobianl * (sigma_xx*etaxl + sigma_yx*etayl + sigma_zx*etazl);
    sh_tempy2[tx] = jacobianl * (sigma_xy*etaxl + sigma_yy*etayl + sigma_zy*etazl);
    sh_tempz2[tx] = jacobianl * (sigma_xz*etaxl + sigma_yz*etayl + sigma_zz*etazl);

    sh_tempx3[tx] = jacobianl * (sigma_xx*gammaxl + sigma_yx*gammayl + sigma_zx*gammazl);
    sh_tempy3[tx] = jacobianl * (sigma_xy*gammaxl + sigma_yy*gammayl + sigma_zy*gammazl);
    sh_tempz3[tx] = jacobianl * (sigma_xz*gammaxl + sigma_yz*gammayl + sigma_zz*gammazl);
  }

  // re-assigns sh_hprime_xx to load hprimewgll
  // note: the sync seems to be necessary, otherwise there is more jitter, not sure why...
  __syncthreads();
  if (tx < NGLL2) {
    // copy hprime from global memory to shared memory
    load_shared_memory_hprimewgll(&tx,d_hprimewgll_xx,sh_hprime_xx);
  }

  // synchronize all the threads (one thread for each of the NGLL grid points of the
  // current spectral element) because we need the whole element to be ready in order
  // to be able to compute the matrix products along cut planes of the 3D element below
  __syncthreads();

// JC JC here we will need to add GPU support for the new C-PML routines

  if (active ){

#ifndef MANUALLY_UNROLLED_LOOPS

    tempx1l = 0.f;
    tempy1l = 0.f;
    tempz1l = 0.f;

    tempx2l = 0.f;
    tempy2l = 0.f;
    tempz2l = 0.f;

    tempx3l = 0.f;
    tempy3l = 0.f;
    tempz3l = 0.f;

    for (l=0;l<NGLLX;l++) {
      // assumes hprimewgll_xx == hprimewgll_yy == hprimewgll_zz
      fac1 = sh_hprime_xx[I*NGLLX+l]; //  d_hprimewgll_xx[I*NGLLX+l];
      tempx1l += sh_tempx1[K*NGLL2+J*NGLLX+l]*fac1;
      tempy1l += sh_tempy1[K*NGLL2+J*NGLLX+l]*fac1;
      tempz1l += sh_tempz1[K*NGLL2+J*NGLLX+l]*fac1;

      fac2 = sh_hprime_xx[J*NGLLX+l]; // d_hprimewgll_xx[J*NGLLX+l];
      tempx2l += sh_tempx2[K*NGLL2+l*NGLLX+I]*fac2;
      tempy2l += sh_tempy2[K*NGLL2+l*NGLLX+I]*fac2;
      tempz2l += sh_tempz2[K*NGLL2+l*NGLLX+I]*fac2;

      fac3 = sh_hprime_xx[K*NGLLX+l]; // d_hprimewgll_xx[K*NGLLX+l];
      tempx3l += sh_tempx3[l*NGLL2+J*NGLLX+I]*fac3;
      tempy3l += sh_tempy3[l*NGLL2+J*NGLLX+I]*fac3;
      tempz3l += sh_tempz3[l*NGLL2+J*NGLLX+I]*fac3;
    }
#else
    tempx1l = sh_tempx1[K*NGLL2+J*NGLLX]*d_hprimewgll_xx[I*NGLLX]
            + sh_tempx1[K*NGLL2+J*NGLLX+1]*d_hprimewgll_xx[I*NGLLX+1]
            + sh_tempx1[K*NGLL2+J*NGLLX+2]*d_hprimewgll_xx[I*NGLLX+2]
            + sh_tempx1[K*NGLL2+J*NGLLX+3]*d_hprimewgll_xx[I*NGLLX+3]
            + sh_tempx1[K*NGLL2+J*NGLLX+4]*d_hprimewgll_xx[I*NGLLX+4];

    tempy1l = sh_tempy1[K*NGLL2+J*NGLLX]*d_hprimewgll_xx[I*NGLLX]
            + sh_tempy1[K*NGLL2+J*NGLLX+1]*d_hprimewgll_xx[I*NGLLX+1]
            + sh_tempy1[K*NGLL2+J*NGLLX+2]*d_hprimewgll_xx[I*NGLLX+2]
            + sh_tempy1[K*NGLL2+J*NGLLX+3]*d_hprimewgll_xx[I*NGLLX+3]
            + sh_tempy1[K*NGLL2+J*NGLLX+4]*d_hprimewgll_xx[I*NGLLX+4];

    tempz1l = sh_tempz1[K*NGLL2+J*NGLLX]*d_hprimewgll_xx[I*NGLLX]
            + sh_tempz1[K*NGLL2+J*NGLLX+1]*d_hprimewgll_xx[I*NGLLX+1]
            + sh_tempz1[K*NGLL2+J*NGLLX+2]*d_hprimewgll_xx[I*NGLLX+2]
            + sh_tempz1[K*NGLL2+J*NGLLX+3]*d_hprimewgll_xx[I*NGLLX+3]
            + sh_tempz1[K*NGLL2+J*NGLLX+4]*d_hprimewgll_xx[I*NGLLX+4];

    tempx2l = sh_tempx2[K*NGLL2+I]*d_hprimewgll_xx[J*NGLLX]
            + sh_tempx2[K*NGLL2+NGLLX+I]*d_hprimewgll_xx[J*NGLLX+1]
            + sh_tempx2[K*NGLL2+2*NGLLX+I]*d_hprimewgll_xx[J*NGLLX+2]
            + sh_tempx2[K*NGLL2+3*NGLLX+I]*d_hprimewgll_xx[J*NGLLX+3]
            + sh_tempx2[K*NGLL2+4*NGLLX+I]*d_hprimewgll_xx[J*NGLLX+4];

    tempy2l = sh_tempy2[K*NGLL2+I]*d_hprimewgll_xx[J*NGLLX]
            + sh_tempy2[K*NGLL2+NGLLX+I]*d_hprimewgll_xx[J*NGLLX+1]
            + sh_tempy2[K*NGLL2+2*NGLLX+I]*d_hprimewgll_xx[J*NGLLX+2]
            + sh_tempy2[K*NGLL2+3*NGLLX+I]*d_hprimewgll_xx[J*NGLLX+3]
            + sh_tempy2[K*NGLL2+4*NGLLX+I]*d_hprimewgll_xx[J*NGLLX+4];

    tempz2l = sh_tempz2[K*NGLL2+I]*d_hprimewgll_xx[J*NGLLX]
            + sh_tempz2[K*NGLL2+NGLLX+I]*d_hprimewgll_xx[J*NGLLX+1]
            + sh_tempz2[K*NGLL2+2*NGLLX+I]*d_hprimewgll_xx[J*NGLLX+2]
            + sh_tempz2[K*NGLL2+3*NGLLX+I]*d_hprimewgll_xx[J*NGLLX+3]
            + sh_tempz2[K*NGLL2+4*NGLLX+I]*d_hprimewgll_xx[J*NGLLX+4];

    tempx3l = sh_tempx3[J*NGLLX+I]*d_hprimewgll_xx[K*NGLLX]
            + sh_tempx3[NGLL2+J*NGLLX+I]*d_hprimewgll_xx[K*NGLLX+1]
            + sh_tempx3[2*NGLL2+J*NGLLX+I]*d_hprimewgll_xx[K*NGLLX+2]
            + sh_tempx3[3*NGLL2+J*NGLLX+I]*d_hprimewgll_xx[K*NGLLX+3]
            + sh_tempx3[4*NGLL2+J*NGLLX+I]*d_hprimewgll_xx[K*NGLLX+4];

    tempy3l = sh_tempy3[J*NGLLX+I]*d_hprimewgll_xx[K*NGLLX]
            + sh_tempy3[NGLL2+J*NGLLX+I]*d_hprimewgll_xx[K*NGLLX+1]
            + sh_tempy3[2*NGLL2+J*NGLLX+I]*d_hprimewgll_xx[K*NGLLX+2]
            + sh_tempy3[3*NGLL2+J*NGLLX+I]*d_hprimewgll_xx[K*NGLLX+3]
            + sh_tempy3[4*NGLL2+J*NGLLX+I]*d_hprimewgll_xx[K*NGLLX+4];

    tempz3l = sh_tempz3[J*NGLLX+I]*d_hprimewgll_xx[K*NGLLX]
            + sh_tempz3[NGLL2+J*NGLLX+I]*d_hprimewgll_xx[K*NGLLX+1]
            + sh_tempz3[2*NGLL2+J*NGLLX+I]*d_hprimewgll_xx[K*NGLLX+2]
            + sh_tempz3[3*NGLL2+J*NGLLX+I]*d_hprimewgll_xx[K*NGLLX+3]
            + sh_tempz3[4*NGLL2+J*NGLLX+I]*d_hprimewgll_xx[K*NGLLX+4];
#endif

    fac1 = d_wgllwgll_yz[K*NGLLX+J];
    fac2 = d_wgllwgll_xz[K*NGLLX+I];
    fac3 = d_wgllwgll_xy[J*NGLLX+I];

    sum_terms1 = - (fac1*tempx1l + fac2*tempx2l + fac3*tempx3l);
    sum_terms2 = - (fac1*tempy1l + fac2*tempy2l + fac3*tempy3l);
    sum_terms3 = - (fac1*tempz1l + fac2*tempz2l + fac3*tempz3l);

    // adds gravity term
    if (gravity ){
      sum_terms1 += rho_s_H1;
      sum_terms2 += rho_s_H2;
      sum_terms3 += rho_s_H3;
    }

#ifdef USE_MESH_COLORING_GPU
    // no atomic operation needed, colors don't share global points between elements

#ifdef USE_TEXTURES_FIELDS
    d_accel[iglob*3]     = texfetch_accel<FORWARD_OR_ADJOINT>(iglob*3) + sum_terms1;
    d_accel[iglob*3 + 1] = texfetch_accel<FORWARD_OR_ADJOINT>(iglob*3 + 1) + sum_terms2;
    d_accel[iglob*3 + 2] = texfetch_accel<FORWARD_OR_ADJOINT>(iglob*3 + 2) + sum_terms3;
#else
    d_accel[iglob*3]     += sum_terms1;
    d_accel[iglob*3 + 1] += sum_terms2;
    d_accel[iglob*3 + 2] += sum_terms3;
#endif // USE_TEXTURES_FIELDS

// JC JC here we will need to add GPU support for the new C-PML routines

#else // MESH_COLORING

    //mesh coloring
    if (use_mesh_coloring_gpu ){

      // no atomic operation needed, colors don't share global points between elements
#ifdef USE_TEXTURES_FIELDS
      d_accel[iglob*3]     = texfetch_accel<FORWARD_OR_ADJOINT>(iglob*3) + sum_terms1;
      d_accel[iglob*3 + 1] = texfetch_accel<FORWARD_OR_ADJOINT>(iglob*3 + 1) + sum_terms2;
      d_accel[iglob*3 + 2] = texfetch_accel<FORWARD_OR_ADJOINT>(iglob*3 + 2) + sum_terms3;
#else
      d_accel[iglob*3]     += sum_terms1;
      d_accel[iglob*3 + 1] += sum_terms2;
      d_accel[iglob*3 + 2] += sum_terms3;
#endif // USE_TEXTURES_FIELDS

    }
    else {
      atomicAdd(&d_accel[iglob*3], sum_terms1);
      atomicAdd(&d_accel[iglob*3+1], sum_terms2);
      atomicAdd(&d_accel[iglob*3+2], sum_terms3);
    } // if (use_mesh_coloring_gpu)

#endif // MESH_COLORING

    // attenuation
    // update memory variables based upon the Runge-Kutta scheme
    compute_element_att_memory(tx,working_element,NSPEC,
                               d_muv,
                               factor_common,alphaval,betaval,gammaval,
                               R_xx,R_yy,R_xy,R_xz,R_yz,
                               epsilondev_xx,epsilondev_yy,epsilondev_xy,epsilondev_xz,epsilondev_yz,
                               epsilondev_xx_loc,epsilondev_yy_loc,epsilondev_xy_loc,epsilondev_xz_loc,epsilondev_yz_loc);

    // save deviatoric strain for Runge-Kutta scheme
    // fortran: epsilondev_xx(:,:,:,ispec) = epsilondev_xx_loc(:,:,:)
    epsilondev_xx[tx + working_element*NGLL3] = epsilondev_xx_loc;
    epsilondev_yy[tx + working_element*NGLL3] = epsilondev_yy_loc;
    epsilondev_xy[tx + working_element*NGLL3] = epsilondev_xy_loc;
    epsilondev_xz[tx + working_element*NGLL3] = epsilondev_xz_loc;
    epsilondev_yz[tx + working_element*NGLL3] = epsilondev_yz_loc;
  } // if (active)

// JC JC here we will need to add GPU support for the new C-PML routines

} // kernel_2_att_impl()


/* ----------------------------------------------------------------------------------------------- */

/*

// original kernel
// please leave it here for reference ...

template<int FORWARD_OR_ADJOINT> __global__ void
#ifdef USE_LAUNCH_BOUNDS
// adds compiler specification
__launch_bounds__(NGLL3_PADDED,LAUNCH_MIN_BLOCKS)
#endif
// main kernel
Kernel_2_att_org_impl(int nb_blocks_to_compute,
                  const int* d_ibool,
                  const int* d_phase_ispec_inner_elastic,const int num_phase_ispec_elastic,
                  const int d_iphase,
                  const int use_mesh_coloring_gpu,
                  realw d_deltat,
                  realw_const_p d_displ,
                  realw_const_p d_veloc,
                  realw_p d_accel,
                  realw_const_p d_xix,realw_const_p d_xiy,realw_const_p d_xiz,
                  realw_const_p d_etax,realw_const_p d_etay,realw_const_p d_etaz,
                  realw_const_p d_gammax,realw_const_p d_gammay,realw_const_p d_gammaz,
                  realw_const_p d_hprime_xx,
                  realw_const_p d_hprimewgll_xx,
                  realw_const_p d_wgllwgll_xy,realw_const_p d_wgllwgll_xz,realw_const_p d_wgllwgll_yz,
                  realw_const_p d_kappav,realw_const_p d_muv,
                  realw_p epsilondev_xx,realw_p epsilondev_yy,realw_p epsilondev_xy,
                  realw_p epsilondev_xz,realw_p epsilondev_yz,
                  realw_p epsilon_trace_over_3,
                  const int SIMULATION_TYPE,
                  const int NSPEC,
                  realw_const_p one_minus_sum_beta,realw_const_p factor_common,
                  realw_p R_xx,realw_p R_yy,realw_p R_xy,realw_p R_xz,realw_p R_yz,
                  realw_const_p alphaval,realw_const_p betaval,realw_const_p gammaval,
                  const int ANISOTROPY,
                  realw_const_p d_c11store,realw_const_p d_c12store,realw_const_p d_c13store,
                  realw_const_p d_c14store,realw_const_p d_c15store,realw_const_p d_c16store,
                  realw_const_p d_c22store,realw_const_p d_c23store,realw_const_p d_c24store,
                  realw_const_p d_c25store,realw_const_p d_c26store,realw_const_p d_c33store,
                  realw_const_p d_c34store,realw_const_p d_c35store,realw_const_p d_c36store,
                  realw_const_p d_c44store,realw_const_p d_c45store,realw_const_p d_c46store,
                  realw_const_p d_c55store,realw_const_p d_c56store,realw_const_p d_c66store,
                  const int gravity,
                  realw_const_p d_minus_g,
                  realw_const_p d_minus_deriv_gravity,
                  realw_const_p d_rhostore,
                  realw_const_p wgll_cube){


// elastic compute kernel with attenuation
// holds for: ATTENUATION = .true.
//            COMPUTE_AND_STORE_STRAIN = .true. (always true for attenuation)

  int bx = blockIdx.y*gridDim.x+blockIdx.x;
  int tx = threadIdx.x;

  int K = (tx/NGLL2);
  int J = ((tx-K*NGLL2)/NGLLX);
  int I = (tx-K*NGLL2-J*NGLLX);

  unsigned short int active;
  int iglob,offset;
  int working_element;

  realw tempx1l,tempx2l,tempx3l,tempy1l,tempy2l,tempy3l,tempz1l,tempz2l,tempz3l;
  realw xixl,xiyl,xizl,etaxl,etayl,etazl,gammaxl,gammayl,gammazl,jacobianl;
  realw duxdxl,duxdyl,duxdzl,duydxl,duydyl,duydzl,duzdxl,duzdyl,duzdzl;
  realw duxdxl_plus_duydyl,duxdxl_plus_duzdzl,duydyl_plus_duzdzl;
  realw duxdyl_plus_duydxl,duzdxl_plus_duxdzl,duzdyl_plus_duydzl;
  realw templ;

  realw tempx1l_att,tempx2l_att,tempx3l_att,tempy1l_att,tempy2l_att,tempy3l_att,tempz1l_att,tempz2l_att,tempz3l_att;
  realw duxdxl_att,duxdyl_att,duxdzl_att,duydxl_att,duydyl_att,duydzl_att,duzdxl_att,duzdyl_att,duzdzl_att;
  realw duxdyl_plus_duydxl_att,duzdxl_plus_duxdzl_att,duzdyl_plus_duydzl_att;

  realw fac1,fac2,fac3;
  realw lambdal,mul,lambdalplus2mul,kappal;

  realw sigma_xx,sigma_yy,sigma_zz,sigma_xy,sigma_xz,sigma_yz;
  realw epsilondev_xx_loc,epsilondev_yy_loc,epsilondev_xy_loc,epsilondev_xz_loc,epsilondev_yz_loc;

  realw c11,c12,c13,c14,c15,c16,c22,c23,c24,c25,c26,c33,c34,c35,c36,c44,c45,c46,c55,c56,c66;
  realw sum_terms1,sum_terms2,sum_terms3;

  // gravity variables
  realw sigma_yx,sigma_zx,sigma_zy;
  realw rho_s_H1,rho_s_H2,rho_s_H3;

#ifndef MANUALLY_UNROLLED_LOOPS
  int l;
#endif

  __shared__ realw sh_tempx1[NGLL3];
  __shared__ realw sh_tempx2[NGLL3];
  __shared__ realw sh_tempx3[NGLL3];

  __shared__ realw sh_tempy1[NGLL3];
  __shared__ realw sh_tempy2[NGLL3];
  __shared__ realw sh_tempy3[NGLL3];

  __shared__ realw sh_tempz1[NGLL3];
  __shared__ realw sh_tempz2[NGLL3];
  __shared__ realw sh_tempz3[NGLL3];

  __shared__ realw sh_hprime_xx[NGLL2];

  __shared__ realw s_dummyx_loc[NGLL3];
  __shared__ realw s_dummyy_loc[NGLL3];
  __shared__ realw s_dummyz_loc[NGLL3];

  // re-assigns shared array to decrease shared memory usage
  // note: this will re-use s_temp arrays from above to save shared memory
  realw* s_dummyx_loc_att = (realw*) sh_tempx1;
  realw* s_dummyy_loc_att = (realw*) sh_tempx2;
  realw* s_dummyz_loc_att = (realw*) sh_tempx3;

  // use only NGLL^3 = 125 active threads, plus 3 inactive/ghost threads,
  // because we used memory padding from NGLL^3 = 125 to 128 to get coalescent memory accesses
  active = (tx < NGLL3 && bx < nb_blocks_to_compute) ? 1:0;

  // copy from global memory to shared memory
  // each thread writes one of the NGLL^3 = 125 data points
  if (active ){

#ifdef USE_MESH_COLORING_GPU
    working_element = bx;
#else
    //mesh coloring
    if (use_mesh_coloring_gpu ){
      working_element = bx;
    }else{
      // iphase-1 and working_element-1 for Fortran->C array conventions
      working_element = d_phase_ispec_inner_elastic[bx + num_phase_ispec_elastic*(d_iphase-1)]-1;
    }
#endif
    // local padded index
    offset = working_element*NGLL3_PADDED + tx;

    // global index
    iglob = d_ibool[offset] - 1;

    // copy displacement from global memory to shared memory
    load_shared_memory_displ<FORWARD_OR_ADJOINT>(&tx,&iglob,d_displ,s_dummyx_loc,s_dummyy_loc,s_dummyz_loc);

// JC JC here we will need to add GPU support for the new C-PML routines

    // attenuation
    // use first order Taylor expansion of displacement for local storage of stresses
    // at this current time step, to fix attenuation in a consistent way
#ifdef USE_TEXTURES_FIELDS
    s_dummyx_loc_att[tx] = s_dummyx_loc[tx] + d_deltat * texfetch_veloc<FORWARD_OR_ADJOINT>(iglob*3);
    s_dummyy_loc_att[tx] = s_dummyy_loc[tx] + d_deltat * texfetch_veloc<FORWARD_OR_ADJOINT>(iglob*3 + 1);
    s_dummyz_loc_att[tx] = s_dummyz_loc[tx] + d_deltat * texfetch_veloc<FORWARD_OR_ADJOINT>(iglob*3 + 2);
#else
    s_dummyx_loc_att[tx] = s_dummyx_loc[tx] + d_deltat * d_veloc[iglob*3];
    s_dummyy_loc_att[tx] = s_dummyy_loc[tx] + d_deltat * d_veloc[iglob*3 + 1];
    s_dummyz_loc_att[tx] = s_dummyz_loc[tx] + d_deltat * d_veloc[iglob*3 + 2];
#endif
  }// active

  // loads hprime's into shared memory
  if (tx < NGLL2) {
    // copy hprime from global memory to shared memory
    load_shared_memory_hprime(&tx,d_hprime_xx,sh_hprime_xx);
  }

  // synchronize all the threads (one thread for each of the NGLL grid points of the
  // current spectral element) because we need the whole element to be ready in order
  // to be able to compute the matrix products along cut planes of the 3D element below
  __syncthreads();

  if (active ){

#ifndef MANUALLY_UNROLLED_LOOPS

    tempx1l = 0.f;
    tempx2l = 0.f;
    tempx3l = 0.f;

    tempy1l = 0.f;
    tempy2l = 0.f;
    tempy3l = 0.f;

    tempz1l = 0.f;
    tempz2l = 0.f;
    tempz3l = 0.f;

    for (l=0;l<NGLLX;l++) {
      //assumes that hprime_xx = hprime_yy = hprime_zz
      fac1 = sh_hprime_xx[l*NGLLX+I];
      tempx1l += s_dummyx_loc[K*NGLL2+J*NGLLX+l]*fac1;
      tempy1l += s_dummyy_loc[K*NGLL2+J*NGLLX+l]*fac1;
      tempz1l += s_dummyz_loc[K*NGLL2+J*NGLLX+l]*fac1;

      fac2 = sh_hprime_xx[l*NGLLX+J];
      tempx2l += s_dummyx_loc[K*NGLL2+l*NGLLX+I]*fac2;
      tempy2l += s_dummyy_loc[K*NGLL2+l*NGLLX+I]*fac2;
      tempz2l += s_dummyz_loc[K*NGLL2+l*NGLLX+I]*fac2;

      fac3 = sh_hprime_xx[l*NGLLX+K];
      tempx3l += s_dummyx_loc[l*NGLL2+J*NGLLX+I]*fac3;
      tempy3l += s_dummyy_loc[l*NGLL2+J*NGLLX+I]*fac3;
      tempz3l += s_dummyz_loc[l*NGLL2+J*NGLLX+I]*fac3;
    }

// JC JC here we will need to add GPU support for the new C-PML routines

    // attenuation
    // temporary variables used for fixing attenuation in a consistent way
    tempx1l_att = 0.f;
    tempx2l_att = 0.f;
    tempx3l_att = 0.f;

    tempy1l_att = 0.f;
    tempy2l_att = 0.f;
    tempy3l_att = 0.f;

    tempz1l_att = 0.f;
    tempz2l_att = 0.f;
    tempz3l_att = 0.f;

    for (l=0;l<NGLLX;l++) {
      fac1 = sh_hprime_xx[l*NGLLX+I];
      tempx1l_att += s_dummyx_loc_att[K*NGLL2+J*NGLLX+l]*fac1;
      tempy1l_att += s_dummyy_loc_att[K*NGLL2+J*NGLLX+l]*fac1;
      tempz1l_att += s_dummyz_loc_att[K*NGLL2+J*NGLLX+l]*fac1;

      fac2 = sh_hprime_xx[l*NGLLX+J];
      tempx2l_att += s_dummyx_loc_att[K*NGLL2+l*NGLLX+I]*fac2;
      tempy2l_att += s_dummyy_loc_att[K*NGLL2+l*NGLLX+I]*fac2;
      tempz2l_att += s_dummyz_loc_att[K*NGLL2+l*NGLLX+I]*fac2;

      fac3 = sh_hprime_xx[l*NGLLX+K];
      tempx3l_att += s_dummyx_loc_att[l*NGLL2+J*NGLLX+I]*fac3;
      tempy3l_att += s_dummyy_loc_att[l*NGLL2+J*NGLLX+I]*fac3;
      tempz3l_att += s_dummyz_loc_att[l*NGLL2+J*NGLLX+I]*fac3;
    }

#else
    tempx1l = s_dummyx_loc[K*NGLL2+J*NGLLX]*d_hprime_xx[I]
            + s_dummyx_loc[K*NGLL2+J*NGLLX+1]*d_hprime_xx[NGLLX+I]
            + s_dummyx_loc[K*NGLL2+J*NGLLX+2]*d_hprime_xx[2*NGLLX+I]
            + s_dummyx_loc[K*NGLL2+J*NGLLX+3]*d_hprime_xx[3*NGLLX+I]
            + s_dummyx_loc[K*NGLL2+J*NGLLX+4]*d_hprime_xx[4*NGLLX+I];

    tempy1l = s_dummyy_loc[K*NGLL2+J*NGLLX]*d_hprime_xx[I]
            + s_dummyy_loc[K*NGLL2+J*NGLLX+1]*d_hprime_xx[NGLLX+I]
            + s_dummyy_loc[K*NGLL2+J*NGLLX+2]*d_hprime_xx[2*NGLLX+I]
            + s_dummyy_loc[K*NGLL2+J*NGLLX+3]*d_hprime_xx[3*NGLLX+I]
            + s_dummyy_loc[K*NGLL2+J*NGLLX+4]*d_hprime_xx[4*NGLLX+I];

    tempz1l = s_dummyz_loc[K*NGLL2+J*NGLLX]*d_hprime_xx[I]
            + s_dummyz_loc[K*NGLL2+J*NGLLX+1]*d_hprime_xx[NGLLX+I]
            + s_dummyz_loc[K*NGLL2+J*NGLLX+2]*d_hprime_xx[2*NGLLX+I]
            + s_dummyz_loc[K*NGLL2+J*NGLLX+3]*d_hprime_xx[3*NGLLX+I]
            + s_dummyz_loc[K*NGLL2+J*NGLLX+4]*d_hprime_xx[4*NGLLX+I];

    tempx2l = s_dummyx_loc[K*NGLL2+I]*d_hprime_xx[J]
            + s_dummyx_loc[K*NGLL2+NGLLX+I]*d_hprime_xx[NGLLX+J]
            + s_dummyx_loc[K*NGLL2+2*NGLLX+I]*d_hprime_xx[2*NGLLX+J]
            + s_dummyx_loc[K*NGLL2+3*NGLLX+I]*d_hprime_xx[3*NGLLX+J]
            + s_dummyx_loc[K*NGLL2+4*NGLLX+I]*d_hprime_xx[4*NGLLX+J];

    tempy2l = s_dummyy_loc[K*NGLL2+I]*d_hprime_xx[J]
            + s_dummyy_loc[K*NGLL2+NGLLX+I]*d_hprime_xx[NGLLX+J]
            + s_dummyy_loc[K*NGLL2+2*NGLLX+I]*d_hprime_xx[2*NGLLX+J]
            + s_dummyy_loc[K*NGLL2+3*NGLLX+I]*d_hprime_xx[3*NGLLX+J]
            + s_dummyy_loc[K*NGLL2+4*NGLLX+I]*d_hprime_xx[4*NGLLX+J];

    tempz2l = s_dummyz_loc[K*NGLL2+I]*d_hprime_xx[J]
            + s_dummyz_loc[K*NGLL2+NGLLX+I]*d_hprime_xx[NGLLX+J]
            + s_dummyz_loc[K*NGLL2+2*NGLLX+I]*d_hprime_xx[2*NGLLX+J]
            + s_dummyz_loc[K*NGLL2+3*NGLLX+I]*d_hprime_xx[3*NGLLX+J]
            + s_dummyz_loc[K*NGLL2+4*NGLLX+I]*d_hprime_xx[4*NGLLX+J];

    tempx3l = s_dummyx_loc[J*NGLLX+I]*d_hprime_xx[K]
            + s_dummyx_loc[NGLL2+J*NGLLX+I]*d_hprime_xx[NGLLX+K]
            + s_dummyx_loc[2*NGLL2+J*NGLLX+I]*d_hprime_xx[2*NGLLX+K]
            + s_dummyx_loc[3*NGLL2+J*NGLLX+I]*d_hprime_xx[3*NGLLX+K]
            + s_dummyx_loc[4*NGLL2+J*NGLLX+I]*d_hprime_xx[4*NGLLX+K];

    tempy3l = s_dummyy_loc[J*NGLLX+I]*d_hprime_xx[K]
            + s_dummyy_loc[NGLL2+J*NGLLX+I]*d_hprime_xx[NGLLX+K]
            + s_dummyy_loc[2*NGLL2+J*NGLLX+I]*d_hprime_xx[2*NGLLX+K]
            + s_dummyy_loc[3*NGLL2+J*NGLLX+I]*d_hprime_xx[3*NGLLX+K]
            + s_dummyy_loc[4*NGLL2+J*NGLLX+I]*d_hprime_xx[4*NGLLX+K];

    tempz3l = s_dummyz_loc[J*NGLLX+I]*d_hprime_xx[K]
            + s_dummyz_loc[NGLL2+J*NGLLX+I]*d_hprime_xx[NGLLX+K]
            + s_dummyz_loc[2*NGLL2+J*NGLLX+I]*d_hprime_xx[2*NGLLX+K]
            + s_dummyz_loc[3*NGLL2+J*NGLLX+I]*d_hprime_xx[3*NGLLX+K]
            + s_dummyz_loc[4*NGLL2+J*NGLLX+I]*d_hprime_xx[4*NGLLX+K];

// JC JC here we will need to add GPU support for the new C-PML routines

    // attenuation
    // temporary variables used for fixing attenuation in a consistent way
    tempx1l_att = s_dummyx_loc_att[K*NGLL2+J*NGLLX]*d_hprime_xx[I]
                  + s_dummyx_loc_att[K*NGLL2+J*NGLLX+1]*d_hprime_xx[NGLLX+I]
                  + s_dummyx_loc_att[K*NGLL2+J*NGLLX+2]*d_hprime_xx[2*NGLLX+I]
                  + s_dummyx_loc_att[K*NGLL2+J*NGLLX+3]*d_hprime_xx[3*NGLLX+I]
                  + s_dummyx_loc_att[K*NGLL2+J*NGLLX+4]*d_hprime_xx[4*NGLLX+I];

    tempy1l_att = s_dummyy_loc_att[K*NGLL2+J*NGLLX]*d_hprime_xx[I]
                  + s_dummyy_loc_att[K*NGLL2+J*NGLLX+1]*d_hprime_xx[NGLLX+I]
                  + s_dummyy_loc_att[K*NGLL2+J*NGLLX+2]*d_hprime_xx[2*NGLLX+I]
                  + s_dummyy_loc_att[K*NGLL2+J*NGLLX+3]*d_hprime_xx[3*NGLLX+I]
                  + s_dummyy_loc_att[K*NGLL2+J*NGLLX+4]*d_hprime_xx[4*NGLLX+I];

    tempz1l_att = s_dummyz_loc_att[K*NGLL2+J*NGLLX]*d_hprime_xx[I]
                  + s_dummyz_loc_att[K*NGLL2+J*NGLLX+1]*d_hprime_xx[NGLLX+I]
                  + s_dummyz_loc_att[K*NGLL2+J*NGLLX+2]*d_hprime_xx[2*NGLLX+I]
                  + s_dummyz_loc_att[K*NGLL2+J*NGLLX+3]*d_hprime_xx[3*NGLLX+I]
                  + s_dummyz_loc_att[K*NGLL2+J*NGLLX+4]*d_hprime_xx[4*NGLLX+I];

    tempx2l_att = s_dummyx_loc_att[K*NGLL2+I]*d_hprime_xx[J]
                  + s_dummyx_loc_att[K*NGLL2+NGLLX+I]*d_hprime_xx[NGLLX+J]
                  + s_dummyx_loc_att[K*NGLL2+2*NGLLX+I]*d_hprime_xx[2*NGLLX+J]
                  + s_dummyx_loc_att[K*NGLL2+3*NGLLX+I]*d_hprime_xx[3*NGLLX+J]
                  + s_dummyx_loc_att[K*NGLL2+4*NGLLX+I]*d_hprime_xx[4*NGLLX+J];

    tempy2l_att = s_dummyy_loc_att[K*NGLL2+I]*d_hprime_xx[J]
                  + s_dummyy_loc_att[K*NGLL2+NGLLX+I]*d_hprime_xx[NGLLX+J]
                  + s_dummyy_loc_att[K*NGLL2+2*NGLLX+I]*d_hprime_xx[2*NGLLX+J]
                  + s_dummyy_loc_att[K*NGLL2+3*NGLLX+I]*d_hprime_xx[3*NGLLX+J]
                  + s_dummyy_loc_att[K*NGLL2+4*NGLLX+I]*d_hprime_xx[4*NGLLX+J];

    tempz2l_att = s_dummyz_loc_att[K*NGLL2+I]*d_hprime_xx[J]
                  + s_dummyz_loc_att[K*NGLL2+NGLLX+I]*d_hprime_xx[NGLLX+J]
                  + s_dummyz_loc_att[K*NGLL2+2*NGLLX+I]*d_hprime_xx[2*NGLLX+J]
                  + s_dummyz_loc_att[K*NGLL2+3*NGLLX+I]*d_hprime_xx[3*NGLLX+J]
                  + s_dummyz_loc_att[K*NGLL2+4*NGLLX+I]*d_hprime_xx[4*NGLLX+J];

    tempx3l_att = s_dummyx_loc_att[J*NGLLX+I]*d_hprime_xx[K]
                  + s_dummyx_loc_att[NGLL2+J*NGLLX+I]*d_hprime_xx[NGLLX+K]
                  + s_dummyx_loc_att[2*NGLL2+J*NGLLX+I]*d_hprime_xx[2*NGLLX+K]
                  + s_dummyx_loc_att[3*NGLL2+J*NGLLX+I]*d_hprime_xx[3*NGLLX+K]
                  + s_dummyx_loc_att[4*NGLL2+J*NGLLX+I]*d_hprime_xx[4*NGLLX+K];

    tempy3l_att = s_dummyy_loc_att[J*NGLLX+I]*d_hprime_xx[K]
                  + s_dummyy_loc_att[NGLL2+J*NGLLX+I]*d_hprime_xx[NGLLX+K]
                  + s_dummyy_loc_att[2*NGLL2+J*NGLLX+I]*d_hprime_xx[2*NGLLX+K]
                  + s_dummyy_loc_att[3*NGLL2+J*NGLLX+I]*d_hprime_xx[3*NGLLX+K]
                  + s_dummyy_loc_att[4*NGLL2+J*NGLLX+I]*d_hprime_xx[4*NGLLX+K];

    tempz3l_att = s_dummyz_loc_att[J*NGLLX+I]*d_hprime_xx[K]
                  + s_dummyz_loc_att[NGLL2+J*NGLLX+I]*d_hprime_xx[NGLLX+K]
                  + s_dummyz_loc_att[2*NGLL2+J*NGLLX+I]*d_hprime_xx[2*NGLLX+K]
                  + s_dummyz_loc_att[3*NGLL2+J*NGLLX+I]*d_hprime_xx[3*NGLLX+K]
                  + s_dummyz_loc_att[4*NGLL2+J*NGLLX+I]*d_hprime_xx[4*NGLLX+K];
#endif

    // compute derivatives of ux, uy and uz with respect to x, y and z
    xixl = get_global_cr( &d_xix[offset] );
    xiyl = get_global_cr( &d_xiy[offset] );
    xizl = get_global_cr( &d_xiz[offset] );
    etaxl = get_global_cr( &d_etax[offset] );
    etayl = get_global_cr( &d_etay[offset] );
    etazl = get_global_cr( &d_etaz[offset] );
    gammaxl = get_global_cr( &d_gammax[offset] );
    gammayl = get_global_cr( &d_gammay[offset] );
    gammazl = get_global_cr( &d_gammaz[offset] );

    jacobianl = 1.0f / (xixl*(etayl*gammazl-etazl*gammayl)
                       -xiyl*(etaxl*gammazl-etazl*gammaxl)
                       +xizl*(etaxl*gammayl-etayl*gammaxl));

    duxdxl = xixl*tempx1l + etaxl*tempx2l + gammaxl*tempx3l;
    duxdyl = xiyl*tempx1l + etayl*tempx2l + gammayl*tempx3l;
    duxdzl = xizl*tempx1l + etazl*tempx2l + gammazl*tempx3l;

    duydxl = xixl*tempy1l + etaxl*tempy2l + gammaxl*tempy3l;
    duydyl = xiyl*tempy1l + etayl*tempy2l + gammayl*tempy3l;
    duydzl = xizl*tempy1l + etazl*tempy2l + gammazl*tempy3l;

    duzdxl = xixl*tempz1l + etaxl*tempz2l + gammaxl*tempz3l;
    duzdyl = xiyl*tempz1l + etayl*tempz2l + gammayl*tempz3l;
    duzdzl = xizl*tempz1l + etazl*tempz2l + gammazl*tempz3l;

// JC JC here we will need to add GPU support for the new C-PML routines

    // precompute some sums to save CPU time
    duxdxl_plus_duydyl = duxdxl + duydyl;
    duxdxl_plus_duzdzl = duxdxl + duzdzl;
    duydyl_plus_duzdzl = duydyl + duzdzl;
    duxdyl_plus_duydxl = duxdyl + duydxl;
    duzdxl_plus_duxdzl = duzdxl + duxdzl;
    duzdyl_plus_duydzl = duzdyl + duydzl;

// JC JC here we will need to add GPU support for the new C-PML routines

    // attenuation
    // temporary variables used for fixing attenuation in a consistent way
    duxdxl_att = xixl*tempx1l_att + etaxl*tempx2l_att + gammaxl*tempx3l_att;
    duxdyl_att = xiyl*tempx1l_att + etayl*tempx2l_att + gammayl*tempx3l_att;
    duxdzl_att = xizl*tempx1l_att + etazl*tempx2l_att + gammazl*tempx3l_att;

    duydxl_att = xixl*tempy1l_att + etaxl*tempy2l_att + gammaxl*tempy3l_att;
    duydyl_att = xiyl*tempy1l_att + etayl*tempy2l_att + gammayl*tempy3l_att;
    duydzl_att = xizl*tempy1l_att + etazl*tempy2l_att + gammazl*tempy3l_att;

    duzdxl_att = xixl*tempz1l_att + etaxl*tempz2l_att + gammaxl*tempz3l_att;
    duzdyl_att = xiyl*tempz1l_att + etayl*tempz2l_att + gammayl*tempz3l_att;
    duzdzl_att = xizl*tempz1l_att + etazl*tempz2l_att + gammazl*tempz3l_att;

    // precompute some sums to save CPU time
    duxdyl_plus_duydxl_att = duxdyl_att + duydxl_att;
    duzdxl_plus_duxdzl_att = duzdxl_att + duxdzl_att;
    duzdyl_plus_duydzl_att = duzdyl_att + duydzl_att;

    // attenuation
    // computes deviatoric strain attenuation and/or for kernel calculations
    templ = 0.33333333333333333333f * (duxdxl_att + duydyl_att + duzdzl_att); // 1./3. = 0.33333
    // local storage: stresses at this current time step
    epsilondev_xx_loc = duxdxl_att - templ;
    epsilondev_yy_loc = duydyl_att - templ;
    epsilondev_xy_loc = 0.5f * duxdyl_plus_duydxl_att;
    epsilondev_xz_loc = 0.5f * duzdxl_plus_duxdzl_att;
    epsilondev_yz_loc = 0.5f * duzdyl_plus_duydzl_att;

    if (SIMULATION_TYPE == 3) {
      epsilon_trace_over_3[tx + working_element*NGLL3] = templ;
    }

    // full anisotropic case, stress calculations
    if (ANISOTROPY){
      c11 = d_c11store[offset];
      c12 = d_c12store[offset];
      c13 = d_c13store[offset];
      c14 = d_c14store[offset];
      c15 = d_c15store[offset];
      c16 = d_c16store[offset];
      c22 = d_c22store[offset];
      c23 = d_c23store[offset];
      c24 = d_c24store[offset];
      c25 = d_c25store[offset];
      c26 = d_c26store[offset];
      c33 = d_c33store[offset];
      c34 = d_c34store[offset];
      c35 = d_c35store[offset];
      c36 = d_c36store[offset];
      c44 = d_c44store[offset];
      c45 = d_c45store[offset];
      c46 = d_c46store[offset];
      c55 = d_c55store[offset];
      c56 = d_c56store[offset];
      c66 = d_c66store[offset];

      sigma_xx = c11*duxdxl + c16*duxdyl_plus_duydxl + c12*duydyl +
                 c15*duzdxl_plus_duxdzl + c14*duzdyl_plus_duydzl + c13*duzdzl;
      sigma_yy = c12*duxdxl + c26*duxdyl_plus_duydxl + c22*duydyl +
                 c25*duzdxl_plus_duxdzl + c24*duzdyl_plus_duydzl + c23*duzdzl;
      sigma_zz = c13*duxdxl + c36*duxdyl_plus_duydxl + c23*duydyl +
                 c35*duzdxl_plus_duxdzl + c34*duzdyl_plus_duydzl + c33*duzdzl;
      sigma_xy = c16*duxdxl + c66*duxdyl_plus_duydxl + c26*duydyl +
                 c56*duzdxl_plus_duxdzl + c46*duzdyl_plus_duydzl + c36*duzdzl;
      sigma_xz = c15*duxdxl + c56*duxdyl_plus_duydxl + c25*duydyl +
                 c55*duzdxl_plus_duxdzl + c45*duzdyl_plus_duydzl + c35*duzdzl;
      sigma_yz = c14*duxdxl + c46*duxdyl_plus_duydxl + c24*duydyl +
                 c45*duzdxl_plus_duxdzl + c44*duzdyl_plus_duydzl + c34*duzdzl;
    }else{

      // isotropic case

      // compute elements with an elastic isotropic rheology
      kappal = get_global_cr( &d_kappav[offset] );
      mul = get_global_cr( &d_muv[offset] );

      // attenuation
      // use unrelaxed parameters if attenuation
      mul  = mul * get_global_cr( &one_minus_sum_beta[tx+working_element*NGLL3] ); // (i,j,k,ispec)

      lambdalplus2mul = kappal + 1.33333333333333333333f * mul;  // 4./3. = 1.3333333
      lambdal = lambdalplus2mul - 2.0f * mul;

      // compute the six components of the stress tensor sigma
      sigma_xx = lambdalplus2mul*duxdxl + lambdal*duydyl_plus_duzdzl;
      sigma_yy = lambdalplus2mul*duydyl + lambdal*duxdxl_plus_duzdzl;
      sigma_zz = lambdalplus2mul*duzdzl + lambdal*duxdxl_plus_duydyl;

      sigma_xy = mul*duxdyl_plus_duydxl;
      sigma_xz = mul*duzdxl_plus_duxdzl;
      sigma_yz = mul*duzdyl_plus_duydzl;
    }

    // attenuation
    // subtracts memory variables if attenuation
    compute_element_att_stress(tx,working_element,NSPEC,
                               R_xx,R_yy,R_xy,R_xz,R_yz,
                               &sigma_xx,&sigma_yy,&sigma_zz,&sigma_xy,&sigma_xz,&sigma_yz);

    // define symmetric components (needed for non-symmetric dot product and sigma for gravity)
    sigma_yx = sigma_xy;
    sigma_zx = sigma_xz;
    sigma_zy = sigma_yz;

    if (gravity ){
      //  computes non-symmetric terms for gravity
      compute_element_gravity(tx,working_element,&iglob,d_minus_g,d_minus_deriv_gravity,
                              d_rhostore,wgll_cube,jacobianl,
                              s_dummyx_loc,s_dummyy_loc,s_dummyz_loc,
                              &sigma_xx,&sigma_yy,&sigma_xz,&sigma_yz,
                              &rho_s_H1,&rho_s_H2,&rho_s_H3);
    }
  } // active

  //note: due to re-assignement of s_dummyx_loc_att,..,we need to sync before updating sh_tempx1...
  __syncthreads();

  if (active ){
    // form dot product with test vector, non-symmetric form
    sh_tempx1[tx] = jacobianl * (sigma_xx*xixl + sigma_yx*xiyl + sigma_zx*xizl);
    sh_tempy1[tx] = jacobianl * (sigma_xy*xixl + sigma_yy*xiyl + sigma_zy*xizl);
    sh_tempz1[tx] = jacobianl * (sigma_xz*xixl + sigma_yz*xiyl + sigma_zz*xizl);

    sh_tempx2[tx] = jacobianl * (sigma_xx*etaxl + sigma_yx*etayl + sigma_zx*etazl);
    sh_tempy2[tx] = jacobianl * (sigma_xy*etaxl + sigma_yy*etayl + sigma_zy*etazl);
    sh_tempz2[tx] = jacobianl * (sigma_xz*etaxl + sigma_yz*etayl + sigma_zz*etazl);

    sh_tempx3[tx] = jacobianl * (sigma_xx*gammaxl + sigma_yx*gammayl + sigma_zx*gammazl);
    sh_tempy3[tx] = jacobianl * (sigma_xy*gammaxl + sigma_yy*gammayl + sigma_zy*gammazl);
    sh_tempz3[tx] = jacobianl * (sigma_xz*gammaxl + sigma_yz*gammayl + sigma_zz*gammazl);
  }

  // re-assigns sh_hprime_xx to load hprimewgll
  // note: the sync seems to be necessary, otherwise there is more jitter, not sure why...
  __syncthreads();
  if (tx < NGLL2) {
    // copy hprime from global memory to shared memory
    load_shared_memory_hprimewgll(&tx,d_hprimewgll_xx,sh_hprime_xx);
  }

  // synchronize all the threads (one thread for each of the NGLL grid points of the
  // current spectral element) because we need the whole element to be ready in order
  // to be able to compute the matrix products along cut planes of the 3D element below
  __syncthreads();

// JC JC here we will need to add GPU support for the new C-PML routines

  if (active ){

#ifndef MANUALLY_UNROLLED_LOOPS

    tempx1l = 0.f;
    tempy1l = 0.f;
    tempz1l = 0.f;

    tempx2l = 0.f;
    tempy2l = 0.f;
    tempz2l = 0.f;

    tempx3l = 0.f;
    tempy3l = 0.f;
    tempz3l = 0.f;

    for (l=0;l<NGLLX;l++) {
      // assumes hprimewgll_xx == hprimewgll_yy == hprimewgll_zz
      fac1 = sh_hprime_xx[I*NGLLX+l]; //  d_hprimewgll_xx[I*NGLLX+l];
      tempx1l += sh_tempx1[K*NGLL2+J*NGLLX+l]*fac1;
      tempy1l += sh_tempy1[K*NGLL2+J*NGLLX+l]*fac1;
      tempz1l += sh_tempz1[K*NGLL2+J*NGLLX+l]*fac1;

      fac2 = sh_hprime_xx[J*NGLLX+l]; // d_hprimewgll_xx[J*NGLLX+l];
      tempx2l += sh_tempx2[K*NGLL2+l*NGLLX+I]*fac2;
      tempy2l += sh_tempy2[K*NGLL2+l*NGLLX+I]*fac2;
      tempz2l += sh_tempz2[K*NGLL2+l*NGLLX+I]*fac2;

      fac3 = sh_hprime_xx[K*NGLLX+l]; // d_hprimewgll_xx[K*NGLLX+l];
      tempx3l += sh_tempx3[l*NGLL2+J*NGLLX+I]*fac3;
      tempy3l += sh_tempy3[l*NGLL2+J*NGLLX+I]*fac3;
      tempz3l += sh_tempz3[l*NGLL2+J*NGLLX+I]*fac3;
    }
#else
    tempx1l = sh_tempx1[K*NGLL2+J*NGLLX]*d_hprimewgll_xx[I*NGLLX]
            + sh_tempx1[K*NGLL2+J*NGLLX+1]*d_hprimewgll_xx[I*NGLLX+1]
            + sh_tempx1[K*NGLL2+J*NGLLX+2]*d_hprimewgll_xx[I*NGLLX+2]
            + sh_tempx1[K*NGLL2+J*NGLLX+3]*d_hprimewgll_xx[I*NGLLX+3]
            + sh_tempx1[K*NGLL2+J*NGLLX+4]*d_hprimewgll_xx[I*NGLLX+4];

    tempy1l = sh_tempy1[K*NGLL2+J*NGLLX]*d_hprimewgll_xx[I*NGLLX]
            + sh_tempy1[K*NGLL2+J*NGLLX+1]*d_hprimewgll_xx[I*NGLLX+1]
            + sh_tempy1[K*NGLL2+J*NGLLX+2]*d_hprimewgll_xx[I*NGLLX+2]
            + sh_tempy1[K*NGLL2+J*NGLLX+3]*d_hprimewgll_xx[I*NGLLX+3]
            + sh_tempy1[K*NGLL2+J*NGLLX+4]*d_hprimewgll_xx[I*NGLLX+4];

    tempz1l = sh_tempz1[K*NGLL2+J*NGLLX]*d_hprimewgll_xx[I*NGLLX]
            + sh_tempz1[K*NGLL2+J*NGLLX+1]*d_hprimewgll_xx[I*NGLLX+1]
            + sh_tempz1[K*NGLL2+J*NGLLX+2]*d_hprimewgll_xx[I*NGLLX+2]
            + sh_tempz1[K*NGLL2+J*NGLLX+3]*d_hprimewgll_xx[I*NGLLX+3]
            + sh_tempz1[K*NGLL2+J*NGLLX+4]*d_hprimewgll_xx[I*NGLLX+4];

    tempx2l = sh_tempx2[K*NGLL2+I]*d_hprimewgll_xx[J*NGLLX]
            + sh_tempx2[K*NGLL2+NGLLX+I]*d_hprimewgll_xx[J*NGLLX+1]
            + sh_tempx2[K*NGLL2+2*NGLLX+I]*d_hprimewgll_xx[J*NGLLX+2]
            + sh_tempx2[K*NGLL2+3*NGLLX+I]*d_hprimewgll_xx[J*NGLLX+3]
            + sh_tempx2[K*NGLL2+4*NGLLX+I]*d_hprimewgll_xx[J*NGLLX+4];

    tempy2l = sh_tempy2[K*NGLL2+I]*d_hprimewgll_xx[J*NGLLX]
            + sh_tempy2[K*NGLL2+NGLLX+I]*d_hprimewgll_xx[J*NGLLX+1]
            + sh_tempy2[K*NGLL2+2*NGLLX+I]*d_hprimewgll_xx[J*NGLLX+2]
            + sh_tempy2[K*NGLL2+3*NGLLX+I]*d_hprimewgll_xx[J*NGLLX+3]
            + sh_tempy2[K*NGLL2+4*NGLLX+I]*d_hprimewgll_xx[J*NGLLX+4];

    tempz2l = sh_tempz2[K*NGLL2+I]*d_hprimewgll_xx[J*NGLLX]
            + sh_tempz2[K*NGLL2+NGLLX+I]*d_hprimewgll_xx[J*NGLLX+1]
            + sh_tempz2[K*NGLL2+2*NGLLX+I]*d_hprimewgll_xx[J*NGLLX+2]
            + sh_tempz2[K*NGLL2+3*NGLLX+I]*d_hprimewgll_xx[J*NGLLX+3]
            + sh_tempz2[K*NGLL2+4*NGLLX+I]*d_hprimewgll_xx[J*NGLLX+4];

    tempx3l = sh_tempx3[J*NGLLX+I]*d_hprimewgll_xx[K*NGLLX]
            + sh_tempx3[NGLL2+J*NGLLX+I]*d_hprimewgll_xx[K*NGLLX+1]
            + sh_tempx3[2*NGLL2+J*NGLLX+I]*d_hprimewgll_xx[K*NGLLX+2]
            + sh_tempx3[3*NGLL2+J*NGLLX+I]*d_hprimewgll_xx[K*NGLLX+3]
            + sh_tempx3[4*NGLL2+J*NGLLX+I]*d_hprimewgll_xx[K*NGLLX+4];

    tempy3l = sh_tempy3[J*NGLLX+I]*d_hprimewgll_xx[K*NGLLX]
            + sh_tempy3[NGLL2+J*NGLLX+I]*d_hprimewgll_xx[K*NGLLX+1]
            + sh_tempy3[2*NGLL2+J*NGLLX+I]*d_hprimewgll_xx[K*NGLLX+2]
            + sh_tempy3[3*NGLL2+J*NGLLX+I]*d_hprimewgll_xx[K*NGLLX+3]
            + sh_tempy3[4*NGLL2+J*NGLLX+I]*d_hprimewgll_xx[K*NGLLX+4];

    tempz3l = sh_tempz3[J*NGLLX+I]*d_hprimewgll_xx[K*NGLLX]
            + sh_tempz3[NGLL2+J*NGLLX+I]*d_hprimewgll_xx[K*NGLLX+1]
            + sh_tempz3[2*NGLL2+J*NGLLX+I]*d_hprimewgll_xx[K*NGLLX+2]
            + sh_tempz3[3*NGLL2+J*NGLLX+I]*d_hprimewgll_xx[K*NGLLX+3]
            + sh_tempz3[4*NGLL2+J*NGLLX+I]*d_hprimewgll_xx[K*NGLLX+4];
#endif

    fac1 = d_wgllwgll_yz[K*NGLLX+J];
    fac2 = d_wgllwgll_xz[K*NGLLX+I];
    fac3 = d_wgllwgll_xy[J*NGLLX+I];

    sum_terms1 = - (fac1*tempx1l + fac2*tempx2l + fac3*tempx3l);
    sum_terms2 = - (fac1*tempy1l + fac2*tempy2l + fac3*tempy3l);
    sum_terms3 = - (fac1*tempz1l + fac2*tempz2l + fac3*tempz3l);

    // adds gravity term
    if (gravity ){
      sum_terms1 += rho_s_H1;
      sum_terms2 += rho_s_H2;
      sum_terms3 += rho_s_H3;
    }

#ifdef USE_MESH_COLORING_GPU
    // no atomic operation needed, colors don't share global points between elements

#ifdef USE_TEXTURES_FIELDS
    d_accel[iglob*3]     = texfetch_accel<FORWARD_OR_ADJOINT>(iglob*3) + sum_terms1;
    d_accel[iglob*3 + 1] = texfetch_accel<FORWARD_OR_ADJOINT>(iglob*3 + 1) + sum_terms2;
    d_accel[iglob*3 + 2] = texfetch_accel<FORWARD_OR_ADJOINT>(iglob*3 + 2) + sum_terms3;
#else
    d_accel[iglob*3]     += sum_terms1;
    d_accel[iglob*3 + 1] += sum_terms2;
    d_accel[iglob*3 + 2] += sum_terms3;
#endif // USE_TEXTURES_FIELDS

// JC JC here we will need to add GPU support for the new C-PML routines

#else // MESH_COLORING

    //mesh coloring
    if (use_mesh_coloring_gpu ){

      // no atomic operation needed, colors don't share global points between elements
#ifdef USE_TEXTURES_FIELDS
      d_accel[iglob*3]     = texfetch_accel<FORWARD_OR_ADJOINT>(iglob*3) + sum_terms1;
      d_accel[iglob*3 + 1] = texfetch_accel<FORWARD_OR_ADJOINT>(iglob*3 + 1) + sum_terms2;
      d_accel[iglob*3 + 2] = texfetch_accel<FORWARD_OR_ADJOINT>(iglob*3 + 2) + sum_terms3;
#else
      d_accel[iglob*3]     += sum_terms1;
      d_accel[iglob*3 + 1] += sum_terms2;
      d_accel[iglob*3 + 2] += sum_terms3;
#endif // USE_TEXTURES_FIELDS

    }
    else {
      atomicAdd(&d_accel[iglob*3], sum_terms1);
      atomicAdd(&d_accel[iglob*3+1], sum_terms2);
      atomicAdd(&d_accel[iglob*3+2], sum_terms3);
    } // if (use_mesh_coloring_gpu)

#endif // MESH_COLORING

    // attenuation
    // update memory variables based upon the Runge-Kutta scheme
    compute_element_att_memory(tx,working_element,NSPEC,
                               d_muv,
                               factor_common,alphaval,betaval,gammaval,
                               R_xx,R_yy,R_xy,R_xz,R_yz,
                               epsilondev_xx,epsilondev_yy,epsilondev_xy,epsilondev_xz,epsilondev_yz,
                               epsilondev_xx_loc,epsilondev_yy_loc,epsilondev_xy_loc,epsilondev_xz_loc,epsilondev_yz_loc);

    // save deviatoric strain for Runge-Kutta scheme
    // fortran: epsilondev_xx(:,:,:,ispec) = epsilondev_xx_loc(:,:,:)
    epsilondev_xx[tx + working_element*NGLL3] = epsilondev_xx_loc;
    epsilondev_yy[tx + working_element*NGLL3] = epsilondev_yy_loc;
    epsilondev_xy[tx + working_element*NGLL3] = epsilondev_xy_loc;
    epsilondev_xz[tx + working_element*NGLL3] = epsilondev_xz_loc;
    epsilondev_yz[tx + working_element*NGLL3] = epsilondev_yz_loc;
  } // if (active)

// JC JC here we will need to add GPU support for the new C-PML routines

} // kernel_2_att_impl()

*/

/* ----------------------------------------------------------------------------------------------- */


template<int FORWARD_OR_ADJOINT> __global__ void
#ifdef USE_LAUNCH_BOUNDS
// adds compiler specification
__launch_bounds__(NGLL3_PADDED,LAUNCH_MIN_BLOCKS)
#endif
// main kernel
Kernel_2_noatt_iso_kelvinvoigt_impl(const int nb_blocks_to_compute,
                        const int* d_ibool,
                        const int* d_phase_ispec_inner_elastic,const int num_phase_ispec_elastic,
                        const int d_iphase,
                        realw* d_kelvin_voigt_eta,
                        realw_const_p d_displ,
                        realw_const_p d_veloc,
                        realw_p d_accel,
                        realw* d_xix,realw* d_xiy,realw* d_xiz,
                        realw* d_etax,realw* d_etay,realw* d_etaz,
                        realw* d_gammax,realw* d_gammay,realw* d_gammaz,
                        realw_const_p d_hprime_xx,
                        realw_const_p d_hprimewgll_xx,
                        realw_const_p d_wgllwgll_xy,realw_const_p d_wgllwgll_xz,realw_const_p d_wgllwgll_yz,
                        realw* d_kappav,realw* d_muv){

// elastic compute kernel without attenuation for isotropic elements with kelvin voigt damping aroung the fault
//
// holds for:
//  ATTENUATION               = .false.
//  ANISOTROPY                = .false.
//  COMPUTE_AND_STORE_STRAIN  = .true. or .false. (true for kernel simulations)
//  gravity                   = .false.
//  use_mesh_coloring_gpu     = .false.
//  COMPUTE_AND_STORE_STRAIN  = .false.
//  mp -> Kelvin_Voigt_damping = .true.

  // block-id == number of local element id in phase_ispec array
  int bx = blockIdx.y*gridDim.x+blockIdx.x;

  // thread-id == GLL node id
  // note: use only NGLL^3 = 125 active threads, plus 3 inactive/ghost threads,
  //       because we used memory padding from NGLL^3 = 125 to 128 to get coalescent memory accesses;
  //       to avoid execution branching and the need of registers to store an active state variable,
  //       the thread ids are put in valid range
  int tx = threadIdx.x;

  int I,J,K;
  int iglob,offset;
  int working_element;

  realw tempx1l,tempx2l,tempx3l,tempy1l,tempy2l,tempy3l,tempz1l,tempz2l,tempz3l;
  realw xixl,xiyl,xizl,etaxl,etayl,etazl,gammaxl,gammayl,gammazl,jacobianl;
  realw duxdxl,duxdyl,duxdzl,duydxl,duydyl,duydzl,duzdxl,duzdyl,duzdzl;
  realw duxdxl_plus_duydyl,duxdxl_plus_duzdzl,duydyl_plus_duzdzl;
  realw duxdyl_plus_duydxl,duzdxl_plus_duxdzl,duzdyl_plus_duydzl;
  realw kelvin_voigt_eta;
  realw fac1,fac2,fac3;
  realw lambdal,mul,lambdalplus2mul,kappal;
  realw sigma_xx,sigma_yy,sigma_zz,sigma_xy,sigma_xz,sigma_yz;
  realw sum_terms1,sum_terms2,sum_terms3;

  // shared memory
  __shared__ realw sh_tempx[NGLL3];
  __shared__ realw sh_tempy[NGLL3];
  __shared__ realw sh_tempz[NGLL3];

  // note: using shared memory for hprime's improves performance
  //       (but could tradeoff with occupancy)
  __shared__ realw sh_hprime_xx[NGLL2];
  __shared__ realw sh_hprimewgll_xx[NGLL2];

// arithmetic intensity: ratio of number-of-arithmetic-operations / number-of-bytes-accessed-on-DRAM
//
// hand-counts on floating-point operations: counts addition/subtraction/multiplication/division
//                                           no counts for operations on indices in for-loops (compiler will likely unrool loops)
//
//                                           counts accesses to global memory, but no shared memory or register loads/stores
//                                           float has 4 bytes

// counts:
// 2 FLOP
  // checks if anything to do
  if (bx >= nb_blocks_to_compute) return;

  // limits thread ids to range [0,125-1]
  if (tx >= NGLL3) tx = NGLL3 - 1;

// counts:
// + 1 FLOP
//
// + 0 BYTE

  // loads hprime's into shared memory
  if (tx < NGLL2) {
    // copy hprime from global memory to shared memory
    load_shared_memory_hprime(&tx,d_hprime_xx,sh_hprime_xx);

    // copy hprimewgll from global memory to shared memory
    load_shared_memory_hprimewgll(&tx,d_hprimewgll_xx,sh_hprimewgll_xx);
  }

// counts:
// + 0 FLOP
//
// 2 * 1 float * 25 threads = 200 BYTE

  // spectral-element id
  // iphase-1 and working_element-1 for Fortran->C array conventions
  working_element = d_phase_ispec_inner_elastic[bx + num_phase_ispec_elastic*(d_iphase-1)] - 1;

  // local padded index
  offset = working_element*NGLL3_PADDED + tx;

  // global index
  iglob = d_ibool[offset] - 1 ;
  // fetch the value of kelvin_voigt eta
  kelvin_voigt_eta = d_kelvin_voigt_eta[bx] ;


// counts:
// + 7 FLOP
//
// + 2 float * 128 threads = 1024 BYTE

  // copy from global memory to shared memory
  // each thread writes one of the NGLL^3 = 125 data points
  if (threadIdx.x < NGLL3 ){
    // copy displacement from global memory to shared memory
    load_shared_memory_displ_visco<FORWARD_OR_ADJOINT>(&tx,&iglob,d_displ,d_veloc,kelvin_voigt_eta,sh_tempx,sh_tempy,sh_tempz);
  }

// counts:
// + 5 FLOP
//
// + 3 float * 125 threads = 1500 BYTE

  kappal = d_kappav[offset];
  mul = d_muv[offset];

// counts:
// + 0 FLOP
//
// + 2 * 1 float * 128 threads = 1024 BYTE

  // loads mesh values here to give compiler possibility to overlap memory fetches with some computations
  // note: arguments defined as realw* instead of const realw* __restrict__ to avoid that the compiler
  //       loads all memory by texture loads
  //       we only use the first loads explicitly by texture loads, all subsequent without. this should lead/trick
  //       the compiler to use global memory loads for all the subsequent accesses.
  //
  // calculates laplacian
  xixl = get_global_cr( &d_xix[offset] ); // first array with texture load
  xiyl = get_global_cr( &d_xiy[offset] ); // first array with texture load
  xizl = get_global_cr( &d_xiz[offset] ); // first array with texture load

//  xixl = d_xix[offset]; // first array with texture load
//  xiyl = d_xiy[offset]; // all subsequent without to avoid over-use of texture for coalescent access
//  xizl = d_xiz[offset];

  etaxl = d_etax[offset];
  etayl = d_etay[offset];
  etazl = d_etaz[offset];

  gammaxl = d_gammax[offset];
  gammayl = d_gammay[offset];
  gammazl = d_gammaz[offset];

  jacobianl = 1.f / (xixl*(etayl*gammazl-etazl*gammayl)
                    -xiyl*(etaxl*gammazl-etazl*gammaxl)
                    +xizl*(etaxl*gammayl-etayl*gammaxl));

// counts:
// + 15 FLOP
//
// + 9 float * 128 threads = 4608 BYTE



  // local index
  K = (tx/NGLL2);
  J = ((tx-K*NGLL2)/NGLLX);
  I = (tx-K*NGLL2-J*NGLLX);

// counts:
// + 8 FLOP
//
// + 0 BYTE

  // synchronize all the threads (one thread for each of the NGLL grid points of the
  // current spectral element) because we need the whole element to be ready in order
  // to be able to compute the matrix products along cut planes of the 3D element below
  __syncthreads();

  // computes first matrix products
  // 1. cut-plane
  sum_hprime_xi(I,J,K,&tempx1l,&tempy1l,&tempz1l,sh_tempx,sh_tempy,sh_tempz,sh_hprime_xx);
  // 2. cut-plane
  sum_hprime_eta(I,J,K,&tempx2l,&tempy2l,&tempz2l,sh_tempx,sh_tempy,sh_tempz,sh_hprime_xx);
  // 3. cut-plane
  sum_hprime_gamma(I,J,K,&tempx3l,&tempy3l,&tempz3l,sh_tempx,sh_tempy,sh_tempz,sh_hprime_xx);

// counts:
// + 3 * 100 FLOP = 300 FLOP
//
// + 0 BYTE

  // compute derivatives of ux, uy and uz with respect to x, y and z
  duxdxl = xixl*tempx1l + etaxl*tempx2l + gammaxl*tempx3l;
  duxdyl = xiyl*tempx1l + etayl*tempx2l + gammayl*tempx3l;
  duxdzl = xizl*tempx1l + etazl*tempx2l + gammazl*tempx3l;

  duydxl = xixl*tempy1l + etaxl*tempy2l + gammaxl*tempy3l;
  duydyl = xiyl*tempy1l + etayl*tempy2l + gammayl*tempy3l;
  duydzl = xizl*tempy1l + etazl*tempy2l + gammazl*tempy3l;

  duzdxl = xixl*tempz1l + etaxl*tempz2l + gammaxl*tempz3l;
  duzdyl = xiyl*tempz1l + etayl*tempz2l + gammayl*tempz3l;
  duzdzl = xizl*tempz1l + etazl*tempz2l + gammazl*tempz3l;

// counts:
// + 9 * 5 FLOP = 45 FLOP
//
// + 0 BYTE

  // precompute some sums to save CPU time
  duxdxl_plus_duydyl = duxdxl + duydyl;
  duxdxl_plus_duzdzl = duxdxl + duzdzl;
  duydyl_plus_duzdzl = duydyl + duzdzl;
  duxdyl_plus_duydxl = duxdyl + duydxl;
  duzdxl_plus_duxdzl = duzdxl + duxdzl;
  duzdyl_plus_duydzl = duzdyl + duydzl;

  // stress calculations

  // isotropic case
  // compute elements with an elastic isotropic rheology

  lambdalplus2mul = kappal + 1.33333333333333333333f * mul;  // 4./3. = 1.3333333
  lambdal = lambdalplus2mul - 2.0f * mul;

  // compute the six components of the stress tensor sigma
  sigma_xx = lambdalplus2mul*duxdxl + lambdal*duydyl_plus_duzdzl;
  sigma_yy = lambdalplus2mul*duydyl + lambdal*duxdxl_plus_duzdzl;
  sigma_zz = lambdalplus2mul*duzdzl + lambdal*duxdxl_plus_duydyl;

  sigma_xy = mul*duxdyl_plus_duydxl;
  sigma_xz = mul*duzdxl_plus_duxdzl;
  sigma_yz = mul*duzdyl_plus_duydzl;

// counts:
// + 22 FLOP
//
// + 0 BYTE

  // form dot product with test vector, symmetric form
  // 1. cut-plane xi
  __syncthreads();
  // fills shared memory arrays
  if (threadIdx.x < NGLL3) {
    sh_tempx[tx] = jacobianl * (sigma_xx*xixl + sigma_xy*xiyl + sigma_xz*xizl); // sh_tempx1
    sh_tempy[tx] = jacobianl * (sigma_xy*xixl + sigma_yy*xiyl + sigma_yz*xizl); // sh_tempy1
    sh_tempz[tx] = jacobianl * (sigma_xz*xixl + sigma_yz*xiyl + sigma_zz*xizl); // sh_tempz1
  }
  __syncthreads();
  // 1. cut-plane xi
  sum_hprimewgll_xi(I,J,K,&tempx1l,&tempy1l,&tempz1l,sh_tempx,sh_tempy,sh_tempz,sh_hprimewgll_xx);

  // 2. cut-plane eta
  __syncthreads();
  // fills shared memory arrays
  if (threadIdx.x < NGLL3) {
    sh_tempx[tx] = jacobianl * (sigma_xx*etaxl + sigma_xy*etayl + sigma_xz*etazl); // sh_tempx2
    sh_tempy[tx] = jacobianl * (sigma_xy*etaxl + sigma_yy*etayl + sigma_yz*etazl); // sh_tempy2
    sh_tempz[tx] = jacobianl * (sigma_xz*etaxl + sigma_yz*etayl + sigma_zz*etazl); // sh_tempz2
  }
  __syncthreads();
  // 2. cut-plane eta
  sum_hprimewgll_eta(I,J,K,&tempx2l,&tempy2l,&tempz2l,sh_tempx,sh_tempy,sh_tempz,sh_hprimewgll_xx);

  // 3. cut-plane gamma
  __syncthreads();
  // fills shared memory arrays
  if (threadIdx.x < NGLL3) {
    sh_tempx[tx] = jacobianl * (sigma_xx*gammaxl + sigma_xy*gammayl + sigma_xz*gammazl); // sh_tempx3
    sh_tempy[tx] = jacobianl * (sigma_xy*gammaxl + sigma_yy*gammayl + sigma_yz*gammazl); // sh_tempy3
    sh_tempz[tx] = jacobianl * (sigma_xz*gammaxl + sigma_yz*gammayl + sigma_zz*gammazl); // sh_tempz3
  }
  __syncthreads();
  // 3. cut-plane gamma
  sum_hprimewgll_gamma(I,J,K,&tempx3l,&tempy3l,&tempz3l,sh_tempx,sh_tempy,sh_tempz,sh_hprimewgll_xx);

// counts:
// + 3 * 3 * 6 FLOP = 54 FLOP
// + 3 * 100 FLOP = 300 FLOP
//
// + 0 BYTE

  // gets double weights
  fac1 = d_wgllwgll_yz[K*NGLLX+J];
  fac2 = d_wgllwgll_xz[K*NGLLX+I];
  fac3 = d_wgllwgll_xy[J*NGLLX+I];

// counts:
// + 3 * 2 FLOP = 6 FLOP
//
// + 3 float * 128 threads = 1536 BYTE

  sum_terms1 = - (fac1*tempx1l + fac2*tempx2l + fac3*tempx3l);
  sum_terms2 = - (fac1*tempy1l + fac2*tempy2l + fac3*tempy3l);
  sum_terms3 = - (fac1*tempz1l + fac2*tempz2l + fac3*tempz3l);

// counts:
// + 3 * 6 FLOP = 18 FLOP
//
// + 0 BYTE

  // assembles acceleration array
  if (threadIdx.x < NGLL3) {
    atomicAdd(&d_accel[iglob*3], sum_terms1);
    atomicAdd(&d_accel[iglob*3+1], sum_terms2);
    atomicAdd(&d_accel[iglob*3+2], sum_terms3);
  }
} // kernel_2_noatt_iso_kelvinvoigt_impl()

/* ----------------------------------------------------------------------------------------------- */






void Kernel_2(int nb_blocks_to_compute,Mesh* mp,int d_iphase,realw d_deltat,
              int COMPUTE_AND_STORE_STRAIN,
              int ATTENUATION,
              int ANISOTROPY,
              int* d_ibool,
              realw* d_xix,realw* d_xiy,realw* d_xiz,
              realw* d_etax,realw* d_etay,realw* d_etaz,
              realw* d_gammax,realw* d_gammay,realw* d_gammaz,
              realw* d_kappav,
              realw* d_muv,
              realw* d_epsilondev_xx,realw* d_epsilondev_yy,realw* d_epsilondev_xy,
              realw* d_epsilondev_xz,realw* d_epsilondev_yz,
              realw* d_epsilon_trace_over_3,
              realw* d_one_minus_sum_beta,
              realw* d_factor_common,
              realw* d_R_xx,realw* d_R_yy,realw* d_R_xy,
              realw* d_R_xz,realw* d_R_yz,
              realw* d_b_epsilondev_xx,realw* d_b_epsilondev_yy,realw* d_b_epsilondev_xy,
              realw* d_b_epsilondev_xz,realw* d_b_epsilondev_yz,
              realw* d_b_epsilon_trace_over_3,
              realw* d_b_R_xx,realw* d_b_R_yy,realw* d_b_R_xy,
              realw* d_b_R_xz,realw* d_b_R_yz,
              realw* d_c11store,realw* d_c12store,realw* d_c13store,
              realw* d_c14store,realw* d_c15store,realw* d_c16store,
              realw* d_c22store,realw* d_c23store,realw* d_c24store,
              realw* d_c25store,realw* d_c26store,realw* d_c33store,
              realw* d_c34store,realw* d_c35store,realw* d_c36store,
              realw* d_c44store,realw* d_c45store,realw* d_c46store,
              realw* d_c55store,realw* d_c56store,realw* d_c66store,
              realw* d_rhostore){

  TRACE("\tKernel_2");

#ifdef ENABLE_VERY_SLOW_ERROR_CHECKING
  exit_on_cuda_error("before kernel Kernel 2");
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
  cudaEvent_t start,stop;
  if (CUDA_TIMING ){
    start_timing_cuda(&start,&stop);
  }

  // cuda kernel call
  if (ATTENUATION ){
    // compute kernels with attenuation
    // forward wavefields -> FORWARD_OR_ADJOINT == 1
    Kernel_2_att_impl<1><<<grid,threads,0,mp->compute_stream>>>(nb_blocks_to_compute,
                                                                d_ibool,
                                                                mp->d_phase_ispec_inner_elastic,
                                                                mp->num_phase_ispec_elastic,
                                                                d_iphase,
                                                                mp->use_mesh_coloring_gpu,
                                                                d_deltat,
                                                                mp->d_displ,mp->d_veloc,mp->d_accel,
                                                                d_xix, d_xiy, d_xiz,
                                                                d_etax, d_etay, d_etaz,
                                                                d_gammax, d_gammay, d_gammaz,
                                                                mp->d_hprime_xx,
                                                                mp->d_hprimewgll_xx,
                                                                mp->d_wgllwgll_xy, mp->d_wgllwgll_xz, mp->d_wgllwgll_yz,
                                                                d_kappav, d_muv,
                                                                d_epsilondev_xx,d_epsilondev_yy,d_epsilondev_xy,
                                                                d_epsilondev_xz,d_epsilondev_yz,
                                                                d_epsilon_trace_over_3,
                                                                mp->simulation_type,
                                                                mp->NSPEC_AB,
                                                                d_one_minus_sum_beta,
                                                                d_factor_common,
                                                                d_R_xx,d_R_yy,d_R_xy,d_R_xz,d_R_yz,
                                                                mp->d_alphaval,mp->d_betaval,mp->d_gammaval,
                                                                ANISOTROPY,
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
                                                                mp->d_wgll_cube);

    if (mp->simulation_type == 3) {
      // backward/reconstructed wavefields -> FORWARD_OR_ADJOINT == 3
      Kernel_2_att_impl<3><<< grid,threads,0,mp->compute_stream>>>(nb_blocks_to_compute,
                                                                   d_ibool,
                                                                   mp->d_phase_ispec_inner_elastic,
                                                                   mp->num_phase_ispec_elastic,
                                                                   d_iphase,
                                                                   mp->use_mesh_coloring_gpu,
                                                                   d_deltat,
                                                                   mp->d_b_displ,mp->d_b_veloc,mp->d_b_accel,
                                                                   d_xix, d_xiy, d_xiz,
                                                                   d_etax, d_etay, d_etaz,
                                                                   d_gammax, d_gammay, d_gammaz,
                                                                   mp->d_hprime_xx,
                                                                   mp->d_hprimewgll_xx,
                                                                   mp->d_wgllwgll_xy, mp->d_wgllwgll_xz, mp->d_wgllwgll_yz,
                                                                   d_kappav, d_muv,
                                                                   d_b_epsilondev_xx,d_b_epsilondev_yy,d_b_epsilondev_xy,
                                                                   d_b_epsilondev_xz,d_b_epsilondev_yz,
                                                                   d_b_epsilon_trace_over_3,
                                                                   mp->simulation_type,
                                                                   mp->NSPEC_AB,
                                                                   d_one_minus_sum_beta,
                                                                   d_factor_common,
                                                                   d_b_R_xx,d_b_R_yy,d_b_R_xy,d_b_R_xz,d_b_R_yz,
                                                                   mp->d_b_alphaval,mp->d_b_betaval,mp->d_b_gammaval,
                                                                   ANISOTROPY,
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
                                                                   mp->d_wgll_cube);
    }
  }else{
    // compute kernels without attenuation
    if (ANISOTROPY ){
      // full anisotropy
      // forward wavefields -> FORWARD_OR_ADJOINT == 1
      Kernel_2_noatt_ani_impl<1><<<grid,threads,0,mp->compute_stream>>>(nb_blocks_to_compute,
                                                                        d_ibool,
                                                                        mp->d_phase_ispec_inner_elastic,mp->num_phase_ispec_elastic,
                                                                        d_iphase,
                                                                        mp->use_mesh_coloring_gpu,
                                                                        mp->d_displ,
                                                                        mp->d_accel,
                                                                        d_xix, d_xiy, d_xiz,
                                                                        d_etax, d_etay, d_etaz,
                                                                        d_gammax, d_gammay, d_gammaz,
                                                                        mp->d_hprime_xx,
                                                                        mp->d_hprimewgll_xx,
                                                                        mp->d_wgllwgll_xy, mp->d_wgllwgll_xz, mp->d_wgllwgll_yz,
                                                                        d_kappav, d_muv,
                                                                        COMPUTE_AND_STORE_STRAIN,
                                                                        d_epsilondev_xx,d_epsilondev_yy,d_epsilondev_xy,
                                                                        d_epsilondev_xz,d_epsilondev_yz,
                                                                        d_epsilon_trace_over_3,
                                                                        mp->simulation_type,
                                                                        ANISOTROPY,
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
                                                                        mp->d_wgll_cube );

      // backward/reconstructed wavefield
      if (mp->simulation_type == 3) {
        // backward/reconstructed wavefields -> FORWARD_OR_ADJOINT == 3
        Kernel_2_noatt_ani_impl<3><<< grid,threads,0,mp->compute_stream>>>(nb_blocks_to_compute,
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
                                                                           ANISOTROPY,
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
                                                                           mp->d_wgll_cube );
      }
    }else{
      // isotropic case
      if (mp->gravity){
        // with gravity
        // forward wavefields -> FORWARD_OR_ADJOINT == 1
        Kernel_2_noatt_iso_grav_impl<1><<<grid,threads,0,mp->compute_stream>>>(nb_blocks_to_compute,
                                                                          d_ibool,
                                                                          mp->d_phase_ispec_inner_elastic,mp->num_phase_ispec_elastic,
                                                                          d_iphase,
                                                                          mp->use_mesh_coloring_gpu,
                                                                          mp->d_displ,
                                                                          mp->d_accel,
                                                                          d_xix, d_xiy, d_xiz,
                                                                          d_etax, d_etay, d_etaz,
                                                                          d_gammax, d_gammay, d_gammaz,
                                                                          mp->d_hprime_xx,
                                                                          mp->d_hprimewgll_xx,
                                                                          mp->d_wgllwgll_xy, mp->d_wgllwgll_xz, mp->d_wgllwgll_yz,
                                                                          d_kappav, d_muv,
                                                                          COMPUTE_AND_STORE_STRAIN,
                                                                          d_epsilondev_xx,d_epsilondev_yy,d_epsilondev_xy,
                                                                          d_epsilondev_xz,d_epsilondev_yz,
                                                                          d_epsilon_trace_over_3,
                                                                          mp->simulation_type,
                                                                          mp->gravity,
                                                                          mp->d_minus_g,
                                                                          mp->d_minus_deriv_gravity,
                                                                          d_rhostore,
                                                                          mp->d_wgll_cube );

        // backward/reconstructed wavefield
        if (mp->simulation_type == 3) {
          // backward/reconstructed wavefields -> FORWARD_OR_ADJOINT == 3
          Kernel_2_noatt_iso_grav_impl<3><<< grid,threads,0,mp->compute_stream>>>(nb_blocks_to_compute,
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
                                                                             mp->gravity,
                                                                             mp->d_minus_g,
                                                                             mp->d_minus_deriv_gravity,
                                                                             d_rhostore,
                                                                             mp->d_wgll_cube );
        }
      }else{
        // without gravity
        if (mp->use_mesh_coloring_gpu) {
          // with mesh coloring
          // forward wavefields -> FORWARD_OR_ADJOINT == 1
          Kernel_2_noatt_iso_col_impl<1><<<grid,threads,0,mp->compute_stream>>>(nb_blocks_to_compute,
                                                                            d_ibool,
                                                                            mp->d_phase_ispec_inner_elastic,mp->num_phase_ispec_elastic,
                                                                            d_iphase,
                                                                            mp->use_mesh_coloring_gpu,
                                                                            mp->d_displ,
                                                                            mp->d_accel,
                                                                            d_xix, d_xiy, d_xiz,
                                                                            d_etax, d_etay, d_etaz,
                                                                            d_gammax, d_gammay, d_gammaz,
                                                                            mp->d_hprime_xx,
                                                                            mp->d_hprimewgll_xx,
                                                                            mp->d_wgllwgll_xy, mp->d_wgllwgll_xz, mp->d_wgllwgll_yz,
                                                                            d_kappav, d_muv,
                                                                            COMPUTE_AND_STORE_STRAIN,
                                                                            d_epsilondev_xx,d_epsilondev_yy,d_epsilondev_xy,
                                                                            d_epsilondev_xz,d_epsilondev_yz,
                                                                            d_epsilon_trace_over_3,
                                                                            mp->simulation_type);

          // backward/reconstructed wavefield
          if (mp->simulation_type == 3) {
            // backward/reconstructed wavefields -> FORWARD_OR_ADJOINT == 3
            Kernel_2_noatt_iso_col_impl<3><<< grid,threads,0,mp->compute_stream>>>(nb_blocks_to_compute,
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
                                                                               mp->simulation_type);
          }
        }else{
          // without mesh coloring
          if (COMPUTE_AND_STORE_STRAIN ){
            // stores strains
            // forward wavefields -> FORWARD_OR_ADJOINT == 1
            Kernel_2_noatt_iso_strain_impl<1><<<grid,threads,0,mp->compute_stream>>>(nb_blocks_to_compute,
                                                                              d_ibool,
                                                                              mp->d_phase_ispec_inner_elastic,mp->num_phase_ispec_elastic,
                                                                              d_iphase,
                                                                              mp->d_displ,
                                                                              mp->d_accel,
                                                                              d_xix, d_xiy, d_xiz,
                                                                              d_etax, d_etay, d_etaz,
                                                                              d_gammax, d_gammay, d_gammaz,
                                                                              mp->d_hprime_xx,
                                                                              mp->d_hprimewgll_xx,
                                                                              mp->d_wgllwgll_xy, mp->d_wgllwgll_xz, mp->d_wgllwgll_yz,
                                                                              d_kappav, d_muv,
                                                                              COMPUTE_AND_STORE_STRAIN,
                                                                              d_epsilondev_xx,d_epsilondev_yy,d_epsilondev_xy,
                                                                              d_epsilondev_xz,d_epsilondev_yz,
                                                                              d_epsilon_trace_over_3,
                                                                              mp->simulation_type);

            // backward/reconstructed wavefield
            if (mp->simulation_type == 3) {
              // backward/reconstructed wavefields -> FORWARD_OR_ADJOINT == 3
              Kernel_2_noatt_iso_strain_impl<3><<< grid,threads,0,mp->compute_stream>>>(nb_blocks_to_compute,
                                                                                 d_ibool,
                                                                                 mp->d_phase_ispec_inner_elastic,mp->num_phase_ispec_elastic,
                                                                                 d_iphase,
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
                                                                                 mp->simulation_type);
            }
          }else{
            if (mp->Kelvin_Voigt_damping) {
                         // Kelvin_Voigt_damping == true means there is fault in this partition
              Kernel_2_noatt_iso_kelvinvoigt_impl<1><<<grid,threads,0,mp->compute_stream>>>(nb_blocks_to_compute,
                                                                                           d_ibool,
                                                                                           mp->d_phase_ispec_inner_elastic,mp->num_phase_ispec_elastic,
                                                                                           d_iphase,
                                                                                           mp->d_Kelvin_Voigt_eta,
                                                                                           mp->d_displ,
                                                                                           mp->d_veloc,
                                                                                           mp->d_accel,
                                                                                           d_xix, d_xiy, d_xiz,
                                                                                           d_etax, d_etay, d_etaz,
                                                                                           d_gammax, d_gammay, d_gammaz,
                                                                                           mp->d_hprime_xx,
                                                                                           mp->d_hprimewgll_xx,
                                                                                           mp->d_wgllwgll_xy, mp->d_wgllwgll_xz, mp->d_wgllwgll_yz,
                                                                                           d_kappav, d_muv);
                       }
            else{
            // without storing strains
            // forward wavefields -> FORWARD_OR_ADJOINT == 1
            Kernel_2_noatt_iso_impl<1><<<grid,threads,0,mp->compute_stream>>>(nb_blocks_to_compute,
                                                                              d_ibool,
                                                                              mp->d_phase_ispec_inner_elastic,mp->num_phase_ispec_elastic,
                                                                              d_iphase,
                                                                              mp->d_displ,
                                                                              mp->d_accel,
                                                                              d_xix, d_xiy, d_xiz,
                                                                              d_etax, d_etay, d_etaz,
                                                                              d_gammax, d_gammay, d_gammaz,
                                                                              mp->d_hprime_xx,
                                                                              mp->d_hprimewgll_xx,
                                                                              mp->d_wgllwgll_xy, mp->d_wgllwgll_xz, mp->d_wgllwgll_yz,
                                                                              d_kappav, d_muv);

            // backward/reconstructed wavefield
            if (mp->simulation_type == 3) {
              // backward/reconstructed wavefields -> FORWARD_OR_ADJOINT == 3
              Kernel_2_noatt_iso_impl<3><<< grid,threads,0,mp->compute_stream>>>(nb_blocks_to_compute,
                                                                                 d_ibool,
                                                                                 mp->d_phase_ispec_inner_elastic,mp->num_phase_ispec_elastic,
                                                                                 d_iphase,
                                                                                 mp->d_b_displ,
                                                                                 mp->d_b_accel,
                                                                                 d_xix, d_xiy, d_xiz,
                                                                                 d_etax, d_etay, d_etaz,
                                                                                 d_gammax, d_gammay, d_gammaz,
                                                                                 mp->d_hprime_xx,
                                                                                 mp->d_hprimewgll_xx,
                                                                                 mp->d_wgllwgll_xy, mp->d_wgllwgll_xz, mp->d_wgllwgll_yz,
                                                                                 d_kappav, d_muv);
            }
            }
          } // COMPUTE_AND_STORE_STRAIN
        } // use_mesh_coloring_gpu
      } // gravity
    } // ANISOTROPY
  } // ATTENUATION

  // Cuda timing
  if (CUDA_TIMING ){
    if (ATTENUATION ){
      stop_timing_cuda(&start,&stop,"Kernel_2_att_impl");
    }else{
      if (ANISOTROPY ){
        stop_timing_cuda(&start,&stop,"Kernel_2_noatt_ani_impl");
      }else{
        if (mp->gravity ){
          stop_timing_cuda(&start,&stop,"Kernel_2_noatt_iso_grav_impl");
        }else{
          if (COMPUTE_AND_STORE_STRAIN ){
            stop_timing_cuda(&start,&stop,"Kernel_2_noatt_iso_strain_impl");
          }else{
            realw time;
            stop_timing_cuda(&start,&stop,"Kernel_2_noatt_iso_impl",&time);
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

#ifdef ENABLE_VERY_SLOW_ERROR_CHECKING
  exit_on_cuda_error("Kernel_2_impl");
#endif
}

/* ----------------------------------------------------------------------------------------------- */


extern "C"
void FC_FUNC_(compute_forces_viscoelastic_cuda,
              COMPUTE_FORCES_VISCOELASTIC_CUDA)(long* Mesh_pointer,
                                                int* iphase,
                                                realw* deltat,
                                                int* nspec_outer_elastic,
                                                int* nspec_inner_elastic,
                                                int* COMPUTE_AND_STORE_STRAIN,
                                                int* ATTENUATION,
                                                int* ANISOTROPY) {

  TRACE("\tcompute_forces_viscoelastic_cuda");
  // EPIK_TRACER("compute_forces_viscoelastic_cuda");
  //printf("Running compute_forces\n");
  //double start_time = get_time();

  Mesh* mp = (Mesh*)(*Mesh_pointer); // get Mesh from fortran integer wrapper

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
               *ATTENUATION,*ANISOTROPY,
               mp->d_ibool + offset,
               mp->d_xix + offset,mp->d_xiy + offset,mp->d_xiz + offset,
               mp->d_etax + offset,mp->d_etay + offset,mp->d_etaz + offset,
               mp->d_gammax + offset,mp->d_gammay + offset,mp->d_gammaz + offset,
               mp->d_kappav + offset,
               mp->d_muv + offset,
               mp->d_epsilondev_xx + offset_nonpadded,mp->d_epsilondev_yy + offset_nonpadded,mp->d_epsilondev_xy + offset_nonpadded,
               mp->d_epsilondev_xz + offset_nonpadded,mp->d_epsilondev_yz + offset_nonpadded,
               mp->d_epsilon_trace_over_3 + offset_nonpadded,
               mp->d_one_minus_sum_beta + offset_nonpadded,
               mp->d_factor_common + offset_nonpadded_att2,
               mp->d_R_xx + offset_nonpadded,mp->d_R_yy + offset_nonpadded,mp->d_R_xy + offset_nonpadded,
               mp->d_R_xz + offset_nonpadded,mp->d_R_yz + offset_nonpadded,
               mp->d_b_epsilondev_xx + offset_nonpadded,mp->d_b_epsilondev_yy + offset_nonpadded,mp->d_b_epsilondev_xy + offset_nonpadded,
               mp->d_b_epsilondev_xz + offset_nonpadded,mp->d_b_epsilondev_yz + offset_nonpadded,
               mp->d_b_epsilon_trace_over_3 + offset_nonpadded,
               mp->d_b_R_xx + offset_nonpadded,mp->d_b_R_yy + offset_nonpadded,mp->d_b_R_xy + offset_nonpadded,
               mp->d_b_R_xz + offset_nonpadded,mp->d_b_R_yz + offset_nonpadded,
               mp->d_c11store + offset,mp->d_c12store + offset,mp->d_c13store + offset,
               mp->d_c14store + offset,mp->d_c15store + offset,mp->d_c16store + offset,
               mp->d_c22store + offset,mp->d_c23store + offset,mp->d_c24store + offset,
               mp->d_c25store + offset,mp->d_c26store + offset,mp->d_c33store + offset,
               mp->d_c34store + offset,mp->d_c35store + offset,mp->d_c36store + offset,
               mp->d_c44store + offset,mp->d_c45store + offset,mp->d_c46store + offset,
               mp->d_c55store + offset,mp->d_c56store + offset,mp->d_c66store + offset,
               mp->d_rhostore + offset);

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
             *ATTENUATION,*ANISOTROPY,
             mp->d_ibool,
             mp->d_xix,mp->d_xiy,mp->d_xiz,
             mp->d_etax,mp->d_etay,mp->d_etaz,
             mp->d_gammax,mp->d_gammay,mp->d_gammaz,
             mp->d_kappav,
             mp->d_muv,
             mp->d_epsilondev_xx,mp->d_epsilondev_yy,mp->d_epsilondev_xy,
             mp->d_epsilondev_xz,mp->d_epsilondev_yz,
             mp->d_epsilon_trace_over_3,
             mp->d_one_minus_sum_beta,
             mp->d_factor_common,
             mp->d_R_xx,mp->d_R_yy,mp->d_R_xy,
             mp->d_R_xz,mp->d_R_yz,
             mp->d_b_epsilondev_xx,mp->d_b_epsilondev_yy,mp->d_b_epsilondev_xy,
             mp->d_b_epsilondev_xz,mp->d_b_epsilondev_yz,
             mp->d_b_epsilon_trace_over_3,
             mp->d_b_R_xx,mp->d_b_R_yy,mp->d_b_R_xy,
             mp->d_b_R_xz,mp->d_b_R_yz,
             mp->d_c11store,mp->d_c12store,mp->d_c13store,
             mp->d_c14store,mp->d_c15store,mp->d_c16store,
             mp->d_c22store,mp->d_c23store,mp->d_c24store,
             mp->d_c25store,mp->d_c26store,mp->d_c33store,
             mp->d_c34store,mp->d_c35store,mp->d_c36store,
             mp->d_c44store,mp->d_c45store,mp->d_c46store,
             mp->d_c55store,mp->d_c56store,mp->d_c66store,
             mp->d_rhostore);
  }
}

