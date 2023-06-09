/*
!=====================================================================
!
!                          S p e c f e m 3 D
!                          -----------------
!
!    Main historical authors: Dimitri Komatitsch and Jeroen Tromp
!                             CNRS, France
!                      and Princeton University, USA
!                (there are currently many more authors!)
!                          (c) October 2017
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


#ifdef USE_TEXTURES_FIELDS
extern realw_texture d_displ_tex;
extern realw_texture d_veloc_tex;
extern realw_texture d_accel_tex;
//backward/reconstructed
extern realw_texture d_b_displ_tex;
extern realw_texture d_b_veloc_tex;
extern realw_texture d_b_accel_tex;

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
extern realw_texture d_hprime_xx_tex;
#endif

/* ----------------------------------------------------------------------------------------------- */

// updates memory variables for convolution & stress

__device__ __forceinline__ void pml_compute_memory_variables_elastic(int tx,int working_element,
                                                                     int I, int J, int K,
                                                                     int ispec_CPML,
                                                                     int NSPEC_CPML,
                                                                     realw_const_p pml_convolution_coef_alpha,
                                                                     realw_const_p pml_convolution_coef_beta,
                                                                     realw_const_p pml_convolution_coef_strain,
                                                                     realw* rmemory_dux_dxl_x,
                                                                     realw* rmemory_duy_dxl_y,
                                                                     realw* rmemory_duz_dxl_z,
                                                                     realw PML_dux_dxl_new,realw PML_dux_dxl_old,
                                                                     realw PML_duy_dxl_new,realw PML_duy_dxl_old,
                                                                     realw PML_duz_dxl_new,realw PML_duz_dxl_old,
                                                                     realw* rmemory_dux_dyl_x,
                                                                     realw* rmemory_duy_dyl_y,
                                                                     realw* rmemory_duz_dyl_z,
                                                                     realw PML_dux_dyl_new,realw PML_dux_dyl_old,
                                                                     realw PML_duy_dyl_new,realw PML_duy_dyl_old,
                                                                     realw PML_duz_dyl_new,realw PML_duz_dyl_old,
                                                                     realw* rmemory_dux_dzl_x,
                                                                     realw* rmemory_duy_dzl_y,
                                                                     realw* rmemory_duz_dzl_z,
                                                                     realw PML_dux_dzl_new,realw PML_dux_dzl_old,
                                                                     realw PML_duy_dzl_new,realw PML_duy_dzl_old,
                                                                     realw PML_duz_dzl_new,realw PML_duz_dzl_old,
                                                                     realw* rmemory_dux_dxl_y,
                                                                     realw* rmemory_dux_dxl_z,
                                                                     realw* rmemory_duy_dxl_x,
                                                                     realw* rmemory_duz_dxl_x,
                                                                     realw* rmemory_dux_dyl_y,
                                                                     realw* rmemory_duy_dyl_x,
                                                                     realw* rmemory_duy_dyl_z,
                                                                     realw* rmemory_duz_dyl_y,
                                                                     realw* rmemory_dux_dzl_z,
                                                                     realw* rmemory_duy_dzl_z,
                                                                     realw* rmemory_duz_dzl_x,
                                                                     realw* rmemory_duz_dzl_y,
                                                                     realw PML_dux_dxl,realw PML_dux_dyl,realw PML_dux_dzl,
                                                                     realw PML_duy_dxl,realw PML_duy_dyl,realw PML_duy_dzl,
                                                                     realw PML_duz_dxl,realw PML_duz_dyl,realw PML_duz_dzl,
                                                                     realw_const_p d_kappav,
                                                                     realw_const_p d_muv,
                                                                     realw* sigma_xx,realw* sigma_yx,realw* sigma_zx,
                                                                     realw* sigma_xy,realw* sigma_yy,realw* sigma_zy,
                                                                     realw* sigma_xz,realw* sigma_yz,realw* sigma_zz) {

  realw A6,A7,A8,A9;
  realw A10,A11,A12,A13;
  realw A14,A15,A16,A17;
  realw A18,A19;
  realw A20,A21;
  realw A22,A23;
  realw coef0_1,coef1_1,coef2_1;
  realw coef0_2,coef1_2,coef2_2;
  realw coef0_3,coef1_3,coef2_3;
  realw duxdxl_x,duxdyl_x,duxdzl_x;
  realw duzdzl_x,duzdxl_x;
  realw duydyl_x,duydxl_x;
  realw duydxl_y,duydyl_y,duydzl_y;
  realw duzdzl_y,duzdyl_y;
  realw duxdyl_y,duxdxl_y;
  realw duzdxl_z,duzdyl_z,duzdzl_z;
  realw duydzl_z,duydyl_z;
  realw duxdzl_z,duxdxl_z;
  realw kappal,mul;
  realw lambdalplus2mul,lambdal;

  int offset_x,offset_y,offset_z;
  int offset;

  //---------------------- A6, A7, A8, A9 --------------------------
  // coefficients
  // alpha_z
  coef0_1 = pml_convolution_coef_alpha[INDEX5(9,NGLLX,NGLLX,NGLLX,6,I,J,K,ispec_CPML)];
  coef1_1 = pml_convolution_coef_alpha[INDEX5(9,NGLLX,NGLLX,NGLLX,7,I,J,K,ispec_CPML)];
  coef2_1 = pml_convolution_coef_alpha[INDEX5(9,NGLLX,NGLLX,NGLLX,8,I,J,K,ispec_CPML)];

  // alpha_y
  coef0_2 = pml_convolution_coef_alpha[INDEX5(9,NGLLX,NGLLX,NGLLX,3,I,J,K,ispec_CPML)];
  coef1_2 = pml_convolution_coef_alpha[INDEX5(9,NGLLX,NGLLX,NGLLX,4,I,J,K,ispec_CPML)];
  coef2_2 = pml_convolution_coef_alpha[INDEX5(9,NGLLX,NGLLX,NGLLX,5,I,J,K,ispec_CPML)];

  // beta_x = alpha_x + d_x / kappa_x
  coef0_3 = pml_convolution_coef_beta[INDEX5(9,NGLLX,NGLLX,NGLLX,0,I,J,K,ispec_CPML)];
  coef1_3 = pml_convolution_coef_beta[INDEX5(9,NGLLX,NGLLX,NGLLX,1,I,J,K,ispec_CPML)];
  coef2_3 = pml_convolution_coef_beta[INDEX5(9,NGLLX,NGLLX,NGLLX,2,I,J,K,ispec_CPML)];

  offset_x = INDEX5(NGLLX,NGLLX,NGLLX,NSPEC_CPML,I,J,K,ispec_CPML,0);
  offset_y = INDEX5(NGLLX,NGLLX,NGLLX,NSPEC_CPML,I,J,K,ispec_CPML,1);
  offset_z = INDEX5(NGLLX,NGLLX,NGLLX,NSPEC_CPML,I,J,K,ispec_CPML,2);

  rmemory_dux_dxl_x[offset_x] = coef0_1 * rmemory_dux_dxl_x[offset_x]
                              + PML_dux_dxl_new * coef1_1 + PML_dux_dxl_old * coef2_1;
  rmemory_duy_dxl_y[offset_x] = coef0_1 * rmemory_duy_dxl_y[offset_x]
                              + PML_duy_dxl_new * coef1_1 + PML_duy_dxl_old * coef2_1;
  rmemory_duz_dxl_z[offset_x] = coef0_1 * rmemory_duz_dxl_z[offset_x]
                              + PML_duz_dxl_new * coef1_1 + PML_duz_dxl_old * coef2_1;

  rmemory_dux_dxl_x[offset_y] = coef0_2 * rmemory_dux_dxl_x[offset_y]
                              + PML_dux_dxl_new * coef1_2 + PML_dux_dxl_old * coef2_2;
  rmemory_duy_dxl_y[offset_y] = coef0_2 * rmemory_duy_dxl_y[offset_y]
                              + PML_duy_dxl_new * coef1_2 + PML_duy_dxl_old * coef2_2;
  rmemory_duz_dxl_z[offset_y] = coef0_2 * rmemory_duz_dxl_z[offset_y]
                              + PML_duz_dxl_new * coef1_2 + PML_duz_dxl_old * coef2_2;

  rmemory_dux_dxl_x[offset_z] = coef0_3 * rmemory_dux_dxl_x[offset_z]
                              + PML_dux_dxl_new * coef1_3 + PML_dux_dxl_old * coef2_3;
  rmemory_duy_dxl_y[offset_z] = coef0_3 * rmemory_duy_dxl_y[offset_z]
                              + PML_duy_dxl_new * coef1_3 + PML_duy_dxl_old * coef2_3;
  rmemory_duz_dxl_z[offset_z] = coef0_3 * rmemory_duz_dxl_z[offset_z]
                              + PML_duz_dxl_new * coef1_3 + PML_duz_dxl_old * coef2_3;

  //---------------------- A10,A11,A12,A13 --------------------------
  // coefficients
  // alpha_x
  coef0_1 = pml_convolution_coef_alpha[INDEX5(9,NGLLX,NGLLX,NGLLX,0,I,J,K,ispec_CPML)];
  coef1_1 = pml_convolution_coef_alpha[INDEX5(9,NGLLX,NGLLX,NGLLX,1,I,J,K,ispec_CPML)];
  coef2_1 = pml_convolution_coef_alpha[INDEX5(9,NGLLX,NGLLX,NGLLX,2,I,J,K,ispec_CPML)];

  // alpha_z
  coef0_2 = pml_convolution_coef_alpha[INDEX5(9,NGLLX,NGLLX,NGLLX,6,I,J,K,ispec_CPML)];
  coef1_2 = pml_convolution_coef_alpha[INDEX5(9,NGLLX,NGLLX,NGLLX,7,I,J,K,ispec_CPML)];
  coef2_2 = pml_convolution_coef_alpha[INDEX5(9,NGLLX,NGLLX,NGLLX,8,I,J,K,ispec_CPML)];

  // beta_y = alpha_y + d_y / kappa_y
  coef0_3 = pml_convolution_coef_beta[INDEX5(9,NGLLX,NGLLX,NGLLX,3,I,J,K,ispec_CPML)];
  coef1_3 = pml_convolution_coef_beta[INDEX5(9,NGLLX,NGLLX,NGLLX,4,I,J,K,ispec_CPML)];
  coef2_3 = pml_convolution_coef_beta[INDEX5(9,NGLLX,NGLLX,NGLLX,5,I,J,K,ispec_CPML)];

  rmemory_dux_dyl_x[offset_x] = coef0_1 * rmemory_dux_dyl_x[offset_x]
                              + PML_dux_dyl_new * coef1_1 + PML_dux_dyl_old * coef2_1;
  rmemory_duy_dyl_y[offset_x] = coef0_1 * rmemory_duy_dyl_y[offset_x]
                              + PML_duy_dyl_new * coef1_1 + PML_duy_dyl_old * coef2_1;
  rmemory_duz_dyl_z[offset_x] = coef0_1 * rmemory_duz_dyl_z[offset_x]
                              + PML_duz_dyl_new * coef1_1 + PML_duz_dyl_old * coef2_1;

  rmemory_dux_dyl_x[offset_y] = coef0_2 * rmemory_dux_dyl_x[offset_y]
                              + PML_dux_dyl_new * coef1_2 + PML_dux_dyl_old * coef2_2;
  rmemory_duy_dyl_y[offset_y] = coef0_2 * rmemory_duy_dyl_y[offset_y]
                              + PML_duy_dyl_new * coef1_2 + PML_duy_dyl_old * coef2_2;
  rmemory_duz_dyl_z[offset_y] = coef0_2 * rmemory_duz_dyl_z[offset_y]
                              + PML_duz_dyl_new * coef1_2 + PML_duz_dyl_old * coef2_2;

  rmemory_dux_dyl_x[offset_z] = coef0_3 * rmemory_dux_dyl_x[offset_z]
                              + PML_dux_dyl_new * coef1_3 + PML_dux_dyl_old * coef2_3;
  rmemory_duy_dyl_y[offset_z] = coef0_3 * rmemory_duy_dyl_y[offset_z]
                              + PML_duy_dyl_new * coef1_3 + PML_duy_dyl_old * coef2_3;
  rmemory_duz_dyl_z[offset_z] = coef0_3 * rmemory_duz_dyl_z[offset_z]
                              + PML_duz_dyl_new * coef1_3 + PML_duz_dyl_old * coef2_3;

  //---------------------- A14,A15,A16,A17 --------------------------
  // coefficients
  // alpha_x
  coef0_1 = pml_convolution_coef_alpha[INDEX5(9,NGLLX,NGLLX,NGLLX,0,I,J,K,ispec_CPML)];
  coef1_1 = pml_convolution_coef_alpha[INDEX5(9,NGLLX,NGLLX,NGLLX,1,I,J,K,ispec_CPML)];
  coef2_1 = pml_convolution_coef_alpha[INDEX5(9,NGLLX,NGLLX,NGLLX,2,I,J,K,ispec_CPML)];

  // alpha_y
  coef0_2 = pml_convolution_coef_alpha[INDEX5(9,NGLLX,NGLLX,NGLLX,3,I,J,K,ispec_CPML)];
  coef1_2 = pml_convolution_coef_alpha[INDEX5(9,NGLLX,NGLLX,NGLLX,4,I,J,K,ispec_CPML)];
  coef2_2 = pml_convolution_coef_alpha[INDEX5(9,NGLLX,NGLLX,NGLLX,5,I,J,K,ispec_CPML)];

  // beta_z = alpha_z + d_z / kappa_z
  coef0_3 = pml_convolution_coef_beta[INDEX5(9,NGLLX,NGLLX,NGLLX,6,I,J,K,ispec_CPML)];
  coef1_3 = pml_convolution_coef_beta[INDEX5(9,NGLLX,NGLLX,NGLLX,7,I,J,K,ispec_CPML)];
  coef2_3 = pml_convolution_coef_beta[INDEX5(9,NGLLX,NGLLX,NGLLX,8,I,J,K,ispec_CPML)];

  rmemory_dux_dzl_x[offset_x] = coef0_1 * rmemory_dux_dzl_x[offset_x]
                              + PML_dux_dzl_new * coef1_1 + PML_dux_dzl_old * coef2_1;
  rmemory_duy_dzl_y[offset_x] = coef0_1 * rmemory_duy_dzl_y[offset_x]
                              + PML_duy_dzl_new * coef1_1 + PML_duy_dzl_old * coef2_1;
  rmemory_duz_dzl_z[offset_x] = coef0_1 * rmemory_duz_dzl_z[offset_x]
                              + PML_duz_dzl_new * coef1_1 + PML_duz_dzl_old * coef2_1;

  rmemory_dux_dzl_x[offset_y] = coef0_2 * rmemory_dux_dzl_x[offset_y]
                              + PML_dux_dzl_new * coef1_2 + PML_dux_dzl_old * coef2_2;
  rmemory_duy_dzl_y[offset_y] = coef0_2 * rmemory_duy_dzl_y[offset_y]
                              + PML_duy_dzl_new * coef1_2 + PML_duy_dzl_old * coef2_2;
  rmemory_duz_dzl_z[offset_y] = coef0_2 * rmemory_duz_dzl_z[offset_y]
                              + PML_duz_dzl_new * coef1_2 + PML_duz_dzl_old * coef2_2;

  rmemory_dux_dzl_x[offset_z] = coef0_3 * rmemory_dux_dzl_x[offset_z]
                              + PML_dux_dzl_new * coef1_3 + PML_dux_dzl_old * coef2_3;
  rmemory_duy_dzl_y[offset_z] = coef0_3 * rmemory_duy_dzl_y[offset_z]
                              + PML_duy_dzl_new * coef1_3 + PML_duy_dzl_old * coef2_3;
  rmemory_duz_dzl_z[offset_z] = coef0_3 * rmemory_duz_dzl_z[offset_z]
                              + PML_duz_dzl_new * coef1_3 + PML_duz_dzl_old * coef2_3;

  //---------------------- A18 and A19 --------------------------
  // coefficients
  // alpha_x
  coef0_1 = pml_convolution_coef_alpha[INDEX5(9,NGLLX,NGLLX,NGLLX,0,I,J,K,ispec_CPML)];
  coef1_1 = pml_convolution_coef_alpha[INDEX5(9,NGLLX,NGLLX,NGLLX,1,I,J,K,ispec_CPML)];
  coef2_1 = pml_convolution_coef_alpha[INDEX5(9,NGLLX,NGLLX,NGLLX,2,I,J,K,ispec_CPML)];

  offset = INDEX4(NGLLX,NGLLX,NGLLX,I,J,K,ispec_CPML);

  rmemory_duz_dzl_y[offset] = coef0_1 * rmemory_duz_dzl_y[offset]
                            + PML_duz_dzl_new * coef1_1 + PML_duz_dzl_old * coef2_1;

  rmemory_duz_dyl_y[offset] = coef0_1 * rmemory_duz_dyl_y[offset]
                            + PML_duz_dyl_new * coef1_1 + PML_duz_dyl_old * coef2_1;

  rmemory_duy_dzl_z[offset] = coef0_1 * rmemory_duy_dzl_z[offset]
                            + PML_duy_dzl_new * coef1_1 + PML_duy_dzl_old * coef2_1;

  rmemory_duy_dyl_z[offset] = coef0_1 * rmemory_duy_dyl_z[offset]
                            + PML_duy_dyl_new * coef1_1 + PML_duy_dyl_old * coef2_1;

  //---------------------- A20 and A21 --------------------------
  // coefficients
  // alpha_y
  coef0_1 = pml_convolution_coef_alpha[INDEX5(9,NGLLX,NGLLX,NGLLX,3,I,J,K,ispec_CPML)];
  coef1_1 = pml_convolution_coef_alpha[INDEX5(9,NGLLX,NGLLX,NGLLX,4,I,J,K,ispec_CPML)];
  coef2_1 = pml_convolution_coef_alpha[INDEX5(9,NGLLX,NGLLX,NGLLX,5,I,J,K,ispec_CPML)];

  rmemory_duz_dzl_x[offset] = coef0_1 * rmemory_duz_dzl_x[offset]
                            + PML_duz_dzl_new * coef1_1 + PML_duz_dzl_old * coef2_1;

  rmemory_duz_dxl_x[offset] = coef0_1 * rmemory_duz_dxl_x[offset]
                            + PML_duz_dxl_new * coef1_1 + PML_duz_dxl_old * coef2_1;

  rmemory_dux_dzl_z[offset] = coef0_1 * rmemory_dux_dzl_z[offset]
                            + PML_dux_dzl_new * coef1_1 + PML_dux_dzl_old * coef2_1;

  rmemory_dux_dxl_z[offset] = coef0_1 * rmemory_dux_dxl_z[offset]
                            + PML_dux_dxl_new * coef1_1 + PML_dux_dxl_old * coef2_1;

  //---------------------- A22 and A23 --------------------------
  // coefficients
  // alpha_z
  coef0_1 = pml_convolution_coef_alpha[INDEX5(9,NGLLX,NGLLX,NGLLX,6,I,J,K,ispec_CPML)];
  coef1_1 = pml_convolution_coef_alpha[INDEX5(9,NGLLX,NGLLX,NGLLX,7,I,J,K,ispec_CPML)];
  coef2_1 = pml_convolution_coef_alpha[INDEX5(9,NGLLX,NGLLX,NGLLX,8,I,J,K,ispec_CPML)];

  rmemory_duy_dyl_x[offset] = coef0_1 * rmemory_duy_dyl_x[offset]
                            + PML_duy_dyl_new * coef1_1 + PML_duy_dyl_old * coef2_1;

  rmemory_duy_dxl_x[offset] = coef0_1 * rmemory_duy_dxl_x[offset]
                            + PML_duy_dxl_new * coef1_1 + PML_duy_dxl_old * coef2_1;

  rmemory_dux_dyl_y[offset] = coef0_1 * rmemory_dux_dyl_y[offset]
                            + PML_dux_dyl_new * coef1_1 + PML_dux_dyl_old * coef2_1;

  rmemory_dux_dxl_y[offset] = coef0_1 * rmemory_dux_dxl_y[offset]
                            + PML_dux_dxl_new * coef1_1 + PML_dux_dxl_old * coef2_1;

  // PML convolution coefficients
  A6 = pml_convolution_coef_strain[INDEX5(18,NGLLX,NGLLX,NGLLX,0,I,J,K,ispec_CPML)];
  A7 = pml_convolution_coef_strain[INDEX5(18,NGLLX,NGLLX,NGLLX,1,I,J,K,ispec_CPML)];
  A8 = pml_convolution_coef_strain[INDEX5(18,NGLLX,NGLLX,NGLLX,2,I,J,K,ispec_CPML)];
  A9 = pml_convolution_coef_strain[INDEX5(18,NGLLX,NGLLX,NGLLX,3,I,J,K,ispec_CPML)];

  A10 = pml_convolution_coef_strain[INDEX5(18,NGLLX,NGLLX,NGLLX,4,I,J,K,ispec_CPML)];
  A11 = pml_convolution_coef_strain[INDEX5(18,NGLLX,NGLLX,NGLLX,5,I,J,K,ispec_CPML)];
  A12 = pml_convolution_coef_strain[INDEX5(18,NGLLX,NGLLX,NGLLX,6,I,J,K,ispec_CPML)];
  A13 = pml_convolution_coef_strain[INDEX5(18,NGLLX,NGLLX,NGLLX,7,I,J,K,ispec_CPML)];

  A14 = pml_convolution_coef_strain[INDEX5(18,NGLLX,NGLLX,NGLLX,8,I,J,K,ispec_CPML)];
  A15 = pml_convolution_coef_strain[INDEX5(18,NGLLX,NGLLX,NGLLX,9,I,J,K,ispec_CPML)];
  A16 = pml_convolution_coef_strain[INDEX5(18,NGLLX,NGLLX,NGLLX,10,I,J,K,ispec_CPML)];
  A17 = pml_convolution_coef_strain[INDEX5(18,NGLLX,NGLLX,NGLLX,11,I,J,K,ispec_CPML)];

  A18 = pml_convolution_coef_strain[INDEX5(18,NGLLX,NGLLX,NGLLX,12,I,J,K,ispec_CPML)];
  A19 = pml_convolution_coef_strain[INDEX5(18,NGLLX,NGLLX,NGLLX,13,I,J,K,ispec_CPML)];

  A20 = pml_convolution_coef_strain[INDEX5(18,NGLLX,NGLLX,NGLLX,14,I,J,K,ispec_CPML)];
  A21 = pml_convolution_coef_strain[INDEX5(18,NGLLX,NGLLX,NGLLX,15,I,J,K,ispec_CPML)];

  A22 = pml_convolution_coef_strain[INDEX5(18,NGLLX,NGLLX,NGLLX,16,I,J,K,ispec_CPML)];
  A23 = pml_convolution_coef_strain[INDEX5(18,NGLLX,NGLLX,NGLLX,17,I,J,K,ispec_CPML)];

  // derivatives
  duxdxl_x = A6 * PML_dux_dxl + A7 * rmemory_dux_dxl_x[offset_x]
           + A8 * rmemory_dux_dxl_x[offset_y] + A9 * rmemory_dux_dxl_x[offset_z];
  duxdyl_x = A10 * PML_dux_dyl + A11 * rmemory_dux_dyl_x[offset_x]
             + A12 * rmemory_dux_dyl_x[offset_y] + A13 * rmemory_dux_dyl_x[offset_z];
  duxdzl_x = A14 * PML_dux_dzl + A15 * rmemory_dux_dzl_x[offset_x]
           + A16 * rmemory_dux_dzl_x[offset_y] + A17 * rmemory_dux_dzl_x[offset_z];

  duzdzl_x = A20 * PML_duz_dzl + A21 * rmemory_duz_dzl_x[offset];
  duzdxl_x = A20 * PML_duz_dxl + A21 * rmemory_duz_dxl_x[offset];
  duydyl_x = A22 * PML_duy_dyl + A23 * rmemory_duy_dyl_x[offset];
  duydxl_x = A22 * PML_duy_dxl + A23 * rmemory_duy_dxl_x[offset];

  duydxl_y = A6 * PML_duy_dxl + A7 * rmemory_duy_dxl_y[offset_x]
           + A8 * rmemory_duy_dxl_y[offset_y] + A9 * rmemory_duy_dxl_y[offset_z];
  duydyl_y = A10 * PML_duy_dyl + A11 * rmemory_duy_dyl_y[offset_x]
           + A12 * rmemory_duy_dyl_y[offset_y] + A13 * rmemory_duy_dyl_y[offset_z];
  duydzl_y = A14 * PML_duy_dzl + A15 * rmemory_duy_dzl_y[offset_x]
           + A16 * rmemory_duy_dzl_y[offset_y] + A17 * rmemory_duy_dzl_y[offset_z];

  duzdzl_y = A18 * PML_duz_dzl + A19 * rmemory_duz_dzl_y[offset];
  duzdyl_y = A18 * PML_duz_dyl + A19 * rmemory_duz_dyl_y[offset];
  duxdyl_y = A22 * PML_dux_dyl + A23 * rmemory_dux_dyl_y[offset];
  duxdxl_y = A22 * PML_dux_dxl + A23 * rmemory_dux_dxl_y[offset];

  duzdxl_z = A6 * PML_duz_dxl + A7 * rmemory_duz_dxl_z[offset_x]
           + A8 * rmemory_duz_dxl_z[offset_y] + A9 * rmemory_duz_dxl_z[offset_z];
  duzdyl_z = A10 * PML_duz_dyl + A11 * rmemory_duz_dyl_z[offset_x]
           + A12 * rmemory_duz_dyl_z[offset_y] + A13 * rmemory_duz_dyl_z[offset_z];
  duzdzl_z = A14 * PML_duz_dzl + A15 * rmemory_duz_dzl_z[offset_x]
           + A16 * rmemory_duz_dzl_z[offset_y] + A17 * rmemory_duz_dzl_z[offset_z];

  duydzl_z = A18 * PML_duy_dzl + A19 * rmemory_duy_dzl_z[offset];
  duydyl_z = A18 * PML_duy_dyl + A19 * rmemory_duy_dyl_z[offset];
  duxdzl_z = A20 * PML_dux_dzl + A21 * rmemory_dux_dzl_z[offset];
  duxdxl_z = A20 * PML_dux_dxl + A21 * rmemory_dux_dxl_z[offset];

  // elastic parameters
  // local padded index
  offset = working_element*NGLL3_PADDED + tx;

  kappal = d_kappav[offset];
  mul = d_muv[offset];

  lambdalplus2mul = kappal + 4.0f/3.0f * mul;
  lambdal = lambdalplus2mul - 2.0f * mul;

  // compute stress sigma (non-symmetric)
  *sigma_xx = lambdalplus2mul*duxdxl_x + lambdal*duydyl_x + lambdal*duzdzl_x;
  *sigma_yx = mul*duxdyl_x + mul*duydxl_x;
  *sigma_zx = mul*duzdxl_x + mul*duxdzl_x;

  *sigma_xy = mul*duxdyl_y + mul*duydxl_y;
  *sigma_yy = lambdal*duxdxl_y + lambdalplus2mul*duydyl_y + lambdal*duzdzl_y;
  *sigma_zy = mul*duzdyl_y + mul*duydzl_y;

  *sigma_xz = mul*duzdxl_z + mul*duxdzl_z;
  *sigma_yz = mul*duzdyl_z + mul*duydzl_z;
  *sigma_zz = lambdal*duxdxl_z + lambdal*duydyl_z + lambdalplus2mul*duzdzl_z;

  return;
}

/* ----------------------------------------------------------------------------------------------- */

// calculates contribution from each C-PML element to update acceleration

__device__ __forceinline__ void pml_compute_accel_contribution_elastic(int tx,int working_element,
                                                                       int iglob, int I, int J, int K,
                                                                       int ispec_CPML,
                                                                       int NSPEC_CPML,
                                                                       realw jacobianl,
                                                                       realw_const_p pml_convolution_coef_alpha,
                                                                       realw_const_p pml_convolution_coef_abar,
                                                                       realw* rmemory_displ_elastic,
                                                                       realw_const_p PML_displ_new,
                                                                       realw_const_p PML_displ_old,
                                                                       realw_const_p d_rhostore,
                                                                       realw_const_p wgll_cube,
                                                                       realw_const_p d_displ,
                                                                       realw_const_p d_veloc,
                                                                       realw* accel_elastic_CPML){

  realw wgllcube;
  realw A_1,A_2,A_3,A_4,A_5;
  realw coef0_x,coef1_x,coef2_x;
  realw coef0_y,coef1_y,coef2_y;
  realw coef0_z,coef1_z,coef2_z;
  realw rhol;

  int offset_x,offset_y,offset_z;
  int offset_xx,offset_yx,offset_zx;
  int offset_xy,offset_yy,offset_zy;
  int offset_xz,offset_yz,offset_zz;
  int offset;

  // coefficients
  // alpha_x
  coef0_x = pml_convolution_coef_alpha[INDEX5(9,NGLLX,NGLLX,NGLLX,0,I,J,K,ispec_CPML)];
  coef1_x = pml_convolution_coef_alpha[INDEX5(9,NGLLX,NGLLX,NGLLX,1,I,J,K,ispec_CPML)];
  coef2_x = pml_convolution_coef_alpha[INDEX5(9,NGLLX,NGLLX,NGLLX,2,I,J,K,ispec_CPML)];

  // alpha_y
  coef0_y = pml_convolution_coef_alpha[INDEX5(9,NGLLX,NGLLX,NGLLX,3,I,J,K,ispec_CPML)];
  coef1_y = pml_convolution_coef_alpha[INDEX5(9,NGLLX,NGLLX,NGLLX,4,I,J,K,ispec_CPML)];
  coef2_y = pml_convolution_coef_alpha[INDEX5(9,NGLLX,NGLLX,NGLLX,5,I,J,K,ispec_CPML)];

  // alpha_z
  coef0_z = pml_convolution_coef_alpha[INDEX5(9,NGLLX,NGLLX,NGLLX,6,I,J,K,ispec_CPML)];
  coef1_z = pml_convolution_coef_alpha[INDEX5(9,NGLLX,NGLLX,NGLLX,7,I,J,K,ispec_CPML)];
  coef2_z = pml_convolution_coef_alpha[INDEX5(9,NGLLX,NGLLX,NGLLX,8,I,J,K,ispec_CPML)];

  offset_x = INDEX5(NDIM,NGLLX,NGLLX,NGLLX,0,I,J,K,ispec_CPML);
  offset_y = INDEX5(NDIM,NGLLX,NGLLX,NGLLX,1,I,J,K,ispec_CPML);
  offset_z = INDEX5(NDIM,NGLLX,NGLLX,NGLLX,2,I,J,K,ispec_CPML);

  offset_xx = INDEX6(NDIM,NGLLX,NGLLX,NGLLX,NSPEC_CPML,0,I,J,K,ispec_CPML,0);
  offset_yx = INDEX6(NDIM,NGLLX,NGLLX,NGLLX,NSPEC_CPML,1,I,J,K,ispec_CPML,0);
  offset_zx = INDEX6(NDIM,NGLLX,NGLLX,NGLLX,NSPEC_CPML,2,I,J,K,ispec_CPML,0);

  // updates memory variables
  rmemory_displ_elastic[offset_xx] = coef0_x * rmemory_displ_elastic[offset_xx]
                                   + PML_displ_new[offset_x] * coef1_x + PML_displ_old[offset_x] * coef2_x;
  rmemory_displ_elastic[offset_yx] = coef0_x * rmemory_displ_elastic[offset_yx]
                                   + PML_displ_new[offset_y] * coef1_x + PML_displ_old[offset_y] * coef2_x;
  rmemory_displ_elastic[offset_zx] = coef0_x * rmemory_displ_elastic[offset_zx]
                                   + PML_displ_new[offset_z] * coef1_x + PML_displ_old[offset_z] * coef2_x;

  offset_xy = INDEX6(NDIM,NGLLX,NGLLX,NGLLX,NSPEC_CPML,0,I,J,K,ispec_CPML,1);
  offset_yy = INDEX6(NDIM,NGLLX,NGLLX,NGLLX,NSPEC_CPML,1,I,J,K,ispec_CPML,1);
  offset_zy = INDEX6(NDIM,NGLLX,NGLLX,NGLLX,NSPEC_CPML,2,I,J,K,ispec_CPML,1);

  rmemory_displ_elastic[offset_xy] = coef0_y * rmemory_displ_elastic[offset_xy]
                                   + PML_displ_new[offset_x] * coef1_y + PML_displ_old[offset_x] * coef2_y;
  rmemory_displ_elastic[offset_yy] = coef0_y * rmemory_displ_elastic[offset_yy]
                                   + PML_displ_new[offset_y] * coef1_y + PML_displ_old[offset_y] * coef2_y;
  rmemory_displ_elastic[offset_zy] = coef0_y * rmemory_displ_elastic[offset_zy]
                                   + PML_displ_new[offset_z] * coef1_y + PML_displ_old[offset_z] * coef2_y;

  offset_xz = INDEX6(NDIM,NGLLX,NGLLX,NGLLX,NSPEC_CPML,0,I,J,K,ispec_CPML,2);
  offset_yz = INDEX6(NDIM,NGLLX,NGLLX,NGLLX,NSPEC_CPML,1,I,J,K,ispec_CPML,2);
  offset_zz = INDEX6(NDIM,NGLLX,NGLLX,NGLLX,NSPEC_CPML,2,I,J,K,ispec_CPML,2);

  rmemory_displ_elastic[offset_xz] = coef0_z * rmemory_displ_elastic[offset_xz]
                                   + PML_displ_new[offset_x] * coef1_z + PML_displ_old[offset_x] * coef2_z;
  rmemory_displ_elastic[offset_yz] = coef0_z * rmemory_displ_elastic[offset_yz]
                                   + PML_displ_new[offset_y] * coef1_z + PML_displ_old[offset_y] * coef2_z;
  rmemory_displ_elastic[offset_zz] = coef0_z * rmemory_displ_elastic[offset_zz]
                                   + PML_displ_new[offset_z] * coef1_z + PML_displ_old[offset_z] * coef2_z;

  // local padded index
  offset = working_element*NGLL3_PADDED + tx;

  rhol = get_global_cr( &d_rhostore[offset] );
  wgllcube = wgll_cube[tx];

  // PML coefficient values
  A_1 = pml_convolution_coef_abar[INDEX5(5,NGLLX,NGLLX,NGLLX,0,I,J,K,ispec_CPML)];
  A_2 = pml_convolution_coef_abar[INDEX5(5,NGLLX,NGLLX,NGLLX,1,I,J,K,ispec_CPML)];
  A_3 = pml_convolution_coef_abar[INDEX5(5,NGLLX,NGLLX,NGLLX,2,I,J,K,ispec_CPML)];
  A_4 = pml_convolution_coef_abar[INDEX5(5,NGLLX,NGLLX,NGLLX,3,I,J,K,ispec_CPML)];
  A_5 = pml_convolution_coef_abar[INDEX5(5,NGLLX,NGLLX,NGLLX,4,I,J,K,ispec_CPML)];

  // updates PML acceleration
  accel_elastic_CPML[0] =  wgllcube * rhol * jacobianl *
       ( A_1 * d_veloc[iglob*3] + A_2 * d_displ[iglob*3] +
         A_3 * rmemory_displ_elastic[offset_xx] +
         A_4 * rmemory_displ_elastic[offset_xy] +
         A_5 * rmemory_displ_elastic[offset_xz]
       );

  accel_elastic_CPML[1] =  wgllcube * rhol * jacobianl *
       ( A_1 * d_veloc[iglob*3+1] + A_2 * d_displ[iglob*3+1] +
         A_3 * rmemory_displ_elastic[offset_yx] +
         A_4 * rmemory_displ_elastic[offset_yy] +
         A_5 * rmemory_displ_elastic[offset_yz]
       );

  accel_elastic_CPML[2] =  wgllcube * rhol * jacobianl *
       ( A_1 * d_veloc[iglob*3+2] + A_2 * d_displ[iglob*3+2] +
         A_3 * rmemory_displ_elastic[offset_zx] +
         A_4 * rmemory_displ_elastic[offset_zy] +
         A_5 * rmemory_displ_elastic[offset_zz]
       );

  return;
}

/* ----------------------------------------------------------------------------------------------- */

// loads displacement into shared memory for element

template<int FORWARD_OR_ADJOINT>
__device__  __forceinline__ void load_shared_memory_displ(const int* tx, const int* iglob,
                                                          realw_const_p d_displ,
                                                          realw* sh_displx,
                                                          realw* sh_disply,
                                                          realw* sh_displz){

  // copy from global memory to shared memory
  // each thread writes one of the NGLL^3 = 125 data points
#ifdef USE_TEXTURES_FIELDS
  sh_displx[(*tx)] = texfetch_displ<FORWARD_OR_ADJOINT>((*iglob)*3);
  sh_disply[(*tx)] = texfetch_displ<FORWARD_OR_ADJOINT>((*iglob)*3 + 1);
  sh_displz[(*tx)] = texfetch_displ<FORWARD_OR_ADJOINT>((*iglob)*3 + 2);
#else
  // changing iglob indexing to match fortran row changes fast style
  sh_displx[(*tx)] = d_displ[(*iglob)*3];
  sh_disply[(*tx)] = d_displ[(*iglob)*3 + 1];
  sh_displz[(*tx)] = d_displ[(*iglob)*3 + 2];
#endif
}

/* ----------------------------------------------------------------------------------------------- */

// loads displacement into shared memory for element

__device__  __forceinline__ void load_shared_memory_PML_displ(const int* tx, const int* offset_x, const int* offset_y, const int* offset_z,
                                                              realw_const_p d_PML_displ,
                                                              realw* sh_pml_displx,
                                                              realw* sh_pml_disply,
                                                              realw* sh_pml_displz){

  // copy from global memory to shared memory
  // each thread writes one of the NGLL^3 = 125 data points
  sh_pml_displx[(*tx)] = d_PML_displ[(*offset_x)];
  sh_pml_disply[(*tx)] = d_PML_displ[(*offset_y)];
  sh_pml_displz[(*tx)] = d_PML_displ[(*offset_z)];
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

// computes the spatial derivatives

__device__  __forceinline__ void get_spatial_derivatives_pml(realw* xixl,realw* xiyl,realw* xizl,
                                                             realw* etaxl,realw* etayl,realw* etazl,
                                                             realw* gammaxl,realw* gammayl,realw* gammazl,
                                                             realw* jacobianl,
                                                             int I,int J,int K,int tx,
                                                             realw* tempx1l,realw* tempy1l,realw* tempz1l,
                                                             realw* tempx2l,realw* tempy2l,realw* tempz2l,
                                                             realw* tempx3l,realw* tempy3l,realw* tempz3l,
                                                             realw* sh_tempx,realw* sh_tempy,realw* sh_tempz,
                                                             realw* sh_hprime_xx,
                                                             realw* duxdxl,realw* duxdyl,realw* duxdzl,
                                                             realw* duydxl,realw* duydyl,realw* duydzl,
                                                             realw* duzdxl,realw* duzdyl,realw* duzdzl,
                                                             realw_const_p d_xix,realw_const_p d_xiy,realw_const_p d_xiz,
                                                             realw_const_p d_etax,realw_const_p d_etay,realw_const_p d_etaz,
                                                             realw_const_p d_gammax,realw_const_p d_gammay,realw_const_p d_gammaz,
                                                             int ispec_irreg, realw xix_regular, realw jacobian_regular,
                                                             int ipass){

  // computes first matrix products

  // determines jacobian
  // this is the only difference to the get_spatial_derivatives() routines in the other kernel_2 implementations.
  // we need `jacobianl` to be valid also for regular elements
  // as it will be needed in the pml_compute_accel_contribution_elastic() routine.
  if (ipass == 0){
    if (ispec_irreg >= 0){
      // irregular_element
      // local padded index
      int offset = ispec_irreg*NGLL3_PADDED + tx;

      *xixl = get_global_cr(&d_xix[offset]);
      *xiyl = get_global_cr(&d_xiy[offset]);
      *xizl = get_global_cr(&d_xiz[offset]);
      *etaxl = get_global_cr(&d_etax[offset]);
      *etayl = get_global_cr(&d_etay[offset]);
      *etazl = get_global_cr(&d_etaz[offset]);
      *gammaxl = get_global_cr(&d_gammax[offset]);
      *gammayl = get_global_cr(&d_gammay[offset]);
      *gammazl = get_global_cr(&d_gammaz[offset]);

      *jacobianl = 1.f / ((*xixl)*((*etayl)*(*gammazl)-(*etazl)*(*gammayl))
                        -(*xiyl)*((*etaxl)*(*gammazl)-(*etazl)*(*gammaxl))
                        +(*xizl)*((*etaxl)*(*gammayl)-(*etayl)*(*gammaxl)));
    }else{
      // regular elements
      *jacobianl = jacobian_regular;
    }
  }

  // 1. cut-plane
  sum_hprime_xi(I,J,K,tempx1l,tempy1l,tempz1l,sh_tempx,sh_tempy,sh_tempz,sh_hprime_xx);
  // 2. cut-plane
  sum_hprime_eta(I,J,K,tempx2l,tempy2l,tempz2l,sh_tempx,sh_tempy,sh_tempz,sh_hprime_xx);
  // 3. cut-plane
  sum_hprime_gamma(I,J,K,tempx3l,tempy3l,tempz3l,sh_tempx,sh_tempy,sh_tempz,sh_hprime_xx);

  // synchronize all the threads (one thread for each of the NGLL grid points of the
  // current spectral element) because we need the whole element to be ready in order
  // to be able to compute the matrix products along cut planes of the 3D element below
  __syncthreads();

  // compute derivatives of ux, uy and uz with respect to x, y and z
  if (ispec_irreg >= 0 ){
    // irregular_element
    (*duxdxl) = (*xixl)*(*tempx1l) + (*etaxl)*(*tempx2l) + (*gammaxl)*(*tempx3l);
    (*duxdyl) = (*xiyl)*(*tempx1l) + (*etayl)*(*tempx2l) + (*gammayl)*(*tempx3l);
    (*duxdzl) = (*xizl)*(*tempx1l) + (*etazl)*(*tempx2l) + (*gammazl)*(*tempx3l);

    (*duydxl) = (*xixl)*(*tempy1l) + (*etaxl)*(*tempy2l) + (*gammaxl)*(*tempy3l);
    (*duydyl) = (*xiyl)*(*tempy1l) + (*etayl)*(*tempy2l) + (*gammayl)*(*tempy3l);
    (*duydzl) = (*xizl)*(*tempy1l) + (*etazl)*(*tempy2l) + (*gammazl)*(*tempy3l);

    (*duzdxl) = (*xixl)*(*tempz1l) + (*etaxl)*(*tempz2l) + (*gammaxl)*(*tempz3l);
    (*duzdyl) = (*xiyl)*(*tempz1l) + (*etayl)*(*tempz2l) + (*gammayl)*(*tempz3l);
    (*duzdzl) = (*xizl)*(*tempz1l) + (*etazl)*(*tempz2l) + (*gammazl)*(*tempz3l);
  }
  else{
    // regular element
    (*duxdxl) = xix_regular*(*tempx1l);
    (*duxdyl) = xix_regular*(*tempx2l);
    (*duxdzl) = xix_regular*(*tempx3l);

    (*duydxl) = xix_regular*(*tempy1l);
    (*duydyl) = xix_regular*(*tempy2l);
    (*duydzl) = xix_regular*(*tempy3l);

    (*duzdxl) = xix_regular*(*tempz1l);
    (*duzdyl) = xix_regular*(*tempz2l);
    (*duzdzl) = xix_regular*(*tempz3l);
  }
}

/* ----------------------------------------------------------------------------------------------- */

// computes dot product between tensor and derivatives

__device__  __forceinline__ void get_dot_product(realw jacobianl,
                                                 realw sigma_xx,realw sigma_xy,realw sigma_yx,
                                                 realw sigma_xz,realw sigma_zx,realw sigma_yy,
                                                 realw sigma_yz,realw sigma_zy,realw sigma_zz,
                                                 realw Dxl,realw Dyl,realw Dzl,
                                                 realw* sh_tempx,realw* sh_tempy,realw* sh_tempz,
                                                 int tx,
                                                 int ispec_irreg,realw xix_regular,realw jacobian_regular,
                                                 int component){

  // fills shared memory arrays
  if (threadIdx.x < NGLL3) {
    if (ispec_irreg >= 0){
      //irregular element
      sh_tempx[tx] = jacobianl * (sigma_xx*Dxl + sigma_yx*Dyl + sigma_zx*Dzl); // sh_tempx1
      sh_tempy[tx] = jacobianl * (sigma_xy*Dxl + sigma_yy*Dyl + sigma_zy*Dzl); // sh_tempy1
      sh_tempz[tx] = jacobianl * (sigma_xz*Dxl + sigma_yz*Dyl + sigma_zz*Dzl); // sh_tempz1
    }
    else if (component == 1){
      sh_tempx[tx] = jacobian_regular * (sigma_xx*xix_regular); // sh_tempx1
      sh_tempy[tx] = jacobian_regular * (sigma_xy*xix_regular); // sh_tempy1
      sh_tempz[tx] = jacobian_regular * (sigma_xz*xix_regular); // sh_tempz1
    }
    else if (component == 2){
      sh_tempx[tx] = jacobian_regular * (sigma_yx*xix_regular); // sh_tempx1
      sh_tempy[tx] = jacobian_regular * (sigma_yy*xix_regular); // sh_tempy1
      sh_tempz[tx] = jacobian_regular * (sigma_yz*xix_regular); // sh_tempz1

    }else{
      sh_tempx[tx] = jacobian_regular * (sigma_zx*xix_regular); // sh_tempx1
      sh_tempy[tx] = jacobian_regular * (sigma_zy*xix_regular); // sh_tempy1
      sh_tempz[tx] = jacobian_regular * (sigma_zz*xix_regular); // sh_tempz1
    }
  }
  __syncthreads();
}


/* ----------------------------------------------------------------------------------------------- */

// compute forces kernel for PML elements

/* ----------------------------------------------------------------------------------------------- */

__global__ void
#ifdef USE_LAUNCH_BOUNDS
// adds compiler specification
__launch_bounds__(NGLL3_PADDED,LAUNCH_MIN_BLOCKS)
#endif
// main kernel
pml_kernel_2_impl(int nb_blocks_to_compute,
                  const int* d_ibool,
                  const int* d_phase_ispec_inner_elastic,const int num_phase_ispec_elastic,
                  const int d_iphase,
                  const int* d_irregular_element_number,
                  realw_const_p d_displ,
                  realw_const_p d_veloc,
                  realw_p d_accel,
                  realw_const_p d_xix,realw_const_p d_xiy,realw_const_p d_xiz,
                  realw_const_p d_etax,realw_const_p d_etay,realw_const_p d_etaz,
                  realw_const_p d_gammax,realw_const_p d_gammay,realw_const_p d_gammaz,
                  const realw xix_regular,const realw jacobian_regular,
                  realw_const_p d_hprime_xx,
                  realw_const_p d_hprimewgll_xx,
                  realw_const_p d_wgllwgll_xy,realw_const_p d_wgllwgll_xz,realw_const_p d_wgllwgll_yz,
                  realw_const_p d_kappav,realw_const_p d_muv,
                  const int COMPUTE_AND_STORE_STRAIN,
                  realw_p epsilondev_xx,realw_p epsilondev_yy,realw_p epsilondev_xy,
                  realw_p epsilondev_xz,realw_p epsilondev_yz,
                  realw_p epsilon_trace_over_3,
                  const int SIMULATION_TYPE,
                  realw_const_p d_rhostore,
                  realw_const_p wgll_cube,
                  const int NSPEC_CPML,
                  const int* d_is_CPML,
                  const int* d_spec_to_CPML,
                  realw_const_p d_PML_displ_new,
                  realw_const_p d_PML_displ_old,
                  realw_p d_rmemory_displ_elastic,
                  realw_p d_rmemory_dux_dxl_x,
                  realw_p d_rmemory_duy_dxl_y,
                  realw_p d_rmemory_duz_dxl_z,
                  realw_p d_rmemory_dux_dyl_x,
                  realw_p d_rmemory_duy_dyl_y,
                  realw_p d_rmemory_duz_dyl_z,
                  realw_p d_rmemory_dux_dzl_x,
                  realw_p d_rmemory_duy_dzl_y,
                  realw_p d_rmemory_duz_dzl_z,
                  realw_p d_rmemory_dux_dxl_y,
                  realw_p d_rmemory_dux_dxl_z,
                  realw_p d_rmemory_duy_dxl_x,
                  realw_p d_rmemory_duz_dxl_x,
                  realw_p d_rmemory_dux_dyl_y,
                  realw_p d_rmemory_duy_dyl_x,
                  realw_p d_rmemory_duy_dyl_z,
                  realw_p d_rmemory_duz_dyl_y,
                  realw_p d_rmemory_dux_dzl_z,
                  realw_p d_rmemory_duy_dzl_z,
                  realw_p d_rmemory_duz_dzl_x,
                  realw_p d_rmemory_duz_dzl_y,
                  realw_const_p pml_convolution_coef_alpha,
                  realw_const_p pml_convolution_coef_beta,
                  realw_const_p pml_convolution_coef_strain,
                  realw_const_p pml_convolution_coef_abar){

// elastic compute kernel for PML elements
//
// holds for: PML_CONDITIONS = .true.
//            COMPUTE_AND_STORE_STRAIN = .false. or .true.

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
  int working_element,ispec_irreg;

  realw tempx1l,tempx2l,tempx3l,tempy1l,tempy2l,tempy3l,tempz1l,tempz2l,tempz3l;
  realw xixl,xiyl,xizl,etaxl,etayl,etazl,gammaxl,gammayl,gammazl,jacobianl;
  realw duxdxl,duxdyl,duxdzl,duydxl,duydyl,duydzl,duzdxl,duzdyl,duzdzl;

  realw sigma_xx,sigma_yy,sigma_zz,sigma_xy,sigma_xz,sigma_yz;
  realw sigma_yx,sigma_zx,sigma_zy;

  realw fac1,fac2,fac3;
  realw sum_terms1,sum_terms2,sum_terms3;

  __shared__ realw sh_hprime_xx[NGLL2];
  __shared__ realw sh_hprimewgll_xx[NGLL2];

  __shared__ realw sh_tempx[NGLL3];
  __shared__ realw sh_tempy[NGLL3];
  __shared__ realw sh_tempz[NGLL3];

  __shared__ realw s_dummyx_loc[NGLL3];
  __shared__ realw s_dummyy_loc[NGLL3];
  __shared__ realw s_dummyz_loc[NGLL3];

  // PML
  __shared__ realw s_dummyx_loc_new[NGLL3];
  __shared__ realw s_dummyy_loc_new[NGLL3];
  __shared__ realw s_dummyz_loc_new[NGLL3];

  __shared__ realw s_dummyx_loc_old[NGLL3];
  __shared__ realw s_dummyy_loc_old[NGLL3];
  __shared__ realw s_dummyz_loc_old[NGLL3];

  // derivatives of ux, uy and uz with respect to x, y and z
  // in PML_du* computation displ at "n" time step is used
  // in PML_du*_old computation we replace displ with displ_old
  realw PML_dux_dxl_old,PML_dux_dyl_old,PML_dux_dzl_old;
  realw PML_duy_dxl_old,PML_duy_dyl_old,PML_duy_dzl_old;
  realw PML_duz_dxl_old,PML_duz_dyl_old,PML_duz_dzl_old;
  // in PML_du*_new computation we replace displ with displ_new
  realw PML_dux_dxl_new,PML_dux_dyl_new,PML_dux_dzl_new;
  realw PML_duy_dxl_new,PML_duy_dyl_new,PML_duy_dzl_new;
  realw PML_duz_dxl_new,PML_duz_dyl_new,PML_duz_dzl_new;
  // PML accel contribution
  realw accel_elastic_CPML[NDIM];

  int ispec_CPML;
  int offset_x,offset_y,offset_z;

  // iphase-1 and working_element-1 for Fortran->C array conventions
  working_element = d_phase_ispec_inner_elastic[bx + num_phase_ispec_elastic*(d_iphase-1)]-1;

  // only PML elements will be computed
  if(d_is_CPML[working_element] == 0) return;

  // loads hprime's into shared memory
  if (tx < NGLL2) {
    // copy hprime from global memory to shared memory
    load_shared_memory_hprime(&tx,d_hprime_xx,sh_hprime_xx);
    // copy hprime from global memory to shared memory
    load_shared_memory_hprimewgll(&tx,d_hprimewgll_xx,sh_hprimewgll_xx);
  }

  ispec_irreg = d_irregular_element_number[working_element] - 1;

  // local padded index
  offset = working_element*NGLL3_PADDED + tx;

  // global index
  iglob = d_ibool[offset] - 1;

  // additional old and new arrays
  ispec_CPML = d_spec_to_CPML[working_element] - 1;

  // copy displacement from global memory to shared memory
  // each thread writes one of the NGLL^3 = 125 data points
  if (threadIdx.x < NGLL3 ){
    // displ
    load_shared_memory_displ<1>(&tx,&iglob,d_displ,s_dummyx_loc,s_dummyy_loc,s_dummyz_loc);

    // index in PML_displ arrays
    offset_x = INDEX5(NDIM,NGLLX,NGLLX,NGLLX,0,I,J,K,ispec_CPML);
    offset_y = INDEX5(NDIM,NGLLX,NGLLX,NGLLX,1,I,J,K,ispec_CPML);
    offset_z = INDEX5(NDIM,NGLLX,NGLLX,NGLLX,2,I,J,K,ispec_CPML);

    // PML_displ_old
    load_shared_memory_PML_displ(&tx,&offset_x,&offset_y,&offset_z,d_PML_displ_old,s_dummyx_loc_old,s_dummyy_loc_old,s_dummyz_loc_old);

    // PML_displ_old
    load_shared_memory_PML_displ(&tx,&offset_x,&offset_y,&offset_z,d_PML_displ_new,s_dummyx_loc_new,s_dummyy_loc_new,s_dummyz_loc_new);
  }

  // synchronize all the threads (one thread for each of the NGLL grid points of the
  // current spectral element) because we need the whole element to be ready in order
  // to be able to compute the matrix products along cut planes of the 3D element below
  __syncthreads();

  // grad(u)
  // computes the spatial derivatives duxdxl
  get_spatial_derivatives_pml(&xixl,&xiyl,&xizl,&etaxl,&etayl,&etazl,
                              &gammaxl,&gammayl,&gammazl,&jacobianl,
                              I,J,K,tx,
                              &tempx1l,&tempy1l,&tempz1l,
                              &tempx2l,&tempy2l,&tempz2l,
                              &tempx3l,&tempy3l,&tempz3l,
                              s_dummyx_loc,s_dummyy_loc,s_dummyz_loc,
                              sh_hprime_xx,
                              &duxdxl,&duxdyl,&duxdzl,&duydxl,&duydyl,&duydzl,&duzdxl,&duzdyl,&duzdzl,
                              d_xix,d_xiy,d_xiz,d_etax,d_etay,d_etaz,d_gammax,d_gammay,d_gammaz,
                              ispec_irreg,xix_regular,jacobian_regular,0); // ipass==0 -> gets xixl,.. and calculates jacobianl; needs to be done only once

  // computes deviatoric strain for kernel calculations
  // (maybe not really needed, but will keep for now based on a "pure" elastic element contribution)
  if (COMPUTE_AND_STORE_STRAIN) {
    // non-attenuation case
    if (threadIdx.x < NGLL3) {
      // local storage: stresses at this current time step
      realw templ = 0.33333333333333333333f * (duxdxl + duydyl + duzdzl); // 1./3. = 0.33333
      // kernel simulations
      if (SIMULATION_TYPE == 3){ epsilon_trace_over_3[tx + working_element*NGLL3] = templ; }
      // fortran: epsilondev_xx(:,:,:,ispec) = epsilondev_xx_loc(:,:,:)
      epsilondev_xx[tx + working_element*NGLL3] = duxdxl - templ; // epsilondev_xx_loc;
      epsilondev_yy[tx + working_element*NGLL3] = duydyl - templ; // epsilondev_yy_loc;
      epsilondev_xy[tx + working_element*NGLL3] = 0.5f * (duxdyl + duydxl); // epsilondev_xy_loc;
      epsilondev_xz[tx + working_element*NGLL3] = 0.5f * (duzdxl + duxdzl); // epsilondev_xz_loc;
      epsilondev_yz[tx + working_element*NGLL3] = 0.5f * (duzdyl + duydzl); //epsilondev_yz_loc;
      //epsilondev_trace[tx + working_element*NGLL3] = 3.0f * templ ! not needed for PML elements
    } // threadIdx.x
  }

  // PML_displ_old
  get_spatial_derivatives_pml(&xixl,&xiyl,&xizl,&etaxl,&etayl,&etazl,
                              &gammaxl,&gammayl,&gammazl,&jacobianl,
                              I,J,K,tx,
                              &tempx1l,&tempy1l,&tempz1l,
                              &tempx2l,&tempy2l,&tempz2l,
                              &tempx3l,&tempy3l,&tempz3l,
                              s_dummyx_loc_old,s_dummyy_loc_old,s_dummyz_loc_old,
                              sh_hprime_xx,
                              &PML_dux_dxl_old,&PML_dux_dyl_old,&PML_dux_dzl_old,
                              &PML_duy_dxl_old,&PML_duy_dyl_old,&PML_duy_dzl_old,
                              &PML_duz_dxl_old,&PML_duz_dyl_old,&PML_duz_dzl_old,
                              d_xix,d_xiy,d_xiz,d_etax,d_etay,d_etaz,d_gammax,d_gammay,d_gammaz,
                              ispec_irreg,xix_regular,jacobian_regular,1);

  // PML_displ_new
  get_spatial_derivatives_pml(&xixl,&xiyl,&xizl,&etaxl,&etayl,&etazl,
                              &gammaxl,&gammayl,&gammazl,&jacobianl,
                              I,J,K,tx,
                              &tempx1l,&tempy1l,&tempz1l,
                              &tempx2l,&tempy2l,&tempz2l,
                              &tempx3l,&tempy3l,&tempz3l,
                              s_dummyx_loc_new,s_dummyy_loc_new,s_dummyz_loc_new,
                              sh_hprime_xx,
                              &PML_dux_dxl_new,&PML_dux_dyl_new,&PML_dux_dzl_new,
                              &PML_duy_dxl_new,&PML_duy_dyl_new,&PML_duy_dzl_new,
                              &PML_duz_dxl_new,&PML_duz_dyl_new,&PML_duz_dzl_new,
                              d_xix,d_xiy,d_xiz,d_etax,d_etay,d_etaz,d_gammax,d_gammay,d_gammaz,
                              ispec_irreg,xix_regular,jacobian_regular,1);

  // updates memory variables for convolution & stress
  pml_compute_memory_variables_elastic(tx,working_element,
                                       I,J,K,
                                       ispec_CPML,NSPEC_CPML,
                                       pml_convolution_coef_alpha,
                                       pml_convolution_coef_beta,
                                       pml_convolution_coef_strain,
                                       d_rmemory_dux_dxl_x,
                                       d_rmemory_duy_dxl_y,
                                       d_rmemory_duz_dxl_z,
                                       PML_dux_dxl_new,PML_dux_dxl_old,
                                       PML_duy_dxl_new,PML_duy_dxl_old,
                                       PML_duz_dxl_new,PML_duz_dxl_old,
                                       d_rmemory_dux_dyl_x,
                                       d_rmemory_duy_dyl_y,
                                       d_rmemory_duz_dyl_z,
                                       PML_dux_dyl_new,PML_dux_dyl_old,
                                       PML_duy_dyl_new,PML_duy_dyl_old,
                                       PML_duz_dyl_new,PML_duz_dyl_old,
                                       d_rmemory_dux_dzl_x,
                                       d_rmemory_duy_dzl_y,
                                       d_rmemory_duz_dzl_z,
                                       PML_dux_dzl_new,PML_dux_dzl_old,
                                       PML_duy_dzl_new,PML_duy_dzl_old,
                                       PML_duz_dzl_new,PML_duz_dzl_old,
                                       d_rmemory_dux_dxl_y,
                                       d_rmemory_dux_dxl_z,
                                       d_rmemory_duy_dxl_x,
                                       d_rmemory_duz_dxl_x,
                                       d_rmemory_dux_dyl_y,
                                       d_rmemory_duy_dyl_x,
                                       d_rmemory_duy_dyl_z,
                                       d_rmemory_duz_dyl_y,
                                       d_rmemory_dux_dzl_z,
                                       d_rmemory_duy_dzl_z,
                                       d_rmemory_duz_dzl_x,
                                       d_rmemory_duz_dzl_y,
                                       duxdxl,duxdyl,duxdzl,
                                       duydxl,duydyl,duydzl,
                                       duzdxl,duzdyl,duzdzl,
                                       d_kappav,
                                       d_muv,
                                       &sigma_xx,&sigma_yx,&sigma_zx,
                                       &sigma_xy,&sigma_yy,&sigma_zy,
                                       &sigma_xz,&sigma_yz,&sigma_zz);

  // calculates contribution from each C-PML element to update acceleration
  pml_compute_accel_contribution_elastic(tx, working_element,
                                         iglob, I, J, K,
                                         ispec_CPML,
                                         NSPEC_CPML,
                                         jacobianl,
                                         pml_convolution_coef_alpha,
                                         pml_convolution_coef_abar,
                                         d_rmemory_displ_elastic,
                                         d_PML_displ_new,
                                         d_PML_displ_old,
                                         d_rhostore,
                                         wgll_cube,
                                         d_displ,
                                         d_veloc,
                                         accel_elastic_CPML);

  // second double-loop over GLL to compute all the terms
  // form dot product with test vector, symmetric form
  // 1. cut-plane xi
  __syncthreads();
  get_dot_product(jacobianl,sigma_xx,sigma_xy,sigma_yx,sigma_xz,sigma_zx,sigma_yy,sigma_yz,sigma_zy,sigma_zz,
                  xixl,xiyl,xizl,
                  sh_tempx,sh_tempy,sh_tempz,
                  tx,ispec_irreg,xix_regular,jacobian_regular,1);

  sum_hprimewgll_xi(I,J,K,&tempx1l,&tempy1l,&tempz1l,
                    sh_tempx,sh_tempy,sh_tempz,sh_hprimewgll_xx);

  // 2. cut-plane eta
  __syncthreads();
  get_dot_product(jacobianl,sigma_xx,sigma_xy,sigma_yx,sigma_xz,sigma_zx,sigma_yy,sigma_yz,sigma_zy,sigma_zz,
                  etaxl,etayl,etazl,
                  sh_tempx,sh_tempy,sh_tempz,
                  tx,ispec_irreg,xix_regular,jacobian_regular,2);

  sum_hprimewgll_eta(I,J,K,&tempx2l,&tempy2l,&tempz2l,
                     sh_tempx,sh_tempy,sh_tempz,sh_hprimewgll_xx);

  // 3. cut-plane gamma
  __syncthreads();
  get_dot_product(jacobianl,sigma_xx,sigma_xy,sigma_yx,sigma_xz,sigma_zx,sigma_yy,sigma_yz,sigma_zy,sigma_zz,
                  gammaxl,gammayl,gammazl,
                  sh_tempx,sh_tempy,sh_tempz,
                  tx,ispec_irreg,xix_regular,jacobian_regular,3);

  sum_hprimewgll_gamma(I,J,K,&tempx3l,&tempy3l,&tempz3l,
                       sh_tempx,sh_tempy,sh_tempz,sh_hprimewgll_xx);

  // gets double weights
  fac1 = d_wgllwgll_yz[K*NGLLX+J];
  fac2 = d_wgllwgll_xz[K*NGLLX+I];
  fac3 = d_wgllwgll_xy[J*NGLLX+I];

  sum_terms1 = - (fac1*tempx1l + fac2*tempx2l + fac3*tempx3l + accel_elastic_CPML[0]);
  sum_terms2 = - (fac1*tempy1l + fac2*tempy2l + fac3*tempy3l + accel_elastic_CPML[1]);
  sum_terms3 = - (fac1*tempz1l + fac2*tempz2l + fac3*tempz3l + accel_elastic_CPML[2]);

  // assembles acceleration array
  if (threadIdx.x < NGLL3) {
    atomicAdd(&d_accel[iglob*3], sum_terms1);
    atomicAdd(&d_accel[iglob*3+1], sum_terms2);
    atomicAdd(&d_accel[iglob*3+2], sum_terms3);
  } // threadIdx.x

}  // pml_kernel_2_impl

