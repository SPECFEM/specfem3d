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

#ifndef COMPUTE_GRADIENT_GPU_H
#define COMPUTE_GRADIENT_GPU_H

// needed in compute_kernels_acoustic_kernel.cu and compute_kernels_hess_ac_cudakernel.cu

__device__ __forceinline__ void compute_gradient_kernel(int ijk,
                                                        int ispec,int ispec_irreg,
                                                        field* scalar_field,
                                                        field* vector_field_loc,
                                                        realw* d_hprime_xx,
                                                        realw* d_xix,realw* d_xiy,realw* d_xiz,
                                                        realw* d_etax,realw* d_etay,realw* d_etaz,
                                                        realw* d_gammax,realw* d_gammay,realw* d_gammaz,
                                                        realw rhol,realw xix_regular,
                                                        int gravity) {

  field temp1l,temp2l,temp3l;
  realw hp1,hp2,hp3;
  realw xixl,xiyl,xizl,etaxl,etayl,etazl,gammaxl,gammayl,gammazl;
  realw rho_invl;

  int K = (ijk/NGLL2);
  int J = ((ijk-K*NGLL2)/NGLLX);
  int I = (ijk-K*NGLL2-J*NGLLX);

  // derivative along x
  temp1l = Make_field(0.f);
  for(int l=0; l<NGLLX; l++){
    hp1 = d_hprime_xx[l*NGLLX+I];
    temp1l += scalar_field[K*NGLL2+J*NGLLX+l]*hp1;
  }

  // derivative along y
  temp2l = Make_field(0.f);
  for(int l=0; l<NGLLX; l++){
    // assumes hprime_xx == hprime_yy
    hp2 = d_hprime_xx[l*NGLLX+J];
    temp2l += scalar_field[K*NGLL2+l*NGLLX+I]*hp2;
  }

  // derivative along z
  temp3l = Make_field(0.f);
  for(int l=0; l<NGLLX; l++){
    // assumes hprime_xx == hprime_zz
    hp3 = d_hprime_xx[l*NGLLX+K];
    temp3l += scalar_field[l*NGLL2+J*NGLLX+I]*hp3;
  }

  //if (gravity){
  //  // daniel: TODO - check gravity case here
  //  rho_invl = 1.0f / rhol;
  //}else{
    rho_invl = 1.0f / rhol;
  //}

  if (ispec_irreg >= 0){
    int offset = ispec_irreg * NGLL3_PADDED + ijk;
    xixl = d_xix[offset];
    xiyl = d_xiy[offset];
    xizl = d_xiz[offset];
    etaxl = d_etax[offset];
    etayl = d_etay[offset];
    etazl = d_etaz[offset];
    gammaxl = d_gammax[offset];
    gammayl = d_gammay[offset];
    gammazl = d_gammaz[offset];

    // derivatives of acoustic scalar potential field on GLL points
    vector_field_loc[0] = (temp1l*xixl + temp2l*etaxl + temp3l*gammaxl) * rho_invl;
    vector_field_loc[1] = (temp1l*xiyl + temp2l*etayl + temp3l*gammayl) * rho_invl;
    vector_field_loc[2] = (temp1l*xizl + temp2l*etazl + temp3l*gammazl) * rho_invl;

  }else{
    // derivatives of acoustic scalar potential field on GLL points
    vector_field_loc[0] = temp1l * xix_regular * rho_invl;
    vector_field_loc[1] = temp2l * xix_regular * rho_invl;
    vector_field_loc[2] = temp3l * xix_regular * rho_invl;
  }
}


#endif   // COMPUTE_GRADIENT_GPU_H
