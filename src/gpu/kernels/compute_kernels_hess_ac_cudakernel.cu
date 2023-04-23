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

// includes device function compute_gradient_kernel()
#include "compute_gradient_kernel.h"


__global__ void compute_kernels_hess_ac_cudakernel(int* ispec_is_acoustic,
                                                   int* d_ibool,
                                                   field* potential_dot_dot_acoustic,
                                                   field* b_potential_dot_dot_acoustic,
                                                   field* b_potential_dot_acoustic,
                                                   realw* rhostore,
                                                   realw* kappastore,
                                                   realw* d_hprime_xx,
                                                   int* d_irregular_element_number,
                                                   realw* d_xix,realw* d_xiy,realw* d_xiz,
                                                   realw* d_etax,realw* d_etay,realw* d_etaz,
                                                   realw* d_gammax,realw* d_gammay,realw* d_gammaz,
                                                   realw xix_regular,
                                                   realw* hess_kl,
                                                   realw* hess_rho_ac_kl,
                                                   realw* hess_kappa_ac_kl,
                                                   realw deltat,
                                                   int NSPEC_AB,
                                                   int gravity) {

  int ispec = blockIdx.x + blockIdx.y*gridDim.x;
  int ijk = threadIdx.x;
  int ijk_ispec_padded = ijk + NGLL3_PADDED*ispec;
  int ijk_ispec = ijk + NGLL3*ispec;
  int iglob;
  int ispec_irreg = d_irregular_element_number[ispec] - 1;

  // shared memory between all threads within this block
  __shared__ field scalar_field_accel[NGLL3];
  __shared__ field scalar_field_b_accel[NGLL3];
  __shared__ field scalar_field_b_veloc[NGLL3];

  int active = 0;

  // handles case when there is 1 extra block (due to rectangular grid)
  if (ispec < NSPEC_AB) {

    // acoustic elements only
    if (ispec_is_acoustic[ispec] ){
      active = 1;

      // global indices
      iglob = d_ibool[ijk_ispec_padded] - 1;

      // copy field values
      scalar_field_accel[ijk] = potential_dot_dot_acoustic[iglob];
      scalar_field_b_accel[ijk] = b_potential_dot_dot_acoustic[iglob];
      scalar_field_b_veloc[ijk] =  b_potential_dot_acoustic[iglob];
    }
  }

  // synchronizes threads
  __syncthreads();

  if (active ){
    field accel_loc[3];
    field b_accel_loc[3];
    field b_veloc_loc[3];
    realw rhol, kappal;

    // gets material parameter
    rhol = rhostore[ijk_ispec_padded];

    // acceleration vector
    compute_gradient_kernel(ijk,ispec,ispec_irreg,
                            scalar_field_accel,accel_loc,
                            d_hprime_xx,
                            d_xix,d_xiy,d_xiz,d_etax,d_etay,d_etaz,d_gammax,d_gammay,d_gammaz,
                            rhol,xix_regular,gravity);

    // acceleration vector from backward field
    compute_gradient_kernel(ijk,ispec,ispec_irreg,
                            scalar_field_b_accel,b_accel_loc,
                            d_hprime_xx,
                            d_xix,d_xiy,d_xiz,d_etax,d_etay,d_etaz,d_gammax,d_gammay,d_gammaz,
                            rhol,xix_regular,gravity);

    // velocity vector from backward field
    compute_gradient_kernel(ijk,ispec,ispec_irreg,
                            scalar_field_b_veloc,b_veloc_loc,
                            d_hprime_xx,
                            d_xix,d_xiy,d_xiz,d_etax,d_etay,d_etaz,d_gammax,d_gammay,d_gammaz,
                            rhol,xix_regular,gravity);

    //sum function is here to sum over wavefields, in case of NB_RUNS_ACOUSTIC_GPU>1

    // approximates hessian
    hess_kl[ijk + NGLL3*ispec] += deltat * sum(accel_loc[0]*b_accel_loc[0] +
                                               accel_loc[1]*b_accel_loc[1] +
                                               accel_loc[2]*b_accel_loc[2]);


    hess_rho_ac_kl[ijk_ispec] += deltat * rhol * sum(b_veloc_loc[0]*b_veloc_loc[0] +
                                                     b_veloc_loc[1]*b_veloc_loc[1] +
                                                     b_veloc_loc[2]*b_veloc_loc[2]);
    kappal = kappastore[ijk_ispec];
    hess_kappa_ac_kl[ijk_ispec] += deltat / kappal * sum(  b_potential_dot_acoustic[iglob]
                                                         * b_potential_dot_acoustic[iglob]);
     //

  } // active
}


