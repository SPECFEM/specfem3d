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


__global__ void compute_kernels_acoustic_kernel(int* ispec_is_acoustic,
                                                int* d_ibool,
                                                realw* rhostore,
                                                realw* d_hprime_xx,
                                                int* d_irregular_element_number,
                                                realw* d_xix,realw* d_xiy,realw* d_xiz,
                                                realw* d_etax,realw* d_etay,realw* d_etaz,
                                                realw* d_gammax,realw* d_gammay,realw* d_gammaz,
                                                realw xix_regular,
                                                field* potential_acoustic,
                                                field* potential_dot_dot_acoustic,
                                                field* b_potential_acoustic,
                                                field* b_potential_dot_dot_acoustic,
                                                realw* rho_ac_kl,
                                                realw* kappa_ac_kl,
                                                realw deltat,
                                                int NSPEC_AB,
                                                int gravity) {

  int ispec = blockIdx.x + blockIdx.y*gridDim.x;
  int ijk = threadIdx.x;

  // local and global indices
  int iglob;

  int ijk_ispec = ijk + NGLL3*ispec;
  int ijk_ispec_padded = ijk + NGLL3_PADDED*ispec;

  int ispec_irreg = d_irregular_element_number[ispec] - 1;

  // shared memory between all threads within this block
  __shared__ field scalar_field_displ[NGLL3];
  __shared__ field scalar_field_accel[NGLL3];

  int active = 0;

  // handles case when there is 1 extra block (due to rectangular grid)
  if (ispec < NSPEC_AB ){
    // acoustic elements only
    if (ispec_is_acoustic[ispec] ){
      active = 1;

      // copy field values
      iglob = d_ibool[ijk_ispec_padded] - 1;
      scalar_field_displ[ijk] = b_potential_acoustic[iglob];
      scalar_field_accel[ijk] = potential_acoustic[iglob];
    }
  }

  // synchronizes threads
  __syncthreads();

  if (active ){
    field accel_loc[3];
    field b_displ_loc[3];
    realw rhol;

    // gets material parameter
    rhol = rhostore[ijk_ispec_padded];

    // displacement vector from backward field
    compute_gradient_kernel(ijk,ispec,ispec_irreg,scalar_field_displ,b_displ_loc,
                            d_hprime_xx,
                            d_xix,d_xiy,d_xiz,d_etax,d_etay,d_etaz,d_gammax,d_gammay,d_gammaz,
                            rhol,xix_regular,gravity);

    // acceleration vector
    compute_gradient_kernel(ijk,ispec,ispec_irreg,scalar_field_accel,accel_loc,
                            d_hprime_xx,
                            d_xix,d_xiy,d_xiz,d_etax,d_etay,d_etaz,d_gammax,d_gammay,d_gammaz,
                            rhol,xix_regular,gravity);

    // note: the gradients grad(..) above will have a 1/rho factor added in routine compute_gradient_kernel()
    //       (as in compute_gradient_in_acoustic() for CPU)
    //
    //       that is, b_displ_elm = 1/rho grad(phi)
    //       and      accel_elm   = 1/rho grad(phi^adj)
    //
    //       however, here we want to compute the contributions for the absolute kernel
    //         K_rho   = - int_0^T [ grad(phi^adj) * grad(phi) ] dt        (A-28)
    //       and thus we multiply the b_displ_elm and accel_elm by rho again to have the time step contribution
    //         b_displ_elem * accel_elm = {rho [1/rho grad(phi)]} * {rho [1/rho grad(phi^adj)]}
    //                                  = grad(phi) * grad(phi^adj) and times dt
    //
    // we could also just set rhol above to 1.0 and call compute_gradient_kernel(), but prefer here to make this more explicit.
    // in case the gradient routine changes, revisit this...
    accel_loc[0] = rhol * accel_loc[0];
    accel_loc[1] = rhol * accel_loc[1];
    accel_loc[2] = rhol * accel_loc[2];

    b_displ_loc[0] = rhol * b_displ_loc[0];
    b_displ_loc[1] = rhol * b_displ_loc[1];
    b_displ_loc[2] = rhol * b_displ_loc[2];

    // the sum function is here to enable summing over wavefields when NB_RUNS_ACOUSTIC_GPU > 1

    // density kernel
    // (multiplication with rho(x) factor will be done after the time loop)
    rho_ac_kl[ijk_ispec] += deltat * sum(accel_loc[0]*b_displ_loc[0] + accel_loc[1]*b_displ_loc[1] + accel_loc[2]*b_displ_loc[2]);

    // bulk modulus kernel
    // (multiplication with 1/kappa(x) factor will be done after the time loop)
    kappa_ac_kl[ijk_ispec] += deltat * sum(potential_acoustic[iglob] * b_potential_dot_dot_acoustic[iglob]);
  } // active
}


