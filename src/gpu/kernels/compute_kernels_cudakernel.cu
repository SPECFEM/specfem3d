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


__global__ void compute_kernels_cudakernel(int* ispec_is_elastic,
                                           int* d_ibool,
                                           realw* accel,
                                           realw* b_displ,
                                           realw* epsilondev_xx,realw* epsilondev_yy,realw* epsilondev_xy,
                                           realw* epsilondev_xz,realw* epsilondev_yz,
                                           realw* b_epsilondev_xx,realw* b_epsilondev_yy,realw* b_epsilondev_xy,
                                           realw* b_epsilondev_xz,realw* b_epsilondev_yz,
                                           realw* rho_kl,
                                           realw deltat,
                                           realw* mu_kl,
                                           realw* kappa_kl,
                                           realw* epsilon_trace_over_3,
                                           realw* b_epsilon_trace_over_3,
                                           int NSPEC_AB) {

  int ispec = blockIdx.x + blockIdx.y*gridDim.x;
  int ijk = threadIdx.x;
  int ijk_ispec = ijk + NGLL3*ispec;

  // handles case when there is 1 extra block (due to rectangular grid)
  if (ispec < NSPEC_AB) {

    // elastic elements only
    if (ispec_is_elastic[ispec]) {
      int iglob = d_ibool[ijk + NGLL3_PADDED*ispec] - 1 ;

      // isotropic kernels:
      // density kernel
      rho_kl[ijk_ispec] += deltat * (accel[3*iglob]*b_displ[3*iglob]+
                                     accel[3*iglob+1]*b_displ[3*iglob+1]+
                                     accel[3*iglob+2]*b_displ[3*iglob+2]);


      // shear modulus kernel
      mu_kl[ijk_ispec] += deltat * (epsilondev_xx[ijk_ispec]*b_epsilondev_xx[ijk_ispec]+
                                    epsilondev_yy[ijk_ispec]*b_epsilondev_yy[ijk_ispec]+
                                    (epsilondev_xx[ijk_ispec]+epsilondev_yy[ijk_ispec])*
                                    (b_epsilondev_xx[ijk_ispec]+b_epsilondev_yy[ijk_ispec])+
                                    2*(epsilondev_xy[ijk_ispec]*b_epsilondev_xy[ijk_ispec]+
                                       epsilondev_xz[ijk_ispec]*b_epsilondev_xz[ijk_ispec]+
                                       epsilondev_yz[ijk_ispec]*b_epsilondev_yz[ijk_ispec]));

      // bulk modulus kernel
      kappa_kl[ijk_ispec] += deltat*(9*epsilon_trace_over_3[ijk_ispec]*b_epsilon_trace_over_3[ijk_ispec]);

      /*
      if (ijk_ispec==100){
        printf(" Kernel, %e  %e \n",b_epsilondev_xx[ijk_ispec], b_epsilondev_yy[ijk_ispec]);
      }
      */
    }
  }
}


