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


__global__ void compute_kernels_hess_el_cudakernel(int* ispec_is_elastic,
                                                   int* d_ibool,
                                                   realw* accel,
                                                   realw* b_accel,
                                                   realw* b_veloc,
                                                   realw* b_epsilondev_xx,realw* b_epsilondev_yy,realw* b_epsilondev_xy,
                                                   realw* b_epsilondev_xz,realw* b_epsilondev_yz,realw* b_epsilon_trace_over_3,
                                                   realw* hess_kl,
                                                   realw* hess_rho_kl,
                                                   realw* hess_kappa_kl,
                                                   realw* hess_mu_kl,
                                                   realw deltat,
                                                   int NSPEC_AB) {

  int ispec = blockIdx.x + blockIdx.y*gridDim.x;
  int ijk = threadIdx.x;
  int ijk_ispec = ijk + NGLL3*ispec;

  // handles case when there is 1 extra block (due to rectangular grid)
  if (ispec < NSPEC_AB) {

    // elastic elements only
    if (ispec_is_elastic[ispec]) {
      int iglob = d_ibool[ijk + NGLL3_PADDED*ispec] - 1;

      // approximate hessian
      hess_kl[ijk + NGLL3*ispec] += deltat * (accel[3*iglob]*b_accel[3*iglob]+
                                              accel[3*iglob+1]*b_accel[3*iglob+1]+
                                              accel[3*iglob+2]*b_accel[3*iglob+2]);

      //
      hess_rho_kl[ijk_ispec] += deltat * (b_veloc[3*iglob]  *b_veloc[3*iglob]+
            b_veloc[3*iglob+1]*b_veloc[3*iglob+1]+
            b_veloc[3*iglob+2]*b_veloc[3*iglob+2]);

      hess_mu_kl[ijk_ispec] += deltat * (b_epsilondev_xx[ijk_ispec]*b_epsilondev_xx[ijk_ispec]+
           b_epsilondev_yy[ijk_ispec]*b_epsilondev_yy[ijk_ispec]+
          (b_epsilondev_xx[ijk_ispec]+b_epsilondev_yy[ijk_ispec])*
          (b_epsilondev_xx[ijk_ispec]+b_epsilondev_yy[ijk_ispec])+
              2*(b_epsilondev_xy[ijk_ispec]*b_epsilondev_xy[ijk_ispec]+
                                         b_epsilondev_xz[ijk_ispec]*b_epsilondev_xz[ijk_ispec]+
                                         b_epsilondev_yz[ijk_ispec]*b_epsilondev_yz[ijk_ispec]));

      hess_kappa_kl[ijk_ispec] += deltat*(9*b_epsilon_trace_over_3[ijk_ispec]*b_epsilon_trace_over_3[ijk_ispec]);

      /*
      if (ijk_ispec==100){
        printf(" Hessian %e  %e \n",b_epsilondev_xx[ijk_ispec], b_epsilondev_yy[ijk_ispec]);
      }
      */
    }
  }
}


