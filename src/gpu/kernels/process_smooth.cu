/*
!=====================================================================
!
!               S p e c f e m 3 D  V e r s i o n  3 . 0
!               ---------------------------------------
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


__global__ void process_smooth(realw_const_p xstore_me,realw_const_p ystore_me,realw_const_p zstore_me,
                               realw_const_p xstore_other,realw_const_p ystore_other,realw_const_p zstore_other,
                               realw_const_p data_other,
                               const realw sigma_h2_inv, const realw sigma_v2_inv,
                               const int iker, const int nspec_me, const int nspec_other,
                               const realw v_criterion, const realw h_criterion,
                               realw_const_p jacobian,
                               const int * irregular_element_number,realw jacobian_regular,
                               realw_p sum_data_smooth,
                               realw_p normalisation,
                               realw_const_p wgll_cube){

  int ispec = blockIdx.x + gridDim.x*blockIdx.y;
  int igll = threadIdx.x;

  int ispec_irreg;
  int gll_other;
  realw x_me,y_me, z_me, x_other,y_other, z_other, coef, normalisation_slice;
  realw dat;
  __shared__ int sh_test[NGLL3];
  __shared__ realw sh_x_other[NGLL3];
  __shared__ realw sh_y_other[NGLL3];
  __shared__ realw sh_z_other[NGLL3];
  __shared__ realw sh_jacobian[NGLL3];
  __shared__ realw sh_wgll_cube[NGLL3];
  __shared__ realw sh_data[NGLL3];

  int n_loop = nspec_other/NGLL3 + 1;

  x_me = xstore_me[NGLL3*ispec + igll ];
  y_me = ystore_me[NGLL3*ispec + igll ];
  z_me = zstore_me[NGLL3*ispec + igll ];

  sh_wgll_cube[igll] = wgll_cube[igll];

  __syncthreads();

  dat = 0.f;
  normalisation_slice = 0.f;

  //We test 125 spectral elements at a time
  for (int i=0;i<n_loop;i++){
    __syncthreads();

    if (NGLL3*i + threadIdx.x < nspec_other){
      x_other = (xstore_other[i*NGLL3*NGLL3 + threadIdx.x*NGLL3 ] + xstore_other[i*NGLL3*NGLL3 + threadIdx.x*NGLL3 + NGLL3 - 1 ])/2;
      y_other = (ystore_other[i*NGLL3*NGLL3 + threadIdx.x*NGLL3 ] + ystore_other[i*NGLL3*NGLL3 + threadIdx.x*NGLL3 + NGLL3 - 1 ])/2;
      z_other = (zstore_other[i*NGLL3*NGLL3 + threadIdx.x*NGLL3 ] + zstore_other[i*NGLL3*NGLL3 + threadIdx.x*NGLL3 + NGLL3 - 1 ])/2;
    }

    sh_test[threadIdx.x] = ( NGLL3*i + threadIdx.x >= nspec_other
                            || ((x_me-x_other)*(x_me-x_other) + (y_me-y_other)*(y_me-y_other)) > h_criterion
                            || (z_me-z_other)*(z_me-z_other) > v_criterion ) ? 1 : 0 ;
    __syncthreads();

    //loop over each spectral element tested
    for (int k=0;k<NGLL3;k++){
      __syncthreads();
      if (sh_test[k]) continue ;

      //Load data from other slice to shared memory
      sh_x_other[igll] = xstore_other[i*NGLL3*NGLL3 + k*NGLL3 + igll ];
      sh_y_other[igll] = ystore_other[i*NGLL3*NGLL3 + k*NGLL3 + igll ];
      sh_z_other[igll] = zstore_other[i*NGLL3*NGLL3 + k*NGLL3 + igll ];

      sh_data[igll] = data_other[i*NGLL3*NGLL3 + k*NGLL3 + igll ];

      ispec_irreg = irregular_element_number[i*NGLL3 + k]-1;
      if (ispec_irreg >= 0){
        sh_jacobian[igll] = jacobian[ispec_irreg*NGLL3 + igll ];
      }else{
        sh_jacobian[igll] = jacobian_regular;
      }

      __syncthreads();

      for (int j=0;j<NGLL3;j++){
        gll_other = (igll + j) % NGLL3;

        x_other = sh_x_other[gll_other];
        y_other = sh_y_other[gll_other];
        z_other = sh_z_other[gll_other];

        coef = expf(- sigma_h2_inv*((x_me-x_other)*(x_me-x_other) + (y_me-y_other)*(y_me-y_other))
                    - sigma_v2_inv*(z_me-z_other)*(z_me-z_other))
               * sh_jacobian[gll_other] * sh_wgll_cube[gll_other];

        normalisation_slice = normalisation_slice + coef;
        dat += sh_data[gll_other]*coef;
      } //loop on each gll_other
    } //loop on each spec_other tested
  } //loop on each serie of 125 spec_other

  sum_data_smooth[NGLL3*nspec_me*iker+NGLL3*ispec + igll] += dat;
  normalisation[NGLL3*ispec + igll] += normalisation_slice;
}

/* ----------------------------------------------------------------------------------------------- */

__global__ void normalize_data(realw_p data_smooth, realw_const_p normalisation,int nker, int nspec_me){
  int ispec = blockIdx.x + gridDim.x*blockIdx.y;

  realw norm = normalisation[NGLL3*ispec + threadIdx.x];
  for (int j=0;j<nker;j++) data_smooth[NGLL3*nspec_me*j + NGLL3*ispec + threadIdx.x] /= norm/nker;
}


