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


__global__ void process_smooth(realw_const_p xstore_me,
                               realw_const_p ystore_me,
                               realw_const_p zstore_me,
                               realw_const_p xstore_other,
                               realw_const_p ystore_other,
                               realw_const_p zstore_other,
                               realw_const_p data_other,
                               const realw sigma_h2_inv, const realw sigma_v2_inv,
                               const int iker,
                               const int nspec_me,
                               const int nspec_other,
                               const realw v_criterion, const realw h_criterion,
                               realw_const_p integ_factor,
                               realw_p data_smooth,
                               realw_p normalisation){

  int ispec = blockIdx.x + gridDim.x*blockIdx.y;
  int igll = threadIdx.x;

  int gll_other;
  realw x_me,y_me,z_me;
  realw x_other,y_other,z_other;
  realw center_x,center_y,center_z;
  realw dist_h,dist_v;
  realw val,val_gaussian;
  realw coef,normalisation_slice;
  realw dat;

  __shared__ int sh_test[NGLL3];
  __shared__ realw sh_x_other[NGLL3];
  __shared__ realw sh_y_other[NGLL3];
  __shared__ realw sh_z_other[NGLL3];
  __shared__ realw sh_integ_factor[NGLL3];
  __shared__ realw sh_data[NGLL3];

  // for each reference GLL point, we can check a block of 125 neighbor elements
  // by convenience, the block size is set to the number of threads 125 of this kernel
  int n_loop = nspec_other/NGLL3 + 1;

  // reference GLL point position
  x_me = xstore_me[NGLL3*ispec + igll ];
  y_me = ystore_me[NGLL3*ispec + igll ];
  z_me = zstore_me[NGLL3*ispec + igll ];

  __syncthreads();

  dat = 0.f;
  normalisation_slice = 0.f;

  //We test 125 spectral elements at a time
  for (int i=0; i < n_loop; i++){
    __syncthreads();

    // each thread helps to test a different element in the other slice (using the center position)
    // number of threads == NGLL3 == 125
    // for i==0: element range [0,124]
    // for i==1: element range [125,(125+124)]
    // ..
    // for i==n_loop-1: element range [NGLL3*(nloop-1),NGLL3*(nloop-1)+124]
    //                  where NGLL3*(nloop-1)+124 is equal to nspec_other (or slightly greater)
    int ispec_other = NGLL3*i + igll;

    if (ispec_other < nspec_other){
      // center position
      center_x = (xstore_other[ispec_other * NGLL3] + xstore_other[ispec_other * NGLL3 + (NGLL3 - 1)]) * 0.5f;
      center_y = (ystore_other[ispec_other * NGLL3] + ystore_other[ispec_other * NGLL3 + (NGLL3 - 1)]) * 0.5f;
      center_z = (zstore_other[ispec_other * NGLL3] + zstore_other[ispec_other * NGLL3 + (NGLL3 - 1)]) * 0.5f;

      // note: instead of distance we use distance squared to avoid too many sqrt() operations

      // Cartesian case
      // distance horizontal = (x-x0)**2 + (y-y0)**2, and vertical = (z-z0)**2
      dist_h = (x_me - center_x)*(x_me - center_x) + (y_me - center_y)*(y_me - center_y);
      dist_v = (z_me - center_z)*(z_me - center_z);
    } else {
      // artificial high values
      dist_v = 99999999.f;
      dist_h = 99999999.f;
    }

    // tests if element is too far away
    sh_test[igll] = ( ispec_other >= nspec_other
                    || dist_h > h_criterion
                    || dist_v > v_criterion ) ? 1 : 0 ;

    __syncthreads();

    // loops over each spectral element tested
    for (int k=0; k < NGLL3; k++){
      __syncthreads();

      // skips element if test was true (too far away)
      if (sh_test[k]) continue ;

      // loads data from other slice to shared memory
      int ispec_test = i*NGLL3 + k;
      sh_x_other[igll] = xstore_other[ispec_test*NGLL3 + igll];
      sh_y_other[igll] = ystore_other[ispec_test*NGLL3 + igll];
      sh_z_other[igll] = zstore_other[ispec_test*NGLL3 + igll];

      sh_data[igll] = data_other[ispec_test*NGLL3 + igll];
      sh_integ_factor[igll] = integ_factor[ispec_test*NGLL3 + igll];

      __syncthreads();

      // loops over gll points
      for (int j=0; j < NGLL3; j++){
        gll_other = (igll + j) % NGLL3;

        x_other = sh_x_other[gll_other];
        y_other = sh_y_other[gll_other];
        z_other = sh_z_other[gll_other];

        // Cartesian case
        // distance horizontal = (x-x0)**2 + (y-y0)**2, and vertical = (z-z0)**2
        dist_h = (x_me - x_other)*(x_me - x_other) + (y_me - y_other)*(y_me - y_other);
        dist_v = (z_me - z_other)*(z_me - z_other);

        // Gaussian function
        val = - dist_h*sigma_h2_inv - dist_v*sigma_v2_inv;

        // limits to single precision
        if (val < - 86.0f){
          // smaller than numerical precision: exp(-86) < 1.e-37
          val_gaussian = 0.0f;
        } else {
          val_gaussian = expf(val);
        }

        coef = val_gaussian * sh_integ_factor[gll_other];

        normalisation_slice = normalisation_slice + coef;
        dat += sh_data[gll_other] * coef;
      } //loop on each gll_other
    } //loop on each spec_other tested
  } //loop on each serie of 125 spec_other

  data_smooth[NGLL3*nspec_me*iker + NGLL3*ispec + igll] += dat;

  // note: normalization coefficient is added nker times
  normalisation[NGLL3*ispec + igll] += normalisation_slice;
}

/* ----------------------------------------------------------------------------------------------- */

__global__ void normalize_data(realw_p data_smooth,
                               realw_const_p normalisation,
                               int nker,
                               int nspec_me){

  int ispec = blockIdx.x + gridDim.x*blockIdx.y;
  int igll = threadIdx.x;

  // note: normalization coefficient is added nker times, thus divide by nker
  realw norm = normalisation[NGLL3*ispec + igll] / nker;

  // avoids division by zero
  if (norm < 1.e-24) norm = 1.0f;

  // normalizes smoothed kernel values
  for (int iker=0; iker<nker; iker++) data_smooth[NGLL3*nspec_me*iker + NGLL3*ispec + igll] /= norm;
}


