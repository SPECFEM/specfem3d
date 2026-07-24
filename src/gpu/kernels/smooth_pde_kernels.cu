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

__global__ void Kernel_2_smooth_pde(const int nb_blocks_to_compute,
                                    const int* d_ibool,
                                    const int* d_irregular_element_number,
                                    const int* d_phase_ispec_inner,
                                    const int num_phase_ispec,
                                    const int d_iphase,
                                    field_const_p d_dat_smooth_glob,
                                    field_p d_ddat_smooth_glob,
                                    realw* d_xix,realw* d_xiy,realw* d_xiz,
                                    realw* d_etax,realw* d_etay,realw* d_etaz,
                                    realw* d_gammax,realw* d_gammay,realw* d_gammaz,
                                    const realw xix_regular, const realw jacobian_regular,
                                    realw_const_p d_hprime_xx,
                                    realw_const_p d_hprimewgll_xx,
                                    realw_const_p d_wgllwgll_xy,
                                    realw_const_p d_wgllwgll_xz,
                                    realw_const_p d_wgllwgll_yz,
                                    const realw cv, const realw ch){
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
  // limits thread ids to range [0,125-1]
  if (tx >= NGLL3) tx = NGLL3-1;

  // local index
  int K = (tx/NGLL2);
  int J = ((tx-K*NGLL2)/NGLLX);
  int I = (tx-K*NGLL2-J*NGLLX);

  int iglob,offset;
  int working_element,ispec_irreg;

  field temp1l,temp2l,temp3l;
  realw xixl,xiyl,xizl,etaxl,etayl,etazl,gammaxl,gammayl,gammazl;
  realw jacobianl;

  field ddxl,ddyl,ddzl;
  realw fac1,fac2,fac3;

  field sum_terms;

  __shared__ field s_dummy_loc[NGLL3];

  __shared__ field s_temp1[NGLL3];
  __shared__ field s_temp2[NGLL3];
  __shared__ field s_temp3[NGLL3];

  __shared__ realw sh_hprime_xx[NGLL2];
  __shared__ realw sh_hprimewgll_xx[NGLL2];

  // spectral-element id
  // iphase-1 and working_element-1 for Fortran->C array conventions
  working_element = d_phase_ispec_inner[bx + num_phase_ispec*(d_iphase-1)]-1;

  ispec_irreg = d_irregular_element_number[working_element] - 1;

  // local padded index
  offset = working_element*NGLL3_PADDED + tx;

  // global index
  iglob = d_ibool[offset] - 1;

  // loads potential values into shared memory
  if (threadIdx.x < NGLL3) {
    // changing iglob indexing to match fortran row changes fast style
    s_dummy_loc[tx] = d_dat_smooth_glob[iglob];
  }

  // calculates laplacian
  if (ispec_irreg >= 0){ //irregular_element
    int offset = ispec_irreg*NGLL3_PADDED + tx;
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
  }

  if (tx < NGLL2) {
    sh_hprime_xx[tx] = d_hprime_xx[tx];
    // loads hprimewgll into shared memory
    sh_hprimewgll_xx[tx] = d_hprimewgll_xx[tx];
  }

  // synchronize all the threads (one thread for each of the NGLL grid points of the
  // current spectral element) because we need the whole element to be ready in order
  // to be able to compute the matrix products along cut planes of the 3D element below
  __syncthreads();

  // summed terms with added gll weights
  fac1 = d_wgllwgll_yz[K*NGLLX+J];
  fac2 = d_wgllwgll_xz[K*NGLLX+I];
  fac3 = d_wgllwgll_xy[J*NGLLX+I];

  // computes first matrix product
  temp1l = Make_field(0.f);
  temp2l = Make_field(0.f);
  temp3l = Make_field(0.f);

  for (int l=0;l<NGLLX;l++) {
    //assumes that hprime_xx = hprime_yy = hprime_zz
    // 1. cut-plane along xi-direction
    temp1l += s_dummy_loc[K*NGLL2+J*NGLLX+l] * sh_hprime_xx[l*NGLLX+I];
    // 2. cut-plane along eta-direction
    temp2l += s_dummy_loc[K*NGLL2+l*NGLLX+I] * sh_hprime_xx[l*NGLLX+J];
    // 3. cut-plane along gamma-direction
    temp3l += s_dummy_loc[l*NGLL2+J*NGLLX+I] * sh_hprime_xx[l*NGLLX+K];
  }

  // compute derivatives of ux, uy and uz with respect to x, y and z
  // derivatives of potential
  if (threadIdx.x < NGLL3) {
    if (ispec_irreg >= 0){ //irregular_element

      ddxl = xixl*temp1l + etaxl*temp2l + gammaxl*temp3l;
      ddyl = xiyl*temp1l + etayl*temp2l + gammayl*temp3l;
      ddzl = xizl*temp1l + etazl*temp2l + gammazl*temp3l;

      // form the dot product with the test vector
      s_temp1[tx] = ((cv-ch) * xizl * ddzl +
              ch * (xixl*ddxl+xiyl*ddyl+xizl*ddzl)) * jacobianl;
      s_temp2[tx] = ((cv-ch) * etazl * ddzl +
              ch * (etaxl*ddxl+etayl*ddyl+etazl*ddzl)) * jacobianl;
      s_temp3[tx] = ((cv-ch) * gammazl * ddzl +
              ch * (gammaxl*ddxl+gammayl*ddyl+gammazl*ddzl)) * jacobianl;
    }
    else{
      s_temp1[tx] = (ch * xix_regular*xix_regular*temp1l) * jacobian_regular;
      s_temp2[tx] = (ch * xix_regular*xix_regular*temp2l) * jacobian_regular;
      s_temp3[tx] = (cv * xix_regular*xix_regular*temp3l) * jacobian_regular;
    }
  }

  // synchronize all the threads (one thread for each of the NGLL grid points of the
  // current spectral element) because we need the whole element to be ready in order
  // to be able to compute the matrix products along cut planes of the 3D element below
  __syncthreads();

  // computes second matrix product
  temp1l = Make_field(0.f);
  temp2l = Make_field(0.f);
  temp3l = Make_field(0.f);

  for (int l=0;l<NGLLX;l++) {
    //assumes hprimewgll_xx = hprimewgll_yy = hprimewgll_zz
    // 1. cut-plane along xi-direction
    temp1l += s_temp1[K*NGLL2+J*NGLLX+l] * sh_hprimewgll_xx[I*NGLLX+l];
    // 2. cut-plane along eta-direction
    temp2l += s_temp2[K*NGLL2+l*NGLLX+I] * sh_hprimewgll_xx[J*NGLLX+l];
    // 3. cut-plane along gamma-direction
    temp3l += s_temp3[l*NGLL2+J*NGLLX+I] * sh_hprimewgll_xx[K*NGLLX+l];
  }

  sum_terms = -(fac1*temp1l + fac2*temp2l + fac3*temp3l);

   __syncthreads();
// assembles potential array
  if (threadIdx.x < NGLL3) {
    atomicAdd(&d_ddat_smooth_glob[iglob], sum_terms);
  }
}


__global__ void kernel_3_smooth_pde_cuda_device(field* d_ddat_smooth_glob,
                                                realw_const_p rvol,
                                                int size) {
  int id = threadIdx.x + (blockIdx.x + blockIdx.y*gridDim.x)*blockDim.x;

  if (id < size) {
    realw r = rvol[id];
    d_ddat_smooth_glob[id] *= r;
  }
}

__global__ void UpdateData_smooth_pde_kernel(field* d_dat_smooth_glob,
                                             field* d_ddat_smooth_glob,
                                             int size) {
  int id = threadIdx.x + (blockIdx.x + blockIdx.y*gridDim.x)*blockDim.x;
  if (id < size) {
    field ddat = d_ddat_smooth_glob[id];
    d_dat_smooth_glob[id] += ddat;
    d_ddat_smooth_glob[id] = Make_field(0.f);
  }
}

__global__ void zero_pml_smooth_pde_kernel(int nb_blocks_to_compute,
                                           field * d_dat_smooth_glob,
                                           const int * d_ibool,
                                           const int * d_CPML_to_spec) {
  // block-id == number of local element id in phase_ispec array
  int bx = blockIdx.y*gridDim.x+blockIdx.x;

  // checks if anything to do
  if (bx >= nb_blocks_to_compute) return;

  int tx = threadIdx.x;
  if (tx >= NGLL3) tx = NGLL3 - 1;

  int working_element = d_CPML_to_spec[bx] - 1;

  // local padded index
  int offset = working_element*NGLL3_PADDED + tx;

  // global index
  int iglob = d_ibool[offset] - 1;

  d_dat_smooth_glob[iglob] = Make_field(0.f);
}
