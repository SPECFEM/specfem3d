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


__global__ void compute_acoustic_seismogram_kernel(int nrec_local,
                                                   field* pressure,
                                                   int* d_ibool,
                                                   realw* hxir_store, realw* hetar_store, realw* hgammar_store,
                                                   field* seismograms,
                                                   int* ispec_selected_rec_loc,
                                                   int it){

  int irec_local = blockIdx.x + blockIdx.y*gridDim.x;
  int tx = threadIdx.x;

  // local index
  int K = (tx/NGLL2);
  int J = ((tx-K*NGLL2)/NGLLX);
  int I = (tx-K*NGLL2-J*NGLLX);

  __shared__ field sh_dxd[NGLL3_PADDED];

  if (irec_local < nrec_local) {

    int ispec = ispec_selected_rec_loc[irec_local]-1;

    sh_dxd[tx] = Make_field(0.f);

    if (tx < NGLL3) {
      realw hxir = hxir_store[INDEX2(NGLLX,I,irec_local)];
      realw hetar = hetar_store[INDEX2(NGLLX,J,irec_local)];
      realw hgammar = hgammar_store[INDEX2(NGLLX,K,irec_local)];

      realw hlagrange = hxir * hetar * hgammar;
      int iglob = d_ibool[INDEX4_PADDED(NGLLX,NGLLX,NGLLX,I,J,K,ispec)]-1;

      sh_dxd[tx] = hlagrange * pressure[iglob];
    }
    __syncthreads();

    for (unsigned int s=1; s<NGLL3_PADDED ; s *= 2) {
      if (tx % (2*s) == 0) {sh_dxd[tx] += sh_dxd[tx + s];}
      __syncthreads();
    }
    //debug
    //if (tx == 0) printf("debug: seismo %i x/y = %f/%f\n",irec_local,sh_dxd[0].x,sh_dxd[0].y);

    int idx = INDEX2(nrec_local,irec_local,it);

    // Signe moins car pression = -potential_dot_dot
    if (tx == 0) seismograms[idx] = -sh_dxd[0];
  }
}

/* ----------------------------------------------------------------------------------------------- */

__global__ void compute_acoustic_vectorial_seismogram_kernel(int nrec_local,
                                                             int*  d_ispec_is_acoustic,
                                                             field* scalar_potential,
                                                             realw* seismograms,
                                                             realw* d_rhostore,
                                                             int* d_ibool,
                                                             int * d_irregular_element_number,
                                                             realw* hxir_store, realw* hetar_store, realw* hgammar_store,
                                                             realw* d_xix, realw* d_xiy, realw* d_xiz,
                                                             realw* d_etax, realw* d_etay, realw* d_etaz,
                                                             realw* d_gammax, realw* d_gammay, realw* d_gammaz,
                                                             realw xix_regular,
                                                             realw* d_hprime_xx,
                                                             realw* nu_rec,
                                                             int* ispec_selected_rec_loc,
                                                             int it){

  int irec_local = blockIdx.x + blockIdx.y*gridDim.x;
  int tx = threadIdx.x;

  // shared memory
  __shared__ realw s_dummy_loc[NGLL3_PADDED];
  __shared__ realw s_temp1[NGLL3_PADDED];
  __shared__ realw s_temp2[NGLL3_PADDED];
  __shared__ realw s_temp3[NGLL3_PADDED];
  __shared__ realw sh_hprime_xx[NGLL2];

  // locals
  realw temp1l, temp2l, temp3l;
  realw rho_invl, hlagrange;
  realw xixl, xiyl, xizl;
  realw etaxl, etayl, etazl;
  realw gammaxl, gammayl, gammazl;
  realw dpotentialdxl, dpotentialdyl, dpotentialdzl;
  int ispec, offset, offset_irreg, iglob, ispec_irreg;

  /*
  // debug
  if (irec_local < nrec_local) {
    ispec = ispec_selected_rec_loc[irec_local] - 1;
    offset = INDEX4_PADDED(NGLLX,NGLLX,NGLLX,I,J,K,ispec);
    iglob = d_ibool[offset]-1;
    rho_invl = 1.f / d_rhostore[offset];
    xixl = d_xix[offset];
    xiyl = d_xiy[offset];
    xizl = d_xiz[offset];
    etaxl = d_etax[offset];
    etayl = d_etay[offset];
    etazl = d_etaz[offset];
    gammaxl = d_gammax[offset];
    gammayl = d_gammay[offset];
    gammazl = d_gammaz[offset];

    realw hxir = hxir_store[INDEX2(NGLLX,I,irec_local)];
    realw hetar = hetar_store[INDEX2(NGLLX,J,irec_local)];
    realw hgammar = hgammar_store[INDEX2(NGLLX,K,irec_local)];

    hlagrange = hxir * hetar * hgammar;

    // loads into shared memory
    if (tx < NGLL2) {
      sh_hprime_xx[tx] = d_hprime_xx[tx];}
    s_dummy_loc[tx] = 1.; //scalar_potential[iglob];
    if (iglob > 0) {
      printf(" iglob =%d, (i,j,k)=(%d,%d,%d), ispec =%d  --- %f \n", iglob, I, J, K, ispec, scalar_potential[iglob]);}
    else{
      printf(" -illegal %d  %d %d %d %d\n", tx, ispec, I, J, K);
    }
  }
  */

  s_temp1[tx] = 0.0f;
  s_temp2[tx] = 0.0f;
  s_temp3[tx] = 0.0f;

  // local index
  int K = (tx/NGLL2);
  int J = ((tx-K*NGLL2)/NGLLX);
  int I = (tx-K*NGLL2-J*NGLLX);

  if (irec_local >= nrec_local) return;

  if (tx < NGLL3) {
    ispec = ispec_selected_rec_loc[irec_local] - 1;
    ispec_irreg = d_irregular_element_number[ispec] - 1;

    // nothing to do if we are in elastic element
    if (d_ispec_is_acoustic[ispec] == 0) {return;}

    offset = INDEX4_PADDED(NGLLX,NGLLX,NGLLX,I,J,K,ispec);
    offset_irreg = INDEX4_PADDED(NGLLX,NGLLX,NGLLX,I,J,K,ispec_irreg);

    iglob = d_ibool[offset]-1;
    rho_invl = 1.f / d_rhostore[offset];

    realw hxir = hxir_store[INDEX2(NGLLX,I,irec_local)];
    realw hetar = hetar_store[INDEX2(NGLLX,J,irec_local)];
    realw hgammar = hgammar_store[INDEX2(NGLLX,K,irec_local)];

    hlagrange = hxir * hetar * hgammar;
  }

  //debug
  //if (tx == 0) printf("thread %d %d %d - %f %f %f\n",ispec,iglob,irec_local,hlagrange,rho_invl, xixl);

  // loads into shared memory
  if (tx < NGLL2) sh_hprime_xx[tx] = d_hprime_xx[tx];
  if (tx < NGLL3) s_dummy_loc[tx] = realw_(scalar_potential[iglob]); // quick fix to convert field to realw


  // synchronize all the threads (one thread for each of the NGLL grid points of the
  // current spectral element) because we need the whole element to be ready
  __syncthreads();

  if (tx < NGLL3) {
    // computes first matrix product
    temp1l = 0.f;
    temp2l = 0.f;
    temp3l = 0.f;

    for (int l=0;l<NGLLX;l++) {
      //assumes that hprime_xx = hprime_yy = hprime_zz
      // 1. cut-plane along xi-direction
      temp1l += s_dummy_loc[K*NGLL2+J*NGLLX+l] * sh_hprime_xx[l*NGLLX+I];
      // 2. cut-plane along eta-direction
      temp2l += s_dummy_loc[K*NGLL2+l*NGLLX+I] * sh_hprime_xx[l*NGLLX+J];
      // 3. cut-plane along gamma-direction
      temp3l += s_dummy_loc[l*NGLL2+J*NGLLX+I] * sh_hprime_xx[l*NGLLX+K];
    }

   if (ispec_irreg >= 0){ //irregular element
      xixl = d_xix[offset_irreg];
      xiyl = d_xiy[offset_irreg];
      xizl = d_xiz[offset_irreg];
      etaxl = d_etax[offset_irreg];
      etayl = d_etay[offset_irreg];
      etazl = d_etaz[offset_irreg];
      gammaxl = d_gammax[offset_irreg];
      gammayl = d_gammay[offset_irreg];
      gammazl = d_gammaz[offset_irreg];
      // compute derivatives of ux, uy and uz with respect to x, y and z
      // derivatives of potential
      dpotentialdxl = xixl*temp1l + etaxl*temp2l + gammaxl*temp3l;
      dpotentialdyl = xiyl*temp1l + etayl*temp2l + gammayl*temp3l;
      dpotentialdzl = xizl*temp1l + etazl*temp2l + gammazl*temp3l;
    }
    else{
      // compute derivatives of ux, uy and uz with respect to x, y and z
      // derivatives of potential
      dpotentialdxl = xix_regular*temp1l;
      dpotentialdyl = xix_regular*temp2l;
      dpotentialdzl = xix_regular*temp3l;
    }

    // store the field in shared memmory
    s_temp1[tx] = hlagrange * dpotentialdxl * rho_invl;
    s_temp2[tx] = hlagrange * dpotentialdyl * rho_invl;
    s_temp3[tx] = hlagrange * dpotentialdzl * rho_invl;
  }

  __syncthreads();

  // reduction
  for (unsigned int s=1; s<NGLL3_PADDED ; s *= 2) {
    if (tx % (2*s) == 0){ s_temp1[tx] += s_temp1[tx + s];
                          s_temp2[tx] += s_temp2[tx + s];
                          s_temp3[tx] += s_temp3[tx + s];}
    __syncthreads();
  }

  int idx = INDEX3(NDIM,nrec_local,0,irec_local,it);

  if (tx == 0) {
    seismograms[0+idx] = nu_rec[0+NDIM*(0+NDIM*irec_local)]*s_temp1[0]
                       + nu_rec[0+NDIM*(1+NDIM*irec_local)]*s_temp2[0]
                       + nu_rec[0+NDIM*(2+NDIM*irec_local)]*s_temp3[0];
  }
  if (tx == 1) {
    seismograms[1+idx] = nu_rec[1+NDIM*(0+NDIM*irec_local)]*s_temp1[0]
                       + nu_rec[1+NDIM*(1+NDIM*irec_local)]*s_temp2[0]
                       + nu_rec[1+NDIM*(2+NDIM*irec_local)]*s_temp3[0];
  }
  if (tx == 2) {
    seismograms[2+idx] = nu_rec[2+NDIM*(0+NDIM*irec_local)]*s_temp1[0]
                       + nu_rec[2+NDIM*(1+NDIM*irec_local)]*s_temp2[0]
                       + nu_rec[2+NDIM*(2+NDIM*irec_local)]*s_temp3[0];
  }
}


