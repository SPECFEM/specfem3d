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
                                                   realw* displ,
                                                   field* potential,
                                                   int* d_ibool,
                                                   realw* hxir_store, realw* hetar_store, realw* hgammar_store,
                                                   field* seismograms,
                                                   int* ispec_selected_rec_loc,
                                                   int* ispec_is_elastic,
                                                   int* ispec_is_acoustic,
                                                   realw* kappastore, realw* mustore,
                                                   realw* d_hprime_xx,
                                                   realw* d_xix,realw* d_xiy,realw* d_xiz,
                                                   realw* d_etax,realw* d_etay,realw* d_etaz,
                                                   realw* d_gammax,realw* d_gammay,realw* d_gammaz,
                                                   int* d_irregular_element_number,
                                                   realw xix_regular,
                                                   const int ANISOTROPY,
                                                   realw* d_c11store,realw* d_c12store,realw* d_c13store,
                                                   realw* d_c14store,realw* d_c15store,realw* d_c16store,
                                                   realw* d_c22store,realw* d_c23store,realw* d_c24store,
                                                   realw* d_c25store,realw* d_c26store,realw* d_c33store,
                                                   realw* d_c34store,realw* d_c35store,realw* d_c36store,
                                                   int it){

  int irec_local = blockIdx.x + blockIdx.y*gridDim.x;
  int tx = threadIdx.x;

  // local index
  int K = (tx/NGLL2);
  int J = ((tx-K*NGLL2)/NGLLX);
  int I = (tx-K*NGLL2-J*NGLLX);

  __shared__ realw sh_tempx[NGLL3_PADDED];
  __shared__ realw sh_tempy[NGLL3_PADDED];
  __shared__ realw sh_tempz[NGLL3_PADDED];

  __shared__ field sh_dxd[NGLL3_PADDED];

  if (irec_local < nrec_local) {
    // initializes
    sh_dxd[tx] = Make_field(0.f);

    // receiver element
    int ispec = ispec_selected_rec_loc[irec_local] - 1;

    // acoustic domains
    if (ispec_is_acoustic[ispec]) {
      if (tx < NGLL3) {
        realw hxir = hxir_store[INDEX2(NGLLX,I,irec_local)];
        realw hetar = hetar_store[INDEX2(NGLLX,J,irec_local)];
        realw hgammar = hgammar_store[INDEX2(NGLLX,K,irec_local)];

        realw hlagrange = hxir * hetar * hgammar;

        int iglob = d_ibool[INDEX4_PADDED(NGLLX,NGLLX,NGLLX,I,J,K,ispec)] - 1;

        // Signe moins car pression = -potential_dot_dot
        sh_dxd[tx] = - hlagrange * potential[iglob];
      }
    }

    // elastic domains
    if (ispec_is_elastic[ispec]) {
      // computes pressure

      // for an elastic element:
      //
      // isostatic stress or pressure: p = - 1/3 trace(sigma)
      // (corresponds to hydrostatic pressure in fluids)
      //
      // from L. S. Bennethum, Compressibility Moduli for Porous Materials Incorporating Volume Fraction,
      // J. Engrg. Mech., vol. 132(11), p. 1205-1214 (2006), below equation (5):
      // for a 3D isotropic solid, pressure is defined in terms of the trace of the stress tensor as
      // p = -1/3 (t11 + t22 + t33) where t is the Cauchy stress tensor.
      //
      // to compute pressure in 3D in an elastic solid, one uses: pressure = - trace(sigma) / 3
      //
      // sigma_ij = lambda delta_ij trace(epsilon) + 2 mu epsilon_ij
      //          = lambda (epsilon_xx + epsilon_yy + epsilon_zz) + 2 mu epsilon_ij
      //
      // sigma_xx = lambda (epsilon_xx + epsilon_yy + epsilon_zz) + 2 mu epsilon_xx
      // sigma_yy = lambda (epsilon_xx + epsilon_yy + epsilon_zz) + 2 mu epsilon_yy
      // sigma_zz = lambda (epsilon_xx + epsilon_yy + epsilon_zz) + 2 mu epsilon_zz
      //
      // pressure = - trace(sigma) / 3 = - (lambda + 2/3 mu) trace(epsilon) = - kappa * trace(epsilon)
      //
      // this routines limits the pressure computations to: non-anisotropic, non-attenuation case
      // todo for the future...

      // loads displ into shared memory
      if (tx < NGLL3) {
        int iglob = d_ibool[INDEX4_PADDED(NGLLX,NGLLX,NGLLX,I,J,K,ispec)] - 1;
        sh_tempx[tx] = displ[NDIM*iglob];
        sh_tempy[tx] = displ[NDIM*iglob+1];
        sh_tempz[tx] = displ[NDIM*iglob+2];
      }

      // synchronizes threads
      __syncthreads();

      if (tx < NGLL3) {
        // derivative along x (dux_xi,duy_xi,duz_xi) (see sum_hprime_xi)
        realw temp1x = 0.0f;
        realw temp1y = 0.0f;
        realw temp1z = 0.0f;
        for(int l=0; l<NGLLX;l++){
          realw hp1 = d_hprime_xx[l*NGLLX+I];
          temp1x += sh_tempx[K*NGLL2+J*NGLLX+l] * hp1;
          temp1y += sh_tempy[K*NGLL2+J*NGLLX+l] * hp1;
          temp1z += sh_tempz[K*NGLL2+J*NGLLX+l] * hp1;
        }
        // derivative along y (duy_eta,duy_eta,duz_eta) (see sum_hprime_eta)
        realw temp2x = 0.0f;
        realw temp2y = 0.0f;
        realw temp2z = 0.0f;
        for(int l=0; l<NGLLX;l++){
          // assumes hprime_yy == hprime_zz
          realw hp2 = d_hprime_xx[l*NGLLX+J];
          temp2x += sh_tempx[K*NGLL2+l*NGLLX+I] * hp2;
          temp2y += sh_tempy[K*NGLL2+l*NGLLX+I] * hp2;
          temp2z += sh_tempz[K*NGLL2+l*NGLLX+I] * hp2;
        }
        // derivative along z (dux_gamma,duy_gamma,duz_gamma) (see sum_hprime_gamma)
        realw temp3x = 0.0f;
        realw temp3y = 0.0f;
        realw temp3z = 0.0f;
        for(int l=0; l<NGLLX;l++){
          // assumes hprime_xx == hprime_zz
          realw hp3 = d_hprime_xx[l*NGLLX+K];
          temp3x += sh_tempx[l*NGLL2+J*NGLLX+I] * hp3;
          temp3y += sh_tempy[l*NGLL2+J*NGLLX+I] * hp3;
          temp3z += sh_tempz[l*NGLL2+J*NGLLX+I] * hp3;
        }

        // derivatives of displ field on GLL points
        realw duxdxl,duydyl,duzdzl;
        realw duxdyl,duxdzl,duydxl,duydzl,duzdxl,duzdyl;

        int ispec_irreg = d_irregular_element_number[ispec] - 1;

        if (ispec_irreg >= 0){
          // irregular element
          int offset_irreg = ispec_irreg * NGLL3_PADDED + tx;
          realw xixl = d_xix[offset_irreg];
          realw xiyl = d_xiy[offset_irreg];
          realw xizl = d_xiz[offset_irreg];
          realw etaxl = d_etax[offset_irreg];
          realw etayl = d_etay[offset_irreg];
          realw etazl = d_etaz[offset_irreg];
          realw gammaxl = d_gammax[offset_irreg];
          realw gammayl = d_gammay[offset_irreg];
          realw gammazl = d_gammaz[offset_irreg];
          // derivatives
          duxdxl = xixl * temp1x + etaxl * temp2x + gammaxl * temp3x;
          duydyl = xiyl * temp1y + etayl * temp2y + gammayl * temp3y;
          duzdzl = xizl * temp1z + etazl * temp2z + gammazl * temp3z;
          if (ANISOTROPY){
            // additional non-symmetric derivatives
            duxdyl = xiyl * temp1x + etayl * temp2x + gammayl * temp3x;
            duxdzl = xizl * temp1x + etazl * temp2x + gammazl * temp3x;
            duydxl = xixl * temp1y + etaxl * temp2y + gammaxl * temp3y;
            duydzl = xizl * temp1y + etazl * temp2y + gammazl * temp3y;
            duzdxl = xixl * temp1z + etaxl * temp2z + gammaxl * temp3z;
            duzdyl = xiyl * temp1z + etayl * temp2z + gammayl * temp3z;
          }
        }else{
          // regular elements
          duxdxl = xix_regular * temp1x;
          duydyl = xix_regular * temp2y;
          duzdzl = xix_regular * temp3z;
          if (ANISOTROPY){
            // additional non-symmetric derivatives
            duxdyl = xix_regular * temp2x;
            duxdzl = xix_regular * temp3x;
            duydxl = xix_regular * temp1y;
            duydzl = xix_regular * temp3y;
            duzdxl = xix_regular * temp1z;
            duzdyl = xix_regular * temp2z;
          }
        }

        // stress trace elements
        realw sigma_xx,sigma_yy,sigma_zz;

        // full anisotropic case, stress calculations
        if (ANISOTROPY){
          int offset = ispec * NGLL3_PADDED + tx;
          realw c11 = d_c11store[offset];
          realw c12 = d_c12store[offset];
          realw c13 = d_c13store[offset];
          realw c14 = d_c14store[offset];
          realw c15 = d_c15store[offset];
          realw c16 = d_c16store[offset];
          realw c22 = d_c22store[offset];
          realw c23 = d_c23store[offset];
          realw c24 = d_c24store[offset];
          realw c25 = d_c25store[offset];
          realw c26 = d_c26store[offset];
          realw c33 = d_c33store[offset];
          realw c34 = d_c34store[offset];
          realw c35 = d_c35store[offset];
          realw c36 = d_c36store[offset];

          sigma_xx = c11*duxdxl + c16*(duxdyl + duydxl) + c12*duydyl +
                     c15*(duzdxl + duxdzl) + c14*(duzdyl + duydzl) + c13*duzdzl;
          sigma_yy = c12*duxdxl + c26*(duxdyl + duydxl) + c22*duydyl +
                     c25*(duzdxl + duxdzl) + c24*(duzdyl + duydzl) + c23*duzdzl;
          sigma_zz = c13*duxdxl + c36*(duxdyl + duydxl) + c23*duydyl +
                     c35*(duzdxl + duxdzl) + c34*(duzdyl + duydzl) + c33*duzdzl;
        }else{
          // isotropic
          int offset = ispec * NGLL3_PADDED + tx;
          realw kappal = kappastore[offset];
          realw mul = mustore[offset];

          realw lambdalplus2mul = kappal + 1.33333333333333333333f * mul;  // 4./3. = 1.3333333
          realw lambdal = lambdalplus2mul - 2.0f * mul;

          // compute diagonal components of the stress tensor
          // isotropic, non-attenuation
          sigma_xx = lambdalplus2mul * duxdxl + lambdal * (duydyl + duzdzl);
          sigma_yy = lambdalplus2mul * duydyl + lambdal * (duxdxl + duzdzl);
          sigma_zz = lambdalplus2mul * duzdzl + lambdal * (duxdxl + duydyl);
        }

        // pressure p = - 1/3 trace(sigma)
        field pressure = Make_field( (sigma_xx + sigma_yy + sigma_zz) / 3.0f);

        realw hxir = hxir_store[INDEX2(NGLLX,I,irec_local)];
        realw hetar = hetar_store[INDEX2(NGLLX,J,irec_local)];
        realw hgammar = hgammar_store[INDEX2(NGLLX,K,irec_local)];

        realw hlagrange = hxir * hetar * hgammar;

        // stores pressure
        sh_dxd[tx] = - hlagrange * pressure;
      }
    }

    // synchronizes threads
    __syncthreads();

    // reduction
    for (unsigned int s=1; s<NGLL3_PADDED ; s *= 2) {
      if (tx % (2*s) == 0) {sh_dxd[tx] += sh_dxd[tx + s];}
      __syncthreads();
    }
    //debug
    //if (tx == 0) printf("debug: seismo %i x/y = %f/%f\n",irec_local,sh_dxd[0].x,sh_dxd[0].y);

    int idx = INDEX2(nrec_local,irec_local,it);

    // seismo
    if (tx == 0) seismograms[idx] = sh_dxd[0];
  }
}

/* ----------------------------------------------------------------------------------------------- */

__global__ void compute_acoustic_vectorial_seismogram_kernel(int nrec_local,
                                                             int*  d_ispec_is_acoustic,
                                                             field* scalar_potential,
                                                             realw* seismograms,
                                                             realw* d_rhostore,
                                                             int* d_ibool,
                                                             int* d_irregular_element_number,
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
  int ispec, iglob, ispec_irreg;

  /*
  // debug
  if (irec_local < nrec_local) {
    ispec = ispec_selected_rec_loc[irec_local] - 1;
    int offset = INDEX4_PADDED(NGLLX,NGLLX,NGLLX,I,J,K,ispec);
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

    int offset = INDEX4_PADDED(NGLLX,NGLLX,NGLLX,I,J,K,ispec);

    iglob = d_ibool[offset]-1;
    rho_invl = 1.f / d_rhostore[offset];

    realw hxir = hxir_store[INDEX2(NGLLX,I,irec_local)];
    realw hetar = hetar_store[INDEX2(NGLLX,J,irec_local)];
    realw hgammar = hgammar_store[INDEX2(NGLLX,K,irec_local)];

    hlagrange = hxir * hetar * hgammar;
  }

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

   if (ispec_irreg >= 0){
      //irregular element
      int offset_irreg = INDEX4_PADDED(NGLLX,NGLLX,NGLLX,I,J,K,ispec_irreg);
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


