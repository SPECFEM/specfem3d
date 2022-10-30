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


__global__ void compute_elastic_seismogram_kernel(int nrec_local,
                                                  realw* displ,
                                                  field* potential,
                                                  int* d_ibool,
                                                  realw* hxir_store, realw* hetar_store, realw* hgammar_store,
                                                  realw* seismograms,
                                                  realw* nu_rec,
                                                  int* ispec_selected_rec_loc,
                                                  int* ispec_is_elastic,
                                                  int* ispec_is_acoustic,
                                                  realw* rhostore,
                                                  realw* d_hprime_xx,
                                                  realw* d_xix,realw* d_xiy,realw* d_xiz,
                                                  realw* d_etax,realw* d_etay,realw* d_etaz,
                                                  realw* d_gammax,realw* d_gammay,realw* d_gammaz,
                                                  int* d_irregular_element_number,
                                                  realw xix_regular,
                                                  int it){

  int irec_local = blockIdx.x + blockIdx.y*gridDim.x;
  int tx = threadIdx.x;

  // local index
  int K = (tx/NGLL2);
  int J = ((tx-K*NGLL2)/NGLLX);
  int I = (tx-K*NGLL2-J*NGLLX);

  __shared__ realw sh_dxd[NGLL3_PADDED];
  __shared__ realw sh_dyd[NGLL3_PADDED];
  __shared__ realw sh_dzd[NGLL3_PADDED];

  __shared__ field scalar_field[NGLL3];

  if (irec_local < nrec_local) {
    // initializes
    sh_dxd[tx] = 0.0f;
    sh_dyd[tx] = 0.0f;
    sh_dzd[tx] = 0.0f;

    // receiver element
    int ispec = ispec_selected_rec_loc[irec_local] - 1;

    // elastic domains
    if (ispec_is_elastic[ispec]) {
      if (tx < NGLL3) {
        realw hxir = hxir_store[INDEX2(NGLLX,I,irec_local)];
        realw hetar = hetar_store[INDEX2(NGLLX,J,irec_local)];
        realw hgammar = hgammar_store[INDEX2(NGLLX,K,irec_local)];

        realw hlagrange = hxir * hetar * hgammar;
        int iglob = d_ibool[INDEX4_PADDED(NGLLX,NGLLX,NGLLX,I,J,K,ispec)] - 1;

        sh_dxd[tx] = hlagrange * displ[NDIM*iglob];
        sh_dyd[tx] = hlagrange * displ[NDIM*iglob+1];
        sh_dzd[tx] = hlagrange * displ[NDIM*iglob+2];

        //debug
        //if (tx == 0) printf("thread %d %d %d - %f %f %f\n",ispec,iglob,irec_local,hlagrange,displ[0 + 2*iglob],displ[1 + 2*iglob]);
      }
    }

    // acoustic domains
    if (ispec_is_acoustic[ispec]) {
      // loads scalar into shared memory
      if (tx < NGLL3) {
        int iglob = d_ibool[INDEX4_PADDED(NGLLX,NGLLX,NGLLX,I,J,K,ispec)] - 1;
        scalar_field[tx] = potential[iglob];
      }

      // synchronizes threads
      __syncthreads();

      if (tx < NGLL3) {
        // compute gradient of potential to calculate vector if acoustic element
        // we then need to divide by density because the potential is a potential of (density * displacement)
        int offset = INDEX4_PADDED(NGLLX,NGLLX,NGLLX,I,J,K,ispec);

        realw rhol = rhostore[offset];

        // gradient vector
        field vec_elem[3];
        const int gravity = 0; // not used yet
        int ispec_irreg = d_irregular_element_number[ispec] - 1;

        compute_gradient_kernel(tx,ispec,ispec_irreg,scalar_field,vec_elem,
                                d_hprime_xx,
                                d_xix,d_xiy,d_xiz,d_etax,d_etay,d_etaz,d_gammax,d_gammay,d_gammaz,
                                rhol,xix_regular,gravity);

        realw hxir = hxir_store[INDEX2(NGLLX,I,irec_local)];
        realw hetar = hetar_store[INDEX2(NGLLX,J,irec_local)];
        realw hgammar = hgammar_store[INDEX2(NGLLX,K,irec_local)];

        realw hlagrange = hxir * hetar * hgammar;

        sh_dxd[tx] = hlagrange * realw_(vec_elem[0]);  // realw_(..) quick fix to convert field to realw
        sh_dyd[tx] = hlagrange * realw_(vec_elem[1]);
        sh_dzd[tx] = hlagrange * realw_(vec_elem[2]);
      }
    }

    // synchronizes threads
    __syncthreads();

    // reduction
    for (unsigned int s=1; s<NGLL3_PADDED ; s *= 2) {
      if (tx % (2*s) == 0){ sh_dxd[tx] += sh_dxd[tx + s];
                            sh_dyd[tx] += sh_dyd[tx + s];
                            sh_dzd[tx] += sh_dzd[tx + s];}
      __syncthreads();
    }

    int idx = INDEX3(NDIM,nrec_local,0,irec_local,it);

    // component rotation
    // seismograms_d(:,irec_local,seismo_current) = rotation_seismo(:,1)*dxd + rotation_seismo(:,2)*dyd + rotation_seismo(:,3)*dzd
    if (tx == 0) {
      seismograms[0+idx] = nu_rec[0+NDIM*(0+NDIM*irec_local)]*sh_dxd[0]
                         + nu_rec[0+NDIM*(1+NDIM*irec_local)]*sh_dyd[0]
                         + nu_rec[0+NDIM*(2+NDIM*irec_local)]*sh_dzd[0];
    }
    if (tx == 1) {
      seismograms[1+idx] = nu_rec[1+NDIM*(0+NDIM*irec_local)]*sh_dxd[0]
                         + nu_rec[1+NDIM*(1+NDIM*irec_local)]*sh_dyd[0]
                         + nu_rec[1+NDIM*(2+NDIM*irec_local)]*sh_dzd[0];
    }
    if (tx == 2) {
      seismograms[2+idx] = nu_rec[2+NDIM*(0+NDIM*irec_local)]*sh_dxd[0]
                         + nu_rec[2+NDIM*(1+NDIM*irec_local)]*sh_dyd[0]
                         + nu_rec[2+NDIM*(2+NDIM*irec_local)]*sh_dzd[0];
    }
  }
}

