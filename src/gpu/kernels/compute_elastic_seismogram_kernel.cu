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


__global__ void compute_elastic_seismogram_kernel(int nrec_local,
                                                  realw* field,
                                                  int* d_ibool,
                                                  realw* hxir_store, realw* hetar_store, realw* hgammar_store,
                                                  realw* seismograms,
                                                  realw* nu_rec,
                                                  int* ispec_selected_rec_loc,
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

  if (irec_local < nrec_local) {

    int ispec = ispec_selected_rec_loc[irec_local] - 1;

    sh_dxd[tx] = 0;
    sh_dyd[tx] = 0;
    sh_dzd[tx] = 0;

    if (tx < NGLL3) {
      realw hxir = hxir_store[INDEX2(NGLLX,I,irec_local)];
      realw hetar = hetar_store[INDEX2(NGLLX,J,irec_local)];
      realw hgammar = hgammar_store[INDEX2(NGLLX,K,irec_local)];

      realw hlagrange = hxir * hetar * hgammar;
      int iglob = d_ibool[INDEX4_PADDED(NGLLX,NGLLX,NGLLX,I,J,K,ispec)]-1;

      sh_dxd[tx] = hlagrange * field[0 + NDIM*iglob];
      sh_dyd[tx] = hlagrange * field[1 + NDIM*iglob];
      sh_dzd[tx] = hlagrange * field[2 + NDIM*iglob];

      //debug
      //if (tx == 0) printf("thread %d %d %d - %f %f %f\n",ispec,iglob,irec_local,hlagrange,field[0 + 2*iglob],field[1 + 2*iglob]);
    }
    __syncthreads();

    // reduction
    for (unsigned int s=1; s<NGLL3_PADDED ; s *= 2) {
      if (tx % (2*s) == 0){ sh_dxd[tx] += sh_dxd[tx + s];
                            sh_dyd[tx] += sh_dyd[tx + s];
                            sh_dzd[tx] += sh_dzd[tx + s];}
      __syncthreads();
    }

    int idx = INDEX3(NDIM,nrec_local,0,irec_local,it);

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

