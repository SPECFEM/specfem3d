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


__global__ void add_sources_el_SIM_TYPE_2_OR_3_kernel(realw* accel,
                                                      int nrec,
                                                      int it,
                                                      int NSTEP_BETWEEN_ADJSRC,
                                                      field* source_adjoint,
                                                      realw* xir_store,
                                                      realw* etar_store,
                                                      realw* gammar_store,
                                                      int* d_ibool,
                                                      int* ispec_is_elastic,
                                                      int* ispec_selected_recloc,
                                                      int nadj_rec_local) {

  int irec_local = blockIdx.x + gridDim.x*blockIdx.y;

  if (irec_local < nadj_rec_local) { // when nrec > 65535, but mod(nspec_top,2) > 0, we end up with an extra block.

    int ispec = ispec_selected_recloc[irec_local]-1;

    if (ispec_is_elastic[ispec]){
      int i = threadIdx.x;
      int j = threadIdx.y;
      int k = threadIdx.z;
      int iglob = d_ibool[INDEX4_PADDED(NGLLX,NGLLX,NGLLX,i,j,k,ispec)]-1;

      realw hxir    = xir_store[INDEX2(NGLLX,i,irec_local)];
      realw hetar   = etar_store[INDEX2(NGLLX,j,irec_local)];
      realw hgammar = gammar_store[INDEX2(NGLLX,k,irec_local)];

      realw lagrange =   hxir * hetar * hgammar ;

      // atomic operations are absolutely necessary for correctness!
      atomicAdd(&accel[0+3*iglob],source_adjoint[INDEX3(NDIM,nadj_rec_local,0,irec_local,it)]*lagrange);
      atomicAdd(&accel[1+3*iglob],source_adjoint[INDEX3(NDIM,nadj_rec_local,1,irec_local,it)]*lagrange);
      atomicAdd(&accel[2+3*iglob],source_adjoint[INDEX3(NDIM,nadj_rec_local,2,irec_local,it)]*lagrange);
    } // ispec_is_elastic
  }

}


