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


__global__ void add_sources_ac_SIM_TYPE_2_OR_3_kernel(field* potential_dot_dot_acoustic,
                                                      int nrec,
                                                      int it,
                                                      int NSTEP_BETWEEN_ADJSRC,
                                                      field* source_adjoint,
                                                      realw* xir_store,
                                                      realw* etar_store,
                                                      realw* gammar_store,
                                                      int* d_ibool,
                                                      int* ispec_is_acoustic,
                                                      int* ispec_selected_recloc,
                                                      int nadj_rec_local,
                                                      realw* kappastore) {

  int irec_local = blockIdx.x + gridDim.x*blockIdx.y;

  // because of grid shape, irec_local can be too big
  if (irec_local < nadj_rec_local) {

    int ispec = ispec_selected_recloc[irec_local]-1;
    if (ispec_is_acoustic[ispec]){
      int i = threadIdx.x;
      int j = threadIdx.y;
      int k = threadIdx.z;

      int iglob = d_ibool[INDEX4_PADDED(NGLLX,NGLLX,NGLLX,i,j,k,ispec)]-1;

      realw hxir    = xir_store[INDEX2(NGLLX,i,irec_local)];
      realw hetar   = etar_store[INDEX2(NGLLX,j,irec_local)];
      realw hgammar = gammar_store[INDEX2(NGLLX,k,irec_local)];

      // note: we take the first component of the adj_sourcearrays
      field source_adj = source_adjoint[INDEX3(NDIM,nadj_rec_local,0,irec_local,it)];


      //realw kappal = kappastore[INDEX4(NGLLX,NGLLY,NGLLZ,i,j,k,ispec)];
      //
      //potential_dot_dot_acoustic[iglob] += adj_sourcearrays[INDEX6(nadj_rec_local,NTSTEP_BETWEEN_ADJSRC,3,5,5,
      //                                            pre_computed_irec_local_index[irec],
      //                                            pre_computed_index,
      //                                            0,
      //                                            i,j,k)]/kappal;

      // beware, for acoustic medium, a pressure source would be taking the negative
      // and divide by Kappa of the fluid;

      //realw stf = - source_adj * hxir * hetar * hgammar / kappal;

      // VM VM : change the adjoint source to be consistent with CPU code
      field stf = source_adj * hxir * hetar * hgammar;

      atomicAdd(&potential_dot_dot_acoustic[iglob],stf);

                //+adj_sourcearrays[INDEX6(nadj_rec_local,NTSTEP_BETWEEN_ADJSRC,3,5,5,
                //                         pre_computed_irec_local_index[irec],pre_computed_index-1,
                //                         0,i,j,k)] // / kappal
                //                         );
    }
  }
}



