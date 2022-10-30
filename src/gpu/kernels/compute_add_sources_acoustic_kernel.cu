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


__global__ void compute_add_sources_acoustic_kernel(field* potential_dot_dot_acoustic,
                                                    int* d_ibool,
                                                    realw* sourcearrays,
                                                    field* stf_pre_compute,
                                                    int myrank,
                                                    int* islice_selected_source,
                                                    int* ispec_selected_source,
                                                    int* ispec_is_acoustic,
                                                    realw* kappastore,
                                                    int NSOURCES) {
  int i = threadIdx.x;
  int j = threadIdx.y;
  int k = threadIdx.z;

  int isource  = blockIdx.x + gridDim.x*blockIdx.y; // bx

  int ispec,iglob;
  field stf;
  realw kappal;

  if (isource < NSOURCES){

    if (myrank == islice_selected_source[isource]) {

      ispec = ispec_selected_source[isource]-1;

      if (ispec_is_acoustic[ispec]) {

        iglob = d_ibool[INDEX4_PADDED(NGLLX,NGLLX,NGLLX,i,j,k,ispec)] - 1;

        stf = stf_pre_compute[isource];
        kappal = kappastore[INDEX4(NGLLX,NGLLX,NGLLX,i,j,k,ispec)];

        atomicAdd(&potential_dot_dot_acoustic[iglob],
                  -sourcearrays[INDEX5(NSOURCES,NDIM,NGLLX,NGLLX,isource, 0,i,j,k)]*stf/kappal);

        // debug: without atomic operation
        //      potential_dot_dot_acoustic[iglob] +=
        //              -sourcearrays[INDEX5(NSOURCES, 3, NGLLX,NGLLX,isource, 0, i,j,k)]*stf/kappal;
      }
    }
  }
}



