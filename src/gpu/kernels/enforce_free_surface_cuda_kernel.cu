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


__global__ void enforce_free_surface_cuda_kernel(field_p potential_acoustic,
                                                 field_p potential_dot_acoustic,
                                                 field_p potential_dot_dot_acoustic,
                                                 const int num_free_surface_faces,
                                                 const int* free_surface_ispec,
                                                 const int* free_surface_ijk,
                                                 const int* d_ibool,
                                                 const int* ispec_is_acoustic) {
  // gets spectral element face id
  int iface = blockIdx.x + gridDim.x*blockIdx.y;

  // for all faces on free surface
  if (iface < num_free_surface_faces ){

    int ispec = free_surface_ispec[iface]-1;

    // checks if element is in acoustic domain
    if (ispec_is_acoustic[ispec] ){

      // gets global point index
      int igll = threadIdx.x + threadIdx.y*blockDim.x;

      int i = free_surface_ijk[INDEX3(NDIM,NGLL2,0,igll,iface)] - 1; // (1,igll,iface)
      int j = free_surface_ijk[INDEX3(NDIM,NGLL2,1,igll,iface)] - 1;
      int k = free_surface_ijk[INDEX3(NDIM,NGLL2,2,igll,iface)] - 1;

      int iglob = d_ibool[INDEX4_PADDED(NGLLX,NGLLX,NGLLX,i,j,k,ispec)] - 1;

      // sets potentials to zero at free surface
      potential_acoustic[iglob] = Make_field(0.f);
      potential_dot_acoustic[iglob] = Make_field(0.f);
      potential_dot_dot_acoustic[iglob] = Make_field(0.f);
    }
  }
}


