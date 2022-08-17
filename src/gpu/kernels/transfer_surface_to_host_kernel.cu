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


__global__ void transfer_surface_to_host_kernel(int* free_surface_ispec,
                                                int* free_surface_ijk,
                                                int num_free_surface_faces,
                                                int* d_ibool,
                                                realw* displ,
                                                realw* noise_surface_movie) {
  int igll = threadIdx.x;
  int iface = blockIdx.x + blockIdx.y*gridDim.x;

  // int id = tx + blockIdx.x*blockDim.x + blockIdx.y*blockDim.x*gridDim.x;

  if (iface < num_free_surface_faces) {
    int ispec = free_surface_ispec[iface]-1; //-1 for C-based indexing

    int i = free_surface_ijk[INDEX3(NDIM,NGLL2,0,igll,iface)]-1;
    int j = free_surface_ijk[INDEX3(NDIM,NGLL2,1,igll,iface)]-1;
    int k = free_surface_ijk[INDEX3(NDIM,NGLL2,2,igll,iface)]-1;

    int iglob = d_ibool[INDEX4_PADDED(NGLLX,NGLLX,NGLLX,i,j,k,ispec)]-1;

    noise_surface_movie[INDEX3(NDIM,NGLL2,0,igll,iface)] = displ[iglob*3];
    noise_surface_movie[INDEX3(NDIM,NGLL2,1,igll,iface)] = displ[iglob*3+1];
    noise_surface_movie[INDEX3(NDIM,NGLL2,2,igll,iface)] = displ[iglob*3+2];
  }
}


