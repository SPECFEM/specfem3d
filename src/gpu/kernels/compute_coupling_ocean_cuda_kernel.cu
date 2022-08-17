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


__global__ void compute_coupling_ocean_cuda_kernel(realw* accel,
                                                   realw* rmassx,realw* rmassy,realw* rmassz,
                                                   realw* rmass_ocean_load,
                                                   int num_free_surface_faces,
                                                   int* free_surface_ispec,
                                                   int* free_surface_ijk,
                                                   realw* free_surface_normal,
                                                   int* d_ibool,
                                                   int* updated_dof_ocean_load) {
  // gets spectral element face id
  int igll = threadIdx.x ;  //  threadIdx.y*blockDim.x will be always = 0 for thread block (25,1,1)
  int iface = blockIdx.x + gridDim.x*blockIdx.y;

  realw nx,ny,nz;
  realw force_normal_comp;

  // for all faces on free surface
  if (iface < num_free_surface_faces ){

    int ispec = free_surface_ispec[iface]-1;

    // gets global point index
    int i = free_surface_ijk[INDEX3(NDIM,NGLL2,0,igll,iface)] - 1; // (1,igll,iface)
    int j = free_surface_ijk[INDEX3(NDIM,NGLL2,1,igll,iface)] - 1;
    int k = free_surface_ijk[INDEX3(NDIM,NGLL2,2,igll,iface)] - 1;

    int iglob = d_ibool[INDEX4_PADDED(NGLLX,NGLLX,NGLLX,i,j,k,ispec)] - 1;

    //if (igll == 0) printf("igll %d %d %d %d\n",igll,i,j,k,iglob);

    // only update this global point once

    // daniel: TODO - there might be better ways to implement a mutex like below,
    //            and find a workaround to not use the temporary update array.
    //            atomicExch: returns the old value, i.e. 0 indicates that we still have to do this point

    if (atomicExch(&updated_dof_ocean_load[iglob],1) == 0){

      // get normal
      nx = free_surface_normal[INDEX3(NDIM,NGLL2,0,igll,iface)]; //(1,igll,iface)
      ny = free_surface_normal[INDEX3(NDIM,NGLL2,1,igll,iface)];
      nz = free_surface_normal[INDEX3(NDIM,NGLL2,2,igll,iface)];

      // make updated component of right-hand side
      // we divide by rmass() which is 1 / M
      // we use the total force which includes the Coriolis term above
      force_normal_comp = accel[iglob*3]*nx / rmassx[iglob]
                          + accel[iglob*3+1]*ny / rmassy[iglob]
                          + accel[iglob*3+2]*nz / rmassz[iglob];

      // probably wouldn't need atomicAdd anymore, but just to be sure...
      atomicAdd(&accel[iglob*3],   + (rmass_ocean_load[iglob] - rmassx[iglob]) * force_normal_comp * nx);
      atomicAdd(&accel[iglob*3+1], + (rmass_ocean_load[iglob] - rmassy[iglob]) * force_normal_comp * ny);
      atomicAdd(&accel[iglob*3+2], + (rmass_ocean_load[iglob] - rmassz[iglob]) * force_normal_comp * nz);
    }
  }
}


