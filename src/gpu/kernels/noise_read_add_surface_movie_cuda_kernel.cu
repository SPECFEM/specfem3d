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


__global__ void noise_read_add_surface_movie_cuda_kernel(realw* accel, int* d_ibool,
                                                         int* free_surface_ispec,
                                                         int* free_surface_ijk,
                                                         int num_free_surface_faces,
                                                         realw* noise_surface_movie,
                                                         realw* normal_x_noise,
                                                         realw* normal_y_noise,
                                                         realw* normal_z_noise,
                                                         realw* mask_noise,
                                                         realw* free_surface_jacobian2Dw) {

  int iface = blockIdx.x + gridDim.x*blockIdx.y; // surface element id

  // when nspec_top > 65535, but mod(nspec_top,2) > 0, we end up with an extra block.
  if (iface < num_free_surface_faces) {
    int ispec = free_surface_ispec[iface]-1;

    int igll = threadIdx.x;

    int ipoin = NGLL2*iface + igll;
    int i = free_surface_ijk[INDEX3(NDIM,NGLL2,0,igll,iface)]-1;
    int j = free_surface_ijk[INDEX3(NDIM,NGLL2,1,igll,iface)]-1;
    int k = free_surface_ijk[INDEX3(NDIM,NGLL2,2,igll,iface)]-1;

    int iglob = d_ibool[INDEX4_PADDED(NGLLX,NGLLX,NGLLX,i,j,k,ispec)]-1;

    realw normal_x = normal_x_noise[ipoin];
    realw normal_y = normal_y_noise[ipoin];
    realw normal_z = normal_z_noise[ipoin];

    // along noise direction
    realw eta = (noise_surface_movie[INDEX3(NDIM,NGLL2,0,igll,iface)]*normal_x +
                 noise_surface_movie[INDEX3(NDIM,NGLL2,1,igll,iface)]*normal_y +
                 noise_surface_movie[INDEX3(NDIM,NGLL2,2,igll,iface)]*normal_z);

    realw val = eta * mask_noise[ipoin] * free_surface_jacobian2Dw[igll+NGLL2*iface];

    // error from cuda-memcheck and ddt seems "incorrect", because we
    // are passing a __constant__ variable pointer around like it was
    // made using cudaMalloc, which *may* be "incorrect", but produces
    // correct results.

    // ========= Invalid __global__ read of size
    // 4 ========= at 0x00000cd8 in
    // compute_add_sources_cuda.cu:260:noise_read_add_surface_movie_cuda_kernel
    // ========= by thread (0,0,0) in block (3443,0) ========= Address
    // 0x203000c8 is out of bounds

    // non atomic version for speed testing -- atomic updates are needed for correctness
    // accel[3*iglob] +=   eta*mask_noise[ipoin] * normal_x * wgllwgll_xy[tx] * free_surface_jacobian2Dw[tx + NGLL2*ispec2D];
    // accel[3*iglob+1] += eta*mask_noise[ipoin] * normal_y * wgllwgll_xy[tx] * free_surface_jacobian2Dw[tx + NGLL2*ispec2D];
    // accel[3*iglob+2] += eta*mask_noise[ipoin] * normal_z * wgllwgll_xy[tx] * free_surface_jacobian2Dw[tx + NGLL2*ispec2D];

    // Fortran version in SVN -- note deletion of wgllwgll_xy?
    // accel(1,iglob) = accel(1,iglob) + eta * mask_noise(ipoin) * normal_x_noise(ipoin) &
    // * free_surface_jacobian2Dw(igll,iface)
    // accel(2,iglob) = accel(2,iglob) + eta * mask_noise(ipoin) * normal_y_noise(ipoin) &
    // * free_surface_jacobian2Dw(igll,iface)
    // accel(3,iglob) = accel(3,iglob) + eta * mask_noise(ipoin) * normal_z_noise(ipoin) &
    // * free_surface_jacobian2Dw(igll,iface) ! wgllwgll_xy(i,j) * jacobian2D_top(i,j,iface)

    // atomicAdd(&accel[iglob*3]  ,eta*mask_noise[ipoin]*normal_x*wgllwgll_xy[tx]*free_surface_jacobian2Dw[igll+NGLL2*iface]);
    // atomicAdd(&accel[iglob*3+1],eta*mask_noise[ipoin]*normal_y*wgllwgll_xy[tx]*free_surface_jacobian2Dw[igll+NGLL2*iface]);
    // atomicAdd(&accel[iglob*3+2],eta*mask_noise[ipoin]*normal_z*wgllwgll_xy[tx]*free_surface_jacobian2Dw[igll+NGLL2*iface]);

    atomicAdd(&accel[iglob*3]  , val * normal_x);
    atomicAdd(&accel[iglob*3+1], val * normal_y);
    atomicAdd(&accel[iglob*3+2], val * normal_z);

  }
}


