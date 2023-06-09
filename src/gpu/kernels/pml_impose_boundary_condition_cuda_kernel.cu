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


__global__ void pml_impose_boundary_condition_cuda_kernel(realw* accel,
                                                          realw* veloc,
                                                          realw* displ,
                                                          realw* PML_displ_old,
                                                          realw* PML_displ_new,
                                                          int* abs_boundary_ispec,
                                                          int* abs_boundary_ijk,
                                                          int num_abs_boundary_faces,
                                                          int* d_ibool,
                                                          int* ispec_is_elastic,
                                                          int* is_CPML,
                                                          int* spec_to_CPML) {

  int igll = threadIdx.x; // tx
  int iface = blockIdx.x + gridDim.x*blockIdx.y; // bx

  int i,j,k,iglob,ispec,ispec_CPML;

  // don't compute points outside NGLLSQUARE==NGLL2==25
  // way 2: no further check needed since blocksize = 25
  if (iface < num_abs_boundary_faces){

  //if (igll < NGLL2 && iface < num_abs_boundary_faces) {

    // "-1" from index values to convert from Fortran-> C indexing
    ispec = abs_boundary_ispec[iface]-1;

    if (ispec_is_elastic[ispec] && is_CPML[ispec]) {

      i = abs_boundary_ijk[INDEX3(NDIM,NGLL2,0,igll,iface)]-1;
      j = abs_boundary_ijk[INDEX3(NDIM,NGLL2,1,igll,iface)]-1;
      k = abs_boundary_ijk[INDEX3(NDIM,NGLL2,2,igll,iface)]-1;

      iglob = d_ibool[INDEX4_PADDED(NGLLX,NGLLX,NGLLX,i,j,k,ispec)]-1;

      // gets associated velocity

      displ[iglob*3+0] = 0.f;
      displ[iglob*3+1] = 0.f;
      displ[iglob*3+2] = 0.f;

      veloc[iglob*3+0] = 0.f;
      veloc[iglob*3+1] = 0.f;
      veloc[iglob*3+2] = 0.f;

      accel[iglob*3+0] = 0.f;
      accel[iglob*3+1] = 0.f;
      accel[iglob*3+2] = 0.f;

      ispec_CPML = spec_to_CPML[ispec]-1;

      PML_displ_old[INDEX5(NDIM,NGLLX,NGLLX,NGLLX,0,i,j,k,ispec_CPML)] = 0.f;
      PML_displ_old[INDEX5(NDIM,NGLLX,NGLLX,NGLLX,1,i,j,k,ispec_CPML)] = 0.f;
      PML_displ_old[INDEX5(NDIM,NGLLX,NGLLX,NGLLX,2,i,j,k,ispec_CPML)] = 0.f;

      PML_displ_new[INDEX5(NDIM,NGLLX,NGLLX,NGLLX,0,i,j,k,ispec_CPML)] = 0.f;
      PML_displ_new[INDEX5(NDIM,NGLLX,NGLLX,NGLLX,1,i,j,k,ispec_CPML)] = 0.f;
      PML_displ_new[INDEX5(NDIM,NGLLX,NGLLX,NGLLX,2,i,j,k,ispec_CPML)] = 0.f;
    }
  } // num_abs_boundary_faces
}
