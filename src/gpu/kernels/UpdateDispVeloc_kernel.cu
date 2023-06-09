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


__global__ void UpdateDispVeloc_kernel(realw* displ,
                                       realw* veloc,
                                       realw* accel,
                                       const int size,
                                       const realw deltat,
                                       const realw deltatsqover2,
                                       const realw deltatover2) {

  // two dimensional array of blocks on grid where each block has one dimensional array of threads
  int id = threadIdx.x + (blockIdx.x + blockIdx.y*gridDim.x)*blockDim.x;

  // because of block and grid sizing problems, there is a small
  // amount of buffer at the end of the calculation
  if (id < size) {
    realw acc = accel[id];
    realw vel = veloc[id];

    displ[id] = displ[id] + deltat * vel + deltatsqover2 * acc;
    veloc[id] = vel + deltatover2 * acc;
    accel[id] = 0.0f; // can do this using memset...not sure if faster,probably not
  }

// -----------------
// total of: 6 FLOP per thread (without int id calculation at beginning)
//
//           8 * 4 BYTE = 32 DRAM accesses per thread
//
// arithmetic intensity: 6 FLOP / 32 BYTES ~ 0.19 FLOP/BYTE
// -----------------
// nvprof: 24599250 flops for 4099875 threads -> 6 FLOP per thread
}

/* ----------------------------------------------------------------------------------------------- */


__global__ void UpdateDispVeloc_PML_kernel(realw* displ,
                                           realw* veloc,
                                           realw* accel,
                                           realw* PML_displ,
                                           const int NSPEC_CPML,
                                           const int* d_CPML_to_spec,
                                           const int* d_ibool,
                                           const realw deltat,
                                           const realw deltatsqover2,
                                           const realw deltatover2) {

  int ispec_cpml = blockIdx.x + blockIdx.y*gridDim.x;
  int ijk = threadIdx.x;

  // because of block and grid sizing problems, there is a small
  // amount of buffer at the end of the calculation
  if (ispec_cpml < NSPEC_CPML) {

    int ispec = d_CPML_to_spec[ispec_cpml] - 1;

    // local and global indices
    int K = ijk/NGLL2;
    int J = (ijk - K*NGLL2)/NGLLX;
    int I = ijk - K*NGLL2 - J*NGLLX;

    int iglob = d_ibool[ijk + NGLL3_PADDED*ispec] - 1;

    const realw theta = 0.125f; // theta = 1.0 / 8.0;

    // updates PML displacement
    PML_displ[INDEX5(NDIM,NGLLX,NGLLX,NGLLX,0,I,J,K,ispec_cpml)] = displ[3*iglob]
                                                                   + deltatover2 * (1.0f - 2.0f * theta) * veloc[3*iglob]
                                                                   + deltatsqover2 * (1.0f - theta) * accel[3*iglob];

    PML_displ[INDEX5(NDIM,NGLLX,NGLLX,NGLLX,1,I,J,K,ispec_cpml)] = displ[3*iglob+1]
                                                                   + deltatover2 * (1.0f - 2.0f * theta) * veloc[3*iglob+1]
                                                                   + deltatsqover2 * (1.0f - theta) * accel[3*iglob+1];

    PML_displ[INDEX5(NDIM,NGLLX,NGLLX,NGLLX,2,I,J,K,ispec_cpml)] = displ[3*iglob+2]
                                                                   + deltatover2 * (1.0f - 2.0f * theta) * veloc[3*iglob+2]
                                                                   + deltatsqover2 * (1.0f - theta) * accel[3*iglob+2];
  }
}
