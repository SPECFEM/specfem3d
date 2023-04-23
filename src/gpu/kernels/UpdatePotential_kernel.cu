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


__global__ void UpdatePotential_kernel(field* potential_acoustic,
                                       field* potential_dot_acoustic,
                                       field* potential_dot_dot_acoustic,
                                       int size,
                                       realw deltat,
                                       realw deltatsqover2,
                                       realw deltatover2) {

  int id = threadIdx.x + (blockIdx.x + blockIdx.y*gridDim.x)*blockDim.x;

  // because of block and grid sizing problems, there is a small
  // amount of buffer at the end of the calculation
  if (id < size) {
    field p_dot = potential_dot_acoustic[id];
    field p_dot_dot = potential_dot_dot_acoustic[id];

    potential_acoustic[id] += deltat*p_dot + deltatsqover2*p_dot_dot;

    potential_dot_acoustic[id] = p_dot + deltatover2*p_dot_dot;

    potential_dot_dot_acoustic[id] = Make_field(0.f);
  }

// -----------------
// total of: 6 FLOP per thread (without id calculation)
//
//           8 * 4 BYTE = 32 DRAM accesses per thread
//
// arithmetic intensity: 6 FLOP / 32 BYTES ~ 0.19 FLOP/BYTE
// -----------------
//
// nvprof: nvprof --metrics flops_sp ./xspecfem3D
//          -> 8199750 FLOPS (Single) floating-point operations for 1366625 threads
//                                    1366625 (NGLOB) -> 10677 * 128 active threads- 31 ghost threads
//          -> 6 FLOP per thread


}

