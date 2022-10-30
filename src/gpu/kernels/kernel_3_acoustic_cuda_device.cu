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


__global__ void kernel_3_acoustic_cuda_device(field* potential_dot_acoustic,
                                                field* potential_dot_dot_acoustic,
                                                field* b_potential_dot_acoustic,
                                                field* b_potential_dot_dot_acoustic,
                                                int simulation_type,
                                                int size,
                                                realw deltatover2,
                                                realw b_deltatover2,
                                                realw_const_p rmass_acoustic) {

  int id = threadIdx.x + (blockIdx.x + blockIdx.y*gridDim.x)*blockDim.x;

  realw rmass;
  field p_dot_dot;
  // because of block and grid sizing problems, there is a small
  // amount of buffer at the end of the calculation
  if (id < size) {
    rmass = rmass_acoustic[id];
    // multiplies pressure with the inverse of the mass matrix
    p_dot_dot = rmass*potential_dot_dot_acoustic[id];
    potential_dot_dot_acoustic[id] = p_dot_dot;
    potential_dot_acoustic[id] += deltatover2*p_dot_dot;
    if (simulation_type==3) {
      p_dot_dot = rmass*b_potential_dot_dot_acoustic[id];
      b_potential_dot_dot_acoustic[id] = p_dot_dot;
      b_potential_dot_acoustic[id] += b_deltatover2*p_dot_dot;
    }
  }
}

/* ----------------------------------------------------------------------------------------------- */

__global__ void kernel_3_acoustic_single_cuda_device(field* potential_dot_acoustic,
                                                     field* potential_dot_dot_acoustic,
                                                     int size,
                                                     realw deltatover2,
                                                     realw_const_p rmass_acoustic) {

  int id = threadIdx.x + (blockIdx.x + blockIdx.y*gridDim.x)*blockDim.x;

  // because of block and grid sizing problems, there is a small
  // amount of buffer at the end of the calculation
  if (id < size) {
    realw rmass = rmass_acoustic[id];
    // multiplies pressure with the inverse of the mass matrix
    field p_dot_dot = rmass * potential_dot_dot_acoustic[id];

    potential_dot_dot_acoustic[id] = p_dot_dot;
    potential_dot_acoustic[id] += deltatover2 * p_dot_dot;
  }
}


