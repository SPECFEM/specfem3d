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


__global__ void get_maximum_field_kernel(field* array, int size, realw* d_max){

  /* simplest version: uses only 1 thread
   realw max;
   max = 0;
   // finds maximum value in array
   if (size > 0){
   max = abs(array[0]);
   for( int i=1; i < size; i++){
   if (abs(array[i]) > max) max = abs(array[i]);
   }
   }
   *d_max = max;
   */

  // reduction example:
  __shared__ realw sdata[BLOCKSIZE_TRANSFER] ;

  // load shared mem
  unsigned int tid = threadIdx.x;
  unsigned int bx = blockIdx.y*gridDim.x+blockIdx.x;
  unsigned int i = tid + bx*blockDim.x;

  // loads absolute values into shared memory
  sdata[tid] = (i < size) ? fabs(array[i]) : 0.0 ;

  __syncthreads();

  // do reduction in shared mem
  for(unsigned int s=blockDim.x/2; s>0; s>>=1)
  {
    if (tid < s){
      // summation:
      //sdata[tid] += sdata[tid + s];
      // maximum:
      if (sdata[tid] < sdata[tid + s]) sdata[tid] = sdata[tid + s];
    }
    __syncthreads();
  }

  // write result for this block to global mem
  if (tid == 0) d_max[bx] = sdata[0];

}

