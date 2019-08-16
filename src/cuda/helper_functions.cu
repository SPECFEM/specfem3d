/*
 !=====================================================================
 !
 !               S p e c f e m 3 D  V e r s i o n  3 . 0
 !               ---------------------------------------
 !
 !     Main historical authors: Dimitri Komatitsch and Jeroen Tromp
 !                              CNRS, France
 !                       and Princeton University, USA
 !                 (there are currently many more authors!)
 !                           (c) October 2017
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

#include "mesh_constants_cuda.h"

/* ----------------------------------------------------------------------------------------------- */

// helper functions

/* ----------------------------------------------------------------------------------------------- */


// copies integer array from CPU host to GPU device
void copy_todevice_int(void** d_array_addr_ptr,int* h_array,int size){
  TRACE("  copy_todevice_int");

  // allocates memory on GPU
  //
  // note: cudaMalloc uses a double-pointer, such that it can return an error code in case it fails
  //          we thus pass the address to the pointer above (as void double-pointer) to have it
  //          pointing to the correct pointer of the array here
  print_CUDA_error_if_any(cudaMalloc((void**)d_array_addr_ptr,size*sizeof(int)),
                          12001);

  // copies values onto GPU
  //
  // note: cudaMemcpy uses the pointer to the array, we thus re-cast the value of
  //          the double-pointer above to have the correct pointer to the array
  print_CUDA_error_if_any(cudaMemcpy((int*) *d_array_addr_ptr,h_array,size*sizeof(int),cudaMemcpyHostToDevice),
                          12002);
}

/* ----------------------------------------------------------------------------------------------- */

// copies integer array from CPU host to GPU device
void copy_todevice_realw(void** d_array_addr_ptr,realw* h_array,int size){
  TRACE("  copy_todevice_realw");

  // allocates memory on GPU
  print_CUDA_error_if_any(cudaMalloc((void**)d_array_addr_ptr,size*sizeof(realw)),
                          22001);
  // copies values onto GPU
  print_CUDA_error_if_any(cudaMemcpy((realw*) *d_array_addr_ptr,h_array,size*sizeof(realw),cudaMemcpyHostToDevice),
                          22002);
}
