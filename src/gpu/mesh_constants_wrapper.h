/*
 !=====================================================================
 !
 !                         S p e c f e m 3 D
 !                         -----------------
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

#ifndef MESH_CONSTANTS_WRAPPER_H
#define MESH_CONSTANTS_WRAPPER_H


/* ----------------------------------------------------------------------------------------------- */

// wrapper functions

/* ----------------------------------------------------------------------------------------------- */


static inline void gpuMemcpy_todevice_realw(realw* d_array,realw* h_array,const size_t size){
  // copies array onto GPU
#ifdef USE_CUDA
  if (run_cuda){
    print_CUDA_error_if_any(cudaMemcpy(d_array,h_array,size*sizeof(realw),cudaMemcpyHostToDevice),1801);
  }
#endif
#ifdef USE_HIP
  if (run_hip){
    print_HIP_error_if_any(hipMemcpy(d_array,h_array,size*sizeof(realw),hipMemcpyHostToDevice),1801);
  }
#endif
}

static inline void gpuMemcpy_todevice_field(field* d_array,field* h_array,const size_t size){
  // copies array onto GPU
#ifdef USE_CUDA
  if (run_cuda){
    print_CUDA_error_if_any(cudaMemcpy(d_array,h_array,size*sizeof(field),cudaMemcpyHostToDevice),1802);
  }
#endif
#ifdef USE_HIP
  if (run_hip){
    print_HIP_error_if_any(hipMemcpy(d_array,h_array,size*sizeof(field),hipMemcpyHostToDevice),1802);
  }
#endif
}

static inline void gpuMemcpy_todevice_int(int* d_array,int* h_array,const size_t size){
  // copies array onto GPU
#ifdef USE_CUDA
  if (run_cuda){
    print_CUDA_error_if_any(cudaMemcpy(d_array,h_array,size*sizeof(int),cudaMemcpyHostToDevice),1803);
  }
#endif
#ifdef USE_HIP
  if (run_hip){
    print_HIP_error_if_any(hipMemcpy(d_array,h_array,size*sizeof(int),hipMemcpyHostToDevice),1803);
  }
#endif
}

static inline void gpuMemcpy_todevice_void(void* d_array,void* h_array,const size_t byte_size){
  // copies array onto GPU
#ifdef USE_CUDA
  if (run_cuda){
    print_CUDA_error_if_any(cudaMemcpy(d_array,h_array,byte_size,cudaMemcpyHostToDevice),1803);
  }
#endif
#ifdef USE_HIP
  if (run_hip){
    print_HIP_error_if_any(hipMemcpy(d_array,h_array,byte_size,hipMemcpyHostToDevice),1803);
  }
#endif
}


/* ----------------------------------------------------------------------------------------------- */


static inline void gpuMemcpy2D_todevice_realw(realw* d_array,const size_t d_size,
                                              realw* h_array,const size_t h_size,
                                              const size_t h_width,const size_t h_height){
  // copies array onto GPU
#ifdef USE_CUDA
  if (run_cuda){
    print_CUDA_error_if_any(cudaMemcpy2D(d_array,d_size*sizeof(realw),h_array,h_size*sizeof(realw),h_width*sizeof(realw),h_height,
                                         cudaMemcpyHostToDevice),1901);
  }
#endif
#ifdef USE_HIP
  if (run_hip){
    print_HIP_error_if_any(hipMemcpy2D(d_array,d_size*sizeof(realw),h_array,h_size*sizeof(realw),h_width*sizeof(realw),h_height,
                                       hipMemcpyHostToDevice),1901);
  }
#endif
}

static inline void gpuMemcpy2D_todevice_int(int* d_array,const size_t d_size,
                                            int* h_array,const size_t h_size,
                                            const size_t h_width,const size_t h_height){
  // copies array onto GPU
#ifdef USE_CUDA
  if (run_cuda){
    print_CUDA_error_if_any(cudaMemcpy2D(d_array,d_size*sizeof(int),h_array,h_size*sizeof(int),h_width*sizeof(int),h_height,
                                         cudaMemcpyHostToDevice),1902);
  }
#endif
#ifdef USE_HIP
  if (run_hip){
    print_HIP_error_if_any(hipMemcpy2D(d_array,d_size*sizeof(int),h_array,h_size*sizeof(int),h_width*sizeof(int),h_height,
                                       hipMemcpyHostToDevice),1902);
  }
#endif
}


/* ----------------------------------------------------------------------------------------------- */


static inline void gpuMemcpy_tohost_realw(realw* h_array,realw* d_array,const size_t size){
  // copies array to CPU
#ifdef USE_CUDA
  if (run_cuda){
    print_CUDA_error_if_any(cudaMemcpy(h_array,d_array,size*sizeof(realw),cudaMemcpyDeviceToHost),2001);
  }
#endif
#ifdef USE_HIP
  if (run_hip){
    print_HIP_error_if_any(hipMemcpy(h_array,d_array,size*sizeof(realw),hipMemcpyDeviceToHost),2001);
  }
#endif
}

static inline void gpuMemcpy_tohost_field(field* h_array,field* d_array,const size_t size){
  // copies array to CPU
#ifdef USE_CUDA
  if (run_cuda){
    print_CUDA_error_if_any(cudaMemcpy(h_array,d_array,size*sizeof(field),cudaMemcpyDeviceToHost),2002);
  }
#endif
#ifdef USE_HIP
  if (run_hip){
    print_HIP_error_if_any(hipMemcpy(h_array,d_array,size*sizeof(field),hipMemcpyDeviceToHost),2002);
  }
#endif
}

/* unused so far...
static inline void gpuMemcpy_tohost_int(int* h_array,int* d_array,const size_t size){
  // copies array to CPU
#ifdef USE_CUDA
  if (run_cuda){
    print_CUDA_error_if_any(cudaMemcpy(h_array,d_array,size*sizeof(int),cudaMemcpyDeviceToHost),2003);
  }
#endif
#ifdef USE_HIP
  if (run_hip){
    print_HIP_error_if_any(hipMemcpy(h_array,d_array,size*sizeof(int),hipMemcpyDeviceToHost),2003);
  }
#endif
}
*/

static inline void gpuMemcpy_tohost_void(void* h_array,void* d_array,const size_t byte_size){
  // copies array onto GPU
#ifdef USE_CUDA
  if (run_cuda){
    print_CUDA_error_if_any(cudaMemcpy(h_array,d_array,byte_size,cudaMemcpyDeviceToHost),2004);
  }
#endif
#ifdef USE_HIP
  if (run_hip){
    print_HIP_error_if_any(hipMemcpy(h_array,d_array,byte_size,hipMemcpyDeviceToHost),2004);
  }
#endif
}

/* ----------------------------------------------------------------------------------------------- */


static inline void gpuMemcpyAsync_todevice_realw(realw* d_array,realw* h_array,const size_t size, gpu_stream stream){
  // asynchronuous copy of array onto GPU
#ifdef USE_CUDA
  if (run_cuda){
    cudaMemcpyAsync(d_array,h_array,size*sizeof(realw),cudaMemcpyHostToDevice,stream);
  }
#endif
#ifdef USE_HIP
  if (run_hip){
    hipMemcpyAsync(d_array,h_array,size*sizeof(realw),hipMemcpyHostToDevice,stream);
  }
#endif
}


static inline void gpuMemcpyAsync_tohost_realw(realw* h_array,realw* d_array,const size_t size, gpu_stream stream){
  // asynchronuous copy of array to CPU
#ifdef USE_CUDA
  if (run_cuda){
    cudaMemcpyAsync(h_array,d_array,size*sizeof(realw),cudaMemcpyDeviceToHost,stream);
  }
#endif
#ifdef USE_HIP
  if (run_hip){
    hipMemcpyAsync(h_array,d_array,size*sizeof(realw),hipMemcpyDeviceToHost,stream);
  }
#endif
}



/* ----------------------------------------------------------------------------------------------- */


static inline void gpuMemset_int(int* d_array, const size_t size, int value){
  // sets value for array on device
#ifdef USE_CUDA
  if (run_cuda){
    print_CUDA_error_if_any(cudaMemset(d_array,value,size*sizeof(int)),2301);
  }
#endif
#ifdef USE_HIP
  if (run_hip){
    print_HIP_error_if_any(hipMemset(d_array,value,size*sizeof(int)),2301);
  }
#endif
}

static inline void gpuMemset_realw(realw* d_array, const size_t size, int value){
  // sets value for array on device
#ifdef USE_CUDA
  if (run_cuda){
    print_CUDA_error_if_any(cudaMemset(d_array, value, size*sizeof(realw)),2302);
  }
#endif
#ifdef USE_HIP
  if (run_hip){
    print_HIP_error_if_any(hipMemset(d_array, value, size*sizeof(realw)),2302);
  }
#endif
}


/* ----------------------------------------------------------------------------------------------- */


static inline void gpuMalloc_int(void** d_array_addr_ptr,const size_t size){
  // allocates array on device
#ifdef USE_CUDA
  if (run_cuda){
    print_CUDA_error_if_any(cudaMalloc((void**)d_array_addr_ptr,size*sizeof(int)),2401);
  }
#endif
#ifdef USE_HIP
  if (run_hip){
    print_HIP_error_if_any(hipMalloc((void**)d_array_addr_ptr,size*sizeof(int)),2401);
  }
#endif
}


static inline void gpuMalloc_realw(void** d_array_addr_ptr,const size_t size){
  // allocates array on device
#ifdef USE_CUDA
  if (run_cuda){
    print_CUDA_error_if_any(cudaMalloc((void**)d_array_addr_ptr,size*sizeof(realw)),2402);
  }
#endif
#ifdef USE_HIP
  if (run_hip){
    print_HIP_error_if_any(hipMalloc((void**)d_array_addr_ptr,size*sizeof(realw)),2402);
  }
#endif
}

static inline void gpuMalloc_field(void** d_array_addr_ptr,const size_t size){
  // allocates array on device
#ifdef USE_CUDA
  if (run_cuda){
    print_CUDA_error_if_any(cudaMalloc((void**)d_array_addr_ptr,size*sizeof(field)),2403);
  }
#endif
#ifdef USE_HIP
  if (run_hip){
    print_HIP_error_if_any(hipMalloc((void**)d_array_addr_ptr,size*sizeof(field)),2403);
  }
#endif
}

static inline void gpuMallocHost_realw(void** h_array_addr_ptr,const size_t size){
  // allocates array on device
#ifdef USE_CUDA
  if (run_cuda){
    print_CUDA_error_if_any(cudaMallocHost((void**)h_array_addr_ptr,size*sizeof(realw)),2404);
  }
#endif
#ifdef USE_HIP
  if (run_hip){
    print_HIP_error_if_any(hipHostMalloc((void**)h_array_addr_ptr,size*sizeof(realw)),2404);
  }
#endif
}



/* ----------------------------------------------------------------------------------------------- */


static inline void gpuFree(realw* d_array){
  // allocates array on device
#ifdef USE_CUDA
  if (run_cuda){ cudaFree(d_array); }
#endif
#ifdef USE_HIP
  if (run_hip){ hipFree(d_array); }
#endif
}

static inline void gpuFree(int* d_array){
  // allocates array on device
#ifdef USE_CUDA
  if (run_cuda){ cudaFree(d_array); }
#endif
#ifdef USE_HIP
  if (run_hip){ hipFree(d_array); }
#endif
}


#endif  // MESH_CONSTANTS_WRAPPER_H
