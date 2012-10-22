/*
 !=====================================================================
 !
 !               S p e c f e m 3 D  V e r s i o n  2 . 1
 !               ---------------------------------------
 !
 !          Main authors: Dimitri Komatitsch and Jeroen Tromp
 !    Princeton University, USA and CNRS / INRIA / University of Pau
 ! (c) Princeton University / California Institute of Technology and CNRS / INRIA / University of Pau
 !                             July 2012
 !
 ! This program is free software; you can redistribute it and/or modify
 ! it under the terms of the GNU General Public License as published by
 ! the Free Software Foundation; either version 2 of the License, or
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

#include <stdio.h>
#include <cuda.h>
#include <cublas.h>

#ifdef WITH_MPI
#include <mpi.h>
#endif

#include <sys/time.h>
#include <sys/resource.h>

#include "config.h"
#include "mesh_constants_cuda.h"


/* ----------------------------------------------------------------------------------------------- */

// GPU initialization

/* ----------------------------------------------------------------------------------------------- */

extern "C"
void FC_FUNC_(initialize_cuda_device,
              INITIALIZE_CUDA_DEVICE)(int* myrank_f,int* ncuda_devices) {
  TRACE("initialize_cuda_device");

  // Gets rank number of MPI process
  int myrank = *myrank_f;

  /*
   // cuda initialization (needs -lcuda library)
   // note:   cuInit initializes the driver API.
   //             it is needed for any following CUDA driver API function call (format cuFUNCTION(..) )
   //             however, for the CUDA runtime API functions (format cudaFUNCTION(..) )
   //             the initialization is implicit, thus cuInit() here would not be needed...
   CUresult status = cuInit(0);
   if ( CUDA_SUCCESS != status ) exit_on_error("CUDA driver API device initialization failed\n");

   // returns a handle to the first cuda compute device
   CUdevice dev;
   status = cuDeviceGet(&dev, 0);
   if ( CUDA_SUCCESS != status ) exit_on_error("CUDA device not found\n");

   // gets device properties
   int major,minor;
   status = cuDeviceComputeCapability(&major,&minor,dev);
   if ( CUDA_SUCCESS != status ) exit_on_error("CUDA device information not found\n");

   // make sure that the device has compute capability >= 1.3
   if (major < 1){
   fprintf(stderr,"Compute capability major number should be at least 1, got: %d \nexiting...\n",major);
   exit_on_error("CUDA Compute capability major number should be at least 1\n");
   }
   if (major == 1 && minor < 3){
   fprintf(stderr,"Compute capability should be at least 1.3, got: %d.%d \nexiting...\n",major,minor);
   exit_on_error("CUDA Compute capability major number should be at least 1.3\n");
   }
   */

  // note: from here on we use the runtime API  ...

  // Gets number of GPU devices
  int device_count = 0;
  cudaGetDeviceCount(&device_count);
  exit_on_cuda_error("CUDA runtime error: cudaGetDeviceCount failed\ncheck if driver and runtime libraries work together\nexiting...\n");

  // returns device count to fortran
  if (device_count == 0) exit_on_error("CUDA runtime error: there is no device supporting CUDA\n");
  *ncuda_devices = device_count;

  // Sets the active device
  if(device_count >= 1) {
    // generalized for more GPUs per node
    // note: without previous context release, cudaSetDevice will complain with the cuda error
    //         "setting the device when a process is active is not allowed"
    // releases previous contexts
    cudaThreadExit();

    //printf("rank %d: cuda device count = %d sets device = %d \n",myrank,device_count,myrank % device_count);
    //MPI_Barrier(MPI_COMM_WORLD);

    // sets active device
    cudaSetDevice( myrank % device_count );
    exit_on_cuda_error("cudaSetDevice");
  }

  // returns a handle to the active device
  int device;
  cudaGetDevice(&device);

  // get device properties
  struct cudaDeviceProp deviceProp;
  cudaGetDeviceProperties(&deviceProp,device);

  // exit if the machine has no CUDA-enabled device
  if (deviceProp.major == 9999 && deviceProp.minor == 9999){
    fprintf(stderr,"No CUDA-enabled device found, exiting...\n\n");
    exit_on_error("CUDA runtime error: there is no CUDA-enabled device found\n");
  }

  // outputs device infos to file
  char filename[BUFSIZ];
  FILE* fp;
  sprintf(filename,"../OUTPUT_FILES/gpu_device_info_proc_%06d.txt",myrank);
  fp = fopen(filename,"a+");
  if (fp != NULL){
    // display device properties
    fprintf(fp,"Device Name = %s\n",deviceProp.name);
    fprintf(fp,"multiProcessorCount: %d\n",deviceProp.multiProcessorCount);
    fprintf(fp,"totalGlobalMem (in MB): %f\n",(unsigned long) deviceProp.totalGlobalMem / (1024.f * 1024.f));
    fprintf(fp,"totalGlobalMem (in GB): %f\n",(unsigned long) deviceProp.totalGlobalMem / (1024.f * 1024.f * 1024.f));
    fprintf(fp,"sharedMemPerBlock (in bytes): %lu\n",(unsigned long) deviceProp.sharedMemPerBlock);
    fprintf(fp,"Maximum number of threads per block: %d\n",deviceProp.maxThreadsPerBlock);
    fprintf(fp,"Maximum size of each dimension of a block: %d x %d x %d\n",
            deviceProp.maxThreadsDim[0],deviceProp.maxThreadsDim[1],deviceProp.maxThreadsDim[2]);
    fprintf(fp,"Maximum sizes of each dimension of a grid: %d x %d x %d\n",
            deviceProp.maxGridSize[0],deviceProp.maxGridSize[1],deviceProp.maxGridSize[2]);
    fprintf(fp,"Compute capability of the device = %d.%d\n", deviceProp.major, deviceProp.minor);
    if(deviceProp.canMapHostMemory){
      fprintf(fp,"canMapHostMemory: TRUE\n");
    }else{
      fprintf(fp,"canMapHostMemory: FALSE\n");
    }
    if(deviceProp.deviceOverlap){
      fprintf(fp,"deviceOverlap: TRUE\n");
    }else{
      fprintf(fp,"deviceOverlap: FALSE\n");
    }

    // outputs initial memory infos via cudaMemGetInfo()
    double free_db,used_db,total_db;
    get_free_memory(&free_db,&used_db,&total_db);
    fprintf(fp,"%d: GPU memory usage: used = %f MB, free = %f MB, total = %f MB\n", myrank,
            used_db/1024.0/1024.0, free_db/1024.0/1024.0, total_db/1024.0/1024.0);

    fclose(fp);
  }

  // make sure that the device has compute capability >= 1.3
  if (deviceProp.major < 1){
    fprintf(stderr,"Compute capability major number should be at least 1, exiting...\n\n");
    exit_on_error("CUDA Compute capability major number should be at least 1\n");
  }
  if (deviceProp.major == 1 && deviceProp.minor < 3){
    fprintf(stderr,"Compute capability should be at least 1.3, exiting...\n");
    exit_on_error("CUDA Compute capability major number should be at least 1.3\n");
  }
  // we use pinned memory for asynchronous copy
  if( ! deviceProp.canMapHostMemory){
    fprintf(stderr,"Device capability should allow to map host memory, exiting...\n");
    exit_on_error("CUDA Device capability canMapHostMemory should be TRUE\n");
  }
}
