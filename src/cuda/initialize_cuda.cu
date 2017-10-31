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

// GPU initialization

/* ----------------------------------------------------------------------------------------------- */

extern "C"
void FC_FUNC_(initialize_cuda_device,
              INITIALIZE_CUDA_DEVICE)(int* myrank_f,int* ncuda_devices) {
  TRACE("initialize_cuda_device");

  int device;
  int device_count;

  // Gets rank number of MPI process
  int myrank = *myrank_f;

  /*
   // cuda initialization (needs -lcuda library)
   // note:   cuInit initializes the driver API.
   //             it is needed for any following CUDA driver API function call (format cuFUNCTION(..) )
   //             however, for the CUDA runtime API functions (format cudaFUNCTION(..) )
   //             the initialization is implicit, thus cuInit() here would not be needed...
   CUresult status = cuInit(0);
   if (CUDA_SUCCESS != status) exit_on_error("CUDA driver API device initialization failed\n");

   // returns a handle to the first cuda compute device
   CUdevice dev;
   status = cuDeviceGet(&dev, 0);
   if (CUDA_SUCCESS != status) exit_on_error("CUDA device not found\n");

   // gets device properties
   int major,minor;
   status = cuDeviceComputeCapability(&major,&minor,dev);
   if (CUDA_SUCCESS != status) exit_on_error("CUDA device information not found\n");

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
  device_count = 0;
  cudaGetDeviceCount(&device_count);
  // Do not check if command failed with `exit_on_cuda_error` since it calls cudaDevice()/ThreadSynchronize():
  // If multiple MPI tasks access multiple GPUs per node, they will try to synchronize
  // GPU 0 and depending on the order of the calls, an error will be raised
  // when setting the device number. If MPS is enabled, some GPUs will silently not be used.
  //
  // being verbose and catches error from first call to CUDA runtime function, without synchronize call
  cudaError_t err = cudaGetLastError();
  if (err != cudaSuccess){
    fprintf(stderr,"Error after cudaGetDeviceCount: %s\n", cudaGetErrorString(err));
    exit_on_error("CUDA runtime error: cudaGetDeviceCount failed\n\nplease check if driver and runtime libraries work together\nor on titan enable environment: CRAY_CUDA_PROXY=1 to use single GPU with multiple MPI processes\n\nexiting...\n");
  }

  // returns device count to fortran
  if (device_count == 0) exit_on_error("CUDA runtime error: there is no device supporting CUDA\n");
  *ncuda_devices = device_count;

  // Sets the active device
  if (device_count >= 1) {
    // generalized for more GPUs per node
    // note: without previous context release, cudaSetDevice will complain with the cuda error
    //         "setting the device when a process is active is not allowed"

    // releases previous contexts
#if CUDA_VERSION < 4000
    cudaThreadExit();
#else
    cudaDeviceReset();
#endif

    //printf("rank %d: cuda device count = %d sets device = %d \n",myrank,device_count,myrank % device_count);
    //MPI_Barrier(MPI_COMM_WORLD);

    // sets active device
#ifdef CUDA_DEVICE_ID
    // uses fixed device id when compile with e.g.: -DCUDA_DEVICE_ID=1
    device = CUDA_DEVICE_ID;
    if (myrank == 0) printf("setting cuda devices with id = %d for all processes by -DCUDA_DEVICE_ID\n\n",device);

    cudaSetDevice( device );
    exit_on_cuda_error("cudaSetDevice has invalid device");

    // double check that device was  properly selected
    cudaGetDevice(&device);
    if (device != CUDA_DEVICE_ID ){
       printf("error rank: %d devices: %d \n",myrank,device_count);
       printf("  cudaSetDevice()=%d\n  cudaGetDevice()=%d\n",CUDA_DEVICE_ID,device);
       exit_on_error("CUDA set/get device error: device id conflict \n");
    }
#else
    // device changes for different mpi processes according to number of device per node
    // (assumes that number of devices per node is the same for different compute nodes)
    device = myrank % device_count;

    cudaSetDevice( device );
    exit_on_cuda_error("cudaSetDevice has invalid device");

    // double check that device was  properly selected
    cudaGetDevice(&device);
    if (device != (myrank % device_count) ){
       printf("error rank: %d devices: %d \n",myrank,device_count);
       printf("  cudaSetDevice()=%d\n  cudaGetDevice()=%d\n",myrank%device_count,device);
       exit_on_error("CUDA set/get device error: device id conflict \n");
    }
#endif
  }

  // returns a handle to the active device
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
  int do_output_info;

  // by default, only master process outputs device infos to avoid file cluttering
  do_output_info = 0;
  if (myrank == 0){
    do_output_info = 1;
    sprintf(filename,OUTPUT_FILES"/gpu_device_info.txt");
  }
  // debugging
  if (DEBUG){
    do_output_info = 1;
    sprintf(filename,OUTPUT_FILES"/gpu_device_info_proc_%06d.txt",myrank);
  }

  // output to file
  if (do_output_info ){
    fp = fopen(filename,"w");
    if (fp != NULL){
      // display device properties
      fprintf(fp,"Device Name = %s\n",deviceProp.name);
      fprintf(fp,"memory:\n");
      fprintf(fp,"  totalGlobalMem (in MB): %f\n",(unsigned long) deviceProp.totalGlobalMem / (1024.f * 1024.f));
      fprintf(fp,"  totalGlobalMem (in GB): %f\n",(unsigned long) deviceProp.totalGlobalMem / (1024.f * 1024.f * 1024.f));
      fprintf(fp,"  totalConstMem (in bytes): %lu\n",(unsigned long) deviceProp.totalConstMem);
      fprintf(fp,"  Maximum 1D texture size (in bytes): %lu\n",(unsigned long) deviceProp.maxTexture1D);
      fprintf(fp,"  sharedMemPerBlock (in bytes): %lu\n",(unsigned long) deviceProp.sharedMemPerBlock);
      fprintf(fp,"  regsPerBlock (in bytes): %lu\n",(unsigned long) deviceProp.regsPerBlock);
      fprintf(fp,"blocks:\n");
      fprintf(fp,"  Maximum number of threads per block: %d\n",deviceProp.maxThreadsPerBlock);
      fprintf(fp,"  Maximum size of each dimension of a block: %d x %d x %d\n",
              deviceProp.maxThreadsDim[0],deviceProp.maxThreadsDim[1],deviceProp.maxThreadsDim[2]);
      fprintf(fp,"  Maximum sizes of each dimension of a grid: %d x %d x %d\n",
              deviceProp.maxGridSize[0],deviceProp.maxGridSize[1],deviceProp.maxGridSize[2]);
      fprintf(fp,"features:\n");
      fprintf(fp,"  Compute capability of the device = %d.%d\n", deviceProp.major, deviceProp.minor);
      fprintf(fp,"  multiProcessorCount: %d\n",deviceProp.multiProcessorCount);
      if (deviceProp.canMapHostMemory){
        fprintf(fp,"  canMapHostMemory: TRUE\n");
      }else{
        fprintf(fp,"  canMapHostMemory: FALSE\n");
      }
      if (deviceProp.deviceOverlap){
        fprintf(fp,"  deviceOverlap: TRUE\n");
      }else{
        fprintf(fp,"  deviceOverlap: FALSE\n");
      }
      if (deviceProp.concurrentKernels){
        fprintf(fp,"  concurrentKernels: TRUE\n");
      }else{
        fprintf(fp,"  concurrentKernels: FALSE\n");
      }
      // outputs initial memory infos via cudaMemGetInfo()
      double free_db,used_db,total_db;
      get_free_memory(&free_db,&used_db,&total_db);
      fprintf(fp,"memory usage:\n");
      fprintf(fp,"  rank %d: GPU memory usage: used = %f MB, free = %f MB, total = %f MB\n",myrank,
              used_db/1024.0/1024.0, free_db/1024.0/1024.0, total_db/1024.0/1024.0);

      // closes output file
      fclose(fp);
    }
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
  if (! deviceProp.canMapHostMemory){
    fprintf(stderr,"Device capability should allow to map host memory, exiting...\n");
    exit_on_error("CUDA Device capability canMapHostMemory should be TRUE\n");
  }

  // checks kernel optimization setting
#ifdef USE_LAUNCH_BOUNDS
  // see: mesh_constants_cuda.h
  // performance statistics: main kernel Kernel_2_**_impl():
  //       shared memory per block = 6200    for Kepler: total = 49152 -> limits active blocks to 7
  //       registers per thread    = 72                                   (limited by LAUNCH_MIN_BLOCKS 7)
  //       registers per block     = 9216                total = 65536    (limited by LAUNCH_MIN_BLOCKS 7)

  // shared memory
  if (deviceProp.sharedMemPerBlock > 49152 && LAUNCH_MIN_BLOCKS <= 7){
    if (myrank == 0){
      printf("GPU non-optimal settings: your setting of using LAUNCH_MIN_BLOCK %i is too low and limits the register usage\n",
             LAUNCH_MIN_BLOCKS);
    }
  }

  // registers
  if (deviceProp.regsPerBlock > 65536 && LAUNCH_MIN_BLOCKS <= 7){
    if (myrank == 0){
      printf("GPU non-optimal settings: your setting of using LAUNCH_MIN_BLOCK %i is too low and limits the register usage\n",
             LAUNCH_MIN_BLOCKS);
    }
  }
#endif

}
