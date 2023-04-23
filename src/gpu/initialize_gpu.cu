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

#include "mesh_constants_gpu.h"

// gpu runtime flags
int run_cuda = 0;
int run_hip = 0;

/* ----------------------------------------------------------------------------------------------- */
// CUDA initialization
/* ----------------------------------------------------------------------------------------------- */

// CUDA version output
#ifdef USE_CUDA

// macros for version output
#define VALUE_TO_STRING(x) #x
#define VALUE(x) VALUE_TO_STRING(x)
#define VAR_NAME_VALUE(var) #var " = "  VALUE(var)

#pragma message ("\n\nCompiling with: " VAR_NAME_VALUE(CUDA_VERSION) "\n")
#if defined(__CUDA_ARCH__)
#pragma message ("\n\nCompiling with: " VAR_NAME_VALUE(__CUDA_ARCH__) "\n")
#endif

// CUDA version >= 4.0 needed for cudaTextureType1D and cudaDeviceSynchronize()
#if CUDA_VERSION < 4000 || (defined (__CUDACC_VER_MAJOR__) && (__CUDACC_VER_MAJOR__ < 4))
#pragma message ("\n\nCompiling for CUDA version < 4.0\n")
#endif


void initialize_cuda_device(int* myrank_f,int* ncuda_devices) {

  TRACE("initialize_cuda_device");

  int device;
  int device_count;

  // sets gpu runtime flag
  run_cuda = 1;

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
  // Do not check if command failed with `exit_on_gpu_error` since it calls cudaDevice()/ThreadSynchronize():
  // If multiple MPI tasks access multiple GPUs per node, they will try to synchronize
  // GPU 0 and depending on the order of the calls, an error will be raised
  // when setting the device number. If MPS is enabled, some GPUs will silently not be used.
  //
  // being verbose and catches error from first call to CUDA runtime function, without synchronize call
  cudaError_t err = cudaGetLastError();

  // adds quick check on versions
  int driverVersion = 0, runtimeVersion = 0;
  cudaDriverGetVersion(&driverVersion);
  cudaRuntimeGetVersion(&runtimeVersion);

  // exit in case first cuda call failed
  if (err != cudaSuccess){
    fprintf(stderr,"Error after cudaGetDeviceCount: %s\n", cudaGetErrorString(err));
    fprintf(stderr,"CUDA Device count: %d\n",device_count);
    fprintf(stderr,"CUDA Driver Version / Runtime Version: %d.%d / %d.%d\n",
            driverVersion / 1000, (driverVersion % 100) / 10,
            runtimeVersion / 1000, (runtimeVersion % 100) / 10);

    exit_on_error("CUDA runtime error: cudaGetDeviceCount failed\n\nplease check if driver and runtime libraries work together\nor on cluster environments enable MPS (Multi-Process Service) to use single GPU with multiple MPI processes\n\nexiting...\n");
  }

  // checks if CUDA devices available
  if (device_count == 0) exit_on_error("CUDA runtime error: there is no device supporting CUDA\n");

  // returns device count to fortran
  *ncuda_devices = device_count;

  // Sets the active device
  if (device_count >= 1) {
    // generalized for more GPUs per node
    // note: without previous context release, cudaSetDevice will complain with the cuda error
    //         "setting the device when a process is active is not allowed"

    // releases previous contexts
    gpuReset();

    //printf("rank %d: cuda device count = %d sets device = %d \n",myrank,device_count,myrank % device_count);
    //MPI_Barrier(MPI_COMM_WORLD);

    // sets active device
#ifdef GPU_DEVICE_ID
    // uses fixed device id when compile with e.g.: -DGPU_DEVICE_ID=1
    device = GPU_DEVICE_ID;
    if (myrank == 0) printf("setting CUDA devices with id = %d for all processes by -DGPU_DEVICE_ID\n\n",device);

    cudaSetDevice( device );
    exit_on_gpu_error("cudaSetDevice has invalid device");

    // double check that device was  properly selected
    cudaGetDevice(&device);
    if (device != GPU_DEVICE_ID ){
       printf("Error rank: %d devices: %d \n",myrank,device_count);
       printf("  cudaSetDevice()=%d\n  cudaGetDevice()=%d\n",GPU_DEVICE_ID,device);
       exit_on_error("CUDA set/get device error: device id conflict \n");
    }
#else
    // device changes for different mpi processes according to number of device per node
    // (assumes that number of devices per node is the same for different compute nodes)
    device = myrank % device_count;

    cudaSetDevice( device );
    exit_on_gpu_error("cudaSetDevice has invalid device");

    // double check that device was  properly selected
    cudaGetDevice(&device);
    if (device != (myrank % device_count) ){
       printf("Error rank: %d devices: %d \n",myrank,device_count);
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

  // memory usage
  double free_db,used_db,total_db;
  get_free_memory(&free_db,&used_db,&total_db);

  // outputs device infos to file
  char filename[BUFSIZ];
  FILE* fp;
  int do_output_info;

  // by default, only main process outputs device infos to avoid file cluttering
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
      fprintf(fp,"CUDA Device count: %d\n",device_count);
      fprintf(fp,"CUDA Driver Version / Runtime Version          %d.%d / %d.%d\n",
              driverVersion / 1000, (driverVersion % 100) / 10,
              runtimeVersion / 1000, (runtimeVersion % 100) / 10);

      // outputs initial memory infos via cudaMemGetInfo()
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
  // see: mesh_constants_gpu.h
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
#endif // USE_CUDA



/* ----------------------------------------------------------------------------------------------- */
// HIP initialization
/* ----------------------------------------------------------------------------------------------- */

#ifdef USE_HIP

void initialize_hip_device(int* myrank_f,int* ncuda_devices) {

  TRACE("initialize_hip_device");

  int device;
  int device_count;

  // Gets rank number of MPI process
  int myrank = *myrank_f;

  // first HIP call
  //
  // explicit initialization
  // (not necessary, most HIP APIs implicitly initialize the HIP runtime)
  //hipError_t status = hipInit(0);
  //if (status != hipSuccess) exit_on_error("HIP initialization failed\n");
  //
  // gets number of devices
  device_count = 0;
  hipGetDeviceCount(&device_count);
  hipError_t err = hipGetLastError();

  // adds quick check on versions
  int driverVersion = 0, runtimeVersion = 0;
  hipDriverGetVersion(&driverVersion);
  hipRuntimeGetVersion(&runtimeVersion);

  // exit in case first HIP call failed
  if (err != hipSuccess){
    fprintf (stderr,"Error after hipGetDeviceCount: %s\n", hipGetErrorString(err));
    fprintf (stderr,"HIP Device count: %d\n",device_count);
    fprintf (stderr,"HIP Driver Version / Runtime Version: %d.%d / %d.%d\n",
                    driverVersion / 1000, (driverVersion % 100) / 10,
                    runtimeVersion / 1000, (runtimeVersion % 100) / 10);

    exit_on_error("HIP runtime error: hipGetDeviceCount failed\n\nPlease check if any HIP devices are available\n\nexiting...\n");
  }

  // checks if HIP devices available
  if (device_count == 0) exit_on_error("HIP runtime error: no HIP devices available\n");

  // returns device count to fortran
  *ncuda_devices = device_count;

  // Sets the active device
  if (device_count >= 1) {
    // generalized for more GPUs per node
    // note: without previous context release, hipSetDevice will complain with the cuda error
    //         "setting the device when a process is active is not allowed"

    // releases previous contexts
#if CUDA_VERSION < 4000
    hipDeviceReset();
#else
    hipDeviceReset();
#endif

    //printf("rank %d: cuda device count = %d sets device = %d \n",myrank,device_count,myrank % device_count);
    //MPI_Barrier(MPI_COMM_WORLD);

    // sets active device
#ifdef GPU_DEVICE_ID
    // uses fixed device id when compile with e.g.: -DGPU_DEVICE_ID=1
    device = GPU_DEVICE_ID;
    if (myrank == 0) printf("setting HIP devices with id = %d for all processes by -DGPU_DEVICE_ID\n\n",device);

    hipSetDevice( device );
    exit_on_gpu_error("hipSetDevice has invalid device");

    // double check that device was  properly selected
    hipGetDevice(&device);
    if (device != GPU_DEVICE_ID ){
       printf("Error rank: %d devices: %d \n",myrank,device_count);
       printf("  hipSetDevice()=%d\n  hipGetDevice()=%d\n",GPU_DEVICE_ID,device);
       exit_on_error("HIP set/get device error: device id conflict \n");
    }
#else
    // device changes for different mpi processes according to number of device per node
    // (assumes that number of devices per node is the same for different compute nodes)
    device = myrank % device_count;

    hipSetDevice( device );
    exit_on_gpu_error("hipSetDevice has invalid device");

    // double check that device was  properly selected
    hipGetDevice(&device);
    if (device != (myrank % device_count) ){
       printf("Error rank: %d devices: %d \n",myrank,device_count);
       printf("  hipSetDevice()=%d\n  hipGetDevice()=%d\n",myrank%device_count,device);
       exit_on_error("HIP set/get device error: device id conflict \n");
    }
#endif
  }

  // returns a handle to the active device
  hipGetDevice(&device);
  exit_on_gpu_error("hipGetDevice failed");

  // get device properties
  struct hipDeviceProp_t deviceProp;
  hipGetDeviceProperties(&deviceProp,device);
  exit_on_gpu_error("hipGetDevicePropoerties failed");

  // memory usage
  double free_db,used_db,total_db;
  get_free_memory(&free_db,&used_db,&total_db);

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
      fprintf (fp, "Device Name = %s\n", deviceProp.name);
      fprintf (fp, "memory:\n");
      fprintf (fp, "  totalGlobalMem (in MB): %f\n",(unsigned long) deviceProp.totalGlobalMem / (1024.f * 1024.f));
      fprintf (fp, "  totalGlobalMem (in GB): %f\n",(unsigned long) deviceProp.totalGlobalMem / (1024.f * 1024.f * 1024.f));
      fprintf (fp, "  totalConstMem (in bytes): %lu\n",(unsigned long) deviceProp.totalConstMem); // seems to be same as GlobalMem
      //fprintf (fp, "  Maximum 1D texture size (in bytes): %lu\n",(unsigned long) deviceProp.maxTexture1D); // not available?
      fprintf (fp, "  sharedMemPerBlock (in bytes): %lu\n",(unsigned long) deviceProp.sharedMemPerBlock);
      fprintf (fp, "  regsPerBlock (in bytes): %lu\n",(unsigned long) deviceProp.regsPerBlock);
      fprintf (fp, "blocks:\n");
      fprintf (fp, "  Maximum number of threads per block: %d\n",deviceProp.maxThreadsPerBlock);
      fprintf (fp, "  Maximum size of each dimension of a block: %d x %d x %d\n",
                       deviceProp.maxThreadsDim[0],deviceProp.maxThreadsDim[1],deviceProp.maxThreadsDim[2]);
      fprintf (fp, "  Maximum sizes of each dimension of a grid: %d x %d x %d\n",
                       deviceProp.maxGridSize[0],deviceProp.maxGridSize[1],deviceProp.maxGridSize[2]);
      fprintf (fp, "features:\n");
      fprintf (fp, "  Compute capability of the device = %d.%d\n", deviceProp.major, deviceProp.minor);
      fprintf (fp, "  multiProcessorCount: %d\n",deviceProp.multiProcessorCount);
      if (deviceProp.canMapHostMemory){
        fprintf (fp, "  canMapHostMemory: TRUE\n");
      }else{
        fprintf (fp, "  canMapHostMemory: FALSE\n");
      }
      if (deviceProp.concurrentKernels){
        fprintf (fp, "  concurrentKernels: TRUE\n");
      }else{
        fprintf (fp, "  concurrentKernels: FALSE\n");
      }

      fprintf(fp,"HIP Device count: %d\n",device_count);
      fprintf(fp,"HIP Driver Version / Runtime Version          %d.%d / %d.%d\n",
              driverVersion / 1000, (driverVersion % 100) / 10,
              runtimeVersion / 1000, (runtimeVersion % 100) / 10);

      // outputs initial memory infos via hipMemGetInfo()
      fprintf(fp,"memory usage:\n");
      fprintf(fp,"  rank %d: GPU memory usage: used = %f MB, free = %f MB, total = %f MB\n",myrank,
              used_db/1024.0/1024.0, free_db/1024.0/1024.0, total_db/1024.0/1024.0);

      // closes output file
      fclose(fp);
    }
  }

  /* daniel todo: check in case this applies...
  // we use pinned memory for asynchronous copy
  if (! deviceProp.canMapHostMemory){
    fprintf(stderr,"Device capability should allow to map host memory, exiting...\n");
    exit_on_error("CUDA Device capability canMapHostMemory should be TRUE\n");
  }
  */

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
#endif // USE_HIP

/* ----------------------------------------------------------------------------------------------- */

// GPU initialization

/* ----------------------------------------------------------------------------------------------- */

extern EXTERN_LANG
void FC_FUNC_(initialize_gpu_device,
              INITIALIZE_GPU_DEVICE)(int* myrank_f,int* ncuda_devices) {

  TRACE("initialize_gpu_device");

  // check if compiled with both CUDA and HIP support
#if defined(USE_CUDA) && defined(USE_HIP)
  if (*myrank_f == 0) {
    printf("Error: GPU version compilation with both USE_CUDA and USE_HIP not supported yet.\nPlease only use one for now...\n\n",);
  }
  exit(1);
#endif

  // initializes gpu cards
#ifdef USE_CUDA
  run_cuda = 1;
  if (run_cuda) {
    initialize_cuda_device(myrank_f, ncuda_devices);
  }
#endif
#ifdef USE_HIP
  run_hip = 1;
  if (run_hip) {
    initialize_hip_device(myrank_f, ncuda_devices);
  }
#endif
}
