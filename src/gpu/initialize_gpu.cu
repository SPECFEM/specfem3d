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

// CUDA-aware support
#ifdef WITH_CUDA_AWARE_MPI
#if defined(OPEN_MPI)
#include <mpi-ext.h> /* extensions */
#endif
#endif

// number of GPU cards (per compute node)
static int number_of_gpu_devices = 0;

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


void initialize_cuda_device(int myrank,int* nb_devices) {

  TRACE("initialize_cuda_device");

  int device;
  int device_count = 0;

  // user error info
  const char* err_info = "Please check GPU settings on your node \n\n";

  // sets gpu runtime flag
  run_cuda = 1;

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

  // stores counts
  number_of_gpu_devices = device_count;

  // returns device count to fortran
  *nb_devices = number_of_gpu_devices;

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

    err = cudaGetLastError();
    if (err != cudaSuccess) {
      fprintf(stderr,"Error cudaGetDevice: %s\n", cudaGetErrorString(err));
      if (err == cudaErrorDevicesUnavailable){ fprintf(stderr,"\n%s\n", err_info); }
      exit_on_error("CUDA runtime error: cudaGetDevice failed\n\n");
    }

    // checks device id
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

    err = cudaGetLastError();
    if (err != cudaSuccess) {
      fprintf(stderr,"Error cudaGetDevice: %s\n", cudaGetErrorString(err));
      if (err == cudaErrorDevicesUnavailable){ fprintf(stderr,"\n%s\n", err_info); }
      exit_on_error("CUDA runtime error: cudaGetDevice failed\n\n");
    }

    // checks device id
    if (device != (myrank % device_count) ){
       printf("Error rank: %d devices: %d \n",myrank,device_count);
       printf("  cudaSetDevice()=%d\n  cudaGetDevice()=%d\n",myrank%device_count,device);
       exit_on_error("CUDA set/get device error: device id conflict \n");
    }
#endif
  }
}

// outputs devices infos

static void output_cuda_device_infos(int myrank){

  struct cudaDeviceProp deviceProp;
  int device;

  // returns a handle to the active device
  cudaGetDevice(&device);
  exit_on_gpu_error("cudaGetDevice failed");

  // Gets number of GPU devices & version infos
  int device_count = 0;
  int driverVersion = 0, runtimeVersion = 0;

  cudaGetDeviceCount(&device_count);
  cudaDriverGetVersion(&driverVersion);
  cudaRuntimeGetVersion(&runtimeVersion);

  // get device properties
  cudaGetDeviceProperties(&deviceProp,device);
  exit_on_gpu_error("cudaGetDeviceProperties failed");

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
  int do_output_info = 0;

  // by default, only main process outputs device infos to avoid file cluttering
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
#if CUDA_VERSION < 13000 || (defined (__CUDACC_VER_MAJOR__) && (__CUDACC_VER_MAJOR__ < 13))
      if (deviceProp.deviceOverlap){
        fprintf(fp,"  deviceOverlap: TRUE\n");
      }else{
        fprintf(fp,"  deviceOverlap: FALSE\n");
      }
#else
      // CUDA version >= 13, deviceOverlap deprecated, replaced by asyncEngineCount
      fprintf(fp,"  asyncEngineCount: %d\n", deviceProp.asyncEngineCount);
#endif
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

void initialize_hip_device(int myrank,int* nb_devices) {

  TRACE("initialize_hip_device");

  int device;
  int device_count = 0;

  // first HIP call
  //
  // explicit initialization
  // (not necessary, most HIP APIs implicitly initialize the HIP runtime)
  //hipError_t status = hipInit(0);
  //if (status != hipSuccess) exit_on_error("HIP initialization failed\n");
  //
  // gets number of devices
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

  // stores counts
  number_of_gpu_devices = device_count;

  // returns device count to fortran
  *nb_devices = number_of_gpu_devices;

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
}

// outputs devices infos

static void output_hip_device_infos(int myrank){

  struct hipDeviceProp_t deviceProp;
  int device;

  // returns a handle to the active device
  hipGetDevice(&device);
  exit_on_gpu_error("hipGetDevice failed");

  // Gets number of GPU devices & version infos
  int device_count = 0;
  int driverVersion = 0, runtimeVersion = 0;

  hipGetDeviceCount(&device_count);
  hipDriverGetVersion(&driverVersion);
  hipRuntimeGetVersion(&runtimeVersion);

  // get device properties
  hipGetDeviceProperties(&deviceProp,device);
  exit_on_gpu_error("hipGetDevicePropoerties failed");

  // memory usage
  double free_db,used_db,total_db;
  get_free_memory(&free_db,&used_db,&total_db);

  // outputs device infos to file
  char filename[BUFSIZ];
  FILE* fp;
  int do_output_info = 0;

  // by default, only master process outputs device infos to avoid file cluttering
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
              INITIALIZE_GPU_DEVICE)(int* myrank_f,int* nb_devices,int* cuda_aware_mpi_init_type) {

  TRACE("initialize_gpu_device");

  // rank
  int myrank = *myrank_f;
  int init_type = *cuda_aware_mpi_init_type;

  // flags to run initialization and output device infos
  int do_init = 1;
  int do_output = 1;

  // check if compiled with both CUDA and HIP support
#if defined(USE_CUDA) && defined(USE_HIP)
  if (myrank == 0) {
    printf("Error: GPU version compilation with both USE_CUDA and USE_HIP not supported yet.\nPlease only use one for now...\n\n",);
  }
  exit(1);
#endif

  // CUDA-aware MPI
  // we need to set the GPU device before MPI_init but should avoid calling further CUDA calls to avoid issues.
  // for example, on Summit the PAMI backend uses "CUDA hooks" and would complain about:
  //    CUDA Hook Library: Failed to find symbol mem_find_dreg_entries, ./bin/xspecfem3D: undefined symbol: __PAMI_Invalidate_region
  // see: https://docs.olcf.ornl.gov/systems/summit_user_guide.html#cuda-hook-error-when-program-uses-cuda-without-first-calling-mpi-init
  //
  // thus, we separate the initialization and the device output (which contains a memory allocation check leading to this problem).
#ifdef WITH_CUDA_AWARE_MPI
  // checks if initialize called by CUDA-aware check
  if (init_type == 1) {
    // initial call to set device
    if (myrank == 0){ printf("using CUDA-aware MPI: initializing GPU devices\n"); }
    // only initialization
    do_init = 1;
    do_output = 0;
  }else if (init_type == 2){
    // called again with Par_file setting
    if (myrank == 0){ printf("using CUDA-aware MPI: returning number of devices = %d\n",number_of_gpu_devices); }
    // already initialized
    *nb_devices = number_of_gpu_devices;
    // only device infos
    do_init = 0;
    do_output = 1;
  }
#endif // WITH_CUDA_AWARE_MPI

  // initializes gpu cards
  if (do_init) {
#ifdef USE_CUDA
    run_cuda = 1;
    if (run_cuda) { initialize_cuda_device(myrank, nb_devices); }
#endif
#ifdef USE_HIP
    run_hip = 1;
    if (run_hip) { initialize_hip_device(myrank, nb_devices); }
#endif
  }

  // outputs device infos
  if (do_output){
#ifdef USE_CUDA
    if (run_cuda) { output_cuda_device_infos(myrank); }
#endif
#ifdef USE_HIP
    if (run_hip) { output_hip_device_infos(myrank); }
#endif
  }
}

/* ----------------------------------------------------------------------------------------------- */

// CUDA-aware MPI
// we need to call cudaSetDevice before MPI_Init to ensure that the same GPU is chosen by MPI and your application

extern EXTERN_LANG
void FC_FUNC_ (check_cuda_aware_mpi,
               CHECK_CUDA_AWARE_MPI) (int* has_cuda_aware_mpi_f) {

  TRACE ("check_cuda_aware_mpi");

  // flags
  int has_cuda_aware_mpi = 0;

#ifdef WITH_CUDA_AWARE_MPI
  // environment variable which allows the reading of the local rank of the current MPI
  // process before the MPI environment gets initialized with MPI_Init().
  //
  // This is necessary when running the CUDA-aware MPI version, which needs this information in order to be able to
  // set the CUDA device for the MPI process before MPI environment initialization.
  //
  // If you are using MVAPICH2, set this constant to "MV2_COMM_WORLD_LOCAL_RANK";
  // for Open MPI, use "OMPI_COMM_WORLD_LOCAL_RANK".
#if defined(OPEN_MPI) && OPEN_MPI
// OpenMPI
#pragma message ("\n\nCompiling with: WITH_CUDA_AWARE_MPI uses OPEN_MPI CUDA-aware local rank\n")
#define ENV_LOCAL_RANK    "OMPI_COMM_WORLD_LOCAL_RANK"

#elif defined(MVAPICH2_NUMVERSION) && (MVAPICH2_NUMVERSION >= 20205300)
// MVAPICH
#pragma message ("\n\nCompiling with: WITH_CUDA_AWARE_MPI uses MVAPICH2 CUDA-aware local rank\n")
#define ENV_LOCAL_RANK    "MV2_COMM_WORLD_LOCAL_RANK"

#else
// unknown
#pragma message ("\n\nCompiling with: unknown CUDA-aware local rank environment, use -DENV_LOCAL_RANK \"<MY_LOCAL_RANK>\" setting\n")
// defines local rank environment variables as unknown if not set by compilation flag, mostly to be able to run getenv() command
#ifndef ENV_LOCAL_RANK
#define ENV_LOCAL_RANK    "UNKNOWN_LOCAL_RANK"
#endif

#endif

  // sets GPU device before MPI initialization
  // MPI will then recognize the setting and take over the GPU device setup

  // determine local rank
  // note: local rank is the rank id per compute node
  //       for example, 4 MPI processes per node -> local rank id = 0,1,2,3 on all cmopute nodes
  //       not the same as the MPI rank which goes from 0 to MPI size-1
  int has_local_rank_info = 0;
  int rank = 0;
  char * localRankStr = NULL;

  // debug output to file
  char filename[BUFSIZ];
  FILE* fp;
  sprintf(filename,OUTPUT_FILES"/gpu_aware_info.txt");

  // local rank info from environment
  if ((localRankStr = getenv(ENV_LOCAL_RANK)) != NULL) {
    // catching OpenMPI environment rank
    rank = atoi(localRankStr);
    has_local_rank_info = 1;
  } else {
    // no OpenMPI environment rank found, initializing myrank to zero
    rank = 0;
    has_local_rank_info = 0;
  }

  // debug
  //printf("debug: [check_cuda_aware_mpi] CUDA-aware check: has_local_rank_info = %d  -  local rank = %d\n",has_local_rank_info,rank);

  // enables CUDA-aware MPI support
  if (has_local_rank_info){
    // user output
    if (rank == 0){ printf("\nchecking CUDA-aware MPI\n\n");}

    // debug
    //printf("debug: compile time check for CUDA-aware MPI - rank %d\n",rank);

#if defined(MPIX_CUDA_AWARE_SUPPORT)
#pragma message ("\n\nCompiling with: WITH_CUDA_AWARE_MPI has MPIX_CUDA_AWARE_SUPPORT\n")
    int ret = MPIX_Query_cuda_support();
    if (ret == 1) {
      // MPI library has CUDA-aware support
      has_cuda_aware_mpi = 1;
    } else {
      // MPI library does not have CUDA-aware support
      has_cuda_aware_mpi = 0;
    }
#else
#pragma message ("\n\nCompiling with: WITH_CUDA_AWARE_MPI has no MPIX_CUDA_AWARE_SUPPORT, please check MPI installation\n")
    // user info
    if (rank == 0){
      printf("\
This version has been compiled with flag WITH_CUDA_AWARE_MPI, but MPI library cannot determine if there is CUDA-aware support.\n \
Please check MPI installation.\n\n");
      // file output
      fp = fopen(filename,"w");
      if (fp != NULL){
        fprintf (fp, "\
This version has been compiled with flag WITH_CUDA_AWARE_MPI, but MPI library cannot determine if there is CUDA-aware support.\n \
Please check MPI installation.\n\n");
        fclose(fp);
      }
    }
    has_cuda_aware_mpi = 0;
#endif  // MPIX_CUDA_AWARE_SUPPORT

    // debug
    //printf("debug: query cuda support: MPI library CUDA-aware support - rank %d has support %d\n\n",rank,has_cuda_aware_mpi);

    // sets local rank's GPU association
    if (has_cuda_aware_mpi){
      // dummy value, not needed at this point
      int dummy_nb_devices;
      int init_type = 1; // type 1 == only initialize, no device info output yet
      // debug
      //printf("debug: setting - rank %d has support %d - running GPU init\n\n",rank,has_cuda_aware_mpi);

      // user info
      if (rank == 0) {
        fp = fopen(filename,"w");
        if (fp != NULL){
          fprintf (fp, "gpu: has CUDA-aware MPI.\n\n");
          fclose(fp);
        }
      }

      // sets device
      FC_FUNC_(initialize_gpu_device,INITIALIZE_GPU_DEVICE)(&rank,&dummy_nb_devices,&init_type);
    }
  }

#endif // WITH_CUDA_AWARE_MPI

  // return value
  *has_cuda_aware_mpi_f = has_cuda_aware_mpi;
}
