/*
**************************

check_cuda_device utility

**************************

this utility program will output GPU device information


for compilation, see the command-line examples given here:

- example without MPI support:

nvcc --gpu-architecture=sm_60 -o check_cuda_device check_cuda_device.cu
./check_cuda_device

- example with MPI support:

nvcc -arch=sm_60 -DWITH_MPI -I/usr/lib/openmpi/include -o check_cuda_device check_cuda_device.cu -lmpi -L/usr/lib/openmpi/lib
mpirun -np 2 ./check_cuda_device


*/

#include <stdio.h>
#include <cuda.h>
#include <unistd.h>

#ifdef WITH_MPI
#include <mpi.h>
#endif

#include <sys/time.h>
#include <sys/resource.h>

/* ----------------------------------------------------------------------------------------------- */

void get_free_memory(double* free_db, double* used_db, double* total_db) {

  // gets memory usage in byte
  size_t free_byte ;
  size_t total_byte ;
  cudaError_t cuda_status = cudaMemGetInfo( &free_byte, &total_byte ) ;
  if ( cudaSuccess != cuda_status ){
    printf("Error: cudaMemGetInfo fails, %s \n", cudaGetErrorString(cuda_status) );
    exit(EXIT_FAILURE);
  }

  *free_db = (double)free_byte ;
  *total_db = (double)total_byte ;
  *used_db = *total_db - *free_db ;
  return;
}

/* ----------------------------------------------------------------------------------------------- */

void exit_on_error(const char* info) {
  printf("\nERROR: %s\n",info);
  fflush(stdout);

  // stops program
#ifdef WITH_MPI
  MPI_Abort(MPI_COMM_WORLD,1);
#endif
  //free(info);
  exit(EXIT_FAILURE);
  return;
}

/* ----------------------------------------------------------------------------------------------- */

void exit_on_gpu_error(const char* kernel_name) {
  // sync and check to catch errors from previous async operations
#if CUDA_VERSION < 4000 || (defined (__CUDACC_VER_MAJOR__) && (__CUDACC_VER_MAJOR__ < 4))
  cudaThreadSynchronize();
#else
  cudaDeviceSynchronize();
#endif

  cudaError_t err = cudaGetLastError();
  if (err != cudaSuccess){
    printf("Error after %s: %s\n", kernel_name, cudaGetErrorString(err));

    // releases previous contexts
#if CUDA_VERSION < 4000 || (defined (__CUDACC_VER_MAJOR__) && (__CUDACC_VER_MAJOR__ < 4))
    cudaThreadExit();
#else
    cudaDeviceReset();
#endif

    // stops program
    //free(kernel_name);
#ifdef WITH_MPI
    MPI_Abort(MPI_COMM_WORLD,1);
#endif
    exit(EXIT_FAILURE);
  }
}



/* ----------------------------------------------------------------------------------------------- */

// GPU initialization

/* ----------------------------------------------------------------------------------------------- */


void initialize_cuda_device(int* myrank_f,int* ncuda_devices) {

  int device;
  int device_count;
  cudaError_t err;

  // Gets rank number of MPI process
  int myrank = *myrank_f;

#ifdef WITH_MPI
  int sizeprocs;
  MPI_Comm_size(MPI_COMM_WORLD,&sizeprocs);
  MPI_Barrier(MPI_COMM_WORLD);
#endif

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
  device_count = 0;
  cudaGetDeviceCount(&device_count);
  // Do not check if command failed with `exit_on_gpu_error` since it calls cudaDevice()/ThreadSynchronize():
  // If multiple MPI tasks access multiple GPUs per node, they will try to synchronize
  // GPU 0 and depending on the order of the calls, an error will be raised
  // when setting the device number. If MPS is enabled, some GPUs will silently not be used.
  //
  // being verbose and catches error from first call to CUDA runtime function, without synchronize call
  err = cudaGetLastError();
  if (err != cudaSuccess) {
    fprintf(stderr,"Error after cudaGetDeviceCount: %s\n", cudaGetErrorString(err));
    exit_on_error("\
CUDA runtime error: cudaGetDeviceCount failed\n\n\
please check if driver and runtime libraries work together\n\n");
  }

  // returns device count to fortran
  if (device_count == 0) {
    exit_on_error("CUDA runtime error: there is no device supporting CUDA\n");
  }

  // output
#ifdef WITH_MPI
  // output infos for mpi ranks ordered
  for(int iproc = 0;iproc < sizeprocs; iproc++){
    if (iproc == myrank){
      printf("process %d found number of CUDA devices = %d\n",myrank,device_count);
      fflush(stdout);
    }
    // synchronizes mpi processes
    MPI_Barrier(MPI_COMM_WORLD);
  }
  sleep(1);
  if (myrank == 0){printf("\n\n");fflush(stdout);}
  // synchronizes mpi processes
  MPI_Barrier(MPI_COMM_WORLD);
#else
  printf("\nfound number of CUDA devices = %d\n",device_count);fflush(stdout);
#endif

  *ncuda_devices = device_count;

  // user error info
  const char* err_info = "\
Please check GPU settings on your node \n\
and/or check CUDA MPS setup to use a single GPU with multiple MPI processes,\n\
e.g., on titan enable environment CRAY_CUDA_MPS=1 to use a single GPU with multiple MPI processes\n\n";

  // Sets the active device
  // generalized for more GPUs per node
  // note: without previous context release, cudaSetDevice will complain with the cuda error
  //         "setting the device when a process is active is not allowed"

  // releases previous contexts
#if CUDA_VERSION < 4000 || (defined (__CUDACC_VER_MAJOR__) && (__CUDACC_VER_MAJOR__ < 4))
  cudaThreadExit();
#else
  cudaDeviceReset();
#endif

#ifdef WITH_MPI
  // check
  if (cudaGetLastError() != cudaSuccess){
    fprintf(stderr,"Error: %s\n", cudaGetErrorString(cudaGetLastError()));
    exit_on_error("CUDA runtime error: cudaDeviceReset failed\n\n");
  }
  // synchronizes mpi processes
  MPI_Barrier(MPI_COMM_WORLD);
  // output infos for mpi ranks ordered
  for(int iproc = 0;iproc < sizeprocs; iproc++){
    if (iproc == myrank){
      printf("process %d: cuda device count = %d (would select device = %d)\n",myrank,device_count,myrank % device_count);
      fflush(stdout);
    }
    // synchronizes mpi processes
    MPI_Barrier(MPI_COMM_WORLD);
  }
  sleep(1);
  if (myrank == 0){printf("\n\n");fflush(stdout);}
  // synchronizes mpi processes
  MPI_Barrier(MPI_COMM_WORLD);
#endif

  // loops over all devices for displaying info
  for(int i=0;i < device_count; i++){

#ifdef WITH_MPI
    // output info
    MPI_Barrier(MPI_COMM_WORLD);
    if (myrank == 0){printf("cuda set device %d\n\n",i);fflush(stdout);}
    for(int iproc = 0;iproc < sizeprocs; iproc++){
      if (iproc == myrank){
        printf("process %d: cudaSetDevice %d\n",myrank,i);fflush(stdout);
#endif

    // sets active device
    //device = myrank % device_count;
    device = i;
    err = cudaSetDevice( device );
    if (err != cudaSuccess) {
      fprintf(stderr,"Error cudaSetDevice: %s\n", cudaGetErrorString(err));
      if (err == cudaErrorDevicesUnavailable){ fprintf(stderr,"\n%s\n", err_info); }
      exit_on_error("CUDA runtime error: cudaSetDevice failed\n\n");
    }

    // checks device execution
    err = cudaDeviceSynchronize();
    if (err != cudaSuccess) {
      fprintf(stderr,"Error cudaDeviceSynchronize: %s\n", cudaGetErrorString(err));
      exit_on_error("CUDA runtime error: cudaDeviceSynchronize failed\n\n");
    }

#ifdef WITH_MPI
      } // if
      MPI_Barrier(MPI_COMM_WORLD);
    }//for
    // double check if CUDA context gets created on multiple processes
    MPI_Barrier(MPI_COMM_WORLD);
    for(int iproc = 0;iproc < sizeprocs; iproc++){
      if (iproc == myrank){
        printf("process %d: create context\n",myrank);
        fflush(stdout);
        // creates context
        err = cudaFree(0);
        if (err != cudaSuccess) {
          printf("Error cudaFree: %s\n", cudaGetErrorString(err));
          exit_on_error("CUDA runtime error: cudaFree for context failed\n\n");
        }
      }
      // synchronizes mpi processes
      MPI_Barrier(MPI_COMM_WORLD);
    }
    sleep(1);
    if (myrank == 0){printf("\n\n");fflush(stdout);}
    MPI_Barrier(MPI_COMM_WORLD);
#endif

    // double check
    exit_on_gpu_error("cudaSetDevice has invalid device");

#ifdef WITH_MPI
    // check
    MPI_Barrier(MPI_COMM_WORLD);
    if (myrank == 0){printf("cuda get device %d\n\n",i);fflush(stdout);}
#endif

    // double check that device was  properly selected
    cudaGetDevice(&device);

    err = cudaGetLastError();
    if (err != cudaSuccess) {
      printf("Error cudaGetDevice: %s\n", cudaGetErrorString(err));
      if (err == cudaErrorDevicesUnavailable){ printf("\n%s\n", err_info); }
      exit_on_error("CUDA runtime error: cudaGetDevice failed\n\n");
    }

#ifdef WITH_MPI
    // output infos for mpi ranks ordered
    MPI_Barrier(MPI_COMM_WORLD);
    for(int iproc = 0;iproc < sizeprocs; iproc++){
      if (iproc == myrank){
        printf("device set/get: rank %d set %d get %d\n - return %s\n",myrank,i,device,cudaGetErrorString(err));
        fflush(stdout);
      }
      // synchronizes mpi processes
      MPI_Barrier(MPI_COMM_WORLD);
    }
    sleep(1);
    if (myrank == 0){printf("\n\n");fflush(stdout);}
    MPI_Barrier(MPI_COMM_WORLD);
#endif

    // checks device id
    //if( device != (myrank % device_count) ){
    if( device != i ){
       printf("error rank: %d devices: %d \n",myrank,device_count);
       printf("  cudaSetDevice()=%d\n  cudaGetDevice()=%d\n",myrank%device_count,device);
       exit_on_error("CUDA set/get device error: device id conflict \n");
    }

    // get device properties
    struct cudaDeviceProp deviceProp;
    cudaGetDeviceProperties(&deviceProp,device);

    exit_on_gpu_error("cudaGetDeviceProperties failed");

    // exit if the machine has no CUDA-enabled device
    if (deviceProp.major == 9999 && deviceProp.minor == 9999){
      printf("No CUDA-enabled device found, exiting...\n\n");
      exit_on_error("CUDA runtime error: there is no CUDA-enabled device found\n");
    }

    // memory infos via cudaMemGetInfo()
    double free_db,used_db,total_db;
    get_free_memory(&free_db,&used_db,&total_db);

    // ordering mpi output
#ifdef WITH_MPI
    // synchronizes mpi processes
    MPI_Barrier(MPI_COMM_WORLD);
    // output infos for mpi ranks ordered
    for(int iproc = 0;iproc < sizeprocs; iproc++){
      if (iproc == myrank){
        //printf("\n\nGPU device for rank: %d - total procs: %d\n\n",myrank,sizeprocs);
#endif

    // outputs device infos to file
    printf("\n\nGPU device id: %d\n\n",i);

    // display device properties
    printf("Device Name = %s\n\n",deviceProp.name);
    printf("memory:\n");
    printf("  totalGlobalMem (in MB, dividing by powers of 1024): %f\n",(unsigned long) deviceProp.totalGlobalMem / (1024.f * 1024.f));
    printf("  totalGlobalMem (in GB, dividing by powers of 1024): %f\n\n",(unsigned long) deviceProp.totalGlobalMem / (1024.f * 1024.f * 1024.f));
    printf("  totalGlobalMem (in MB, dividing by powers of 1000): %f\n",(unsigned long) deviceProp.totalGlobalMem / (1000.f * 1000.f));
    printf("  totalGlobalMem (in GB, dividing by powers of 1000): %f\n\n",(unsigned long) deviceProp.totalGlobalMem / (1000.f * 1000.f * 1000.f));
    printf("  sharedMemPerBlock (in bytes): %lu\n\n",(unsigned long) deviceProp.sharedMemPerBlock);
    printf("blocks:\n");
    printf("  Maximum number of registers per block: %d\n",deviceProp.regsPerBlock);
    printf("  Maximum number of threads per block: %d\n",deviceProp.maxThreadsPerBlock);
    printf("  Maximum size of each dimension of a block: %d x %d x %d\n",
            deviceProp.maxThreadsDim[0],deviceProp.maxThreadsDim[1],deviceProp.maxThreadsDim[2]);
    printf("  Maximum sizes of each dimension of a grid: %d x %d x %d\n\n",
            deviceProp.maxGridSize[0],deviceProp.maxGridSize[1],deviceProp.maxGridSize[2]);
    printf("features:\n");
    printf("  Compute capability of the device = %d.%d\n", deviceProp.major, deviceProp.minor);
    printf("  multiProcessorCount: %d\n",deviceProp.multiProcessorCount);
    if(deviceProp.canMapHostMemory){
      printf("  canMapHostMemory: TRUE\n");
    }else{
      printf("  canMapHostMemory: FALSE\n");
    }
    if(deviceProp.deviceOverlap){
      printf("  deviceOverlap: TRUE\n");
    }else{
      printf("  deviceOverlap: FALSE\n");
    }
    printf("  Compute Mode: %d\n", deviceProp.computeMode);
    fflush(stdout);


    // outputs initial memory infos via cudaMemGetInfo()
    printf("\n%d: GPU memory usage (dividing by powers of 1024): used = %f MB, free = %f MB, total = %f MB",myrank,
            used_db/1024.0/1024.0, free_db/1024.0/1024.0, total_db/1024.0/1024.0);
    printf("\n%d: GPU memory usage (dividing by powers of 1000): used = %f MB, free = %f MB, total = %f MB\n\n",myrank,
            used_db/1000.0/1000.0, free_db/1000.0/1000.0, total_db/1000.0/1000.0);

    // ordering mpi output
#ifdef WITH_MPI
      } //if
      // synchronizes mpi processes
      MPI_Barrier(MPI_COMM_WORLD);
    } //for
#endif

    // make sure that the device has compute capability >= 1.3
    if (deviceProp.major < 1){
      printf("Compute capability major number should be at least 1, exiting...\n\n");
      exit_on_error("CUDA Compute capability major number should be at least 1\n");
    }
    if (deviceProp.major == 1 && deviceProp.minor < 3){
      printf("Compute capability should be at least 1.3, exiting...\n");
      exit_on_error("CUDA Compute capability major number should be at least 1.3\n");
    }
    // we use pinned memory for asynchronous copy
    if( ! deviceProp.canMapHostMemory){
      printf("Device capability should allow to map host memory, exiting...\n");
      exit_on_error("CUDA Device capability canMapHostMemory should be TRUE\n");
    }

    // tests the device with a small memory allocation
    int size = 128;
    int* d_array;
    err = cudaMalloc((void**)&d_array,size*sizeof(int));
    if (err != cudaSuccess) {
      printf("Error testing memory allocation on device failed\n");
      printf("Error rank %d: cudaMalloc failed: %s\n", myrank,cudaGetErrorString(err));
      if (err == cudaErrorDevicesUnavailable){ printf("\n%s\n", err_info); }
      exit_on_error("CUDA runtime error: cudaMalloc failed\n\n");
    }
    err = cudaFree(d_array);
    if (err != cudaSuccess) {
      printf("Error cudaFree failed: %s\n", cudaGetErrorString(err));
      if (err == cudaErrorDevicesUnavailable){ printf("\n%s\n", err_info); }
      exit_on_error("CUDA runtime error: cudaFree failed\n\n");
    }

    // double check
    exit_on_gpu_error("cuda Malloc/Free test failed");

    // synchronizes GPU
#if CUDA_VERSION < 4000 || (defined (__CUDACC_VER_MAJOR__) && (__CUDACC_VER_MAJOR__ < 4))
    cudaThreadSynchronize();
#else
    cudaDeviceSynchronize();
#endif

    // synchronizes mpi processes
#ifdef WITH_MPI
    // output infos for mpi ranks ordered
    MPI_Barrier(MPI_COMM_WORLD);
    for(int iproc = 0;iproc < sizeprocs; iproc++){
      if (iproc == myrank){
        printf("rank %d on device %d okay\n",myrank,device);
        fflush(stdout);
      }
      // synchronizes mpi processes
      MPI_Barrier(MPI_COMM_WORLD);
    }
    sleep(1);
    if (myrank == 0){printf("\n\n");fflush(stdout);}
    MPI_Barrier(MPI_COMM_WORLD);
#endif
  } // loop device_count

}


int main(int argc, char **argv)
{
  int myrank,ndevices;

  // initialize
#ifdef WITH_MPI
  int size;
  int rc;
  rc = MPI_Init(&argc,&argv);
  if (rc != MPI_SUCCESS) {
    printf ("Error starting MPI program. Terminating.\n");
    MPI_Abort(MPI_COMM_WORLD, rc);
  }
  MPI_Comm_size(MPI_COMM_WORLD,&size);
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  if( myrank == 0 ){ printf ("Number of MPI processes = %d\n\n",size);fflush(stdout); }
  MPI_Barrier(MPI_COMM_WORLD);
#else
  myrank = 0;
#endif

  ndevices = 0;

  // initializes cuda devices
  initialize_cuda_device(&myrank,&ndevices);

  // releases previous contexts
#if CUDA_VERSION < 4000 || (defined (__CUDACC_VER_MAJOR__) && (__CUDACC_VER_MAJOR__ < 4))
  cudaThreadExit();
#else
  cudaDeviceReset();
#endif

  printf("number of total devices: %d\n\n",ndevices);

#ifdef WITH_MPI
  MPI_Finalize();
#endif
  return 0;
}

