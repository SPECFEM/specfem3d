/*
**************************

check_hip_device utility

**************************

this utility program will output GPU device informations helpful for debugging HIP.


for compilation, see the command-line examples given here:

hipcc -o check_hip_device check_hip_device.cpp

execute by:

./check_hip_device

*/


// includes

#include <stdio.h>
#include <stdlib.h>

// HIP header file
#include <hip/hip_runtime.h>

/*----------------------------------------------------------------------------------------*/

//helper functions

void exit_on_error(const char* info) {
  printf("\nERROR: %s\n",info);
  fflush(stdout);

  // stops program
  exit(EXIT_FAILURE);
  return;
}


void check_hip_error(const char* info) {
  // checks error
  hipError_t err = hipGetLastError();
  if (err != hipSuccess){
    fprintf (stderr,"HIP error: %s\n", hipGetErrorString(err));

    printf("\nERROR: %s\n",info);
    fflush(stdout);

    // stops program
    exit(EXIT_FAILURE);
    return;
  }
}


void get_free_memory(double* free_db, double* used_db, double* total_db) {
  // gets memory usage in byte
  size_t free_byte ;
  size_t total_byte ;
  hipError_t status = hipMemGetInfo( &free_byte, &total_byte ) ;
  if (status != hipSuccess){
    printf("Error: hipMemGetInfo failed, %s \n", hipGetErrorString(status) );
    exit(EXIT_FAILURE);
  }

  *free_db = (double)free_byte ;
  *total_db = (double)total_byte ;
  *used_db = *total_db - *free_db ;
  return;
}


/*----------------------------------------------------------------------------------------*/

int main(int argc, char* const argv[]) {

  // gets platform info
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

  // output info
  fprintf (stdout,"HIP Driver Version / Runtime Version          %d.%d / %d.%d\n",
                  driverVersion / 1000, (driverVersion % 100) / 10,
                  runtimeVersion / 1000, (runtimeVersion % 100) / 10);
  fprintf (stdout,"HIP Device count: %d\n\n",device_count);

  // devices
  for (int i = 0; i < device_count; i++) {
    fprintf (stdout, "-----------\n");
    fprintf (stdout, "Device %i:\n",i);
    fprintf (stdout, "-----------\n");

    // sets device
    int device = i;
    hipSetDevice( device );
    check_hip_error("hipSetDevice failed");

    // double check that device was properly selected
    hipGetDevice(&device);
    check_hip_error("hipGetDevice failed");

    if (device != i ){
      printf("Error: hipSetDevice()=%d\n  hipGetDevice()=%d\n",i,device);
      exit_on_error("HIP set/get device error: device id conflict \n");
    }

    // get device properties
    struct hipDeviceProp_t deviceProp;

    hipGetDeviceProperties(&deviceProp,device);
    check_hip_error("hipGetDevicePropoerties failed");

    // display device properties
    fprintf (stdout, "Device Name = %s\n", deviceProp.name);
    fprintf (stdout, "memory:\n");
    fprintf (stdout, "  totalGlobalMem (in MB): %f\n",(unsigned long) deviceProp.totalGlobalMem / (1024.f * 1024.f));
    fprintf (stdout, "  totalGlobalMem (in GB): %f\n",(unsigned long) deviceProp.totalGlobalMem / (1024.f * 1024.f * 1024.f));
    fprintf (stdout, "  totalConstMem (in bytes): %lu\n",(unsigned long) deviceProp.totalConstMem); // seems to be same as GlobalMem
    //fprintf (stdout, "  Maximum 1D texture size (in bytes): %lu\n",(unsigned long) deviceProp.maxTexture1D); // not available?
    fprintf (stdout, "  sharedMemPerBlock (in bytes): %lu\n",(unsigned long) deviceProp.sharedMemPerBlock);
    fprintf (stdout, "  regsPerBlock (in bytes): %lu\n",(unsigned long) deviceProp.regsPerBlock);
    fprintf (stdout, "blocks:\n");
    fprintf (stdout, "  Maximum number of threads per block: %d\n",deviceProp.maxThreadsPerBlock);
    fprintf (stdout, "  Maximum size of each dimension of a block: %d x %d x %d\n",
                     deviceProp.maxThreadsDim[0],deviceProp.maxThreadsDim[1],deviceProp.maxThreadsDim[2]);
    fprintf (stdout, "  Maximum sizes of each dimension of a grid: %d x %d x %d\n",
                     deviceProp.maxGridSize[0],deviceProp.maxGridSize[1],deviceProp.maxGridSize[2]);
    fprintf (stdout, "features:\n");
    fprintf (stdout, "  Compute capability of the device = %d.%d\n", deviceProp.major, deviceProp.minor);
    fprintf (stdout, "  multiProcessorCount: %d\n",deviceProp.multiProcessorCount);
    if (deviceProp.canMapHostMemory){
      fprintf (stdout, "  canMapHostMemory: TRUE\n");
    }else{
      fprintf (stdout, "  canMapHostMemory: FALSE\n");
    }
    if (deviceProp.concurrentKernels){
      fprintf (stdout, "  concurrentKernels: TRUE\n");
    }else{
      fprintf (stdout, "  concurrentKernels: FALSE\n");
    }

    // outputs initial memory infos via hipMemGetInfo()
    // memory usage
    double free_db,used_db,total_db;
    get_free_memory(&free_db,&used_db,&total_db);

    fprintf(stdout, "memory used:\n");
    fprintf(stdout, "  GPU memory usage: used = %f MB, free = %f MB, total = %f MB\n",
                    used_db/1024.0/1024.0, free_db/1024.0/1024.0, total_db/1024.0/1024.0);
    fprintf (stdout, "\n\n");

    // free device
    hipDeviceReset();
  }
}

