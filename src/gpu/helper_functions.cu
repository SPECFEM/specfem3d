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

#include "mesh_constants_gpu.h"

/* ----------------------------------------------------------------------------------------------- */

// Helper functions

/* ----------------------------------------------------------------------------------------------- */

// copies integer array from CPU host to GPU device
void gpuCopy_todevice_int(void** d_array_addr_ptr,int* h_array,int size){

  TRACE("gpuCopy_todevice_int");

#ifdef USE_CUDA
  if (run_cuda) {
    // allocates memory on GPU
    //
    // note: cudaMalloc uses a double-pointer, such that it can return an error code in case it fails
    //          we thus pass the address to the pointer above (as void double-pointer) to have it
    //          pointing to the correct pointer of the array here
    print_CUDA_error_if_any(cudaMalloc((void**)d_array_addr_ptr,size*sizeof(int)),12001);

    // copies values onto GPU
    //
    // note: cudaMemcpy uses the pointer to the array, we thus re-cast the value of
    //          the double-pointer above to have the correct pointer to the array
    print_CUDA_error_if_any(cudaMemcpy((int*) *d_array_addr_ptr,h_array,size*sizeof(int),cudaMemcpyHostToDevice),12002);
  }
#endif
#ifdef USE_HIP
  if (run_hip) {
    // allocates memory on GPU
    print_HIP_error_if_any(hipMalloc((void**)d_array_addr_ptr,size*sizeof(int)),12001);
    // copies values onto GPU
    int *tmp_ptr = (int*) *d_array_addr_ptr;
    print_HIP_error_if_any(hipMemcpy( (void*) tmp_ptr, (void*) h_array,size*sizeof(int), hipMemcpyHostToDevice),12002);
  }
#endif
}

/* ----------------------------------------------------------------------------------------------- */

// copies integer array from CPU host to GPU device
void gpuCopy_todevice_realw(void** d_array_addr_ptr,realw* h_array,int size){

  TRACE("gpuCopy_todevice_realw");

#ifdef USE_CUDA
  if (run_cuda) {
    // allocates memory on GPU
    print_CUDA_error_if_any(cudaMalloc((void**)d_array_addr_ptr,size*sizeof(realw)),22001);

    // copies values onto GPU
    print_CUDA_error_if_any(cudaMemcpy((realw*) *d_array_addr_ptr,h_array,size*sizeof(realw),cudaMemcpyHostToDevice),22002);
  }
#endif
#ifdef USE_HIP
  if (run_hip) {
    // allocates memory on GPU
    print_HIP_error_if_any(hipMalloc((void**)d_array_addr_ptr,size*sizeof(realw)),22001);

    // copies values onto GPU
    print_HIP_error_if_any(hipMemcpy((realw*) *d_array_addr_ptr,h_array,size*sizeof(realw),hipMemcpyHostToDevice),22002);
  }
#endif
}

/*----------------------------------------------------------------------------------------------- */

void gpuFree (void *d_ptr) {

  TRACE ("gpuFree");

  // frees memory on GPU
#ifdef USE_CUDA
  if (run_cuda) { cudaFree(d_ptr); }
#endif
#ifdef USE_HIP
  if (run_hip) { hipFree(d_ptr); }
#endif
}

/* ----------------------------------------------------------------------------------------------- */
// GPU reset
/* ----------------------------------------------------------------------------------------------- */

void gpuReset() {
  // releases previous contexts

  // cuda version
#ifdef USE_CUDA
  if (run_cuda) {
#if CUDA_VERSION < 4000 || (defined (__CUDACC_VER_MAJOR__) && (__CUDACC_VER_MAJOR__ < 4))
    cudaThreadExit();
#else
    cudaDeviceReset();
#endif
  }
#endif

  // hip version
#ifdef USE_HIP
  if (run_hip) { hipDeviceReset(); }
#endif
}


/* ----------------------------------------------------------------------------------------------- */
// GPU synchronization
/* ----------------------------------------------------------------------------------------------- */

void gpuSynchronize() {
  // synchronizes device

  // cuda version
#ifdef USE_CUDA
  if (run_cuda) {
#if CUDA_VERSION < 4000 || (defined (__CUDACC_VER_MAJOR__) && (__CUDACC_VER_MAJOR__ < 4))
    cudaThreadSynchronize();
#else
    cudaDeviceSynchronize();
#endif
  }
#endif

  // hip version
#ifdef USE_HIP
  if (run_hip){ hipDeviceSynchronize(); }
#endif
}

/* ----------------------------------------------------------------------------------------------- */

void gpuStreamCreate(gpu_stream* stream_ptr) {
  // creates stream
#ifdef USE_CUDA
  if (run_cuda){ cudaStreamCreate(stream_ptr); }
#endif
#ifdef USE_HIP
  if (run_hip){ hipStreamCreate(stream_ptr); }
#endif
}


void gpuStreamSynchronize(gpu_stream stream) {
  // synchronizes stream
#ifdef USE_CUDA
  if (run_cuda){ cudaStreamSynchronize(stream); }
#endif
#ifdef USE_HIP
  if (run_hip){ hipStreamSynchronize(stream); }
#endif
}


/* ----------------------------------------------------------------------------------------------- */

// Timing helper functions

/* ----------------------------------------------------------------------------------------------- */

void start_timing_gpu(gpu_event* start,gpu_event* stop){

#ifdef USE_CUDA
  // creates & starts event
  cudaEventCreate(start);
  cudaEventCreate(stop);
  cudaEventRecord( *start, 0 );
#endif
#ifdef USE_HIP
  // creates & starts event
  hipEventCreate(start);
  hipEventCreate(stop);
  hipEventRecord( *start, 0);
#endif
}

/* ----------------------------------------------------------------------------------------------- */

void stop_timing_gpu(gpu_event* start,gpu_event* stop, const char* info_str){
  realw time;

#ifdef USE_CUDA
  // stops events
  cudaEventRecord( *stop, 0 );
  cudaEventSynchronize( *stop );
  cudaEventElapsedTime( &time, *start, *stop );
  cudaEventDestroy( *start );
  cudaEventDestroy( *stop );
#endif
#ifdef USE_HIP
  // stops events
  hipEventRecord( *stop, 0);
  hipEventSynchronize( *stop );
  hipEventElapsedTime( &time, *start, *stop );
  hipEventDestroy( *start );
  hipEventDestroy( *stop );
#endif

  // user output
  printf("%s: Execution Time = %f ms\n",info_str,time);
}

/* ----------------------------------------------------------------------------------------------- */

void stop_timing_gpu(gpu_event* start,gpu_event* stop, const char* info_str, realw* t){
  realw time;

#ifdef USE_CUDA
  // stops events
  cudaEventRecord( *stop, 0);
  cudaEventSynchronize( *stop );
  cudaEventElapsedTime( &time, *start, *stop );
  cudaEventDestroy( *start );
  cudaEventDestroy( *stop );
#endif
#ifdef USE_HIP
  // stops events
  hipEventRecord( *stop, 0);
  hipEventSynchronize( *stop );
  hipEventElapsedTime( &time, *start, *stop );
  hipEventDestroy( *start );
  hipEventDestroy( *stop );
#endif

  // user output
  printf("%s: Execution Time = %f ms\n",info_str,time);

  // returns time
  *t = time;
}


/* ----------------------------------------------------------------------------------------------- */

void exit_on_gpu_error(const char* kernel_name) {
  //check to catch errors from previous operations
  int error = 0;
  const char *strerr = NULL;

  // synchronizes GPU
  gpuSynchronize();

#ifdef USE_CUDA
  if (run_cuda) {
    cudaError_t err = cudaGetLastError();
    error = err != cudaSuccess;
    strerr = cudaGetErrorString(err);
  }
#endif
#ifdef USE_HIP
  if (run_hip) {
    hipError_t err = hipGetLastError();
    error = err != hipSuccess;
    strerr = hipGetErrorString(err);
  }
#endif

  if (error){
    fprintf(stderr,"GPU Error: after %s: %s\n", kernel_name, strerr);

    //debugging
    //pause_for_debugger(0);

    // outputs error file
    FILE* fp;
    int myrank;
    char filename[BUFSIZ];
#ifdef WITH_MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
#else
    myrank = 0;
#endif
    sprintf(filename,OUTPUT_FILES"/error_message_%06d.txt",myrank);
    fp = fopen(filename,"a+");
    if (fp != NULL){
      fprintf(fp,"GPU Error: after %s: %s\n", kernel_name, strerr);
      fclose(fp);
    }

    // stops program
    //free(kernel_name);
#ifdef WITH_MPI
    MPI_Abort(MPI_COMM_WORLD,1);
#endif
    exit(EXIT_FAILURE);
  }
}

/* ----------------------------------------------------------------------------------------------- */

void exit_on_error(const char* info) {
  printf("\nERROR: %s\n",info);
  fflush(stdout);

  // outputs error file
  FILE* fp;
  int myrank;
  char filename[BUFSIZ];
#ifdef WITH_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
#else
  myrank = 0;
#endif
  sprintf(filename,OUTPUT_FILES"/error_message_%06d.txt",myrank);
  fp = fopen(filename,"a+");
  if (fp != NULL){
    fprintf(fp,"ERROR: %s\n",info);
    fclose(fp);
  }

  // stops program
#ifdef WITH_MPI
  MPI_Abort(MPI_COMM_WORLD,1);
#endif
  //free(info);
  exit(EXIT_FAILURE);
}


/*----------------------------------------------------------------------------------------------- */

// additional helper functions

/*----------------------------------------------------------------------------------------------- */

double get_time_val() {
  struct timeval t;
  gettimeofday(&t, NULL);
  return t.tv_sec + t.tv_usec*1e-6;
}

/* ----------------------------------------------------------------------------------------------- */

void pause_for_debugger(int pause) {
  if (pause) {
    int myrank;
#ifdef WITH_MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
#else
    myrank = 0;
#endif
    printf("I'm rank %d\n",myrank);
    int i = 0;
    char hostname[256];
    gethostname(hostname, sizeof(hostname));
    printf("PID %d on %s:%d ready for attach\n", getpid(), hostname,myrank);

    FILE *file = fopen("./attach_gdb.txt","w+");
    if (file != NULL){
      fprintf(file,"PID %d on %s:%d ready for attach\n", getpid(), hostname,myrank);
      fclose(file);
    }

    fflush(stdout);
    while (0 == i)
      sleep(5);
  }
}

/* ----------------------------------------------------------------------------------------------- */

extern EXTERN_LANG
void FC_FUNC_(pause_for_debug,
              PAUSE_FOR_DEBUG)() {
  TRACE("pause_for_debug");

  pause_for_debugger(1);
}


/* ----------------------------------------------------------------------------------------------- */

// MPI synchronization

/* ----------------------------------------------------------------------------------------------- */

void synchronize_mpi () {
#ifdef WITH_MPI
  MPI_Barrier (MPI_COMM_WORLD);
#endif
}


/* ----------------------------------------------------------------------------------------------- */

// for debugging purposes, unused so far...

/* ----------------------------------------------------------------------------------------------- */

//extern EXTERN_LANG
//void FC_FUNC_(fortranflush,FORTRANFLUSH)(int* rank){
//TRACE("fortranflush");
//
//  fflush(stdout);
//  fflush(stderr);
//  printf("Flushing proc %d!\n",*rank);
//}
//
//extern EXTERN_LANG
//void FC_FUNC_(fortranprint,FORTRANPRINT)(int* id) {
//TRACE("fortranprint");
//
//  int procid;
//#ifdef WITH_MPI
//  MPI_Comm_rank(MPI_COMM_WORLD,&procid);
//#else
//  procid = 0;
//#endif
//  printf("%d: sends msg_id %d\n",procid,*id);
//}
//
//extern EXTERN_LANG
//void FC_FUNC_(fortranprintf,FORTRANPRINTF)(realw* val) {
//TRACE("fortranprintf");
//
//  int procid;
//#ifdef WITH_MPI
//  MPI_Comm_rank(MPI_COMM_WORLD,&procid);
//#else
//  procid = 0;
//#endif
//  printf("%d: sends val %e\n",procid,*val);
//}
//
//extern EXTERN_LANG
//void FC_FUNC_(fortranprintd,FORTRANPRINTD)(double* val) {
//TRACE("fortranprintd");
//
//  int procid;
//#ifdef WITH_MPI
//  MPI_Comm_rank(MPI_COMM_WORLD,&procid);
//#else
//  procid = 0;
//#endif
//  printf("%d: sends val %e\n",procid,*val);
//}
