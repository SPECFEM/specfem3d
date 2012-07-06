/*
 !=====================================================================
 !
 !               S p e c f e m 3 D  V e r s i o n  2 . 0
 !               ---------------------------------------
 !
 !          Main authors: Dimitri Komatitsch and Jeroen Tromp
 !    Princeton University, USA and University of Pau / CNRS / INRIA
 ! (c) Princeton University / California Institute of Technology and University of Pau / CNRS / INRIA
 !                            April 2011
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
#include "prepare_constants_cuda.h"


/* ----------------------------------------------------------------------------------------------- */

// Check functions

/* ----------------------------------------------------------------------------------------------- */

extern "C"
void FC_FUNC_(check_max_norm_displ_gpu,
              CHECK_MAX_NORM_DISPL_GPU)(int* size, realw* displ,long* Mesh_pointer_f,int* announceID) {

TRACE("check_max_norm_displ_gpu");

  Mesh* mp = (Mesh*)(*Mesh_pointer_f); //get mesh pointer out of fortran integer container

  cudaMemcpy(displ, mp->d_displ,*size*sizeof(realw),cudaMemcpyDeviceToHost);
  realw maxnorm=0;

  for(int i=0;i<*size;i++) {
    maxnorm = MAX(maxnorm,fabsf(displ[i]));
  }
  printf("%d: maxnorm of forward displ = %e\n",*announceID,maxnorm);
}

/* ----------------------------------------------------------------------------------------------- */

extern "C"
void FC_FUNC_(check_max_norm_vector,
              CHECK_MAX_NORM_VECTOR)(int* size, realw* vector1, int* announceID) {

TRACE("check_max_norm_vector");

  int procid;
#ifdef WITH_MPI
  MPI_Comm_rank(MPI_COMM_WORLD,&procid);
#else
  procid = 0;
#endif
  realw maxnorm=0;
  int maxloc;
  for(int i=0;i<*size;i++) {
    if(maxnorm<fabsf(vector1[i])) {
      maxnorm = vector1[i];
      maxloc = i;
    }
  }
  printf("%d:maxnorm of vector %d [%d] = %e\n",procid,*announceID,maxloc,maxnorm);
}

/* ----------------------------------------------------------------------------------------------- */

extern "C"
void FC_FUNC_(check_max_norm_displ,
              CHECK_MAX_NORM_DISPL)(int* size, realw* displ, int* announceID) {

TRACE("check_max_norm_displ");

  realw maxnorm=0;

  for(int i=0;i<*size;i++) {
    maxnorm = MAX(maxnorm,fabsf(displ[i]));
  }
  printf("%d: maxnorm of forward displ = %e\n",*announceID,maxnorm);
}

/* ----------------------------------------------------------------------------------------------- */

extern "C"
void FC_FUNC_(check_max_norm_b_displ_gpu,
              CHECK_MAX_NORM_B_DISPL_GPU)(int* size, realw* b_displ,long* Mesh_pointer_f,int* announceID) {

TRACE("check_max_norm_b_displ_gpu");

  Mesh* mp = (Mesh*)(*Mesh_pointer_f); //get mesh pointer out of fortran integer container

  realw* b_accel = (realw*)malloc(*size*sizeof(realw));

  cudaMemcpy(b_displ, mp->d_b_displ,*size*sizeof(realw),cudaMemcpyDeviceToHost);
  cudaMemcpy(b_accel, mp->d_b_accel,*size*sizeof(realw),cudaMemcpyDeviceToHost);

  realw maxnorm=0;
  realw maxnorm_accel=0;

  for(int i=0;i<*size;i++) {
    maxnorm = MAX(maxnorm,fabsf(b_displ[i]));
    maxnorm_accel = MAX(maxnorm,fabsf(b_accel[i]));
  }
  free(b_accel);
  printf("%d: maxnorm of backward displ = %e\n",*announceID,maxnorm);
  printf("%d: maxnorm of backward accel = %e\n",*announceID,maxnorm_accel);
}

/* ----------------------------------------------------------------------------------------------- */

extern "C"
void FC_FUNC_(check_max_norm_b_accel_gpu,
              CHECK_MAX_NORM_B_ACCEL_GPU)(int* size, realw* b_accel,long* Mesh_pointer_f,int* announceID) {

TRACE("check_max_norm_b_accel_gpu");

  Mesh* mp = (Mesh*)(*Mesh_pointer_f); //get mesh pointer out of fortran integer container

  cudaMemcpy(b_accel, mp->d_b_accel,*size*sizeof(realw),cudaMemcpyDeviceToHost);

  realw maxnorm=0;

  for(int i=0;i<*size;i++) {
    maxnorm = MAX(maxnorm,fabsf(b_accel[i]));
  }
  printf("%d: maxnorm of backward accel = %e\n",*announceID,maxnorm);
}

/* ----------------------------------------------------------------------------------------------- */

extern "C"
void FC_FUNC_(check_max_norm_b_veloc_gpu,
              CHECK_MAX_NORM_B_VELOC_GPU)(int* size, realw* b_veloc,long* Mesh_pointer_f,int* announceID) {

TRACE("check_max_norm_b_veloc_gpu");

  Mesh* mp = (Mesh*)(*Mesh_pointer_f); //get mesh pointer out of fortran integer container

  cudaMemcpy(b_veloc, mp->d_b_veloc,*size*sizeof(realw),cudaMemcpyDeviceToHost);

  realw maxnorm=0;

  for(int i=0;i<*size;i++) {
    maxnorm = MAX(maxnorm,fabsf(b_veloc[i]));
  }
  printf("%d: maxnorm of backward veloc = %e\n",*announceID,maxnorm);
}

/* ----------------------------------------------------------------------------------------------- */

extern "C"
void FC_FUNC_(check_max_norm_b_displ,
              CHECK_MAX_NORM_B_DISPL)(int* size, realw* b_displ,int* announceID) {

TRACE("check_max_norm_b_displ");

  realw maxnorm=0;

  for(int i=0;i<*size;i++) {
    maxnorm = MAX(maxnorm,fabsf(b_displ[i]));
  }
  printf("%d:maxnorm of backward displ = %e\n",*announceID,maxnorm);
}

/* ----------------------------------------------------------------------------------------------- */

extern "C"
void FC_FUNC_(check_max_norm_b_accel,
              CHECK_MAX_NORM_B_ACCEL)(int* size, realw* b_accel,int* announceID) {

TRACE("check_max_norm_b_accel");

  realw maxnorm=0;

  for(int i=0;i<*size;i++) {
    maxnorm = MAX(maxnorm,fabsf(b_accel[i]));
  }
  printf("%d:maxnorm of backward accel = %e\n",*announceID,maxnorm);
}

/* ----------------------------------------------------------------------------------------------- */

extern "C"
void FC_FUNC_(check_error_vectors,
              CHECK_ERROR_VECTORS)(int* sizef, realw* vector1,realw* vector2) {

TRACE("check_error_vectors");

  int size = *sizef;

  double diff2 = 0;
  double sum = 0;
  double temp;
  double maxerr=0;
  int maxerrorloc;

  for(int i=0;i<size;++i) {
    temp = vector1[i]-vector2[i];
    diff2 += temp*temp;
    sum += vector1[i]*vector1[i];
    if(maxerr < fabsf(temp)) {
      maxerr = abs(temp);
      maxerrorloc = i;
    }
  }

  printf("rel error = %f, maxerr = %e @ %d\n",diff2/sum,maxerr,maxerrorloc);
  int myrank;
#ifdef WITH_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
#else
  myrank = 0;
#endif
  if(myrank == 0) {
    for(int i=maxerrorloc;i>maxerrorloc-5;i--) {
      printf("[%d]: %e vs. %e\n",i,vector1[i],vector2[i]);
    }
  }

}


/* ----------------------------------------------------------------------------------------------- */

// Auxiliary functions

/* ----------------------------------------------------------------------------------------------- */


/* ----------------------------------------------------------------------------------------------- */

extern "C"
void FC_FUNC_(get_max_accel,
              GET_MAX_ACCEL)(int* itf,int* sizef,long* Mesh_pointer) {

TRACE("get_max_accel");

  Mesh* mp = (Mesh*)(*Mesh_pointer);
  int procid;
#ifdef WITH_MPI
  MPI_Comm_rank(MPI_COMM_WORLD,&procid);
#else
  procid = 0;
#endif
  int size = *sizef;
  int it = *itf;
  realw* accel_cpy = (realw*)malloc(size*sizeof(realw));
  cudaMemcpy(accel_cpy,mp->d_accel,size*sizeof(realw),cudaMemcpyDeviceToHost);
  realw maxval=0;
  for(int i=0;i<size;++i) {
    maxval = MAX(maxval,accel_cpy[i]);
  }
  printf("%d/%d: max=%e\n",it,procid,maxval);
  free(accel_cpy);
}

/* ----------------------------------------------------------------------------------------------- */

// ACOUSTIC simulations

/* ----------------------------------------------------------------------------------------------- */

__global__ void get_maximum_kernel(realw* array, int size, realw* d_max){

  /* simplest version: uses only 1 thread
   realw max;
   max = 0;
   // finds maximum value in array
   if( size > 0 ){
   max = abs(array[0]);
   for( int i=1; i < size; i++){
   if( abs(array[i]) > max ) max = abs(array[i]);
   }
   }
   *d_max = max;
   */

  // reduction example:
  __shared__ realw sdata[256] ;

  // load shared mem
  unsigned int tid = threadIdx.x;
  unsigned int i = blockIdx.x*blockDim.x + threadIdx.x;

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
      if( sdata[tid] < sdata[tid + s] ) sdata[tid] = sdata[tid + s];
    }
    __syncthreads();
  }

  // write result for this block to global mem
  if (tid == 0) d_max[blockIdx.x] = sdata[0];

}

/* ----------------------------------------------------------------------------------------------- */

extern "C"
void FC_FUNC_(get_norm_acoustic_from_device,
              GET_NORM_ACOUSTIC_FROM_DEVICE)(realw* norm,
                                                  long* Mesh_pointer_f,
                                                  int* SIMULATION_TYPE) {

TRACE("get_norm_acoustic_from_device");
  //double start_time = get_time();

  Mesh* mp = (Mesh*)(*Mesh_pointer_f); //get mesh pointer out of fortran integer container
  realw max;
  realw *d_max;

  max = 0;

  /* way 1 : timing Elapsed time: 8.464813e-03
   realw* h_array;
   h_array = (realw*)calloc(mp->NGLOB_AB,sizeof(realw));

   print_CUDA_error_if_any(cudaMemcpy(h_array,mp->d_potential_dot_dot_acoustic,
   sizeof(realw)*(mp->NGLOB_AB),cudaMemcpyDeviceToHost),131);

   // finds maximum value in array
   max = h_array[0];
   for( int i=1; i < mp->NGLOB_AB; i++){
   if( abs(h_array[i]) > max ) max = abs(h_array[i]);
   }
   free(h_array);
   */

  /* way 2: timing Elapsed time: 8.818102e-02
   // launch simple kernel
   cudaMalloc((void**)&d_max,sizeof(realw));

   dim3 grid(1,1);
   dim3 threads(1,1,1);

   get_maximum_kernel<<<grid,threads>>>(mp->d_potential_dot_dot_acoustic,
   mp->NGLOB_AB,
   d_max);
   print_CUDA_error_if_any(cudaMemcpy(&max,d_max, sizeof(realw), cudaMemcpyDeviceToHost),222);

   cudaFree(d_max);
   */

  // way 2 b: timing Elapsed time: 1.236916e-03
  // launch simple reduction kernel
  realw* h_max;
  int blocksize = 256;

  int num_blocks_x = (int) ceil(mp->NGLOB_AB/blocksize);
  //printf("num_blocks_x %i \n",num_blocks_x);

  h_max = (realw*) calloc(num_blocks_x,sizeof(realw));
  cudaMalloc((void**)&d_max,num_blocks_x*sizeof(realw));

  dim3 grid(num_blocks_x,1);
  dim3 threads(blocksize,1,1);

  if(*SIMULATION_TYPE == 1 ){
    get_maximum_kernel<<<grid,threads>>>(mp->d_potential_dot_dot_acoustic,
                                         mp->NGLOB_AB,
                                         d_max);
  }

  if(*SIMULATION_TYPE == 3 ){
    get_maximum_kernel<<<grid,threads>>>(mp->d_b_potential_dot_dot_acoustic,
                                         mp->NGLOB_AB,
                                         d_max);
  }

  print_CUDA_error_if_any(cudaMemcpy(h_max,d_max,num_blocks_x*sizeof(realw),cudaMemcpyDeviceToHost),222);

  // determines max for all blocks
  max = h_max[0];
  for(int i=1;i<num_blocks_x;i++) {
    if( max < h_max[i]) max = h_max[i];
  }

  cudaFree(d_max);
  free(h_max);

  /* way 3: doesn't work properly...
   cublasStatus status;

   // Initialize CUBLAS
   status = cublasInit();
   if (status != CUBLAS_STATUS_SUCCESS) {
   fprintf (stderr, "!!!! CUBLAS initialization error\n");
   exit(1);
   }

   // cublas function: cublasIsamax
   //       finds the smallest index of the maximum magnitude element of single 
   //      precision vector x
   int incr = 1;
   int imax = 0;
   imax = cublasIsamax(mp->NGLOB_AB,(realw*)mp->d_potential_dot_dot_acoustic, incr);
   status= cublasGetError();
   if (status != CUBLAS_STATUS_SUCCESS) {
   fprintf (stderr, "!!!! CUBLAS error in cublasIsamax\n");
   exit(1);
   }

   print_CUDA_error_if_any(cudaMemcpy(&max,&(mp->d_potential_dot_dot_acoustic[imax]),
                      sizeof(realw), cudaMemcpyDeviceToHost),222);

   printf("maximum %i %i %f \n",mp->NGLOB_AB,imax,max);

   // Shutdown
   status = cublasShutdown();
   if (status != CUBLAS_STATUS_SUCCESS) {
   fprintf (stderr, "!!!! shutdown error (A)\n");
   exit(1);
   }

   */

  // return result
  *norm = max;

#ifdef ENABLE_VERY_SLOW_ERROR_CHECKING
  //double end_time = get_time();
  //printf("Elapsed time: %e\n",end_time-start_time);
  exit_on_cuda_error("after get_norm_acoustic_from_device");
#endif
}

/* ----------------------------------------------------------------------------------------------- */

// ELASTIC simulations

/* ----------------------------------------------------------------------------------------------- */

__global__ void get_maximum_vector_kernel(realw* array, int size, realw* d_max){

  // reduction example:
  __shared__ realw sdata[256] ;

  // load shared mem
  unsigned int tid = threadIdx.x;
  unsigned int i = blockIdx.x*blockDim.x + threadIdx.x;

  // loads values into shared memory: assume array is a vector array
  sdata[tid] = (i < size) ? sqrt(array[i*3]*array[i*3]
                                 + array[i*3+1]*array[i*3+1]
                                 + array[i*3+2]*array[i*3+2]) : 0.0 ;

  __syncthreads();

  // do reduction in shared mem
  for(unsigned int s=blockDim.x/2; s>0; s>>=1)
  {
    if (tid < s){
      // summation:
      //sdata[tid] += sdata[tid + s];
      // maximum:
      if( sdata[tid] < sdata[tid + s] ) sdata[tid] = sdata[tid + s];
    }
    __syncthreads();
  }

  // write result for this block to global mem
  if (tid == 0) d_max[blockIdx.x] = sdata[0];

}

/* ----------------------------------------------------------------------------------------------- */

extern "C"
void FC_FUNC_(get_norm_elastic_from_device,
              GET_NORM_ELASTIC_FROM_DEVICE)(realw* norm,
                                                 long* Mesh_pointer_f,
                                                 int* SIMULATION_TYPE) {

  TRACE("get_norm_elastic_from_device");
  //double start_time = get_time();

  Mesh* mp = (Mesh*)(*Mesh_pointer_f); //get mesh pointer out of fortran integer container
  realw max;
  realw *d_max;

  max = 0;

  // launch simple reduction kernel
  realw* h_max;
  int blocksize = 256;

  int num_blocks_x = (int) ceil(mp->NGLOB_AB/blocksize);
  //printf("num_blocks_x %i \n",num_blocks_x);

  h_max = (realw*) calloc(num_blocks_x,sizeof(realw));
  cudaMalloc((void**)&d_max,num_blocks_x*sizeof(realw));

  dim3 grid(num_blocks_x,1);
  dim3 threads(blocksize,1,1);

  if(*SIMULATION_TYPE == 1 ){
    get_maximum_vector_kernel<<<grid,threads>>>(mp->d_displ,
                                                mp->NGLOB_AB,
                                                d_max);
  }

  if(*SIMULATION_TYPE == 3 ){
    get_maximum_vector_kernel<<<grid,threads>>>(mp->d_b_displ,
                                                mp->NGLOB_AB,
                                                d_max);
  }

  print_CUDA_error_if_any(cudaMemcpy(h_max,d_max,num_blocks_x*sizeof(realw),cudaMemcpyDeviceToHost),222);

  // determines max for all blocks
  max = h_max[0];
  for(int i=1;i<num_blocks_x;i++) {
    if( max < h_max[i]) max = h_max[i];
  }

  cudaFree(d_max);
  free(h_max);

  // return result
  *norm = max;

#ifdef ENABLE_VERY_SLOW_ERROR_CHECKING
  //double end_time = get_time();
  //printf("Elapsed time: %e\n",end_time-start_time);
  exit_on_cuda_error("after get_norm_elastic_from_device");
#endif
}


