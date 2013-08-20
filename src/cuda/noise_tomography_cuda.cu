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

#include <sys/types.h>
#include <unistd.h>

#include "config.h"
#include "mesh_constants_cuda.h"
// #include "epik_user.h"


/* ----------------------------------------------------------------------------------------------- */

extern "C"
void FC_FUNC_(fortranflush,FORTRANFLUSH)(int* rank){
TRACE("fortranflush");

  fflush(stdout);
  fflush(stderr);
  printf("Flushing proc %d!\n",*rank);
}

/* ----------------------------------------------------------------------------------------------- */

extern "C"
void FC_FUNC_(fortranprint,FORTRANPRINT)(int* id) {
TRACE("fortranprint");

  int procid;
#ifdef WITH_MPI
  MPI_Comm_rank(MPI_COMM_WORLD,&procid);
#else
  procid = 0;
#endif
  printf("%d: sends msg_id %d\n",procid,*id);
}

/* ----------------------------------------------------------------------------------------------- */

extern "C"
void FC_FUNC_(fortranprintf,FORTRANPRINTF)(realw* val) {
TRACE("fortranprintf");

  int procid;
#ifdef WITH_MPI
  MPI_Comm_rank(MPI_COMM_WORLD,&procid);
#else
  procid = 0;
#endif
  printf("%d: sends val %e\n",procid,*val);
}

/* ----------------------------------------------------------------------------------------------- */

extern "C"
void FC_FUNC_(fortranprintd,FORTRANPRINTD)(double* val) {
TRACE("fortranprintd");

  int procid;
#ifdef WITH_MPI
  MPI_Comm_rank(MPI_COMM_WORLD,&procid);
#else
  procid = 0;
#endif
  printf("%d: sends val %e\n",procid,*val);
}

/* ----------------------------------------------------------------------------------------------- */

// randomize displ for testing
extern "C"
void FC_FUNC_(make_displ_rand,MAKE_DISPL_RAND)(long* Mesh_pointer,realw* h_displ) {
TRACE("make_displ_rand");

  Mesh* mp = (Mesh*)(*Mesh_pointer); // get Mesh from fortran integer wrapper
  // realw* displ_rnd = (realw*)malloc(mp->NGLOB_AB*3*sizeof(realw));
  for(int i=0;i<mp->NGLOB_AB*3;i++) {
    h_displ[i] = rand();
  }
  print_CUDA_error_if_any(cudaMemcpy(mp->d_displ,h_displ,mp->NGLOB_AB*3*sizeof(realw),cudaMemcpyHostToDevice),44001);
}

/* ----------------------------------------------------------------------------------------------- */

__global__ void transfer_surface_to_host_kernel(int* free_surface_ispec,
                                                int* free_surface_ijk,
                                                int num_free_surface_faces,
                                                int* ibool,
                                                realw* displ,
                                                realw* noise_surface_movie) {
  int igll = threadIdx.x;
  int iface = blockIdx.x + blockIdx.y*gridDim.x;

  // int id = tx + blockIdx.x*blockDim.x + blockIdx.y*blockDim.x*gridDim.x;

  if(iface < num_free_surface_faces) {
    int ispec = free_surface_ispec[iface]-1; //-1 for C-based indexing

    int i = free_surface_ijk[INDEX3(NDIM,NGLL2,0,igll,iface)]-1;
    int j = free_surface_ijk[INDEX3(NDIM,NGLL2,1,igll,iface)]-1;
    int k = free_surface_ijk[INDEX3(NDIM,NGLL2,2,igll,iface)]-1;

    int iglob = ibool[INDEX4(5,5,5,i,j,k,ispec)]-1;

    noise_surface_movie[INDEX3(NDIM,NGLL2,0,igll,iface)] = displ[iglob*3];
    noise_surface_movie[INDEX3(NDIM,NGLL2,1,igll,iface)] = displ[iglob*3+1];
    noise_surface_movie[INDEX3(NDIM,NGLL2,2,igll,iface)] = displ[iglob*3+2];
  }
}

/* ----------------------------------------------------------------------------------------------- */

extern "C"
void FC_FUNC_(transfer_surface_to_host,
              TRANSFER_SURFACE_TO_HOST)(long* Mesh_pointer,
                                        realw* h_noise_surface_movie) {
TRACE("transfer_surface_to_host");

  Mesh* mp = (Mesh*)(*Mesh_pointer); // get Mesh from fortran integer wrapper

  int num_blocks_x, num_blocks_y;
  get_blocks_xy(mp->num_free_surface_faces,&num_blocks_x,&num_blocks_y);

  dim3 grid(num_blocks_x,num_blocks_y,1);
  dim3 threads(NGLL2,1,1);

  transfer_surface_to_host_kernel<<<grid,threads>>>(mp->d_free_surface_ispec,
                                                    mp->d_free_surface_ijk,
                                                    mp->num_free_surface_faces,
                                                    mp->d_ibool,
                                                    mp->d_displ,
                                                    mp->d_noise_surface_movie);

  print_CUDA_error_if_any(cudaMemcpy(h_noise_surface_movie,mp->d_noise_surface_movie,
                                     3*NGLL2*(mp->num_free_surface_faces)*sizeof(realw),cudaMemcpyDeviceToHost),44002);

#ifdef ENABLE_VERY_SLOW_ERROR_CHECKING
  exit_on_cuda_error("transfer_surface_to_host");
#endif
}

/* ----------------------------------------------------------------------------------------------- */

__global__ void noise_read_add_surface_movie_cuda_kernel(realw* accel, int* ibool,
                                                         int* free_surface_ispec,
                                                         int* free_surface_ijk,
                                                         int num_free_surface_faces,
                                                         realw* noise_surface_movie,
                                                         realw* normal_x_noise,
                                                         realw* normal_y_noise,
                                                         realw* normal_z_noise,
                                                         realw* mask_noise,
                                                         realw* free_surface_jacobian2Dw) {

  int iface = blockIdx.x + gridDim.x*blockIdx.y; // surface element id

  // when nspec_top > 65535, but mod(nspec_top,2) > 0, we end up with an extra block.
  if(iface < num_free_surface_faces) {
    int ispec = free_surface_ispec[iface]-1;

    int igll = threadIdx.x;

    int ipoin = NGLL2*iface + igll;
    int i=free_surface_ijk[INDEX3(NDIM,NGLL2,0,igll,iface)]-1;
    int j=free_surface_ijk[INDEX3(NDIM,NGLL2,1,igll,iface)]-1;
    int k=free_surface_ijk[INDEX3(NDIM,NGLL2,2,igll,iface)]-1;

    int iglob = ibool[INDEX4(5,5,5,i,j,k,ispec)]-1;

    realw normal_x = normal_x_noise[ipoin];
    realw normal_y = normal_y_noise[ipoin];
    realw normal_z = normal_z_noise[ipoin];

    realw eta = (noise_surface_movie[INDEX3(NDIM,NGLL2,0,igll,iface)]*normal_x +
                noise_surface_movie[INDEX3(NDIM,NGLL2,1,igll,iface)]*normal_y +
                noise_surface_movie[INDEX3(NDIM,NGLL2,2,igll,iface)]*normal_z);

    // error from cuda-memcheck and ddt seems "incorrect", because we
    // are passing a __constant__ variable pointer around like it was
    // made using cudaMalloc, which *may* be "incorrect", but produces
    // correct results.

    // ========= Invalid __global__ read of size
    // 4 ========= at 0x00000cd8 in
    // compute_add_sources_cuda.cu:260:noise_read_add_surface_movie_cuda_kernel
    // ========= by thread (0,0,0) in block (3443,0) ========= Address
    // 0x203000c8 is out of bounds

    // non atomic version for speed testing -- atomic updates are needed for correctness
    // accel[3*iglob] +=   eta*mask_noise[ipoin] * normal_x * wgllwgll_xy[tx] * free_surface_jacobian2Dw[tx + NGLL2*ispec2D];
    // accel[3*iglob+1] += eta*mask_noise[ipoin] * normal_y * wgllwgll_xy[tx] * free_surface_jacobian2Dw[tx + NGLL2*ispec2D];
    // accel[3*iglob+2] += eta*mask_noise[ipoin] * normal_z * wgllwgll_xy[tx] * free_surface_jacobian2Dw[tx + NGLL2*ispec2D];

    // Fortran version in SVN -- note deletion of wgllwgll_xy?
    // accel(1,iglob) = accel(1,iglob) + eta * mask_noise(ipoin) * normal_x_noise(ipoin) &
    // * free_surface_jacobian2Dw(igll,iface)
    // accel(2,iglob) = accel(2,iglob) + eta * mask_noise(ipoin) * normal_y_noise(ipoin) &
    // * free_surface_jacobian2Dw(igll,iface)
    // accel(3,iglob) = accel(3,iglob) + eta * mask_noise(ipoin) * normal_z_noise(ipoin) &
    // * free_surface_jacobian2Dw(igll,iface) ! wgllwgll_xy(i,j) * jacobian2D_top(i,j,iface)

    // atomicAdd(&accel[iglob*3]  ,eta*mask_noise[ipoin]*normal_x*wgllwgll_xy[tx]*free_surface_jacobian2Dw[igll+NGLL2*iface]);
    // atomicAdd(&accel[iglob*3+1],eta*mask_noise[ipoin]*normal_y*wgllwgll_xy[tx]*free_surface_jacobian2Dw[igll+NGLL2*iface]);
    // atomicAdd(&accel[iglob*3+2],eta*mask_noise[ipoin]*normal_z*wgllwgll_xy[tx]*free_surface_jacobian2Dw[igll+NGLL2*iface]);

    atomicAdd(&accel[iglob*3]  ,eta*mask_noise[ipoin]*normal_x*free_surface_jacobian2Dw[igll+NGLL2*iface]);
    atomicAdd(&accel[iglob*3+1],eta*mask_noise[ipoin]*normal_y*free_surface_jacobian2Dw[igll+NGLL2*iface]);
    atomicAdd(&accel[iglob*3+2],eta*mask_noise[ipoin]*normal_z*free_surface_jacobian2Dw[igll+NGLL2*iface]);

  }
}

/* ----------------------------------------------------------------------------------------------- */

extern "C"
void FC_FUNC_(noise_read_add_surface_movie_cu,
              NOISE_READ_ADD_SURFACE_MOVIE_CU)(long* Mesh_pointer,
                                               realw* h_noise_surface_movie,
                                               int* NOISE_TOMOGRAPHYf) {
  TRACE("noise_read_add_surface_movie_cu");

  // EPIK_TRACER("noise_read_add_surface_movie_cu");

  Mesh* mp = (Mesh*)(*Mesh_pointer); //get mesh pointer out of fortran integer container
  int NOISE_TOMOGRAPHY = *NOISE_TOMOGRAPHYf;

  print_CUDA_error_if_any(cudaMemcpy(mp->d_noise_surface_movie,h_noise_surface_movie,
                                     3*NGLL2*(mp->num_free_surface_faces)*sizeof(realw),cudaMemcpyHostToDevice),44003);

  int num_blocks_x, num_blocks_y;
  get_blocks_xy(mp->num_free_surface_faces,&num_blocks_x,&num_blocks_y);

  dim3 grid(num_blocks_x,num_blocks_y,1);
  dim3 threads(NGLL2,1,1);

  if(NOISE_TOMOGRAPHY == 2) { // add surface source to forward field
    noise_read_add_surface_movie_cuda_kernel<<<grid,threads>>>(mp->d_accel,
                                                               mp->d_ibool,
                                                               mp->d_free_surface_ispec,
                                                               mp->d_free_surface_ijk,
                                                               mp->num_free_surface_faces,
                                                               mp->d_noise_surface_movie,
                                                               mp->d_normal_x_noise,
                                                               mp->d_normal_y_noise,
                                                               mp->d_normal_z_noise,
                                                               mp->d_mask_noise,
                                                               mp->d_free_surface_jacobian2Dw);
  }
  else if(NOISE_TOMOGRAPHY == 3) { // add surface source to adjoint (backward) field
    noise_read_add_surface_movie_cuda_kernel<<<grid,threads>>>(mp->d_b_accel,
                                                               mp->d_ibool,
                                                               mp->d_free_surface_ispec,
                                                               mp->d_free_surface_ijk,
                                                               mp->num_free_surface_faces,
                                                               mp->d_noise_surface_movie,
                                                               mp->d_normal_x_noise,
                                                               mp->d_normal_y_noise,
                                                               mp->d_normal_z_noise,
                                                               mp->d_mask_noise,
                                                               mp->d_free_surface_jacobian2Dw);
  }

#ifdef ENABLE_VERY_SLOW_ERROR_CHECKING
  exit_on_cuda_error("noise_read_add_surface_movie_cuda_kernel");
#endif
}
