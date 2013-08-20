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

#include <sys/time.h>
#include <sys/resource.h>

#include "config.h"
#include "mesh_constants_cuda.h"


/* ----------------------------------------------------------------------------------------------- */

__global__ void compute_stacey_elastic_kernel(realw* veloc,
                                              realw* accel,
                                              int* abs_boundary_ispec,
                                              int* abs_boundary_ijk,
                                              realw* abs_boundary_normal,
                                              realw* abs_boundary_jacobian2Dw,
                                              int* ibool,
                                              realw* rho_vp,
                                              realw* rho_vs,
                                              int* ispec_is_inner,
                                              int* ispec_is_elastic,
                                              int phase_is_inner,
                                              int SIMULATION_TYPE,
                                              int SAVE_FORWARD,
                                              int num_abs_boundary_faces,
                                              realw* b_absorb_field) {

  int igll = threadIdx.x; // tx
  int iface = blockIdx.x + gridDim.x*blockIdx.y; // bx

  int i,j,k,iglob,ispec;
  realw vx,vy,vz,vn;
  realw nx,ny,nz;
  realw rho_vp_temp,rho_vs_temp;
  realw tx,ty,tz;
  realw jacobianw;

  // don't compute points outside NGLLSQUARE==NGLL2==25
  // way 2: no further check needed since blocksize = 25
  if( iface < num_abs_boundary_faces){

  //if(igll < NGLL2 && iface < num_abs_boundary_faces) {

    // "-1" from index values to convert from Fortran-> C indexing
    ispec = abs_boundary_ispec[iface]-1;

    if(ispec_is_inner[ispec] == phase_is_inner && ispec_is_elastic[ispec] ) {

      i = abs_boundary_ijk[INDEX3(NDIM,NGLL2,0,igll,iface)]-1;
      j = abs_boundary_ijk[INDEX3(NDIM,NGLL2,1,igll,iface)]-1;
      k = abs_boundary_ijk[INDEX3(NDIM,NGLL2,2,igll,iface)]-1;
      iglob = ibool[INDEX4(NGLLX,NGLLX,NGLLX,i,j,k,ispec)]-1;

      // gets associated velocity

      vx = veloc[iglob*3+0];
      vy = veloc[iglob*3+1];
      vz = veloc[iglob*3+2];

      // gets associated normal
      nx = abs_boundary_normal[INDEX3(NDIM,NGLL2,0,igll,iface)];
      ny = abs_boundary_normal[INDEX3(NDIM,NGLL2,1,igll,iface)];
      nz = abs_boundary_normal[INDEX3(NDIM,NGLL2,2,igll,iface)];

      // // velocity component in normal direction (normal points out of element)
      vn = vx*nx + vy*ny + vz*nz;

      rho_vp_temp = rho_vp[INDEX4(NGLLX,NGLLX,NGLLX,i,j,k,ispec)];
      rho_vs_temp = rho_vs[INDEX4(NGLLX,NGLLX,NGLLX,i,j,k,ispec)];

      tx = rho_vp_temp*vn*nx + rho_vs_temp*(vx-vn*nx);
      ty = rho_vp_temp*vn*ny + rho_vs_temp*(vy-vn*ny);
      tz = rho_vp_temp*vn*nz + rho_vs_temp*(vz-vn*nz);

      jacobianw = abs_boundary_jacobian2Dw[INDEX2(NGLL2,igll,iface)];

      atomicAdd(&accel[iglob*3],-tx*jacobianw);
      atomicAdd(&accel[iglob*3+1],-ty*jacobianw);
      atomicAdd(&accel[iglob*3+2],-tz*jacobianw);

      if(SAVE_FORWARD && SIMULATION_TYPE == 1) {
        b_absorb_field[INDEX3(NDIM,NGLL2,0,igll,iface)] = tx*jacobianw;
        b_absorb_field[INDEX3(NDIM,NGLL2,1,igll,iface)] = ty*jacobianw;
        b_absorb_field[INDEX3(NDIM,NGLL2,2,igll,iface)] = tz*jacobianw;
      } // SIMULATION_TYPE
    }
  } // num_abs_boundary_faces

}

/* ----------------------------------------------------------------------------------------------- */

__global__ void compute_stacey_elastic_sim3_kernel(int* abs_boundary_ispec,
                                                   int* abs_boundary_ijk,
                                                   int* ibool,
                                                   int* ispec_is_inner,
                                                   int* ispec_is_elastic,
                                                   int phase_is_inner,
                                                   int num_abs_boundary_faces,
                                                   realw* b_accel,
                                                   realw* b_absorb_field) {

  int igll = threadIdx.x; // tx
  int iface = blockIdx.x + gridDim.x*blockIdx.y; // bx

  int i,j,k,iglob,ispec;

  // don't compute points outside NGLLSQUARE==NGLL2==25
  // way 2: no further check needed since blocksize = 25
  if( iface < num_abs_boundary_faces){

  //if(igll < NGLL2 && iface < num_abs_boundary_faces) {

    // "-1" from index values to convert from Fortran-> C indexing
    ispec = abs_boundary_ispec[iface]-1;

    if(ispec_is_inner[ispec] == phase_is_inner && ispec_is_elastic[ispec] ) {

      i = abs_boundary_ijk[INDEX3(NDIM,NGLL2,0,igll,iface)]-1;
      j = abs_boundary_ijk[INDEX3(NDIM,NGLL2,1,igll,iface)]-1;
      k = abs_boundary_ijk[INDEX3(NDIM,NGLL2,2,igll,iface)]-1;

      iglob = ibool[INDEX4(NGLLX,NGLLX,NGLLX,i,j,k,ispec)]-1;

      atomicAdd(&b_accel[iglob*3  ],-b_absorb_field[INDEX3(NDIM,NGLL2,0,igll,iface)]);
      atomicAdd(&b_accel[iglob*3+1],-b_absorb_field[INDEX3(NDIM,NGLL2,1,igll,iface)]);
      atomicAdd(&b_accel[iglob*3+2],-b_absorb_field[INDEX3(NDIM,NGLL2,2,igll,iface)]);
    }
  } // num_abs_boundary_faces

}

/* ----------------------------------------------------------------------------------------------- */


extern "C"
void FC_FUNC_(compute_stacey_viscoelastic_cuda,
              COMPUTE_STACEY_VISCOELASTIC_CUDA)(long* Mesh_pointer,
                                           int* phase_is_innerf,
                                           realw* b_absorb_field) {

  TRACE("\tcompute_stacey_viscoelastic_cuda");

  Mesh* mp = (Mesh*)(*Mesh_pointer); //get mesh pointer out of fortran integer container

  // checks if anything to do
  if( mp->d_num_abs_boundary_faces == 0 ) return;

  int phase_is_inner    = *phase_is_innerf;

  // way 1
  // > NGLLSQUARE==NGLL2==25, but we handle this inside kernel
  //int blocksize = 32;

  // way 2: seems sligthly faster
  // > NGLLSQUARE==NGLL2==25, no further check inside kernel
  int blocksize = NGLL2;

  int num_blocks_x, num_blocks_y;
  get_blocks_xy(mp->d_num_abs_boundary_faces,&num_blocks_x,&num_blocks_y);

  dim3 grid(num_blocks_x,num_blocks_y);
  dim3 threads(blocksize,1,1);

  if(mp->simulation_type == 3 ) {
    // reading is done in fortran routine
    print_CUDA_error_if_any(cudaMemcpy(mp->d_b_absorb_field,b_absorb_field,
                                       mp->d_b_reclen_field,cudaMemcpyHostToDevice),7700);
  }

#ifdef ENABLE_VERY_SLOW_ERROR_CHECKING
  exit_on_cuda_error("between cudamemcpy and compute_stacey_elastic_kernel");
#endif

  compute_stacey_elastic_kernel<<<grid,threads,0,mp->compute_stream>>>(mp->d_veloc,
                                                  mp->d_accel,
                                                  mp->d_abs_boundary_ispec,
                                                  mp->d_abs_boundary_ijk,
                                                  mp->d_abs_boundary_normal,
                                                  mp->d_abs_boundary_jacobian2Dw,
                                                  mp->d_ibool,
                                                  mp->d_rho_vp,
                                                  mp->d_rho_vs,
                                                  mp->d_ispec_is_inner,
                                                  mp->d_ispec_is_elastic,
                                                  phase_is_inner,
                                                  mp->simulation_type,
                                                  mp->save_forward,
                                                  mp->d_num_abs_boundary_faces,
                                                  mp->d_b_absorb_field);

  // adjoint simulations
  if(mp->simulation_type == 3 ){
    compute_stacey_elastic_sim3_kernel<<<grid,threads,0,mp->compute_stream>>>(mp->d_abs_boundary_ispec,
                                                         mp->d_abs_boundary_ijk,
                                                         mp->d_ibool,
                                                         mp->d_ispec_is_inner,
                                                         mp->d_ispec_is_elastic,
                                                         phase_is_inner,
                                                         mp->d_num_abs_boundary_faces,
                                                         mp->d_b_accel,
                                                         mp->d_b_absorb_field);
  }

#ifdef ENABLE_VERY_SLOW_ERROR_CHECKING
  exit_on_cuda_error("compute_stacey_elastic_kernel");
#endif

  // ! adjoint simulations: stores absorbed wavefield part
  // if (mp->simulation_type == 1 .and. SAVE_FORWARD .and. num_abs_boundary_faces > 0 ) &
  //   write(IOABS,rec=it) b_reclen_field,b_absorb_field,b_reclen_field

  if(mp->simulation_type == 1 && mp->save_forward ) {
    // explicitly wait until compute stream is done
    // (cudaMemcpy implicitly synchronizes all other cuda operations)
    cudaStreamSynchronize(mp->compute_stream);

    // copies absorb_field values to CPU
    print_CUDA_error_if_any(cudaMemcpy(b_absorb_field,mp->d_b_absorb_field,
                                       mp->d_b_reclen_field,cudaMemcpyDeviceToHost),7701);
    // writing is done in fortran routine
  }

#ifdef ENABLE_VERY_SLOW_ERROR_CHECKING
  exit_on_cuda_error("after compute_stacey_elastic after cudamemcpy");
#endif
}

