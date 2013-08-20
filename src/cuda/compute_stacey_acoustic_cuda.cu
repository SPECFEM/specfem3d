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

__global__ void compute_stacey_acoustic_kernel(realw* potential_dot_acoustic,
                                               realw* potential_dot_dot_acoustic,
                                               int* abs_boundary_ispec,
                                               int* abs_boundary_ijk,
                                               realw* abs_boundary_jacobian2Dw,
                                               int* ibool,
                                               realw* rhostore,
                                               realw* kappastore,
                                               int* ispec_is_inner,
                                               int* ispec_is_acoustic,
                                               int phase_is_inner,
                                               int SIMULATION_TYPE,
                                               int SAVE_FORWARD,
                                               int num_abs_boundary_faces,
                                               realw* b_potential_dot_acoustic,
                                               realw* b_potential_dot_dot_acoustic,
                                               realw* b_absorb_potential,
                                               int gravity) {

  int igll = threadIdx.x;
  int iface = blockIdx.x + gridDim.x*blockIdx.y;

  int i,j,k,iglob,ispec;
  realw rhol,kappal,cpl;
  realw jacobianw;
  realw vel;

  // don't compute points outside NGLLSQUARE==NGLL2==25
  // way 2: no further check needed since blocksize = 25
  if( iface < num_abs_boundary_faces){

  //  if(igll<NGLL2 && iface < num_abs_boundary_faces) {

    // "-1" from index values to convert from Fortran-> C indexing
    ispec = abs_boundary_ispec[iface]-1;

    if(ispec_is_inner[ispec] == phase_is_inner && ispec_is_acoustic[ispec] ) {

      i = abs_boundary_ijk[INDEX3(NDIM,NGLL2,0,igll,iface)]-1;
      j = abs_boundary_ijk[INDEX3(NDIM,NGLL2,1,igll,iface)]-1;
      k = abs_boundary_ijk[INDEX3(NDIM,NGLL2,2,igll,iface)]-1;
      iglob = ibool[INDEX4(5,5,5,i,j,k,ispec)]-1;

      // determines bulk sound speed
      rhol = rhostore[INDEX4_PADDED(NGLLX,NGLLX,NGLLX,i,j,k,ispec)];

      kappal = kappastore[INDEX4(5,5,5,i,j,k,ispec)];

      cpl = sqrt( kappal / rhol );

      // velocity
      if( gravity ){
        // daniel: TODO - check gravity and stacey condition here...
        // uses a potential definition of: s = grad(chi)
        vel = potential_dot_acoustic[iglob] / rhol ;
      }else{
        // uses a potential definition of: s = 1/rho grad(chi)
        vel = potential_dot_acoustic[iglob] / rhol;
      }

      // gets associated, weighted jacobian
      jacobianw = abs_boundary_jacobian2Dw[INDEX2(NGLL2,igll,iface)];

      // Sommerfeld condition
      atomicAdd(&potential_dot_dot_acoustic[iglob],-vel*jacobianw/cpl);

      // adjoint simulations
      if( SIMULATION_TYPE == 3 ){
        // Sommerfeld condition
        atomicAdd(&b_potential_dot_dot_acoustic[iglob],-b_absorb_potential[INDEX2(NGLL2,igll,iface)]);
      }else if( SIMULATION_TYPE == 1 && SAVE_FORWARD ){
        // saves boundary values
        b_absorb_potential[INDEX2(NGLL2,igll,iface)] = vel*jacobianw/cpl;
      }
    }
//  }
  }
}

/* ----------------------------------------------------------------------------------------------- */

extern "C"
void FC_FUNC_(compute_stacey_acoustic_cuda,
              COMPUTE_STACEY_ACOUSTIC_CUDA)(long* Mesh_pointer,
                                            int* phase_is_innerf,
                                            realw* h_b_absorb_potential) {
TRACE("compute_stacey_acoustic_cuda");
  //double start_time = get_time();

  Mesh* mp = (Mesh*)(*Mesh_pointer); //get mesh pointer out of fortran integer container

  // checks if anything to do
  if( mp->d_num_abs_boundary_faces == 0 ) return;

  int phase_is_inner          = *phase_is_innerf;

  // way 1: Elapsed time: 4.385948e-03
  // > NGLLSQUARE==NGLL2==25, but we handle this inside kernel
  //  int blocksize = 32;

  // way 2: Elapsed time: 4.379034e-03
  // > NGLLSQUARE==NGLL2==25, no further check inside kernel
  int blocksize = NGLL2;

  int num_blocks_x, num_blocks_y;
  get_blocks_xy(mp->d_num_abs_boundary_faces,&num_blocks_x,&num_blocks_y);

  dim3 grid(num_blocks_x,num_blocks_y);
  dim3 threads(blocksize,1,1);

  //  adjoint simulations: reads in absorbing boundary
  if (mp->simulation_type == 3 ){
    // copies array to GPU
    print_CUDA_error_if_any(cudaMemcpy(mp->d_b_absorb_potential,h_b_absorb_potential,
                                       mp->d_b_reclen_potential,cudaMemcpyHostToDevice),7700);
  }

  compute_stacey_acoustic_kernel<<<grid,threads>>>(mp->d_potential_dot_acoustic,
                                                   mp->d_potential_dot_dot_acoustic,
                                                   mp->d_abs_boundary_ispec,
                                                   mp->d_abs_boundary_ijk,
                                                   mp->d_abs_boundary_jacobian2Dw,
                                                   mp->d_ibool,
                                                   mp->d_rhostore,
                                                   mp->d_kappastore,
                                                   mp->d_ispec_is_inner,
                                                   mp->d_ispec_is_acoustic,
                                                   phase_is_inner,
                                                   mp->simulation_type,
                                                   mp->save_forward,
                                                   mp->d_num_abs_boundary_faces,
                                                   mp->d_b_potential_dot_acoustic,
                                                   mp->d_b_potential_dot_dot_acoustic,
                                                   mp->d_b_absorb_potential,
                                                   mp->gravity);

  //  adjoint simulations: stores absorbed wavefield part
  if (mp->simulation_type == 1 && mp->save_forward ){
    // (cudaMemcpy implicitly synchronizes all other cuda operations)
    // copies array to CPU
    print_CUDA_error_if_any(cudaMemcpy(h_b_absorb_potential,mp->d_b_absorb_potential,
                                       mp->d_b_reclen_potential,cudaMemcpyDeviceToHost),7701);
  }

#ifdef ENABLE_VERY_SLOW_ERROR_CHECKING
  exit_on_cuda_error("compute_stacey_acoustic_kernel");
#endif
}

