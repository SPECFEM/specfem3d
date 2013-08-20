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

#include <sys/types.h>
#include <unistd.h>
#include <sys/time.h>
#include <sys/resource.h>

#include "config.h"
#include "mesh_constants_cuda.h"


/* ----------------------------------------------------------------------------------------------- */

// ELASTIC SIMULATIONS

/* ----------------------------------------------------------------------------------------------- */

__global__ void compute_kernels_ani_cudakernel(int* ispec_is_elastic,
                                           int* ibool,
                                           realw* accel,
                                           realw* b_displ,
                                           realw* epsilondev_xx,realw* epsilondev_yy,realw* epsilondev_xy,
                                           realw* epsilondev_xz,realw* epsilondev_yz,
                                           realw* b_epsilondev_xx,realw* b_epsilondev_yy,realw* b_epsilondev_xy,
                                           realw* b_epsilondev_xz,realw* b_epsilondev_yz,
                                           realw* rho_kl,
                                           realw deltat,
                                           realw* cijkl_kl,
                                           realw* epsilon_trace_over_3,
                                           realw* b_epsilon_trace_over_3,
                                           int NSPEC_AB) {

  int ispec = blockIdx.x + blockIdx.y*gridDim.x;
  int ijk = threadIdx.x;
  int ijk_ispec = ijk + NGLL3*ispec;
  int ijk21_ispec = ijk + 21*NGLL3*ispec;

  realw prod[21];
  realw eps[6];
  realw b_eps[6];
  realw epsdev[6];
  realw b_epsdev[6];
  realw eps_trace_over_3,b_eps_trace_over_3;
  int i,j;

  // handles case when there is 1 extra block (due to rectangular grid)
  if(ispec < NSPEC_AB) {

    // elastic elements only
    if( ispec_is_elastic[ispec] ) {
      int iglob = ibool[ijk_ispec] - 1 ;

      // anisotropic kernels:
      // density kernel
      rho_kl[ijk_ispec] += deltat * (accel[3*iglob]*b_displ[3*iglob]+
                                     accel[3*iglob+1]*b_displ[3*iglob+1]+
                                     accel[3*iglob+2]*b_displ[3*iglob+2]);


      // anisotropic kernel
      epsdev[0] = epsilondev_xx[ijk_ispec];
      epsdev[1] = epsilondev_yy[ijk_ispec];
      epsdev[2] = epsilondev_xy[ijk_ispec];
      epsdev[3] = epsilondev_xz[ijk_ispec];
      epsdev[4] = epsilondev_yz[ijk_ispec];

      b_epsdev[0] = b_epsilondev_xx[ijk_ispec];
      b_epsdev[1] = b_epsilondev_yy[ijk_ispec];
      b_epsdev[2] = b_epsilondev_xy[ijk_ispec];
      b_epsdev[3] = b_epsilondev_xz[ijk_ispec];
      b_epsdev[4] = b_epsilondev_yz[ijk_ispec];

      eps_trace_over_3 = epsilon_trace_over_3[ijk_ispec];
      b_eps_trace_over_3 = b_epsilon_trace_over_3[ijk_ispec];

      //! Building of the local matrix of the strain tensor
      //! for the adjoint field and the regular backward field
      //!eps11 et eps22
      eps[0] = epsdev[0] + eps_trace_over_3;
      eps[1] = epsdev[1] + eps_trace_over_3;
      //!eps33
      eps[2] = -(eps[0]+eps[1])+3*eps_trace_over_3;
      //!eps23
      eps[3] = epsdev[4];
      //!eps13
      eps[4] = epsdev[3];
      //!eps12
      eps[5] = epsdev[2];

      // backward arrays
      b_eps[0] = b_epsdev[0] + b_eps_trace_over_3;
      b_eps[1] = b_epsdev[1] + b_eps_trace_over_3;
      b_eps[2] = -(b_eps[0]+b_eps[1])+3*b_eps_trace_over_3;
      b_eps[3] = b_epsdev[4];
      b_eps[4] = b_epsdev[3];
      b_eps[5] = b_epsdev[2];

      //! Computing the 21 strain products without assuming eps(i)*b_eps(j) = eps(j)*b_eps(i)
      int p = 0;
      for( i=0; i<6; i++){
        for( j=i; j<6; j++){
          prod[p] = eps[i] * b_eps[j];
          if( j > i ){
            prod[p] = prod[p] + eps[j]*b_eps[i];
            if( j > 2 && i < 3 ){ prod[p] = prod[p]*2; }
          }
          if(i > 2 ){ prod[p] = prod[p]*4; }
          p++;
        }
      }

      // all 21 anisotropic coefficients
      for( i=0; i<21; i++){
        cijkl_kl[i+ijk21_ispec] += deltat * prod[i];
      }

    }
  }
}

/* ----------------------------------------------------------------------------------------------- */

__global__ void compute_kernels_cudakernel(int* ispec_is_elastic,
                                           int* ibool,
                                           realw* accel,
                                           realw* b_displ,
                                           realw* epsilondev_xx,realw* epsilondev_yy,realw* epsilondev_xy,
                                           realw* epsilondev_xz,realw* epsilondev_yz,
                                           realw* b_epsilondev_xx,realw* b_epsilondev_yy,realw* b_epsilondev_xy,
                                           realw* b_epsilondev_xz,realw* b_epsilondev_yz,
                                           realw* rho_kl,
                                           realw deltat,
                                           realw* mu_kl,
                                           realw* kappa_kl,
                                           realw* epsilon_trace_over_3,
                                           realw* b_epsilon_trace_over_3,
                                           int NSPEC_AB) {

  int ispec = blockIdx.x + blockIdx.y*gridDim.x;
  int ijk = threadIdx.x;
  int ijk_ispec = ijk + NGLL3*ispec;

  // handles case when there is 1 extra block (due to rectangular grid)
  if(ispec < NSPEC_AB) {

    // elastic elements only
    if( ispec_is_elastic[ispec] ) {
      int iglob = ibool[ijk_ispec] - 1 ;

      // isotropic kernels:
      // density kernel
      rho_kl[ijk_ispec] += deltat * (accel[3*iglob]*b_displ[3*iglob]+
                                     accel[3*iglob+1]*b_displ[3*iglob+1]+
                                     accel[3*iglob+2]*b_displ[3*iglob+2]);


      // shear modulus kernel
      mu_kl[ijk_ispec] += deltat * (epsilondev_xx[ijk_ispec]*b_epsilondev_xx[ijk_ispec]+
                                    epsilondev_yy[ijk_ispec]*b_epsilondev_yy[ijk_ispec]+
                                    (epsilondev_xx[ijk_ispec]+epsilondev_yy[ijk_ispec])*
                                    (b_epsilondev_xx[ijk_ispec]+b_epsilondev_yy[ijk_ispec])+
                                    2*(epsilondev_xy[ijk_ispec]*b_epsilondev_xy[ijk_ispec]+
                                       epsilondev_xz[ijk_ispec]*b_epsilondev_xz[ijk_ispec]+
                                       epsilondev_yz[ijk_ispec]*b_epsilondev_yz[ijk_ispec]));

      // bulk modulus kernel
      kappa_kl[ijk_ispec] += deltat*(9*epsilon_trace_over_3[ijk_ispec]*
                                     b_epsilon_trace_over_3[ijk_ispec]);

    }
  }
}

/* ----------------------------------------------------------------------------------------------- */

extern "C"
void FC_FUNC_(compute_kernels_elastic_cuda,
              COMPUTE_KERNELS_ELASTIC_CUDA)(long* Mesh_pointer,
                                            realw* deltat_f) {

  TRACE("compute_kernels_elastic_cuda");

  Mesh* mp = (Mesh*)(*Mesh_pointer); //get mesh pointer out of fortran integer container

  int blocksize = NGLL3; // NGLLX*NGLLY*NGLLZ
  realw deltat = *deltat_f;

  int num_blocks_x, num_blocks_y;
  get_blocks_xy(mp->NSPEC_AB,&num_blocks_x,&num_blocks_y);

  dim3 grid(num_blocks_x,num_blocks_y);
  dim3 threads(blocksize,1,1);

  if( mp->anisotropic_kl ){
    compute_kernels_ani_cudakernel<<<grid,threads>>>(mp->d_ispec_is_elastic,mp->d_ibool,
                                                     mp->d_accel, mp->d_b_displ,
                                                     mp->d_epsilondev_xx,
                                                     mp->d_epsilondev_yy,
                                                     mp->d_epsilondev_xy,
                                                     mp->d_epsilondev_xz,
                                                     mp->d_epsilondev_yz,
                                                     mp->d_b_epsilondev_xx,
                                                     mp->d_b_epsilondev_yy,
                                                     mp->d_b_epsilondev_xy,
                                                     mp->d_b_epsilondev_xz,
                                                     mp->d_b_epsilondev_yz,
                                                     mp->d_rho_kl,
                                                     deltat,
                                                     mp->d_cijkl_kl,
                                                     mp->d_epsilon_trace_over_3,
                                                     mp->d_b_epsilon_trace_over_3,
                                                     mp->NSPEC_AB);

  }else{
    compute_kernels_cudakernel<<<grid,threads>>>(mp->d_ispec_is_elastic,mp->d_ibool,
                                                 mp->d_accel, mp->d_b_displ,
                                                 mp->d_epsilondev_xx,
                                                 mp->d_epsilondev_yy,
                                                 mp->d_epsilondev_xy,
                                                 mp->d_epsilondev_xz,
                                                 mp->d_epsilondev_yz,
                                                 mp->d_b_epsilondev_xx,
                                                 mp->d_b_epsilondev_yy,
                                                 mp->d_b_epsilondev_xy,
                                                 mp->d_b_epsilondev_xz,
                                                 mp->d_b_epsilondev_yz,
                                                 mp->d_rho_kl,
                                                 deltat,
                                                 mp->d_mu_kl,
                                                 mp->d_kappa_kl,
                                                 mp->d_epsilon_trace_over_3,
                                                 mp->d_b_epsilon_trace_over_3,
                                                 mp->NSPEC_AB);
  }

#ifdef ENABLE_VERY_SLOW_ERROR_CHECKING
  exit_on_cuda_error("compute_kernels_elastic_cuda");
#endif
}


/* ----------------------------------------------------------------------------------------------- */

// NOISE SIMULATIONS

/* ----------------------------------------------------------------------------------------------- */


__global__ void compute_kernels_strength_noise_cuda_kernel(realw* displ,
                                                           int* free_surface_ispec,
                                                           int* free_surface_ijk,
                                                           int* ibool,
                                                           realw* noise_surface_movie,
                                                           realw* normal_x_noise,
                                                           realw* normal_y_noise,
                                                           realw* normal_z_noise,
                                                           realw* Sigma_kl,
                                                           realw deltat,
                                                           int num_free_surface_faces) {
  int iface = blockIdx.x + blockIdx.y*gridDim.x;
  int igll = threadIdx.x;
  int ipoin = igll + NGLL2*iface;

  if(iface < num_free_surface_faces) {

    int ispec = free_surface_ispec[iface]-1;

    int i = free_surface_ijk[INDEX3(NDIM,NGLL2,0,igll,iface)] - 1 ;
    int j = free_surface_ijk[INDEX3(NDIM,NGLL2,1,igll,iface)] - 1;
    int k = free_surface_ijk[INDEX3(NDIM,NGLL2,2,igll,iface)] - 1;

    int iglob = ibool[INDEX4(5,5,5,i,j,k,ispec)] - 1 ;

    realw eta = ( noise_surface_movie[INDEX3(NDIM,NGLL2,0,igll,iface)]*normal_x_noise[ipoin]+
                 noise_surface_movie[INDEX3(NDIM,NGLL2,1,igll,iface)]*normal_y_noise[ipoin]+
                 noise_surface_movie[INDEX3(NDIM,NGLL2,2,igll,iface)]*normal_z_noise[ipoin]);

    Sigma_kl[INDEX4(5,5,5,i,j,k,ispec)] += deltat*eta*(normal_x_noise[ipoin]*displ[3*iglob]+
                                                       normal_y_noise[ipoin]*displ[1+3*iglob]+
                                                       normal_z_noise[ipoin]*displ[2+3*iglob]);
  }

}

/* ----------------------------------------------------------------------------------------------- */

extern "C"
void FC_FUNC_(compute_kernels_strgth_noise_cu,
              COMPUTE_KERNELS_STRGTH_NOISE_CU)(long* Mesh_pointer,
                                                    realw* h_noise_surface_movie,
                                                    realw* deltat) {

TRACE("compute_kernels_strgth_noise_cu");

  Mesh* mp = (Mesh*)(*Mesh_pointer); //get mesh pointer out of fortran integer container

  // checks if anything to do
  if( mp->num_free_surface_faces == 0 ) return;

  int num_blocks_x, num_blocks_y;
  get_blocks_xy(mp->num_free_surface_faces,&num_blocks_x,&num_blocks_y);

  dim3 grid(num_blocks_x,num_blocks_y);
  dim3 threads(NGLL2,1,1);

  print_CUDA_error_if_any(cudaMemcpy(mp->d_noise_surface_movie,h_noise_surface_movie,
                          NDIM*NGLL2*(mp->num_free_surface_faces)*sizeof(realw),cudaMemcpyHostToDevice),81000);

  compute_kernels_strength_noise_cuda_kernel<<<grid,threads>>>(mp->d_displ,
                                                               mp->d_free_surface_ispec,
                                                               mp->d_free_surface_ijk,
                                                               mp->d_ibool,
                                                               mp->d_noise_surface_movie,
                                                               mp->d_normal_x_noise,
                                                               mp->d_normal_y_noise,
                                                               mp->d_normal_z_noise,
                                                               mp->d_Sigma_kl,*deltat,
                                                               mp->num_free_surface_faces);

#ifdef ENABLE_VERY_SLOW_ERROR_CHECKING
  exit_on_cuda_error("compute_kernels_strength_noise_cuda_kernel");
#endif
}



/* ----------------------------------------------------------------------------------------------- */

// ACOUSTIC SIMULATIONS

/* ----------------------------------------------------------------------------------------------- */


__device__ void compute_gradient_kernel(int ijk,
                                        int ispec,
                                        realw* scalar_field,
                                        realw* vector_field_element,
                                        realw* d_hprime_xx,
                                        realw* d_xix,realw* d_xiy,realw* d_xiz,
                                        realw* d_etax,realw* d_etay,realw* d_etaz,
                                        realw* d_gammax,realw* d_gammay,realw* d_gammaz,
                                        realw rhol,
                                        int gravity) {

  realw temp1l,temp2l,temp3l;
  realw hp1,hp2,hp3;
  realw xixl,xiyl,xizl,etaxl,etayl,etazl,gammaxl,gammayl,gammazl;
  realw rho_invl;
  int l,offset,offset1,offset2,offset3;

  const int NGLL3_ALIGN = NGLL3_PADDED;

  int K = (ijk/NGLL2);
  int J = ((ijk-K*NGLL2)/NGLLX);
  int I = (ijk-K*NGLL2-J*NGLLX);

  // derivative along x
  temp1l = 0.f;
  for( l=0; l<NGLLX;l++){
    hp1 = d_hprime_xx[l*NGLLX+I];
    offset1 = K*NGLL2+J*NGLLX+l;
    temp1l += scalar_field[offset1]*hp1;
  }

  // derivative along y
  temp2l = 0.f;
  for( l=0; l<NGLLX;l++){
    // assumes hprime_xx == hprime_yy
    hp2 = d_hprime_xx[l*NGLLX+J];
    offset2 = K*NGLL2+l*NGLLX+I;
    temp2l += scalar_field[offset2]*hp2;
  }

  // derivative along z
  temp3l = 0.f;
  for( l=0; l<NGLLX;l++){
    // assumes hprime_xx == hprime_zz
    hp3 = d_hprime_xx[l*NGLLX+K];
    offset3 = l*NGLL2+J*NGLLX+I;
    temp3l += scalar_field[offset3]*hp3;
  }

  offset = ispec*NGLL3_ALIGN + ijk;

  xixl = d_xix[offset];
  xiyl = d_xiy[offset];
  xizl = d_xiz[offset];
  etaxl = d_etax[offset];
  etayl = d_etay[offset];
  etazl = d_etaz[offset];
  gammaxl = d_gammax[offset];
  gammayl = d_gammay[offset];
  gammazl = d_gammaz[offset];

  if( gravity ){
    // daniel: TODO - check gravity case here
    rho_invl = 1.0f / rhol;
  }else{
    rho_invl = 1.0f / rhol;
  }
  // derivatives of acoustic scalar potential field on GLL points
  vector_field_element[0] = (temp1l*xixl + temp2l*etaxl + temp3l*gammaxl) * rho_invl;
  vector_field_element[1] = (temp1l*xiyl + temp2l*etayl + temp3l*gammayl) * rho_invl;
  vector_field_element[2] = (temp1l*xizl + temp2l*etazl + temp3l*gammazl) * rho_invl;

}

/* ----------------------------------------------------------------------------------------------- */


__global__ void compute_kernels_acoustic_kernel(int* ispec_is_acoustic,
                                                int* ibool,
                                                realw* rhostore,
                                                realw* kappastore,
                                                realw* d_hprime_xx,
                                                realw* d_xix,realw* d_xiy,realw* d_xiz,
                                                realw* d_etax,realw* d_etay,realw* d_etaz,
                                                realw* d_gammax,realw* d_gammay,realw* d_gammaz,
                                                realw* potential_dot_dot_acoustic,
                                                realw* b_potential_acoustic,
                                                realw* b_potential_dot_dot_acoustic,
                                                realw* rho_ac_kl,
                                                realw* kappa_ac_kl,
                                                realw deltat,
                                                int NSPEC_AB,
                                                int gravity) {

  int ispec = blockIdx.x + blockIdx.y*gridDim.x;
  int ijk = threadIdx.x;

  // local and global indices
  int ijk_ispec = ijk + NGLL3*ispec;
  int ijk_ispec_padded = ijk + NGLL3_PADDED*ispec;
  int iglob;

  // shared memory between all threads within this block
  __shared__ realw scalar_field_displ[NGLL3];
  __shared__ realw scalar_field_accel[NGLL3];

  int active = 0;

  // handles case when there is 1 extra block (due to rectangular grid)
  if( ispec < NSPEC_AB ){
    // acoustic elements only
    if( ispec_is_acoustic[ispec] ){
      active = 1;

      // copy field values
      iglob = ibool[ijk_ispec] - 1;
      scalar_field_displ[ijk] = b_potential_acoustic[iglob];
      scalar_field_accel[ijk] = potential_dot_dot_acoustic[iglob];
    }
  }

  // synchronizes threads
  __syncthreads();

  if( active ){
    realw accel_elm[3];
    realw b_displ_elm[3];
    realw rhol,kappal;

    // gets material parameter
    rhol = rhostore[ijk_ispec_padded];

    // displacement vector from backward field
    compute_gradient_kernel(ijk,ispec,scalar_field_displ,b_displ_elm,
                            d_hprime_xx,
                            d_xix,d_xiy,d_xiz,d_etax,d_etay,d_etaz,d_gammax,d_gammay,d_gammaz,
                            rhol,gravity);

    // acceleration vector
    compute_gradient_kernel(ijk,ispec,scalar_field_accel,accel_elm,
                            d_hprime_xx,
                            d_xix,d_xiy,d_xiz,d_etax,d_etay,d_etaz,d_gammax,d_gammay,d_gammaz,
                            rhol,gravity);

    // density kernel
    rho_ac_kl[ijk_ispec] -= deltat * rhol * (accel_elm[0]*b_displ_elm[0] +
                                             accel_elm[1]*b_displ_elm[1] +
                                             accel_elm[2]*b_displ_elm[2]);

    // bulk modulus kernel
    kappal = kappastore[ijk_ispec];
    kappa_ac_kl[ijk_ispec] -= deltat / kappal * potential_dot_dot_acoustic[iglob]
                                              * b_potential_dot_dot_acoustic[iglob];
  } // active
}

/* ----------------------------------------------------------------------------------------------- */


extern "C"
void FC_FUNC_(compute_kernels_acoustic_cuda,
              COMPUTE_KERNELS_ACOUSTIC_CUDA)(long* Mesh_pointer,
                                             realw* deltat_f) {

TRACE("compute_kernels_acoustic_cuda");

  Mesh* mp = (Mesh*)(*Mesh_pointer); //get mesh pointer out of fortran integer container

  int blocksize = NGLL3; // NGLLX*NGLLY*NGLLZ
  realw deltat = *deltat_f;

  int num_blocks_x, num_blocks_y;
  get_blocks_xy(mp->NSPEC_AB,&num_blocks_x,&num_blocks_y);

  dim3 grid(num_blocks_x,num_blocks_y);
  dim3 threads(blocksize,1,1);

  compute_kernels_acoustic_kernel<<<grid,threads>>>(mp->d_ispec_is_acoustic,
                                                    mp->d_ibool,
                                                    mp->d_rhostore,
                                                    mp->d_kappastore,
                                                    mp->d_hprime_xx,
                                                    mp->d_xix,mp->d_xiy,mp->d_xiz,
                                                    mp->d_etax,mp->d_etay,mp->d_etaz,
                                                    mp->d_gammax,mp->d_gammay,mp->d_gammaz,
                                                    mp->d_potential_dot_dot_acoustic,
                                                    mp->d_b_potential_acoustic,
                                                    mp->d_b_potential_dot_dot_acoustic,
                                                    mp->d_rho_ac_kl,
                                                    mp->d_kappa_ac_kl,
                                                    deltat,
                                                    mp->NSPEC_AB,
                                                    mp->gravity);

#ifdef ENABLE_VERY_SLOW_ERROR_CHECKING
  exit_on_cuda_error("compute_kernels_acoustic_kernel");
#endif
}

/* ----------------------------------------------------------------------------------------------- */

// preconditioner (approximate Hessian kernel)

/* ----------------------------------------------------------------------------------------------- */

__global__ void compute_kernels_hess_el_cudakernel(int* ispec_is_elastic,
                                                   int* ibool,
                                                   realw* accel,
                                                   realw* b_accel,
                                                   realw* hess_kl,
                                                   realw deltat,
                                                   int NSPEC_AB) {

  int ispec = blockIdx.x + blockIdx.y*gridDim.x;
  int ijk = threadIdx.x;
  int ijk_ispec = ijk + NGLL3*ispec;

  // handles case when there is 1 extra block (due to rectangular grid)
  if(ispec < NSPEC_AB) {

    // elastic elements only
    if( ispec_is_elastic[ispec] ) {
      int iglob = ibool[ijk_ispec] - 1 ;

      // approximate hessian
      hess_kl[ijk_ispec] += deltat * (accel[3*iglob]*b_accel[3*iglob]+
                                      accel[3*iglob+1]*b_accel[3*iglob+1]+
                                      accel[3*iglob+2]*b_accel[3*iglob+2]);
    }
  }
}

/* ----------------------------------------------------------------------------------------------- */

__global__ void compute_kernels_hess_ac_cudakernel(int* ispec_is_acoustic,
                                                   int* ibool,
                                                   realw* potential_dot_dot_acoustic,
                                                   realw* b_potential_dot_dot_acoustic,
                                                   realw* rhostore,
                                                   realw* d_hprime_xx,
                                                   realw* d_xix,realw* d_xiy,realw* d_xiz,
                                                   realw* d_etax,realw* d_etay,realw* d_etaz,
                                                   realw* d_gammax,realw* d_gammay,realw* d_gammaz,
                                                   realw* hess_kl,
                                                   realw deltat,
                                                   int NSPEC_AB,
                                                   int gravity) {

  int ispec = blockIdx.x + blockIdx.y*gridDim.x;
  int ijk = threadIdx.x;
  int ijk_ispec = ijk + NGLL3*ispec;
  int ijk_ispec_padded = ijk + NGLL3_PADDED*ispec;
  int iglob;

  // shared memory between all threads within this block
  __shared__ realw scalar_field_accel[NGLL3];
  __shared__ realw scalar_field_b_accel[NGLL3];

  int active = 0;

  // handles case when there is 1 extra block (due to rectangular grid)
  if(ispec < NSPEC_AB) {

    // acoustic elements only
    if( ispec_is_acoustic[ispec] ){
      active = 1;

      // global indices
      iglob = ibool[ijk_ispec] - 1 ;

      // copy field values
      scalar_field_accel[ijk] = potential_dot_dot_acoustic[iglob];
      scalar_field_b_accel[ijk] = b_potential_dot_dot_acoustic[iglob];
    }
  }

  // synchronizes threads
  __syncthreads();

  if( active ){
    realw accel_elm[3];
    realw b_accel_elm[3];
    realw rhol;

    // gets material parameter
    rhol = rhostore[ijk_ispec_padded];

    // acceleration vector
    compute_gradient_kernel(ijk,ispec,
                            scalar_field_accel,accel_elm,
                            d_hprime_xx,
                            d_xix,d_xiy,d_xiz,d_etax,d_etay,d_etaz,d_gammax,d_gammay,d_gammaz,
                            rhol,gravity);

    // acceleration vector from backward field
    compute_gradient_kernel(ijk,ispec,
                            scalar_field_b_accel,b_accel_elm,
                            d_hprime_xx,
                            d_xix,d_xiy,d_xiz,d_etax,d_etay,d_etaz,d_gammax,d_gammay,d_gammaz,
                            rhol,gravity);
    // approximates hessian
    hess_kl[ijk_ispec] += deltat * (accel_elm[0]*b_accel_elm[0] +
                                    accel_elm[1]*b_accel_elm[1] +
                                    accel_elm[2]*b_accel_elm[2]);

  } // active
}

/* ----------------------------------------------------------------------------------------------- */

extern "C"
void FC_FUNC_(compute_kernels_hess_cuda,
              COMPUTE_KERNELS_HESS_CUDA)(long* Mesh_pointer,
                                         realw* deltat_f,
                                         int* ELASTIC_SIMULATION,
                                         int* ACOUSTIC_SIMULATION) {
  TRACE("compute_kernels_hess_cuda");

  Mesh* mp = (Mesh*)(*Mesh_pointer); //get mesh pointer out of fortran integer container

  int blocksize = NGLL3; // NGLLX*NGLLY*NGLLZ
  realw deltat = *deltat_f;

  int num_blocks_x, num_blocks_y;
  get_blocks_xy(mp->NSPEC_AB,&num_blocks_x,&num_blocks_y);

  dim3 grid(num_blocks_x,num_blocks_y);
  dim3 threads(blocksize,1,1);

  if( *ELASTIC_SIMULATION ) {
    compute_kernels_hess_el_cudakernel<<<grid,threads>>>(mp->d_ispec_is_elastic,
                                                         mp->d_ibool,
                                                         mp->d_accel,
                                                         mp->d_b_accel,
                                                         mp->d_hess_el_kl,
                                                         deltat,
                                                         mp->NSPEC_AB);
  }

  if( *ACOUSTIC_SIMULATION ) {
    compute_kernels_hess_ac_cudakernel<<<grid,threads>>>(mp->d_ispec_is_acoustic,
                                                         mp->d_ibool,
                                                         mp->d_potential_dot_dot_acoustic,
                                                         mp->d_b_potential_dot_dot_acoustic,
                                                         mp->d_rhostore,
                                                         mp->d_hprime_xx,
                                                         mp->d_xix,mp->d_xiy,mp->d_xiz,
                                                         mp->d_etax,mp->d_etay,mp->d_etaz,
                                                         mp->d_gammax,mp->d_gammay,mp->d_gammaz,
                                                         mp->d_hess_ac_kl,
                                                         deltat,
                                                         mp->NSPEC_AB,
                                                         mp->gravity);
  }


#ifdef ENABLE_VERY_SLOW_ERROR_CHECKING
  exit_on_cuda_error("compute_kernels_hess_cuda");
#endif
}

