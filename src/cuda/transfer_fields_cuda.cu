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
//#include "prepare_constants_cuda.h"

/* ----------------------------------------------------------------------------------------------- */

// Transfer functions

/* ----------------------------------------------------------------------------------------------- */



/* ----------------------------------------------------------------------------------------------- */

// for ELASTIC simulations

/* ----------------------------------------------------------------------------------------------- */

extern "C"
void FC_FUNC_(transfer_fields_el_to_device,
              TRANSFER_FIELDS_EL_TO_DEVICE)(int* size, realw* displ, realw* veloc, realw* accel,long* Mesh_pointer) {

  TRACE("transfer_fields_el_to_device");

  Mesh* mp = (Mesh*)(*Mesh_pointer); //get mesh pointer out of fortran integer container

  print_CUDA_error_if_any(cudaMemcpy(mp->d_displ,displ,sizeof(realw)*(*size),cudaMemcpyHostToDevice),40003);
  print_CUDA_error_if_any(cudaMemcpy(mp->d_veloc,veloc,sizeof(realw)*(*size),cudaMemcpyHostToDevice),40004);
  print_CUDA_error_if_any(cudaMemcpy(mp->d_accel,accel,sizeof(realw)*(*size),cudaMemcpyHostToDevice),40005);

}

/* ----------------------------------------------------------------------------------------------- */

extern "C"
void FC_FUNC_(transfer_fields_el_from_device,
              TRANSFER_FIELDS_EL_FROM_DEVICE)(int* size, realw* displ, realw* veloc, realw* accel,long* Mesh_pointer) {

  TRACE("transfer_fields_el_from_device");

  Mesh* mp = (Mesh*)(*Mesh_pointer); //get mesh pointer out of fortran integer container

  print_CUDA_error_if_any(cudaMemcpy(displ,mp->d_displ,sizeof(realw)*(*size),cudaMemcpyDeviceToHost),40006);
  print_CUDA_error_if_any(cudaMemcpy(veloc,mp->d_veloc,sizeof(realw)*(*size),cudaMemcpyDeviceToHost),40007);
  print_CUDA_error_if_any(cudaMemcpy(accel,mp->d_accel,sizeof(realw)*(*size),cudaMemcpyDeviceToHost),40008);

}

/* ----------------------------------------------------------------------------------------------- */

extern "C"
void FC_FUNC_(transfer_b_fields_to_device,
              TRANSFER_B_FIELDS_TO_DEVICE)(int* size, realw* b_displ, realw* b_veloc, realw* b_accel,
                                           long* Mesh_pointer) {

  TRACE("transfer_b_fields_to_device");

  Mesh* mp = (Mesh*)(*Mesh_pointer); //get mesh pointer out of fortran integer container

  print_CUDA_error_if_any(cudaMemcpy(mp->d_b_displ,b_displ,sizeof(realw)*(*size),cudaMemcpyHostToDevice),41006);
  print_CUDA_error_if_any(cudaMemcpy(mp->d_b_veloc,b_veloc,sizeof(realw)*(*size),cudaMemcpyHostToDevice),41007);
  print_CUDA_error_if_any(cudaMemcpy(mp->d_b_accel,b_accel,sizeof(realw)*(*size),cudaMemcpyHostToDevice),41008);

}

/* ----------------------------------------------------------------------------------------------- */

extern "C"
void FC_FUNC_(transfer_b_fields_from_device,
              TRANSFER_B_FIELDS_FROM_DEVICE)(int* size, realw* b_displ, realw* b_veloc, realw* b_accel,long* Mesh_pointer) {

  TRACE("transfer_b_fields_from_device");

  Mesh* mp = (Mesh*)(*Mesh_pointer); //get mesh pointer out of fortran integer container

  print_CUDA_error_if_any(cudaMemcpy(b_displ,mp->d_b_displ,sizeof(realw)*(*size),cudaMemcpyDeviceToHost),42006);
  print_CUDA_error_if_any(cudaMemcpy(b_veloc,mp->d_b_veloc,sizeof(realw)*(*size),cudaMemcpyDeviceToHost),42007);
  print_CUDA_error_if_any(cudaMemcpy(b_accel,mp->d_b_accel,sizeof(realw)*(*size),cudaMemcpyDeviceToHost),42008);

}


/* ----------------------------------------------------------------------------------------------- */

extern "C"
void FC_FUNC_(transfer_accel_to_device,
              TRNASFER_ACCEL_TO_DEVICE)(int* size, realw* accel,long* Mesh_pointer) {

  TRACE("transfer_accel_to_device");

  Mesh* mp = (Mesh*)(*Mesh_pointer); //get mesh pointer out of fortran integer container

  print_CUDA_error_if_any(cudaMemcpy(mp->d_accel,accel,sizeof(realw)*(*size),cudaMemcpyHostToDevice),40016);

}

/* ----------------------------------------------------------------------------------------------- */

extern "C"
void FC_FUNC_(transfer_accel_from_device,
              TRANSFER_ACCEL_FROM_DEVICE)(int* size, realw* accel,long* Mesh_pointer) {

  TRACE("transfer_accel_from_device");

  Mesh* mp = (Mesh*)(*Mesh_pointer); //get mesh pointer out of fortran integer container

  print_CUDA_error_if_any(cudaMemcpy(accel,mp->d_accel,sizeof(realw)*(*size),cudaMemcpyDeviceToHost),40026);

}

/* ----------------------------------------------------------------------------------------------- */

extern "C"
void FC_FUNC_(transfer_b_accel_from_device,
              TRNASFER_B_ACCEL_FROM_DEVICE)(int* size, realw* b_accel,long* Mesh_pointer) {

  TRACE("transfer_b_accel_from_device");

  Mesh* mp = (Mesh*)(*Mesh_pointer); //get mesh pointer out of fortran integer container

  print_CUDA_error_if_any(cudaMemcpy(b_accel,mp->d_b_accel,sizeof(realw)*(*size),cudaMemcpyDeviceToHost),40036);

}

/* ----------------------------------------------------------------------------------------------- */

extern "C"
void FC_FUNC_(transfer_sigma_from_device,
              TRANSFER_SIGMA_FROM_DEVICE)(int* size, realw* sigma_kl,long* Mesh_pointer) {

  TRACE("transfer_sigma_from_device");

  Mesh* mp = (Mesh*)(*Mesh_pointer); //get mesh pointer out of fortran integer container

  print_CUDA_error_if_any(cudaMemcpy(sigma_kl,mp->d_Sigma_kl,sizeof(realw)*(*size),cudaMemcpyDeviceToHost),40046);

}

/* ----------------------------------------------------------------------------------------------- */

extern "C"
void FC_FUNC_(transfer_b_displ_from_device,
              TRANSFER_B_DISPL_FROM_DEVICE)(int* size, realw* displ,long* Mesh_pointer) {

  TRACE("transfer_b_displ_from_device");

  Mesh* mp = (Mesh*)(*Mesh_pointer); //get mesh pointer out of fortran integer container

  print_CUDA_error_if_any(cudaMemcpy(displ,mp->d_b_displ,sizeof(realw)*(*size),cudaMemcpyDeviceToHost),40056);

}

/* ----------------------------------------------------------------------------------------------- */

extern "C"
void FC_FUNC_(transfer_displ_from_device,
              TRANSFER_DISPL_FROM_DEVICE)(int* size, realw* displ,long* Mesh_pointer) {

  TRACE("transfer_displ_from_device");

  Mesh* mp = (Mesh*)(*Mesh_pointer); //get mesh pointer out of fortran integer container

  print_CUDA_error_if_any(cudaMemcpy(displ,mp->d_displ,sizeof(realw)*(*size),cudaMemcpyDeviceToHost),40066);

}

/* ----------------------------------------------------------------------------------------------- */

// attenuation fields

extern "C"
void FC_FUNC_(transfer_b_fields_att_to_device,
              TRANSFER_B_FIELDS_ATT_TO_DEVICE)(long* Mesh_pointer,
                                               realw* b_R_xx,realw* b_R_yy,realw* b_R_xy,
                                               realw* b_R_xz,realw* b_R_yz,
                                               int* size_R,
                                               realw* b_epsilondev_xx,realw* b_epsilondev_yy,realw* b_epsilondev_xy,
                                               realw* b_epsilondev_xz,realw* b_epsilondev_yz,
                                               int* size_epsilondev) {

  TRACE("transfer_b_fields_att_to_device");

  //get mesh pointer out of fortran integer container
  Mesh* mp = (Mesh*)(*Mesh_pointer);

  print_CUDA_error_if_any(cudaMemcpy(mp->d_b_R_xx,b_R_xx,*size_R*sizeof(realw),cudaMemcpyHostToDevice),43011);
  print_CUDA_error_if_any(cudaMemcpy(mp->d_b_R_yy,b_R_yy,*size_R*sizeof(realw),cudaMemcpyHostToDevice),43012);
  print_CUDA_error_if_any(cudaMemcpy(mp->d_b_R_xy,b_R_xy,*size_R*sizeof(realw),cudaMemcpyHostToDevice),43013);
  print_CUDA_error_if_any(cudaMemcpy(mp->d_b_R_xz,b_R_xz,*size_R*sizeof(realw),cudaMemcpyHostToDevice),43014);
  print_CUDA_error_if_any(cudaMemcpy(mp->d_b_R_yz,b_R_yz,*size_R*sizeof(realw),cudaMemcpyHostToDevice),43015);

  print_CUDA_error_if_any(cudaMemcpy(mp->d_b_epsilondev_xx,b_epsilondev_xx,*size_epsilondev*sizeof(realw),cudaMemcpyHostToDevice),43016);
  print_CUDA_error_if_any(cudaMemcpy(mp->d_b_epsilondev_yy,b_epsilondev_yy,*size_epsilondev*sizeof(realw),cudaMemcpyHostToDevice),43017);
  print_CUDA_error_if_any(cudaMemcpy(mp->d_b_epsilondev_xy,b_epsilondev_xy,*size_epsilondev*sizeof(realw),cudaMemcpyHostToDevice),43018);
  print_CUDA_error_if_any(cudaMemcpy(mp->d_b_epsilondev_xz,b_epsilondev_xz,*size_epsilondev*sizeof(realw),cudaMemcpyHostToDevice),43019);
  print_CUDA_error_if_any(cudaMemcpy(mp->d_b_epsilondev_yz,b_epsilondev_yz,*size_epsilondev*sizeof(realw),cudaMemcpyHostToDevice),43020);

#ifdef ENABLE_VERY_SLOW_ERROR_CHECKING
  exit_on_cuda_error("after transfer_b_fields_att_to_device");
#endif
}

/* ----------------------------------------------------------------------------------------------- */

// attenuation fields

extern "C"
void FC_FUNC_(transfer_fields_att_from_device,
              TRANSFER_FIELDS_ATT_FROM_DEVICE)(long* Mesh_pointer,
                                               realw* R_xx,realw* R_yy,realw* R_xy,realw* R_xz,realw* R_yz,
                                               int* size_R,
                                               realw* epsilondev_xx,realw* epsilondev_yy,realw* epsilondev_xy,
                                               realw* epsilondev_xz,realw* epsilondev_yz,
                                               int* size_epsilondev) {
  TRACE("transfer_fields_att_from_device");

  //get mesh pointer out of fortran integer container
  Mesh* mp = (Mesh*)(*Mesh_pointer);

  print_CUDA_error_if_any(cudaMemcpy(R_xx,mp->d_R_xx,*size_R*sizeof(realw),cudaMemcpyDeviceToHost),43021);
  print_CUDA_error_if_any(cudaMemcpy(R_yy,mp->d_R_yy,*size_R*sizeof(realw),cudaMemcpyDeviceToHost),43022);
  print_CUDA_error_if_any(cudaMemcpy(R_xy,mp->d_R_xy,*size_R*sizeof(realw),cudaMemcpyDeviceToHost),43023);
  print_CUDA_error_if_any(cudaMemcpy(R_xz,mp->d_R_xz,*size_R*sizeof(realw),cudaMemcpyDeviceToHost),43024);
  print_CUDA_error_if_any(cudaMemcpy(R_yz,mp->d_R_yz,*size_R*sizeof(realw),cudaMemcpyDeviceToHost),43025);

  print_CUDA_error_if_any(cudaMemcpy(epsilondev_xx,mp->d_epsilondev_xx,*size_epsilondev*sizeof(realw),cudaMemcpyDeviceToHost),43026);
  print_CUDA_error_if_any(cudaMemcpy(epsilondev_yy,mp->d_epsilondev_yy,*size_epsilondev*sizeof(realw),cudaMemcpyDeviceToHost),43027);
  print_CUDA_error_if_any(cudaMemcpy(epsilondev_xy,mp->d_epsilondev_xy,*size_epsilondev*sizeof(realw),cudaMemcpyDeviceToHost),43028);
  print_CUDA_error_if_any(cudaMemcpy(epsilondev_xz,mp->d_epsilondev_xz,*size_epsilondev*sizeof(realw),cudaMemcpyDeviceToHost),43029);
  print_CUDA_error_if_any(cudaMemcpy(epsilondev_yz,mp->d_epsilondev_yz,*size_epsilondev*sizeof(realw),cudaMemcpyDeviceToHost),43030);


#ifdef ENABLE_VERY_SLOW_ERROR_CHECKING
  exit_on_cuda_error("after transfer_fields_att_from_device");
#endif
}

// JC JC here we will need to add GPU support for the new C-PML routines

/* ----------------------------------------------------------------------------------------------- */

extern "C"
void FC_FUNC_(transfer_kernels_el_to_host,
              TRANSFER_KERNELS_EL_TO_HOST)(long* Mesh_pointer,
                                            realw* h_rho_kl,
                                            realw* h_mu_kl,
                                            realw* h_kappa_kl,
                                            realw* h_cijkl_kl,
                                            int* NSPEC_AB) {
  TRACE("transfer_kernels_el_to_host");

  //get mesh pointer out of fortran integer container
  Mesh* mp = (Mesh*)(*Mesh_pointer);

  print_CUDA_error_if_any(cudaMemcpy(h_rho_kl,mp->d_rho_kl,*NSPEC_AB*NGLL3*sizeof(realw),
                                     cudaMemcpyDeviceToHost),40101);
  if( mp->anisotropic_kl ){
    print_CUDA_error_if_any(cudaMemcpy(h_cijkl_kl,mp->d_cijkl_kl,*NSPEC_AB*21*NGLL3*sizeof(realw),
                                       cudaMemcpyDeviceToHost),40102);
  }else{
    print_CUDA_error_if_any(cudaMemcpy(h_mu_kl,mp->d_mu_kl,*NSPEC_AB*NGLL3*sizeof(realw),
                                       cudaMemcpyDeviceToHost),40102);
    print_CUDA_error_if_any(cudaMemcpy(h_kappa_kl,mp->d_kappa_kl,*NSPEC_AB*NGLL3*sizeof(realw),
                                       cudaMemcpyDeviceToHost),40103);
  }
}

/* ----------------------------------------------------------------------------------------------- */

// for NOISE simulations

/* ----------------------------------------------------------------------------------------------- */

extern "C"
void FC_FUNC_(transfer_kernels_noise_to_host,
              TRANSFER_KERNELS_NOISE_TO_HOST)(long* Mesh_pointer,
                                              realw* h_Sigma_kl,
                                              int* NSPEC_AB) {
  TRACE("transfer_kernels_noise_to_host");

  //get mesh pointer out of fortran integer container
  Mesh* mp = (Mesh*)(*Mesh_pointer);

  print_CUDA_error_if_any(cudaMemcpy(h_Sigma_kl,mp->d_Sigma_kl,NGLL3*(*NSPEC_AB)*sizeof(realw),
                                     cudaMemcpyDeviceToHost),40201);

}


/* ----------------------------------------------------------------------------------------------- */

// for ACOUSTIC simulations

/* ----------------------------------------------------------------------------------------------- */

extern "C"
void FC_FUNC_(transfer_fields_ac_to_device,
              TRANSFER_FIELDS_AC_TO_DEVICE)(int* size,
                                            realw* potential_acoustic,
                                            realw* potential_dot_acoustic,
                                            realw* potential_dot_dot_acoustic,
                                            long* Mesh_pointer) {

  TRACE("transfer_fields_ac_to_device");

  //get mesh pointer out of fortran integer container
  Mesh* mp = (Mesh*)(*Mesh_pointer);

  print_CUDA_error_if_any(cudaMemcpy(mp->d_potential_acoustic,potential_acoustic,
                                     sizeof(realw)*(*size),cudaMemcpyHostToDevice),50110);
  print_CUDA_error_if_any(cudaMemcpy(mp->d_potential_dot_acoustic,potential_dot_acoustic,
                                     sizeof(realw)*(*size),cudaMemcpyHostToDevice),50120);
  print_CUDA_error_if_any(cudaMemcpy(mp->d_potential_dot_dot_acoustic,potential_dot_dot_acoustic,
                                     sizeof(realw)*(*size),cudaMemcpyHostToDevice),50130);

#ifdef ENABLE_VERY_SLOW_ERROR_CHECKING
  exit_on_cuda_error("after transfer_fields_ac_to_device");
#endif
}

/* ----------------------------------------------------------------------------------------------- */

extern "C"
void FC_FUNC_(transfer_b_fields_ac_to_device,
              TRANSFER_B_FIELDS_AC_TO_DEVICE)(int* size,
                                              realw* b_potential_acoustic,
                                              realw* b_potential_dot_acoustic,
                                              realw* b_potential_dot_dot_acoustic,
                                              long* Mesh_pointer) {

  TRACE("transfer_b_fields_ac_to_device");

  //get mesh pointer out of fortran integer container
  Mesh* mp = (Mesh*)(*Mesh_pointer);

  print_CUDA_error_if_any(cudaMemcpy(mp->d_b_potential_acoustic,b_potential_acoustic,
                                     sizeof(realw)*(*size),cudaMemcpyHostToDevice),51110);
  print_CUDA_error_if_any(cudaMemcpy(mp->d_b_potential_dot_acoustic,b_potential_dot_acoustic,
                                     sizeof(realw)*(*size),cudaMemcpyHostToDevice),51120);
  print_CUDA_error_if_any(cudaMemcpy(mp->d_b_potential_dot_dot_acoustic,b_potential_dot_dot_acoustic,
                                     sizeof(realw)*(*size),cudaMemcpyHostToDevice),51130);

#ifdef ENABLE_VERY_SLOW_ERROR_CHECKING
  exit_on_cuda_error("after transfer_b_fields_ac_to_device");
#endif
}


/* ----------------------------------------------------------------------------------------------- */

extern "C"
void FC_FUNC_(transfer_fields_ac_from_device,
              TRANSFER_FIELDS_AC_FROM_DEVICE)(int* size,
                                              realw* potential_acoustic,
                                              realw* potential_dot_acoustic,
                                              realw* potential_dot_dot_acoustic,
                                              long* Mesh_pointer) {
  TRACE("transfer_fields_ac_from_device");

  //get mesh pointer out of fortran integer container
  Mesh* mp = (Mesh*)(*Mesh_pointer);

  print_CUDA_error_if_any(cudaMemcpy(potential_acoustic,mp->d_potential_acoustic,
                                     sizeof(realw)*(*size),cudaMemcpyDeviceToHost),52111);
  print_CUDA_error_if_any(cudaMemcpy(potential_dot_acoustic,mp->d_potential_dot_acoustic,
                                     sizeof(realw)*(*size),cudaMemcpyDeviceToHost),52121);
  print_CUDA_error_if_any(cudaMemcpy(potential_dot_dot_acoustic,mp->d_potential_dot_dot_acoustic,
                                     sizeof(realw)*(*size),cudaMemcpyDeviceToHost),52131);

#ifdef ENABLE_VERY_SLOW_ERROR_CHECKING
  exit_on_cuda_error("after transfer_fields_ac_from_device");
#endif
}

/* ----------------------------------------------------------------------------------------------- */

extern "C"
void FC_FUNC_(transfer_b_fields_ac_from_device,
              TRANSFER_B_FIELDS_AC_FROM_DEVICE)(int* size,
                                                realw* b_potential_acoustic,
                                                realw* b_potential_dot_acoustic,
                                                realw* b_potential_dot_dot_acoustic,
                                                long* Mesh_pointer) {
  TRACE("transfer_b_fields_ac_from_device");

  //get mesh pointer out of fortran integer container
  Mesh* mp = (Mesh*)(*Mesh_pointer);

  print_CUDA_error_if_any(cudaMemcpy(b_potential_acoustic,mp->d_b_potential_acoustic,
                                     sizeof(realw)*(*size),cudaMemcpyDeviceToHost),53111);
  print_CUDA_error_if_any(cudaMemcpy(b_potential_dot_acoustic,mp->d_b_potential_dot_acoustic,
                                     sizeof(realw)*(*size),cudaMemcpyDeviceToHost),53121);
  print_CUDA_error_if_any(cudaMemcpy(b_potential_dot_dot_acoustic,mp->d_b_potential_dot_dot_acoustic,
                                     sizeof(realw)*(*size),cudaMemcpyDeviceToHost),53131);

#ifdef ENABLE_VERY_SLOW_ERROR_CHECKING
  exit_on_cuda_error("after transfer_b_fields_ac_from_device");
#endif
}

/* ----------------------------------------------------------------------------------------------- */

extern "C"
void FC_FUNC_(transfer_dot_dot_from_device,
              TRNASFER_DOT_DOT_FROM_DEVICE)(int* size, realw* potential_dot_dot_acoustic,long* Mesh_pointer) {

  TRACE("transfer_dot_dot_from_device");

  //get mesh pointer out of fortran integer container
  Mesh* mp = (Mesh*)(*Mesh_pointer);

  print_CUDA_error_if_any(cudaMemcpy(potential_dot_dot_acoustic,mp->d_potential_dot_dot_acoustic,
                                     sizeof(realw)*(*size),cudaMemcpyDeviceToHost),50041);

}

/* ----------------------------------------------------------------------------------------------- */

extern "C"
void FC_FUNC_(transfer_b_dot_dot_from_device,
              TRNASFER_B_DOT_DOT_FROM_DEVICE)(int* size, realw* b_potential_dot_dot_acoustic,long* Mesh_pointer) {

  TRACE("transfer_b_dot_dot_from_device");

  //get mesh pointer out of fortran integer container
  Mesh* mp = (Mesh*)(*Mesh_pointer);

  print_CUDA_error_if_any(cudaMemcpy(b_potential_dot_dot_acoustic,mp->d_b_potential_dot_dot_acoustic,
                                     sizeof(realw)*(*size),cudaMemcpyDeviceToHost),50042);

}


/* ----------------------------------------------------------------------------------------------- */

extern "C"
void FC_FUNC_(transfer_kernels_ac_to_host,
              TRANSFER_KERNELS_AC_TO_HOST)(long* Mesh_pointer,realw* h_rho_ac_kl,realw* h_kappa_ac_kl,int* NSPEC_AB) {

  TRACE("transfer_kernels_ac_to_host");

  //get mesh pointer out of fortran integer container
  Mesh* mp = (Mesh*)(*Mesh_pointer);

  int size = *NSPEC_AB*NGLL3;

  // copies kernel values over to CPU host
  print_CUDA_error_if_any(cudaMemcpy(h_rho_ac_kl,mp->d_rho_ac_kl,size*sizeof(realw),
                                     cudaMemcpyDeviceToHost),54101);
  print_CUDA_error_if_any(cudaMemcpy(h_kappa_ac_kl,mp->d_kappa_ac_kl,size*sizeof(realw),
                                     cudaMemcpyDeviceToHost),54102);
}

/* ----------------------------------------------------------------------------------------------- */

// for Hess kernel calculations

/* ----------------------------------------------------------------------------------------------- */

extern "C"
void FC_FUNC_(transfer_kernels_hess_el_tohost,
              TRANSFER_KERNELS_HESS_EL_TOHOST)(long* Mesh_pointer,realw* h_hess_kl,int* NSPEC_AB) {

  TRACE("transfer_kernels_hess_el_tohost");

  //get mesh pointer out of fortran integer container
  Mesh* mp = (Mesh*)(*Mesh_pointer);

  print_CUDA_error_if_any(cudaMemcpy(h_hess_kl,mp->d_hess_el_kl,NGLL3*(*NSPEC_AB)*sizeof(realw),
                                     cudaMemcpyDeviceToHost),70201);
}

/* ----------------------------------------------------------------------------------------------- */

extern "C"
void FC_FUNC_(transfer_kernels_hess_ac_tohost,
              TRANSFER_KERNELS_HESS_AC_TOHOST)(long* Mesh_pointer,realw* h_hess_ac_kl,int* NSPEC_AB) {

  TRACE("transfer_kernels_hess_ac_tohost");

  //get mesh pointer out of fortran integer container
  Mesh* mp = (Mesh*)(*Mesh_pointer);

  print_CUDA_error_if_any(cudaMemcpy(h_hess_ac_kl,mp->d_hess_ac_kl,NGLL3*(*NSPEC_AB)*sizeof(realw),
                                     cudaMemcpyDeviceToHost),70202);
}

// unused...

/* ----------------------------------------------------------------------------------------------- */
/*
extern "C"
void FC_FUNC_(transfer_compute_kernel_answers_from_device,
              TRANSFER_COMPUTE_KERNEL_ANSWERS_FROM_DEVICE)(long* Mesh_pointer,
                                                           realw* rho_kl,int* size_rho,
                                                           realw* mu_kl, int* size_mu,
                                                           realw* kappa_kl, int* size_kappa) {
TRACE("transfer_compute_kernel_answers_from_device");

  Mesh* mp = (Mesh*)(*Mesh_pointer); //get mesh pointer out of fortran integer container
  cudaMemcpy(rho_kl,mp->d_rho_kl,*size_rho*sizeof(realw),cudaMemcpyDeviceToHost);
  if( ! mp->anisotropic_kl ){
    cudaMemcpy(mu_kl,mp->d_mu_kl,*size_mu*sizeof(realw),cudaMemcpyDeviceToHost);
    cudaMemcpy(kappa_kl,mp->d_kappa_kl,*size_kappa*sizeof(realw),cudaMemcpyDeviceToHost);
  }
}
*/

/* ----------------------------------------------------------------------------------------------- */
/*
extern "C"
void FC_FUNC_(transfer_compute_kernel_fields_from_device,
              TRANSFER_COMPUTE_KERNEL_FIELDS_FROM_DEVICE)(long* Mesh_pointer,
                                                          realw* accel, int* size_accel,
                                                          realw* b_displ, int* size_b_displ,
                                                          realw* epsilondev_xx,realw* epsilondev_yy,realw* epsilondev_xy,
                                                          realw* epsilondev_xz,realw* epsilondev_yz,
                                                          int* size_epsilondev,
                                                          realw* b_epsilondev_xx,realw* b_epsilondev_yy,realw* b_epsilondev_xy,
                                                          realw* b_epsilondev_xz,realw* b_epsilondev_yz,
                                                          int* size_b_epsilondev,
                                                          realw* rho_kl,int* size_rho,
                                                          realw* mu_kl, int* size_mu,
                                                          realw* kappa_kl, int* size_kappa,
                                                          realw* epsilon_trace_over_3,
                                                          realw* b_epsilon_trace_over_3,
                                                          int* size_epsilon_trace_over_3) {
TRACE("transfer_compute_kernel_fields_from_device");

  Mesh* mp = (Mesh*)(*Mesh_pointer); //get mesh pointer out of fortran integer container
  cudaMemcpy(accel,mp->d_accel,*size_accel*sizeof(realw),cudaMemcpyDeviceToHost);
  cudaMemcpy(b_displ,mp->d_b_displ,*size_b_displ*sizeof(realw),cudaMemcpyDeviceToHost);
  cudaMemcpy(epsilondev_xx,mp->d_epsilondev_xx,*size_epsilondev*sizeof(realw),cudaMemcpyDeviceToHost);
  cudaMemcpy(epsilondev_yy,mp->d_epsilondev_yy,*size_epsilondev*sizeof(realw),cudaMemcpyDeviceToHost);
  cudaMemcpy(epsilondev_xy,mp->d_epsilondev_xy,*size_epsilondev*sizeof(realw),cudaMemcpyDeviceToHost);
  cudaMemcpy(epsilondev_xz,mp->d_epsilondev_xz,*size_epsilondev*sizeof(realw),cudaMemcpyDeviceToHost);
  cudaMemcpy(epsilondev_yz,mp->d_epsilondev_yz,*size_epsilondev*sizeof(realw),cudaMemcpyDeviceToHost);
  cudaMemcpy(b_epsilondev_xx,mp->d_b_epsilondev_xx,*size_b_epsilondev*sizeof(realw),cudaMemcpyDeviceToHost);
  cudaMemcpy(b_epsilondev_yy,mp->d_b_epsilondev_yy,*size_b_epsilondev*sizeof(realw),cudaMemcpyDeviceToHost);
  cudaMemcpy(b_epsilondev_xy,mp->d_b_epsilondev_xy,*size_b_epsilondev*sizeof(realw),cudaMemcpyDeviceToHost);
  cudaMemcpy(b_epsilondev_xz,mp->d_b_epsilondev_xz,*size_b_epsilondev*sizeof(realw),cudaMemcpyDeviceToHost);
  cudaMemcpy(b_epsilondev_yz,mp->d_b_epsilondev_yz,*size_b_epsilondev*sizeof(realw),cudaMemcpyDeviceToHost);
  cudaMemcpy(rho_kl,mp->d_rho_kl,*size_rho*sizeof(realw),cudaMemcpyDeviceToHost);

  if( ! mp->anisotropic_kl ){
    cudaMemcpy(mu_kl,mp->d_mu_kl,*size_mu*sizeof(realw),cudaMemcpyDeviceToHost);
    cudaMemcpy(kappa_kl,mp->d_kappa_kl,*size_kappa*sizeof(realw),cudaMemcpyDeviceToHost);
  }

  cudaMemcpy(epsilon_trace_over_3,mp->d_epsilon_trace_over_3,*size_epsilon_trace_over_3*sizeof(realw),
       cudaMemcpyDeviceToHost);
  cudaMemcpy(b_epsilon_trace_over_3,mp->d_b_epsilon_trace_over_3,*size_epsilon_trace_over_3*sizeof(realw),
       cudaMemcpyDeviceToHost);

#ifdef ENABLE_VERY_SLOW_ERROR_CHECKING
  exit_on_cuda_error("after transfer_compute_kernel_fields_from_device");
#endif
}
*/

