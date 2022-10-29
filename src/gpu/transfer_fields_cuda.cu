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

/* ----------------------------------------------------------------------------------------------- */

// Transfer functions

/* ----------------------------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------------------------------- */

// for ELASTIC simulations

/* ----------------------------------------------------------------------------------------------- */

extern EXTERN_LANG
void FC_FUNC_(transfer_fields_el_to_device,
              TRANSFER_FIELDS_EL_TO_DEVICE)(int* size, realw* displ, realw* veloc, realw* accel,long* Mesh_pointer) {

  TRACE("transfer_fields_el_to_device");

  Mesh* mp = (Mesh*)(*Mesh_pointer); //get mesh pointer out of fortran integer container

  gpuMemcpy_todevice_realw(mp->d_displ,displ,(*size));
  gpuMemcpy_todevice_realw(mp->d_veloc,veloc,(*size));
  gpuMemcpy_todevice_realw(mp->d_accel,accel,(*size));

}

/* ----------------------------------------------------------------------------------------------- */

extern EXTERN_LANG
void FC_FUNC_(transfer_fields_el_from_device,
              TRANSFER_FIELDS_EL_FROM_DEVICE)(int* size, realw* displ, realw* veloc, realw* accel,long* Mesh_pointer) {

  TRACE("transfer_fields_el_from_device");

  Mesh* mp = (Mesh*)(*Mesh_pointer); //get mesh pointer out of fortran integer container

  gpuMemcpy_tohost_realw(displ,mp->d_displ,(*size));
  gpuMemcpy_tohost_realw(veloc,mp->d_veloc,(*size));
  gpuMemcpy_tohost_realw(accel,mp->d_accel,(*size));

}

/* ----------------------------------------------------------------------------------------------- */

extern EXTERN_LANG
void FC_FUNC_(transfer_b_fields_to_device,
              TRANSFER_B_FIELDS_TO_DEVICE)(int* size, realw* b_displ, realw* b_veloc, realw* b_accel,
                                           long* Mesh_pointer) {

  TRACE("transfer_b_fields_to_device");

  Mesh* mp = (Mesh*)(*Mesh_pointer); //get mesh pointer out of fortran integer container

  gpuMemcpy_todevice_realw(mp->d_b_displ,b_displ,(*size));
  gpuMemcpy_todevice_realw(mp->d_b_veloc,b_veloc,(*size));
  gpuMemcpy_todevice_realw(mp->d_b_accel,b_accel,(*size));

}

/* ----------------------------------------------------------------------------------------------- */

extern EXTERN_LANG
void FC_FUNC_(transfer_b_fields_from_device,
              TRANSFER_B_FIELDS_FROM_DEVICE)(int* size, realw* b_displ, realw* b_veloc, realw* b_accel,long* Mesh_pointer) {

  TRACE("transfer_b_fields_from_device");

  Mesh* mp = (Mesh*)(*Mesh_pointer); //get mesh pointer out of fortran integer container

  gpuMemcpy_tohost_realw(b_displ,mp->d_b_displ,(*size));
  gpuMemcpy_tohost_realw(b_veloc,mp->d_b_veloc,(*size));
  gpuMemcpy_tohost_realw(b_accel,mp->d_b_accel,(*size));

}

/* ----------------------------------------------------------------------------------------------- */

extern EXTERN_LANG
void FC_FUNC_(transfer_veloc_from_device,
              TRANSFER_VELOC_FROM_DEVICE)(int* size, realw* veloc, long* Mesh_pointer) {

  TRACE("transfer_veloc_from_device");

  Mesh* mp = (Mesh*)(*Mesh_pointer); //get mesh pointer out of fortran integer container

  gpuMemcpy_tohost_realw(veloc,mp->d_veloc,(*size));

}

/* ----------------------------------------------------------------------------------------------- */

extern EXTERN_LANG
void FC_FUNC_(transfer_veloc_to_device,
              TRANSFER_VELOC_TO_DEVICE)(int* size, realw* veloc, long* Mesh_pointer) {

  TRACE("transfer_veloc_to_device");

  Mesh* mp = (Mesh*)(*Mesh_pointer); //get mesh pointer out of fortran integer container

  gpuMemcpy_todevice_realw(mp->d_veloc,veloc,(*size));

}

/* ----------------------------------------------------------------------------------------------- */

extern EXTERN_LANG
void FC_FUNC_(transfer_b_veloc_from_device,
              TRNASFER_B_VELOC_FROM_DEVICE)(int* size, realw* b_veloc,long* Mesh_pointer) {

  TRACE("transfer_b_veloc_from_device");

  Mesh* mp = (Mesh*)(*Mesh_pointer); //get mesh pointer out of fortran integer container

  gpuMemcpy_tohost_realw(b_veloc,mp->d_b_veloc,(*size));

}

/* ----------------------------------------------------------------------------------------------- */

extern EXTERN_LANG
void FC_FUNC_(transfer_b_veloc_to_device,
              TRANSFER_B_VELOC_TO_DEVICE)(int* size, realw* b_veloc,long* Mesh_pointer) {

  TRACE("transfer_b_veloc_to_device");

  Mesh* mp = (Mesh*)(*Mesh_pointer); //get mesh pointer out of fortran integer container

  gpuMemcpy_todevice_realw(mp->d_b_veloc,b_veloc,(*size));

}

/* ----------------------------------------------------------------------------------------------- */

extern EXTERN_LANG
void FC_FUNC_(transfer_accel_to_device,
              TRNASFER_ACCEL_TO_DEVICE)(int* size, realw* accel,long* Mesh_pointer) {

  TRACE("transfer_accel_to_device");

  Mesh* mp = (Mesh*)(*Mesh_pointer); //get mesh pointer out of fortran integer container

  gpuMemcpy_todevice_realw(mp->d_accel,accel,(*size));

}

/* ----------------------------------------------------------------------------------------------- */

extern EXTERN_LANG
void FC_FUNC_(transfer_accel_from_device,
              TRANSFER_ACCEL_FROM_DEVICE)(int* size, realw* accel,long* Mesh_pointer) {

  TRACE("transfer_accel_from_device");

  Mesh* mp = (Mesh*)(*Mesh_pointer); //get mesh pointer out of fortran integer container

  gpuMemcpy_tohost_realw(accel,mp->d_accel,(*size));

}

/* ----------------------------------------------------------------------------------------------- */

extern EXTERN_LANG
void FC_FUNC_(transfer_b_accel_from_device,
              TRNASFER_B_ACCEL_FROM_DEVICE)(int* size, realw* b_accel,long* Mesh_pointer) {

  TRACE("transfer_b_accel_from_device");

  Mesh* mp = (Mesh*)(*Mesh_pointer); //get mesh pointer out of fortran integer container

  gpuMemcpy_tohost_realw(b_accel,mp->d_b_accel,(*size));

}

/* ----------------------------------------------------------------------------------------------- */

extern EXTERN_LANG
void FC_FUNC_(transfer_b_accel_to_device,
              TRANSFER_B_accel_to_DEVICE)(int* size, realw* b_accel,long* Mesh_pointer) {

  TRACE("transfer_b_accel_to_device");

  Mesh* mp = (Mesh*)(*Mesh_pointer); //get mesh pointer out of fortran integer container

  gpuMemcpy_todevice_realw(mp->d_b_accel,b_accel,(*size));

}

/* ----------------------------------------------------------------------------------------------- */

extern EXTERN_LANG
void FC_FUNC_(transfer_sigma_from_device,
              TRANSFER_SIGMA_FROM_DEVICE)(int* size, realw* sigma_kl,long* Mesh_pointer) {

  TRACE("transfer_sigma_from_device");

  Mesh* mp = (Mesh*)(*Mesh_pointer); //get mesh pointer out of fortran integer container

  gpuMemcpy_tohost_realw(sigma_kl,mp->d_sigma_kl,(*size));

}

/* ----------------------------------------------------------------------------------------------- */

extern EXTERN_LANG
void FC_FUNC_(transfer_b_displ_from_device,
              TRANSFER_B_DISPL_FROM_DEVICE)(int* size, realw* b_displ,long* Mesh_pointer) {

  TRACE("transfer_b_displ_from_device");

  Mesh* mp = (Mesh*)(*Mesh_pointer); //get mesh pointer out of fortran integer container

  gpuMemcpy_tohost_realw(b_displ,mp->d_b_displ,(*size));

}

/* ----------------------------------------------------------------------------------------------- */

extern EXTERN_LANG
void FC_FUNC_(transfer_b_displ_to_device,
              TRANSFER_B_DISPL_to_DEVICE)(int* size, realw* b_displ,long* Mesh_pointer) {

  TRACE("transfer_b_displ_to_device");

  Mesh* mp = (Mesh*)(*Mesh_pointer); //get mesh pointer out of fortran integer container

  gpuMemcpy_todevice_realw(mp->d_b_displ,b_displ,(*size));

}

/* ----------------------------------------------------------------------------------------------- */

extern EXTERN_LANG
void FC_FUNC_(transfer_displ_from_device,
              TRANSFER_DISPL_FROM_DEVICE)(int* size, realw* displ,long* Mesh_pointer) {

  TRACE("transfer_displ_from_device");

  Mesh* mp = (Mesh*)(*Mesh_pointer); //get mesh pointer out of fortran integer container

  gpuMemcpy_tohost_realw(displ,mp->d_displ,(*size));

}

/* ----------------------------------------------------------------------------------------------- */

extern EXTERN_LANG
void FC_FUNC_(transfer_displ_to_device,
              TRANSFER_DISPL_TO_DEVICE)(int* size, realw* displ, long* Mesh_pointer) {

  TRACE("transfer_displ_to_device");

  Mesh* mp = (Mesh*)(*Mesh_pointer); //get mesh pointer out of fortran integer container

  gpuMemcpy_todevice_realw(mp->d_displ,displ,(*size));

}

/* ----------------------------------------------------------------------------------------------- */

// attenuation fields

/* ----------------------------------------------------------------------------------------------- */

// backward/reconstructed wavefields

extern EXTERN_LANG
void FC_FUNC_(transfer_b_rmemory_to_device,
              TRANSFER_B_RMEMORY_TO_DEVICE)(long* Mesh_pointer,
                                            realw* b_R_xx,realw* b_R_yy,realw* b_R_xy,
                                            realw* b_R_xz,realw* b_R_yz,
                                            realw* b_R_trace,
                                            int* size_R) {

  TRACE("transfer_b_rmemory_to_device");

  //get mesh pointer out of fortran integer container
  Mesh* mp = (Mesh*)(*Mesh_pointer);

  gpuMemcpy_todevice_realw(mp->d_b_R_xx,b_R_xx,*size_R);
  gpuMemcpy_todevice_realw(mp->d_b_R_yy,b_R_yy,*size_R);
  gpuMemcpy_todevice_realw(mp->d_b_R_xy,b_R_xy,*size_R);
  gpuMemcpy_todevice_realw(mp->d_b_R_xz,b_R_xz,*size_R);
  gpuMemcpy_todevice_realw(mp->d_b_R_yz,b_R_yz,*size_R);
  gpuMemcpy_todevice_realw(mp->d_b_R_trace,b_R_trace,*size_R);

  GPU_ERROR_CHECKING("after transfer_b_rmemory_to_device");
}

/* ----------------------------------------------------------------------------------------------- */

extern EXTERN_LANG
void FC_FUNC_(transfer_b_strain_to_device,
              TRANSFER_B_strain_TO_DEVICE)(long* Mesh_pointer,
                                           realw* b_epsilondev_xx,realw* b_epsilondev_yy,realw* b_epsilondev_xy,
                                           realw* b_epsilondev_xz,realw* b_epsilondev_yz,
                                           realw* b_epsilondev_trace,
                                           int* size_epsilondev) {

  TRACE("transfer_b_strain_to_device");

  //get mesh pointer out of fortran integer container
  Mesh* mp = (Mesh*)(*Mesh_pointer);

  gpuMemcpy_todevice_realw(mp->d_b_epsilondev_xx,b_epsilondev_xx,*size_epsilondev);
  gpuMemcpy_todevice_realw(mp->d_b_epsilondev_yy,b_epsilondev_yy,*size_epsilondev);
  gpuMemcpy_todevice_realw(mp->d_b_epsilondev_xy,b_epsilondev_xy,*size_epsilondev);
  gpuMemcpy_todevice_realw(mp->d_b_epsilondev_xz,b_epsilondev_xz,*size_epsilondev);
  gpuMemcpy_todevice_realw(mp->d_b_epsilondev_yz,b_epsilondev_yz,*size_epsilondev);
  gpuMemcpy_todevice_realw(mp->d_b_epsilondev_trace,b_epsilondev_trace,*size_epsilondev);

  GPU_ERROR_CHECKING("after transfer_b_strain_to_device");
}

/* ----------------------------------------------------------------------------------------------- */

// forward/adjoint wavefields

extern EXTERN_LANG
void FC_FUNC_(transfer_rmemory_from_device,
              TRANSFER_RMEMORY_FROM_DEVICE)(long* Mesh_pointer,
                                            realw* R_xx,realw* R_yy,realw* R_xy,realw* R_xz,realw* R_yz,
                                            realw* R_trace,
                                            int* size_R) {
  TRACE("transfer_rmemory_from_device");

  //get mesh pointer out of fortran integer container
  Mesh* mp = (Mesh*)(*Mesh_pointer);

  gpuMemcpy_tohost_realw(R_xx,mp->d_R_xx,*size_R);
  gpuMemcpy_tohost_realw(R_yy,mp->d_R_yy,*size_R);
  gpuMemcpy_tohost_realw(R_xy,mp->d_R_xy,*size_R);
  gpuMemcpy_tohost_realw(R_xz,mp->d_R_xz,*size_R);
  gpuMemcpy_tohost_realw(R_yz,mp->d_R_yz,*size_R);
  gpuMemcpy_tohost_realw(R_trace,mp->d_R_trace,*size_R);

  GPU_ERROR_CHECKING("after transfer_rmemory_from_device");
}

/* ----------------------------------------------------------------------------------------------- */

extern EXTERN_LANG
void FC_FUNC_(transfer_strain_from_device,
              TRANSFER_STRAIN_FROM_DEVICE)(long* Mesh_pointer,
                                           realw* epsilondev_xx,realw* epsilondev_yy,realw* epsilondev_xy,
                                           realw* epsilondev_xz,realw* epsilondev_yz,
                                           realw* epsilondev_trace,
                                           int* size_epsilondev) {
  TRACE("transfer_strain_from_device");

  //get mesh pointer out of fortran integer container
  Mesh* mp = (Mesh*)(*Mesh_pointer);

  gpuMemcpy_tohost_realw(epsilondev_xx,mp->d_epsilondev_xx,*size_epsilondev);
  gpuMemcpy_tohost_realw(epsilondev_yy,mp->d_epsilondev_yy,*size_epsilondev);
  gpuMemcpy_tohost_realw(epsilondev_xy,mp->d_epsilondev_xy,*size_epsilondev);
  gpuMemcpy_tohost_realw(epsilondev_xz,mp->d_epsilondev_xz,*size_epsilondev);
  gpuMemcpy_tohost_realw(epsilondev_yz,mp->d_epsilondev_yz,*size_epsilondev);
  gpuMemcpy_tohost_realw(epsilondev_trace,mp->d_epsilondev_trace,*size_epsilondev);

  GPU_ERROR_CHECKING("after transfer_strain_from_device");
}

/* ----------------------------------------------------------------------------------------------- */

// kernels

/* ----------------------------------------------------------------------------------------------- */

extern EXTERN_LANG
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

  gpuMemcpy_tohost_realw(h_rho_kl,mp->d_rho_kl,*NSPEC_AB*NGLL3);

  if (mp->anisotropic_kl ){
    gpuMemcpy_tohost_realw(h_cijkl_kl,mp->d_cijkl_kl,*NSPEC_AB*21*NGLL3);
  }else{
    gpuMemcpy_tohost_realw(h_mu_kl,mp->d_mu_kl,*NSPEC_AB*NGLL3);
    gpuMemcpy_tohost_realw(h_kappa_kl,mp->d_kappa_kl,*NSPEC_AB*NGLL3);
  }
}

/* ----------------------------------------------------------------------------------------------- */

// for NOISE simulations

/* ----------------------------------------------------------------------------------------------- */

extern EXTERN_LANG
void FC_FUNC_(transfer_kernels_noise_to_host,
              TRANSFER_KERNELS_NOISE_TO_HOST)(long* Mesh_pointer,
                                              realw* h_sigma_kl,
                                              int* NSPEC_AB) {
  TRACE("transfer_kernels_noise_to_host");

  //get mesh pointer out of fortran integer container
  Mesh* mp = (Mesh*)(*Mesh_pointer);

  gpuMemcpy_tohost_realw(h_sigma_kl,mp->d_sigma_kl,NGLL3*(*NSPEC_AB));
}


/* ----------------------------------------------------------------------------------------------- */

// for ACOUSTIC simulations

/* ----------------------------------------------------------------------------------------------- */

extern EXTERN_LANG
void FC_FUNC_(transfer_fields_ac_to_device,
              TRANSFER_FIELDS_AC_TO_DEVICE)(int* size,
                                            field* potential_acoustic,
                                            field* potential_dot_acoustic,
                                            field* potential_dot_dot_acoustic,
                                            long* Mesh_pointer) {

  TRACE("transfer_fields_ac_to_device");

  //get mesh pointer out of fortran integer container
  Mesh* mp = (Mesh*)(*Mesh_pointer);

  gpuMemcpy_todevice_field(mp->d_potential_acoustic,potential_acoustic,(*size));
  gpuMemcpy_todevice_field(mp->d_potential_dot_acoustic,potential_dot_acoustic,(*size));
  gpuMemcpy_todevice_field(mp->d_potential_dot_dot_acoustic,potential_dot_dot_acoustic,(*size));

  GPU_ERROR_CHECKING("after transfer_fields_ac_to_device");
}

/* ----------------------------------------------------------------------------------------------- */

extern EXTERN_LANG
void FC_FUNC_(transfer_b_fields_ac_to_device,
              TRANSFER_B_FIELDS_AC_TO_DEVICE)(int* size,
                                              field* b_potential_acoustic,
                                              field* b_potential_dot_acoustic,
                                              field* b_potential_dot_dot_acoustic,
                                              long* Mesh_pointer) {

  TRACE("transfer_b_fields_ac_to_device");

  //get mesh pointer out of fortran integer container
  Mesh* mp = (Mesh*)(*Mesh_pointer);

  gpuMemcpy_todevice_field(mp->d_b_potential_acoustic,b_potential_acoustic,(*size));
  gpuMemcpy_todevice_field(mp->d_b_potential_dot_acoustic,b_potential_dot_acoustic,(*size));
  gpuMemcpy_todevice_field(mp->d_b_potential_dot_dot_acoustic,b_potential_dot_dot_acoustic,(*size));

  GPU_ERROR_CHECKING("after transfer_b_fields_ac_to_device");
}


/* ----------------------------------------------------------------------------------------------- */

extern EXTERN_LANG
void FC_FUNC_(transfer_fields_ac_from_device,
              TRANSFER_FIELDS_AC_FROM_DEVICE)(int* size,
                                              field* potential_acoustic,
                                              field* potential_dot_acoustic,
                                              field* potential_dot_dot_acoustic,
                                              long* Mesh_pointer) {
  TRACE("transfer_fields_ac_from_device");

  //get mesh pointer out of fortran integer container
  Mesh* mp = (Mesh*)(*Mesh_pointer);

  gpuMemcpy_tohost_field(potential_acoustic,mp->d_potential_acoustic,(*size));
  gpuMemcpy_tohost_field(potential_dot_acoustic,mp->d_potential_dot_acoustic,(*size));
  gpuMemcpy_tohost_field(potential_dot_dot_acoustic,mp->d_potential_dot_dot_acoustic,(*size));

  GPU_ERROR_CHECKING("after transfer_fields_ac_from_device");
}

/* ----------------------------------------------------------------------------------------------- */

extern EXTERN_LANG
void FC_FUNC_(transfer_b_fields_ac_from_device,
              TRANSFER_B_FIELDS_AC_FROM_DEVICE)(int* size,
                                                field* b_potential_acoustic,
                                                field* b_potential_dot_acoustic,
                                                field* b_potential_dot_dot_acoustic,
                                                long* Mesh_pointer) {
  TRACE("transfer_b_fields_ac_from_device");

  //get mesh pointer out of fortran integer container
  Mesh* mp = (Mesh*)(*Mesh_pointer);

  gpuMemcpy_tohost_field(b_potential_acoustic,mp->d_b_potential_acoustic,(*size));
  gpuMemcpy_tohost_field(b_potential_dot_acoustic,mp->d_b_potential_dot_acoustic,(*size));
  gpuMemcpy_tohost_field(b_potential_dot_dot_acoustic,mp->d_b_potential_dot_dot_acoustic,(*size));

  GPU_ERROR_CHECKING("after transfer_b_fields_ac_from_device");
}

/* ----------------------------------------------------------------------------------------------- */

extern EXTERN_LANG
void FC_FUNC_(transfer_b_potential_ac_from_device,
              TRANSFER_B_potentical_AC_FROM_DEVICE)(int* size,
                                                    field* b_potential_acoustic,
                                                    long* Mesh_pointer) {
  TRACE("transfer_b_potential_ac_from_device");

  //get mesh pointer out of fortran integer container
  Mesh* mp = (Mesh*)(*Mesh_pointer);

  gpuMemcpy_tohost_field(b_potential_acoustic,mp->d_b_potential_acoustic,(*size));

  GPU_ERROR_CHECKING("after transfer_b_potential_ac_from_device");
}

/* ----------------------------------------------------------------------------------------------- */

extern EXTERN_LANG
void FC_FUNC_(transfer_b_potential_dot_dot_ac_from_device,
              TRANSFER_B_potentical_DOT_DOT_AC_FROM_DEVICE)(int* size,
                                                            field* b_potential_dot_dot_acoustic,
                                                            long* Mesh_pointer) {
  TRACE("transfer_b_potential_dot_dot_ac_from_device");

  //get mesh pointer out of fortran integer container
  Mesh* mp = (Mesh*)(*Mesh_pointer);

  gpuMemcpy_tohost_field(b_potential_dot_dot_acoustic,mp->d_b_potential_dot_dot_acoustic,(*size));

  GPU_ERROR_CHECKING("after transfer_b_potential_dot_dot_ac_from_device");
}

/* ----------------------------------------------------------------------------------------------- */

extern EXTERN_LANG
void FC_FUNC_(transfer_b_potential_ac_to_device,
              TRANSFER_B_potentical_AC_TO_DEVICE)(int* size,
                                                  field* b_potential_acoustic,
                                                  long* Mesh_pointer) {
  TRACE("transfer_b_potential_ac_to_device");

  //get mesh pointer out of fortran integer container
  Mesh* mp = (Mesh*)(*Mesh_pointer);

  gpuMemcpy_todevice_field(mp->d_b_potential_acoustic,b_potential_acoustic,(*size));

  GPU_ERROR_CHECKING("after transfer_b_potential_ac_to_device");
}

/* ----------------------------------------------------------------------------------------------- */

extern EXTERN_LANG
void FC_FUNC_(transfer_b_potential_dot_dot_ac_to_device,
              TRANSFER_B_potentical_DOT_DOT_AC_TO_DEVICE)(int* size,
                                                          field* b_potential_dot_dot_acoustic,
                                                          long* Mesh_pointer) {
  TRACE("transfer_b_potential_ac_to_device");

  //get mesh pointer out of fortran integer container
  Mesh* mp = (Mesh*)(*Mesh_pointer);

  gpuMemcpy_todevice_field(mp->d_b_potential_dot_dot_acoustic,b_potential_dot_dot_acoustic,(*size));

  GPU_ERROR_CHECKING("after transfer_b_potential_dot_dot_ac_to_device");
}

/* ----------------------------------------------------------------------------------------------- */

extern EXTERN_LANG
void FC_FUNC_(transfer_dot_dot_from_device,
              TRNASFER_DOT_DOT_FROM_DEVICE)(int* size, field* potential_dot_dot_acoustic,long* Mesh_pointer) {

  TRACE("transfer_dot_dot_from_device");

  //get mesh pointer out of fortran integer container
  Mesh* mp = (Mesh*)(*Mesh_pointer);

  gpuMemcpy_tohost_field(potential_dot_dot_acoustic,mp->d_potential_dot_dot_acoustic,(*size));

}

/* ----------------------------------------------------------------------------------------------- */

extern EXTERN_LANG
void FC_FUNC_(transfer_b_dot_dot_from_device,
              TRNASFER_B_DOT_DOT_FROM_DEVICE)(int* size, field* b_potential_dot_dot_acoustic,long* Mesh_pointer) {

  TRACE("transfer_b_dot_dot_from_device");

  //get mesh pointer out of fortran integer container
  Mesh* mp = (Mesh*)(*Mesh_pointer);

  gpuMemcpy_tohost_field(b_potential_dot_dot_acoustic,mp->d_b_potential_dot_dot_acoustic,(*size));

}


/* ----------------------------------------------------------------------------------------------- */

extern EXTERN_LANG
void FC_FUNC_(transfer_kernels_ac_to_host,
              TRANSFER_KERNELS_AC_TO_HOST)(long* Mesh_pointer,realw* h_rho_ac_kl,realw* h_kappa_ac_kl,int* NSPEC_AB) {

  TRACE("transfer_kernels_ac_to_host");

  //get mesh pointer out of fortran integer container
  Mesh* mp = (Mesh*)(*Mesh_pointer);

  int size = *NSPEC_AB*NGLL3;

  // copies kernel values over to CPU host
  gpuMemcpy_tohost_realw(h_rho_ac_kl,mp->d_rho_ac_kl,size);
  gpuMemcpy_tohost_realw(h_kappa_ac_kl,mp->d_kappa_ac_kl,size);
}

/* ----------------------------------------------------------------------------------------------- */

// for Hess kernel calculations

/* ----------------------------------------------------------------------------------------------- */

extern EXTERN_LANG
void FC_FUNC_(transfer_kernels_hess_el_tohost,
              TRANSFER_KERNELS_HESS_EL_TOHOST)(long* Mesh_pointer,
                 realw* h_hess_kl,
                 realw* h_hess_rho_kl,
                 realw* h_hess_kappa_kl,
                 realw* h_hess_mu_kl,
                 int* NSPEC_AB) {

  TRACE("transfer_kernels_hess_el_tohost");

  //get mesh pointer out of fortran integer container
  Mesh* mp = (Mesh*)(*Mesh_pointer);

  gpuMemcpy_tohost_realw(h_hess_kl,mp->d_hess_el_kl,NGLL3*(*NSPEC_AB));
  gpuMemcpy_tohost_realw(h_hess_rho_kl,mp->d_hess_rho_el_kl,NGLL3*(*NSPEC_AB));

  gpuMemcpy_tohost_realw(h_hess_kappa_kl,mp->d_hess_kappa_el_kl,NGLL3*(*NSPEC_AB));
  gpuMemcpy_tohost_realw(h_hess_mu_kl,mp->d_hess_mu_el_kl,NGLL3*(*NSPEC_AB));

  //printf("%e %e \n",h_hess_kappa_kl[125], h_hess_mu_kl[125]);
}

/* ----------------------------------------------------------------------------------------------- */

extern EXTERN_LANG
void FC_FUNC_(transfer_kernels_hess_ac_tohost,
              TRANSFER_KERNELS_HESS_AC_TOHOST)(long* Mesh_pointer,
                 realw* h_hess_ac_kl,
                 realw* h_hess_rho_ac_kl,
                 realw* h_hess_kappa_ac_kl,
                 int* NSPEC_AB) {

  TRACE("transfer_kernels_hess_ac_tohost");

  //get mesh pointer out of fortran integer container
  Mesh* mp = (Mesh*)(*Mesh_pointer);

  gpuMemcpy_tohost_realw(h_hess_ac_kl,mp->d_hess_ac_kl,NGLL3*(*NSPEC_AB));
  gpuMemcpy_tohost_realw(h_hess_rho_ac_kl,mp->d_hess_rho_ac_kl,NGLL3*(*NSPEC_AB));
  gpuMemcpy_tohost_realw(h_hess_kappa_ac_kl,mp->d_hess_kappa_ac_kl,NGLL3*(*NSPEC_AB));
}

// unused...

/* ----------------------------------------------------------------------------------------------- */
/*
extern EXTERN_LANG
void FC_FUNC_(transfer_compute_kernel_answers_from_device,
              TRANSFER_COMPUTE_KERNEL_ANSWERS_FROM_DEVICE)(long* Mesh_pointer,
                                                           realw* rho_kl,int* size_rho,
                                                           realw* mu_kl, int* size_mu,
                                                           realw* kappa_kl, int* size_kappa) {
TRACE("transfer_compute_kernel_answers_from_device");

  Mesh* mp = (Mesh*)(*Mesh_pointer); //get mesh pointer out of fortran integer container
  gpuMemcpy_tohost_realw(rho_kl,mp->d_rho_kl,*size_rho);
  if (! mp->anisotropic_kl ){
    gpuMemcpy_tohost_realw(mu_kl,mp->d_mu_kl,*size_mu);
    gpuMemcpy_tohost_realw(kappa_kl,mp->d_kappa_kl,*size_kappa);
  }
}
*/

/* ----------------------------------------------------------------------------------------------- */
/*
extern EXTERN_LANG
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
  gpuMemcpy_tohost_realw(accel,mp->d_accel,*size_accel);
  gpuMemcpy_tohost_realw(b_displ,mp->d_b_displ,*size_b_displ);
  gpuMemcpy_tohost_realw(epsilondev_xx,mp->d_epsilondev_xx,*size_epsilondev);
  gpuMemcpy_tohost_realw(epsilondev_yy,mp->d_epsilondev_yy,*size_epsilondev);
  gpuMemcpy_tohost_realw(epsilondev_xy,mp->d_epsilondev_xy,*size_epsilondev);
  gpuMemcpy_tohost_realw(epsilondev_xz,mp->d_epsilondev_xz,*size_epsilondev);
  gpuMemcpy_tohost_realw(epsilondev_yz,mp->d_epsilondev_yz,*size_epsilondev);
  gpuMemcpy_tohost_realw(b_epsilondev_xx,mp->d_b_epsilondev_xx,*size_b_epsilondev);
  gpuMemcpy_tohost_realw(b_epsilondev_yy,mp->d_b_epsilondev_yy,*size_b_epsilondev);
  gpuMemcpy_tohost_realw(b_epsilondev_xy,mp->d_b_epsilondev_xy,*size_b_epsilondev);
  gpuMemcpy_tohost_realw(b_epsilondev_xz,mp->d_b_epsilondev_xz,*size_b_epsilondev);
  gpuMemcpy_tohost_realw(b_epsilondev_yz,mp->d_b_epsilondev_yz,*size_b_epsilondev);
  gpuMemcpy_tohost_realw(rho_kl,mp->d_rho_kl,*size_rho);

  if (! mp->anisotropic_kl ){
    gpuMemcpy_tohost_realw(mu_kl,mp->d_mu_kl,*size_mu);
    gpuMemcpy_tohost_realw(kappa_kl,mp->d_kappa_kl,*size_kappa);
  }

  gpuMemcpy_tohost_realw(epsilon_trace_over_3,mp->d_epsilon_trace_over_3,*size_epsilon_trace_over_3);
  gpuMemcpy_tohost_realw(b_epsilon_trace_over_3,mp->d_b_epsilon_trace_over_3,*size_epsilon_trace_over_3);

  GPU_ERROR_CHECKING("after transfer_compute_kernel_fields_from_device");
}
*/

/* ----------------------------------------------------------------------------------------------- */
// register host array for pinned memory

extern EXTERN_LANG
void FC_FUNC_(register_host_array,
              REGISTER_HOST_ARRAY)(int *size, realw *h_array) {

  TRACE("register_host_array");

  // page-locks the memory to automatically accelerate calls to functions such as cudaMemcpy()
  // since the memory can be accessed directly by the device, it can be read or written with
  // much higher bandwidth than pageable memory that has not been registered.
  // Page-locking excessive amounts of memory may degrade system performance,
  // since it reduces the amount of memory available to the system for paging.
  // As a result, this function is best used sparingly to register staging areas for data exchange between host and device.

#ifdef USE_CUDA
  if (run_cuda){
    print_CUDA_error_if_any(cudaHostRegister(h_array, (*size)*sizeof(realw), 0),55001);
  }
#endif
#ifdef USE_HIP
  if (run_hip){
    print_HIP_error_if_any(hipHostRegister(h_array, (*size)*sizeof(realw), 0),55001);
  }
#endif

  GPU_ERROR_CHECKING ("after register_host_array");
}


extern EXTERN_LANG
void FC_FUNC_(unregister_host_array,
              UNREGISTER_HOST_ARRAY)(realw *h_array) {

  TRACE("unregister_host_array");

#ifdef USE_CUDA
  if (run_cuda){
    print_CUDA_error_if_any(cudaHostUnregister(h_array),55002);
  }
#endif
#ifdef USE_HIP
  if (run_hip){
    print_HIP_error_if_any(hipHostUnregister(h_array),55002);
  }
#endif

  GPU_ERROR_CHECKING ("after unregister_host_array");
}
