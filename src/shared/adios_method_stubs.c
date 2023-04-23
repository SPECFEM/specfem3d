/*
!=====================================================================
!
!                          S p e c f e m 3 D
!                          -----------------
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

#include <stdio.h>
#include <stdlib.h>

#include "config.h"

typedef float realw;


// for xmeshfem3D compilation

void FC_FUNC_(save_databases_adios,SAVE_DATABASES_ADIOS)(char* LOCAL_PATH, int* sizeprocs,
                                                         int* nspec, int* nglob,
                                                         int* iMPIcut_xi, int* iMPIcut_eta,
                                                         double* nodes_coords, int* ispec_material_id,
                                                         int* nspec2D_xmin, int* nspec2D_xmax,
                                                         int* nspec2D_ymin, int* nspec2D_ymax,
                                                         int* ibelm_xmin, int* ibelm_xmax,
                                                         int* ibelm_ymin, int* ibelm_ymax,
                                                         int* ibelm_bottom, int* ibelm_top) {}


// for xgenerate_databases compilation

// subroutines from model_gll_adios.F90
void FC_FUNC_(model_gll_adios,MODEL_GLL_ADIOS)(int* myrank, int* nspec, char* LOCAL_PATH) {}

// subroutines from model_ipati_adios.F90
void FC_FUNC_(model_ipati_adios,MODEL_IPATI_ADIOS)(int* myrank, int* nspec, char* LOCAL_PATH) {}
void FC_FUNC_(model_ipati_water_adios,MODEL_IPATI_WATER_ADIOS)(int* myrank, int* nspec, char* LOCAL_PATH) {}

// subroutines from read_partition_files_adios.F90
void FC_FUNC_(read_partition_files_adios,READ_PARTITION_FILES_ADIOS)(void) {}

// subroutines from save_arrays_solver_adios.F90
void FC_FUNC_(save_arrays_solver_mesh_adios,SAVE_ARRAYS_SOLVER_MESH_ADIOS)(void) {}

// subroutines from save_moho_adios.F90
void FC_FUNC_(crm_save_moho_adios,CRM_SAVE_MOHO_ADIOS)(void) {}


// for xspecfem3D compilation

void FC_FUNC_(read_mesh_for_init_adios,READ_MESH_FOR_INIT_ADIOS)(void) {}

void FC_FUNC_(read_mesh_databases_adios,READ_MESH_DATABASES_ADIOS)(void) {}

void FC_FUNC_(read_mesh_databases_moho_adios,READ_MESH_DATABASES_MOHO_ADIOS)(void) {}

void FC_FUNC_(define_kernel_adios_variables,DEFINE_KERNEL_ADIOS_VARIABLES)(void) {}

void FC_FUNC_(perform_write_adios_kernels,PERFORM_WRITE_ADIOS_KERNELS)(void) {}

void FC_FUNC_(read_forward_arrays_adios,READ_FORWARD_ARRAYS_ADIOS)(void) {}

void FC_FUNC_(read_forward_arrays_undoatt_adios,READ_FORWARD_ARRAYS_UNDOATT_ADIOS)(int* iteration_on_subset_tmp) {}

void FC_FUNC_(save_forward_arrays_adios,SAVE_FORWARD_ARRAYS_ADIOS)(void) {}

void FC_FUNC_(save_forward_arrays_undoatt_adios,SAVE_FORWARD_UNDOATT_ARRAYS_ADIOS)(void) {}

void FC_FUNC_(save_kernels_acoustic_adios,SAVE_KERNELS_ACOUSTIC_ADIOS)(void) {}

void FC_FUNC_(save_kernels_elastic_iso_adios,SAVE_KERNELS_ELASTIC_ISO_ADIOS)(realw* rhop_kl, realw* alpha_kl, realw* beta_kl) {}

void FC_FUNC_(save_kernels_elastic_aniso_adios,SAVE_KERNELS_ELASTIC_ANISO_ADIOS)(realw* alphav_kl, realw* alphah_kl,
                                                                                 realw* betav_kl, realw* betah_kl, realw* eta_kl,
                                                                                 realw* c11_kl,realw* c12_kl,realw* c13_kl,realw* c14_kl,realw* c15_kl,realw* c16_kl,
                                                                                 realw* c22_kl,realw* c23_kl,realw* c24_kl,realw* c25_kl,realw* c26_kl,
                                                                                 realw* c33_kl,realw* c34_kl,realw* c35_kl,realw* c36_kl,
                                                                                 realw* c44_kl,realw* c45_kl,realw* c46_kl,
                                                                                 realw* c55_kl,realw* c56_kl,
                                                                                 realw* c66_kl) {}

void FC_FUNC_(save_kernels_poroelastic_adios,SAVE_KERNELS_POROELASTIC_ADIOS)(void) {}

void FC_FUNC_(save_kernels_moho_adios,SAVE_KERNELS_MOHO_ADIOS)(void) {}

void FC_FUNC_(save_kernels_hessian_adios,SAVE_KERNELS_HESSIAN_ADIOS)(void) {}

