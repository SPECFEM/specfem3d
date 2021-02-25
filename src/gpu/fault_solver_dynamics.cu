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

#include "mesh_constants_cuda.h"
#include "fault_struct_cuda.h"

/* ----------------------------------------------------------------------------------------------- */

extern EXTERN_LANG
void FC_FUNC_(initialize_fault_solver_gpu,
              INITIALIZE_FAULT_SOLVER_GPU)(long* Fault_solver,
                                           int* num_of_faults,
                                           realw* v_healing,
                                           realw* v_rupt) {

  TRACE("initialize_fault_solver_gpu");

  // allocates fault parameter structure
  Fault_solver_dynamics *Fdyn = (Fault_solver_dynamics*) malloc(sizeof(Fault_solver_dynamics));
  if (Fdyn == NULL) exit_on_error("error allocating fault_solver pointer");

  *Fault_solver = (long) Fdyn;

  // initializes
  Fdyn->Nbfaults = *num_of_faults;
  Fdyn->faults = (Fault*) malloc((*num_of_faults)*sizeof(Fault));

  Fdyn->v_healing = *v_healing;
  Fdyn->v_rupt = *v_rupt;
}

/* ----------------------------------------------------------------------------------------------- */

extern EXTERN_LANG
void FC_FUNC_(initialize_fault_data_gpu,
              INITIALIZE_FAULT_DATA_GPU)(long* Fault_solver,
                                         int* iglob,
                                         int* num_of_records,
                                         int* nt,
                                         int* ifault) {
  TRACE("initialize_fault_data_gpu");

  Fault_solver_dynamics* Fsolver = (Fault_solver_dynamics*)(*Fault_solver);

  // fault data recorder
  Fault_data *fault_data_recorder = (Fault_data*) malloc(sizeof(Fault_data));
  fault_data_recorder->NRECORD = *num_of_records;
  fault_data_recorder->NT = *nt;

  int recordlength = 7; //store 7 different quantities

  // allocates arrays on GPU
  print_CUDA_error_if_any(cudaMalloc((void**)&(fault_data_recorder->iglob),
                                     fault_data_recorder->NRECORD*sizeof(int)),60001);

  print_CUDA_error_if_any(cudaMemcpy(fault_data_recorder->iglob, iglob,
                                     fault_data_recorder->NRECORD*sizeof(int), cudaMemcpyHostToDevice),60002);

  print_CUDA_error_if_any(cudaMalloc((void**)&(fault_data_recorder->dataT),
                                     recordlength*fault_data_recorder->NRECORD*fault_data_recorder->NT*sizeof(realw)),60003);

  // stores structure pointers
  Fsolver->faults[*ifault].output_dataT = fault_data_recorder;

  GPU_ERROR_CHECKING("initialize_fault_data_gpu");
}


/* ----------------------------------------------------------------------------------------------- */

// copies realw array from CPU host to GPU device
void copy_todevice_realw_test(void** d_array_addr_ptr,realw* h_array,int size) {

  // allocates memory on GPU
  cudaMalloc((void**)d_array_addr_ptr,size*sizeof(realw));
  // copies values onto GPU
  cudaMemcpy((realw*) *d_array_addr_ptr,h_array,size*sizeof(realw),cudaMemcpyHostToDevice);
}

/* ----------------------------------------------------------------------------------------------- */

void copy_tohost_realw_test(void** d_array_addr_ptr,realw* h_array,int size) {

  // copies values onto GPU
  cudaMemcpy(h_array, (realw*) *d_array_addr_ptr,size*sizeof(realw),cudaMemcpyDeviceToHost);
}

/* ----------------------------------------------------------------------------------------------- */

// copies integer array from CPU host to GPU device
void copy_todevice_int_test(void** d_array_addr_ptr,int* h_array,int size) {

  // allocates memory on GPU
  cudaMalloc((void**)d_array_addr_ptr,size*sizeof(int));
  // copies values onto GPU
  cudaMemcpy((realw*) *d_array_addr_ptr,h_array,size*sizeof(int),cudaMemcpyHostToDevice);
}

/* ----------------------------------------------------------------------------------------------- */

void copy_tohost_int_test(void** d_array_addr_ptr,int* h_array,int size) {

  // allocates memory on GPU
  cudaMemcpy(h_array,(realw*) *d_array_addr_ptr,size*sizeof(int),cudaMemcpyDeviceToHost);
}

/* ----------------------------------------------------------------------------------------------- */

void allocate_cuda_memory_test(void** d_array_addr_ptr,int size) {

  // allocates memory on GPU
  cudaMalloc((void**)d_array_addr_ptr,size*sizeof(int));
}

/* ----------------------------------------------------------------------------------------------- */

extern EXTERN_LANG
void FC_FUNC_(transfer_fault_data_to_device,
              TRANSFER_FAULT_DATA_TO_DEVICE)(long* Fault_pointer,
                                             int* fault_index,
                                             int* NSPEC_AB,
                                             int* NGLOB_AB,
                                             realw* D,
                                             realw* T0,
                                             realw* T,
                                             realw* B,
                                             realw* R,
                                             realw* V0,
                                             realw* Z,
                                             realw* invM1,
                                             realw* invM2,
                                             int* ibulk1,
                                             int* ibulk2) {
  TRACE("transfer_fault_data_to_device");

  Fault_solver_dynamics* Fsolver = (Fault_solver_dynamics*)(*Fault_pointer);
  Fault* Flt = &(Fsolver->faults[*fault_index]);

  Flt->NSPEC_AB = *NSPEC_AB;
  Flt->NGLOB_AB = *NGLOB_AB;

  if (*NGLOB_AB > 0){
    copy_todevice_realw_test((void **)&(Flt->B),B,*NGLOB_AB);
    copy_todevice_realw_test((void **)&(Flt->R),R,(*NGLOB_AB)*9);
    copy_todevice_realw_test((void **)&(Flt->Z),Z,(*NGLOB_AB));

    copy_todevice_realw_test((void **)&(Flt->D),D,(*NGLOB_AB)*3);
    copy_todevice_realw_test((void **)&(Flt->V),V0,(*NGLOB_AB)*3);

    copy_todevice_realw_test((void **)&(Flt->T0),T0,(*NGLOB_AB)*3);
    copy_todevice_realw_test((void **)&(Flt->T),T,(*NGLOB_AB)*3);

    copy_todevice_realw_test((void **)&(Flt->invM1),invM1,*NGLOB_AB);
    copy_todevice_realw_test((void **)&(Flt->invM2),invM2,*NGLOB_AB);

    copy_todevice_int_test((void **)&(Flt->ibulk1),ibulk1,(*NGLOB_AB));
    copy_todevice_int_test((void **)&(Flt->ibulk2),ibulk2,(*NGLOB_AB));

  }

  GPU_ERROR_CHECKING("transfer_fault_data_to_device");
}

/* ----------------------------------------------------------------------------------------------- */

extern EXTERN_LANG
void FC_FUNC_(transfer_fault_data_to_host,
              TRANSFER_FAULT_DATA_TO_HOST)(long* Fault_pointer,
                                           int* fault_index,
                                           int* NSPEC_AB,
                                           int* NGLOB_AB,
                                           realw* D,
                                           realw* V,
                                           realw* T) {

  TRACE("transfer_fault_data_to_host");

  Fault_solver_dynamics* Fsolver = (Fault_solver_dynamics*)(*Fault_pointer);
  Fault* Flt = &(Fsolver->faults[*fault_index]);

  if (*NGLOB_AB > 0){
    copy_tohost_realw_test((void **)&(Flt->V),V,(*NGLOB_AB)*3);
    copy_tohost_realw_test((void **)&(Flt->D),D,(*NGLOB_AB)*3);
    copy_tohost_realw_test((void **)&(Flt->T),T,(*NGLOB_AB)*3);
  }

  GPU_ERROR_CHECKING("transfer_fault_data_to_host");
}

/* ----------------------------------------------------------------------------------------------- */

extern EXTERN_LANG
void FC_FUNC_(transfer_datat_to_host,
              TRANSFER_DATAT_TO_HOST)(long* Fault_pointer,
                                      realw* h_dataT,
                                      int* it,
                                      int* ifault) {

  TRACE("transfer_dataT_to_host");

  Fault_solver_dynamics* Fsolver = (Fault_solver_dynamics*)(*Fault_pointer);
  Fault* Flt = &(Fsolver->faults[*ifault]);

  int recordlength = 7; //store 7 different quantities

  if (Flt->output_dataT->NRECORD > 0){
    copy_tohost_realw_test((void **)&(Flt->output_dataT->dataT),
                           h_dataT + (*it - Flt->output_dataT->NT) * Flt->output_dataT->NRECORD * recordlength,
                           Flt->output_dataT->NRECORD * recordlength * Flt->output_dataT->NT);
  }

  GPU_ERROR_CHECKING("transfer_dataT_to_host");
}

/*-------------------------------------------------------------------------------------------------*/

extern EXTERN_LANG
void FC_FUNC_(transfer_rsf_data_todevice,
              TRANSFER_RSF_DATA_TODEVICE)(long* Fault_pointer,
                                          int *NGLOB_AB,
                                          int* fault_index,
                                          realw* V0,
                                          realw* f0,
                                          realw* V_init,
                                          realw* a,
                                          realw* b,
                                          realw* L,
                                          realw* theta,
                                          realw* T,
                                          realw* C,
                                          realw* fw,
                                          realw* Vw) {

  TRACE("transfer_rsf_data_todevice");

  Fault_solver_dynamics* Fsolver = (Fault_solver_dynamics*)(*Fault_pointer);
  Rsf_type* Rsf  = &((Fsolver->faults[*fault_index]).rsf);

  if (*NGLOB_AB > 0){
    copy_todevice_realw_test((void **)&(Rsf->V0),V0,*NGLOB_AB);
    copy_todevice_realw_test((void **)&(Rsf->f0),f0,*NGLOB_AB);
    copy_todevice_realw_test((void **)&(Rsf->V_init),V_init,*NGLOB_AB);
    copy_todevice_realw_test((void **)&(Rsf->a),a,*NGLOB_AB);
    copy_todevice_realw_test((void **)&(Rsf->b),b,*NGLOB_AB);
    copy_todevice_realw_test((void **)&(Rsf->L),L,*NGLOB_AB);
    copy_todevice_realw_test((void **)&(Rsf->theta),theta,*NGLOB_AB);
    copy_todevice_realw_test((void **)&(Rsf->T),T,*NGLOB_AB);
    copy_todevice_realw_test((void **)&(Rsf->C),C,*NGLOB_AB);
    copy_todevice_realw_test((void **)&(Rsf->fw),fw,*NGLOB_AB);
    copy_todevice_realw_test((void **)&(Rsf->Vw),Vw,*NGLOB_AB);
  }

  GPU_ERROR_CHECKING("transfer_rsf_data_todevice");
}

/* ----------------------------------------------------------------------------------------------- */

extern EXTERN_LANG
void FC_FUNC_(transfer_swf_data_todevice,
              TRANSFER_SWF_DATA_TODEVICE)(long* Fault_pointer,
                                          int *NGLOB_AB,
                                          int *fault_index,
                                          realw* Dc,
                                          realw* mus,
                                          realw* mud,
                                          realw* T,
                                          realw* C,
                                          realw* theta) {

  TRACE("transfer_swf_data_todevice");

  Fault_solver_dynamics* Fsolver = (Fault_solver_dynamics*)(*Fault_pointer);
  Swf_type* Swf  = &((Fsolver->faults[*fault_index]).swf);

  if (*NGLOB_AB > 0){
    copy_todevice_realw_test((void **)&(Swf->Dc),Dc,*NGLOB_AB);
    copy_todevice_realw_test((void **)&(Swf->mus),mus,*NGLOB_AB);
    copy_todevice_realw_test((void **)&(Swf->mud),mud,*NGLOB_AB);
    copy_todevice_realw_test((void **)&(Swf->Coh),C,*NGLOB_AB);
    copy_todevice_realw_test((void **)&(Swf->T),T,*NGLOB_AB);
    copy_todevice_realw_test((void **)&(Swf->theta),theta,*NGLOB_AB);
  }

  GPU_ERROR_CHECKING("transfer_swf_data_todevice");
}

/* ----------------------------------------------------------------------------------------------- */

extern EXTERN_LANG
void FC_FUNC_(transfer_rsf_data_tohost,
              TRANSFER_RSF_DATA_TOHOST)(long* Fault_pointer,
                                        int *NGLOB_AB,
                                        int *fault_index,
                                        realw* V0,
                                        realw* f0,
                                        realw* V_init,
                                        realw* a,
                                        realw* b,
                                        realw* L,
                                        realw* theta,
                                        realw* T,
                                        realw* C,
                                        realw* fw,
                                        realw* Vw) {

  TRACE("transfer_rsf_data_tohost");

  Fault_solver_dynamics* Fsolver = (Fault_solver_dynamics*)(*Fault_pointer);
  Rsf_type* Rsf  = &((Fsolver->faults[*fault_index]).rsf);

  if (*NGLOB_AB > 0){
    copy_tohost_realw_test((void **)&(Rsf->V0),V0,*NGLOB_AB);
    copy_tohost_realw_test((void **)&(Rsf->f0),f0,*NGLOB_AB);
    copy_tohost_realw_test((void **)&(Rsf->V_init),V_init,*NGLOB_AB);
    copy_tohost_realw_test((void **)&(Rsf->a),a,*NGLOB_AB);
    copy_tohost_realw_test((void **)&(Rsf->b),b,*NGLOB_AB);
    copy_tohost_realw_test((void **)&(Rsf->L),L,*NGLOB_AB);
    copy_tohost_realw_test((void **)&(Rsf->theta),theta,*NGLOB_AB);
    copy_tohost_realw_test((void **)&(Rsf->T),T,*NGLOB_AB);
    copy_tohost_realw_test((void **)&(Rsf->C),C,*NGLOB_AB);
    copy_tohost_realw_test((void **)&(Rsf->fw),fw,*NGLOB_AB);
    copy_tohost_realw_test((void **)&(Rsf->Vw),Vw,*NGLOB_AB);
  }

  GPU_ERROR_CHECKING("transfer_rsf_data_tohost");
}

/* ----------------------------------------------------------------------------------------------- */

extern EXTERN_LANG
void FC_FUNC_(transfer_swf_data_tohost,
              TRANSFER_SWF_DATA_TOHOST)(long* Fault_pointer,
                                        int *NGLOB_AB,
                                        int *fault_index,
                                        realw* Dc,
                                        realw* mus,
                                        realw* mud,
                                        realw* T,
                                        realw* C,
                                        realw* theta) {

  TRACE("transfer_swf_data_tohost");

  Fault_solver_dynamics* Fsolver = (Fault_solver_dynamics*)(*Fault_pointer);
  Swf_type *Swf = &((Fsolver -> faults[*fault_index]).swf);

  if (*NGLOB_AB > 0){
    copy_tohost_realw_test((void **)&(Swf->Dc),Dc,*NGLOB_AB);
    copy_tohost_realw_test((void **)&(Swf->mus),mus,*NGLOB_AB);
    copy_tohost_realw_test((void **)&(Swf->mud),mud,*NGLOB_AB);
    copy_tohost_realw_test((void **)&(Swf->Coh),C,*NGLOB_AB);
    copy_tohost_realw_test((void **)&(Swf->T),T,*NGLOB_AB);
    copy_tohost_realw_test((void **)&(Swf->theta),theta,*NGLOB_AB);
  }

  GPU_ERROR_CHECKING("transfer_swf_data_tohost");
}


/* ----------------------------------------------------------------------------------------------- */

extern EXTERN_LANG
void FC_FUNC_(fault_solver_gpu,
              FAULT_SOLVER_GPU)(long* Mesh_pointer,
                                long* Fault_pointer,
                                realw* dt,
                                int* myrank,
                                int* it) {
  TRACE("fault_solver_gpu");

  Mesh*  mp = (Mesh*)(*Mesh_pointer);

  Fault_solver_dynamics* Fsolver = (Fault_solver_dynamics*)(*Fault_pointer);

  for(int ifault = 0; ifault < (Fsolver->Nbfaults); ifault++){
    Fault* Flt = &(Fsolver->faults[ifault]);
    Rsf_type* rsf = &(Flt->rsf);
    Swf_type* swf = &(Flt->swf);

    int size = Flt->NGLOB_AB;

    if(size > 0){
      int blocksize = 128;
      int size_padded = ((int)ceil(((double)size)/((double)blocksize)))*blocksize;

      int num_blocks_x, num_blocks_y;
      get_blocks_xy(size_padded/blocksize,&num_blocks_x,&num_blocks_y);

      dim3 grid(num_blocks_x,num_blocks_y);
      dim3 threads(blocksize,1,1);

      /*
      // this is dirty implementation
      // fault kernel for rate and state friction
      compute_dynamic_fault_cuda_rsf<<<grid,threads>>>( mp->d_displ,
                                                        mp->d_veloc,
                                                        mp->d_accel,
                                                        Flt->NGLOB_AB,
                                                        Flt->invM1,
                                                        Flt->invM2,
                                                        Flt->B,
                                                        Flt->Z,
                                                        Flt->R,
                                                        Flt->T0,
                                                        Flt->T,
                                                        rsf->C,   // rate and state friction quantities
                                                        rsf->a,
                                                        rsf->b,
                                                        rsf->L,
                                                        rsf->f0,
                                                        rsf->V0,
                                                        rsf->V_init,
                                                        rsf->theta,
                                                        rsf->Vw,
                                                        rsf->fw,
                                                        Flt->V,
                                                        Flt->D,
                                                        Flt->ibulk1,
                                                        Flt->ibulk2,
                                                        *dt,
                                                        *myrank);
      }
      */

      // fault kernel for slip weakening friction
      compute_dynamic_fault_cuda_swf<<<grid,threads>>>( mp->d_displ,
                                                        mp->d_veloc,
                                                        mp->d_accel,
                                                        Flt->NGLOB_AB,
                                                        Flt->invM1,
                                                        Flt->invM2,
                                                        Flt->B,
                                                        Flt->Z,
                                                        Flt->R,
                                                        Flt->T0,
                                                        Flt->T,
                                                        swf->Dc,
                                                        swf->theta,
                                                        swf->mus,
                                                        swf->mud,
                                                        swf->Coh,
                                                        swf->T,
                                                        Flt->V,
                                                        Flt->D,
                                                        Flt->ibulk1,
                                                        Flt->ibulk2,
                                                        *dt,
                                                        *myrank);

      // fills output dataT array
      int num_of_block2 = (int) (Flt->output_dataT->NRECORD/128)+1;

      store_dataT<<<num_of_block2,128>>>( Flt->output_dataT->dataT,
                                          Flt->V,
                                          Flt->D,
                                          Flt->T,
                                          Flt->output_dataT->iglob,
                                          *it,
                                          Flt->output_dataT->NRECORD,
                                          Flt->output_dataT->NT);
    }
  }

  GPU_ERROR_CHECKING("fault_solver_gpu");
}

