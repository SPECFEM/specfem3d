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
#include "fault_struct_cuda.h"

/* ----------------------------------------------------------------------------------------------- */

extern EXTERN_LANG
void FC_FUNC_(initialize_fault_solver_gpu,
              INITIALIZE_FAULT_SOLVER_GPU)(long* Fault_solver,
                                           int* num_of_faults,
                                           realw* v_healing,
                                           realw* v_rupt,
                                           int* RATE_AND_STATE) {

  TRACE("initialize_fault_solver_gpu");

  // allocates fault parameter structure
  Fault_solver_dynamics *Fsolver = (Fault_solver_dynamics*) malloc(sizeof(Fault_solver_dynamics));
  if (Fsolver == NULL) exit_on_error("error allocating fault_solver pointer");

  *Fault_solver = (long) Fsolver;

  // initializes
  Fsolver->Nbfaults = *num_of_faults;
  Fsolver->faults = (Fault*) malloc((*num_of_faults)*sizeof(Fault));

  Fsolver->v_healing = *v_healing;
  Fsolver->v_rupt = *v_rupt;

  // flag for rate and state friction
  Fsolver->RATE_AND_STATE = *RATE_AND_STATE; // 0 == false, 1 == true
}

/* ----------------------------------------------------------------------------------------------- */

extern EXTERN_LANG
void FC_FUNC_(initialize_fault_data_gpu,
              INITIALIZE_FAULT_DATA_GPU)(long* Fault_solver,
                                         int* fault_index,
                                         int* iglob,
                                         int* num_of_records,
                                         int* ndat,
                                         int* NT_RECORD_LENGTH) {
  TRACE("initialize_fault_data_gpu");

  Fault_solver_dynamics* Fsolver = (Fault_solver_dynamics*)(*Fault_solver);

  // fault data recorder
  Fault_data *fault_data_recorder = (Fault_data*) malloc(sizeof(Fault_data));

  // initializes faults
  fault_data_recorder->NRECORD = *num_of_records;  // number of points npoin
  fault_data_recorder->NT = *NT_RECORD_LENGTH;     // NT_RECORD_LENGTH : set in fault_solver_dynamic.f90

  int recordlength = *ndat; // data stores for each record/point: default == 7 or when RATE_AND_STATE == 8
  fault_data_recorder->recordlength = recordlength;

  // checks
  if (Fsolver->RATE_AND_STATE){
    if (fault_data_recorder->recordlength != 8){exit_on_error("Error: invalid recordlength for fault_data_recorder, should be 8");}
  }else{
    if (fault_data_recorder->recordlength != 7){exit_on_error("Error: invalid recordlength for fault_data_recorder, should be 7");}
  }

  // allocates arrays on GPU
  gpuCreateCopy_todevice_int((void**)&(fault_data_recorder->iglob),iglob,fault_data_recorder->NRECORD);

  gpuMalloc_realw((void**)&(fault_data_recorder->dataT),recordlength*fault_data_recorder->NRECORD*fault_data_recorder->NT);

  // stores structure pointers
  Fsolver->faults[*fault_index].output_dataT = fault_data_recorder;

  GPU_ERROR_CHECKING("initialize_fault_data_gpu");
}


/* ----------------------------------------------------------------------------------------------- */

extern EXTERN_LANG
void FC_FUNC_(transfer_fault_data_to_device,
              TRANSFER_FAULT_DATA_TO_DEVICE)(long* Fault_pointer,
                                             int* fault_index,
                                             int* NSPEC_FLT,
                                             int* NGLOB_FLT,
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
                                             int* ibulk2,
                                             int* allow_opening) {
  TRACE("transfer_fault_data_to_device");

  Fault_solver_dynamics* Fsolver = (Fault_solver_dynamics*)(*Fault_pointer);
  Fault* Flt = &(Fsolver->faults[*fault_index]);

  Flt->NSPEC_FLT = *NSPEC_FLT;
  Flt->NGLOB_FLT = *NGLOB_FLT;

  // default : do not allow opening (see setting in fault_solver_common.f90)
  Flt->allow_opening = *allow_opening;
  // checks
  if (Flt->allow_opening){ exit_on_error("Fault opening not implemented yet on GPU\n"); }

  // copies data to GPU
  if (*NGLOB_FLT > 0){
    gpuCreateCopy_todevice_realw((void **)&(Flt->B),B,*NGLOB_FLT);
    gpuCreateCopy_todevice_realw((void **)&(Flt->R),R,(*NGLOB_FLT)*9);
    gpuCreateCopy_todevice_realw((void **)&(Flt->Z),Z,(*NGLOB_FLT));

    gpuCreateCopy_todevice_realw((void **)&(Flt->D),D,(*NGLOB_FLT)*3);
    gpuCreateCopy_todevice_realw((void **)&(Flt->V),V0,(*NGLOB_FLT)*3);

    gpuCreateCopy_todevice_realw((void **)&(Flt->T0),T0,(*NGLOB_FLT)*3);
    gpuCreateCopy_todevice_realw((void **)&(Flt->T),T,(*NGLOB_FLT)*3);

    gpuCreateCopy_todevice_realw((void **)&(Flt->invM1),invM1,*NGLOB_FLT);
    gpuCreateCopy_todevice_realw((void **)&(Flt->invM2),invM2,*NGLOB_FLT);

    gpuCreateCopy_todevice_int((void **)&(Flt->ibulk1),ibulk1,(*NGLOB_FLT));
    gpuCreateCopy_todevice_int((void **)&(Flt->ibulk2),ibulk2,(*NGLOB_FLT));
  }

  GPU_ERROR_CHECKING("transfer_fault_data_to_device");
}

/* ----------------------------------------------------------------------------------------------- */

extern EXTERN_LANG
void FC_FUNC_(transfer_fault_data_to_host,
              TRANSFER_FAULT_DATA_TO_HOST)(long* Fault_pointer,
                                           int* fault_index,
                                           int* NSPEC_FLT,
                                           int* NGLOB_FLT,
                                           realw* D,
                                           realw* V,
                                           realw* T) {

  TRACE("transfer_fault_data_to_host");

  Fault_solver_dynamics* Fsolver = (Fault_solver_dynamics*)(*Fault_pointer);
  Fault* Flt = &(Fsolver->faults[*fault_index]);

  if (*NGLOB_FLT > 0){
    gpuMemcpy_tohost_realw(V,Flt->V,(*NGLOB_FLT)*3);
    gpuMemcpy_tohost_realw(D,Flt->D,(*NGLOB_FLT)*3);
    gpuMemcpy_tohost_realw(T,Flt->T,(*NGLOB_FLT)*3);
  }

  GPU_ERROR_CHECKING("transfer_fault_data_to_host");
}

/* ----------------------------------------------------------------------------------------------- */

extern EXTERN_LANG
void FC_FUNC_(transfer_datat_to_host,
              TRANSFER_DATAT_TO_HOST)(long* Fault_pointer,
                                      int* fault_index,
                                      realw* h_dataT,
                                      int* it_in) {

  TRACE("transfer_dataT_to_host");

  Fault_solver_dynamics* Fsolver = (Fault_solver_dynamics*)(*Fault_pointer);
  Fault* Flt = &(Fsolver->faults[*fault_index]);

  int recordlength = Flt->output_dataT->recordlength; // stores default == 7 different quantities or when RATE_AND_STATE == 8
  int it = *it_in;

  if (Flt->output_dataT->NRECORD > 0){
    int it_index, size;

    // determines dataT array length
    int length_left = (it % Flt->output_dataT->NT);

    if (length_left == 0){
      // multiple of NT_RECORD_LENGTH == Flt->output_dataT->NT
      // copies full dataT array
      size = Flt->output_dataT->NRECORD * recordlength * Flt->output_dataT->NT;
      it_index = (it - Flt->output_dataT->NT) * Flt->output_dataT->NRECORD * recordlength;
    }else{
      // output at last step, might not be a multiple of NT_RECORD_LENGTH
      // need to copy only array of length_left
      it_index = (it - length_left) * Flt->output_dataT->NRECORD * recordlength;
      size = Flt->output_dataT->NRECORD * recordlength * length_left;
    }

    // copies dataT record to CPU
    gpuMemcpy_tohost_realw(&h_dataT[it_index],Flt->output_dataT->dataT,size);
  }

  GPU_ERROR_CHECKING("transfer_dataT_to_host");
}

/*-------------------------------------------------------------------------------------------------*/

extern EXTERN_LANG
void FC_FUNC_(transfer_rsf_data_todevice,
              TRANSFER_RSF_DATA_TODEVICE)(long* Fault_pointer,
                                          int* fault_index,
                                          int *NGLOB_FLT,
                                          realw* Fload,
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
                                          realw* Vw,
                                          int* StateLaw) {

  TRACE("transfer_rsf_data_todevice");

  Fault_solver_dynamics* Fsolver = (Fault_solver_dynamics*)(*Fault_pointer);
  Rsf_type* rsf  = &((Fsolver->faults[*fault_index]).rsf);

  // checks with rsf flag
  if (! Fsolver->RATE_AND_STATE){ exit_on_error("Error with RSF setup, RATE_AND_STATE flag is off; please check fault setup and rerun\n");}

  // set StateLaw type: 1 == ageing law, 2 == slip law
  rsf->StateLaw = *StateLaw;

  // copies arrays onto GPU
  if (*NGLOB_FLT > 0){
    gpuCreateCopy_todevice_realw((void **)&(rsf->V0),V0,*NGLOB_FLT);
    gpuCreateCopy_todevice_realw((void **)&(rsf->f0),f0,*NGLOB_FLT);
    gpuCreateCopy_todevice_realw((void **)&(rsf->V_init),V_init,*NGLOB_FLT);
    gpuCreateCopy_todevice_realw((void **)&(rsf->a),a,*NGLOB_FLT);
    gpuCreateCopy_todevice_realw((void **)&(rsf->b),b,*NGLOB_FLT);
    gpuCreateCopy_todevice_realw((void **)&(rsf->L),L,*NGLOB_FLT);
    gpuCreateCopy_todevice_realw((void **)&(rsf->theta),theta,*NGLOB_FLT);
    gpuCreateCopy_todevice_realw((void **)&(rsf->T),T,*NGLOB_FLT);
    gpuCreateCopy_todevice_realw((void **)&(rsf->Coh),C,*NGLOB_FLT);
    gpuCreateCopy_todevice_realw((void **)&(rsf->fw),fw,*NGLOB_FLT);
    gpuCreateCopy_todevice_realw((void **)&(rsf->Vw),Vw,*NGLOB_FLT);
    gpuCreateCopy_todevice_realw((void **)&(rsf->Fload),Fload,*NGLOB_FLT);
  }

  GPU_ERROR_CHECKING("transfer_rsf_data_todevice");
}

/* ----------------------------------------------------------------------------------------------- */

extern EXTERN_LANG
void FC_FUNC_(transfer_swf_data_todevice,
              TRANSFER_SWF_DATA_TODEVICE)(long* Fault_pointer,
                                          int *fault_index,
                                          int *NGLOB_FLT,
                                          realw* Dc,
                                          realw* mus,
                                          realw* mud,
                                          realw* T,
                                          realw* C,
                                          realw* theta) {

  TRACE("transfer_swf_data_todevice");

  Fault_solver_dynamics* Fsolver = (Fault_solver_dynamics*)(*Fault_pointer);
  Swf_type* swf  = &((Fsolver->faults[*fault_index]).swf);

  // checks with rsf flag
  if (Fsolver->RATE_AND_STATE){ exit_on_error("Error with SWF setup, RATE_AND_STATE flag is on; please check fault setup and rerun\n");}

  if (*NGLOB_FLT > 0){
    gpuCreateCopy_todevice_realw((void **)&(swf->Dc),Dc,*NGLOB_FLT);
    gpuCreateCopy_todevice_realw((void **)&(swf->mus),mus,*NGLOB_FLT);
    gpuCreateCopy_todevice_realw((void **)&(swf->mud),mud,*NGLOB_FLT);
    gpuCreateCopy_todevice_realw((void **)&(swf->Coh),C,*NGLOB_FLT);
    gpuCreateCopy_todevice_realw((void **)&(swf->T),T,*NGLOB_FLT);
    gpuCreateCopy_todevice_realw((void **)&(swf->theta),theta,*NGLOB_FLT);
  }

  GPU_ERROR_CHECKING("transfer_swf_data_todevice");
}

/* ----------------------------------------------------------------------------------------------- */

// not used yet...

extern EXTERN_LANG
void FC_FUNC_(transfer_rsf_data_tohost,
              TRANSFER_RSF_DATA_TOHOST)(long* Fault_pointer,
                                        int *fault_index,
                                        int *NGLOB_FLT,
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
  Rsf_type* rsf  = &((Fsolver->faults[*fault_index]).rsf);

  if (*NGLOB_FLT > 0){
    gpuMemcpy_tohost_realw(V0,rsf->V0,*NGLOB_FLT);
    gpuMemcpy_tohost_realw(f0,rsf->f0,*NGLOB_FLT);
    gpuMemcpy_tohost_realw(V_init,rsf->V_init,*NGLOB_FLT);
    gpuMemcpy_tohost_realw(a,rsf->a,*NGLOB_FLT);
    gpuMemcpy_tohost_realw(b,rsf->b,*NGLOB_FLT);
    gpuMemcpy_tohost_realw(L,rsf->L,*NGLOB_FLT);
    gpuMemcpy_tohost_realw(theta,rsf->theta,*NGLOB_FLT);
    gpuMemcpy_tohost_realw(T,rsf->T,*NGLOB_FLT);
    gpuMemcpy_tohost_realw(C,rsf->Coh,*NGLOB_FLT);
    gpuMemcpy_tohost_realw(fw,rsf->fw,*NGLOB_FLT);
    gpuMemcpy_tohost_realw(Vw,rsf->Vw,*NGLOB_FLT);
  }

  GPU_ERROR_CHECKING("transfer_rsf_data_tohost");
}

/* ----------------------------------------------------------------------------------------------- */

// not used yet...

extern EXTERN_LANG
void FC_FUNC_(transfer_swf_data_tohost,
              TRANSFER_SWF_DATA_TOHOST)(long* Fault_pointer,
                                        int *fault_index,
                                        int *NGLOB_FLT,
                                        realw* Dc,
                                        realw* mus,
                                        realw* mud,
                                        realw* T,
                                        realw* C,
                                        realw* theta) {

  TRACE("transfer_swf_data_tohost");

  Fault_solver_dynamics* Fsolver = (Fault_solver_dynamics*)(*Fault_pointer);
  Swf_type *swf = &((Fsolver -> faults[*fault_index]).swf);

  if (*NGLOB_FLT > 0){
    gpuMemcpy_tohost_realw(Dc,swf->Dc,*NGLOB_FLT);
    gpuMemcpy_tohost_realw(mus,swf->mus,*NGLOB_FLT);
    gpuMemcpy_tohost_realw(mud,swf->mud,*NGLOB_FLT);
    gpuMemcpy_tohost_realw(C,swf->Coh,*NGLOB_FLT);
    gpuMemcpy_tohost_realw(T,swf->T,*NGLOB_FLT);
    gpuMemcpy_tohost_realw(theta,swf->theta,*NGLOB_FLT);
  }

  GPU_ERROR_CHECKING("transfer_swf_data_tohost");
}


/* ----------------------------------------------------------------------------------------------- */

extern EXTERN_LANG
void FC_FUNC_(fault_solver_gpu,
              FAULT_SOLVER_GPU)(long* Mesh_pointer,
                                long* Fault_pointer,
                                realw* dt,
                                int* it) {
  TRACE("fault_solver_gpu");

  Mesh*  mp = (Mesh*)(*Mesh_pointer);

  Fault_solver_dynamics* Fsolver = (Fault_solver_dynamics*)(*Fault_pointer);

  for(int ifault = 0; ifault < (Fsolver->Nbfaults); ifault++){
    Fault* Flt = &(Fsolver->faults[ifault]);

    Rsf_type* rsf = &(Flt->rsf);
    Swf_type* swf = &(Flt->swf);

    int size = Flt->NGLOB_FLT;

    if(size > 0){
      int blocksize = 128;
      int size_padded = ((int)ceil(((double)size)/((double)blocksize)))*blocksize;

      int num_blocks_x, num_blocks_y;
      get_blocks_xy(size_padded/blocksize,&num_blocks_x,&num_blocks_y);

      dim3 grid(num_blocks_x,num_blocks_y);
      dim3 threads(blocksize,1,1);

      if (Fsolver->RATE_AND_STATE){
        // this is dirty implementation
        // fault kernel for rate and state friction
#ifdef USE_CUDA
        if (run_cuda){
          compute_dynamic_fault_cuda_rsf<<<grid,threads>>>( mp->d_displ,
                                                            mp->d_veloc,
                                                            mp->d_accel,
                                                            Flt->NGLOB_FLT,
                                                            Flt->invM1,
                                                            Flt->invM2,
                                                            Flt->B,
                                                            Flt->Z,
                                                            Flt->R,
                                                            Flt->T0,
                                                            Flt->T,
                                                            rsf->Coh,       // rate and state friction quantities
                                                            rsf->a,
                                                            rsf->b,
                                                            rsf->L,
                                                            rsf->f0,
                                                            rsf->V0,
                                                            rsf->V_init,
                                                            rsf->theta,
                                                            rsf->Vw,
                                                            rsf->fw,
                                                            rsf->Fload,
                                                            rsf->StateLaw,
                                                            Flt->V,
                                                            Flt->D,
                                                            Flt->ibulk1,
                                                            Flt->ibulk2,
                                                            *dt,
                                                            *it);
        }
#endif
#ifdef USE_HIP
        if (run_hip){
          hipLaunchKernelGGL(compute_dynamic_fault_cuda_rsf, dim3(grid), dim3(threads), 0, 0,
                                                             mp->d_displ,
                                                             mp->d_veloc,
                                                             mp->d_accel,
                                                             Flt->NGLOB_FLT,
                                                             Flt->invM1,
                                                             Flt->invM2,
                                                             Flt->B,
                                                             Flt->Z,
                                                             Flt->R,
                                                             Flt->T0,
                                                             Flt->T,
                                                             rsf->Coh,       // rate and state friction quantities
                                                             rsf->a,
                                                             rsf->b,
                                                             rsf->L,
                                                             rsf->f0,
                                                             rsf->V0,
                                                             rsf->V_init,
                                                             rsf->theta,
                                                             rsf->Vw,
                                                             rsf->fw,
                                                             rsf->Fload,
                                                             rsf->StateLaw,
                                                             Flt->V,
                                                             Flt->D,
                                                             Flt->ibulk1,
                                                             Flt->ibulk2,
                                                             *dt,
                                                             *it);
        }
#endif

      } else {
        // fault kernel for slip weakening friction
#ifdef USE_CUDA
        if (run_cuda){
          compute_dynamic_fault_cuda_swf<<<grid,threads>>>( mp->d_displ,
                                                            mp->d_veloc,
                                                            mp->d_accel,
                                                            Flt->NGLOB_FLT,
                                                            Flt->invM1,
                                                            Flt->invM2,
                                                            Flt->B,
                                                            Flt->Z,
                                                            Flt->R,
                                                            Flt->T0,
                                                            Flt->T,
                                                            swf->Dc,        // slip weakening friction quantities
                                                            swf->theta,
                                                            swf->mus,
                                                            swf->mud,
                                                            swf->Coh,
                                                            swf->T,
                                                            Flt->V,
                                                            Flt->D,
                                                            Flt->ibulk1,
                                                            Flt->ibulk2,
                                                            *dt);
        }
#endif
#ifdef USE_HIP
        if (run_hip){
          hipLaunchKernelGGL(compute_dynamic_fault_cuda_swf, dim3(grid), dim3(threads), 0, 0,
                                                             mp->d_displ,
                                                             mp->d_veloc,
                                                             mp->d_accel,
                                                             Flt->NGLOB_FLT,
                                                             Flt->invM1,
                                                             Flt->invM2,
                                                             Flt->B,
                                                             Flt->Z,
                                                             Flt->R,
                                                             Flt->T0,
                                                             Flt->T,
                                                             swf->Dc,        // slip weakening friction quantities
                                                             swf->theta,
                                                             swf->mus,
                                                             swf->mud,
                                                             swf->Coh,
                                                             swf->T,
                                                             Flt->V,
                                                             Flt->D,
                                                             Flt->ibulk1,
                                                             Flt->ibulk2,
                                                             *dt);
        }
#endif

      }

      // output dataT array
      int num_of_block2 = (int) (Flt->output_dataT->NRECORD/128)+1;

      // case for rate and state friction
      int StateLaw = 0;
      realw* theta = NULL;
      if (Fsolver->RATE_AND_STATE){
        StateLaw = rsf->StateLaw;
        theta = rsf->theta;
      }

      // fills dataT arrays on GPU
#ifdef USE_CUDA
      if (run_cuda){
        store_dataT<<<num_of_block2,128>>>( Flt->output_dataT->dataT,
                                            Flt->V,
                                            Flt->D,
                                            Flt->T,
                                            Fsolver->RATE_AND_STATE,
                                            StateLaw,
                                            theta,
                                            Flt->output_dataT->iglob,
                                            *it,
                                            Flt->output_dataT->NRECORD,
                                            Flt->output_dataT->NT);
      }
#endif
#ifdef USE_HIP
      if (run_hip){
        hipLaunchKernelGGL(store_dataT, dim3(num_of_block2), dim3(128), 0, 0,
                                        Flt->output_dataT->dataT,
                                        Flt->V,
                                        Flt->D,
                                        Flt->T,
                                        Fsolver->RATE_AND_STATE,
                                        StateLaw,
                                        theta,
                                        Flt->output_dataT->iglob,
                                        *it,
                                        Flt->output_dataT->NRECORD,
                                        Flt->output_dataT->NT);
      }
#endif

    }
  }

  GPU_ERROR_CHECKING("fault_solver_gpu");
}

