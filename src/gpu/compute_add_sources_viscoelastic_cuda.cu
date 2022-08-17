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

#include "mesh_constants_gpu.h"

/* ----------------------------------------------------------------------------------------------- */

// elastic domain sources

/* ----------------------------------------------------------------------------------------------- */

extern EXTERN_LANG
void FC_FUNC_(compute_add_sources_el_cuda,
              COMPUTE_ADD_SOURCES_EL_CUDA)(long* Mesh_pointer,
                                           double* h_stf_pre_compute,
                                           int* h_NSOURCES) {

  TRACE("compute_add_sources_el_cuda");

  Mesh* mp = (Mesh*)(*Mesh_pointer); //get mesh pointer out of fortran integer container

  // check if anything to do
  if (mp->nsources_local == 0) return;

  int NSOURCES = *h_NSOURCES;

  // convert to GPU precision
  realw* stf_pre_compute = (realw*) malloc(NSOURCES * sizeof(realw));
  for (int i_source=0;i_source < NSOURCES;i_source++) stf_pre_compute[i_source] = (realw) h_stf_pre_compute[i_source];

  // note: d_stf_pre_compute has field type, as it is used also for acoustic simulations.
  //       to use this also for elastic simulations, one would need to do:
  //          field* stf_pre_compute = (field*) malloc(NSOURCES * sizeof(field));
  //          get_stf_for_gpu(stf_pre_compute,h_stf_pre_compute,run_number_of_the_source,NSOURCES);
  //       however, NB_RUNS_ACOUSTIC_GPU > 1 is not supported by elastic sources, thus field declaration is realw by default.
  //       to avoid compilation issues we use the copy function for (void*) and exact byte size.
  gpuMemcpy_todevice_void((void*)mp->d_stf_pre_compute,(void*)stf_pre_compute,NSOURCES*sizeof(realw));

  free(stf_pre_compute);

  GPU_ERROR_CHECKING("compute_add_sources_el_cuda copy");

  int num_blocks_x, num_blocks_y;
  get_blocks_xy(NSOURCES,&num_blocks_x,&num_blocks_y);

  dim3 grid(num_blocks_x,num_blocks_y);
  dim3 threads(NGLLX,NGLLY,NGLLZ);

#ifdef USE_CUDA
  if (run_cuda){
    compute_add_sources_kernel<<<grid,threads,0,mp->compute_stream>>>(mp->d_accel,mp->d_ibool,
                                                                      mp->d_sourcearrays,
                                                                      mp->d_stf_pre_compute,
                                                                      mp->myrank,
                                                                      mp->d_islice_selected_source,mp->d_ispec_selected_source,
                                                                      mp->d_ispec_is_elastic,
                                                                      NSOURCES);
  }
#endif
#ifdef USE_HIP
  if (run_hip){
    hipLaunchKernelGGL(compute_add_sources_kernel, dim3(grid), dim3(threads), 0, mp->compute_stream,
                                                   mp->d_accel,mp->d_ibool,
                                                   mp->d_sourcearrays,
                                                   mp->d_stf_pre_compute,
                                                   mp->myrank,
                                                   mp->d_islice_selected_source,mp->d_ispec_selected_source,
                                                   mp->d_ispec_is_elastic,
                                                   NSOURCES);
  }
#endif

  GPU_ERROR_CHECKING("compute_add_sources_el_cuda");
}

/* ----------------------------------------------------------------------------------------------- */

extern EXTERN_LANG
void FC_FUNC_(compute_add_sources_el_s3_cuda,
              COMPUTE_ADD_SOURCES_EL_S3_CUDA)(long* Mesh_pointer,
                                              double* h_stf_pre_compute,
                                              int* h_NSOURCES) {

  TRACE("compute_add_sources_el_s3_cuda");
  // EPIK_TRACER("compute_add_sources_el_s3_cuda");

  Mesh* mp = (Mesh*)(*Mesh_pointer); //get mesh pointer out of fortran integer container

  int NSOURCES = *h_NSOURCES;

  // convert to GPU precision
  realw* stf_pre_compute = (realw*) malloc(NSOURCES * sizeof(realw));
  for (int i_source=0;i_source < NSOURCES;i_source++) stf_pre_compute[i_source] = (realw) h_stf_pre_compute[i_source];

  // note: d_stf_pre_compute has field type, as it is used also for acoustic simulations.
  //       to use this also for elastic simulations, one would need to do:
  //          field* stf_pre_compute = (field*) malloc(NSOURCES * sizeof(field));
  //          get_stf_for_gpu(stf_pre_compute,h_stf_pre_compute,run_number_of_the_source,NSOURCES);
  //       however, NB_RUNS_ACOUSTIC_GPU > 1 is not supported by elastic sources, thus field declaration is realw by default.
  //       to avoid compilation issues we use the copy function for (void*) and exact byte size.
  gpuMemcpy_todevice_void((void*)mp->d_stf_pre_compute,(void*)stf_pre_compute,NSOURCES*sizeof(realw));

  free(stf_pre_compute);

  GPU_ERROR_CHECKING("compute_add_sources_el_s3_cuda copy");

  int num_blocks_x, num_blocks_y;
  get_blocks_xy(NSOURCES,&num_blocks_x,&num_blocks_y);

  dim3 grid(num_blocks_x,num_blocks_y);
  dim3 threads(NGLLX,NGLLY,NGLLZ);

#ifdef USE_CUDA
  if (run_cuda){
    compute_add_sources_kernel<<<grid,threads,0,mp->compute_stream>>>(mp->d_b_accel,mp->d_ibool,
                                                                      mp->d_sourcearrays,
                                                                      mp->d_stf_pre_compute,
                                                                      mp->myrank,
                                                                      mp->d_islice_selected_source,mp->d_ispec_selected_source,
                                                                      mp->d_ispec_is_elastic,
                                                                      NSOURCES);
  }
#endif
#ifdef USE_HIP
  if (run_hip){
    hipLaunchKernelGGL(compute_add_sources_kernel, dim3(grid), dim3(threads), 0, mp->compute_stream,
                                                   mp->d_accel,mp->d_ibool,
                                                   mp->d_sourcearrays,
                                                   mp->d_stf_pre_compute,
                                                   mp->myrank,
                                                   mp->d_islice_selected_source,mp->d_ispec_selected_source,
                                                   mp->d_ispec_is_elastic,
                                                   NSOURCES);
  }
#endif

  GPU_ERROR_CHECKING("compute_add_sources_el_s3_cuda");
}

/* ----------------------------------------------------------------------------------------------- */

// NOISE sources

/* ----------------------------------------------------------------------------------------------- */

extern EXTERN_LANG
void FC_FUNC_(add_source_main_rec_noise_cu,
              ADD_SOURCE_MAIN_REC_NOISE_CU)(long* Mesh_pointer,
                                            int* it_f,
                                            int* irec_main_noise_f,
                                            int* islice_selected_rec) {

TRACE("add_source_main_rec_noise_cu");

  Mesh* mp = (Mesh*)(*Mesh_pointer); //get mesh pointer out of fortran integer container

  int it = *it_f-1; // -1 for Fortran -> C indexing differences
  int irec_main_noise = *irec_main_noise_f;

  dim3 grid(1,1,1);
  dim3 threads(NGLL3,1,1);

  if (mp->myrank == islice_selected_rec[irec_main_noise-1]) {

#ifdef USE_CUDA
    if (run_cuda){
      add_source_main_rec_noise_cuda_kernel<<<grid,threads,0,mp->compute_stream>>>(mp->d_ibool,
                                                                                   mp->d_ispec_selected_rec,
                                                                                   irec_main_noise,
                                                                                   mp->d_accel,
                                                                                   mp->d_noise_sourcearray,
                                                                                   it);
    }
#endif
#ifdef USE_HIP
    if (run_hip){
      hipLaunchKernelGGL(add_source_main_rec_noise_cuda_kernel, dim3(grid), dim3(threads), 0, mp->compute_stream,
                                                                mp->d_ibool,
                                                                mp->d_ispec_selected_rec,
                                                                irec_main_noise,
                                                                mp->d_accel,
                                                                mp->d_noise_sourcearray,
                                                                it);
    }
#endif

  GPU_ERROR_CHECKING("add_source_main_rec_noise_cuda_kernel");
  }
}

/* ----------------------------------------------------------------------------------------------- */

// ADJOINT sources

/* ----------------------------------------------------------------------------------------------- */

extern EXTERN_LANG
void FC_FUNC_(add_sources_el_sim_type_2_or_3,
              ADD_SOURCES_EL_SIM_TYPE_2_OR_3)(long* Mesh_pointer,
                                              realw* h_source_adjoint,
                                              int* nrec,
                                              int* nadj_rec_local,
                                              int* NTSTEP_BETWEEN_READ_ADJSRC,
                                              int* it) {

  TRACE("add_sources_el_sim_type_2_or_3");

  Mesh* mp = (Mesh*)(*Mesh_pointer); //get mesh pointer out of fortran integer container

  // checks
  if (*nadj_rec_local != mp->nadj_rec_local) exit_on_error("add_sources_el_sim_type_2_or_3: nadj_rec_local not equal\n");

  // checks if anything to do
  if (mp->nadj_rec_local == 0) return;

  int num_blocks_x, num_blocks_y;
  get_blocks_xy(mp->nadj_rec_local,&num_blocks_x,&num_blocks_y);

  dim3 grid(num_blocks_x,num_blocks_y,1);
  dim3 threads(NGLLX,NGLLY,NGLLZ);

  int it_index = *NTSTEP_BETWEEN_READ_ADJSRC - (*it-1) % *NTSTEP_BETWEEN_READ_ADJSRC - 1 ;

  // copies extracted array values onto GPU
  if ( (*it-1) % *NTSTEP_BETWEEN_READ_ADJSRC==0){
    // note: field declaration is only equal to realw if NB_RUNS_ACOUSTIC_GPU == 1.
    //       for any other setting of NB_RUNS_ACOUSTIC_GPU, the compilation would fail for
    //          gpuMemcpy_todevice_field(mp->d_source_adjoint,h_source_adjoint,mp->nadj_rec_local*NDIM*(*NTSTEP_BETWEEN_READ_ADJSRC));
    //       since the host array is still defined as realw. in that case, we will need to construct a field array first.
    //       however, the case with NB_RUNS_ACOUSTIC_GPU > 1 is not fully implemented yet for adjoint/kernels simulations and elastic cases,
    //       and on the todo for future ...
    // copies adjoint source array onto GPU using (void*) as variable and actual byte size to avoid compilation errors
    gpuMemcpy_todevice_void((void*)mp->d_source_adjoint,(void*)h_source_adjoint,
                             mp->nadj_rec_local*NDIM*(*NTSTEP_BETWEEN_READ_ADJSRC)*sizeof(realw));
  }

#ifdef USE_CUDA
  if (run_cuda){
    add_sources_el_SIM_TYPE_2_OR_3_kernel<<<grid,threads,0,mp->compute_stream>>>(mp->d_accel,
                                                                                 *nrec,it_index,*NTSTEP_BETWEEN_READ_ADJSRC,
                                                                                 mp->d_source_adjoint,
                                                                                 mp->d_hxir_adj,
                                                                                 mp->d_hetar_adj,
                                                                                 mp->d_hgammar_adj,
                                                                                 mp->d_ibool,
                                                                                 mp->d_ispec_is_elastic,
                                                                                 mp->d_ispec_selected_adjrec_loc,
                                                                                 mp->nadj_rec_local);
  }
#endif
#ifdef USE_HIP
  if (run_hip){
    hipLaunchKernelGGL(add_sources_el_SIM_TYPE_2_OR_3_kernel, dim3(grid), dim3(threads), 0, mp->compute_stream,
                                                              mp->d_accel,
                                                              *nrec,it_index,*NTSTEP_BETWEEN_READ_ADJSRC,
                                                              mp->d_source_adjoint,
                                                              mp->d_hxir_adj,
                                                              mp->d_hetar_adj,
                                                              mp->d_hgammar_adj,
                                                              mp->d_ibool,
                                                              mp->d_ispec_is_elastic,
                                                              mp->d_ispec_selected_rec_loc,
                                                              mp->nadj_rec_local);
  }
#endif

  GPU_ERROR_CHECKING("add_sources_SIM_TYPE_2_OR_3_kernel");
}

