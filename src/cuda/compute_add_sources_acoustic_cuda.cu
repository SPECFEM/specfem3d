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

/* ----------------------------------------------------------------------------------------------- */

// acoustic sources

/* ----------------------------------------------------------------------------------------------- */

__global__ void compute_add_sources_acoustic_kernel(field* potential_dot_dot_acoustic,
                                                    int* d_ibool,
                                                    realw* sourcearrays,
                                                    field* stf_pre_compute,
                                                    int myrank,
                                                    int* islice_selected_source,
                                                    int* ispec_selected_source,
                                                    int* ispec_is_acoustic,
                                                    realw* kappastore,
                                                    int NSOURCES) {
  int i = threadIdx.x;
  int j = threadIdx.y;
  int k = threadIdx.z;

  int isource  = blockIdx.x + gridDim.x*blockIdx.y; // bx

  int ispec,iglob;
  field stf;
  realw kappal;

  if (isource < NSOURCES){

    if (myrank == islice_selected_source[isource]) {

      ispec = ispec_selected_source[isource]-1;

      if (ispec_is_acoustic[ispec]) {

        iglob = d_ibool[INDEX4_PADDED(NGLLX,NGLLX,NGLLX,i,j,k,ispec)] - 1;

        stf = stf_pre_compute[isource];
        kappal = kappastore[INDEX4(NGLLX,NGLLX,NGLLX,i,j,k,ispec)];

        atomicAdd(&potential_dot_dot_acoustic[iglob],
                  -sourcearrays[INDEX5(NSOURCES,NDIM,NGLLX,NGLLX,isource, 0,i,j,k)]*stf/kappal);

        // debug: without atomic operation
        //      potential_dot_dot_acoustic[iglob] +=
        //              -sourcearrays[INDEX5(NSOURCES, 3, NGLLX,NGLLX,isource, 0, i,j,k)]*stf/kappal;
      }
    }
  }
}

// Converts source time function to the correct GPU precision, and adapts format for NB_RUNS_ON_ACOUSTIC_GPU option
void get_stf_for_gpu(field* stf_pre_compute, double* h_stf_pre_compute, int * run_number_of_the_source, int NSOURCES) {

  TRACE("get_stf_for_gpu");
  realw realw_to_field[NB_RUNS_ACOUSTIC_GPU];
  //Conversion to GPU precision
  //Converts source time function to the field format. The stf value is saved only into its corresponding run. For other runs, a zero will be added


  for (int i_source=0;i_source < NSOURCES;i_source++){
    for (int i_run=0;i_run < NB_RUNS_ACOUSTIC_GPU;i_run++)
      if (run_number_of_the_source[i_source] == i_run){
        realw_to_field[i_run]= (realw)h_stf_pre_compute[i_source];
      }
      else{
        realw_to_field[i_run] = 0.0f;
      }
      //function Make_field is overloaded to convert array of realw into field structure
      stf_pre_compute[i_source] = Make_field(realw_to_field);
  }
}
/* ----------------------------------------------------------------------------------------------- */

extern "C"
void FC_FUNC_(compute_add_sources_ac_cuda,
              COMPUTE_ADD_SOURCES_AC_CUDA)(long* Mesh_pointer,
                                           int* NSOURCESf,
                                           double* h_stf_pre_compute,int* run_number_of_the_source) {

  TRACE("compute_add_sources_ac_cuda");

  Mesh* mp = (Mesh*)(*Mesh_pointer); //get mesh pointer out of fortran integer container

  // check if anything to do
  if (mp->nsources_local == 0) return;

  int NSOURCES = *NSOURCESf;


  field* stf_pre_compute = (field*)malloc(NSOURCES * sizeof(field));
  get_stf_for_gpu(stf_pre_compute,h_stf_pre_compute,run_number_of_the_source,NSOURCES);

  // copies pre-computed source time factors onto GPU
  print_CUDA_error_if_any(cudaMemcpy(mp->d_stf_pre_compute,stf_pre_compute,
                                     NSOURCES*sizeof(field),cudaMemcpyHostToDevice),1877);
  free(stf_pre_compute);

  int num_blocks_x, num_blocks_y;
  get_blocks_xy(NSOURCES,&num_blocks_x,&num_blocks_y);

  dim3 grid(num_blocks_x,num_blocks_y);
  dim3 threads(NGLLX,NGLLY,NGLLZ);

  compute_add_sources_acoustic_kernel<<<grid,threads,0,mp->compute_stream>>>(mp->d_potential_dot_dot_acoustic,
                                                                              mp->d_ibool,
                                                                              mp->d_sourcearrays,
                                                                              mp->d_stf_pre_compute,
                                                                              mp->myrank,
                                                                              mp->d_islice_selected_source,
                                                                              mp->d_ispec_selected_source,
                                                                              mp->d_ispec_is_acoustic,
                                                                              mp->d_kappastore,
                                                                              NSOURCES);


#ifdef ENABLE_VERY_SLOW_ERROR_CHECKING
  exit_on_cuda_error("compute_add_sources_ac_cuda");
#endif
}

/* ----------------------------------------------------------------------------------------------- */

extern "C"
void FC_FUNC_(compute_add_sources_ac_s3_cuda,
              COMPUTE_ADD_SOURCES_AC_s3_CUDA)(long* Mesh_pointer,
                                              int* NSOURCESf,
                                              double* h_stf_pre_compute,int* run_number_of_the_source) {

  TRACE("compute_add_sources_ac_s3_cuda");

  Mesh* mp = (Mesh*)(*Mesh_pointer); //get mesh pointer out of fortran integer container

  // check if anything to do
  if (mp->nsources_local == 0) return;

  int NSOURCES = *NSOURCESf;

  field* stf_pre_compute = (field*)malloc(NSOURCES * sizeof(field));
  get_stf_for_gpu(stf_pre_compute,h_stf_pre_compute,run_number_of_the_source,NSOURCES);

  // copies source time factors onto GPU
  print_CUDA_error_if_any(cudaMemcpy(mp->d_stf_pre_compute,stf_pre_compute,
                                     NSOURCES*sizeof(field),cudaMemcpyHostToDevice),55);

  free(stf_pre_compute);

  int num_blocks_x, num_blocks_y;
  get_blocks_xy(NSOURCES,&num_blocks_x,&num_blocks_y);

  dim3 grid(num_blocks_x,num_blocks_y);
  dim3 threads(NGLLX,NGLLY,NGLLZ);

  compute_add_sources_acoustic_kernel<<<grid,threads,0,mp->compute_stream>>>(mp->d_b_potential_dot_dot_acoustic,
                                                                              mp->d_ibool,
                                                                              mp->d_sourcearrays,
                                                                              mp->d_stf_pre_compute,
                                                                              mp->myrank,
                                                                              mp->d_islice_selected_source,
                                                                              mp->d_ispec_selected_source,
                                                                              mp->d_ispec_is_acoustic,
                                                                              mp->d_kappastore,
                                                                              NSOURCES);

#ifdef ENABLE_VERY_SLOW_ERROR_CHECKING
  exit_on_cuda_error("compute_add_sources_ac_s3_cuda");
#endif
}


/* ----------------------------------------------------------------------------------------------- */

// acoustic adjoint sources

/* ----------------------------------------------------------------------------------------------- */

__global__ void add_sources_ac_SIM_TYPE_2_OR_3_kernel(field* potential_dot_dot_acoustic,
                                                      int nrec,
                                                      int it,
                                                      int NSTEP_BETWEEN_ADJSRC,
                                                      field* source_adjoint,
                                                      realw* xir_store,
                                                      realw* etar_store,
                                                      realw* gammar_store,
                                                      int* d_ibool,
                                                      int* ispec_is_acoustic,
                                                      int* ispec_selected_recloc,
                                                      int nadj_rec_local,
                                                      realw* kappastore) {

  int irec_local = blockIdx.x + gridDim.x*blockIdx.y;

  // because of grid shape, irec_local can be too big
  if (irec_local < nadj_rec_local) {

    int ispec = ispec_selected_recloc[irec_local]-1;
    if (ispec_is_acoustic[ispec]){
      int i = threadIdx.x;
      int j = threadIdx.y;
      int k = threadIdx.z;

      int iglob = d_ibool[INDEX4_PADDED(NGLLX,NGLLX,NGLLX,i,j,k,ispec)]-1;

      realw xir    = xir_store[INDEX2(nadj_rec_local,irec_local,i)];
      realw etar   = etar_store[INDEX2(nadj_rec_local,irec_local,j)];
      realw gammar = gammar_store[INDEX2(nadj_rec_local,irec_local,k)];
      field source_adj = source_adjoint[INDEX3(NDIM,nadj_rec_local,0,irec_local,it)];
      //realw kappal = kappastore[INDEX4(NGLLX,NGLLY,NGLLZ,i,j,k,ispec)];

      //potential_dot_dot_acoustic[iglob] += adj_sourcearrays[INDEX6(nadj_rec_local,NTSTEP_BETWEEN_ADJSRC,3,5,5,
      //                                            pre_computed_irec_local_index[irec],
      //                                            pre_computed_index,
      //                                            0,
      //                                            i,j,k)]/kappal;

      // beware, for acoustic medium, a pressure source would be taking the negative
      // and divide by Kappa of the fluid;
      //
      // note: we take the first component of the adj_sourcearrays

      //realw stf = - source_adj * xir * etar * gammar / kappal;

      // VM VM : change the adjoint source to be consistent with CPU code
      field stf = source_adj * xir * etar * gammar;
      atomicAdd(&potential_dot_dot_acoustic[iglob],stf);

                //+adj_sourcearrays[INDEX6(nadj_rec_local,NTSTEP_BETWEEN_ADJSRC,3,5,5,
                //                         pre_computed_irec_local_index[irec],pre_computed_index-1,
                //                         0,i,j,k)] // / kappal
                //                         );
    }
  }
}

/* ----------------------------------------------------------------------------------------------- */


extern "C"
void FC_FUNC_(add_sources_ac_sim_2_or_3_cuda,
              ADD_SOURCES_AC_SIM_2_OR_3_CUDA)(long* Mesh_pointer,
                                              realw* h_source_adjoint,
                                              int* nrec,
                                              int* nadj_rec_local,
                                              int* NTSTEP_BETWEEN_READ_ADJSRC,
                                              int* it) {

  TRACE("add_sources_ac_sim_2_or_3_cuda");

  Mesh* mp = (Mesh*)(*Mesh_pointer); //get mesh pointer out of fortran integer container

  // checks
  if (*nadj_rec_local/NB_RUNS_ACOUSTIC_GPU != mp->nadj_rec_local) exit_on_cuda_error("add_sources_ac_sim_type_2_or_3: nadj_rec_local not equal\n");

  int num_blocks_x, num_blocks_y;
  get_blocks_xy(mp->nadj_rec_local,&num_blocks_x,&num_blocks_y);

  dim3 grid(num_blocks_x,num_blocks_y,1);
  dim3 threads(NGLLX,NGLLY,NGLLZ);
  int it_index = *NTSTEP_BETWEEN_READ_ADJSRC - (*it-1) % *NTSTEP_BETWEEN_READ_ADJSRC - 1 ;
  // copies extracted array values onto GPU
  if ( (*it-1) % *NTSTEP_BETWEEN_READ_ADJSRC==0) print_CUDA_error_if_any(cudaMemcpy(mp->d_source_adjoint,h_source_adjoint,
                                                                mp->nadj_rec_local*3*sizeof(field)*(*NTSTEP_BETWEEN_READ_ADJSRC),cudaMemcpyHostToDevice),99099);

  // launches cuda kernel for acoustic adjoint sources
  add_sources_ac_SIM_TYPE_2_OR_3_kernel<<<grid,threads,0,mp->compute_stream>>>(mp->d_potential_dot_dot_acoustic,
                                                                                *nrec,it_index,*NTSTEP_BETWEEN_READ_ADJSRC,
                                                                                mp->d_source_adjoint,
                                                                                mp->d_hxir,
                                                                                mp->d_hetar,
                                                                                mp->d_hgammar,
                                                                                mp->d_ibool,
                                                                                mp->d_ispec_is_acoustic,
                                                                                mp->d_ispec_selected_rec_loc,
                                                                                mp->nadj_rec_local,
                                                                                mp->d_kappastore);


#ifdef ENABLE_VERY_SLOW_ERROR_CHECKING
  exit_on_cuda_error("add_sources_acoustic_SIM_TYPE_2_OR_3_kernel");
#endif
}
