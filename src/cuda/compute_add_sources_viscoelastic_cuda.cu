/*
 !=====================================================================
 !
 !               S p e c f e m 3 D  V e r s i o n  3 . 0
 !               ---------------------------------------
 !
 !     Main historical authors: Dimitri Komatitsch and Jeroen Tromp
 !                        Princeton University, USA
 !                and CNRS / University of Marseille, France
 !                 (there are currently many more authors!)
 ! (c) Princeton University and CNRS / University of Marseille, July 2012
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

#include "mesh_constants_cuda.h"

/* ----------------------------------------------------------------------------------------------- */

// elastic domain sources

/* ----------------------------------------------------------------------------------------------- */

__global__ void compute_add_sources_kernel(realw* accel,
                                           int* d_ibool,
                                           realw* sourcearrays,
                                           double* stf_pre_compute,
                                           int myrank,
                                           int* islice_selected_source,
                                           int* ispec_selected_source,
                                           int* ispec_is_elastic,
                                           int NSOURCES) {
  int i = threadIdx.x;
  int j = threadIdx.y;
  int k = threadIdx.z;

  int isource  = blockIdx.x + gridDim.x*blockIdx.y; // bx

  int ispec,iglob;
  realw stf;

  if (isource < NSOURCES) { // when NSOURCES > 65535, but mod(nspec_top,2) > 0, we end up with an extra block.

    if (myrank == islice_selected_source[isource]) {

      ispec = ispec_selected_source[isource]-1;

      if (ispec_is_elastic[ispec]) {

        stf = (realw) stf_pre_compute[isource];
        iglob = d_ibool[INDEX4_PADDED(NGLLX,NGLLX,NGLLX,i,j,k,ispec)]-1;

        atomicAdd(&accel[iglob*3],sourcearrays[INDEX5(NSOURCES,NDIM,NGLLX,NGLLX,isource, 0,i,j,k)]*stf);
        atomicAdd(&accel[iglob*3+1],sourcearrays[INDEX5(NSOURCES,NDIM,NGLLX,NGLLX,isource, 1,i,j,k)]*stf);
        atomicAdd(&accel[iglob*3+2],sourcearrays[INDEX5(NSOURCES,NDIM,NGLLX,NGLLX,isource, 2,i,j,k)]*stf);
      }
    }
  }

}


/* ----------------------------------------------------------------------------------------------- */

extern "C"
void FC_FUNC_(compute_add_sources_el_cuda,
              COMPUTE_ADD_SOURCES_EL_CUDA)(long* Mesh_pointer,
                                           double* h_stf_pre_compute,
                                           int* h_NSOURCES) {

  TRACE("\tcompute_add_sources_el_cuda");

  Mesh* mp = (Mesh*)(*Mesh_pointer); //get mesh pointer out of fortran integer container

  // check if anything to do
  if (mp->nsources_local == 0) return;

  int NSOURCES = *h_NSOURCES;

  print_CUDA_error_if_any(cudaMemcpy(mp->d_stf_pre_compute,h_stf_pre_compute,
                                     NSOURCES*sizeof(double),cudaMemcpyHostToDevice),18);

#ifdef ENABLE_VERY_SLOW_ERROR_CHECKING
  exit_on_cuda_error("compute_add_sources_el_cuda copy");
#endif

  int num_blocks_x, num_blocks_y;
  get_blocks_xy(NSOURCES,&num_blocks_x,&num_blocks_y);

  dim3 grid(num_blocks_x,num_blocks_y);
  dim3 threads(5,5,5);

  compute_add_sources_kernel<<<grid,threads,0,mp->compute_stream>>>(mp->d_accel,mp->d_ibool,
                                                                    mp->d_sourcearrays,
                                                                    mp->d_stf_pre_compute,
                                                                    mp->myrank,
                                                                    mp->d_islice_selected_source,mp->d_ispec_selected_source,
                                                                    mp->d_ispec_is_elastic,
                                                                    NSOURCES);

#ifdef ENABLE_VERY_SLOW_ERROR_CHECKING
  exit_on_cuda_error("compute_add_sources_el_cuda");
#endif
}

/* ----------------------------------------------------------------------------------------------- */

extern "C"
void FC_FUNC_(compute_add_sources_el_s3_cuda,
              COMPUTE_ADD_SOURCES_EL_S3_CUDA)(long* Mesh_pointer,
                                              double* h_stf_pre_compute,
                                              int* h_NSOURCES) {

  TRACE("\tcompute_add_sources_el_s3_cuda");
  // EPIK_TRACER("compute_add_sources_el_s3_cuda");

  Mesh* mp = (Mesh*)(*Mesh_pointer); //get mesh pointer out of fortran integer container

  int NSOURCES = *h_NSOURCES;

  print_CUDA_error_if_any(cudaMemcpy(mp->d_stf_pre_compute,h_stf_pre_compute,
                                     NSOURCES*sizeof(double),cudaMemcpyHostToDevice),18);

#ifdef ENABLE_VERY_SLOW_ERROR_CHECKING
  exit_on_cuda_error("compute_add_sources_el_s3_cuda copy");
#endif

  int num_blocks_x, num_blocks_y;
  get_blocks_xy(NSOURCES,&num_blocks_x,&num_blocks_y);

  dim3 grid(num_blocks_x,num_blocks_y);
  dim3 threads(5,5,5);

  compute_add_sources_kernel<<<grid,threads,0,mp->compute_stream>>>(mp->d_b_accel,mp->d_ibool,
                                                                    mp->d_sourcearrays,
                                                                    mp->d_stf_pre_compute,
                                                                    mp->myrank,
                                                                    mp->d_islice_selected_source,mp->d_ispec_selected_source,
                                                                    mp->d_ispec_is_elastic,
                                                                    NSOURCES);

#ifdef ENABLE_VERY_SLOW_ERROR_CHECKING
  exit_on_cuda_error("compute_add_sources_el_s3_cuda");
#endif
}

/* ----------------------------------------------------------------------------------------------- */

// NOISE sources

/* ----------------------------------------------------------------------------------------------- */

__global__ void add_source_master_rec_noise_cuda_kernel(int* d_ibool,
                                                        int* ispec_selected_rec,
                                                        int irec_master_noise,
                                                        realw* accel,
                                                        realw* noise_sourcearray,
                                                        int it) {
  int tx = threadIdx.x;
  int iglob = d_ibool[tx + NGLL3_PADDED*(ispec_selected_rec[irec_master_noise-1]-1)]-1;

  // not sure if we need atomic operations but just in case...
  // accel[3*iglob] += noise_sourcearray[3*tx + 3*125*it];
  // accel[1+3*iglob] += noise_sourcearray[1+3*tx + 3*125*it];
  // accel[2+3*iglob] += noise_sourcearray[2+3*tx + 3*125*it];

  atomicAdd(&accel[iglob*3],noise_sourcearray[3*tx + 3*NGLL3*it]);
  atomicAdd(&accel[iglob*3+1],noise_sourcearray[1+3*tx + 3*NGLL3*it]);
  atomicAdd(&accel[iglob*3+2],noise_sourcearray[2+3*tx + 3*NGLL3*it]);

}

/* ----------------------------------------------------------------------------------------------- */

extern "C"
void FC_FUNC_(add_source_master_rec_noise_cu,
              ADD_SOURCE_MASTER_REC_NOISE_CU)(long* Mesh_pointer,
                                              int* it_f,
                                              int* irec_master_noise_f,
                                              int* islice_selected_rec) {

TRACE("\tadd_source_master_rec_noise_cu");

  Mesh* mp = (Mesh*)(*Mesh_pointer); //get mesh pointer out of fortran integer container

  int it = *it_f-1; // -1 for Fortran -> C indexing differences
  int irec_master_noise = *irec_master_noise_f;

  dim3 grid(1,1,1);
  dim3 threads(NGLL3,1,1);

  if (mp->myrank == islice_selected_rec[irec_master_noise-1]) {
    add_source_master_rec_noise_cuda_kernel<<<grid,threads,0,mp->compute_stream>>>(mp->d_ibool,
                                                                                    mp->d_ispec_selected_rec,
                                                                                    irec_master_noise,
                                                                                    mp->d_accel,
                                                                                    mp->d_noise_sourcearray,
                                                                                    it);

#ifdef ENABLE_VERY_SLOW_ERROR_CHECKING
  exit_on_cuda_error("add_source_master_rec_noise_cuda_kernel");
#endif
  }
}

/* ----------------------------------------------------------------------------------------------- */

// ADJOINT sources

/* ----------------------------------------------------------------------------------------------- */

__global__ void add_sources_el_SIM_TYPE_2_OR_3_kernel(realw* accel,
                                                     int nrec,
                                                     realw* adj_sourcearrays,
                                                     int* d_ibool,
                                                     int* ispec_is_elastic,
                                                     int* ispec_selected_rec,
                                                     int* pre_computed_irec,
                                                     int nadj_rec_local) {

  int irec_local = blockIdx.x + gridDim.x*blockIdx.y;

  if (irec_local < nadj_rec_local) { // when nrec > 65535, but mod(nspec_top,2) > 0, we end up with an extra block.

    int irec = pre_computed_irec[irec_local];

    int ispec = ispec_selected_rec[irec]-1;

    if (ispec_is_elastic[ispec]){
      int i = threadIdx.x;
      int j = threadIdx.y;
      int k = threadIdx.z;
      int iglob = d_ibool[INDEX4_PADDED(NGLLX,NGLLX,NGLLX,i,j,k,ispec)]-1;

      // atomic operations are absolutely necessary for correctness!
      atomicAdd(&accel[3*iglob],adj_sourcearrays[INDEX5(NGLLX,NGLLX,NGLLX,NDIM,i,j,k,0,irec_local)]);
      atomicAdd(&accel[1+3*iglob], adj_sourcearrays[INDEX5(NGLLX,NGLLX,NGLLX,NDIM,i,j,k,1,irec_local)]);
      atomicAdd(&accel[2+3*iglob],adj_sourcearrays[INDEX5(NGLLX,NGLLX,NGLLX,NDIM,i,j,k,2,irec_local)]);
    } // ispec_is_elastic
  }

}

/* ----------------------------------------------------------------------------------------------- */

extern "C"
void FC_FUNC_(add_sources_el_sim_type_2_or_3,
              ADD_SOURCES_EL_SIM_TYPE_2_OR_3)(long* Mesh_pointer,
                                               realw* h_adj_sourcearrays,
                                               int* h_ispec_is_elastic,
                                               int* h_ispec_selected_rec,
                                               int* nrec,
                                               int* time_index,
                                               int* h_islice_selected_rec,
                                               int* nadj_rec_local,
                                               int* NTSTEP_BETWEEN_READ_ADJSRC) {

  TRACE("\tadd_sources_el_sim_type_2_or_3");

  Mesh* mp = (Mesh*)(*Mesh_pointer); //get mesh pointer out of fortran integer container

  // checks
  if (*nadj_rec_local != mp->nadj_rec_local) exit_on_error("add_sources_el_sim_type_2_or_3: nadj_rec_local not equal\n");

  int num_blocks_x, num_blocks_y;
  get_blocks_xy(mp->nadj_rec_local,&num_blocks_x,&num_blocks_y);

  dim3 grid(num_blocks_x,num_blocks_y,1);
  dim3 threads(5,5,5);

  // build slice of adj_sourcearrays because full array is *very* large.
  // note: this extracts array values for local adjoint sources at given time step "time_index"
  //          from large adj_sourcearrays array into h_adj_sourcearrays_slice
  int ispec,i,j,k,irec_local,it_index;

  it_index = (*time_index) - 1;
  irec_local = 0;

  for(int irec = 0; irec < *nrec; irec++) {
    if (mp->myrank == h_islice_selected_rec[irec]) {
      // takes only elastic sources
      ispec = h_ispec_selected_rec[irec] - 1;

      // only for elastic elements
      if (h_ispec_is_elastic[ispec]){

        for(k=0;k<NGLLX;k++) {
          for(j=0;j<NGLLX;j++) {
            for(i=0;i<NGLLX;i++) {

              mp->h_adj_sourcearrays_slice[INDEX5(NGLLX,NGLLX,NGLLX,NDIM,i,j,k,0,irec_local)]
                      = h_adj_sourcearrays[INDEX6(*nadj_rec_local,
                                                  *NTSTEP_BETWEEN_READ_ADJSRC,NDIM,NGLLX,NGLLX,
                                                  irec_local,it_index,0,i,j,k)];

              mp->h_adj_sourcearrays_slice[INDEX5(NGLLX,NGLLX,NGLLX,NDIM,i,j,k,1,irec_local)]
                      = h_adj_sourcearrays[INDEX6(*nadj_rec_local,
                                                  *NTSTEP_BETWEEN_READ_ADJSRC,NDIM,NGLLX,NGLLX,
                                                  irec_local,it_index,1,i,j,k)];

              mp->h_adj_sourcearrays_slice[INDEX5(NGLLX,NGLLX,NGLLX,NDIM,i,j,k,2,irec_local)]
                      = h_adj_sourcearrays[INDEX6(*nadj_rec_local,
                                                  *NTSTEP_BETWEEN_READ_ADJSRC,NDIM,NGLLX,NGLLX,
                                                  irec_local,it_index,2,i,j,k)];
            }
          }
        }
      } // h_ispec_is_elastic

      // increases local receivers counter
      irec_local++;
    }
  }
  // check all local sources were added
  if (irec_local != mp->nadj_rec_local) exit_on_error("irec_local not equal to nadj_rec_local\n");

  // copies extracted array values onto GPU
  print_CUDA_error_if_any(cudaMemcpy(mp->d_adj_sourcearrays, mp->h_adj_sourcearrays_slice,
                                    (mp->nadj_rec_local)*3*NGLL3*sizeof(realw),cudaMemcpyHostToDevice),98001);


  // the irec_local variable needs to be precomputed (as
  // h_pre_comp..), because normally it is in the loop updating accel,
  // and due to how it's incremented, it cannot be parallelized

  add_sources_el_SIM_TYPE_2_OR_3_kernel<<<grid,threads,0,mp->compute_stream>>>(mp->d_accel,
                                                                               *nrec,
                                                                               mp->d_adj_sourcearrays,
                                                                               mp->d_ibool,
                                                                               mp->d_ispec_is_elastic,
                                                                               mp->d_ispec_selected_rec,
                                                                               mp->d_pre_computed_irec,
                                                                               mp->nadj_rec_local);

#ifdef ENABLE_VERY_SLOW_ERROR_CHECKING
  exit_on_cuda_error("add_sources_SIM_TYPE_2_OR_3_kernel");
#endif
}

