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

// elastic domain sources

/* ----------------------------------------------------------------------------------------------- */

__global__ void compute_add_sources_kernel(realw* accel,
                                           int* d_ibool,
                                           realw* sourcearrays,
                                           field* stf_pre_compute,
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
  field stf;

  if (isource < NSOURCES) { // when NSOURCES > 65535, but mod(nspec_top,2) > 0, we end up with an extra block.

    if (myrank == islice_selected_source[isource]) {

      ispec = ispec_selected_source[isource]-1;

      if (ispec_is_elastic[ispec]) {

        stf = stf_pre_compute[isource];
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

  // convert to GPU precision
  realw* stf_pre_compute;
  stf_pre_compute = (realw*)malloc(NSOURCES * sizeof(realw));
  for (int i_source=0;i_source < NSOURCES;i_source++) stf_pre_compute[i_source] = (realw)h_stf_pre_compute[i_source];

  print_CUDA_error_if_any(cudaMemcpy(mp->d_stf_pre_compute,stf_pre_compute,
                                     NSOURCES*sizeof(realw),cudaMemcpyHostToDevice),18);
  free(stf_pre_compute);

#ifdef ENABLE_VERY_SLOW_ERROR_CHECKING
  exit_on_cuda_error("compute_add_sources_el_cuda copy");
#endif

  int num_blocks_x, num_blocks_y;
  get_blocks_xy(NSOURCES,&num_blocks_x,&num_blocks_y);

  dim3 grid(num_blocks_x,num_blocks_y);
  dim3 threads(NGLLX,NGLLY,NGLLZ);

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

  // convert to GPU precision
  realw* stf_pre_compute;
  stf_pre_compute = (realw*)malloc(NSOURCES * sizeof(realw));
  for (int i_source=0;i_source < NSOURCES;i_source++) stf_pre_compute[i_source] = (realw)h_stf_pre_compute[i_source];

  print_CUDA_error_if_any(cudaMemcpy(mp->d_stf_pre_compute,stf_pre_compute,
                                     NSOURCES*sizeof(realw),cudaMemcpyHostToDevice),18);
  free(stf_pre_compute);

#ifdef ENABLE_VERY_SLOW_ERROR_CHECKING
  exit_on_cuda_error("compute_add_sources_el_s3_cuda copy");
#endif

  int num_blocks_x, num_blocks_y;
  get_blocks_xy(NSOURCES,&num_blocks_x,&num_blocks_y);

  dim3 grid(num_blocks_x,num_blocks_y);
  dim3 threads(NGLLX,NGLLY,NGLLZ);

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
                                                      int it,
                                                      int NSTEP_BETWEEN_ADJSRC,
                                                      field* source_adjoint,
                                                      realw* xir_store,
                                                      realw* etar_store,
                                                      realw* gammar_store,
                                                      int* d_ibool,
                                                      int* ispec_is_elastic,
                                                      int* ispec_selected_recloc,
                                                      int nadj_rec_local) {

  int irec_local = blockIdx.x + gridDim.x*blockIdx.y;

  if (irec_local < nadj_rec_local) { // when nrec > 65535, but mod(nspec_top,2) > 0, we end up with an extra block.

    int ispec = ispec_selected_recloc[irec_local]-1;

    if (ispec_is_elastic[ispec]){
      int i = threadIdx.x;
      int j = threadIdx.y;
      int k = threadIdx.z;
      int iglob = d_ibool[INDEX4_PADDED(NGLLX,NGLLX,NGLLX,i,j,k,ispec)]-1;
      realw xir    = xir_store[INDEX2(nadj_rec_local,irec_local,i)];
      realw etar   = etar_store[INDEX2(nadj_rec_local,irec_local,j)];
      realw gammar = gammar_store[INDEX2(nadj_rec_local,irec_local,k)];

      realw lagrange =   xir * etar * gammar ;

      // atomic operations are absolutely necessary for correctness!
      atomicAdd(&accel[3*iglob],source_adjoint[INDEX3(NDIM,nadj_rec_local,0,irec_local,it)]*lagrange);
      atomicAdd(&accel[1+3*iglob],source_adjoint[INDEX3(NDIM,nadj_rec_local,1,irec_local,it)]*lagrange);
      atomicAdd(&accel[2+3*iglob],source_adjoint[INDEX3(NDIM,nadj_rec_local,2,irec_local,it)]*lagrange);
    } // ispec_is_elastic
  }

}

/* ----------------------------------------------------------------------------------------------- */

extern "C"
void FC_FUNC_(add_sources_el_sim_type_2_or_3,
              ADD_SOURCES_EL_SIM_TYPE_2_OR_3)(long* Mesh_pointer,
                                              realw* h_source_adjoint,
                                              int* nrec,
                                              int* nadj_rec_local,
                                              int* NTSTEP_BETWEEN_READ_ADJSRC,
                                              int* it) {

  TRACE("\tadd_sources_el_sim_type_2_or_3");

  Mesh* mp = (Mesh*)(*Mesh_pointer); //get mesh pointer out of fortran integer container

  // checks
  if (*nadj_rec_local != mp->nadj_rec_local) exit_on_error("add_sources_el_sim_type_2_or_3: nadj_rec_local not equal\n");

  int num_blocks_x, num_blocks_y;
  get_blocks_xy(mp->nadj_rec_local,&num_blocks_x,&num_blocks_y);

  dim3 grid(num_blocks_x,num_blocks_y,1);
  dim3 threads(NGLLX,NGLLY,NGLLZ);

  int it_index = *NTSTEP_BETWEEN_READ_ADJSRC - (*it-1) % *NTSTEP_BETWEEN_READ_ADJSRC - 1 ;
  // copies extracted array values onto GPU
  if ( (*it-1) % *NTSTEP_BETWEEN_READ_ADJSRC==0) print_CUDA_error_if_any(cudaMemcpy(mp->d_source_adjoint,h_source_adjoint,
                                                                mp->nadj_rec_local*3*sizeof(realw)*(*NTSTEP_BETWEEN_READ_ADJSRC),cudaMemcpyHostToDevice),99099);




  add_sources_el_SIM_TYPE_2_OR_3_kernel<<<grid,threads,0,mp->compute_stream>>>(mp->d_accel,
                                                                               *nrec,it_index,*NTSTEP_BETWEEN_READ_ADJSRC,
                                                                               mp->d_source_adjoint,
                                                                               mp->d_hxir,
                                                                               mp->d_hetar,
                                                                               mp->d_hgammar,
                                                                               mp->d_ibool,
                                                                               mp->d_ispec_is_elastic,
                                                                               mp->d_ispec_selected_rec_loc,
                                                                               mp->nadj_rec_local);

#ifdef ENABLE_VERY_SLOW_ERROR_CHECKING
  exit_on_cuda_error("add_sources_SIM_TYPE_2_OR_3_kernel");
#endif
}

