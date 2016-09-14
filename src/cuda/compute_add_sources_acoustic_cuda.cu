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

// acoustic sources

/* ----------------------------------------------------------------------------------------------- */

__global__ void compute_add_sources_acoustic_kernel(realw* potential_dot_dot_acoustic,
                                                    int* d_ibool,
                                                    realw* sourcearrays,
                                                    double* stf_pre_compute,
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
  realw stf,kappal;

  if (isource < NSOURCES){

    if (myrank == islice_selected_source[isource]) {

      ispec = ispec_selected_source[isource]-1;

      if (ispec_is_acoustic[ispec]) {

        iglob = d_ibool[INDEX4_PADDED(NGLLX,NGLLX,NGLLX,i,j,k,ispec)] - 1;

        stf = (realw) stf_pre_compute[isource];
        kappal = kappastore[INDEX4(NGLLX,NGLLX,NGLLX,i,j,k,ispec)];

        atomicAdd(&potential_dot_dot_acoustic[iglob],
                  -sourcearrays[INDEX5(NSOURCES,NDIM,NGLLX,NGLLX,isource, 0,i,j,k)]*stf/kappal);

        // debug: without atomic operation
        //      potential_dot_dot_acoustic[iglob] +=
        //                -sourcearrays[INDEX5(NSOURCES, 3, 5, 5,isource, 0, i,j,k)]*stf/kappal;
      }
    }
  }
}


/* ----------------------------------------------------------------------------------------------- */

extern "C"
void FC_FUNC_(compute_add_sources_ac_cuda,
              COMPUTE_ADD_SOURCES_AC_CUDA)(long* Mesh_pointer,
                                           int* NSOURCESf,
                                           double* h_stf_pre_compute) {

  TRACE("compute_add_sources_ac_cuda");

  Mesh* mp = (Mesh*)(*Mesh_pointer); //get mesh pointer out of fortran integer container

  // check if anything to do
  if (mp->nsources_local == 0) return;

  int NSOURCES = *NSOURCESf;

  // copies pre-computed source time factors onto GPU
  print_CUDA_error_if_any(cudaMemcpy(mp->d_stf_pre_compute,h_stf_pre_compute,
                                     NSOURCES*sizeof(double),cudaMemcpyHostToDevice),18);

  int num_blocks_x, num_blocks_y;
  get_blocks_xy(NSOURCES,&num_blocks_x,&num_blocks_y);

  dim3 grid(num_blocks_x,num_blocks_y);
  dim3 threads(5,5,5);

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
                                              double* h_stf_pre_compute) {

  TRACE("compute_add_sources_ac_s3_cuda");

  Mesh* mp = (Mesh*)(*Mesh_pointer); //get mesh pointer out of fortran integer container

  // check if anything to do
  if (mp->nsources_local == 0) return;

  int NSOURCES = *NSOURCESf;

  // copies source time factors onto GPU
  print_CUDA_error_if_any(cudaMemcpy(mp->d_stf_pre_compute,h_stf_pre_compute,
                                     NSOURCES*sizeof(double),cudaMemcpyHostToDevice),18);

  int num_blocks_x, num_blocks_y;
  get_blocks_xy(NSOURCES,&num_blocks_x,&num_blocks_y);

  dim3 grid(num_blocks_x,num_blocks_y);
  dim3 threads(5,5,5);

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

__global__ void add_sources_ac_SIM_TYPE_2_OR_3_kernel(realw* potential_dot_dot_acoustic,
                                                      int nrec,
                                                      realw* adj_sourcearrays,
                                                      int* d_ibool,
                                                      int* ispec_is_acoustic,
                                                      int* ispec_selected_rec,
                                                      int* pre_computed_irec,
                                                      int nadj_rec_local,
                                                      realw* kappastore) {

  int irec_local = blockIdx.x + gridDim.x*blockIdx.y;

  // because of grid shape, irec_local can be too big
  if (irec_local < nadj_rec_local) {

    int irec = pre_computed_irec[irec_local];

    int ispec = ispec_selected_rec[irec]-1;
    if (ispec_is_acoustic[ispec]){
      int i = threadIdx.x;
      int j = threadIdx.y;
      int k = threadIdx.z;

      int iglob = d_ibool[INDEX4_PADDED(NGLLX,NGLLX,NGLLX,i,j,k,ispec)]-1;

      //kappal = kappastore[INDEX4(5,5,5,i,j,k,ispec)];

      //potential_dot_dot_acoustic[iglob] += adj_sourcearrays[INDEX6(nadj_rec_local,NTSTEP_BETWEEN_ADJSRC,3,5,5,
      //                                            pre_computed_irec_local_index[irec],
      //                                            pre_computed_index,
      //                                            0,
      //                                            i,j,k)]/kappal;

      // beware, for acoustic medium, a pressure source would be taking the negative
      // and divide by Kappa of the fluid;
      // this would have to be done when constructing the adjoint source.
      //
      // note: we take the first component of the adj_sourcearrays
      //          the idea is to have e.g. a pressure source, where all 3 components would be the same
      realw stf = adj_sourcearrays[INDEX5(NGLLX,NGLLX,NGLLX,NDIM,i,j,k,0,irec_local)]; // / kappal

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
                                               realw* h_adj_sourcearrays,
                                               int* h_ispec_is_acoustic,
                                               int* h_ispec_selected_rec,
                                               int* nrec,
                                               int* time_index,
                                               int* h_islice_selected_rec,
                                               int* nadj_rec_local,
                                               int* NTSTEP_BETWEEN_READ_ADJSRC) {

  TRACE("add_sources_ac_sim_2_or_3_cuda");

  Mesh* mp = (Mesh*)(*Mesh_pointer); //get mesh pointer out of fortran integer container

  // checks
  if (*nadj_rec_local != mp->nadj_rec_local) exit_on_cuda_error("add_sources_ac_sim_type_2_or_3: nadj_rec_local not equal\n");

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
      // takes only acoustic sources
      ispec = h_ispec_selected_rec[irec] - 1;

      // only for acoustic elements
      if (h_ispec_is_acoustic[ispec]){
        for(k=0;k<5;k++) {
          for(j=0;j<5;j++) {
            for(i=0;i<5;i++) {

              mp->h_adj_sourcearrays_slice[INDEX5(NGLLX,NGLLX,NGLLX,NDIM,i,j,k,0,irec_local)]
                = h_adj_sourcearrays[INDEX6(mp->nadj_rec_local,
                                          *NTSTEP_BETWEEN_READ_ADJSRC,NDIM,NGLLX,NGLLX,
                                          irec_local,it_index,0,i,j,k)];

              mp->h_adj_sourcearrays_slice[INDEX5(NGLLX,NGLLX,NGLLX,NDIM,i,j,k,1,irec_local)]
                = h_adj_sourcearrays[INDEX6(mp->nadj_rec_local,
                                          *NTSTEP_BETWEEN_READ_ADJSRC,NDIM,NGLLX,NGLLX,
                                          irec_local,it_index,1,i,j,k)];

              mp->h_adj_sourcearrays_slice[INDEX5(NGLLX,NGLLX,NGLLX,NDIM,i,j,k,2,irec_local)]
                = h_adj_sourcearrays[INDEX6(mp->nadj_rec_local,
                                          *NTSTEP_BETWEEN_READ_ADJSRC,NDIM,NGLLX,NGLLX,
                                          irec_local,it_index,2,i,j,k)];
            }
          }
        }
      } // h_ispec_is_acoustic

      // increases local receivers counter
      irec_local++;
    }
  }
  // check all local sources were added
  if (irec_local != mp->nadj_rec_local) exit_on_error("irec_local not equal to nadj_rec_local\n");

  // copies extracted array values onto GPU
  print_CUDA_error_if_any(cudaMemcpy(mp->d_adj_sourcearrays, mp->h_adj_sourcearrays_slice,
                                    (mp->nadj_rec_local)*3*NGLL3*sizeof(realw),cudaMemcpyHostToDevice),99099);

  // launches cuda kernel for acoustic adjoint sources
  add_sources_ac_SIM_TYPE_2_OR_3_kernel<<<grid,threads,0,mp->compute_stream>>>(mp->d_potential_dot_dot_acoustic,
                                                                                *nrec,
                                                                                mp->d_adj_sourcearrays,
                                                                                mp->d_ibool,
                                                                                mp->d_ispec_is_acoustic,
                                                                                mp->d_ispec_selected_rec,
                                                                                mp->d_pre_computed_irec,
                                                                                mp->nadj_rec_local,
                                                                                mp->d_kappastore);

#ifdef ENABLE_VERY_SLOW_ERROR_CHECKING
  exit_on_cuda_error("add_sources_acoustic_SIM_TYPE_2_OR_3_kernel");
#endif
}
