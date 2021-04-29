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


extern EXTERN_LANG
void FC_FUNC_(compute_seismograms_cuda,
              COMPUTE_SEISMOGRAMS_CUDA)(long* Mesh_pointer_f,
                                        realw* seismograms_d,
                                        realw* seismograms_v,
                                        realw* seismograms_a,
                                        realw* seismograms_p,
                                        int* seismo_currentf,
                                        int* nlength_seismogramf,
                                        int* itf, int* it_endf,
                                        int* ACOUSTIC_SIMULATION,
                                        int* ELASTIC_SIMULATION,
                                        int* USE_TRICK_FOR_BETTER_PRESSURE) {

// compute_seismograms
  TRACE("compute_seismograms_cuda");

  Mesh* mp = (Mesh*)(*Mesh_pointer_f); // get Mesh from fortran integer wrapper

  //checks if anything to do
  if (mp->nrec_local == 0) return;

  int num_blocks_x, num_blocks_y;
  get_blocks_xy(mp->nrec_local,&num_blocks_x,&num_blocks_y);

  dim3 grid(num_blocks_x,num_blocks_y);
  dim3 threads(NGLL3_PADDED,1,1);

  int seismo_current = *seismo_currentf - 1;  // starts indexing from 0 for cuda arrays
  int nlength_seismogram = *nlength_seismogramf;

  int it = *itf;
  int it_end = *it_endf;

  // sets seismogram types array to loop over (to treat the calls below similar to 2D version)
  int h_seismo[4];
  int nseismos = 0;
  if (mp->save_seismograms_d){ h_seismo[nseismos] = 1; nseismos += 1; }
  if (mp->save_seismograms_v){ h_seismo[nseismos] = 2; nseismos += 1; }
  if (mp->save_seismograms_a){ h_seismo[nseismos] = 3; nseismos += 1; }
  if (mp->save_seismograms_p){ h_seismo[nseismos] = 4; nseismos += 1; }

  // note: mp->d_ispec_selected_rec_loc is the array holding spectral elements in which the local receivers are located
  //       for "pure" adjoint simulation (SIMULATION_TYPE == 2), adjoint "receivers" are located at CMT source positions,
  //       otherwise receivers are located at station positions.
  //       the array mp->d_ispec_selected_rec_loc is setup accordingly in prepare_constants_device() routine.

  // computes seismograms
  for(int i = 0; i < nseismos; i++){
    int seismotype = h_seismo[i];

    // skip if zero
    if (seismotype == 0){ continue; }

    // selects wavefields (see corresponding handling in compute_seismograms.f90)
    realw* displ = NULL;
    field* potential = NULL;
    realw* d_seismo = NULL;
    field* d_seismo_p = NULL;

    if (seismotype == 1){
      // deplacement
      if (mp->simulation_type == 1 || mp->simulation_type == 2){
        displ = mp->d_displ;
        potential = mp->d_potential_acoustic;
      }else{
        // kernel simulations
        // reconstructed forward wavefield stored in b_displ, b_veloc, b_accel
        displ = mp->d_b_displ;
        potential = mp->d_b_potential_acoustic;
      }
      d_seismo = mp->d_seismograms_d;

    }else if (seismotype == 2){
      // vitesse
      if (mp->simulation_type == 1 || mp->simulation_type == 2){
        displ = mp->d_veloc;
        potential = mp->d_potential_dot_acoustic;
      }else{
        // kernel simulations
        displ = mp->d_b_veloc;
        potential = mp->d_b_potential_dot_acoustic;
      }
      d_seismo = mp->d_seismograms_v;

    }else if (seismotype == 3){
      // acceleration
      if (mp->simulation_type == 1 || mp->simulation_type == 2){
        displ = mp->d_accel;
        potential = mp->d_potential_dot_dot_acoustic;
      }else{
        // kernel simulations
        displ = mp->d_b_accel;
        potential = mp->d_b_potential_dot_dot_acoustic;
      }
      d_seismo = mp->d_seismograms_a;

    }else if (seismotype == 4){
      // pression
      if (mp->simulation_type == 1 || mp->simulation_type == 2){
        displ = mp->d_displ;
        if (*USE_TRICK_FOR_BETTER_PRESSURE){
          potential = mp->d_potential_acoustic;
        }else{
          potential = mp->d_potential_dot_dot_acoustic;
        }
      }else{
        // kernel simulations
        displ = mp->d_b_displ;
        if (*USE_TRICK_FOR_BETTER_PRESSURE){
          potential = mp->d_b_potential_acoustic;
        }else{
          potential = mp->d_b_potential_dot_dot_acoustic;
        }
      }
      d_seismo_p = mp->d_seismograms_p;
    }

    // computes current seismograms value
    switch (seismotype){
    case 1 :
    case 2 :
    case 3 :
      //Displ/Veloc/Accel
#ifdef USE_CUDA
      if (run_cuda){
        compute_elastic_seismogram_kernel<<<grid,threads,0,mp->compute_stream>>>(mp->nrec_local,
                                                                                 displ,
                                                                                 potential,
                                                                                 mp->d_ibool,
                                                                                 mp->d_hxir,mp->d_hetar,mp->d_hgammar,
                                                                                 d_seismo,
                                                                                 mp->d_nu_rec,
                                                                                 mp->d_ispec_selected_rec_loc,
                                                                                 mp->d_ispec_is_elastic,
                                                                                 mp->d_ispec_is_acoustic,
                                                                                 mp->d_rhostore,
                                                                                 mp->d_hprime_xx,
                                                                                 mp->d_xix,mp->d_xiy,mp->d_xiz,
                                                                                 mp->d_etax,mp->d_etay,mp->d_etaz,
                                                                                 mp->d_gammax,mp->d_gammay,mp->d_gammaz,
                                                                                 mp->d_irregular_element_number,
                                                                                 mp->xix_regular,
                                                                                 seismo_current);
      }
#endif
#ifdef USE_HIP
      if (run_hip){
        hipLaunchKernelGGL(compute_elastic_seismogram_kernel, dim3(grid), dim3(threads), 0, mp->compute_stream,
                                                                 mp->nrec_local,
                                                                 displ,
                                                                 potential,
                                                                 mp->d_ibool,
                                                                 mp->d_hxir,mp->d_hetar,mp->d_hgammar,
                                                                 d_seismo,
                                                                 mp->d_nu_rec,
                                                                 mp->d_ispec_selected_rec_loc,
                                                                 mp->d_ispec_is_elastic,
                                                                 mp->d_ispec_is_acoustic,
                                                                 mp->d_rhostore,
                                                                 mp->d_hprime_xx,
                                                                 mp->d_xix,mp->d_xiy,mp->d_xiz,
                                                                 mp->d_etax,mp->d_etay,mp->d_etaz,
                                                                 mp->d_gammax,mp->d_gammay,mp->d_gammaz,
                                                                 mp->d_irregular_element_number,
                                                                 mp->xix_regular,
                                                                 seismo_current);
      }
#endif
      break;

    case 4 :
      //Pression
#ifdef USE_CUDA
      if (run_cuda){
        compute_acoustic_seismogram_kernel<<<grid,threads,0,mp->compute_stream>>>(mp->nrec_local,
                                                                                  displ,
                                                                                  potential,
                                                                                  mp->d_ibool,
                                                                                  mp->d_hxir,mp->d_hetar,mp->d_hgammar,
                                                                                  d_seismo_p,
                                                                                  mp->d_ispec_selected_rec_loc,
                                                                                  mp->d_ispec_is_elastic,
                                                                                  mp->d_ispec_is_acoustic,
                                                                                  mp->d_kappav,mp->d_muv,
                                                                                  mp->d_hprime_xx,
                                                                                  mp->d_xix,mp->d_xiy,mp->d_xiz,
                                                                                  mp->d_etax,mp->d_etay,mp->d_etaz,
                                                                                  mp->d_gammax,mp->d_gammay,mp->d_gammaz,
                                                                                  mp->d_irregular_element_number,
                                                                                  mp->xix_regular,
                                                                                  mp->ANISOTROPY,
                                                                                  mp->d_c11store,mp->d_c12store,mp->d_c13store,
                                                                                  mp->d_c14store,mp->d_c15store,mp->d_c16store,
                                                                                  mp->d_c22store,mp->d_c23store,mp->d_c24store,
                                                                                  mp->d_c25store,mp->d_c26store,mp->d_c33store,
                                                                                  mp->d_c34store,mp->d_c35store,mp->d_c36store,
                                                                                  seismo_current);
      }
#endif
#ifdef USE_HIP
      if (run_hip){
        hipLaunchKernelGGL(compute_acoustic_seismogram_kernel, dim3(grid), dim3(threads), 0, mp->compute_stream,
                                                               mp->nrec_local,
                                                               displ,
                                                               potential,
                                                               mp->d_ibool,
                                                               mp->d_hxir,mp->d_hetar,mp->d_hgammar,
                                                               d_seismo_p,
                                                               mp->d_ispec_selected_rec_loc,
                                                               mp->d_ispec_is_elastic,
                                                               mp->d_ispec_is_acoustic,
                                                               mp->d_kappav,mp->d_muv,
                                                               mp->d_hprime_xx,
                                                               mp->d_xix,mp->d_xiy,mp->d_xiz,
                                                               mp->d_etax,mp->d_etay,mp->d_etaz,
                                                               mp->d_gammax,mp->d_gammay,mp->d_gammaz,
                                                               mp->d_irregular_element_number,
                                                               mp->xix_regular,
                                                               mp->ANISOTROPY,
                                                               mp->d_c11store,mp->d_c12store,mp->d_c13store,
                                                               mp->d_c14store,mp->d_c15store,mp->d_c16store,
                                                               mp->d_c22store,mp->d_c23store,mp->d_c24store,
                                                               mp->d_c25store,mp->d_c26store,mp->d_c33store,
                                                               mp->d_c34store,mp->d_c35store,mp->d_c36store,
                                                               seismo_current);
      }
#endif
      break;
    } //switch
  }//for

  // note: due to subsampling, the last time step it == it_end might not be reached,
  //       but computing seismogram entries might end before.
  //       thus, both checks
  //         it%NTSTEP_BETWEEN_OUTPUT_SEISMOS == 0 || it == it_end
  //       might not be reached. instead we test if the seismogram array is full by
  //         seismo_current == nlength_seismogram - 1
  //       and copy it back whenever.
  //printf("debug: gpu seismo: seismo current/lenght %i/%i - it/it_end %i/%i\n",seismo_current,nlength_seismogram,it,it_end);

  // copies array to CPU host
  if (seismo_current == nlength_seismogram - 1 || it == it_end){
    // explicitly waits until previous compute stream finishes
    // (cudaMemcpy implicitly synchronizes all other cuda operations)
    gpuStreamSynchronize(mp->compute_stream);

    int size = mp->nrec_local * nlength_seismogram;

    // (cudaMemcpy implicitly synchronizes all other cuda operations)
    if (mp->save_seismograms_d && ! *USE_TRICK_FOR_BETTER_PRESSURE)
      gpuMemcpy_tohost_realw(seismograms_d,mp->d_seismograms_d,NDIM * size);
    if (mp->save_seismograms_v && ! *USE_TRICK_FOR_BETTER_PRESSURE)
      gpuMemcpy_tohost_realw(seismograms_v,mp->d_seismograms_v,NDIM * size);
    if (mp->save_seismograms_a && ! *USE_TRICK_FOR_BETTER_PRESSURE)
      gpuMemcpy_tohost_realw(seismograms_a,mp->d_seismograms_a,NDIM * size);

    // EB EB Temporary solution : in the future we will also declare host pressure seismograms as
    //                            (1,nrec_local,NTSTEP_BETWEEN_OUTPUT_SEISMOS)
    if (mp->save_seismograms_p){
      // EB EB We need to reorganize data to match host array shape :
      //       if NB_RUNS_ACOUSTIC_GPU = 1: from fortran shape
      //          (1,nrec_local,NTSTEP_BETWEEN_OUTPUT_SEISMOS)
      //          to (NDIM,nrec_local,NTSTEP_BETWEEN_OUTPUT_SEISMOS)
      //       if NB_RUNS_ACOUSTIC_GPU > 1: from fortran shape
      //          (NB_RUNS_ACOUSTIC_GPU,nrec_local,NTSTEP_BETWEEN_OUTPUT_SEISMOS) to
      //          to (NDIM,nrec_local*NB_RUNS_ACOUSTIC_GPU,NTSTEP_BETWEEN_OUTPUT_SEISMOS)
      realw *seismo_temp = (realw*) malloc(size * NB_RUNS_ACOUSTIC_GPU * sizeof(realw));
      gpuMemcpy_tohost_void((void*)seismo_temp,(void*)mp->d_seismograms_p,size * NB_RUNS_ACOUSTIC_GPU * sizeof(realw));

      for (int j = 0; j < nlength_seismogram; j++){
        for (int i_recloc = 0; i_recloc < mp->nrec_local; i_recloc++){
          for (int i_run = 0; i_run < NB_RUNS_ACOUSTIC_GPU; i_run++){
            seismograms_p[INDEX4(NDIM,mp->nrec_local,NB_RUNS_ACOUSTIC_GPU,0,i_recloc,i_run,j)] =
                    seismo_temp[INDEX3(NB_RUNS_ACOUSTIC_GPU,mp->nrec_local,i_run,i_recloc,j)];
            seismograms_p[INDEX4(NDIM,mp->nrec_local,NB_RUNS_ACOUSTIC_GPU,1,i_recloc,i_run,j)] = 0.f;
            seismograms_p[INDEX4(NDIM,mp->nrec_local,NB_RUNS_ACOUSTIC_GPU,2,i_recloc,i_run,j)] = 0.f;
          }
        }
      }

      free(seismo_temp);

      // debug - checks min/max
      /*
      for (int i_recloc=0; i_recloc<mp->nrec_local; i_recloc++){
        for (int i_run=0; i_run<NB_RUNS_ACOUSTIC_GPU; i_run++){
          float xmin,xmax;
          xmin = 0.f;
          xmax = 0.f;
          for (int j = 0; j< nlength_seismogram; j++){
            int idx = INDEX4(NDIM,mp->nrec_local,NB_RUNS_ACOUSTIC_GPU,0,i_recloc,i_run,j);
            xmin = min(xmin,seismograms_p[idx]);
            xmax = max(xmax,seismograms_p[idx]);
          }
          printf("debug: gpu seismo: run %i receiver %i min/max = %f/%f\n",i_run,i_recloc,xmin,xmax);
        }
      }
      */
    }
  }

  GPU_ERROR_CHECKING("after compute_seismograms_cuda");
}

