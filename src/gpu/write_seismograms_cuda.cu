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

  // selects wavefields (see corresponding handling in compute_seismograms.f90)
  realw* displ, *veloc, *accel;
  field* potential_acoustic, *potential_dot_acoustic, *potential_dot_dot_acoustic;
  if (mp->simulation_type == 1 || mp->simulation_type == 2){
    // forward simulations & pure adjoint simulations
    // wavefields stored in displ,veloc,accel
    displ = mp->d_displ;
    veloc = mp->d_veloc;
    accel = mp->d_accel;
    potential_acoustic = mp->d_potential_acoustic;
    potential_dot_acoustic = mp->d_potential_dot_acoustic;
    potential_dot_dot_acoustic = mp->d_potential_dot_dot_acoustic;
  }else{
    // kernel simulations
    // reconstructed forward wavefield stored in b_displ, b_veloc, b_accel
    displ = mp->d_b_displ;
    veloc = mp->d_b_veloc;
    accel = mp->d_b_accel;
    potential_acoustic = mp->d_b_potential_acoustic;
    potential_dot_acoustic = mp->d_b_potential_dot_acoustic;
    potential_dot_dot_acoustic = mp->d_b_potential_dot_dot_acoustic;
  }

  // note: mp->d_ispec_selected_rec_loc is the array holding spectral elements in which the local receivers are located
  //       for "pure" adjoint simulation (SIMULATION_TYPE == 2), adjoint "receivers" are located at CMT source positions,
  //       otherwise receivers are located at station positions.
  //       the array mp->d_ispec_selected_rec_loc is setup accordingly in prepare_constants_device() routine.

  // warning: put in fortran routine prepare_GPU()
  /*
  if (it == 0){
    if (mp->save_seismograms_d || mp->save_seismograms_v || mp->save_seismograms_a){
      // warnings
      if (! *ELASTIC_SIMULATION)
        printf("\nWarning: Wrong type of seismogram for a pure fluid simulation, use pressure in seismotype\n");
      if (*ELASTIC_SIMULATION && *ACOUSTIC_SIMULATION)
        printf("\nWarning: Coupled elastic/fluid simulation has only valid displacement seismograms in elastic domain for GPU simulation\n\n");
    }
    if (mp->save_seismograms_p){
      if (! *ACOUSTIC_SIMULATION)
        printf("\nWarning: Wrong type of seismogram for a pure elastic simulation, use displ veloc or accel in seismotype\n");
      if (*ELASTIC_SIMULATION && *ACOUSTIC_SIMULATION)
        printf("\nWarning: Coupled elastic/fluid simulation has only valid pressure seismograms in fluid domain for GPU simulation\n\n");
    }
  }
  */

  // todo: for coupled simulations, one should check in which domain the receiver lies to output displacement
  //       similar to what routine compute_vector_one_element(..) is doing

  // computes current seismograms value

  // elastic wavefield
  // acoustic wavefield
  if (*ELASTIC_SIMULATION){
    if (mp->save_seismograms_d){
#ifdef USE_CUDA
      if (run_cuda){
        compute_elastic_seismogram_kernel<<<grid,threads,0,mp->compute_stream>>>(mp->nrec_local,
                                                                                 displ,
                                                                                 mp->d_ibool,
                                                                                 mp->d_hxir,mp->d_hetar,mp->d_hgammar,
                                                                                 mp->d_seismograms_d,
                                                                                 mp->d_nu_rec,
                                                                                 mp->d_ispec_selected_rec_loc,
                                                                                 seismo_current);
      }
#endif
#ifdef USE_HIP
      if (run_hip){
        hipLaunchKernelGGL(compute_elastic_seismogram_kernel, dim3(grid), dim3(threads), 0, mp->compute_stream,
                                                                 mp->nrec_local,
                                                                 displ,
                                                                 mp->d_ibool,
                                                                 mp->d_hxir,mp->d_hetar,mp->d_hgammar,
                                                                 mp->d_seismograms_d,
                                                                 mp->d_nu_rec,
                                                                 mp->d_ispec_selected_rec_loc,
                                                                 seismo_current);
      }
#endif
    }

    if (mp->save_seismograms_v){
#ifdef USE_CUDA
      if (run_cuda){
        compute_elastic_seismogram_kernel<<<grid,threads,0,mp->compute_stream>>>(mp->nrec_local,
                                                                                 veloc,
                                                                                 mp->d_ibool,
                                                                                 mp->d_hxir,mp->d_hetar,mp->d_hgammar,
                                                                                 mp->d_seismograms_v,
                                                                                 mp->d_nu_rec,
                                                                                 mp->d_ispec_selected_rec_loc,
                                                                                 seismo_current);
      }
#endif
#ifdef USE_HIP
      if (run_hip){
        hipLaunchKernelGGL(compute_elastic_seismogram_kernel, dim3(grid), dim3(threads), 0, mp->compute_stream,
                                                              mp->nrec_local,
                                                              veloc,
                                                              mp->d_ibool,
                                                              mp->d_hxir,mp->d_hetar,mp->d_hgammar,
                                                              mp->d_seismograms_v,
                                                              mp->d_nu_rec,
                                                              mp->d_ispec_selected_rec_loc,
                                                              seismo_current);
      }
#endif
    }

    if (mp->save_seismograms_a){
#ifdef USE_CUDA
      if (run_cuda){
        compute_elastic_seismogram_kernel<<<grid,threads,0,mp->compute_stream>>>(mp->nrec_local,
                                                                                 accel,
                                                                                 mp->d_ibool,
                                                                                 mp->d_hxir,mp->d_hetar,mp->d_hgammar,
                                                                                 mp->d_seismograms_a,
                                                                                 mp->d_nu_rec,
                                                                                 mp->d_ispec_selected_rec_loc,
                                                                                 seismo_current);
      }
#endif
#ifdef USE_HIP
      if (run_hip){
        hipLaunchKernelGGL(compute_elastic_seismogram_kernel, dim3(grid), dim3(threads), 0, mp->compute_stream,
                                                              mp->nrec_local,
                                                              accel,
                                                              mp->d_ibool,
                                                              mp->d_hxir,mp->d_hetar,mp->d_hgammar,
                                                              mp->d_seismograms_a,
                                                              mp->d_nu_rec,
                                                              mp->d_ispec_selected_rec_loc,
                                                              seismo_current);
      }
#endif
    }
  } // elastic

  // acoustic wavefield
  if (*ACOUSTIC_SIMULATION){
    if (mp->save_seismograms_p){
      if (*USE_TRICK_FOR_BETTER_PRESSURE){
#ifdef USE_CUDA
        if (run_cuda){
          compute_acoustic_seismogram_kernel<<<grid,threads,0,mp->compute_stream>>>(mp->nrec_local,
                                                                                    potential_acoustic,
                                                                                    mp->d_ibool,
                                                                                    mp->d_hxir,mp->d_hetar,mp->d_hgammar,
                                                                                    mp->d_seismograms_p,
                                                                                    mp->d_ispec_selected_rec_loc,
                                                                                    seismo_current);
        }
#endif
#ifdef USE_HIP
        if (run_hip){
          hipLaunchKernelGGL(compute_acoustic_seismogram_kernel, dim3(grid), dim3(threads), 0, mp->compute_stream,
                                                                 mp->nrec_local,
                                                                 potential_acoustic,
                                                                 mp->d_ibool,
                                                                 mp->d_hxir,mp->d_hetar,mp->d_hgammar,
                                                                 mp->d_seismograms_p,
                                                                 mp->d_ispec_selected_rec_loc,
                                                                 seismo_current);
        }
#endif
      }else{
#ifdef USE_CUDA
        if (run_cuda){
          compute_acoustic_seismogram_kernel<<<grid,threads,0,mp->compute_stream>>>(mp->nrec_local,
                                                                                    potential_dot_dot_acoustic,
                                                                                    mp->d_ibool,
                                                                                    mp->d_hxir,mp->d_hetar,mp->d_hgammar,
                                                                                    mp->d_seismograms_p,
                                                                                    mp->d_ispec_selected_rec_loc,
                                                                                    seismo_current);
        }
#endif
#ifdef USE_HIP
        if (run_hip){
          hipLaunchKernelGGL(compute_acoustic_seismogram_kernel, dim3(grid), dim3(threads), 0, mp->compute_stream,
                                                                 mp->nrec_local,
                                                                 potential_dot_dot_acoustic,
                                                                 mp->d_ibool,
                                                                 mp->d_hxir,mp->d_hetar,mp->d_hgammar,
                                                                 mp->d_seismograms_p,
                                                                 mp->d_ispec_selected_rec_loc,
                                                                 seismo_current);
        }
#endif
      }
    }

// VM VM add computation of vectorial field in fluids ----------------------------------------------------------------
    if (mp->save_seismograms_d){
#ifdef USE_CUDA
      if (run_cuda){
        compute_acoustic_vectorial_seismogram_kernel<<<grid,threads,0,mp->compute_stream>>>(mp->nrec_local,
                                                                                            mp->d_ispec_is_acoustic,
                                                                                            potential_acoustic,
                                                                                            mp->d_seismograms_d,
                                                                                            mp->d_rhostore,
                                                                                            mp->d_ibool,
                                                                                            mp->d_irregular_element_number,
                                                                                            mp->d_hxir,mp->d_hetar,mp->d_hgammar,
                                                                                            mp->d_xix,mp->d_xiy,mp->d_xiz,
                                                                                            mp->d_etax,mp->d_etay,mp->d_etaz,
                                                                                            mp->d_gammax,mp->d_gammay,mp->d_gammaz,
                                                                                            mp->xix_regular,
                                                                                            mp->d_hprime_xx,
                                                                                            mp->d_nu_rec,
                                                                                            mp->d_ispec_selected_rec_loc,
                                                                                            seismo_current);
      }
#endif
#ifdef USE_HIP
      if (run_hip){
        hipLaunchKernelGGL(compute_acoustic_vectorial_seismogram_kernel, dim3(grid), dim3(threads), 0, mp->compute_stream,
                                                                         mp->nrec_local,
                                                                         mp->d_ispec_is_acoustic,
                                                                         potential_acoustic,
                                                                         mp->d_seismograms_d,
                                                                         mp->d_rhostore,
                                                                         mp->d_ibool,
                                                                         mp->d_irregular_element_number,
                                                                         mp->d_hxir,mp->d_hetar,mp->d_hgammar,
                                                                         mp->d_xix,mp->d_xiy,mp->d_xiz,
                                                                         mp->d_etax,mp->d_etay,mp->d_etaz,
                                                                         mp->d_gammax,mp->d_gammay,mp->d_gammaz,
                                                                         mp->xix_regular,
                                                                         mp->d_hprime_xx,
                                                                         mp->d_nu_rec,
                                                                         mp->d_ispec_selected_rec_loc,
                                                                         seismo_current);
      }
#endif
    }

    if (mp->save_seismograms_v){
#ifdef USE_CUDA
      if (run_cuda){
        compute_acoustic_vectorial_seismogram_kernel<<<grid,threads,0,mp->compute_stream>>>(mp->nrec_local,
                                                                                            mp->d_ispec_is_acoustic,
                                                                                            potential_dot_acoustic,
                                                                                            mp->d_seismograms_v,
                                                                                            mp->d_rhostore,
                                                                                            mp->d_ibool,
                                                                                            mp->d_irregular_element_number,
                                                                                            mp->d_hxir,mp->d_hetar,mp->d_hgammar,
                                                                                            mp->d_xix,mp->d_xiy,mp->d_xiz,
                                                                                            mp->d_etax,mp->d_etay,mp->d_etaz,
                                                                                            mp->d_gammax,mp->d_gammay,mp->d_gammaz,
                                                                                            mp->xix_regular,
                                                                                            mp->d_hprime_xx,
                                                                                            mp->d_nu_rec,
                                                                                            mp->d_ispec_selected_rec_loc,
                                                                                            seismo_current);
      }
#endif
#ifdef USE_HIP
      if (run_hip){
        hipLaunchKernelGGL(compute_acoustic_vectorial_seismogram_kernel, dim3(grid), dim3(threads), 0, mp->compute_stream,
                                                                         mp->nrec_local,
                                                                         mp->d_ispec_is_acoustic,
                                                                         potential_dot_acoustic,
                                                                         mp->d_seismograms_v,
                                                                         mp->d_rhostore,
                                                                         mp->d_ibool,
                                                                         mp->d_irregular_element_number,
                                                                         mp->d_hxir,mp->d_hetar,mp->d_hgammar,
                                                                         mp->d_xix,mp->d_xiy,mp->d_xiz,
                                                                         mp->d_etax,mp->d_etay,mp->d_etaz,
                                                                         mp->d_gammax,mp->d_gammay,mp->d_gammaz,
                                                                         mp->xix_regular,
                                                                         mp->d_hprime_xx,
                                                                         mp->d_nu_rec,
                                                                         mp->d_ispec_selected_rec_loc,
                                                                         seismo_current);
      }
#endif
    }

    if (mp->save_seismograms_a){
#ifdef USE_CUDA
      if (run_cuda){
        compute_acoustic_vectorial_seismogram_kernel<<<grid,threads,0,mp->compute_stream>>>(mp->nrec_local,
                                                                                            mp->d_ispec_is_acoustic,
                                                                                            potential_dot_dot_acoustic,
                                                                                            mp->d_seismograms_a,
                                                                                            mp->d_rhostore,
                                                                                            mp->d_ibool,
                                                                                            mp->d_irregular_element_number,
                                                                                            mp->d_hxir,mp->d_hetar,mp->d_hgammar,
                                                                                            mp->d_xix,mp->d_xiy,mp->d_xiz,
                                                                                            mp->d_etax,mp->d_etay,mp->d_etaz,
                                                                                            mp->d_gammax,mp->d_gammay,mp->d_gammaz,
                                                                                            mp->xix_regular,
                                                                                            mp->d_hprime_xx,
                                                                                            mp->d_nu_rec,
                                                                                            mp->d_ispec_selected_rec_loc,
                                                                                            seismo_current);
      }
#endif
#ifdef USE_HIP
      if (run_hip){
        hipLaunchKernelGGL(compute_acoustic_vectorial_seismogram_kernel, dim3(grid), dim3(threads), 0, mp->compute_stream,
                                                                         mp->nrec_local,
                                                                         mp->d_ispec_is_acoustic,
                                                                         potential_dot_dot_acoustic,
                                                                         mp->d_seismograms_a,
                                                                         mp->d_rhostore,
                                                                         mp->d_ibool,
                                                                         mp->d_irregular_element_number,
                                                                         mp->d_hxir,mp->d_hetar,mp->d_hgammar,
                                                                         mp->d_xix,mp->d_xiy,mp->d_xiz,
                                                                         mp->d_etax,mp->d_etay,mp->d_etaz,
                                                                         mp->d_gammax,mp->d_gammay,mp->d_gammaz,
                                                                         mp->xix_regular,
                                                                         mp->d_hprime_xx,
                                                                         mp->d_nu_rec,
                                                                         mp->d_ispec_selected_rec_loc,
                                                                         seismo_current);
      }
#endif
    }
  } // ACOUSTIC_SIMULATION

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
    int size = mp->nrec_local * nlength_seismogram;

    // (cudaMemcpy implicitly synchronizes all other cuda operations)
    if (mp->save_seismograms_d)
      gpuMemcpy_tohost_realw(seismograms_d,mp->d_seismograms_d,NDIM * size);
    if (mp->save_seismograms_v)
      gpuMemcpy_tohost_realw(seismograms_v,mp->d_seismograms_v,NDIM * size);
    if (mp->save_seismograms_a)
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
      gpuMemcpy_tohost_realw(seismo_temp,mp->d_seismograms_p,size * NB_RUNS_ACOUSTIC_GPU);

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

