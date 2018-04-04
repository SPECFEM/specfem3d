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

__global__ void compute_elastic_seismogram_kernel(int nrec_local,
                                                  realw* field,
                                                  int* d_ibool,
                                                  realw* hxir, realw* hetar, realw* hgammar,
                                                  realw* seismograms,
                                                  realw* nu,
                                                  int* ispec_selected_rec_loc,
                                                  int it){

  int irec_local = blockIdx.x + blockIdx.y*gridDim.x;
  int tx = threadIdx.x;

  // local index
  int K = (tx/NGLL2);
  int J = ((tx-K*NGLL2)/NGLLX);
  int I = (tx-K*NGLL2-J*NGLLX);

  __shared__ realw sh_dxd[NGLL3_PADDED];
  __shared__ realw sh_dyd[NGLL3_PADDED];
  __shared__ realw sh_dzd[NGLL3_PADDED];

  if (irec_local < nrec_local) {

    int ispec = ispec_selected_rec_loc[irec_local] - 1;

    sh_dxd[tx] = 0;
    sh_dyd[tx] = 0;
    sh_dzd[tx] = 0;

    if (tx < NGLL3) {
      realw hlagrange = hxir[irec_local + nrec_local*I]*hetar[irec_local + nrec_local*J]*hgammar[irec_local + nrec_local*K];
      int iglob = iglob = d_ibool[INDEX4_PADDED(NGLLX,NGLLX,NGLLX,I,J,K,ispec)]-1;

      sh_dxd[tx] = hlagrange * field[0 + 3*iglob];
      sh_dyd[tx] = hlagrange * field[1 + 3*iglob];
      sh_dzd[tx] = hlagrange * field[2 + 3*iglob];

      //debug
      //if (tx == 0) printf("thread %d %d %d - %f %f %f\n",ispec,iglob,irec_local,hlagrange,field[0 + 2*iglob],field[1 + 2*iglob]);
    }
    __syncthreads();

    // reduction
    for (unsigned int s=1; s<NGLL3_PADDED ; s *= 2) {
      if (tx % (2*s) == 0){ sh_dxd[tx] += sh_dxd[tx + s];
                            sh_dyd[tx] += sh_dyd[tx + s];
                            sh_dzd[tx] += sh_dzd[tx + s];}
      __syncthreads();
    }

    int idx = INDEX3(NDIM,nrec_local,0,irec_local,it);

    if (tx == 0) {
      seismograms[0+idx] = nu[0+3*(0+3*irec_local)]*sh_dxd[0] + nu[0+3*(1+3*irec_local)]*sh_dyd[0] + nu[0+3*(2+3*irec_local)]*sh_dzd[0];
    }
    if (tx == 1) {
      seismograms[1+idx] = nu[1+3*(0+3*irec_local)]*sh_dxd[0] + nu[1+3*(1+3*irec_local)]*sh_dyd[0] + nu[1+3*(2+3*irec_local)]*sh_dzd[0];
    }
    if (tx == 2) {
      seismograms[2+idx] = nu[2+3*(0+3*irec_local)]*sh_dxd[0] + nu[2+3*(1+3*irec_local)]*sh_dyd[0] + nu[2+3*(2+3*irec_local)]*sh_dzd[0];
    }
  }
}

/* ----------------------------------------------------------------------------------------------- */

__global__ void compute_acoustic_seismogram_kernel(int nrec_local,
                                                   field* pressure,
                                                   int* d_ibool,
                                                   realw* hxir, realw* hetar, realw* hgammar,
                                                   field* seismograms,
                                                   int* ispec_selected_rec_loc,
                                                   int it){

  int irec_local = blockIdx.x + blockIdx.y*gridDim.x;
  int tx = threadIdx.x;

  // local index
  int K = (tx/NGLL2);
  int J = ((tx-K*NGLL2)/NGLLX);
  int I = (tx-K*NGLL2-J*NGLLX);

  __shared__ field sh_dxd[NGLL3_PADDED];

  if (irec_local < nrec_local) {

    int ispec = ispec_selected_rec_loc[irec_local]-1;

    sh_dxd[tx] = Make_field(0.f);

    if (tx < NGLL3) {

      realw hlagrange = hxir[irec_local + nrec_local*I]*hetar[irec_local + nrec_local*J]*hgammar[irec_local + nrec_local*K];
      int iglob = iglob = d_ibool[INDEX4_PADDED(NGLLX,NGLLX,NGLLX,I,J,K,ispec)]-1;

      sh_dxd[tx] = hlagrange*pressure[iglob];
    }
    __syncthreads();

    for (unsigned int s=1; s<NGLL3_PADDED ; s *= 2) {
      if (tx % (2*s) == 0) {sh_dxd[tx] += sh_dxd[tx + s];}
      __syncthreads();
    }

    int idx = INDEX2(nrec_local,irec_local,it);

    // Signe moins car pression = -potential_dot_dot
   if (tx == 0) seismograms[idx] = -sh_dxd[0];
  }
}

/* ----------------------------------------------------------------------------------------------- */

__global__ void compute_acoustic_vectorial_seismogram_kernel(int nrec_local,
                   int*  d_ispec_is_acoustic,
                   field* scalar_potential,
                   realw* seismograms,
                   realw* d_rhostore,
                   int* d_ibool,
                   int * d_irregular_element_number,
                   realw* hxir, realw* hetar, realw* hgammar,
                   realw* d_xix, realw* d_xiy, realw* d_xiz,
                   realw* d_etax, realw* d_etay, realw* d_etaz,
                   realw* d_gammax, realw* d_gammay, realw* d_gammaz,
                   realw xix_regular,
                   realw* d_hprime_xx,
                   realw* nu,
                   int* ispec_selected_rec_loc,
                   int it){

  int irec_local = blockIdx.x + blockIdx.y*gridDim.x;
  int tx = threadIdx.x;

  // shared memory
  __shared__ realw s_dummy_loc[NGLL3_PADDED];
  __shared__ realw s_temp1[NGLL3_PADDED];
  __shared__ realw s_temp2[NGLL3_PADDED];
  __shared__ realw s_temp3[NGLL3_PADDED];
  __shared__ realw sh_hprime_xx[NGLL2];

  // locals
  realw temp1l, temp2l, temp3l;
  realw rho_invl, hlagrange;
  realw xixl, xiyl, xizl;
  realw etaxl, etayl, etazl;
  realw gammaxl, gammayl, gammazl;
  realw dpotentialdxl, dpotentialdyl, dpotentialdzl;
  int ispec, offset, offset_irreg, iglob, ispec_irreg;

  /*
  // debug
  if (irec_local < nrec_local) {
    ispec = ispec_selected_rec_loc[irec_local] - 1;
    offset = INDEX4_PADDED(NGLLX,NGLLX,NGLLX,I,J,K,ispec);
    iglob = d_ibool[offset]-1;
    rho_invl = 1.f / d_rhostore[offset];
    xixl = d_xix[offset];
    xiyl = d_xiy[offset];
    xizl = d_xiz[offset];
    etaxl = d_etax[offset];
    etayl = d_etay[offset];
    etazl = d_etaz[offset];
    gammaxl = d_gammax[offset];
    gammayl = d_gammay[offset];
    gammazl = d_gammaz[offset];

    hlagrange = hxir[irec_local + nrec_local*I]*hetar[irec_local + nrec_local*J]*hgammar[irec_local + nrec_local*K];
    // loads into shared memory
    if (tx < NGLL2) {
      sh_hprime_xx[tx] = d_hprime_xx[tx];}
    s_dummy_loc[tx] = 1.; //scalar_potential[iglob];
    if (iglob > 0) {
      printf(" iglob =%d, (i,j,k)=(%d,%d,%d), ispec =%d  --- %f \n", iglob, I, J, K, ispec, scalar_potential[iglob]);}
    else{
      printf(" -illegal %d  %d %d %d %d\n", tx, ispec, I, J, K);
    }
  }
  */

  s_temp1[tx] = 0.0f;
  s_temp2[tx] = 0.0f;
  s_temp3[tx] = 0.0f;

  // local index
  int K = (tx/NGLL2);
  int J = ((tx-K*NGLL2)/NGLLX);
  int I = (tx-K*NGLL2-J*NGLLX);


  if (irec_local >= nrec_local) return;

  if (tx < NGLL3) {
    ispec = ispec_selected_rec_loc[irec_local] - 1;
    ispec_irreg = d_irregular_element_number[ispec] - 1;
    // nothing to do if we are in elastic element
    if (d_ispec_is_acoustic[ispec] == 0) {return;}

    offset = INDEX4_PADDED(NGLLX,NGLLX,NGLLX,I,J,K,ispec);
    offset_irreg = INDEX4_PADDED(NGLLX,NGLLX,NGLLX,I,J,K,ispec_irreg);

    iglob = d_ibool[offset]-1;
    rho_invl = 1.f / d_rhostore[offset];
    hlagrange = hxir[irec_local + nrec_local*I]*hetar[irec_local + nrec_local*J]*hgammar[irec_local + nrec_local*K];
  }

  //debug
  //if (tx == 0) printf("thread %d %d %d - %f %f %f\n",ispec,iglob,irec_local,hlagrange,rho_invl, xixl);

  // loads into shared memory
  if (tx < NGLL2) sh_hprime_xx[tx] = d_hprime_xx[tx];
  if (tx < NGLL3) s_dummy_loc[tx] = (realw)scalar_potential[iglob];


  // synchronize all the threads (one thread for each of the NGLL grid points of the
  // current spectral element) because we need the whole element to be ready
  __syncthreads();

  if (tx < NGLL3) {
    // computes first matrix product
    temp1l = 0.f;
    temp2l = 0.f;
    temp3l = 0.f;

    for (int l=0;l<NGLLX;l++) {
      //assumes that hprime_xx = hprime_yy = hprime_zz
      // 1. cut-plane along xi-direction
      temp1l += s_dummy_loc[K*NGLL2+J*NGLLX+l] * sh_hprime_xx[l*NGLLX+I];
      // 2. cut-plane along eta-direction
      temp2l += s_dummy_loc[K*NGLL2+l*NGLLX+I] * sh_hprime_xx[l*NGLLX+J];
      // 3. cut-plane along gamma-direction
      temp3l += s_dummy_loc[l*NGLL2+J*NGLLX+I] * sh_hprime_xx[l*NGLLX+K];
    }

   if (ispec_irreg >= 0){ //irregular element
      xixl = d_xix[offset_irreg];
      xiyl = d_xiy[offset_irreg];
      xizl = d_xiz[offset_irreg];
      etaxl = d_etax[offset_irreg];
      etayl = d_etay[offset_irreg];
      etazl = d_etaz[offset_irreg];
      gammaxl = d_gammax[offset_irreg];
      gammayl = d_gammay[offset_irreg];
      gammazl = d_gammaz[offset_irreg];
      // compute derivatives of ux, uy and uz with respect to x, y and z
      // derivatives of potential
      dpotentialdxl = xixl*temp1l + etaxl*temp2l + gammaxl*temp3l;
      dpotentialdyl = xiyl*temp1l + etayl*temp2l + gammayl*temp3l;
      dpotentialdzl = xizl*temp1l + etazl*temp2l + gammazl*temp3l;
    }
    else{
      // compute derivatives of ux, uy and uz with respect to x, y and z
      // derivatives of potential
      dpotentialdxl = xix_regular*temp1l;
      dpotentialdyl = xix_regular*temp2l;
      dpotentialdzl = xix_regular*temp3l;
    }

    // store the field in shared memmory
    s_temp1[tx] = hlagrange *dpotentialdxl * rho_invl;
    s_temp2[tx] = hlagrange *dpotentialdyl * rho_invl;
    s_temp3[tx] = hlagrange *dpotentialdzl * rho_invl;
  }

  __syncthreads();

  // reduction
  for (unsigned int s=1; s<NGLL3_PADDED ; s *= 2) {
    if (tx % (2*s) == 0){ s_temp1[tx] += s_temp1[tx + s];
                          s_temp2[tx] += s_temp2[tx + s];
                          s_temp3[tx] += s_temp3[tx + s];}
    __syncthreads();
  }

  int idx = INDEX3(NDIM,nrec_local,0,irec_local,it);

  if (tx == 0) {
    seismograms[0+idx] = nu[0+3*(0+3*irec_local)]*s_temp1[0] + nu[0+3*(1+3*irec_local)]*s_temp2[0] + nu[0+3*(2+3*irec_local)]*s_temp3[0];
  }
  if (tx == 1) {
    seismograms[1+idx] = nu[1+3*(0+3*irec_local)]*s_temp1[0] + nu[1+3*(1+3*irec_local)]*s_temp2[0] + nu[1+3*(2+3*irec_local)]*s_temp3[0];
  }
  if (tx == 2) {
    seismograms[2+idx] = nu[2+3*(0+3*irec_local)]*s_temp1[0] + nu[2+3*(1+3*irec_local)]*s_temp2[0] + nu[2+3*(2+3*irec_local)]*s_temp3[0];
  }
}

/* ----------------------------------------------------------------------------------------------- */

extern "C"
void FC_FUNC_(compute_seismograms_cuda,
              COMPUTE_SEISMOGRAMS_CUDA)(long* Mesh_pointer_f,
                                        realw* seismograms_d,
                                        realw* seismograms_v,
                                        realw* seismograms_a,
                                        realw* seismograms_p,
                                        int* seismo_currentf,
                                        int* NTSTEP_BETWEEN_OUTPUT_SEISMOSf,
                                        int* it, int* it_end,
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

  int seismo_current = *seismo_currentf - 1 ;
  int NTSTEP_BETWEEN_OUTPUT_SEISMOS = *NTSTEP_BETWEEN_OUTPUT_SEISMOSf;

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
    if (mp->save_seismograms_d)
        compute_elastic_seismogram_kernel<<<grid,threads,0,mp->compute_stream>>>(mp->nrec_local,
                                                                                 mp->d_displ,
                                                                                 mp->d_ibool,
                                                                                 mp->d_hxir,mp->d_hetar,mp->d_hgammar,
                                                                                 mp->d_seismograms_d,
                                                                                 mp->d_nu,
                                                                                 mp->d_ispec_selected_rec_loc,
                                                                                 seismo_current);


    if (mp->save_seismograms_v)
        compute_elastic_seismogram_kernel<<<grid,threads,0,mp->compute_stream>>>(mp->nrec_local,
                                                                                 mp->d_veloc,
                                                                                 mp->d_ibool,
                                                                                 mp->d_hxir,mp->d_hetar,mp->d_hgammar,
                                                                                 mp->d_seismograms_v,
                                                                                 mp->d_nu,
                                                                                 mp->d_ispec_selected_rec_loc,
                                                                                 seismo_current);

    if (mp->save_seismograms_a)
        compute_elastic_seismogram_kernel<<<grid,threads,0,mp->compute_stream>>>(mp->nrec_local,
                                                                                 mp->d_accel,
                                                                                 mp->d_ibool,
                                                                                 mp->d_hxir,mp->d_hetar,mp->d_hgammar,
                                                                                 mp->d_seismograms_a,
                                                                                 mp->d_nu,
                                                                                 mp->d_ispec_selected_rec_loc,
                                                                                 seismo_current);

  }
  // acoustic wavefield
  if (*ACOUSTIC_SIMULATION){
    if (mp->save_seismograms_p){
        if (*USE_TRICK_FOR_BETTER_PRESSURE){
          compute_acoustic_seismogram_kernel<<<grid,threads,0,mp->compute_stream>>>(mp->nrec_local,
                                                                                    mp->d_potential_acoustic,
                                                                                    mp->d_ibool,
                                                                                    mp->d_hxir,mp->d_hetar,mp->d_hgammar,
                                                                                    mp->d_seismograms_p,
                                                                                    mp->d_ispec_selected_rec_loc,
                                                                                    seismo_current);
        }else{
          compute_acoustic_seismogram_kernel<<<grid,threads,0,mp->compute_stream>>>(mp->nrec_local,
                                                                                    mp->d_potential_dot_dot_acoustic,
                                                                                    mp->d_ibool,
                                                                                    mp->d_hxir,mp->d_hetar,mp->d_hgammar,
                                                                                    mp->d_seismograms_p,
                                                                                    mp->d_ispec_selected_rec_loc,
                                                                                    seismo_current);
        }
    }

// VM VM add computation of vectorial field in fluids ----------------------------------------------------------------
    if (mp->save_seismograms_d)
      compute_acoustic_vectorial_seismogram_kernel<<<grid,threads,0,mp->compute_stream>>>(mp->nrec_local,
                                                                                          mp->d_ispec_is_acoustic,
                                                                                          mp->d_potential_acoustic,
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
                                                                                          mp->d_nu,
                                                                                          mp->d_ispec_selected_rec_loc,
                                                                                          seismo_current);


    if (mp->save_seismograms_v)
      compute_acoustic_vectorial_seismogram_kernel<<<grid,threads,0,mp->compute_stream>>>(mp->nrec_local,
                                                                                          mp->d_ispec_is_acoustic,
                                                                                          mp->d_potential_dot_acoustic,
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
                                                                                          mp->d_nu,
                                                                                          mp->d_ispec_selected_rec_loc,
                                                                                          seismo_current);


    if (mp->save_seismograms_a)
      compute_acoustic_vectorial_seismogram_kernel<<<grid,threads,0,mp->compute_stream>>>(mp->nrec_local,
                                                                                          mp->d_ispec_is_acoustic,
                                                                                          mp->d_potential_dot_dot_acoustic,
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
                                                                                          mp->d_nu,
                                                                                          mp->d_ispec_selected_rec_loc,
                                                                                          seismo_current);
  } // ACOUSTIC_SIMULATION

  if (seismo_current == NTSTEP_BETWEEN_OUTPUT_SEISMOS || *it == *it_end ){
    int size = mp->nrec_local * NTSTEP_BETWEEN_OUTPUT_SEISMOS * sizeof(realw);

    // (cudaMemcpy implicitly synchronizes all other cuda operations)
    if (mp->save_seismograms_d)
      print_CUDA_error_if_any(cudaMemcpy(seismograms_d,mp->d_seismograms_d,NDIM * size,cudaMemcpyDeviceToHost),72001);
    if (mp->save_seismograms_v)
      print_CUDA_error_if_any(cudaMemcpy(seismograms_v,mp->d_seismograms_v,NDIM * size,cudaMemcpyDeviceToHost),72002);
    if (mp->save_seismograms_a)
      print_CUDA_error_if_any(cudaMemcpy(seismograms_a,mp->d_seismograms_a,NDIM * size,cudaMemcpyDeviceToHost),72003);
    // EB EB Temporary solution : in the future we will also declare host pressure seismograms as (1,nrec_local,NTSTEP_BETWEEN_OUTPUT_SEISMOS)
    realw * seismo_temp;
    if (mp->save_seismograms_p){
      // EB EB We need to reorganize data to match host array shape :
      // if NB_RUNS_ACOUSTIC_GPU = 1 from fortran shape (1,nrec_local,NTSTEP_BETWEEN_OUTPUT_SEISMOS) to (NDIM,nrec_local,NTSTEP_BETWEEN_OUTPUT_SEISMOS)
      // if NB_RUNS_ACOUSTIC_GPU > 1 from fortran shape (NB_RUNS_ACOUSTIC_GPU,nrec_local,NTSTEP_BETWEEN_OUTPUT_SEISMOS) to (NDIM,nrec_local*NB_RUNS_ACOUSTIC_GPU,NTSTEP_BETWEEN_OUTPUT_SEISMOS)
      seismo_temp = (realw*)malloc(size*NB_RUNS_ACOUSTIC_GPU);
      print_CUDA_error_if_any(cudaMemcpy(seismo_temp,mp->d_seismograms_p,size*NB_RUNS_ACOUSTIC_GPU,cudaMemcpyDeviceToHost),72004);
      for (int it = 0; it<NTSTEP_BETWEEN_OUTPUT_SEISMOS; it++)
        for (int i_recloc=0; i_recloc<mp->nrec_local; i_recloc++)
          for (int i_run=0; i_run<NB_RUNS_ACOUSTIC_GPU; i_run++){
          seismograms_p[INDEX4(NDIM,mp->nrec_local,NB_RUNS_ACOUSTIC_GPU,0,i_recloc,i_run,it)] = seismo_temp[INDEX3(NB_RUNS_ACOUSTIC_GPU,mp->nrec_local,i_run,i_recloc,it)];
          seismograms_p[INDEX4(NDIM,mp->nrec_local,NB_RUNS_ACOUSTIC_GPU,1,i_recloc,i_run,it)] = 0.f;
          seismograms_p[INDEX4(NDIM,mp->nrec_local,NB_RUNS_ACOUSTIC_GPU,2,i_recloc,i_run,it)] = 0.f;
          }
      free(seismo_temp);
    }
  }

#ifdef ENABLE_VERY_SLOW_ERROR_CHECKING
  exit_on_cuda_error("after compute_seismograms_cuda");
#endif
}

