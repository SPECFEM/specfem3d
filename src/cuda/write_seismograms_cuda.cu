/*
 !=====================================================================
 !
 !               S p e c f e m 3 D  V e r s i o n  2 . 1
 !               ---------------------------------------
 !
 !          Main authors: Dimitri Komatitsch and Jeroen Tromp
 !    Princeton University, USA and CNRS / INRIA / University of Pau
 ! (c) Princeton University / California Institute of Technology and CNRS / INRIA / University of Pau
 !                             July 2012
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

#include <stdio.h>
#include <cuda.h>
#include <cublas.h>

#include <sys/types.h>
#include <unistd.h>

#include "config.h"
#include "mesh_constants_cuda.h"


/* ----------------------------------------------------------------------------------------------- */

// ELASTIC simulations

/* ----------------------------------------------------------------------------------------------- */

//fortran code snippet...
/*
  ! gets global number of that receiver
  irec = number_receiver_global(irec_local)

  ! gets local receiver interpolators
  ! (1-D Lagrange interpolators)
  hxir(:) = hxir_store(irec_local,:)
  hetar(:) = hetar_store(irec_local,:)
  hgammar(:) = hgammar_store(irec_local,:)
*/

/* ----------------------------------------------------------------------------------------------- */

// unused...
/*
__device__ double my_atomicAdd(double* address, double val) {

    unsigned long long int* address_as_ull = (unsigned long long int*)address;
    unsigned long long int old = *address_as_ull, assumed;
    do{
      assumed = old;
      old = atomicCAS(address_as_ull, assumed, __double_as_longlong(val + __longlong_as_double(assumed)));
    } while (assumed != old);
    return __longlong_as_double(old);
}
*/

/* ----------------------------------------------------------------------------------------------- */

__global__ void compute_interpolated_dva_plus_seismogram(int nrec_local,
                                                         realw* displ, realw* veloc, realw* accel,
                                                         int* ibool,
                                                         double* hxir, double* hetar, double* hgammar,
                                                         realw* seismograms_d, realw* seismograms_v, realw* seismograms_a,
                                                         double* nu,
                                                         int* number_receiver_global,
                                                         int* ispec_selected_rec) {
  int irec_local = blockIdx.x + blockIdx.y*gridDim.x;
  int i = threadIdx.x;
  int j = threadIdx.y;
  int k = threadIdx.z;
  int ijk = i+5*(j+5*(k));

  // we do the **d variable reduction in shared memory, because the
  // atomicAdd() should be faster on the shared memory registers
  // according to
  // http://supercomputingblog.com/cuda/cuda-tutorial-4-atomic-operations/
  __shared__ double sh_dxd[NGLL3];
  __shared__ double sh_dyd[NGLL3];
  __shared__ double sh_dzd[NGLL3];
  __shared__ double sh_vxd[NGLL3];
  __shared__ double sh_vyd[NGLL3];
  __shared__ double sh_vzd[NGLL3];
  __shared__ double sh_axd[NGLL3];
  __shared__ double sh_ayd[NGLL3];
  __shared__ double sh_azd[NGLL3];

  if(irec_local < nrec_local) {
    int irec = number_receiver_global[irec_local]-1;
    int ispec = ispec_selected_rec[irec]-1;
    int iglob = ibool[ijk+125*ispec]-1;
    double hlagrange = hxir[irec_local + nrec_local*i]*hetar[irec_local + nrec_local*j]*hgammar[irec_local + nrec_local*k];
    sh_dxd[ijk] = hlagrange*displ[0+3*iglob];
    sh_dyd[ijk] = hlagrange*displ[1+3*iglob];
    sh_dzd[ijk] = hlagrange*displ[2+3*iglob];

    sh_vxd[ijk] = hlagrange*veloc[0+3*iglob];
    sh_vyd[ijk] = hlagrange*veloc[1+3*iglob];
    sh_vzd[ijk] = hlagrange*veloc[2+3*iglob];

    sh_axd[ijk] = hlagrange*accel[0+3*iglob];
    sh_ayd[ijk] = hlagrange*accel[1+3*iglob];
    sh_azd[ijk] = hlagrange*accel[2+3*iglob];

    // the reduction has to skip the first element (we don't need to
    // add element 0 to itself) This reduction serializes the code,
    // but it should be fast enough --- it can be made faster with a
    // proper reduction algorithm.
    __syncthreads();

    // if(ijk>0) {
    // reduction needs to be done atomically to avoid race conditions
      // atomicAdd(&sh_dxd[0],sh_dxd[ijk]);
      // atomicAdd(&sh_dyd[0],sh_dyd[ijk]);
      // atomicAdd(&sh_dzd[0],sh_dzd[ijk]);

      // atomicAdd(&sh_vxd[0],sh_vxd[ijk]);
      // atomicAdd(&sh_vyd[0],sh_vyd[ijk]);
      // atomicAdd(&sh_vzd[0],sh_vzd[ijk]);

      // atomicAdd(&sh_axd[0],sh_axd[ijk]);
      // atomicAdd(&sh_ayd[0],sh_ayd[ijk]);
      // atomicAdd(&sh_azd[0],sh_azd[ijk]);
    // }
    // __syncthreads();
    if(ijk==0) {
      // a loop in thread 0 is 4 times faster than atomic operations
      for(int i=1;i<125;i++) {
        sh_dxd[0] += sh_dxd[i];
        sh_dyd[0] += sh_dyd[i];
        sh_dzd[0] += sh_dzd[i];

        sh_vxd[0] += sh_vxd[i];
        sh_vyd[0] += sh_vyd[i];
        sh_vzd[0] += sh_vzd[i];

        sh_axd[0] += sh_axd[i];
        sh_ayd[0] += sh_ayd[i];
        sh_azd[0] += sh_azd[i];

      }

      seismograms_d[0+3*irec_local] = nu[0+3*(0+3*irec)]*sh_dxd[0] + nu[0+3*(1+3*irec)]*sh_dyd[0] + nu[0+3*(2+3*irec)]*sh_dzd[0];
      seismograms_d[1+3*irec_local] = nu[1+3*(0+3*irec)]*sh_dxd[0] + nu[1+3*(1+3*irec)]*sh_dyd[0] + nu[1+3*(2+3*irec)]*sh_dzd[0];
      seismograms_d[2+3*irec_local] = nu[2+3*(0+3*irec)]*sh_dxd[0] + nu[2+3*(1+3*irec)]*sh_dyd[0] + nu[2+3*(2+3*irec)]*sh_dzd[0];

      seismograms_v[0+3*irec_local] = nu[0+3*(0+3*irec)]*sh_vxd[0] + nu[0+3*(1+3*irec)]*sh_vyd[0] + nu[0+3*(2+3*irec)]*sh_vzd[0];
      seismograms_v[1+3*irec_local] = nu[1+3*(0+3*irec)]*sh_vxd[0] + nu[1+3*(1+3*irec)]*sh_vyd[0] + nu[1+3*(2+3*irec)]*sh_vzd[0];
      seismograms_v[2+3*irec_local] = nu[2+3*(0+3*irec)]*sh_vxd[0] + nu[2+3*(1+3*irec)]*sh_vyd[0] + nu[2+3*(2+3*irec)]*sh_vzd[0];

      seismograms_a[0+3*irec_local] = nu[0+3*(0+3*irec)]*sh_axd[0] + nu[0+3*(1+3*irec)]*sh_ayd[0] + nu[0+3*(2+3*irec)]*sh_azd[0];
      seismograms_a[1+3*irec_local] = nu[1+3*(0+3*irec)]*sh_axd[0] + nu[1+3*(1+3*irec)]*sh_ayd[0] + nu[1+3*(2+3*irec)]*sh_azd[0];
      seismograms_a[2+3*irec_local] = nu[2+3*(0+3*irec)]*sh_axd[0] + nu[2+3*(1+3*irec)]*sh_ayd[0] + nu[2+3*(2+3*irec)]*sh_azd[0];

    }
  }
}


/* ----------------------------------------------------------------------------------------------- */

extern "C"
void FC_FUNC_(transfer_seismograms_el_from_d,
              TRANSFER_SEISMOGRAMS_EL_FROM_D)(int* nrec_local,
                                              long* Mesh_pointer_f,
                                              realw* seismograms_d,
                                              realw* seismograms_v,
                                              realw* seismograms_a,
                                              int* it) {

// transfers seismograms from device to host

  TRACE("\ttransfer_seismograms_el_from_d");

  Mesh* mp = (Mesh*)(*Mesh_pointer_f); // get Mesh from fortran integer wrapper

  int num_blocks_x, num_blocks_y;
  get_blocks_xy(*nrec_local,&num_blocks_x,&num_blocks_y);

  dim3 grid(num_blocks_x,num_blocks_y);
  dim3 threads(5,5,5);

  // double h_debug[125]; for(int i=0;i<125;i++) h_debug[i] = 0;
  // double* d_debug; cudaMalloc((void**)&d_debug,125*sizeof(double));
  // cudaMemcpy(d_debug,h_debug,125*sizeof(double),cudaMemcpyHostToDevice);
  // Cuda timing
  // cudaEvent_t start, stop;
  // realw time;
  // cudaEventCreate(&start);
  // cudaEventCreate(&stop);
  // cudaEventRecord( start, 0 );

  compute_interpolated_dva_plus_seismogram<<<grid,threads,0,mp->compute_stream>>>(*nrec_local,
                                                                                  mp->d_displ,mp->d_veloc,mp->d_accel,
                                                                                  mp->d_ibool,
                                                                                  mp->d_hxir, mp->d_hetar, mp->d_hgammar,
                                                                                  mp->d_seismograms_d,
                                                                                  mp->d_seismograms_v,
                                                                                  mp->d_seismograms_a,
                                                                                  mp->d_nu,
                                                                                  mp->d_number_receiver_global,
                                                                                  mp->d_ispec_selected_rec
                                                                                  );

  // cudaMemcpy(h_debug,d_debug,125*sizeof(double),cudaMemcpyDeviceToHost);

  // (cudaMemcpy implicitly synchronizes all other cuda operations)
  print_CUDA_error_if_any(cudaMemcpy(mp->h_seismograms_d_it,mp->d_seismograms_d,sizeof(realw)*3* *nrec_local,cudaMemcpyDeviceToHost),72001);
  print_CUDA_error_if_any(cudaMemcpy(mp->h_seismograms_v_it,mp->d_seismograms_v,sizeof(realw)*3* *nrec_local,cudaMemcpyDeviceToHost),72002);
  print_CUDA_error_if_any(cudaMemcpy(mp->h_seismograms_a_it,mp->d_seismograms_a,sizeof(realw)*3* *nrec_local,cudaMemcpyDeviceToHost),72003);

  // cudaEventRecord( stop, 0 );
  // cudaEventSynchronize( stop );
  // cudaEventElapsedTime( &time, start, stop );
  // cudaEventDestroy( start );
  // cudaEventDestroy( stop );
  // printf("seismogram Execution Time: %f ms\n",time);

  // if(abs(mp->h_seismograms_d_it[0]) < 1e-25) printf("seismo1_x=%e\n",mp->h_seismograms_d_it[0]);
  // if(abs(mp->h_seismograms_d_it[1]) < 1e-25) printf("seismo1_y=%e\n",mp->h_seismograms_d_it[1]);

  // if(abs(mp->h_seismograms_d_it[2]) < 1e-25) {

  // printf("%d:seismo1_z=%e\n",*it,mp->h_seismograms_d_it[2]);

  // }


  memcpy(&seismograms_d[3**nrec_local*(*it-1)],mp->h_seismograms_d_it,3* *nrec_local*sizeof(realw));
  memcpy(&seismograms_v[3**nrec_local*(*it-1)],mp->h_seismograms_v_it,3* *nrec_local*sizeof(realw));
  memcpy(&seismograms_a[3**nrec_local*(*it-1)],mp->h_seismograms_a_it,3* *nrec_local*sizeof(realw));

}

/* ----------------------------------------------------------------------------------------------- */

__global__ void transfer_stations_fields_from_device_kernel(int* number_receiver_global,
                                                            int* ispec_selected_rec,
                                                            int* ibool,
                                                            realw* station_seismo_field,
                                                            realw* desired_field,
                                                            int nrec_local) {
  int blockID = blockIdx.x + blockIdx.y*gridDim.x;
  if(blockID<nrec_local) {
    int irec = number_receiver_global[blockID]-1;
    int ispec = ispec_selected_rec[irec]-1;
    int iglob = ibool[threadIdx.x + NGLL3*ispec]-1;

    station_seismo_field[3*NGLL3*blockID + 3*threadIdx.x+0] = desired_field[3*iglob];
    station_seismo_field[3*NGLL3*blockID + 3*threadIdx.x+1] = desired_field[3*iglob+1];
    station_seismo_field[3*NGLL3*blockID + 3*threadIdx.x+2] = desired_field[3*iglob+2];
  }
}


/* ----------------------------------------------------------------------------------------------- */

void transfer_field_from_device(Mesh* mp, realw* d_field,realw* h_field,
                                int* number_receiver_global,
                                int* d_ispec_selected,
                                int* h_ispec_selected,
                                int* ibool) {

TRACE("\ttransfer_field_from_device");

  // checks if anything to do
  if( mp->nrec_local == 0 ) return;

  int blocksize = NGLL3;

  int num_blocks_x, num_blocks_y;
  get_blocks_xy(mp->nrec_local,&num_blocks_x,&num_blocks_y);

  dim3 grid(num_blocks_x,num_blocks_y);
  dim3 threads(blocksize,1,1);

  // prepare field transfer array on device
  transfer_stations_fields_from_device_kernel<<<grid,threads,0,mp->compute_stream>>>(mp->d_number_receiver_global,
                                                                                      d_ispec_selected,
                                                                                      mp->d_ibool,
                                                                                      mp->d_station_seismo_field,
                                                                                      d_field,
                                                                                      mp->nrec_local);

  // (cudaMemcpy implicitly synchronizes all other cuda operations)
  print_CUDA_error_if_any(cudaMemcpy(mp->h_station_seismo_field,mp->d_station_seismo_field,
                                    (3*NGLL3)*(mp->nrec_local)*sizeof(realw),cudaMemcpyDeviceToHost),71001);

  int irec_local;
  for(irec_local=0;irec_local<mp->nrec_local;irec_local++) {
    int irec = number_receiver_global[irec_local] - 1;
    int ispec = h_ispec_selected[irec] - 1;

    for(int i=0;i<NGLL3;i++) {
      int iglob = ibool[i+NGLL3*ispec] - 1;
      h_field[0+3*iglob] = mp->h_station_seismo_field[0+3*i+irec_local*NGLL3*3];
      h_field[1+3*iglob] = mp->h_station_seismo_field[1+3*i+irec_local*NGLL3*3];
      h_field[2+3*iglob] = mp->h_station_seismo_field[2+3*i+irec_local*NGLL3*3];
    }

  }
#ifdef ENABLE_VERY_SLOW_ERROR_CHECKING
  exit_on_cuda_error("transfer_field_from_device");
#endif
}

/* ----------------------------------------------------------------------------------------------- */

extern "C"
void FC_FUNC_(transfer_station_el_from_device,
              TRANSFER_STATION_EL_FROM_DEVICE)(realw* displ,realw* veloc,realw* accel,
                                                   realw* b_displ, realw* b_veloc, realw* b_accel,
                                                   long* Mesh_pointer_f,int* number_receiver_global,
                                                   int* ispec_selected_rec,int* ispec_selected_source,
                                                   int* ibool) {
TRACE("transfer_station_el_from_device");

  Mesh* mp = (Mesh*)(*Mesh_pointer_f); // get Mesh from fortran integer wrapper

  // checks if anything to do
  if( mp->nrec_local == 0 ) return;

  if(mp->simulation_type == 1) {
    transfer_field_from_device(mp,mp->d_displ,displ, number_receiver_global,
             mp->d_ispec_selected_rec, ispec_selected_rec, ibool);
    transfer_field_from_device(mp,mp->d_veloc,veloc, number_receiver_global,
             mp->d_ispec_selected_rec, ispec_selected_rec, ibool);
    transfer_field_from_device(mp,mp->d_accel,accel, number_receiver_global,
             mp->d_ispec_selected_rec, ispec_selected_rec, ibool);
  }
  else if(mp->simulation_type == 2) {
    transfer_field_from_device(mp,mp->d_displ,displ, number_receiver_global,
             mp->d_ispec_selected_source, ispec_selected_source, ibool);
    transfer_field_from_device(mp,mp->d_veloc,veloc, number_receiver_global,
             mp->d_ispec_selected_source, ispec_selected_source, ibool);
    transfer_field_from_device(mp,mp->d_accel,accel, number_receiver_global,
             mp->d_ispec_selected_source, ispec_selected_source, ibool);
  }
  else if(mp->simulation_type == 3) {
    transfer_field_from_device(mp,mp->d_b_displ,b_displ, number_receiver_global,
             mp->d_ispec_selected_rec, ispec_selected_rec, ibool);
    transfer_field_from_device(mp,mp->d_b_veloc,b_veloc, number_receiver_global,
             mp->d_ispec_selected_rec, ispec_selected_rec, ibool);
    transfer_field_from_device(mp,mp->d_b_accel,b_accel, number_receiver_global,
             mp->d_ispec_selected_rec, ispec_selected_rec, ibool);
  }

}

/* ----------------------------------------------------------------------------------------------- */

// ACOUSTIC simulations

/* ----------------------------------------------------------------------------------------------- */

__global__ void transfer_stations_fields_acoustic_from_device_kernel(int* number_receiver_global,
                                                                     int* ispec_selected_rec,
                                                                     int* ibool,
                                                                     realw* station_seismo_potential,
                                                                     realw* desired_potential) {

  int blockID = blockIdx.x + blockIdx.y*gridDim.x;
  int nodeID = threadIdx.x + blockID*blockDim.x;

  int irec = number_receiver_global[blockID]-1;
  int ispec = ispec_selected_rec[irec]-1;
  int iglob = ibool[threadIdx.x + NGLL3*ispec]-1;

  //if(threadIdx.x == 0 ) printf("node acoustic: %i %i %i %i %i %e \n",blockID,nodeID,irec,ispec,iglob,desired_potential[iglob]);

  station_seismo_potential[nodeID] = desired_potential[iglob];
}

/* ----------------------------------------------------------------------------------------------- */

void transfer_field_acoustic_from_device(Mesh* mp,
                                         realw* d_potential,
                                         realw* h_potential,
                                         int* number_receiver_global,
                                         int* d_ispec_selected,
                                         int* h_ispec_selected,
                                         int* ibool) {

TRACE("transfer_field_acoustic_from_device");

  int irec_local,irec,ispec,iglob,j;

  // checks if anything to do
  if( mp->nrec_local < 1 ) return;

  // sets up kernel dimensions
  int blocksize = NGLL3;

  int num_blocks_x, num_blocks_y;
  get_blocks_xy(mp->nrec_local,&num_blocks_x,&num_blocks_y);

  dim3 grid(num_blocks_x,num_blocks_y);
  dim3 threads(blocksize,1,1);

  // prepare field transfer array on device
  transfer_stations_fields_acoustic_from_device_kernel<<<grid,threads>>>(mp->d_number_receiver_global,
                                                                         d_ispec_selected,
                                                                         mp->d_ibool,
                                                                         mp->d_station_seismo_potential,
                                                                         d_potential);


#ifdef ENABLE_VERY_SLOW_ERROR_CHECKING
  exit_on_cuda_error("transfer_field_acoustic_from_device kernel");
#endif

  // (cudaMemcpy implicitly synchronizes all other cuda operations)
  print_CUDA_error_if_any(cudaMemcpy(mp->h_station_seismo_potential,mp->d_station_seismo_potential,
                                     mp->nrec_local*NGLL3*sizeof(realw),cudaMemcpyDeviceToHost),55000);

  //printf("copy local receivers: %i \n",mp->nrec_local);

  for(irec_local=0; irec_local < mp->nrec_local; irec_local++) {
    irec = number_receiver_global[irec_local]-1;
    ispec = h_ispec_selected[irec]-1;

    // copy element values
    // note: iglob may vary and can be irregularly accessing the h_potential array
    for(j=0; j < NGLL3; j++){
      iglob = ibool[j+NGLL3*ispec]-1;
      h_potential[iglob] = mp->h_station_seismo_potential[j+irec_local*NGLL3];
    }

    // copy each station element's points to working array
    // note: this works if iglob values would be all aligned...
    //memcpy(&(h_potential[iglob]),&(mp->h_station_seismo_potential[irec_local*NGLL3]),NGLL3*sizeof(realw));

  }
#ifdef ENABLE_VERY_SLOW_ERROR_CHECKING
  exit_on_cuda_error("transfer_field_acoustic_from_device");
#endif
}

/* ----------------------------------------------------------------------------------------------- */

extern "C"
void FC_FUNC_(transfer_station_ac_from_device,
              TRANSFER_STATION_AC_FROM_DEVICE)(realw* potential_acoustic,
                                                realw* potential_dot_acoustic,
                                                realw* potential_dot_dot_acoustic,
                                                realw* b_potential_acoustic,
                                                realw* b_potential_dot_acoustic,
                                                realw* b_potential_dot_dot_acoustic,
                                                long* Mesh_pointer_f,
                                                int* number_receiver_global,
                                                int* ispec_selected_rec,
                                                int* ispec_selected_source,
                                                int* ibool) {

TRACE("transfer_station_ac_from_device");
  //double start_time = get_time();

  Mesh* mp = (Mesh*)(*Mesh_pointer_f); // get Mesh from fortran integer wrapper

  // checks if anything to do
  if( mp->nrec_local == 0 ) return;

  if(mp->simulation_type == 1) {
    transfer_field_acoustic_from_device(mp,mp->d_potential_acoustic,potential_acoustic,
                                        number_receiver_global,
                                        mp->d_ispec_selected_rec, ispec_selected_rec, ibool);
    transfer_field_acoustic_from_device(mp,mp->d_potential_dot_acoustic,potential_dot_acoustic,
                                        number_receiver_global,
                                        mp->d_ispec_selected_rec, ispec_selected_rec, ibool);
    transfer_field_acoustic_from_device(mp,mp->d_potential_dot_dot_acoustic,potential_dot_dot_acoustic,
                                        number_receiver_global,
                                        mp->d_ispec_selected_rec, ispec_selected_rec, ibool);
  }
  else if(mp->simulation_type == 2) {
    transfer_field_acoustic_from_device(mp,mp->d_potential_acoustic,potential_acoustic,
                                        number_receiver_global,
                                        mp->d_ispec_selected_source, ispec_selected_source, ibool);
    transfer_field_acoustic_from_device(mp,mp->d_potential_dot_acoustic,potential_dot_acoustic,
                                        number_receiver_global,
                                        mp->d_ispec_selected_source, ispec_selected_source, ibool);
    transfer_field_acoustic_from_device(mp,mp->d_potential_dot_dot_acoustic,potential_dot_dot_acoustic,
                                        number_receiver_global,
                                        mp->d_ispec_selected_source, ispec_selected_source, ibool);
  }
  else if(mp->simulation_type == 3) {
    transfer_field_acoustic_from_device(mp,mp->d_b_potential_acoustic,b_potential_acoustic,
                                        number_receiver_global,
                                        mp->d_ispec_selected_rec, ispec_selected_rec, ibool);
    transfer_field_acoustic_from_device(mp,mp->d_b_potential_dot_acoustic,b_potential_dot_acoustic,
                                        number_receiver_global,
                                        mp->d_ispec_selected_rec, ispec_selected_rec, ibool);
    transfer_field_acoustic_from_device(mp,mp->d_b_potential_dot_dot_acoustic,b_potential_dot_dot_acoustic,
                                        number_receiver_global,
                                        mp->d_ispec_selected_rec, ispec_selected_rec, ibool);
  }

#ifdef ENABLE_VERY_SLOW_ERROR_CHECKING
  //double end_time = get_time();
  //printf("Elapsed time: %e\n",end_time-start_time);
  exit_on_cuda_error("transfer_station_ac_from_device");
#endif
}

