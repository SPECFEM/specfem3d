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

#ifdef WITH_MPI
#include <mpi.h>
#endif

#include <sys/time.h>
#include <sys/resource.h>

#include "config.h"
#include "mesh_constants_cuda.h"
#include "prepare_constants_cuda.h"

#ifdef USE_OLDER_CUDA4_GPU
#else
  #ifdef USE_TEXTURES_FIELDS
extern texture<realw, cudaTextureType1D, cudaReadModeElementType> d_displ_tex;
extern texture<realw, cudaTextureType1D, cudaReadModeElementType> d_veloc_tex;
extern texture<realw, cudaTextureType1D, cudaReadModeElementType> d_accel_tex;
  #endif

  #ifdef USE_TEXTURES_CONSTANTS
extern texture<realw, cudaTextureType1D, cudaReadModeElementType> d_hprime_xx_tex;
  #endif
#endif


/* ----------------------------------------------------------------------------------------------- */

// helper functions

/* ----------------------------------------------------------------------------------------------- */


// copies integer array from CPU host to GPU device
void copy_todevice_int(void** d_array_addr_ptr,int* h_array,int size){
  TRACE("copy_todevice_int");

  // allocates memory on GPU
  //
  // note: cudaMalloc uses a double-pointer, such that it can return an error code in case it fails
  //          we thus pass the address to the pointer above (as void double-pointer) to have it
  //          pointing to the correct pointer of the array here
  print_CUDA_error_if_any(cudaMalloc((void**)d_array_addr_ptr,size*sizeof(int)),
                          12001);

  // copies values onto GPU
  //
  // note: cudaMemcpy uses the pointer to the array, we thus re-cast the value of
  //          the double-pointer above to have the correct pointer to the array
  print_CUDA_error_if_any(cudaMemcpy((int*) *d_array_addr_ptr,h_array,size*sizeof(int),cudaMemcpyHostToDevice),
                          12002);
}

/* ----------------------------------------------------------------------------------------------- */

// copies integer array from CPU host to GPU device
void copy_todevice_realw(void** d_array_addr_ptr,realw* h_array,int size){
  TRACE("copy_todevice_realw");

  // allocates memory on GPU
  print_CUDA_error_if_any(cudaMalloc((void**)d_array_addr_ptr,size*sizeof(realw)),
                          22001);
  // copies values onto GPU
  print_CUDA_error_if_any(cudaMemcpy((realw*) *d_array_addr_ptr,h_array,size*sizeof(realw),cudaMemcpyHostToDevice),
                          22002);
}


/* ----------------------------------------------------------------------------------------------- */

// GPU preparation

/* ----------------------------------------------------------------------------------------------- */

extern "C"
void FC_FUNC_(prepare_constants_device,
              PREPARE_CONSTANTS_DEVICE)(long* Mesh_pointer,
                                        int* h_NGLLX,
                                        int* NSPEC_AB, int* NGLOB_AB,
                                        realw* h_xix, realw* h_xiy, realw* h_xiz,
                                        realw* h_etax, realw* h_etay, realw* h_etaz,
                                        realw* h_gammax, realw* h_gammay, realw* h_gammaz,
                                        realw* h_kappav, realw* h_muv,
                                        int* h_ibool,
                                        int* num_interfaces_ext_mesh,
                                        int* max_nibool_interfaces_ext_mesh,
                                        int* h_nibool_interfaces_ext_mesh,
                                        int* h_ibool_interfaces_ext_mesh,
                                        realw* h_hprime_xx,
                                        realw* h_hprimewgll_xx,
                                        realw* h_wgllwgll_xy,realw* h_wgllwgll_xz,realw* h_wgllwgll_yz,
                                        int* ABSORBING_CONDITIONS,
                                        int* h_abs_boundary_ispec, int* h_abs_boundary_ijk,
                                        realw* h_abs_boundary_normal,
                                        realw* h_abs_boundary_jacobian2Dw,
                                        int* h_num_abs_boundary_faces,
                                        int* h_ispec_is_inner,
                                        int* NSOURCES,
                                        int* nsources_local_f,
                                        realw* h_sourcearrays,
                                        int* h_islice_selected_source,
                                        int* h_ispec_selected_source,
                                        int* h_number_receiver_global,
                                        int* h_ispec_selected_rec,
                                        int* nrec,int* nrec_local,
                                        int* SIMULATION_TYPE,
                                        int* USE_MESH_COLORING_GPU_f,
                                        int* nspec_acoustic,int* nspec_elastic,
                                        int* my_neighbours_ext_mesh,
                                        int* request_send_vector_ext_mesh,
                                        int* request_recv_vector_ext_mesh,
                                        realw* buffer_recv_vector_ext_mesh
                                        ) {

TRACE("prepare_constants_device");

  // allocates mesh parameter structure
  Mesh* mp = (Mesh*) malloc( sizeof(Mesh) );
  if (mp == NULL) exit_on_error("error allocating mesh pointer");
  *Mesh_pointer = (long)mp;

  // checks if NGLLX == 5
  if( *h_NGLLX != NGLLX ){
    exit_on_error("NGLLX must be 5 for CUDA devices");
  }

  // sets number of processes
#ifdef WITH_MPI
  int nproc;
  MPI_Comm_size(MPI_COMM_WORLD,&nproc);
  mp->NPROC = nproc;
#else
  mp->NPROC = 1;
#endif

  // sets global parameters
  mp->NSPEC_AB = *NSPEC_AB;
  mp->NGLOB_AB = *NGLOB_AB;

  // constants
  mp->simulation_type = *SIMULATION_TYPE;
  mp->absorbing_conditions = *ABSORBING_CONDITIONS;

  // sets constant arrays
  setConst_hprime_xx(h_hprime_xx,mp);
  // setConst_hprime_yy(h_hprime_yy,mp); // only needed if NGLLX != NGLLY != NGLLZ
  // setConst_hprime_zz(h_hprime_zz,mp); // only needed if NGLLX != NGLLY != NGLLZ

  setConst_hprimewgll_xx(h_hprimewgll_xx,mp);
  //setConst_hprimewgll_yy(h_hprimewgll_yy,mp); // only needed if NGLLX != NGLLY != NGLLZ
  //setConst_hprimewgll_zz(h_hprimewgll_zz,mp); // only needed if NGLLX != NGLLY != NGLLZ

  setConst_wgllwgll_xy(h_wgllwgll_xy,mp);
  setConst_wgllwgll_xz(h_wgllwgll_xz,mp);
  setConst_wgllwgll_yz(h_wgllwgll_yz,mp);

  // Using texture memory for the hprime-style constants is slower on
  // Fermi generation hardware, but *may* be faster on Kepler
  // generation. We will reevaluate this again, so might as well leave
  // in the code with with #USE_TEXTURES_FIELDS not-defined.
  #ifdef USE_TEXTURES_CONSTANTS
  {
    #ifdef USE_OLDER_CUDA4_GPU
      const textureReference* d_hprime_xx_tex_ptr;
      print_CUDA_error_if_any(cudaGetTextureReference(&d_hprime_xx_tex_ptr, "d_hprime_xx_tex"), 4101);
      cudaChannelFormatDesc channelDesc = cudaCreateChannelDesc<float>();
      print_CUDA_error_if_any(cudaBindTexture(0, d_hprime_xx_tex_ptr, mp->d_hprime_xx, &channelDesc, sizeof(realw)*(NGLL2)), 4001);
   #else
      cudaChannelFormatDesc channelDesc = cudaCreateChannelDesc<float>();
      print_CUDA_error_if_any(cudaBindTexture(0, &d_hprime_xx, mp->d_hprime_xx, &channelDesc, sizeof(realw)*(NGLL2)), 4001);
   #endif
  }
  #endif

  /* Assuming NGLLX=5. Padded is then 128 (5^3+3) */
  int size_padded = NGLL3_PADDED * (mp->NSPEC_AB);

  // mesh
  print_CUDA_error_if_any(cudaMalloc((void**) &mp->d_xix, size_padded*sizeof(realw)),1001);
  print_CUDA_error_if_any(cudaMalloc((void**) &mp->d_xiy, size_padded*sizeof(realw)),1002);
  print_CUDA_error_if_any(cudaMalloc((void**) &mp->d_xiz, size_padded*sizeof(realw)),1003);
  print_CUDA_error_if_any(cudaMalloc((void**) &mp->d_etax, size_padded*sizeof(realw)),1004);
  print_CUDA_error_if_any(cudaMalloc((void**) &mp->d_etay, size_padded*sizeof(realw)),1005);
  print_CUDA_error_if_any(cudaMalloc((void**) &mp->d_etaz, size_padded*sizeof(realw)),1006);
  print_CUDA_error_if_any(cudaMalloc((void**) &mp->d_gammax, size_padded*sizeof(realw)),1007);
  print_CUDA_error_if_any(cudaMalloc((void**) &mp->d_gammay, size_padded*sizeof(realw)),1008);
  print_CUDA_error_if_any(cudaMalloc((void**) &mp->d_gammaz, size_padded*sizeof(realw)),1009);
  print_CUDA_error_if_any(cudaMalloc((void**) &mp->d_kappav, size_padded*sizeof(realw)),1010);
  print_CUDA_error_if_any(cudaMalloc((void**) &mp->d_muv, size_padded*sizeof(realw)),1011);

  // transfer constant element data with padding
  for(int i=0;i < mp->NSPEC_AB;i++) {
    print_CUDA_error_if_any(cudaMemcpy(mp->d_xix + i*NGLL3_PADDED, &h_xix[i*NGLL3],
                                       NGLL3*sizeof(realw),cudaMemcpyHostToDevice),1501);
    print_CUDA_error_if_any(cudaMemcpy(mp->d_xiy+i*NGLL3_PADDED,   &h_xiy[i*NGLL3],
                                       NGLL3*sizeof(realw),cudaMemcpyHostToDevice),1502);
    print_CUDA_error_if_any(cudaMemcpy(mp->d_xiz+i*NGLL3_PADDED,   &h_xiz[i*NGLL3],
                                       NGLL3*sizeof(realw),cudaMemcpyHostToDevice),1503);
    print_CUDA_error_if_any(cudaMemcpy(mp->d_etax+i*NGLL3_PADDED,  &h_etax[i*NGLL3],
                                       NGLL3*sizeof(realw),cudaMemcpyHostToDevice),1504);
    print_CUDA_error_if_any(cudaMemcpy(mp->d_etay+i*NGLL3_PADDED,  &h_etay[i*NGLL3],
                                       NGLL3*sizeof(realw),cudaMemcpyHostToDevice),1505);
    print_CUDA_error_if_any(cudaMemcpy(mp->d_etaz+i*NGLL3_PADDED,  &h_etaz[i*NGLL3],
                                       NGLL3*sizeof(realw),cudaMemcpyHostToDevice),1506);
    print_CUDA_error_if_any(cudaMemcpy(mp->d_gammax+i*NGLL3_PADDED,&h_gammax[i*NGLL3],
                                       NGLL3*sizeof(realw),cudaMemcpyHostToDevice),1507);
    print_CUDA_error_if_any(cudaMemcpy(mp->d_gammay+i*NGLL3_PADDED,&h_gammay[i*NGLL3],
                                       NGLL3*sizeof(realw),cudaMemcpyHostToDevice),1508);
    print_CUDA_error_if_any(cudaMemcpy(mp->d_gammaz+i*NGLL3_PADDED,&h_gammaz[i*NGLL3],
                                       NGLL3*sizeof(realw),cudaMemcpyHostToDevice),1509);
    print_CUDA_error_if_any(cudaMemcpy(mp->d_kappav+i*NGLL3_PADDED,&h_kappav[i*NGLL3],
                                       NGLL3*sizeof(realw),cudaMemcpyHostToDevice),1510);
    print_CUDA_error_if_any(cudaMemcpy(mp->d_muv+i*NGLL3_PADDED,   &h_muv[i*NGLL3],
                                       NGLL3*sizeof(realw),cudaMemcpyHostToDevice),1511);
  }

  // global indexing
  copy_todevice_int((void**)&mp->d_ibool,h_ibool,NGLL3*(mp->NSPEC_AB));

  //print_CUDA_error_if_any(cudaMalloc((void**) &mp->d_ibool,size_padded*sizeof(int)),1021);
  //print_CUDA_error_if_any(cudaMemcpy(mp->d_ibool, h_ibool,
  //                                  NGLL3*(mp->NSPEC_AB)*sizeof(int),cudaMemcpyHostToDevice),1022);


  // prepare interprocess-edge exchange information
  mp->num_interfaces_ext_mesh = *num_interfaces_ext_mesh;
  mp->max_nibool_interfaces_ext_mesh = *max_nibool_interfaces_ext_mesh;
  if( mp->num_interfaces_ext_mesh > 0 ){

    copy_todevice_int((void**)&mp->d_nibool_interfaces_ext_mesh,h_nibool_interfaces_ext_mesh,
                      mp->num_interfaces_ext_mesh);

    //print_CUDA_error_if_any(cudaMalloc((void**) &mp->d_nibool_interfaces_ext_mesh,
    //                                   (mp->num_interfaces_ext_mesh)*sizeof(int)),1201);
    //print_CUDA_error_if_any(cudaMemcpy(mp->d_nibool_interfaces_ext_mesh,h_nibool_interfaces_ext_mesh,
    //                                   (mp->num_interfaces_ext_mesh)*sizeof(int),cudaMemcpyHostToDevice),1202);

    copy_todevice_int((void**)&mp->d_ibool_interfaces_ext_mesh,h_ibool_interfaces_ext_mesh,
                      (mp->num_interfaces_ext_mesh)*(mp->max_nibool_interfaces_ext_mesh));

    //print_CUDA_error_if_any(cudaMalloc((void**) &mp->d_ibool_interfaces_ext_mesh,
    //                                   (mp->num_interfaces_ext_mesh)*(mp->max_nibool_interfaces_ext_mesh)*sizeof(int)),1203);
    //print_CUDA_error_if_any(cudaMemcpy(mp->d_ibool_interfaces_ext_mesh,h_ibool_interfaces_ext_mesh,
    //                                   (mp->num_interfaces_ext_mesh)*(mp->max_nibool_interfaces_ext_mesh)*sizeof(int),
    //                                   cudaMemcpyHostToDevice),1204);
  }

  // Allocate pinned mpi-buffers.
  // MPI buffers use pinned memory allocated by cudaMallocHost, which
  // enables the use of asynchronous memory copies from host <->
  // device
  int size_mpi_buffer = 3 * (mp->num_interfaces_ext_mesh) * (mp->max_nibool_interfaces_ext_mesh);
  // send buffer
  mp->size_mpi_buffer = size_mpi_buffer;
  if( mp->size_mpi_buffer > 0 ){
    print_CUDA_error_if_any(cudaMallocHost((void**)&(mp->h_send_accel_buffer),sizeof(float)*(mp->size_mpi_buffer)),8004);
    mp->send_buffer = (float*)malloc((mp->size_mpi_buffer)*sizeof(float));
    // adjoint
    print_CUDA_error_if_any(cudaMallocHost((void**)&(mp->h_send_b_accel_buffer),sizeof(float)*(size_mpi_buffer)),8004);
    // mp->b_send_buffer = (float*)malloc((size_mpi_buffer)*sizeof(float));

    // receive buffer
    print_CUDA_error_if_any(cudaMallocHost((void**)&(mp->h_recv_accel_buffer),sizeof(float)*(mp->size_mpi_buffer)),8004);
    mp->recv_buffer = (float*)malloc((mp->size_mpi_buffer)*sizeof(float));
  }

  //daniel: check if needed
  //mp->nibool_interfaces_ext_mesh = h_nibool_interfaces_ext_mesh;
  //mp->my_neighbours_ext_mesh = my_neighbours_ext_mesh;
  //mp->request_send_vector_ext_mesh = request_send_vector_ext_mesh;
  //mp->request_recv_vector_ext_mesh = request_recv_vector_ext_mesh;
  //mp->buffer_recv_vector_ext_mesh = buffer_recv_vector_ext_mesh;

  // setup two streams, one for compute and one for host<->device memory copies
  // compute stream
  cudaStreamCreate(&mp->compute_stream);
  // copy stream (needed to transfer mpi buffers)
  if( mp->size_mpi_buffer > 0 ){
    cudaStreamCreate(&mp->copy_stream);
    //cudaStreamCreate(&mp->b_copy_stream);
  }

  // inner elements
  copy_todevice_int((void**)&mp->d_ispec_is_inner,h_ispec_is_inner,mp->NSPEC_AB);

//  print_CUDA_error_if_any(cudaMalloc((void**) &mp->d_ispec_is_inner,mp->NSPEC_AB*sizeof(int)),1205);
//  print_CUDA_error_if_any(cudaMemcpy(mp->d_ispec_is_inner, h_ispec_is_inner,
 //                                    mp->NSPEC_AB*sizeof(int),cudaMemcpyHostToDevice),1206);

  // absorbing boundaries
  mp->d_num_abs_boundary_faces = *h_num_abs_boundary_faces;
  if( mp->absorbing_conditions && mp->d_num_abs_boundary_faces > 0 ){
    copy_todevice_int((void**)&mp->d_abs_boundary_ispec,h_abs_boundary_ispec,mp->d_num_abs_boundary_faces);

//    print_CUDA_error_if_any(cudaMalloc((void**) &(mp->d_abs_boundary_ispec),
//                                       (mp->d_num_abs_boundary_faces)*sizeof(int)),1101);
//    print_CUDA_error_if_any(cudaMemcpy(mp->d_abs_boundary_ispec, h_abs_boundary_ispec,
//                                       (mp->d_num_abs_boundary_faces)*sizeof(int),
//                                       cudaMemcpyHostToDevice),1102);

    copy_todevice_int((void**)&mp->d_abs_boundary_ijk,h_abs_boundary_ijk,
                      3*NGLL2*(mp->d_num_abs_boundary_faces));

//    print_CUDA_error_if_any(cudaMalloc((void**) &(mp->d_abs_boundary_ijk),
//                                       3*NGLL2*(mp->d_num_abs_boundary_faces)*sizeof(int)),1103);
//    print_CUDA_error_if_any(cudaMemcpy(mp->d_abs_boundary_ijk, h_abs_boundary_ijk,
//                                       3*NGLL2*(mp->d_num_abs_boundary_faces)*sizeof(int),
//                                       cudaMemcpyHostToDevice),1104);

    copy_todevice_realw((void**)&mp->d_abs_boundary_normal,h_abs_boundary_normal,
                        NDIM*NGLL2*(mp->d_num_abs_boundary_faces));

//    print_CUDA_error_if_any(cudaMalloc((void**) &(mp->d_abs_boundary_normal),
//                                       3*NGLL2*(mp->d_num_abs_boundary_faces)*sizeof(realw)),1105);
//    print_CUDA_error_if_any(cudaMemcpy(mp->d_abs_boundary_normal, h_abs_boundary_normal,
//                                       3*NGLL2*(mp->d_num_abs_boundary_faces)*sizeof(realw),
//                                       cudaMemcpyHostToDevice),1106);

    copy_todevice_realw((void**)&mp->d_abs_boundary_jacobian2Dw,h_abs_boundary_jacobian2Dw,
                        NGLL2*(mp->d_num_abs_boundary_faces));

//    print_CUDA_error_if_any(cudaMalloc((void**) &(mp->d_abs_boundary_jacobian2Dw),
//                                       NGLL2*(mp->d_num_abs_boundary_faces)*sizeof(realw)),1107);
//    print_CUDA_error_if_any(cudaMemcpy(mp->d_abs_boundary_jacobian2Dw, h_abs_boundary_jacobian2Dw,
//                                       NGLL2*(mp->d_num_abs_boundary_faces)*sizeof(realw),
 //                                      cudaMemcpyHostToDevice),1108);
  }

  // sources
  mp->nsources_local = *nsources_local_f;
  if (mp->simulation_type == 1  || mp->simulation_type == 3){
    // not needed in case of pure adjoint simulations (SIMULATION_TYPE == 2)
    copy_todevice_realw((void**)&mp->d_sourcearrays,h_sourcearrays,(*NSOURCES)*NDIM*NGLL3);

    // buffer for source time function values
    print_CUDA_error_if_any(cudaMalloc((void**)&mp->d_stf_pre_compute,
                                       *NSOURCES*sizeof(double)),1303);
  }
  copy_todevice_int((void**)&mp->d_islice_selected_source,h_islice_selected_source,(*NSOURCES));
  copy_todevice_int((void**)&mp->d_ispec_selected_source,h_ispec_selected_source,(*NSOURCES));


  // receiver stations
  mp->nrec_local = *nrec_local; // number of receiver located in this partition
  // note that:
  // size(number_receiver_global) = nrec_local
  // size(ispec_selected_rec) = nrec
  if( mp->nrec_local > 0 ){
    copy_todevice_int((void**)&mp->d_number_receiver_global,h_number_receiver_global,mp->nrec_local);
  }
  copy_todevice_int((void**)&mp->d_ispec_selected_rec,h_ispec_selected_rec,(*nrec));


#ifdef USE_MESH_COLORING_GPU
  mp->use_mesh_coloring_gpu = 1;
  if( ! *USE_MESH_COLORING_GPU_f ) exit_on_error("error with USE_MESH_COLORING_GPU constant; please re-compile\n");
#else
  // mesh coloring
  // note: this here passes the coloring as an option to the kernel routines
  //          the performance seems to be the same if one uses the pre-processing directives above or not
  mp->use_mesh_coloring_gpu = *USE_MESH_COLORING_GPU_f;
#endif

  // number of elements per domain
  mp->nspec_acoustic = *nspec_acoustic;
  mp->nspec_elastic = *nspec_elastic;

  // gravity flag initialization
  mp->gravity = 0;

  // JC JC here we will need to add GPU support for the new C-PML routines

#ifdef ENABLE_VERY_SLOW_ERROR_CHECKING
  exit_on_cuda_error("prepare_constants_device");
#endif
}


/* ----------------------------------------------------------------------------------------------- */

// for ACOUSTIC simulations

/* ----------------------------------------------------------------------------------------------- */

extern "C"
void FC_FUNC_(prepare_fields_acoustic_device,
              PREPARE_FIELDS_ACOUSTIC_DEVICE)(long* Mesh_pointer_f,
                                              realw* rmass_acoustic,
                                              realw* rhostore,
                                              realw* kappastore,
                                              int* num_phase_ispec_acoustic,
                                              int* phase_ispec_inner_acoustic,
                                              int* ispec_is_acoustic,
                                              int* NOISE_TOMOGRAPHY,
                                              int* num_free_surface_faces,
                                              int* free_surface_ispec,
                                              int* free_surface_ijk,
                                              int* b_reclen_potential,
                                              realw* b_absorb_potential,
                                              int* ELASTIC_SIMULATION,
                                              int* num_coupling_ac_el_faces,
                                              int* coupling_ac_el_ispec,
                                              int* coupling_ac_el_ijk,
                                              realw* coupling_ac_el_normal,
                                              realw* coupling_ac_el_jacobian2Dw,
                                              int* num_colors_outer_acoustic,
                                              int* num_colors_inner_acoustic,
                                              int* num_elem_colors_acoustic) {

  TRACE("prepare_fields_acoustic_device");

  Mesh* mp = (Mesh*)(*Mesh_pointer_f);
  /* Assuming NGLLX==5. Padded is then 128 (5^3+3) */
  int size_padded = NGLL3_PADDED * mp->NSPEC_AB;
//  int size_nonpadded = NGLL3 * mp->NSPEC_AB;
  int size_glob = mp->NGLOB_AB;

  // allocates arrays on device (GPU)
  print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_potential_acoustic),sizeof(realw)*size_glob),2001);
  print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_potential_dot_acoustic),sizeof(realw)*size_glob),2002);
  print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_potential_dot_dot_acoustic),sizeof(realw)*size_glob),2003);

  // mpi buffer
  if( mp->num_interfaces_ext_mesh > 0 ){
    print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_send_potential_dot_dot_buffer),
                      (mp->max_nibool_interfaces_ext_mesh)*(mp->num_interfaces_ext_mesh)*sizeof(realw)),2004);
  }

  // mass matrix
  copy_todevice_realw((void**)&mp->d_rmass_acoustic,rmass_acoustic,mp->NGLOB_AB);

//  print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_rmass_acoustic),sizeof(realw)*size_glob),2005);
//  print_CUDA_error_if_any(cudaMemcpy(mp->d_rmass_acoustic,rmass_acoustic,
//                                     sizeof(realw)*size_glob,cudaMemcpyHostToDevice),2100);

  // padded array
  print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_rhostore),size_padded*sizeof(realw)),2006);
  // transfer constant element data with padding
  for(int i=0; i < mp->NSPEC_AB; i++) {
    print_CUDA_error_if_any(cudaMemcpy(mp->d_rhostore+i*NGLL3_PADDED, &rhostore[i*NGLL3],
                                       NGLL3*sizeof(realw),cudaMemcpyHostToDevice),2106);
  }

  // non-padded array
  copy_todevice_realw((void**)&mp->d_kappastore,kappastore,NGLL3*mp->NSPEC_AB);

//  print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_kappastore),size_nonpadded*sizeof(realw)),2007);
//  print_CUDA_error_if_any(cudaMemcpy(mp->d_kappastore,kappastore,
//                                     NGLL3*mp->NSPEC_AB*sizeof(realw),cudaMemcpyHostToDevice),2105);

  // phase elements
  mp->num_phase_ispec_acoustic = *num_phase_ispec_acoustic;

  copy_todevice_int((void**)&mp->d_phase_ispec_inner_acoustic,phase_ispec_inner_acoustic,
                    2*mp->num_phase_ispec_acoustic);

//  print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_phase_ispec_inner_acoustic),
//                                      mp->num_phase_ispec_acoustic*2*sizeof(int)),2008);
//  print_CUDA_error_if_any(cudaMemcpy(mp->d_phase_ispec_inner_acoustic,phase_ispec_inner_acoustic,
//                                     mp->num_phase_ispec_acoustic*2*sizeof(int),cudaMemcpyHostToDevice),2101);

  copy_todevice_int((void**)&mp->d_ispec_is_acoustic,ispec_is_acoustic,mp->NSPEC_AB);

//  print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_ispec_is_acoustic),
//                                     mp->NSPEC_AB*sizeof(int)),2009);
//  print_CUDA_error_if_any(cudaMemcpy(mp->d_ispec_is_acoustic,ispec_is_acoustic,
//                                     mp->NSPEC_AB*sizeof(int),cudaMemcpyHostToDevice),2102);

  // free surface
  if( *NOISE_TOMOGRAPHY == 0 ){
    // allocate surface arrays
    mp->num_free_surface_faces = *num_free_surface_faces;
    if( mp->num_free_surface_faces > 0 ){

      copy_todevice_int((void**)&mp->d_free_surface_ispec,free_surface_ispec,mp->num_free_surface_faces);

//      print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_free_surface_ispec),
//                                       mp->num_free_surface_faces*sizeof(int)),2201);
//      print_CUDA_error_if_any(cudaMemcpy(mp->d_free_surface_ispec,free_surface_ispec,
//                                       mp->num_free_surface_faces*sizeof(int),cudaMemcpyHostToDevice),2203);


      copy_todevice_int((void**)&mp->d_free_surface_ijk,free_surface_ijk,
                        3*NGLL2*mp->num_free_surface_faces);

//      print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_free_surface_ijk),
//                                       3*NGLL2*mp->num_free_surface_faces*sizeof(int)),2202);
//      print_CUDA_error_if_any(cudaMemcpy(mp->d_free_surface_ijk,free_surface_ijk,
//                                       3*NGLL2*mp->num_free_surface_faces*sizeof(int),cudaMemcpyHostToDevice),2204);
    }
  }

  // absorbing boundaries
  if( mp->absorbing_conditions ){
    mp->d_b_reclen_potential = *b_reclen_potential;

    copy_todevice_realw((void**)&mp->d_b_absorb_potential,b_absorb_potential,mp->d_b_reclen_potential);

//    print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_b_absorb_potential),mp->d_b_reclen_potential),2301);
//    print_CUDA_error_if_any(cudaMemcpy(mp->d_b_absorb_potential,b_absorb_potential,
//                                       mp->d_b_reclen_potential,cudaMemcpyHostToDevice),2302);
  }


  // for seismograms
  if( mp->nrec_local > 0 ){
    print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_station_seismo_potential),
                                       mp->nrec_local*NGLL3*sizeof(realw)),2107);

    mp->h_station_seismo_potential = (realw*) malloc( mp->nrec_local*NGLL3*sizeof(realw) );
    if( mp->h_station_seismo_potential == NULL) exit_on_error("error allocating h_station_seismo_potential");
  }


  // coupling with elastic parts
  if( *ELASTIC_SIMULATION && *num_coupling_ac_el_faces > 0 ){

    copy_todevice_int((void**)&mp->d_coupling_ac_el_ispec,coupling_ac_el_ispec,(*num_coupling_ac_el_faces));

//    print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_coupling_ac_el_ispec),
//                                       (*num_coupling_ac_el_faces)*sizeof(int)),2601);
//    print_CUDA_error_if_any(cudaMemcpy(mp->d_coupling_ac_el_ispec,coupling_ac_el_ispec,
//                                       (*num_coupling_ac_el_faces)*sizeof(int),cudaMemcpyHostToDevice),2602);

    copy_todevice_int((void**)&mp->d_coupling_ac_el_ijk,coupling_ac_el_ijk,3*NGLL2*(*num_coupling_ac_el_faces));

//    print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_coupling_ac_el_ijk),
//                                       3*NGLL2*(*num_coupling_ac_el_faces)*sizeof(int)),2603);
//    print_CUDA_error_if_any(cudaMemcpy(mp->d_coupling_ac_el_ijk,coupling_ac_el_ijk,
//                                       3*NGLL2*(*num_coupling_ac_el_faces)*sizeof(int),cudaMemcpyHostToDevice),2604);

    copy_todevice_realw((void**)&mp->d_coupling_ac_el_normal,coupling_ac_el_normal,
                        3*NGLL2*(*num_coupling_ac_el_faces));

//    print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_coupling_ac_el_normal),
//                                        3*NGLL2*(*num_coupling_ac_el_faces)*sizeof(realw)),2605);
//    print_CUDA_error_if_any(cudaMemcpy(mp->d_coupling_ac_el_normal,coupling_ac_el_normal,
//                                        3*NGLL2*(*num_coupling_ac_el_faces)*sizeof(realw),cudaMemcpyHostToDevice),2606);

    copy_todevice_realw((void**)&mp->d_coupling_ac_el_jacobian2Dw,coupling_ac_el_jacobian2Dw,
                        NGLL2*(*num_coupling_ac_el_faces));

//    print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_coupling_ac_el_jacobian2Dw),
//                                        NGLL2*(*num_coupling_ac_el_faces)*sizeof(realw)),2607);
//    print_CUDA_error_if_any(cudaMemcpy(mp->d_coupling_ac_el_jacobian2Dw,coupling_ac_el_jacobian2Dw,
//                                        NGLL2*(*num_coupling_ac_el_faces)*sizeof(realw),cudaMemcpyHostToDevice),2608);

  }

  // mesh coloring
  if( mp->use_mesh_coloring_gpu ){
    mp->num_colors_outer_acoustic = *num_colors_outer_acoustic;
    mp->num_colors_inner_acoustic = *num_colors_inner_acoustic;
    mp->h_num_elem_colors_acoustic = (int*) num_elem_colors_acoustic;
  }

#ifdef ENABLE_VERY_SLOW_ERROR_CHECKING
  exit_on_cuda_error("prepare_fields_acoustic_device");
#endif
}


/* ----------------------------------------------------------------------------------------------- */

extern "C"
void FC_FUNC_(prepare_fields_acoustic_adj_dev,
              PREPARE_FIELDS_ACOUSTIC_ADJ_DEV)(long* Mesh_pointer_f,
                                              int* APPROXIMATE_HESS_KL) {

  TRACE("prepare_fields_acoustic_adj_dev");

  Mesh* mp = (Mesh*)(*Mesh_pointer_f);

  int size_glob = mp->NGLOB_AB;

  // kernel simulations
  if( mp->simulation_type != 3 ) return;

  // allocates backward/reconstructed arrays on device (GPU)
  print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_b_potential_acoustic),sizeof(realw)*size_glob),3014);
  print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_b_potential_dot_acoustic),sizeof(realw)*size_glob),3015);
  print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_b_potential_dot_dot_acoustic),sizeof(realw)*size_glob),3016);

  // allocates kernels
  print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_rho_ac_kl),NGLL3*mp->NSPEC_AB*sizeof(realw)),3017);
  print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_kappa_ac_kl),NGLL3*mp->NSPEC_AB*sizeof(realw)),3018);

  // initializes kernel values to zero
  print_CUDA_error_if_any(cudaMemset(mp->d_rho_ac_kl,0,
                                     NGLL3*mp->NSPEC_AB*sizeof(realw)),3019);
  print_CUDA_error_if_any(cudaMemset(mp->d_kappa_ac_kl,0,
                                     NGLL3*mp->NSPEC_AB*sizeof(realw)),3020);

  // preconditioner
  if( *APPROXIMATE_HESS_KL ){
    print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_hess_ac_kl),NGLL3*mp->NSPEC_AB*sizeof(realw)),3030);
    // initializes with zeros
    print_CUDA_error_if_any(cudaMemset(mp->d_hess_ac_kl,0,
                                       NGLL3*mp->NSPEC_AB*sizeof(realw)),3031);
  }

#ifdef ENABLE_VERY_SLOW_ERROR_CHECKING
  exit_on_cuda_error("prepare_fields_acoustic_adj_dev");
#endif
}


/* ----------------------------------------------------------------------------------------------- */

// for ELASTIC simulations

/* ----------------------------------------------------------------------------------------------- */

extern "C"
void FC_FUNC_(prepare_fields_elastic_device,
              PREPARE_FIELDS_ELASTIC_DEVICE)(long* Mesh_pointer_f,
                                             int* size,
                                             realw* rmassx,
                                             realw* rmassy,
                                             realw* rmassz,
                                             realw* rho_vp,
                                             realw* rho_vs,
                                             int* num_phase_ispec_elastic,
                                             int* phase_ispec_inner_elastic,
                                             int* ispec_is_elastic,
                                             realw* h_b_absorb_field,
                                             int* h_b_reclen_field,
                                             int* SAVE_FORWARD,
                                             int* COMPUTE_AND_STORE_STRAIN,
                                             realw* epsilondev_xx,realw* epsilondev_yy,realw* epsilondev_xy,
                                             realw* epsilondev_xz,realw* epsilondev_yz,
                                             int* ATTENUATION,
                                             int* R_size,
                                             realw* R_xx,realw* R_yy,realw* R_xy,realw* R_xz,realw* R_yz,
                                             realw* one_minus_sum_beta,realw* factor_common,
                                             realw* alphaval,realw* betaval,realw* gammaval,
                                             int* APPROXIMATE_OCEAN_LOAD,
                                             realw* rmass_ocean_load,
                                             int* NOISE_TOMOGRAPHY,
                                             realw* free_surface_normal,
                                             int* free_surface_ispec,
                                             int* free_surface_ijk,
                                             int* num_free_surface_faces,
                                             int* ACOUSTIC_SIMULATION,
                                             int* num_colors_outer_elastic,
                                             int* num_colors_inner_elastic,
                                             int* num_elem_colors_elastic,
                                             int* ANISOTROPY,
                                             realw *c11store,
                                             realw *c12store,
                                             realw *c13store,
                                             realw *c14store,
                                             realw *c15store,
                                             realw *c16store,
                                             realw *c22store,
                                             realw *c23store,
                                             realw *c24store,
                                             realw *c25store,
                                             realw *c26store,
                                             realw *c33store,
                                             realw *c34store,
                                             realw *c35store,
                                             realw *c36store,
                                             realw *c44store,
                                             realw *c45store,
                                             realw *c46store,
                                             realw *c55store,
                                             realw *c56store,
                                             realw *c66store){

TRACE("prepare_fields_elastic_device");

  Mesh* mp = (Mesh*)(*Mesh_pointer_f);
  /* Assuming NGLLX==5. Padded is then 128 (5^3+3) */
  int size_padded = NGLL3_PADDED * (mp->NSPEC_AB);
//  int size_nonpadded = NGLL3 * (mp->NSPEC_AB);

  print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_displ),sizeof(realw)*(*size)),4001);
  print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_veloc),sizeof(realw)*(*size)),4002);
  print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_accel),sizeof(realw)*(*size)),4003);

  #ifdef USE_TEXTURES_FIELDS
  {
    #ifdef USE_OLDER_CUDA4_GPU
      print_CUDA_error_if_any(cudaGetTextureReference(&mp->d_displ_tex_ref_ptr, "d_displ_tex"), 4001);
      cudaChannelFormatDesc channelDesc = cudaCreateChannelDesc<float>();
      print_CUDA_error_if_any(cudaBindTexture(0, mp->d_displ_tex_ref_ptr, mp->d_displ, &channelDesc, sizeof(realw)*(*size)), 4001);
    #else
      cudaChannelFormatDesc channelDesc = cudaCreateChannelDesc<float>();
      print_CUDA_error_if_any(cudaBindTexture(0, &d_displ_tex, mp->d_displ, &channelDesc, sizeof(realw)*(*size)), 4001);
    #endif
  }

  {
    #ifdef USE_OLDER_CUDA4_GPU
      print_CUDA_error_if_any(cudaGetTextureReference(&mp->d_veloc_tex_ref_ptr, "d_veloc_tex"), 4002);
      cudaChannelFormatDesc channelDesc = cudaCreateChannelDesc<float>();
      print_CUDA_error_if_any(cudaBindTexture(0, mp->d_veloc_tex_ref_ptr, mp->d_veloc, &channelDesc, sizeof(realw)*(*size)), 4002);
    #else
      cudaChannelFormatDesc channelDesc = cudaCreateChannelDesc<float>();
      print_CUDA_error_if_any(cudaBindTexture(0, &d_veloc_tex, mp->d_veloc, &channelDesc, sizeof(realw)*(*size)), 4002);
    #endif
  }

  {
    #ifdef USE_OLDER_CUDA4_GPU
      print_CUDA_error_if_any(cudaGetTextureReference(&mp->d_accel_tex_ref_ptr, "d_accel_tex"), 4003);
      cudaChannelFormatDesc channelDesc = cudaCreateChannelDesc<float>();
      print_CUDA_error_if_any(cudaBindTexture(0, mp->d_accel_tex_ref_ptr, mp->d_accel, &channelDesc, sizeof(realw)*(*size)), 4003);
    #else
      cudaChannelFormatDesc channelDesc = cudaCreateChannelDesc<float>();
      print_CUDA_error_if_any(cudaBindTexture(0, &d_accel_tex, mp->d_accel, &channelDesc, sizeof(realw)*(*size)), 4003);
    #endif
  }
  #endif

  // mpi buffer
  if( mp->size_mpi_buffer > 0 ){
    print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_send_accel_buffer),
                                      mp->size_mpi_buffer*sizeof(realw)),4004);
  }

  // mass matrix
  copy_todevice_realw((void**)&mp->d_rmassx,rmassx,mp->NGLOB_AB);
  copy_todevice_realw((void**)&mp->d_rmassy,rmassy,mp->NGLOB_AB);
  copy_todevice_realw((void**)&mp->d_rmassz,rmassz,mp->NGLOB_AB);

//  print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_rmass),sizeof(realw)*mp->NGLOB_AB),4005);
//  // transfer element data
//  print_CUDA_error_if_any(cudaMemcpy(mp->d_rmass,rmass,
//                                     sizeof(realw)*mp->NGLOB_AB,cudaMemcpyHostToDevice),4010);


  // element indices
  copy_todevice_int((void**)&mp->d_ispec_is_elastic,ispec_is_elastic,mp->NSPEC_AB);

//  print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_ispec_is_elastic),mp->NSPEC_AB*sizeof(int)),4009);
//  print_CUDA_error_if_any(cudaMemcpy(mp->d_ispec_is_elastic,ispec_is_elastic,
//                                     mp->NSPEC_AB*sizeof(int),cudaMemcpyHostToDevice),4012);

  // phase elements
  mp->num_phase_ispec_elastic = *num_phase_ispec_elastic;

  copy_todevice_int((void**)&mp->d_phase_ispec_inner_elastic,phase_ispec_inner_elastic,2*mp->num_phase_ispec_elastic);

//  print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_phase_ispec_inner_elastic),
//                                     mp->num_phase_ispec_elastic*2*sizeof(int)),4008);
//  print_CUDA_error_if_any(cudaMemcpy(mp->d_phase_ispec_inner_elastic,phase_ispec_inner_elastic,
//                                     mp->num_phase_ispec_elastic*2*sizeof(int),cudaMemcpyHostToDevice),4011);

  // for seismograms
  if( mp->nrec_local > 0 ){
    print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_station_seismo_field),
                                     3*NGLL3*(mp->nrec_local)*sizeof(realw)),4015);

    mp->h_station_seismo_field = (realw*) malloc( 3*NGLL3*(mp->nrec_local)*sizeof(realw) );
    if( mp->h_station_seismo_field == NULL) exit_on_error("h_station_seismo_field not allocated \n");
  }

  // absorbing conditions
  if( mp->absorbing_conditions && mp->d_num_abs_boundary_faces > 0){
    // non-padded arrays
    copy_todevice_realw((void**)&mp->d_rho_vp,rho_vp,NGLL3*mp->NSPEC_AB);
    copy_todevice_realw((void**)&mp->d_rho_vs,rho_vs,NGLL3*mp->NSPEC_AB);

//    print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_rho_vp),size_nonpadded*sizeof(realw)),4006);
//    print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_rho_vs),size_nonpadded*sizeof(realw)),4007);

    // rho_vp, rho_vs non-padded; they are needed for stacey boundary condition
//    print_CUDA_error_if_any(cudaMemcpy(mp->d_rho_vp, rho_vp,
//                                       NGLL3*mp->NSPEC_AB*sizeof(realw),cudaMemcpyHostToDevice),4013);
//    print_CUDA_error_if_any(cudaMemcpy(mp->d_rho_vs, rho_vs,
//                                       NGLL3*mp->NSPEC_AB*sizeof(realw),cudaMemcpyHostToDevice),4014);

    // absorb_field array used for file i/o
    if(mp->simulation_type == 3 || ( mp->simulation_type == 1 && *SAVE_FORWARD )){
      mp->d_b_reclen_field = *h_b_reclen_field;

      copy_todevice_realw((void**)&mp->d_b_absorb_field,h_b_absorb_field,mp->d_b_reclen_field);

//      print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_b_absorb_field),
//                                       mp->d_b_reclen_field),4016);
//      print_CUDA_error_if_any(cudaMemcpy(mp->d_b_absorb_field, h_b_absorb_field,
//                                       mp->d_b_reclen_field,cudaMemcpyHostToDevice),4017);
    }
  }

  // strains used for attenuation and kernel simulations
  if( *COMPUTE_AND_STORE_STRAIN ){
    // strains
    int epsilondev_size = NGLL3*mp->NSPEC_AB; // note: non-aligned; if align, check memcpy below and indexing

    copy_todevice_realw((void**)&mp->d_epsilondev_xx,epsilondev_xx,epsilondev_size);

//    print_CUDA_error_if_any(cudaMalloc((void**)&mp->d_epsilondev_xx,
//                                       epsilondev_size*sizeof(realw)),4301);
//    print_CUDA_error_if_any(cudaMemcpy(mp->d_epsilondev_xx,epsilondev_xx,epsilondev_size*sizeof(realw),
//                                       cudaMemcpyHostToDevice),4302);

    copy_todevice_realw((void**)&mp->d_epsilondev_yy,epsilondev_yy,epsilondev_size);

//    print_CUDA_error_if_any(cudaMalloc((void**)&mp->d_epsilondev_yy,
//                                       epsilondev_size*sizeof(realw)),4302);
//    print_CUDA_error_if_any(cudaMemcpy(mp->d_epsilondev_yy,epsilondev_yy,epsilondev_size*sizeof(realw),
//                                       cudaMemcpyHostToDevice),4303);

    copy_todevice_realw((void**)&mp->d_epsilondev_xy,epsilondev_xy,epsilondev_size);

//    print_CUDA_error_if_any(cudaMalloc((void**)&mp->d_epsilondev_xy,
//                                       epsilondev_size*sizeof(realw)),4304);
//    print_CUDA_error_if_any(cudaMemcpy(mp->d_epsilondev_xy,epsilondev_xy,epsilondev_size*sizeof(realw),
//                                       cudaMemcpyHostToDevice),4305);

    copy_todevice_realw((void**)&mp->d_epsilondev_xz,epsilondev_xz,epsilondev_size);

//    print_CUDA_error_if_any(cudaMalloc((void**)&mp->d_epsilondev_xz,
//                                       epsilondev_size*sizeof(realw)),4306);
//    print_CUDA_error_if_any(cudaMemcpy(mp->d_epsilondev_xz,epsilondev_xz,epsilondev_size*sizeof(realw),
//                                       cudaMemcpyHostToDevice),4307);

    copy_todevice_realw((void**)&mp->d_epsilondev_yz,epsilondev_yz,epsilondev_size);

//    print_CUDA_error_if_any(cudaMalloc((void**)&mp->d_epsilondev_yz,
//                                       epsilondev_size*sizeof(realw)),4308);
//    print_CUDA_error_if_any(cudaMemcpy(mp->d_epsilondev_yz,epsilondev_yz,epsilondev_size*sizeof(realw),
//                                       cudaMemcpyHostToDevice),4309);

  }

  // attenuation memory variables
  if( *ATTENUATION ){
    // memory arrays
    copy_todevice_realw((void**)&mp->d_R_xx,R_xx,(*R_size));

//    print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_R_xx),
//                                       (*R_size)*sizeof(realw)),4401);
//    print_CUDA_error_if_any(cudaMemcpy(mp->d_R_xx,R_xx,(*R_size)*sizeof(realw),
//                                       cudaMemcpyHostToDevice),4402);

    copy_todevice_realw((void**)&mp->d_R_yy,R_yy,(*R_size));

//    print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_R_yy),
//                                       (*R_size)*sizeof(realw)),4403);
//    print_CUDA_error_if_any(cudaMemcpy(mp->d_R_yy,R_yy,(*R_size)*sizeof(realw),
//                                       cudaMemcpyHostToDevice),4404);

    copy_todevice_realw((void**)&mp->d_R_xy,R_xy,(*R_size));

//    print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_R_xy),
//                                       (*R_size)*sizeof(realw)),4405);
//    print_CUDA_error_if_any(cudaMemcpy(mp->d_R_xy,R_xy,(*R_size)*sizeof(realw),
//                                       cudaMemcpyHostToDevice),4406);

    copy_todevice_realw((void**)&mp->d_R_xz,R_xz,(*R_size));

//    print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_R_xz),
//                                       (*R_size)*sizeof(realw)),4407);
//    print_CUDA_error_if_any(cudaMemcpy(mp->d_R_xz,R_xz,(*R_size)*sizeof(realw),
//                                       cudaMemcpyHostToDevice),4408);

    copy_todevice_realw((void**)&mp->d_R_yz,R_yz,(*R_size));

//    print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_R_yz),
//                                       (*R_size)*sizeof(realw)),4409);
//    print_CUDA_error_if_any(cudaMemcpy(mp->d_R_yz,R_yz,(*R_size)*sizeof(realw),
//                                       cudaMemcpyHostToDevice),4410);

    // attenuation factors
    copy_todevice_realw((void**)&mp->d_one_minus_sum_beta,one_minus_sum_beta,NGLL3*mp->NSPEC_AB);

//    print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_one_minus_sum_beta),
//                                       NGLL3*mp->NSPEC_AB*sizeof(realw)),4430);
//    print_CUDA_error_if_any(cudaMemcpy(mp->d_one_minus_sum_beta ,one_minus_sum_beta,
//                                       NGLL3*mp->NSPEC_AB*sizeof(realw),cudaMemcpyHostToDevice),4431);

    copy_todevice_realw((void**)&mp->d_factor_common,factor_common,N_SLS*NGLL3*mp->NSPEC_AB);

//    print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_factor_common),
//                                       N_SLS*NGLL3*mp->NSPEC_AB*sizeof(realw)),4432);
//    print_CUDA_error_if_any(cudaMemcpy(mp->d_factor_common ,factor_common,
//                                       N_SLS*NGLL3*mp->NSPEC_AB*sizeof(realw),cudaMemcpyHostToDevice),4433);

    // alpha,beta,gamma factors
    copy_todevice_realw((void**)&mp->d_alphaval,alphaval,N_SLS);

//    print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_alphaval),
//                                       N_SLS*sizeof(realw)),4434);
//    print_CUDA_error_if_any(cudaMemcpy(mp->d_alphaval ,alphaval,
//                                       N_SLS*sizeof(realw),cudaMemcpyHostToDevice),4435);

    copy_todevice_realw((void**)&mp->d_betaval,betaval,N_SLS);

//    print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_betaval),
//                                       N_SLS*sizeof(realw)),4436);
//    print_CUDA_error_if_any(cudaMemcpy(mp->d_betaval ,betaval,
//                                       N_SLS*sizeof(realw),cudaMemcpyHostToDevice),4437);

    copy_todevice_realw((void**)&mp->d_gammaval,gammaval,N_SLS);

//    print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_gammaval),
//                                       N_SLS*sizeof(realw)),4438);
//    print_CUDA_error_if_any(cudaMemcpy(mp->d_gammaval ,gammaval,
//                                       N_SLS*sizeof(realw),cudaMemcpyHostToDevice),4439);

  }

  // anisotropy
  if( *ANISOTROPY ){
    // allocates memory on GPU
    print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_c11store),
                                       size_padded*sizeof(realw)),4700);
    print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_c12store),
                                       size_padded*sizeof(realw)),4701);
    print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_c13store),
                                       size_padded*sizeof(realw)),4702);
    print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_c14store),
                                       size_padded*sizeof(realw)),4703);
    print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_c15store),
                                       size_padded*sizeof(realw)),4704);
    print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_c16store),
                                       size_padded*sizeof(realw)),4705);
    print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_c22store),
                                       size_padded*sizeof(realw)),4706);
    print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_c23store),
                                       size_padded*sizeof(realw)),4707);
    print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_c24store),
                                       size_padded*sizeof(realw)),4708);
    print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_c25store),
                                       size_padded*sizeof(realw)),4709);
    print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_c26store),
                                       size_padded*sizeof(realw)),4710);
    print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_c33store),
                                       size_padded*sizeof(realw)),4711);
    print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_c34store),
                                       size_padded*sizeof(realw)),4712);
    print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_c35store),
                                       size_padded*sizeof(realw)),4713);
    print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_c36store),
                                       size_padded*sizeof(realw)),4714);
    print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_c44store),
                                       size_padded*sizeof(realw)),4715);
    print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_c45store),
                                       size_padded*sizeof(realw)),4716);
    print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_c46store),
                                       size_padded*sizeof(realw)),4717);
    print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_c55store),
                                       size_padded*sizeof(realw)),4718);
    print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_c56store),
                                       size_padded*sizeof(realw)),4719);
    print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_c66store),
                                       size_padded*sizeof(realw)),4720);

    // transfer constant element data with padding
    for(int i=0;i < mp->NSPEC_AB;i++) {
      print_CUDA_error_if_any(cudaMemcpy(mp->d_c11store + i*NGLL3_PADDED, &c11store[i*NGLL3],
                                         NGLL3*sizeof(realw),cudaMemcpyHostToDevice),4800);
      print_CUDA_error_if_any(cudaMemcpy(mp->d_c12store + i*NGLL3_PADDED, &c12store[i*NGLL3],
                                         NGLL3*sizeof(realw),cudaMemcpyHostToDevice),4801);
      print_CUDA_error_if_any(cudaMemcpy(mp->d_c13store + i*NGLL3_PADDED, &c13store[i*NGLL3],
                                         NGLL3*sizeof(realw),cudaMemcpyHostToDevice),4802);
      print_CUDA_error_if_any(cudaMemcpy(mp->d_c14store + i*NGLL3_PADDED, &c14store[i*NGLL3],
                                         NGLL3*sizeof(realw),cudaMemcpyHostToDevice),4803);
      print_CUDA_error_if_any(cudaMemcpy(mp->d_c15store + i*NGLL3_PADDED, &c15store[i*NGLL3],
                                         NGLL3*sizeof(realw),cudaMemcpyHostToDevice),4804);
      print_CUDA_error_if_any(cudaMemcpy(mp->d_c16store + i*NGLL3_PADDED, &c16store[i*NGLL3],
                                         NGLL3*sizeof(realw),cudaMemcpyHostToDevice),4805);
      print_CUDA_error_if_any(cudaMemcpy(mp->d_c22store + i*NGLL3_PADDED, &c22store[i*NGLL3],
                                         NGLL3*sizeof(realw),cudaMemcpyHostToDevice),4806);
      print_CUDA_error_if_any(cudaMemcpy(mp->d_c23store + i*NGLL3_PADDED, &c23store[i*NGLL3],
                                         NGLL3*sizeof(realw),cudaMemcpyHostToDevice),4807);
      print_CUDA_error_if_any(cudaMemcpy(mp->d_c24store + i*NGLL3_PADDED, &c24store[i*NGLL3],
                                         NGLL3*sizeof(realw),cudaMemcpyHostToDevice),4808);
      print_CUDA_error_if_any(cudaMemcpy(mp->d_c25store + i*NGLL3_PADDED, &c25store[i*NGLL3],
                                         NGLL3*sizeof(realw),cudaMemcpyHostToDevice),4809);
      print_CUDA_error_if_any(cudaMemcpy(mp->d_c26store + i*NGLL3_PADDED, &c26store[i*NGLL3],
                                         NGLL3*sizeof(realw),cudaMemcpyHostToDevice),4810);
      print_CUDA_error_if_any(cudaMemcpy(mp->d_c33store + i*NGLL3_PADDED, &c33store[i*NGLL3],
                                         NGLL3*sizeof(realw),cudaMemcpyHostToDevice),4811);
      print_CUDA_error_if_any(cudaMemcpy(mp->d_c34store + i*NGLL3_PADDED, &c34store[i*NGLL3],
                                         NGLL3*sizeof(realw),cudaMemcpyHostToDevice),4812);
      print_CUDA_error_if_any(cudaMemcpy(mp->d_c35store + i*NGLL3_PADDED, &c35store[i*NGLL3],
                                         NGLL3*sizeof(realw),cudaMemcpyHostToDevice),4813);
      print_CUDA_error_if_any(cudaMemcpy(mp->d_c36store + i*NGLL3_PADDED, &c36store[i*NGLL3],
                                         NGLL3*sizeof(realw),cudaMemcpyHostToDevice),4814);
      print_CUDA_error_if_any(cudaMemcpy(mp->d_c44store + i*NGLL3_PADDED, &c44store[i*NGLL3],
                                         NGLL3*sizeof(realw),cudaMemcpyHostToDevice),4815);
      print_CUDA_error_if_any(cudaMemcpy(mp->d_c45store + i*NGLL3_PADDED, &c45store[i*NGLL3],
                                         NGLL3*sizeof(realw),cudaMemcpyHostToDevice),4816);
      print_CUDA_error_if_any(cudaMemcpy(mp->d_c46store + i*NGLL3_PADDED, &c46store[i*NGLL3],
                                         NGLL3*sizeof(realw),cudaMemcpyHostToDevice),4817);
      print_CUDA_error_if_any(cudaMemcpy(mp->d_c55store + i*NGLL3_PADDED, &c55store[i*NGLL3],
                                         NGLL3*sizeof(realw),cudaMemcpyHostToDevice),4818);
      print_CUDA_error_if_any(cudaMemcpy(mp->d_c56store + i*NGLL3_PADDED, &c56store[i*NGLL3],
                                         NGLL3*sizeof(realw),cudaMemcpyHostToDevice),4819);
      print_CUDA_error_if_any(cudaMemcpy(mp->d_c66store + i*NGLL3_PADDED, &c66store[i*NGLL3],
                                         NGLL3*sizeof(realw),cudaMemcpyHostToDevice),4820);
    }
  }

  // ocean load approximation
  if( *APPROXIMATE_OCEAN_LOAD ){
    // oceans needs a free surface
    mp->num_free_surface_faces = *num_free_surface_faces;
    if( mp->num_free_surface_faces > 0 ){
      // mass matrix
      copy_todevice_realw((void**)&mp->d_rmass_ocean_load,rmass_ocean_load,mp->NGLOB_AB);

//      print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_rmass_ocean_load),
//                                         sizeof(realw)*mp->NGLOB_AB),4501);
//      print_CUDA_error_if_any(cudaMemcpy(mp->d_rmass_ocean_load,rmass_ocean_load,
//                                         sizeof(realw)*mp->NGLOB_AB,cudaMemcpyHostToDevice),4502);
      // surface normal
      copy_todevice_realw((void**)&mp->d_free_surface_normal,free_surface_normal,
                          3*NGLL2*(mp->num_free_surface_faces));

//      print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_free_surface_normal),
//                                         3*NGLL2*(mp->num_free_surface_faces)*sizeof(realw)),4503);
//      print_CUDA_error_if_any(cudaMemcpy(mp->d_free_surface_normal,free_surface_normal,
//                                         3*NGLL2*(mp->num_free_surface_faces)*sizeof(realw),cudaMemcpyHostToDevice),4504);

      // temporary global array: used to synchronize updates on global accel array
      print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_updated_dof_ocean_load),
                                         sizeof(int)*mp->NGLOB_AB),4505);

      if( *NOISE_TOMOGRAPHY == 0 && *ACOUSTIC_SIMULATION == 0 ){

        copy_todevice_int((void**)&mp->d_free_surface_ispec,free_surface_ispec,mp->num_free_surface_faces);

//        print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_free_surface_ispec),
//                                          mp->num_free_surface_faces*sizeof(int)),4601);
//        print_CUDA_error_if_any(cudaMemcpy(mp->d_free_surface_ispec,free_surface_ispec,
//                                          mp->num_free_surface_faces*sizeof(int),cudaMemcpyHostToDevice),4603);

        copy_todevice_int((void**)&mp->d_free_surface_ijk,free_surface_ijk,
                          3*NGLL2*mp->num_free_surface_faces);

//        print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_free_surface_ijk),
//                                          3*NGLL2*mp->num_free_surface_faces*sizeof(int)),4602);
//        print_CUDA_error_if_any(cudaMemcpy(mp->d_free_surface_ijk,free_surface_ijk,
//                                          3*NGLL2*mp->num_free_surface_faces*sizeof(int),cudaMemcpyHostToDevice),4604);
      }
    }
  }

  // mesh coloring
  if( mp->use_mesh_coloring_gpu ){
    mp->num_colors_outer_elastic = *num_colors_outer_elastic;
    mp->num_colors_inner_elastic = *num_colors_inner_elastic;
    mp->h_num_elem_colors_elastic = (int*) num_elem_colors_elastic;
  }

  // JC JC here we will need to add GPU support for the new C-PML routines

#ifdef ENABLE_VERY_SLOW_ERROR_CHECKING
  exit_on_cuda_error("prepare_fields_elastic_device");
#endif
}

/* ----------------------------------------------------------------------------------------------- */

extern "C"
void FC_FUNC_(prepare_fields_elastic_adj_dev,
              PREPARE_FIELDS_ELASTIC_ADJ_DEV)(long* Mesh_pointer_f,
                                             int* size,
                                             int* COMPUTE_AND_STORE_STRAIN,
                                             realw* epsilon_trace_over_3,
                                             realw* b_epsilondev_xx,realw* b_epsilondev_yy,realw* b_epsilondev_xy,
                                             realw* b_epsilondev_xz,realw* b_epsilondev_yz,
                                             realw* b_epsilon_trace_over_3,
                                             int* ATTENUATION,
                                             int* R_size,
                                             realw* b_R_xx,realw* b_R_yy,realw* b_R_xy,realw* b_R_xz,realw* b_R_yz,
                                             realw* b_alphaval,realw* b_betaval,realw* b_gammaval,
                                             int* APPROXIMATE_HESS_KL){

  TRACE("prepare_fields_elastic_adj_dev");

  Mesh* mp = (Mesh*)(*Mesh_pointer_f);

  // checks if kernel simulation
  if( mp->simulation_type != 3 ) return;

  // kernel simulations
  // allocates backward/reconstructed arrays on device (GPU)
  print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_b_displ),sizeof(realw)*(*size)),5201);
  print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_b_veloc),sizeof(realw)*(*size)),5202);
  print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_b_accel),sizeof(realw)*(*size)),5203);

  // allocates kernels
  print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_rho_kl),NGLL3*mp->NSPEC_AB*sizeof(realw)),5204);
  print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_mu_kl),NGLL3*mp->NSPEC_AB*sizeof(realw)),5205);
  print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_kappa_kl),NGLL3*mp->NSPEC_AB*sizeof(realw)),5206);

  // initializes kernel values to zero
  print_CUDA_error_if_any(cudaMemset(mp->d_rho_kl,0,
                                     NGLL3*mp->NSPEC_AB*sizeof(realw)),5207);
  print_CUDA_error_if_any(cudaMemset(mp->d_mu_kl,0,
                                     NGLL3*mp->NSPEC_AB*sizeof(realw)),5208);
  print_CUDA_error_if_any(cudaMemset(mp->d_kappa_kl,0,
                                     NGLL3*mp->NSPEC_AB*sizeof(realw)),5209);

  // strains used for attenuation and kernel simulations
  if( *COMPUTE_AND_STORE_STRAIN ){
    // strains
    int epsilondev_size = NGLL3*mp->NSPEC_AB; // note: non-aligned; if align, check memcpy below and indexing

    // solid pressure
    copy_todevice_realw((void**)&mp->d_epsilon_trace_over_3,epsilon_trace_over_3,NGLL3*mp->NSPEC_AB);

//    print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_epsilon_trace_over_3),
//                                       NGLL3*mp->NSPEC_AB*sizeof(realw)),5310);
//    print_CUDA_error_if_any(cudaMemcpy(mp->d_epsilon_trace_over_3,epsilon_trace_over_3,
//                                       NGLL3*mp->NSPEC_AB*sizeof(realw),cudaMemcpyHostToDevice),5311);
    // backward solid pressure

    copy_todevice_realw((void**)&mp->d_b_epsilon_trace_over_3,b_epsilon_trace_over_3,NGLL3*mp->NSPEC_AB);

//    print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_b_epsilon_trace_over_3),
//                                       NGLL3*mp->NSPEC_AB*sizeof(realw)),5312);
//    print_CUDA_error_if_any(cudaMemcpy(mp->d_b_epsilon_trace_over_3 ,b_epsilon_trace_over_3,
//                                       NGLL3*mp->NSPEC_AB*sizeof(realw),cudaMemcpyHostToDevice),5313);

    // prepares backward strains

    copy_todevice_realw((void**)&mp->d_b_epsilondev_xx,b_epsilondev_xx,epsilondev_size);
    copy_todevice_realw((void**)&mp->d_b_epsilondev_yy,b_epsilondev_yy,epsilondev_size);
    copy_todevice_realw((void**)&mp->d_b_epsilondev_xy,b_epsilondev_xy,epsilondev_size);
    copy_todevice_realw((void**)&mp->d_b_epsilondev_xz,b_epsilondev_xz,epsilondev_size);
    copy_todevice_realw((void**)&mp->d_b_epsilondev_yz,b_epsilondev_yz,epsilondev_size);

//    print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_b_epsilondev_xx),
//                                       epsilondev_size*sizeof(realw)),5321);
//    print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_b_epsilondev_yy),
//                                       epsilondev_size*sizeof(realw)),5322);
//    print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_b_epsilondev_xy),
//                                       epsilondev_size*sizeof(realw)),5323);
//    print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_b_epsilondev_xz),
//                                       epsilondev_size*sizeof(realw)),5324);
//    print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_b_epsilondev_yz),
//                                       epsilondev_size*sizeof(realw)),5325);

//    print_CUDA_error_if_any(cudaMemcpy(mp->d_b_epsilondev_xx,b_epsilondev_xx,
//                                       epsilondev_size*sizeof(realw),cudaMemcpyHostToDevice),5326);
//    print_CUDA_error_if_any(cudaMemcpy(mp->d_b_epsilondev_yy,b_epsilondev_yy,
//                                       epsilondev_size*sizeof(realw),cudaMemcpyHostToDevice),5327);
//    print_CUDA_error_if_any(cudaMemcpy(mp->d_b_epsilondev_xy,b_epsilondev_xy,
//                                       epsilondev_size*sizeof(realw),cudaMemcpyHostToDevice),5328);
//    print_CUDA_error_if_any(cudaMemcpy(mp->d_b_epsilondev_xz,b_epsilondev_xz,
//                                       epsilondev_size*sizeof(realw),cudaMemcpyHostToDevice),5329);
//    print_CUDA_error_if_any(cudaMemcpy(mp->d_b_epsilondev_yz,b_epsilondev_yz,
//                                       epsilondev_size*sizeof(realw),cudaMemcpyHostToDevice),5330);
  }

  // attenuation memory variables
  if( *ATTENUATION ){
    copy_todevice_realw((void**)&mp->d_b_R_xx,b_R_xx,(*R_size));
    copy_todevice_realw((void**)&mp->d_b_R_yy,b_R_yy,(*R_size));
    copy_todevice_realw((void**)&mp->d_b_R_xy,b_R_xy,(*R_size));
    copy_todevice_realw((void**)&mp->d_b_R_xz,b_R_xz,(*R_size));
    copy_todevice_realw((void**)&mp->d_b_R_yz,b_R_yz,(*R_size));

//    print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_b_R_xx),
//                                       (*R_size)*sizeof(realw)),5421);
//    print_CUDA_error_if_any(cudaMemcpy(mp->d_b_R_xx,b_R_xx,(*R_size)*sizeof(realw),
//                                       cudaMemcpyHostToDevice),5422);

//    print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_b_R_yy),
//                                       (*R_size)*sizeof(realw)),5423);
//    print_CUDA_error_if_any(cudaMemcpy(mp->d_b_R_yy,b_R_yy,(*R_size)*sizeof(realw),
//                                       cudaMemcpyHostToDevice),5424);

//    print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_b_R_xy),
//                                       (*R_size)*sizeof(realw)),5425);
//    print_CUDA_error_if_any(cudaMemcpy(mp->d_b_R_xy,b_R_xy,(*R_size)*sizeof(realw),
//                                       cudaMemcpyHostToDevice),5426);

//    print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_b_R_xz),
//                                       (*R_size)*sizeof(realw)),5427);
//    print_CUDA_error_if_any(cudaMemcpy(mp->d_b_R_xz,b_R_xz,(*R_size)*sizeof(realw),
//                                       cudaMemcpyHostToDevice),5428);

//    print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_b_R_yz),
//                                       (*R_size)*sizeof(realw)),5429);
//   print_CUDA_error_if_any(cudaMemcpy(mp->d_b_R_yz,b_R_yz,(*R_size)*sizeof(realw),
//                                       cudaMemcpyHostToDevice),5420);

    // alpha,beta,gamma factors for backward fields
    copy_todevice_realw((void**)&mp->d_b_alphaval,b_alphaval,N_SLS);
    copy_todevice_realw((void**)&mp->d_b_betaval,b_betaval,N_SLS);
    copy_todevice_realw((void**)&mp->d_b_gammaval,b_gammaval,N_SLS);

//    print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_b_alphaval),
//                                       N_SLS*sizeof(realw)),5434);
//    print_CUDA_error_if_any(cudaMemcpy(mp->d_b_alphaval ,b_alphaval,
//                                       N_SLS*sizeof(realw),cudaMemcpyHostToDevice),5435);

//    print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_b_betaval),
//                                       N_SLS*sizeof(realw)),5436);
//    print_CUDA_error_if_any(cudaMemcpy(mp->d_b_betaval ,b_betaval,
//                                       N_SLS*sizeof(realw),cudaMemcpyHostToDevice),5437);

//    print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_b_gammaval),
//                                       N_SLS*sizeof(realw)),5438);
//    print_CUDA_error_if_any(cudaMemcpy(mp->d_b_gammaval ,b_gammaval,
//                                       N_SLS*sizeof(realw),cudaMemcpyHostToDevice),5439);
  }

  if( *APPROXIMATE_HESS_KL ){
    print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_hess_el_kl),NGLL3*mp->NSPEC_AB*sizeof(realw)),5450);
    // initializes with zeros
    print_CUDA_error_if_any(cudaMemset(mp->d_hess_el_kl,0,
                                       NGLL3*mp->NSPEC_AB*sizeof(realw)),5451);
  }

#ifdef ENABLE_VERY_SLOW_ERROR_CHECKING
  exit_on_cuda_error("prepare_fields_elastic_adj_dev");
#endif
}

/* ----------------------------------------------------------------------------------------------- */

// purely adjoint & kernel simulations

/* ----------------------------------------------------------------------------------------------- */

extern "C"
void FC_FUNC_(prepare_sim2_or_3_const_device,
              PREPARE_SIM2_OR_3_CONST_DEVICE)(
                                              long* Mesh_pointer_f,
                                              int* islice_selected_rec,
                                              int* islice_selected_rec_size,
                                              int* nadj_rec_local,
                                              int* nrec,
                                              int* myrank) {

  TRACE("prepare_sim2_or_3_const_device");

  Mesh* mp = (Mesh*)(*Mesh_pointer_f);

  // adjoint source arrays
  mp->nadj_rec_local = *nadj_rec_local;
  if( mp->nadj_rec_local > 0 ){
    print_CUDA_error_if_any(cudaMalloc((void**)&mp->d_adj_sourcearrays,
                                       (mp->nadj_rec_local)*3*NGLL3*sizeof(realw)),6003);

    print_CUDA_error_if_any(cudaMalloc((void**)&mp->d_pre_computed_irec,
                                       (mp->nadj_rec_local)*sizeof(int)),6004);

    // prepares local irec array:
    // the irec_local variable needs to be precomputed (as
    // h_pre_comp..), because normally it is in the loop updating accel,
    // and due to how it's incremented, it cannot be parallelized
    int* h_pre_computed_irec = (int*) malloc( (mp->nadj_rec_local)*sizeof(int) );
    if( h_pre_computed_irec == NULL ) exit_on_error("prepare_sim2_or_3_const_device: h_pre_computed_irec not allocated\n");

    int irec_local = 0;
    for(int irec = 0; irec < *nrec; irec++) {
      if(*myrank == islice_selected_rec[irec]) {
        irec_local++;
        h_pre_computed_irec[irec_local-1] = irec;
      }
    }
    if( irec_local != mp->nadj_rec_local ) exit_on_error("prepare_sim2_or_3_const_device: irec_local not equal\n");
    // copies values onto GPU
    print_CUDA_error_if_any(cudaMemcpy(mp->d_pre_computed_irec,h_pre_computed_irec,
                                       (mp->nadj_rec_local)*sizeof(int),cudaMemcpyHostToDevice),6010);
    free(h_pre_computed_irec);

    // temporary array to prepare extracted source array values
    mp->h_adj_sourcearrays_slice = (realw*) malloc( (mp->nadj_rec_local)*3*NGLL3*sizeof(realw) );
    if( mp->h_adj_sourcearrays_slice == NULL ) exit_on_error("h_adj_sourcearrays_slice not allocated\n");

  }

#ifdef ENABLE_VERY_SLOW_ERROR_CHECKING
  exit_on_cuda_error("prepare_sim2_or_3_const_device");
#endif
}


/* ----------------------------------------------------------------------------------------------- */

// for NOISE simulations

/* ----------------------------------------------------------------------------------------------- */

extern "C"
void FC_FUNC_(prepare_fields_noise_device,
              PREPARE_FIELDS_NOISE_DEVICE)(long* Mesh_pointer_f,
                                           int* NSPEC_AB, int* NGLOB_AB,
                                           int* free_surface_ispec,
                                           int* free_surface_ijk,
                                           int* num_free_surface_faces,
                                           int* NOISE_TOMOGRAPHY,
                                           int* NSTEP,
                                           realw* noise_sourcearray,
                                           realw* normal_x_noise,
                                           realw* normal_y_noise,
                                           realw* normal_z_noise,
                                           realw* mask_noise,
                                           realw* free_surface_jacobian2Dw) {

  TRACE("prepare_fields_noise_device");

  Mesh* mp = (Mesh*)(*Mesh_pointer_f);

  // free surface
  mp->num_free_surface_faces = *num_free_surface_faces;

  copy_todevice_int((void**)&mp->d_free_surface_ispec,free_surface_ispec,mp->num_free_surface_faces);

//  print_CUDA_error_if_any(cudaMalloc((void**) &mp->d_free_surface_ispec,
//                                     mp->num_free_surface_faces*sizeof(int)),7001);
//  print_CUDA_error_if_any(cudaMemcpy(mp->d_free_surface_ispec, free_surface_ispec,
//                                     mp->num_free_surface_faces*sizeof(int),cudaMemcpyHostToDevice),7002);

  copy_todevice_int((void**)&mp->d_free_surface_ijk,free_surface_ijk,
                    3*NGLL2*mp->num_free_surface_faces);

//  print_CUDA_error_if_any(cudaMalloc((void**) &mp->d_free_surface_ijk,
//                                     3*NGLL2*mp->num_free_surface_faces*sizeof(int)),7003);
//  print_CUDA_error_if_any(cudaMemcpy(mp->d_free_surface_ijk,free_surface_ijk,
//                                     3*NGLL2*mp->num_free_surface_faces*sizeof(int),cudaMemcpyHostToDevice),7004);

  // alloc storage for the surface buffer to be copied
  print_CUDA_error_if_any(cudaMalloc((void**) &mp->d_noise_surface_movie,
                                     3*NGLL2*mp->num_free_surface_faces*sizeof(realw)),7005);

  // prepares noise source array
  if( *NOISE_TOMOGRAPHY == 1 ){
    copy_todevice_realw((void**)&mp->d_noise_sourcearray,noise_sourcearray,
                        3*NGLL3*(*NSTEP));

//    print_CUDA_error_if_any(cudaMalloc((void**)&mp->d_noise_sourcearray,
//                                       3*NGLL3*(*NSTEP)*sizeof(realw)),7101);
//    print_CUDA_error_if_any(cudaMemcpy(mp->d_noise_sourcearray, noise_sourcearray,
//                                       3*NGLL3*(*NSTEP)*sizeof(realw),cudaMemcpyHostToDevice),7102);
  }

  // prepares noise directions
  if( *NOISE_TOMOGRAPHY > 1 ){
    int nface_size = NGLL2*(*num_free_surface_faces);
    // allocates memory on GPU
    copy_todevice_realw((void**)&mp->d_normal_x_noise,normal_x_noise,nface_size);
    copy_todevice_realw((void**)&mp->d_normal_y_noise,normal_y_noise,nface_size);
    copy_todevice_realw((void**)&mp->d_normal_z_noise,normal_z_noise,nface_size);

    copy_todevice_realw((void**)&mp->d_mask_noise,mask_noise,nface_size);
    copy_todevice_realw((void**)&mp->d_free_surface_jacobian2Dw,free_surface_jacobian2Dw,nface_size);

//    print_CUDA_error_if_any(cudaMalloc((void**)&mp->d_normal_x_noise,
//                                       nface_size*sizeof(realw)),7301);
//    print_CUDA_error_if_any(cudaMalloc((void**)&mp->d_normal_y_noise,
//                                       nface_size*sizeof(realw)),7302);
//    print_CUDA_error_if_any(cudaMalloc((void**)&mp->d_normal_z_noise,
//                                       nface_size*sizeof(realw)),7303);

//    print_CUDA_error_if_any(cudaMalloc((void**)&mp->d_mask_noise,
//                                       nface_size*sizeof(realw)),7304);
//    print_CUDA_error_if_any(cudaMalloc((void**)&mp->d_free_surface_jacobian2Dw,
//                                       nface_size*sizeof(realw)),7305);
    // transfers data onto GPU
//    print_CUDA_error_if_any(cudaMemcpy(mp->d_normal_x_noise, normal_x_noise,
//                                       nface_size*sizeof(realw),cudaMemcpyHostToDevice),7306);
//    print_CUDA_error_if_any(cudaMemcpy(mp->d_normal_y_noise, normal_y_noise,
//                                       nface_size*sizeof(realw),cudaMemcpyHostToDevice),7307);
//    print_CUDA_error_if_any(cudaMemcpy(mp->d_normal_z_noise, normal_z_noise,
//                                       nface_size*sizeof(realw),cudaMemcpyHostToDevice),7308);
//    print_CUDA_error_if_any(cudaMemcpy(mp->d_mask_noise, mask_noise,
//                                       nface_size*sizeof(realw),cudaMemcpyHostToDevice),7309);
//    print_CUDA_error_if_any(cudaMemcpy(mp->d_free_surface_jacobian2Dw, free_surface_jacobian2Dw,
//                                       nface_size*sizeof(realw),cudaMemcpyHostToDevice),7310);
  }

  // prepares noise strength kernel
  if( *NOISE_TOMOGRAPHY == 3 ){
    print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_Sigma_kl),
                                       NGLL3*(mp->NSPEC_AB)*sizeof(realw)),7401);
    // initializes kernel values to zero
    print_CUDA_error_if_any(cudaMemset(mp->d_Sigma_kl,0,
                                       NGLL3*mp->NSPEC_AB*sizeof(realw)),7403);

  }

#ifdef ENABLE_VERY_SLOW_ERROR_CHECKING
  //printf("jacobian_size = %d\n",25*(*num_free_surface_faces));
  exit_on_cuda_error("prepare_fields_noise_device");
#endif
}

/* ----------------------------------------------------------------------------------------------- */

// GRAVITY simulations

/* ----------------------------------------------------------------------------------------------- */

extern "C"
void FC_FUNC_(prepare_fields_gravity_device,
              PREPARE_FIELDS_gravity_DEVICE)(long* Mesh_pointer_f,
                                             int* GRAVITY,
                                             realw* minus_deriv_gravity,
                                             realw* minus_g,
                                             realw* h_wgll_cube,
                                             int* ACOUSTIC_SIMULATION,
                                             realw* rhostore) {

  TRACE("prepare_fields_gravity_device");

  Mesh* mp = (Mesh*)(*Mesh_pointer_f);

  setConst_wgll_cube(h_wgll_cube,mp);

  mp->gravity = *GRAVITY;
  if( mp->gravity ){

    copy_todevice_realw((void**)&mp->d_minus_deriv_gravity,minus_deriv_gravity,mp->NGLOB_AB);

//    print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_minus_deriv_gravity),
//                                       (mp->NGLOB_AB)*sizeof(realw)),8000);
//    print_CUDA_error_if_any(cudaMemcpy(mp->d_minus_deriv_gravity, minus_deriv_gravity,
//                                       (mp->NGLOB_AB)*sizeof(realw),cudaMemcpyHostToDevice),8001);

    copy_todevice_realw((void**)&mp->d_minus_g,minus_g,mp->NGLOB_AB);

//    print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_minus_g),
//                                       (mp->NGLOB_AB)*sizeof(realw)),8002);
//    print_CUDA_error_if_any(cudaMemcpy(mp->d_minus_g, minus_g,
//                                       (mp->NGLOB_AB)*sizeof(realw),cudaMemcpyHostToDevice),8003);


    if( *ACOUSTIC_SIMULATION == 0 ){
      // rhostore not allocated yet
      int size_padded = NGLL3_PADDED * (mp->NSPEC_AB);
      // padded array
      print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_rhostore),size_padded*sizeof(realw)),8006);
      // transfer constant element data with padding
      for(int i=0; i < mp->NSPEC_AB; i++) {
        print_CUDA_error_if_any(cudaMemcpy(mp->d_rhostore+i*NGLL3_PADDED, &rhostore[i*NGLL3],
                                           NGLL3*sizeof(realw),cudaMemcpyHostToDevice),8007);
      }
    }
  }

}

/* ----------------------------------------------------------------------------------------------- */

extern "C"
void FC_FUNC_(prepare_seismogram_fields,
              PREPARE_SEISMOGRAM_FIELDS)(long* Mesh_pointer,int* nrec_local, double* nu, double* hxir, double* hetar, double* hgammar) {

  TRACE("prepare_constants_device");
  Mesh* mp = (Mesh*)(*Mesh_pointer);

  print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_nu),3*3*(*nrec_local)*sizeof(double)),8100);
  print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_hxir),5*(*nrec_local)*sizeof(double)),8100);
  print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_hetar),5*(*nrec_local)*sizeof(double)),8100);
  print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_hgammar),5*(*nrec_local)*sizeof(double)),8100);

  print_CUDA_error_if_any(cudaMalloc((void**)&mp->d_seismograms_d,3*(*nrec_local)*sizeof(realw)),8101);
  print_CUDA_error_if_any(cudaMalloc((void**)&mp->d_seismograms_v,3*(*nrec_local)*sizeof(realw)),8101);
  print_CUDA_error_if_any(cudaMalloc((void**)&mp->d_seismograms_a,3*(*nrec_local)*sizeof(realw)),8101);

  print_CUDA_error_if_any(cudaMemcpy(mp->d_nu,nu,3*3*(*nrec_local)*sizeof(double),cudaMemcpyHostToDevice),8101);
  print_CUDA_error_if_any(cudaMemcpy(mp->d_hxir,hxir,5*(*nrec_local)*sizeof(double),cudaMemcpyHostToDevice),8101);
  print_CUDA_error_if_any(cudaMemcpy(mp->d_hetar,hetar,5*(*nrec_local)*sizeof(double),cudaMemcpyHostToDevice),8101);
  print_CUDA_error_if_any(cudaMemcpy(mp->d_hgammar,hgammar,5*(*nrec_local)*sizeof(double),cudaMemcpyHostToDevice),8101);

  cudaMallocHost((void**)&mp->h_seismograms_d_it,3**nrec_local*sizeof(realw));
  cudaMallocHost((void**)&mp->h_seismograms_v_it,3**nrec_local*sizeof(realw));
  cudaMallocHost((void**)&mp->h_seismograms_a_it,3**nrec_local*sizeof(realw));
}

/* ----------------------------------------------------------------------------------------------- */

// cleanup

/* ----------------------------------------------------------------------------------------------- */

extern "C"
void FC_FUNC_(prepare_cleanup_device,
              PREPARE_CLEANUP_DEVICE)(long* Mesh_pointer_f,
                                      int* SAVE_FORWARD,
                                      int* ACOUSTIC_SIMULATION,
                                      int* ELASTIC_SIMULATION,
                                      int* ABSORBING_CONDITIONS,
                                      int* NOISE_TOMOGRAPHY,
                                      int* COMPUTE_AND_STORE_STRAIN,
                                      int* ATTENUATION,
                                      int* ANISOTROPY,
                                      int* APPROXIMATE_OCEAN_LOAD,
                                      int* APPROXIMATE_HESS_KL) {

TRACE("prepare_cleanup_device");

  // frees allocated memory arrays
  Mesh* mp = (Mesh*)(*Mesh_pointer_f);

  // frees memory on GPU
  // mesh
  cudaFree(mp->d_xix);
  cudaFree(mp->d_xiy);
  cudaFree(mp->d_xiz);
  cudaFree(mp->d_etax);
  cudaFree(mp->d_etay);
  cudaFree(mp->d_etaz);
  cudaFree(mp->d_gammax);
  cudaFree(mp->d_gammay);
  cudaFree(mp->d_gammaz);
  cudaFree(mp->d_muv);

  // absorbing boundaries
  if( *ABSORBING_CONDITIONS && mp->d_num_abs_boundary_faces > 0 ){
    cudaFree(mp->d_abs_boundary_ispec);
    cudaFree(mp->d_abs_boundary_ijk);
    cudaFree(mp->d_abs_boundary_normal);
    cudaFree(mp->d_abs_boundary_jacobian2Dw);
  }

  // interfaces
  cudaFree(mp->d_nibool_interfaces_ext_mesh);
  cudaFree(mp->d_ibool_interfaces_ext_mesh);

  // global indexing
  cudaFree(mp->d_ispec_is_inner);
  cudaFree(mp->d_ibool);

  // sources
  if (mp->simulation_type == 1  || mp->simulation_type == 3){
    cudaFree(mp->d_sourcearrays);
    cudaFree(mp->d_stf_pre_compute);
  }

  cudaFree(mp->d_islice_selected_source);
  cudaFree(mp->d_ispec_selected_source);

  // receivers
  if( mp->nrec_local > 0 ) cudaFree(mp->d_number_receiver_global);
  cudaFree(mp->d_ispec_selected_rec);

  // ACOUSTIC arrays
  if( *ACOUSTIC_SIMULATION ){
    cudaFree(mp->d_potential_acoustic);
    cudaFree(mp->d_potential_dot_acoustic);
    cudaFree(mp->d_potential_dot_dot_acoustic);
    cudaFree(mp->d_send_potential_dot_dot_buffer);
    cudaFree(mp->d_rmass_acoustic);
    cudaFree(mp->d_rhostore);
    cudaFree(mp->d_kappastore);
    cudaFree(mp->d_phase_ispec_inner_acoustic);
    cudaFree(mp->d_ispec_is_acoustic);

    if( *NOISE_TOMOGRAPHY == 0 ){
      cudaFree(mp->d_free_surface_ispec);
      cudaFree(mp->d_free_surface_ijk);
    }

    if( *ABSORBING_CONDITIONS ) cudaFree(mp->d_b_absorb_potential);

    if( mp->simulation_type == 3 ) {
      cudaFree(mp->d_b_potential_acoustic);
      cudaFree(mp->d_b_potential_dot_acoustic);
      cudaFree(mp->d_b_potential_dot_dot_acoustic);
      cudaFree(mp->d_rho_ac_kl);
      cudaFree(mp->d_kappa_ac_kl);
      if( *APPROXIMATE_HESS_KL) cudaFree(mp->d_hess_ac_kl);
    }


    if(mp->nrec_local > 0 ){
      cudaFree(mp->d_station_seismo_potential);
      free(mp->h_station_seismo_potential);
    }

  } // ACOUSTIC_SIMULATION

  // ELASTIC arrays
  if( *ELASTIC_SIMULATION ){
    cudaFree(mp->d_displ);
    cudaFree(mp->d_veloc);
    cudaFree(mp->d_accel);
    cudaFree(mp->d_send_accel_buffer);

    cudaFree(mp->d_rmassx);
    cudaFree(mp->d_rmassy);
    cudaFree(mp->d_rmassz);

    cudaFree(mp->d_phase_ispec_inner_elastic);
    cudaFree(mp->d_ispec_is_elastic);

    if( mp->nrec_local > 0 ){
      cudaFree(mp->d_station_seismo_field);
      free(mp->h_station_seismo_field);
    }

    if( *ABSORBING_CONDITIONS && mp->d_num_abs_boundary_faces > 0){
      cudaFree(mp->d_rho_vp);
      cudaFree(mp->d_rho_vs);

      if(mp->simulation_type == 3 || ( mp->simulation_type == 1 && *SAVE_FORWARD ))
          cudaFree(mp->d_b_absorb_field);
    }

    if( mp->simulation_type == 3 ) {
      cudaFree(mp->d_b_displ);
      cudaFree(mp->d_b_veloc);
      cudaFree(mp->d_b_accel);
      cudaFree(mp->d_rho_kl);
      cudaFree(mp->d_mu_kl);
      cudaFree(mp->d_kappa_kl);
      if( *APPROXIMATE_HESS_KL ) cudaFree(mp->d_hess_el_kl);
    }

    if( *COMPUTE_AND_STORE_STRAIN ){
      cudaFree(mp->d_epsilondev_xx);
      cudaFree(mp->d_epsilondev_yy);
      cudaFree(mp->d_epsilondev_xy);
      cudaFree(mp->d_epsilondev_xz);
      cudaFree(mp->d_epsilondev_yz);
      if( mp->simulation_type == 3 ){
        cudaFree(mp->d_epsilon_trace_over_3);
        cudaFree(mp->d_b_epsilon_trace_over_3);
        cudaFree(mp->d_b_epsilondev_xx);
        cudaFree(mp->d_b_epsilondev_yy);
        cudaFree(mp->d_b_epsilondev_xy);
        cudaFree(mp->d_b_epsilondev_xz);
        cudaFree(mp->d_b_epsilondev_yz);
      }
    }

    if( *ATTENUATION ){
      cudaFree(mp->d_factor_common);
      cudaFree(mp->d_one_minus_sum_beta);
      cudaFree(mp->d_alphaval);
      cudaFree(mp->d_betaval);
      cudaFree(mp->d_gammaval);
      cudaFree(mp->d_R_xx);
      cudaFree(mp->d_R_yy);
      cudaFree(mp->d_R_xy);
      cudaFree(mp->d_R_xz);
      cudaFree(mp->d_R_yz);
      if( mp->simulation_type == 3){
        cudaFree(mp->d_b_R_xx);
        cudaFree(mp->d_b_R_yy);
        cudaFree(mp->d_b_R_xy);
        cudaFree(mp->d_b_R_xz);
        cudaFree(mp->d_b_R_yz);
        cudaFree(mp->d_b_alphaval);
        cudaFree(mp->d_b_betaval);
        cudaFree(mp->d_b_gammaval);
      }
    }

    if( *ANISOTROPY ){
      cudaFree(mp->d_c11store);
      cudaFree(mp->d_c12store);
      cudaFree(mp->d_c13store);
      cudaFree(mp->d_c14store);
      cudaFree(mp->d_c15store);
      cudaFree(mp->d_c16store);
      cudaFree(mp->d_c22store);
      cudaFree(mp->d_c23store);
      cudaFree(mp->d_c24store);
      cudaFree(mp->d_c25store);
      cudaFree(mp->d_c26store);
      cudaFree(mp->d_c33store);
      cudaFree(mp->d_c34store);
      cudaFree(mp->d_c35store);
      cudaFree(mp->d_c36store);
      cudaFree(mp->d_c44store);
      cudaFree(mp->d_c45store);
      cudaFree(mp->d_c46store);
      cudaFree(mp->d_c55store);
      cudaFree(mp->d_c56store);
      cudaFree(mp->d_c66store);
    }

    if( *APPROXIMATE_OCEAN_LOAD ){
      if( mp->num_free_surface_faces > 0 ){
        cudaFree(mp->d_rmass_ocean_load);
        cudaFree(mp->d_free_surface_normal);
        cudaFree(mp->d_updated_dof_ocean_load);
        if( *NOISE_TOMOGRAPHY == 0){
          cudaFree(mp->d_free_surface_ispec);
          cudaFree(mp->d_free_surface_ijk);
        }
      }
    }
  } // ELASTIC_SIMULATION

  // purely adjoint & kernel array
  if( mp->simulation_type == 2 || mp->simulation_type == 3 ){
    if(mp->nadj_rec_local > 0 ){
      cudaFree(mp->d_adj_sourcearrays);
      cudaFree(mp->d_pre_computed_irec);
      free(mp->h_adj_sourcearrays_slice);
    }
  }

  // NOISE arrays
  if( *NOISE_TOMOGRAPHY > 0 ){
    cudaFree(mp->d_free_surface_ispec);
    cudaFree(mp->d_free_surface_ijk);
    cudaFree(mp->d_noise_surface_movie);
    if( *NOISE_TOMOGRAPHY == 1 ) cudaFree(mp->d_noise_sourcearray);
    if( *NOISE_TOMOGRAPHY > 1 ){
      cudaFree(mp->d_normal_x_noise);
      cudaFree(mp->d_normal_y_noise);
      cudaFree(mp->d_normal_z_noise);
      cudaFree(mp->d_mask_noise);
      cudaFree(mp->d_free_surface_jacobian2Dw);
    }
    if( *NOISE_TOMOGRAPHY == 3 ) cudaFree(mp->d_Sigma_kl);
  }

  // mesh pointer - not needed anymore
  free(mp);
}
