/*
 !=====================================================================
 !
 !                         S p e c f e m 3 D
 !                         -----------------
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
#include "prepare_constants_cuda.h"

#ifdef USE_CUDA
#ifdef USE_OLDER_CUDA4_GPU
#pragma message ("\nCompiling with: USE_OLDER_CUDA4_GPU enabled\n")
#endif
#endif

#ifdef USE_OLDER_CUDA4_GPU
#else
  #ifdef USE_TEXTURES_FIELDS
    // elastic
    extern realw_texture d_displ_tex;
    extern realw_texture d_veloc_tex;
    extern realw_texture d_accel_tex;
    // backward/reconstructed
    extern realw_texture d_b_displ_tex;
    extern realw_texture d_b_veloc_tex;
    extern realw_texture d_b_accel_tex;
    // acoustic
    extern realw_texture d_potential_tex;
    extern realw_texture d_potential_dot_dot_tex;
    // backward/reconstructed
    extern realw_texture d_b_potential_tex;
    extern realw_texture d_b_potential_dot_dot_tex;
  #endif
  #ifdef USE_TEXTURES_CONSTANTS
    extern realw_texture d_hprime_xx_tex;
  #endif
#endif


/* ----------------------------------------------------------------------------------------------- */

// GPU preparation

/* ----------------------------------------------------------------------------------------------- */

extern EXTERN_LANG
void FC_FUNC_(prepare_constants_device,
              PREPARE_CONSTANTS_DEVICE)(long* Mesh_pointer,
                                        int* h_NGLLX, int* NSPEC_AB, int* NGLOB_AB,
                                        int* NSPEC_IRREGULAR,int* h_irregular_element_number,
                                        realw* h_xix, realw* h_xiy, realw* h_xiz,
                                        realw* h_etax, realw* h_etay, realw* h_etaz,
                                        realw* h_gammax, realw* h_gammay, realw* h_gammaz,
                                        realw* xix_regular, realw* jacobian_regular,
                                        int* h_ibool,
                                        int* num_interfaces_ext_mesh, int* max_nibool_interfaces_ext_mesh,
                                        int* h_nibool_interfaces_ext_mesh, int* h_ibool_interfaces_ext_mesh,
                                        realw* h_hprime_xx, realw* h_hprimewgll_xx,
                                        realw* h_wgllwgll_xy,realw* h_wgllwgll_xz,realw* h_wgllwgll_yz,
                                        int* ABSORBING_CONDITIONS,
                                        int* h_abs_boundary_ispec, int* h_abs_boundary_ijk,
                                        realw* h_abs_boundary_normal,
                                        realw* h_abs_boundary_jacobian2Dw,
                                        int* h_num_abs_boundary_faces,
                                        int* h_ispec_is_inner,
                                        int* NSOURCES, int* nsources_local_f,
                                        realw* h_sourcearrays,
                                        int* h_islice_selected_source, int* h_ispec_selected_source,
                                        int* h_ispec_selected_rec,
                                        int* nrec, int* nrec_local,
                                        int* SIMULATION_TYPE,
                                        int* USE_MESH_COLORING_GPU_f,
                                        int* nspec_acoustic,int* nspec_elastic,
                                        int* ispec_is_acoustic,
                                        int* ispec_is_elastic,
                                        int* h_myrank,
                                        int* SAVE_FORWARD,
                                        realw* h_xir,realw* h_etar, realw* h_gammar,
                                        double* nu_rec, double* nu_source,
                                        int* h_islice_selected_rec,
                                        int* nlength_seismogram,
                                        int* SAVE_SEISMOGRAMS_DISPLACEMENT,int* SAVE_SEISMOGRAMS_VELOCITY,
                                        int* SAVE_SEISMOGRAMS_ACCELERATION,int* SAVE_SEISMOGRAMS_PRESSURE,
                                        int* h_NB_RUNS_ACOUSTIC_GPU,
                                        int* FAULT_SIMULATION,
                                        int* UNDO_ATTENUATION_AND_OR_PML) {

  TRACE("prepare_constants_device");

  // allocates mesh parameter structure
  Mesh* mp = (Mesh*) malloc( sizeof(Mesh) );
  if (mp == NULL) exit_on_error("error allocating mesh pointer");
  *Mesh_pointer = (long)mp;

  // sets processes mpi rank
  mp->myrank = *h_myrank;

  // sets global parameters
  mp->NSPEC_AB = *NSPEC_AB;
  mp->NSPEC_IRREGULAR = *NSPEC_IRREGULAR;
  mp->NGLOB_AB = *NGLOB_AB;

  // simulation flags
  mp->simulation_type = *SIMULATION_TYPE;
  mp->absorbing_conditions = *ABSORBING_CONDITIONS;  // STACEY_ABSORBING_CONDITIONS
  mp->save_forward = *SAVE_FORWARD;
  mp->undo_attenuation = *UNDO_ATTENUATION_AND_OR_PML;

  // checks setup
// DK DK August 2018: adding this test, following a suggestion by Etienne Bachmann
  if (*h_NGLLX != NGLLX) {
    exit_on_error("make sure that the NGLL constants are equal in the two files:\n" \
                  "  setup/constants.h and src/gpu/mesh_constants_gpu.h\n" \
                  "and then please re-compile; also make sure that the value of NGLL3_PADDED " \
                  "is consistent with the value of NGLL\n");
  }
  if (*h_NB_RUNS_ACOUSTIC_GPU != NB_RUNS_ACOUSTIC_GPU){
    exit_on_error("make sure that the NB_RUNS_ACOUSTIC_GPU constants are equal in the two files:\n" \
                  "  setup/constants.h and src/gpu/mesh_constants_gpu.h\n" \
                  "and then please re-compile...\n");
  }

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
#ifdef USE_CUDA
    #ifdef USE_OLDER_CUDA4_GPU
      cudaChannelFormatDesc channelDesc = cudaCreateChannelDesc<float>();
      const textureReference* d_hprime_xx_tex_ptr;
      print_CUDA_error_if_any(cudaGetTextureReference(&d_hprime_xx_tex_ptr, "d_hprime_xx_tex"), 4101);
      print_CUDA_error_if_any(cudaBindTexture(0, d_hprime_xx_tex_ptr, mp->d_hprime_xx, &channelDesc, sizeof(realw)*(NGLL2)), 4001);
   #else
      cudaChannelFormatDesc channelDesc = cudaCreateChannelDesc<float>();
      print_CUDA_error_if_any(cudaBindTexture(0, &d_hprime_xx_tex, mp->d_hprime_xx, &channelDesc, sizeof(realw)*(NGLL2)), 4001);
   #endif
#endif
#ifdef USE_HIP
    hipChannelFormatDesc channelDesc = hipCreateChannelDesc<float>();
    print_HIP_error_if_any(hipBindTexture(0, &d_hprime_xx_tex, mp->d_hprime_xx, &channelDesc, sizeof(realw)*(NGLL2)), 4001);
#endif
  }
  #endif

  gpuCreateCopy_todevice_int((void**)&mp->d_irregular_element_number,h_irregular_element_number,mp->NSPEC_AB);
  mp->xix_regular = *xix_regular;
  mp->jacobian_regular = *jacobian_regular;

  // mesh
  // Assuming NGLLX=5. Padded is then 128 (5^3+3)
  int size_padded = NGLL3_PADDED * (mp->NSPEC_IRREGULAR > 0 ? *NSPEC_IRREGULAR : 1);

  gpuMalloc_realw((void**) &mp->d_xix, size_padded);
  gpuMalloc_realw((void**) &mp->d_xiy, size_padded);
  gpuMalloc_realw((void**) &mp->d_xiz, size_padded);
  gpuMalloc_realw((void**) &mp->d_etax, size_padded);
  gpuMalloc_realw((void**) &mp->d_etay, size_padded);
  gpuMalloc_realw((void**) &mp->d_etaz, size_padded);
  gpuMalloc_realw((void**) &mp->d_gammax, size_padded);
  gpuMalloc_realw((void**) &mp->d_gammay, size_padded);
  gpuMalloc_realw((void**) &mp->d_gammaz, size_padded);

  // transfer constant element data with padding
  /*
  // way 1: slow...
  for(int i=0;i < mp->NSPEC_IRREGULAR;i++) {
    gpuMemcpy_todevice_realw(mp->d_xix + i*NGLL3_PADDED, &h_xix[i*NGLL3],NGLL3);
    gpuMemcpy_todevice_realw(mp->d_xiy+i*NGLL3_PADDED,   &h_xiy[i*NGLL3],NGLL3);
    gpuMemcpy_todevice_realw(mp->d_xiz+i*NGLL3_PADDED,   &h_xiz[i*NGLL3],NGLL3);
    gpuMemcpy_todevice_realw(mp->d_etax+i*NGLL3_PADDED,  &h_etax[i*NGLL3],NGLL3);
    gpuMemcpy_todevice_realw(mp->d_etay+i*NGLL3_PADDED,  &h_etay[i*NGLL3],NGLL3);
    gpuMemcpy_todevice_realw(mp->d_etaz+i*NGLL3_PADDED,  &h_etaz[i*NGLL3],NGLL3);
    gpuMemcpy_todevice_realw(mp->d_gammax+i*NGLL3_PADDED,&h_gammax[i*NGLL3],NGLL3);
    gpuMemcpy_todevice_realw(mp->d_gammay+i*NGLL3_PADDED,&h_gammay[i*NGLL3],NGLL3);
    gpuMemcpy_todevice_realw(mp->d_gammaz+i*NGLL3_PADDED,&h_gammaz[i*NGLL3],NGLL3);
    gpuMemcpy_todevice_realw(mp->d_kappav+i*NGLL3_PADDED,&h_kappav[i*NGLL3],NGLL3);
    gpuMemcpy_todevice_realw(mp->d_muv+i*NGLL3_PADDED,   &h_muv[i*NGLL3],NGLL3);
  }
  */
  // way 2: faster ....
  if (*NSPEC_IRREGULAR > 0 ){
    gpuMemcpy2D_todevice_realw(mp->d_xix, NGLL3_PADDED, h_xix, NGLL3, NGLL3, mp->NSPEC_IRREGULAR);
    gpuMemcpy2D_todevice_realw(mp->d_xiy, NGLL3_PADDED, h_xiy, NGLL3, NGLL3, mp->NSPEC_IRREGULAR);
    gpuMemcpy2D_todevice_realw(mp->d_xiz, NGLL3_PADDED, h_xiz, NGLL3, NGLL3, mp->NSPEC_IRREGULAR);
    gpuMemcpy2D_todevice_realw(mp->d_etax, NGLL3_PADDED, h_etax, NGLL3, NGLL3, mp->NSPEC_IRREGULAR);
    gpuMemcpy2D_todevice_realw(mp->d_etay, NGLL3_PADDED, h_etay, NGLL3, NGLL3, mp->NSPEC_IRREGULAR);
    gpuMemcpy2D_todevice_realw(mp->d_etaz, NGLL3_PADDED, h_etaz, NGLL3, NGLL3, mp->NSPEC_IRREGULAR);
    gpuMemcpy2D_todevice_realw(mp->d_gammax, NGLL3_PADDED, h_gammax, NGLL3, NGLL3, mp->NSPEC_IRREGULAR);
    gpuMemcpy2D_todevice_realw(mp->d_gammay, NGLL3_PADDED, h_gammay, NGLL3, NGLL3, mp->NSPEC_IRREGULAR);
    gpuMemcpy2D_todevice_realw(mp->d_gammaz, NGLL3_PADDED, h_gammaz, NGLL3, NGLL3, mp->NSPEC_IRREGULAR);
  }

  size_padded = NGLL3_PADDED * (mp->NSPEC_AB);

  // global indexing (padded)
  gpuMalloc_int((void**) &mp->d_ibool, size_padded);
  gpuMemcpy2D_todevice_int(mp->d_ibool, NGLL3_PADDED, h_ibool, NGLL3, NGLL3, mp->NSPEC_AB);


  // prepare interprocess-edge exchange information
  mp->num_interfaces_ext_mesh = *num_interfaces_ext_mesh;
  mp->max_nibool_interfaces_ext_mesh = *max_nibool_interfaces_ext_mesh;
  if (mp->num_interfaces_ext_mesh > 0){
    gpuCreateCopy_todevice_int((void**)&mp->d_nibool_interfaces_ext_mesh,h_nibool_interfaces_ext_mesh,mp->num_interfaces_ext_mesh);
    gpuCreateCopy_todevice_int((void**)&mp->d_ibool_interfaces_ext_mesh,h_ibool_interfaces_ext_mesh,
                         (mp->num_interfaces_ext_mesh)*(mp->max_nibool_interfaces_ext_mesh));
  }

  // setup two streams, one for compute and one for host<->device memory copies
  // compute stream
  gpuStreamCreate(&mp->compute_stream);
  // copy stream (needed to transfer mpi buffers)
  if (mp->num_interfaces_ext_mesh * mp->max_nibool_interfaces_ext_mesh > 0){
    gpuStreamCreate(&mp->copy_stream);
  }

  // inner elements
  gpuCreateCopy_todevice_int((void**)&mp->d_ispec_is_inner,h_ispec_is_inner,mp->NSPEC_AB);

  // absorbing boundaries
  mp->d_num_abs_boundary_faces = *h_num_abs_boundary_faces;
  if (mp->absorbing_conditions && mp->d_num_abs_boundary_faces > 0){
    gpuCreateCopy_todevice_int((void**)&mp->d_abs_boundary_ispec,h_abs_boundary_ispec,mp->d_num_abs_boundary_faces);
    gpuCreateCopy_todevice_int((void**)&mp->d_abs_boundary_ijk,h_abs_boundary_ijk,3*NGLL2*(mp->d_num_abs_boundary_faces));
    gpuCreateCopy_todevice_realw((void**)&mp->d_abs_boundary_normal,h_abs_boundary_normal,NDIM*NGLL2*(mp->d_num_abs_boundary_faces));
    gpuCreateCopy_todevice_realw((void**)&mp->d_abs_boundary_jacobian2Dw,h_abs_boundary_jacobian2Dw,NGLL2*(mp->d_num_abs_boundary_faces));
  }

  // sources
  mp->nsources_local = *nsources_local_f;
  if (mp->simulation_type == 1  || mp->simulation_type == 3){
    // not needed in case of pure adjoint simulations (SIMULATION_TYPE == 2)
    gpuCreateCopy_todevice_realw((void**)&mp->d_sourcearrays,h_sourcearrays,(*NSOURCES)*NDIM*NGLL3);

    // buffer for source time function values
    gpuMalloc_field((void**)&mp->d_stf_pre_compute,(*NSOURCES));
  }
  gpuCreateCopy_todevice_int((void**)&mp->d_islice_selected_source,h_islice_selected_source,(*NSOURCES));
  gpuCreateCopy_todevice_int((void**)&mp->d_ispec_selected_source,h_ispec_selected_source,(*NSOURCES));

  // seismogram outputs
  mp->save_seismograms_d = *SAVE_SEISMOGRAMS_DISPLACEMENT;
  mp->save_seismograms_v = *SAVE_SEISMOGRAMS_VELOCITY;
  mp->save_seismograms_a = *SAVE_SEISMOGRAMS_ACCELERATION;
  mp->save_seismograms_p = *SAVE_SEISMOGRAMS_PRESSURE;

  // receiver stations
  mp->nrec_local = *nrec_local; // number of receiver located in this partition

  // note: for adjoint simulations (SIMULATION_TYPE == 2),
  //         nrec_local     - is set to the number of sources (CMTSOLUTIONs), which act as "adjoint receiver" locations
  //                          for storing seismograms or strains
  //
  //         nadj_rec_local - determines the number of "adjoint sources", i.e., number of station locations (STATIONS_ADJOINT),
  //                          which act as sources to drive the adjoint wavefield
  //
  //         still, size(ispec_selected_rec) = nrec
  //
  //         hxir,.. arrays are interpolators for: - receiver locations (STATIONS) in case SIMULATION_TYPE == 1 or 3,
  //                                               - "adjoint receiver" locations (CMTSOLUTIONs) in case SIMULATION_TYPE == 2
  if (mp->nrec_local > 0){
    gpuCreateCopy_todevice_realw((void**)&mp->d_hxir,h_xir,NGLLX*mp->nrec_local);
    gpuCreateCopy_todevice_realw((void**)&mp->d_hetar,h_etar,NGLLY*mp->nrec_local);
    gpuCreateCopy_todevice_realw((void**)&mp->d_hgammar,h_gammar,NGLLZ*mp->nrec_local);

    // seismograms
    int size =  (*nlength_seismogram) * (*nrec_local);

    if (mp->save_seismograms_d)
      gpuMalloc_realw((void**)&mp->d_seismograms_d,NDIM * size);
    if (mp->save_seismograms_v)
      gpuMalloc_realw((void**)&mp->d_seismograms_v,NDIM * size);
    if (mp->save_seismograms_a)
      gpuMalloc_realw((void**)&mp->d_seismograms_a,NDIM * size);
    if (mp->save_seismograms_p)
      gpuMalloc_field((void**)&mp->d_seismograms_p,size);

    // stores only local receiver rotations in d_nu_rec
    realw* h_nu_rec;
    h_nu_rec = (realw*) calloc(NDIM * NDIM * mp->nrec_local, sizeof(realw));
    int irec_loc = 0;
    if (mp->simulation_type == 1 || mp->simulation_type == 3){
      // forward/kernel simulations: receiver positions at STATIONS locations
      //                             nu_rec,.. are for actual receiver locations
      //                             seismograms are taken at these receiver locations (specified by STATIONS)
      for (int irec=0; irec < (*nrec); irec++){
        if (mp->myrank == h_islice_selected_rec[irec]){
          for (int j = 0; j < 9; j++) h_nu_rec[j + NDIM * NDIM * irec_loc] = (realw)nu_rec[j + NDIM * NDIM * irec];
          irec_loc = irec_loc + 1;
        }
      }
    }else{
      // "pure" adjoint simulation: "adjoint receivers" are located at CMTSOLUTION source locations
      //                            nu_source,.. are for "adjoint receivers" locations
      //                            seismograms are taken at the "adjoint receivers" location (specified by CMTSOLUTION)
      for(int irec=0; irec < (*NSOURCES); irec++) {
        if (mp->myrank == h_islice_selected_source[irec]){
          for (int j = 0; j < 9; j++) h_nu_rec[j + NDIM * NDIM * irec_loc] = (realw)nu_source[j + NDIM * NDIM * irec];
          irec_loc = irec_loc + 1;
        }
      }
    }
    // checks
    if (irec_loc != mp->nrec_local) exit_on_error("prepare_constants_device: nrec_local not equal for d_nu_rec\n");
    // allocates on device
    gpuCreateCopy_todevice_realw((void**)&mp->d_nu_rec,h_nu_rec,NDIM * NDIM * mp->nrec_local);
    free(h_nu_rec);

    // stores only local receiver array
    int *ispec_selected_rec_loc;
    ispec_selected_rec_loc = (int*) calloc(mp->nrec_local, sizeof(int));
    irec_loc = 0;
    if (mp->simulation_type == 1 || mp->simulation_type == 3){
      // forward/kernel simulations: receiver positions at STATIONS locations
      //                             xir_store,.. are for actual receiver locations
      //                             seismograms are taken at these receiver locations (specified by STATIONS)
      for(int i=0; i < (*nrec); i++) {
        if (mp->myrank == h_islice_selected_rec[i]){
          ispec_selected_rec_loc[irec_loc] = h_ispec_selected_rec[i];
          irec_loc = irec_loc+1;
        }
      }
    }else{
      // "pure" adjoint simulation: "adjoint receivers" are located at CMTSOLUTION source locations
      //                            xir_store,.. are for "adjoint receivers" locations
      //                            seismograms are taken at the "adjoint receivers" location (specified by CMTSOLUTION)
      for(int i=0; i < (*NSOURCES); i++) {
        if (mp->myrank == h_islice_selected_source[i]){
          ispec_selected_rec_loc[irec_loc] = h_ispec_selected_source[i];
          irec_loc = irec_loc+1;
        }
      }
    }
    // checks
    if (irec_loc != mp->nrec_local) exit_on_error("prepare_constants_device: nrec_local not equal for d_ispec_selected_rec_loc\n");
    // allocates on device
    gpuCreateCopy_todevice_int((void**)&mp->d_ispec_selected_rec_loc,ispec_selected_rec_loc,mp->nrec_local);
    free(ispec_selected_rec_loc);
  }
  gpuCreateCopy_todevice_int((void**)&mp->d_ispec_selected_rec,h_ispec_selected_rec,(*nrec));

#ifdef USE_MESH_COLORING_GPUX
  mp->use_mesh_coloring_gpu = 1;
  if (! *USE_MESH_COLORING_GPU_f) exit_on_error("error with USE_MESH_COLORING_GPU constant; please re-compile\n");
#else
  // mesh coloring
  // note: this here passes the coloring as an option to the kernel routines
  //          the performance seems to be the same if one uses the pre-processing directives above or not
  mp->use_mesh_coloring_gpu = *USE_MESH_COLORING_GPU_f;
#endif

  // number of elements per domain
  mp->nspec_acoustic = *nspec_acoustic;
  mp->nspec_elastic = *nspec_elastic;

  // element flags (always needed for seismogram routines)
  gpuCreateCopy_todevice_int((void**)&mp->d_ispec_is_acoustic,ispec_is_acoustic,mp->NSPEC_AB);
  gpuCreateCopy_todevice_int((void**)&mp->d_ispec_is_elastic,ispec_is_elastic,mp->NSPEC_AB);

  // gravity flag initialization
  mp->gravity = 0;

  // fault simulation
  mp->fault_simulation = *FAULT_SIMULATION;
  // Kelvin_voigt initialization
  mp->use_Kelvin_Voigt_damping = 0;

  // JC JC here we will need to add GPU support for the new C-PML routines

  GPU_ERROR_CHECKING("prepare_constants_device");
}


/* ----------------------------------------------------------------------------------------------- */

// for ACOUSTIC simulations

/* ----------------------------------------------------------------------------------------------- */

extern EXTERN_LANG
void FC_FUNC_(prepare_fields_acoustic_device,
              PREPARE_FIELDS_ACOUSTIC_DEVICE)(long* Mesh_pointer,
                                              realw* rmass_acoustic, realw* rhostore, realw* kappastore,
                                              int* num_phase_ispec_acoustic, int* phase_ispec_inner_acoustic,
                                              int* NOISE_TOMOGRAPHY,
                                              int* num_free_surface_faces,
                                              int* free_surface_ispec,
                                              int* free_surface_ijk,
                                              int* b_reclen_potential, realw* b_absorb_potential,
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

  Mesh* mp = (Mesh*)(*Mesh_pointer);

  // allocates arrays on device (GPU)
  int size = mp->NGLOB_AB;
  gpuMalloc_field((void**)&(mp->d_potential_acoustic),size);
  gpuMalloc_field((void**)&(mp->d_potential_dot_acoustic),size);
  gpuMalloc_field((void**)&(mp->d_potential_dot_dot_acoustic),size);
  // initializes values to zero
  //gpuMemset_field(mp->d_potential_acoustic,size,0);
  //gpuMemset_field(mp->d_potential_dot_acoustic,size,0);
  //gpuMemset_field(mp->d_potential_dot_dot_acoustic,size,0);

  #ifdef USE_TEXTURES_FIELDS
  {
#ifdef USE_CUDA
    #ifdef USE_OLDER_CUDA4_GPU
      cudaChannelFormatDesc channelDesc = cudaCreateChannelDesc<float>();
      const textureReference* d_potential_tex_ref_ptr;
      print_CUDA_error_if_any(cudaGetTextureReference(&d_potential_tex_ref_ptr, "d_potential_tex"), 2221);
      print_CUDA_error_if_any(cudaBindTexture(0, d_potential_tex_ref_ptr, mp->d_potential_acoustic, &channelDesc, sizeof(realw)*size), 2001);

      const textureReference* d_potential_dot_dot_tex_ref_ptr;
      print_CUDA_error_if_any(cudaGetTextureReference(&d_potential_dot_dot_tex_ref_ptr, "d_potential_dot_dot_tex"), 2222);
      print_CUDA_error_if_any(cudaBindTexture(0, d_potential_dot_dot_tex_ref_ptr, mp->d_potential_dot_dot_acoustic, &channelDesc, sizeof(realw)*size), 2003);
    #else
      cudaChannelFormatDesc channelDesc = cudaCreateChannelDesc<float>();
      print_CUDA_error_if_any(cudaBindTexture(0, &d_potential_tex, mp->d_potential_acoustic, &channelDesc, sizeof(realw)*size), 2221);
      print_CUDA_error_if_any(cudaBindTexture(0, &d_potential_dot_dot_tex, mp->d_potential_dot_dot_acoustic, &channelDesc, sizeof(realw)*size), 2003);
    #endif
#endif
#ifdef USE_HIP
      hipChannelFormatDesc channelDesc = hipCreateChannelDesc<float>();
      print_HIP_error_if_any(hipBindTexture(0, &d_potential_tex, mp->d_potential_acoustic, &channelDesc, sizeof(realw)*size), 2221);
      print_HIP_error_if_any(hipBindTexture(0, &d_potential_dot_dot_tex, mp->d_potential_dot_dot_acoustic, &channelDesc, sizeof(realw)*size), 2003);
#endif
  }
  #endif

  // mpi buffer
  mp->size_mpi_buffer_potential = (mp->num_interfaces_ext_mesh) * (mp->max_nibool_interfaces_ext_mesh);
  if (mp->size_mpi_buffer_potential > 0){
    gpuMalloc_field((void**)&(mp->d_send_potential_dot_dot_buffer),mp->size_mpi_buffer_potential);
  }

  // mass matrix
  gpuCreateCopy_todevice_realw((void**)&mp->d_rmass_acoustic,rmass_acoustic,mp->NGLOB_AB);

  // density
  // padded array
  // Assuming NGLLX==5. Padded is then 128 (5^3+3)
  int size_padded = NGLL3_PADDED * mp->NSPEC_AB;
  gpuMalloc_realw((void**)&(mp->d_rhostore),size_padded);
  // transfer constant element data with padding
  /*
  // way 1: slow...
  for(int i=0; i < mp->NSPEC_AB; i++) {
    gpuMemcpy_todevice_realw(mp->d_rhostore+i*NGLL3_PADDED, &rhostore[i*NGLL3],NGLL3);
  }
  */
  // way 2: faster ...
  gpuMemcpy2D_todevice_realw(mp->d_rhostore, NGLL3_PADDED, rhostore, NGLL3, NGLL3, mp->NSPEC_AB);

  // non-padded array
  gpuCreateCopy_todevice_realw((void**)&mp->d_kappastore,kappastore,NGLL3*mp->NSPEC_AB);

  // phase elements
  mp->num_phase_ispec_acoustic = *num_phase_ispec_acoustic;
  gpuCreateCopy_todevice_int((void**)&mp->d_phase_ispec_inner_acoustic,phase_ispec_inner_acoustic,
                       2*mp->num_phase_ispec_acoustic);

  // free surface
  if (*NOISE_TOMOGRAPHY == 0){
    // allocate surface arrays
    mp->num_free_surface_faces = *num_free_surface_faces;
    if (mp->num_free_surface_faces > 0){
      gpuCreateCopy_todevice_int((void**)&mp->d_free_surface_ispec,free_surface_ispec,mp->num_free_surface_faces);
      gpuCreateCopy_todevice_int((void**)&mp->d_free_surface_ijk,free_surface_ijk,3*NGLL2*mp->num_free_surface_faces);
    }
  }

  // absorbing boundaries
  if (mp->absorbing_conditions && mp->d_num_abs_boundary_faces > 0){
    // absorb_field array used for file i/o
    if (mp->simulation_type == 3 || ( mp->simulation_type == 1 && mp->save_forward )){
      // note: b_reclen_potential is record length in bytes ( CUSTOM_REAL * NGLLSQUARE * num_abs_boundary_faces )
      mp->d_b_reclen_potential = *b_reclen_potential;

      // note: reclen includes sizeof(realw), thus for gpuMalloc_realw the size needs to be divided
      gpuMalloc_realw((void**)&mp->d_b_absorb_potential,mp->d_b_reclen_potential/sizeof(realw));
      // note: for copying with gpuMemcpy_**_void the actualy byte size is used, thus no need to divide here by sizeof(realw)
      gpuMemcpy_todevice_void((void*)mp->d_b_absorb_potential,(void*)b_absorb_potential,mp->d_b_reclen_potential);
    }
  }

  // coupling with elastic parts
  if (*ELASTIC_SIMULATION && *num_coupling_ac_el_faces > 0){
    gpuCreateCopy_todevice_int((void**)&mp->d_coupling_ac_el_ispec,coupling_ac_el_ispec,(*num_coupling_ac_el_faces));
    gpuCreateCopy_todevice_int((void**)&mp->d_coupling_ac_el_ijk,coupling_ac_el_ijk,3*NGLL2*(*num_coupling_ac_el_faces));
    gpuCreateCopy_todevice_realw((void**)&mp->d_coupling_ac_el_normal,coupling_ac_el_normal,3*NGLL2*(*num_coupling_ac_el_faces));
    gpuCreateCopy_todevice_realw((void**)&mp->d_coupling_ac_el_jacobian2Dw,coupling_ac_el_jacobian2Dw,NGLL2*(*num_coupling_ac_el_faces));
  }

  // mesh coloring
  if (mp->use_mesh_coloring_gpu ){
    mp->num_colors_outer_acoustic = *num_colors_outer_acoustic;
    mp->num_colors_inner_acoustic = *num_colors_inner_acoustic;
    mp->h_num_elem_colors_acoustic = (int*) num_elem_colors_acoustic;
  }

  GPU_ERROR_CHECKING("prepare_fields_acoustic_device");
}


/* ----------------------------------------------------------------------------------------------- */

extern EXTERN_LANG
void FC_FUNC_(prepare_fields_acoustic_adj_dev,
              PREPARE_FIELDS_ACOUSTIC_ADJ_DEV)(long* Mesh_pointer,
                                               int* APPROXIMATE_HESS_KL) {

  TRACE("prepare_fields_acoustic_adj_dev");

  Mesh* mp = (Mesh*)(*Mesh_pointer);

  // kernel simulations
  if (mp->simulation_type != 3) return;

  // allocates backward/reconstructed arrays on device (GPU)
  int size = mp->NGLOB_AB;
  gpuMalloc_field((void**)&(mp->d_b_potential_acoustic),size);
  gpuMalloc_field((void**)&(mp->d_b_potential_dot_acoustic),size);
  gpuMalloc_field((void**)&(mp->d_b_potential_dot_dot_acoustic),size);
  // initializes values to zero
  //gpuMemset_field(mp->d_b_potential_acoustic,size,0);
  //gpuMemset_field(mp->d_b_potential_dot_acoustic,size,0);
  //gpuMemset_field(mp->d_b_potential_dot_dot_acoustic,size,0);

  #ifdef USE_TEXTURES_FIELDS
  {
#ifdef USE_CUDA
    #ifdef USE_OLDER_CUDA4_GPU
      cudaChannelFormatDesc channelDesc = cudaCreateChannelDesc<float>();
      const textureReference* d_b_potential_tex_ref_ptr;
      print_CUDA_error_if_any(cudaGetTextureReference(&d_b_potential_tex_ref_ptr, "d_b_potential_tex"), 3001);
      print_CUDA_error_if_any(cudaBindTexture(0, d_b_potential_tex_ref_ptr, mp->d_b_potential_acoustic, &channelDesc, sizeof(realw)*size), 3001);

      const textureReference* d_b_potential_dot_dot_tex_ref_ptr;
      print_CUDA_error_if_any(cudaGetTextureReference(&d_b_potential_dot_dot_tex_ref_ptr, "d_b_potential_dot_dot_tex"),3003);
      print_CUDA_error_if_any(cudaBindTexture(0, d_b_potential_dot_dot_tex_ref_ptr, mp->d_b_potential_dot_dot_acoustic, &channelDesc, sizeof(realw)*size), 3003);
    #else
      cudaChannelFormatDesc channelDesc = cudaCreateChannelDesc<float>();
      print_CUDA_error_if_any(cudaBindTexture(0, &d_b_potential_tex, mp->d_b_potential_acoustic, &channelDesc, sizeof(realw)*size), 3001);
      print_CUDA_error_if_any(cudaBindTexture(0, &d_b_potential_dot_dot_tex, mp->d_b_potential_dot_dot_acoustic, &channelDesc, sizeof(realw)*size), 3003);
    #endif
#endif
#ifdef USE_HIP
      hipChannelFormatDesc channelDesc = hipCreateChannelDesc<float>();
      print_HIP_error_if_any(hipBindTexture(0, &d_b_potential_tex, mp->d_b_potential_acoustic, &channelDesc, sizeof(realw)*size), 3001);
      print_HIP_error_if_any(hipBindTexture(0, &d_b_potential_dot_dot_tex, mp->d_b_potential_dot_dot_acoustic, &channelDesc, sizeof(realw)*size), 3003);
#endif
  }
  #endif

  // allocates kernels
  size = NGLL3*mp->NSPEC_AB;
  gpuMalloc_realw((void**)&(mp->d_rho_ac_kl),size);
  gpuMalloc_realw((void**)&(mp->d_kappa_ac_kl),size);
  // initializes kernel values to zero
  gpuMemset_realw(mp->d_rho_ac_kl,size,0);
  gpuMemset_realw(mp->d_kappa_ac_kl,size,0);

  // preconditioner
  mp->approximate_hess_kl = *APPROXIMATE_HESS_KL;
  if (mp->approximate_hess_kl){
    gpuMalloc_realw((void**)&(mp->d_hess_ac_kl),size);
    gpuMalloc_realw((void**)&(mp->d_hess_rho_ac_kl),size);
    gpuMalloc_realw((void**)&(mp->d_hess_kappa_ac_kl),size);

    // initializes with zeros
    gpuMemset_realw(mp->d_hess_ac_kl,size,0);
    gpuMemset_realw(mp->d_hess_rho_ac_kl,size,0);
    gpuMemset_realw(mp->d_hess_kappa_ac_kl,size,0);
  }

  // mpi buffer
  if (mp->size_mpi_buffer_potential > 0){
    gpuMalloc_field((void**)&(mp->d_b_send_potential_dot_dot_buffer),mp->size_mpi_buffer_potential);
  }

  GPU_ERROR_CHECKING("prepare_fields_acoustic_adj_dev");
}


/* ----------------------------------------------------------------------------------------------- */

// for ELASTIC simulations

/* ----------------------------------------------------------------------------------------------- */

extern EXTERN_LANG
void FC_FUNC_(prepare_fields_elastic_device,
              PREPARE_FIELDS_ELASTIC_DEVICE)(long* Mesh_pointer,
                                             realw* rmassx, realw* rmassy, realw* rmassz,
                                             realw* rho_vp, realw* rho_vs,
                                             realw* h_kappav, realw* h_muv,
                                             int* num_phase_ispec_elastic,
                                             int* phase_ispec_inner_elastic,
                                             realw* b_absorb_field, int* b_reclen_field,
                                             int* COMPUTE_AND_STORE_STRAIN,
                                             realw* epsilondev_xx,realw* epsilondev_yy,realw* epsilondev_xy,
                                             realw* epsilondev_xz,realw* epsilondev_yz,
                                             int* ATTENUATION,
                                             int* R_size,
                                             realw* R_xx,realw* R_yy,realw* R_xy,realw* R_xz,realw* R_yz,
                                             realw* factor_common,
                                             realw* R_trace,realw* epsilondev_trace,
                                             realw* factor_common_kappa,
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
                                             realw *c11store,realw *c12store,realw *c13store,
                                             realw *c14store,realw *c15store,realw *c16store,
                                             realw *c22store,realw *c23store,realw *c24store,
                                             realw *c25store,realw *c26store,realw *c33store,
                                             realw *c34store,realw *c35store,realw *c36store,
                                             realw *c44store,realw *c45store,realw *c46store,
                                             realw *c55store,realw *c56store,realw *c66store ){

  TRACE("prepare_fields_elastic_device");

  Mesh* mp = (Mesh*)(*Mesh_pointer);
  int size;

  // debug
  //printf("prepare_fields_elastic_device: rank %d - wavefield setup\n",mp->myrank);
  //synchronize_mpi();

  // elastic wavefields
  size = NDIM * mp->NGLOB_AB;
  gpuMalloc_realw((void**)&(mp->d_displ),size);
  gpuMalloc_realw((void**)&(mp->d_veloc),size);
  gpuMalloc_realw((void**)&(mp->d_accel),size);
  // initializes values to zero
  //gpuMemset_realw(mp->d_displ,size,0);
  //gpuMemset_realw(mp->d_veloc,size,0);
  //gpuMemset_realw(mp->d_accel,size,0);

  #ifdef USE_TEXTURES_FIELDS
  {
#ifdef USE_CUDA
    #ifdef USE_OLDER_CUDA4_GPU
      cudaChannelFormatDesc channelDesc = cudaCreateChannelDesc<float>();
      const textureReference* d_displ_tex_ref_ptr;
      print_CUDA_error_if_any(cudaGetTextureReference(&d_displ_tex_ref_ptr, "d_displ_tex"), 4001);
      print_CUDA_error_if_any(cudaBindTexture(0, d_displ_tex_ref_ptr, mp->d_displ, &channelDesc, sizeof(realw)*size), 4001);

      const textureReference* d_veloc_tex_ref_ptr;
      print_CUDA_error_if_any(cudaGetTextureReference(&d_veloc_tex_ref_ptr, "d_veloc_tex"), 4002);
      print_CUDA_error_if_any(cudaBindTexture(0, d_veloc_tex_ref_ptr, mp->d_veloc, &channelDesc, sizeof(realw)*size), 4002);

      const textureReference* d_accel_tex_ref_ptr;
      print_CUDA_error_if_any(cudaGetTextureReference(&d_accel_tex_ref_ptr, "d_accel_tex"), 4003);
      print_CUDA_error_if_any(cudaBindTexture(0, d_accel_tex_ref_ptr, mp->d_accel, &channelDesc, sizeof(realw)*size), 4003);
    #else
      cudaChannelFormatDesc channelDesc = cudaCreateChannelDesc<float>();
      print_CUDA_error_if_any(cudaBindTexture(0, &d_displ_tex, mp->d_displ, &channelDesc, sizeof(realw)*size), 4001);
      print_CUDA_error_if_any(cudaBindTexture(0, &d_veloc_tex, mp->d_veloc, &channelDesc, sizeof(realw)*size), 4002);
      print_CUDA_error_if_any(cudaBindTexture(0, &d_accel_tex, mp->d_accel, &channelDesc, sizeof(realw)*size), 4003);
    #endif
#endif
#ifdef USE_HIP
      hipChannelFormatDesc channelDesc = hipCreateChannelDesc<float>();
      print_HIP_error_if_any(hipBindTexture(0, &d_displ_tex, mp->d_displ, &channelDesc, sizeof(realw)*size), 4001);
      print_HIP_error_if_any(hipBindTexture(0, &d_veloc_tex, mp->d_veloc, &channelDesc, sizeof(realw)*size), 4002);
      print_HIP_error_if_any(hipBindTexture(0, &d_accel_tex, mp->d_accel, &channelDesc, sizeof(realw)*size), 4003);
#endif
  }
  #endif

  // debug
  //synchronize_mpi();

  // MPI buffer
  mp->size_mpi_buffer = NDIM * (mp->num_interfaces_ext_mesh) * (mp->max_nibool_interfaces_ext_mesh);
  if (mp->size_mpi_buffer > 0){
    // note: Allocate pinned mpi-buffers.
    //       MPI buffers use pinned memory allocated by cudaMallocHost, which
    //       enables the use of asynchronous memory copies from host <-> device

    // send buffer
    gpuMallocHost_realw((void**)&(mp->h_send_accel_buffer),mp->size_mpi_buffer);

    // unused so far..
    //mp->send_buffer = (realw*)malloc((mp->size_mpi_buffer)*sizeof(realw));
    // extra buffer for adjoint, not needed so far..., can use the same buffer for both forward/adjoint mpi exchanges
    //gpuMallocHost_realw((void**)&(mp->h_send_b_accel_buffer),mp->size_mpi_buffer);
    //mp->b_send_buffer = (realw*)malloc((size_mpi_buffer)*sizeof(realw));

    // receive buffer
    gpuMallocHost_realw((void**)&(mp->h_recv_accel_buffer),mp->size_mpi_buffer);

    // unused so far..
    //mp->recv_buffer = (realw*) malloc((mp->size_mpi_buffer)*sizeof(realw));

    // non-pinned buffer
    gpuMalloc_realw((void**)&(mp->d_send_accel_buffer),mp->size_mpi_buffer);
    // adjoint
    if (mp->simulation_type == 3){
      gpuMalloc_realw((void**)&(mp->d_b_send_accel_buffer),mp->size_mpi_buffer);
    }
  }

  // debug
  //printf("prepare_fields_elastic_device: rank %d - mass matrix\n",mp->myrank);
  //synchronize_mpi();

  // mass matrix
  gpuCreateCopy_todevice_realw((void**)&mp->d_rmassx,rmassx,mp->NGLOB_AB);
  gpuCreateCopy_todevice_realw((void**)&mp->d_rmassy,rmassy,mp->NGLOB_AB);
  gpuCreateCopy_todevice_realw((void**)&mp->d_rmassz,rmassz,mp->NGLOB_AB);

  // phase elements
  mp->num_phase_ispec_elastic = *num_phase_ispec_elastic;

  gpuCreateCopy_todevice_int((void**)&mp->d_phase_ispec_inner_elastic,phase_ispec_inner_elastic,2*mp->num_phase_ispec_elastic);

  // debug
  //synchronize_mpi();

  // absorbing conditions
  if (mp->absorbing_conditions && mp->d_num_abs_boundary_faces > 0){

    // debug
    //printf("prepare_fields_elastic_device: rank %d - absorbing boundary setup\n",mp->myrank);

    // non-padded arrays
    // rho_vp, rho_vs non-padded; they are needed for stacey boundary condition
    gpuCreateCopy_todevice_realw((void**)&mp->d_rho_vp,rho_vp,NGLL3*mp->NSPEC_AB);
    gpuCreateCopy_todevice_realw((void**)&mp->d_rho_vs,rho_vs,NGLL3*mp->NSPEC_AB);

    // absorb_field array used for file i/o
    if (mp->simulation_type == 3 || ( mp->simulation_type == 1 && mp->save_forward )){
      // note: b_reclen_field is length in bytes already (CUSTOM_REAL * NDIM * NGLLSQUARE * num_abs_boundary_faces )
      mp->d_b_reclen_field = *b_reclen_field;

      // debug
      //printf("prepare_fields_elastic_device: rank %d - absorbing boundary i/o %d\n",mp->myrank,mp->d_b_reclen_field);

      // note: reclen includes sizeof(realw), thus for gpuMalloc_realw the size needs to be divided
      gpuMalloc_realw((void**)&mp->d_b_absorb_field,mp->d_b_reclen_field/sizeof(realw));
      // note: for copying with gpuMemcpy_**_void the actualy byte size is used, thus no need to divide here by sizeof(realw)
      gpuMemcpy_todevice_void((void*)mp->d_b_absorb_field,(void*)b_absorb_field,mp->d_b_reclen_field);
    }
  }
  int size_padded = NGLL3_PADDED*mp->NSPEC_AB;
  gpuMalloc_realw((void**) &mp->d_kappav, size_padded);
  gpuMalloc_realw((void**) &mp->d_muv, size_padded);
  gpuMemcpy2D_todevice_realw(mp->d_kappav, NGLL3_PADDED, h_kappav, NGLL3, NGLL3, mp->NSPEC_AB);
  gpuMemcpy2D_todevice_realw(mp->d_muv, NGLL3_PADDED, h_muv, NGLL3, NGLL3, mp->NSPEC_AB);

  // debug
  //synchronize_mpi();

  // strains used for attenuation and kernel simulations
  mp->compute_and_store_strain = *COMPUTE_AND_STORE_STRAIN;
  if (mp->compute_and_store_strain){
    // debug
    //printf("prepare_fields_elastic_device: rank %d - strain setup\n",mp->myrank);
    //synchronize_mpi();

    // strains
    size = NGLL3 * mp->NSPEC_AB; // note: non-aligned; if align, check memcpy below and indexing
    gpuCreateCopy_todevice_realw((void**)&mp->d_epsilondev_xx,epsilondev_xx,size);
    gpuCreateCopy_todevice_realw((void**)&mp->d_epsilondev_yy,epsilondev_yy,size);
    gpuCreateCopy_todevice_realw((void**)&mp->d_epsilondev_xy,epsilondev_xy,size);
    gpuCreateCopy_todevice_realw((void**)&mp->d_epsilondev_xz,epsilondev_xz,size);
    gpuCreateCopy_todevice_realw((void**)&mp->d_epsilondev_yz,epsilondev_yz,size);
    gpuCreateCopy_todevice_realw((void**)&mp->d_epsilondev_trace,epsilondev_trace,size);
  }

  // attenuation memory variables
  mp->attenuation = *ATTENUATION;
  if (mp->attenuation){
    // debug
    //printf("prepare_fields_elastic_device: rank %d - attenuation setup\n",mp->myrank);
    //synchronize_mpi();

    // memory arrays
    size = *R_size;
    gpuCreateCopy_todevice_realw((void**)&mp->d_R_xx,R_xx,size);
    gpuCreateCopy_todevice_realw((void**)&mp->d_R_yy,R_yy,size);
    gpuCreateCopy_todevice_realw((void**)&mp->d_R_xy,R_xy,size);
    gpuCreateCopy_todevice_realw((void**)&mp->d_R_xz,R_xz,size);
    gpuCreateCopy_todevice_realw((void**)&mp->d_R_yz,R_yz,size);
    gpuCreateCopy_todevice_realw((void**)&mp->d_R_trace,R_trace,size);

    // attenuation factors
    gpuCreateCopy_todevice_realw((void**)&mp->d_factor_common,factor_common,N_SLS*NGLL3*mp->NSPEC_AB);
    gpuCreateCopy_todevice_realw((void**)&mp->d_factor_common_kappa,factor_common_kappa,N_SLS*NGLL3*mp->NSPEC_AB);

    // alpha,beta,gamma factors
    gpuCreateCopy_todevice_realw((void**)&mp->d_alphaval,alphaval,N_SLS);
    gpuCreateCopy_todevice_realw((void**)&mp->d_betaval,betaval,N_SLS);
    gpuCreateCopy_todevice_realw((void**)&mp->d_gammaval,gammaval,N_SLS);
  }

  // anisotropy
  mp->ANISOTROPY = *ANISOTROPY;
  if (mp->ANISOTROPY){
    // debug
    //printf("prepare_fields_elastic_device: rank %d - attenuation setup\n",mp->myrank);
    //synchronize_mpi();

    // Assuming NGLLX==5. Padded is then 128 (5^3+3)
    size_padded = NGLL3_PADDED * (mp->NSPEC_AB);

    // allocates memory on GPU
    gpuMalloc_realw((void**)&(mp->d_c11store),size_padded);
    gpuMalloc_realw((void**)&(mp->d_c12store),size_padded);
    gpuMalloc_realw((void**)&(mp->d_c13store),size_padded);
    gpuMalloc_realw((void**)&(mp->d_c14store),size_padded);
    gpuMalloc_realw((void**)&(mp->d_c15store),size_padded);
    gpuMalloc_realw((void**)&(mp->d_c16store),size_padded);
    gpuMalloc_realw((void**)&(mp->d_c22store),size_padded);
    gpuMalloc_realw((void**)&(mp->d_c23store),size_padded);
    gpuMalloc_realw((void**)&(mp->d_c24store),size_padded);
    gpuMalloc_realw((void**)&(mp->d_c25store),size_padded);
    gpuMalloc_realw((void**)&(mp->d_c26store),size_padded);
    gpuMalloc_realw((void**)&(mp->d_c33store),size_padded);
    gpuMalloc_realw((void**)&(mp->d_c34store),size_padded);
    gpuMalloc_realw((void**)&(mp->d_c35store),size_padded);
    gpuMalloc_realw((void**)&(mp->d_c36store),size_padded);
    gpuMalloc_realw((void**)&(mp->d_c44store),size_padded);
    gpuMalloc_realw((void**)&(mp->d_c45store),size_padded);
    gpuMalloc_realw((void**)&(mp->d_c46store),size_padded);
    gpuMalloc_realw((void**)&(mp->d_c55store),size_padded);
    gpuMalloc_realw((void**)&(mp->d_c56store),size_padded);
    gpuMalloc_realw((void**)&(mp->d_c66store),size_padded);

    // transfer constant element data with padding
    /*
    // way 1: slower ...
    for(int i=0;i < mp->NSPEC_AB;i++) {
      gpuMemcpy_todevice_realw(mp->d_c11store + i*NGLL3_PADDED, &c11store[i*NGLL3], NGLL3);
      gpuMemcpy_todevice_realw(mp->d_c12store + i*NGLL3_PADDED, &c12store[i*NGLL3], NGLL3);
      gpuMemcpy_todevice_realw(mp->d_c13store + i*NGLL3_PADDED, &c13store[i*NGLL3], NGLL3);
      gpuMemcpy_todevice_realw(mp->d_c14store + i*NGLL3_PADDED, &c14store[i*NGLL3], NGLL3);
      gpuMemcpy_todevice_realw(mp->d_c15store + i*NGLL3_PADDED, &c15store[i*NGLL3], NGLL3);
      gpuMemcpy_todevice_realw(mp->d_c16store + i*NGLL3_PADDED, &c16store[i*NGLL3], NGLL3);
      gpuMemcpy_todevice_realw(mp->d_c22store + i*NGLL3_PADDED, &c22store[i*NGLL3], NGLL3);
      gpuMemcpy_todevice_realw(mp->d_c23store + i*NGLL3_PADDED, &c23store[i*NGLL3], NGLL3);
      gpuMemcpy_todevice_realw(mp->d_c24store + i*NGLL3_PADDED, &c24store[i*NGLL3], NGLL3);
      gpuMemcpy_todevice_realw(mp->d_c25store + i*NGLL3_PADDED, &c25store[i*NGLL3], NGLL3);
      gpuMemcpy_todevice_realw(mp->d_c26store + i*NGLL3_PADDED, &c26store[i*NGLL3], NGLL3);
      gpuMemcpy_todevice_realw(mp->d_c33store + i*NGLL3_PADDED, &c33store[i*NGLL3], NGLL3);
      gpuMemcpy_todevice_realw(mp->d_c34store + i*NGLL3_PADDED, &c34store[i*NGLL3], NGLL3);
      gpuMemcpy_todevice_realw(mp->d_c35store + i*NGLL3_PADDED, &c35store[i*NGLL3], NGLL3);
      gpuMemcpy_todevice_realw(mp->d_c36store + i*NGLL3_PADDED, &c36store[i*NGLL3], NGLL3);
      gpuMemcpy_todevice_realw(mp->d_c44store + i*NGLL3_PADDED, &c44store[i*NGLL3], NGLL3);
      gpuMemcpy_todevice_realw(mp->d_c45store + i*NGLL3_PADDED, &c45store[i*NGLL3], NGLL3);
      gpuMemcpy_todevice_realw(mp->d_c46store + i*NGLL3_PADDED, &c46store[i*NGLL3], NGLL3);
      gpuMemcpy_todevice_realw(mp->d_c55store + i*NGLL3_PADDED, &c55store[i*NGLL3], NGLL3);
      gpuMemcpy_todevice_realw(mp->d_c56store + i*NGLL3_PADDED, &c56store[i*NGLL3], NGLL3);
      gpuMemcpy_todevice_realw(mp->d_c66store + i*NGLL3_PADDED, &c66store[i*NGLL3], NGLL3);
    }
    */
    // way 2: faster ...
    gpuMemcpy2D_todevice_realw(mp->d_c11store, NGLL3_PADDED, c11store, NGLL3, NGLL3, mp->NSPEC_AB);
    gpuMemcpy2D_todevice_realw(mp->d_c12store, NGLL3_PADDED, c12store, NGLL3, NGLL3, mp->NSPEC_AB);
    gpuMemcpy2D_todevice_realw(mp->d_c13store, NGLL3_PADDED, c13store, NGLL3, NGLL3, mp->NSPEC_AB);
    gpuMemcpy2D_todevice_realw(mp->d_c14store, NGLL3_PADDED, c14store, NGLL3, NGLL3, mp->NSPEC_AB);
    gpuMemcpy2D_todevice_realw(mp->d_c15store, NGLL3_PADDED, c15store, NGLL3, NGLL3, mp->NSPEC_AB);
    gpuMemcpy2D_todevice_realw(mp->d_c16store, NGLL3_PADDED, c16store, NGLL3, NGLL3, mp->NSPEC_AB);
    gpuMemcpy2D_todevice_realw(mp->d_c22store, NGLL3_PADDED, c22store, NGLL3, NGLL3, mp->NSPEC_AB);
    gpuMemcpy2D_todevice_realw(mp->d_c23store, NGLL3_PADDED, c23store, NGLL3, NGLL3, mp->NSPEC_AB);
    gpuMemcpy2D_todevice_realw(mp->d_c24store, NGLL3_PADDED, c24store, NGLL3, NGLL3, mp->NSPEC_AB);
    gpuMemcpy2D_todevice_realw(mp->d_c25store, NGLL3_PADDED, c25store, NGLL3, NGLL3, mp->NSPEC_AB);
    gpuMemcpy2D_todevice_realw(mp->d_c26store, NGLL3_PADDED, c26store, NGLL3, NGLL3, mp->NSPEC_AB);
    gpuMemcpy2D_todevice_realw(mp->d_c33store, NGLL3_PADDED, c33store, NGLL3, NGLL3, mp->NSPEC_AB);
    gpuMemcpy2D_todevice_realw(mp->d_c34store, NGLL3_PADDED, c34store, NGLL3, NGLL3, mp->NSPEC_AB);
    gpuMemcpy2D_todevice_realw(mp->d_c35store, NGLL3_PADDED, c35store, NGLL3, NGLL3, mp->NSPEC_AB);
    gpuMemcpy2D_todevice_realw(mp->d_c36store, NGLL3_PADDED, c36store, NGLL3, NGLL3, mp->NSPEC_AB);
    gpuMemcpy2D_todevice_realw(mp->d_c44store, NGLL3_PADDED, c44store, NGLL3, NGLL3, mp->NSPEC_AB);
    gpuMemcpy2D_todevice_realw(mp->d_c45store, NGLL3_PADDED, c45store, NGLL3, NGLL3, mp->NSPEC_AB);
    gpuMemcpy2D_todevice_realw(mp->d_c46store, NGLL3_PADDED, c46store, NGLL3, NGLL3, mp->NSPEC_AB);
    gpuMemcpy2D_todevice_realw(mp->d_c55store, NGLL3_PADDED, c55store, NGLL3, NGLL3, mp->NSPEC_AB);
    gpuMemcpy2D_todevice_realw(mp->d_c56store, NGLL3_PADDED, c56store, NGLL3, NGLL3, mp->NSPEC_AB);
    gpuMemcpy2D_todevice_realw(mp->d_c66store, NGLL3_PADDED, c66store, NGLL3, NGLL3, mp->NSPEC_AB);
  }

  // ocean load approximation
  mp->approximate_ocean_load = *APPROXIMATE_OCEAN_LOAD;
  if (mp->approximate_ocean_load){
    // debug
    //printf("prepare_fields_elastic_device: rank %d - ocean load setup\n",mp->myrank);
    //synchronize_mpi();

    // oceans needs a free surface
    mp->num_free_surface_faces = *num_free_surface_faces;
    if (mp->num_free_surface_faces > 0){
      // mass matrix
      gpuCreateCopy_todevice_realw((void**)&mp->d_rmass_ocean_load,rmass_ocean_load,mp->NGLOB_AB);
      // surface normal
      gpuCreateCopy_todevice_realw((void**)&mp->d_free_surface_normal,free_surface_normal,3*NGLL2*(mp->num_free_surface_faces));
      // temporary global array: used to synchronize updates on global accel array
      gpuMalloc_int((void**)&(mp->d_updated_dof_ocean_load),mp->NGLOB_AB);

      if (*NOISE_TOMOGRAPHY == 0 && *ACOUSTIC_SIMULATION == 0){
        gpuCreateCopy_todevice_int((void**)&mp->d_free_surface_ispec,free_surface_ispec,mp->num_free_surface_faces);
        gpuCreateCopy_todevice_int((void**)&mp->d_free_surface_ijk,free_surface_ijk,
                             3*NGLL2*mp->num_free_surface_faces);
      }
    }
  }

  // mesh coloring
  if (mp->use_mesh_coloring_gpu ){
    mp->num_colors_outer_elastic = *num_colors_outer_elastic;
    mp->num_colors_inner_elastic = *num_colors_inner_elastic;
    mp->h_num_elem_colors_elastic = (int*) num_elem_colors_elastic;
  }

  // JC JC here we will need to add GPU support for the new C-PML routines

  // debug
  //printf("prepare_fields_elastic_device: rank %d - done\n",mp->myrank);
  //synchronize_mpi();

  GPU_ERROR_CHECKING("prepare_fields_elastic_device");
}

/* ----------------------------------------------------------------------------------------------- */

extern EXTERN_LANG
void FC_FUNC_(prepare_fields_elastic_adj_dev,
              PREPARE_FIELDS_ELASTIC_ADJ_DEV)(long* Mesh_pointer,
                                             int* size_f,
                                             realw* epsilon_trace_over_3,
                                             realw* b_epsilondev_xx,realw* b_epsilondev_yy,realw* b_epsilondev_xy,
                                             realw* b_epsilondev_xz,realw* b_epsilondev_yz,
                                             realw* b_epsilon_trace_over_3,
                                             int* R_size,
                                             realw* b_R_xx,realw* b_R_yy,realw* b_R_xy,realw* b_R_xz,realw* b_R_yz,
                                             realw* b_R_trace,realw* b_epsilondev_trace,
                                             realw* b_alphaval,realw* b_betaval,realw* b_gammaval,
                                             int* ANISOTROPIC_KL,
                                             int* APPROXIMATE_HESS_KL){

  TRACE("prepare_fields_elastic_adj_dev");

  Mesh* mp = (Mesh*)(*Mesh_pointer);
  int size;

  // checks if kernel simulation
  if (mp->simulation_type != 3) return;

  // kernel simulations
  // debug
  //printf("prepare_fields_elastic_adj_dev: rank %d - kernel setup\n",mp->myrank);
  //synchronize_mpi();

  // backward/reconstructed wavefields
  // allocates backward/reconstructed arrays on device (GPU)
  size = *size_f;
  gpuMalloc_realw((void**)&(mp->d_b_displ),size);
  gpuMalloc_realw((void**)&(mp->d_b_veloc),size);
  gpuMalloc_realw((void**)&(mp->d_b_accel),size);
  // initializes values to zero
  //gpuMemset_realw(mp->d_b_displ,size,0);
  //gpuMemset_realw(mp->d_b_veloc,size,0);
  //gpuMemset_realw(mp->d_b_accel,size,0);

  #ifdef USE_TEXTURES_FIELDS
  {
#ifdef USE_CUDA
    #ifdef USE_OLDER_CUDA4_GPU
      cudaChannelFormatDesc channelDesc = cudaCreateChannelDesc<float>();
      const textureReference* d_b_displ_tex_ref_ptr;
      print_CUDA_error_if_any(cudaGetTextureReference(&d_b_displ_tex_ref_ptr, "d_b_displ_tex"), 4001);
      print_CUDA_error_if_any(cudaBindTexture(0, d_b_displ_tex_ref_ptr, mp->d_b_displ, &channelDesc, sizeof(realw)*size), 4001);

      const textureReference* d_b_veloc_tex_ref_ptr;
      print_CUDA_error_if_any(cudaGetTextureReference(&d_b_veloc_tex_ref_ptr, "d_b_veloc_tex"), 4002);
      print_CUDA_error_if_any(cudaBindTexture(0, d_b_veloc_tex_ref_ptr, mp->d_b_veloc, &channelDesc, sizeof(realw)*size), 4002);

      const textureReference* d_b_accel_tex_ref_ptr;
      print_CUDA_error_if_any(cudaGetTextureReference(&d_b_accel_tex_ref_ptr, "d_b_accel_tex"), 4003);
      print_CUDA_error_if_any(cudaBindTexture(0, d_b_accel_tex_ref_ptr, mp->d_b_accel, &channelDesc, sizeof(realw)*size), 4003);
    #else
      cudaChannelFormatDesc channelDesc = cudaCreateChannelDesc<float>();
      print_CUDA_error_if_any(cudaBindTexture(0, &d_b_displ_tex, mp->d_b_displ, &channelDesc, sizeof(realw)*size), 4001);
      print_CUDA_error_if_any(cudaBindTexture(0, &d_b_veloc_tex, mp->d_b_veloc, &channelDesc, sizeof(realw)*size), 4002);
      print_CUDA_error_if_any(cudaBindTexture(0, &d_b_accel_tex, mp->d_b_accel, &channelDesc, sizeof(realw)*size), 4003);
    #endif
#endif
#ifdef USE_HIP
      hipChannelFormatDesc channelDesc = hipCreateChannelDesc<float>();
      print_HIP_error_if_any(hipBindTexture(0, &d_b_displ_tex, mp->d_b_displ, &channelDesc, sizeof(realw)*size), 4001);
      print_HIP_error_if_any(hipBindTexture(0, &d_b_veloc_tex, mp->d_b_veloc, &channelDesc, sizeof(realw)*size), 4002);
      print_HIP_error_if_any(hipBindTexture(0, &d_b_accel_tex, mp->d_b_accel, &channelDesc, sizeof(realw)*size), 4003);
#endif
  }
  #endif


  // anisotropic kernel flag
  mp->anisotropic_kl = *ANISOTROPIC_KL;

  // anisotropic/isotropic kernels
  // debug
  //printf("prepare_fields_elastic_adj_dev: rank %d -  anisotropic/isotropic kernels\n",mp->myrank);
  //synchronize_mpi();

  // allocates kernels
  size = NGLL3 * mp->NSPEC_AB; // note: non-aligned; if align, check memcpy below and indexing
  // density kernel
  gpuMalloc_realw((void**)&(mp->d_rho_kl),size);
  // initializes kernel values to zero
  gpuMemset_realw(mp->d_rho_kl,size,0);

  if (mp->anisotropic_kl ){
    // anisotropic kernels
    gpuMalloc_realw((void**)&(mp->d_cijkl_kl),21*size);
    gpuMemset_realw(mp->d_cijkl_kl,21*size,0);

  }else{
    // isotropic kernels
    gpuMalloc_realw((void**)&(mp->d_mu_kl),size);
    gpuMalloc_realw((void**)&(mp->d_kappa_kl),size);
    gpuMemset_realw(mp->d_mu_kl,size,0);
    gpuMemset_realw(mp->d_kappa_kl,size,0);
  }

  // strains used for attenuation and kernel simulations
  if (mp->compute_and_store_strain){
    // strains
    // debug
    //printf("prepare_fields_elastic_adj_dev: rank %d - strains\n",mp->myrank);
    //synchronize_mpi();

    size = NGLL3 * mp->NSPEC_AB; // note: non-aligned; if align, check memcpy below and indexing

    // solid pressure
    gpuCreateCopy_todevice_realw((void**)&mp->d_epsilon_trace_over_3,epsilon_trace_over_3,size);

    // kernel simulations
    if (mp->simulation_type == 3){
      // backward solid pressure
      gpuCreateCopy_todevice_realw((void**)&mp->d_b_epsilon_trace_over_3,b_epsilon_trace_over_3,size);
      // prepares backward strains
      gpuCreateCopy_todevice_realw((void**)&mp->d_b_epsilondev_xx,b_epsilondev_xx,size);
      gpuCreateCopy_todevice_realw((void**)&mp->d_b_epsilondev_yy,b_epsilondev_yy,size);
      gpuCreateCopy_todevice_realw((void**)&mp->d_b_epsilondev_xy,b_epsilondev_xy,size);
      gpuCreateCopy_todevice_realw((void**)&mp->d_b_epsilondev_xz,b_epsilondev_xz,size);
      gpuCreateCopy_todevice_realw((void**)&mp->d_b_epsilondev_yz,b_epsilondev_yz,size);
      gpuCreateCopy_todevice_realw((void**)&mp->d_b_epsilondev_trace,b_epsilondev_trace,size);
    }
  }

  // attenuation memory variables
  if (mp->attenuation){
    // debug
    //printf("prepare_fields_elastic_adj_dev: rank %d - attenuation\n",mp->myrank);
    //synchronize_mpi();

    size = *R_size;

    gpuCreateCopy_todevice_realw((void**)&mp->d_b_R_xx,b_R_xx,size);
    gpuCreateCopy_todevice_realw((void**)&mp->d_b_R_yy,b_R_yy,size);
    gpuCreateCopy_todevice_realw((void**)&mp->d_b_R_xy,b_R_xy,size);
    gpuCreateCopy_todevice_realw((void**)&mp->d_b_R_xz,b_R_xz,size);
    gpuCreateCopy_todevice_realw((void**)&mp->d_b_R_yz,b_R_yz,size);
    gpuCreateCopy_todevice_realw((void**)&mp->d_b_R_trace,b_R_trace,size);

    // alpha,beta,gamma factors for backward fields
    gpuCreateCopy_todevice_realw((void**)&mp->d_b_alphaval,b_alphaval,N_SLS);
    gpuCreateCopy_todevice_realw((void**)&mp->d_b_betaval,b_betaval,N_SLS);
    gpuCreateCopy_todevice_realw((void**)&mp->d_b_gammaval,b_gammaval,N_SLS);
  }

  // approximate hessian kernel
  mp->approximate_hess_kl = *APPROXIMATE_HESS_KL;
  if (mp->approximate_hess_kl){
    // debug
    //printf("prepare_fields_elastic_adj_dev: rank %d - hessian kernel\n",mp->myrank);
    //synchronize_mpi();

    size = NGLL3 * mp->NSPEC_AB; // note: non-aligned; if align, check memcpy below and indexing
    gpuMalloc_realw((void**)&(mp->d_hess_el_kl),size);
    gpuMemset_realw(mp->d_hess_el_kl,size,0);

    gpuMalloc_realw((void**)&(mp->d_hess_rho_el_kl),size);
    gpuMemset_realw(mp->d_hess_rho_el_kl,size,0);

    gpuMalloc_realw((void**)&(mp->d_hess_kappa_el_kl),size);
    gpuMemset_realw(mp->d_hess_kappa_el_kl,size,0);

    gpuMalloc_realw((void**)&(mp->d_hess_mu_el_kl),size);
    gpuMemset_realw(mp->d_hess_mu_el_kl,size,0);
  }

  // debug
  //printf("prepare_fields_elastic_adj_dev: rank %d - done\n",mp->myrank);
  //synchronize_mpi();

  GPU_ERROR_CHECKING("prepare_fields_elastic_adj_dev");
}

/* ----------------------------------------------------------------------------------------------- */

// purely adjoint & kernel simulations

/* ----------------------------------------------------------------------------------------------- */

extern EXTERN_LANG
void FC_FUNC_(prepare_sim2_or_3_const_device,
              PREPARE_SIM2_OR_3_CONST_DEVICE)(long* Mesh_pointer,int *nadj_rec_local, int* NTSTEP_BETWEEN_READ_ADJSRC,
                                              realw* hxir_adjstore, realw* hetar_adjstore, realw* hgammar_adjstore,
                                              int* nrec, int* h_islice_selected_rec, int* h_ispec_selected_rec) {

  TRACE("prepare_sim2_or_3_const_device");

  Mesh* mp = (Mesh*)(*Mesh_pointer);

  // for SIMULATION_TYPE == 2 or 3:
  //         nadj_rec_local - determines the number of "adjoint sources", i.e., number of station locations (STATIONS_ADJOINT),
  //                          which act as sources to drive the adjoint wavefield

  // adjoint source arrays
  mp->nadj_rec_local = *nadj_rec_local;
  if (mp->nadj_rec_local > 0){
    gpuMalloc_field((void**)&mp->d_source_adjoint,(mp->nadj_rec_local)* NDIM * (*NTSTEP_BETWEEN_READ_ADJSRC));

    // adjoint simulations
    if (mp->simulation_type == 2){
      // "pure" adjoint simulation: "adjoint sources" are located at receiver STATIONS locations
      //                            however, xir_store,.. are for "adjoint receiver" locations (specified by CMTSOLUTION)
      //                            thus, we need to store separate arrays xir_adj,.. for "adjoint sources"

      // hxir for "adjoint receivers" are at CMT locations (needed for seismograms),
      // hxir_adj for "adjoint sources" are at receiver STATIONS locations (needed to add adjoint sources)
      // sets local adjoint source positions
      gpuCreateCopy_todevice_realw((void**)&mp->d_hxir_adj,hxir_adjstore,NGLLX*mp->nadj_rec_local);
      gpuCreateCopy_todevice_realw((void**)&mp->d_hetar_adj,hetar_adjstore,NGLLX*mp->nadj_rec_local);
      gpuCreateCopy_todevice_realw((void**)&mp->d_hgammar_adj,hgammar_adjstore,NGLLX*mp->nadj_rec_local);

      // stores only local "adjoint sources" array
      int *ispec_selected_adjrec_loc;
      ispec_selected_adjrec_loc = (int*) calloc(mp->nadj_rec_local, sizeof(int));
      int iadjrec_loc = 0;
      for(int i=0; i < (*nrec); i++) {
        if (mp->myrank == h_islice_selected_rec[i]){
          ispec_selected_adjrec_loc[iadjrec_loc] = h_ispec_selected_rec[i];
          iadjrec_loc = iadjrec_loc+1;
        }
      }
      // checks
      if (iadjrec_loc != mp->nadj_rec_local) exit_on_error("prepare_sim2_or_3_const_device: nadj_rec_local not equal\n");
      // allocates on GPU
      gpuCreateCopy_todevice_int((void**)&mp->d_ispec_selected_adjrec_loc,ispec_selected_adjrec_loc,mp->nadj_rec_local);
      free(ispec_selected_adjrec_loc);

    }else{
      // kernel simulations (SIMULATION_TYPE == 3)
      // adjoint source arrays and receiver arrays are the same, no need to allocate new arrays, just point to the existing ones
      mp->d_hxir_adj = mp->d_hxir;
      mp->d_hetar_adj = mp->d_hetar;
      mp->d_hgammar_adj = mp->d_hgammar;
      // "adjoint source" locations and receiver location arrays are the same.
      mp->d_ispec_selected_adjrec_loc = mp->d_ispec_selected_rec_loc;
    }
  }

  GPU_ERROR_CHECKING("prepare_sim2_or_3_const_device");
}


/* ----------------------------------------------------------------------------------------------- */

// for NOISE simulations

/* ----------------------------------------------------------------------------------------------- */

extern EXTERN_LANG
void FC_FUNC_(prepare_fields_noise_device,
              PREPARE_FIELDS_NOISE_DEVICE)(long* Mesh_pointer,
                                           int* NSPEC_AB, int* NGLOB_AB,
                                           int* free_surface_ispec,
                                           int* free_surface_ijk,
                                           int* num_free_surface_faces,
                                           int* NOISE_TOMOGRAPHY,
                                           int* NSTEP,
                                           realw* noise_sourcearray,
                                           realw* normal_x_noise, realw* normal_y_noise, realw* normal_z_noise,
                                           realw* mask_noise,
                                           realw* free_surface_jacobian2Dw) {

  TRACE("prepare_fields_noise_device");

  Mesh* mp = (Mesh*)(*Mesh_pointer);

  // free surface
  mp->num_free_surface_faces = *num_free_surface_faces;

  gpuCreateCopy_todevice_int((void**)&mp->d_free_surface_ispec,free_surface_ispec,mp->num_free_surface_faces);
  gpuCreateCopy_todevice_int((void**)&mp->d_free_surface_ijk,free_surface_ijk,NDIM*NGLL2*mp->num_free_surface_faces);

  // alloc storage for the surface buffer to be copied
  gpuMalloc_realw((void**) &mp->d_noise_surface_movie,NDIM*NGLL2*mp->num_free_surface_faces);

  // prepares noise source array
  if (*NOISE_TOMOGRAPHY == 1){
    gpuCreateCopy_todevice_realw((void**)&mp->d_noise_sourcearray,noise_sourcearray,NDIM*NGLL3*(*NSTEP));
  }

  // prepares noise directions
  if (*NOISE_TOMOGRAPHY > 1){
    int nface_size = NGLL2*(*num_free_surface_faces);
    // allocates memory on GPU
    gpuCreateCopy_todevice_realw((void**)&mp->d_normal_x_noise,normal_x_noise,nface_size);
    gpuCreateCopy_todevice_realw((void**)&mp->d_normal_y_noise,normal_y_noise,nface_size);
    gpuCreateCopy_todevice_realw((void**)&mp->d_normal_z_noise,normal_z_noise,nface_size);
    gpuCreateCopy_todevice_realw((void**)&mp->d_mask_noise,mask_noise,nface_size);
    gpuCreateCopy_todevice_realw((void**)&mp->d_free_surface_jacobian2Dw,free_surface_jacobian2Dw,nface_size);
  }

  // prepares noise strength kernel
  if (*NOISE_TOMOGRAPHY == 3){
    gpuMalloc_realw((void**)&(mp->d_sigma_kl),NGLL3*(mp->NSPEC_AB));
    // initializes kernel values to zero
    gpuMemset_realw(mp->d_sigma_kl,NGLL3*mp->NSPEC_AB,0);
  }

  //printf("jacobian_size = %d\n",25*(*num_free_surface_faces));
  GPU_ERROR_CHECKING("prepare_fields_noise_device");
}

/* ----------------------------------------------------------------------------------------------- */

// GRAVITY simulations

/* ----------------------------------------------------------------------------------------------- */

extern EXTERN_LANG
void FC_FUNC_(prepare_fields_gravity_device,
              PREPARE_FIELDS_gravity_DEVICE)(long* Mesh_pointer,
                                             int* GRAVITY,
                                             realw* minus_deriv_gravity,
                                             realw* minus_g,
                                             realw* h_wgll_cube,
                                             int* ACOUSTIC_SIMULATION,
                                             realw* rhostore) {

  TRACE("prepare_fields_gravity_device");

  Mesh* mp = (Mesh*)(*Mesh_pointer);

  setConst_wgll_cube(h_wgll_cube,mp);

  mp->gravity = *GRAVITY;
  if (mp->gravity ){

    gpuCreateCopy_todevice_realw((void**)&mp->d_minus_deriv_gravity,minus_deriv_gravity,mp->NGLOB_AB);
    gpuCreateCopy_todevice_realw((void**)&mp->d_minus_g,minus_g,mp->NGLOB_AB);

    if (*ACOUSTIC_SIMULATION == 0){
      // density
      // rhostore not allocated yet
      int size_padded = NGLL3_PADDED * (mp->NSPEC_AB);
      // padded array
      gpuMalloc_realw((void**)&(mp->d_rhostore),size_padded);
      // transfer constant element data with padding
      /*
      // way 1: slower ...
      for(int i=0; i < mp->NSPEC_AB; i++) {
        gpuMemcpy_todevice_realw(mp->d_rhostore+i*NGLL3_PADDED, &rhostore[i*NGLL3],NGLL3);
      }
      */
      // way 2: faster...
      gpuMemcpy2D_todevice_realw(mp->d_rhostore, NGLL3_PADDED, rhostore, NGLL3, NGLL3, mp->NSPEC_AB);
    }
  }

  GPU_ERROR_CHECKING("prepare_fields_gravity_device");
}

/* ----------------------------------------------------------------------------------------------- */

// unused yet...

/*
//extern EXTERN_LANG
void FC_FUNC_(prepare_seismogram_fields,
              PREPARE_SEISMOGRAM_FIELDS)(long* Mesh_pointer,
                                         int* nrec_local,
                                         double* nu_rec,
                                         double* hxir,
                                         double* hetar,
                                         double* hgammar) {

  TRACE("prepare_constants_device");
  Mesh* mp = (Mesh*)(*Mesh_pointer);

  gpuMalloc_double((void**)&(mp->d_nu_rec),3*3*(*nrec_local));
  gpuMalloc_double((void**)&(mp->d_hxir),5*(*nrec_local));
  gpuMalloc_double((void**)&(mp->d_hetar),5*(*nrec_local));
  gpuMalloc_double((void**)&(mp->d_hgammar),5*(*nrec_local));

  gpuMalloc_realw((void**)&mp->d_seismograms_d,3*(*nrec_local));
  gpuMalloc_realw((void**)&mp->d_seismograms_v,3*(*nrec_local));
  gpuMalloc_realw((void**)&mp->d_seismograms_a,3*(*nrec_local));

  gpuMemcpy_todevice_double(mp->d_nu_rec,nu_rec,3*3*(*nrec_local));
  gpuMemcpy_todevice_double(mp->d_hxir,hxir,5*(*nrec_local));
  gpuMemcpy_todevice_double(mp->d_hetar,hetar,5*(*nrec_local));
  gpuMemcpy_todevice_double(mp->d_hgammar,hgammar,5*(*nrec_local));

  gpuMallocHost_realw((void**)&mp->h_seismograms_d_it,3**nrec_local);
  gpuMallocHost_realw((void**)&mp->h_seismograms_v_it,3**nrec_local);
  gpuMallocHost_realw((void**)&mp->h_seismograms_a_it,3**nrec_local);
}
*/





/* ----------------------------------------------------------------------------------------------- */

// FAULT simulations

/* ----------------------------------------------------------------------------------------------- */

extern EXTERN_LANG
void FC_FUNC_(prepare_fault_device,
              PREPARE_FAULT_DEVICE)(long* Mesh_pointer,
                                    int* KELVIN_VOIGT_DAMPING,
                                    realw* Kelvin_Voigt_eta) {

  TRACE("prepare_fault_device");

  Mesh* mp = (Mesh*)(*Mesh_pointer);

  // uses Kelvin-Voigt damping if Kelvin_Voigt_eta array has been allocated
  mp->use_Kelvin_Voigt_damping = *KELVIN_VOIGT_DAMPING;

  //debug
  //printf("debug: prepare_fault_device: myrank = %d , isAllocated = %d\n",mp->myrank,mp->use_Kelvin_Voigt_damping);

  // allocates & copies damping array onto GPU
  if (mp->use_Kelvin_Voigt_damping ){
    gpuCreateCopy_todevice_realw((void**)&mp->d_Kelvin_Voigt_eta,Kelvin_Voigt_eta,mp->NSPEC_AB);
  }
}


/* ----------------------------------------------------------------------------------------------- */

// cleanup

/* ----------------------------------------------------------------------------------------------- */

extern EXTERN_LANG
void FC_FUNC_(prepare_cleanup_device,
              PREPARE_CLEANUP_DEVICE)(long* Mesh_pointer,
                                      int* ACOUSTIC_SIMULATION,
                                      int* ELASTIC_SIMULATION,
                                      int* NOISE_TOMOGRAPHY) {

TRACE("prepare_cleanup_device");

  // frees allocated memory arrays
  Mesh* mp = (Mesh*)(*Mesh_pointer);

  gpuFree(mp->d_irregular_element_number);

  // frees memory on GPU
  // mesh
  gpuFree(mp->d_xix);
  gpuFree(mp->d_xiy);
  gpuFree(mp->d_xiz);
  gpuFree(mp->d_etax);
  gpuFree(mp->d_etay);
  gpuFree(mp->d_etaz);
  gpuFree(mp->d_gammax);
  gpuFree(mp->d_gammay);
  gpuFree(mp->d_gammaz);

  // absorbing boundaries
  if (mp->absorbing_conditions && mp->d_num_abs_boundary_faces > 0){
    gpuFree(mp->d_abs_boundary_ispec);
    gpuFree(mp->d_abs_boundary_ijk);
    gpuFree(mp->d_abs_boundary_normal);
    gpuFree(mp->d_abs_boundary_jacobian2Dw);
  }

  // interfaces
  if (mp->num_interfaces_ext_mesh > 0){
    gpuFree(mp->d_nibool_interfaces_ext_mesh);
    gpuFree(mp->d_ibool_interfaces_ext_mesh);
  }

  // global indexing
  gpuFree(mp->d_ispec_is_inner);
  gpuFree(mp->d_ibool);
  // element flags
  gpuFree(mp->d_ispec_is_acoustic);
  gpuFree(mp->d_ispec_is_elastic);

  // sources
  if (mp->simulation_type == 1  || mp->simulation_type == 3){
    gpuFree(mp->d_sourcearrays);
    gpuFree(mp->d_stf_pre_compute);
  }

  gpuFree(mp->d_islice_selected_source);
  gpuFree(mp->d_ispec_selected_source);

  // receivers
  if (mp->nrec_local > 0){
    gpuFree(mp->d_hxir);
    gpuFree(mp->d_hetar);
    gpuFree(mp->d_hgammar);
    if (mp->save_seismograms_d) gpuFree(mp->d_seismograms_d);
    if (mp->save_seismograms_v) gpuFree(mp->d_seismograms_v);
    if (mp->save_seismograms_a) gpuFree(mp->d_seismograms_a);
    if (mp->save_seismograms_p) gpuFree(mp->d_seismograms_p);
    gpuFree(mp->d_nu_rec);
    gpuFree(mp->d_ispec_selected_rec_loc);
    }
    gpuFree(mp->d_ispec_selected_rec);

  // ACOUSTIC arrays
  if (*ACOUSTIC_SIMULATION ){
    gpuFree(mp->d_potential_acoustic);
    gpuFree(mp->d_potential_dot_acoustic);
    gpuFree(mp->d_potential_dot_dot_acoustic);
    gpuFree(mp->d_send_potential_dot_dot_buffer);
    gpuFree(mp->d_rmass_acoustic);
    gpuFree(mp->d_rhostore);
    gpuFree(mp->d_kappastore);
    gpuFree(mp->d_phase_ispec_inner_acoustic);
    if (*NOISE_TOMOGRAPHY == 0){
      gpuFree(mp->d_free_surface_ispec);
      gpuFree(mp->d_free_surface_ijk);
    }
    if (mp->absorbing_conditions) gpuFree(mp->d_b_absorb_potential);
    if (mp->simulation_type == 3) {
      gpuFree(mp->d_b_potential_acoustic);
      gpuFree(mp->d_b_potential_dot_acoustic);
      gpuFree(mp->d_b_potential_dot_dot_acoustic);
      gpuFree(mp->d_rho_ac_kl);
      gpuFree(mp->d_kappa_ac_kl);
      if (mp->approximate_hess_kl) {
        gpuFree(mp->d_hess_ac_kl);
        gpuFree(mp->d_hess_rho_ac_kl);
        gpuFree(mp->d_hess_kappa_ac_kl);
      }
    }
  } // ACOUSTIC_SIMULATION

  // ELASTIC arrays
  if (*ELASTIC_SIMULATION ){
    gpuFree(mp->d_displ);
    gpuFree(mp->d_veloc);
    gpuFree(mp->d_accel);
    gpuFree(mp->d_send_accel_buffer);
    if (mp->simulation_type == 3) gpuFree(mp->d_b_send_accel_buffer);
    gpuFree(mp->d_rmassx);
    gpuFree(mp->d_rmassy);
    gpuFree(mp->d_rmassz);
    gpuFree(mp->d_phase_ispec_inner_elastic);
    if (mp->absorbing_conditions && mp->d_num_abs_boundary_faces > 0){
      gpuFree(mp->d_rho_vp);
      gpuFree(mp->d_rho_vs);
      if (mp->simulation_type == 3 || ( mp->simulation_type == 1 && mp->save_forward ))
          gpuFree(mp->d_b_absorb_field);
    }
    gpuFree(mp->d_kappav);
    gpuFree(mp->d_muv);
    if (mp->simulation_type == 3) {
      gpuFree(mp->d_b_displ);
      gpuFree(mp->d_b_veloc);
      gpuFree(mp->d_b_accel);
      gpuFree(mp->d_rho_kl);
      if (mp->anisotropic_kl ){
        gpuFree(mp->d_cijkl_kl);
      }else{
        gpuFree(mp->d_mu_kl);
        gpuFree(mp->d_kappa_kl);
      }
      if (mp->approximate_hess_kl) {
        gpuFree(mp->d_hess_el_kl);
        gpuFree(mp->d_hess_rho_el_kl);
        gpuFree(mp->d_hess_kappa_el_kl);
        gpuFree(mp->d_hess_mu_el_kl);
      }
    }
    if (mp->compute_and_store_strain){
      gpuFree(mp->d_epsilondev_xx);
      gpuFree(mp->d_epsilondev_yy);
      gpuFree(mp->d_epsilondev_xy);
      gpuFree(mp->d_epsilondev_xz);
      gpuFree(mp->d_epsilondev_yz);
      gpuFree(mp->d_epsilondev_trace);
      if (mp->simulation_type == 3){
        gpuFree(mp->d_epsilon_trace_over_3);
        gpuFree(mp->d_b_epsilon_trace_over_3);
        gpuFree(mp->d_b_epsilondev_xx);
        gpuFree(mp->d_b_epsilondev_yy);
        gpuFree(mp->d_b_epsilondev_xy);
        gpuFree(mp->d_b_epsilondev_xz);
        gpuFree(mp->d_b_epsilondev_yz);
        gpuFree(mp->d_b_epsilondev_trace);
      }
    }
    if (mp->attenuation){
      gpuFree(mp->d_factor_common);
      gpuFree(mp->d_factor_common_kappa);
      gpuFree(mp->d_alphaval);
      gpuFree(mp->d_betaval);
      gpuFree(mp->d_gammaval);
      gpuFree(mp->d_R_xx);
      gpuFree(mp->d_R_yy);
      gpuFree(mp->d_R_xy);
      gpuFree(mp->d_R_xz);
      gpuFree(mp->d_R_yz);
      gpuFree(mp->d_R_trace);
      if (mp->simulation_type == 3){
        gpuFree(mp->d_b_R_xx);
        gpuFree(mp->d_b_R_yy);
        gpuFree(mp->d_b_R_xy);
        gpuFree(mp->d_b_R_xz);
        gpuFree(mp->d_b_R_yz);
        gpuFree(mp->d_b_R_trace);
        gpuFree(mp->d_b_alphaval);
        gpuFree(mp->d_b_betaval);
        gpuFree(mp->d_b_gammaval);
      }
    }
    if (mp->ANISOTROPY){
      gpuFree(mp->d_c11store);
      gpuFree(mp->d_c12store);
      gpuFree(mp->d_c13store);
      gpuFree(mp->d_c14store);
      gpuFree(mp->d_c15store);
      gpuFree(mp->d_c16store);
      gpuFree(mp->d_c22store);
      gpuFree(mp->d_c23store);
      gpuFree(mp->d_c24store);
      gpuFree(mp->d_c25store);
      gpuFree(mp->d_c26store);
      gpuFree(mp->d_c33store);
      gpuFree(mp->d_c34store);
      gpuFree(mp->d_c35store);
      gpuFree(mp->d_c36store);
      gpuFree(mp->d_c44store);
      gpuFree(mp->d_c45store);
      gpuFree(mp->d_c46store);
      gpuFree(mp->d_c55store);
      gpuFree(mp->d_c56store);
      gpuFree(mp->d_c66store);
    }
    if (mp->approximate_ocean_load){
      if (mp->num_free_surface_faces > 0){
        gpuFree(mp->d_rmass_ocean_load);
        gpuFree(mp->d_free_surface_normal);
        gpuFree(mp->d_updated_dof_ocean_load);
        if (*NOISE_TOMOGRAPHY == 0){
          gpuFree(mp->d_free_surface_ispec);
          gpuFree(mp->d_free_surface_ijk);
        }
      }
    }
  } // ELASTIC_SIMULATION

  // purely adjoint & kernel array
  if (mp->simulation_type == 2 || mp->simulation_type == 3){
    if (mp->nadj_rec_local > 0){
      gpuFree(mp->d_source_adjoint);
    }
  }

  // NOISE arrays
  if (*NOISE_TOMOGRAPHY > 0){
    gpuFree(mp->d_free_surface_ispec);
    gpuFree(mp->d_free_surface_ijk);
    gpuFree(mp->d_noise_surface_movie);
    if (*NOISE_TOMOGRAPHY == 1) gpuFree(mp->d_noise_sourcearray);
    if (*NOISE_TOMOGRAPHY > 1){
      gpuFree(mp->d_normal_x_noise);
      gpuFree(mp->d_normal_y_noise);
      gpuFree(mp->d_normal_z_noise);
      gpuFree(mp->d_mask_noise);
      gpuFree(mp->d_free_surface_jacobian2Dw);
    }
    if (*NOISE_TOMOGRAPHY == 3) gpuFree(mp->d_sigma_kl);
  }

  // releases previous contexts
  gpuReset();

  // mesh pointer - not needed anymore
  free(mp);
}


