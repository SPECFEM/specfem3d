#include "mesh_constants_gpu.h"
#include "smooth_pde_cuda.h"

extern EXTERN_LANG
void FC_FUNC_(prepare_smooth_pde_gpu,
              PREPARE_SMOOTH_PDE_GPU)(long * Mesh_pointer,
                                      long * Container_smooth_pde,
                                      field* dat_smooth_glob,
                                      realw* rvol,
                                      int* phase_ispec_inner,
                                      int* num_phase_ispec,
                                      int* CPML_to_spec,
                                      int* NSPEC_CPML,
                                      realw* cv,
                                      realw* ch){
  TRACE("prepare_smooth_pde_gpu");
  Mesh* mp = (Mesh*)(*Mesh_pointer);
  Smooth_pde_data* sp = (Smooth_pde_data*) malloc(sizeof(Smooth_pde_data));
  *Container_smooth_pde = (long)sp;

  gpuCreateCopy_todevice_realw((void**)&sp->d_dat_smooth_glob, dat_smooth_glob, mp->NGLOB_AB);
  gpuMalloc_field((void**)&sp->d_ddat_smooth_glob, mp->NGLOB_AB);

  // initializes values to zero
  gpuMemset_realw(sp->d_ddat_smooth_glob,mp->NGLOB_AB,0);

  sp->size_mpi_buffer_smooth = (mp->num_interfaces_ext_mesh) * (mp->max_nibool_interfaces_ext_mesh);
  gpuMalloc_field((void**)&sp->d_send_buffer, sp->size_mpi_buffer_smooth);
  gpuCreateCopy_todevice_realw((void**)&sp->d_rvol, rvol, mp->NGLOB_AB);

  sp->num_phase_ispec = *num_phase_ispec;
  gpuCreateCopy_todevice_int((void**)&sp->d_phase_ispec_inner, phase_ispec_inner, sp->num_phase_ispec*2);

  sp->NSPEC_CPML = *NSPEC_CPML;
  if (sp->NSPEC_CPML > 0) {
    gpuCreateCopy_todevice_int((void**)&sp->d_CPML_to_spec, CPML_to_spec, sp->NSPEC_CPML);
  }
  sp->cv = *cv;
  sp->ch = *ch;

  GPU_ERROR_CHECKING("prepare_smooth_pde_gpu");
}

/* ----------------------------------------------------------------------------------------------- */

extern EXTERN_LANG
void FC_FUNC_(compute_update_element_smooth_pde_gpu,
              COMPUTE_UPDATE_ELEMENT_SMOOTH_PDE_GPU)(long * Mesh_pointer,
                                                     long * Container_smooth_pde,
                                                     int * iphase,
                                                     int * num_elements){
  TRACE("compute_update_element_smooth_pde_gpu");
  Mesh* mp = (Mesh*)(*Mesh_pointer);
  Smooth_pde_data* sp = (Smooth_pde_data*)(*Container_smooth_pde);

  int blocksize = NGLL3_PADDED;

  int num_blocks_x, num_blocks_y;
  get_blocks_xy(*num_elements,&num_blocks_x,&num_blocks_y);

  dim3 grid(num_blocks_x,num_blocks_y);
  dim3 threads(blocksize,1,1);

#ifdef USE_CUDA
  if (run_cuda) {
    Kernel_2_smooth_pde<<<grid,threads,0,mp->compute_stream>>>(*num_elements,
                                                               mp->d_ibool,
                                                               mp->d_irregular_element_number,
                                                               sp->d_phase_ispec_inner,
                                                               sp->num_phase_ispec,
                                                               *iphase,
                                                               sp->d_dat_smooth_glob,
                                                               sp->d_ddat_smooth_glob,
                                                               mp->d_xix, mp->d_xiy, mp->d_xiz,
                                                               mp->d_etax,mp->d_etay,mp->d_etaz,
                                                               mp->d_gammax,mp->d_gammay,mp->d_gammaz,
                                                               mp->xix_regular,mp->jacobian_regular,
                                                               mp->d_hprime_xx,
                                                               mp->d_hprimewgll_xx,
                                                               mp->d_wgllwgll_xy,
                                                               mp->d_wgllwgll_xz,
                                                               mp->d_wgllwgll_yz,
                                                               sp->cv, sp->ch
                                                               );
  }
#endif
#ifdef USE_HIP
  if (run_hip) {
    hipLaunchKernelGGL(HIP_KERNEL_NAME(Kernel_2_smooth_pde), dim3(grid), dim3(threads), 0, mp->compute_stream,
                                                               *num_elements,
                                                               mp->d_ibool,
                                                               mp->d_irregular_element_number,
                                                               sp->d_phase_ispec_inner,
                                                               sp->num_phase_ispec,
                                                               *iphase,
                                                               sp->d_dat_smooth_glob,
                                                               sp->d_ddat_smooth_glob,
                                                               mp->d_xix, mp->d_xiy, mp->d_xiz,
                                                               mp->d_etax,mp->d_etay,mp->d_etaz,
                                                               mp->d_gammax,mp->d_gammay,mp->d_gammaz,
                                                               mp->xix_regular,mp->jacobian_regular,
                                                               mp->d_hprime_xx,
                                                               mp->d_hprimewgll_xx,
                                                               mp->d_wgllwgll_xy,
                                                               mp->d_wgllwgll_xz,
                                                               mp->d_wgllwgll_yz,
                                                               sp->cv, sp->ch
                                                               );
  }
#endif

  GPU_ERROR_CHECKING("compute_update_element_smooth_pde_gpu");
}

/* ----------------------------------------------------------------------------------------------- */

extern EXTERN_LANG
void FC_FUNC_(transfer_boun_dat_smooth_pde_from_device,
              TRANSFER_BOUN_DAT_SMOOTH_PDE_FROM_DEVICE)(long * Mesh_pointer,
                                                        long * Container_smooth_pde,
                                                        field * send_buffer_smooth){
  TRACE("transfer_boun_dat_smooth_pde_from_device");
  Mesh* mp = (Mesh*)(*Mesh_pointer);
  Smooth_pde_data* sp = (Smooth_pde_data*)(*Container_smooth_pde);

  if (sp->size_mpi_buffer_smooth > 0) {
    int blocksize = BLOCKSIZE_TRANSFER;
    int size_padded = ((int)ceil(((double)(mp->max_nibool_interfaces_ext_mesh))/((double)blocksize)))*blocksize;

    int num_blocks_x, num_blocks_y;
    get_blocks_xy(size_padded/blocksize,&num_blocks_x,&num_blocks_y);

    dim3 grid(num_blocks_x,num_blocks_y);
    dim3 threads(blocksize,1,1);
#ifdef USE_CUDA
    if (run_cuda) {
      prepare_boundary_potential_on_device<<<grid,threads,0,mp->compute_stream>>>(sp->d_ddat_smooth_glob,
                                                                                  sp->d_send_buffer,
                                                                                  mp->num_interfaces_ext_mesh,
                                                                                  mp->max_nibool_interfaces_ext_mesh,
                                                                                  mp->d_nibool_interfaces_ext_mesh,
                                                                                  mp->d_ibool_interfaces_ext_mesh
                                                                                  );
    }
#endif
#ifdef USE_HIP
    if (run_hip) {
       hipLaunchKernelGGL(prepare_boundary_potential_on_device, dim3(grid), dim3(threads), 0, mp->compute_stream,
                                                                                  sp->d_ddat_smooth_glob,
                                                                                  sp->d_send_buffer,
                                                                                  mp->num_interfaces_ext_mesh,
                                                                                  mp->max_nibool_interfaces_ext_mesh,
                                                                                  mp->d_nibool_interfaces_ext_mesh,
                                                                                  mp->d_ibool_interfaces_ext_mesh
                                                                                  );
    }
#endif
    gpuStreamSynchronize(mp->compute_stream);

    // copies buffer to CPU
    gpuMemcpy_tohost_field(send_buffer_smooth,sp->d_send_buffer,sp->size_mpi_buffer_smooth);
  }

  GPU_ERROR_CHECKING("transfer_boun_dat_smooth_pde_from_device");
}

/* ----------------------------------------------------------------------------------------------- */

extern EXTERN_LANG
void FC_FUNC_(transfer_asmbl_dat_smooth_pde_from_device,
              TRANSFER_ASMBL_DAT_SMOOTH_PDE_FROM_DEVICE)(long * Mesh_pointer,
                                                         long * Container_smooth_pde,
                                                         field * recv_buffer_smooth){
  TRACE("transfer_asmbl_dat_smooth_pde_from_device");

  Mesh* mp = (Mesh*)(*Mesh_pointer);
  Smooth_pde_data* sp = (Smooth_pde_data*)(*Container_smooth_pde);

  if (sp->size_mpi_buffer_smooth > 0) {
    int blocksize = BLOCKSIZE_TRANSFER;
    int size_padded = ((int)ceil(((double)(mp->max_nibool_interfaces_ext_mesh))/((double)blocksize)))*blocksize;

    int num_blocks_x, num_blocks_y;
    get_blocks_xy(size_padded/blocksize,&num_blocks_x,&num_blocks_y);

    dim3 grid(num_blocks_x,num_blocks_y);
    dim3 threads(blocksize,1,1);

    // synchronizes
    gpuSynchronize();

    // copies buffer onto GPU
    gpuMemcpy_todevice_field(sp->d_send_buffer,recv_buffer_smooth,sp->size_mpi_buffer_smooth);

#ifdef USE_CUDA
    if (run_cuda) {
      assemble_boundary_potential_on_device<<<grid,threads,0,mp->compute_stream>>>(sp->d_ddat_smooth_glob,
                                                                                  sp->d_send_buffer,
                                                                                  mp->num_interfaces_ext_mesh,
                                                                                  mp->max_nibool_interfaces_ext_mesh,
                                                                                  mp->d_nibool_interfaces_ext_mesh,
                                                                                  mp->d_ibool_interfaces_ext_mesh
                                                                                  );
    }
#endif
#ifdef USE_HIP
    if (run_hip) {
       hipLaunchKernelGGL(assemble_boundary_potential_on_device, dim3(grid), dim3(threads), 0, mp->compute_stream,
                                                                                  sp->d_ddat_smooth_glob,
                                                                                  sp->d_send_buffer,
                                                                                  mp->num_interfaces_ext_mesh,
                                                                                  mp->max_nibool_interfaces_ext_mesh,
                                                                                  mp->d_nibool_interfaces_ext_mesh,
                                                                                  mp->d_ibool_interfaces_ext_mesh
                                                                                  );
    }
#endif
  }

  GPU_ERROR_CHECKING("transfer_asmbl_dat_smooth_pde_from_device");
}

/* ----------------------------------------------------------------------------------------------- */

extern EXTERN_LANG
void FC_FUNC_(kernel_3_smooth_pde_cuda,
              KERNEL_3_SMOOTH_PDE_CUDA)(long * Mesh_pointer,
                                        long * Container_smooth_pde) {
  TRACE("kernel_3_smooth_pde_cuda");

  Mesh* mp = (Mesh*)(*Mesh_pointer);
  Smooth_pde_data* sp = (Smooth_pde_data*)(*Container_smooth_pde);

  int size = mp->NGLOB_AB;

  int blocksize = BLOCKSIZE_KERNEL3;
  int size_padded = ((int)ceil(((double)size)/((double)blocksize)))*blocksize;

  int num_blocks_x, num_blocks_y;
  get_blocks_xy(size_padded/blocksize,&num_blocks_x,&num_blocks_y);

  dim3 grid(num_blocks_x,num_blocks_y);
  dim3 threads(blocksize,1,1);

#ifdef USE_CUDA
  if (run_cuda) {
    kernel_3_smooth_pde_cuda_device<<<grid, threads, 0, mp->compute_stream>>>(sp->d_ddat_smooth_glob,
                                                                              sp->d_rvol,
                                                                              size);
  }
#endif
#ifdef USE_HIP
  if (run_hip) {
    hipLaunchKernelGGL(kernel_3_smooth_pde_cuda_device, dim3(grid), dim3(threads), 0, mp->compute_stream,
                                                                              sp->d_ddat_smooth_glob,
                                                                              sp->d_rvol,
                                                                              size);
  }
#endif

  GPU_ERROR_CHECKING("kernel_3_smooth_pde_cuda");
}

/* ----------------------------------------------------------------------------------------------- */

extern EXTERN_LANG
void FC_FUNC_(update_dat_smooth_pde_cuda,
              UPDATE_DAT_SMOOTH_PDE_CUDA)(long * Mesh_pointer,
                                          long * Container_smooth_pde) {
  TRACE("update_dat_smooth_pde_cuda");

  Mesh* mp = (Mesh*)(*Mesh_pointer);
  Smooth_pde_data* sp = (Smooth_pde_data*)(*Container_smooth_pde);

  int size = mp->NGLOB_AB;

  int blocksize = BLOCKSIZE_KERNEL1;
  int size_padded = ((int)ceil(((double)size)/((double)blocksize)))*blocksize;

  int num_blocks_x, num_blocks_y;
  get_blocks_xy(size_padded/blocksize,&num_blocks_x,&num_blocks_y);

  dim3 grid(num_blocks_x,num_blocks_y);
  dim3 threads(blocksize,1,1);

#ifdef USE_CUDA
  if (run_cuda) {
    UpdateData_smooth_pde_kernel<<<grid, threads, 0, mp->compute_stream>>>(sp->d_dat_smooth_glob,
                                                                           sp->d_ddat_smooth_glob,
                                                                           size
                                                                           );
  }
#endif
#ifdef USE_HIP
  if (run_hip) {
    hipLaunchKernelGGL(UpdateData_smooth_pde_kernel, dim3(grid), dim3(threads), 0, mp->compute_stream,
                                                                           sp->d_dat_smooth_glob,
                                                                           sp->d_ddat_smooth_glob,
                                                                           size
                                                                           );
  }
#endif

  GPU_ERROR_CHECKING("update_dat_smooth_pde_cuda");
}

/* ----------------------------------------------------------------------------------------------- */

extern EXTERN_LANG
void FC_FUNC_(zero_pml_smooth_pde_cuda,
              ZERO_PML_SMOOTH_PDE_CUDA)(long * Mesh_pointer,
                                        long * Container_smooth_pde) {
  TRACE("zero_pml_smooth_pde_cuda");

  Mesh* mp = (Mesh*)(*Mesh_pointer);
  Smooth_pde_data* sp = (Smooth_pde_data*)(*Container_smooth_pde);

  int num_elements = sp->NSPEC_CPML;

  int blocksize = NGLL3_PADDED;

  int num_blocks_x, num_blocks_y;
  get_blocks_xy(num_elements,&num_blocks_x,&num_blocks_y);

  dim3 grid(num_blocks_x,num_blocks_y);
  dim3 threads(blocksize,1,1);

#ifdef USE_CUDA
  if (run_cuda) {
    zero_pml_smooth_pde_kernel<<<grid, threads, 0, mp->compute_stream>>>(num_elements,
                                                                         sp->d_dat_smooth_glob,
                                                                         mp->d_ibool,
                                                                         sp->d_CPML_to_spec);
  }
#endif
#ifdef USE_HIP
  if (run_hip) {
    hipLaunchKernelGGL(zero_pml_smooth_pde_kernel, dim3(grid), dim3(threads), 0, mp->compute_stream,
                                                                         num_elements,
                                                                         sp->d_dat_smooth_glob,
                                                                         mp->d_ibool,
                                                                         sp->d_CPML_to_spec);
  }
#endif

  GPU_ERROR_CHECKING("zero_pml_smooth_pde_cuda");
}

/* ----------------------------------------------------------------------------------------------- */

extern EXTERN_LANG
void FC_FUNC_(get_norm_smooth_pde_from_device,
              GET_NORM_SMOOTH_PDE_FROM_DEVICE)(long * Mesh_pointer,
                                               long * Container_smooth_pde,
                                               realw * norm, int * ind_val){
  TRACE("get_norm_smooth_pde_from_device");

  Mesh* mp = (Mesh*)(*Mesh_pointer);
  Smooth_pde_data* sp = (Smooth_pde_data*)(*Container_smooth_pde);

  realw max = 0.0;
  realw *d_max;

  //initializes
  *norm = 0.0f;

  // way 2 b: timing Elapsed time: 1.236916e-03
  // launch simple reduction kernel
  realw* h_max;
  int blocksize = BLOCKSIZE_TRANSFER;

  int size = mp->NGLOB_AB;
  int size_padded = ((int)ceil(((double)size)/((double)blocksize)))*blocksize;

  int num_blocks_x, num_blocks_y;
  get_blocks_xy(size_padded/blocksize,&num_blocks_x,&num_blocks_y);

  dim3 grid(num_blocks_x,num_blocks_y);
  dim3 threads(blocksize,1,1);

  // on host (allocates & initializes to zero)
  h_max = (realw*) calloc(num_blocks_x*num_blocks_y,sizeof(realw));

  // allocates memory on device
  gpuMalloc_realw((void**)&d_max,num_blocks_x*num_blocks_y);
  // initializes values to zero
  gpuMemset_realw(d_max,num_blocks_x*num_blocks_y,0);

  gpuStreamSynchronize(mp->compute_stream);
#ifdef USE_CUDA
  if (run_cuda){
    if ((*ind_val) == 0) {
      get_maximum_field_kernel<<<grid,threads,0,mp->compute_stream>>>(sp->d_dat_smooth_glob,size,d_max);
    }
    else if ((*ind_val) == 1) {
      get_maximum_field_kernel<<<grid,threads,0,mp->compute_stream>>>(sp->d_ddat_smooth_glob,size,d_max);
    }
  }
#endif
#ifdef USE_HIP
  if (run_hip){
    if ((*ind_val) == 0) {
      hipLaunchKernelGGL(get_maximum_field_kernel, dim3(grid), dim3(threads), 0, mp->compute_stream,
                                                                    sp->d_dat_smooth_glob,size,d_max);
    }
    else if ((*ind_val) == 1) {
      hipLaunchKernelGGL(get_maximum_field_kernel, dim3(grid), dim3(threads), 0, mp->compute_stream,
                                                                    sp->d_ddat_smooth_glob,size,d_max);
    }
  }
#endif

  // synchronizes
  //gpuSynchronize();
  // explicitly waits for stream to finish
  // (cudaMemcpy implicitly synchronizes all other cuda operations)
  gpuStreamSynchronize(mp->compute_stream);

  gpuMemcpy_tohost_realw(h_max,d_max,num_blocks_x*num_blocks_y);

  // determines max for all blocks
  max = h_max[0];
  for(int i=1;i<num_blocks_x*num_blocks_y;i++) {
    if (max < h_max[i]) max = h_max[i];
  }

  gpuFree(d_max);
  free(h_max);

  // return result
  *norm = max;

  GPU_ERROR_CHECKING("get_norm_smooth_pde_from_device");
}

/* ----------------------------------------------------------------------------------------------- */

extern EXTERN_LANG
void FC_FUNC_(transfer_dat_smooth_pde_from_device,
              TRANSFER_DAT_SMOOTH_PDE_FROM_DEVICE)(long * Mesh_pointer,
                                                   long * Container_smooth_pde,
                                                   realw* dat_smooth_glob){

  TRACE("transfer_dat_smooth_pde_from_device");

  Mesh* mp = (Mesh*)(*Mesh_pointer);
  Smooth_pde_data* sp = (Smooth_pde_data*)(*Container_smooth_pde);

  gpuMemcpy_tohost_realw(dat_smooth_glob, sp->d_dat_smooth_glob, mp->NGLOB_AB);

  GPU_ERROR_CHECKING("transfer_dat_smooth_pde_from_device");
}
