/*
!=====================================================================
!
!                          S p e c f e m 3 D
!                          -----------------
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
#include "smooth_cuda.h"


/* ----------------------------------------------------------------------------------------------- */

// smoothing routines

/* ----------------------------------------------------------------------------------------------- */

extern EXTERN_LANG
void FC_FUNC_(prepare_smooth_gpu,
              PREPARE_SMOOTH_GPU)(long * Container,
                                  realw * xstore_me,
                                  realw * ystore_me,
                                  realw * zstore_me,
                                  realw * sigma_h2,
                                  realw * sigma_v2,
                                  realw * sigma_h3,
                                  realw * sigma_v3,
                                  int * nspec_me,
                                  int * nker){

  TRACE("prepare_smooth_gpu");

  // allocates structure
  Smooth_data* sp = (Smooth_data*) malloc( sizeof(Smooth_data) );

  // checks pointer
  if (! sp) exit_on_error ("Error allocating smooth data pointer");

  // sets fortran pointer
  *Container = (long)sp;

  gpuCreateCopy_todevice_realw((void**)&sp->x_me,xstore_me, NGLL3*(*nspec_me));
  gpuCreateCopy_todevice_realw((void**)&sp->y_me,ystore_me, NGLL3*(*nspec_me));
  gpuCreateCopy_todevice_realw((void**)&sp->z_me,zstore_me, NGLL3*(*nspec_me));

  // sets variable values
  realw sigma_h2_inv = 1.0f / (*sigma_h2);
  realw sigma_v2_inv = 1.0f / (*sigma_v2);

  realw sigma_h3_sq = (*sigma_h3) * (*sigma_h3);  // squared
  realw sigma_v3_sq = (*sigma_v3) * (*sigma_v3);

  sp->sigma_h2_inv = sigma_h2_inv;
  sp->sigma_v2_inv = sigma_v2_inv;

  sp->h_criterion = sigma_h3_sq;
  sp->v_criterion = sigma_v3_sq;

  sp->nspec_me =  *nspec_me;
  sp->nker = *nker;

  // smoothed kernel
  gpuMalloc_realw((void**)&sp->data_smooth,NGLL3*(*nspec_me)*(*nker));
  gpuMemset_realw(sp->data_smooth,NGLL3*(*nspec_me)*(*nker),0);

  // normalization factor
  gpuMalloc_realw((void**)&sp->normalisation,NGLL3*(*nspec_me));
  gpuMemset_realw(sp->normalisation,NGLL3*(*nspec_me),0);

  GPU_ERROR_CHECKING ("prepare_smooth_gpu");
}

/* ----------------------------------------------------------------------------------------------- */

extern EXTERN_LANG
void FC_FUNC_(compute_smooth_gpu,
              COMPUTE_SMOOTH_GPU)(long * smooth_pointer,
                                  realw * data_other,
                                  realw * integ_factor,
                                  realw * xstore_other,
                                  realw * ystore_other,
                                  realw * zstore_other,
                                  const int * nspec_other_f){

  TRACE("compute_smooth_gpu");

  Smooth_data * sp = (Smooth_data*) *smooth_pointer;

  int nspec_other = *nspec_other_f;

  realw * x_other;
  realw * y_other;
  realw * z_other;

  realw * d_data_other;
  realw * d_integ_factor;

  gpuCreateCopy_todevice_realw((void**)&x_other,xstore_other,NGLL3 * nspec_other);
  gpuCreateCopy_todevice_realw((void**)&y_other,ystore_other,NGLL3 * nspec_other);
  gpuCreateCopy_todevice_realw((void**)&z_other,zstore_other,NGLL3 * nspec_other);

  gpuCreateCopy_todevice_realw((void**)&d_integ_factor,integ_factor,NGLL3 * nspec_other);

  //dim3 grid(sp->nspec_me,1);
  //dim3 threads(NGLL3,1,1);

  int num_blocks_x, num_blocks_y;
  get_blocks_xy (sp->nspec_me, &num_blocks_x, &num_blocks_y);

  dim3 grid(num_blocks_x,num_blocks_y);
  dim3 threads(NGLL3,1,1);

  for (int iker=0; iker < sp->nker; iker++){
    // pointer to corresponding kernel section
    realw * data_p = &data_other[NGLL3 * nspec_other * iker];

    // copy single kernel
    gpuCreateCopy_todevice_realw((void**)&d_data_other,data_p,NGLL3 * nspec_other);

#ifdef USE_CUDA
    if (run_cuda){
      process_smooth<<<grid,threads>>>(sp->x_me,
                                       sp->y_me,
                                       sp->z_me,
                                       x_other,
                                       y_other,
                                       z_other,
                                       d_data_other,
                                       sp->sigma_h2_inv,
                                       sp->sigma_v2_inv,
                                       iker,
                                       sp->nspec_me,
                                       nspec_other,
                                       sp->v_criterion,
                                       sp->h_criterion,
                                       d_integ_factor,
                                       sp->data_smooth,
                                       sp->normalisation);
    }
#endif
#ifdef USE_HIP
    if (run_hip){
      hipLaunchKernelGGL(process_smooth, dim3(grid), dim3(threads), 0, 0,
                                         sp->x_me,
                                         sp->y_me,
                                         sp->z_me,
                                         x_other,
                                         y_other,
                                         z_other,
                                         d_data_other,
                                         sp->sigma_h2_inv,
                                         sp->sigma_v2_inv,
                                         iker,
                                         sp->nspec_me,
                                         nspec_other,
                                         sp->v_criterion,
                                         sp->h_criterion,
                                         d_integ_factor,
                                         sp->data_smooth,
                                         sp->normalisation);
    }
#endif

    gpuFree(d_data_other);
  }
  gpuSynchronize();

  gpuFree(x_other);
  gpuFree(y_other);
  gpuFree(z_other);
  gpuFree(d_integ_factor);

  GPU_ERROR_CHECKING ("compute_smooth_gpu");
}

/* ----------------------------------------------------------------------------------------------- */

extern EXTERN_LANG
void FC_FUNC_(get_smooth_gpu,
              GET_SMOOTH_gpu)(long * smooth_pointer,
                              realw * data_smooth) {

  TRACE("get_smooth_gpu");

  Smooth_data * sp = (Smooth_data*) *smooth_pointer;

  //dim3 grid(sp->nspec_me,1);
  //dim3 threads(NGLL3,1,1);

  int num_blocks_x, num_blocks_y;
  get_blocks_xy (sp->nspec_me, &num_blocks_x, &num_blocks_y);

  dim3 grid(num_blocks_x,num_blocks_y);
  dim3 threads(NGLL3,1,1);

#ifdef USE_CUDA
  if (run_cuda){
    normalize_data<<<grid,threads>>>(sp->data_smooth,
                                     sp->normalisation,
                                     sp->nker,
                                     sp->nspec_me);
  }
#endif
#ifdef USE_HIP
  if (run_hip){
    hipLaunchKernelGGL(normalize_data, dim3(grid), dim3(threads), 0, 0,
                                       sp->data_smooth,
                                       sp->normalisation,
                                       sp->nker,
                                       sp->nspec_me);
  }
#endif

  gpuMemcpy_tohost_realw(data_smooth, sp->data_smooth, NGLL3 * sp->nspec_me * sp->nker);

  gpuFree(sp->x_me);
  gpuFree(sp->y_me);
  gpuFree(sp->z_me);
  gpuFree(sp->data_smooth);
  gpuFree(sp->normalisation);

  free(sp);

  GPU_ERROR_CHECKING ("get_smooth_gpu");
}

