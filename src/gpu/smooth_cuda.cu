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
#include "smooth_cuda.h"


/* ----------------------------------------------------------------------------------------------- */

// smoothing routines

/* ----------------------------------------------------------------------------------------------- */

extern EXTERN_LANG
void FC_FUNC_(prepare_gpu_smooth,
              PREPARE_GPU_SMOOTH)(long * Container,
                                  realw * xstore_me,
                                  realw * ystore_me,
                                  realw * zstore_me,
                                  realw * sigma_h2_inv,
                                  realw * sigma_v2_inv,
                                  realw * h_criterion,
                                  realw * v_criterion,
                                  int * nspec_me,
                                  int * nker,
                                  realw * wgll_cube){

  TRACE("prepare_gpu_smooth");

  Smooth_data* sp = (Smooth_data*) malloc( sizeof(Smooth_data) );
  *Container = (long)sp;

  gpuCopy_todevice_realw((void**)&sp->x_me,xstore_me, NGLL3*(*nspec_me));
  gpuCopy_todevice_realw((void**)&sp->y_me,ystore_me, NGLL3*(*nspec_me));
  gpuCopy_todevice_realw((void**)&sp->z_me,zstore_me, NGLL3*(*nspec_me));
  gpuCopy_todevice_realw((void**)&sp->wgll_cube,wgll_cube, NGLL3);

  sp->sigma_h2_inv= *sigma_h2_inv;
  sp->sigma_v2_inv= *sigma_v2_inv;
  sp->h_criterion = *h_criterion;
  sp->v_criterion = *v_criterion;
  sp->nspec_me =  *nspec_me;
  sp->nker = *nker;

  gpuMalloc_realw((void**)&sp->data_smooth,NGLL3*(*nspec_me)*(*nker));
  gpuMemset_realw(sp->data_smooth,0,NGLL3*(*nspec_me)*(*nker));

  gpuMalloc_realw((void**)&sp->normalisation,NGLL3*(*nspec_me));
  gpuMemset_realw(sp->normalisation,0,NGLL3*(*nspec_me));
}

/* ----------------------------------------------------------------------------------------------- */

extern EXTERN_LANG
void FC_FUNC_(compute_smooth,
              COMPUTE_SMOOTH)(long * smooth_pointer,
                              realw * jacobian,
                              realw * jacobian_regular,
                              int * irregular_element_number,
                              realw * xstore_other,
                              realw * ystore_other,
                              realw * zstore_other,
                              realw * data_other,
                              const int * nspec_other,
                              const int * nspec_irregular_other){

  TRACE("compute_smooth");

  realw * x_other;
  realw * y_other;
  realw * z_other;
  realw * d_data_other;
  realw * d_jacobian;
  int * d_irregular_element_number;

  Smooth_data * sp = (Smooth_data*)*smooth_pointer;

  gpuCopy_todevice_realw((void**)&x_other,xstore_other,NGLL3*(*nspec_other));
  gpuCopy_todevice_realw((void**)&y_other,ystore_other,NGLL3*(*nspec_other));
  gpuCopy_todevice_realw((void**)&z_other,zstore_other,NGLL3*(*nspec_other));

  if (*nspec_irregular_other > 0){
    gpuCopy_todevice_realw((void**)&d_jacobian,jacobian,NGLL3*(*nspec_irregular_other));
  } else {
    gpuCopy_todevice_realw((void**)&d_jacobian,jacobian,1);
  }
  gpuCopy_todevice_int((void**)&d_irregular_element_number,irregular_element_number,(*nspec_other));

  dim3 grid(sp->nspec_me,1);
  dim3 threads(NGLL3,1,1);

  for (int i=0;i<sp->nker;i++){
    gpuCopy_todevice_realw((void**)&d_data_other,&data_other[NGLL3*(*nspec_other)*i],NGLL3*(*nspec_other));

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
                                       i,
                                       sp->nspec_me,
                                       *nspec_other,
                                       sp->v_criterion,
                                       sp->h_criterion,
                                       d_jacobian,
                                       d_irregular_element_number,
                                       *jacobian_regular,
                                       sp->data_smooth,
                                       sp->normalisation,
                                       sp->wgll_cube);
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
                                         i,
                                         sp->nspec_me,
                                         *nspec_other,
                                         sp->v_criterion,
                                         sp->h_criterion,
                                         d_jacobian,
                                         d_irregular_element_number,
                                         *jacobian_regular,
                                         sp->data_smooth,
                                         sp->normalisation,
                                         sp->wgll_cube);
    }
#endif

    gpuFree(d_data_other);
  }
  gpuSynchronize();

  gpuFree(x_other);
  gpuFree(y_other);
  gpuFree(z_other);
  gpuFree(d_jacobian);
  gpuFree(d_irregular_element_number);
}

/* ----------------------------------------------------------------------------------------------- */

extern EXTERN_LANG
void FC_FUNC_(get_smooth,
              GET_SMOOTH)(long * smooth_pointer,
                          realw * data_smooth) {

  TRACE("get_smooth");

  Smooth_data * sp = (Smooth_data*)*smooth_pointer;

  dim3 grid(sp->nspec_me,1);
  dim3 threads(NGLL3,1,1);

#ifdef USE_CUDA
  if (run_cuda){
    normalize_data<<<grid,threads>>>(sp->data_smooth,sp->normalisation,sp->nker,sp->nspec_me);
  }
#endif
#ifdef USE_HIP
  if (run_hip){
    hipLaunchKernelGGL(normalize_data, dim3(grid), dim3(threads), 0, 0,
                                       sp->data_smooth,sp->normalisation,sp->nker,sp->nspec_me);
  }
#endif

  gpuMemcpy_tohost_realw(data_smooth, sp->data_smooth,NGLL3*sp->nspec_me*sp->nker);

  gpuFree(sp->x_me);
  gpuFree(sp->y_me);
  gpuFree(sp->z_me);
  gpuFree(sp->data_smooth);
  gpuFree(sp->wgll_cube);
  gpuFree(sp->normalisation);

  free(sp);
}

