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

#include "smooth_cuda.h"
#include "config.h"
#include <stdio.h>
// copies integer array from CPU host to GPU device
void copy_todevice_int(void** d_array_addr_ptr,int* h_array,int size){
   cudaMalloc((void**)d_array_addr_ptr,size*sizeof(int));
   cudaMemcpy((int*) *d_array_addr_ptr,h_array,size*sizeof(int),cudaMemcpyHostToDevice);
}

void copy_todevice_realw(void** d_array_addr_ptr,realw* h_array,int size){
   cudaMalloc((void**)d_array_addr_ptr,size*sizeof(realw));
   cudaMemcpy((realw*) *d_array_addr_ptr,h_array,size*sizeof(realw),cudaMemcpyHostToDevice);
}

__global__ void process_smooth(realw_const_p xstore_me,realw_const_p ystore_me,realw_const_p zstore_me,realw_const_p xstore_other,realw_const_p ystore_other,realw_const_p zstore_other, realw_const_p data_other, const realw sigma_h2_inv, const realw sigma_v2_inv, const int iker, const int nspec_me, const int nspec_other, const realw v_criterion, const realw h_criterion, realw_const_p jacobian, const int * irregular_element_number,realw jacobian_regular,realw_p sum_data_smooth, realw_p normalisation,realw_const_p wgll_cube){

int ispec = blockIdx.x + gridDim.x*blockIdx.y;
int ispec_irreg;
int igll = threadIdx.x;
int gll_other;
realw x_me,y_me, z_me, x_other,y_other, z_other, coef, normalisation_slice;
realw dat;
__shared__ int sh_test[NGLL3];
__shared__ realw sh_x_other[NGLL3];
__shared__ realw sh_y_other[NGLL3];
__shared__ realw sh_z_other[NGLL3];
__shared__ realw sh_jacobian[NGLL3];
__shared__ realw sh_wgll_cube[NGLL3];
__shared__ realw sh_data[NGLL3];

int n_loop = nspec_other/NGLL3 + 1;
x_me = xstore_me[NGLL3*ispec + igll ];
y_me = ystore_me[NGLL3*ispec + igll ];
z_me = zstore_me[NGLL3*ispec + igll ];
sh_wgll_cube[igll]=wgll_cube[igll];
__syncthreads();

dat=0.f;
normalisation_slice=0.f;
//We test 125 spectral elements at a time
for (int i=0;i<n_loop;i++)
{
__syncthreads();
if (NGLL3*i + threadIdx.x < nspec_other){
x_other = (xstore_other[i*NGLL3*NGLL3 + threadIdx.x*NGLL3 ] + xstore_other[i*NGLL3*NGLL3 + threadIdx.x*NGLL3 + NGLL3 - 1 ])/2;
y_other = (ystore_other[i*NGLL3*NGLL3 + threadIdx.x*NGLL3 ] + ystore_other[i*NGLL3*NGLL3 + threadIdx.x*NGLL3 + NGLL3 - 1 ])/2;
z_other = (zstore_other[i*NGLL3*NGLL3 + threadIdx.x*NGLL3 ] + zstore_other[i*NGLL3*NGLL3 + threadIdx.x*NGLL3 + NGLL3 - 1 ])/2;

}

sh_test[threadIdx.x] = ( NGLL3*i + threadIdx.x >= nspec_other || ((x_me-x_other)*(x_me-x_other)+ (y_me-y_other)*(y_me-y_other)) > h_criterion || (z_me-z_other)*(z_me-z_other) > v_criterion ) ? 1 : 0 ;
__syncthreads();

//loop over each spectral element tested
for (int k=0;k<NGLL3;k++)
{
__syncthreads();
if (sh_test[k]) continue ;
//Load data from other slice to shared memory
sh_x_other[igll] = xstore_other[i*NGLL3*NGLL3 + k*NGLL3 + igll ];
sh_y_other[igll] = ystore_other[i*NGLL3*NGLL3 + k*NGLL3 + igll ];
sh_z_other[igll] = zstore_other[i*NGLL3*NGLL3 + k*NGLL3 + igll ];

sh_data[igll] = data_other[i*NGLL3*NGLL3 + k*NGLL3 + igll ];
ispec_irreg = irregular_element_number[i*NGLL3 + k]-1;
if (ispec_irreg >= 0){sh_jacobian[igll] = jacobian[ispec_irreg*NGLL3 + igll ];}
else{sh_jacobian[igll] = jacobian_regular;}
__syncthreads();

for (int j=0;j<NGLL3;j++){

gll_other = (igll + j) % NGLL3;

x_other = sh_x_other[gll_other];
y_other = sh_y_other[gll_other];
z_other = sh_z_other[gll_other];
coef = expf(- sigma_h2_inv*((x_me-x_other)*(x_me-x_other) + (y_me-y_other)*(y_me-y_other))- sigma_v2_inv*(z_me-z_other)*(z_me-z_other))*sh_jacobian[gll_other]*sh_wgll_cube[gll_other];
normalisation_slice = normalisation_slice + coef;
dat += sh_data[gll_other]*coef;
} //loop on each gll_other
} //loop on each spec_other tested
} //loop on each serie of 125 spec_other

sum_data_smooth[NGLL3*nspec_me*iker+NGLL3*ispec + igll] += dat;
normalisation[NGLL3*ispec + igll] += normalisation_slice;
}

__global__ void normalize_data(realw_p data_smooth, realw_const_p normalisation,int nker, int nspec_me){
int ispec = blockIdx.x + gridDim.x*blockIdx.y;
realw norm = normalisation[NGLL3*ispec + threadIdx.x];
for (int j=0;j<nker;j++) data_smooth[NGLL3*nspec_me*j + NGLL3*ispec + threadIdx.x] /= norm/nker;
}

extern "C"
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

  Smooth_data* sp = (Smooth_data*)malloc( sizeof(Smooth_data) );
  *Container = (long)sp;

  copy_todevice_realw((void**)&sp->x_me,xstore_me, NGLL3*(*nspec_me));
  copy_todevice_realw((void**)&sp->y_me,ystore_me, NGLL3*(*nspec_me));
  copy_todevice_realw((void**)&sp->z_me,zstore_me, NGLL3*(*nspec_me));
  copy_todevice_realw((void**)&sp->wgll_cube,wgll_cube, NGLL3);

  sp->sigma_h2_inv= *sigma_h2_inv;
  sp->sigma_v2_inv= *sigma_v2_inv;
  sp->h_criterion = *h_criterion;
  sp->v_criterion = *v_criterion;
  sp->nspec_me =  *nspec_me;
  sp->nker = *nker;

  print_CUDA_error_if_any(cudaMalloc((void**)&sp->data_smooth,NGLL3*(*nspec_me)*(*nker)*sizeof(realw)),2000);
  print_CUDA_error_if_any(cudaMemset(sp->data_smooth,0,NGLL3*(*nspec_me)*(*nker)*sizeof(realw)),2001);

  print_CUDA_error_if_any(cudaMalloc((void**)&sp->normalisation,NGLL3*(*nspec_me)*sizeof(realw)),2002);
  print_CUDA_error_if_any(cudaMemset(sp->normalisation,0,NGLL3*(*nspec_me)*sizeof(realw)),2003);
}

extern "C"
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
realw * x_other;
realw * y_other;
realw * z_other;
realw * d_data_other;
realw * d_jacobian;
int * d_irregular_element_number;

Smooth_data * sp = (Smooth_data*)*smooth_pointer;

copy_todevice_realw((void**)&x_other,xstore_other,NGLL3*(*nspec_other));
copy_todevice_realw((void**)&y_other,ystore_other,NGLL3*(*nspec_other));
copy_todevice_realw((void**)&z_other,zstore_other,NGLL3*(*nspec_other));
if (*nspec_irregular_other > 0){copy_todevice_realw((void**)&d_jacobian,jacobian,NGLL3*(*nspec_irregular_other));}
else {copy_todevice_realw((void**)&d_jacobian,jacobian,1);}
copy_todevice_int((void**)&d_irregular_element_number,irregular_element_number,(*nspec_other));

dim3 grid(sp->nspec_me,1);
dim3 threads(NGLL3,1,1);

for (int i=0;i<sp->nker;i++)
{
copy_todevice_realw((void**)&d_data_other,&data_other[NGLL3*(*nspec_other)*i],NGLL3*(*nspec_other));

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
cudaFree(d_data_other);
}
synchronize_cuda();
cudaFree(x_other);
cudaFree(y_other);
cudaFree(z_other);
cudaFree(d_jacobian);
cudaFree(d_irregular_element_number);
}

extern "C"
void FC_FUNC_(get_smooth,
              GET_SMOOTH)(long * smooth_pointer,realw * data_smooth){

Smooth_data * sp = (Smooth_data*)*smooth_pointer;
dim3 grid(sp->nspec_me,1);
dim3 threads(NGLL3,1,1);

normalize_data<<<grid,threads>>>(sp->data_smooth,sp->normalisation,sp->nker,sp->nspec_me);

print_CUDA_error_if_any(cudaMemcpy(data_smooth, sp->data_smooth,
                                     NGLL3*sp->nspec_me*sizeof(int)*sp->nker, cudaMemcpyDeviceToHost),98010);

cudaFree(sp->x_me);
cudaFree(sp->y_me);
cudaFree(sp->z_me);
cudaFree(sp->data_smooth);
cudaFree(sp->wgll_cube);
cudaFree(sp->normalisation);
free(sp);
}

