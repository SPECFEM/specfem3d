/*
 !=====================================================================
 !
 !               S p e c f e m 3 D  V e r s i o n  3 . 0
 !               ---------------------------------------
 !
 !     Main historical authors: Dimitri Komatitsch and Jeroen Tromp
 !                        Princeton University, USA
 !                and CNRS / University of Marseille, France
 !                 (there are currently many more authors!)
 ! (c) Princeton University and CNRS / University of Marseille, July 2012
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

/* trivia

- for most working arrays we use now "realw" instead of "float" type declarations to make it easier to switch
  between a real or double precision simulation
  (matching CUSTOM_REAL == 4 or 8 in fortran routines).

- instead of boolean "logical" declared in fortran routines, in C (or Cuda-C) we have to use "int" variables.
  ifort / gfortran caveat:
    to check whether it is true or false, do not check for == 1 to test for true values since ifort just uses
    non-zero values for true (e.g. can be -1 for true). however, false will be always == 0.
  thus, rather use: if (var ) {...}  for testing if true instead of if (var == 1){...} (alternative: one could use if (var != 0){...}

*/

#ifndef MESH_CONSTANTS_CUDA_H
#define MESH_CONSTANTS_CUDA_H

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <sys/types.h>
#include <sys/time.h>
#include <sys/resource.h>

#include "config.h"

#ifdef WITH_MPI
#include <mpi.h>
#endif

#include <cuda.h>
#include <cuda_runtime.h>


/* ----------------------------------------------------------------------------------------------- */

// for debugging and benchmarking

/* ----------------------------------------------------------------------------------------------- */

#define DEBUG 0
#if DEBUG == 1
#define TRACE(x) printf("%s\n",x);
#else
#define TRACE(x) // printf("%s\n",x);
#endif

#define MAXDEBUG 0
#if MAXDEBUG == 1
#define LOG(x) printf("%s\n",x)
#define PRINT5(var,offset) for(;print_count<5;print_count++) printf("var(%d)=%2.20f\n",print_count,var[offset+print_count]);
#define PRINT10(var) if (print_count<10) { printf("var=%1.20e\n",var); print_count++; }
#define PRINT10i(var) if (print_count<10) { printf("var=%d\n",var); print_count++; }
#else
#define LOG(x) // printf("%s\n",x);
#define PRINT5(var,offset) // for(i=0;i<10;i++) printf("var(%d)=%f\n",i,var[offset+i]);
#endif

// performance timers
#define CUDA_TIMING 0
#define CUDA_TIMING_UPDATE 0

// error checking after cuda function calls
// (note: this synchronizes many calls, thus e.g. no asynchronuous memcpy possible)
//#define ENABLE_VERY_SLOW_ERROR_CHECKING

// maximum function
#define MAX(x,y)     (((x) < (y)) ? (y) : (x))
// minimum function
#define MIN(a,b)     (((a) > (b)) ? (b) : (a))

/* ----------------------------------------------------------------------------------------------- */

// cuda constant arrays

/* ----------------------------------------------------------------------------------------------- */

// dimensions
#define NDIM 3

// Gauss-Lobatto-Legendre
#define NGLLX 5
#define NGLLY NGLLX
#define NGLLZ NGLLX
#define NGLL2 25
#define NGLL3 125  // no padding: requires same size as in fortran for NGLLX * NGLLY * NGLLZ

// padding: 128 == 2**7 might improve on older graphics cards w/ coalescent memory accesses:
#define NGLL3_PADDED 128
// no padding: 125 == 5*5*5 to avoid allocation of extra memory
//#define NGLL3_PADDED 125

// number of standard linear solids
#define N_SLS 3

/* ----------------------------------------------------------------------------------------------- */

// Output paths, see setup/constants.h
#define OUTPUT_FILES "./OUTPUT_FILES/"

/* ----------------------------------------------------------------------------------------------- */

// (optional) pre-processing directive used in kernels: if defined check that it is also set in setup/constants.h:
// leads up to ~ 5% performance increase
//#define USE_MESH_COLORING_GPU

/* ----------------------------------------------------------------------------------------------- */

// Texture memory usage:
// requires CUDA version >= 4.0, see check below
// Use textures for d_displ and d_accel -- 10% performance boost
//#define USE_TEXTURES_FIELDS

// Using texture memory for the hprime-style constants is slower on
// Fermi generation hardware, but *may* be faster on Kepler
// generation.
// Use textures for hprime_xx
//#define USE_TEXTURES_CONSTANTS

// CUDA version >= 4.0 needed for cudaTextureType1D and cudaDeviceSynchronize()
#if CUDA_VERSION < 4000
#undef USE_TEXTURES_FIELDS
#undef USE_TEXTURES_CONSTANTS
#endif

#ifdef USE_TEXTURES_FIELDS
#pragma message ("\nCompiling with: USE_TEXTURES_FIELDS enabled\n")
#endif
#ifdef USE_TEXTURES_CONSTANTS
#pragma message ("\nCompiling with: USE_TEXTURES_CONSTANTS enabled\n")
#endif

// (optional) unrolling loops
// leads up to ~1% performance increase
//#define MANUALLY_UNROLLED_LOOPS

// compiler specifications
// (optional) use launch_bounds specification to increase compiler optimization
// (depending on GPU type, register spilling might slow down the performance)
// (uncomment if desired)
//#define USE_LAUNCH_BOUNDS

// elastic kernel
// note: main kernel is Kernel_2_***_impl() which is limited by shared memory usage to 8 active blocks
//       while register usage might use up to 9 blocks
//
// performance statistics: kernel Kernel_2_noatt_impl():
//       shared memory per block = 1700    for Kepler: total = 49152 -> limits active blocks to 16
//       registers per thread    = 48
//       registers per block     = 6144                total = 65536 -> limits active blocks to 10
//
// performance statistics: kernel Kernel_2_att_impl():
//       shared memory per block = 6100    for Kepler: total = 49152 -> limits active blocks to 8
//       registers per thread    = 59
//       registers per block     = 8192                total = 65536 -> limits active blocks to 8
#define LAUNCH_MIN_BLOCKS 10

// acoustic kernel
// performance statistics: kernel Kernel_2_acoustic_impl():
//       shared memory per block = 2200    for Kepler: -> limits active blocks to 16 (maximum possible)
//       registers per thread    = 40
//       registers per block     = 5120                -> limits active blocks to 12
// note: for K20x, using a minimum of 16 blocks leads to register spilling.
//       this slows down the kernel by ~ 4%
#define LAUNCH_MIN_BLOCKS_ACOUSTIC 16

/* ----------------------------------------------------------------------------------------------- */

// cuda kernel block size for updating displacements/potential (newmark time scheme)
// current hardware: 128 is slightly faster than 256 ( ~ 4%)
#define BLOCKSIZE_KERNEL1 128
#define BLOCKSIZE_KERNEL3 128
#define BLOCKSIZE_TRANSFER 256

// maximum grid dimension in one direction of GPU
#define MAXIMUM_GRID_DIM 65535

/* ----------------------------------------------------------------------------------------------- */

// indexing
#define INDEX2(xsize,x,y) x + (y)*xsize
#define INDEX3(xsize,ysize,x,y,z) x + xsize*(y + ysize*z)
#define INDEX4(xsize,ysize,zsize,x,y,z,i) x + xsize*(y + ysize*(z + zsize*i))
#define INDEX5(xsize,ysize,zsize,isize,x,y,z,i,j) x + xsize*(y + ysize*(z + zsize*(i + isize*(j))))
#define INDEX6(xsize,ysize,zsize,isize,jsize,x,y,z,i,j,k) x + xsize*(y + ysize*(z + zsize*(i + isize*(j + jsize*k))))

#define INDEX4_PADDED(xsize,ysize,zsize,x,y,z,i) x + xsize*(y + ysize*z) + (i)*NGLL3_PADDED

/* ----------------------------------------------------------------------------------------------- */

// custom type declarations

/* ----------------------------------------------------------------------------------------------- */

// type of "working" variables: see also CUSTOM_REAL
// double precision temporary variables leads to 10% performance decrease
// in Kernel_2_impl (not very much..)
typedef float realw;
// textures
typedef texture<float, cudaTextureType1D, cudaReadModeElementType> realw_texture;

// pointer declarations
// restricted pointers: can improve performance on Kepler ~ 10%
//   see: http://docs.nvidia.com/cuda/cuda-c-programming-guide/index.html#restrict
//   however, compiler tends to use texture loads for restricted memory arrays, which might slow down performance
//
// non-restricted (default)
//typedef const realw* realw_const_p;
// restricted
typedef const realw* __restrict__ realw_const_p;
//
// non-restricted (default)
//typedef realw* realw_p;
// restricted
typedef realw* __restrict__ realw_p;

// wrapper for global memory load function
// usage:  val = get_global_cr( &A[index] );
#if __CUDA_ARCH__ >= 350
// Device has ldg
__device__ __forceinline__ realw get_global_cr(realw_const_p ptr) { return __ldg(ptr); }
#else
//Device does not, fall back.
__device__ __forceinline__ realw get_global_cr(realw_const_p ptr) { return (*ptr); }
#endif


/* ----------------------------------------------------------------------------------------------- */

// Specific definitions and overloads for option NB_RUNS_ON_ACOUSTIC_GPU

/* ----------------------------------------------------------------------------------------------- */

// option to run several events within the same MPI slice
#define NB_RUNS_ACOUSTIC_GPU 1

// Functions and operators are overloaded according to the number of runs to perform
#if NB_RUNS_ACOUSTIC_GPU == 1
typedef realw field;
inline __host__ __device__ field Make_field(realw* b){ return b[0];}
inline __host__ __device__ field Make_field(realw b){ return b;}
inline __host__ __device__ field max(field b){ return b;}
inline __host__ __device__ realw sum(field b){ return b;}
#endif

#if NB_RUNS_ACOUSTIC_GPU == 2
typedef float2 field;
inline __host__ __device__ field Make_field(realw* b){ return make_float2(b[0],b[1]);}
inline __host__ __device__ field Make_field(realw b){ return make_float2(b,b);}
inline __host__ __device__ realw fabs(field b){ return max(fabs(b.x),fabs(b.y));}
inline __host__ __device__ void operator+=(field &a, field b){a.x += b.x;a.y += b.y;}
inline __host__ __device__ void operator-=(field &a, field b){a.x -= b.x;a.y -= b.y;}
inline __host__ __device__ field operator*(field a, realw b){ return make_float2(a.x * b, a.y * b);}
inline __host__ __device__ field operator*(realw b, field a){ return make_float2(a.x * b, a.y * b);}
inline __host__ __device__ field operator*(field a, field b){ return make_float2(a.x * b.x, a.y * b.y);}
inline __host__ __device__ field operator+(field a, field b){ return make_float2(a.x + b.x, a.y + b.y);}
inline __host__ __device__ field operator/(field a, realw b){ return make_float2(a.x / b, a.y / b);}
inline __host__ __device__ field operator-(field a, field b){ return make_float2(a.x - b.x, a.y - b.y);}
inline __host__ __device__ field operator-(field a){ return make_float2(-a.x,-a.y);}
inline __host__ __device__ realw sum(field b){ return b.x+b.y;}
inline __device__ void atomicAdd(field* address, field val){atomicAdd(&(address->x),val.x); atomicAdd(&(address->y),val.y);}

//dummy overloads, just to enable compilation
inline __host__ __device__ field operator+(field a, realw b){ return a;}
inline __device__ void atomicAdd(field* address, float val){}
inline __device__ void atomicAdd(realw* address, field val){}
inline __host__ __device__ void operator+=(realw &a, field b){}

//texture are not compatible for the moment with floatN data
#undef USE_TEXTURES_FIELDS
#endif

#if NB_RUNS_ACOUSTIC_GPU == 4
typedef float4 field;
inline __host__ __device__ field Make_field(realw* b){ return make_float4(b[0],b[1],b[2],b[3]);}
inline __host__ __device__ field Make_field(realw b){ return make_float4(b,b,b,b);}
inline __host__ __device__ realw fabs(field b){ return max(max(fabs(b.x),fabs(b.y)),max(fabs(b.z),fabs(b.w)));}
inline __host__ __device__ void operator+=(field &a, field b){a.x += b.x;a.y += b.y;a.z += b.z;a.w += b.w;}
inline __host__ __device__ void operator-=(field &a, field b){a.x -= b.x;a.y -= b.y;a.z -= b.z;a.w -= b.w;}
inline __host__ __device__ field operator*(field a, realw b){ return make_float4(a.x * b, a.y * b,a.z * b, a.w * b);}
inline __host__ __device__ field operator*(realw b, field a){ return make_float4(a.x * b, a.y * b,a.z * b, a.w * b);}
inline __host__ __device__ field operator*(field a, field b){ return make_float4(a.x * b.x, a.y * b.y,a.z * b.z, a.w * b.w);}
inline __host__ __device__ field operator+(field a, field b){ return make_float4(a.x + b.x, a.y + b.y,a.z + b.z, a.w + b.w);}
inline __host__ __device__ field operator/(field a, realw b){ return make_float4(a.x / b, a.y / b,a.z / b, a.w / b);}
inline __host__ __device__ field operator-(field a, field b){ return make_float4(a.x - b.x, a.y - b.y,a.z - b.z, a.w - b.w);}
inline __host__ __device__ field operator-(field a){ return make_float4(-a.x,-a.y,-a.z,-a.w);}
inline __host__ __device__ realw sum(field b){ return b.x+b.y+b.z+b.w;}
inline __device__ void atomicAdd(field* address, field val){atomicAdd(&(address->x),val.x); atomicAdd(&(address->y),val.y);atomicAdd(&(address->z),val.z); atomicAdd(&(address->w),val.w);}

inline __host__ __device__ int operator>(field &a, realw b){return (a.x>b || a.y>b || a.z>b || a.w >b ) ? 1:0;}


//dummy overloads, just to enable compilation
inline __host__ __device__ field operator+(field a, realw b){ return a;}
inline __device__ void atomicAdd(field* address, float val){}
inline __device__ void atomicAdd(realw* address, field val){}
inline __host__ __device__ void operator+=(realw &a, field b){}

//texture are not compatible for the moment with floatN data
#undef USE_TEXTURES_FIELDS
#endif
typedef const field* __restrict__ field_const_p;
typedef field* __restrict__ field_p;


/* ----------------------------------------------------------------------------------------------- */

// utility functions: defined in check_fields_cuda.cu

/* ----------------------------------------------------------------------------------------------- */

double get_time();
void get_free_memory(double* free_db, double* used_db, double* total_db);
void print_CUDA_error_if_any(cudaError_t err, int num);
void pause_for_debugger(int pause);
void exit_on_cuda_error(const char* kernel_name);
void exit_on_error(const char* info);
void synchronize_cuda();
void synchronize_mpi();
void start_timing_cuda(cudaEvent_t* start,cudaEvent_t* stop);
void stop_timing_cuda(cudaEvent_t* start,cudaEvent_t* stop, const char* info_str);
void stop_timing_cuda(cudaEvent_t* start,cudaEvent_t* stop, const char* info_str,realw* t);
void get_blocks_xy(int num_blocks,int* num_blocks_x,int* num_blocks_y);
realw get_device_array_maximum_value(realw* array,int size);


/* ----------------------------------------------------------------------------------------------- */

// mesh pointer wrapper structure

/* ----------------------------------------------------------------------------------------------- */

typedef struct mesh_ {

  // mesh resolution
  int NSPEC_AB;
  int NGLOB_AB;

  // mpi process
  int myrank;

  // constants
  int simulation_type;
  int save_forward;
  int use_mesh_coloring_gpu;
  int absorbing_conditions;
  int gravity;


  // ------------------------------------------------------------------ //
  // GLL points & weights
  // ------------------------------------------------------------------ //

  int* d_irregular_element_number;
  realw xix_regular,jacobian_regular;

  // interpolators
  realw* d_xix; realw* d_xiy; realw* d_xiz;
  realw* d_etax; realw* d_etay; realw* d_etaz;
  realw* d_gammax; realw* d_gammay; realw* d_gammaz;

  // model parameters
  realw* d_kappav; realw* d_muv;

  // global indexing
  int* d_ibool;

  // inner / outer elements
  int* d_ispec_is_inner;

  // pointers to constant memory arrays
  realw* d_hprime_xx;
  //realw* d_hprime_yy; // only needed if NGLLX != NGLLY != NGLLZ
  //realw* d_hprime_zz; // only needed if NGLLX != NGLLY != NGLLZ

  realw* d_hprimewgll_xx;
  //realw* d_hprimewgll_yy; // only needed if NGLLX != NGLLY != NGLLZ
  //realw* d_hprimewgll_zz; // only needed if NGLLX != NGLLY != NGLLZ

  realw* d_wgllwgll_xy; realw* d_wgllwgll_xz; realw* d_wgllwgll_yz;
  realw* d_wgll_cube;

  // A buffer for mpi-send/recv, which is duplicated in fortran but is
  // allocated with pinned memory to facilitate asynchronus device <->
  // host memory transfers
  float* h_send_accel_buffer;
  float* h_send_b_accel_buffer;

  float* send_buffer;
  float* h_recv_accel_buffer;
  float* h_recv_b_accel_buffer;
  float* recv_buffer;

  int size_mpi_buffer;
  int size_mpi_buffer_potential;

  // mpi interfaces
  int num_interfaces_ext_mesh;
  int max_nibool_interfaces_ext_mesh;
  int* d_nibool_interfaces_ext_mesh;
  int* d_ibool_interfaces_ext_mesh;

  // overlapped memcpy streams
  cudaStream_t compute_stream;
  cudaStream_t copy_stream;
  //cudaStream_t b_copy_stream;

  // sources
  int nsources_local;
  realw* d_sourcearrays;
  field* d_stf_pre_compute;
  int* d_islice_selected_source;
  int* d_ispec_selected_source;

  // receivers
  int* d_ispec_selected_rec_loc;
  int* d_ispec_selected_rec;
  int nrec_local;

  realw* d_hxir, *d_hetar, *d_hgammar;
  realw* d_seismograms_d, *d_seismograms_v, *d_seismograms_a, * d_nu;
  field* d_seismograms_p;

  int save_seismograms_d;
  int save_seismograms_v;
  int save_seismograms_a;
  int save_seismograms_p;

  //realw* h_seismograms_d_it;
  //realw* h_seismograms_v_it;
  //realw* h_seismograms_a_it;

  // adjoint receivers/sources
  int nadj_rec_local;
  field* d_source_adjoint;

  // ------------------------------------------------------------------ //
  // elastic wavefield parameters
  // ------------------------------------------------------------------ //

  // displacement, velocity, acceleration
  realw* d_displ; realw* d_veloc; realw* d_accel;
  // backward/reconstructed elastic wavefield
  realw* d_b_displ; realw* d_b_veloc; realw* d_b_accel;

  // elastic elements
  int* d_ispec_is_elastic;

  // elastic domain parameters
  int* d_phase_ispec_inner_elastic;
  int num_phase_ispec_elastic;

  // mesh coloring
  int* h_num_elem_colors_elastic;
  int num_colors_outer_elastic,num_colors_inner_elastic;
  int nspec_elastic;

  realw* d_rmassx;
  realw* d_rmassy;
  realw* d_rmassz;

  // mpi buffer
  realw* d_send_accel_buffer;
  realw* d_b_send_accel_buffer;

  //used for absorbing stacey boundaries
  int d_num_abs_boundary_faces;
  int* d_abs_boundary_ispec;
  int* d_abs_boundary_ijk;
  realw* d_abs_boundary_normal;
  realw* d_abs_boundary_jacobian2Dw;

  realw* d_b_absorb_field;
  int d_b_reclen_field;

  realw* d_rho_vp;
  realw* d_rho_vs;

  // surface elements (to save for noise tomography and acoustic simulations)
  int* d_free_surface_ispec;
  int* d_free_surface_ijk;
  int num_free_surface_faces;

  // surface movie elements to save for noise tomography
  realw* d_noise_surface_movie;

  // attenuation
  realw* d_R_xx;
  realw* d_R_yy;
  realw* d_R_xy;
  realw* d_R_xz;
  realw* d_R_yz;

  realw* d_factor_common;

  realw* d_alphaval;
  realw* d_betaval;
  realw* d_gammaval;

  // attenuation & kernel
  realw* d_epsilondev_xx;
  realw* d_epsilondev_yy;
  realw* d_epsilondev_xy;
  realw* d_epsilondev_xz;
  realw* d_epsilondev_yz;
  realw* d_epsilon_trace_over_3;

  // anisotropy
  realw* d_c11store;
  realw* d_c12store;
  realw* d_c13store;
  realw* d_c14store;
  realw* d_c15store;
  realw* d_c16store;
  realw* d_c22store;
  realw* d_c23store;
  realw* d_c24store;
  realw* d_c25store;
  realw* d_c26store;
  realw* d_c33store;
  realw* d_c34store;
  realw* d_c35store;
  realw* d_c36store;
  realw* d_c44store;
  realw* d_c45store;
  realw* d_c46store;
  realw* d_c55store;
  realw* d_c56store;
  realw* d_c66store;

  // noise
  realw* d_normal_x_noise;
  realw* d_normal_y_noise;
  realw* d_normal_z_noise;
  realw* d_mask_noise;
  realw* d_free_surface_jacobian2Dw;

  realw* d_noise_sourcearray;

  // attenuation & kernel backward fields
  realw* d_b_R_xx;
  realw* d_b_R_yy;
  realw* d_b_R_xy;
  realw* d_b_R_xz;
  realw* d_b_R_yz;

  realw* d_b_epsilondev_xx;
  realw* d_b_epsilondev_yy;
  realw* d_b_epsilondev_xy;
  realw* d_b_epsilondev_xz;
  realw* d_b_epsilondev_yz;
  realw* d_b_epsilon_trace_over_3;

  realw* d_b_alphaval;
  realw* d_b_betaval;
  realw* d_b_gammaval;

  // sensitivity kernels
  int anisotropic_kl;
  realw* d_rho_kl;
  realw* d_mu_kl;
  realw* d_kappa_kl;
  realw* d_cijkl_kl;

  // noise sensitivity kernel
  realw* d_sigma_kl;

  // approximative hessian for preconditioning kernels
  realw* d_hess_el_kl, *d_hess_rho_el_kl, *d_hess_mu_el_kl, *d_hess_kappa_el_kl;

  // oceans
  realw* d_rmass_ocean_load;
  realw* d_free_surface_normal;
  int* d_updated_dof_ocean_load;

  // JC JC here we will need to add GPU support for the new C-PML routines

  // ------------------------------------------------------------------ //
  // acoustic wavefield
  // ------------------------------------------------------------------ //
  // potential and first and second time derivative
  field* d_potential_acoustic, *d_potential_dot_acoustic, *d_potential_dot_dot_acoustic;
  // backward/reconstructed wavefield
  field* d_b_potential_acoustic, *d_b_potential_dot_acoustic, *d_b_potential_dot_dot_acoustic;

  // acoustic domain parameters
  int* d_ispec_is_acoustic;

  int* d_phase_ispec_inner_acoustic;
  int num_phase_ispec_acoustic;

  // mesh coloring
  int* h_num_elem_colors_acoustic;
  int num_colors_outer_acoustic,num_colors_inner_acoustic;
  int nspec_acoustic;

  realw* d_rhostore;
  realw* d_kappastore;
  realw* d_rmass_acoustic;

  // mpi buffer
  field* d_send_potential_dot_dot_buffer;
  field* d_b_send_potential_dot_dot_buffer;

  field* d_b_absorb_potential;
  int d_b_reclen_potential;

  // sensitivity kernels
  realw* d_rho_ac_kl;
  realw* d_kappa_ac_kl;

  // approximative hessian for preconditioning kernels
  realw* d_hess_ac_kl, *d_hess_rho_ac_kl, *d_hess_kappa_ac_kl;

  // coupling acoustic-elastic
  int* d_coupling_ac_el_ispec;
  int* d_coupling_ac_el_ijk;
  realw* d_coupling_ac_el_normal;
  realw* d_coupling_ac_el_jacobian2Dw;

  // gravity
  realw* d_minus_deriv_gravity;
  realw* d_minus_g;
  // FAULT
  int Kelvin_Voigt_damping;
  realw* d_Kelvin_Voigt_eta;

  // for option NB_RUNS_FOR_ACOUSTIC_GPU
  int* run_number_of_the_source;

} Mesh;


#endif
