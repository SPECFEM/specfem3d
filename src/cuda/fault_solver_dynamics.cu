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

#include "fault_struct_cuda.h"
#include "mesh_constants_cuda.h"

// asserts
#include <assert.h>

/* ----------------------------------------------------------------------------------------------- */

extern "C"
void FC_FUNC_(initialize_fault_solver,
              INITIALIZE_FAULT_SOLVER)(long* Fault_solver,
                                       int* num_of_faults,
                                       realw* v_healing,
                                       realw* v_rupt) {

  Fault_solver_dynamics *Fdyn = (Fault_solver_dynamics*)malloc(sizeof(Fault_solver_dynamics));
  Fdyn->faults = (Fault*)malloc((*num_of_faults)*sizeof(Fault));
  Fdyn->v_healing = *v_healing;
  Fdyn->v_rupt = *v_rupt;
  Fdyn->Nbfaults = *num_of_faults;
  *Fault_solver = (long) Fdyn;
}

extern "C"
void FC_FUNC_(initialize_fault_data,
              INITIALIZE_FAULT_DATA)(long* Fault_solver,
                                     int* iglob,
                                     int* num_of_records,
                                     int* nt,
                                     int* ifault) {

  Fault_data *fault_data_recorder = (Fault_data*)malloc(sizeof(Fault_data));
  fault_data_recorder->NRECORD = *num_of_records;
  fault_data_recorder->NT = *nt;

  int recordlength = 7; //store 7 different quantities

  print_CUDA_error_if_any(cudaMalloc((void**)&(fault_data_recorder->iglob), fault_data_recorder->NRECORD*sizeof(int)),60001);

  print_CUDA_error_if_any(cudaMemcpy(fault_data_recorder->iglob, iglob, fault_data_recorder->NRECORD*sizeof(int), cudaMemcpyHostToDevice),60002);

  print_CUDA_error_if_any(cudaMalloc((void**)&(fault_data_recorder->dataT), recordlength*fault_data_recorder->NRECORD*fault_data_recorder->NT*sizeof(realw)),60003);

  Fault_solver_dynamics* Fsolver = (Fault_solver_dynamics*)(*Fault_solver);

  Fsolver->faults[*ifault].output_dataT = fault_data_recorder;
}


/* ----------------------------------------------------------------------------------------------- */

// copies integer array from CPU host to GPU device
void copy_todevice_realw_test(void** d_array_addr_ptr,realw* h_array,int size) {

  // allocates memory on GPU
  cudaMalloc((void**)d_array_addr_ptr,size*sizeof(realw));
  // copies values onto GPU
  cudaMemcpy((realw*) *d_array_addr_ptr,h_array,size*sizeof(realw),cudaMemcpyHostToDevice);
}

/* ----------------------------------------------------------------------------------------------- */

void copy_tohost_realw_test(void** d_array_addr_ptr,realw* h_array,int size) {

  // copies values onto GPU
  cudaMemcpy(h_array, (realw*) *d_array_addr_ptr,size*sizeof(realw),cudaMemcpyDeviceToHost);
}

/* ----------------------------------------------------------------------------------------------- */

// copies integer array from CPU host to GPU device
void copy_todevice_int_test(void** d_array_addr_ptr,int* h_array,int size) {

  // allocates memory on GPU
  cudaMalloc((void**)d_array_addr_ptr,size*sizeof(int));
  // copies values onto GPU
  cudaMemcpy((realw*) *d_array_addr_ptr,h_array,size*sizeof(int),cudaMemcpyHostToDevice);
}

/* ----------------------------------------------------------------------------------------------- */

void copy_tohost_int_test(void** d_array_addr_ptr,int* h_array,int size) {

  // allocates memory on GPU
  cudaMemcpy(h_array,(realw*) *d_array_addr_ptr,size*sizeof(int),cudaMemcpyDeviceToHost);
}

/* ----------------------------------------------------------------------------------------------- */

void allocate_cuda_memory_test(void** d_array_addr_ptr,int size) {

  // allocates memory on GPU
  cudaMalloc((void**)d_array_addr_ptr,size*sizeof(int));
}

/* ----------------------------------------------------------------------------------------------- */

extern "C"
void FC_FUNC_(transfer_todevice_fault_data,
              TRANSFER_TODEVICE_FAULT_DATA)(long* Fault_pointer,
                                            int* fault_index,
                                            int* NSPEC_AB,
                                            int* NGLOB_AB,
                                            realw* D,
                                            realw* T0,
                                            realw* T,
                                            realw* B,
                                            realw* R,
                                            realw* V0,
                                            realw* Z,
                                            realw* invM1,
                                            realw* invM2,
                                            int* ibulk1,
                                            int* ibulk2) {

  Fault_solver_dynamics* Fsolver = (Fault_solver_dynamics*)(*Fault_pointer);
  Fault* Flt = &(Fsolver->faults[*fault_index]);
  Flt->NSPEC_AB = *NSPEC_AB;
  Flt->NGLOB_AB = *NGLOB_AB;
  if (*NGLOB_AB > 0){
    copy_todevice_realw_test((void **)&(Flt->B),B,*NGLOB_AB);
    copy_todevice_realw_test((void **)&(Flt->R),R,(*NGLOB_AB)*9);
    copy_todevice_realw_test((void **)&(Flt->Z),Z,(*NGLOB_AB));
    copy_todevice_realw_test((void **)&(Flt->T0),T0,(*NGLOB_AB)*3);
    copy_todevice_int_test((void **)&(Flt->ibulk1),ibulk1,(*NGLOB_AB));
    copy_todevice_int_test((void **)&(Flt->ibulk2),ibulk2,(*NGLOB_AB));
    copy_todevice_realw_test((void **)&(Flt->V),V0,(*NGLOB_AB)*3);
    copy_todevice_realw_test((void **)&(Flt->D),D,(*NGLOB_AB)*3);
    copy_todevice_realw_test((void **)&(Flt->invM1),invM1,*NGLOB_AB);
    copy_todevice_realw_test((void **)&(Flt->invM2),invM2,*NGLOB_AB);
    copy_todevice_realw_test((void **)&(Flt->T),T,(*NGLOB_AB)*3);
  }

}

/* ----------------------------------------------------------------------------------------------- */

extern "C"
void FC_FUNC_(transfer_tohost_fault_data,
              TRANSFER_TOHOST_FAULT_DATA)(long* Fault_pointer,
                                          int* fault_index,
                                          int* NSPEC_AB,
                                          int* NGLOB_AB,
                                          realw* D,
                                          realw* V,
                                          realw* T) {

  Fault_solver_dynamics* Fsolver = (Fault_solver_dynamics*)(*Fault_pointer);
  Fault* Flt = &(Fsolver->faults[*fault_index]);
  if (*NGLOB_AB > 0){
    copy_tohost_realw_test((void **)&(Flt->V),V,(*NGLOB_AB)*3);
    copy_tohost_realw_test((void **)&(Flt->D),D,(*NGLOB_AB)*3);
    copy_tohost_realw_test((void **)&(Flt->T),T,(*NGLOB_AB)*3);
  }
}

/* ----------------------------------------------------------------------------------------------- */

extern "C"
void FC_FUNC_(transfer_tohost_datat,
              TRANSFER_TOHOST_DATAT)(long* Fault_pointer,
                                     realw* h_dataT,
                                     int* it,
                                     int* ifault) {

  Fault_solver_dynamics* Fsolver = (Fault_solver_dynamics*)(*Fault_pointer);
  Fault* fault_pointer = Fsolver->faults + *ifault;
  if (fault_pointer->output_dataT->NRECORD > 0){
    copy_tohost_realw_test((void **)&(fault_pointer->output_dataT->dataT),h_dataT + (*it -fault_pointer->output_dataT->NT) * fault_pointer->output_dataT->NRECORD*7
      ,fault_pointer->output_dataT->NRECORD*7*fault_pointer->output_dataT->NT);
  }
}

/*-------------------------------------------------------------------------------------------------*/

extern "C"
void FC_FUNC_(transfer_todevice_rsf_data,
              TRANSFER_TODEVICE_RSF_DATA)(long* Fault_pointer,
                                          int *NGLOB_AB,
                                          int* fault_index,
                                          realw* V0,
                                          realw* f0,
                                          realw* V_init,
                                          realw* a,
                                          realw* b,
                                          realw* L,
                                          realw* theta,
                                          realw* T,
                                          realw* C,
                                          realw* fw,
                                          realw* Vw) {

  Fault_solver_dynamics* Fsolver = (Fault_solver_dynamics*)(*Fault_pointer);
  Rsf_type* Rsf  = &((Fsolver->faults[*fault_index]).rsf);
  if (*NGLOB_AB > 0){
    copy_todevice_realw_test((void **)&(Rsf->V0),V0,*NGLOB_AB);
    copy_todevice_realw_test((void **)&(Rsf->f0),f0,*NGLOB_AB);
    copy_todevice_realw_test((void **)&(Rsf->V_init),V_init,*NGLOB_AB);
    copy_todevice_realw_test((void **)&(Rsf->a),a,*NGLOB_AB);
    copy_todevice_realw_test((void **)&(Rsf->b),b,*NGLOB_AB);
    copy_todevice_realw_test((void **)&(Rsf->L),L,*NGLOB_AB);
    copy_todevice_realw_test((void **)&(Rsf->theta),theta,*NGLOB_AB);
    copy_todevice_realw_test((void **)&(Rsf->T),T,*NGLOB_AB);
    copy_todevice_realw_test((void **)&(Rsf->C),C,*NGLOB_AB);
    copy_todevice_realw_test((void **)&(Rsf->fw),fw,*NGLOB_AB);
    copy_todevice_realw_test((void **)&(Rsf->Vw),Vw,*NGLOB_AB);
  }
}

/* ----------------------------------------------------------------------------------------------- */

extern "C"
void FC_FUNC_(transfer_todevice_swf_data,
              TRANSFER_TODEVICE_SWF_DATA)(long* Fault_pointer,
                                          int *NGLOB_AB,
                                          int *fault_index,
                                          realw* Dc,
                                          realw* mus,
                                          realw* mud,
                                          realw* T,
                                          realw* C,
                                          realw* theta) {

  Fault_solver_dynamics* Fsolver = (Fault_solver_dynamics*)(*Fault_pointer);
  Swf_type* Swf  = &((Fsolver->faults[*fault_index]).swf);
  if (*NGLOB_AB > 0){
    copy_todevice_realw_test((void **)&(Swf->Dc),Dc,*NGLOB_AB);
    copy_todevice_realw_test((void **)&(Swf->mus),mus,*NGLOB_AB);
    copy_todevice_realw_test((void **)&(Swf->mud),mud,*NGLOB_AB);
    copy_todevice_realw_test((void **)&(Swf->Coh),C,*NGLOB_AB);
    copy_todevice_realw_test((void **)&(Swf->T),T,*NGLOB_AB);
    copy_todevice_realw_test((void **)&(Swf->theta),theta,*NGLOB_AB);
  }
}

/* ----------------------------------------------------------------------------------------------- */

extern "C"
void FC_FUNC_(transfer_tohost_rsf_data,
              TRANSFER_TOHOST_RSF_DATA)(long* Fault_pointer,
                                        int *NGLOB_AB,
                                        int *fault_index,
                                        realw* V0,
                                        realw* f0,
                                        realw* V_init,
                                        realw* a,
                                        realw* b,
                                        realw* L,
                                        realw* theta,
                                        realw* T,
                                        realw* C,
                                        realw* fw,
                                        realw* Vw) {

  Fault_solver_dynamics* Fsolver = (Fault_solver_dynamics*)(*Fault_pointer);
  Rsf_type* Rsf  = &((Fsolver->faults[*fault_index]).rsf);
  if (*NGLOB_AB > 0){
    copy_tohost_realw_test((void **)&(Rsf->V0),V0,*NGLOB_AB);
    copy_tohost_realw_test((void **)&(Rsf->f0),f0,*NGLOB_AB);
    copy_tohost_realw_test((void **)&(Rsf->V_init),V_init,*NGLOB_AB);
    copy_tohost_realw_test((void **)&(Rsf->a),a,*NGLOB_AB);
    copy_tohost_realw_test((void **)&(Rsf->b),b,*NGLOB_AB);
    copy_tohost_realw_test((void **)&(Rsf->L),L,*NGLOB_AB);
    copy_tohost_realw_test((void **)&(Rsf->theta),theta,*NGLOB_AB);
    copy_tohost_realw_test((void **)&(Rsf->T),T,*NGLOB_AB);
    copy_tohost_realw_test((void **)&(Rsf->C),C,*NGLOB_AB);
    copy_tohost_realw_test((void **)&(Rsf->fw),fw,*NGLOB_AB);
    copy_tohost_realw_test((void **)&(Rsf->Vw),Vw,*NGLOB_AB);
  }
}

/* ----------------------------------------------------------------------------------------------- */

extern "C"
void FC_FUNC_(transfer_tohost_swf_data,
              TRANSFER_TOHOST_SWF_DATA)(long* Fault_pointer,
                                        int *NGLOB_AB,
                                        int *fault_index,
                                        realw* Dc,
                                        realw* mus,
                                        realw* mud,
                                        realw* T,
                                        realw* C,
                                        realw* theta) {

  Fault_solver_dynamics* Fsolver = (Fault_solver_dynamics*)(*Fault_pointer);
  Swf_type *Swf = &((Fsolver -> faults[*fault_index]).swf);
  if (*NGLOB_AB > 0){
    copy_tohost_realw_test((void **)&(Swf->Dc),Dc,*NGLOB_AB);
    copy_tohost_realw_test((void **)&(Swf->mus),mus,*NGLOB_AB);
    copy_tohost_realw_test((void **)&(Swf->mud),mud,*NGLOB_AB);
    copy_tohost_realw_test((void **)&(Swf->Coh),C,*NGLOB_AB);
    copy_tohost_realw_test((void **)&(Swf->T),T,*NGLOB_AB);
    copy_tohost_realw_test((void **)&(Swf->theta),theta,*NGLOB_AB);

  }
}

/* ----------------------------------------------------------------------------------------------- */

__device__ __forceinline__ double csevl(const double x,const double* cs,int n) {

  int i, ni ;
  double  b0, b1, b2, twox ,result;

  if (n < 1) return -1.0;

  if (n > 1000) return -1.0;

  if (x < -1.1e0 || x > 1.1e0) return -1.0;

  b1 = 0.E0;
  b0 = 0.E0;
  twox = 2.E0 * x;

  for(i=1; i<=n; i++)
  {

      b2 = b1;
      b1 = b0;
      ni = n  - i + 1;
      b0 = twox*b1 - b2 + cs[ni-1];
  }

  result = 0.5E0 * (b0 - b2);
  return result;
}

/* ----------------------------------------------------------------------------------------------- */

__device__ __forceinline__ int  inits(const double* os,int nos,double eta) {
  int i, ii;
  double   err;


  if (nos < 1) return -1.0;


  err = 0.E0;

  for(ii=1; ii<=nos; ii++)
  {
      i = nos  - ii + 1;
      err = err + fabs(os[i-1]);
      if (err > eta) break;
  }

  //  if (i == nos) print *,'warning: Math::inits: eta may be too small'

  return i;
}

/* ----------------------------------------------------------------------------------------------- */

__device__ __forceinline__ double  asinh_slatec(realw x) {

  const double asnhcs[39]= {
      -.12820039911738186343372127359268E+0,  -.58811761189951767565211757138362E-1,
      +.47274654322124815640725249756029E-2,  -.49383631626536172101360174790273E-3,
      +.58506207058557412287494835259321E-4,  -.74669983289313681354755069217188E-5,
      +.10011693583558199265966192015812E-5,  -.13903543858708333608616472258886E-6,
      +.19823169483172793547317360237148E-7,  -.28847468417848843612747272800317E-8,
      +.42672965467159937953457514995907E-9,  -.63976084654366357868752632309681E-10,
      +.96991686089064704147878293131179E-11, -.14844276972043770830246658365696E-11,
      +.22903737939027447988040184378983E-12, -.35588395132732645159978942651310E-13,
      +.55639694080056789953374539088554E-14, -.87462509599624678045666593520162E-15,
      +.13815248844526692155868802298129E-15, -.21916688282900363984955142264149E-16,
      +.34904658524827565638313923706880E-17, -.55785788400895742439630157032106E-18,
      +.89445146617134012551050882798933E-19, -.14383426346571317305551845239466E-19,
      +.23191811872169963036326144682666E-20, -.37487007953314343674570604543999E-21,
      +.60732109822064279404549242880000E-22, -.98599402764633583177370173440000E-23,
      +.16039217452788496315232638293333E-23, -.26138847350287686596716134399999E-24,
      +.42670849606857390833358165333333E-25, -.69770217039185243299730773333333E-26,
      +.11425088336806858659812693333333E-26, -.18735292078860968933021013333333E-27,
      +.30763584414464922794065920000000E-28, -.50577364031639824787046399999999E-29,
      +.83250754712689142224213333333333E-30, -.13718457282501044163925333333333E-30,
      +.22629868426552784104106666666666E-31
  };

  double  aln2 = 0.69314718055994530941723212145818E0;


// series for asnh       on the interval  0.          to  1.00000d+00
//                                        with weighted error   2.19e-17
//                                         log weighted error  16.66
//                               significant figures required  15.60
//                                    decimal places required  17.31
//

  int nterms = 0;
  double  xmax = 0.E0, sqeps = 0.E0;
  double asinh_slatec=0.0E0;

  // taken from http://people.sc.fsu.edu/~jburkardt/f_src/machine/machine.f90
  double d1mach_3 = 1.110223024625157E-016;

  double  y;

  if (nterms == 0){
      nterms = inits(asnhcs, 39, 0.1E0*d1mach_3);
  //nterms = 39;
      sqeps = sqrt(d1mach_3);
      xmax = 1.E0/sqeps;
  }
  y = fabs(x);

  if (y <= 1.E0){
      asinh_slatec = x;
      if (y > sqeps) asinh_slatec = x*(1.E0 )+csevl(2.E0*x*x-1.E0, asnhcs, nterms);
      return asinh_slatec;
  }
  if (y < xmax ) asinh_slatec = log(y + sqrt(y*y + 1.E0F));
  if (y >= xmax) asinh_slatec = aln2 + log(y);
  asinh_slatec = x>0.0 ? fabs(asinh_slatec):-fabs(asinh_slatec);

  return asinh_slatec;
  /*
   April 1977 version.  W. Fullerton, C3, Los Alamos Scientific Lab.
   Evaluate the n-term Chebyshev series cs at x.  Adapted from
   R. Broucke, Algorithm 446, C.A.C.M., 16, 254 (1973).  Also see Fox
   and Parker, Chebyshev polynomials in numerical analysis, Oxford Press, p.56.

               input arguments --
   x      value at which the series is to be evaluated.
   cs     array of n terms of a Chebyshev series.
          in evaluating cs, only half the first coefficient is summed.
   n      number of terms in array cs.*/
}

/* ----------------------------------------------------------------------------------------------- */

__device__ __forceinline__ void funcd(double x,double *fn,double *df,realw tStick,realw Seff,
                                      realw Z,realw f0,realw V0,realw a,realw b,realw L,realw theta,realw cohesion,int statelaw) {
  /*real(kind=CUSTOM_REAL) :: tStick,Seff,Z,f0,V0,a,b,L,theta
  double precision :: arg,fn,df,x
  integer :: statelaw*/
  double arg,xarg;

  if (statelaw == 1){
      arg = exp((f0+b*log(V0*theta/L))/a)/2.0/V0;
  }else{
      arg = exp(theta/a)/2.0E0/V0;
  }
  xarg = x*arg;
  *fn = tStick - Z*x - a*Seff*asinh_slatec(xarg) - cohesion;
  *df = -Z - a*Seff/sqrt(1.0E0 + pow((x*arg),2.0))*arg;
}

/* ----------------------------------------------------------------------------------------------- */

/*// April 1977 version.  W. Fullerton, C3, Los Alamos Scientific Lab.
//
// Initialize the orthogonal series so that inits is the number of terms
// needed to ensure that the error is no larger than eta. Ordinarily, eta
will be chosen to be one-tenth machine precision.
!
!             input arguments --
! os     array of nos coefficients in an orthogonal series.
! nos    number of coefficients in os.
! eta    requested accuracy of series.*/

__device__ __forceinline__ double rtsafe(realw x1,realw x2,realw xacc,realw tStick,realw Seff,realw Z,
        realw f0,realw V0,realw a,realw b,realw L,realw theta,realw cohesion,int statelaw) {

  const int  MAXIT=200;
  int j;
  double   df,dx,dxold,f,fh,fl,temp,xh,xl,rtsafe;
  funcd((double)x1,&fl,&df,tStick,Seff,Z,f0,V0,a,b,L,theta,cohesion,statelaw);
  funcd((double)x2,&fh,&df,tStick,Seff,Z,f0,V0,a,b,L,theta,cohesion,statelaw);
  if ((fl>0. && fh>0.) || (fl<0. && fh<0.) ) return -1.0;
  if (fl==0.){
      rtsafe=x1;
      return rtsafe;
  } else if (fh==0.){
      rtsafe=x2;
      return rtsafe;
  } else if (fl<0.){
      xl=x1;
      xh=x2;
  } else{
      xh=x1;
      xl=x2;
  }
  rtsafe = 0.5E0*(x1+x2);
  dxold = fabsf(x2-x1);
  dx = dxold;
  funcd(rtsafe,&f,&df,tStick,Seff,Z,f0,V0,a,b,L,theta,cohesion,statelaw);
  for(j=1; j<MAXIT; j++) {
    if (((rtsafe-xh)*df-f)*((rtsafe-xl)*df-f)>0 || fabsf(2.0F*f)>fabsf(dxold*df)){
      dxold=dx;
      dx=0.5E0*(xh-xl);
      rtsafe=xl+dx;
      if (xl==rtsafe) return rtsafe;
    }else{
      dxold=dx;
      dx=f/df;
      temp=rtsafe;
      rtsafe=rtsafe-dx;
      if (temp==rtsafe) return rtsafe;
    }
    if (fabsf(dx)<xacc) return rtsafe;
    funcd(rtsafe,&f,&df,tStick,Seff,Z,f0,V0,a,b,L,theta,cohesion,statelaw);
    if (f<0.){
      xl=rtsafe;
    } else {
      xh=rtsafe;
    }
  }
  return -2.0;
}

/* ----------------------------------------------------------------------------------------------- */

__device__ __forceinline__ realw update_state_rsf(realw Ll,
                                                  realw theta,
                                                  realw Vslip,
                                                  realw dt){
  double vDtL;

  realw theta_r;
  vDtL = Vslip*dt/Ll;
  if(vDtL > 1.0e-5){
    theta_r = theta*exp(-vDtL) +
                Ll/Vslip*(1.0 - exp(-vDtL));
  }else{
    theta_r = theta*exp(-vDtL) +
                dt*(1.0 - 0.5*vDtL);
  }
  return theta_r;
}

/* ----------------------------------------------------------------------------------------------- */

__device__ __forceinline__ realw update_state_swf(realw  Dx,
                                                  realw  Dy,
                                                  realw* D_slip,
                                                  int index,
                                                  realw theta_old) {
  realw theta_r;
  theta_r = theta_old + sqrt((Dx-D_slip[index*3])*(Dx-D_slip[index*3])+(Dy-D_slip[index*3+1])*(Dy-D_slip[index*3+1]));
  return theta_r;
}

/* ----------------------------------------------------------------------------------------------- */

__device__ __forceinline__ realw swf_mu(realw Dcl,
                                        realw musl,
                                        realw mudl,
                                        realw thetal) {
  realw mul,tmp;
  tmp = MIN(thetal/Dcl,1.00);
  mul = musl - (musl - mudl)*tmp;
  return mul;
}

/* ----------------------------------------------------------------------------------------------- */

__device__ __forceinline__ void rotate(realw* R,
                                       realw* vrx,
                                       realw* vry,
                                       realw* vrz,
                                       int tx,
                                       int isForward) {
  realw vx,vy,vz;
  vx = *vrx;
  vy = *vry;
  vz = *vrz;
  if(isForward){
// Percy, tangential direction Vt, equation 7 of Pablo's notes in agreement with SPECFEM3D
// forward rotation
    *vrx = vx*R[0+9*tx]+vy*R[3+9*tx]+vz*R[6+9*tx];  //vx
    *vry = vx*R[1+9*tx]+vy*R[4+9*tx]+vz*R[7+9*tx];  //vy
    *vrz = vx*R[2+9*tx]+vy*R[5+9*tx]+vz*R[8+9*tx];  //vz
  }else {
  // backward rotation
    *vrx = vx*R[0+9*tx]+vy*R[1+9*tx]+vz*R[2+9*tx];  //vx
    *vry = vx*R[3+9*tx]+vy*R[4+9*tx]+vz*R[5+9*tx];  //vy
    *vrz = vx*R[6+9*tx]+vy*R[7+9*tx]+vz*R[8+9*tx];  //vz
  }
}

/* ----------------------------------------------------------------------------------------------- */

__device__ __forceinline__ void get_jump(const realw* Vector,realw* Dx, realw* Dy, realw* Dz,int index1,int index2) {
  *Dx = Vector[3*index2] - Vector[3*index1];
  *Dy = Vector[3*index2 + 1] - Vector[3*index1 + 1];
  *Dz = Vector[3*index2 + 2] - Vector[3*index1 + 2];
  return;
}

/* ----------------------------------------------------------------------------------------------- */

__device__ __forceinline__ void get_weighted_jump(const realw* Vector,const realw Weigh1,const realw Weigh2, realw* Dx, realw* Dy, realw* Dz, int index1, int index2){
  *Dx = Vector[3*index2] * Weigh2 - Vector[3*index1] * Weigh1;
  *Dy = Vector[3*index2 + 1] * Weigh2 - Vector[3*index1 + 1] * Weigh1;
  *Dz = Vector[3*index2 + 2] * Weigh2 - Vector[3*index1 + 2] * Weigh1;
  return;
}

/* ----------------------------------------------------------------------------------------------- */

__global__  void compute_dynamic_fault_cuda_swf(realw* Displ, //this is a mesh vector
                                                realw* Veloc,
                                                realw* MxAccel,
                                                int NGLOB_AB,
                                                realw* invM1,  // this is a fault vector
                                                realw* invM2,
                                                realw* B,
                                                realw* Z,
                                                realw* R,
                                                realw* T0,
                                                realw* T,     //for output
                                                realw* Dc,
                                                realw* theta,
                                                realw* mus,
                                                realw* mud,
                                                realw* Coh,
                                                realw* RT,
                                                realw* V_slip,
                                                realw* D_slip,
                                                int* ibulk1,
                                                int* ibulk2,
                                                realw dt,
                                                int myrank) {
  int tx,iglob1,iglob2;
  realw Dx,Dy,Dz,Vx,Vy,Vz,Ax,Ay,Az;
  realw Tx,Ty,Tz,T0xl,T0yl,T0zl;
  realw Tstick;
  realw Zl,mudl,musl,Dcl,thetal,Cohl;
  //realw RTl;
  realw strength;
  realw mul;
  realw thetaold;
  realw Tnew;

  tx = blockDim.x * blockIdx.x + threadIdx.x;  /*calculate thread id*/
  if(tx>=NGLOB_AB) return;

  Zl = Z[tx];

  thetal = theta[tx];
  mudl = mud[tx];
  musl = mus[tx];
  Dcl = Dc[tx];
  Cohl = Coh[tx];
  //RTl = RT[tx];
  T0xl = T0[tx*3];
  T0yl = T0[tx*3+1];
  T0zl = T0[tx*3+2];

  iglob1 = ibulk1[tx]-1;
  iglob2 = ibulk2[tx]-1;

  get_jump(Displ, &Dx, &Dy, &Dz, iglob1, iglob2);
  get_jump(Veloc, &Vx, &Vy, &Vz, iglob1, iglob2);
  get_weighted_jump(MxAccel, invM1[tx], invM2[tx], &Ax, &Ay, &Az,iglob1,iglob2);

  rotate(R,&Dx,&Dy,&Dz,tx,1);
  rotate(R,&Vx,&Vy,&Vz,tx,1);
  rotate(R,&Ax,&Ay,&Az,tx,1);

  Tx = Zl*(Vx + 0.50*dt*Ax);
  Ty = Zl*(Vy + 0.50*dt*Ay);
  Tz = Zl*(Vz + 0.50*dt*Az);

  /*rotate to fault frame*/
  Tx = Tx + T0xl;
  Ty = Ty + T0yl;
  Tz = Tz + T0zl;

  Tstick = sqrt(Tx * Tx + Ty * Ty);

  thetaold = thetal;

  thetal = update_state_swf(Dx,Dy,D_slip,tx,thetaold);

  theta[tx] = thetal;

  mul = swf_mu(Dcl,musl,mudl,thetal);

  strength = -mul * (MIN(Tz,0.00)) + Cohl;

  Tnew = MIN(Tstick,strength);

  Tstick = MAX(Tstick,1.0E0);

  Tx = Tnew * Tx/Tstick;
  Ty = Tnew * Ty/Tstick;

  T[tx*3]   = Tx;
  T[tx*3+1] = Ty;
  T[tx*3+2] = Tz;

  Tx = Tx - T0xl;
  Ty = Ty - T0yl;
  Tz = Tz - T0zl;
  Ax = Ax - Tx/(Zl*0.5E0*dt);
  Ay = Ay - Ty/(Zl*0.5E0*dt);
  Az = Az - Tz/(Zl*0.5E0*dt);

  D_slip[tx*3]   = Dx;
  D_slip[tx*3+1] = Dy;
  D_slip[tx*3+2] = Dz;

  V_slip[tx*3]   = Vx + 0.5E0*dt*Ax;
  V_slip[tx*3+1] = Vy + 0.5E0*dt*Ay;
  V_slip[tx*3+2] = Vz + 0.5E0*dt*Az;
  rotate(R,&Tx,&Ty,&Tz,tx,0);

  MxAccel[3*iglob1] = MxAccel[3*iglob1] + B[tx]*Tx;
  MxAccel[3*iglob2] = MxAccel[3*iglob2] - B[tx]*Tx;

  MxAccel[3*iglob1+1] = MxAccel[3*iglob1+1] + B[tx]*Ty;
  MxAccel[3*iglob2+1] = MxAccel[3*iglob2+1] - B[tx]*Ty;

  MxAccel[3*iglob1+2] = MxAccel[3*iglob1+2] + B[tx]*Tz;
  MxAccel[3*iglob2+2] = MxAccel[3*iglob2+2] - B[tx]*Tz;
}

/* ----------------------------------------------------------------------------------------------- */

__global__  void compute_dynamic_fault_cuda(realw* Displ, /*mesh quantities*/
                                            realw* Veloc,
                                            realw* MxAccel,
                                            int NGLOB_AB,
                                            realw* invM1,  /* fault quantities*/
                                            realw* invM2,
                                            realw* B,
                                            realw* Z,
                                            realw* R,
                                            realw* T0,
                                            realw* T,
                                            realw* Coh,
                                            realw* a,
                                            realw* b,
                                            realw* L,
                                            realw* f0,
                                            realw* V0,    /*frictional quantities*/
                                            realw* V_init,
                                            realw* theta,
                                            realw* Vw,
                                            realw* fw,
                                            realw* V_slip,
                                            realw* D_slip,
                                            int* ibulk1,
                                            int* ibulk2,
                                            realw dt,
                                            int myrank) {
  int tx,iglob1,iglob2;
  realw Dx,Dy,Dz,Vx,Vy,Vz,Ax,Ay,Az;
  realw Tx,Ty,Tz,T0xl,T0yl,T0zl;
  realw Tstick;
  realw Zl,al,bl,Ll,f0l,V0l,thetal;
  //realw V_initl,Vwl,fwl;
  realw thetaold;
  realw Vf_oldl,Vf_newl,Vf_tmp;
  realw Ztmp;
  realw Tnew;
  realw Cohl;
  realw netTstick;

  tx = blockDim.x * blockIdx.x + threadIdx.x;  /*calculate thread id*/
  if(tx>=NGLOB_AB) return;

  Vf_oldl = sqrt(V_slip[3*tx]*V_slip[3*tx]+V_slip[3*tx+1]*V_slip[3*tx+1]);
  Zl = Z[tx];
  al = a[tx];
  bl = b[tx];
  thetal = theta[tx];
  f0l = f0[tx];
  Ll = L[tx];
  //Vwl = Vw[tx];
  //fwl = fw[tx];
  V0l = V0[tx];
  //V_initl=V_init[tx];
  Cohl = Coh[tx];

  T0xl = T0[tx*3];
  T0yl = T0[tx*3+1];
  T0zl = T0[tx*3+2];

  iglob1 = ibulk1[tx]-1;
  iglob2 = ibulk2[tx]-1;

  get_jump(Displ, &Dx, &Dy, &Dz, iglob1, iglob2);
  get_jump(Veloc, &Vx, &Vy, &Vz, iglob1, iglob2);
  get_weighted_jump(MxAccel, invM1[tx], invM2[tx], &Ax, &Ay, &Az,iglob1,iglob2);
  Ztmp = Z[tx];

  rotate(R,&Dx,&Dy,&Dz,tx,1);
  rotate(R,&Vx,&Vy,&Vz,tx,1);
  rotate(R,&Ax,&Ay,&Az,tx,1);

  Tx = Ztmp*(Vx + 0.50*dt*Ax);
  Ty = Ztmp*(Vy + 0.50*dt*Ay);
  Tz = Ztmp*(Vz + 0.50*dt*Az);

  /*rotate back to fault frame*/
  Tx = Tx + T0xl;
  Ty = Ty + T0yl;
  Tz = Tz + T0zl;

  Tstick = sqrt(Tx * Tx + Ty * Ty);

  netTstick = Tstick - Cohl; /** add cohesion into simulation*/

  netTstick = MAX(netTstick, 0.0E0);  /** prevent Tstick from being negative */

  thetaold = thetal;

  thetal = update_state_rsf(Ll ,thetaold , Vf_oldl, dt );

  Vf_newl = (realw)rtsafe(0.0E0, Vf_oldl+5.0E0, 1.0E-5, netTstick, -Tz, Zl, f0l, V0l, al, bl, Ll, thetal, 0.0, 1);

  Vf_tmp = 0.5E0*(Vf_oldl + Vf_newl);

  thetal = update_state_rsf(Ll ,thetaold , Vf_tmp, dt );

  theta[tx] = thetal;



  Vf_newl = (realw)rtsafe(0.0E0, Vf_oldl+5.0E0, 1.0E-5, netTstick, -Tz, Zl, f0l, V0l, al, bl, Ll, thetal, 0.0, 1);

  assert(Vf_newl > -1.0E-6);      /** Double precisio to single precision conversion may cause an error of about 1e-6*/

  Vf_newl = MAX(Vf_newl,0.0E0);

  Tstick = MAX(Tstick,1.0E0);

  Tnew = Tstick - Zl*Vf_newl;


  Tx = Tnew * Tx/Tstick;
  Ty = Tnew * Ty/Tstick;
  T[tx*3]   = Tx;
  T[tx*3+1] = Ty;
  T[tx*3+2] = Tz;

  Tx = Tx - T0xl;
  Ty = Ty - T0yl;
  Tz = Tz - T0zl;

  Ax = Ax - Tx/(Zl*0.5E0*dt);
  Ay = Ay - Ty/(Zl*0.5E0*dt);
  Az = Az - Tz/(Zl*0.5E0*dt);

  D_slip[tx*3]   = Dx;
  D_slip[tx*3+1] = Dy;
  D_slip[tx*3+2] = Dz;

  V_slip[tx*3]   = Vx + 0.5E0*dt*Ax;
  V_slip[tx*3+1] = Vy + 0.5E0*dt*Ay;
  V_slip[tx*3+2] = Vz + 0.5E0*dt*Az;
  rotate(R,&Tx,&Ty,&Tz,tx,0);

  MxAccel[3*iglob1] = MxAccel[3*iglob1] + B[tx]*Tx;
  MxAccel[3*iglob2] = MxAccel[3*iglob2] - B[tx]*Tx;

  MxAccel[3*iglob1+1] = MxAccel[3*iglob1+1] + B[tx]*Ty;
  MxAccel[3*iglob2+1] = MxAccel[3*iglob2+1] - B[tx]*Ty;

  MxAccel[3*iglob1+2] = MxAccel[3*iglob1+2] + B[tx]*Tz;
  MxAccel[3*iglob2+2] = MxAccel[3*iglob2+2] - B[tx]*Tz;
}

/* ----------------------------------------------------------------------------------------------- */
__global__ void store_dataT(realw* store_dataT,
                            realw* V_slip,
                            realw* D_slip,
                            realw* T,
                            int* iglob,
                            int istep,
                            int n_record,
                            int nt) {
  int tx = blockDim.x * blockIdx.x + threadIdx.x;  /*calculate thread id*/
  if(tx>=n_record) return;
  int it = (istep - 1)%nt ;
  int irec = iglob[tx];

  store_dataT[it*n_record*7 + tx*7 + 0] = D_slip[3*irec + 0];
  store_dataT[it*n_record*7 + tx*7 + 1] = V_slip[3*irec + 0];
  store_dataT[it*n_record*7 + tx*7 + 2] = T[3*irec + 0];
  store_dataT[it*n_record*7 + tx*7 + 3] = -D_slip[3*irec + 1];
  store_dataT[it*n_record*7 + tx*7 + 4] = -V_slip[3*irec + 1];
  store_dataT[it*n_record*7 + tx*7 + 5] = -T[3*irec + 1];
  store_dataT[it*n_record*7 + tx*7 + 6] = T[3*irec + 2];
}


extern "C"
void FC_FUNC_(fault_solver_gpu,
              FAULT_SOLVER_GPU)(long* Mesh_pointer,
                                long* Fault_pointer,
                                realw* dt,
                                int* myrank,
                                int* it) {
  Fault_solver_dynamics* Fsolver = (Fault_solver_dynamics*)(*Fault_pointer);
  Fault* Flt = Fsolver->faults;
  Mesh*  mp = (Mesh*)(*Mesh_pointer);
  int num_of_block;

  int num_of_block2;
  for(int ifault = 0; ifault < (Fsolver->Nbfaults); ifault++){
    Flt = &(Fsolver->faults[ifault]);
    Rsf_type* rsf = &(Flt->rsf);
    Swf_type* swf = &(Flt->swf);
    if(Flt->NGLOB_AB > 0){
      num_of_block = (int) (Flt->NGLOB_AB/128)+1;
      num_of_block2 = (int) (Flt->output_dataT->NRECORD/128)+1;
      if (false){ // this is dirty implementation

        compute_dynamic_fault_cuda<<<num_of_block,128>>>( mp->d_displ,
                                                          mp->d_veloc,
                                                          mp->d_accel,
                                                          Flt->NGLOB_AB,
                                                          Flt->invM1,
                                                          Flt->invM2,
                                                          Flt->B,
                                                          Flt->Z,
                                                          Flt->R,
                                                          Flt->T0,
                                                          Flt->T,
                                                          rsf->C,   /**rate and state friction quantities*/
                                                          rsf->a,
                                                          rsf->b,
                                                          rsf->L,
                                                          rsf->f0,
                                                          rsf->V0,
                                                          rsf->V_init,
                                                          rsf->theta,
                                                          rsf->Vw,
                                                          rsf->fw,
                                                          Flt->V,
                                                          Flt->D,
                                                          Flt->ibulk1,
                                                          Flt->ibulk2,
                                                          *dt,
                                                          *myrank);
      }else{
        compute_dynamic_fault_cuda_swf<<<num_of_block,128>>>( mp->d_displ,
                                                              mp->d_veloc,
                                                              mp->d_accel,
                                                              Flt->NGLOB_AB,
                                                              Flt->invM1,
                                                              Flt->invM2,
                                                              Flt->B,
                                                              Flt->Z,
                                                              Flt->R,
                                                              Flt->T0,
                                                              Flt->T,
                                                              swf->Dc,
                                                              swf->theta,
                                                              swf->mus,
                                                              swf->mud,
                                                              swf->Coh,
                                                              swf->T,
                                                              Flt->V,
                                                              Flt->D,
                                                              Flt->ibulk1,
                                                              Flt->ibulk2,
                                                              *dt,
                                                              *myrank);

        store_dataT<<<num_of_block2,128>>>( Flt->output_dataT->dataT,
                                            Flt->V,
                                            Flt->D,
                                            Flt->T,
                                            Flt->output_dataT->iglob,
                                            *it,
                                            Flt->output_dataT->NRECORD,
                                            Flt->output_dataT->NT);

      }
    }
  }
}

