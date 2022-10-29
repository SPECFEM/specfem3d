/*
!=====================================================================
!
!                          S p e c f e m 3 D
!                          -----------------
!
!    Main historical authors: Dimitri Komatitsch and Jeroen Tromp
!                             CNRS, France
!                      and Princeton University, USA
!                (there are currently many more authors!)
!                          (c) October 2017
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


// asserts
#include <assert.h>

/* ----------------------------------------------------------------------------------------------- */

__device__ __forceinline__ double csevl(const double x,const double* cs,int n) {

/*
! April 1977 version.  W. Fullerton, C3, Los Alamos Scientific Lab.
! Evaluate the n-term Chebyshev series cs at x.  Adapted from
! R. Broucke, Algorithm 446, C.A.C.M., 16, 254 (1973).  Also see Fox
! and Parker, Chebyshev polynomials in numerical analysis, Oxford Press, p.56.
!
! input arguments --
! x      value at which the series is to be evaluated.
! cs     array of n terms of a Chebyshev series.
!        in evaluating cs, only half the first coefficient is summed.
! n      number of terms in array cs.
 */

  int i, ni ;
  double  b0, b1, b2, twox ,result;

  if (n < 1)   {assert(n < 1);    return -1.0; }  // 'Math::csevl: number of terms <= 0'
  if (n > 1000){assert(n > 1000); return -1.0; } // 'Math::csevl: number of terms > 1000'
  if (x < -1.1 || x > 1.1){ assert(x < -1.1 || x > 1.1); return -1.0; } // 'Math::csevl: x outside (-1,+1)'

  b1 = 0.0;
  b0 = 0.0;
  twox = 2.0 * x;

  for(i=1; i<=n; i++) {
      b2 = b1;
      b1 = b0;
      ni = n + 1 - i;
      b0 = twox * b1 - b2 + cs[ni-1];
  }

  result = 0.5 * (b0 - b2);
  return result;
}

/* ----------------------------------------------------------------------------------------------- */

__device__ __forceinline__ int  inits(const double* os,int nos,double eta) {

/*
 ! April 1977 version.  W. Fullerton, C3, Los Alamos Scientific Lab.
 !
 ! Initialize the orthogonal series so that inits is the number of terms
 ! needed to ensure that the error is no larger than eta. Ordinarily, eta
 ! will be chosen to be one-tenth machine precision.
 !
 !             input arguments --
 ! os     array of nos coefficients in an orthogonal series.
 ! nos    number of coefficients in os.
 ! eta    requested accuracy of series.
*/

  int i, ii;
  double err;

  if (nos < 1){
    assert(nos < 1);  // stop 'Math::inits: number of terms <= 0'
    return -1.0;
  }

  err = 0.0;

  for(ii=1; ii<=nos; ii++) {
      i = nos + 1 - ii;
      err = err + fabs(os[i-1]);
      if (err > eta) break;
  }

  //  if (i == nos) print *,'warning: Math::inits: eta may be too small'

  return i;
}

/* ----------------------------------------------------------------------------------------------- */

__device__ __forceinline__ double  asinh_slatec(double x) {

  /*
  ! asinh() function taken from Netlib
  ! April 1977 edition.  W. Fullerton, C3, Los Alamos Scientific Lab.

  ! taken from http://www.tddft.org/trac/octopus/browser/trunk/src/asinh.F90?rev=2

  ! and modified by Dimitri Komatitsch in December 2012 for portability
   */

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
  double xmax = 0.0, sqeps = 0.0;
  double asinh_slatec = 0.0;

  // taken from http://people.sc.fsu.edu/~jburkardt/f_src/machine/machine.f90
  double d1mach_3 = 1.110223024625157E-016;

  double y;

  if (nterms == 0){
      nterms = inits(asnhcs, 39, 0.1*d1mach_3);
      sqeps = sqrt(d1mach_3);
      xmax = 1.0/sqeps;
  }
  y = fabs(x);

  if (y <= 1.0){
      asinh_slatec = x;
      if (y > sqeps) asinh_slatec = x + csevl(2.0 * x*x - 1.0, asnhcs, nterms);
      return asinh_slatec;
  }
  if (y < xmax ) asinh_slatec = log(y + sqrt(y*y + 1.0));
  if (y >= xmax) asinh_slatec = aln2 + log(y);

  asinh_slatec = x > 0.0 ? fabs(asinh_slatec) : -fabs(asinh_slatec);

  return asinh_slatec;
}

/* ----------------------------------------------------------------------------------------------- */

__device__ __forceinline__ void funcd(double x,double *fn,double *df,
                                      realw Tstick_in,
                                      realw Seff_in,
                                      realw Z_in,
                                      realw f0_in,
                                      realw V0_in,
                                      realw a_in,
                                      realw b_in,
                                      realw L_in,
                                      realw theta_in,
                                      realw cohesion_in,
                                      int statelaw) {

  double arg,xarg;

  // converts calculations to double precision
  // todo: using double precision calculation for now,
  //       but might be okay to go with realw, i.e., float - to check...
  double Tstick = (double) Tstick_in;
  double Seff = (double) Seff_in;
  double Z = (double) Z_in;
  double f0 = (double) f0_in;
  double V0 = (double) V0_in;
  double a = (double) a_in;
  double b = (double) b_in;
  double L = (double) L_in;
  double theta = (double) theta_in;
  double cohesion = (double) cohesion_in;

  if (statelaw == 1){
    // ageing law
    arg = (double) exp((f0+b*log(V0*theta/L))/a)/2.0/V0;
  }else{
    // slip law
    arg = (double) exp(theta/a)/2.0/V0;
  }
  xarg = x * arg;

  // traction Tau(x)

  // using netlib's asinh_slatec() implementation
  // (not sure why this explicit asinh function is used, seems to be slower...)
  //*fn = Tstick - Z*x - a * Seff * asinh_slatec(xarg) - cohesion;
  //
  // using intrinsic <math.h> asinh() function
  *fn = Tstick - Z*x - a * Seff * asinh(xarg) - cohesion;

  // derivative of traction Tau'(x)
  *df = -Z - a * Seff / sqrt(1.0 + xarg*xarg)*arg;
}

/* ----------------------------------------------------------------------------------------------- */

/*
// April 1977 version.  W. Fullerton, C3, Los Alamos Scientific Lab.
//
// Initialize the orthogonal series so that inits is the number of terms
// needed to ensure that the error is no larger than eta. Ordinarily, eta
will be chosen to be one-tenth machine precision.
!
!             input arguments --
! os     array of nos coefficients in an orthogonal series.
! nos    number of coefficients in os.
! eta    requested accuracy of series.
 */

__device__ __forceinline__ double rtsafe(realw x1_in,
                                         realw x2_in,
                                         realw xacc_in,
                                         realw Tstick,
                                         realw Seff,
                                         realw Z, realw f0,realw V0,
                                         realw a,realw b,
                                         realw L,
                                         realw theta,
                                         realw cohesion,
                                         int statelaw) {

  const int  MAXIT=200; // maximum number of iterations
  int j;
  double df,dx,dxold,f,fh,fl,temp,xh,xl,rtsafe;

  double x1 = (double) x1_in;
  double x2 = (double) x2_in;
  double xacc = (double) xacc_in;

  funcd(x1,&fl,&df,Tstick,Seff,Z,f0,V0,a,b,L,theta,cohesion,statelaw);
  funcd(x2,&fh,&df,Tstick,Seff,Z,f0,V0,a,b,L,theta,cohesion,statelaw);

  if ((fl > 0.0 && fh > 0.0) || (fl < 0.0 && fh < 0.0) ) {
    // case should not occur
    assert((fl > 0.0 && fh > 0.0));
    assert((fl < 0.0 && fh < 0.0));
    rtsafe = -1.0;
    assert(rtsafe == -1.0);
    return rtsafe;
  }

  if (fl == 0.0){               // todo: comparison of float against zero, should add numerical tolerance
      rtsafe = x1;              // todo check: return x2 or x1?
      return rtsafe;
  } else if (fh == 0.0){        // todo: comparison of float against zero, should add numerical tolerance
      rtsafe = x2;
      return rtsafe;
  } else if (fl < 0.0){
      xl = x1;
      xh = x2;
  } else{
      xh = x1;
      xl = x2;
  }

  rtsafe = 0.5 * (x1+x2);
  dxold = fabs(x2-x1);
  dx = dxold;

  funcd(rtsafe,&f,&df,Tstick,Seff,Z,f0,V0,a,b,L,theta,cohesion,statelaw);

  for(j=1; j<=MAXIT; j++) {
    if (((rtsafe-xh)*df-f)*((rtsafe-xl)*df-f) > 0.0 || fabs(2.0 * f) > fabs(dxold*df)){
      dxold = dx;
      dx = 0.5 * (xh-xl);
      rtsafe = xl + dx;
      if (xl == rtsafe) return rtsafe;    // todo: comparison of float against float, should add numerical tolerance
    }else{
      dxold = dx;
      dx = f/df;
      temp = rtsafe;
      rtsafe = rtsafe - dx;
      if (temp == rtsafe) return rtsafe;  // todo: comparison of float against float, should add numerical tolerance
    }

    // check if solution within accuracy xacc
    if (fabs(dx) < xacc) return rtsafe;

    funcd(rtsafe,&f,&df,Tstick,Seff,Z,f0,V0,a,b,L,theta,cohesion,statelaw);

    if (f < 0.0){
      xl = rtsafe;
    } else {
      xh = rtsafe;
    }
  }

  // case should not occur, might need higher number of iterations?
  rtsafe = -2.0;
  assert(rtsafe == -2.0);
  return rtsafe;
}

/* ----------------------------------------------------------------------------------------------- */

__device__ __forceinline__ realw update_state_rsf(realw L_in,
                                                  realw theta_in,
                                                  realw Vslip_in,
                                                  realw dt_in,
                                                  realw f0_in,
                                                  realw fw_in,
                                                  realw a_in,
                                                  realw b_in,
                                                  realw V0_in,
                                                  realw Vw_in,
                                                  int StateLaw){

  realw theta_ret;

  // converts calculations to double precision
  // todo: using double precision calculation for now,
  //       but might be okay to go with realw, i.e., float - to check...
  double L = (double) L_in;
  double theta = (double) theta_in;
  double Vslip = (double) Vslip_in;
  double dt = (double) dt_in;
  double a = (double) a_in;
  double b = (double) b_in;
  double f0 = (double) f0_in;
  double fw = (double) fw_in;
  double V0 = (double) V0_in;
  double Vw = (double) Vw_in;

  double vDtL = Vslip * dt / L;

  // state update
  if (StateLaw == 1){
    // ageing law
    if (vDtL > 1.0e-5){
      theta_ret = theta * exp(-vDtL) + L / Vslip * (1.0 - exp(-vDtL));
    }else{
      theta_ret = theta * exp(-vDtL) + dt * (1.0 - 0.5 * vDtL);
    }
  }else{
    // slip law
    if (Vslip != 0.0){ // todo: comparison of float against zero, should add numerical tolerance
      double fLV = f0 - (b - a) * log(Vslip/V0);
      double x = 1.0 + pow(Vslip/Vw,8.0);
      double f_ss = fw + (fLV - fw) / pow(x,0.125);
      double xi_ss = a * log( 2.0 * V0/Vslip * sinh(f_ss/a) );
      theta_ret = xi_ss + (theta - xi_ss) * exp(-vDtL);
    }else{
      theta_ret = theta;
    }
  }
  return theta_ret;
}

/* ----------------------------------------------------------------------------------------------- */

__device__ __forceinline__ realw update_state_swf(realw  Dx,
                                                  realw  Dy,
                                                  realw* D_slip,
                                                  int index,
                                                  realw theta_old) {
  realw theta_ret;

  // fault state variable theta (magnitude of accumulated slip on fault)
  theta_ret = theta_old + sqrt((Dx-D_slip[index*3])*(Dx-D_slip[index*3])+(Dy-D_slip[index*3+1])*(Dy-D_slip[index*3+1]));

  return theta_ret;
}

/* ----------------------------------------------------------------------------------------------- */

__device__ __forceinline__ realw swf_mu(realw Dcl,
                                        realw musl,
                                        realw mudl,
                                        realw thetal) {
  realw mul,tmp;

  // slip weakening friction law
  tmp = MIN(thetal/Dcl,1.0f);
  mul = musl - (musl - mudl)*tmp;

  return mul;
}

/* ----------------------------------------------------------------------------------------------- */

__device__ __forceinline__ void rotate(realw* R,
                                       realw* vrx,
                                       realw* vry,
                                       realw* vrz,
                                       int id,
                                       int isForward) {
  realw vx,vy,vz;
  vx = *vrx;
  vy = *vry;
  vz = *vrz;

  if(isForward){
    // Percy, tangential direction Vt, equation 7 of Pablo's notes in agreement with SPECFEM3D
    // forward rotation
    *vrx = vx*R[0+9*id]+vy*R[3+9*id]+vz*R[6+9*id];  //vs strike
    *vry = vx*R[1+9*id]+vy*R[4+9*id]+vz*R[7+9*id];  //vd dip
    *vrz = vx*R[2+9*id]+vy*R[5+9*id]+vz*R[8+9*id];  //vn normal direction
  }else {
    // backward rotation
    *vrx = vx*R[0+9*id]+vy*R[1+9*id]+vz*R[2+9*id];  //vx
    *vry = vx*R[3+9*id]+vy*R[4+9*id]+vz*R[5+9*id];  //vy
    *vrz = vx*R[6+9*id]+vy*R[7+9*id]+vz*R[8+9*id];  //vz
  }
}

/* ----------------------------------------------------------------------------------------------- */

__device__ __forceinline__ void get_jump(const realw* Vector,realw* Dx, realw* Dy, realw* Dz,int index1,int index2) {

  *Dx = Vector[3*index2]     - Vector[3*index1];
  *Dy = Vector[3*index2 + 1] - Vector[3*index1 + 1];
  *Dz = Vector[3*index2 + 2] - Vector[3*index1 + 2];

  return;
}

/* ----------------------------------------------------------------------------------------------- */

__device__ __forceinline__ void get_weighted_jump(const realw* Vector,const realw Weigh1,const realw Weigh2,
                                                  realw* dAx, realw* dAy, realw* dAz, int index1, int index2){

  *dAx = Vector[3*index2]     * Weigh2 - Vector[3*index1]     * Weigh1;
  *dAy = Vector[3*index2 + 1] * Weigh2 - Vector[3*index1 + 1] * Weigh1;
  *dAz = Vector[3*index2 + 2] * Weigh2 - Vector[3*index1 + 2] * Weigh1;

  return;
}

/* ----------------------------------------------------------------------------------------------- */

__global__  void compute_dynamic_fault_cuda_swf(realw* Displ,   // this is a mesh vector
                                                realw* Veloc,
                                                realw* MxAccel,
                                                int NGLOB_FLT,
                                                realw* invM1,   // this is a fault vector
                                                realw* invM2,
                                                realw* B,
                                                realw* Z,
                                                realw* R,
                                                realw* T0,
                                                realw* T,       // for output
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
                                                realw dt) {

  int iglob1,iglob2;
  realw Dx,Dy,Dz,Vx,Vy,Vz,Ax,Ay,Az;
  realw Tx,Ty,Tz,T0x,T0y,T0z;
  realw Tstick;
  realw Zl,mudl,musl,Dcl,thetal,Cohl;
  //realw RTl;
  realw strength;
  realw mul;
  realw theta_old;
  realw Tnew;

  // calculate thread id
  int id = threadIdx.x + (blockIdx.x + blockIdx.y*gridDim.x)*blockDim.x;

  // check if anything to do
  if (id >= NGLOB_FLT) return;

  // gets local values for iglob/id point on fault
  Zl = Z[id];

  Dcl = Dc[id];
  thetal = theta[id];
  mudl = mud[id];
  musl = mus[id];
  Cohl = Coh[id];
  //RTl = RT[id];

  // initial stress
  T0x = T0[id*3];
  T0y = T0[id*3+1];
  T0z = T0[id*3+2];

  iglob1 = ibulk1[id]-1;
  iglob2 = ibulk2[id]-1;

  // get predicted values
  get_jump(Displ, &Dx, &Dy, &Dz, iglob1, iglob2);
  get_jump(Veloc, &Vx, &Vy, &Vz, iglob1, iglob2);
  get_weighted_jump(MxAccel, invM1[id], invM2[id], &Ax, &Ay, &Az, iglob1, iglob2);

  // rotate to fault frame
  rotate(R,&Dx,&Dy,&Dz,id,1);
  rotate(R,&Vx,&Vy,&Vz,id,1);
  rotate(R,&Ax,&Ay,&Az,id,1);

  // T_stick
  Tx = Zl*(Vx + 0.5f * dt * Ax);
  Ty = Zl*(Vy + 0.5f * dt * Ay);
  Tz = Zl*(Vz + 0.5f * dt * Az);

  // add initial stress
  Tx = Tx + T0x;
  Ty = Ty + T0y;
  Tz = Tz + T0z;

  // Opening implies free stress
  //if (bc%allow_opening) T(3,:) = min(T(3,:),0.0_CUSTOM_REAL)

  Tstick = sqrt(Tx * Tx + Ty * Ty);

  // slip weakening friction
  theta_old = thetal;
  thetal = update_state_swf(Dx,Dy,D_slip,id,theta_old);

  mul = swf_mu(Dcl,musl,mudl,thetal);

  theta[id] = thetal;

  // update strength
  strength = -mul * (MIN(Tz,0.0f)) + Cohl;

  // solve for shear stress
  Tnew = MIN(Tstick,strength);

  // to avoid division by zero
  Tstick = MAX(Tstick,1.0f);

  Tx = Tnew * Tx/Tstick;
  Ty = Tnew * Ty/Tstick;

  // save total traction
  T[id*3]   = Tx;
  T[id*3+1] = Ty;
  T[id*3+2] = Tz;

  // subtract initial stress
  Tx = Tx - T0x;
  Ty = Ty - T0y;
  Tz = Tz - T0z;

  // update slip acceleration
  Ax = Ax - Tx/(Zl * 0.5f * dt);
  Ay = Ay - Ty/(Zl * 0.5f * dt);
  Az = Az - Tz/(Zl * 0.5f * dt);

  // Update slip and slip rate, in fault frame
  D_slip[id*3]   = Dx;
  D_slip[id*3+1] = Dy;
  D_slip[id*3+2] = Dz; // unused, done for completeness

  V_slip[id*3]   = Vx + 0.5f * dt * Ax;
  V_slip[id*3+1] = Vy + 0.5f * dt * Ay;
  V_slip[id*3+2] = Vz + 0.5f * dt * Az; // unused, done for completeness

  // Rotate tractions back to (x,y,z) frame
  rotate(R,&Tx,&Ty,&Tz,id,0);

  // Add boundary term B*T to M*a
  MxAccel[3*iglob1]   = MxAccel[3*iglob1]   + B[id]*Tx;
  MxAccel[3*iglob1+1] = MxAccel[3*iglob1+1] + B[id]*Ty;
  MxAccel[3*iglob1+2] = MxAccel[3*iglob1+2] + B[id]*Tz;

  MxAccel[3*iglob2]   = MxAccel[3*iglob2]   - B[id]*Tx;
  MxAccel[3*iglob2+1] = MxAccel[3*iglob2+1] - B[id]*Ty;
  MxAccel[3*iglob2+2] = MxAccel[3*iglob2+2] - B[id]*Tz;
}

/* ----------------------------------------------------------------------------------------------- */

__global__  void compute_dynamic_fault_cuda_rsf(realw* Displ,   // mesh quantities
                                                realw* Veloc,
                                                realw* MxAccel,
                                                int NGLOB_FLT,
                                                realw* invM1,   // fault quantities
                                                realw* invM2,
                                                realw* B,
                                                realw* Z,
                                                realw* R,
                                                realw* T0,
                                                realw* T,       // for output
                                                realw* Coh,     // cohesion
                                                realw* a,
                                                realw* b,
                                                realw* L,
                                                realw* f0,
                                                realw* V0,      // frictional quantities
                                                realw* V_init,
                                                realw* theta,
                                                realw* Vw,
                                                realw* fw,
                                                realw* Fload,
                                                int StateLaw,
                                                realw* V_slip,
                                                realw* D_slip,
                                                int* ibulk1,
                                                int* ibulk2,
                                                realw dt,
                                                int it) {

  int iglob1,iglob2;
  realw Dx,Dy,Dz,Vx,Vy,Vz,Ax,Ay,Az;
  realw Tx,Ty,Tz,T0x,T0y,T0z;
  realw Tstick;
  realw Zl,al,bl,Ll,f0l,fwl,V0l,Vwl,thetal;
  //realw V_initl;
  realw theta_old;
  realw Vf_oldl,Vf_newl,Vf_tmp;
  realw Tnew;
  realw Cohl;
  realw netTstick;
  realw TxExt;

  // calculate thread id
  int id = threadIdx.x + (blockIdx.x + blockIdx.y*gridDim.x)*blockDim.x;

  // check if anything to do
  if(id >= NGLOB_FLT) return;

  // gets local values for iglob/id point on fault
  Zl = Z[id];

  al = a[id];
  bl = b[id];
  Ll = L[id];
  thetal = theta[id];

  f0l = f0[id];
  fwl = fw[id];
  V0l = V0[id];
  Vwl = Vw[id];
  //V_initl=V_init[id];

  // note: CPU version doesn't implement cohesion, this is a modification for GPU-only version
  Cohl = Coh[id];

  Vf_oldl = sqrt(V_slip[3*id]*V_slip[3*id] + V_slip[3*id+1]*V_slip[3*id+1]);

  // initial stress
  T0x = T0[id*3];
  T0y = T0[id*3+1];
  T0z = T0[id*3+2];

  iglob1 = ibulk1[id]-1;
  iglob2 = ibulk2[id]-1;

  // get predicted values
  get_jump(Displ, &Dx, &Dy, &Dz, iglob1, iglob2);
  get_jump(Veloc, &Vx, &Vy, &Vz, iglob1, iglob2);
  get_weighted_jump(MxAccel, invM1[id], invM2[id], &Ax, &Ay, &Az, iglob1, iglob2);

  // rotate to fault frame
  rotate(R,&Dx,&Dy,&Dz,id,1);
  rotate(R,&Vx,&Vy,&Vz,id,1);
  rotate(R,&Ax,&Ay,&Az,id,1);

  // T_stick
  Tx = Zl*(Vx + 0.5f * dt * Ax);
  Ty = Zl*(Vy + 0.5f * dt * Ay);
  Tz = Zl*(Vz + 0.5f * dt * Az);

  // add initial stress
  Tx = Tx + T0x;
  Ty = Ty + T0y;
  Tz = Tz + T0z;

  // Opening implies free stress
  //if (bc%allow_opening) T(3,:) = min(T(3,:),0.0_CUSTOM_REAL)

  // smooth loading within nucleation patch
  //WARNING : ad hoc for SCEC benchmark TPV10x
  {
    realw TLoad = 1.0f;
    realw DTau0 = 1.0f;
    realw GLoad = 1.0f;
    realw timeval = it * dt; // time will never be zero. it starts from 1
    if (timeval <= TLoad){
      GLoad = exp( (timeval-TLoad)*(timeval-TLoad) / (timeval*(timeval - 2.0f * TLoad)) );
    }
    realw Floadl = Fload[id];
    TxExt = DTau0 * Floadl * GLoad;
    Tx = Tx + TxExt;
  }

  Tstick = sqrt(Tx * Tx + Ty * Ty);

  // add cohesion into simulation
  if (1 == 0){
    // GPU version modification
    // adds cohesion
    netTstick = Tstick - Cohl;
    // prevent Tstick from being negative
    netTstick = MAX(netTstick, 0.0E0);
  }else{
    // corresponds to CPU-version
    netTstick = Tstick;
  }

  // rate and state friction

  // input for root finding
  realw x1l = 0.0f;             // lower bound
  realw x2l = Vf_oldl + 5.0f;   // upper bound for finding root
  realw xaccl = 1.0E-5;         // numerical accuracy
  realw loc_coh = 0.0f;         // no cohesion case

  // the solver below can be refactored into a loop with two passes
  // first pass
  theta_old = thetal;
  thetal = update_state_rsf(Ll ,theta_old , Vf_oldl, dt, f0l, fwl, al, bl, V0l, Vwl, StateLaw);

  Vf_newl = (realw) rtsafe(x1l, x2l, xaccl, netTstick, -Tz, Zl, f0l, V0l, al, bl, Ll, thetal, loc_coh, StateLaw);

  // second pass
  Vf_tmp = 0.5f * (Vf_oldl + Vf_newl);
  thetal = update_state_rsf(Ll ,theta_old , Vf_tmp, dt, f0l, fwl, al, bl, V0l, Vwl, StateLaw);

  Vf_newl = (realw) rtsafe(x1l, x2l, xaccl, netTstick, -Tz, Zl, f0l, V0l, al, bl, Ll, thetal, loc_coh, StateLaw);

  // save state
  theta[id] = thetal;

  // Double precision to single precision conversion may cause an error of about 1e-6
  //printf("debug: compute_dynamic_fault_cuda_rsf: thread id %i %i Vf_newl %f\n",id,NGLOB_FLT,Vf_newl);
  //assert(Vf_newl > -1.0E-6);

  // updates stress
  if (1 == 0){
    // GPU version modification
    // prevent from being negative
    Vf_newl = MAX(Vf_newl,0.0f);
    Tnew = Tstick - Zl * Vf_newl;   // todo: check if Tstick or netTstick?
  }else{
    // corresponds to CPU-version
    Tnew = Tstick - Zl * Vf_newl;
  }

  // to avoid division by zero
  Tstick = MAX(Tstick,1.0f);

  Tx = Tnew * Tx/Tstick;
  Ty = Tnew * Ty/Tstick;

  // save total traction
  T[id*3]   = Tx;
  T[id*3+1] = Ty;
  T[id*3+2] = Tz;

  // subtract initial stress
  Tx = Tx - T0x;
  Ty = Ty - T0y;
  Tz = Tz - T0z;

  Tx = Tx - TxExt;
  //JPA: this eliminates the effect of TxExt on the equations of motion. Why is it needed?

  // update slip acceleration
  Ax = Ax - Tx/(Zl * 0.5f * dt);
  Ay = Ay - Ty/(Zl * 0.5f * dt);
  Az = Az - Tz/(Zl * 0.5f * dt);

  // Update slip and slip rate, in fault frame
  D_slip[id*3]   = Dx;
  D_slip[id*3+1] = Dy;
  D_slip[id*3+2] = Dz; // unused, done for completeness

  V_slip[id*3]   = Vx + 0.5f * dt * Ax;
  V_slip[id*3+1] = Vy + 0.5f * dt * Ay;
  V_slip[id*3+2] = Vz + 0.5f * dt * Az; // unused, done for completeness

  // Rotate tractions back to (x,y,z) frame
  rotate(R,&Tx,&Ty,&Tz,id,0);

  // Add boundary term B*T to M*a
  MxAccel[3*iglob1]   = MxAccel[3*iglob1]   + B[id]*Tx;
  MxAccel[3*iglob1+1] = MxAccel[3*iglob1+1] + B[id]*Ty;
  MxAccel[3*iglob1+2] = MxAccel[3*iglob1+2] + B[id]*Tz;

  MxAccel[3*iglob2]   = MxAccel[3*iglob2]   - B[id]*Tx;
  MxAccel[3*iglob2+1] = MxAccel[3*iglob2+1] - B[id]*Ty;
  MxAccel[3*iglob2+2] = MxAccel[3*iglob2+2] - B[id]*Tz;
}

/* ----------------------------------------------------------------------------------------------- */

__global__ void store_dataT(realw* dataT,
                            realw* V_slip,
                            realw* D_slip,
                            realw* T,
                            int RATE_AND_STATE,
                            int StateLaw,
                            realw* theta,
                            int* iglob,
                            int it_step,
                            int n_record,
                            int nt) {

  // calculate thread record
  int irec = threadIdx.x + (blockIdx.x + blockIdx.y*gridDim.x)*blockDim.x;

  if(irec >= n_record) return;

  int it = (it_step - 1)%nt ;
  int id = iglob[irec] - 1; // fortran to C array indexing -> iglob[irec]-1

  int recordlength;
  if (RATE_AND_STATE){
    recordlength = 8;
  }else{
    recordlength = 7;
  }

  dataT[it*n_record*recordlength + irec*recordlength + 0] = D_slip[3*id];           // horizontal right-lateral slip (m)
  dataT[it*n_record*recordlength + irec*recordlength + 1] = V_slip[3*id];           // horizontal right-lateral slip rate (m/s)
  dataT[it*n_record*recordlength + irec*recordlength + 2] = T[3*id] / 1.0e6;        // horizontal right-lateral shear stress (MPa)

  dataT[it*n_record*recordlength + irec*recordlength + 3] = D_slip[3*id + 1];      // vertical up-dip slip (m)
  dataT[it*n_record*recordlength + irec*recordlength + 4] = V_slip[3*id + 1];      // vertical up-dip slip rate (m/s)
  dataT[it*n_record*recordlength + irec*recordlength + 5] = T[3*id + 1] / 1.0e6;   // vertical up-dip shear stress (MPa)

  dataT[it*n_record*recordlength + irec*recordlength + 6] = T[3*id + 2] / 1.0e6;    // normal stress (MPa)

  // for RATE_AND_STATE, last storage is be theta or log(theta)
  //
  //if (bc%rsf%StateLaw == 1) then
  //  ! ageing law
  //  bc%dataT%dat(8,ipoin,it) = log10(theta_new(iglob))
  //else
  //  ! slip law
  //  bc%dataT%dat(8,ipoin,it) = theta_new(iglob)
  //endif
  if (RATE_AND_STATE){
    realw theta_new;
    if (StateLaw == 1){
      // ageing law
      theta_new = log10(theta[id]);
    }else{
      // slip law
      theta_new = theta[id];
    }
    dataT[it*n_record*recordlength + irec*recordlength + 7] = theta_new;    // log10 of state variable (log-seconds)
  }
}


