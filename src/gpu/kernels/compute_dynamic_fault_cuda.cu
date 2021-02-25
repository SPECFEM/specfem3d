/*
!=====================================================================
!
!               S p e c f e m 3 D  V e r s i o n  3 . 0
!               ---------------------------------------
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

  int i, ni ;
  double  b0, b1, b2, twox ,result;

  if (n < 1) return -1.0;
  if (n > 1000) return -1.0;
  if (x < -1.1e0 || x > 1.1e0) return -1.0;

  b1 = 0.E0;
  b0 = 0.E0;
  twox = 2.E0 * x;

  for(i=1; i<=n; i++) {
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

  for(ii=1; ii<=nos; ii++) {
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

__device__ __forceinline__ void funcd(double x,double *fn,double *df,
                                      realw tStick,realw Seff,
                                      realw Z,realw f0,realw V0,
                                      realw a,realw b,
                                      realw L,realw theta,
                                      realw cohesion,
                                      int statelaw) {
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

__device__ __forceinline__ double rtsafe(realw x1,realw x2,
                                         realw xacc,
                                         realw tStick,realw Seff,
                                         realw Z, realw f0,realw V0,
                                         realw a,realw b,
                                         realw L,realw theta,
                                         realw cohesion,
                                         int statelaw) {

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
    theta_r = theta*exp(-vDtL) + Ll/Vslip*(1.0 - exp(-vDtL));
  }else{
    theta_r = theta*exp(-vDtL) + dt*(1.0 - 0.5*vDtL);
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
                                       int id,
                                       int isForward) {
  realw vx,vy,vz;
  vx = *vrx;
  vy = *vry;
  vz = *vrz;

  if(isForward){
// Percy, tangential direction Vt, equation 7 of Pablo's notes in agreement with SPECFEM3D
// forward rotation
    *vrx = vx*R[0+9*id]+vy*R[3+9*id]+vz*R[6+9*id];  //vx
    *vry = vx*R[1+9*id]+vy*R[4+9*id]+vz*R[7+9*id];  //vy
    *vrz = vx*R[2+9*id]+vy*R[5+9*id]+vz*R[8+9*id];  //vz
  }else {
  // backward rotation
    *vrx = vx*R[0+9*id]+vy*R[1+9*id]+vz*R[2+9*id];  //vx
    *vry = vx*R[3+9*id]+vy*R[4+9*id]+vz*R[5+9*id];  //vy
    *vrz = vx*R[6+9*id]+vy*R[7+9*id]+vz*R[8+9*id];  //vz
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

__global__  void compute_dynamic_fault_cuda_swf(realw* Displ,   // this is a mesh vector
                                                realw* Veloc,
                                                realw* MxAccel,
                                                int NGLOB_AB,
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
                                                realw dt,
                                                int myrank) {

  int iglob1,iglob2;
  realw Dx,Dy,Dz,Vx,Vy,Vz,Ax,Ay,Az;
  realw Tx,Ty,Tz,T0xl,T0yl,T0zl;
  realw Tstick;
  realw Zl,mudl,musl,Dcl,thetal,Cohl;
  //realw RTl;
  realw strength;
  realw mul;
  realw thetaold;
  realw Tnew;

  // calculate thread id
  int id = threadIdx.x + (blockIdx.x + blockIdx.y*gridDim.x)*blockDim.x;

  // check if anything to do
  if (id >= NGLOB_AB) return;

  Zl = Z[id];

  thetal = theta[id];
  mudl = mud[id];
  musl = mus[id];
  Dcl = Dc[id];
  Cohl = Coh[id];
  //RTl = RT[id];

  T0xl = T0[id*3];
  T0yl = T0[id*3+1];
  T0zl = T0[id*3+2];

  iglob1 = ibulk1[id]-1;
  iglob2 = ibulk2[id]-1;

  get_jump(Displ, &Dx, &Dy, &Dz, iglob1, iglob2);
  get_jump(Veloc, &Vx, &Vy, &Vz, iglob1, iglob2);
  get_weighted_jump(MxAccel, invM1[id], invM2[id], &Ax, &Ay, &Az,iglob1,iglob2);

  // rotate to fault frame
  rotate(R,&Dx,&Dy,&Dz,id,1);
  rotate(R,&Vx,&Vy,&Vz,id,1);
  rotate(R,&Ax,&Ay,&Az,id,1);

  // T_stick
  Tx = Zl*(Vx + 0.50*dt*Ax);
  Ty = Zl*(Vy + 0.50*dt*Ay);
  Tz = Zl*(Vz + 0.50*dt*Az);

  Tx = Tx + T0xl;
  Ty = Ty + T0yl;
  Tz = Tz + T0zl;

  Tstick = sqrt(Tx * Tx + Ty * Ty);

  // slip weakening friction
  thetaold = thetal;

  thetal = update_state_swf(Dx,Dy,D_slip,id,thetaold);

  theta[id] = thetal;

  mul = swf_mu(Dcl,musl,mudl,thetal);

  // update strength
  strength = -mul * (MIN(Tz,0.00)) + Cohl;

  // solve for shear stress
  Tnew = MIN(Tstick,strength);

  // to avoid division by zero
  Tstick = MAX(Tstick,1.0E0);

  Tx = Tnew * Tx/Tstick;
  Ty = Tnew * Ty/Tstick;

  // save total traction
  T[id*3]   = Tx;
  T[id*3+1] = Ty;
  T[id*3+2] = Tz;

  Tx = Tx - T0xl;
  Ty = Ty - T0yl;
  Tz = Tz - T0zl;

  // update slip acceleration
  Ax = Ax - Tx/(Zl*0.5E0*dt);
  Ay = Ay - Ty/(Zl*0.5E0*dt);
  Az = Az - Tz/(Zl*0.5E0*dt);

  // Update slip and slip rate, in fault frame
  D_slip[id*3]   = Dx;
  D_slip[id*3+1] = Dy;
  D_slip[id*3+2] = Dz;

  V_slip[id*3]   = Vx + 0.5E0*dt*Ax;
  V_slip[id*3+1] = Vy + 0.5E0*dt*Ay;
  V_slip[id*3+2] = Vz + 0.5E0*dt*Az;

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
                                                int NGLOB_AB,
                                                realw* invM1,   // fault quantities
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
                                                realw* V0,      // frictional quantities
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

  int iglob1,iglob2;
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

  // calculate thread id
  int id = threadIdx.x + (blockIdx.x + blockIdx.y*gridDim.x)*blockDim.x;

  // check if anything to do
  if(id >= NGLOB_AB) return;

  Vf_oldl = sqrt(V_slip[3*id]*V_slip[3*id] + V_slip[3*id+1]*V_slip[3*id+1]);

  Zl = Z[id];
  al = a[id];
  bl = b[id];
  thetal = theta[id];
  f0l = f0[id];
  Ll = L[id];
  //Vwl = Vw[id];
  //fwl = fw[id];
  V0l = V0[id];
  //V_initl=V_init[id];
  Cohl = Coh[id];

  T0xl = T0[id*3];
  T0yl = T0[id*3+1];
  T0zl = T0[id*3+2];

  iglob1 = ibulk1[id]-1;
  iglob2 = ibulk2[id]-1;

  get_jump(Displ, &Dx, &Dy, &Dz, iglob1, iglob2);
  get_jump(Veloc, &Vx, &Vy, &Vz, iglob1, iglob2);
  get_weighted_jump(MxAccel, invM1[id], invM2[id], &Ax, &Ay, &Az,iglob1,iglob2);

  Ztmp = Z[id];

  // rotate to fault frame
  rotate(R,&Dx,&Dy,&Dz,id,1);
  rotate(R,&Vx,&Vy,&Vz,id,1);
  rotate(R,&Ax,&Ay,&Az,id,1);

  // T_stick
  Tx = Ztmp*(Vx + 0.50*dt*Ax);
  Ty = Ztmp*(Vy + 0.50*dt*Ay);
  Tz = Ztmp*(Vz + 0.50*dt*Az);

  Tx = Tx + T0xl;
  Ty = Ty + T0yl;
  Tz = Tz + T0zl;

  Tstick = sqrt(Tx * Tx + Ty * Ty);

  // add cohesion into simulation
  netTstick = Tstick - Cohl;
  // prevent Tstick from being negative
  netTstick = MAX(netTstick, 0.0E0);

  // rate and state friction
  thetaold = thetal;

  // first pass
  thetal = update_state_rsf(Ll ,thetaold , Vf_oldl, dt );

  Vf_newl = (realw)rtsafe(0.0E0, Vf_oldl+5.0E0, 1.0E-5, netTstick, -Tz, Zl, f0l, V0l, al, bl, Ll, thetal, 0.0, 1);

  // second pass
  Vf_tmp = 0.5E0*(Vf_oldl + Vf_newl);

  thetal = update_state_rsf(Ll ,thetaold , Vf_tmp, dt );

  theta[id] = thetal;

  Vf_newl = (realw)rtsafe(0.0E0, Vf_oldl+5.0E0, 1.0E-5, netTstick, -Tz, Zl, f0l, V0l, al, bl, Ll, thetal, 0.0, 1);

  // Double precision to single precision conversion may cause an error of about 1e-6
  assert(Vf_newl > -1.0E-6);

  // prevent from being negative
  Vf_newl = MAX(Vf_newl,0.0E0);

  // to avoid division by zero
  Tstick = MAX(Tstick,1.0E0);

  Tnew = Tstick - Zl*Vf_newl;

  Tx = Tnew * Tx/Tstick;
  Ty = Tnew * Ty/Tstick;

  // save total traction
  T[id*3]   = Tx;
  T[id*3+1] = Ty;
  T[id*3+2] = Tz;

  Tx = Tx - T0xl;
  Ty = Ty - T0yl;
  Tz = Tz - T0zl;

  // update slip acceleration
  Ax = Ax - Tx/(Zl*0.5E0*dt);
  Ay = Ay - Ty/(Zl*0.5E0*dt);
  Az = Az - Tz/(Zl*0.5E0*dt);

  // Update slip and slip rate, in fault frame
  D_slip[id*3]   = Dx;
  D_slip[id*3+1] = Dy;
  D_slip[id*3+2] = Dz;

  V_slip[id*3]   = Vx + 0.5E0*dt*Ax;
  V_slip[id*3+1] = Vy + 0.5E0*dt*Ay;
  V_slip[id*3+2] = Vz + 0.5E0*dt*Az;

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

__global__ void store_dataT(realw* store_dataT,
                            realw* V_slip,
                            realw* D_slip,
                            realw* T,
                            int* iglob,
                            int istep,
                            int n_record,
                            int nt) {


  // calculate thread id
  int id = threadIdx.x + (blockIdx.x + blockIdx.y*gridDim.x)*blockDim.x;

  if(id >= n_record) return;

  int it = (istep - 1)%nt ;
  int irec = iglob[id];

  store_dataT[it*n_record*7 + id*7 + 0] = D_slip[3*irec + 0];
  store_dataT[it*n_record*7 + id*7 + 1] = V_slip[3*irec + 0];
  store_dataT[it*n_record*7 + id*7 + 2] = T[3*irec + 0];
  store_dataT[it*n_record*7 + id*7 + 3] = -D_slip[3*irec + 1];
  store_dataT[it*n_record*7 + id*7 + 4] = -V_slip[3*irec + 1];
  store_dataT[it*n_record*7 + id*7 + 5] = -T[3*irec + 1];
  store_dataT[it*n_record*7 + id*7 + 6] = T[3*irec + 2];
}


