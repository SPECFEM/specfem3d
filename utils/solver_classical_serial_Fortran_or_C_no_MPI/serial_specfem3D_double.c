
/*
!=====================================================================
!
!          S p e c f e m 3 D  G l o b e  V e r s i o n  4 . 0
!          --------------------------------------------------
!
!          Main authors: Dimitri Komatitsch and Jeroen Tromp
!    Seismological Laboratory, California Institute of Technology, USA
!             and University of Pau / CNRS / INRIA, France
! (c) California Institute of Technology and University of Pau / CNRS / INRIA
!                            February 2008
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

//
// All the arrays below use static memory allocation,
// using constant sizes defined in values_from_mesher.h.
// This is done purposely to improve performance (Fortran compilers
// can optimize much more when the size of the loops and arrays
// is known at compile time).
// NGLLX, NGLLY and NGLLZ are set equal to 5,
// therefore each element contains NGLLX * NGLLY * NGLLZ = 125 points.
//

//
// All the calculations are done in single precision.
// We do not need double precision in SPECFEM3D.
//

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>

// include values created by the mesher
// done for performance only using static allocation to allow for loop unrolling
#include "DATABASES_FOR_SOLVER/values_from_mesher_C.h"

// constant value of the time step in the main time loop
#define deltatover2 0.5*deltat
#define deltatsqover2 0.5*deltat*deltat

// for the source time function
#define pi 3.14159265
#define f0 (1. / 50.)
#define t0 (1.2 / f0)
#define a pi*pi*f0*f0

// number of GLL integration points in each direction of an element (degree plus one)
#define NGLLX 5
#define NGLLY 5
#define NGLLZ 5

// 3-D simulation
#define NDIM 3

// displacement threshold above which we consider that the code became unstable
#define STABILITY_THRESHOLD 1.e+25

#define VERYSMALLVAL 1.e-24
#define NTSTEP_BETWEEN_OUTPUT_INFO 100

// approximate density of the geophysical medium in which the source is located
// this value is only a constant scaling factor therefore it does not really matter
#define rho 4500.

int main(int argc, char *argv[])
{

// global displacement, velocity and acceleration vectors
 static double displ[NGLOB][NDIM];
 static double veloc[NGLOB][NDIM];
 static double accel[NGLOB][NDIM];

// global diagonal mass matrix
 static double rmass_inverse[NGLOB];

// record a seismogram to check that the simulation went well
 static double seismogram[NSTEP];

// arrays with mesh parameters per slice
 static int ibool[NSPEC][NGLLZ][NGLLY][NGLLX];

 static double xix[NSPEC][NGLLZ][NGLLY][NGLLX];
 static double xiy[NSPEC][NGLLZ][NGLLY][NGLLX];
 static double xiz[NSPEC][NGLLZ][NGLLY][NGLLX];

 static double etax[NSPEC][NGLLZ][NGLLY][NGLLX];
 static double etay[NSPEC][NGLLZ][NGLLY][NGLLX];
 static double etaz[NSPEC][NGLLZ][NGLLY][NGLLX];

 static double gammax[NSPEC][NGLLZ][NGLLY][NGLLX];
 static double gammay[NSPEC][NGLLZ][NGLLY][NGLLX];
 static double gammaz[NSPEC][NGLLZ][NGLLY][NGLLX];

 static double kappav[NSPEC][NGLLZ][NGLLY][NGLLX];
 static double muv[NSPEC][NGLLZ][NGLLY][NGLLX];

 static double dummyx_loc[NGLLZ][NGLLY][NGLLX];
 static double dummyy_loc[NGLLZ][NGLLY][NGLLX];
 static double dummyz_loc[NGLLZ][NGLLY][NGLLX];

// array with derivatives of Lagrange polynomials and precalculated products
 static double hprime_xx[NGLLX][NGLLX];
 static double hprimewgll_xx[NGLLX][NGLLX];
 static double wgllwgll_xy[NGLLY][NGLLX];
 static double wgllwgll_xz[NGLLZ][NGLLX];
 static double wgllwgll_yz[NGLLZ][NGLLY];

 static double tempx1[NGLLZ][NGLLY][NGLLX];
 static double tempx2[NGLLZ][NGLLY][NGLLX];
 static double tempx3[NGLLZ][NGLLY][NGLLX];
 static double tempy1[NGLLZ][NGLLY][NGLLX];
 static double tempy2[NGLLZ][NGLLY][NGLLX];
 static double tempy3[NGLLZ][NGLLY][NGLLX];
 static double tempz1[NGLLZ][NGLLY][NGLLX];
 static double tempz2[NGLLZ][NGLLY][NGLLX];
 static double tempz3[NGLLZ][NGLLY][NGLLX];

// time step
 int it;

 clock_t timeloop_begin;
 float timeloop_total;

 int ispec,iglob,i,j,k,l;

 double xixl,xiyl,xizl,etaxl,etayl,etazl,gammaxl,gammayl,gammazl,jacobianl;
 double duxdxl,duxdyl,duxdzl,duydxl,duydyl,duydzl,duzdxl,duzdyl,duzdzl;
 double duxdxl_plus_duydyl,duxdxl_plus_duzdzl,duydyl_plus_duzdzl;
 double duxdyl_plus_duydxl,duzdxl_plus_duxdzl,duzdyl_plus_duydzl;
 double sigma_xx,sigma_yy,sigma_zz,sigma_xy,sigma_xz,sigma_yz;
 double hp1,hp2,hp3,fac1,fac2,fac3,lambdal,mul,lambdalplus2mul,kappal;
 double tempx1l,tempx2l,tempx3l,tempy1l,tempy2l,tempy3l,tempz1l,tempz2l,tempz3l;

 double Usolidnorm,current_value,time,memory_size;

// to read external files
 FILE *IIN;

// estimate of total memory size used
 printf("\nNSPEC = %d\n",NSPEC);
 printf("NGLOB = %d\n\n",NGLOB);
 printf("NSTEP = %d\n",NSTEP);
 printf("deltat = %lf\n\n",deltat);

// estimate total memory size (the size of a real number is 4 bytes)
// we perform the calculation in single precision rather than integer
// to avoid integer overflow in the case of very large meshes
 memory_size = 4. * ((3.*NDIM + 1.) * NGLOB + 12. * (double)(NGLLX*NGLLY*NGLLZ)*(double)(NSPEC));
 printf("approximate total memory size used = %lf Mb\n\n",memory_size/1024./1024.);

// make sure the source element number is an integer
 if(NSPEC % 2 != 0) {
         fprintf(stderr,"source element number is not an integer, exiting...\n");
         exit(1);
       }

 printf("reading file DATABASES_FOR_SOLVER/proc000000_reg1_database.dat\n");
// read the mesh from external file
 if((IIN=fopen("DATABASES_FOR_SOLVER/proc000000_reg1_database.dat","r"))==NULL) {
         fprintf(stderr,"Cannot open file DATABASES_FOR_SOLVER/proc000000_reg1_database.dat, exiting...\n");
         exit(1);
       }

 for (ispec=0;ispec<NSPEC;ispec++) {
   for (k=0;k<NGLLZ;k++) {
     for (j=0;j<NGLLY;j++) {
       for (i=0;i<NGLLX;i++) {
// read real numbers here
         fscanf(IIN, "%le\n", &xix[ispec][k][j][i]);
         fscanf(IIN, "%le\n", &xiy[ispec][k][j][i]);
         fscanf(IIN, "%le\n", &xiz[ispec][k][j][i]);
         fscanf(IIN, "%le\n", &etax[ispec][k][j][i]);
         fscanf(IIN, "%le\n", &etay[ispec][k][j][i]);
         fscanf(IIN, "%le\n", &etaz[ispec][k][j][i]);
         fscanf(IIN, "%le\n", &gammax[ispec][k][j][i]);
         fscanf(IIN, "%le\n", &gammay[ispec][k][j][i]);
         fscanf(IIN, "%le\n", &gammaz[ispec][k][j][i]);
         fscanf(IIN, "%le\n", &kappav[ispec][k][j][i]);
         fscanf(IIN, "%le\n", &muv[ispec][k][j][i]);

// read an integer here
         fscanf(IIN, "%d\n", &ibool[ispec][k][j][i]);
// subtract one because indices start at zero in C but this array was created by a Fortran
// program and therefore starts at one in the file stored on the disk
         ibool[ispec][k][j][i]--;
       }
     }
   }
 }
 for (i=0;i<NGLOB;i++) {
   fscanf(IIN, "%le\n", &rmass_inverse[i]);
// the real exactly diagonal mass matrix is read (not its inverse)
// therefore invert it here once and for all
   rmass_inverse[i] = 1.f / rmass_inverse[i];
 }
 fclose(IIN);

 printf("reading file DATABASES_FOR_SOLVER/matrices.dat\n");
// read the derivation matrices from external file
 if((IIN=fopen("DATABASES_FOR_SOLVER/matrices.dat","r"))==NULL) {
         fprintf(stderr,"Cannot open file DATABASES_FOR_SOLVER/matrices.dat, exiting...\n");
         exit(1);
       }

 for (j=0;j<NGLLY;j++) {
   for (i=0;i<NGLLX;i++) {
     fscanf(IIN, "%le\n", &hprime_xx[j][i]);
     fscanf(IIN, "%le\n", &hprimewgll_xx[j][i]);
     fscanf(IIN, "%le\n", &wgllwgll_yz[j][i]);
     fscanf(IIN, "%le\n", &wgllwgll_xz[j][i]);
     fscanf(IIN, "%le\n", &wgllwgll_xy[j][i]);
   }
 }
 fclose(IIN);

// clear initial vectors before starting the time loop
// (can remain serial because done only once before entering the time loop)
 for (i=0;i<NGLOB;i++) {
   displ[i][0] = 1; // VERYSMALLVAL;
   displ[i][1] = 1; // VERYSMALLVAL;
   displ[i][2] = 1; // VERYSMALLVAL;

   veloc[i][0] = 1; // 0.;
   veloc[i][1] = 1; // 0.;
   veloc[i][2] = 1; // 0.;

   accel[i][0] = 1; // 0.;
   accel[i][1] = 1; // 0.;
   accel[i][2] = 1; // 0.;
 }

 printf("starting the time loop\n");

 timeloop_begin = clock();

// start of the time loop (which must remain serial obviously)
 for (it=1;it<=NSTEP;it++) {

// compute maximum of norm of displacement from time to time and display it
// in order to monitor the simulation
// this can remain serial because it is done only every 200 time steps
   if((it % NTSTEP_BETWEEN_OUTPUT_INFO) == 0 || it == 5 || it == NSTEP) {

     Usolidnorm = -1.;

     for (iglob = 0; iglob < NGLOB; iglob++) {
       current_value = sqrt(displ[iglob][0]*displ[iglob][0] + displ[iglob][1]*displ[iglob][1] + displ[iglob][2]*displ[iglob][2]);
       if(current_value > Usolidnorm) { Usolidnorm = current_value; }
     }

     printf("\nTime step # %d out of %d\n",it,NSTEP);
// compute current time
     time = (it-1)*deltat;
     printf("Max norm displacement vector U in the solid (m) = %.8g\n",Usolidnorm);
     timeloop_total = ((clock()-timeloop_begin)/(float)CLOCKS_PER_SEC);
     printf("Total elapsed time so far: %f\n",timeloop_total);
     if (it>100) {
         printf("Average elapsed time per time step: %f\n",timeloop_total/(float)(it-1));
     }
// check stability of the code, exit if unstable
     if(Usolidnorm > STABILITY_THRESHOLD || Usolidnorm < 0) {
         fprintf(stderr,"code became unstable and blew up\n");
         exit(1);
       }
   }

// big loop over all the global points (not elements) in the mesh to update
// the displacement and velocity vectors and clear the acceleration vector
 for (i=0;i<NGLOB;i++) {
   displ[i][0] += deltat*veloc[i][0] + deltatsqover2*accel[i][0];
   displ[i][1] += deltat*veloc[i][1] + deltatsqover2*accel[i][1];
   displ[i][2] += deltat*veloc[i][2] + deltatsqover2*accel[i][2];

   veloc[i][0] += deltatover2*accel[i][0];
   veloc[i][1] += deltatover2*accel[i][1];
   veloc[i][2] += deltatover2*accel[i][2];

   accel[i][0] = 0.;
   accel[i][1] = 0.;
   accel[i][2] = 0.;
 }

// big loop over all the elements in the mesh to localize data
// from the global vectors to the local mesh
// using indirect addressing (contained in array ibool)
// and then to compute the elemental contribution
// to the acceleration vector of each element of the finite-element mesh
 for (ispec=0;ispec<NSPEC;ispec++) {

   for (k=0;k<NGLLZ;k++) {
     for (j=0;j<NGLLY;j++) {
       for (i=0;i<NGLLX;i++) {
           iglob = ibool[ispec][k][j][i];
           dummyx_loc[k][j][i] = displ[iglob][0];
           dummyy_loc[k][j][i] = displ[iglob][1];
           dummyz_loc[k][j][i] = displ[iglob][2];
       }
     }
   }

// big loop over all the elements in the mesh to compute the elemental contribution
// to the acceleration vector of each element of the finite-element mesh

   for (k=0;k<NGLLZ;k++) {
     for (j=0;j<NGLLY;j++) {
       for (i=0;i<NGLLX;i++) {

         tempx1l = 0.;
         tempx2l = 0.;
         tempx3l = 0.;

         tempy1l = 0.;
         tempy2l = 0.;
         tempy3l = 0.;

         tempz1l = 0.;
         tempz2l = 0.;
         tempz3l = 0.;

         for (l=0;l<NGLLX;l++) {
           hp1 = hprime_xx[l][i];
           tempx1l += dummyx_loc[k][j][l]*hp1;
           tempy1l += dummyy_loc[k][j][l]*hp1;
           tempz1l += dummyz_loc[k][j][l]*hp1;

           hp2 = hprime_xx[l][j];
           tempx2l += dummyx_loc[k][l][i]*hp2;
           tempy2l += dummyy_loc[k][l][i]*hp2;
           tempz2l += dummyz_loc[k][l][i]*hp2;

           hp3 = hprime_xx[l][k];
           tempx3l += dummyx_loc[l][j][i]*hp3;
           tempy3l += dummyy_loc[l][j][i]*hp3;
           tempz3l += dummyz_loc[l][j][i]*hp3;
         }

// compute derivatives of ux, uy and uz with respect to x, y and z
         xixl = xix[ispec][k][j][i];
         xiyl = xiy[ispec][k][j][i];
         xizl = xiz[ispec][k][j][i];
         etaxl = etax[ispec][k][j][i];
         etayl = etay[ispec][k][j][i];
         etazl = etaz[ispec][k][j][i];
         gammaxl = gammax[ispec][k][j][i];
         gammayl = gammay[ispec][k][j][i];
         gammazl = gammaz[ispec][k][j][i];
         jacobianl = 1. / (xixl*(etayl*gammazl-etazl*gammayl)-xiyl*(etaxl*gammazl-etazl*gammaxl)+xizl*(etaxl*gammayl-etayl*gammaxl));

         duxdxl = xixl*tempx1l + etaxl*tempx2l + gammaxl*tempx3l;
         duxdyl = xiyl*tempx1l + etayl*tempx2l + gammayl*tempx3l;
         duxdzl = xizl*tempx1l + etazl*tempx2l + gammazl*tempx3l;

         duydxl = xixl*tempy1l + etaxl*tempy2l + gammaxl*tempy3l;
         duydyl = xiyl*tempy1l + etayl*tempy2l + gammayl*tempy3l;
         duydzl = xizl*tempy1l + etazl*tempy2l + gammazl*tempy3l;

         duzdxl = xixl*tempz1l + etaxl*tempz2l + gammaxl*tempz3l;
         duzdyl = xiyl*tempz1l + etayl*tempz2l + gammayl*tempz3l;
         duzdzl = xizl*tempz1l + etazl*tempz2l + gammazl*tempz3l;

// precompute some sums to save CPU time
         duxdxl_plus_duydyl = duxdxl + duydyl;
         duxdxl_plus_duzdzl = duxdxl + duzdzl;
         duydyl_plus_duzdzl = duydyl + duzdzl;
         duxdyl_plus_duydxl = duxdyl + duydxl;
         duzdxl_plus_duxdzl = duzdxl + duxdzl;
         duzdyl_plus_duydzl = duzdyl + duydzl;

// compute isotropic elements
         kappal = kappav[ispec][k][j][i];
         mul = muv[ispec][k][j][i];

         lambdalplus2mul = kappal + (4./3.) * mul;
         lambdal = lambdalplus2mul - 2.*mul;

// compute stress sigma
         sigma_xx = lambdalplus2mul*duxdxl + lambdal*duydyl_plus_duzdzl;
         sigma_yy = lambdalplus2mul*duydyl + lambdal*duxdxl_plus_duzdzl;
         sigma_zz = lambdalplus2mul*duzdzl + lambdal*duxdxl_plus_duydyl;

         sigma_xy = mul*duxdyl_plus_duydxl;
         sigma_xz = mul*duzdxl_plus_duxdzl;
         sigma_yz = mul*duzdyl_plus_duydzl;

// form dot product with test vector
     tempx1[k][j][i] = jacobianl * (sigma_xx*xixl + sigma_xy*xiyl + sigma_xz*xizl);
     tempy1[k][j][i] = jacobianl * (sigma_xy*xixl + sigma_yy*xiyl + sigma_yz*xizl);
     tempz1[k][j][i] = jacobianl * (sigma_xz*xixl + sigma_yz*xiyl + sigma_zz*xizl);

     tempx2[k][j][i] = jacobianl * (sigma_xx*etaxl + sigma_xy*etayl + sigma_xz*etazl);
     tempy2[k][j][i] = jacobianl * (sigma_xy*etaxl + sigma_yy*etayl + sigma_yz*etazl);
     tempz2[k][j][i] = jacobianl * (sigma_xz*etaxl + sigma_yz*etayl + sigma_zz*etazl);

     tempx3[k][j][i] = jacobianl * (sigma_xx*gammaxl + sigma_xy*gammayl + sigma_xz*gammazl);
     tempy3[k][j][i] = jacobianl * (sigma_xy*gammaxl + sigma_yy*gammayl + sigma_yz*gammazl);
     tempz3[k][j][i] = jacobianl * (sigma_xz*gammaxl + sigma_yz*gammayl + sigma_zz*gammazl);

         }
       }
     }

   for (k=0;k<NGLLZ;k++) {
     for (j=0;j<NGLLY;j++) {
       for (i=0;i<NGLLX;i++) {

         tempx1l = 0.;
         tempy1l = 0.;
         tempz1l = 0.;

         tempx2l = 0.;
         tempy2l = 0.;
         tempz2l = 0.;

         tempx3l = 0.;
         tempy3l = 0.;
         tempz3l = 0.;

         for (l=0;l<NGLLX;l++) {
           fac1 = hprimewgll_xx[i][l];
           tempx1l += tempx1[k][j][l]*fac1;
           tempy1l += tempy1[k][j][l]*fac1;
           tempz1l += tempz1[k][j][l]*fac1;

           fac2 = hprimewgll_xx[j][l];
           tempx2l += tempx2[k][l][i]*fac2;
           tempy2l += tempy2[k][l][i]*fac2;
           tempz2l += tempz2[k][l][i]*fac2;

           fac3 = hprimewgll_xx[k][l];
           tempx3l += tempx3[l][j][i]*fac3;
           tempy3l += tempy3[l][j][i]*fac3;
           tempz3l += tempz3[l][j][i]*fac3;
         }

         fac1 = wgllwgll_yz[k][j];
         fac2 = wgllwgll_xz[k][i];
         fac3 = wgllwgll_xy[j][i];

// sum contributions from each element to the global mesh using indirect addressing
         iglob = ibool[ispec][k][j][i];
         accel[iglob][0] -= (fac1*tempx1l + fac2*tempx2l + fac3*tempx3l);
         accel[iglob][1] -= (fac1*tempy1l + fac2*tempy2l + fac3*tempy3l);
         accel[iglob][2] -= (fac1*tempz1l + fac2*tempz2l + fac3*tempz3l);

       }
     }
   }

 }   // end of main loop on all the elements

// big loop over all the global points (not elements) in the mesh to update
// the acceleration and velocity vectors
 for (i=0;i<NGLOB;i++) {
   accel[i][0] *= rmass_inverse[i];
   accel[i][1] *= rmass_inverse[i];
   accel[i][2] *= rmass_inverse[i];
 }

// add the earthquake source at a given grid point
// this is negligible and is intrinsically serial because it is done by only
// one grid point out of several millions typically
// we subtract one to the element number of the source because arrays start at 0 in C
// compute current time
 time = (it-1)*deltat;
 accel[ibool[NSPEC_SOURCE-1][1][1][1]][2] += 1.e4 * (1.-2.*a*(time-t0)*(time-t0)) * exp(-a*(time-t0)*(time-t0)) / rho;

 for (i=0;i<NGLOB;i++) {
   veloc[i][0] += deltatover2*accel[i][0];
   veloc[i][1] += deltatover2*accel[i][1];
   veloc[i][2] += deltatover2*accel[i][2];
 }

// record a seismogram to check that the simulation went well
// we subtract one to the element number of the receiver because arrays start at 0 in C
   seismogram[it-1] = displ[ibool[NSPEC_STATION-1][1][1][1]][2];

 } // end of the serial time loop

// save the seismogram at the end of the run
 if((IIN = fopen("seismogram_C_double.txt","w")) == NULL) {
         fprintf(stderr,"Cannot open file seismogram_C_double.txt, exiting...\n");
         exit(1);
       }
 for (it=0;it<NSTEP;it++)
 {  fprintf(IIN,"%e     %e\n",(float)(it*deltat),(float)seismogram[it]);
 }
 fclose(IIN);

 printf("\nEnd of the program\n\n");

 }

