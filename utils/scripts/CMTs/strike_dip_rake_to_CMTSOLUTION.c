/*---------------------------------------------------------------------------
 *
 *      Copyright (c) 2000-2004 by Onur TAN
 *      See COPYING file for copying and redistribution conditions.
 *
 *      This program is free software; you can redistribute it and/or modify
 *      it under the terms of the GNU General Public License as published by
 *      the Free Software Foundation; version 2 of the License.
 *
 *      This program is distributed in the hope that it will be useful,
 *      but WITHOUT ANY WARRANTY; without even the implied warranty of
 *      MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *      GNU General Public License for more details.
 *
 *      Contact info:   Onur TAN,
 *      Istanbul Technical University, Faculty of Mines
 *      Department of Geophysics, Maslak, Istanbul-TURKEY
 *
 *      New address (in 2017):  onur.tan AT mam.gov.tr    https://en.wikipedia.org/wiki/T%C3%9CB%C4%B0TAK_Marmara_Research_Center
 *
 *--------------------------------------------------------------------------*/

/*  dc2mt.c  : Calculates moment tensor elements from strike/dip/rake
  by Onur TAN
  ITU, Dept. of Geophysics, Istanbul, Turkey, 02 Jan 2004
               21 Apr 2004 bug fix for Mxy
*/

// Modified by Dimitri Komatitsch, University of Pau, France, September 2007
// compile with "  gcc -o strike_dip_rake_to_CMTSOLUTION strike_dip_rake_to_CMTSOLUTION.c -lm " to include the math library

// Message from Onur TAN to Dimitri Komatitsch, about the reverse conversion program, Feb 2016:
// I am sorry, I have no reverse conversion program.
// But I found these web sites for Matlab codes:
// https://github.com/g2e/seizmo/blob/master/cmt/sdr2mt.m
// https://github.com/g2e/seizmo/blob/master/cmt/mt2sdr.m
// Note from Dimitri: I have downloaded them and put them in file seizmo_master_containing_Matlab_codes_for_mt2dc_conversion.zip in this directory.
// Note from Elliott Sales de Andrade: ObsPy https://github.com/obspy also provides conversion for moment tensors from/to several different bases and strike/dip/rake.
// Note from Dimitri: we also have a Matlab code from Carl Tape that does that now, in this directory.

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

int main ( int argc, char* argv[] )
{
float S,D,R,d2r;
float Mxx, Myy, Mzz, Mxy, Mxz, Myz;
float Mrr, Mtt, Mpp, Mrt, Mrp, Mtp;

if ( argc != 4 )  {
printf("\nDouble Couple - Moment Tensor Converter (by Onur TAN)\n");
printf("Usage : dc2mt Strike Dip Rake (in degrees)\n" );
printf("Output Aki&Richards1980:  Mxx  Myy  Mzz  Mxy  Mxz  Myz \n\n");
exit(1);
}

// PI / 180 to convert degrees to radians
d2r =  0.017453293;

S  = atof ( argv[1] );
D  = atof ( argv[2] );
R  = atof ( argv[3] );

printf("\nStrike = %9.5f degrees\n",S);
printf("Dip = %9.5f degrees\n",D);
printf("Rake/Slip = %9.5f degrees\n",R);

// convert to radians
S *= d2r;
D *= d2r;
R *= d2r;

Mxx = -1.0 * ( sin(D) * cos(R) * sin (2*S) + sin(2*D) * sin(R) * sin(S)*sin(S) ); /* Mxx */
Myy =        ( sin(D) * cos(R) * sin (2*S) - sin(2*D) * sin(R) * cos(S)*cos(S) ); /* Myy */
Mzz = -1.0 * ( Mxx + Myy);                /* Mzz */
Mxy =        ( sin(D) * cos(R) * cos (2*S) + 0.5 * sin(2*D) * sin(R) * sin(2*S) );  /* Mxy */
Mxz = -1.0 * ( cos(D) * cos(R) * cos (S)   + cos(2*D) * sin(R) * sin(S) );    /* Mxz */
Myz = -1.0 * ( cos(D) * cos(R) * sin (S)   - cos(2*D) * sin(R) * cos(S) );    /* Myz */

printf("\nOutput Aki&Richards1980:  Mxx  Myy  Mzz  Mxy  Mxz  Myz \n");
printf("%9.5f %9.5f %9.5f %9.5f %9.5f %9.5f\n",Mxx, Myy, Mzz, Mxy, Mxz, Myz);

// also convert to Harvard CMTSOLUTION format
Mtt = Mxx;
Mpp = Myy;
Mrr = Mzz;
Mtp = Mxy * -1. ;
Mrt = Mxz;
Mrp = Myz * -1. ;

/* Harvard CMTSOLUTION format is for instance
Mrr:      -7.600000e+27
Mtt:       7.700000e+27
Mpp:      -2.000000e+26
Mrt:      -2.500000e+28
Mrp:       4.000000e+26
Mtp:      -2.500000e+27
*/

printf("\nOutput Harvard CMTSOLUTION:  Mrr Mtt Mpp Mrt Mrp Mtp\n");
printf("Mrr: %9.5f\n",Mrr);
printf("Mtt: %9.5f\n",Mtt);
printf("Mpp: %9.5f\n",Mpp);
printf("Mrt: %9.5f\n",Mrt);
printf("Mrp: %9.5f\n",Mrp);
printf("Mtp: %9.5f\n\n",Mtp);

}

