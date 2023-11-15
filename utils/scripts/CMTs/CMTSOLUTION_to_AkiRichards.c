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

/*  cmt2aki.c  : Converts Harvard-CMT moment moment tensor elements to Aki&Richard 1980 convention
  by Onur TAN
  ITU, Dept. of Geophysics, Istanbul, Turkey, 10 Jan 2004
*/

// Modified by Dimitri Komatitsch, University of Pau, France, September 2007

/* Harvard CMTSOLUTION format is for instance
Mrr:      -7.600000e+27
Mtt:       7.700000e+27
Mpp:      -2.000000e+26
Mrt:      -2.500000e+28
Mrp:       4.000000e+26
Mtp:      -2.500000e+27
*/

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

int main ( int argc, char* argv[] )
{
float Mrr, Mtt, Mpp, Mrt, Mrp, Mtp;
float Mxx, Myy, Mzz, Mxy, Mxz, Myz;

if ( argc != 7 )  {
printf("\nHRV-CMT - Aki&Richards1980 Moment Tensor Element Converter (by Onur TAN)\n");
printf("Usage : cmt2aki  Mrr Mtt Mpp Mrt Mrp Mtp\n") ;
printf("Input Harvard CMTSOLUTION:  Mrr Mtt Mpp Mrt Mrp Mtp\n") ;
printf("Output Aki&Richards1980:  Mxx  Myy  Mzz  Mxy  Mxz  Myz \n\n");
exit(1);
}

Mrr  = atof ( argv[1] );
Mtt  = atof ( argv[2] );
Mpp  = atof ( argv[3] );
Mrt  = atof ( argv[4] );
Mrp  = atof ( argv[5] );
Mtp  = atof ( argv[6] );

Mxx  = Mtt ;
Myy  = Mpp ;
Mzz  = Mrr ;
Mxy  = Mtp * -1. ;
Mxz  = Mrt ;
Myz  = Mrp * -1. ;

printf("Output Aki&Richards1980:  Mxx  Myy  Mzz  Mxy  Mxz  Myz \n");
printf("%9.5f %9.5f %9.5f %9.5f %9.5f %9.5f\n",Mxx, Myy, Mzz, Mxy, Mxz, Myz);

}

