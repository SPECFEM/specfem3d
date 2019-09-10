#!/usr/bin/env python
#
# based on C-code routine strike_dip_rake_to_CMTSOLUTION.c
#
#/*---------------------------------------------------------------------------
# *
# *      Copyright (c) 2000-2004 by Onur TAN
# *      See COPYING file for copying and redistribution conditions.
# *
# *      This program is free software; you can redistribute it and/or modify
# *      it under the terms of the GNU General Public License as published by
# *      the Free Software Foundation; version 2 of the License.
# *
# *      This program is distributed in the hope that it will be useful,
# *      but WITHOUT ANY WARRANTY; without even the implied warranty of
# *      MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# *      GNU General Public License for more details.
# *
# *      Contact info:   Onur TAN,
# *      Istanbul Technical University, Faculty of Mines
# *      Department of Geophysics, Maslak, Istanbul-TURKEY
# *
# *      New address (in 2017):  onur.tan AT mam.gov.tr    https://en.wikipedia.org/wiki/T%C3%9CB%C4%B0TAK_Marmara_Research_Center
# *
# *--------------------------------------------------------------------------*/
#
#
#/*  dc2mt.c  : Calculates moment tensor elements from strike/dip/rake
#  by Onur TAN
#  ITU, Dept. of Geophysics, Istanbul, Turkey, 02 Jan 2004
#               21 Apr 2004 bug fix for Mxy
#*/
#
#// Modified by Dimitri Komatitsch, University of Pau, France, September 2007
#// compile with "  gcc -o strike_dip_rake_to_CMTSOLUTION strike_dip_rake_to_CMTSOLUTION.c -lm " to include the math library
#
#// Message from Onur TAN to Dimitri Komatitsch, about the reverse conversion program, Feb 2016:
#// I am sorry, I have no reverse conversion program.
#// But I found these web sites for Matlab codes:
#// https://github.com/g2e/seizmo/blob/master/cmt/sdr2mt.m
#// https://github.com/g2e/seizmo/blob/master/cmt/mt2sdr.m
#// Note from Dimitri: I have downloaded them and put them in file seizmo_master_containing_Matlab_codes_for_mt2dc_conversion.zip in this directory.
#// Note from Elliott Sales de Andrade: ObsPy https://github.com/obspy also provides conversion for moment tensors from/to several different bases and strike/dip/rake.
#// Note from Dimitri: we also have a Matlab code from Carl Tape that does that now, in this directory.
#
##include <stdio.h>
##include <string.h>
##include <stdlib.h>
##include <math.h>
#
from __future__ import print_function

import sys
import os
from math import sin,cos


def convert_strike_dip_rake_to_CMTSOLUTION(strike,dip,rake):
    """
    converts given strike/dip/rake to moment tensor
    """
    S = strike
    D = dip
    R = rake

    # PI / 180 to convert degrees to radians
    d2r =  0.017453293

    print("Strike    = %9.5f degrees" % S)
    print("Dip       = %9.5f degrees" % D)
    print("Rake/Slip = %9.5f degrees" % R)
    print("")

    # convert to radians
    S *= d2r
    D *= d2r
    R *= d2r

    # Aki & Richards
    Mxx = -1.0 * ( sin(D) * cos(R) * sin (2*S) + sin(2*D) * sin(R) * sin(S)*sin(S) )
    Myy =        ( sin(D) * cos(R) * sin (2*S) - sin(2*D) * sin(R) * cos(S)*cos(S) )
    Mzz = -1.0 * ( Mxx + Myy)
    Mxy =        ( sin(D) * cos(R) * cos (2*S) + 0.5 * sin(2*D) * sin(R) * sin(2*S) )
    Mxz = -1.0 * ( cos(D) * cos(R) * cos (S)   + cos(2*D) * sin(R) * sin(S) )
    Myz = -1.0 * ( cos(D) * cos(R) * sin (S)   - cos(2*D) * sin(R) * cos(S) )

    print("Output Aki&Richards1980:  Mxx  Myy  Mzz  Mxy  Mxz  Myz")
    print("%9.5f %9.5f %9.5f %9.5f %9.5f %9.5f\n" %(Mxx, Myy, Mzz, Mxy, Mxz, Myz))
    print("")

    # also convert to Harvard CMTSOLUTION format
    Mtt = Mxx
    Mpp = Myy
    Mrr = Mzz
    Mtp = Mxy * -1.0  # sign change
    Mrt = Mxz
    Mrp = Myz * -1.0  # sign change

    # Harvard CMTSOLUTION format is for instance
    # Mrr:      -7.600000e+27
    # Mtt:       7.700000e+27
    # Mpp:      -2.000000e+26
    # Mrt:      -2.500000e+28
    # Mrp:       4.000000e+26
    # Mtp:      -2.500000e+27

    print("Output Harvard CMTSOLUTION:  Mrr Mtt Mpp Mrt Mrp Mtp")
    print("Mrr: %9.5f" % Mrr)
    print("Mtt: %9.5f" % Mtt)
    print("Mpp: %9.5f" % Mpp)
    print("Mrt: %9.5f" % Mrt)
    print("Mrp: %9.5f" % Mrp)
    print("Mtp: %9.5f" % Mtp)
    print("")

    return


if __name__ == '__main__':
    # gets arguments
    if len(sys.argv) != 4:
        print("\nDouble Couple - Moment Tensor Converter (original dc2mt routine by Onur TAN)\n")
        print("Usage : ./strike_dip_rake_to_CMTSOLUTION.py Strike Dip Rake (in degrees)" )
        print("  outputs Aki & Richards 1980:  Mxx  Myy  Mzz  Mxy  Mxz  Myz")
        print("          Harvard CMTSOLUTION:  Mrr  Mtt  Mpp  Mrt  Mrp  Mtp")
        sys.exit(1)
    else:
        strike = float(sys.argv[1])
        dip = float(sys.argv[2])
        rake = float(sys.argv[3])

    convert_strike_dip_rake_to_CMTSOLUTION(strike,dip,rake)


