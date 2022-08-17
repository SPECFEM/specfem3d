#!/usr/bin/env python
#
# creates adjoint sources f^adj = (syn - data)
#
#########################################
from __future__ import print_function

import sys
import numpy as np
import array
import os

from helper_functions import helper

#########################################

## globals
# takes 2nd-derivative of pressure for adjoint source
use_derivative_of_pressure = True

#########################################



def adj_seismogram(filename_syn,filename_dat):
    """
    creates adjoint seismograms
    """
    global use_derivative_of_pressure

    print("")
    print("creating adjoint seismograms:")
    print("  input syn : ",filename_syn)
    print("  input data: ",filename_dat)
    print("")

    # helper functions
    hf = helper()

    # basename
    name = os.path.basename(filename_syn)

    # reads in seismograms
    syn = hf.read_SU_file(filename_syn)
    dat = hf.read_SU_file(filename_dat)

    # adjoint source f^adj = (s - d)
    adj = syn - dat

    # misfit values
    print("misfit:")
    diff_max = np.abs(adj).max()
    print("  maximum waveform difference (syn - dat) = ",diff_max)

    # sampling rate given in microseconds
    DT = hf.sampling_DT * 1.e-6
    print("  trace time steps: DT = ",DT,"(s)")

    # total misfit
    total_misfit = 0.0
    for irec in range(len(adj)):
        # single receiver trace
        adj_trace = adj[irec]
        # inner product
        total_misfit += np.sum(adj_trace * adj_trace) * DT

    print("")
    print("  total misfit: sum(s - d)^2 = {:e}".format(total_misfit))
    print("")

    # number of receivers/traces (SU files have only single components)
    num_receivers = len(adj)

    print("adjoint source:")
    print("  number of traces = ",num_receivers)

    # checks
    if num_receivers == 0:
        print("Did find no receivers or traces, please check...")
        sys.exit(1)

    # for acoustic FWI, L2 adjoint source is the second derivative of pressure difference
    # (e.g., see Peter et al. 2011, GJI, eq. (A8))
    if "_p_SU" in name:
        # pressure output
        # note: pressure in fluid is defined as p = - \partial_t^2 phi
        #       thus, if potential phi is chosen as output, there is a minus sign and time derivative difference.
        #
        #       assuming a pressure output for syn and dat, the adjoint source expression is given by (A8) in Peter et al. (2011)
        #       note the negative sign in the definition.
        #       the adjoint source for pressure is: f^adj = - \partial_t^2 p_syn - \partial_t^2 p_obs
        #                                                 = - \partial_t^2 ( p_syn - p_obs )
        #
        if use_derivative_of_pressure:
            print("  creating adjoint sources for pressure (taking second-derivative of pressure differences)...")
            # takes second-derivative
            adj_new = adj.copy()
            fac = 1.0 / DT**2
            for irec in range(num_receivers):
                # single receiver trace
                adj_trace = adj[irec]

                # 2nd-order scheme
                #for i in range(1,len(adj_trace)-1):
                #    # central finite difference (2nd-order scheme)
                #    val = (adj_trace[i+1] - 2.0 * adj_trace[i] + adj_trace[i-1]) * fac
                #    # adding negative sign
                #    adj_new[irec][i] = - val

                # 4th-order scheme
                for i in range(2,len(adj_trace)-2):
                    # central finite difference (4th-order scheme)
                    val = ( - 1.0/12.0 * adj_trace[i+2] + 4.0/3.0 * adj_trace[i+1] - 5.0/2.0 * adj_trace[i]
                            + 4.0/3.0 * adj_trace[i-1] - 1.0/12.0 * adj_trace[i-2] ) * fac
                    # adding negative sign
                    adj_new[irec][i] = - val

            # saves as adjoint source
            adj = adj_new.copy()

    # statistics
    amp_max = np.abs(adj).max()
    print("  maximum amplitude |f^adj| = ",amp_max)

    print("")

    # SEM output directory
    os.system("mkdir -p SEM/")

    # adjoint trace name
    adjname = name + ".adj"

    filename_adj = "SEM/" + adjname
    hf.write_SU_file(adj,filename_adj)

    # user output
    print("  receivers: ",len(adj))
    print("  adjoint sources written to: ",filename_adj)
    print("")

#
#------------------------------------------------------------------------------------------
#

def usage():
    print("Usage: ./adj_seismogram.py synthetics-SU-file data-SU-file")
    print("")
    sys.exit(1)


#
#------------------------------------------------------------------------------------------
#

if __name__ == '__main__':
    # gets arguments
    if len(sys.argv) < 3:
        usage()

    filename_syn = sys.argv[1]
    filename_dat = sys.argv[2]

    adj_seismogram(filename_syn,filename_dat)
