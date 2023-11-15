#!/usr/bin/env python
#
# plot_SU_file_with_obspy.py
#
from __future__ import print_function
import os,sys

## matplotlib
#import matplotlib as mpl
#print("matplotlib version: ",mpl.__version__)

import matplotlib.pyplot as plt

plt.style.use('ggplot')
plt.rcParams['figure.figsize'] = 20, 5
plt.rcParams['lines.linewidth'] = 0.5

## numpy
import numpy as np
print("numpy version: ",np.__version__)

## do not show scipy warnings
#import warnings
#warnings.filterwarnings('ignore')

## obspy
import obspy
print("obspy version: ",obspy.__version__)


## Seismic Unix
#
# needed on mcmillan
#suoldtonew < $file | suxwigb #perc=99 &

# reads header info:
#   ns = number of samples
#   dt = time step
#   sx,sy = source coordinates
#   gx,gy = receiver coordinates
#sugethw <$sufile key=dt,ns,sx,sy,gx,gy verbose=1
#
# removes headers
#sustrip <$sufile > $sufile.raw
#
# wiggle plot
#xwigb <$sufile n1=$NSTEP n2=$ntraces &
#
# image plot
#ximage <$sufile n1=$NSTEP n2=$ntraces &
#
# ps image
#suoldtonew <$sufile | supsimage label1="Time (s)" label2="Offset (m)" title="$sufile"  perc=99 labelsize=22 verbose=1 > tmp.eps
#
# convert to png
#ps2pdf -r300 -sDEVICE=png16m -sOutputFile=$sufile.png tmp.eps



def plot_SU_file(file,show_plot):
    print("")
    print("getting station information...")
    print("")

    st = obspy.read(file)

    print("traces: ")
    print(st)

    print("")
    print("number of traces: ",len(st))
    print("")

    npts = 0
    duration = 0.0

    # trace statistics
    if len(st) > 0:
        tr = st[0]
        #print(tr.stats)
        npts = tr.stats.npts      # number of samples
        dt = tr.stats.delta       # time step
        duration = npts * dt

        freq_Ny = 1.0/(2.0 * dt)  # Nyquist frequency
        freq_min = 1./duration    # minimal possible frequency for length of trace

        print("trace info:")
        print("  number of samples = ",npts)
        print("  time step         = ",dt,"(s)")
        print("  duration          = ",duration,"(s)")
        print("  Nyquist frequency = ",freq_Ny)
        print("  minimum frequency = ",freq_min)
        print("")
    else:
        print("no trace information, exiting...")
        sys.exit(1)

    # plotting
    if show_plot:
        # plots all traces
        st.plot()

        # time axis for plotting (in s)
        t_d = np.linspace(0,duration,npts)

        # vertical trace
        for i,tr in enumerate(st):
            print("Trace "+str(i+1)+":")
            print(tr)
            #print(tr.stats)
            # plotting
            plt.title("Trace " + str(i+1))
            # vertical component
            plt.plot(t_d, tr.data, color='black', linewidth=1.0,label="raw data trace "+str(i+1))
            plt.legend()
            plt.xlabel("Time (s)")
            plt.ylabel("Counts")
            plt.show()

        # Wait for the user input to terminate the program
        input("Press any key to terminate the program")

    print("")
    print("all done")
    print("")


def usage():
    print("usage: ./plot_SU_file_with_obspy.py su_file [show]")
    print(" with")
    print("   su_file - seismogram file e.g. OUTPUT_FILES/0_dz_SU")
    print("   show   - (optional) plot waveforms")
    sys.exit(1)


if __name__ == "__main__":
    # gets arguments
    do_plot_waveforms = False
    if len(sys.argv) == 2:
        file = sys.argv[1]
    elif len(sys.argv) == 3:
        file = sys.argv[1]
        if sys.argv[2] == "show": do_plot_waveforms = True
    else:
        usage()
        sys.exit(1)

    plot_SU_file(file,do_plot_waveforms)




