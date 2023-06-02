#!/usr/bin/env python
#
# plot_binary_seismo_file_with_obspy.py
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

#########################################################################################

# sets default custom_real type
# defines float 'f' or double precision 'd' for binary values
custom_type = 'f'

#########################################################################################

def read_binary_seismo_file(filename,num_receivers,num_samples,NDIM,verbose=False,show=False):
    """
    reads binary seismo file

    in fortran code, (gather) file must have been written with access='direct'
    and a record length of recl=(CUSTOM_REAL * num_receivers * num_samples)
    """
    global custom_type

    # user output
    if verbose:
        print("binary file: ",filename)
        print("number of receivers     : ",num_receivers)
        print("number of trace samples : ",num_samples)
        print("")

    # checks if file exists
    if not os.path.isfile(filename):
        print("file %s does not exist, please check..." % filename)
        sys.exit(1)

    if custom_type == 'f':
        # float (4 bytes) for single precision
        bintype = 'float32'
        binlength = 4
    else:
        # double precision
        bintype = 'float64'
        binlength = 8

    # empty array
    data_all = []

    # open file
    with open(filename,'rb') as f:
        for idim in range(NDIM):
            #debug
            #print("read: dim ",idim," - NDIM ",NDIM)

            # reads data slice
            # positions file pointer
            f.seek(idim * binlength * num_receivers * num_samples)
            # reads data slice
            array = np.fromfile(f,dtype=bintype,count=num_receivers * num_samples)

            # checks
            if array.size == 0:
                print("Error: no more data in file: idim ",idim)
                print("Please check if NDIM, number of stations and samples are correct for given binary file")
                sys.exit(1)

            # reshapes array
            array = np.reshape(array,(num_samples,num_receivers))

            #debug
            #print("read: dim ",idim," - array: ",array[0])

            # plots traces
            if show:
                for irec in range(num_receivers):
                    plt.title("station {}".format(irec))
                    plt.plot(np.arange(num_samples),array[:,irec])
                    plt.show()

            # adds to data
            data_all.append(array)

    # converts to numpy array
    data = np.array(data_all)

    return data

def compare_data_misfit(st_ref,st_compare):
    # vertical trace
    max_err = 0.0
    max_rel_err = 0.0
    max_name = ""
    total_err = 0.0

    for i,tr in enumerate(st_ref):
        #info
        #print("Trace "+str(i+1)+":")
        #print(tr)
        #print(tr.stats)

        # trace
        name = "{}.{}.{}".format(tr.stats.network,tr.stats.station,tr.stats.channel)
        title = "Trace " + str(i+1) + " : " + name

        # reference trace
        dat0 = tr.data

        # comparison trace
        tr1 = st_compare[i]
        dat1 = tr1.data

        # cuts common length
        length = min(len(dat0),len(dat1))
        if length <= 1: continue

        # length warning
        if len(dat0) != len(dat1):
          print("** warning: mismatch of trace length in both files = %d / %d" %(len(dat0),len(dat1)))

        # least square test
        err = np.linalg.norm(dat0 - dat1)

        # relative error in percent
        norm0 = np.linalg.norm(dat0)
        if norm0 > 0.0:
            rel_err = err / norm0 * 100.0
        else:
            rel_err = err

        # stats
        if rel_err > max_rel_err:
            max_rel_err = rel_err
            max_err = err
            max_name = name

        # total error
        total_err += err

        # print results to screen
        if i < 5 or i >= len(st_ref) - 5:
            print("  %-30s| %13.5e    (%6.2f %%)" % (name, err, rel_err))
        elif i == 6:
            print("  ..")
        else:
            continue

    print("")
    print("  %-30s| %13.5e" % ("total error",total_err))
    print("")
    print("  maximum error:\n  %-30s| %13.5e    (%6.2f %%)" % (max_name,max_err,max_rel_err))
    print("")


def plot_binary_seismo_file(file,NSTA,Nt,NDIM,DT,cfile1,cfile2,normalize_data=False,show_plot=False,show_single_trace=""):

    # check if cfile1 available
    if len(cfile1) > 0:
        compare_file_1 = True
    else:
        compare_file_1 = False

    # check if cfile2 available
    if len(cfile2) > 0:
        compare_file_2 = True
    else:
        compare_file_2 = False

    # user output
    print("")
    print("getting station information...")
    print("")
    print("file : ",file)
    print("number of stations (NSTA)   : ",NSTA)
    print("number of time samples (Nt) : ",Nt)
    print("number of components (NDIM) : ",NDIM)
    print("DT         : ",DT)
    if compare_file_1: print("file1 to compare : ",cfile1)
    if compare_file_2: print("file2 to compare : ",cfile2)
    if normalize_data: print("normalize amplitudes for comparing traces: ",normalize_data)
    print("")
    print("show plots : ",show_plot)
    if show_single_trace:
        print("show single trace : ",show_single_trace)
    print("")

    # helper function
    data = read_binary_seismo_file(file,NSTA,Nt,NDIM,verbose=True)

    print("read: data shape = ",data.shape)
    print("      data size  = ",data.size)
    print("")

    if data.size == 0:
        print("no data")
        sys.exit(1)

    # obspy stream
    st = obspy.Stream()

    # file comparison
    if compare_file_1:
        data1 = read_binary_seismo_file(cfile1,NSTA,Nt,NDIM,verbose=True)
        print("read: data1 shape = ",data1.shape)
        print("      data1 size  = ",data1.size)
        print("")
        # new stream
        st1 = obspy.Stream()

    # file comparison
    if compare_file_2:
        data2 = read_binary_seismo_file(cfile2,NSTA,Nt,NDIM,verbose=True)
        print("read: data2 shape = ",data2.shape)
        print("      data2 size  = ",data2.size)
        print("")
        # new stream
        st2 = obspy.Stream()

    # sets traces
    for idim in range(NDIM):
        data_comp = data[idim]

        #debug
        #print("data component: ",data_comp.shape)

        for irec in range(NSTA):
            #print("receiver: ",irec+1)
            # sets name
            sta = "S{:02d}".format(irec+1)
            net = "DB"
            cha = [ "MUX", "MUY", "MUZ" ][idim]
            # artificial distance to plot stream as 'section'
            dist = irec * 1.0

            # gets trace data
            dat = data_comp[:,irec]

            # normalize for comparison
            if normalize_data:
                maxval = np.abs(dat).max()
                if maxval > 0.0: dat = dat / maxval

            # creates new trace
            tr = obspy.Trace()

            tr.stats.station = sta
            tr.stats.network = net
            tr.stats.channel = cha
            tr.stats.distance = dist
            tr.stats.delta = DT
            tr.data = dat

            st.append(tr)

            # file comparisons
            if compare_file_1 > 0:
                dat = data1[idim,:,irec]

                # normalize for comparison
                if normalize_data:
                    maxval = np.abs(dat).max()
                    if maxval > 0.0: dat = dat / maxval

                tr1 = obspy.Trace()
                tr1.stats.station = sta
                tr1.stats.network = net
                tr1.stats.channel = cha + "-COMPARE1"
                tr1.stats.distance = dist
                tr1.stats.delta = DT
                tr1.data = dat

                st1.append(tr1)

            if compare_file_2 > 0:
                dat = data2[idim,:,irec]

                # normalize for comparison
                if normalize_data:
                    maxval = np.abs(dat).max()
                    if maxval > 0.0: dat = dat / maxval

                tr2 = obspy.Trace()
                tr2.stats.station = sta
                tr2.stats.network = net
                tr2.stats.channel = cha + "-COMPARE2"
                tr2.stats.distance = dist
                tr2.stats.delta = DT
                tr2.data = dat

                st2.append(tr2)

    print("traces: ")
    print("  number of traces = ",len(st))
    print("")
    #print(st)
    #print("")

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

    # data comparison
    if compare_file_1:
        print("trace comparison 1:")
        print("  file1: ",cfile1)
        # data misfits
        compare_data_misfit(st,st1)

    if compare_file_2:
        print("trace comparison 2:")
        print("  file2: ",cfile2)
        # data misfits
        compare_data_misfit(st,st2)

    # plotting
    if show_plot:
        print("plotting traces...")

        # plots all traces
        if not show_single_trace:
            st.plot(type='default',size=(1800, 900))
            #st.plot(type='section',size=(1800, 900))

        # time axis for plotting (in s)
        t_d = np.linspace(0,duration,npts)

        # vertical trace
        for i,tr in enumerate(st):
            #info
            #print("Trace "+str(i+1)+":")
            #print(tr)
            #print(tr.stats)

            # plotting
            name = "{}.{}.{}".format(tr.stats.network,tr.stats.station,tr.stats.channel)

            if show_single_trace and name not in show_single_trace:
                continue

            # main title
            title = "Trace " + str(i+1) + " : " + name
            plt.suptitle(title)

            # sub-title
            title = "{}".format(file)
            if compare_file_1: title += "\ncomparison 1: {} (red)".format(cfile1)
            if compare_file_2: title += "\ncomparison 2: {} (blue)".format(cfile2)
            plt.title(title,fontsize=8)

            # plot trace
            plt.plot(t_d, tr.data, color='black', linewidth=1.0,label=name)

            # comparison
            if compare_file_1:
                tr1 = st1[i]
                name1 = "{}.{}.{}".format(tr1.stats.network,tr1.stats.station,tr1.stats.channel)
                plt.plot(t_d, tr1.data, color='red', linewidth=1.0,label=name1)

            # comparison
            if compare_file_2:
                tr2 = st2[i]
                name2 = "{}.{}.{}".format(tr2.stats.network,tr2.stats.station,tr2.stats.channel)
                plt.plot(t_d, tr2.data, color='blue', linestyle='dashed', linewidth=1.0,label=name2)

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
    print("usage: ./plot_binary_seismo_file_with_obspy.py file NSTA Nt NDIM [DT]")
    print("                                               [compare1=file1] [compare2=file2] [normalize]")
    print("                                               [show] [show_single_trace=name]")
    print(" with")
    print("   file   - binary seismogram file, e.g. RAW_DATA/event_01.bin")
    print("   NSTA   - number of stations (per component)")
    print("   Nt     - number of time samples (per trace)")
    print("   NDIM   - number of components")
    print("   DT     - (optional) time step size")
    print("   compare1=file1 - (optional) to compare against another binary file, e.g. compare1=RAW_DATA/event_01.bin_initial_syn")
    print("   compare2=file2 - (optional) to compare against another binary file, e.g. compare2=RAW_DATA/event_01.bin_iter0000_syn")
    print("   normalize      - (optional) normalize amplitudes for comparison")
    print("   show                    - (optional) plot waveforms")
    print("   show_single_trace=name  - (optional) plot only single trace matching name, e.g., show_single_trace=DB.S485.MUX")
    sys.exit(1)


if __name__ == "__main__":
    # gets arguments
    NSTA = 0
    Nt = 0
    NDIM = 0
    DT = 1.0
    cfile1 = ""
    cfile2 = ""
    normalize_data = False
    do_plot_waveforms = False
    show_single_trace = ""

    # reads arguments
    if len(sys.argv) < 5:
        usage()
    else:
        file = sys.argv[1]
        NSTA = int(sys.argv[2])
        Nt = int(sys.argv[3])
        NDIM = int(sys.argv[4])

    if len(sys.argv) >= 6:
        DT = float(sys.argv[5])

    # optional arguments
    for arg in sys.argv:
        if "--help" in arg:
            usage()

        elif "compare1" in arg:
            cfile1 = arg.split('=')[1]

        elif "compare2" in arg:
            cfile2 = arg.split('=')[1]

        elif "normalize" in arg:
            normalize_data = True

        elif "show" in arg:
            do_plot_waveforms = True
            if "show_single_trace" in arg:
                show_single_trace = arg.split('=')[1]

    plot_binary_seismo_file(file,NSTA,Nt,NDIM,DT,cfile1,cfile2,normalize_data,do_plot_waveforms,show_single_trace)




