#!/usr/bin/env python
#
# reads in an HDF5-format seismogram file
# and outputs infos about content and checks validity
#
from __future__ import print_function

import sys,os

try:
    import h5py
except:
    print("Error importing h5py, check if h5py module is installed")
    print("")
    # python version
    print("python version:")
    print(sys.version)
    print("")
    # import module paths
    print("module paths:")
    for path in sys.path:
        print(path)
    print("")
    sys.exit(1)

try:
    import numpy as np
except:
    print("Error importing numpy, check if numpy module is installed")
    print("")
    # python version
    print("python version:")
    print(sys.version)
    print("")
    # import module paths
    print("module paths:")
    for path in sys.path:
        print(path)
    print("")
    sys.exit(1)

# visualization
import matplotlib.pyplot as plt

######################################################################
## USER PARAMETERS

# compare traces with reference solution
# works only with REF_SEIS/ folder in EXAMPLES/homogeneous_halfspace/
do_compare_ref = False

######################################################################

# read and plot the seismograms in hdf5 format
def plot_HDF5_data_info(file_in,show_plot,convert_waveforms):
    print("")
    print("reading file: ",file_in)
    print("")

    # open
    with h5py.File(file_in, "r") as f:
        # contents
        # coords                   Dataset {3, 4}
        # displ                    Dataset {4, 5000, 3}
        # network                  Dataset {4}
        # station                  Dataset {4}
        # time                     Dataset {5000}
        keys = list(f.keys())

        # file info
        print("file info...")
        print("  keys: ",keys)
        print("")

        if "AuxiliaryData" in keys:
            aux = f["AuxiliaryData"]
            for item in aux: print("  AuxiliaryData: ",item)

        if "Provenance" in keys:
            prov = f["Provenance"]
            for item in prov: print("  Provenance: ",item)

        # load data
        print("basic data:")
        if "network" in keys:
            network = f["network"][:]
            print("  network: shape = ",network.shape)
        if "station" in keys:
            station = f["station"][:]
            print("  station: shape = ",station.shape)
        if "coords" in keys:
            coords = f["coords"][:]
            print("  coords : shape = ",coords.shape)
        if "time" in keys:
            time = f["time"][:]
            print("  time   : shape = ",time.shape)
        if "channel_press" in keys:
            channel_press = f["channel_press"][:]
            print("  channel_press: shape = ",channel_press.shape)
        if "channel" in keys:
            channel = f["channel"][:]
            print("  channel: shape = ",channel.shape)
        print("")

        # info
        nrec = len(station)
        print("number of stations  : ",nrec)
        ntime = len(time)
        print("number of time steps: ",ntime)
        print("")

        # checks lengths
        if len(network) != nrec:
            print("Error: invalid data: station and network arrays have different different lengths",nrec,len(network))
            sys.exit(1)
        if len(coords[0,:]) != nrec:
            print("Error: invalid data: station and coords arrays have different different lengths",nrec,len(coords))
            sys.exit(1)

        print("seismogram data:")
        data = []
        data_label = []
        if "displ" in keys:
            displ = f["displ"][:]
            print("  displacement: shape = ",displ.shape)
            data.append(displ)
            data_label.append("displacement (m)")
        if "veloc" in keys:
            veloc = f["veloc"][:]
            print("  velocity    : shape = ",veloc.shape)
            data.append(veloc)
            data_label.append("velocity (m/s)")
        if "accel" in keys:
            accel = f["accel"][:]
            print("  acceleration: shape = ",accel.shape)
            data.append(accel)
            data_label.append("acceleration (m/s^2)")
        if "press" in keys:
            press = f["press"][:]
            print("  pressure    : shape = ",press.shape)
            data.append(press)
            data_label.append("pressure (Pa)")
        print("")

    # read reference seismograms
    if do_compare_ref:
        fx20 = np.loadtxt("./REF_SEIS/DB.X20.BXZ.semd")
        fx30 = np.loadtxt("./REF_SEIS/DB.X30.BXZ.semd")
        fx40 = np.loadtxt("./REF_SEIS/DB.X40.BXZ.semd")
        fx50 = np.loadtxt("./REF_SEIS/DB.X50.BXZ.semd")
        fx_ref = np.array([fx20, fx30, fx40, fx50])

    # plot
    if show_plot:
        # figure
        for id,data in enumerate(data):
            # label
            label = data_label[id]
            print("showing: ",label)
            # plot 4 stations
            if nrec % 4 == 0:
                npages = int(nrec/4)
            else:
                npages = int(nrec/4) + 1
            for ipage in range(npages):
                # single page - plot with 4 station maximum
                if nrec % 4 == 0:
                    nstas = 4
                else:
                    if nrec - (ipage+1)*4 > 0:
                        nstas = 4
                    else:
                        nstas = (nrec - ipage*4) % 4
                fig, ax = plt.subplots(2, 2, figsize=(10, 10))
                for i in range(nstas):
                    # station index
                    ista = i + ipage*4

                    # station name
                    sta = station[ista]
                    # avoid bytes string issues with strings like b'Hello', converts to text string
                    if isinstance(sta, (bytes, bytearray)): sta = sta.decode("utf-8")
                    # network name
                    net = network[ista]
                    # avoid bytes string issues with strings like b'Hello', converts to text string
                    if isinstance(net, (bytes, bytearray)): net = net.decode("utf-8")
                    # channel
                    if "press" in label:
                        cha = channel_press[0]
                    else:
                        cha = channel[2]
                    # avoid bytes string issues with strings like b'Hello', converts to text string
                    if isinstance(cha, (bytes, bytearray)): cha = cha.decode("utf-8")

                    name = "{}.{}.{}".format(net,sta,cha)
                    print("  station: ",name)

                    if "press" in label:
                        # pressure single component
                        ax[i // 2, i % 2].plot(time, data[ista, :, 2], "b")
                    else:
                        # plot z component
                        ax[i // 2, i % 2].plot(time, data[ista, :, 2], "b")

                    # trace comparisons
                    if do_compare_ref:
                        ax[i // 2, i % 2].plot(fx_ref[i][:, 0], fx_ref[i][:, 1], "k--")

                    # labels
                    ax[i // 2, i % 2].set_title(name)
                    ax[i // 2, i % 2].set_xlabel("time (s)")
                    ax[i // 2, i % 2].set_ylabel(label)

                plt.tight_layout()
                plt.show()

        # Wait for the user input to terminate the program
        input("Press any key to terminate the program")

    if convert_waveforms:
        print("converting to ASCII:")

        # takes same output directory as input file
        # for example: ./plot_HDF5_data_info.py OUTPUT_FILES/seismograms.h5
        #              -> output in OUTPUT_FILES/ directory
        out_dir = os.path.dirname(file_in)
        if out_dir == '':
            out_dir = "./"  # -> ./ + DB.STA.XXX.semd
        else:
            out_dir = out_dir + "/"  # OUTPUT_FILES -> OUTPUT_FILES/ + DB.STA.XXX.semd

        for id,data in enumerate(data):
            # label
            label = data_label[id]
            for irec in range(nrec):
                # station name
                sta = station[irec]
                # avoid bytes string issues with strings like b'Hello', converts to text string
                if isinstance(sta, (bytes, bytearray)): sta = sta.decode("utf-8")
                # network name
                net = network[irec]
                # avoid bytes string issues with strings like b'Hello', converts to text string
                if isinstance(net, (bytes, bytearray)): net = net.decode("utf-8")

                # components
                if "press" in label:
                    ncomp = 1
                else:
                    ncomp = 3
                for icomp in range(ncomp):
                    # trace data
                    if "press" in label:
                        # pressure single-component w/ shape (nrec,nstep)
                        trace = data[irec, :]
                    else:
                        # displ/veloc/accel w/ shape (nrec,nstep,ndim)
                        trace = data[irec, :, icomp]

                    length = len(trace)

                    # check if time and trace have same length
                    if length != len(time):
                        print("Error: different trace / time lengths: {} / {}".format(length,len(time)))
                        sys.exit(1)

                    # channel
                    if "press" in label:
                        cha = channel_press[icomp]
                    else:
                        cha = channel[icomp]
                    # avoid bytes string issues with strings like b'Hello', converts to text string
                    if isinstance(cha, (bytes, bytearray)): cha = cha.decode("utf-8")

                    # file ending
                    if "displ" in label:
                        ending = ".semd"
                    elif "veloc" in label:
                        ending = ".semv"
                    elif "accel" in label:
                        ending = ".sema"
                    elif "press" in label:
                        ending = ".semp"
                    else:
                        print("unknown data type")
                        sys.exit(1)

                    name = "{}.{}.{}".format(net,sta,cha)
                    print("  station: ",name)

                    filename = out_dir + name + ending

                    # single file per station and component
                    # file output
                    f = open(filename, "w")

                    t0 = time[0]
                    # file header
                    f.write("# STATION %s\n" % (sta))
                    f.write("# CHANNEL %s\n" % (cha))
                    f.write("# START_TIME %s\n" % (str(t0)))
                    #f.write("# SAMP_FREQ %f\n" % (tr.stats.sampling_rate))
                    f.write("# NDAT %d\n" % (length))

                    # data section
                    # fills numpy array
                    xy = np.empty(shape=(0,2))
                    for ii in range(length):
                        t = time[ii]
                        val = trace[ii]
                        xy = np.append(xy,[[t,val]],axis=0)

                    # data column
                    #print(xy[:,1])
                    #print(xy[0,1],xy[1,1],xy[2,1],xy[3,1])

                    # saves as ascii
                    np.savetxt(f, xy, fmt="%f\t%e")
                    f.close()

                    # user output
                    print("  written file: ",filename)

    print("")
    print("all done")
    print("")


def usage():
    print("usage: ./plot_HDF5_data_info.py filename[e.g. seismograms.h5] [show] [convert]")
    print(" with")
    print("   filename  - e.g. OUTPUT_FILES/seismograms.h5")
    print("   show      - (optional) plot waveforms")
    print("   convert   - (optional) convert seismograms to ASCII traces")
    sys.exit(1)


if __name__ == "__main__":
    # gets arguments
    do_plot_waveforms = False
    do_convert_waveforms = False

    if len(sys.argv) == 2:
        filename = sys.argv[1]
    elif len(sys.argv) == 3:
        filename = sys.argv[1]
        if sys.argv[2] == "show": do_plot_waveforms = True
    elif len(sys.argv) == 4:
        filename = sys.argv[1]
        if sys.argv[2] == "show": do_plot_waveforms = True
        if sys.argv[3] == "convert": do_convert_waveforms = True
    else:
        usage()
        sys.exit(1)

    plot_HDF5_data_info(filename,do_plot_waveforms,do_convert_waveforms)


