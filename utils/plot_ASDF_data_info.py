#!/usr/bin/env python
#
# reads in an ASDF-format file and outputs infos about content and checks validity
# see: http://seismicdata.github.io/pyasdf/tutorial.html
#
from __future__ import print_function

import sys

try:
    import pyasdf
except:
    print("Error importing pyasdf, check if pyasdf module is installed")
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


def plot_ASDF_data_info(file,show_plot):
    print("")
    print("reading file: ",filename)
    print("")

    ds = pyasdf.ASDFDataSet(filename)

    # file info
    print("file info...")
    print(ds)
    print("")

    # validate
    print("validate...")
    ds.validate()
    print("")

    # event info
    print("event info...")
    print("  number of events: ",len(ds.events))
    print("")
    type(ds.events)
    print(ds.events)
    for event in ds.events:
        print(event)
    print("")

    # station info
    print("")
    print("station list...")
    print("  number of stations: ",len(ds.waveforms))
    print("")
    print(ds.waveforms.list())
    print("")

    for station in ds.waveforms.list():
        print("")
        print("station:")
        print("")
        print(station)

        # how to get station name and associated waveform?
        # waveforms
        sta = ds.waveforms[station]
        type(sta)
        print(sta.synthetic)
        print("")
        print("waveform list:",sta.list())
        print("waveform tags:",sta.get_waveform_tags())
        print("")

        # plotting
        if show_plot:
            print("plotting stream...")
            st = sta.synthetic
            st.plot()

    if show_plot:
        # Wait for the user input to terminate the program
        input("Press any key to terminate the program")

    print("")
    print("all done")
    print("")

def usage():
    print("usage: ./plot_ASDF_data_info.py filename[e.g. synthetic.h5] [show]")
    print(" with")
    print("   filename - e.g. OUTPUT_FILES/synthetic.h5")
    print("   show   - (optional) plot waveforms")
    sys.exit(1)


if __name__ == "__main__":
    # gets arguments
    do_plot_waveforms = False

    if len(sys.argv) == 2:
        filename = sys.argv[1]
    elif len(sys.argv) == 3:
        filename = sys.argv[1]
        if sys.argv[2] == "show": do_plot_waveforms = True
    else:
        usage()
        sys.exit(1)

    plot_ASDF_data_info(filename,do_plot_waveforms)
