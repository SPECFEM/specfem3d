#!/usr/bin/env python
#
# see: http://seismicdata.github.io/pyasdf/tutorial.html
#
from __future__ import print_function
import sys
import pyasdf

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


def usage():
    print("usage: ./read_asdf_data_info.py filename[e.g. synthetics.h5] [show]")
    print(" with")
    print("   filename - e.g. synthetics.h5")
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
