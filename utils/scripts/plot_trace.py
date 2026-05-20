#!/usr/bin/env python3
"""
Plot seismograms from SPECFEM3D output
"""

import sys,os
import matplotlib.pyplot as plt
from mpl_toolkits.axisartist.axislines import Axes
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes 
from mpl_toolkits.axes_grid1.inset_locator import mark_inset

import numpy as np
from scipy import signal


def apply_filter(data, dt, filter_type='bandpass', lowcut=None, highcut=None, order=4):
    """
    Apply Butterworth filter to seismogram data
    
    Parameters:
        data: input signal
        dt: sampling interval
        filter_type: 'lowpass', 'highpass', or 'bandpass'
        lowcut: low cutoff frequency (Hz)
        highcut: high cutoff frequency (Hz)
        order: filter order
    
    Returns:
        filtered data
    """
    # Calculate sampling frequency
    fs = 1.0 / dt
    nyquist = 0.5 * fs
    print(f"filter: Nyquist = {nyquist:.4f} Hz")

    if filter_type == 'lowpass':
        if highcut is None:
            raise ValueError("highcut frequency must be specified for lowpass filter")
        # Normalize the frequencies by the Nyquist frequency
        normal_cutoff = highcut / nyquist
        #b, a = signal.butter(order, normal_cutoff, btype='low', analog=False)
        # Use second-order sections format (output='sos' instead of the default 'ba')
        sos = signal.butter(order, normal_cutoff, btype='low', analog=False, output='sos')

    elif filter_type == 'highpass':
        if lowcut is None:
            raise ValueError("lowcut frequency must be specified for highpass filter")
        # Normalize the frequencies by the Nyquist frequency
        normal_cutoff = lowcut / nyquist
        #b, a = signal.butter(order, normal_cutoff, btype='high', analog=False)
        # Use second-order sections format (output='sos' instead of the default 'ba')
        sos = signal.butter(order, normal_cutoff, btype='high', analog=False, output='sos')

    elif filter_type == 'bandpass':
        if lowcut is None or highcut is None:
            raise ValueError("Both lowcut and highcut frequencies must be specified for bandpass filter")
        # Normalize the frequencies by the Nyquist frequency
        low = lowcut / nyquist
        high = highcut / nyquist
        # check low < high
        if low > high:
            raise ValueError("Invalid lowcut > highcut frequency, make sure lowcut is smaller than highcut")
        #b, a = signal.butter(order, [low, high], btype='band', analog=False)
        # Use second-order sections format (output='sos' instead of the default 'ba')
        sos = signal.butter(order, [low, high], btype='band', analog=False, output='sos')
    else:
        raise ValueError(f"Unknown filter type: {filter_type}")
    
    # Apply the filter forward and backward (zero-phase shift)
    #filtered_data = signal.filtfilt(b, a, data)
    # Use second-order sections format
    filtered_data = signal.sosfiltfilt(sos, data)

    return filtered_data


def plot_trace(files=[]):
    print("Plot trace:")
    print(f"files: {' '.join(files)}")
    print("")

    # checks
    if len(files) == 0: usage()
    for file in files:
        if not os.path.exists(file):
            print(f"File {file} not found.")
            sys.exit(1)

    times = []
    values = []

    # Read seismogram data
    for file in files:
        print(f"file: {file}")

        # data
        data = np.loadtxt(file)
        print(f"data: length           = {len(data)}")

        time = data[:, 0]
        trace = data[:, 1]

        if len(trace) <= 1:
            print("not enough data.")
            sys.exit(1)

        # filter options
        use_filter = False
        filter_type = 'bandpass'
        lowcut = 0.001   # 1 mHz
        highcut = 1.0    # 1 Hz
        order = 4

        if "--filter" in sys.argv:
            use_filter = True
        if "--filter-type" in sys.argv:
            use_filter = True
            idx = sys.argv.index("--filter-type")
            filter_type = sys.argv[idx + 1]
        if "--lowcut" in sys.argv:
            use_filter = True
            idx = sys.argv.index("--lowcut")
            lowcut = float(sys.argv[idx + 1])
        if "--highcut" in sys.argv:
            use_filter = True
            idx = sys.argv.index("--highcut")
            highcut = float(sys.argv[idx + 1])
        if "--filter-order" in sys.argv:
            use_filter = True
            idx = sys.argv.index("--filter-order")
            order = int(sys.argv[idx + 1])

        # Apply filter if requested
        if use_filter:
            # sampling interval
            dt = time[1] - time[0]

            print(f"filter: {filter_type} - lowcut {lowcut} (Hz) / highcut {highcut} (Hz) / order {order}")
            print(f"filter: sampling interval = {dt}")

            trace = apply_filter(trace, dt, filter_type, lowcut, highcut, order)

        print(f"      time  start/end = {time[0]} / {time[-1]}")
        print(f"      trace min/max   = {np.min(trace)} / {np.max(trace)}")
        print("")

        times.append(time)
        values.append(trace)

    # Determine figure layout
    # color-blind friendly colors
    if "--color-blind" in sys.argv:
        # see: https://matplotlib.org/stable/gallery/style_sheets/style_sheets_reference.html
        plt.style.use('tableau-colorblind10')
    else:
        # default, classic, fast, ..
        plt.style.use('classic')

    #ax = fig.add_subplot(axes_class=Axes)
    #ax.axis['right'].set_visible(False)
    #ax.axis['top'].set_visible(False)

    fig,ax = plt.subplots(figsize=(20, 8))

    # turn axis off
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.tick_params(top=False, right=False)

    # Plot
    ylabel = 'amplitude'
    for i,file in enumerate(files):
        basename = os.path.basename(file)
    
        time = times[i]
        trace = values[i]

        # line plot
        if len(files) == 1:
            if "--label" in sys.argv:
                idx = sys.argv.index("--label")
                label = sys.argv[idx + 1] if idx + 1 < len(sys.argv) else basename
            else:
                label = basename
            ax.plot(time, trace, 'b-', linewidth=1.0, label=f"{label}")
        else:
            labelN = f"--label{i+1}"
            if labelN in sys.argv:
                idx = sys.argv.index(labelN)
                label = sys.argv[idx + 1] if idx + 1 < len(sys.argv) else basename
            else:
                label = basename
            ax.plot(time, trace, '-', linewidth=1.0, label=f"{label}")

        # axis label
        if "semd" in basename.split("."):
            ylabel = 'displacement (m)'
        elif "semv" in basename.split("."):
            ylabel = 'velocity (m/s)'
        elif "sema" in basename.split("."):
            ylabel = 'acceleration (m/s^2)'

    # Only show x-axis label on bottom row
    ax.set_xlabel('time (s)', fontsize=10)
    ax.set_ylabel(ylabel, fontsize=10)

    # zoom-in window
    if "--zoom" in sys.argv:
        # get input x-range
        idx = sys.argv.index("--zoom")
        x1 = float(sys.argv[idx + 1]) if idx + 1 < len(sys.argv) else time[0]
        x2 = float(sys.argv[idx + 2]) if idx + 2 < len(sys.argv) else time[-1]
        print(f"  zoom window: x-range {x1} to {x2}")
        # select y range for zoomin
        y1 = float('inf')
        y2 = -float('inf')
        for i,_ in enumerate(files):
            time = times[i]
            trace = values[i]
            mask = (time >= x1) & (time <= x2)
            y1 = min(y1,np.min(trace[mask]))
            y2 = max(y2,np.max(trace[mask]))
        y1 = y1 * 1.05 if y1 < 0 else y1 * 0.95 # add margin
        y2 = y2 * 1.05 if y2 > 0 else y2 * 0.95
        print(f"               y min/max = {y1} / {y2}\n")
        # Plot a zoom-in graph
        axins = ax.inset_axes([0.7, 0.8, 0.25, 0.25])   # x0,y0,width,height
        #axins = zoomed_inset_axes(ax, 2, loc=1) # zoom = 2
        for i,_ in enumerate(files):
            time = times[i]
            trace = values[i]
            axins.plot(time, trace, linestyle='-', linewidth=1.0)
        # setting of a zoomed graph
        axins.set_xlim(x1, x2)
        axins.set_ylim(y1, y2)
        axins.set_xticks([])
        axins.set_yticks([])
        #ax.indicate_inset_zoom(axins, edgecolor="black", alpha=0.2)
        mark_inset(ax, axins, loc1=2, loc2=4, fc="none", ec="0.5", alpha=0.4)

    # Set overall title
    if "--title" in sys.argv:
        idx = sys.argv.index("--title")
        title = sys.argv[idx + 1]
        print(f"  title: {title}")
        plt.suptitle(f'{title}',fontsize=14, fontweight='bold')
    else:
        title = f"{'\n'.join(files)}"
        if use_filter:
            if filter_type == 'bandpass':
                title += f' - Bandpass: {lowcut}-{highcut} Hz'
            elif filter_type == 'lowpass':
                title += f' - Lowpass: <{highcut} Hz'
            elif filter_type == 'highpass':
                title += f' - Highpass: >{lowcut} Hz'
        plt.suptitle(title, fontsize=14, fontweight='bold')
    #plt.title(f'{file}',fontsize=11, fontweight='bold')

    # axis limits
    if "--xlim" in sys.argv:
        # get input x-range
        idx = sys.argv.index("--xlim")
        x1 = float(sys.argv[idx + 1]) if idx + 1 < len(sys.argv) else time[0]
        x2 = float(sys.argv[idx + 2]) if idx + 2 < len(sys.argv) else time[-1]
        print(f"  xlim: x-range {x1} to {x2}")
        ax.set_xlim(x1, x2)
    if "--ylim" in sys.argv:
        # get input x-range
        idx = sys.argv.index("--ylim")
        y1 = float(sys.argv[idx + 1]) if idx + 1 < len(sys.argv) else trace.min()
        y2 = float(sys.argv[idx + 2]) if idx + 2 < len(sys.argv) else trace.max()
        print(f"  ylim: y-range {y1} to {y2}")
        ax.set_ylim(y1, y2)

    #ax.grid(True, linestyle=':', alpha=0.5)
    #ax.ticklabel_format(style='scientific', axis='both', scilimits=(-3,3))
    ax.ticklabel_format(style='scientific', axis='y', scilimits=(-1,1))
    ax.legend(loc='upper left',frameon=False)  # or 'lower right'
    #plt.tight_layout()

    # Save
    output_file = "out.jpg"
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    #plt.savefig(output_file, dpi=300)
    print(f"\nsaved to: {output_file}\n")

    # show
    if "--show" in sys.argv: plt.show()


def usage():
    print("usage:")
    print("    ./plot_trace.py <file,e.g., DB.A1.FXZ.semd> [--title 'text'] [--label 'my trace' or --label1 'trace1', --label2 'trace2', ..]")
    print("                                                [--filter] [--filter-type type] [--lowcut val] [--highcut val] [--filter-order val]")
    print("                                                [--show] [--color-blind]")
    print("                                                [--xlim x1 x2] [--ylim y1 y2]")
    print("                                                [--zoom x1 x2]")
    print("  with")
    print("    title              - title text")
    print("    label              - trace label text")
    print("    filter             - Apply frequency filter to seismograms")
    print("    filter-type type   - Filter type, e.g., lowpass, highpass, bandpass (default bandpass)")
    print("    lowcut val         - Low cutoff frequency, val in Hz (for highpass/bandpass)")
    print("    highcut val        - High cutoff frequency, val in Hz (for lowpass/bandpass)")
    print("    filter-order val   - Butterworth filter order (default: 4)")
    print("    show               - show figure (default just plots to file")
    print("    color-blind        - use color-blind friendly palette")
    print("    xlim x1 x2         - limit x-axis range")
    print("    ylim y1 y2         - limit y-axis range")
    print("    zoom x1 x2         - add zoom-in figure for x-range x1 to x2")
    sys.exit(1)

if __name__ == "__main__":
    files = []
    if len(sys.argv) < 2: usage()

    file = sys.argv[1]
    files.append(file)

    # check for more files
    for i in range(2,len(sys.argv)):
        if ".semd" in sys.argv[i] or ".semv" in sys.argv[i] or ".sema" in sys.argv[i]:
            file = sys.argv[i]
            files.append(file)

    # Plot trace
    plot_trace(files)

