#!/usr/bin/env python
#
# plot the cross-correlation and L2-norm between reference and output seismograms
#
import sys
import glob
import os
import numpy as np

# tolerance values
TOL_CORR = 0.8
TOL_ERR = 0.01
TOL_SHIFT = 0.01

###############################################################
# USER PARAMETERS

# computes correlations within a moving window
# (set either to False or True)
USE_SUB_WINDOW_CORR = False

# for moving window correlations:
# apprixomate minimum period of simulation
# (default: NEX = 80 -> T_min = 256/80 * 17 s = 54.4 s)
TMIN = 54.4
###############################################################

def get_cross_correlation_timeshift(x,y,dt):
    """
    computes the time shift of the maximum cross-correlation of signal x with respect to y
    """
    # checks signals
    if len(x) != len(y):
        print "Error: lengths in cross-correlation don't match"
        return 1.e30

    # cross-correlation length
    signal_length = len(x)
    length = 2 * signal_length - 1

    # cross-correlation array
    crosscorrelation = np.correlate(x, y, mode='full')

    # index of maximum (between [0,2 * signal_length - 1]
    indexmax = np.argmax(crosscorrelation)

    # position (negative -> signal shifted to right, positive -> signal shifted to left)
    # time lag (will have steps of dt)
    lag = (indexmax + 1) - signal_length

    # subsample precision
    maxval = crosscorrelation[indexmax]

    #debug
    #print "xcorr: ",indexmax,maxval,len(crosscorrelation),length

    # gets values left/right from maximum value
    if indexmax >= 1 and indexmax < length-1:
        val_left = crosscorrelation[indexmax-1]
        val_right = crosscorrelation[indexmax+1]
    elif indexmax == 0:
        # wrap-around values
        val_left = crosscorrelation[length-1]
        val_right = crosscorrelation[1]
    elif indexmax == length-1:
        # wrap-around values
        val_left = crosscorrelation[indexmax-1]
        val_right = crosscorrelation[0]

    # quadratic interpolation will give this maximum
    # see: http://www.dsprelated.com/freebooks/sasp/Peak_Detection_Steps_3.html
    if (val_left - 2.0*maxval + val_right) != 0.0:
        peak_shift = 0.5 * (val_left - val_right) / (val_left - 2.0*maxval + val_right)
    else:
        peak_shift = 0.0

    # adds subsample shift
    lag += peak_shift

    # cross-correlation time lag
    time_shift = lag * dt

    # debug
    #print "cross-correlation:",length,signal_length,"shift = ",indexmax,lag,time_shift

    return time_shift


def plot_correlations(out_dir,ref_dir):
    """
    plots correlation and L2-norm values between reference and output seismograms
    """
    print('comparing seismograms')
    print('  reference directory: %s' % ref_dir)
    print('  output directory   : %s\n' % out_dir)

    # seismogram file ending
    ending = '.sem*' # .semd, .semv, .sema, .semp, ..

    # gets seismograms
    files = glob.glob(out_dir + '/*' + ending)
    if len(files) == 0:
        print "no seismogram files with ending ",ending," found"
        print "Please check directory: ",out_dir
        sys.exit(1)

    files.sort()

    corr_min = 1.0
    err_max = 0.0
    shift_max = 0.0

    # gets time step size from first file
    syn_file = files[0]
    print "  time step: reading from first file ",syn_file
    syn_time = np.loadtxt(syn_file)[:, 0]
    dt = syn_time[1] - syn_time[0]
    print "  time step: size = ",dt
    # start time
    dt_start = syn_time[0]

    # warning
    if dt <= 0.0:
        print "warning: invalid time step size for file ",files[0]

    # determines window length
    if USE_SUB_WINDOW_CORR:
        # moving window
        print "  using correlations in moving sub-windows"
        print "  minimum period: ",TMIN
        # checks
        if dt <= 0.0:
            # use no moving window
            window_length = len(syn_time) - 1
        else:
            # window length for minimum period
            window_length = int(TMIN/dt)

        print "  moving window length: ",window_length


    print ""
    print "comparing ",len(files),"seismograms"
    print ""

    # outputs table header
    print("|%-30s| %13s| %13s| %13s|" % ('file name', 'corr', 'err', 'time shift'))

    # counter
    n = 0

    for f in files:
        # build reference and synthetics file names
        # specfem file: **network**.**station**.**comp**.sem.ascii
        fname = os.path.basename(f)
        names = str.split(fname,".")

        # trace
        net = names[0]
        sta = names[1]
        cha = names[2]

        # filenames
        # old format
        #fname_old = sta + '.' + net + '.' + cha + '.sem.ascii'
        #ref_file = ref_dir + '/' + fname_old
        #syn_file = out_dir + '/' + fname_old
        # new format
        ref_file = ref_dir + '/' + fname
        syn_file = out_dir + '/' + fname

        # makes sure files are both available
        if not os.path.isfile(ref_file):
            print "  file " + ref_file + " not found"
            continue
        if not os.path.isfile(syn_file):
            print "  file " + syn_file + " not found"
            continue

        # numpy: reads in file data
        ref0 = np.loadtxt(ref_file)[:, 1]
        syn0 = np.loadtxt(syn_file)[:, 1]

        #debug
        #print "  seismogram: ", fname, "vs", fname_old,"  lengths: ",len(ref0),len(syn0)

        # cuts common length
        length = min(len(ref0),len(syn0))
        if length <= 1: continue

        # length warning
        if len(ref0) != len(syn0):
          print("** warning: mismatch of file length in both files syn/ref = %d / %d" %(len(syn0),len(ref0)))
          #print("** warning: using smaller length %d" % length)

        # time step size in reference file
        ref_time = np.loadtxt(ref_file)[:, 0]
        dt_ref = ref_time[1] - ref_time[0]
        # start time
        dt_ref_start = ref_time[0]

        # mismatch warning
        if abs(dt - dt_ref)/dt > 1.e-5:
          print("** warning: mismatch of time step size in both files syn/ref = %e / %e" %(dt,dt_ref))
          #print("** warning: using time step size %e from file %s" %(dt,syn_file))

        #debug
        #print "common length: ",length

        ref = ref0[0:length]
        syn = syn0[0:length]

        # least square test
        norm = np.linalg.norm
        sqrt = np.sqrt

        # normalized by power in reference solution
        fac_norm = norm(ref)
        # or normalized by power in (ref*syn)
        #fac_norm = sqrt(norm(ref)*norm(syn))

        if fac_norm > 0.0:
            err = norm(ref-syn)/fac_norm
        else:
            err = norm(ref-syn)

        #debug
        #print('norm syn = %e norm ref = %e' % (norm(syn),fac_norm))

        # correlation test
        # total length
        if fac_norm > 0.0:
            corr_mat = np.corrcoef(ref, syn)
        else:
            if norm(ref-syn) > 0.0:
                corr_mat = np.cov(ref-syn)
            else:
                # both zero traces
                print("** warning: comparing zero traces")
                corr_mat = 1.0
        corr = np.min(corr_mat)

        # time shift
        if fac_norm > 0.0:
          # shift (in s) by cross correlation
          shift = get_cross_correlation_timeshift(ref,syn,dt)
        else:
          # no correlation with zero trace
          shift = 0.0

        # correlation in moving window
        if USE_SUB_WINDOW_CORR:
            # moves window through seismogram
            for i in range(0,length-window_length):
                # windowed signals
                x = ref[i:i+window_length]
                y = syn[i:i+window_length]

                # correlations
                corr_win = np.corrcoef(x, y)
                corr_w = np.min(corr_win)
                corr = min(corr, corr_w)

                # cross-correlation array
                shift_w = get_cross_correlation_timeshift(x,y,dt)
                if abs(shift) < abs(shift_w): shift = shift_w

        # adding shift in start times
        shift += (dt_ref_start - dt_start)

        # statistics
        corr_min = min(corr, corr_min)
        err_max = max(err, err_max)
        if abs(shift_max) < abs(shift): shift_max = shift

        # info string
        info = ""
        if corr < TOL_CORR:   info += "  poor correlation"
        if err > TOL_ERR:     info += "      poor match"
        if abs(shift) > TOL_SHIFT: info += "      significant shift"

        # print results to screen
        print("|%-30s| %13.5f| %13.5le| %13.5le| %s" % (fname, corr, err, shift, info))

        # counter
        n += 1


    # check if any comparison done
    if n == 0:
        # values indicating failure
        corr_min = 0.0
        err_max = 1.e9
        shift_max = 1.e9

    # print min(coor) max(err)
    print("|---------------------------------------------------------------------------|")
    print("|%30s| %13.5f| %13.5le| %13.5le|" % ('min/max', corr_min, err_max, shift_max))

    # output summary
    print("\nsummary:")
    print("%d seismograms compared\n" % n)
    if n == 0:
        print("\nno seismograms found for comparison!\n\n")

    print("correlations: values 1.0 perfect, < %.1f poor correlation" % TOL_CORR)
    if corr_min < TOL_CORR:
        print("              poor correlation seismograms found")
    else:
        print("              no poor correlations found")
    print ""

    print("L2-error    : values 0.0 perfect, > %.2f poor match" % TOL_ERR)
    if err_max > TOL_ERR:
        print("              poor matching seismograms found")
    else:
        print("              no poor matches found")
    print ""

    print("Time shift  : values 0.0 perfect, > %.2f significant shift" % TOL_SHIFT)
    if abs(shift_max) > TOL_SHIFT:
        print("              significant time shift in seismograms found")
    else:
        print("              no significant time shifts found")
    print ""


def usage():
    print "usage: ./compare_seismogram_correlations.py directory1/ directory2/"
    print "  with"
    print "     directory1 - directory holding seismogram files (***.sem.ascii),"
    print "                    e.g. OUTPUT_FILES/"
    print "     directory2 - directory holding corresponding reference seismogram files,"
    print "                    e.g. OUTPUT_FILES_reference_OK/"

if __name__ == '__main__':
    # gets arguments
    if len(sys.argv) != 3:
        usage()
        sys.exit(1)
    else:
        out_dir = sys.argv[1]
        ref_dir = sys.argv[2]

    plot_correlations(out_dir,ref_dir)

