#!/usr/bin/env python
#
# creates an external source time function file
#
# example:
# ./create_external_source_time_function.py OUTPUT_FILES/plot_source_time_function.txt
#
# requires access to: DATA/Par_file
#
import sys
import os
import array

import numpy as np


def write_ascii_file(filename,data):
    """
    writes data array to ascii file
    """
    with open(filename, mode='w') as f:
        # user output
        print("  filename: ",filename)
        print("  number of points in array = ",len(data))

        # comment lines
        f.write("# Source time function\n")
        f.write("#   number of time steps: {}\n".format(len(data)))
        # data values
        f.write("# STF values\n")
        for i in range(len(data)):
            val = data[i]
            f.write("{}\n".format(val))

        print("  file written")
        print("")


def write_binary_file(filename,data):
    """
    writes data array to binary file
    """
    # float 'f' by default (otherwise 'd' for double precision)
    custom_type = 'f'

    with open(filename, mode='wb') as f:
        # gets array length in bytes
        # marker
        binlength = array.array('i')
        num_points = data.size
        if custom_type == 'f':
            # float (4 bytes) for single precision
            binlength.fromlist([num_points * 4])
        else:
            # double precision
            binlength.fromlist([num_points * 8])

        # user output
        print("  filename: ",filename)
        print("  array length = ",binlength," Bytes")
        print("  number of points in array = ",num_points)

        # writes array data
        binvalues = array.array(custom_type)

        data = np.reshape(data, (num_points), order='F') # fortran-like index ordering
        #print("debug: data ",data.tolist())

        binvalues.fromlist(data.tolist()) #.tolist())

        # fortran binary file: file starts with array-size again
        binlength.tofile(f)
        # data
        binvalues.tofile(f)
        # fortran binary file: file ends with array-size again
        binlength.tofile(f)

        print("  file written")
        print("")


def create_external_source_time_function(input_file,use_ascii=True):
    """
    creates an external STF file in DATA/
    """

    print("")
    print("creating external source time function file: DATA/stf.dat")
    print("  input file  : ",input_file)
    print("")

    if len(input_file) == 0: usage()

    # checks if file exists
    if not os.path.isfile(input_file):
        print("Please check if file exists: ",input_file)
        sys.exit(1)

    # reads in input file
    with open(input_file, 'r') as f:
        content = f.readlines()

    # checks number of lines
    if len(content) == 0:
        print("File is empty, please check input file")
        sys.exit(1)

    # get data entries
    stf_data = np.array([])
    for line in content:
        # count valid data lines
        if "#" in line or "!" in line:
            # comment line
            continue
        else:
            # assumes this is a data line
            # cleans line from whitespace and extra characters
            line = line.strip()
            # check for empty line
            if len(line) == 0: continue

            # reads in value(s)
            numbers = line.split()

            if len(numbers) == 0:
                # nothing to read
                continue
            if len(numbers) == 1:
                # single data entry per line
                val = float(numbers[0])
            elif len(numbers) == 2:
                # line contains two numbers (like plot_source_time_function.txt)
                # we'll take the second
                val = float(numbers[1])
            else:
                # line contains multiple numbers
                # we'll take the last
                val = float(numbers[-1])

            # store data value
            stf_data = np.append(stf_data,val)

    print("number of STF time steps found: ",len(stf_data))
    print("")

    # checks with Par_file setting
    filename = "DATA/Par_file"
    with open(filename,'r') as f:
        content_par = f.readlines()

    NSTEP = None
    for line in content_par:
        # skip comments
        if "#" in line or "!" in line: continue

        # search for NSTEP
        if "NSTEP" in line and "=" in line:
            NSTEP = int(line.split('=')[1])

    print("Par_file:")
    print("  simulation NSTEP = ",NSTEP)
    print("")

    # checks
    if NSTEP is None:
        print("Couldn't find NSTEP in DATA/Par_file")
        sys.exit(1)

    # checks lengths
    # must have at least NSTEP data entries in input file
    if NSTEP > len(stf_data):
        print("Error: simulation length NSTEP = {} is longer than the available STF data {}".format(NSTEP,len(stf_data)))
        print("       Please use another input file...")
        sys.exit(1)

    # warning if data length is too long
    if NSTEP < len(stf_data):
        print("WARNING: simulation length NSTEP == {} is shorter than the available STF data {}".format(NSTEP,len(stf_data)))
        print("         Only the first {} data entries will be written to the output STF file, all the others will be ignored".format(NSTEP))
        print("")

    # cuts STF data to simulation length
    stf_data = stf_data[0:NSTEP]

    # stats
    stf_min = stf_data.min()
    stf_max = stf_data.max()
    print("stats:")
    print("  source time function values: min/max = {} / {}".format(stf_min,stf_max))
    print("")

    # output STF file
    if use_ascii:
        filename = "DATA/stf.dat"
        print("writing ascii output:")
        write_ascii_file(filename,stf_data)
    else:
        filename ="DATA/stf.dat.bin"
        print("writing binary output:")
        write_binary_file(filename,stf_data)

    print("written external source time function to: ",filename)
    print("")
    print("all done")
    print("")


def usage():
    print("usage: ./create_external_source_time_function.py input_file [--ascii or --binary]")
    print("  with")
    print("     input_file              - ascii file with trace, e.g., OUTPUT_FILES/plot_source_time_function.txt")
    print("     --ascii                 - (optional) writes out STF as ascii file (default)")
    print("     --binary                - (optional) writes out STF as binary file")
    sys.exit(1)


if __name__ == '__main__':
    # init
    input_file = ""
    use_ascii = True
    # reads arguments
    #print("\nnumber of arguments: " + str(len(sys.argv)))
    if len(sys.argv) <= 1: usage()

    input_file = sys.argv[1]
    if len(sys.argv) == 3:
        if sys.argv[2] == "--ascii":
            use_ascii = True
        elif sys.argv[2] == "--binary":
            use_ascii = False

    # main routine
    create_external_source_time_function(input_file,use_ascii)
