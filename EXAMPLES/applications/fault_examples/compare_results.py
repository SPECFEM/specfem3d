#!/usr/bin/env python
#
# compares results in a multi-column file, e.g. faultst-120dp000.dat
#
#########################################################################
from __future__ import print_function

import numpy as np
import numpy.linalg as alg
import sys

def usage():
    print("compare_results filename1 filename2")

def diff_measure(x1,x2):
    p = alg.norm(x1 + x2)
    if p != 0.0:
        measure = alg.norm(x1 - x2)/p
    else:
        measure = alg.norm(x1 - x2)
    return measure


def main(argv):
    #print(len(sys.argv))
    if len(argv) != 2:
        usage()
        exit()

    # reads in data, skipping first rows (comment lines)
    dat1 = np.loadtxt(argv[0], skiprows = 21)
    dat2 = np.loadtxt(argv[1], skiprows = 21)

    _,n = dat1.shape
    for index in range(n):
        col1 = dat1[:,index]
        col2 = dat2[:,index]
        diff = diff_measure(col1, col2);
        print("the difference between the 2 files for column ", index, " is : ",diff)

if __name__=="__main__":
    main(sys.argv[1:])
