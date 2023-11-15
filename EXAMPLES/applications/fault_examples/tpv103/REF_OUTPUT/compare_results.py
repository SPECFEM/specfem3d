import numpy as np
import sys

def usage():
    print "compare_results filename1 filename2"

#print len(sys.argv)
if len(sys.argv) != 3:
    usage()
dat = np.loadtxt(sys.argv[1], skiprows = 21)
dat2= np.loadtxt(sys.argv[2], skiprows = 21)


