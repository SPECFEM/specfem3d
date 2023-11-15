from __future__ import print_function

import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.pyplot import figure, show, axes, sci
from matplotlib import cm, colors
from matplotlib.font_manager import FontProperties
from numpy import amin, amax, ravel
from numpy.random import rand

#*********************************************************************
nrcvr_x = 91
nrcvr_y = 4

dl = 1

x_init = 2.0
y_init = 8.0
dx = 1
dy = 8

rcvr_elev = 0
rcvr_bur = 0.2
rcvr_file = open('STATIONS', "w")

rcvr_ind = 0
for rcvr_y in range(0, nrcvr_y):
    for rcvr_x in range(0, nrcvr_x):
	rcvr_ind += 1
	rcvr_x_loc = (x_init) + rcvr_x*dx*dl
	rcvr_y_loc = (y_init) + rcvr_y*dy*dl

        rcvr_file.write('X' + str(rcvr_ind) + '\t DB \t' + str(rcvr_y_loc) + '\t \t' + str(rcvr_x_loc) + '\t \t' + str(rcvr_elev) + '\t' + str(rcvr_bur) + '\n')


rcvr_file.close()
