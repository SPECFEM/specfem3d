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
nsrc_x = 6
nsrc_y = 3

dl = 18.75

x_init = 9
y_init = 21
#ix_fin = 119
#y_fin = 78
dx = 22
dy = 28


src_ind=0
for src_y in range(0, nsrc_y):
    for src_x in range(0, nsrc_x):

	src_x_loc = (x_init-1)*dl + src_x*dx*dl
	src_y_loc = (y_init-1)*dl + src_y*dy*dl

        src_file = open("FORCESOLUTION_" + str(src_ind).zfill(6), "w")
        src_file.write('FORCE  001 \n')
        src_file.write('time shift: \t \t 0.000 \n')
        src_file.write('hdurorf0: \t \t 5.0 \n')
        src_file.write('latorUTM: \t \t ' + str(src_y_loc) + '\n')
        src_file.write('longorUTM: \t \t ' + str(src_x_loc) + '\n')
        src_file.write('depth: \t \t \t -0.01 \n')
        src_file.write('source time function: \t \t 0 \n')
        src_file.write('factor force source: \t \t 1.d15 \n')
        src_file.write('component dir vect source E: \t 0.d0 \n')
        src_file.write('component dir vect source N: \t 0.d0 \n')
        src_file.write('component dir vect source Z_UP: -1.d0 \n')
        src_file.close()
	src_ind += 1
