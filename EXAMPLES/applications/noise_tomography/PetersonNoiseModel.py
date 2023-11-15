# based on
# PetersonNoiseModel.m
from __future__ import print_function

#from math import pi,log10

#import numpy as np
from numpy import pi,log10,array

def PetersonNoiseModel(Period,s=''):

    #debug
    #print "PetersonNoiseModel: ",Period,s

    # initialize
    accel = 0.0
    veloc = 0.0
    displ = 0.0

    if s == 'NLNM':
        # Peterson's New Low Noise Model
        P = array([0.1,0.17,0.4,0.8,1.24,2.4,4.3,5.0,6.0,10.0,12.0,15.6,21.9,31.6,45.0,70.0,101.0,154.0,328.0,600.0,10000.0])
        A = array([-162.36,-166.7,-170.0,-166.4,-168.6,-159.98,-141.1,-71.36,-97.26,-132.18,-205.27,-37.65,-114.37,-160.58,-187.5,-216.47,-185.0,-168.34,-217.43,-258.28,-346.88])
        B = array([5.64,0.0,-8.3,28.9,52.48,29.81,0.0,-99.77,-66.49,-31.57,36.16,-104.33,-47.1,-16.28,0.0,15.7,0.0,-7.61,11.9,26.6,48.75])
    elif s == 'NHNM':
        # Peterson's New High Noise Model
        P = array([0.1,0.22,0.32,0.8,3.8,4.6,6.3,7.9,15.4,20.0,354.8])
        A = array([-108.73,-150.34,-122.31,-116.85,-108.48,-74.66,0.66,-93.37,73.54,-151.52,-206.66])
        B = array([-17.23,-80.5,-23.87,32.51,18.08,-32.95,-127.18,-22.42,-162.98,10.01,31.63])

    N = len(P)
    l = 1

    # Period is less than P(1) (minimum period defined in the model)
    if Period < P[0]:
        accel = A[0] + B[0] * log10(P[0])
        veloc = accel + 20.0 * log10(P[0] / 2.0 / pi)
        displ = accel + 20.0 * log10(P[0]**2 / 4.0 / pi**2)

    # Period is between P(1) and P(N)
    while l <= N - 2:
        if Period >= P[l] and Period < P[l + 1]:
            accel = A[l] + B[l] * log10(Period)
            veloc = accel + 20.0 * log10(Period / 2.0 / pi)
            displ = accel + 20.0 * log10(Period**2 / 4.0 / pi**2)
            break
        else:
            l = l + 1

    # Period is larger than P(N) and less than 1e5 (maximum period defined in the model)
    if Period >= P[N-1] and Period < 100000.0:
        accel = A[N-1] + B[N-1] * log10(Period)
        veloc = accel + 20.0 * log10(Period / 2.0 / pi)
        displ = accel + 20.0 * log10(Period**2 / 4.0 / pi**2)
    elif Period > 100000.0:
        accel = A[N-1] + B[N-1] * log10(100000.0)
        veloc = accel + 20.0 * log10(100000.0 / 2.0 / pi)
        displ = accel + 20.0 * log10(100000.0**2 / 4.0 / pi ** 2)

    #debug
    #print "PetersonNoiseModel: ",accel,veloc,displ

    return accel,veloc,displ

