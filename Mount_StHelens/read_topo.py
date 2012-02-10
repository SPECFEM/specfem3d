#!/usr/bin/env python

import cubit

import os
import sys
import fileinput
import string
import math

print sys.path

#############################################################
# USER PARAMETERS

# topography file, data points per lon-increment
inputFile = 'ptopo.mean.utm'

# X coordinate in topography file repeats after line
nstep = 34

#############################################################

# converts xyz to utm
#os.system('./convert_lonlat2utm.pl ptopo.mean.xyz 10 > ptopo.mean.utm ')

# reads in points
print '#reading from file: ',inputFile

cubit.cmd('reset')
cubit.cmd('echo off')
cubit.cmd('Graphics Pause')


# creates point vertices
print '#creating points...'
count = 0
for line in fileinput.input( inputFile ):
  count = count + 1
  #print '# '+str(count)+' point: '+line
  lineitems = line.split()
  x = lineitems[0]
  y = lineitems[1]
  z = lineitems[2]
  xyz = str(x) + ' ' + str(y) + ' ' + str(z)
  #print '#point: ',xyz
  cubit.cmd('create vertex '+ xyz )
fileinput.close()
print '#done points: '+str(count)
print ''

cubit.cmd('Display')    

# creates smooth spline curves for surface
print '#creating curves...'
countcurves = 0
for i in range(1,count+1):
  if i > 1 :
    if i % nstep == 0 :
      countcurves = countcurves + 1
      cubit.cmd('create curve spline vertex '+str(i-nstep+1) + ' to ' + str(i) + ' delete' )
print '#done curves: '+str(countcurves)
print ''

cubit.cmd('Display')    
cubit.cmd('pause')

# creates surface
print '#creating skin surface...'
cubit.cmd('create surface skin curve all')


cubit.cmd('Display')    
cubit.cmd('pause')

# cleans up
cubit.cmd('merge all ')
cubit.cmd('delete vertex all')
cubit.cmd('delete curve all')

print '#done cleaning up'
cubit.cmd('Display')
cubit.cmd('echo on')

# saves and exports surface
# cubit file (uses ACIS file format)
cubit.cmd('save as "topo.cub" overwrite')
# export surface as STL file
cubit.cmd('export stl ascii "topo.stl" overwrite')

print '#exporting done'
print '#finished'
