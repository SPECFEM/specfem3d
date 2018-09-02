#! /usr/bin/env pvpython
#
#
# usage: ./paraviewpython-example.py alpha_kernel.pvsm [2]
#
# creates jpg: image*.jpg

import os
import sys
import fileinput
import string

from paraview import servermanager

## input:
if len(sys.argv) == 2:
    filename = str(sys.argv[1])
    number = ""
elif len(sys.argv) == 3:
    filename = str(sys.argv[1])
    number = str(sys.argv[2])
else:
    print "usage: ./paraviewpython-example.py state-file [counter]"
    print "  with"
    print "    state-file - paraview state file, e.g. alpha_kernel.pvsm"
    print "    counter    - (optional) integer counter appended to filename, e.g. 0"
    sys.exit(1)

#outfile = "paraview_movie." + number
outfile = "image" + number
print "file root: ",outfile

## paraview
servermanager.Connect()
view = servermanager.CreateRenderView()
servermanager.LoadState(filename)
view = servermanager.GetRenderView()

# turn off axis visibility
view.CenterAxesVisibility = 0
view.OrientationAxesVisibility = 0

# to avoid segmentation fault
if servermanager.vtkSMProxyManager.GetVersionMajor() <= 5 and servermanager.vtkSMProxyManager.GetVersionMinor() < 5:
    view.UseOffscreenRenderingForScreenshots = 0

# sets view size for display
#view.ViewSize = [1920,1080] # width x height
print "view size: " + str(view.ViewSize)


## save as jpeg
jpegfilename = outfile + ".jpg"
print "plot to: " + jpegfilename
view.WriteImage(jpegfilename, "vtkJPEGWriter", 1)
print

## save as png
#print "plot to: " + "image.png"
#view.WriteImage("image.png", "vtkPNGWriter",1)

## save as bmp
#print "plot to: " + "image.bmp"
#view.WriteImage("image.bmp", "vtkBMPWriter",1)

