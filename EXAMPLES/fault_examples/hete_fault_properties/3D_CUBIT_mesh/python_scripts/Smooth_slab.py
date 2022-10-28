import sys
import numpy
from   scipy.interpolate     import griddata
from   scipy.ndimage.filters import gaussian_filter

for i in range(len(sys.argv)):                        # read parameter from shell(contain the prefix)
  if(sys.argv[i].find("-p")==0):
    P1=sys.argv[i+1]
    P2=sys.argv[i+2]
    P3=sys.argv[i+3]
sigma       =   float(P2)
inc         =   float(P3)
input_file  =   P1
output_file =   P1 + ".smooth_" + P2

read_data   =   numpy.loadtxt(input_file,skiprows=0)
DeltaX      =   read_data[1,0] - read_data[0,0]
Xnum        =   int(abs(read_data[-1,0]-read_data[0,0]) / DeltaX) + 1
Ynum        =   read_data.shape[0] / Xnum

grid_x,grid_y = numpy.mgrid[read_data[0,0]:read_data[-1,0]:(Xnum)*1j, read_data[0,1]:read_data[-1,1]:(Ynum)*1j]
Grid_data     = griddata(read_data[:,0:2], read_data[:,2], (grid_x,grid_y), method='linear')

smooth_data   = gaussian_filter(Grid_data, sigma, mode='nearest')

for y in range(Ynum):
  for x in range(Xnum):
    if(x%inc == 0 and y%inc == 0 and not numpy.math.isnan(smooth_data[x,y])):
        print grid_x[x,y], grid_y[x,y], smooth_data[x,y]


