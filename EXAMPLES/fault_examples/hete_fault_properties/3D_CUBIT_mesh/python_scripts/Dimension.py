import sys
import numpy

for i in range(len(sys.argv)):                        # read parameter from shell(contain the prefix)
  if(sys.argv[i].find("-p")==0):
    P1=sys.argv[i+1]
input_file = P1

read_data = numpy.loadtxt(input_file,skiprows=0)
Total_num = read_data.shape[0]
DeltaX    = read_data[1,0] - read_data[0,0]
Xnum      = int(abs(read_data[-1,0]-read_data[0,0]) / DeltaX) + 1
Ynum      = Total_num / Xnum
if(Xnum*Ynum != Total_num):
    Xnum = Xnum + 1
    Ynum = Total_num / Xnum
if(Xnum*Ynum != Total_num):
    Xnum = Xnum - 2
    Ynum = Total_num / Xnum

print Xnum, Ynum

