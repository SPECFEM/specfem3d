#!/usr/bin/env python
# similar to Carl's prepare_meas_all.pl script
# Qinya Liu, Oct 2009, UofT

import os,sys,glob

if (len(sys.argv) < 3):
  sys.exit("write_flexwin_out.py measure_dirs out_filename ")

dirs=sys.argv[1:-1]; outfile=sys.argv[-1]
f=open(outfile, 'w')

nfiles=0
output=''
for dir in dirs:
  if (not os.path.isdir(dir)):
    sys.exit('check if '+dir+' exists or not')
    
  files=glob.glob(dir+'/*mt*')
  if len(files) <= 0:
    sys.exit('No such file as '+dir+'/*mt*')

  for file in files:
    nfiles=nfiles+1
    output=output+''.join(os.popen('cat '+file).readlines())


output=str(nfiles)+'\n'+output
f.write(output)
f.close()

