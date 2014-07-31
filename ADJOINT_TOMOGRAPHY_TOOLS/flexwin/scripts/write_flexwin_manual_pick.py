#!/usr/bin/env python
# similar to Carl's prepare_meas_all.pl script
# Qinya Liu, Oct 2009, UofT

# adding capability to apply a manual pick file, Nov 2009
# manual pick file format
# input_flexwin  --- line 1
# - KIEV.IU.LHZ 1      --- delete window 1 pick for KIEV.IU.LHZ
# + UCH.KN.LHR 400 725 --- add window pick between [400,725] for UCH.KN.LHR
# possibly can be used for automation after manual picking in a GUI interface

import os,sys,glob

if (len(sys.argv) != 3 and len(sys.argv) != 4):
  sys.exit("write_flexwin_out.py measure_dirs out_filename [manual-pick-file]")

dir=sys.argv[1]; outfile=sys.argv[2]
if (not os.path.isdir(dir)):
  sys.exit('check if '+dir+' exists or not')

manual_pick=False
if (len(sys.argv) == 4):
  manual_pick=True
  manual_file=sys.argv[3]
  if (not os.path.isfile(manual_file)):
    sys.exit('Check if '+manual_file+' exists or not')
  manual_lines=open(manual_file,'r').readlines()
  flexwin_input = manual_lines[0].rstrip('\n');
  print '*** Adding manual pick info from '+ manual_file+ ' ***'
  
output=''
files=glob.glob(dir+'/*mt*')
if len(files) <= 0:
  sys.exit('No such file as '+dir+'/*mt*')
nfiles=0

for file in files:
  if (not manual_pick):
    nfiles=nfiles+1
    output=output+''.join(os.popen('cat '+file).readlines())
   
  else:  # manual pick file
    basefile='.'.join(os.path.basename(file).split('.')[0:3])
#    print 'processing '+basefile+' ...'
    out=os.popen('grep '+basefile+' '+manual_file).readlines()
    if (len(out) == 0):
      output=output+''.join(os.popen('cat '+file).readlines())
      nfiles=nfiles+1
    elif (len(out) == 1):
      out2=out[0].rstrip('\n').split(' ') # manual pick info
      mtfile=os.popen('cat '+file).readlines()
      nwin=int(mtfile[2])
      if (out2[0] == '-'):  # delete window(s)
        nwin_subtract=len(out2)-2
        print 'Deleting windows from '+file, nwin, nwin_subtract
        if (nwin > nwin_subtract): # delete some windows
          nfiles=nfiles+1
          output=output+mtfile[0]+mtfile[1]+'           '+str(nwin-nwin_subtract)+'\n'
          for i in range(0,nwin):
            if not (str(i+1) in out2[2:]):
              output=output+mtfile[i+3]
      elif (out2[0] == '+'): # add windows
        nwin_add=(len(out2)-2)/2
        nfiles=nfiles+1
        print 'Adding windows for '+file, nwin, nwin_add
        output=output+mtfile[0]+mtfile[1]+'           '+str(nwin+nwin_add)+'\n'
        output=output+''.join(mtfile[3:])
        for i in range(0,nwin_add):
          if (float(out2[i*2+2]) < float(out2[i*2+3])):
            output=output+'  '+out2[i*2+2]+'   '+out2[i*2+3]+'\n'
          else:
            print 'End time < start time in line: '+ out[0]; exit()
      else:
        print 'Format error in line: '+out[0]; exit()
    else:
      print 'Multiple lines for '+basefile; exit()

# now add files that did not have mt_input
if manual_pick:
  print ' *** Adding missing data and syn pairs *** '
  f=open(outfile, 'w')
  f.write(output)
  f.close()
  for i in range(1,len(manual_lines)):
    out2=manual_lines[i].split(' ')
    if (os.system('cat '+outfile+' | grep '+out2[1]+'>/dev/null') != 0 and out2[0] == '+'): # not an existing window
      out3=os.popen('grep -n "'+dir+'.*'+out2[1]+'" '+flexwin_input).readlines()
      if (len(out3) != 1):
        print 'Error grep '+out2[1]+' from '+flexwin_input+' : '+synfile; exit()
      else:
        line_num=int(out3[0].split(':')[0])
        datafile=os.popen("awk 'NR=="+str(line_num-2)+" {print $0}' "+flexwin_input).readline()
        synfile=os.popen("awk 'NR=="+str(line_num-1)+" {print $0}' "+flexwin_input).readline()
        print 'Adding windows to '+datafile+' and '+synfile
        nfiles=nfiles+1
        nwin_add=(len(out2)-2)/2
        output=output+datafile+synfile+str(nwin_add)+'\n'
        for i in range(0,nwin_add):
          output=output+'  '+out2[i*2+2]+'   '+out2[i*2+3]+'\n'
    
output=str(nfiles)+'\n'+output
f=open(outfile, 'w')
f.write(output)
f.close()

