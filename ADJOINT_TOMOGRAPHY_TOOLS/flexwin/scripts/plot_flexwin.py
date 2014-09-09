#!/usr/bin/env python
# Similar to Carl's plot_windows_all.pl script
# Qinya Liu, Oct 2009, UofT
# plots are sort by distance and then components

import sys,os

if (len(sys.argv) != 2):
  sys.exit('plot_flexwin.py measure_dir(MEASURE)')

dir=sys.argv[1]
if not os.path.isdir(dir):
  sys.exit("no such dir "+dir)
script1='scripts/plot_seismos_gmt.sh'
script2='scripts/extract_event_windowing_stats.sh'
if not os.path.isfile(script1) or not os.path.isfile(script2):
  sys.exit('no '+script1+' or '+script2)

if (os.system(script2+' '+dir+' > /dev/null') != 0):
  sys.exit('Error executing '+ script2)

ps=dir+'/event_winstats.pdf'+' '+dir+'/event_recordsection.pdf'

for basename in os.popen('grep DDG '+dir+"/*.info | awk '{print $1,$4}'| sort -k 2 -g | awk '{print $1}' | sed 's/\.info:#//g'").readlines():
  input = basename.rstrip()  
  if (os.system(script1+' '+input) != 0):
    sys.exit("Error plotting individual seismograms")
  ps = ps + ' '+input+'.seis.pdf'

if (os.system('pdcat -r '+ps+' flexwin_seis.pdf') != 0):
  sys.exit("Error concatenate all seismogram plots")


