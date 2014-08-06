#!/usr/bin/env python
# note: geometrical spreading factor needs to be significantly improved,
# right now only scale by 1/dist**exp
from numpy import *
import os,sys,commands,glob,getopt,re

# command-line options
usage="""
write_flexwin_shift_weight.py  [-C]
     [-c Z/R/T  -d B/R/L -a az_weight -r corr_weight -t ]
     -D data_dir -s measure_dir  [-o adj-dir] 
     flexwin_summary_file new_flexwin_summary_file
   -c Z/R/T component weights; -d B/R/L body/rayleigh/love wave weights
   -a azimuthal weight -r correlation weight -t apply time shift
   -C socal (regional) weight scheme in which -d pnl/rayleigh/love waves
   -D and -s gives data and measure dir prefixes
"""

try:
  opts,args=getopt.getopt(sys.argv[1:],'c:d:a:D:s:r:o:tC')
except getopt.GetoptError:
  sys.exit(usage)
if len(args) != 2:
  sys.exit(usage)
  
# default weight values
comp_z_weight=1.; comp_r_weight=1.; comp_t_weight=1.
pnl_dist_weigth=0.; rayleigh_dist_weight=0.; love_dist_weight=0.
az_exp_weigth=0.; corr_exp_weight = 0.; data_dir='data'
adj_dir='.'; write_tshift=False; write_adj_dir=False; socal=False

#opts
for o,a in opts:
  if o == '-c':
    tmp=a.split('/')
    if len(tmp)!=3:
      sys.exit('-c Z/R/T '+a)
    comp_z_weight=float(tmp[0]); comp_r_weight=float(tmp[1]); comp_t_weight=float(tmp[2])
  if o == '-d':
    tmp=a.split('/')
    if len(tmp)!=3:
      sys.exit('-d P/R/L '+a)
    pnl_dist_weight =float(tmp[0]); rayleigh_dist_weight=float(tmp[1]); love_dist_weight=float(tmp[2])
  if o == '-a':
    az_exp_weight=float(a)
  if o == '-r':
    corr_exp_weight=float(a)
  if o == '-s':
    measure_dir=a
#    if not os.path.isdir(measure_dir): # measure_dir is just a prefix
#      sys.exit('Check if '+measure_dir+' exist or not')
  if o == '-o':
    write_adj_dir=True
    adj_dir=a.rstrip('/')
    if not os.path.isdir(adj_dir):
      os.mkdir(adj_dir) 
  if o == '-t':
    write_tshift=True
  if o == '-C':
    socal=True
  if o == '-D':
    data_dir=a # data_dir is a prefix 

infile=args[0]; outfile=args[1]

if not os.path.isfile(infile):
  sys.exit('Check if input flexwin summary file '+infile+' exist or not')
          
print 'Component linear weight = ', comp_z_weight, comp_r_weight, comp_t_weight
print 'Distance exp weight = ', pnl_dist_weight, rayleigh_dist_weight, love_dist_weight
print 'Azimuth exp weight = ', az_exp_weight
print 'correlation exp weight = ', corr_exp_weight

if socal:
  ref_dist=100 # km
else:
  ref_dist=60 # degree
  vr=3.1; vl=4.1
naz=10; daz=360/naz;
eps=1.0e-2; pi=3.1416926;
naz_array=ones((naz))

alldata=commands.getoutput('grep '+data_dir+'  '+infile).split('\n')
if (len(alldata) == 0):
  sys.exit('No selected files of '+data_dir+' in '+infile)
output=commands.getoutput('saclst az f '+' '.join(alldata)+"|awk '{print $2}'").split('\n')
nfile_az=len(output)
for i in range(0,nfile_az):
  az=float(output[i])
  #azimuth counts
  k = int(divmod(az,daz)[0])
  if (k < 0 or k >= naz):
    sys.exit('Error binning azimuth')
  naz_array[k] = naz_array[k] + 1
  
print 'azimuth = '+str(naz_array)

f=open(infile,'r')
nfile=int(f.readline())
print 'nfile = '+str(nfile)
if nfile != nfile_az:
  sys.exit('Check nfile and nfile_az')
g=open(outfile,'w')
g.write(str(nfile)+'\n')

for i in range(0,nfile):
  datfile=f.readline().rstrip('\n')
  synfile=f.readline().rstrip('\n')
  g.write(datfile+'\n')
  g.write(synfile+'\n')

  # read data and synthetics, check headers
  [status1,out1]=commands.getstatusoutput('saclst kstnm knetwk kcmpnm az dist gcarc f '+datfile)
  [status2,out2]=commands.getstatusoutput('saclst kstnm knetwk kcmpnm az dist gcarc f '+synfile)
  if status1 != 0 or status2 != 0:
    sys.exit('Error taking kcmpnm az dist from '+datfile+' and '+synfile)
  [sta,net,comp,az,dist,gcarc]=out1.split()[1:7]
  az=float(az); dist=float(dist); gcarc=float(gcarc)
  [sta2,net2,comp2,az2,dist2,gcarc2]=out2.split()[1:7]
  az2=float(az2); dist2=float(dist2); gcarc2=float(gcarc2)
  if (sta != sta2 or net != net2 or comp != comp2 or abs(az-az2) >eps or abs(dist-dist2)/dist > eps):
#    sys.exit('Check if kstnm knetwk kcmpnm az dist are the same for '+datfile+' and '+synfile)
    print 'Unmatching sta, net, comp, az, dist for =====\n',out1, '\n', out2

  # componennt
  if comp[2] == 'Z':
    comp_weight=comp_z_weight
  elif comp[2] == 'R':
    comp_weight=comp_r_weight
  elif comp[2] == 'T':
    comp_weight=comp_t_weight
  else:
    sys.exit('Only deal with Z/R/T components right now')

  filename=re.sub(data_dir,measure_dir,os.path.dirname(datfile))+'/'+ sta+'.'+net+'.'+comp+'.win.qual'
  if write_adj_dir:
    g.write(adj_dir+'/'+sta+'.'+net+'.'+comp+'.adj.sac\n')

  nwin=int(f.readline())
  g.write(str(nwin)+'\n')

  for j in range(0,nwin):
    (tstart,tend)=f.readline().split()

    # distance
    if socal: 
      if comp[2] == 'T':
        dist_exp_weight=love_dist_weight
      elif nwin > 1 and j == 0:
        dist_exp_weight=pnl_dist_weight
      else:
        dist_exp_weight=rayleigh_dist_weight
      dist_weight=(dist/ref_dist)**dist_exp_weight
    else:
      if comp[2] == 'T':
        if float(tend) < dist/vl:
          dist_weight=(sin(gcarc*pi/180)/sin(ref_dist*pi/180)) ** pnl_dist_weight
        else:
          dist_weight=(sin(gcarc*pi/180)/sin(ref_dist*pi/180)) ** love_dist_weight
      else:
        if float(tend) < dist/vr:
          dist_weight=(sin(gcarc*pi/180)/sin(ref_dist*pi/180)) ** pnl_dist_weight
        else:
          dist_weight=(sin(gcarc*pi/180)/sin(ref_dist*pi/180)) ** rayleigh_dist_weight

    # azimuth
    nn=3+j;
    [istat,out]=commands.getstatusoutput("awk 'NF == 6 {print $2,$3,$4,$5}' "+filename)
    if (istat != 0  or len(out)== 0):
      sys.exit('Error with '+"awk 'NF == 6 {print $2,$3,$4,$5}' "+measure_dir+'/'+filename)

    set=False
    for line in out.split('\n'):
      [t1,t2,tt,cc]=line.split()
      if abs(float(t1)-float(tstart)) < eps and abs(float(t2)-float(tend)) < eps:
        tshift=tt; corr=cc; set=True; break
    if not set:
      sys.exit('can not set tshift for '+filename)
 
    # correlation
    k = int(divmod(az,daz)[0])
    corr_weight = float(corr)**corr_exp_weight
    
    weight=comp_weight*dist_weight*corr_weight/(naz_array[k]**az_exp_weight)
#    print k, weight, comp_weight, dist_weight, corr_weight, naz_array[k]**az_exp_weight

    if write_tshift:    
      g.write('%s  %s  %s  %f\n' % (tstart, tend, tshift, weight))
    else:
      g.write('%s  %s  %s  %f\n' % (tstart, tend, '0.0', weight))
    
#    print weight

    
f.close()
g.close()
