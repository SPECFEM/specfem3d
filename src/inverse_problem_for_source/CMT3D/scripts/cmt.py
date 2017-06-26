import os,commands,sys,math
from rotate import *
from numpy import *

# -------------------------------------------------------
# utm_center (in meter) is only used for utm coordinates (utm in [1,60]
# read cmt_par[0:npar] from cmt file
def read_cmt(cmt_file,cmt_par,npar,utm=0,nevent=1,scale_par=array([1.0e22,1.0e22,1.0e22,1.0e22,1.0e22,1.0e22,1,1,1,1,1]),utm_center=[0.,0.]):

# check inputs
  if npar not in [6,7,9,10,11]:
    print 'npar=6,7,9,10,11';   return -1
  if size(cmt_par) != npar*nevent:
    print 'Error size of cmt_par'; return -1

  if not os.path.isfile(cmt_file):
    print 'No such file '+ cmt_file;  return -2
  ref_lines=file(cmt_file,'r').readlines() # read cmt file
  if mod(len(ref_lines),13) != 0:
    print 'Error number of lines in cmt file '+cmt_file; return -2
  ne = len(ref_lines)/13
  if ne != nevent:
    print 'Error number of events'; return -2

  if utm == 0: # check utm zone
    print 'Reading cmt file '+ cmt_file+' into local R,T,P coord sys'
  elif utm == -1:
    print 'Reading cmt file '+ cmt_file+' into global X(0 lon), Y(90 lon), Z(N.P.) coord sys'
    if npar < 9:
      print 'npar >= 9 for utm = -1 (global coord)'; return -3
  elif utm <= 60 and utm >= 1:
    print 'Reading cmt file '+ cmt_file+' into local UTM coordinates '+str(utm)
  else:
    print 'Error utm zone'; return -3


  for i in range(0,nevent):
    j=i*13 # start line of i'th subevent
    mrr=float(ref_lines[j+7].rstrip('\n').split()[1])
    mtt=float(ref_lines[j+8].rstrip('\n').split()[1])
    mpp=float(ref_lines[j+9].rstrip('\n').split()[1])
    mrt=float(ref_lines[j+10].rstrip('\n').split()[1])
    mrp=float(ref_lines[j+11].rstrip('\n').split()[1])
    mtp=float(ref_lines[j+12].rstrip('\n').split()[1])
    depth=float(ref_lines[j+6].rstrip('\n').split()[1])
    lat=float(ref_lines[j+4].rstrip('\n').split()[1])
    lon=float(ref_lines[j+5].rstrip('\n').split()[1])
    tshift=float(ref_lines[j+2].rstrip('\n').split()[2])
    hdur=float(ref_lines[j+3].rstrip('\n').split()[2])

# rotate to appropriate coordinates
    if (utm >= 1 and utm <= 60): # utm
      if os.system('which sph2utm > /dev/null') != 0:
        print'Error finding sph2utm program'; return 1
      [ist,out]=commands.getstatusoutput('sph2utm '+str(lon)+' '+str(lat)+' '+str(utm))
      if ist != 0:
        print'Error running sph2utm program'; return 1
      lon=(float(out.split()[2])-utm_center[0])/1000.
      lat=(float(out.split()[5])-utm_center[1])/1000.    # E and N components (km)
    elif utm == -1: # global code, global coord [X,Y,Z]
      [loc_x,loc_y,loc_z]=rotate_loc(lon,lat,R_EARTH-depth,0,0,local2global=True)
      [mrr,mtt,mpp,mrt,mrp,mtp]=rotate_cmt(lon,lat,mrr,mtt,mpp,mrt,mrp,mtp,local2global=True)
      depth=loc_x; lon=loc_y; lat=loc_z
#    else: # utm=0, global code, local [R, T, P]
#      depth=R_EARTH-depth; lon=0.; lat=0.

    a=array([mrr,mtt,mpp,mrt,mrp,mtp,depth,lon,lat,tshift,hdur])

    if nevent == 1:
      cmt_par[0:npar]=a[0:npar]/scale_par[0:npar]
    else:
      cmt_par[0:npar,i]=a[0:npar]/scale_par[0:npar]

  return 0

# -------------------------------------------------------
# write cmt file with moment tensor and location pars cmt_par[0:npar] scaled by scale_par[0:npar]
# input cmt_par order:
# utm [1,60] -- M(rtp), depth, UTM_X, UTM_Y (km)
# utm=0 --      M(rtp), R, T, P (km)
# utm=-1 --     M(xyz), X, Y, Z (km)

def write_cmt(cmt_file,cmt_par,npar,ref_cmt_file,utm=0,nevent=1,utm_center=[0.,0.],scale_par=array([1.0e22,1.0e22,1.0e22,1.0e22,1.0e22,1.0e22,1,1,1,1,1])):

# check pars
  if npar not in [6,7,9,10,11]:
    print ' npar = 6 or 7 or 9 ';    return -1
  if size(cmt_par) != npar*nevent:
    print 'Error size of cmt_par'; return -1

  if not os.path.isfile(ref_cmt_file): # read ref cmt
    print 'no reference cmt available: '+ref_cmt_file; return -2
  else:
    ref_lines=file(ref_cmt_file,'r').readlines()
    ne=len(ref_lines)/13
    if mod(len(ref_lines),13) != 0 or ne != nevent:
      print 'Error nevent for '+ref_cmt_file; return -2

  if utm < -1 or utm > 60 or (utm==-1 and npar < 9):
    print 'Check these pararameters: utm,npar =', utm, npar; return -3

  mloc=zeros((6)); mgl=zeros((6)); loc=zeros((3));  cmt_par_unscaled=zeros((npar))
  f=open(cmt_file,'w')

  NSIG_DIGIT_CMT=6
  if nevent > 1:
    exp_largest = floor(log10(max(abs(dot(diag(scale_par[0:6]),cmt_par[0:6,:])).reshape(nevent*6))))
  else:
    exp_largest = floor(log10(max(abs(dot(diag(scale_par[0:6]),cmt_par[0:6])).reshape(nevent*6))))
  nn = exp_largest - NSIG_DIGIT_CMT
  exponent = 1.0e1 ** nn

  for i in range(0,nevent):
    j = i*13 # start of i'th subevent
    if nevent == 1:
      cmt_par_unscaled=cmt_par[:] * scale_par[0:npar]
    else:
      cmt_par_unscaled=cmt_par[:,i] * scale_par[0:npar]
    mloc=cmt_par_unscaled[0:6]
    depth=float(ref_lines[j+6].rstrip('\n').split()[1])
    lat=float(ref_lines[j+4].rstrip('\n').split()[1])
    lon=float(ref_lines[j+5].rstrip('\n').split()[1])
    tshift=float(ref_lines[j+2].rstrip('\n').split()[2])
    hdur=float(ref_lines[j+3].rstrip('\n').split()[2])

    f.write(''.join(ref_lines[j:j+2]))
    if npar > 10:
      hdur=cmt_par_unscaled[10]
    if npar > 9:
      tshift=cmt_par_unscaled[9]

    if utm > 0 and utm <= 60: # basin code
      if npar >= 9:  # lon, lat, depth
        lon=cmt_par_unscaled[7]*1000.+utm_center[0] # E (m)
        lat=cmt_par_unscaled[8]*1000.+utm_center[1] # N (m)
        if os.system('which utm2sph > /dev/null') != 0:
          print 'Error finding utm2sph program'; return -3
        [ist,out]=commands.getstatusoutput('utm2sph '+str(lon)+' '+str(lat)+' '+str(utm))
        if ist != 0:
          print 'Error running utm2sph program'; return -4
        lon=float(out.split()[2]); lat=float(out.split()[5])  #lon,lat (in degrees)
      if npar >= 7: # depth
        depth=cmt_par_unscaled[6]
    elif utm == 0: # global code, [R,T,P]
      if npar >= 9:
        [lon,lat]=cmt_par_unscaled[7:9]
#        loc=cmt_par_unscaled[6:9]
#        [x,y,z]=rotate_loc(lon,lat,loc[0],loc[1],loc[2],local2global=True)
#        mgl=rotate_cmt(lon,lat,mloc[0],mloc[1],mloc[2],mloc[3],mloc[4],mloc[5],local2global=True)
#        [lon, lat, r]=xyz2sph(x,y,z); depth=R_EARTH-r
#        mloc=rotate_cmt(lon,lat,mgl[0],mgl[1],mgl[2],mgl[3],mgl[4],mgl[5],local2global=False)
      if npar >= 7:
        depth=cmt_par_unscaled[6]
    elif utm == -1: # global code, [X,Y,Z]
      if npar >= 9:
        loc=cmt_par_unscaled[6:9]
        [lon, lat, r]=xyz2sph(loc[0],loc[1],loc[2])
        depth=R_EARTH-r
        mloc=rotate_cmt(lon,lat,mloc[0],mloc[1],mloc[2],mloc[3],mloc[4],mloc[5],local2global=False)

    mloc = floor(mloc/exponent)*exponent

    f.write('time shift:%12.4f\n' % (tshift))
    f.write('half duration:%9.4f\n' % (hdur))
    f.write('latitude:%14.4f\n' % (lat))
    f.write('longitude:%13.4f\n' % (lon))
    f.write('depth:%17.4f\n' % (depth))
    f.write("Mrr:     %14.9g\n" % (mloc[0]))
    f.write("Mtt:     %14.9g\n" % (mloc[1]))
    f.write("Mpp:     %14.9g\n" % (mloc[2]))
    f.write("Mrt:     %14.9g\n" % (mloc[3]))
    f.write("Mrp:     %14.9g\n" % (mloc[4]))
    f.write("Mtp:     %14.9g\n" % (mloc[5]))

  f.close()
  print 'Writing cmt file '+cmt_file
  return 0

# -------------------------------------------------------

def diff_cmt_files(cmt1,cmt2,npar,utm=0):

  if npar not in [6,7,9,10,11]:
    print 'Error npar: 6,7,9,10,11'
    return -1

  [istat,out]=commands.getstatusoutput('wc '+cmt1)
  if istat != 0:
    print 'Error reading '+cmt1; return 1
  nline=int(out.split()[0])
  if mod(nline,13) != 0:
    print nline, mod(nline,13)
    print 'Check number of events '+str(nline); return 2
  nevent=nline/13

  if nevent == 1:
    cpar1=zeros((npar)); cpar2=zeros((npar))
  else:
    cpar1=zeros((npar,nevent)); cpar2=zeros((npar,nevent))

  if read_cmt(cmt1,cpar1,npar,nevent=nevent,utm=utm) != 0:
    print 'Error reading cmtsolution '+ cmt1; return -2

  if read_cmt(cmt2,cpar2,npar,nevent=nevent,utm=utm) != 0:
    print 'Error reading cmtsolution '+ cmt2
    return -3

  dval=cpar1-cpar2
  for k in range(0,nevent):
    print str(k)+'th event:'
    if nevent == 1:
      print dval[:]
    else:
      print dval[:,k]

  return 0

# ------------------------------------

def get_cmt(cmt_file,info='veryshort'):

  if not os.path.isfile(cmt_file):
    sys.exit('Check if '+cmt_file+' exists or not')

  lines=file(cmt_file,'r').readlines();
  M=[]

  event_info=lines[0].split()
  (pde,year,month,day,hr,min,sec)=event_info[0:7]
  mw=event_info[10]
  evnm=''+'_'.join(event_info[12:])
  (sec,msec)=sec.split('.'); msec=str(int(msec)*10)
  evid=lines[1].split()[2]
  tshift=lines[2].split()[2]
  hdur=lines[3].split()[2]
  elat=lines[4].split()[1]
  elon=lines[5].split()[1]
  edep=lines[6].split()[1]
  for i in range(0,6):
    M.append(lines[7+i].split()[1])

  if info == 'M':
    return M
  elif info == 'location':
    return elat,elon,edep
  elif info == 'evnm':
    return evnm
  if info == 'veryshort':
    return elat,elon,edep,M
  elif info == 'short':
    return elat,elon,edep,M,evnm
  elif info == 'long':
    return elat,elon,dep,M,tshift,hdur,evnm
  elif info == 'verylong':
    return elat,elon,edep,M,mw,tshift,hdur,evnm,year,month,day,hr,min,sec,msec,evid
  elif info == 'centroid':
    return elat,elon,edep,M,hdur,evnm,year1,month1,day1,hr1,min1,sec1,msec1
  else:
    return evnm

# ==================================

def mij2dc(M):

  if os.system('which mij2dc > /dev/null') != 0:
    sys.exit('Error which mij2dc')

  mij2dc_input= '1 '+' '.join(M)+'\nEOF\n'
  [istat,out]=commands.getstatusoutput('echo "'+mij2dc_input+'" | mij2dc')
  if istat != 0:
    sys.exit('Error executing mij2dc '+' '.join(M))

  for line in out.split('\n'):
    if 'output1' in line:
      [s1,d1,r1]=line.split()[1:4]
    elif 'output2' in line:
      [s2,d2,r2]=line.split()[1:4]
    elif 'epsilon' in line:
      eps = line.split()[1]
    elif 'scalar moment' in line:
      moment=line.split()[2]
      mw=math.log10(float(moment))/1.5-10.73;

  return s1,d1,r1,s2,d2,r2,mw,moment

# ** FUTURE IMPROVEMENTS: add tshift

# ====================================
# extension name array is just a name holder, the deriatives may differ from the names they are called.
# dlocation mean different things for utm = 0, and utm in [1,60], and dlocation should be set to the same value as ddepth for utm=-1
def gen_cmt_der(cmt,npar=9,dmoment=1e22,ddepth=1,dlocation=0.05,utm=0,nevent=1,ext=['Mrr','Mtt','Mpp','Mrt','Mrp','Mtp','dep','lon','lat']):

 # ddepth and dlocation are in terms of kms
  if not os.path.isfile(cmt):
    print 'Check if '+cmt_file+' exists or not'; return -1

  if utm > 60 or utm < -1:
    print 'Check utm in [-1,60]'; return -2
#  elif utm == 0:
#    ext[6]='rrr'; ext[7]='ttt'; ext[8]='ppp'
 #   ext=['Mrr','Mtt','Mpp','Mrt','Mrp','Mtp','rrr','sss','eee']
  elif utm == -1:
    ext[0]='Mxx'; ext[1]='Myy'; ext[2]='Mzz'; ext[3]='Mxy'; ext[4]='Mxz'; ext[5]='Myz';
    ext[6]='xxx'; ext[7]='yyy'; ext[8]='zzz'
    if npar < 9:
      print 'utm == -1 requires npar >= 9'; return -2

  npar_f=9
  if nevent==1:
    cmt_par=zeros((npar_f))
  else:
    cmt_par=zeros((npar_f,nevent))

  scale_par=array([dmoment,dmoment,dmoment,dmoment,dmoment,dmoment,1.,1.,1.,1.,1.,1.])

  if read_cmt(cmt,cmt_par,npar_f,utm=utm,nevent=nevent,scale_par=scale_par) != 0:
    print 'Error reading cmt_file'; return 1

  # initialization
  cmt_par_out=zeros((6))
  for k in range(0,nevent):

    if nevent == 1:
      depth=cmt_par[6]; x=cmt_par[7];  y=cmt_par[8]
    else:
      depth=cmt_par[6,k]; x=cmt_par[7,k];  y=cmt_par[8,k]
    dnew=depth+ddepth; xnew=x+dlocation; ynew=y+dlocation # name holders again
    print depth, ddepth, dnew

    # moment derivatives
    for i in range(0,6):
      cmt_par_out=cmt_par.copy()
      if nevent == 1:
        cmt_par_out[0:6] = 0.
        cmt_par_out[i] = dmoment/scale_par[i]
      else:
        for j in range(0,nevent):
          cmt_par_out[0:6,j] = 0.
        cmt_par_out[i,k]=dmoment/scale_par[i]
      cmt_file=cmt+'_'+str(ext[i])
      if nevent > 1:
        cmt_file+='.'+('%1.1d' % k)
      if write_cmt(cmt_file,cmt_par_out,npar_f,cmt,utm=utm,nevent=nevent,scale_par=scale_par) != 0:
        print 'Error writing cmt file for M_i, '+str(i); return 2

    if npar > 6:
      # depth derivatives
      cmt_file=cmt+'_'+ext[6]
      if nevent > 1:
        cmt_file+='.'+('%1.1d' % k)
      cmt_par_out=cmt_par.copy()
      if nevent > 1:
        cmt_par_out[6,k]=dnew
      else:
        cmt_par_out[6]=dnew
      print cmt_file, cmt_par_out, npar_f
      if write_cmt(cmt_file,cmt_par_out,npar_f,cmt,utm=utm,nevent=nevent,scale_par=scale_par) != 0:
        print 'Error writing cmt file for depth, '; return 3
    if npar >= 9:
      # longitude
      cmt_file=cmt+'_'+ext[7]
      if nevent > 1:
        cmt_file+='.'+('%1.1d' % k)
      cmt_par_out=cmt_par.copy()
      if nevent > 1:
        cmt_par_out[7,k]=xnew
      else:
        cmt_par_out[7]=xnew
      if write_cmt(cmt_file,cmt_par_out,npar_f,cmt,utm=utm,nevent=nevent,scale_par=scale_par) != 0:
        print 'Error writing cmt file for longitude, '; return 4
      # latitude derivatives
      cmt_file=cmt+'_'+ext[8]
      if nevent > 1:
        cmt_file+='.'+('%1.1d' % k)
      cmt_par_out=cmt_par.copy()
      if nevent > 1:
        cmt_par_out[8,k]=ynew
      else:
        cmt_par_out[8]=ynew
      if write_cmt(cmt_file,cmt_par_out,npar_f,cmt,utm=utm,nevent=nevent,scale_par=scale_par) != 0:
        print 'Error writing cmt file for latitude, '; return 5
  return 0

#-------------------------------
def sem_der_script(cmt,run_script,npar,out_dir,utm=0,nevent=1,ext=['Mrr','Mtt','Mpp','Mrt','Mrp','Mtp','dep','lon','lat'],check_sem_finish=True):

  # check existence of cmt files
  cmt=os.path.abspath(cmt)
  out_dir=os.path.abspath(out_dir)
  if not os.path.isfile(run_script):
    print 'check if '+run_script+' exists or not'; return -1
  if not os.path.isdir(out_dir):
    os.system('mkdir -p '+out_dir)

  if utm > 60 or utm < -1:
    print 'Check utm in [-1,60]'; return -2
  elif utm == 0:
    ext=['Mrr','Mtt','Mpp','Mrt','Mrp','Mtp','rrr','ttt','ppp']
  elif utm == -1:
    ext=['Mxx','Myy','Mzz','Mxy','Mxz','Myz','xxx','yyy','zzz']
    if npar < 9:
      print 'utm == -1 requires npar >= 9'; return -2
  if npar > 9:
    print 'SEM simulations only need to be done for npar up to 9'; return -2

  # check existence of derivative files
  for i in range(0,nevent):
    for k in range(0,npar):
      if nevent > 1:
        cmt_new=cmt+'_'+ext[k]+'.'+('%1.1d' % i)
      else:
        cmt_new=cmt+'_'+ext[k]
      if not os.path.isfile(cmt_new):
        print 'Error finding '+cmt_new; return 1

  out=file(run_script).readlines()
  lines=''.join(out[0:17])+'\n'
  if nevent > 1:
    lines+='for i in `seq 0 '+str(nevent-1)+'`; do \n'
    ext1='${ext}.$i'
  else:
    ext1='${ext}'
  cmt_new=cmt+'_'+ext1
  out_new=out_dir+'/'+ext1
  lines+='for ext in '+' '.join(ext[0:npar])+'; do\n'
  lines+='  cp '+cmt_new+' DATA/CMTSOLUTION\n  mkdir -p '+out_new+'\n'
  lines+='  '+out[17]+'  '+' '.join(out[18].split()[0:6])+' '+out_new+'; '+' '.join(out[18].split()[7:9])+'\n'
  lines+='  '+out[19]
  lines+='  '+' '.join(out[20].split()[0:3])+' '+out_new+'\n'
  lines+='  cp -f output_solver.txt '+out_new+'\n'
  lines+='  '+out[21]

  if check_sem_finish:
    lines+='  /scratch/liuqy/aix_bin/check_sem_finish.py OUTPUT_FILES/output_solver.txt\n  if [ $? != 0 ]; then echo "Error executing job '+ext1+' "; exit; fi\n'
  if nevent > 1:
    lines+='done\n'
  lines+='done\n'
  if len(out) > 22:
    lines+=''.join(out[22:])

  file(run_script,'w').write(lines)
  print '  finished writing run script '+run_script+'\n'
  return 0


def sem_der_ext(utm,ext):
  if utm == 0:
    ext[0]='Mrr'; ext[1]='Mtt'; ext[2]='Mpp'; ext[3]='Mrt'; ext[4]='Mrp';
    ext[5]='Mtp'; ext[6]='rrr'; ext[7]='ttt'; ext[8]='ppp'
  elif utm == -1:
    ext[0]='Mxx'; ext[1]='Myy'; ext[2]='Mzz'; ext[3]='Mxy'; ext[4]='Mxz';
    ext[5]='Myz'; ext[6]='xxx'; ext[7]='yyy'; ext[8]='zzz'
  elif utm <= 60 and utm >= 1:
    ext[0]='Mrr'; ext[1]='Mtt'; ext[2]='Mpp'; ext[3]='Mrt'; ext[4]='Mrp';
    ext[5]='Mtp'; ext[6]='dep'; ext[7]='lon'; ext[8]='lat'
  else:
    print 'Error utm [-1,60]'; return -1
  return 0
