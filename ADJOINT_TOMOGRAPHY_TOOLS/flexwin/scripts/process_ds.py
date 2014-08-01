#!/usr/bin/env python

# find matching data and synthetics from data and syn directories, process them
# and write input to flexwin package
# similar to Carl's process_data_and_syn.pl
# Qinya Liu, Oct 2009

# remaining question: do we need -h for process_syn?
# back up your raw data and synthetics before processing
# channels may become a command-line option in the future
# when number of matching files is large, we may face issues of
# 'argument list too long' (especially when derivatives are included)
# a. data/syn_tmin_tmax are wiped out before processing
#
# need process_data_new.pl, process_syn_new.pl, pad_zeros.pl and rotate.pl
# as well as saclst sacch (?), and sac_merge

# data_list,syn_list -- input data/syn
# all_syn_list -- input syn+derivative syn (.???)
# new_data_list,new_syn_list -- new_datadir/data.c new_datadir/syn.c(hole)

import os,sys,glob,getopt

##########
def calc_snr(dfile):
  asc_file=dfile+'.a'
  error=os.system("sac2ascii.pl "+dfile+" >/dev/null")
  if (error != 0):
    sys.exit("Error converting sac file")
  t1= str(float(os.popen("saclst t1 f "+dfile).readline().split()[1])-5)
  mm_nl=os.popen("awk '$1< "+ t1 + " {print $0}' " +asc_file+ " | minmax -C").readline().split()
  if (len(mm_nl) != 4):
    sys.exit("Error computing noise level")
  max_nl=max(abs(float(mm_nl[2])),abs(float(mm_nl[3])))
  mm_sn=os.popen("minmax -C "+asc_file).readline().split()
  max_sn=max(abs(float(mm_sn[2])),abs(float(mm_sn[3])))
  data_snr=max_sn/max_nl
  os.system("rm -f ascfile")
  return(data_snr)
###########

usage= 'Usage: process_ds.py\n        -d ddir,dext,chan -s sdir,sext,chan (required)\n        -p cmt,sta -t tmin,tmax -m -w -D (optional)\n\n        -m merge sac files from the same (sta,net,cmp,hole)\n        -w generate flexwin input\n       -D allow simultaneous processing of derivative synthetics \n        -r (for regional data) -M measure_dir\n'

preprocess=False; bandpass=False; merge=False; windowing=False; der_syn=False
regional_data=False; model_name='global'; eps=0.001; pzdir=''; rot_reg='-d'; chan='LH'; synchan=chan; meas_dir='MEASURE'

# command-line options
try:
  opts,arg=getopt.getopt(sys.argv[1:],'d:s:p:t:mwDrM:')
except getopt.GetoptError:
  sys.exit(usage)
if (len(opts) == 0):
  sys.exit(usage)

for o,a in opts:
  if o == '-d':
    (datadir,dataext,chan)=a.split(',')
  if o == '-s':
    aa=a.split(',')
    syndir=aa[0]
    if (len(aa) >= 2):
      synext=aa[1]
    if (len(aa) == 3):
      synchan=aa[2]
    if (synext != ''):
      synext='.'+synext
    if (len(aa) == 2):
      synchan=chan
  if o == '-p':
    preprocess=True
    (cmt,st)=a.split(',')
    if (cmt == '' or st== '' or not os.path.isfile(cmt) or not os.path.isfile(st)):
      sys.exit('check cmt and sta file')
  if o == '-t':
    bandpass=True
    (tmin,tmax)=a.split(',')
    junk=float(tmin); junk=float(tmax) #?????
    ext='_'+tmin+'_'+tmax
    new_datadir=datadir+ext; new_syndir=syndir+ext
    if (not os.path.isdir(new_datadir)):
      os.mkdir(new_datadir)
    if (not os.path.isdir(new_syndir)):
      os.mkdir(new_syndir)
    os.system('rm -rf '+new_datadir+'/*')
    os.system('rm -rf '+new_syndir+'/*')
  if o == '-m':
    merge=True
  if o == '-w':
    windowing=True
  if o == '-D':
    der_syn=True
  if o == '-r':
    regional_data=True
    model_name='socal'
    pzdir='/opt/seismo/datalib/CALIFORNIA_PZS'; rot_reg=''
  if o == '-M':
    meas_dir=a

if (not os.path.isdir(datadir) or not os.path.isdir(syndir)):
  sys.exit(datadir+' '+syndir+' exist or not')
if merge:
  merge_dir=datadir+'/merge_dir/'
  if (not os.path.isdir(merge_dir)):
    os.mkdir(merge_dir)
if (windowing and not bandpass):
  sys.exit('pass bands -t has to be set before windowing (-w)')
if (windowing):
  meas_dir=meas_dir+ext
  input_flexwin='INPUT_FLEXWIN'+ext 
#if (preprocess and bandpass): # only clean up directories if preprocessing is done
#  os.system('rm -rf '+new_datadir+'/*')
#  os.system('rm -rf '+new_syndir+'/*')
#channel name
if chan == 'LH':
  sps='20'
else: 
  sps='1'
  resp_dir='/data2/Datalib/global_resp/'
  
#padding zeros
if (regional_data):
  bb=-40
else:
  bb=-100 
if (bandpass):
  bb=-float(tmin)*3.0
bb2=bb+1.0 # to make sure cut range is inside [b,e]
  
#########  obtain data/syn list (sort by station) ##################
nfile=0; outfile='out_ds.tmp'; os.system('rm -f '+outfile)
data_list=[]; syn_list=[]; all_syn_list=[]; hole_list=[]; comp_list=[]

print '**** Data/syn list for the same (sta,net,'+chan+'?,khole) ****'

for dd in os.popen('saclst kstnm knetwk kcmpnm khole f '+datadir+'/*'+chan+'[ZEN12]*'+dataext+"| awk '{print $2,$3,$4,$5}' | sort -k 1 | uniq").readlines():
  tmp=dd.split()
  if (len(tmp) > 3):
    [sta,net,cmp,hole]=tmp
    if abs(float(hole)+12345)<0.1:
      hole=''
  else:
    [sta,net,cmp]=tmp; hole=''
  if (list(cmp)[2] == '1'):
    scmp=cmp[0:2]+'E'
  elif (list(cmp)[2] == '2'):
    scmp=cmp[0:2]+'N'
  else:
    scmp=cmp
  if regional_data: # the naming convention is subject to change
    dat=glob.glob(datadir+'/*'+net+'.'+sta+'.'+cmp+'*'+dataext)
  else:
    dat=glob.glob(datadir+'/*'+net+'.'+sta+'.'+hole+'.'+cmp+'*'+dataext)
  syn=glob.glob(syndir+'/'+sta+'.'+net+'.'+synchan+list(scmp)[2]+synext)
  if len(syn) != 1:
   print 'no matching syn. for ', dat[0], ':', syndir+'/'+sta+'.'+net+'.'+synchan+list(scmp)[2]+synext
  else:
    if len(dat) > 1:
      print 'multiple files for merging ===', dd.rstrip()
      if merge:
        dat.sort()
        os.system('process_data_new.pl -m '+cmt+' '+' '.join(dat)+'>>'+outfile)
        if os.system('sac_merge '+' '.join(dat)+' merge_tmp.sac') != 0:
          print 'cannot merge sac files ', dat
          os.system('mv -f '+' '.join(dat[1:])+' '+merge_dir)
        else:
          print 'successfully merging sac files ', dat
          os.system('mv -f '+' '.join(dat)+' '+merge_dir+';mv -f merge_tmp.sac '+dat[0]);
      else:
        sys.exit('sac files need to be merged or resolved', ' '.join(dat), ' use -m');
    data_list.append(dat[0]); syn_list.append(syn[0]) # ENZ components
    hole_list.append(hole); comp_list.append(cmp)
    nfile=nfile+1
      
if der_syn:
  for i in range(0,len(syn_list)):
    all_syn_list.append(syn_list[i]+' '+' '.join(glob.glob(syn_list[i]+'.???')))
else:
  all_syn_list=syn_list
print "**** Total number of matching traces ", nfile, " *****"

############ preprocessing data and synthetics ################
if preprocess:
    
  print '**** pre-processing data and synthetics ****'
  data_cmd='process_data_new.pl -m '+cmt+' -p --model='+model_name
  print data_cmd+' data_files'
  error=os.system(data_cmd+' '+' '.join(data_list)+'>>'+outfile)
  if error != 0:
    sys.exit('Error pre-processing data')

  # why can't we run this separately, because of the limitation on arg length?
  f=open('syn.command','w')
  syn_cmd='process_syn_new.pl -S -m '+cmt+' -a '+st+' -p --model='+model_name
  print syn_cmd+' all_syn_files'
  f.write(syn_cmd+' '+' '.join(all_syn_list)+' >> '+outfile+'\n')
  pad_cmd='pad_zeros.pl -s -l '+str(bb)
  print pad_cmd+' all_syn_files'
  f.write(pad_cmd+' '+' '.join(all_syn_list)+' >>'+outfile+'\n')
  f.close()
  if os.system('chmod a+x syn.command; syn.command') != 0:
    sys.exit('Error pre-processing syn and padding synthetics')

########### cutting and band-passing data and synthetics #############
cdata_list=[];
if bandpass:

  print '**** cutting data and synthetics ****'
  sac_input='cuterr fillz\n' # this should be redundant
  for i in range(0,nfile):
    data=data_list[i]; syn=syn_list[i]; allsyn=all_syn_list[i].rstrip()   

    [db,de,ddt]=os.popen('saclst b e delta f '+data).readline().split()[1:]
    [sb,se,sdt]=os.popen('saclst b e delta f '+syn).readline().split()[1:]
    
    if (float(ddt)-float(sdt)) > eps:
      sys.exit('check dt for '+data+' and '+syn)
    else:
      b=max(float(db),float(sb),bb2); e=max(min(float(de),float(se)),b+30)
      if (os.system('sacch t3 '+str(b)+' t4 '+str(e)+' f '+data+' '+allsyn+'>>'+outfile)!=0):
        sys.exit('error setting start and end time to t3 and t4')
      cdata_list.append(new_datadir+'/'+os.path.basename(data));

  # further band-pass filtering data and all syn
  print '**** band-passing data and synthetics ****'
  data_cmd= 'process_data_new.pl -l t3/t4 -d '+new_datadir+' -t '+tmin+'/'+tmax+' -i '+pzdir+' -f -s '+sps+' --model='+model_name
  print data_cmd+' data_files'
  error=os.system(data_cmd+' '+' '.join(data_list)+'>>'+outfile)
  if (error != 0):
    sys.exit('Error bandpass filtering data '+str(error1))
    
  f=open('syn.command','w')
  syn_cmd='process_syn_new.pl -S -l t3/t4 -d '+new_syndir+' -t '+tmin+'/'+tmax+' -f -s '+sps+' --model='+model_name
  print syn_cmd+' all_syn_files'
  f.write(syn_cmd+' '+' '.join(all_syn_list)+'>>'+ outfile+'\n')
  f.close()
  if os.system('chmod a+x syn.command; syn.command') != 0:
    sys.exit('Error filtering syn ')

  # rotating LHE and LH1 components
  print '**** rotating processed data and synthetics ****'
  rot_cmd='rotate.pl '+rot_reg+' '+new_datadir+'/*'+chan+'[E1]*'+dataext
  print rot_cmd
  error1=os.system(rot_cmd+' >>'+outfile)

  rot_cmd='rotate.pl '+new_syndir+'/*'+synchan+'[E1]'+synext+'*'# '+new_syndir+'/*'+synchan+'[E1]'+synext+'.???'
  print rot_cmd
  error2=os.system(rot_cmd+' >>'+outfile)

  if (error1 != 0 or error2 != 0):
    sys.exit('Error rotating raw data and syn '+str(error1)+' '+str(error2))


# ############### flexwin input file ##########################
if windowing:
  if (not bandpass):
    sys.exit("make sure -t tmin,tmax is used before -w (flexwin input)")
  print '**** writing flexwin input file ****'
  string=[]
  j=-1
  flexwin_dict={}
  for i in range(0,nfile):
    # flexwin input file
    [sta,net,cmp]=os.popen('saclst kstnm knetwk kcmpnm f '+cdata_list[i]).readline().split()[1:]
    if (list(cmp)[2] == 'E' or list(cmp)[2] == '1'):
      cmp=''.join(list(cmp[0:2])+['T'])
    elif (list(cmp)[2] == 'N' or list(cmp)[2] == '2'):
      cmp=''.join(list(cmp[0:2])+['R'])
    if (not regional_data):
      new_data=glob.glob(new_datadir+'/*'+net+'.'+sta+'.'+hole_list[i]+'.'+cmp+'*'+dataext)
    else:
      new_data=glob.glob(new_datadir+'/*'+net+'.'+sta+'.'+cmp+'*'+dataext)
    new_syn=glob.glob(new_syndir+'/'+sta+'.'+net+'.'+cmp+synext)
    if (len(new_data) == 1 and len(new_syn) == 1):
      string.append(new_data[0]+'\n'+new_syn[0]+'\n'+meas_dir+'/'+sta+'.'+net+'.'+cmp+'\n')
      j=j+1
      key=sta+'.'+net+'.'+cmp
      if (key in flexwin_dict.keys()):
        print 'multiple files for key ' + key + '===='
        # compare data quality to decide data from which hole
        data1=string[flexwin_dict[key]].split('\n')[0]; data2=string[j].split('\n')[0]
        snr1=calc_snr(data1); snr2=calc_snr(data2)
        if (snr1 < snr2):
          flexwin_dict[key]=j
          print 'Keeping '+data2+' instead of '+data1+', determined by SNR', snr2, snr1
        else:
          print 'Keeping '+data1+' over '+data2+', determined by SNR', snr1, snr2 
      else:
        flexwin_dict[key]=j
    else:
      print "No exact matching for ", sta, net, cmp, hole_list[i]
# check the beginning and end times

  nkeys=len(flexwin_dict.keys())
  outstring=str(nkeys)+'\n'
  for key in sorted(flexwin_dict,key=flexwin_dict.__getitem__):
    outstring=outstring+string[flexwin_dict[key]]

  f=open(input_flexwin,'w')
  f.write(outstring)
  f.close()
  
  if (not os.path.isdir(meas_dir)):
    os.mkdir(meas_dir)

#os.system('rm -f sac.cmd syn.command')
print '**** Done ****'


