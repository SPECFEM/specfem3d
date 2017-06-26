#!/scinet/gpc/tools/Python/Python262/bin/python -u

# find matching data and synthetics from data and syn directories, process them
# and write input to flexwin package
# similar to Carl's process_data_and_syn.pl
# Qinya Liu, Oct 2009

# remaining question: do we need -h for process_syn?
# back up your raw data and synthetics before processing
# channels may become a command-line option in the future
# when number of matching files is large, we may face issues of
# 'argument list too long' (especially when derivatives are includede)
#
# need process_data_new.pl, process_syn_new.pl, pad_zeros.pl and rotate.pl as well as sac_merge

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

usage= 'Processing synthetic data and synthetics (including derivative syn)\nUsage: process_ds_syn.py\n        -d ddir,dext -s sdir,sext (required)\n        -p cmt,sta -t tmin,tmax -m -w -D -r(optional, merge not needed)\n\n        -w generate flexwin input\n        -D allow simultaneous processing of derivative synthetics \n        -r regional data [BH/LH] components\n'

preprocess=False; bandpass=False; merge=False; windowing=False; der_syn=False
regional_data=False
# command-line options
try:
  opts,arg=getopt.getopt(sys.argv[1:],'d:s:p:t:mwDr')
except getopt.GetoptError:
  sys.exit(usage)
if (len(opts) == 0):
  sys.exit(usage)

for o,a in opts:
  if o == '-d':
    (datadir,dataext)=a.split(',')
  if o == '-s':
    (syndir,synext)=a.split(',')
    if (synext != ''):
      synext='.'+synext
  if o == '-p':
    preprocess=True
    (cmt,st)=a.split(',')
    if (cmt == '' or st== '' or not os.path.isfile(cmt) or not os.path.isfile(st)):
      sys.exit('check cmt and sta file')
  if o == '-t':
    bandpass=True
    (tmin,tmax)=a.split(',')
    junk=float(tmin); junk=float(tmax)
    ext='_'+tmin+'_'+tmax
    new_datadir=datadir+ext; new_syndir=syndir+ext
    if (not os.path.isdir(new_datadir)):
      os.mkdir(new_datadir)
    else:
      os.system('rm -rf '+new_datadir+'/*')
    if (not os.path.isdir(new_syndir)):
      os.mkdir(new_syndir)
    else:
      os.system('rm -rf '+new_syndir+'/*')
  if o == '-m':
    merge=True
  if o == '-w':
    windowing=True
  if o == '-D':
    der_syn=True
  if o == '-r':
    regional_data=True

if (not os.path.isdir(datadir) or not os.path.isdir(syndir)):
  sys.exit(datadir+' '+syndir+' exist or not')

# clean tmp files
os.system('rm -f data.tmp '+datadir+'/*.c '+syndir+'/*.c '+syndir+'/*.c??')

#channel name
if regional_data:
  chan='BH'; sps='20'
else:
  chan='LH'; sps='1.0'
#padding zeros
bb=-40
if (bandpass):
  bb=-float(tmin)*3.0

#########  obtain data/syn list##################
nsta=0
data_list=[]; syn_list=[]; all_syn_list=[]; hole_list=[]; comp_list=[]

print '**** Data/syn list for the same (sta,net,LH?,khole) ****'

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
  dat=glob.glob(datadir+'/'+sta+'.'+net+'.'+cmp+'*'+dataext)
  syn=glob.glob(syndir+'/'+sta+'.'+net+'.'+scmp+synext)
  if len(syn) != 1:
    print 'no matching syn. for ', dat[0]
  else:
    if len(dat) > 1:
      print 'multiple files for merging ===', dd.rstrip()
      sys.exit()
    data_list.append(dat[0]); syn_list.append(syn[0])
    hole_list.append(hole); comp_list.append(cmp)
    nsta=nsta+1

if der_syn:
  for i in range(0,len(syn_list)):
    all_syn_list.append(syn_list[i]+' '+' '.join(glob.glob(syn_list[i]+'.???')))
else:
  all_syn_list=syn_list
print "**** Total number of matching traces ", nsta, " *****"

############ preprocessing data and synthetics ################
if preprocess:

  print '**** pre-processing data and synthetics ****'
  error1=os.system('process_syn_new.pl  -S -m '+cmt+' -p -s '+sps+' '+' '.join(data_list)+'>/dev/null')
#  print 'process_syn_new.pl -S -m '+cmt+' -p -s '+sps+' '+' '.join(data_list)
  if error1 != 0:
    sys.exit('Error pre-processing data')

  f=open('syn.command','w')
  f.write('process_syn_new.pl  -S -m '+cmt+' -a '+st+' -p -s '+sps+' '+' '.join(all_syn_list)+' >/dev/null\n')
  f.write('pad_zeros.pl -s -l '+str(bb)+' '+' '.join(all_syn_list)+'>/dev/null\n')
  f.close()
  if os.system('chmod a+x syn.command; syn.command') != 0:
    sys.exit('Error pre-processing syn and padding synthetics')

########### cutting and band-passing data and synthetics #############
new_data_list=[]; new_syn_list=[]
cdata_list=''; csyn_list=''; rdata_list=''; rsyn_list=''

if bandpass:

  print '**** cutting data and synthetics ****'
  sac_input='cuterr fillz\n'
  for i in range(0,nsta):
    data=data_list[i]; syn=syn_list[i]; allsyn=all_syn_list[i].rstrip()
    new_data_list.append(new_datadir+'/'+os.path.basename(data)+'.c')
    new_syn_list.append(new_syndir+'/'+os.path.basename(syn)+'.c')

    [db,de,ddt]=os.popen('saclst b e delta f '+data).readline().split()[1:]
    [sb,se,sdt]=os.popen('saclst b e delta f '+syn).readline().split()[1:]
    if (float(ddt)-float(sdt)) > 0.001:
      sys.exit('check dt for '+data+' and '+syn)
    else:
#      b=min(float(db),float(sb),bb); e=max(float(de),float(se))
      b=max(float(db),float(sb),bb); e=max(min(float(de),float(se)),b+30)
      sac_input=sac_input + 'cut '+str(b)+' '+str(e)+'\n' \
                 +'r '+data+' '+allsyn+'\n cut off\n w ' \
                 + data_list[i]+'.c '+'.c '.join(allsyn.split(' ')+[''])+'\n'
      cdata_list=cdata_list+data_list[i]+'.c '
      csyn_list=csyn_list+' '+'.c '.join(allsyn.split(' ')+[''])
      if (list(comp_list[i])[2] == 'E' or list(comp_list[i])[2] == '1'):
        rdata_list=rdata_list+new_data_list[i]+' '
        rsyn_list=rsyn_list+new_syn_list[i]+' '
        if der_syn:
          rsyn_list=rsyn_list+' '+new_syndir+'/'+os.path.basename(syn)+'.???.c '
  sac_input=sac_input+'quit\n'
  f = open('sac.cmd', "w"); f.write('#!/bin/bash\nsac<<EOF\n')
  f.write(''.join(sac_input)); f.write('EOF\n'); f.close()
  if (os.system('chmod a+rx sac.cmd; sac.cmd') != 0):
    sys.exit('error cutting data and synthetics ')

  # further band-pass filtering data and all syn
  print '**** band-passing data and synthetics ****'

  error1=os.system('process_syn_new.pl -S -d '+new_datadir+' -t '+tmin+'/'+tmax+' '+cdata_list+'>/dev/null')
#  print 'process_syn_new.pl -S -d '+new_datadir+' -t '+tmin+'/'+tmax+' '+cdata_list+'\n'
  if (error1 != 0):
    sys.exit('Error bandpass filtering data and syn '+str(error1))
  f=open('syn.command','w')
  f.write('process_syn_new.pl -S -d '+new_syndir+' -t'+tmin+'/'+tmax+' '+csyn_list+'>/dev/null\n')
  f.close()
  if os.system('chmod a+x syn.command; syn.command') != 0:
    sys.exit('Error filtering syn ')

  # rotating LHE and LH1 components
  print '**** rotating processed data and synthetics ****'
  error1=os.system('rotate.pl '+rdata_list+' >/dev/null ')
  error2=os.system('rotate.pl '+rsyn_list+' >/dev/null')
  if (error1 != 0 or error2 != 0):
    sys.exit('Error rotating raw data and syn '+str(error1)+' '+str(error2))

  # correct derivative extension
  if (der_syn):
    for ext in ['Mrr','Mtt','Mpp','Mrt','Mrp','Mtp','dep','lat','lon','ttt','rrr','ppp']:
      for file in glob.glob(new_syndir+'/*'+ext+'.c')+glob.glob(new_syndir+'/*'+ext+'.c??'):
        file_list=os.path.basename(file).split('.')
        os.rename(file,new_syndir+'/'+'.'.join(file_list[:-2]+[file_list[-1],file_list[-2]]))

# ############### flexwin input file ##########################
if windowing:
  if (not bandpass):
    sys.exit("use -t tmin,tmax before -w (flexwin input)")
  print '**** writing flexwin input file ****'
  string=[]
  j=-1
  flexwin_dict={}
  for i in range(0,nsta):
    # flexwin input file
    [sta,net,cmp]=os.popen('saclst kstnm knetwk kcmpnm f '+new_data_list[i]).readline().split()[1:]
    if (list(cmp)[2] == 'E' or list(cmp)[2] == '1'):
      cmp=''.join(list(cmp[0:2])+['T'])
    elif (list(cmp)[2] == 'N' or list(cmp)[2] == '2'):
      cmp=''.join(list(cmp[0:2])+['R'])
    new_data=glob.glob(new_datadir+'/'+sta+'.'+net+'.'+cmp+'*'+dataext+'.c')
    new_syn=glob.glob(new_syndir+'/'+sta+'.'+net+'.'+cmp+synext+'.c')
    if (len(new_data) == 1 and len(new_syn) == 1):
      string.append(new_data[0]+'\n'+new_syn[0]+'\n'+'MEASURE/'+sta+'.'+net+'.'+cmp+'\n')
      j=j+1
      key=sta+'.'+net+'.'+cmp
      flexwin_dict[key]=j
    else:
      print "No exact matching for ", sta, net, cmp, hole_list[i]


  nkeys=len(flexwin_dict.keys())
  outstring=str(nkeys)+'\n'
  for key in sorted(flexwin_dict,key=flexwin_dict.__getitem__):
    outstring=outstring+string[flexwin_dict[key]]

  f=open('INPUT_flexwin','w')
  f.write(outstring.rstrip('\n'))
  f.close()

  if (not os.path.isdir('MEASURE')):
    os.mkdir('MEASURE')

os.system('rm -f data.tmp sac.cmd '+datadir+'/*.c '+syndir+'/*.c '+syndir+'/*.c?? syn.command')
print '**** Done ****'


