#!/usr/bin/perl

use Time::Local;
use Getopt::Std;
use POSIX;
use lib '/opt/seismo-util/lib/perl';
use SACLST;
use File::Basename;
sub Usage{
  print STDERR <<END;

  plot_mtm.pl mtmfiles.dlnA
  EXAMPLE (execute from OUTPUT_FILES):
     ../scripts_meas/plot_mtm.pl *.dlnA

  original coding by Vala, adapted by Qinya for the MT_MEASURE_ADJ package

END
  exit(1);
}

if (@ARGV == 0) { Usage(); }
if (!getopts('i:')) {die('Check input arguments\n');}
#if($opt_i){$inputdir = $opt_i;}else{$inputdir = ".";}
#$tp = "mtm";

#@data_files=`ls -1 $inputdir/*.*.*.mtm.dlnA`;
@dlnA_files = @ARGV;

# USER MUST CHANGE THIS PATH
#$progdir = "/home/lqy/SVN/mt_measure_adj/scripts";
#$progdir = "/opt/seismo-util/source/multitaper";
$progdir = "/net/denali/scratch1/carltape/svn/cig/seismo/3D/ADJOINT_TOMO/measure_adj_work/scripts_meas";

my(@inputdirs);

#open(MATLAB,"> plot1.m");
# NOTE: Type "which matlab" to make sure that matlab is not alised.
open(MATLAB,"|/opt/matlab/bin/matlab -nosplash -nojvm") || die("Error opening matlab program\n");
#open(MATLAB,"|/usr/local/bin/matlab -nosplash -nojvm") || die("Error opening matlab program\n");

print MATLAB "P=path; path(P,'$progdir');\n";
print MATLAB "amps=[];phs=[];dlnAs=[];dts=[];gcarcs=[];azs=[];stats=[];\n";

print MATLAB "stats=char('first line');\n";

$iall = 1;  # 1: plot full transfer functions; 0: plot only average measurements

for ($i=0;$i< @dlnA_files;$i++) {

  ($filebase) = split("dlnA",`basename $dlnA_files[$i]`); chop($filebase);
  ($inputdir) = split(" ",`dirname $dlnA_files[$i]`);
  `sac2ascii.pl $inputdir/$filebase.obs.sac $inputdir/$filebase.syn.sac $inputdir/$filebase.recon_syn.sac $inputdir/$filebase.recon_syn_dt.sac $inputdir/$filebase.recon_syn_cc.sac $inputdir/$filebase.recon_syn_cc_dt.sac`;
  $epsfile = "$inputdir/${filebase}_measure.eps";

  (undef,$stname,$stnetwk,$comp,$az,$gcarc)=split(" ",`saclst kstnm knetwk kcmpnm az gcarc f $inputdir/${filebase}.obs.sac`);

  print MATLAB "disp(['$filebase'])\n";
  print MATLAB "data=dlmread('$inputdir/$filebase.obs');\n";
  print MATLAB "syn=dlmread('$inputdir/$filebase.syn');\n";
  print MATLAB "dataf=dlmread('$inputdir/$filebase.obs.power');\n";
  print MATLAB "synf=dlmread('$inputdir/$filebase.syn.power');\n";
  print MATLAB "new=dlmread('$inputdir/$filebase.recon_syn');\n";
  print MATLAB "new_dt=dlmread('$inputdir/$filebase.recon_syn_dt');\n";
  if ($iall == 1) {
    # read sub-sampled error curves
    print MATLAB "dlnA=dlmread('$inputdir/$filebase.err_dlnA');\n";
    print MATLAB "dt  =dlmread('$inputdir/$filebase.err_dt');\n";

    # read NON-sub-sampled error curves
    print MATLAB "dlnA_full=dlmread('$inputdir/$filebase.err_dlnA_full');\n";
    print MATLAB "dt_full  =dlmread('$inputdir/$filebase.err_dt_full');\n";

  } else {
    print MATLAB "dlnA=0; dlnA_full=0; dt=0; dt_full=0;\n";
  }

  # read frequency limits
  print MATLAB "flims=dlmread('$inputdir/$filebase.freq_limits');\n";

  # read multitaper average measurements
  print MATLAB "dlnA_ave=dlmread('$inputdir/$filebase.dlnA_average');\n";
  print MATLAB "dt_ave  =dlmread('$inputdir/$filebase.dt_average');\n";

  # read cross-correlation measurements
  print MATLAB "dlnA_cc=dlmread('$inputdir/$filebase.dlnA_cc');\n";
  print MATLAB "dt_cc  =dlmread('$inputdir/$filebase.dt_cc');\n";

  # read cross-correlation reconstruction
  print MATLAB "new_cc=dlmread('$inputdir/$filebase.recon_syn_cc');\n";
  print MATLAB "new_cc_dt=dlmread('$inputdir/$filebase.recon_syn_cc_dt');\n";

  print MATLAB "title = '$stname.$stnetwk.$comp \Delta = $gcarc^\\circ Az = $az^\\circ';\n";
  print MATLAB "figure;\n";

  # KEY COMMAND (plot_mtm.m)
  print MATLAB "plot_mtm(data,syn,dataf,synf,new,new_dt,dlnA,dt,dlnA_full,dt_full,dlnA_ave,dt_ave,new_cc,new_cc_dt,dlnA_cc,dt_cc,flims,$iall);\n";
  print MATLAB "orient tall; suptitle(title);\n";
  #print MATLAB "print ('-depsc2','$inputdir/${filebase}_measure.eps');\n";
  print MATLAB "print ('-depsc2','$epsfile');\n";
  print MATLAB "close;\n";
  #  print MATLAB "plot_mtm_quality('$inputdir');\n";

  print "\nplotted to: $epsfile \n";

}

print MATLAB "\n";
close(MATLAB);

