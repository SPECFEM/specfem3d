#!/usr/bin/perl

use Time::Local;
use Getopt::Std;
use POSIX;
use lib '/opt/seismo-util/lib/perl';
use SACLST;
use File::Basename;
sub Usage{
  print STDERR <<END;

  plot_cc.pl ccfiles.dlnA
  Example (from OUTPUT_FILES): ../scripts/plot_cc.pl *.dlnA_cc

  original coding by Vala, adapted by Qinya for the MT_MEASURE_ADJ package

END
  exit(1);
}

if (@ARGV == 0) { Usage(); }
if (!getopts('i:')) {die('Check input arguments\n');}
#if($opt_i){$inputdir = $opt_i;}else{$inputdir = ".";}

@dlnA_files = @ARGV;

# USER MUST CHANGE THIS PATH
#$progdir = "/home/lqy/SVN/mt_measure_adj/scripts";
#$progdir = "/opt/seismo-util/source/multitaper";
$progdir = "/home/carltape/svn/specfem/mt_measure_adj_work/scripts";

my(@inputdirs);

#open(MATLAB,"> plot1.m");
# NOTE: Type "which matlab" to make sure that matlab is not alised.
open(MATLAB,"|/opt/matlab/bin/matlab -nosplash -nojvm") || die("Error opening matlab program\n");
print MATLAB "P=path; path(P,'$progdir');\n";
print MATLAB "amps=[];phs=[];dlnAs=[];dts=[];gcarcs=[];azs=[];stats=[];\n";

print MATLAB "stats=char('first line');\n";

for ($i=0;$i< @dlnA_files;$i++) {

  ($filebase) = split("dlnA",`basename $dlnA_files[$i]`); chop($filebase);
  ($inputdir) = split(" ",`dirname $dlnA_files[$i]`);
  `sac2ascii.pl $inputdir/$filebase.obs.sac $inputdir/$filebase.syn.sac $inputdir/$filebase.recon_syn_cc.sac $inputdir/$filebase.recon_syn_cc_dt.sac`; 
  $psfile = "$inputdir/${filebase}_measure_cc.eps";

  (undef,$stname,$stnetwk,$comp,$az,$gcarc)=split(" ",`saclst kstnm knetwk kcmpnm az gcarc f $inputdir/${filebase}.obs.sac`);

  print MATLAB "disp(['$filebase'])\n";
  print MATLAB "data=dlmread('$inputdir/$filebase.obs');\n";
  print MATLAB "syn=dlmread('$inputdir/$filebase.syn');\n";
  print MATLAB "new=dlmread('$inputdir/$filebase.recon_syn_cc');\n";
  print MATLAB "new_dt=dlmread('$inputdir/$filebase.recon_syn_cc_dt');\n";

  # read cross-correlation measurements
  print MATLAB "dlnA_cc=dlmread('$inputdir/$filebase.dlnA_cc');\n";
  print MATLAB "dt_cc  =dlmread('$inputdir/$filebase.dt_cc');\n";  

  print MATLAB "title = '$stname.$stnetwk.$comp \Delta = $gcarc^\\circ Az = $az^\\circ';\n";
  print MATLAB "figure;\n";

  # KEY COMMAND (plot_cc.m)
  print MATLAB "plot_cc(data,syn,new,new_dt,dlnA_cc,dt_cc);\n";
  print MATLAB "orient tall; suptitle(title);\n";
  #print MATLAB "print ('-depsc2','$inputdir/${filebase}_measure.eps');\n";
  print MATLAB "print ('-depsc2','$psfile');\n";
  print MATLAB "close;\n";
}

print MATLAB "\n";
close(MATLAB);
print "\ndone with $psfile \n\n";
