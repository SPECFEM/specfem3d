#!/usr/bin/perl

use Time::Local;
use Getopt::Std;
use POSIX;
use lib '/opt/seismo-util/lib/perl';
use SACLST;
use File::Basename;
sub Usage{
  print STDERR <<END;

  plot_adj.pl *.mtm.adj
  Example (from OUTPUT_FILES): ../scripts/plot_adj.pl *.mtm.adj

  Adapted by Qinya from plot_mtm.pl for the MT_MEASURE_ADJ package

END
  exit(1);
}

if (@ARGV == 0) { Usage(); }
if (!getopts('i:')) {die('Check input arguments\n');}
#if($opt_i){$inputdir = $opt_i;}else{$inputdir = ".";}

#@data_files=`ls -1 $inputdir/*.*.*.mtm.adj`;
@mtm_adj_files = @ARGV;

# USER MUST CHANGE THIS PATH
$progdir = "/home/lqy/SVN/mt_measure_adj/scripts";
#$progdir = "/opt/seismo-util/source/multitaper";
#$progdir = "/home/carltape/svn/specfem/mt_measure_adj_work/scripts";

my(@inputdirs);

#open(MATLAB,"> plot1.m");
# NOTE: Type "which matlab" to make sure that matlab is not alised.
open(MATLAB,"|/opt/matlab/bin/matlab -nosplash -nojvm") || die("Error opening matlab program\n");
print MATLAB "P=path; path(P,'$progdir');\n";
print MATLAB "amps=[];phs=[];dlnAs=[];dts=[];gcarcs=[];azs=[];stats=[];\n";

print MATLAB "stats=char('first line');\n";

for($i=0;$i< @mtm_adj_files;$i++){
  ($filebase) = split("mtm",`basename $mtm_adj_files[$i]`); chop($filebase); 
  $filebase2 = $filebase; $filebase = "$filebase.01";
  ($inputdir) = split(" ",`dirname $mtm_adj_files[$i]`);
  `sac2ascii.pl $inputdir/$filebase.mtm.obs.sac $inputdir/$filebase.mtm.syn.sac $inputdir/$filebase.mtm.recon_syn.sac $inputdir/$filebase.?tp.obs.sac $inputdir/$filebase.?tp.syn.sac $inputdir/$filebase.?tp.recon_syn.sac`;
  (undef,$stname,$stnetwk,$comp,$az,$gcarc,$b,$e)=split(" ",`saclst kstnm knetwk kcmpnm az gcarc b e f $inputdir/${filebase}.mtm.obs.sac`);

  print MATLAB "disp(['$filebase'])\n";
  print MATLAB "data=dlmread('$inputdir/$filebase.mtm.obs');\n";
  print MATLAB "syn=dlmread('$inputdir/$filebase.mtm.syn');\n";
  print MATLAB "new1=dlmread('$inputdir/$filebase.mtm.recon_syn');\n";
  print MATLAB "new2=dlmread('$inputdir/$filebase.ctp.recon_syn');\n";

  print MATLAB "mtm_adj = dlmread('$inputdir/$filebase2.mtm.adj');\n";
  print MATLAB "ctp_adj = dlmread('$inputdir/$filebase2.cc.adj');\n";
  print MATLAB "title = '$stname.$stnetwk.$comp \Delta = $gcarc^\\circ Az = $az^\\circ';\n";
  print MATLAB "figure;\n";
  print MATLAB "plot_adj(data,syn,new1,new2,mtm_adj,ctp_adj,$b,$e);\n";
  print MATLAB "orient tall; suptitle(title);\n";
  print MATLAB "print ('-depsc2','$inputdir/${filebase}_adj.eps');\n";
  print MATLAB "close;\n";

}


print MATLAB "\n";
close(MATLAB);
print "\n";
