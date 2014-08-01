#!/usr/bin/perl

use Time::Local;
use Getopt::Std;
use POSIX;
use lib '/opt/seismo-util/lib/perl';
use SACLST;
use SAC_TOOLS;

sub Usage{
  print STDERR <<END;

  mtm.pl -s syn_dir,yn_suffix -o outputdir -w wave -m ntapers/waterlevel/npi -V [-B -C] sac_files

       -s syn_dir,syn_suffix -- plot synthetics corresponding to datafiles
                           with names: syn_dir/sta.netw.comp.syn_suffix
       -w wave -- wave is one of r, l, v1/v2, or t?/dt1/dt2 
                           (r:v1=3.4,v2=3.8 l: v1=4.2, v2=4.6)
                      the last option cuts from header t?-dt1 to t?+dt2

       -m ntapers/waterlevel/npi -- multitaper parameters [5/0.02/2.5]
       -V -- verbose
       -B -- make boxcar taper measurements
       -C -- make cosine taper measurements

  originally written by Vala, adapted by Qinya for the MT_MEASURE_ADJ package

END
  exit(1);
}

if (@ARGV == 0) { Usage(); }
if (!getopts('s:w:o:m:VBC')) {die('Check input arguments\n');}
if (@ARGV == 0) { Usage(); }

if($opt_s){
  ($syn_dir,$syn_suff)=split(/\,/,$opt_s);
  (-d $syn_dir) or die "No such directory: $syn_dir\n";
  $plot_syn=1;}

if(!$opt_w){die "You have to specify a wave-velocity \n";}
else{$wave = $opt_w;}
if($wave =~ /^r/){$vel = 3.6;$vel1=3.4;$vel2=3.8;$pick=0;}
elsif($wave =~ /^l/){$vel = 4.4;$vel1=4.2;$vel2=4.6;$pick=0;}
elsif($wave =~/^\d/){($vel1,$vel2)=split(/\//,$wave);$pick=0;}
elsif($wave =~/^t/){($head,$dt1,$dt2)=split(/\//,$wave);$pick=1;}
else{die "Don't recognize wave $wave, use r or l \n";}
if($pick==1){print "cutting from $head-$dt1 to $head+$dt2\n";}
else{print "cutting using v1 $v1 and v2 $v2\n";}

if ($opt_o) {$outputdir = $opt_o;} else {$outputdir = "OUTPUT_FILES/";}
if (not -d $outputdir) {system("mkdir $outputdir");}

if($opt_m){($ntapers,$waterlevel,$npi)=split(/\//,$opt_m);}
else{$ntapers=5;$waterlevel=0.02;$npi=2.5;}

if($opt_V){$verbose=1; print "Verbose is on\n";}
else{$verbose=0;}

if($opt_B) {$taper = 3; $adj_src = 2;} elsif ($opt_C) {$taper = 2; $adj_src = 2;}
else {$taper = 1; $adj_src = 1;}


$mt_dir = "/home/lqy/SVN/mt_measure_adj";
#print "using ntapers: $ntapers,waterlevel:$waterlevel,npi:$npi\n";

# prepare data and matching synthetics
@data_files=@ARGV;
foreach $data_file (@data_files){
  chomp($data_file);
  if (not -f $data_file) {die("Check if $data_file exist or not\n");}
}
@gcarcs=sac_files_get_values("gcarc",@data_files);
@az=sac_files_get_values("az",@data_files);
@stnames=sac_files_get_values("kstnm",@data_files);
@stnetwk=sac_files_get_values("knetwk",@data_files);
@comps=sac_files_get_values("kcmpnm",@data_files);

# find matching syn files
($data_file, $syn_file) = match_data_and_syn("${syn_dir},${syn_suff}","@data_files");
(@used_data_files) = split(" ",$data_file);
(@used_syn_files) = split(" ",$syn_file);

$nsyn=@used_syn_files;
print "Total: $nsyn \n";

# prepare input files
open(MTM,"> MEASUREMENT.WINDOWS");
print MTM "$nsyn\n";

for($i=0;$i<=$#used_data_files;$i++){
  $datafile=$used_data_files[$i];
  $synfile =$used_syn_files[$i];
  $gcarc=@gcarcs[$i];
  if($pick==1){
    $temp=`saclst $head f ${datafile}`;
    chomp($temp);
    (undef,$t)=split(/ +/,$temp);
#    print "header is: $t \n";
    $cut1=$t-$dt1;
    $cut2=$t+$dt2;
  }
  else{
    $cut1 = ($gcarc*111/$vel2);
    $cut2 = ($gcarc*111/$vel1);
    }
  $cut1f = sprintf("%8.3f",$cut1);
  $cut2f = sprintf("%8.3f",$cut2);

#run multitaper
#  print "gcarc: $gcarc cut: $cut1 $cut2 \n";

  print MTM "$datafile\n";
  print MTM "$synfile\n";
  print MTM "1\n";  # only one window per trace
  print MTM "$cut1f  $cut2f\n";

#  #saving outputfiles
#  $filebase="$stnames[$i].$stnetwk[$i].$comps[$i].mtm";
#  print "saving $outputdir/$filebase.* \n";
#  `echo $stnames[$i] $stnetwk[$i] $comps[$i] $az[$i] $gcarcs[$i] $data_files[$i] $syn_files[$i] >  $outputdir/$filebase.am  `;
#  `echo $stnames[$i] $stnetwk[$i] $comps[$i] $az[$i] $gcarcs[$i] $data_files[$i] $syn_files[$i] >  $outputdir/$filebase.dlnA`;
#  `echo $stnames[$i] $stnetwk[$i] $comps[$i] $az[$i] $gcarcs[$i] $data_files[$i] $syn_files[$i] >  $outputdir/$filebase.ph  `;
#  `echo $stnames[$i] $stnetwk[$i] $comps[$i] $az[$i] $gcarcs[$i] $data_files[$i] $syn_files[$i] >  $outputdir/$filebase.dt  `;

#  `paste trans.am.mtm     err.am   | awk '{print \$1, \$2, \$4}' >> $outputdir/$filebase.am  `;
#  `paste trans.dlnA.mtm   err.dlnA | awk '{print \$1, \$2, \$4}' >> $outputdir/$filebase.dlnA`;
#  `paste trans.ph.cor.mtm err.phi  | awk '{print \$1, \$2, \$4}' >> $outputdir/$filebase.ph  `;
#  `paste trans.dt.mtm     err.dt   | awk '{print \$1, \$2, \$4}' >> $outputdir/$filebase.dt  `;

#  `mv syn.win $outputdir/$filebase.syn.sac`;
#  `mv obs.win $outputdir/$filebase.obs.sac`;
#  `mv obs.new.sac.mtm $outputdir/$filebase.new.sac`;
#  `mv rh_quality.dat $outputdir/$filebase.quality`;
}
 close(MTM);

open(MTM,">MEASUREMENT.PAR");
print MTM "$taper\n$ntapers $waterlevel $npi\n$adj_src\n$outputdir\n.true.\n.true.\n";
close(MTM);

# ** run the program **
if ($verbose)  {system("$mt_dir/mt_measure_adj");}
else {system("$mt_dir/mt_measrue_adj > log.mt_measure");}

