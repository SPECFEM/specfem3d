#!/usr/bin/perl -w
#
#----------------------------------------------------
# plot_win_adj_all.pl
#
# This calls plot_win_adj.pl to make plots of data, syn, windows, measurements, and adjoint sources.
#
# EXAMPLE:
#   plot_win_adj_all.pl -l -10/200 -m ../CMTSOLUTION_9818433 -n BH -b 0 -k 7/1 -a STATIONS_ADJOINT -d DATA -s SYN -c RECON -w MEASUREMENT.WINDOWS -i m16 -j 6/30
#
#----------------------------------------------------

use Getopt::Std;
use POSIX;

sub Usage{
  print STDERR <<END;

  plot_win_adj_all.pl -m CMTFILE -a STATION_FILE -n chan -b iboth -l tmin/tmax -k imeas/iadj -d data_dir -s syn_dir -c recon_dir -w winfile -i smodel -j Tmin/Tmax

END
  exit(1);
}

if (@ARGV == 0) { Usage(); }
if (!getopts('m:a:n:b:d:k:l:s:c:w:i:j:')) {die('Check input arguments\n');}
if($opt_m) {$cmtfile=$opt_m;}
if(!$opt_b) {$iboth = 0} else {$iboth = $opt_b}
#if($opt_k) {$imeas = $opt_k;}     # does not work for imeas = 1
if($opt_l) {$tcuts=$opt_l;}
if($opt_a) {$station_list=$opt_a;}
if($opt_n) {$chan=$opt_n;}
if($opt_w) {$winfile=$opt_w;}
if($opt_d) {$data_dir=$opt_d;}
if($opt_s) {$syn_dir=$opt_s;}
if($opt_c) {$recon_dir=$opt_c;}
if($opt_i) {$smodel=$opt_i;}
if($opt_j) {$Ts=$opt_j;}

print "\nplot_win_adj_all.pl:\n";

# loop over STATIONS -- first line is the number of stations
$ns = `sed -n 1p $opt_a`;
$nline = 1;
for($i = 0; $i < $ns; $i++){
  $nline = $nline+1;
  $info = `sed -n ${nline}p $opt_a`;
  chomp($info);
  ($sta,$net) = split(" ",$info);
  print "$sta $net\n";

  # include the -r command here if you want the traces scaled by the max Z-R-T value
  print "\n plot_win_adj.pl -m $cmtfile -n $sta/$net/$chan -b $iboth -l $tcuts -k $opt_k -a $station_list -d $data_dir -s $syn_dir -c $recon_dir -w $winfile -i $smodel -j $Ts\n";
  `plot_win_adj.pl -m $cmtfile -n $sta/$net/$chan -b $iboth -l $tcuts -k $opt_k -a $station_list -d $data_dir -s $syn_dir -c $recon_dir -w $winfile -i $smodel -j $Ts`;

}

# sleep command to allow the last file to convert to PDF
`sleep 1s`;
#----------------
