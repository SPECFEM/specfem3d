#!/usr/bin/perl -w

#-----------------------------------
# Carl Tape, 27-March-2009
# plot_seis_all.pl
#
# This script reads a list of desired source-station pairs
# and plots the data and synthetics by calling plot_seis.pl
#
# EXAMPLE:
#    plot_seis_all.pl -f INPUT/plot_seis_input -i 26/26 -t 6/30 m16 m99
#    plot_seis_all.pl -t 6/30 m16 m00
#    plot_seis_all.pl -t 6/30 m00 m99
#
# Using multi:
#    plot_seis_all.pl -t 6/30 m16 m04 m01 m00 m99        # 9983429 - DAN (new)
#    plot_seis_all.pl -t 6/30 m16 m08 m04 m01 m00        # 9983429 - DAN (old)
#    plot_seis_all.pl -t 6/30 m16 m08 m04 m01 m00 m99    # 9983429 - DAN (full)
#    plot_seis_all.pl -t 6/30 m16 m00 m99                # 9983429 - DAN
#    plot_seis_all.pl -t 6/30 m16                        # 14179736
#
# For GJI paper:
#    plot_seis_all.pl -f INPUT/plot_seis_input -i 1/1 -t 6/30 m16 m04 m01 m00 m99 # b
#    plot_seis_all.pl -f INPUT/plot_seis_input -i 5/5 -t 3/30 m16 m04 m01 m00 m99 # b
#    plot_seis_all.pl -f INPUT/plot_seis_input -i 8/8 -t 6/30 m16 m00 m99 # BLANK
#
#    plot_seis_all.pl -f INPUT/plot_seis_input -i 7/7 -t 6/30 m16 m99 # b
#    plot_seis_all.pl -f INPUT/plot_seis_input -i 26/26 -t 6/30 m16 # b
#    plot_seis_all.pl -f INPUT/plot_seis_input -i 26/26 -t 3/30 m16 # c
#    plot_seis_all.pl -f INPUT/plot_seis_input -i 26/26 -t 2/30 m16 # d
#    plot_seis_all.pl -f INPUT/plot_seis_input -i 26/26 -t 2/30 m99 # e
#    plot_seis_all.pl -f INPUT/plot_seis_input -i 8/8 -t 6/30 m16 m00 m99 # b
#    plot_seis_all.pl -f INPUT/plot_seis_input -i 21/21 -t 2/30 m16 m00 m99 # b
#    plot_seis_all.pl -f INPUT/plot_seis_input -i 13/13 -t 6/30 m16 m00 m99 # b
#    plot_seis_all.pl -f INPUT/plot_seis_input -i 17/17 -t 6/30 m16 m00 m99 # b
#    plot_seis_all.pl -f INPUT/plot_seis_input -i 87/87 -t 6/30 m16 # NONE
#
#    plot_seis_all.pl -f INPUT/plot_seis_input -i 5/5 -t 6/30 m16 m00 m99 # b - win
#    plot_seis_all.pl -f INPUT/plot_seis_input -i 9/9 -t 6/30 m16 m00 m99 # b - win
#    plot_seis_all.pl -f INPUT/plot_seis_input -i 83/83 -t 6/30 m16 m00 m99 # b - win
#    plot_seis_all.pl -f INPUT/plot_seis_input -i 27/27 -t 6/30 m16 m00 m99 # b - win
#    plot_seis_all.pl -f INPUT/plot_seis_input_reflection -i 2/2 -t 3/30 m16 m00 m99 # b - win
#    plot_seis_all.pl -f INPUT/plot_seis_input_reflection -i 3/3 -t 3/30 m16 # b - win
#    plot_seis_all.pl -f INPUT/plot_seis_input -i 25/25 -t 2/30 m16 m00 m99 # b - win
#    plot_seis_all.pl -f INPUT/plot_seis_input -i 3/3 -t 6/30 m16 m00 m99 # b - win
#
# Worst VR in dataset:
#    plot_seis_all.pl -f INPUT/plot_seis_input -i 88/88 -t 6/30 m16 m00
#    plot_seis_all.pl -f INPUT/plot_seis_input -i 89/89 -t 6/30 m16 m00
#
# All reflections:
#    plot_seis_all.pl -f INPUT/plot_seis_input_reflection -i 1/20 -t 3/30 m16 m00 m99
#
#-----------------------------------

use Getopt::Std;

sub Usage{
  print STDERR <<END;
  plot_seis_all.pl -f file -i imin/imax -t Tmin/Tmax smodels
END
  exit(1);
}

if (@ARGV == 0) { Usage(); }
if (!getopts('f:t:i:')) {die('Check input arguments\n');}
if ($opt_f) {$fileinput = $opt_f;}
if ($opt_i) {($imin,$imax) = split(/\//,$opt_i);}
if ($opt_t) {($Tmin,$Tmax) = split(/\//,$opt_t);}

@smodels = @ARGV;
$nmod = @smodels;

# labels for bandpass-filtered data and synthetics
$Ts = $opt_t;
$sTmin = sprintf("T%3.3i",$Tmin);
$sTmax = sprintf("T%3.3i",$Tmax);
$Ttag = "${sTmin}_${sTmax}";

$pwd = $ENV{PWD};
$imulti = 0;   # KEY COMMAND

print "-- comparing seismograms from $nmod models --> @smodels --\n";
print "-- bandpass $Ttag --\n";
print "\n$pwd\n";

#-------------------------------------------
# USER INPUT

$dir0 = "/home/carltape/SOCAL_ADJOINT";

# directory containing all CMTSOLUTION files
$slab = "16";
$dir_src = "${dir0}/results/SOURCES/socal_${slab}";
$dir_cmt = "${dir_src}/v${slab}_files";
if (not -e $dir_src) {die("check if dir_src $dir_src exist or not\n")}
if (not -e $dir_cmt) {die("check if dir_cmt $dir_cmt exist or not\n")}

# list of event IDs to use
$file_eid0 = "${dir0}/results/EID_LISTS/syn_run_m16";
if (not -f ${file_eid0}) {die("check if file_eid ${file_eid0} exist or not\n")}

# FULL stations file
$file_stations0 = "/home/carltape/gmt/stations/seismic/Matlab_output/STATIONS_CALIFORNIA_TOMO_INNER_specfem";
if (not -f ${file_stations0}) {die("check if file_stations ${file_stations0} exist or not\n")}

# data and synthetics directories
$dir_data0  = "${dir0}/DATA/FINAL";
$dir_syn0   = "${dir0}/SYN";
if (not -e ${dir_data0}) {die("check if dir_data ${dir_data0} exist or not\n")}
if (not -e ${dir_syn0}) {die("check if dir_syn ${dir_syn0} exist or not\n")}
for ($kk = 1; $kk <= $nmod; $kk++) {
   $smodel = $smodels[$kk-1];
   $dirsyn = "${dir_syn0}/model_${smodel}";
   if (not -e $dirsyn) {die("check if dirsyn $dirsyn exist or not\n")}
}

# directory containing the output for windows, measurements, adjoint sources, etc
$dir_run = "${dir0}/RUNS";
if (not -e ${dir_run}) {die("check if dir_run ${dir_run} exist or not\n")}

#-------------------------------------------

# open EID list
open(IN,"${file_eid0}"); @eids = <IN>; close(IN);
$nevent = @eids;
print "$nevent events in the list\n";

# open STATIONS file
open(IN,"${file_stations0}"); @slines = <IN>; close(IN);
$nrec = @slines - 1;
print "$nrec stations in the list\n";

#-------------------------------------------

if (0==1) {
  for ($ievent = 1; $ievent <= $nevent; $ievent++) {
    $eid = $eids[$ievent-1]; chomp($eid);
    print "\n $ievent : $eid";
  }
  die("testing");
}

if (0==1) {
  for ($irec = 1; $irec <= $nrec; $irec++) {
    ($sta,$net,$slat,$slon,undef,undef) = split(" ",$slines[$irec]);
    $station = "$sta.$net";
    print "\n $irec : $sta $net";
  }
  die("testing");
}

#-------------------------------------------

# KEY: open INPUT file
# EXAMPLE:
# 9983429    DAN CI    80/150
# 9983429    LDF CI    80/160
# 14096196   WGR CI    20/100
# 14096196   SCI2 CI   90/180
# 14096196   SMM CI    0/70

#$finput = "plot_seis_input";
#$finput = "plot_seis_input_reflection";
#$finput = "plot_seis_input_science";
#$finput = "plot_seis_input_anisotropy";

#$fileinput = "INPUT/$finput";
if (not -f $fileinput) {die("check if fileinput $fileinput exist or not\n")}
open(IN,"$fileinput"); @flines = <IN>; close(IN);
$nseis = @flines;
print "$nseis seismograms in input file\n";

for ($i = 1; $i <= $nseis; $i++) {
  ($eid,$sta,$net,$lcut) = split(" ",@flines[$i-1]);
  $stanet = "${sta}/${net}";
  print "-- $i -- $eid -- $stanet -- $lcut --\n";

}
#die("TESTING");

#=============================================
# GENERATE GMT FIGURES

#$imin = 1; $imax = $nseis;
#$imin = 1; $imax = 10;
#$imin = 1; $imax = $imin;

if($imin > $nseis) {print "imin is larger than nseis ($imin, $nseis)\n"; $imin = $nseis;}
if($imax > $nseis) {print "imax is larger than nseis ($imax, $nseis)\n"; $imax = $nseis;}

$cshfile = "plot_seis_all.csh";
open(CSH,">$cshfile");

# remove previous links
`rm DATA/* SYN/*`;

for ($i = $imin; $i <= $imax; $i++) {

  ($eid,$sta,$net,$lcut) = split(" ",@flines[$i-1]);
  $stanet = "${sta}/${net}";
  print "-- $i -- $eid -- $stanet -- $lcut --\n";

  $smodel1 = $smodels[0];
  $cmtfile = "${dir_cmt}/CMTSOLUTION_${eid}";
  $stafile = "${dir_run}/${eid}/${smodel1}/MEASURE_${Ttag}/ADJOINT_${Ttag}/STATIONS_ADJOINT";
  $stafile = "${dir_run}/${eid}/${smodel1}/MEASURE_${Ttag}/ADJOINT_${Ttag}/STATIONS_ADJOINT";
  if (not -f $cmtfile) {die("check if cmtfile $cmtfile exist or not\n")}
  if (not -f $stafile) {die("check if stafile $stafile exist or not\n")}

  for ($kk = 1; $kk <= $nmod; $kk++) {
    # check window and measurement files
    $smodel = $smodels[$kk-1];
    $mdir = "${smodel}";

    #if($smodel eq "m00") {$mdir = "${smodel}_SAVE";}   # KEY FOR NOW
    if(1==1) {
       #if($smodel eq "m16") {$mdir = "${smodel}_DAN_all";}
       #$mdir = "${smodel}_RVR_all";
       $mdir = "${smodel}_SMS_S";
       #$mdir = "${smodel}_SMM_reflection";
       #$mdir = "${smodel}_WGR_rayleigh";
       #$mdir = "${smodel}_EDW2_all";

       #if($smodel eq "m16") {$mdir = "${smodel}_SMM_reflection";}
    }
    print "$kk -- $smodel -- $mdir --\n";

    $winfile = "${dir_run}/${eid}/${mdir}/WINDOW_${Ttag}/MEASUREMENT_WINDOWS_${eid}_${Ttag}_${smodel}";
    $measfile = "${dir_run}/${eid}/${mdir}/MEASURE_${Ttag}/${eid}_${Ttag}_${smodel}_window_chi";
    #if (not -f $winfile) {die("check if winfile $winfile exist or not\n")}
    #if (not -f $measfile) {die("check if measfile $measfile exist or not\n")}
    $winfiles[$kk-1]  = $winfile;
    $measfiles[$kk-1] = $measfile;

    # link the synthetics into the temporary directory
    $syndir = "${dir_syn0}/model_${smodel}/${eid}/PROCESSED/PROCESSED_${Ttag}";
    if (not -e $syndir) {die("check if syndir $syndir exist or not\n")}
    `ln -s $syndir/*$sta*BHZ* SYN`; `ln -s $syndir/*$sta*BHR* SYN`; `ln -s $syndir/*$sta*BHT* SYN`;
  }
  print "-- @smodels --\n-- @winfiles --\n-- @measfiles --\n";
  #die("TESTING");

  # link the data into the temporary directory
  $datadir = "${dir_data0}/${eid}/PROCESSED/PROCESSED_${Ttag}";
  if (not -e $datadir) {die("check if datadir $datadir exist or not\n")}
  `ln -s $datadir/*$sta*BHZ* DATA`; `ln -s $datadir/*$sta*BHR* DATA`; `ln -s $datadir/*$sta*BHT* DATA`;
  `ln -s $datadir/*$sta*HHZ* DATA`; `ln -s $datadir/*$sta*HHR* DATA`; `ln -s $datadir/*$sta*HHT* DATA`;

  #$lcut = "0/200";  # KEY
  #$lcut = "0/150";

  if ($imulti==0) {
    # KEY: plot 1, 2, or 3 sets of ZRT seismograms with map
    #`plot_seis.pl -m $cmtfile -n $stanet -l $lcut -a $stafile -d DATA -s SYN -i ${eid}/${Ts} @smodels @winfiles @measfiles`;
    print CSH "plot_seis.pl -m $cmtfile -n $stanet -l $lcut -a $stafile -d DATA -s SYN -i ${eid}/${Ts} @smodels @winfiles @measfiles\n";
  } else {
    # KEY: plot any number of sets of ZRT seismograms (no map)
    `plot_seis_multi.pl -m $cmtfile -n $stanet -l $lcut -d DATA -s SYN -i ${eid}/${Ts} @smodels @winfiles @measfiles`;
    #print CSH "plot_seis_multi.pl -m $cmtfile -n $stanet -l $lcut -d DATA -s SYN -i ${eid}/${Ts} @smodels @winfiles @measfiles\n";
  }

}

close(CSH);
#print "csh -f $cshfile\n";
system("csh -f $cshfile");

#================================================
print "\n done with plot_seis_all.pl\n\n";
#================================================

#=================================================================
