#!/usr/bin/perl
#
#--------------------------------------------------------------------
# prepare_adj_src.pl
#
# This script reads in a set of adjoint sources made from Z-R-T records
# and outputs a set of adjoint sources in Z-E-N that can be used in SPECFEM3D.
# The key is that it will create an all-zeros record if no measurement is made on
# a particular component (say, Z), but IS made on another component (say, R or T).
# The output STATIONS file is used within plotting routines (plot_win_adj.pl).
#
# EXAMPLES:
#   prepare_adj_src.pl -m PLOTS/CMTSOLUTION_9818433 -z BH -s STATIONS_CI -o ADJOINT_SOURCES -i OUTPUT_FILES/*adj
#   prepare_adj_src.pl -m CMTSOLUTION_9703873 -s PLOTS/STATIONS_TOMO -o ADJOINT_SOURCES OUTPUT_FILES/*recon.cc.sac
#
#--------------------------------------------------------------------

use File::Basename;
use Getopt::Std;
use lib 'UTIL/perl';
use CMT_TOOLS;
use DELAZ5;

#Reading input arguments:
@ARGV>0 || die("prepare_adj_src.pl -m CMT -z BH -s STATION -o OUTDIR -i  all_adj_RTZ_files\n");

if (!getopts('m:z:s:o:i')) {die("Check arguments\n");}

if ($opt_m) {$cmt = $opt_m;} else {$cmt = "CMTSOLUTION";}
if ($opt_z) {$chan = $opt_z;} else {$chan = "BH";}
if ($opt_s) {$stafile = $opt_s;} else {$stafile = "STATIONS";}
if ($opt_o) {$outdir = $opt_o;} else {$outdir = "ADJ_OUT";}

if (not -f $cmt) {die("Check for $cmt\n");}
if (not -f $stafile) {die("Check for $stafile\n");}
if (not -d $outdir) {die("Check for $outdir\n");}

($elat,$elon) = get_cmt_location($cmt);

# fill up the hash table based on adjoint sources
print "Processing adjoint source files ...\n";
foreach $file (@ARGV) {
  if (not -f $file) {die("Check file $file\n");}
  ($basefile) = basename($file);
  ($sta,$net,$cmp) = split(/\./,$basefile);
  $stanet = "$sta.$net";
  #$stanet = "$sta";
  $adj{$stanet}{$cmp} = $file;
  if (defined $adj{$stanet}{all}) {$adj{$stanet}{all} ++ ;}
  else {$adj{$stanet}{all} = 1;}
  print "  $basefile -- $stanet $cmp\n";
}

#die("testing");

$cmpZ = "${chan}Z";
$cmpR = "${chan}R";
$cmpT = "${chan}T";
$cmpE = "${chan}E";
$cmpN = "${chan}N";

# NOTE: channel is either BH or LH (from synthetics)
$stafile_temp = "station.tmp";
if (-f ${stafile_temp}) {system("rm -f ${stafile_temp}");}
$nstation = 0;
print "Copying station adjoint source files to $outdir ...\n";
foreach $stanet (keys(%adj)) {
  ($sta,$net) = split(/\./,$stanet);

  if (defined $adj{$stanet}{$cmpT} or defined $adj{$stanet}{$cmpR} or defined $adj{$stanet}{$cmpZ}) {
    if (defined $adj{$stanet}{$cmpT} ) {
      $tcomp = $adj{$stanet}{$cmpT}; $rcomp = $tcomp; $zcomp = $tcomp;
      $ncomp = $rcomp; $ecomp = $tcomp; $zcomp=~s/$cmpT/$cmpZ/;
      $rcomp =~s/$cmpT/$cmpR/; $ncomp=~s/$cmpT/$cmpN/; $ecomp=~s/$cmpT/$cmpE/;}
    elsif (defined $adj{$stanet}{$cmpR}) {
      $rcomp = $adj{$stanet}{$cmpR}; $tcomp = $rcomp; $zcomp = $rcomp;
      $ncomp = $rcomp; $ecomp = $rcomp; $zcomp=~s/$cmpR/$cmpZ/;
      $tcomp =~s/$cmpR/$cmpT/; $ncomp=~s/$cmpR/$cmpN/; $ecomp=~s/$cmpR/$cmpE/;}
    else {
      $zcomp = $adj{$stanet}{$cmpZ}; $tcomp = $zcomp; $rcomp = $zcomp;
      $ncomp = $rcomp; $ecomp = $rcomp; $rcomp=~s/$cmpZ/$cmpR/;
      $tcomp =~s/$cmpZ/$cmpT/; $ncomp=~s/$cmpZ/$cmpN/; $ecomp=~s/$cmpZ/$cmpE/;}

    # get station latitude and longitude from stations file
    # NOTE: There could be stations in the list with the same name
    #       but different networks.
    @lines = `grep "$sta " $stafile`;
    if (@lines == 0) {die("station not listed in $stafile\n")}
    foreach $line (@lines) {
       ($sta1,$net1,$slat1,$slon1,$selev,$sburial) = split(" ",$line);
       #print "\n -- $sta1,$net1,$slat1,$slon1 -- \n";
       if("$sta1.$net1" eq $stanet) {
          $slon = $slon1; $slat = $slat1;
          print "  *** $stanet ***\n";
          $text = sprintf("%8s %4s %.6f %.6f %.6f %.6f",$sta1,$net1,$slat1,$slon1,$selev,$sburial);
          `echo $text >> ${stafile_temp}`;
          $nstation ++ ;
       }
    }

    #(undef,undef,$slat,$slon) = split(" ",`grep "$sta " $stafile | head -n 1`);
    #print "\n $slat,$slon \n";
    #print "*** $stanet ***\n";
    #`grep "$sta " $stafile | head -n 1 >> ${stafile_temp}`; $nstation ++ ;

    if ($opt_i) {
      # get back-azimuth and apply rotation using rotate_adj_src.f90,
      # which will WRITE OUT the adjoint sources on the east and north components
      (undef,undef,undef,$baz) = delaz5($slat,$slon,$elat,$elon,0); # station - event loc in radians
      #print "***$tcomp***\nbin/rotate_adj_src $baz $zcomp $tcomp $rcomp $ecomp $ncomp\n\n";
      system("rotate_adj_src $baz $zcomp $tcomp $rcomp $ecomp $ncomp");
      if ($? != 0) {
  die("Error: rotate_adj_src $baz $tcomp $rcomp $ecomp $ncomp\n");
      }
      system("cp -f $ecomp $ncomp $zcomp $outdir");
    }

  } else {die("check hash table for adj{$stanet}{$cmpT,$cmpR,$cmpZ}\n");}
}

# write out STATIONS_ADJOINT file
system("echo $nstation > STATIONS_ADJOINT; cat ${stafile_temp} >> STATIONS_ADJOINT; rm -f ${stafile_temp}");
