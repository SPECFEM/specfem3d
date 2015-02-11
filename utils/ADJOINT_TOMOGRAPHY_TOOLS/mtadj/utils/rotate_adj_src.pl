#!/usr/bin/perl -w

use File::Basename;
use Getopt::Std;
use lib '/opt/seismo/lib/perl';
use CMT_TOOLS;
use DELAZ5;

# works on chanel name given by opt_n
# note that rotate_adj_src works only on sac files
# and the at end of this script, all adjoint sources are converted to ASCII.
# this used to be prepare_adj_src.pl in the old version.
# $bin has to be changed if not in the standard PATH

#Reading input arguments:
@ARGV>0 || die("rotate_adj_src.pl  -m CMT -s STATION -o OUTDIR [-n BH/LH] adj-files(*H[RTZ]*adj*)\n");

if (!getopts('m:s:o:n:')) {die("Check arguments\n");}

if ($opt_m) {$cmt = $opt_m;} else {$cmt = "CMTSOLUTION";}
if ($opt_s) {$stafile = $opt_s;} else {$stafile = "STATIONS";}
if ($opt_o) {$outdir = $opt_o;} else {$outdir = "ADJ_SRC";}
if (not -d $outdir) {mkdir("$outdir", 0777);}
if ($opt_n) {$name = $opt_n;} else {$name="BH";}
$tname="${name}T"; $rname="${name}R"; $zname="${name}Z";
$ename="${name}E"; $nname="${name}N";

system("rm -f $outdir/*");
if (not -f $cmt or not -f $stafile or not -d $outdir) {die("Check files/dirs\n");}

($elat,$elon) = get_cmt_location($cmt);
$bin="";

# fill up the hash table that stores the adj src existence information
foreach $file (@ARGV) {
  if (not -f $file) {die("Check file $file\n");}
  ($basefile) = basename($file);
  ($sta,$net,$cmp) = split(/\./,$basefile);
  $stanet="$sta.$net";
  $adj{$stanet}{$cmp} = $file;
  if (defined $adj{$stanet}{all}) {$adj{$stanet}{all} ++ ;}
  else {$adj{$stanet}{all} = 1;}
#  print "$basefile -- $stanet\n";
}

$stafile_temp="station.tmp"; system("rm -f $stafile_temp");
$nstation = 0;
print "\nCopying stations adjoint source file to: $outdir ...\n";
foreach $stanet (keys(%adj)) {
  ($sta,$net) = split(/\./,$stanet);
  if (defined $adj{$stanet}{$tname} or defined $adj{$stanet}{$rname} or defined $adj{$stanet}{$zname}) {
    if (defined  $adj{$stanet}{$tname} ) {
      $tcomp = $adj{$stanet}{$tname}; $rcomp = $tcomp; $zcomp = $tcomp;
      $ncomp = $rcomp; $ecomp = $tcomp; $zcomp=~s/$tname/$zname/;
      $rcomp =~s/$tname/$rname/; $ncomp=~s/$tname/$nname/;
      $ecomp=~s/$tname/$ename/;}
    elsif (defined $adj{$stanet}{$rname}) {
      $rcomp = $adj{$stanet}{$rname}; $tcomp = $rcomp; $zcomp = $rcomp;
      $ncomp = $rcomp; $ecomp = $rcomp; $zcomp=~s/$rname/$zname/;
      $tcomp =~s/$rname/$tname/; $ncomp=~s/$rname/$nname/;
      $ecomp=~s/$rname/$ename/;}
    else {
      $zcomp = $adj{$stanet}{$zname}; $tcomp = $zcomp; $rcomp = $zcomp;
      $ncomp = $rcomp; $ecomp = $rcomp; $rcomp=~s/$zname/$rname/;
      $tcomp =~s/$zname/$tname/; $ncomp=~s/$zname/$nname/;
      $ecomp=~s/$zname/$ename/;}

    # write station info
    (undef,undef,$slat,$slon) = split(" ",`grep "$sta .*$net" $stafile | head -n 1`);
    print "*** $stanet: ($slat, $slon) ***\n";
    `grep "$sta .*$net" $stafile | head -n 1 >> $stafile_temp`; $nstation ++ ;

    # station/event locs in degrees (input i = 0) and output baz in radian
    (undef,undef,undef,$baz) = delaz5($slat,$slon,$elat,$elon,0);
    print "rotate_adj_src $baz $zcomp $tcomp $rcomp $ecomp $ncomp\n";
    system("${bin}rotate_adj_src $baz $zcomp $tcomp $rcomp $ecomp $ncomp");
    if ($? != 0) {die("Error in rotate_adj_src \n");}
    system("cp -f $ecomp $ncomp $zcomp $outdir");
  }
}
system("sac2ascii.pl $outdir/*");
system("echo $nstation > STATIONS_ADJOINT; cat $stafile_temp >> STATIONS_ADJOINT; mv STATIONS_ADJOINT $outdir/");

