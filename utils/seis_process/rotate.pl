#!/usr/bin/perl

use Getopt::Std;
use POSIX;

sub Usage {
  print STDERR <<EOF;

 Usage: rotate.pl -lStart -LEnd -c -d  [East|1]_Component_file_name

    rotate.pl rotates iris data(-d), trinet data and synthetics (with
    out -d option)
    -l -L specifies the start and end point of the traces to rotate
    -c check SAC output on screen
       ex. rotate.pl -d *.LHE.SAC for iris data
           rotate.pl PAS.*.MXE.sac for synthetics
    For iris data (with -d), the timing part of the name will be ignored
    and only the component part of the name will be changed correspondingly

    Make sure you have sac, saclst in the PATH before execution

    Qinya Liu, Caltech, May 2007

EOF
exit(1)
}


if (!getopts('l:L:cdts')) {die('Check input arguments\n');}
@ARGV > 0 or Usage();
if (!$opt_l) {$opt_l = 0;}
$saclst = "/opt/seismo-util/bin/saclst";
if (not -f $saclst) {die("No such file as $saclst\n");}
$undef=-12345.0;
$eps=0.1;


foreach $file (@ARGV) {
  print "processing $file\n";
  if (! -f $file) {die (" check to see if $file exists or not\n");}
  (undef,$comp)=split(" ",`$saclst kcmpnm f $file`);
  if ($comp eq "-12345") {die("No component name defined in the file\n");}
  if (not ($comp=~/E/ or $comp=~/1/)) {die("Please input only E/1 comp\n");}
  $ecomp = $comp; $ncomp = $ecomp; $rcomp = $ecomp; $tcomp = $ecomp;
  $ncomp=~s/E/N/; $ncomp =~s/1/2/;
  $tcomp=~s/E/T/; $tcomp =~s/1/T/;
  $rcomp=~s/E/R/; $rcomp =~s/1/R/;
  ($dir) = split(" ",`dirname $file`);$east = $file;
  if ($opt_d) {
    (undef,undef,undef,undef,undef,undef,$network,undef)=split(/\./,`basename $file`);
    (undef,$east1) = split(/\.$network\./,$file);
    $east1 = "$network.$east1";
    $north1 = $east1; $north1=~s/$ecomp/$ncomp/;
    ($north) = split(" ",`ls -1 $dir/*$north1`);
    $tang1=$east1;$radial1=$east1;
    $tang1=~s/$ecomp/$tcomp/; $radial1=~s/$ecomp/$rcomp/;
    $tang = "$dir/$tang1"; $radial = "$dir/$radial1";
  } else {
    $north = $east; $north =~s/$ecomp/$ncomp/;
    $tang = $east; $tang =~s/$ecomp/$tcomp/;
    $radial = $east; $radial =~s/$ecomp/$rcomp/;
}

#  print "   rotate $north $east \n   to $tang $radial\n";
  if (!-f $north or !-f $east ) {
    print " no such files: $north and $east, check if -d option is missing\n";
    next;}

# check some header variables
  @tmp=`$saclst b kcmpnm npts gcarc cmpinc cmpaz nzyear nzjday nzhour nzmin nzsec nzmsec f $north $east`;
  (undef,$nb,$ncomp,$npts_n,$gcarc_n,$ninc,$naz,$ny,$nj,$nh,$nm,$ns,$nms)=split(" ",$tmp[0]);
  (undef,$eb,$ecomp,$npts_e,$gcarc_e,$einc,$eaz,$ey,$ej,$eh,$em,$es,$ems)=split(" ",$tmp[1]);

  if (abs($ninc-$undef)<$eps or abs($einc-$undef)<$eps or abs($naz-$undef)<$eps or abs($eaz-$undef)<$eps ) {die("Check header cmpinc and cmpaz for $north and $east\n");}

  if ($ny != $ey or $nj != $ej or $nh != $eh or $nm != $em or $ns != $es or $nms != $ems)
    {die(" Not same reference time for $north and $east\n");}

  if (abs($gcarc_e-$undef)<$eps or abs($gcarc_n-$undef)<$eps )
    {die(" Check to see if GCARC is defined in the header\n");}

  # check if the reference time is the same or not!

  if ($opt_c) {open(SAC,"|sac2000");}
  else {open(SAC,"|sac2000 > /dev/null");}
  print SAC "echo on\n";
  print SAC "r $north $east\n ";
  if ($opt_L or $npts_n != $npts_e or abs($nb-$eb) > $eps) { #cut properly
    if(!$opt_L){
      print SAC "setbb begin ( max &1,b &2,b ) \n";
      print SAC "setbb end   ( min &1,e &2,e ) \n";}
    else{
      print SAC "setbb begin ( max &1,b &2,b ( &1,o + $opt_l ) ( &2,o + $opt_l ) ) \n";
      print SAC "setbb end   ( min &1,e &2,e ( &1,o + $opt_L ) ( &2,o + $opt_L ) ) \n"; }
    print SAC "cut %begin% %end%\n r $north $east\n cut off\n";

  }
  print SAC "rot\n";           # this is the default rotate to normal
  print SAC "ch file 1 kcmpnm $rcomp\n";
  print SAC "ch file 2 kcmpnm $tcomp\n";
  print SAC "w $radial $tang\nquit\n";
  close(SAC);
}

print "   Done !\n";
