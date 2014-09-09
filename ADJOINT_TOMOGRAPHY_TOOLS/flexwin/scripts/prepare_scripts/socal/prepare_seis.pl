#!/usr/bin/perl -w

#===============================
#  prepare_seis.pl
#  Carl Tape
#  20-Oct-2007
#
#  This script synchronizes the data and synthetics and cuts them
#  to equivalent lengths.  If the data begin before the synthetics,
#  then the synthetics are filled with zeros prior to the origin time.
#  The pre-event interval is important to quantify the signal-to-noise ratio.
#  The end time of the records is determined by the end time of the data
#  or synthetics, whichever comes first.
#
#  This is adapted from prepare_seis, a shell script,
#  which is listed below (in comments).
#
#  Check that things are lined up:
#     saclst b e delta npts f SYN_TEST/* DATA_TEST/*
#
#  EXAMPLE:
#     prepare_seis.pl
#
#===============================

print "\n Running prepare_seis.pl";

# testing -- copy the UN-prepared files here
if (0==1) {
  `rm DATA_TEST/* SYN_TEST/*`;
  @stations = ("PAS");
  #@stations = ("CLC","PDR","SRN","PAS");
  foreach $sta (@stations) { 
    print "\n $sta";
    `cp DATA/*${sta}* DATA_TEST`;
    `cp SYN/*${sta}* SYN_TEST`;
  }
  #die("\ncopying seismograms to test directories...\n");
}

# USER PARAMETERS
$dt = 0.05;       # interpolation value for data and synthetics
$tfac = 1.25;     # factor to extend lengths of records (should be > 1.0)
$itest = 0;       # test directories or not
$syn_suffix = "semd.sac.d";     # suffix for synthetic files

if($itest == 1) {
   $dirdat = "DATA_TEST"; $dirsyn = "SYN_TEST";
} else {
   $dirdat = "DATA"; $dirsyn = "SYN";
} 

# grab all data files for the windowing code
@files = glob("${dirdat}/*");
$nfile = @files;
print "\n $nfile data files to line up with synthetics\n";

#------------------------

`echo echo on > sac.mac`;
`echo readerr badfile fatal >> sac.mac`;

foreach $file (@files) { 
  # read the sac headers
  (undef,$net,$sta,$chan) = split(" ",`saclst knetwk kstnm kcmpnm f $file`);
  $comp = `echo $chan | awk '{print substr(\$1,3,1)}'`;
  chomp($comp);
  #print "\n $net $sta $chan $comp\n";

  $synt = "${dirsyn}/${sta}.${net}.BH${comp}.${syn_suffix}";
  #print "\n $file $synt ";

  if (-f $synt) { 
    print "\n $file $synt";

    # process data (interpolate)
    `echo r $file >> sac.mac`;
    `echo rtrend >> sac.mac`;
    `echo taper >> sac.mac`;
    `echo interpolate delta $dt >> sac.mac`;
    `echo w over >> sac.mac`;

    # process synthetics (interpolate)
    # NOTE: The synthetics from SPECFEM may have some static trend that
    #       one might want to remove.  However, I found that this should
    #       NOT be done, because it dramatically affects the amplitudes
    #       of the close stations (where the pulses are near the ends).
    `echo r $synt >> sac.mac`;
    #`echo rtrend >> sac.mac`;   # NO!
    #`echo taper >> sac.mac`;    # NO!
    `echo interpolate delta $dt >> sac.mac`;
    `echo w over >> sac.mac`;

    # get info on data and synthetics
    (undef,$bd,$ed,$deltad,$nptd) = split(" ",`saclst b e delta npts f $file`);
    (undef,$bs,$es,$deltas,$npts) = split(" ",`saclst b e delta npts f $synt`);
    $tlend = $ed - $bd;
    $tlens = $es - $bs;

    # cut the records
    # b = earliest start time
    # e = earliest end time, multiplied by some factor
    if($bd < $bs) {$b0 = $bd} else {$b0 = $bs}
    if($ed < $es) {$e0 = $ed} else {$e0 = $es}

    if (1==1) {
      $b = $b0;
      $tlen0 = $e0 - $b;
      $tlen = $tfac * $tlen0;	# extend record length
      $e = $b0 + $tlen;

    } else {
      $e0 = 130;   # cut the records to avoid a later, larger event (14179288)
      $b = $b0;
      $tlen0 = $e0 - $b;
      $tlen = $tlen0;
      $e = $b0 + $tlen;
    }

    $npt = int( ($e-$b)/$dt );

    if (0==1) {
      print "\n Data : $bd $ed $deltad $nptd -- $tlend";
      print "\n Syn  : $bs $es $deltas $npts -- $tlens";
      print "\n b0, e0, tlen0 : $b0, $e0, $tlen0 ";
      print "\n   b : $bd, $bs, $b ";
      print "\n   e : $ed, $es, $e ";
      print "\n npt : $nptd, $npts, $npt ";
      print "\n $tlen = $tfac * ($e0 - $b)";
      print "\n $e = $b0 + $tlen \n";
    }

    # cut records and fill zeros
    `echo r $file $synt >> sac.mac`;
    `echo cuterr fillz >> sac.mac`;
    `echo "cut $b n $npt" >> sac.mac`;
    `echo r $file $synt >> sac.mac`;      # cut both
    `echo cut off >> sac.mac`;
    `echo w over >> sac.mac`;
  }

} 
`echo quit >> sac.mac`;
`sac sac.mac`;              # KEY: EXECUTE the SAC commands
#`rm sac.mac`;

print "\n Done with prepare_seis.pl\n\n";

#----------------------------

##!/bin/sh

#echo echo on > sac.mac
#echo readerr badfile fatal >> sac.mac

#for file in DATA_TEST/* ; do

#  net=`saclst knetwk f $file | awk '{print $2}'`
#  sta=`saclst kstnm f $file | awk '{print $2}'`
#  chan=`saclst kcmpnm f $file | awk '{print $2}'`
#  #net=`echo $file | awk -F"." '{print $2}'`
#  #sta=`echo $file | awk -F"." '{print $3}'`
#  #chan=`echo $file | awk -F"." '{print $4}'`
#  comp=`echo $chan | awk '{print substr($1,3,1)}'`
#  #echo $net $sta $chan $comp

#  synt=SYN_TEST/${sta}.${net}.BH${comp}.semd.sac.d
#  #echo $synt

#  if [ -e $synt ] ; then
#    # process data (interpolate)
#    echo r $file >> sac.mac
#    echo rtrend >> sac.mac
#    echo taper >> sac.mac
#    echo interpolate delta 0.05 >> sac.mac
#    echo w over >> sac.mac

#    # process synthetics (interpolate)
#    echo r $synt >> sac.mac
#    echo "setbb etime &1,e" >> sac.mac        # end time syn
#    echo interpolate delta 0.05 >> sac.mac
#    echo w over >> sac.mac

#    # line up b and npts
#    # Check: saclst npts b e delta f files
#    echo r $file $synt >> sac.mac             # read data, syn
#    echo "setbb btime1 &1,b" >> sac.mac        # begin time data
#    echo "setbb btime2 &2,b" >> sac.mac        # begin time syn
#    echo "setbb npoints &2,npts" >> sac.mac   # number of points syn
#    echo cuterr fillz >> sac.mac
#    echo "cut %btime n %npoints" >> sac.mac   # cut the records
#    echo r $synt >> sac.mac
#    echo cut off >> sac.mac
#    echo w over >> sac.mac
#  fi
#done
#echo quit >> sac.mac

#sac sac.mac
#rm sac.mac

#----------------------------
