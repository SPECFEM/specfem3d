#!/usr/bin/perl -w

#-----------------------------------
# Carl Tape, 03-Oct-2008
# combine_3_adj_src_all.pl
#
# This script calls combine_adj_src.pl to plot make composite adjoint sources
# for a set of events.
# 
# EXAMPLE:
#   combine_3_adj_src_all.pl 2/30 3/30 6/30 mtm mtm mtm m15
#   combine_3_adj_src_all.pl 2/30 3/30 6/30 cc cc cc m14
#   combine_3_adj_src_all.pl 2/30 3/30 6/30 bdcc bdcc bdcc m16
#
#-----------------------------------

if (@ARGV < 7) {die("Usage: combine_3_adj_src_all.pl Tmin Tmax tag1 tag2 tag3 smodel\n")}
($Ts1,$Ts2,$Ts3,$tag1,$tag2,$tag3,$smodel) = @ARGV;

$pwd = $ENV{PWD};

($Tmin1,$Tmax1) = split("/",$Ts1);
($Tmin2,$Tmax2) = split("/",$Ts2);
($Tmin3,$Tmax3) = split("/",$Ts3);
$Ttag1 = sprintf("T%3.3i_T%3.3i",$Tmin1,$Tmax1);
$Ttag2 = sprintf("T%3.3i_T%3.3i",$Tmin2,$Tmax2);
$Ttag3 = sprintf("T%3.3i_T%3.3i",$Tmin3,$Tmax3);

# list of event IDs
$file_eids = "/net/sierra/raid1/carltape/results/EID_LISTS/kernels_run_${smodel}";
if (not -f $file_eids) {die("\n check if $file_eids exists\n")}
open(IN,$file_eids); @eids = <IN>; $nevent = @eids;

# run directory
$dir_run = "/net/sierra/raid1/carltape/socal/socal_3D/RUNS";
if (not -e $dir_run) {die("check if $dir_run exist or not\n")}

if (0==1) {
  for ($i = 1; $i <= $nevent; $i = $i+1) {
    $eid = $eids[$i-1]; chomp($eid);
    print "$i -- $eid\n";
  }
  die("testing");
}

#----------------------------------------------------

#$nevent = 5;

# write the C-shell script to file
$cshfile = "combine_3_adj_src_all.csh";
print "\nWriting to $cshfile ...\n";
open(CSH,">$cshfile");

# loop over all events
$k = 0;
$imin = 1; $imax = $nevent;   # default
#$imin = 1; $imax = 128;
#$imin = 130; $imax = $imin;

for ($i = $imin; $i <= $imax; $i = $i+1) {

  $eid = $eids[$i-1]; chomp($eid);
  print "------------------------------\n";
  print "$i, $eid\n";
  print CSH "echo $i, $eid\n";

  # directories
  $dir_run_eid = "${dir_run}/${eid}/${smodel}";        # note EID
  $dir_meas1 = "${dir_run_eid}/MEASURE_${Ttag1}";
  $dir_meas2 = "${dir_run_eid}/MEASURE_${Ttag2}";
  $dir_meas3 = "${dir_run_eid}/MEASURE_${Ttag3}";
  $dir_adj1 = "${dir_meas1}/ADJOINT_${Ttag1}";
  $dir_adj2 = "${dir_meas2}/ADJOINT_${Ttag2}";
  $dir_adj3 = "${dir_meas3}/ADJOINT_${Ttag3}";

  # name of stations file
  $sta_file1 = "${dir_adj1}/STATIONS_ADJOINT";
  $sta_file2 = "${dir_adj2}/STATIONS_ADJOINT";
  $sta_file3 = "${dir_adj3}/STATIONS_ADJOINT";
  
  # combined output directory
  $odir = "${dir_run_eid}/ADJOINT_all";

  # check if the combined adjoint sources have already been run
  if ( -e $odir ) {
    print " --> $odir already exists\n";

  } else {
    if ( (-f $sta_file1 && -f $sta_file2) && (-f $sta_file3) ) {
      print " --> all three STATIONS_ADJOINT files exist\n";
      print CSH "echo running combine_adj_src.pl...\n";
      print CSH "combine_3_adj_src.pl $dir_adj1 $dir_adj2 $dir_adj3 $odir $tag1 $tag2 $tag3\n";

    } else {
      print " --> all three STATIONS_ADJOINT files do NOT exist\n";
      if (not -f $sta_file1) {print "sta_file1 --> $sta_file1\n";}
      if (not -f $sta_file2) {print "sta_file2 --> $sta_file2\n";}
      if (not -f $sta_file3) {print "sta_file3 --> $sta_file3\n";}
      #die("STOPPING: all three STATIONS_ADJOINT files do NOT exist\n");
    }
  }
}

#-----------------------------------------
close(CSH);
#system("csh -f $cshfile");

print "\n--> done with $cshfile\n\n";

#================================================
