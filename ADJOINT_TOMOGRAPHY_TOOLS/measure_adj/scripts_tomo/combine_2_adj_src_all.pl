#!/usr/bin/perl -w

#-----------------------------------
# Carl Tape, 03-Oct-2008
# combine_2_adj_src_all.pl
#
# This script calls combine_adj_src.pl to plot make composite adjoint sources.
# for a set of events.
# 
# EXAMPLE:
#    combine_2_adj_src_all.pl 2/30 6/30 cc cc m07
#
#-----------------------------------

if (@ARGV < 5) {die("Usage: combine_2_adj_src_all.pl Tmin Tmax tag1 tag2 model iplot\n")}
($Ts1,$Ts2,$tag1,$tag2,$smodel) = @ARGV;

$pwd = $ENV{PWD};

($Tmin1,$Tmax1) = split("/",$Ts1);
($Tmin2,$Tmax2) = split("/",$Ts2);
$Ttag1 = sprintf("T%3.3i_T%3.3i",$Tmin1,$Tmax1);
$Ttag2 = sprintf("T%3.3i_T%3.3i",$Tmin2,$Tmax2);

# list of event IDs
#$file_eids = "/net/sierra/raid1/carltape/results/SOURCES/eid_run";     # favorites
$file_eids = "/net/sierra/raid1/carltape/results/SOURCES/socal_11/EIDs_only_eid";
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
$cshfile = "combine_2_adj_src_all.csh";
print "\nWriting to $cshfile ...\n";
open(CSH,">$cshfile");

# loop over all events
$k = 0;
$imin = 1; $imax = $nevent;   # default
#$imin = 1; $imax = 100;
#$imin = 177; $imax = $imin;

for ($i = $imin; $i <= $imax; $i = $i+1) {

  $eid = $eids[$i-1]; chomp($eid);
  $ftag = "${eid}_${sTmin1}";

  print "\n------------------------------";
  print "\n $i, $eid";
  print CSH "echo $i, $eid\n";

  # directories
  $dir_run_eid = "${dir_run}/${eid}/${smodel}";        # note EID
  $dir_meas1 = "${dir_run_eid}/MEASURE_${Ttag1}";
  $dir_meas2 = "${dir_run_eid}/MEASURE_${Ttag2}";
  $dir_adj1 = "${dir_meas1}/ADJOINT_${Ttag1}";
  $dir_adj2 = "${dir_meas2}/ADJOINT_${Ttag2}";

  # name of stations file
  $sta_file1 = "${dir_adj1}/STATIONS_ADJOINT";
  $sta_file2 = "${dir_adj2}/STATIONS_ADJOINT";
  
  # combined output directory
  $odir = "${dir_run_eid}/ADJOINT_all";

  # check if the combined adjoint sources have already been run
  if ( -e $odir ) {
    print " --> $odir already exists\n";

  } else {

    if ( (-f $sta_file1) && (-f $sta_file2) ) {
      print " --> both STATIONS_ADJOINT files exist\n";

    #if ( -f $sta_file1 || -f $sta_file2) {
    #  print " --> at least one STATIONS_ADJOINT file exists";

      # make output directory
      print CSH "rm -rf $odir\n mkdir $odir\n";

      if ( (-f $sta_file1) && (-f $sta_file2) ) {
	# Case 1: both STATIONS_ADJOINT files exist: execute combine_adj_src.pl
	print CSH "echo running combine_adj_src.pl...\n";
	print CSH "combine_adj_src.pl $dir_adj1 $dir_adj2 $odir $tag1 $tag2\n";

      } else {

	# Case 2: one STATIONS_ADJOINT file: move files into output directory

        # NOTE: mv_files.pl DOES NOT SEEM TO WORK FOR THE LONG FILE NAMES (06-FEB-2008)
	print CSH "echo copying adjoint sources from one directory...\n";
	print CSH "echo CHECK THAT THE FILE NAME EXTENSIONS HAVE BEEN REMOVED...\n";
	if ( -f $sta_file1 ) {
	  print CSH "cp ${dir_adj1}/* $odir\n";
	  print CSH "mv_files.pl -x \"$odir/*.cc.adj\" \"$odir/*.adj\"\n";

	} elsif ( -f $sta_file2 ) {
	  print CSH "cp ${dir_adj2}/* $odir\n";
	  print CSH "mv_files.pl -x \"$odir/*.mtm.adj\" \"$odir/*.adj\"\n";
	}
      }

    } else {
      # Case 3: neither STATIONS_ADJOINT exists
      print " --> neither STATIONS_ADJOINT file exists";
    }

  }

}

#-----------------------------------------
close(CSH);
#system("csh -f $cshfile");

print "\n\n";

#================================================
