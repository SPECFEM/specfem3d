#!/usr/bin/perl -w

#
# EXAMPLE: setup_forward_dir.pl m16 1 10
#

if (@ARGV < 3) {die("Usage: setup_forward_dir.pl smodel imin imax\n")}
($smodel,$imin,$imax) = @ARGV;

#-----------------------------------
# USER INPUT

$main_dir = "/ibrixfs1/home/carltape/BASIN_FORWARD/Main_Dir";
$scripts_dir = "/ibrixfs1/home/carltape/BASIN_FORWARD/Scripts";

# directory containing CMTSOLUTION files for all events
$cmt_dir = "/ibrixfs1/home/carltape/CMT/CMT_files";

# directory containing finished synthetics
$done_dir = "/net/sierra/raid1/carltape/socal/socal_3D/SYN/model_${smodel}";

# list of ALL event IDs and simulation durations
$event_list = "/ibrixfs1/home/carltape/CMT/EIDs_simulation_duration";

# list of EIDs of events for which you want synthetics
$eid_list = "/net/sierra/raid1/carltape/results/EID_LISTS/syn_run_${smodel}";

#-----------------------------------

# read in the list of ALL event IDs and simulation durations
$event_list = "/ibrixfs1/home/carltape/CMT/EIDs_simulation_duration";
open(IN,"${event_list}"); @lines = <IN>; close(IN);
$neid = @lines;

# read in list of synthetics for which you want synthetics
if (not -f ${eid_list}) {print "\n check if eid_list ${eid_list} exists\n";}
open(IN,"${eid_list}"); @keids = <IN>; close(IN);
$nrun = @keids;

# half-duration for the simulation
$hdur = 0;

if ($imax > $neid) {$imax = $neid;}

# loop over all events
for ($i = $imin; $i <= $imax; $i++) {

  ($index,$eid,$sdur_round,$npt0,$hdur0,$mag0) = split(" ",$lines[$i-1]);
  print "\n $i : $eid";

  # find the matching EID in the subset list
  $imatch = 0;
  for ($j = 1; $j <= $nrun; $j++) {
    $keid = $keids[$j-1];
    if ($eid == $keid) {
      $imatch = $imatch+1;
    }
  }

  if ($imatch == 0) {
    print "\n --> no match in EID list";

  } else {

    $dir1 = $eid;
    $dir2 = "DONE/${eid}";
    $dir3 = "DONE_SINGLE/${eid}";

    if ( (-e "${done_dir}/$eid" || -e $dir3) || (-e $dir1 || -e $dir2 ) ) {
      print "\n --> event has already been setup or done";

    } else {

      $cmtfile = "${cmt_dir}/CMTSOLUTION_$eid";
      if (not -f $cmtfile) {
  print "\n --> cmtfile does not exist";

      } else {

  # copy in a new CMTSOLUTION
  $cmtfile1 = "$eid/DATA/CMTSOLUTION";
  `${scripts_dir}/copy_basin_sem_dir.bash ${main_dir} $eid`;
  `cp $cmtfile $cmtfile1`;

  # change the half-duration
  $hdur_st = sprintf("%8.4f",$hdur);
  $ctemp   = "${cmtfile}_temp";
  `sed '/^half duration:/s/^.*\$/half duration: ${hdur_st}/' $cmtfile1 > $ctemp`;
  `mv $ctemp $cmtfile1`;

  # change the simulation time in the Par_file
  $sdur = sprintf("%.1f",$sdur_round);
  $parfile = "$eid/DATA/Par_file";
  $ptemp   = "$eid/DATA/Par_file_temp";
  `sed '/^RECORD_LENGTH_IN_SECONDS        =/s/^.*\$/RECORD_LENGTH_IN_SECONDS        = $sdur/' $parfile > $ptemp`;
  `mv $ptemp $parfile`;

  # change the label for the run
  `sed "/socal_earthquakes/s/socal_earthquakes/${eid}_earthquake/" $eid/go_mesher_solver_lsf_basin.bash > $eid/run.tmp`;
  `mv -f $eid/run.tmp $eid/go_mesher_solver_lsf_basin.bash`;
      }
    }
  }
}


#--------------------
#Main_dir=/ibrixfs1/home/carltape/BASIN_FORWARD/Main_Dir
#Scripts_dir=/ibrixfs1/home/carltape/BASIN_FORWARD/Scripts
#
#for eid in `cat /ibrixfs1/home/carltape/CMT/test_events`; do
#
# mkdir -p $eid
#
# echo "==== $eid ====="
#
# ${Scripts_dir}/copy_basin_sem_dir.bash ${Main_dir} $eid
# cp /ibrixfs1/home/carltape/CMT/CMT_files/CMTSOLUTION_$eid $eid/DATA/CMTSOLUTION
#
#done
#--------------------

#=================================================================
