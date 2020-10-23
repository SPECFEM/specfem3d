#!/usr/bin/perl -w

#-----------------------------------
# Carl Tape, 03-Oct-2008
# setup_kernel_dir.pl
#
# This script ...
#
# INPUT
#    smodel    model within the CG iteration (m0, m0t, m1, etc)
#    hdur_lab  label for the measurements (T02, T06, all)
#    nstat     minimum number of stations to use
#
# EXAMPLE:
#    setup_kernel_dir.pl m00 all 10
#    setup_kernel_dir.pl m12 all 10
#
#-----------------------------------

if (@ARGV < 3) {die("Usage: setup_kernel_dir.pl model dir_suffix nstat\n")}
($smodel,$hdur_lab,$nstat) = @ARGV;

$main_dir = "/ibrixfs1/home/carltape/BASIN_KERNEL/Main_Dir";
$scripts_dir = "/ibrixfs1/home/carltape/BASIN_KERNEL/Scripts";

# directory containing SEM-inverted CMTSOLUTION files
$cmt_dir = "/ibrixfs1/home/carltape/CMT/CMT_files";

# directory containing window files, measurement files, adjoint sources, etc
$run_dir = "/net/sierra/raid1/carltape/socal/socal_3D/RUNS";

# read in the list of ALL event IDs and simulation durations
$event_list = "/home/carltape/CMT/EIDs_simulation_duration";
open(IN,"${event_list}"); @lines = <IN>; close(IN);
$neid = @lines;

# read in list of kernels that will be used for this model update
$ker_list = "/net/sierra/raid1/carltape/results/EID_LISTS/kernels_run_${smodel}";
if (not -f ${ker_list}) {die("check if ker_list ${ker_list} exists\n");}
open(IN,"${ker_list}"); @keids = <IN>; close(IN);
$nker = @keids;

if (0==1) {
  print "\n ALL EVENTS \n";
  for ($i = 1; $i <= $neid; $i++) {
    ($index,$eid,$sdur_round,$nstep,$hdur,$dur,$mw) = split(" ",$lines[$i-1]);
    print "\n $i : $eid";

  }
  die("testing");
}

$imin = 1; $imax = $neid;    # default
#$imin = 1; $imax = 40;
#$imin = 19; $imax = $imin;

#===================================

for ($i = $imin; $i <= $imax; $i++) {

  ($index,$eid,$sdur_round,$nstep,$hdur,$dur,$mw) = split(" ",$lines[$i-1]);
  print "$i : $eid\n";

  # setup only the EIDs that are in the subset list
  $imatch = 0;
  for ($j = 1; $j <= $nker; $j++) {
    $keid = $keids[$j-1];
    if ($eid == $keid) {
      $imatch = $imatch+1;
    }
  }

  if ($imatch == 0) {
    print "--> no match in kernel list\n";

  } else {

    $run_dir_eid = "${run_dir}/${eid}/${smodel}";
    $cmtfile = "${cmt_dir}/CMTSOLUTION_$eid";

    if (not -f $cmtfile) {
      print "--> cmtfile does not exist\n";
    } else {

      if ( -e $eid || (-e "done_dirs/${eid}" || -e "PROCESS/${eid}")  ) {
  print "--> kernel simulation has been setup or has been done\n";

      } else {

  # output directories on denali
  $output_dir = "${run_dir_eid}/OUTPUT_${hdur_lab}";
  $adj_dir = "${run_dir_eid}/ADJOINT_${hdur_lab}";

  if (-e $output_dir) {
    print "--> kernel simulation has been run\n";

  } else {

    if (not -e "${adj_dir}") {
      print "--> measurement code has not been run\n";

    } else {
      $nadj = `ls -1 ${adj_dir}/*adj | wc | awk '{print \$1}'`; chomp($nadj);
      $nstat_adj = $nadj/3;
      print "--> number of adjoint sources/stations : $nadj/$nstat_adj\n";

      if ($nstat_adj < $nstat) {
        print "--> fewer than $nstat stations\n";
      } else {
        print "--> RUN THIS EVENT\n";

        #----------------------------

        # replicate the main directory
        `${scripts_dir}/copy_basin_sem_dir.bash ${main_dir} $eid`;

              # modify the EID label for the post-processing script
              # eid=9818433
              $processfile = "$eid/go_process.bash";
              $ptemp = "$eid/go_process_temp.bash";
        `sed '/^eid=/s/^.*\$/eid=$eid/' $processfile > $ptemp`;
        `mv $ptemp $processfile`;
              `chmod +x $processfile`;

        # copy in a new CMTSOLUTION
        $cmtfile0 = "$cmtfile";
        $cmtfile1 = "$eid/DATA/CMTSOLUTION";
        `cp $cmtfile0 $cmtfile1`;

        # change the half-duration
        #$hdur_st = sprintf("%8.4f",$hdur);
        #$ctemp   = "${cmtfile1}_temp";
        #`sed '/^half duration:/s/^.*\$/half duration: ${hdur_st}/' $cmtfile1 > $ctemp`;
        #`mv $ctemp $cmtfile1`;

        # change the simulation TYPE in the Par_file
        `change_simulation_type2.pl -F`;

        # change the simulation TIME in the Par_file
        #$sdur_round = 5;   # testing only
        $sdur = sprintf("%.1f",$sdur_round);
        $parfile = "$eid/DATA/Par_file";
        $ptemp   = "$eid/DATA/Par_file_temp";
        `sed '/^RECORD_LENGTH_IN_SECONDS        =/s/^.*\$/RECORD_LENGTH_IN_SECONDS        = $sdur/' $parfile > $ptemp`;
        `mv $ptemp $parfile`;

        # make sure attenuation is OFF
        `sed '/^ATTENUATION                     =/s/^.*\$/ATTENUATION                     = .false./' $parfile > $ptemp`;
        `mv $ptemp $parfile`;

        # copy the STATIONS file
        $stafile = "${adj_dir}/STATIONS_ADJOINT";
        `cp $stafile $eid/DATA/STATIONS`;
        `cp $stafile $eid/DATA/STATIONS_ADJOINT`;

        # number of stations
        $ns = `sed -n 1p $eid/DATA/STATIONS`;
        $nadj1 = 3*$ns;

        # copy adjoint sources
        `rm $eid/SEM/*adj`;
        print "\n copying the adjoint sources onto pangu...\n";
        `cp ${adj_dir}/*adj $eid/SEM`;

        # check that all the adjoint sources are there
        $nadj2 = `ls -1 $eid/SEM/*adj | wc | awk '{print \$1}'`; chomp($nadj2);
        print "\n NUMBER OF RECORDS -- $nadj1 $nadj2 --";
        if ($nadj1 != $nadj2) {
    die("Check adjoint sources: $nadj1 does not equal $nadj2");
        }

        # make directories (and remove old ones)
        $sem_dir = "$eid/SEM/SEM_${hdur_lab}";
        $output_dir = "$eid/SEM/OUTPUT_${hdur_lab}";
        $mesh_dir = "$eid/SEM/MESH_${hdur_lab}";
        `rm -r $eid/SEM/SEM_*    ; mkdir ${sem_dir}`;
        `rm -r $eid/SEM/OUTPUT_* ; mkdir ${output_dir}`;
        `rm -r $eid/SEM/MESH_*   ; mkdir ${mesh_dir}`;

        # change the label for the run
        `sed "/socal_kernel/s/socal_kernel/${eid}_kernel/" $eid/go_kernel.bash > $eid/run.tmp`;
        `mv -f $eid/run.tmp $eid/go_kernel.bash`;

        #----------------------------

      }
    }
  }
      }
    }
  }
}       # for loop

#=================================================================
