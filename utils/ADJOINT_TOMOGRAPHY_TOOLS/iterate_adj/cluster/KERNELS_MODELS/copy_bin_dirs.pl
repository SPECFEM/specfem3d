#!/usr/bin/perl -w

#-----------------------------------
# Carl Tape, 04-Oct-2008
# copy_bin_dirs.pl
#
# This script copies directories with binary files of a specified parameter
# from one place to a local directory on pangu.
#
# INPUT
#    smodel    model within the CG iteration (m00, m00t, m01, etc)
#    hdur_lab  label for the measurements (T02, T06, all)
#    nproc     number of processors
#    ftag      type of files you want to copy
#
# EXAMPLE:
#    copy_bin_dirs.pl m10 all 168 mu_kernel ; copy_bin_dirs.pl m10 all 168 kappa_kernel
#
#-----------------------------------

if (@ARGV < 2) {die("Usage: copy_bin_dirs.pl xxx\n")}
($smodel,$hdur_lab,$nproc,$ftag) = @ARGV;

# directory containing window files, measurement files, adjoint sources, etc
$run_dir = "/net/sierra/raid1/carltape/socal/socal_3D/RUNS";

# read in the list of kernels that you want to sum
#$kernel_list = "kernels_run";
$kernel_list = "kernel_lists/kernels_use_${smodel}";
if (not -e ${kernel_list}) {die("check if kernel_list $kernel_list exist or not\n")}
open(IN,"${kernel_list}"); @lines = <IN>; close(IN);
$neid = @lines;

# copy link the kernel list to the generic file name
`ln -s ${kernel_list} kernels_run`;

# list all events
if (0==1) {
  for ($i = 1; $i <= $neid; $i++) {
    $eid = $lines[$i-1]; chomp($eid);
    print "$i : $eid\n";
  }
  die("testing");
}

$imin = 1; $imax = $neid; # default
#$imin = 1; $imax = 30;
#$imin = 89; $imax = $imin;

#===================================

for ($i = $imin; $i <= $imax; $i++) {

  $eid = $lines[$i-1]; chomp($eid);
  print "$i : $eid\n";

  # directory containing the binary kernel files
  $run_dir_bin = "${run_dir}/${eid}/${smodel}/MESH_${hdur_lab}/collect";
  $bin_files = "${run_dir_bin}/*$ftag*.bin";

  # number of files
  $nfile_in = `ls -1 ${bin_files} | wc | awk '{print \$1}'`; chomp($nfile_in);

  if ($nfile_in != $nproc) {
    print "\n -- nfile_in = $nfile_in -- nproc = $nproc --\n";
    print "mismatch of expected vs actual bin files\n";
    #die("mismatch of expected vs actual bin files\n");

  } else {

    # number of files in output directory
    $odir = "event_kernels/kernel_${smodel}/$eid";
    $nfile_out = `ls -1 $odir/*$ftag*.bin | wc | awk '{print \$1}'`; chomp($nfile_out);

    #-------------

    if ($nfile_out == $nproc) {
      print "--> already have the kernel files\n";

    } else {

      `mkdir -p $odir`;
      `cp $bin_files $odir`;
      `sleep 10s`;

      $nfile_out = `ls -1 $odir/*$ftag*.bin | wc | awk '{print \$1}'`; chomp($nfile_out);

      if ($nfile_out != $nproc) {
  print "\n -- nfile_out = $nfile_out -- nproc = $nproc --\n";
  die("mismatch of expected vs actual bin files\n");
      }
    }

  }

}  # for loop

#=================================================================
