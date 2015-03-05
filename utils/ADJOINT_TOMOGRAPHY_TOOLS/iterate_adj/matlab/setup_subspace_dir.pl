#!/usr/bin/perl -w

#-----------------------------------
# Carl Tape, 03-Sept-2008
# setup_subspace_dir.pl
#
# This script creates a set of directories for computing misfit
# functions and mu-vectors for the subspace method for different
# choices of the data covariance matrix.
#
# The user should keep a copy of the output files in OUTPUT_SUBSPACE
# somewhere outside of the working SVN directory, since these output
# files will not be checked into SVN.
#
# This should be run prior to running subspace_specfem.m.
#
# EXAMPLE:
#    setup_subspace_dir.pl 0 10
#    setup_subspace_dir.pl 11 11
#
#-----------------------------------

if (@ARGV < 2) {die("Usage: setup_subspace_dir.pl smodel_min smodel_max\n")}
($minm,$maxm) = @ARGV;

# possible labels for choice of data covariance
@dcov_tags = ("event","window","none");

# possible labels for choice of kernels
@kers = ("mu_kernel","kappa_kernel","mu_kernel_smooth","kappa_kernel_smooth","both","both_smooth");

$ndcov = @dcov_tags;
$nker = @kers;

$dir0 = "OUTPUT_SUBSPACE";
`mkdir -p $dir0`;
#if(not -e $dir0) {die("dir dir0 $dir0 does not exist\n");}

#$dir1 = "$dir0/$smodel";
#`mkdir -p $dir1`;

for ($h = $minm; $h <= $maxm; $h = $h+1) {

  $smodel = sprintf("m%2.2i",$h);
  $dir1 = "${dir0}/${smodel}";
  `mkdir -p $dir1`;

  for ($j = 1; $j <= $ndcov; $j = $j+1) {

  $dcov = $dcov_tags[$j-1];
  $dir2 = "${dir1}/${dcov}";
  `mkdir -p $dir2`;

    for ($i = 1; $i <= $nker; $i = $i+1) {

      $ker  = $kers[$i-1];
      $dir  = "$dir2/${ker}";
      print "$dir \n";

      if (-e $dir) {
  #print "--> dir exists -- now deleting and remaking\n"; `rm -rf $dir`;
        print "--> dir exists\n";
      } else {
  print "--> dir does not exist -- now making\n";
      }
      `mkdir -p $dir`;
      `mkdir -p $dir/mu_all`;
    }
  }
}

#================================================
