#!/usr/bin/perl -w

#-----------------------------------
# Carl Tape, 01-Feb-2009
# setup_model_vp_vs.pl
#
# Start with one vp and vs model, and then add a perturbation and compute new vp and vs models.
# Note the minmax wavespeed specifications within the fortran code.
#
# EXAMPLE:
#    setup_model_vp_vs.pl 168 m15 dm15 vp_m15 vs_m15 dm_kappa_kernel_smooth_p090 dm_mu_kernel_smooth_p090
#    setup_model_vp_vs.pl 168 m14 dm14 vp_m14 vs_m14 dm_kappa_kernel_smooth_p080 dm_mu_kernel_smooth_p080
#    setup_model_vp_vs.pl 168 m13 dm13 vp_m13 vs_m13 dm_kappa_kernel_smooth_p078 dm_mu_kernel_smooth_p078
#    setup_model_vp_vs.pl 168 m12 dm12 vp_m12 vs_m12 dm_kappa_kernel_smooth_p084 dm_mu_kernel_smooth_p084
#    setup_model_vp_vs.pl 168 m11 dm11 vp_m11 vs_m11 dm_kappa_kernel_smooth_p040 dm_mu_kernel_smooth_p040
#    setup_model_vp_vs.pl 168 m10 dm10 vp_m10 vs_m10 dm_kappa_kernel_smooth_p040 dm_mu_kernel_smooth_p040
#    setup_model_vp_vs.pl 168 m09 dm09 vp_m09 vs_m09 dm_kappa_kernel_smooth_p041 dm_mu_kernel_smooth_p041
#    setup_model_vp_vs.pl 168 m08 dm08_smooth vp_m08 vs_m08 dm_kappa_kernel_smooth_p048_smooth dm_mu_kernel_smooth_p048_smooth
#    setup_model_vp_vs.pl 168 m07 dm07_smooth vp_m07 vs_m07 dm_kappa_kernel_smooth_p044_smooth dm_mu_kernel_smooth_p044_smooth
#    setup_model_vp_vs.pl 168 m06 dm06_smooth vp_m06 vs_m06 dm_kappa_kernel_smooth_p048_smooth dm_mu_kernel_smooth_p048_smooth
#    setup_model_vp_vs.pl 168 m05 dm05_smooth vp_m05 vs_m05 dm_kappa_kernel_smooth_p044_smooth dm_mu_kernel_smooth_p044_smooth
#    setup_model_vp_vs.pl 168 m04 dm04_smooth vp_m04 vs_m04 dm_kappa_kernel_smooth_p040_smooth dm_mu_kernel_smooth_p040_smooth
#    setup_model_vp_vs.pl 168 m03 dm03_smooth vp_m03 vs_m03 dm_kappa_kernel_smooth_p025_smooth dm_mu_kernel_smooth_p025_smooth
#    setup_model_vp_vs.pl 168 m02 dm02_smooth vp_m02 vs_m02 dm_kappa_kernel_smooth_p031_smooth dm_mu_kernel_smooth_p031_smooth
#    setup_model_vp_vs.pl 168 m01 dm01_smooth vp_m01 vs_m01 dm_kappa_kernel_smooth_p060_smooth dm_mu_kernel_smooth_p060_smooth
#    setup_model_vp_vs.pl 168 m00 dm00_smooth vp_m00 vs_m00 dm_kappa_kernel_smooth_p068_smooth dm_mu_kernel_smooth_p068_smooth
#
#-----------------------------------

if (@ARGV < 7) {die("Usage: setup_model_vp_vs.pl xxx\n")}
($nproc,$mtag1,$mtag2,$ftag1,$ftag2,$ftag3,$ftag4) = @ARGV;

$basedir = "/ibrixfs1/home/carltape/ADJOINT_TOMO_EXTRA/KERNELS_MODELS/models";
if(not -e $basedir) {die("basedir $basedir does not exist");}

@ftags = ($ftag1,$ftag2,$ftag3,$ftag4);

$olam = "INPUT/ftags";
open(OUT,">$olam");

$count = 0;
for ($i = 1; $i <= 4; $i++) {

  if($i == 1) {$mdir = "$basedir/$mtag1";}
  if($i == 2) {$mdir = "$basedir/$mtag1";}
  if($i == 3) {$mdir = "$basedir/$mtag2/bulk";}
  if($i == 4) {$mdir = "$basedir/$mtag2/beta";}
  if(not -e $mdir) {die("mdir $mdir does not exist");}

  $ftag = $ftags[$i-1];
  $efiles_in  = "$mdir/*${ftag}*";
  $efiles_out = "INPUT_MODEL/*${ftag}*";
  $nfiles_in  = `ls -1 $efiles_in | wc | awk '{print \$1}'`; chomp($nfiles_in);
  $nfiles_out = `ls -1 $efiles_out | wc | awk '{print \$1}'`; chomp($nfiles_out);
  print "nfiles = ${nfiles_in}, nfiles_out = ${nfiles_out}\n";

  if (${nfiles_in} != $nproc ) {
    die("no $nproc files in ${efiles_in}\n");
  } else {
    print "--> ${nfiles_in} bin files exist\n";
  }

  if ( ${nfiles_out} == $nproc ) {
    print "--> ${nfiles_out} bin files exist\n";
    $count = $count + 1;
  } else {
    `ln -s ${efiles_in} INPUT_MODEL`;
  }

  # write out file tags to read in
  print OUT "$ftag\n";

}
close(OUT);

if($count == 4) {print "\n===> YOU ARE READY TO GO\n\n"}

#=================================================================
