#!/usr/bin/perl -w

#-----------------------------------
# Carl Tape, 30-Jan-2009
# setup_model_pert.pl
#
# Start with one vp and vs model, and then add a perturbation to new vp and vs models.
#
# Removed linked files in INPUT_MODEL prior to running.
#
# EXAMPLE:
#    setup_model_pert.pl 168 m16 m00 vs_m16 vs_m00 vs_m16_m00
#    setup_model_pert.pl 168 m16 m15 vs_m16 vs_m15 vs_m16_m15
#    setup_model_pert.pl 168 m16 m00 vp_m16 vp_m00 vp_m16_m00
#    setup_model_pert.pl 168 m16 m00 vb_m16 vb_m00 vb_m16_m00
#    setup_model_pert.pl 168 m16 m00 poisson_m16 poisson_m00 poisson_m16_m00
#
#    setup_model_pert.pl 168 m09 m00 vp_m09 vp_m00 vp_m09_m00
#    setup_model_pert.pl 168 m09 m08 vp_m09 vp_m08 vp_m09_m08
#-----------------------------------

if (@ARGV < 6) {die("Usage: setup_model_pert.pl nproc mtag1 mtag2 ftag1 ftag2 ftag3\n")}
($nproc,$mtag1,$mtag2,$ftag1,$ftag2,$ftag3) = @ARGV;

$basedir = "/ibrixfs1/home/carltape/ADJOINT_TOMO_EXTRA/KERNELS_MODELS/models";
if(not -e $basedir) {die("basedir $basedir does not exist");}

# write out tags to be read in by model_pert.f90
@ftags = ($ftag1,$ftag2,$ftag3);
$olam = "INPUT/ftags";
open(OUT,">$olam");
for ($i = 1; $i <= 3; $i++) {
  $ftag = $ftags[$i-1];
  print OUT "$ftag\n";
}
close(OUT);

# link input files into INPUT_MODEL
$count = 0;
for ($i = 1; $i <= 2; $i++) {

  if($i == 1) {$mdir = "$basedir/$mtag1";}
  if($i == 2) {$mdir = "$basedir/$mtag2";}
  if(not -e $mdir) {die("mdir $mdir does not exist");}

  $ftag = $ftags[$i-1];
  $efiles_in  = "$mdir/*${ftag}.bin";
  $efiles_out = "INPUT_MODEL/*${ftag}.bin";
  $nfiles_in  = `ls -1 $efiles_in | wc | awk '{print \$1}'`; chomp($nfiles_in);
  $nfiles_out = `ls -1 $efiles_out | wc | awk '{print \$1}'`; chomp($nfiles_out);

  print "ftag $ftag\n";
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
}

if($count == 2) {
   print "\n===> YOU ARE READY TO GO\n\n";
} else {
   print "\n===> you are NOT ready to go (try executing again)\n\n";
}

#=================================================================
