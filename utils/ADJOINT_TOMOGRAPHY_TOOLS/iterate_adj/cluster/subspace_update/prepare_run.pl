#!/usr/bin/perl -w

#-----------------------------------
# Carl Tape, 11-Oct-2008
# prepare_run.pl
#
# Prepare run to compute multiple possible model updates from
# a set of mu vectors computed in Matlab.
#
# After you view the possible models, you pick one.
# Make sure you have both the mu and the kappa update.
#
# EXAMPLE:
#    prepare_run.pl 60/4/80 m14 mu_kernel both 1 window 1
#    prepare_run.pl 20/4/60 m10 mu_kernel both 1 window 1
#    prepare_run.pl 44/1/44 m10 kappa_kernel both 1 window 1
#
#    prepare_run.pl 1/5/100 m01 mu_kernel mu_kernel 1 none 1
#    prepare_run.pl 1/5/134 m01 mu_kernel mu_kernel 1 none 1
#    prepare_run.pl 2/5/80 m01 kappa_kernel kappa_kernel 1 none 1
#
#-----------------------------------

if (@ARGV < 7) {die("Usage: prepare_run.pl pvals smodel parm0 plam0 ismooth ftag idelete\n")}
($pvals,$smodel,$parm0,$plab0,$ismooth,$ftag,$idelete) = @ARGV;

# base directory with event kernels
$kdir0 = "/ibrixfs1/home/carltape/ADJOINT_TOMO/iterate_adj/pangu_work/KERNELS_MODELS";
if(not -e $kdir0) {die("dir kdir0 $kdir0 does not exist")}

($pstart,$pinc,$pend) = split("/",$pvals);
if($ismooth == 1) {
   $plab = "${plab0}_smooth";
   $parm = "${parm0}_smooth";
} else {
   $plab = "$plab0";
   $parm = "$parm0";
}

if($ismooth == 1) {
   $kdir1 = "$kdir0/event_kernels_smooth";
} else {
   $kdir1 = "$kdir0/event_kernels";
}
$kdir2 = "$kdir1/kernel_${smodel}";
if(not -e $kdir2) {die("dir kdir2 $kdir2 does not exist")}

# link the desired kernels into INPUT_KERNELS
$kdir3 = "INPUT_KERNELS";
print "linking from $kdir2 to $kdir3\n";
`rm $kdir3/*`;
`ln -s $kdir2/[0-9]* $kdir3`;

# TO DO: check that all the files are there

# check directories
$idir = "INPUT";
$mudir0 = "$idir/OUTPUT_SUBSPACE/$smodel/$ftag";
$mudir1 = "$mudir0/$plab";
$mudir = "$mudir1/mu_all";
if(not -e $mudir) {die("dir mudir $mudir does not exist")}

# check the number of pmax files
($nump,undef,undef) = split(" ",`ls -1 $mudir/mu* | wc`);
print "$nump files in $mudir\n";
if($pend > $nump) {die("$pend > $nump")}
if($nump == 0) {die("no files in $mudir")}

# copy files for the run
$rdir = "$idir/mu_run";
`mkdir -p $rdir`;
`rm $rdir/*`;

#-----------------------

# write the list of pmax values to try
$pmax_file = "$rdir/pmax";
`seq $pstart $pinc $pend > $pmax_file`;

# read in the list of pmax values you want to use
open(IN,"${pmax_file}"); @lines = <IN>; close(IN);
$npmax = @lines;
if($npmax <= 0) {die("npmax $npmax less than or equal to zero\n");}

# read in the list of pmax values you want to use
$lambda_file = "$mudir/lambda_${ftag}";
open(IN,"${lambda_file}"); @lamlines = <IN>; close(IN);
$nlam = @lamlines;
if($nlam != $nump) {die("nlam ($nlam) does not equal nump ($nump)");}

$olam = "$rdir/lams";
open(LAM,">$olam");

for ($i = 1; $i <= $npmax; $i++) {

  $pmax = $lines[$i-1]; chomp($pmax);
  $stp = sprintf("p%3.3i",$pmax);

  # copy mu files into run directory
  @files = glob("$mudir/*${stp}*");
  $nfile = @files;
  if($nfile != 1) {die("there is not one file in files @files\n")}
  $ifile = $files[0];
  print "$i -- $pmax -- $ifile\n";
  $ofile = "$rdir/mu_${stp}";
  `cp $ifile $ofile`;

  # create a lambda file for reading in
  ($lamval,$invlamval,$index) = split(" ",$lamlines[$pmax-1]); chomp($lamval);
  print LAM "$lamval\n";

}

close(LAM);

($nmu,undef,undef) = split(" ",`ls -1 $rdir/mu* | wc`);
if($nmu != $npmax) {die("nmu not equal npmax, $nmu, $npmax\n");}

#-----------------------

$din = "$mudir0/${smodel}_${ftag}_data_norm";
$dout = "$rdir/data_norm";
if(not -f $din) {die("file din $din does not exist")}
`cp $din $dout`;

$din = "$mudir0/${smodel}_${ftag}_dcov_fac";
$dout = "$rdir/dcov_fac";
if(not -f $din) {die("file din $din does not exist")}
`cp $din $dout`;

# print the dcov label to be read by subspace_update.f90
open(OUT,">$rdir/dcov_tag"); print OUT "$ftag"; close(OUT);

# print the kernel label to be read by subspace_update.f90
open(OUT,">$rdir/ker_tag"); print OUT "$parm"; close(OUT);

#=================================================================
