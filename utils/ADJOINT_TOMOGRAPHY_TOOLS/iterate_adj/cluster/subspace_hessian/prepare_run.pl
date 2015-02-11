#!/usr/bin/perl -w

#-----------------------------------
# Carl Tape, 03-July-2008
# prepare_run.pl
#
# This prepares a run for computing the subspace hessian matrix
# from a set of event kernels.
#
# EXAMPLE:
#    prepare_run.pl m14 1
#
#-----------------------------------

if (@ARGV < 2) {die("Usage: prepare_run.pl xxx\n")}
($smodel,$ismooth) = @ARGV;

# base directory with event kernels
$kdir0 = "/ibrixfs1/home/carltape/ADJOINT_TOMO/iterate_adj/pangu_work/KERNELS_MODELS";
#$kdir0 = "../KERNELS_MODELS";
if(not -e $kdir0) {die("dir kdir0 $kdir0 does not exist")}

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
`rm -rf $kdir3/*`;
`ln -s $kdir2/[0-9]* $kdir3`;

# print labels to be read by subspace_hessian.f90
open(OUT,">INPUT/ismooth"); print OUT "$ismooth"; close(OUT);
open(OUT,">INPUT/smodel"); print OUT "$smodel"; close(OUT);

#=================================================================
