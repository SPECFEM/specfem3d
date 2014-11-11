#!/usr/bin/perl -w

#-----------------------------------
# Carl Tape, 12-March-2008
# copy_files.pl
#
# INPUT:
#    parm    -- mu_kernel, mu_kernel_smooth, kappa_kernel, kappa_kernel_smooth
#    ftag    -- event, window
#    idelete -- 0, 1
#
# EXAMPLE:
#    copy_files.pl m0 both event
#    copy_files.pl m0 both_smooth event
#
#    copy_files.pl m0 mu_kernel_smooth event
#    copy_files.pl m0 mu_kernel_smooth window
#    copy_files.pl m0 mu_kernel_smooth none
#
#-----------------------------------

#if (@ARGV < 3) {die("Usage: copy_files.pl xxx\n")}
#($smodel,$parm,$ftag,$idelete) = @ARGV;

# copy all mu files onto pangu
$dir0 = "/net/denali/scratch1/carltape/svn/cig/seismo/3D/ADJOINT_TOMO/iterate_adj_work/matlab";
$odir = "${dir0}/OUTPUT";
if (not -e $odir) {die("check if odir $odir exist or not\n")}

`cp -r $odir .`;

die("stopping here");

#--------------------------------------

#if($idelete == 1) {
#  `rm -rf $ftag`;
#  `mkdir $ftag`;
#}

$file1 = "$dir0/data_norm_${ftag}";
$file2 = "$dir0/dcov_fac_${ftag}";
$file3 = "$dir0/nwin_tot";
$file4 = "$dir0/dmisfit";
#$dir1 = "$dir0/mu_dirs/$smodel/mu_${parm}_${ftag}";    # also includes lambda

if (not -f $file1) {die("check if file1 $file1 exist or not\n")}
if (not -f $file2) {die("check if file2 $file2 exist or not\n")}
if (not -f $file3) {die("check if file3 $file3 exist or not\n")}
if (not -f $file4) {die("check if file4 $file4 exist or not\n")}
#if (not -e $dir1) {die("check if dir1 $dir1 exist or not\n")}

print "\n copying files into local directory $ftag\n";
`cp -r $file1 $ftag`;
`cp -r $file2 $ftag`;
`cp -r $file3 .`;
`cp -r $file4 .`;
`cp -r $dir1 $ftag`;

#=================================================================
