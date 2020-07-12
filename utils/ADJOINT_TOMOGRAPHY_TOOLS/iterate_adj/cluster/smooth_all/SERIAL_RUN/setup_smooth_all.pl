#!/usr/bin/perl -w

#-----------------------------------
# Carl Tape, 22-Feb-2009
# setup_smooth_all.pl
#
# This script copies directories with binary files of a specified parameter
# from one place to a local directory on pangu.
#
#
# EXAMPLE:
#    setup_smooth_all.pl m16 168 mu_kernel ; setup_smooth_all.pl m16 168 kappa_kernel
#
#-----------------------------------

if (@ARGV < 3) {die("Usage: setup_smooth_all.pl smodel nproc ftag\n")}
($smodel,$nproc,$ftag) = @ARGV;

$basedir = "/ibrixfs1/home/carltape/ADJOINT_TOMO/iterate_adj/pangu_work";
$edir = "$basedir/KERNELS_MODELS/event_kernels/kernel_${smodel}";
$edir_smooth = "$basedir/KERNELS_MODELS/event_kernels_smooth/kernel_${smodel}";
$dir_done = "/net/sierra/raid1/carltape/socal/socal_3D/RUNS/SMOOTH_EVENT_KERNELS/$smodel";

$maindir = "$basedir/smooth";

if(not -e $edir) {die("edir $edir does not exist");}
if(not -e $edir_smooth) {`mkdir $edir_smooth`;}
if(not -e $maindir) {die("maindir $maindir does not exist");}
if(not -e ${dir_done}) {die("dir_done $dir_done does not exist");}

# read in the list of kernels that you want to sum
$kernel_list = "$basedir/KERNELS_MODELS/kernels_run";
if (not -e ${kernel_list}) {die("check if kernel_list $kernel_list exist or not\n")}
open(IN,"${kernel_list}"); @eids = <IN>; close(IN);
$neid = @eids;

print "number of event kernels : $neid\n";

$imin = 1; $imax = $neid; # default
#$imin = 1; $imax = 20;
#$imin = 10; $imax = $imin;

for ($i = $imin; $i <= $imax; $i++) {

   $eid = $eids[$i-1]; chomp($eid);
   print "======== $eid -- $i out of $imax ========\n";
   $edir_in   = "$edir/$eid";
   ${efiles_in}  = "${edir_in}/*${ftag}.bin";
   $nfiles_in  = `ls -1 $efiles_in | wc | awk '{print \$1}'`; chomp($nfiles_in);
   print "nfiles_in = ${nfiles_in}\n";

   if (${nfiles_in} != $nproc ) {
     #die("no $nproc files in ${efiles_in}");
     print "no $nproc files in ${efiles_in}\n";

   } else {

     $edir_done = "${dir_done}/$eid";
     $edir_done1 = "${basedir}/smooth_all/TO_COMBINE/${eid}";
     $edir_done2 = "${basedir}/smooth_all/TO_COMBINE/DONE/${eid}";
     #if ( (-e ${edir_done}) || ( (-e ${eid_done1}) || (-e ${eid_done2}) ) {
     if ( (-e ${edir_done}) || (-e ${eid_done1}) ) {
       print "--> kernel is already smoothed\n";

     } else {

       $edir_copy = "$eid";
       $edir_out  = "${edir_copy}/inout_smooth";
       $efiles_out = "${edir_out}/*${ftag}.bin";
       $nfiles_out = `ls -1 ${efiles_out} | wc | awk '{print \$1}'`; chomp($nfiles_out);
       print "nfiles_out = ${nfiles_out}\n";

       if ( ${nfiles_out} == $nproc ) {
   print "--> $nfiles_out bin files exist -- ready to go\n";

       } else {

   `mkdir -p $edir_smooth/$eid`;
   `mkdir -p $edir_copy`;
   `mkdir -p $eid/OUTPUT_FILES`;
   `mkdir -p $eid/inout_smooth`;

   $sfiles = "${edir_out}/*${ftag}_smooth.bin";
   $nsfiles  = `ls -1 $sfiles | wc | awk '{print \$1}'`; chomp($nsfiles);
   if ( $nsfiles == $nproc ) {
     print "--> $nsfiles bin files already exist -- no need to run\n";

   } else {
     # link input event kernels
     `ln -s $efiles_in $edir_out`;

     # copy run scripts and link executable file
     `cp $maindir/run.lsf $edir_copy`;
     `cp $maindir/go.bash $edir_copy`;
     `ln -s $maindir/topo $edir_copy`;
     `ln -s $maindir/smooth_sem_fun $edir_copy`;
     `ln -s $maindir/xcombine_vol_data $edir_copy`;
     `ln -s $maindir/slice_file $edir_copy`;
     `ln -s $maindir/combine.bash $edir_copy`;

     # change the label for the run
     `sed "/smooth_model/s/smooth_model/${eid}_smooth/" ${edir_copy}/go.bash > ${edir_copy}/run.tmp`;
     `mv -f ${edir_copy}/run.tmp ${edir_copy}/go.bash`;
   }
       }
     }
   }

}

#=================================================================
