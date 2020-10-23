#!/usr/bin/perl -w

#-----------------------------------
# Carl Tape, 10-Nov-2009
# setup_smooth_all.pl
#
# This script takes a template smoothing directory and sets up many copies.
#
# EXAMPLE:
#    setup_smooth_all.pl m16 168 rho_cbr_kernel
#
#-----------------------------------

if (@ARGV < 3) {die("Usage: setup_smooth_all.pl smodel nproc ftag\n")}
($smodel,$nproc,$ftag) = @ARGV;

#==========================================================
# USER INPUT: will depend on how files are organized

$basedir = "/n/home/ctape/ADJOINT_TOMO/iterate_adj/pangu_work";
$edirall = "/n/nobackup1/shaw_lab/ctape/SEM_RUNS/kernels_${smodel}_imaging";
$edirall = "${edirall}_target1";

$dtag1 = "SEM/MESH_T003_T010";
$dtag2 = "collect";

#==========================================================

$maindir = "$basedir/smooth";
if(not -e $edirall) {die("edirall $edirall does not exist");}
if(not -e $maindir) {die("maindir $maindir does not exist");}

# read in the list of kernels that you want to sum
$kernel_list = "$basedir/KERNELS_MODELS/kernels_run";
if (not -e ${kernel_list}) {die("check if kernel_list $kernel_list exist or not\n")}
open(IN,"${kernel_list}"); @eids = <IN>; close(IN);
$neid = @eids;

print "number of event kernels : $neid\n";

$imin = 1; $imax = $neid; # default
#$imin = 1; $imax = 2;
#$imin = 14; $imax = $imin;

for ($i = $imin; $i <= $imax; $i++) {

   $eid = $eids[$i-1]; chomp($eid);
   print "======== $eid -- $i out of $imax ========\n";

   # unsmoothed files
   $edir0      = "$edirall/$eid/$dtag1";
   $edir       = "$edir0/$dtag2";
   $efiles_in  = "$edir/*${ftag}.bin";
   $nfiles_in  = `ls -1 ${efiles_in} | wc | awk '{print \$1}'`; chomp($nfiles_in);
   if ($nfiles_in == $nproc) {$boolin = 1;} else {$boolin = 0;}
   print "nfiles_in = ${nfiles_in}\n";

   # smoothed files
   $efiles_out = "$edir/*${ftag}_smooth.bin";
   $nfiles_out = `ls -1 ${efiles_out} | wc | awk '{print \$1}'`; chomp($nfiles_out);
   if ($nfiles_out == $nproc) {$boolout = 1;} else {$boolout = 0;}
   print "nfiles_out = ${nfiles_out}\n";

   # local run directory
   $edir_local = "$eid";

   # Four possible cases:
   #  1. both smoothed and unsmoothed are done
   #  2. smoothed done, but missing unsmoothed
   #  3. neither unsmoothed nor smoothed done
   #  4. unsmoothed done, but not smoothed ==> SETUP

   if ($boolin==1 && $boolout==1) {

     # check if mesh files have been made
     $mesh_local = "${edir_local}/${ftag}_smooth*";
     $mesh_out   = "$edir0/${ftag}_smooth*";
     $nmesh_local = `ls -1 ${mesh_local} | wc | awk '{print \$1}'`; chomp($nmesh_local);
     $nmesh_out = `ls -1 ${mesh_out} | wc | awk '{print \$1}'`; chomp($nmesh_out);
     if ($nmesh_out == 1) {
       print "--> DONE\n";
     } elsif ( $nmesh_out==0 && $nmesh_local==1 ) {
       `cp $mesh_local $edir0`; `sleep 1s`;
     } elsif ( $nmesh_out==0 && $nmesh_local==0 ) {
       print "--> STILL NEED TO COMBINE BIN FILES\n";
     }

   } elsif ($boolin==0 && $boolout==1) {
     print "missing unsmoothed files, but smoothed files are present\n";

   } elsif ($boolin==0 && $boolout==0) {
     print "no unsmoothed files or smoothed files\n";

   } else {
     print "unsmoothed files present but no smoothed files\n";

     if (-e $edir_local) {
       print "--> READY TO EXECUTE SMOOTHING\n";
     } else {

       `mkdir -p ${edir_local}`;
       `mkdir -p ${edir_local}/OUTPUT_FILES`;
       `ln -s $edir ${edir_local}/inout_smooth`; # KEY LINK
       # NOTE: the link leaves the danger of over-writing input files

       # copy run scripts and link executable file
       `cp $maindir/run* $edir_local`;
       `cp $maindir/go.bash $edir_local`;
       `ln -s $maindir/topo $edir_local`;
       `ln -s $maindir/smooth_sem_fun $edir_local`;
       `ln -s $maindir/xcombine_vol_data $edir_local`;
       `ln -s $maindir/slice_file $edir_local`;
       `ln -s $maindir/combine.bash $edir_local`;

       # change the label for the run
       `sed "/smooth_model/s/smooth_model/${eid}S/" ${edir_local}/go.bash > ${edir_local}/run.tmp`;
       `mv -f ${edir_local}/run.tmp ${edir_local}/go.bash`;
     }
   }

}

#=================================================================
