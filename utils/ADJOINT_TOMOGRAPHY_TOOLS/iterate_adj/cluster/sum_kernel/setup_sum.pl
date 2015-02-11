#!/usr/bin/perl -w

#-----------------------------------
# Carl Tape, 10-Nov-2009
# setup_sum.pl
#
# Example script for setting up this directory; this will vary for each user.
#
# EXAMPLE:
#    setup_sum.pl 168 rho_cbr_kernel
#
#-----------------------------------

if (@ARGV < 2) {die("Usage: setup_sum.pl xxx\n")}
($nproc,$ktag) = @ARGV;

#-------------------------------------
# USER INPUT

# directory containing the SEM runs
$basedir = "/n/nobackup1/shaw_lab/ctape/SEM_RUNS/kernels_m16_imaging";
if(not -e $basedir) {die("basedir $basedir does not exist");}

$dtag = "SEM/MESH_T003_T010/collect";

# defaults
$idir = "INPUT_KERNELS";
$kernel_list = "kernels_run";
$kernel_list_sub = "${kernel_list}_sub";
$otag = "ftags";

#-------------------------------------

# write kernel tag to be read in
open(OUT,">$otag"); print OUT "$ktag"; close(OUT);

# delete all symbolic links in input directory
`rm -rf $idir/*`; `sleep 1s`;

# read in the list of kernels that you want to sum
if (not -e ${kernel_list}) {die("check if kernel_list $kernel_list exist or not\n")}
open(IN,"${kernel_list}"); @eids = <IN>; close(IN);
$neid = @eids;

# open list for subset of EIDs
open(KEIDS,">${kernel_list_sub}");

$imin = 1; $imax = $neid; # default
#$imin = 1; $imax = 2;
#$imin = 18; $imax = $imin;

for ($i = $imin; $i <= $imax; $i++) {

   $eid = $eids[$i-1]; chomp($eid);
   print "======== $eid -- $i out of $imax ========\n";
   $edir_in = "$basedir/$eid/$dtag";
   if(not -e $edir_in) {die("edir_in $edir_in does not exist");}

   # print EID to file
   print KEIDS "$eid\n";

   # link directory containing .bin files to the local directory
   `ln -s $edir_in $idir/$eid`;

   # check that all the linked files exist
   $efiles = "$idir/$eid/*${ktag}*";
   $nfiles = `ls -1 $efiles | wc | awk '{print \$1}'`; chomp($nfiles);
   print "nfiles = ${nfiles} for $ktag\n";
   if ($nfiles != $nproc) {die("you do not have $nproc $ktag files setup");}

}

close(KEIDS);

# instructions to make mesh file
print "to combine into a low-res mesh, use this command:\n";
print "   xcombine_vol_data ../topo_input/slice_file $ktag topo OUTPUT_SUM . 0\n";

#=================================================================
