#!/usr/bin/perl -w

# this script collects the mesh database file needed by Paraview from the $jobid
# subdirectory of the local node scratch disks corresponding to a list of slices
# the scp command is:
# /scratch/$ENV{USER}/DATABASES_MPI[.$jobid]/proc%04d_$files[$j].$exts[$j]
#  Qinya Liu, May 2007, Caltech
#  (version without collecting AVS_DX files)

use POSIX;

if (@ARGV != 4 and @ARGV != 5) {die("copy_basin_database.pl slice_file lsf_machine_file filename dir [jobid]\n");}
$sfile = $ARGV[0];
$machine = $ARGV[1];
$filename = $ARGV[2];
$dir = $ARGV[3];
if (@ARGV == 5) {$jobid = ".$ARGV[4]";} else {$jobid="";}

if (not -d $dir) {die("No such dir $dir\n");}

open(FILE,"$sfile") or die("Error opening file $sfile\n");
(@slices) = <FILE>;
close(FILE);
open(FILE,"$machine") or die("Error opening file $machine\n");
($machines) = <FILE>;
close(FILE);

# obtain the hash array of slice to node mapping
(@mac) = split(" ",$machines);
$sl = 0;
for ($i=0;$i<@mac;$i=$i+2) {
  $sl = $sl + $mac[$i+1];
}

for ($i=0;$i<@mac;$i=$i+2) {
  $num_slices = $mac[$i+1];
  for ($j=0;$j < $num_slices;$j++) {
    $sl = $sl  - 1;
    $slice_to_node{$sl} = $mac[$i];
  }
}

for($i=0;$i<@slices;$i++) {
  ($slice[$i]) = split(" ",$slices[$i]);
  $node[$i] = $slice_to_node{$slice[$i]};
  print "$slice[$i], $node[$i]\n";
}

@files = ("$filename","ibool","x","y","z");
@exts = ("txt","txt","bin","bin","bin","bin","bin");


for($i=0;$i<@slices;$i++) {
  for($j=0;$j<@files;$j++){
  $string = sprintf("scp $node[$i]:/scratch/$ENV{USER}/DATABASES_MPI$jobid/proc*%04d_$files[$j].$exts[$j] $dir", $slice[$i]);
  print "$string\n";
  system("$string");
}}

