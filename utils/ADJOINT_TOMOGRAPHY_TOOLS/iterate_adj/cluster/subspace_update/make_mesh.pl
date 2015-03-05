#!/usr/bin/perl -w

#-----------------------------------
# Carl Tape, 12-April-2008
# make_mesh.pl
#
# This script makes all the .mesh files for VTK plotting.
#
#-----------------------------------

$slicefile = "../topo_input/slice_file";

$idir = "INPUT/mu_run";
$odir = "OUTPUT_MODEL";
$pwd = $ENV{PWD};

# read the label that designates the normalization in the data covariance
open(IN,"$idir/dcov_tag"); $ftag = <IN>; chomp($ftag); close(IN);

# read the label that designates the type of model variable
open(IN,"$idir/ker_tag"); $parm = <IN>; chomp($parm); close(IN);

$odir = "OUTPUT_MODEL";

# read in the list of pmax values you want to use
open(IN,"$idir/pmax"); @lines = <IN>; close(IN);
$npmax = @lines;

#===================================

# one-line run file
$rfile = "run_file";
open(OUT,">$rfile");

# launch command for parallel batch
$nminutes = 30;
$queue = "normal";
$efile = "run_file_launch.bash";
open(RUN,">$efile");
print RUN "#!/bin/bash\n";
print RUN "bsub -W $nminutes -n $npmax -o ./%J.out -a mpich_gm -q $queue -J parallel_batch.py mpirun.lsf `which mpipython.exe` /usr/bin/parallel_batch.py $rfile\n";
`chmod +x $efile`;
close(RUN);

#------------------------

$imin = 1; $imax = $npmax;
#$imin = 2; $imax = $imin;

for ($i = $imin; $i <= $imax; $i++) {

  $pmax = $lines[$i-1]; chomp($pmax);
  $dmtag = sprintf("dm_${parm}_p%3.3i",$pmax);
  $ifiles = "$odir/*$dmtag*";
  $nfiles = `ls -1 $ifiles | wc | awk '{print \$1}'`; chomp($nfiles);

  print "------------------------\n";
  print "$i out of $npmax -- pmax = $pmax -- making the mesh file $dmtag...\n";
  print "nfiles ($ifiles) = $nfiles\n";
  print "$dmtag\n";

  # this will run on the head node in serial
  #print "xcombine_vol_data $slicefile $dmtag topo $odir . 0\n";
  #`xcombine_vol_data $slicefile $dmtag topo $odir . 0`;

  #-------------------------------
  # instead we setup a parallel batch script to run on the cluster nodes

  # commands to make low-res mesh
  $mfile = sprintf("make_mesh_%2.2i.bash",$i);
  open(MFILE,">$mfile");
  print MFILE "#!/bin/bash\n";
  print MFILE "cd $pwd\n";
  print MFILE "$pwd/xcombine_vol_data $slicefile $dmtag topo $odir . 0\n";
  close(MFILE);
  `chmod +x $mfile`;

  # one-line run file
  print OUT "cd $pwd ; ./$mfile\n";
}

close(OUT);

#=================================================================

