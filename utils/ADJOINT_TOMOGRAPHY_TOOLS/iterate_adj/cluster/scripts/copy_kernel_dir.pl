#!/usr/bin/perl -w

#-------------------------
# copies kernel output files from cluster to disk
#
# EXAMPLE: copy_kernel_dir.pl m14 all
#-------------------------

if (@ARGV < 2) {die("Usage: copy_kernel_dir.pl smodel ftag\n")}
($smodel,$ftag) = @ARGV;

$basedir = "/net/sierra/raid1/carltape/socal/socal_3D/RUNS";
if (not -e $basedir) {die("check if basedir $basedir exist or not\n")}

@dirs = `ls -1 -d [1-9]*`;
$neid = @dirs;

#print "\n -- @dirs --\n";

#$neid = 1;

# loop over all directories
for ($i = 1; $i <= $neid; $i++) {

  $eid = $dirs[$i-1]; chomp($eid);
  print "==== $eid =====\n";

  $odir = "$basedir/$eid/$smodel";
  if (not -e $odir) {die("odir $odir does not exist\n")}

  # make sure not to overwrite any directories
  $dir_output = "$odir/OUTPUT_${ftag}";
  $dir_sem = "$odir/SEM_${ftag}";
  $dir_mesh = "$odir/MESH_${ftag}";
  if (-e ${dir_output}) {die("dir_output ${dir_output} already exists\n")}
  if (-e ${dir_sem}) {die("dir_sem ${dir_sem} already exists\n")}
  if (-e ${dir_mesh}) {die("dir_mesh ${dir_mesh} already exists\n")}

  # copy files from cluster to raid disk
  #`cp -r $eid/SEM/MESH_${ftag} ${dir_mesh}`;
  `mkdir $dir_mesh`; `sleep 2s`;
  `mkdir ${dir_mesh}/collect`; `sleep 2s`;
  `cp $eid/SEM/MESH_${ftag}/*mesh ${dir_mesh}`;              # mesh
  `cp $eid/SEM/MESH_${ftag}/collect/*kernel* ${dir_mesh}/collect`;   # only get the kernel files
  `sleep 2s`;

  `cp -r $eid/SEM/SEM_${ftag} ${dir_sem}`;
  `sleep 2s`;
  `cp -r $eid/SEM/OUTPUT_${ftag} ${dir_output}`;
  `sleep 2s`;
  `cp $eid/OUTPUT_FILES/*.o ${dir_output}`;     # output file from run
}

#=================================================================
