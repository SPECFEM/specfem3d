#!/usr/bin/perl -w

if (@ARGV < 1) {die("Usage: copy_SEM_dir.pl smodel\n")}
($smodel) = @ARGV;

#$smodel = "m1";

$basedir = "/net/sierra/raid1/carltape/socal/socal_3D/SYN";
$mdir = "$basedir/model_${smodel}";
if (not -e $mdir) {die("check if mdir $mdir exist or not\n")}

@dirs = `ls -1 -d [1-9]*`;
$neid = @dirs;

#print "\n @dirs\n";

#$neid = 1;

# loop over all directories
for ($i = 1; $i <= $neid; $i++) {

  $eid = $dirs[$i-1]; chomp($eid);
  print "==== $eid =====\n";

  $odir = "$mdir/$eid";
  if (-e $odir) {
     die("--> $odir already exists\n");
  } else {
     #`mkdir $odir`;
     #`cp $eid/SEM/* $odir`;
     `cp -r $eid/SEM $odir`;
  }
}

#=================================================================
