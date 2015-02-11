#!/usr/bin/perl -w

#-------------------------
# Setup post-processing for parallel batch launch.
#
# EXAMPLE:
#   setup_batch.pl 1/100 10 1   # check that all files are there
#   setup_batch.pl 1/100 10 0   # prepare parallel batch scripts
#-------------------------

if (@ARGV < 3) {die("Usage: setup_batch.pl ncore\n")}
($inds,$ncore,$icheck) = @ARGV;

$pwd = $ENV{PWD};

#$queue = "test"; $nminutes = 60;
$queue = "normal"; $nminutes = 60;
$nmesh0 = 2;
$nbin0 = 336;

`rm *.out`;

@dirs = `ls -1 -d [1-9]*`;
$neid = @dirs;

($imin,$imax) = split("/",$inds);
if($imin < 1) {$imin = 1}
if($imin > $neid) {$imin = $neid}
if($imax > $neid) {$imax = $neid}

#print "\n -- @dirs --\n";

#$neid = 1;

if (0==1) {
  $k = 0; $p = 0;
  for ($i = 1; $i <= $neid; $i++) {
    $k = $k+1;
    if ($k % $ncore == 1) {
      $p = $p+1;
    }
    print "i = $i, k = $k, p = $p\n";
  }
  die("TESTING");
}

#----------------------

$k = 0; $p = 0;
for ($i = $imin; $i <= $imax; $i++) {

  cd_directory($pwd);
  $eid = $dirs[$i-1]; chomp($eid);
  print "==== $eid =====\n";

  if ($icheck == 0) {

    $k = $k+1;
    if ( ($k % $ncore == 1) || ($imin == $imax) ) {
      $p = $p+1;
      $rfile = sprintf("run_file_%3.3i",$p);
      open(OUT,">$rfile");

      # launch command -- NEEDS TO BE MODIFIED FOR THE LAST BATCH
      $efile = sprintf("run_file_launch_%3.3i.bash",$p);
      open(RUN,">$efile");
      print RUN "#!/bin/bash\n";
      print RUN "bsub -W $nminutes -n $ncore -o ./%J.out -a mpich_gm -q $queue -J parallel_batch.py mpirun.lsf `which mpipython.exe` /usr/bin/parallel_batch.py $rfile\n";
      close(RUN);
      `chmod +x $efile`;
    }

    # each event gets its own line
    print OUT "cd /ibrixfs1/home/carltape/ADJOINT_TOMO/iterate_adj/pangu_work/smooth_all/TO_COMBINE ; ./${eid}/go_process.bash\n";

  } else {

    # check for mesh files
    @meshfiles = glob("$eid/*.mesh");
    $nmesh = @meshfiles;

    if ($nmesh == $nmesh0) {
      print "--> $nmesh mesh files already exist\n";
    } else {
      print "$nmesh mesh files\n";
      @sfiles = `ls -1 -d $eid/inout_smooth/*smooth.bin`;
      $nbin = @sfiles;
      if ($nbin != $nbin0) {
  print "--> re-run the smoothing, since you have $nbin not $nbin0 files\n";
      } else {
        print "--> ready to set this up\n";
        `cp go_process.bash $eid`;
        #`chmod +x $eid/go_process.bash`;

        # modify the EID label for the post-processing script
        # eid=9818433
  $processfile = "$eid/go_process.bash";
  $ptemp = "$eid/go_process_temp.bash";
  `sed '/^eid=/s/^.*\$/eid=$eid/' $processfile > $ptemp`;
  `mv $ptemp $processfile`;
  `chmod +x $processfile`;
      }
    }
  }   # icheck
}

#================================================

sub cd_directory {
    my($new) = @_;
    my $old = `pwd`;
    chomp($old);
    check_directory($new);
    #print "$prog: cd $new\n";
    chdir $new;
    return($old);
}

sub check_directory {
    if(! -e $_[0] ) {
        print "Directory not found: $_[0]\n";
        exit(-1);
    }
}

#=================================================================
