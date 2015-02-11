#!/usr/bin/perl -w

#-----------------------------------
# Carl Tape, 29-Jan-2009
# finish_smooth_all.pl
#
# This script copies directories with binary files of a specified parameter
# from one place to a local directory on pangu.
#
# EXAMPLE:
#    finish_smooth_all.pl m15 168 mu_kernel 6 1
#    finish_smooth_all.pl m15 168 kappa_kernel 6 1
#
#-----------------------------------

if (@ARGV < 5) {die("Usage: setup_smooth_all.pl xxx\n")}
($smodel,$nproc,$ftag,$hsmooth,$vsmooth) = @ARGV;

$pwd = $ENV{PWD};

$odir0 = "/net/sierra/raid1/carltape/socal/socal_3D/RUNS/SMOOTH_EVENT_KERNELS";
$odir1 = "$odir0/$smodel";

#$basedir = "..";
$basedir = "/ibrixfs1/home/carltape/ADJOINT_TOMO/iterate_adj/pangu_work";
$pdir0 = "$basedir/KERNELS_MODELS/event_kernels_smooth/kernel_${smodel}";
$sdir0 = "$basedir/smooth_all/TO_COMBINE/DONE";

if(not -e $pdir0) {die("pdir0 $pdir0 does not exist");}
if(not -e $sdir0) {die("sdir0 $sdir0 does not exist");}

# read in the list of kernels that you want to sum
$kernel_list = "$basedir/KERNELS_MODELS/kernels_run";
if (not -e ${kernel_list}) {die("check if kernel_list $kernel_list exist or not\n")}
open(IN,"${kernel_list}"); @eids = <IN>; close(IN);
$neid = @eids;

# label denoting horizontal and vertical smoothing
$stg = sprintf("h%3.3ikm_v%3.3ikm",$hsmooth,$vsmooth);

print "number of event kernels : $neid\n";
print "smoothing label : $stg\n";

if (0==1) {
  for ($i = 1; $i <= $neid; $i++) {
    $eid = $eids[$i-1]; chomp($eid);
    print "$i, $eid\n";
  }
  die("testing");
}

#===================================

$imin = 1; $imax = $neid; # default
#$imin = 1; $imax = 10;
#$imin = 1; $imax = $imin;

for ($i = $imin; $i <= $imax; $i++) {

  $eid = $eids[$i-1]; chomp($eid);
  print "======== $eid -- $i to $imax ========\n";
  $pdir = "$pdir0/$eid";
  $sdir = "$sdir0/$eid";
  $odir = "$odir1/$eid";

  # (1) check if job is done
  # (2) check if mesh files are off pangu
  # (3) check if mesh files are in smooth_all

  $mfiles = "${ftag}_smooth_${stg}.mesh";
  $bfiles = "*${ftag}_smooth.bin";

  $mfile_out = `ls -1 $odir/$mfiles | wc | awk '{print \$1}'`; chomp($mfile_out);
  #$bfile_out = `ls -1 $pdir/$bfiles | wc | awk '{print \$1}'`; chomp($bfile_out);

  $ofile = `ls -1 $eid/OUTPUT_FILES/*.o`; chomp($ofile);

  if ( 0==1 ) {
  #if ( not -f $ofile ) {
    print "--> job has not been run, or is running, or has finished\n";

  } else {

    if ( ${mfile_out} == 1 ) {
      print "--> ${mfile_out} mesh file exists -- $odir/$mfiles\n";

      # check that there are 168 files there
      $ncheck = `ls -1 $pdir/$bfiles | wc | awk '{print \$1}'`; chomp($ncheck);
      if ($ncheck != $nproc) {
  die("ncheck ($ncheck) not equal to nproc ($nproc)");
      }

    } else {

      # low-res mesh files
      $mfiles2 = "$sdir/$mfiles";
      $mfile_out2 = `ls -1 $mfiles2 | wc | awk '{print \$1}'`; chomp($mfile_out2);

      # high-res binary files
      $bfiles2 = "$sdir/inout_smooth/$bfiles";
      $bfile_out2 = `ls -1 $bfiles2 | wc | awk '{print \$1}'`; chomp($bfile_out2);

      if ( (${mfile_out2} == 1) && (${bfile_out2} == $nproc) ) {
  # (1) copy the mesh files off pangu
  # (2) copy the bin files to event_kernels_smooth
        # NOTE: cp is a safer operation than mv -- especially for the sensitive IBRIX system
  print "--> copying the mesh and binary files\n";
  `mkdir -p $odir`;
  `cp $mfiles2 $odir`;
        `sleep 3s`;
        `mkdir -p $pdir`;
        if (0==1) {
    `cp $bfiles2 $pdir`;
  } else {
    `mv $bfiles2 $pdir`;
  }
        `sleep 3s`;

        # check that all the binary files made it
        # (And copy any that were left behind.)
        # (This absurd section of code is due to the faulty filesystem IBRIX at Caltech.)
        $bfiles3 = "$pdir/$bfiles";
        $bfile_out3 = `ls -1 $bfiles3 | wc | awk '{print \$1}'`; chomp($bfile_out3);
        if ( ${bfile_out3} != $nproc ) {
    #die("only ${bfile_out3} out of $nproc made it into $pdir");
    print "only ${bfile_out3} out of $nproc made it into $pdir\n";
    @leftfiles = glob("$bfiles2");
          print "leftfiles are @leftfiles\n";
    `cp @leftfiles $pdir`;
    `sleep 3s`;

          # check the numbers again
    $bfiles3 = "$pdir/$bfiles";
    $bfile_out3 = `ls -1 $bfiles3 | wc | awk '{print \$1}'`; chomp($bfile_out3);
    if ( ${bfile_out3} != $nproc ) {
      die("only ${bfile_out3} out of $nproc made it into $pdir");
    }
  }

      } else {
        if( ${mfile_out2} == 1 ) {
          die("mesh file was made but only ${bfile_out2} out of $nproc binary files");
        }
  if ( ${bfile_out2} == $nproc ) {
    die("$nproc binary files but no $mesh file\n");
  }
        if ( (${mfile_out2} == 0) && (${bfile_out2} == 0) ) {
    print "--> this event has not been setup yet\n";
  }
      }
    }
  }
}

#================================================

#sub cd_directory {
#    my($new) = @_;
#    my $old = `pwd`;
#    chomp($old);
#    check_directory($new);
#    #print "$prog: cd $new\n";
#    chdir $new;
#    return($old);
#}

#sub check_directory {
#    if(! -e $_[0] ) {
#        print "Directory not found: $_[0]\n";
#        exit(-1);
#    }
#}

#=================================================================

