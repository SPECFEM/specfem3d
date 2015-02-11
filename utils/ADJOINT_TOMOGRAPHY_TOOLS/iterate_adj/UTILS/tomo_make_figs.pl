#!/usr/bin/perl -w

#-----------------------------------
# Carl Tape, 04-Oct-2008
# tomo_make_figs.pl
#
# This script makes the event kernel plots by calling vtk scripts.
#
# EXAMPLE (unsmoothed kernels):
#    ~/UTILS/tomo_make_figs.pl 1 1 all m16 0 0
#    ~/UTILS/tomo_make_figs.pl 1 0 all m16 0 0 (make vtu files only)
#
# EXAMPLE (smoothed kernels):
#    ~/UTILS/tomo_make_figs.pl 1 1 all m16 6 1
#    ~/UTILS/tomo_make_figs.pl 1 0 all m16 6 1 (make vtu files only)
#
#-----------------------------------

if (@ARGV < 6) {die("Usage: tomo_make_figs.pl imesh ifig ftag smodel hsmooth vsmooth\n")}
($imesh,$ifig,$ftag,$smodel,$hsmooth,$vsmooth) = @ARGV;

# to smooth or not to smooth
$ismooth = 1;
if($hsmooth==0 && $vsmooth==0) {$ismooth = 0};
$stg = sprintf("h%3.3ikm_v%3.3ikm",$hsmooth,$vsmooth);

$pwd = $ENV{PWD};
$dir_run = "/net/sierra/raid1/carltape/socal/socal_3D/RUNS";
if (not -e ${dir_run}) {die("check if dir_run ${dir_run} exist or not\n")}

# directories
$dir_vtk = "/net/denali/raid1/carltape/vtk/carl_new";      # dir to run vtk script
$dir_ker_out = "/net/sierra/raid1/carltape/results/KERNELS/kernel_${smodel}";   # dir for output figures
if (not -e ${dir_vtk}) {die("check if dir_vtk ${dir_vtk} exist or not\n")}

# directory containing smoothed event kernels
if($ismooth == 1) {
   $dir_ker_out = "/net/sierra/raid1/carltape/results/KERNELS_SMOOTH/kernel_${smodel}";
   $sdir = "${dir_run}/SMOOTH_EVENT_KERNELS/$smodel";
   if (not -e $sdir) {die("check if sdir $sdir exist or not\n")}
}
if (not -e $dir_ker_out) {die("check if dir_ker_out $dir_ker_out exist or not\n")}

# list of event IDs
#$file_eids = "/net/sierra/raid1/carltape/results/EID_LISTS/kernels_run_${smodel}";
$file_eids = "/net/sierra/raid1/carltape/results/SOURCES/socal_16/SOCAL_FINAL_CMT_v16_eid";
#$file_eids = "/net/sierra/raid1/carltape/results/WINDOWS/EIDs_pass_10_plus";
if (not -f $file_eids) {die("\n check if $file_eids exists\n")}
open(IN,$file_eids); @eids = <IN>; $nevent = @eids;

# COLOR SCALE for kernels
$cmax0 = 2.5e-11; # default value 5.0e-11
#$cmax0 = 15.0e-11;

#========================================================

# # read in list of kernels TO EXCLUDE
# $kexclude = "/net/sierra/raid1/carltape/results/KERNELS/kernel_${smodel}/kernels_exclude_${smodel}";
# open(IN,"$kexclude"); @kex = <IN>; close(IN);
# $nex = @kex;

if (0==1) {
  for ($i = 1; $i <= $nevent; $i = $i+1) {
    $eid = $eids[$i-1]; chomp($eid);
    print "$i, $eid \n";
  }
  die("testing");
}

# write the C-shell script to file
$ofile1 = "${dir_ker_out}/kernels_done_nrec";
$ofile2 = "${dir_ker_out}/kernels_done";
print "\nWriting to $ofile1 ...\n";
open(DONE,">$ofile1");

# loop over all events
$imin = 1; $imax = $nevent;   # default
#$imin = 21; $imax = $nevent;
#$imin = 29; $imax = $imin;

for ($i = $imin; $i <= $imax; $i = $i+1) {

  $eid = $eids[$i-1]; chomp($eid);
  $elab = sprintf("%3.3i_$eid",$i);

  print "------------------------------\n";
  print "$imin to $imax, $i, $eid --\n";

#    # check if event has been excluded
    $imatch = 1;
#    for ($j = 1; $j <= $nex; $j++) {
#      $keid = $kex[$j-1]; chomp($keid);
#      if ($eid == $keid) {
#        $imatch = $imatch*0.0;
#      }
#    }

  if($imatch == 0) {
    print "--> event has been excluded\n";

  } else {

  # directories
  $dir_run_eid = "${dir_run}/${eid}/${smodel}";
  #$dir_run_eid = "${dir_run}/${eid}/${smodel}_BDK_SAVE";   # TEMPORARY
  $dir_mesh    = "${dir_run_eid}/MESH_${ftag}";
  $dir_output  = "${dir_run_eid}/OUTPUT_${ftag}";

  # number of vtu files in MESH directory
  if ($ismooth == 0) {
    # unsmoothed event kernels
    $mesh_files = "${dir_mesh}/*low*.mesh"; # low-res only for now
    $vtu_files  = "${dir_mesh}/*vtu";
    $nmesh0 = 6;
    $nvtu0 = 6;
    $nmesh = `ls -1 ${mesh_files} | wc | awk '{print \$1}'`; chomp($nmesh);
    $nvtu  = `ls -1 ${vtu_files} | wc | awk '{print \$1}'`; chomp($nvtu);

  } else {
    # smoothed event kernels
    $mesh_files = "$sdir/$eid/*${stg}.mesh";
    $vtu_files  = "$sdir/$eid/*${stg}.vtu";
    $nmesh0 = 2;
    $nvtu0 = 2;
    $nmesh = `ls -1 ${mesh_files} | wc | awk '{print \$1}'`; chomp($nmesh);
    $nvtu  = `ls -1 ${vtu_files} | wc | awk '{print \$1}'`; chomp($nvtu);
  }

  # make the vtu files, if desired
  if ($imesh == 1) {
    if ( $nvtu == $nvtu0 ) {
      print "--> $nvtu VTU files exist\n";
    } else {
      print "--> $nmesh MESH files exist\n";

      # get the mesh files
      if ( $nmesh == $nmesh0) {
        print "--> making vtu files\n";
        @files = glob("${mesh_files}");
  foreach $ifile (@files) {
          $ofile = `echo $ifile | sed s/mesh/vtu/`; chomp($ofile);
          `mesh2vtu.pl -i $ifile -o $ofile`;
    print "mesh2vtu.pl -i $ifile -o $ofile\n";
  }
  #die("testing");
      }
    }
  }

  # update the number of vtu files
  $nvtu  = `ls -1 ${vtu_files} | wc | awk '{print \$1}'`; chomp($nvtu);

  # make the vtu files, if desired
  if ($ifig == 1) {
    if ( $nvtu != $nvtu0 ) {
      print "--> VTU files do not exist -- try imesh = 1\n";
    } else {

      # output directories
      $dir_fig1 = "${dir_run_eid}/FIGURES_${ftag}";
      $dir_fig2 = $dir_ker_out;
      if($ismooth == 1) {$dir_fig1 = $dir_fig2;}

      # check if figure has already been made
      @figs = glob("${dir_fig2}/*${eid}*set.pdf");
      $nfig = @figs;

      $srvtk = "${dir_output}/sr.vtk";
      $nrec = `wc $srvtk | awk '{print \$1}'` - 6;

      #if (0==1) {
      if ($nfig >= 1) {
        print "--> PDF files are already done\n";

      } else {

  # RUN output directory
        if ($ismooth == 0) {

    # first, generate the source and receiver vtk files
    #print CSH "\n/opt/seismo-util/source/perl/split_sr_vtk/split_sr_vtk.pl $srvtk ${dir_output}";
    `split_sr_vtk.pl $srvtk ${dir_output}`;

    # color scale value
    $cfile = "${dir_mesh}/cmax_kernel";

    `echo $cmax0 > $cfile`;
    #if(not -f $cfile) {
    #  $cmax0 = 1.0e-11;   # default value
    #  `echo $cmax0 > $cfile`;
    #}

    `mkdir -p ${dir_fig1}`;
    `rm ${dir_fig1}/*`; # remove all figures (for now)

    # second, make figures
    # view_kernel.pl deletes all figures in its local dir each time
    cd_directory($dir_vtk);
    `rm ${dir_vtk}/*.vtk`;
    `cp ${dir_output}/*.vtk ${dir_vtk}`;

    #`${dir_vtk}/view_kernel.pl $eid 1 $smodel $ftag $elab`;   # kappa kernel
    #`cp ${dir_vtk}/*ps ${dir_vtk}/*pdf ${dir_fig1}`;          # kappa kernel

    `${dir_vtk}/view_kernel.pl $eid 2 $smodel $ftag $elab`; # mu kernel
    `cp ${dir_vtk}/*ps ${dir_vtk}/*pdf ${dir_fig1}`; # mu kernel

    `cp ${dir_vtk}/view_kernel*pl ${dir_vtk}/view_kernel*tcl ${dir_fig1}`;

    # copy some files to the extra directory
    # the index denotes which depth cross-section to take
    $index = 3; $sti = sprintf("%2.2i",$index);
    `cp ${dir_fig1}/*set.pdf ${dir_fig1}/*${sti}.pdf ${dir_fig1}/*${sti}.ps ${dir_fig2}`;

        } else {

    cd_directory($dir_vtk);
    `rm ${dir_vtk}/*.vtk`;
    `cp ${dir_output}/*.vtk ${dir_vtk}`;

    #`${dir_vtk}/view_kernel_smooth.pl $eid 1 $smodel $ftag $elab $hsmooth $vsmooth`;   # kappa kernel
    #`cp ${dir_vtk}/*ps ${dir_vtk}/*pdf ${dir_fig1}`;                          # kappa kernel

    `${dir_vtk}/view_kernel_smooth.pl $eid 2 $smodel $ftag $elab $hsmooth $vsmooth`; # mu kernel
    `cp ${dir_vtk}/*ps ${dir_vtk}/*pdf ${dir_fig1}`;                        # mu kernel

    `cp ${dir_vtk}/view_kernel_smooth*pl ${dir_vtk}/view_kernel_smooth*tcl ${dir_fig1}`;
  }
      }

      # write the event to an output file
      print DONE "$eid $nrec\n";

    }
  }
}

}

#-----------------------------------------
close(DONE);
#system("csh -f $cshfile");

# make a file with ONLY the event IDs
`awk '{print \$1}' $ofile1 > $ofile2`;

#pdcat -r *mu_all_kernel_03.pdf mu_all_kernels_surface.pdf
#pdcat -r *kappa_all_kernel_03.pdf kappa_all_kernels_surface.pdf

print "\n";

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

#================================================
