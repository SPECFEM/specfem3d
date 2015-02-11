#!/usr/bin/perl -w

#-----------------------------------
# Carl Tape, 15-Oct-2008
# tomo_make_figs_pmax.pl
#
# This script plots model updates from the subspace method
# by calling vtk scripts.
#
# EXAMPLE:
#    ~/UTILS/tomo_make_figs_pmax.pl 1 0 dm15 0.80 0.0 mu_kernel_smooth beta_window
#
#    ~/UTILS/tomo_make_figs_pmax.pl 1 0 dm08 0.80 0.0 mu_kernel_smooth beta
#    ~/UTILS/tomo_make_figs_pmax.pl 1 0 dm08 0.40 0.0 kappa_kernel_smooth bulk
#
#    ~/UTILS/tomo_make_figs_pmax.pl 1 0 dm00 0.1 0.0 mu_kernel_smooth beta
#    ~/UTILS/tomo_make_figs_pmax.pl 1 1 dm00 0.1 0.0 mu_kernel_smooth beta
#
#    ~/UTILS/tomo_make_figs_pmax.pl 1 0 dm00 0.05 0.0 kappa_kernel_smooth bulk
#    ~/UTILS/tomo_make_figs_pmax.pl 1 1 dm00 0.05 0.0 kappa_kernel_smooth bulk
#
#-----------------------------------

if (@ARGV < 7) {die("Usage: tomo_make_figs_pmax.pl imesh smodel\n")}
($imesh,$ifig,$smodel,$cmax,$gsmooth,$parm,$fdir) = @ARGV;

$pwd = $ENV{PWD};

# directories
# dir_vtk -- to run vtk script
# dir_ker_out -- collect copies of the figures
$dir_vtk = "/net/denali/scratch1/carltape/vtk/carl_new";
$dir_ker_out = "/net/sierra/raid1/carltape/socal/socal_3D/RUNS/MODELS/${smodel}/${fdir}";
$dir_mesh = $dir_ker_out;
if (not -e $dir_vtk) {die("check if dir_vtk $dir_vtk exist or not\n")}
if (not -e $dir_ker_out) {die("check if dir_ker_out $dir_ker_out exist or not\n")}

# read in the list of pmax values you want to use
$pfile = "${dir_ker_out}/pmax";
open(IN,$pfile); @plines = <IN>; close(IN);
if (not -f $pfile) {die("check if pfile $pfile exist or not\n")}
$npmax = @plines;

# read in the list of cmax values you want to use
#$pfile = "${dir_ker_out}/cmax";
#open(IN,$pfile); @clines = <IN>; close(IN);
#$npmax = @clines;

$imin = 1; $imax = $npmax;
#$imin = 10; $imax = $npmax;
#$imin = 5; $imax = $imin;

for ($i = $imin; $i <= $imax; $i++) {

  $pmax = $plines[$i-1]; chomp($pmax);
  $dmtag = sprintf("dm_${parm}_p%3.3i",$pmax);
  #($cmax,undef) = split(" ",$clines[$i-1]); chomp($cmax);
  #$cmax = 0.5;
  print "--------------------------------\n";
  print "$i out of $npmax -- pmax = $pmax -- cmax = $cmax \n";

  $fname = $dmtag;

  # files
  $mesh_files = "${dir_mesh}/${dmtag}*mesh"; # low-res only for now
  $vtu_files  = "${dir_mesh}/${dmtag}*vtu";

  # number of vtu files in MESH directory
  $nmesh0 = 1;
  $nvtu0 = 1;
  $nmesh = `ls -1 ${mesh_files} | wc | awk '{print \$1}'`; chomp($nmesh);
  $nvtu  = `ls -1 ${vtu_files} | wc | awk '{print \$1}'`; chomp($nvtu);

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
      $dir_fig1 = "${dir_ker_out}/FIGURES";
      #`rm -rf $dir_fig1/*`;

      # second, make figures
      # view_update.pl deletes all figures in its local dir each time
      cd_directory($dir_vtk);

      `${dir_vtk}/view_update.pl 2 $smodel $pmax $cmax $gsmooth $fname $fdir`;  # beta update
      `cp ${dir_vtk}/*ps ${dir_vtk}/*pdf ${dir_fig1}`;    # beta update
      `cp ${dir_vtk}/view_update* ${dir_fig1}`;

      #`${dir_vtk}/view_update.pl 1 $smodel $pmax $cmax $gsmooth $fname $fdir`;  # bulk update
      #`cp ${dir_vtk}/*ps ${dir_vtk}/*pdf ${dir_fig1}`;    # bulk update
      #`cp ${dir_vtk}/view_update* ${dir_fig1}`;

      # copy some files to the extra directory
      # the index denotes which depth cross-section to take
      #$index = 3; $sti = sprintf("%2.2i",$index);
      #`cp ${dir_fig1}/*set.pdf ${dir_fig1}/*${sti}.pdf ${dir_fig1}/*${sti}.ps ${dir_fig2}`;
    }
  }
}

#-----------------------------------------

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
