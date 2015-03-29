#!/usr/bin/perl -w

#---------------------------------
#  view_update.pl
#  Carl Tape and Qinya Liu
#  July 29, 2008
#
#  This script creates a PDF of cross-sections of a kernel.
#
#  CALLED BY: tomo_make_figs_pmax.pl
#
#  EXAMPLE:
#     view_update.pl 1 dm00 68 0.10  0 dm_kappa_kernel_smooth_p068 bulk
#     view_update.pl 1 dm00 68 0.10  6 dm_kappa_kernel_smooth_p068_smooth_06km bulk_smooth_06km
#     view_update.pl 1 dm00 68 0.10 10 dm_kappa_kernel_smooth_p068_smooth_10km bulk_smooth_10km
#     view_update.pl 1 dm01 60 0.03 10 dm_kappa_kernel_smooth_p060_smooth_10km bulk_smooth_10km
#     view_update.pl 1 dm02 31 0.40  6 dm_kappa_kernel_smooth_p031_smooth_06km bulk_smooth_06km
#     view_update.pl 1 dm02 31 0.40 10 dm_kappa_kernel_smooth_p031_smooth_10km bulk_smooth_10km
#     view_update.pl 1 dm02 31 0.40 16 dm_kappa_kernel_smooth_p031_smooth_16km bulk_smooth_16km
#     view_update.pl 1 dm03 25 0.80  0 dm_kappa_kernel_smooth_p025 bulk
#     view_update.pl 1 dm03 25 0.80 10 dm_kappa_kernel_smooth_p025_smooth_10km bulk_smooth_10km
#     view_update.pl 1 dm04 40 0.80  0 dm_kappa_kernel_smooth_p040 bulk
#     view_update.pl 1 dm05 44 0.50  0 dm_kappa_kernel_smooth_p044 bulk
#     view_update.pl 1 dm05 44 0.50 10 dm_kappa_kernel_smooth_p044_smooth_10km bulk_smooth_10km
#     view_update.pl 1 dm06 48 0.50  0 dm_kappa_kernel_smooth_p048 bulk
#     view_update.pl 1 dm07 44 0.50 10 dm_kappa_kernel_smooth_p044_smooth_h006km_v001km bulk_smooth_h006km_v001km
#
#     view_update.pl 2 dm00 68 0.10  0 dm_mu_kernel_smooth_p068 beta
#     view_update.pl 2 dm00 68 0.10  6 dm_mu_kernel_smooth_p068_smooth_06km beta_smooth_06km
#     view_update.pl 2 dm00 68 0.10 10 dm_mu_kernel_smooth_p068_smooth_10km beta_smooth_10km
#     view_update.pl 2 dm01 60 0.10 10 dm_mu_kernel_smooth_p060_smooth_10km beta_smooth_10km
#     view_update.pl 2 dm02 31 0.80  6 dm_mu_kernel_smooth_p031_smooth_06km beta_smooth_06km
#     view_update.pl 2 dm02 31 0.80 10 dm_mu_kernel_smooth_p031_smooth_10km beta_smooth_10km
#     view_update.pl 2 dm02 31 0.80 16 dm_mu_kernel_smooth_p031_smooth_16km beta_smooth_16km
#     view_update.pl 2 dm03 25 0.80  0 dm_mu_kernel_smooth_p025 beta
#     view_update.pl 2 dm03 25 0.80 10 dm_mu_kernel_smooth_p025_smooth_10km beta_smooth_10km
#     view_update.pl 2 dm03 25 2.00 10 dm_mu_kernel_smooth_p025_smooth_10km beta_smooth_10km
#     view_update.pl 2 dm04 40 1.40  0 dm_mu_kernel_smooth_p040 beta
#     view_update.pl 2 dm04 40 2.00 10 dm_mu_kernel_smooth_p040_smooth_10km beta_smooth_10km
#     view_update.pl 2 dm05 44 1.00  0 dm_mu_kernel_smooth_p044 beta
#     view_update.pl 2 dm05 44 0.80 10 dm_mu_kernel_smooth_p044_smooth_10km beta_smooth_10km
#     view_update.pl 2 dm06 48 1.00  0 dm_mu_kernel_smooth_p048 beta
#     view_update.pl 2 dm06 48 0.80 10 dm_mu_kernel_smooth_p048_smooth_h006km_v001km beta_smooth_h006km_v001km
#     view_update.pl 2 dm07 44 1.00  0 dm_mu_kernel_smooth_p044 beta
#     view_update.pl 2 dm10 28 1.60  0 dm_mu_kernel_smooth_p028 beta
#     view_update.pl 2 dm11 40 1.60  0 dm_mu_kernel_smooth_p040 beta
#     view_update.pl 2 dm12 84 1.20  0 dm_mu_kernel_smooth_p084 beta
#     view_update.pl 2 dm13 82 1.00  0 dm_mu_kernel_smooth_p082 beta_window
#     view_update.pl 2 dm13 74 2.80  0 dm_mu_kernel_smooth_p074 beta_window
#     view_update.pl 2 dm13 78 2.40  0 dm_mu_kernel_smooth_p074 beta_window

#     view_update.pl 2 dm14 60 3.40  0 dm_mu_kernel_smooth_p060 beta_window  # 123
#     view_update.pl 2 dm14 80 2.00  0 dm_mu_kernel_smooth_p080 beta_window  # 47
#     view_update.pl 2 dm14 72 2.00  0 dm_mu_kernel_smooth_p072 beta_window  # 80
#     view_update.pl 2 dm14 80 1.80  0 dm_mu_kernel_smooth_p080 beta_window  # 40
#
#     view_update.pl 2 dm15 80 1.80  0 dm_mu_kernel_smooth_p080 beta_window
#     view_update.pl 2 dm15 90 0.80  0 dm_mu_kernel_smooth_p090 beta_window
#
#     view_update.pl 2 dm00  1 0.05  6 mu_kernel_smooth_dm beta_cg_smooth_06km
#     view_update.pl 1 dm00  1 0.05  6 kappa_kernel_smooth_dm bulk_cg_smooth_06km
#
#---------------------------------

use File::Basename;

if (@ARGV < 7) {die("Usage: view_update.pl xxx\n")}
($iker,$smodel,$pmax,$cmax,$gsmooth,$fname,$fdir) = @ARGV;

$stpmax = sprintf("%3.3i",$pmax);
$dmtag = $fname;

$tcl_tag = "view_update";
$tcl_file = "${tcl_tag}.tcl";
$tcl_file_local ="${tcl_tag}_local.tcl";
if (not -f $tcl_file) {die("Check if sample file $tcl_file exists or not\n");}

# remove figures and vtk files in local directory
`rm *.pdf *.ps`;

$clabel = " ";
$stg = sprintf("%2.2i",$gsmooth);

# directories
$dir_mesh = "/net/sierra/raid1/carltape/socal/socal_3D/RUNS/MODELS/${smodel}/${fdir}";
$dir_ker_lab = "\\/net\\/sierra\\/raid1\\/carltape\\/socal\\/socal_3D\\/RUNS\\/MODELS\\/${smodel}\\/${fdir}";

# check for vtu files
$nvtu0 = 1;
$nvtu  = `ls -1 ${dir_mesh}/${dmtag}*vtu | wc | awk '{print \$1}'`; chomp($nvtu);
#if ( $nvtu != $nvtu0 ) {die("You have $nvtu files, not $nvtu0")}

# smoothing options
if ( $gsmooth == 0 ) {$stag = ""} else {$stag = "_smooth_${stg}km"}

# kernel options
@klabs = ("bulk","beta");
@ktitles = ("BULK WAVE-SPEED","SHEAR-WAVE-SPEED");
$klab = $klabs[$iker-1];
$ktitle = $ktitles[$iker-1];
#$file_tag = "${dir_ker_lab}\\/${dmtag}${stag}";
$file_tag = "${dir_ker_lab}\\/${fname}";
$tlab = "Subspace update for $ktitle -- model $pmax";
$ftag = "${klab}_update_${fname}";

$pwd = basename($ENV{PWD});

# corners of the UTM mesh
$utm_xmin = 0.06623919273678*1e6;
$utm_xmax = 0.70520037473958*1e6;
$utm_ymin = 3.57170704093370*1e6;
$utm_ymax = 4.07495136244577*1e6;
$utm_zmax = 0.;
$utm_zmin = -60000.;
$xran = $utm_xmax - $utm_xmin;
$yran = $utm_ymax - $utm_ymin;
$zran = $utm_zmax - $utm_zmin;

# center of the mesh
$x_center = ($utm_xmin + $utm_xmax)/2;
$y_center = ($utm_ymin + $utm_ymax)/2;
$z_center = ($utm_zmin + $utm_zmax)/2;

# increment between cross sections
$dxinc = 50000;     # meters
$dyinc = 50000;     # meters
$dzinc = 5000;      # meters

# number of cuts
$Nx = int($xran/$dxinc) + 1;
$Ny = int($yran/$dyinc) + 1;
$Nz = int($zran/$dzinc) + 1;

print "\n UTM-bounds for mesh:\n";
print "$utm_xmin, $utm_xmax, $utm_ymin, $utm_ymax, $utm_zmin, $utm_zmax\n";
print "DX increment between cross-sections: $dxinc\n";
print "DY increment between cross-sections: $dyinc\n";
print "DZ increment between cross-sections: $dzinc\n";
print "Number of cross-sections: $Ny\n";

if ($Nx < 1 || $Ny < 1 || $Nz < 1) {
  die("Number of cross-sections must be at least 1\n");
}

# volumetric VTU files
$file1 = "${file_tag}.vtu";

$xcen = $x_center;
$ycen = $y_center;

# normal for plane
$nx = 0; $ny = 0; $nz = 1;

# SetViewUp
$ux = 0; $uy = 1; $uz = 0;

# depth layers
# Why does depth z = 0 not work, so we must use z = 0.001 ?
#@dlayers = (-0.25,0.001,5,10,15,20,25,30,35,40);
@dlayers = (-0.25,0.001,2,4,6,8,10,15,20,25,30,35,40);
$Nz = @dlayers;

$imin = 1; $imax = $Nz;   # default
#$imin = 2; $imax = $imin;  # testing

for ($i = $imin; $i <= $imax; $i++) {
  #$zcen = $utm_zmax - ($i-1)*$dzinc;
  $zcen = -1000 * $dlayers[$i-1];
  $zcen_km = $zcen/1000;
  printf ("%2.2i : %.1f, %.1f, %.1f\n",$i,$xcen,$ycen,$zcen);

  # title for plot
  $title = sprintf("%s -- Cut at z = %.1f km",$tlab,-$zcen_km);

  # open file for replacing the cross-section line
  open(SED,">sed.txt");

  # file names
  print SED "/kReader1 SetFileName/s/SetFileName.*\$/SetFileName $file1/ \n";
  #print SED "/kReader2 SetFileName/s/SetFileName.*\$/SetFileName $file2/ \n";
  #print SED "/kReader3 SetFileName/s/SetFileName.*\$/SetFileName $file3/ \n";

  # color scale
  print SED "/scalarBar SetTitle/s/SetTitle.*\$/SetTitle \"$clabel\"/ \n";

  $scalar_low  = sprintf("%.4e",-$cmax);
  $scalar_high = sprintf("%.4e",$cmax);

  print SED "/kMapper1 SetScalarRange/s/SetScalarRange.*\$/SetScalarRange $scalar_low $scalar_high/ \n";
  print SED "/hrMapper1 SetScalarRange/s/SetScalarRange.*\$/SetScalarRange $scalar_low $scalar_high/ \n";

  # cross-section
  printf SED ("/hrPlane1 SetOrigin/s/SetOrigin.*\$/SetOrigin %.1f %.1f %.1f / \n",$xcen,$ycen,$zcen);
  print SED "/hrPlane1 SetNormal/s/SetNormal.*\$/SetNormal $nx $ny $nz / \n";

#  print SED "/kMapper1 SetScalarRange/s/SetScalarRange.*\$/SetScalarRange $scalar_low $scalar_high/ \n";
#  print SED "/hrMapper1 SetScalarRange/s/SetScalarRange.*\$/SetScalarRange $scalar_low $scalar_high/ \n";

#  print SED "/kMapper2 SetScalarRange/s/SetScalarRange.*\$/SetScalarRange $scalar_low $scalar_high/ \n";
#  print SED "/hrMapper2 SetScalarRange/s/SetScalarRange.*\$/SetScalarRange $scalar_low $scalar_high/ \n";

#  print SED "/kMapper3 SetScalarRange/s/SetScalarRange.*\$/SetScalarRange $scalar_low $scalar_high/ \n";
#  print SED "/hrMapper3 SetScalarRange/s/SetScalarRange.*\$/SetScalarRange $scalar_low $scalar_high/ \n";

#  # cross-section
#  printf SED ("/hrPlane1 SetOrigin/s/SetOrigin.*\$/SetOrigin %.1f %.1f %.1f / \n",$xcen,$ycen,$zcen);
#  print SED "/hrPlane1 SetNormal/s/SetNormal.*\$/SetNormal $nx $ny $nz / \n";

#  printf SED ("/hrPlane2 SetOrigin/s/SetOrigin.*\$/SetOrigin %.1f %.1f %.1f / \n",$xcen,$ycen,$zcen);
#  print SED "/hrPlane2 SetNormal/s/SetNormal.*\$/SetNormal $nx $ny $nz / \n";

#  printf SED ("/hrPlane3 SetOrigin/s/SetOrigin.*\$/SetOrigin %.1f %.1f %.1f / \n",$xcen,$ycen,$zcen);
#  print SED "/hrPlane3 SetNormal/s/SetNormal.*\$/SetNormal $nx $ny $nz / \n";

  #    # color scale
  #    $scalar_low  = $cmean[$i-1] * 1000 * (1 - $pert/100);
  #    $scalar_high = $cmean[$i-1] * 1000 * (1 + $pert/100);
  #    print SED "/kMapper SetScalarRange/s/SetScalarRange.*\$/SetScalarRange $scalar_low $scalar_high/ \n";
  #    print SED "/hrMapper SetScalarRange/s/SetScalarRange.*\$/SetScalarRange $scalar_low $scalar_high/ \n";

  #    # horizontal cross-section
  #    printf SED ("/hrPlane1 SetOrigin/s/SetOrigin.*\$/SetOrigin %.1f %.1f %.1f / \n",$xcen,$ycen,$zcen);
  #    print SED "/hrPlane1 SetNormal/s/SetNormal.*\$/SetNormal $nx $ny $nz / \n";

  # SetPosition
  $cam_dist = 50000;
  $nxd = $xcen;
  $nyd = $ycen;
  $nzd = $zcen + $cam_dist;

  # orientation
  printf SED ("/cam1 SetFocalPoint/s/SetFocalPoint.*\$/SetFocalPoint %.1f %.1f %.1f / \n",$xcen,$ycen,$zcen);
  printf SED ("/cam1 SetPosition/s/SetPosition.*\$/SetPosition %.1f %.1f %.1f / \n",$nxd,$nyd,$nzd);
  print SED "/cam1 SetViewUp/s/SetViewUp.*\$/SetViewUp $ux $uy $uz/ \n";
  print SED "/titleActor SetInput/s/SetInput.*\$/SetInput \"$title\"/ \n";

  # file name
  #$filename = sprintf("${ftag}_%2.2i.ps",$i);
  $filename = sprintf("${fname}_%2.2i.ps",$i);
  print SED "/writer SetFileName/s/SetFileName.*\$/SetFileName \"$filename\"/ \n";

  # close SED file
  close(SED);

  # make a local executable file and run it to generate a PS figure
  #print "\n ${tcl_file_local} \n";
  system("sed -f sed.txt ${tcl_file} > ${tcl_file_local}");
  system("vtk ${tcl_file_local}");
  system("ps2pdf $filename");

}

#---------------------------

$ofile = "${ftag}_set.pdf";
#system("\\rm ${ofile}");
system("/home/carltape/bin/pdcat -r *.pdf ${ofile}");

#---------------------------
