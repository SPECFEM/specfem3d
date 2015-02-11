#!/usr/bin/perl -w

#---------------------------------
#  view_model.pl
#  Carl Tape and Qinya Liu
#  18-June-2008
#
#  This script creates a PDF of cross-sections of a kernel.
#
#  CALLED BY: tomo_make_figs.pl
#
#  EXAMPLES:
#     view_model.pl 5 vs_m61
#
#     view_model.pl 4 vp_m16
#     view_model.pl 4 vp_m00
#     view_model.pl 5 vs_m16
#     view_model.pl 5 vs_m00
#     view_model.pl 6 vb_m16
#     view_model.pl 6 vb_m00
#     view_model.pl 7 poisson_m16
#     view_model.pl 7 poisson_m00
#     view_model.pl 8 vp_m16_m15
#     view_model.pl 8 vp_m16_m00
#     view_model.pl 9 vs_m16_m15
#     view_model.pl 9 vs_m16_m00
#     view_model.pl 10 vb_m16_m00
#
#---------------------------------

use File::Basename;

if (@ARGV < 2) {die("Usage: view_model.pl klab ktitle \n")}
($iker,$kfile) = @ARGV;

# copy files into local directory
#`cp vtu_files/lin_model_2/*vtu .`;
#`cp /media/raid/carltape/SOCAL_ADJOINT/RUNS/MODELS/m00/vtu_old/*vtu .`;
#`cp /media/raid/carltape/SOCAL_ADJOINT/RUNS/MODELS/m00/*vtu .`; $title_tag = "Model m00 for";
#`cp /media/raid/carltape/SOCAL_ADJOINT/RUNS/MODELS/m01/*vtu .`; $title_tag = "Model m01 for";
#`cp /media/raid/carltape/SOCAL_ADJOINT/RUNS/MODELS/m02/*vtu .`; $title_tag = "Model m02 for";
#`cp /media/raid/carltape/SOCAL_ADJOINT/RUNS/MODELS/m03/*vtu .`; $title_tag = "Model m03 for";
#`cp /media/raid/carltape/SOCAL_ADJOINT/RUNS/MODELS/m04/*vtu .`; $title_tag = "Model m04 for";
#`cp /media/raid/carltape/SOCAL_ADJOINT/RUNS/MODELS/m05/*vtu .`; $title_tag = "Model m05 for";
#`cp /media/raid/carltape/SOCAL_ADJOINT/RUNS/MODELS/m06/*vtu .`; $title_tag = "Model m06 for";
#`cp /media/raid/carltape/SOCAL_ADJOINT/RUNS/MODELS/m07/*vtu .`; $title_tag = "Model m07 for";
#`cp /media/raid/carltape/SOCAL_ADJOINT/RUNS/MODELS/m08/*vtu .`; $title_tag = "Model m08 for";
#`cp /media/raid/carltape/SOCAL_ADJOINT/RUNS/MODELS/m09/*vtu .`; $title_tag = "Model m09 for";
#`cp /media/raid/carltape/SOCAL_ADJOINT/RUNS/MODELS/m10/*vtu .`; $title_tag = "Model m10 for";
#`cp /media/raid/carltape/SOCAL_ADJOINT/RUNS/MODELS/m11/*vtu .`; $title_tag = "Model m11 for";
#`cp /media/raid/carltape/SOCAL_ADJOINT/RUNS/MODELS/m12/*vtu .`; $title_tag = "Model m12 for";
#`cp /media/raid/carltape/SOCAL_ADJOINT/RUNS/MODELS/m13/*vtu .`; $title_tag = "Model m13 for";
#`cp /media/raid/carltape/SOCAL_ADJOINT/RUNS/MODELS/m14/*vtu .`; $title_tag = "Model m14 for";
#`cp /media/raid/carltape/SOCAL_ADJOINT/RUNS/MODELS/m15/*vtu .`; $title_tag = "Model m15 for";
#`cp /media/raid/carltape/SOCAL_ADJOINT/RUNS/MODELS/m16/*vtu .`; $title_tag = "Model m16 for";

`cp /media/raid/carltape/SOCAL_ADJOINT/RUNS/MODELS/m61/*vtu .`; $title_tag = "Model m16-smooth for";

#$title_tag = "LN(m13\\/m12) for"; $cmax = 0.05;
#$title_tag = "LN(m16\\/m00) for"; $cmax = 0.20;

$tcl_tag = "view_model";
$tcl_file = "${tcl_tag}.tcl";
$tcl_file_local ="${tcl_tag}_local.tcl";
if (not -f $tcl_file) {die("Check if sample file $tcl_file exists or not\n");}

# remove figures and vtk files in local directory
`rm *.pdf *.ps`;

# directories
#$dir_output = "/net/sierra/raid1/carltape/socal/socal_3D/RUNS/${eid}/${smodel}/OUTPUT_${ftag}";
#$dir_mesh = "/net/sierra/raid1/carltape/socal/socal_3D/RUNS/${eid}/${smodel}/MESH_${ftag}";
#$dir_ker_lab = "\\/net\\/sierra\\/raid1\\/carltape\\/socal\\/socal_3D\\/RUNS\\/${eid}\\/${smodel}\\/MESH_${ftag}";

# kernel options
@klabs = ("kappa","mu","rho","vp","vs","vb","nu","vp","vs","vb");
@ktitles = ("BULK MODULUS","SHEAR MODULUS","DENSITY","P WAVE-SPEED","S WAVE-SPEED","BULK WAVE-SPEED","POISSON RATIO","P WAVE-SPEED","S WAVE-SPEED","BULK WAVE-SPEED");
$klab = $klabs[$iker-1];
$ktitle = $ktitles[$iker-1];
#$file_tag = "${klab}_low"; $ftag = "${klab}_low";
$file_tag = $kfile; $ftag = $kfile;
$tlab = "$title_tag $ktitle";

$pwd = basename($ENV{PWD});

# OUTPUT from hauksson_tomo_model.m
#1D averaged model:
#     depth  thickness   vp-min   vp-mean    vp-max    vs-min   vs-mean    vs-max    rho-min  rho-mean  rho-max
#         0    1.5000    2.8109    5.2083    6.2289    1.4284    2.9194    3.7731    2.1797    2.5678    2.7677
#    3.0000    3.0000    3.7365    5.6315    6.7879    1.8581    3.1706    3.8945    2.3553    2.6424    2.9085
#    6.0000    3.5000    4.4321    6.1101    7.3782    2.5197    3.4845    4.1260    2.4528    2.7407    3.0822
#   10.0000    4.5000    5.3617    6.2208    6.8249    3.0108    3.5586    4.3172    2.5936    2.7658    2.9187
#   15.0000    3.5000    6.0794    6.4695    7.1689    3.3740    3.6642    4.1778    2.7339    2.8254    3.0178
#   17.0000    3.5000    6.2422    6.5064    6.8832    3.5313    3.6656    4.1101    2.7707    2.8347    2.9349
#   22.0000    7.0000    6.4608    6.5516    6.8584    3.5813    3.6844    3.9630    2.8233    2.8461    2.9279
#   31.0000   33.5000    6.6050    6.6604    6.7609    3.7213    3.7419    3.8484    2.8599    2.8743    2.9012

# depth layers
# Why does depth z = 0 not work, so we must use z = 0.001 ?

#@dlayers = (-0.25,0.001,5,10,15,20,25,30,35,40);
@dlayers = (-0.25,0.001,2,4,6,8,10,15,20,25,30,35,40);
#@dlayers = (-0.25,0.001,1,2,3,4,5,6,7,8,9,10);    # density plots
#@dlayers = (-0.25,0.001,2.75,10.75,24,46);
$Nz = @dlayers;

# percent perturbation for each layer
@players = (15,15,15,15,15,15,10,10,10,10,10,10,5);

# mean values for each depth slice
@cmean_vp = (5.20,5.20,5.50,5.9,6.11,6.11,6.22,6.46,6.50,6.55,6.7,7.30,7.80);
@cmean_vs = (2.95,2.95,3.1,3.25,3.4,3.45,3.55,3.65,3.7,3.7,3.75,4.2,4.5);
#@cmean_vb = (3.92,3.92,4.17,4.55,4.68,4.63,4.67,4.89,4.89,4.96,5.11,5.45,5.81);

$pert = 10;
if($iker == 1) {
  $cmin = 0e10; $cmax = 10e10; $cunits = "Pa";
} elsif($iker == 2) {
  $cmin = 0e10; $cmax = 6e10; $cunits = "Pa";
} elsif($iker == 3) {
  $cmin = 1800; $cmax = 3200; $cunits = "kg m\\^-3";

} elsif($iker == 4) {
  $cmin = 1000; $cmax = 6000; $cunits = "m\\/s";
  @cmean = @cmean_vp;

} elsif($iker == 5) {
  $cmin = 1000; $cmax = 5000; $cunits = "m\\/s";
  @cmean = @cmean_vs;

} elsif($iker == 6) {
  $cmin = 1000; $cmax = 5000; $cunits = "m\\/s";
  #@cmean = @cmean_vb;
  for ($i = 1; $i <= $Nz; $i++) {
     $vp = $cmean_vp[$i-1];
     $vs = $cmean_vs[$i-1];
     $cmean[$i-1] = sqrt( $vp**2 - (4/3)*$vs**2 );
  }
  #print "\n @cmean_vp \n @cmean_vs \n @cmean \n"; die("TESTING");

} elsif($iker == 7) {
  $cmin = 0.1; $cmax = 0.4; $cunits = " ";

} elsif($iker == 8) {
  $cmin = -$cmax; $cunits = " ";

} elsif($iker == 9) {
  $cmin = -$cmax; $cunits = " ";

} elsif($iker == 10) {
  $cmin = -$cmax; $cunits = " ";

}

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

$Nz = @dlayers;
$imin = 1; $imax = $Nz;   # default
#$imin = 9; $imax = 13;  # testing
#$imin = 2; $imax = $imin;  # testing

for ($i = $imin; $i <= $imax; $i++) {

  #$zcen = $utm_zmax - ($i-1)*$dzinc;
  $zcen = -1000 * $dlayers[$i-1];
  $zcen_km = $zcen/1000;
  printf ("%2.2i : %.1f, %.1f, %.1f\n",$i,$xcen,$ycen,$zcen);
  $pert = $players[$i-1];
  $cmid = $cmean[$i-1];

  # title for plot
  $title = sprintf("%s -- Cut at z = %.1f km",$tlab,-$zcen_km);

  if($iker <= 6) {
    $clabel = sprintf("%.0f $cunits  +-  %.0f percent",$cmid*1000,$pert);
    $scalar_low  = $cmid * 1000 * (1 - $pert/100);
    $scalar_high = $cmid * 1000 * (1 + $pert/100);

  } elsif($iker == 7) {
    $clabel = "  ";
    $scalar_low  = sprintf("%.4e",$cmin);
    $scalar_high = sprintf("%.4e",$cmax);

  } else {
    $clabel = "Percent Change";
    $scalar_low  = sprintf("%.4e",$cmin);
    $scalar_high = sprintf("%.4e",$cmax);
  }

  # open file for replacing the cross-section line
  open(SED,">sed.txt");

  # file names
  print SED "/kReader1 SetFileName/s/SetFileName.*\$/SetFileName $file1/ \n";

  # color scale
  print SED "/scalarBar SetTitle/s/SetTitle.*\$/SetTitle \"$clabel\"/ \n";

  print SED "/kMapper1 SetScalarRange/s/SetScalarRange.*\$/SetScalarRange $scalar_low $scalar_high/ \n";
  print SED "/hrMapper1 SetScalarRange/s/SetScalarRange.*\$/SetScalarRange $scalar_low $scalar_high/ \n";

  # cross-section
  printf SED ("/hrPlane1 SetOrigin/s/SetOrigin.*\$/SetOrigin %.1f %.1f %.1f / \n",$xcen,$ycen,$zcen);
  print SED "/hrPlane1 SetNormal/s/SetNormal.*\$/SetNormal $nx $ny $nz / \n";

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
  $filename = sprintf("${ftag}_%2.2i.ps",$i);
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
