#!/usr/bin/perl -w

#---------------------------------
#  view_kernel.pl
#  Carl Tape, after Qinya Liu
#  02-Nov-2009
#
#  This script creates a PDF of cross-sections of a kernel.
#
#  Originally, this was designed to use a three-part vtu file,
#  since the high-res files can be too large for the memory of
#  some computers.  Now, we assume that there is a single
#  vtu file (note: can be a low-res version).
#
#  CALLED BY: tomo_make_figs.pl
#
#  EXAMPLE:
#     view_kernel.pl 9818433 2 m00 all 133_9818433
#     view_kernel.pl 9818433 2 m00 T02 133_9818433
#     view_kernel.pl 9818433 2 m00 T06 133_9818433
#
#     view_kernel.pl 14236768 7 m16 T003_T010 14236768
#     view_kernel.pl 14236768 8 m16 T003_T010 14236768
#     view_kernel.pl 14236768 9 m16 T003_T010 14236768
#
#     view_kernel.pl 10992159 9 m16 T003_T010 10992159
#
#---------------------------------

use File::Basename;

if (@ARGV < 5) {die("Usage: view_kernel.pl eid iker smodel elab\n")}
($eid,$iker,$smodel,$ftag,$elab) = @ARGV;

$tcl_tag = "view_kernel";
$tcl_tag = "view_kernel_mesa";
$tcl_file = "${tcl_tag}.tcl";
$tcl_file_local ="${tcl_tag}_local.tcl";
if (not -f $tcl_file) {die("Check if sample file $tcl_file exists or not\n");}

# remove figures and vtk files in local directory
`rm *.pdf *.ps`;

# KEY COMMANDS
#$eid = 9818433;
#$smodel = "m0";
#$iker = 2;      # kappa, mu, rho, alpha, beta, rhop
#$cmax = 6e-11;  # color scale
#$cmax = 3e-10;  # color scale
#$clabel = "s\\^2  m\\^-3";
$clabel = "m\\^-3";

# directories
$dir0 = "/media/raid/carltape/SOCAL_ADJOINT/RUNS";

$dir_smodel = "$dir0/${eid}/${smodel}";
$dir_smodel = "$dir0/${eid}/${smodel}_BDK_SAVE";    # TEMPORARY
$dir_smodel = "$dir0/${eid}/${smodel}/KERNEL_ORIG";    # TEMPORARY

$dir_output = "${dir_smodel}/OUTPUT_${ftag}";
$dir_mesh = "${dir_smodel}/MESH_${ftag}";
#$dir_ker_lab = "\\/net\\/sierra\\/raid1\\/carltape\\/socal\\/socal_3D\\/RUNS\\/${eid}\\/${smodel}\\/MESH_${ftag}";
$dir_ker_lab = "\\/media\\/raid\\/carltape\\/SOCAL_ADJOINT\\/RUNS\\/${eid}\\/${smodel}\\/KERNEL_ORIG\\/MESH_${ftag}";

# check for vtu files
$nvtu0 = 1;
$nvtu  = `ls -1 ${dir_mesh}/*vtu | wc | awk '{print \$1}'`; chomp($nvtu);
if ( $nvtu < $nvtu0 ) {die("You have $nvtu files, not $nvtu0")}

# kernel options
#@klabs = ("kappa","mu","rho","alpha","beta","rhop");
#@ktitles = ("BULK MODULUS","SHEAR MODULUS","DENSITY","P-WAVE-SPEED","S-WAVE-SPEED","DENSITY");
@klabs = ("kappa_kmr","mu_kmr","rho_kmr","alpha_abr","beta_abr","rho_abr","c_cbr","beta_cbr","rho_cbr");
@ktitles = ("BULK MODULUS","SHEAR MODULUS","DENSITY","P-WAVE-SPEED","S-WAVE-SPEED","IMPEDANCE","BULK-SOUND-SPEED","SHEAR-WAVE-SPEED","IMPEDANCE");
$klab = $klabs[$iker-1];
$ktitle = $ktitles[$iker-1];
$file_tag = "${dir_ker_lab}\\/${klab}_kernel_low";
$tlab = "$eid : $smodel kernel for $ktitle";
$ftag = "${elab}_${klab}_${ftag}_${smodel}_kernel";

$pwd = basename($ENV{PWD});

# copy vtk files here
#`split_sr_vtk.pl sr.vtk`;
#`cp ${dir_output}/*vtk .`;

# number of receivers
$nrec = `wc sr.vtk | awk '{print \$1}'` - 6;
#print "\n -- $nrec -- \n"; die("testing");

# copy cmax value here
$cmax_file = "${dir_mesh}/cmax_kernel";
if(-f ${cmax_file}) {
  `cp $cmax_file .`;
  open(IN,"cmax_kernel"); $cmax = <IN>; close(IN);
} else {
  $cmax = 5e-11;   # defalt
}

$cmax = 3e-13;
#$cmax = 5e-12;   # 9875657 example

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
$dxinc = 50000;			# meters
$dyinc = 50000;			# meters
$dzinc = 5000;			# meters

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
#$file1 = "${file_tag}_1.vtu";
#$file2 = "${file_tag}_2.vtu";
#$file3 = "${file_tag}_3.vtu";

$xcen = $x_center;
$ycen = $y_center;

# normal for plane
$nx = 0; $ny = 0; $nz = 1;

# SetViewUp
$ux = 0; $uy = 1; $uz = 0;

#         0    2.8109    5.2083    6.2289    1.4284    2.9194    3.7731
#    3.0000    3.7365    5.6315    6.7879    1.8581    3.1706    3.8945
#    6.0000    4.4321    6.1101    7.3782    2.5197    3.4845    4.1260
#   10.0000    5.3617    6.2208    6.8249    3.0108    3.5586    4.3172
#   15.0000    6.0794    6.4695    7.1689    3.3740    3.6642    4.1778
#   17.0000    6.2422    6.5064    6.8832    3.5313    3.6656    4.1101
#   22.0000    6.4608    6.5516    6.8584    3.5813    3.6844    3.9630
#   31.0000    6.6050    6.6604    6.7609    3.7213    3.7419    3.8484

# depth layers
# Why does depth z = 0 not work, so we must use z = 0.001 ?
#@dlayers = (-0.25,0.001,5,10,15,20,25,30,35,40);
#@dlayers = (-0.25,0.001,2,4,6,8,10,15,20,25,30,35,40);
@dlayers = (0.001,1,2,3,4,6,8,10,15);
$Nz = @dlayers;

$imin = 1; $imax = $Nz;   # default
#$imin = 2; $imax = 9;  # testing
$imin = 4; $imax = $imin;  # testing

for ($i = $imin; $i <= $imax; $i++) {
  #$zcen = $utm_zmax - ($i-1)*$dzinc;
  $zcen = -1000 * $dlayers[$i-1];
  $zcen_km = $zcen/1000;
  printf ("%2.2i : %.1f, %.1f, %.1f\n",$i,$xcen,$ycen,$zcen);

  # title for plot
  $title = sprintf("%s -- %i stations -- Cut at z = %.1f km",$tlab,$nrec,-$zcen_km);

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

  #print SED "/kMapper2 SetScalarRange/s/SetScalarRange.*\$/SetScalarRange $scalar_low $scalar_high/ \n";
  #print SED "/hrMapper2 SetScalarRange/s/SetScalarRange.*\$/SetScalarRange $scalar_low $scalar_high/ \n";

  #print SED "/kMapper3 SetScalarRange/s/SetScalarRange.*\$/SetScalarRange $scalar_low $scalar_high/ \n";
  #print SED "/hrMapper3 SetScalarRange/s/SetScalarRange.*\$/SetScalarRange $scalar_low $scalar_high/ \n";

  # cross-section
  printf SED ("/hrPlane1 SetOrigin/s/SetOrigin.*\$/SetOrigin %.1f %.1f %.1f / \n",$xcen,$ycen,$zcen);
  print SED "/hrPlane1 SetNormal/s/SetNormal.*\$/SetNormal $nx $ny $nz / \n";

  #printf SED ("/hrPlane2 SetOrigin/s/SetOrigin.*\$/SetOrigin %.1f %.1f %.1f / \n",$xcen,$ycen,$zcen);
  #print SED "/hrPlane2 SetNormal/s/SetNormal.*\$/SetNormal $nx $ny $nz / \n";

  #printf SED ("/hrPlane3 SetOrigin/s/SetOrigin.*\$/SetOrigin %.1f %.1f %.1f / \n",$xcen,$ycen,$zcen);
  #print SED "/hrPlane3 SetNormal/s/SetNormal.*\$/SetNormal $nx $ny $nz / \n";

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
  #system("ps2pdf $filename");

}

die("TESTING");

#---------------------------

$ofile = "${ftag}_set.pdf";
#system("\\rm ${ofile}");
system("/home/carltape/bin/pdcat -r *.pdf ${ofile}");

#---------------------------
