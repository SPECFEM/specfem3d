#!/usr/bin/perl -w

use File::Basename;

if (@ARGV != 0 and @ARGV != 1) {die("view_model.pl tcl_file_tag\n");}

if (@ARGV == 0) {$tcl_tag = "xc_ew";}
else {$tcl_tag=$ARGV[0];}

$tcl_file = "${tcl_tag}.tcl";
$tcl_file_local ="${tcl_tag}_local.tcl";

if (not -f $tcl_file) {die("Check if sample file $tcl_file exists or not\n");}

$cmt = "CMTSOLUTION";
$par = "Par_file";
$sta = "STATIONS_ADJOINT";
$scalar_low = "-3e-7";
$scalar_high = "3e-7";

$pwd = basename($ENV{PWD});

#@temp = `global_slice_number.pl $cmt $sta $par`;

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

# color range for models (m/s)
$scalar_low  = 2000;
$scalar_high = 8000;

$iew = 0;
$ihr = 1;

if ( $iew == 1 ) {

  $xcen = $x_center;
  $zcen = $z_center;

  # SetViewUp
  $ux = 0; $uy = 0; $uz = 1;

  # normal for plane
  $nx = 0; $ny = 1; $nz = 0;

  for ($i = 1; $i <= $Ny; $i++) {
    $ycen = $utm_ymin + ($i-1)*$dyinc;
    $ycen_km = ($ycen-$utm_ymin)/1000;
    printf ("%2.2i : %.1f, %.1f, %.1f\n",$i,$xcen,$ycen,$zcen);

    # title for plot
    $title = sprintf("Lin model with Harvard model -- EW cross-section at UTM-Y = %.0f km",$ycen_km);

    # open file for replacing the cross-section line
    open(SED,">sed.txt");

    # color scale
    print SED "/kMapper SetScalarRange/s/SetScalarRange.*\$/SetScalarRange $scalar_low $scalar_high/ \n";
    print SED "/ewMapper SetScalarRange/s/SetScalarRange.*\$/SetScalarRange $scalar_low $scalar_high/ \n";

    # cross-section
    printf SED ("/ewPlane SetOrigin/s/SetOrigin.*\$/SetOrigin %.1f %.1f %.1f / \n",$xcen,$ycen,$zcen);
    print SED "/ewPlane SetNormal/s/SetNormal.*\$/SetNormal $nx $ny $nz / \n";

    # SetPosition
    $cam_dist = 50000;
    $nxd = $xcen;
    $nyd = $ycen - $cam_dist;
    $nzd = $zcen;

    # orientation
    printf SED ("/cam1 SetFocalPoint/s/SetFocalPoint.*\$/SetFocalPoint %.1f %.1f %.1f / \n",$xcen,$ycen,$zcen);
    printf SED ("/cam1 SetPosition/s/SetPosition.*\$/SetPosition %.1f %.1f %.1f / \n",$nxd,$nyd,$nzd);
    print SED "/cam1 SetViewUp/s/SetViewUp.*\$/SetViewUp $ux $uy $uz/ \n";
    print SED "/titleActor SetInput/s/SetInput.*\$/SetInput \"$title\"/ \n";

    # file name
    $filename = sprintf("${tcl_tag}_%2.2i.ps",$i);
    print SED "/writer SetFileName/s/SetFileName.*\$/SetFileName \"$filename\"/ \n";  

    # close SED file
    close(SED);

    # make a local executable file and run it to generate a PS figure
    system("sed -f sed.txt ${tcl_file} > ${tcl_file_local}");
    system("vtk ${tcl_file_local}");
    system("ps2pdf $filename");
  }
}

#---------------------------

if ( $ihr == 1 ) {

  $tcl_tag = "xc_hr_multi";
  $tcl_file = "${tcl_tag}.tcl";
  $tcl_file_local ="${tcl_tag}_local.tcl";

  # volumetric VTU files
  $file_tag = "vtu_files\\/lin_model\\/vp_low";  # KEY
  $file1 = "${file_tag}_1.vtu";
  $file2 = "${file_tag}_2.vtu";
  $file3 = "${file_tag}_3.vtu";

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
  # Why does depth z = 0 not work, so we must use z = 0.01 ?
  @dlayers = (-0.25,0.001,3,6,10,15,17,22,31,40);
  $Nz = @dlayers;

  # color bar
  $pert = 10;
  @cmean = (5.2083,5.2083,5.6315,6.1101,6.2208,6.4695,6.5064,6.5516,6.6604,7);

  $Nz = 2;
  $imin = 1;

  for ($i = $imin; $i <= $Nz; $i++) {
    #$zcen = $utm_zmax - ($i-1)*$dzinc;
    $zcen = -1000 * $dlayers[$i-1];
    $zcen_km = $zcen/1000;
    printf ("%2.2i : %.1f, %.1f, %.1f\n",$i,$xcen,$ycen,$zcen);

    # title for plot
    $title = sprintf("Lin model with Harvard model -- Depth = %.1f km",-$zcen_km);

    # open file for replacing the cross-section line
    open(SED,">sed.txt");
 
    # file names
    print SED "/kReader1 SetFileName/s/SetFileName.*\$/SetFileName $file1/ \n";
    print SED "/kReader2 SetFileName/s/SetFileName.*\$/SetFileName $file2/ \n";
    print SED "/kReader3 SetFileName/s/SetFileName.*\$/SetFileName $file3/ \n";

    # color scale
    $scalar_low  = $cmean[$i-1] * 1000 * (1 - $pert/100);
    $scalar_high = $cmean[$i-1] * 1000 * (1 + $pert/100);
    print SED "/kMapper1 SetScalarRange/s/SetScalarRange.*\$/SetScalarRange $scalar_low $scalar_high/ \n";
    print SED "/hrMapper1 SetScalarRange/s/SetScalarRange.*\$/SetScalarRange $scalar_low $scalar_high/ \n";

    print SED "/kMapper2 SetScalarRange/s/SetScalarRange.*\$/SetScalarRange $scalar_low $scalar_high/ \n";
    print SED "/hrMapper2 SetScalarRange/s/SetScalarRange.*\$/SetScalarRange $scalar_low $scalar_high/ \n";

    print SED "/kMapper3 SetScalarRange/s/SetScalarRange.*\$/SetScalarRange $scalar_low $scalar_high/ \n";
    print SED "/hrMapper3 SetScalarRange/s/SetScalarRange.*\$/SetScalarRange $scalar_low $scalar_high/ \n";

    # horizontal cross-section
    printf SED ("/hrPlane1 SetOrigin/s/SetOrigin.*\$/SetOrigin %.1f %.1f %.1f / \n",$xcen,$ycen,$zcen);
    print SED "/hrPlane1 SetNormal/s/SetNormal.*\$/SetNormal $nx $ny $nz / \n";

    printf SED ("/hrPlane2 SetOrigin/s/SetOrigin.*\$/SetOrigin %.1f %.1f %.1f / \n",$xcen,$ycen,$zcen);
    print SED "/hrPlane2 SetNormal/s/SetNormal.*\$/SetNormal $nx $ny $nz / \n";

    printf SED ("/hrPlane3 SetOrigin/s/SetOrigin.*\$/SetOrigin %.1f %.1f %.1f / \n",$xcen,$ycen,$zcen);
    print SED "/hrPlane3 SetNormal/s/SetNormal.*\$/SetNormal $nx $ny $nz / \n";

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
    $filename = sprintf("${tcl_tag}_%2.2i.ps",$i);
    print SED "/writer SetFileName/s/SetFileName.*\$/SetFileName \"$filename\"/ \n";  

    # close SED file
    close(SED);

    # make a local executable file and run it to generate a PS figure
    #print "\n ${tcl_file_local} \n";
    system("sed -f sed.txt ${tcl_file} > ${tcl_file_local}");
    system("vtk ${tcl_file_local}");
    system("ps2pdf $filename");

  }
}

#---------------------------

system("/home/wagholi2/lqy/local/PDF-CLE/bin/pdcat -r *.pdf ${tcl_tag}_all.pdf");

#---------------------------

