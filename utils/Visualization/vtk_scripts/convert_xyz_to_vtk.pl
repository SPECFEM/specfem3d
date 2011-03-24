#!/usr/bin/perl -w

# convert_xyz_to_vtk.pl
# This program converts an xy file or xyz file to a vtk file.
# It will also work with segmented files with ">" as the delimiter.
# Qinya Liu, May 2007, Caltech

if (@ARGV == 0) { die("Usage: convert_xyz_to_vtk.pl xyzfile\n");}

$xyzfile=$ARGV[0];
$vtkfile="${xyzfile}.vtk";

open(FILE,">$vtkfile");

print FILE "# vtk DataFile Version 2.0\n";
print FILE "XYZ data\n";
print FILE "ASCII\n";
print FILE "DATASET POLYDATA\n";

# figure out the total number of points
open(XYZ,"$xyzfile") || die("Check $xyzfile\n");
@xyz = <XYZ>;
close(XYZ);
# record the points
$total_points = 0; $points="";
# record the lines
$nlines = 0;  $npoints = 0; $nnum = 0; $lines=""; $lines2="";

for ($i=0;$i<@xyz;$i++) {
  @tmp = split(" ",$xyz[$i]);
  if ($tmp[0] ne ">") {
    # if input has 2 columns, set 3rd column to zero (flat earth scenario, z up, e.g. UTM)
    if (@tmp == 2) {
      $points.="$tmp[0] $tmp[1] 0.0\n";}
    elsif (@tmp == 3) {
      $points.="$tmp[0] $tmp[1] $tmp[2]\n";}
    else {die("Check this line @tmp\n");}
    $lines.="$total_points ";
    $total_points++; $npoints++;
  }
  if ($tmp[0] eq ">" or $i == @xyz-1) {
    # change to a new line
    $nlines ++ ;
    $lines2.="$npoints $lines\n";
    $nnum += $npoints+1;
    $npoints=0; $lines="";
  }
}
print FILE "POINTS $total_points float\n";
print FILE "$points";
print FILE "LINES $nlines $nnum\n";
print FILE "$lines2";
close(FILE);
