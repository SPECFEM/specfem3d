#!/usr/bin/perl -w

# this program converts the xy file (format as in psxy, can be segmented
# by >) to the vtk polydata format to be viewed in Paraview.
# Qinya Liu, May 2007, Caltech

if (@ARGV == 0) { die("Usage: convert_xy_to_vtk.pl xyfile\n");}

$xyfile=$ARGV[0];
$vtkfile="${xyfile}.vtk";

open(FILE,">$vtkfile");

print FILE "# vtk DataFile Version 2.0\n";
print FILE "Really cool data\n";
print FILE "ASCII\n";
print FILE "DATASET POLYDATA\n";

# figure out the total number of points
open(XY,"$xyfile") || die("Check $xyfile\n");
@xy = <XY>;
close(XY);
# record the points
$total_points = 0; $points="";
# record the lines
$nlines = 0;  $npoints = 0; $nnum = 0; $lines=""; $lines2="";

for ($i=0;$i<@xy;$i++) {
  @tmp = split(" ",$xy[$i]);
  if ($tmp[0] ne ">") {
    if (@tmp == 2) {
      $points.="$tmp[0] $tmp[1] 0.0\n";}
    elsif (@tmp == 3) {
      $points.="$tmp[0] $tmp[1] $tmp[2]\n";}
    else {die("Check this line @tmp\n");}
    $lines.="$total_points ";
    $total_points++; $npoints++;
  }
  if ($tmp[0] eq ">" or $i == @xy-1) {
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
