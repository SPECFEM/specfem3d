#!/usr/bin/perl

use Getopt::Std;
use File::Basename;
use POSIX;

$progname = basename($0);
$ugrid    = "/opt/seismo-util/source/mesh2vtu/ugrid_pts";
$REQLIBS  = "env LD_LIBRARY_PATH=/opt/seismo-util/lib/vtk";

sub Usage {
    print STDERR <<END;
Usage: $progname -i input-file -o output-file
    Takes an input file (ascii) with a number of points
    and transforms them into an unstructured grid file

    -i input-file (ascii file)
    -o output-file (XML Unstructured Grid File)
    -L - input is in lon lat depth scalar

    Input ascii files have this structure:
      number_of_points          
      x_1 y_1 z_1 scalar_1   
      ...
      x_n y_n z_n scalar_n   

    This is a wrapper around $ugrid

    Brian Savage 6/26/2004
    
END
    exit(-1);
}

if(@ARGV == 0) {
    Usage ();
}

if(!getopts('i:o:L')){die "Check input paramaters \n";}

if(!defined($opt_i)) {
    die "$progname: Must specify input file -i input-file\n";
}
if(!defined($opt_o)) {
    die "$progname: Must specify output file -o output-file\n";
}
#print "$REQLIBS $ugrid $opt_i $opt_o\n";
if($opt_L){
  open(FILE, "$opt_i");
  @lines = <FILE>;
  close(FILE);
  print "Converting to xyzm\n";
  open(OUT,">points2vtu.tempfile");
  #  print @lines . "\n";
  print OUT "@lines[0]";
  foreach $line (@lines[1..$#lines]) {
    ($lon,$lat,$depth,$mag) = split(/\s+/,$line);
    @pt = lat_lon_depth_2_xyz($lat, $lon, $depth);
    print OUT "@pt $mag\n";
  }
  print "Running ugrid\n";
  system("$REQLIBS $ugrid points2vtu.tempfile $opt_o");
}
else{
  system("$REQLIBS $ugrid $opt_i $opt_o");
}

1;
sub lat_lon_depth_2_xyz {
    my($lat, $lon, $depth) = @_;
    my($PI, $D2R, $theta, $phi, $r0, $r, $x, $y, $z);

    $R_EARTH_KM = 6371.0;
    $PI = 3.141592653589793;
    $D2R = $PI/180.0;

    $theta = ($PI/2.0) - atan(0.99329534*tan($lat*$D2R));
    $phi = $lon * $D2R;
    $r0 = 1.0;
    
    $r = ($R_EARTH_KM - $depth) / $R_EARTH_KM;
    $x = $r * sin($theta) * cos($phi);
    $y = $r * sin($theta) * sin($phi);
    $z = $r * cos($theta);
    return($x, $y, $z);
}
