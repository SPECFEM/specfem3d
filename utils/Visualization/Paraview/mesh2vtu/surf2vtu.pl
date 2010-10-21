#!/usr/bin/perl

use Getopt::Std;
use File::Basename;

$progname = basename($0);
$surf2vtu = "/opt/seismo-util/source/mesh2vtu/surf2vtu";
$REQLIBS  = "env LD_LIBRARY_PATH=/opt/seismo-util/lib/vtk";

sub Usage {
    print STDERR <<END;
Usage: $progname -i input-file -o output-file
    Takes an input file (binary) with a number of points and a number of cells
    and transforms them into an unstructured grid file

    -i input-file (Binary file)
    -o output-file (XML Unstructured Grid File)

    Input Binary files have this structure:
      number_of_points          integer (4 bytes)
      x_1, y_1, z_1, scalar_1   4 reals (4 bytes each)
      ...
      x_n, y_n, z_n, scalar_n   4 reals (4 bytes each)
      number_of_cells           integer (4 bytes)
      cell_1 (four points)      4 integers (4 bytes each)
      ...
      cell_n                    4 integers (4 bytes each)

    This is a wrapper around surf2vtu

    Brian Savage 6/26/2004
    Qinya Liu 9/29/2005 modified to deal with quad elements
    
END
    exit(-1);
}

if(@ARGV == 0) {
    Usage ();
}

if(!getopts('i:o:')){die "Check input paramaters \n";}

if(!defined($opt_i)) {
    die "$progname: Must specify input file -i input-file\n";
}
if(!defined($opt_o)) {
    die "$progname: Must specify output file -o output-file\n";
}
#print "$REQLIBS $surf2vtu $opt_i $opt_o\n";
system("$REQLIBS $surf2vtu $opt_i $opt_o");

1;
