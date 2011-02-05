#!/usr/bin/perl -w

#-------------------
# split_sr_vtk.pl
#
# this program input a sr.vtk file (from SPECFEM3D) and outputs three files
#   source.vtk
#   epicenter.vtk
#   receiver.vtk
#
# NOTE: assumes there is a SINGLE SOURCE. (This can be modified later.)
#
# EXAMPLE:
#   split_sr_vtk.pl sr.vtk .
#   split_sr_vtk.pl sr.vtk
#
#-------------------

if (@ARGV == 0) { die("Usage: split_sr_vtk.pl srfile odir\n");}

$srfile = $ARGV[0];
if($ARGV[1]) {$odir = $ARGV[1]} else {$odir = "."}
#$odir = $ARGV[0];
if(not -e $odir) {$odir = "."}

# read in sr.vtk
#$srfile = "sr.vtk";
open(XY,"$srfile") || die("Check $srfile\n");
@lines = <XY>;
close(XY);

#print "\n -- $lines[0] -- $lines[1] -- $lines[2] -- $lines[3] -- $lines[4] -- $lines[5] \n";

($xsrc,$ysrc,undef) = split(" ",$lines[5]);
$nlines = @lines;
$nrec = $nlines - 6;

# write source file (hypocenter)
open(FILE1,">$odir/source.vtk");
print FILE1 "$lines[0]";
print FILE1 "Source VTK file\n";
print FILE1 "$lines[2]";
print FILE1 "$lines[3]";
print FILE1 "POINTS      1 float\n";
print FILE1 "$lines[5]";
close(FILE1);

# write source file (epicenter)
$zepi = "0.0";
open(FILE2,">$odir/epicenter.vtk");
print FILE2 "$lines[0]";
print FILE2 "Epicenter VTK file\n";
print FILE2 "$lines[2]";
print FILE2 "$lines[3]";
print FILE2 "POINTS      1 float\n";
print FILE2 "$xsrc $ysrc $zepi";
close(FILE2);

# write receiver file
$zepi = "0.0";
open(FILE3,">$odir/receiver.vtk");
print FILE3 "$lines[0]";
print FILE3 "Receiver VTK file\n";
print FILE3 "$lines[2]";
print FILE3 "$lines[3]";
print FILE3 "POINTS     $nrec float\n";
print FILE3 "@lines[6..$nlines-1]";
close(FILE3);

print "\n-- writing out source.vtk, epicenter.vtk, receiver.vtk in directory $odir --\n";

#-----------------------------------
