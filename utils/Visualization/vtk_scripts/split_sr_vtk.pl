#!/usr/bin/perl -w

#-------------------
# split_sr_vtk.pl
#
# this program input a sr.vtk file (from SPECFEM3D) and outputs three files
#   source.vtk
#   epicenter.vtk
#   receiver.vtk
#
# future modifications:
#   compute epicenter for spherical/elliptical model (assumes FLAT model only),
#   or simply have SPECFEM3D output the three files instead of sr.vtk
#
# EXAMPLE:
#   split_sr_vtk.pl sr_tmp.vtk . 1
#
#-------------------

if (@ARGV < 3) { die("Usage: split_sr_vtk.pl srfile odir nsrc\n");}

$srfile = $ARGV[0];
$odir = $ARGV[1];
$nsrc = $ARGV[2];

print "src-rec file is $srfile\n";
print "output directory for files is $odir\n";
print "number of sources is $nsrc (see CMTSOLUTION)\n";

if($ARGV[1]) {$odir = $ARGV[1]} else {$odir = "."}
#$odir = $ARGV[0];
if(not -e $odir) {$odir = "."}

# read in sr.vtk
#$srfile = "sr.vtk";
open(XY,"$srfile") || die("Check $srfile\n");
@lines = <XY>;
close(XY);

#print "\n -- $lines[0] -- $lines[1] -- $lines[2] -- $lines[3] -- $lines[4] -- $lines[5] \n";

#($xsrc,$ysrc,undef) = split(" ",$lines[5]);
$nlines = @lines;
$nrec = $nlines - $nsrc - 6;

print "nlines = $nlines, nsrc = $nsrc , nrec = $nrec\n";

# write receiver file
open(FILE1,">$odir/receiver.vtk");
print FILE1 "$lines[0]";
print FILE1 "Receiver VTK file\n";
print FILE1 "$lines[2]";
print FILE1 "$lines[3]";
print FILE1 "POINTS     $nrec float\n";
print FILE1 "@lines[$nsrc+5..$nrec+$nsrc+4]";
close(FILE1);

# write source file (hypocenter)
open(FILE2,">$odir/source.vtk");
print FILE2 "$lines[0]";
print FILE2 "Source VTK file\n";
print FILE2 "$lines[2]";
print FILE2 "$lines[3]";
print FILE2 "POINTS     $nsrc float\n";
print FILE2 "@lines[5..$nsrc+4]";
close(FILE2);

# write source file (epicenter)
# note: needs to be revised for non-planar surfaces (e.g., SPECFEM3D_GLOBE)
$zepi = "0.0";
open(FILE3,">$odir/epicenter.vtk");
print FILE3 "$lines[0]";
print FILE3 "Epicenter VTK file\n";
print FILE3 "$lines[2]";
print FILE3 "$lines[3]";
print FILE3 "POINTS     $nsrc float\n";
for ($k=5; $k <= $nsrc+4; $k=$k+1) {
   ($x,$y,undef) = split("\ ",$lines[$k]);
   print FILE3 "$x $y $zepi\n";
}
close(FILE3);

print "\n-- writing out source.vtk, epicenter.vtk, receiver.vtk in directory $odir --\n";

#-----------------------------------
