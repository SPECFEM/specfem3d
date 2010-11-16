#!/usr/bin/perl

# this script figures out the minimum set of slices to collect
# in order to create a source-receiver cross-section in Paraview.
# Use information from Par_file, output_solver.txt, and write
# to slice_file
# Qinya Liu, Caltech, May 2007

use POSIX;

if (@ARGV != 3) {die("usage: slice_number.pl Mesh_Par_file output_solver.txt slice_file\n");}

$par_file = $ARGV[0];
$output = $ARGV[1];
$slice_file = $ARGV[2];

if (not -f $par_file) {die("check if $par_file exists or not\n");}
if (not -f $output) {die("check if $output exists or not\n");}

# utm conversion
# for installation: see utils/Visualization/sph2utm_Qinya_Liu.tar.gz
$sph2utm="/opt/seismo-util/bin/sph2utm";
if (! -e $sph2utm) {die(" No $sph2utm file\n");}

# gets parameters
($nproc) = split(" ",`grep NPROC_XI $par_file | cut -d = -f 2 `);
if ( $nproc < 1 ) {die ("check if NPROC_XI is > 1 in Mesh_Par_file\n");}

($lat1) = split(" ",`grep LATITUDE_MIN $par_file | cut -d = -f 2 `);
($lon1) = split(" ",`grep LONGITUDE_MIN $par_file | cut -d = -f 2 `);
($lat2) = split(" ",`grep LATITUDE_MAX $par_file | cut -d = -f 2 `);
($lon2) = split(" ",`grep LONGITUDE_MAX $par_file | cut -d = -f 2 `);

# UTM coordinates of the corners of the model
(undef,undef,$x1,undef,undef,$y1) = split(" ",`echo $lon1 $lat1 | $sph2utm | grep x`);
(undef,undef,$x2,undef,undef,$y2) = split(" ",`echo $lon2 $lat2 | $sph2utm | grep x`);

# source and receiver UTM coordinates
(undef,undef,$xs) = split(" ",`grep "^          UTM x" $output`);
(undef,undef,$ys) = split(" ",`grep "^          UTM y" $output`);
(undef,undef,undef,$xr) = split(" ",`grep "^         original UTM x:" $output`);
(undef,undef,undef,$yr) = split(" ",`grep "^         original UTM y:" $output`);

#print "$x1,$xs,$x2\n$y1,$ys,$y2\n";
if( "$xs" == "" ){ die("source UTM x not found\n");}
if( "$ys" == "" ){ die("source UTM y not found\n");}
if( "$xr" == "" ){ die("receiver UTM x not found\n");}
if( "$yr" == "" ){ die("receiver UTM y not found\n");}

# approximate slice location of the source sslice (nxs, nys), and the receiver rslice(nxr,nyr)
$dx = ($x2-$x1)/$nproc;
$dy = ($y2-$y1)/$nproc;
$nxs = floor(($xs-$x1)/$dx); $nys = floor(($ys-$y1)/$dy); $sslice = $nys * $nproc + $nxs;
$nxr = floor(($xr-$x1)/$dx); $nyr = floor(($yr-$y1)/$dy); $rslice = $nyr * $nproc + $nxr;

# double check with output files
($s_slice) = split(" ",`grep "source located in slice" $output | cut -c 36- `);
(undef,undef,$r_slice) = split(" ",`grep "^  in slice" $output`);
if ($sslice != $s_slice || $rslice != $r_slice) {die("Check source and receiver slice numbers:\n  source: $s");}

print "NPROC = $nproc;  source slice = $s_slice; receiver slice = $r_slice\n";

# interpolate along the source-receiver line and assign the corresponding slice number
$ninc = 100;
for ($i=0;$i<$ninc+1;$i++) {
  $ddx = ($xr-$xs)/$ninc; $x = $xs + $ddx * $i;
  if ($xr != $xs) {$y = $ys + ($yr-$ys)/($xr-$xs) * ($x-$xs);} else {$y = $ys + ($yr-$ys)/$ninc * $i;}
  $nx = floor(($x-$x1)/$dx); $ny = floor(($y-$y1)/$dy);
  $nslice = $ny * $nproc + $nx;
  push(@nslice,$nslice);
}

# sort and compact the slice array
@nslice = sort(@nslice);

for ($i=0;$i<@nslice;$i++) {
  $all{$nslice[$i]} = $nslice[$i];
}

@keys = sort{$a <=> $b}(keys(%all));
print "slice numbers are @keys\n";

# write the slice number file
open(SLICE,">$slice_file");
foreach $key (@keys) {
  print SLICE "$key\n";
}
close(SLICE);
