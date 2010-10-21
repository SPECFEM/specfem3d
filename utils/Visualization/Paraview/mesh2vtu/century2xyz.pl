#!/usr/bin/perl

use POSIX;

open(FILE, "century.llrm");
@lines = <FILE>;
close(FILE);
print @lines . "\n";
foreach $line (@lines) {
    ($lat,$lon,$depth,$mag) = split(/\s+/,$line);
    @pt = lat_lon_depth_2_xyz($lat, $lon, $depth);
    print "@pt $depth\n";
}

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

