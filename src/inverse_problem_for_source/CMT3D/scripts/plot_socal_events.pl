#!/usr/bin/perl

# this script plots given station and cmt files on
# a southern california topography map.
# Special thanks to Carl's main script and his effort
# of putting together the topo-bathymetry grid file for
# this region

use lib '/opt/seismo/lib/perl';
use GMT_PLOT;
use GMT_PLOT_SC;
use Getopt::Std;

if (@ARGV == 0) {die("plot_socal_events.pl  sta-file cmt-solution-files\n");}

$station=$ARGV[0];
if (not -f $station) {die("Check $station file\n");}
@cmt_files=@ARGV[1..$#ARGV];
foreach $file (@cmt_files) {
  if (not -f $file) {die("Check cmt file $file\n");}
}

# this controls the paper orientation;

$file="/data2/Datalib/SC/w140n40.Bathmetry.srtm.swap";
$R="-122/-114.5/32/37";
$JM="-JM7.5i";
$psfile="socal_map.ps";

open(BASH,">socal_map.bash");

print BASH "gmtset PAPER_MEDIA letter BASEMAP_TYPE plain PLOT_DEGREE_FORMAT D TICK_LENGTH 0.3c LABEL_FONT_SIZE 12 ANOT_FONT_SIZE 10  HEADER_FONT 1 ANOT_FONT 1 LABEL_FONT 1 HEADER_FONT_SIZE 18 FRAME_PEN 2p TICK_PEN 2p MEASURE_UNIT inch\n";

plot_sc_all(\*BASH,$psfile,"$JM -R$R -W1.5 -w0.5 -B2/2WesN -Na/1.0p,255/255/255 -I$file.int -C$file.cpt $file.grd -K","");

plot_psxy(\*BASH,$psfile,"-JM -R -St0.07 -W0.3 -G255/0/0","$station");

plot_psmeca(\*BASH,$psfile,"-JM -R -Sm0.20 ",@cmt_files);

plot_pstext(\*BASH,$psfile,"-JM -R -B -G250 -S1p,0","/data2/Datalib/SC/socal_topo_labs.xyz");

$evtext="";
foreach $cmt (@cmt_files) {
  ($evnm)=split(" ",`grep event $cmt | awk '{print \$3}'`);
  ($lat)=split(" ",`grep latitude $cmt | awk '{print \$2}'`);
  ($lon)=split(" ",`grep longitude $cmt | awk '{print \$2}'`);
  $evtext.="$lon $lat 7 0 4 CM $evnm\n";
}
chomp($evtext);
plot_pstext(\*BASH,$psfile,"-JM -R -B -G0/0/255 -D0/0.12","$evtext");

plot_psxy(\*BASH,$psfile,"-JX2 -R0/1/0/1 -O","");
close(BASH);
system("bash socal_map.bash 2> /dev/null");
#system("gv socal_map.ps");
