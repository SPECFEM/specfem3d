#!/usr/bin/perl -w
#
# usage:
#
#  ./plot_seismos.gmt.pl OUTPUT_FILES/Y*.BHZ.semd
#
# GMT package must be installed...
use POSIX;
use Getopt::Std;


#---------------------------------------------------------------------------
## PARAMETERS

# min/max range scaling
$SCALE = 10. ;

#---------------------------------------------------------------------------

sub Usage{
  print STDERR <<END;

Usage: e.g.  ./plot_seismos.gmt.pl OUTPUT_FILES/Y*.BXZ.semd

END
exit(1);
}

@ARGV > 0 or Usage();

# find start and end time and a reasonable step
$narg = @ARGV;
$mid = int($narg/2.0);
$trace = $ARGV[$mid];
#print "trace: $trace\n";

# set region
$minmax=`minmax $trace -C  `;
chomp($minmax);

($t_start,$t_end,$min,$max) = split(" ",$minmax);

$min = $SCALE*$min;
$max = $SCALE*$max;

$region="$t_start/$t_end/$min/$max";

#print "region: $region\n";

$proj="X6/1.5";
$color="0/0/200";

open(GMT,">plot_gmtseismos.sh");
print GMT "gmtset PAPER_MEDIA letter MEASURE_UNIT inch HEADER_FONT_SIZE 14p LABEL_FONT_SIZE 16p\n";


# set output filename
$out="seis.ps";
print GMT "psbasemap -R$region -J$proj -B::.:'Time (s)':/S -K -P -Y1 > $out \n";


#################################
# plot seismograms
#################################

$offset = 8./$narg;

$counter=0;
$xoff = 0;
$yoff = 0;
foreach $file (@ARGV) {

$counter++;

$xoff=0;
$yoff=$offset;

# plots
print GMT "psxy $file -R$region -J$proj -W2/$color -X$xoff -Y$yoff -O -K >> $out \n";

}

# finishes plot with annotations
#print GMT "pstext -R$region -J$proj -N -O -K <<  END >> $out \n";
#print GMT "$t_start $max 12 0 0 LT Seismograms \n";
#print GMT "END \n";

# end ps-file
print GMT "psxy -J -R -O -P -V <<EOF >>$out\nEOF\n";

print GMT "convert $out seis.pdf \n";
print GMT "rm -f $out\n";

close(GMT);

system("sh plot_gmtseismos.sh");
system("rm -f plot_gmtseismos.sh");

print "plotted to: seis.pdf \n";

