#!/usr/bin/perl -w
use lib '/opt/seismo-util/lib/perl';
use CMT_TOOLS;

@frames = (6..6);

$dx = 0.01458333333333;
$dy = 0.01197916666667;
$R = "-R-120.3/-114.7/32.2/36.8";
$JM = "-JM5i";
$B = " -B2/2wesn ";
$null = 127.5000;

$datalib = "/opt/seismo-util/data/datalib_SC";
$int_file = "$datalib/topo_cal.int";
$fault_file = "$datalib/jennings.xy";
$cmt_file = "CMTSOLUTION";
$cpt_file = "disp.cpt";

($elat,$elon) = get_cmt_location($cmt_file);

open(CSH,">convert_tiff.csh");

print CSH "gmtset BASEMAP_TYPE plain ANOT_FONT_SIZE 9 HEADER_FONT_SIZE 15\n";
print CSH "makecpt -Cpolar -T0/255/2 -Z -V > $cpt_file\n";
foreach $frame (@frames) {
  if ($frame < 10) {$frame = "0$frame";}
  $file = "gmt_movie_0000$frame";
  $time = $frame * 200 * 0.009;
  $xyz_file = "$file.xyz";
  $grd_file = "$file.grd";
  $sample_grd_file = "$grd_file.sample";
  $cut_int_file = "$file.int";
  $ps_file = "$file.ps";
  $tiff_file = "$file.tiff";
  print CSH "\n# frame $frame\n";
  print CSH "awk '\$3 != $null {print \$0}' $xyz_file | xyz2grd -I$dx/$dy $R -G$grd_file -N128 -V \n";
  print CSH "grdsample $grd_file -G$sample_grd_file -I0.00833333/0.00833333 -F -V\n";
  print CSH "grdcut $int_file -G$cut_int_file $R -V \n";
  print CSH "grdimage $sample_grd_file $JM $R -I$cut_int_file -C$cpt_file $B -K -V -P > $ps_file\n";
  print CSH "pscoast $JM $R -W4 -Na -Dh -K -O -P -V >> $ps_file \n";
  print CSH "psxy $JM $R -M -W2 -K -O -P -V $fault_file >> $ps_file\n";
  print CSH "psxy $JM $R -Sa0.10 -W1 -G0/0/0 -K -O -P -V <<EOF >> $ps_file\n$elon $elat\nEOF\n";
  print CSH "pstext $JM $R -N -W255/255/255 -O -P -V <<EOF >>$ps_file \n";
  print CSH "-114.8 36.5 12 0 0 RT time = $time s \nEOF\n";
  print CSH "convert $ps_file $tiff_file\n";
}
close(CSH);

