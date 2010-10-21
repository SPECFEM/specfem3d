#!/usr/bin/perl -w

use POSIX;

$cmt_file = "CMTSOLUTION";
if (not -f $cmt_file) {die("No such file as $cmt_file");}
(undef,undef,undef,,undef,undef,undef,undef,undef,undef,$elat,$elon) = &get_cmt($cmt_file);

$dx = 0.0145833333332916;
$dy = 0.01197916666667;
@frame = ("04","08","12","16","20","24","28","32","36","40","44","48");
#@frame = ("000800","001600","002400","003200","004000","004800","005600","006400","007200","008000","008800","009600");
$nframe = @frame;
$ncolumn = 3;
$nrow = ceil($nframe/$ncolumn);

$xlen = 7; $ylen = 9;
$drow = $ylen/$nrow; $dcolumn = $xlen/$ncolumn;
$dcolumn1 = $dcolumn * 0.9;
$x0 = 0; $y0 = 0;

$R = "-R-120.3/-114.7/32.2/36.8";
$JM = "-JM${dcolumn1}i";
$ps_file = "movie.ps";
$null = 127.5000;

open(CSH,">plot_movie.csh");
print CSH "gmtset BASEMAP_TYPE plain ANOT_FONT_SIZE 9 HEADER_FONT_SIZE 15\n";
print CSH "makecpt -Cpolar -T0/255/2 -Z -V > grd.cpt\n";
print CSH "psxy -JX1/1 -R0/1/0/1 -K -P -V -X1 -Y1 <<EOF >$ps_file\nEOF\n";
print CSH "psxy -JX1/1 -R0/1/0/1 -X$x0 -Y$y0 -K -O -P -V  <<EOF >>$ps_file\nEOF\n";

column: for ($j=0; $j<$ncolumn; $j++) {
  for ($i=0;$i<$nrow; $i++) {
    if ($i + $j*$nrow + 1 > $nframe) {last column;}
    $frame = $frame[$i+$nrow* $j]; $time = $frame * 200 * 0.009 - 1.0;
    $file = "gmt_movie_0000${frame}.xyz";
    if (not -f $file) {die("No $file\n");}
    print "Processing frame $file...\n";
    $grdfile = "gmt_movie_0000${frame}.grd";
    $xpos = $j * $dcolumn; $ypos = $drow * ($nrow - $i-1);
    $XY = " -X$xpos -Y$ypos "; $NXY = " -X-$xpos -Y-$ypos ";
    $B = " -B2/2wesn ";
    print CSH "awk '\$3 != $null {print \$0}' $file | xyz2grd -I$dx/$dy $R -G$grdfile -N128 -V \n";
    print CSH "grdsample $grdfile -G$grdfile.1 -I0.00833333/0.00833333 -F -V\n";
    print CSH "grdimage $grdfile.1 $JM $R $XY -Itopo.int -Cgrd.cpt $B -K -O -V -P >> $ps_file\n";
    print CSH "pscoast $JM $R -W4 -Na -Dh -K -O -P -V >> $ps_file \n";
    print CSH "psxy $JM $R -M -W2 -K -O -P -V /home/lqy/datalib/jennings.xy >> $ps_file\n";
    print CSH "psxy $JM $R -Sa0.10 -W1 -G0/0/0 -K -O -P -V <<EOF >> $ps_file\n$elon $elat\nEOF\n";
    print CSH "pstext $JM $R -N -W255/255/255 -K -O -P -V <<EOF >>$ps_file \n";
    print CSH "-114.8 36.5 12 0 0 RT time = $time s \nEOF\n";
    print  CSH "psxy -JX1/1 -R0/1/0/1 $NXY -K -O -P -V <<EOF >>$ps_file\nEOF\n";
  }
}
print  CSH "psxy -JX1/1 -R0/1/0/1 -O -P -V <<EOF >>$ps_file\nEOF\n";
close(CSH);

system("csh -fv plot_movie.csh; xv $ps_file");


#######################################
sub get_cmt {
  my ($cmt_file)=@_;
  open(CMT, "$cmt_file") or die("Error opening $cmt_file\n");
  @cmt = <CMT>;
  my(undef,$oyear,$omonth,$oday,$ohr,$omin,$osec1,undef,undef,undef,$mw,undef,@evnm)=split(" ",$cmt[0]);
  my($osec,$omsec)=split(/\./,$osec1); $omsec=$omsec*10;
  my(undef,undef,$evid) = split(" ",$cmt[1]);
  my(undef,undef,$tshift)=split(" ",$cmt[2]);
  my(undef,undef,$hdur)=split(" ",$cmt[3]);
  my(undef,$elat)=split(" ",$cmt[4]);
  my(undef,$elon)=split(" ",$cmt[5]);
  my(undef,$edep)=split(" ",$cmt[6]);
  my(undef,$Mrr)=split(" ",$cmt[7]);
  my(undef,$Mtt)=split(" ",$cmt[8]);
  my(undef,$Mpp)=split(" ",$cmt[9]);
  my(undef,$Mrt)=split(" ",$cmt[10]);
  my(undef,$Mrp)=split(" ",$cmt[11]);
  my(undef,$Mtp)=split(" ",$cmt[12]);
  close(CMT);
  return ($oyear,$omonth,$oday,$ohr,$omin,$osec,$omsec,$tshift,$hdur,$elat,$elon,$edep,$Mrr,$Mtt,$Mpp,$Mrt,$Mrp,$Mtp,$mw,@evnm);
}
