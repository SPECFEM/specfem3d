#!/usr/bin/perl -w

#==========================================================
#
#  plot_geometry.pl
#  Carl Tape
#  18-Nov-2009
#
#  This plots the source-receiver geometry and source time function (STF) for two possible cases
#     1.  PSV wavefield, view from side, three components of STF
#     2.  SH wavefield, top view from surface, one component of STF
#
#  EXAMPLES (from /ADJOINT_TOMO/iterate_adj/SEM2D_iterate/PLOTTING/FIGURES):
#    ../plot_geometry.pl 0 100 events_dat_lonlat.dat recs_lonlat.dat stf_00001 structure_syn.dat 10.0
#    ../plot_geometry.pl 0 200 events_dat_lonlat.dat recs_lonlat.dat stf_00001 structure_syn.dat 10.0
#    ../plot_geometry.pl 0 300 events_dat_lonlat.dat recs_lonlat.dat stf_00001 structure_syn.dat 10.0
#    ../plot_geometry.pl 0 400 events_dat_lonlat.dat recs_lonlat.dat stf_00001 structure_syn.dat 10.0
#    ../plot_geometry.pl 0 500 events_dat_lonlat.dat recs_lonlat.dat stf_00001 structure_syn.dat 10.0
#
#==========================================================

if (@ARGV < 7) {die("Usage: plot_geometry.pl xxx \n");}
($ibody,$irun,$efile_syn,$rfile_syn,$stf_tag,$mfile_syn,$hdur) = @ARGV;

# base directory
$pwd = $ENV{PWD};
$dirplot = `dirname $pwd`; chomp($dirplot);
$basedir = `dirname $dirplot`; chomp($basedir);
#$basedir = "/data1/cig/seismo/3D/ADJOINT_TOMO/iterate_adj/SEM2D_iterate_work";

# directories
$plotdir = "${basedir}/PLOTTING";
$odir = "${basedir}/OUTPUT";
$figdir = "${plotdir}/FIGURES";
if (not -e $figdir) {die("Check if figdir $figdir exist or not\n");}

$strun = sprintf("%4.4i",$irun);
$rdir = "run_${strun}";
$basedir = "${odir}/${rdir}";
$evefile = "$basedir/${efile_syn}";
$recfile = "$basedir/${rfile_syn}";
$file1syn = "$basedir/${mfile_syn}";

# event to choose for getting source time function
$ievent = 5;
$stev  = sprintf("%3.3i",$ievent);
$edir = "${basedir}/event_${stev}";
$stffile1 = "$edir/${stf_tag}_1";
$stffile2 = "$edir/${stf_tag}_2";
$stffile3 = "$edir/${stf_tag}_3";

if (not -f $evefile) {die("Check if evefile $evefile exist or not\n");}
if (not -f $recfile) {die("Check if recfile $recfile exist or not\n");}
if (not -f $file1syn) {die("Check if file1syn $file1syn exist or not\n");}
#if (not -f $stffile1) { die("Check if $stffile1 exist or not\n") }
#if (not -f $stffile2) { die("Check if $stffile2 exist or not\n") }
#if (not -f $stffile3) { die("Check if $stffile3 exist or not\n") }

#$edir  = sprintf("$edir/event_%3.3i",$event);
#print "\n $dir,$edir,$mfile_dat,$mfile_syn,$c0,$per,$ifinite,$iopt \n";

$plabel = "${plotdir}/plot_geometry.pl";

#die("\ntesting\n");

# plotting specifications
$fsize0 = "14";
$fsize1 = "11";
$fsize2 = "9";
$fsize3 = "5";
$fontno = "1";      # 1 or 4
$tick   = "0.2c";

# plot symbols for sources, receivers, and shelf
$size = 15; $hsize = $size/2;
$src = "-Sa${size}p -W0.5p,0/0/0 -G255/255/255";
#$rec = "-Si${size}p -W0.5p -D0/${hsize}p -G0/255/0";
$rec = "-Si${size}p -W0.5p -G0/255/0";

# write plotting scripts

if ($ibody == 1) {
  $wid = 7;
  $dx0 = $wid;
  $dy0 = $dx0*($zran/$xran);
  $J = "-JX${dx0}i/${dy0}i";
  $origin = "-X1.0 -Y6.5";
  $nstf = 3;

  ($xmin,$xmax,$zmin,$zmax) = split(" ",`minmax -C -I2 $file1syn`);
  $xmin = $xmin/1000; $xmax = $xmax/1000;
  $zmin = $zmin/1000; $zmax = $zmax/1000;
  $pinc = 0.0;      # buffer around min/max
  $B = "-B20:.\" \":WESn";
  @stflabs = ("x, h (t)","y, h (t)","z, h (t)");
  $xtick1 = 1; $xtick2 = 1;
  $ytick1 = $xtick1; $ytick2 = $ytick1;
  $Bopt = "WESn";

} else {
  $wid = 6;
  $J = "-JM${wid}i";
  $origin = "-X1.25 -Y3.5";
  $nstf = 1;
  $pinc = 0.;
  ($xmin,$xmax,$zmin,$zmax) = split(" ",`minmax -C $file1syn`);
  @stflabs = ("h (t)","h (t)","h (t)");
  $xtick1 = 1; $xtick2 = 0.25;
  $ytick1 = $xtick1; $ytick2 = $ytick1;
  $Bopt = "WESN";
}
$B = "-Ba${xtick1}f${xtick2}:.\" \":${Bopt}";

# bounds for plotting
$xran = $xmax - $xmin; $zran = $zmax - $zmin;
$xmin = $xmin-$pinc*$xran;  $xmax = $xmax+$pinc*$xran;
$zmin = $zmin-$pinc*$zran;  $zmax = $zmax+$pinc*$zran;
$R = "-R$xmin/$xmax/$zmin/$zmax";

print "\n $R \n";

# plot title
$xtx = $xmin+0.5*($xmax-$xmin);
$ztx = $zmin+1.10*($zmax-$zmin);

$name = "geometry_${strun}";
$psfile  = "${name}.ps";
$jpgfile = "${name}.jpg";

#===============================================
print "\nWriting CSH file...\n";
$cshfile = "plot_geometry.csh";
open(CSH,">$cshfile");
print CSH "gmtset BASEMAP_TYPE plain PAPER_MEDIA letter MEASURE_UNIT inch TICK_LENGTH $tick LABEL_FONT_SIZE $fsize2 ANOT_FONT_SIZE $fsize2 PLOT_DEGREE_FORMAT D HEADER_FONT $fontno ANOT_FONT $fontno LABEL_FONT $fontno HEADER_FONT_SIZE $fsize1\n";
#===============================================

#$origin = "-X1.0 -Y6.5";

# make text files for plotting source and receiver
$just = "CM"; $angle = 0;
if ($ibody==1) {
  print CSH "awk '{print \$1/1000,\$2/1000,$fsize3,$angle,$fontno,\"$just\",\$3}' $recfile > rectext\n";
  print CSH "awk '{print \$1/1000,\$2/1000,$fsize3,$angle,$fontno,\"$just\",\$3}' $evefile > evetext\n";
} else {
  print CSH "awk '{print \$1,\$2,$fsize3,$angle,$fontno,\"$just\",\$3}' $recfile > rectext\n";
  print CSH "awk '{print \$1,\$2,$fsize3,$angle,$fontno,\"$just\",\$3}' $evefile > evetext\n";
}

# phase velocity map
print CSH "psbasemap $J $R $B -K -P -V $origin > $psfile\n"; # START

# plot discontinutities
if ($ibody==1) {
  print CSH "psxy -N $J $R -K -O -V -W0.5p >> $psfile<<EOF\n 0 45 \n 400 45 \nEOF\n";
  print CSH "psxy -N $J $R -K -O -V -W0.5p >> $psfile<<EOF\n 0 65 \n 400 65 \nEOF\n";
  print CSH "psxy -N $J $R -K -O -V -W0.5p >> $psfile<<EOF\n 0 75 \n 400 75 \nEOF\n";
} else {
  print CSH "pscoast $J $R $B -W1p -Na/1p -Dh -K -O -V >> $psfile\n";
}

# plot sources
if ($ibody==1) {
  print CSH "awk '{print \$1/1000,\$2/1000}' $recfile |psxy -N $J $R -K -O -V $rec >> $psfile\n";
} else {
  print CSH "awk '{print \$1,\$2}' $recfile |psxy -N $J $R -K -O -V $rec >> $psfile\n";
}
print CSH "pstext rectext -N $J $R -K -O -V >> $psfile\n";

# plot receivers
if ($ibody==1) {
  print CSH "awk '{print \$1/1000,\$2/1000}' $evefile |psxy -N $J $R -K -O -V $src >> $psfile\n";
} else {
  print CSH "awk '{print \$1,\$2}' $evefile |psxy -N $J $R -K -O -V $src >> $psfile\n";
}
print CSH "pstext evetext -N $J $R -K -O -V >> $psfile\n";

# plot title and GMT header
$Utag = "-U/0/0.25/$plabel";
if($ibody==1) {$ytemp = $dy0 * 1.2;} else {$ytemp = $wid*1.1}
$dX = 0;
$shift = "-Xa-$dX -Ya$ytemp";
$title = "Source receiver geometry ($rdir)";
print CSH "pstext -N $J $R $Utag -K -O -V $shift >>$psfile<<EOF\n $xmin $zmin $fsize0 0 $fontno LM $title \nEOF\n";

#===============================================

# source time function
print "\n Source time function is $stffile1\n";

# determine the max value among the source time functions
for ($k = 1; $k <= $nstf; $k ++ ) {
  $sfile = "${edir}/${stf_tag}_${k}";
  if (not -f $sfile) {die("Check if sfile $sfile exist or not\n");}
  ($tmin,$tmax,$smin,$smax) = split(" ",`minmax -C $sfile`);
  ($ss) = sort{$b <=> $a} (abs($smin),abs($smax));
}
if ($ss==0) {
  $ss = 1;
}       # do not normalize by zero
print "tmin = $tmin\ntmax = $tmax\nss = $ss\n";
$tinc = 20;

#$stf_data = "temp";
$Bs1 = sprintf("-Ba%3.3ff10:\" \":/%3.3f:\"%s\"::.\"Source time function (hdur = %3.1f s)\":WesN",$tinc,0.5,$stflabs[0],$hdur);
$Bs2 = sprintf("-Ba%3.3ff10:\"Time  (s)\":/%3.3f:\"%s\"::.\" \":Wesn",$tinc,0.5,$stflabs[1]);
$Bs3 = sprintf("-Ba%3.3ff10:\"Time  (s)\":/%3.3f:\"%s\"::.\" \":WeSn",$tinc,0.5,$stflabs[2]);
$Rs = "-R$tmin/$tmax/-1/1";
$ywid = 1.25;
$Js = "-JX$wid/$ywid";

$shift1 = "-X0 -Y-2.75";
$shift2 = "-X0 -Y-$ywid";

print CSH "awk '{print \$1,\$2/$ss}' $stffile1 | psxy $Rs $Js $Bs1 -W1p -K -O -V $shift1 >> $psfile\n";
if ($ibody==1) {
  print CSH "awk '{print \$1,\$2/$ss}' $stffile2 | psxy $Rs $Js $Bs2 -W1p -K -O -V $shift2 >> $psfile\n";
  print CSH "awk '{print \$1,\$2/$ss}' $stffile3 | psxy $Rs $Js $Bs3 -W1p -K -O -V $shift2 >> $psfile\n";
}

$plot_title = "Velocity model and source time function";
print CSH "pstext -N $J -R0/1/0/1 -O -V >>$psfile<<EOF\n 10 10 $fsize0 0 $fontno CM $plot_title\nEOF\n"; # FINISH

#-----------------------------
print CSH "convert $psfile $jpgfile\n";
print CSH "echo done with $psfile\n";

close (CSH);
system("csh -f $cshfile");
system("gv $psfile &")

#=================================================================
