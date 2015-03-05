#!/usr/bin/perl -w

#==========================================================
#
#  plot_body_compare_seis.pl
#  Carl Tape
#  07-Aug-2006
#
#
#  EXAMPLES:
#    scripts/plot_body_compare_seis.pl OUTPUT_body/run_5011 OUTPUT_body/run_5011_5012 events_xy.dat recs_xz.dat structure_syn.dat 3
#    scripts/plot_body_compare_seis.pl OUTPUT_body/run_5013 OUTPUT_body/run_5013_5017 events_xz.dat recs_xz.dat structure_syn.dat 3
#    scripts/plot_body_compare_seis.pl OUTPUT_body/run_5013 OUTPUT_body/run_5013_5017b events_xz.dat recs_xz.dat structure_syn.dat 3
#
#
#==========================================================

if (@ARGV < 3) {die("Usage: plot_body_compare_seis.pl basedir iopt xxx \n")}
($dir1,$dir,$efile_syn,$rfile_syn,$mfile_syn,$ncomp) = @ARGV;

$strun = sprintf("%4.4i",$irun);
$stev  = sprintf("%3.3i",$ievent);

$evefile = "$dir1/${efile_syn}";
$recfile = "$dir1/${rfile_syn}";
$mfile = "$dir1/${mfile_syn}";

if (not -f $evefile)  { die("Check if $evefile exist or not\n") }
if (not -f $recfile)  { die("Check if $recfile exist or not\n") }
if (not -f $mfile) { die("Check if $mfile exist or not\n") }

#$edir  = sprintf("$dir/event_%3.3i",$event);
#print "\n $dir,$edir,$mfile_dat,$mfile_syn,$c0,$per,$ifinite,$iopt \n";

$plabel = "/home/denali2/carltape/wave2d/2d_adjoint/scripts/plot_body_compare_seis.pl";

#die("\ntesting\n");

# plotting specifications
$fsize0 = "14";
$fsize1 = "11";
$fsize2 = "9";
$fsize3 = "5";
$fontno = "1";    # 1 or 4
$tick   = "0.2c";

# plot symbols for sources, receivers, and shelf
$size = 15; $hsize = $size/2;
$src = "-Sa${size}p -W0.5p,0/0/0 -G255/255/255";
#$rec = "-Si${size}p -W0.5p -D0/${hsize}p -G0/255/0";
$rec = "-Si${size}p -W0.5p -G255/255/255";

# bounds for plotting (we could also use the file global_mesh.dat)
($xmin,$xmax,$zmin,$zmax) = split(" ",`minmax -C $mfile`);
$xmin = $xmin/1000; $xmax = $xmax/1000;
$zmin = $zmin/1000; $zmax = $zmax/1000;
$pinc = 0.0;  # buffer around min/max
$xran = $xmax - $xmin; $zran = $zmax - $zmin;
$xmin = $xmin-$pinc*$xran;  $xmax = $xmax+$pinc*$xran;
$zmin = $zmin-$pinc*$zran;  $zmax = $zmax+$pinc*$zran;
$R = "-R$xmin/$xmax/$zmin/$zmax";

print "\n $R \n";

# write plotting scripts
  $wid = 7;
  $dx0 = $wid;                       # width of each subplot
  $dy0 = $dx0*($zran/$xran);
  $J = "-JX${dx0}i/${dy0}i";

$B = "-B20:.\" \":WESn";

# plot title
$xtx = $xmin+0.5*($xmax-$xmin);
$ztx = $zmin+1.10*($zmax-$zmin);

# get sources
if (not -f $evefile) {die("Check if $evefile exist or not\n")}
open(IN,$evefile); @srclines = <IN>; $nsrc = @srclines;

# get receivers
open(IN,$recfile); @reclines = <IN>; $nrec = @reclines;

#$nsrc = 1; $nrec = 1; $ncomp = 1;

  #===============================================
  print "\nWriting CSH file...\n";
  $cshfile = "plot_body_compare_seis.csh";
  open(CSH,">$cshfile");
  print CSH "gmtset BASEMAP_TYPE plain MEASURE_UNIT inch TICK_LENGTH $tick LABEL_FONT_SIZE $fsize2 ANOT_FONT_SIZE $fsize2 PLOT_DEGREE_FORMAT D HEADER_FONT $fontno ANOT_FONT $fontno LABEL_FONT $fontno HEADER_FONT_SIZE $fsize1\n";
  #===============================================

# loop over components
for ($icomp = 1; $icomp <= $ncomp; $icomp++) {

print "\n component is $icomp \n";

for ($isrc = 1; $isrc <= $nsrc; $isrc++) {

for ($irec = 1; $irec <= $nrec; $irec++) {

($xsrc,$zsrc,$junk1,$namesrc) = split(" ",$srclines[$isrc-1]);
print "\n source $namesrc is at (x, z) = ($xsrc, $zsrc) \n";

($xrec,$zrec,$namerec) = split(" ",$reclines[$irec-1]);
print "\n rec $namerec is at (x, z) = ($xrec, $zrec) \n";

$ftag = sprintf("%3.3i_%3.3i_%1.1i",$isrc,$irec,$icomp);

$name = "$dir/seis_$ftag";
$psfile  = "$name.ps";
$jpgfile = "$name.jpg";

  #$plot_title = "Model for synthetics";
  $origin = "-X1.0 -Y8.0";

  # phase velocity map
  print CSH "psbasemap $J $R $B -K -P -V $origin > $psfile\n";  # START

  # plot receivers with numbered label
  #print CSH "awk '\$1 == \"R\" {print \$2,\$3}' $rec_file |psxy -N $J $R -K -O -P -V $rec0 >> $psfile\n";
  #$rec_file2 = text_rec; $angle = 0; $just = "CM";
  #print CSH "awk '\$1 == \"R\" {print \$2,\$3,$fsize3,$angle,$fontno,\"$just\",\$4}' $rec_file > $rec_file2\n";

  $xs = $xsrc/1000; $zs = $zsrc/1000;
  $xr = $xrec/1000; $zr = $zrec/1000;
  print CSH "psxy -N $J $R -K -O -P -V $src >> $psfile<<EOF\n $xs $zs \nEOF\n";
  print CSH "psxy -N $J $R -K -O -P -V $rec >> $psfile<<EOF\n $xr $zr \nEOF\n";

  # plot discontinutities
  print CSH "psxy -N $J $R -K -O -P -V -W0.5p >> $psfile<<EOF\n 0 45 \n 400 45 \nEOF\n";
  print CSH "psxy -N $J $R -K -O -P -V -W0.5p >> $psfile<<EOF\n 0 65 \n 400 65 \nEOF\n";
  print CSH "psxy -N $J $R -K -O -P -V -W0.5p >> $psfile<<EOF\n 0 75 \n 400 75 \nEOF\n";

  $just = "CM"; $angle = 0;
  print CSH "pstext -N $J $R -K -O -P -V >> $psfile<<EOF\n $xs $zs $fsize3 $angle $fontno $just $namesrc\nEOF\n";
  print CSH "pstext -N $J $R -K -O -P -V >> $psfile<<EOF\n $xr $zr $fsize3 $angle $fontno $just $namerec\nEOF\n";

  # plot title and GMT header
  $Utag = "-U/0/0.25/$plabel";
  $ytemp = $dy0 * 1.3;
  $shift = "-Xa-$dX -Ya$ytemp";
  $title = "$dir/seis_compare_${ftag}.dat";
  print CSH "pstext -N $J $R $Utag -K -O -P -V $shift >>$psfile<<EOF\n $xmin $zmin $fsize0 0 $fontno LM $title \nEOF\n";

  #===============================================

  # file with seismograms
  $seisfile = "$dir/seis_compare_${ftag}.dat";
  if (not -f $seisfile)  { die("Check if $seisfile exist or not\n") }

  # source time function
  print "\n seismogram file is $seisfile\n";

  # determine the max value among the source time functions
  ($tmin,$tmax,$junk1,$junk2,$junk3,$junk4,$rmin,$rmax) = split(" ",`minmax -C $seisfile`);
  ($ss) = sort{$b <=> $a} (abs($rmin),abs($rmax));
  $ss = 0.25;   # fixed residual amplitude

  print "smax = $ss \n tmin = $tmin tmax = $tmax\n";
  $ttick1 = 10;
  $ttick2 = $ttick1/2;
  #$ytick = (1.1*$ss)/4;
  $rtick1 = 0.10;
  $rtick2 = 0.05;

  #@titles = ("s0 (NGLL = 7)","s0 - s1","s0 - s2","s0 - s3");
  @titles = ("s1 (NGLL = 5)","s1 - s2","s1 - s3");

  #$stf_data = "temp";
  $Bs0 = sprintf("-Ba%3.3ff%3.3fg%3.3f:\" \":/%3.3f:\"$titles[0]\"::.\" \":WesN",$ttick1,$ttick2,$ttick2,0.5);
  $Bs1 = sprintf("-Ba{%3.3f}f{%3.3f}g{%3.3f}:\" \":/a{%3.3f}f{%3.3f}g{%3.3f}:\"$titles[1]\"::.\" \":Wesn",
         $ttick1,$ttick2,$ttick2,$rtick1,$rtick2,$rtick2);
  $Bs2 = sprintf("-Ba{%3.3f}f{%3.3f}g{%3.3f}:\" \":/a{%3.3f}f{%3.3f}g{%3.3f}:\"$titles[2]\"::.\" \":Wesn",
         $ttick1,$ttick2,$ttick2,$rtick1,$rtick2,$rtick2);
  $Bs3 = sprintf("-Ba{%3.3f}f{%3.3f}g{%3.3f}:\"Time  (s)\":/a{%3.3f}f{%3.3f}g{%3.3f}:\"$titles[3]\"::.\" \":WeSn",
         $ttick1,$ttick2,$ttick2,$rtick1,$rtick2,$rtick2);
  $Rs1 = "-R$tmin/$tmax/-1/1";
  $Rs2 = "-R$tmin/$tmax/-$ss/$ss";
  $ywid = 1.5;
  $Js = "-JX$wid/$ywid";

  #print "\n $Bs1 $Bs2 $Rs $Js \n";

  $shift1 = "-X0 -Y-2.5";
  $shift2 = "-X0 -Y-$ywid";
  $ytxt = 1.5; $txtinfo = "-C3p -W255o";

if(0==1) {

  # plot seismograms (normalized)
  print CSH "awk '{print \$1,\$2}' $seisfile |psxy $Rs1 $Js $Bs0 -W0.5p -P -K -O -V $shift1 >> $psfile\n";
  #print CSH "awk '{print \$1,\$3}' $seisfile |psxy $Rs1 $Js $Bs0 -W0.5p,-- -P -K -O -V >> $psfile\n";

  @titles = ("s0 (NGLL=7, 32 elements)",
             "s0 - s1a (NGLL=6, 32 el), s0 - s1b (NGLL=5, 32 el)",
             "s0 - s2 (NGLL=5, 33 el, NOT honoring mesh)",
             "s0 - s3 (NGLL=5, 34 el, NOT honoring mesh, but GLL pts equidistant from discontinuities)");
  print CSH "pstext -N -JX1 -R0/1/0/1 -K -O -P -V $txtinfo >>$psfile<<EOF\n $wid $ytxt $fsize2 0 $fontno RM $titles[0]\nEOF\n";

  # plot residuals
  print CSH "awk '{print \$1,\$3}' $seisfile |psxy $Rs2 $Js $Bs1 -W0.5p,0/0/255,-- -P -K -O -V $shift2 >> $psfile\n";
  print CSH "awk '{print \$1,\$4}' $seisfile |psxy $Rs2 $Js $Bs1 -W0.5p,0/0/255 -P -K -O -V >> $psfile\n";
  print CSH "pstext -N -JX1 -R0/1/0/1 -K -O -P -V $txtinfo >>$psfile<<EOF\n $wid $ytxt $fsize2 0 $fontno RM $titles[1]\nEOF\n";

  print CSH "awk '{print \$1,\$5}' $seisfile |psxy $Rs2 $Js $Bs2 -W0.5p,255/0/0 -P -K -O -V $shift2 >> $psfile\n";
  print CSH "pstext -N -JX1 -R0/1/0/1 -K -O -P -V $txtinfo >>$psfile<<EOF\n $wid $ytxt $fsize2 0 $fontno RM $titles[2]\nEOF\n";

  print CSH "awk '{print \$1,\$6}' $seisfile |psxy $Rs2 $Js $Bs3 -W0.5p,0/255/0 -P -K -O -V $shift2 >> $psfile\n";
  print CSH "pstext -N -JX1 -R0/1/0/1 -K -O -P -V $txtinfo >>$psfile<<EOF\n $wid $ytxt $fsize2 0 $fontno RM $titles[3]\nEOF\n";

} else {

  # plot seismograms (normalized)
  print CSH "awk '{print \$1,\$2}' $seisfile |psxy $Rs1 $Js $Bs0 -W0.5p,0/0/255 -P -K -O -V $shift1 >> $psfile\n";
  print CSH "awk '{print \$1,\$3}' $seisfile |psxy $Rs1 $Js $Bs0 -W0.5p,255/0/0,-- -P -K -O -V >> $psfile\n";

  @titles = ("s1 (NGLL=5, honoring mesh, 32 elements)",
             "s1 - s2 (NGLL=5, NOT honoring mesh, 33 el)",
             "s1 - s3 (NGLL=5, NOT honoring mesh, but GLL pts equidistant from discontinuities, 34 el)");
  print CSH "pstext -N -JX1 -R0/1/0/1 -K -O -P -V $txtinfo >>$psfile<<EOF\n $wid $ytxt $fsize2 0 $fontno RM $titles[0]\nEOF\n";

  # plot residuals
  print CSH "awk '{print \$1,\$2-\$3}' $seisfile |psxy $Rs2 $Js $Bs1 -W0.5p,255/0/0 -P -K -O -V $shift2 >> $psfile\n";
  print CSH "pstext -N -JX1 -R0/1/0/1 -K -O -P -V $txtinfo >>$psfile<<EOF\n $wid $ytxt $fsize2 0 $fontno RM $titles[1]\nEOF\n";

  print CSH "awk '{print \$1,\$2-\$4}' $seisfile |psxy $Rs2 $Js $Bs2 -W0.5p,0/255/0 -P -K -O -V $shift2 >> $psfile\n";
  print CSH "pstext -N -JX1 -R0/1/0/1 -K -O -P -V $txtinfo >>$psfile<<EOF\n $wid $ytxt $fsize2 0 $fontno RM $titles[2]\nEOF\n";

}

#-----------------------------
  print CSH "pstext -N $J -R0/1/0/1 -O -P -V >>$psfile<<EOF\n 10 10 $fsize0 0 $fontno CM junk\nEOF\n";  # FINISH

  print CSH "convert $psfile $jpgfile\n";
  print CSH "echo done with $psfile\n";

}  # comp
}  # src
}  # rec

  close (CSH);
  system("csh -f $cshfile");
  system("xv $jpgfile &")

#=================================================================
