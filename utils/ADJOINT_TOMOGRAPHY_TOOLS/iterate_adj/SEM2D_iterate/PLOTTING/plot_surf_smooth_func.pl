#!/usr/bin/perl -w

#==========================================================
#
#  plot_surf_smooth_func.pl
#  Carl Tape
#  28-Nov-2006
#
# This script plots the unsmoothed and smoothed functions.
# Smoothing is computed by convolving a function with a 2D Gaussian.
# Copied from plot_model.pl on 16-Oct-2005.
#
# EXAMPLE
#   scripts/plot_surf_smooth_func.pl . fun_smooth.dat 10000
#
#   scripts/plot_surf_smooth_func.pl OUTPUT/run_0160 fun_smooth.dat 30000
#   scripts/plot_surf_smooth_func.pl OUTPUT/run_0220 fun_smooth.dat 25000
#   scripts/plot_surf_smooth_func.pl OUTPUT/run_0100 fun_smooth.dat 20000
#   scripts/plot_surf_smooth_func.pl OUTPUT/run_0200 fun_smooth.dat 15000
#   scripts/plot_surf_smooth_func.pl OUTPUT/run_0040 fun_smooth.dat 10000
#   scripts/plot_surf_smooth_func.pl OUTPUT/run_0240 fun_smooth.dat  5000
#
#   scripts/plot_surf_smooth_func.pl OUTPUT/run_0700 fun_smooth.dat 30000
#
#   scripts/plot_surf_smooth_func.pl OUTPUT/run_0460 fun_smooth.dat 60000
#   scripts/plot_surf_smooth_func.pl . fun_smooth.dat 30000
#
#==========================================================

if (@ARGV < 3) {die("Usage: plot_surf_smooth_func.pl basedir xxx \n");}
($dir,$file_smooth,$gamma) = @ARGV;

# plotting specifications
$fsize0 = "14";
$fsize1 = "11";
$fsize2 = "9";
$fsize3 = "5";
$fontno = "1";    # 1 or 4
$tick   = "0.1c";

# plot symbols for sources, receivers, and shelf
if($ifinite==0) {$src = "-W1.0p -Sa0.20 -G255/0/0";} else {$src = "-Sc0.05";}
$rec = "-W0.5p/0/0/0 -St0.10";
$rec0 = "-Sc10p -W0.5p";
$Wshelf = "-W1.0/0/0/0tap";

$plabel = "/home/carltape/sem2d/2d_adjoint/scripts/plot_surf_smooth_func.pl";

#-------------------------
# color

# KEY: scaling for color
$scale_color = 21.0;
$colorbar = "seis";

#-------------------------

# write plotting scripts
$wid = 3.4;
#$J1 = "-JM${wid}i";      # in lat-lon
$J1 = "-JX${wid}i";
$origin = "-X0.65 -Y6.0";

$Bl = "-B50:.\" \":WeSN";
$Br = "-B50:.\" \":wESN";
$Bb = "-B50:.\" \":WESN";

$Dlen = 2.0;
$Dx = $wid/2;
$Dy = -0.35;

$Dscale1 = "-D$Dx/$Dy/$Dlen/0.10h";

$cmax = 3e-6;
$cmax = 5e-7;

$cmax2 = 4/(3.14159*$gamma*$gamma);
#$cmax2 = 1/(2*3.14159*$sigma*$sigma);

print "\n cmax for gaussian is $cmax2\n";

$bs1 = 1;
$bs2 = 1;
$bs4 = 1000;
$Bscale1  = sprintf("-B%2.2e:\" Rough function, K(x, y)\":",$bs1);
$Bscale2  = sprintf("-B%2.2e:\" Smooth function, K'(x, y)\":",$bs1);
$Bscale3  = sprintf("-B%2.2e:\" Residual value, K - K' \":",$bs1);
$BscaleG  = sprintf("-B%2.2e:\" Gaussian filter,  \@~\147\@~ = %3.3f km \":",$bs2,$gamma/1000);
$Bscale4  = sprintf("-B%2.2e:\" Smooth function, g\@+\@~\141\@~\@~\142\@~\@+  \":",$bs4);

#-------------------------

$file = "$dir/$file_smooth";
if (not -f $file)   { die("Check if $file exist or not\n") }

# set bounds for the plotting
$xmin = 0; $xmax = 480; $zmin = 0; $zmax = 480;
$fac = 40;
$xmin = $xmin - $fac; $xmax = $xmax + $fac;
$zmin = $zmin - $fac; $zmax = $zmax + $fac;
$R = "-R$xmin/$xmax/$zmin/$zmax";

# plot title
$xtx = $xmin+0.5*($xmax-$xmin);
$ztx = $zmin+1.10*($zmax-$zmin);

$name    = "smooth";
$psfile  = "$dir/$name.ps";
$jpgfile = "$dir/$name.jpg";

  #===============================================
  print "\nWriting CSH file...\n";
  $cshfile = "plot_smoothed_func.csh";
  open(CSH,">$cshfile");
  print CSH "gmtset BASEMAP_TYPE plain TICK_LENGTH $tick LABEL_FONT_SIZE $fsize2 ANOT_FONT_SIZE $fsize2 PLOT_DEGREE_FORMAT D HEADER_FONT $fontno ANOT_FONT $fontno LABEL_FONT $fontno HEADER_FONT_SIZE $fsize1\n";
  #===============================================

  # make colorpoint file
  $dc = $cmax/10;
  $T1 = "-T-${cmax}/${cmax}/$dc";
  $cptfile1 = "color1.cpt";
  print CSH "makecpt -C$colorbar $T1 -D > $cptfile1\n";
  #print CSH "sed 's/^B.*/B       170     0      0  /' temp1  >  temp2\n";
  #print CSH "sed 's/^F.*/F         0     0    255  /' temp2 > $cptfile1\n";

  # make colorpoint file
  $dc2 = $cmax2/10;
  $T2  = "-T-${cmax2}/${cmax2}/$dc2";
  $cptfile2 = "color2.cpt";
  print CSH "makecpt -C$colorbar $T2 -D > $cptfile2\n";
  #print CSH "sed 's/^B.*/B       170     0      0  /' temp1  >  temp2\n";
  #print CSH "sed 's/^F.*/F         0     0    255  /' temp2 > $cptfile2\n";

  # model for data (optional)
  $plot_title = "Gaussian smoothing of a function";
  $shift = "-Y4.0 -X-0.9";

  # rough function
  print CSH "psbasemap $Bl $J1 $R -P -K -V $origin > $psfile\n";  # START
  print CSH "awk '{print \$3/1000,\$4/1000,\$5}' $file | pscontour $R $J1 -A- -C$cptfile1 -I -P -K -O -V >> $psfile\n";

  #$grdfile = "temp.grd"; $interp = "$I $S";
  #print CSH "nearneighbor $dfile -G$grdfile $R $interp \n";
  #print CSH "grdimage $grdfile -C$cptfile $J -K -O -V >> $psfile\n";   # -T for no interpolation

  #print CSH "pscoast $J1 $R $B1dat -W1p -Na/1p -Dl -P -K -O -V >> $psfile\n";
  print CSH "awk '{print \$2,\$1}' INPUT/oms_shelf |psxy $J1 $R $Wshelf -K -O -P -V >> $psfile\n";
  print CSH "psscale -C$cptfile1 $Dscale1 $Bscale1 -P -K -O -V >> $psfile \n";

  # gaussian function
  $dX = 1.1*$wid; $dY = 0;
  $shift = "-X$dX -Y$dY";
  print CSH "psbasemap $Br $J1 $R -P -K -O -V $shift >> $psfile\n";
  print CSH "awk '{print \$3/1000,\$4/1000,\$6}' $file | pscontour $R $J1 -A- -C$cptfile2 -I -P -K -O -V >> $psfile\n";
  print CSH "awk '{print \$2,\$1}' INPUT/oms_shelf |psxy $J1 $R $Wshelf -K -O -P -V >> $psfile\n";
  print CSH "psscale -C$cptfile2 $Dscale1 $BscaleG -P -K -O -V >> $psfile \n";

  # smooth function
  $dX = -1.1*$wid; $dY = -1.4*$wid;
  $shift = "-X$dX -Y$dY";
  print CSH "psbasemap $Bl $J1 $R -P -K -O -V $shift >> $psfile\n";
  print CSH "awk '{print \$3/1000,\$4/1000,\$7}' $file | pscontour $R $J1 -A- -C$cptfile1 -I -P -O -K -V >> $psfile\n";
  print CSH "awk '{print \$2,\$1}' INPUT/oms_shelf |psxy $J1 $R $Wshelf -K -O -P -V >> $psfile\n";
  print CSH "psscale -C$cptfile1 $Dscale1 $Bscale2 -K -P -O -V >> $psfile \n";

  # residual function
  $dX = 1.1*$wid; $dY = 0;
  $shift = "-X$dX -Y$dY";
  print CSH "psbasemap $Br $J1 $R -P -K -O -V $shift >> $psfile\n";
  print CSH "awk '{print \$3/1000,\$4/1000,\$8}' $file | pscontour $R $J1 -A- -C$cptfile1 -I -P -O -K -V >> $psfile\n";
  print CSH "awk '{print \$2,\$1}' INPUT/oms_shelf |psxy $J1 $R $Wshelf -K -O -P -V >> $psfile\n";
  print CSH "psscale -C$cptfile1 $Dscale1 $Bscale3 -K -P -O -V >> $psfile \n";

  # plot title and GMT header
  $Utag = "-U/0/0.25/$plabel";
  $shift = "-Xa-$dX -Ya8.75";
  print CSH "pstext -N $J1 $R $Utag -O -P -V $shift >>$psfile<<EOF\n $xmin $zmin $fsize0 0 $fontno LM $plot_title\nEOF\n";  # FINISH

#-----------------------------
  print CSH "convert $psfile $jpgfile\n";
  print CSH "echo done with $psfile\n";

  close (CSH);
  system("csh -f $cshfile");
  system("xv $jpgfile &");

#===============================================

#$name    = "smooth2";
#$psfile  = "$dir/$name.ps";
#$jpgfile = "$dir/$name.jpg";

#  #===============================================
#  print "\nWriting CSH file...\n";
#  $cshfile = "plot_smoothed_func.csh";
#  open(CSH,">$cshfile");
#  print CSH "gmtset BASEMAP_TYPE plain TICK_LENGTH $tick LABEL_FONT_SIZE $fsize2 ANOT_FONT_SIZE $fsize2 PLOT_DEGREE_FORMAT D HEADER_FONT $fontno ANOT_FONT $fontno LABEL_FONT $fontno HEADER_FONT_SIZE $fsize1\n";
#  #===============================================

#  # make colorpoint file
#  $cmax3 = 20;
#  $dc = $cmax3/10;
#  $T2 = "-T-${cmax3}/${cmax3}/$dc";
#  $cptfile3 = "color3.cpt";
#  print CSH "makecpt -C$colorbar $T2 > temp1\n";
#  print CSH "sed 's/^B.*/B       170     0      0  /' temp1  >  temp2\n";
#  print CSH "sed 's/^F.*/F         0     0    255  /' temp2 > $cptfile3\n";

#  # model for data (optional)
#  $plot_title = "Two representations of the misfit gradient";
#  $shift = "-Y4.0 -X-0.9";

#  # smooth function
#  $dX = -1.1*$wid; $dY = -1.4*$wid;
#  $shift = "-X$dX -Y$dY";
#  print CSH "psbasemap $Bl $J1 $R -P -K -V $origin > $psfile\n";  # START
#  print CSH "awk '{print \$3/1000,\$4/1000,\$5}' $file | pscontour $R $J1 -A- -C$cptfile1 -I -P -O -K -V >> $psfile\n";
#  print CSH "awk '{print \$2,\$1}' INPUT/oms_shelf |psxy $J1 $R $Wshelf -K -O -P -V >> $psfile\n";
#  print CSH "psscale -C$cptfile1 $Dscale1 $Bscale2 -K -P -O -V >> $psfile \n";

#  # residual function
#  $dX = 1.1*$wid; $dY = 0;
#  $shift = "-X$dX -Y$dY";
#  print CSH "psbasemap $Br $J1 $R -P -K -O -V $shift >> $psfile\n";
#  print CSH "awk '{print \$3/1000,\$4/1000,\$7}' $file | pscontour $R $J1 -A- -C$cptfile3 -I -P -O -K -V >> $psfile\n";
#  print CSH "awk '{print \$2,\$1}' INPUT/oms_shelf |psxy $J1 $R $Wshelf -K -O -P -V >> $psfile\n";
#  print CSH "psscale -C$cptfile3 $Dscale1 $Bscale4 -K -P -O -V >> $psfile \n";

#  # plot title and GMT header
#  $Utag = "-U/0/0.25/$plabel";
#  $shift = "-Xa-$dX -Ya8.75";
#  print CSH "pstext -N $J1 $R $Utag -O -P -V $shift >>$psfile<<EOF\n $xmin $zmin $fsize0 0 $fontno LM $plot_title\nEOF\n";  # FINISH

##-----------------------------
#  print CSH "convert $psfile $jpgfile\n";

#  close (CSH);
#  system("csh -f $cshfile");
#  system("xv $jpgfile &");

#=================================================================
