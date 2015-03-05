#!/usr/bin/perl -w

#==========================================================
#
#  plot_image.pl
#  Carl Tape
#  20-Oct-2006
#
#  This function inputs a lon-lat-val file and plots the image.
#  It is designed to either take a lat-lon plot ($isurf = 1) or a depth section ($isurf = 0).
#
#  Examples:
#    scripts/plot_image.pl 0 OUTPUT_2/run_5013 structure_dat.dat 1/2/3 5.5/8.0/0.50/3 SoCal_1D P_wave_speed socal_P -I1 -S4
#    scripts/plot_image.pl 0 OUTPUT_2/run_5013 structure_dat.dat 1/2/4 3.0/4.5/0.25/3 SoCal_1D S_wave_speed socal_S -I1 -S4
#    scripts/plot_image.pl 0 OUTPUT_2/run_5013 structure_dat.dat 1/2/5 2.4/3.0/0.10/3 SoCal_1D Density socal_d -I1 -S4
#
#    scripts/plot_image.pl 0 OUTPUT_2/run_5000 structure_dat.dat 1/2/5 2.4/3.0/0.10/3 SoCal_1D Density socal_d -I1 -S4
#
#    scripts/plot_image.pl 0 OUTPUT_2/run_5100 structure_dat.dat 1/2/3 5.5/8.0/0.50/3 SoCal_1D P_wave_speed socal_P 1 -I1 -S4
#    scripts/plot_image.pl 0 OUTPUT_2/run_5100 structure_dat.dat 1/2/4 3.0/4.5/0.25/3 SoCal_1D S_wave_speed socal_S 1 -I1 -S4
#    scripts/plot_image.pl 0 OUTPUT_2/run_5100 structure_syn.dat 1/2/3 5.5/8.0/0.50/3 Homogeneous P_wave_speed socal_P_m0 1 -I1 -S4
#    scripts/plot_image.pl 0 OUTPUT_2/run_5100 structure_syn.dat 1/2/4 3.0/4.5/0.25/3 Homogeneous S_wave_speed socal_S_m0 1 -I1 -S4
#    scripts/plot_image.pl 0 OUTPUT_2/run_5101 structure_syn.dat 1/2/3 5.5/8.0/0.50/3 Test_model P_wave_speed socal_P_m0t 1 -I1 -S4
#    scripts/plot_image.pl 0 OUTPUT_2/run_5101 structure_syn.dat 1/2/4 3.0/4.5/0.25/3 Test_model S_wave_speed socal_S_m0t 1 -I1 -S4
#    scripts/plot_image.pl 0 OUTPUT_2/run_5102 structure_syn.dat 1/2/3 5.5/8.0/0.50/3 model_1 P_wave_speed socal_P_m1 1 -I1 -S4
#    scripts/plot_image.pl 0 OUTPUT_2/run_5102 structure_syn.dat 1/2/4 3.0/4.5/0.25/3 model_1 S_wave_speed socal_S_m1 1 -I1 -S4
#
#    scripts/plot_image.pl 0 OUTPUT_2/run_5250 structure_dat.dat 1/2/3 5.7/6.3/0.2/3 SoCal_1D P_wave_speed socal_P 1 -I1 -S4
#    scripts/plot_image.pl 0 OUTPUT_2/run_5100 structure_dat.dat 1/2/4 3.3/3.7/0.2/3 SoCal_1D S_wave_speed socal_S 1 -I1 -S4
#    scripts/plot_image.pl 0 OUTPUT_2/run_5250 structure_syn.dat 1/2/3 5.7/6.3/0.2/3 Homogeneous P_wave_speed socal_P_m0 1 -I1 -S4
#    scripts/plot_image.pl 0 OUTPUT_2/run_5100 structure_syn.dat 1/2/4 3.3/3.7/0.2/3 Homogeneous S_wave_speed socal_S_m0 1 -I1 -S4
#    scripts/plot_image.pl 0 OUTPUT_2/run_5106 structure_syn.dat 1/2/3 5.7/6.3/0.2/3 Homogeneous P_wave_speed socal_P_m3 1 -I1 -S4
#    scripts/plot_image.pl 0 OUTPUT_2/run_5106 structure_syn.dat 1/2/4 3.3/3.7/0.2/3 Homogeneous S_wave_speed socal_S_m3 1 -I1 -S4
#    scripts/plot_image.pl 0 OUTPUT_2/run_5118 structure_syn.dat 1/2/3 5.7/6.3/0.2/3 Homogeneous P_wave_speed socal_P_m9 1 -I1 -S4
#    scripts/plot_image.pl 0 OUTPUT_2/run_5118 structure_syn.dat 1/2/4 3.3/3.7/0.2/3 Homogeneous S_wave_speed socal_S_m9 1 -I1 -S4
#
#    scripts/plot_image.pl 0 OUTPUT_2/run_5200 summed_ker.dat 1/2/3 -5/5/1/-11 SoCal_1D alpha_kernel atest 1 -I1 -S4
#    scripts/plot_image.pl 0 OUTPUT_2/run_5200 summed_ker.dat 1/2/4 -5/5/1/-10 SoCal_1D beta_kernel btest 1 -I1 -S4
#    scripts/plot_image.pl 0 . fun_smooth.dat 1/2/5 -5/5/1/-10 Unsmoothed_event_kernel_beta K_beta unsmoothed 1 -I1 -S4
#    scripts/plot_image.pl 0 . fun_smooth.dat 1/2/7 -5/5/1/-10 Smoothed_event_kernel_beta K_beta smoothed 1 -I1 -S4
#
#    scripts/plot_image.pl 0 OUTPUT_2/run_5100/event_001 kernel_basis 1/2/3 -3/3/1/-12 Kalpha Kalpha Kalpha_5100 1 -I1 -S4
#    scripts/plot_image.pl 0 OUTPUT_2/run_5100/event_001 kernel_basis 1/2/4 -3/3/1/-11 Kbeta Kbeta Kbeta_5100 1 -I1 -S4
#
#    scripts/plot_image.pl 0 . fun_smooth.dat 1/2/5 -1/1/1/-6 Unsmoothed_event_kernel_beta K_beta unsmoothed 1 -I1 -S4
#    scripts/plot_image.pl 0 . fun_smooth.dat 1/2/7 -1/1/1/-6 Smoothed_event_kernel_beta K_beta smoothed 1 -I1 -S4
#
#    scripts/plot_image.pl 0 OUTPUT_2/run_5300 structure_dat.dat 1/2/3 5.5/7.8/0.2/3 SoCal_1D P_wave_speed socal_P 1 -I1 -S4
#    scripts/plot_image.pl 0 OUTPUT_2/run_5300 structure_dat.dat 1/2/4 3.2/4.5/0.2/3 SoCal_1D S_wave_speed socal_S 1 -I1 -S4
#    scripts/plot_image.pl 0 OUTPUT_2/run_5300 structure_syn.dat 1/2/3 5.5/7.8/0.2/3 SoCal_1D P_wave_speed socal_P_m0 1 -I1 -S4
#    scripts/plot_image.pl 0 OUTPUT_2/run_5300 structure_syn.dat 1/2/4 3.2/4.5/0.2/3 SoCal_1D S_wave_speed socal_S_m0 1 -I1 -S4
#
#-------------------------------
#
#    scripts/plot_image.pl 1 OUTPUT_1/run_4086 socal_vel_dat.dat 1/2/4 -20/20/5/-2 Phase_velocity_map Percent_Perturbation -I1m -S4m
#    scripts/plot_image.pl 1 OUTPUT_1/run_0460 summed_ker.dat 1/2/3 -1.5/1.5/0.5/-6 Misfit_kernel Misfit_kernel -I1m -S4m
#    scripts/plot_image.pl 1 OUTPUT_1/run_0460/event_005 kernel_basis 1/2/3 -0.5/0.5/0.25/-7 Event_kernel Event_kernel -I1m -S4m
#
#    scripts/plot_image.pl 1 OUTPUT_1/run_4480/event_005 kernel_basis 1/2/3 -5/5/2/-8 Event_kernel_TT_xcorr Event_kernel Event_kernel_TT_xcorr 1 -I1m -S4m
#    scripts/plot_image.pl 1 OUTPUT_1/run_4512 socal_vel_syn.dat 1/2/4 -10/10/2/-2 Model_16 Percent_Perturbation fun 1 -I0.5m/0.5m -S4m
#    scripts/plot_image.pl 1 OUTPUT_1/run_4512 socal_vel_dat.dat 1/2/4 -10/10/2/-2 Model_16 Percent_Perturbation fun 1 -I0.5m/0.5m -S4m
#
#    scripts/plot_image.pl 1 OUTPUT_1/run_4600/event_005 kernel_basis 1/2/3 -5/5/2/-5 Event_kernel_TT_mtm Event_kernel Event_kernel_TT_mtm 1 -I1m -S4m
#
#    scripts/plot_image.pl 1 OUTPUT_1/run_4300 summed_ker.dat 1/2/3 -1/1/1/-14 Event_kernel_IKER0 Event_kernel Event_kernel_IKER0 1 -I1m -S4m
#    scripts/plot_image.pl 1 OUTPUT_1/run_4350 summed_ker.dat 1/2/3 -5/5/2/-8 Event_kernel_IKER1 Event_kernel Event_kernel_IKER1 1 -I1m -S4m
#    scripts/plot_image.pl 1 OUTPUT_1/run_4400 summed_ker.dat 1/2/3 -5/5/2/-10 Event_kernel_IKER2 Event_kernel Event_kernel_IKER2 1 -I1m -S4m
#    scripts/plot_image.pl 1 OUTPUT_1/run_4450 summed_ker.dat 1/2/3 -5/5/2/-5 Event_kernel_IKER3 Event_kernel Event_kernel_IKER3 1 -I1m -S4m
#    scripts/plot_image.pl 1 OUTPUT_1/run_4500 summed_ker.dat 1/2/3 -5/5/2/-7 Event_kernel_IKER4 Event_kernel Event_kernel_IKER4 1 -I1m -S4m
#
#    scripts/plot_image.pl 1 OUTPUT_1/run_4600/event_019 kernel_basis 1/2/3 -1/1/0.5/-9 BD_kernel_IKER6 BD_kernel BD_kernel_IKER6 1 -I1m -S4m
#
#    scripts/plot_image.pl 1 OUTPUT_1/run_0460 fun_smooth.dat 1/2/5 -1/1/0.5/-6 Unsmoothed_event_kernel kernel test0 1 -I1m -S4m
#    scripts/plot_image.pl 1 OUTPUT_1/run_0460 fun_smooth.dat 1/2/6 -3/3/1.0/-10 Gaussian_smoothing_function amplitude test1 1 -I1m -S4m
#    scripts/plot_image.pl 1 OUTPUT_1/run_0460 fun_smooth.dat 1/2/7 -1/1/0.5/-6 Smoothed_event_kernel kernel test2  1 -I1m -S4m
#    scripts/plot_image.pl 1 OUTPUT_1/run_0460 fun_smooth.dat 1/2/8 -1/1/0.5/-6 Residual residual test3 1 -I1m -S4m
#
#    scripts/plot_image.pl 1 OUTPUT_1/run_4950/event_001 kernel_basis 1/2/3 -5/5/1/-8 test test test_4950 1 -I1m -S4m
#    scripts/plot_image.pl 1 OUTPUT_1/run_4900/event_001 kernel_basis 1/2/3 -5/5/1/-8 test test test_4900 1 -I1 -S4
#
#    scripts/plot_image.pl 1 OUTPUT_1/run_5102 structure_syn.dat 1/2/4 -10/10/5/-2 Phase_velocity_map Percent_Perturbation test 1 -I1m -S4m
#    scripts/plot_image.pl 1 OUTPUT_1/run_5101 structure_syn.dat 1/2/4 -10/10/5/-2 Phase_velocity_map Percent_Perturbation test_model 1 -I1m -S4m
#
#    scripts/plot_image.pl 1 OUTPUT_2/run_5700 structure_dat.dat 1/2/4 3/4/0.5/3 Target_model S_wave_speed test 1 -I1m -S4m
#    scripts/plot_image.pl 1 OUTPUT_2/run_5751 structure_syn.dat 1/2/4 3/4/0.5/3 Reference_model S_wave_speed test 1 -I1m -S4m
#    scripts/plot_image.pl 1 OUTPUT_2/run_5752 structure_syn.dat 1/2/4 3/4/0.5/3 Reference_model S_wave_speed test 1 -I1m -S4m
#    scripts/plot_image.pl 1 OUTPUT_2/run_5703 structure_syn.dat 1/2/4 3/4/0.5/3 Reference_model S_wave_speed test 1 -I1m -S4m
#
#    scripts/plot_image.pl 1 OUTPUT_2/run_5800 structure_dat.dat 1/2/4 3/4/0.5/3 Target_model S_wave_speed test 1 -I1m -S4m
#    scripts/plot_image.pl 1 OUTPUT_2/run_5801 structure_syn.dat 1/2/4 3/4/0.5/3 Reference_model S_wave_speed test 1 -I1m -S4m
#    scripts/plot_image.pl 1 OUTPUT_2/run_5806 structure_syn.dat 1/2/4 3/4/0.5/3 Reference_model S_wave_speed test 1 -I1m -S4m
#    scripts/plot_image.pl 1 OUTPUT_2/run_5803 structure_syn.dat 1/2/4 3/4/0.5/3 Reference_model S_wave_speed test 1 -I1m -S4m
#
#    scripts/plot_image.pl 1 OUTPUT_2/run_5750 summed_ker.dat 1/2/4 -1/1/1/-6 Event_kernel_IKER0 Event_kernel Event_kernel_IKER0 1 -I1m -S4m
#
#==========================================================

if (@ARGV < 3) {die("Usage: plot_image.pl basedir xxx \n");}
($isurf,$idir,$infile,$dcols,$colors,$title,$scale_label,$filename,$iportrait,$I,$S) = @ARGV;

#print "\n $isurf\n $idir\n $infile\n $dcols\n $colors\n $title\n $scale_label\n $filename\n $I\n $S\n"; die("testing");
#print "\n $filename\n"; die("testing");

# color specifications
@temp  = split("/",$colors);
$cmin  = $temp[0];
$cmax  = $temp[1];
$ctick = $temp[2];
$cpwr  = $temp[3];
$cnorm = "1e$cpwr";
print "\n $cmin $cmax $cpwr $cnorm \n";

@icol = split("/",$dcols);
@plot_title = split("_",$title);
@scale_title = split("_",$scale_label);

$datafile = "$idir/$infile";
if (not -f $datafile) { die("Check if $datafile exist or not\n") }
print "\n data file is $datafile\n";

$nwrd = @plot_title;
$plot_title[$nwrd] = " ($datafile)";

$plabel = "/home/denali2/carltape/wave2d/2d_adjoint/scripts/plot_image.pl";

#die("\ntesting\n");

# plotting specifications
$fsize0 = "14";
$fsize1 = "11";
$fsize2 = "9";
$fsize3 = "5";
$fontno = "1";    # 1 or 4
$tick   = "0.1c";

# plot symbols for sources, receivers, and shelf
$src = "-Sc0.05";
$rec = "-W0.5p/0/0/0 -St0.10";
$rec0 = "-Sc10p -W0.5p";
$Wshelf = "-W1.0/0/0/0tap";

$coastinfo = "-W1p -N1/1p -N2/1p -Df -A100";   # why doesn't the border pen-width work?

#-------------------------
# color

# KEY: scaling for color
$colorinc = 21.0;
$colorbar = "seis";

#-------------------------

# write plotting scripts
if($isurf == 1) {$btick1 = 1; $btick2 = 0.25}
if($isurf == 0) {$btick1 = 40; $btick2 = 10}
$B = "-Ba${btick1}f${btick2}:.\" \":WESN";

#$Bscale  = sprintf("-B%2.2e:\" Percent Perturbation from %3.3f  km s\@+-1\@+ \": -E10p",$ptick,$c0/1000);
$Bscale  = sprintf("-B%2.2e:\" @scale_title  ( 10\@+%2.2i\@+ )\": -E10p",$ctick,$cpwr);

#-------------------------

# bounds for the plotting
#($xmin,$xmax,$zmin,$zmax,$smin,$smax,$tmin,$tmax) = split(" ",`minmax -C $datafile`);
#($tt) = sort{$b <=> $a} ($tt,abs($tmin),abs($tmax));
($xmin,$xmax,$zmin,$zmax) = split(" ",`minmax -C $datafile`);
if($isurf==0) {
  $xmin = $xmin/1000; $xmax = $xmax/1000;
  $zmin = $zmin/1000; $zmax = $zmax/1000;
}
if($isurf==1){$pinc = 0.05}
if($isurf==0){$pinc = 0.0}
$xran = $xmax - $xmin; $zran = $zmax - $zmin;
$xmin = $xmin-$pinc*$xran;  $xmax = $xmax+$pinc*$xran;
$zmin = $zmin-$pinc*$zran;  $zmax = $zmax+$pinc*$zran;
$R = "-R$xmin/$xmax/$zmin/$zmax";
print "\n bounds are $R \n";

# projection
$xwid = 6.0;
$ywid = $xwid*($zran/$xran);
if($isurf == 1) {$J = "-JM${xwid}i"; $ywid = $xwid}
if($isurf == 0) {$J = "-JX${xwid}i/${ywid}i"}
$origin = "-X1.0i -Y3.0i";
print "\n projection is $J \n";

# page orientation
if($iportrait == 1){$orient = "-P"} else {$orient = " "}

$Dlen = 0.7*$xwid;
$Dx = $xwid/2;
$Dy = -0.35;
$Dscale = "-D$Dx/$Dy/$Dlen/0.15h";
print "\n color bar is $Dscale \n $Bscale \n";

#die("testing");

# plot title
#$xtx = $xmin+0.5*($xmax-$xmin);
#$ztx = $zmin+1.10*($zmax-$zmin);

# plot title (see plot_body_geometry.pl ?)
$R_title = "-R0/1/0/1";
$J_title = "-JX${xwid}i/${ywid}i";
$x_title = 0.5;
$y_title = 1.5;
$p_title = 1.75 * $ywid;

#$filename    = "image";
$psfile  = "$filename.ps";
$jpgfile = "$filename.jpg";

  #===============================================
  print "\nWriting CSH file...\n";
  $cshfile = "plot_image.csh";
  open(CSH,">$cshfile");
  print CSH "gmtset BASEMAP_TYPE plain MEASURE_UNIT inch COLOR_NAN 255/255/255 TICK_LENGTH $tick LABEL_FONT_SIZE $fsize2 ANOT_FONT_SIZE $fsize2 PLOT_DEGREE_FORMAT D HEADER_FONT $fontno ANOT_FONT $fontno LABEL_FONT $fontno HEADER_FONT_SIZE $fsize1\n";
  #===============================================

  # make colorpoint file
#  $T1 = "-T3/4/0.1";
#  $pmax = 10;
  $cran = $cmax - $cmin;
  $dc = $cran/$colorinc;
  $T1 = "-T$cmin/$cmax/$dc";
  $cptfile = "color0.cpt";
  print CSH "makecpt -C$colorbar -D -Z $T1 > $cptfile\n";
  print "\n $T1\n";

  #------------------------

  print CSH "psbasemap $J $R $B -K $orient -V $origin > $psfile\n";  # START

  # generate datafile
  $dfile = dtemp;
  if($isurf == 1) {
     print CSH "awk '{print \$($icol[0]),\$($icol[1]),\$($icol[2])/$cnorm }' $datafile > $dfile\n";
  } else {
     print CSH "awk '{print \$($icol[0])/1000,\$($icol[1])/1000,\$($icol[2])/$cnorm }' $datafile > $dfile\n";
  }


  #print CSH "awk '{print \$($icol[0]),\$($icol[1]),\$($icol[2])/$cnorm }' $datafile | pscontour $R $J -A- -C$cptfile -I -K -O -V >> $psfile\n";
  #print CSH "pscontour $dfile $R $J -A- -C$cptfile -I -K -O -V >> $psfile\n";

  $grdfile = "temp.grd"; $interp = "$I $S";
  print CSH "nearneighbor $dfile -G$grdfile $R $interp\n";
  print CSH "grdimage $grdfile -C$cptfile $J -K -O -V >> $psfile\n";   # -T for no interpolation

  # coast-lines
  if($isurf==1) {print CSH "pscoast $J $R $B $coastinfo -K -O -V >> $psfile\n"}

  print CSH "psscale -C$cptfile $Dscale $Bscale -K -O -V >> $psfile\n";

  # plot the element boundaries
  if(0==1) {
     $elementfile = "$idir/elements.dat";
     if (not -f $elementfile) { die("Check if $elementfile exist or not\n") }
     print CSH "awk '{print \$4,\$6}' $elementfile |psxy -N -Sc1p -G0/0/0 $J $R -K -O -V >> $psfile\n";
     print CSH "awk '{print \$5,\$6}' $elementfile |psxy -N -Sc1p -G0/0/0 $J $R -K -O -V >> $psfile\n";
     print CSH "awk '{print \$5,\$7}' $elementfile |psxy -N -Sc1p -G0/0/0 $J $R -K -O -V >> $psfile\n";
     print CSH "awk '{print \$4,\$7}' $elementfile |psxy -N -Sc1p -G0/0/0 $J $R -K -O -V >> $psfile\n";
  }

  # plot receivers with numbered label
  #print CSH "awk '\$1 == \"R\" {print \$2,\$3}' $rec_file |psxy -N $J $R -K -O -V $rec0 >> $psfile\n";
  #$rec_file2 = text_rec; $angle = 0; $just = "CM";
  #print CSH "awk '\$1 == \"R\" {print \$2,\$3,$fsize3,$angle,$fontno,\"$just\",\$4}' $rec_file > $rec_file2\n";
  #print CSH "pstext $rec_file2 -N $J $R -K -O -V >> $psfile\n";
  #print CSH "awk '\$1 == \"S\" {print \$2,\$3}' $rec_file |psxy -N $J $R -K -O -V $src >> $psfile\n";

  # plot title and GMT header
  $Utag = "-U/0/${p_title}/$plabel";
  #$ytemp = $ywid * 1.5;
  #$shift = "-Xa-$dX -Ya$ytemp";
  print CSH "pstext -N $J_title $R_title $Utag -O -V >>$psfile<<EOF\n $x_title $y_title $fsize0 0 $fontno CM @plot_title\nEOF\n";  # FINISH

#-----------------------------
  if($iportrait == 0){print CSH "convert $psfile -rotate 90 $jpgfile\n"}
  if($iportrait == 1){print CSH "convert $psfile $jpgfile\n"}
  print CSH "echo done with $psfile\n";

  close (CSH);
  system("csh -f $cshfile");
  system("xv $jpgfile &")
  #if($iopt <= 2) {system("xv $jpgfile &")}

#=================================================================
