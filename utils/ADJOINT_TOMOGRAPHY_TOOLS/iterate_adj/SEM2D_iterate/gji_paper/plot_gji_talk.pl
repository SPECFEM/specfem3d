#!/usr/bin/perl -w
#==========================================================
#
#  plot_gji_talk.pl
#  Carl Tape
#  01-March-2006
#
#  ISTRUCTURE:
#    1. checker, Nfac=3
#    2. checker, Nfac=2
#    3. checker, Nfac=1
#    4. rayleigh wave (T=20s, smoothed with sig=10km)
#
#    plot_gji_talk.pl 25 60000 -6/3.0/0/80/8 1 460 0 1 1 1 0   # big checker
#    plot_gji_talk.pl 25 30000 -6/3.0/0/20/8 0 280 0 0 1 4 0   # rayleigh
#
#    plot_gji_talk.pl 25 60000 -9/4.0/0/80/8 1 460 0 1 1 1 0   # ifig03 (princeton kernels)
#    plot_gji_talk.pl 25 60000 -9/4.0/0/80/8 1 460 0 1 1 1 0   # ifig03p (kernel, ray)
#
#    plot_gji_talk.pl  1 30000 -7/1.0/0/80/8 1 000 0 1 0 1 0   # ifig04,05 irun_000 (132 recs)
#    plot_gji_talk.pl  1 30000 -8/4.0/0/80/8 1 001 0 1 0 1 0   # ifig04,05 irun_001 (single src-rec)
#    plot_gji_talk.pl  1 30000 -7/0.5/0/80/8 5 002 0 1 0 1 0   # ifig05 irun_002 (single src-rec) -- no data
#
#    plot_gji_talk.pl   # ifig05p (seismograms)
#    plot_gji_talk.pl   # ifig06 (hessian)
#
#    plot_gji_talk.pl 25 30000 -6/3.0/0/20/8 0 280 0 0 1 4 0   # ifig07 (hessian and gradient for kernels)
#
#
#==========================================================

#if (@ARGV < 7) {die("Usage: plot_kernels.pl out_dir_pre event_dir_pre modelfile_dat modelfile_syn kernelfile kpwr/kmax/oprw/omax IKER irun0 iter model \n");}
#($odir,$edir,$mfile_dat,$mfile_syn,$kerfile,$colors,$iker,$irun0,$iter,$ichecker) = @ARGV;
($nevent,$gamma,$colors,$iker,$irun0,$iter,$ichecker,$ipoly,$istructure,$ibanana) = @ARGV;

if($ibanana==1) {$odir      = "../../../2d_adjoint_banana/OUTPUT_banana/run_"}
else            {$odir      = "../../OUTPUT/run_"}
$edir      = "event_";
$mfile_dat = "socal_vel_dat.dat";
$mfile_syn = "socal_vel_syn.dat";
$kerfile   = "kernel_basis";

$cshfile = "plot_gji_talk.csh";

if($istructure==1){$Nfac=3}
if($istructure==2){$Nfac=2}
if($istructure==3){$Nfac=1}
$sNfac = sprintf("%1i",$Nfac);

# boolean commands for plotting
$icolor = 1;   # ccc
$ixv    = 1;
$ieps   = 0;

$ifig01  = 0;   # CG sequence 1
$ifig02  = 0;   # CG sequence 2
$ifig03  = 0;   # banana-doughnut kernels
$ifig03p = 0;   # banana-doughnut kernels
$ifig04  = 0;   # mod, mod, ker
$ifig04p = 0;   #
$ifig05  = 0;   # 4x1 forward wavefield
$ifig05p = 0;   # seismograms

$ifig06  = 0;   # Hessian
$ifig07  = 1;   # Hessian and gradient

$ifig10  = 0;   # 3 x 3 CG summary

$irun = $irun0 + $iter;
@mods = ("0","0t","1","1t","2","2t","3","3t","4","4t","5","5t","6","6t","7","7t","8","8t");
$mod = $mods[$iter];
$smod = "m\@+$mod\@+";

# label for the type of map used to generate the data
$smap = sprintf("m%2.2i",$istructure);

# wave2d run number
$strun0 = sprintf("%4.4i",$irun0);
$strun = sprintf("%4.4i",$irun);
$dir0 = "${odir}$strun0";
$dir = "${odir}$strun";

#print "\n $nevent,$gamma,$colors,$iker,$irun0,$iter,$ichecker,$ipoly,$istructure,$ibanana -- $dir0, $dir\n"; die("testing");

$stgam  = sprintf("\@~\147\@~ = %.1f km",$gamma/1000);
$stgam2 = sprintf("%3.3i",$gamma/1000);

print "\n $dir, $edir, $mfile_dat, $mfile_syn, $kerfile";
print "\n $colors, $iker, $mod, $stgam2 \n";
#die("testing");

# colors for the kernel and misfit function
@cols = split("/",$colors);
$kpwr = $cols[0];  # power for kernels
$kmax = $cols[1];  # value for kernels
$opwr = $cols[2];  # power for misfit maps
$omax = $cols[3];  # value for misfit maps
$cmax = $cols[4];  # value for phase speed maps (percent pert)

#@files = glob("$dir/$kerfile");
#$numk = @files;

@titles = ("Waveform","Traveltime (xcorr), misfit","Amplitude (xcorr), misfit","Traveltime (MT), misfit","Amplitude (MT), misfit","Traveltime (xcorr), sampling","Amplitude (xcorr), sampling");
@units = ("m\@+2\@+ s","s\@+2\@+","xxx","xxx","xxx","xxx","xxx");
$ktype = $titles[$iker];
$utype = $units[$iker];

$plabel = "/home/carltape/sem2d/2d_adjoint/scripts/plot_ker_mod.pl";

# data files
$ishelf = 0;
$shelf_file = "../../INPUT/oms_shelf";

# plotting specifications
$fsize0 = "18";
$fsize1 = "14";
$fsize2 = "12";
$fsize3 = "8";
$fontno = "4";    # 1 (helvetica) or 4 (roman)
$tick   = "0.15c";
$fpen   = "1.5p";
$tpen   = "1.0p";

# plot symbols for sources, receivers, and shelf
$src          = "-Sa12p -W1p";
$src_r        = "-Sa12p -W1p/255/0/0 ";
$src_b        = "-Sa12p -W1p/0/0/255 ";
$rec          = "-Sc5p -W0.5p/0/0/0";
$rec_r        = "-Sc5p -W0.5p/255/0/0";
$rec_w        = "-Sc5p -W0.5p/255";

$src_fr       = "$src -G255/0/0";
$src_fw       = "$src -G255";
$src_fk       = "$src -G0";

$rec_fr       = "$rec -G255/0/0";
$rec_fw       = "$rec -G255";
$rec_fk       = "$rec -G0";

$src1   = $src;
$rec1   = "-St12p -W1p";

$Wshelf = "-W0.5/0/0/0tap";
$coast_info = "-W1p -Na/1p -Df";
$coast_info2 = "-W0.5p -Na/0.5p -Df";

#-------------------------
# color for kernels

# KEY: scaling for color
$scale_color = 41.0;
$colorbar = "seis";
#open(IN,"$dir/$cfile");
#@colors = <IN>;

$norm2 = "1e$kpwr";
$ss = $kmax;
$kmin = -$kmax;
$kran = $kmax-$kmin;
$ds = $kran / ${scale_color};
#$ds = 2*$ss/$scale_color;
$bs2 = sprintf("%2.2e",0.9*$ss); # colorbar
$T2 = sprintf("-T%3.3e/%3.3e/%3.3e",-$ss,$ss,$ds);
print "T2 = $T2\n";

#-------------------------
# color for chi (objective function)

$norm3 = "1e$opwr";
$ss = $omax;
$ds = $ss/$scale_color;
$bs3 = sprintf("%3.3e",0.9*$ss/3); # colorbar
$T3 = sprintf("-T%3.3e/%3.3e/%3.3e",0,$ss,$ds);
print "T3 = $T3\n";

#-------------------------
# get reference phase speed and period

$c0 = 3500; $per = 20;

open(IN,"${dir0}/socal_vel_c0.dat");
@vals = <IN>;
$c0       = $vals[0];
$per      = $vals[1];
$lam      = $c0*$per;
$per_lab  = sprintf("T = %.1f s",$per);
$c0_lab   =  sprintf("c0 = %.2f km/s",$c0/1000);
$lam_lab  = sprintf("\@~\154\@~ = %.0f km",$lam/1000);
print "\n$per_lab, $c0_lab, $lam_lab \n";

#-------------------------

# write plotting scripts
$Jwid = 3;
$J = "-JM${Jwid}i";      # in lat-lon
$origin = "-X0.5 -Y4";

# which borders to plot the lat-lon
# 1 four sides, 4 single sides, 6 two sides, 4 three sides, 1 zero sides
@Bopts = ("WESN","Wesn","wesN","wEsn","weSn","WesN","wEsN","wESn","WeSn","WEsn","weSN","wESN","WESn","WeSN","WEsN","wesn");
$B0 = "-B1:.\" \":";

# axes scale for phase speed maps: c(th, ph)

#$bs1 = 0.5;
#$Bscale1d  = sprintf("-B%2.2e:\" Phase Speed for data ( km s\@+-1\@+ )\":",$bs1);
#$Bscale1s  = sprintf("-B%2.2e:\" Phase Speed for model $smod ( km s\@+-1\@+ )\":",$bs1);
$bs1 = 3;
$Bscale1  = sprintf("-B%2.2ef1:\" \@~\045\@~ pert. from %2.2f km/s\": -E10p",$bs1,$c0/1000);
$Bscale1b = sprintf("-B%2.2ef1:\" \": -E10p",$bs1);

$Bscale_spl  = "-B0.5:\" B ( \@~\161\@~, \@~\146\@~ )\":";

# axes scale for kernels: K(th, ph)
$tp = "\@~\146\@~, \@~\161\@~";
$Bscale2  = sprintf("-B${bs2}:\" K ( $tp )  ( 10\@+%2.2i\@+  m\@+-2\@+ [\@~\143\@~] )\": -E10p",$kpwr);
$Bscale2b = sprintf("-B${bs2}:\" K ( $tp )  ( 10\@+%2.2i\@+  m\@+-2\@+ s )\": -E10p",$kpwr);
$Bscale2c = "-B${bs2}:\" \": -E10p";

# axes scale for chi_plots: chi(th_r, ph_r)
$Bscale3  = sprintf("-B%2.2e:\" \@~\143\@~ ( \@~\146\@~\@-r\@- , \@~\161\@~\@-r\@- )  ( 10\@+%2.2i\@+ )\": -Ef10p",$bs3,$opwr);

#-------------------------
# phase speed model
$file1dat = "${dir}/$mfile_dat";
$file1syn = "${dir}/$mfile_syn";

# set bounds for the plotting
#$xmin = -121; $xmax = -115.0; $zmin = 31.75; $zmax = 36.5;
($xmin,$xmax,$zmin,$zmax) = split(" ",`minmax -C $file1dat`);
$dinc = 0.25;  # buffer around min/max
$xmin = $xmin-$dinc;  $xmax = $xmax+$dinc;
$zmin = $zmin-$dinc;  $zmax = $zmax+$dinc;
$R = "-R$xmin/$xmax/$zmin/$zmax";

# plot title
$J_title = "-JX${Jwid}";
$R_title = "-R0/1/0/1";
$x_title = 0.5;
$z_title = 1.08;
$fsize_title = $fsize1;

#=================================================================
# files and parameters related to polynomial plots and misfit-vs-iteration plots

  $p_info_k   = "-N -Sc10p -W1.0p -G0";
  $p_info_w   = "-N -Sc10p -W1.0p -G255";
  $p_info_r   = "-N -Sc10p -W1.0p -G255/0/0";
  $c_info_ks  = "-W1.0p";
  $c_info_ks2 = "-W2.0p";
  $c_info_rd  = "-W1.0/255/0/0tap";
  $c_info_bd  = "-W1.0/0/0/255tap";
  $c_info_kd  = "-W0.5/0/0/0tap";

if ($ipoly==1) {

  # axes limits for polynomial plots
  $axes_file = "axes_${strun0}.dat";
  if (not -f $axes_file)   { die("Check if $axes_file exist or not\n") }
  open(IN,"$axes_file");
  @ax_lims = <IN>;
  ($a1x1,$a1x2,$a1y1,$a1y2) = split(" ",$ax_lims[0]);  # polynomial
  ($a2x1,$a2x2,$a2y1,$a2y2) = split(" ",$ax_lims[1]);  # chi vs iteration
  ($a3x1,$a3x2,$a3y1,$a3y2) = split(" ",$ax_lims[2]);  # var-red vs iteration

  # scale the axes
  #$lpwr_x = 3; $xtick = 4;
  #$lpwr_y = 2; $ytick = 2;
  $lpwr_x = int( log($a1x2) / log(10) ); $xtick = 1;
  $lpwr_y = int( log($a1y2) / log(10) ); $ytick = 1;
  $normgx = "1e$lpwr_x";
  $normgy = "1e$lpwr_y";

  #die("testing");

  $xming = $a1x1/$normgx;
  $xmaxg = $a1x2/$normgx;
  $yming = $a1y1/$normgy;
  $ymaxg = $a1y2/$normgy;
  $R_poly = "-R$xming/$xmaxg/$yming/$ymaxg";
  $J_poly_wid = 0.7*$Jwid;
  $J_poly = "-JX${J_poly_wid}/${Jwid}";

  # polynomial lines and curves
  $poly_curve = "poly_curve_${strun0}.dat";
  $poly_pts   = "poly_points_${strun0}.dat";
  if (not -f $poly_curve) { die("Check if $poly_curve exist or not\n") }
  if (not -f $poly_pts)   { die("Check if $poly_pts exist or not\n") }
  open(IN,"$poly_pts");
  @poly_pts = <IN>;
  ($x1t,$y1t) = split(" ",$poly_pts[0]); $x1 = $x1t/$normgx; $y1 = $y1t/$normgy;
  ($x2t,$y2t) = split(" ",$poly_pts[1]); $x2 = $x2t/$normgx; $y2 = $y2t/$normgy;
  ($x3t,$y3t) = split(" ",$poly_pts[2]); $x3 = $x3t/$normgx; $y3 = $y3t/$normgy;
  ($x4t,$y4t) = split(" ",$poly_pts[3]); $x4 = $x4t/$normgx; $y4 = $y4t/$normgy;
  ($x5t,$y5t) = split(" ",$poly_pts[4]); $x5 = $x5t/$normgx; $y5 = $y5t/$normgy;
  ($x6t,$y6t) = split(" ",$poly_pts[5]); $x6 = $x6t/$normgx; $y6 = $y6t/$normgy;
  ($x7t,$y7t) = split(" ",$poly_pts[6]); $x7 = $x7t/$normgx; $y7 = $y7t/$normgy;

  # chi-vs-iteration
  $xming = $a2x1;
  $xmaxg = $a2x2;
  $yming = $a2y1;
  $ymaxg = $a2y2;
  $R_chi = "-R$xming/$xmaxg/$yming/$ymaxg";
  $J_chi = "-JX${Jwid}i/${Jwid}il";        # note log scale on y-axis

  $schi0 = "\@~\143\@~ ( m\@+0\@+ )";
  $schi1 = "\@~\143\@~ ( m\@+1\@+ )";

  $chi_curve = "chi_curve_${strun0}.dat";
  if (not -f $chi_curve)   { die("Check if $chi_curve exist or not\n") }

  # var.red-vs-iteration
  $xming = $a3x1;
  $xmaxg = $a3x2;
  $yming = $a3y1;
  $ymaxg = $a3y2;
  $R_var = "-R$xming/$xmaxg/$yming/$ymaxg";
  $J_var = "-JX$Jwid";

  $var_pts = "var_${strun0}.dat";
  $var_fit = "var_fit_${strun0}.dat";
  if (not -f $var_fit)   { die("Check if $var_fit exist or not\n") }
  if (not -f $var_pts)   { die("Check if $var_pts exist or not\n") }

  # scale for polynomial plots (chi0, chi1)
  $B1 = sprintf("-B${xtick}:\"\@~\154\@~   ( 10\@+%2.2i\@+ ) \":/${ytick}:\" \@~\143\@~\@+0\@+ [ m (\@~\154\@~) ]   ( 10\@+%2.2i\@+  $utype ) \":",$lpwr_x,$lpwr_y);
  $B2 = sprintf("-B${xtick}:\"\@~\154\@~   ( 10\@+%2.2i\@+ ) \":/${ytick}:\" \@~\143\@~\@+1\@+ [ m (\@~\154\@~) ]   ( 10\@+%2.2i\@+  $utype ) \":",$lpwr_x,$lpwr_y);

  # scale for chi-vs-m plots
  $B3a    = "-B1:\" model number (iteration) \":/a1f2p:\" \@~\143\@~ ( m )   ( $utype ) \":";
  $B3b    = "-B1:\" model number (iteration) \":/a1f2p:\" \":";
  $B4     = "-B1:\" model number (iteration) \":/20:\" \":";

}

#===========================================================================
# create colorpoint files

  open(CSH,">$cshfile");

  $cpt_vel_map = "../../model_files/socal_color.cpt";

  # phase speed model
  if($ichecker==0) {$cpt_vel = $cpt_vel_map}
  else {
     # make colorpoint file
     #$T1 = "-T3/4/0.1";
     $cmin = -$cmax;
     $cran = $cmax - $cmin;
     $dc = $cran/${scale_color};
     $T1 = "-T$cmin/$cmax/$dc";
     $cpt_vel = "color0.cpt";
     print CSH "makecpt -C$colorbar $T1 > temp1\n";
     print CSH "sed 's/^B.*/B       170     0      0  /' temp1 >  temp2\n";
     print CSH "sed 's/^F.*/F         0     0    255  /' temp2 > $cpt_vel\n";
  }

  # kernel
  $cpt_ker = "color1.cpt";
  print CSH "makecpt -C$colorbar $T2 > temp1\n";
  print CSH "sed 's/^B.*/B       170     0      0  /' temp1 >  temp2\n";
  print CSH "sed 's/^F.*/F         0     0    255  /' temp2 > $cpt_ker\n";

  # spline
  $cpt_spl = "color5.cpt";
  $cmin = -1; $cmax = 1;
  $cran = $cmax - $cmin;
  $dc = $cran/${scale_color};
  $T_spline = "-T$cmin/$cmax/$dc";
  print CSH "makecpt -C$colorbar $T_spline > $cpt_spl \n";

  $cpt_spl2 = "/home/carltape/sem2d/2d_adjoint_banana/mat_SAVE/Gik_examples/color_spline.cpt";

  # misfit (as a function of receiver)
  $cpt_chi = "color2.cpt";
  print CSH "makecpt -Chot $T3 -I > $cpt_chi\n";

  close (CSH); system("csh -f $cshfile");
  #die("testing");

  # colorbar
  $Dlen = 0.75*${Jwid}; $Dx = $Jwid/2; $Dy = -0.10*$Jwid;
  $Dscale = "-D$Dx/$Dy/$Dlen/0.20h";

#===========================================================================
if ($ifig01 == 1) {

  # shifting the subplots
  $xfac = 1.25;
  $yfac = 1.6;
  $dX = 1.15*$Jwid; $dY = 0; $shift1 = "-X$dX -Y$dY";
  $dX = 1.25*$Jwid; $dY = 0; $shift2 = "-X$dX -Y$dY";

  #===============================================
  print "\nWriting CSH file...\n";

  open(CSH,">$cshfile");
  print CSH "gmtset BASEMAP_TYPE plain PLOT_DEGREE_FORMAT D TICK_LENGTH $tick LABEL_FONT_SIZE $fsize2 ANOT_FONT_SIZE $fsize2  HEADER_FONT $fontno ANOT_FONT $fontno LABEL_FONT $fontno HEADER_FONT_SIZE $fsize1 FRAME_PEN $fpen TICK_PEN $tpen \n";
  #===============================================

  # get the summed kernel and chi files
  #$file2 = "${dir}/summed_ker.dat";
  $file2 = "${dir}/fun_smooth.dat";       # smoothed kernel
  $file3 = "${dir}/summed_chi_r.dat";
  $recfile = "${dir}/${edir}001/sr.txt";  # src-rec for first event
  $evefile = "${dir0}/events_lonlat.dat";
  $chifile = "${dir}/summed_chi_all.dat";

  if (not -f $file2)   { die("Check if $file2 exist or not\n") }
  if (not -f $file3)   { die("Check if $file3 exist or not\n") }
  if (not -f $recfile) { die("Check if $recfile exist or not\n") }
  if (not -f $chifile) { die("Check if $chifile exist or not\n") }

  #=======================================================================
  # file names for figures
  $name = "talk_cg01_1"; $psfile = "$name.ps"; $jpgfile = "$name.jpg"; $epsfile = "$name.eps";

  # number of receivers
  open(IN,$file3); @temp = <IN>; $nrec = @temp;
  # $npath = $nrec * $nevent;

  #=============================================
  # model for the data

  $B = "$B0".$Bopts[15];
  $title = "Phase speed for data";

  # phase speed map
  print CSH "psbasemap $B $R $J -K -V $origin > $psfile\n"; # START
  #if($icolor==1) {print CSH "awk '{print \$1,\$2,\$3/1000}' $file1dat | pscontour $R $J -A- -C$cpt_vel -I -K -O -V >> $psfile\n"}
  if($icolor==1) {print CSH "awk '{print \$1,\$2,\$4*100}' $file1dat | pscontour $R $J -A- -C$cpt_vel -I -K -O -V  >> $psfile\n"}
  print CSH "pscoast $J $R $coast_info -K -O -V >> $psfile\n";
  if($ishelf) {print CSH "awk '{print \$2,\$1}' $shelf_file |psxy $J $R $Wshelf -K -O -V >> $psfile\n"}
  print CSH "awk '\$1 == \"R\" {print \$2,\$3}' $recfile |psxy -N $J $R -K -O -V $rec >> $psfile\n";
  print CSH "awk '{print \$1,\$2}' $evefile |psxy -N $J $R -K -O -V $src >> $psfile\n";
  print CSH "psscale -C$cpt_vel $Dscale $Bscale1 -K -O -V >> $psfile \n";
  print CSH "pstext -N $J_title $R_title -K -O -V >>$psfile<<EOF\n $x_title $z_title $fsize_title 0 $fontno CM $title \nEOF\n";

  #=============================================
  # model for synthetics

  $title = "Phase speed model $smod";
  $shift = $shift1;

  # phase speed map
  print CSH "psbasemap $B $R $J -K -O -V $shift1 >> $psfile\n";
  #if($icolor==1) {print CSH "awk '{print \$1,\$2,\$3/1000}' $file1syn | pscontour $R $J -A- -C$cpt_vel -I -O -K -V >> $psfile\n"}
  if($icolor==1) {print CSH "awk '{print \$1,\$2,\$4*100}' $file1syn | pscontour $R $J -A- -C$cpt_vel -I -O -K -V >> $psfile\n"}
  print CSH "pscoast $J $R $coast_info -K -O -V >> $psfile\n";
  if($ishelf) {print CSH "awk '{print \$2,\$1}' $shelf_file |psxy $J $R $Wshelf -K -O -V >> $psfile\n"}
  print CSH "awk '{print \$1,\$2}' $evefile |psxy -N $J $R -K -O -V $src >> $psfile\n";
  print CSH "awk '\$1 == \"R\" {print \$2,\$3}' $recfile |psxy -N $J $R -K -O -V $rec >> $psfile\n";
  print CSH "psscale -C$cpt_vel $Dscale $Bscale1 -K -O -V >> $psfile \n";
  print CSH "pstext -N $J_title $R_title -K -O -V >>$psfile<<EOF\n $x_title $z_title $fsize_title 0 $fontno CM $title \nEOF\n";

  # blank template
  $shift = $shift2;
  $title = "CG line search";
  $Bg = "$B1".$Bopts[8];

  print CSH "psbasemap $Bg $R_poly $J_poly -K -O -V $shift >> $psfile\n";
  #print CSH "psxy $J_poly $R_poly $p_info_k -K -O -V >>$psfile<<EOF\n$x1 $y1 \nEOF\n";
  #print CSH "psxy $J_poly $R_poly $c_info_ks -K -O -V >>$psfile<<EOF\n$x1 $a1y1 \n$x1 $y1 \n$a1x1 $y1 \nEOF\n";
  print CSH "pstext -N $J_poly $R_title -K -O -V >>$psfile<<EOF\n $x_title $z_title $fsize_title 0 $fontno CM $title \nEOF\n";

  #-----------------------------
  print CSH "pstext -N $J_title $R_title -O -V >>$psfile<<EOF\n 10 10 $fsize_title 0 $fontno CM $title \nEOF\n"; # FINISH
  print CSH "convert $psfile -rotate 90 $jpgfile\n";
  if($ieps==1){print CSH "convert $psfile $epsfile\n"}

  #=======================================================================
  # file names for figures
  $name = "talk_cg01_2"; $psfile = "$name.ps"; $jpgfile = "$name.jpg"; $epsfile = "$name.eps";

  #=============================================
  # model for the data

  $B = "$B0".$Bopts[15];
  $title = "Phase speed for data";

  # phase speed map
  print CSH "psbasemap $B $R $J -K -V $origin > $psfile\n"; # START
  if($icolor==1) {print CSH "awk '{print \$1,\$2,\$4*100}' $file1dat | pscontour $R $J -A- -C$cpt_vel -I -K -O -V  >> $psfile\n"}
  print CSH "pscoast $J $R $coast_info -K -O -V >> $psfile\n";
  if($ishelf) {print CSH "awk '{print \$2,\$1}' $shelf_file |psxy $J $R $Wshelf -K -O -V >> $psfile\n"}
  print CSH "awk '\$1 == \"R\" {print \$2,\$3}' $recfile |psxy -N $J $R -K -O -V $rec >> $psfile\n";
  print CSH "awk '{print \$1,\$2}' $evefile |psxy -N $J $R -K -O -V $src >> $psfile\n";
  print CSH "psscale -C$cpt_vel $Dscale $Bscale1 -K -O -V >> $psfile \n";
  print CSH "pstext -N $J_title $R_title -K -O -V >>$psfile<<EOF\n $x_title $z_title $fsize_title 0 $fontno CM $title \nEOF\n";

  #=============================================
  # chi as a function of receiver

  # total misfit for the model
  open(IN,"$chifile");
  $chi = <IN>;
  $schia = "\@~\143\@~ ( $smod )";
  $schi = sprintf("${schia}  =  %3.3e",$chi);
  $title = "Summed misfit : $schi";
  $B = "$B0".$Bopts[15];
  $shift = $shift1;

  print CSH "psbasemap $B $R $J -K -O -V $shift >> $psfile\n";
  #print CSH "awk '{print \$1,\$2,\$3 / $norm3 }' $file > temp\n";
  print CSH "awk '{print \$1,\$2,\$3 / $norm3 }' $file3 | pscontour $J $R -A- -C$cpt_chi -I -K -O -V >> $psfile\n";
  print CSH "pscoast $J $R $coast_info -K -O -V >> $psfile\n";
  if($ishelf) {print CSH "awk '{print \$2,\$1}' $shelf_file |psxy $J $R $Wshelf -K -O -V >> $psfile\n"}
  print CSH "psscale -C$cpt_chi $Dscale $Bscale3 -K -O -V >> $psfile \n";
  print CSH "awk '{print \$1,\$2}' $evefile |psxy -N $J $R -K -O -V $src >> $psfile\n";
  print CSH "awk '\$1 == \"R\" {print \$2,\$3}' $recfile |psxy -N $J $R -K -O -V $rec >> $psfile\n";
  print CSH "pstext -N $J_title $R_title -K -O -V >>$psfile<<EOF\n $x_title $z_title $fsize_title 0 $fontno CM $title \nEOF\n";

  #=============================================
  # graph

  $shift = $shift2;
  $title = "Plot  $schia";
  $Bg = "$B1".$Bopts[8];

  print CSH "psbasemap $Bg $R_poly $J_poly -K -O -V $shift >> $psfile\n";
  print CSH "psxy $J_poly $R_poly $p_info_k -K -O -V >>$psfile<<EOF\n$x1 $y1 \nEOF\n";
  print CSH "psxy $J_poly $R_poly $c_info_ks -K -O -V >>$psfile<<EOF\n$x1 $a1y1 \n$x1 $y1 \n$a1x1 $y1 \nEOF\n";
  print CSH "pstext -N $J_poly $R_title -K -O -V >>$psfile<<EOF\n $x_title $z_title $fsize_title 0 $fontno CM $title \nEOF\n";
  #-----------------------------
  print CSH "pstext -N $J_title $R_title -O -V >>$psfile<<EOF\n 10 10 $fsize_title 0 $fontno CM $title \nEOF\n"; # FINISH
  print CSH "convert $psfile -rotate 90 $jpgfile\n";
  if($ieps==1){print CSH "convert $psfile $epsfile\n"}

  #=======================================================================
  # file names for figures
  $name = "talk_cg01_3"; $psfile  = "$name.ps"; $jpgfile = "$name.jpg"; $epsfile = "$name.eps";

  #=============================================
  # model for the data

  $B = "$B0".$Bopts[15];
  $title = "Phase speed for data";

  # phase speed map
  print CSH "psbasemap $B $R $J -K -V $origin > $psfile\n"; # START
  if($icolor==1) {print CSH "awk '{print \$1,\$2,\$4*100}' $file1dat | pscontour $R $J -A- -C$cpt_vel -I -K -O -V  >> $psfile\n"}
  print CSH "pscoast $J $R $coast_info -K -O -V >> $psfile\n";
  if($ishelf) {print CSH "awk '{print \$2,\$1}' $shelf_file |psxy $J $R $Wshelf -K -O -V >> $psfile\n"}
  print CSH "awk '\$1 == \"R\" {print \$2,\$3}' $recfile |psxy -N $J $R -K -O -V $rec >> $psfile\n";
  print CSH "awk '{print \$1,\$2}' $evefile |psxy -N $J $R -K -O -V $src >> $psfile\n";
  print CSH "psscale -C$cpt_vel $Dscale $Bscale1 -K -O -V >> $psfile \n";
  print CSH "pstext -N $J_title $R_title -K -O -V >>$psfile<<EOF\n $x_title $z_title $fsize_title 0 $fontno CM $title \nEOF\n";

  #=============================================
  # kernel

  $shift = $shift1;
  $title = "Gradient for $smod";
  $B = "$B0".$Bopts[15];

  print CSH "psbasemap $B $R $J -K -O -V $shift >> $psfile\n";
  #if($icolor==1) {print CSH "awk '{print \$1,\$2,\$3 / $norm2 }' $file2 | pscontour $J $R -A- -C$cpt_ker -I -K -O -V >> $psfile\n"}
  if($icolor==1) {print CSH "awk '{print \$1,\$2,\$7 / $norm2 }' $file2 | pscontour $J $R -A- -C$cpt_ker -I -K -O -V >> $psfile\n"}
  print CSH "pscoast $J $R $coast_info -K -O -V >> $psfile\n";
  if($ishelf) {print CSH "awk '{print \$2,\$1}' $shelf_file |psxy $J $R $Wshelf -K -O -V >> $psfile\n"}
  print CSH "psscale -C$cpt_ker $Dscale $Bscale2 -K -O -V >> $psfile \n";
  print CSH "awk '{print \$1,\$2}' $evefile |psxy -N $J $R -K -O -V $src >> $psfile\n";
  print CSH "pstext -N $J_title $R_title -K -O -V >>$psfile<<EOF\n $x_title $z_title $fsize_title 0 $fontno CM $title \nEOF\n";
  print CSH "awk '\$1 == \"R\" {print \$2,\$3}' $recfile |psxy -N $J $R -K -O -V $rec >> $psfile\n";

  #=============================================
  # graph

  $shift = $shift2;
  $Bg = "$B1".$Bopts[8];
  $title = "Plot  g\@+0\@+ = d\@~\143\@~ / dm($smod)";

  print CSH "psbasemap $Bg $R_poly $J_poly -K -O -V $shift >> $psfile\n";
  print CSH "awk '{print \$1 / $normgx,\$2 / $normgy}' $poly_curve | psxy $c_info_rd $J_poly $R_poly -K -O -V >> $psfile\n";
  print CSH "psxy $J_poly $R_poly $c_info_ks -K -O -V >>$psfile<<EOF\n$x1 $a1y1 \n$x1 $y1 \n$a1x1 $y1 \nEOF\n";
  print CSH "psxy $J_poly $R_poly $p_info_k -K -O -V >>$psfile<<EOF\n$x1 $y1 \nEOF\n";
  print CSH "pstext -N $J_poly $R_title -K -O -V >>$psfile<<EOF\n $x_title $z_title $fsize_title 0 $fontno CM $title \nEOF\n";

  #-----------------------------
  print CSH "pstext -N $J_title $R_title -O -V >>$psfile<<EOF\n 10 10 $fsize_title 0 $fontno CM $title \nEOF\n"; # FINISH
  print CSH "convert $psfile -rotate 90 $jpgfile\n";
  if($ieps==1){print CSH "convert $psfile $epsfile\n"}

  #=======================================================================
  # file names for figures
  $name = "talk_cg01_4"; $psfile  = "$name.ps"; $jpgfile = "$name.jpg"; $epsfile = "$name.eps";

  #=============================================
  # model for the data

  $B = "$B0".$Bopts[15];
  $title = "Phase speed for data";

  # phase speed map
  print CSH "psbasemap $B $R $J -K -V $origin > $psfile\n"; # START
  if($icolor==1) {print CSH "awk '{print \$1,\$2,\$4*100}' $file1dat | pscontour $R $J -A- -C$cpt_vel -I -K -O -V  >> $psfile\n"}
  print CSH "pscoast $J $R $coast_info -K -O -V >> $psfile\n";
  if($ishelf) {print CSH "awk '{print \$2,\$1}' $shelf_file |psxy $J $R $Wshelf -K -O -V >> $psfile\n"}
  print CSH "awk '\$1 == \"R\" {print \$2,\$3}' $recfile |psxy -N $J $R -K -O -V $rec >> $psfile\n";
  print CSH "awk '{print \$1,\$2}' $evefile |psxy -N $J $R -K -O -V $src >> $psfile\n";
  print CSH "psscale -C$cpt_vel $Dscale $Bscale1 -K -O -V >> $psfile \n";
  print CSH "pstext -N $J_title $R_title -K -O -V >>$psfile<<EOF\n $x_title $z_title $fsize_title 0 $fontno CM $title \nEOF\n";

#=================================================================
# NEW MODEL

  # wave2d run number
  $strun = sprintf("%4.4i",$irun + 1);
  $dir = "${odir}$strun";

  $mod = $mods[1];
  $smod = "m\@+$mod\@+";

  # get the summed kernel and chi files
  $file1syn = "${dir}/$mfile_syn";
  $file2 = "${dir}/fun_smooth.dat";       # smoothed kernel
  #$file2 = "$dir/summed_ker.dat";
  $file3 = "${dir}/summed_chi_r.dat";
  $recfile = "${dir}/${edir}001/sr.txt";  # src-rec for first event
  $chifile = "${dir}/summed_chi_all.dat";

  if (not -f $file1syn){ die("Check if $file1syn exist or not\n") }
  if (not -f $file2)   { die("Check if $file2 exist or not\n") }
  if (not -f $file3)   { die("Check if $file3 exist or not\n") }
  if (not -f $recfile) { die("Check if $recfile exist or not\n") }
  if (not -f $chifile) { die("Check if $chifile exist or not\n") }
#=================================================================

  #=============================================
  # model for synthetics

  $title = "Phase speed model $smod";
  $shift = $shift1;

  # phase speed map
  print CSH "psbasemap $B $R $J -K -O -V $shift >> $psfile\n";
  if($icolor==1) {print CSH "awk '{print \$1,\$2,\$4*100}' $file1syn | pscontour $R $J -A- -C$cpt_vel -I -O -K -V >> $psfile\n"}
  print CSH "pscoast $J $R $coast_info -K -O -V >> $psfile\n";
  if($ishelf) {print CSH "awk '{print \$2,\$1}' $shelf_file |psxy $J $R $Wshelf -K -O -V >> $psfile\n"}
  print CSH "awk '{print \$1,\$2}' $evefile |psxy -N $J $R -K -O -V $src >> $psfile\n";
  print CSH "awk '\$1 == \"R\" {print \$2,\$3}' $recfile |psxy -N $J $R -K -O -V $rec >> $psfile\n";
  print CSH "psscale -C$cpt_vel $Dscale $Bscale1 -K -O -V >> $psfile \n";
  print CSH "pstext -N $J_title $R_title -K -O -V >>$psfile<<EOF\n $x_title $z_title $fsize_title 0 $fontno CM $title \nEOF\n";

  #=============================================
  # graph

  $shift = $shift2;
  $Bg = "$B1".$Bopts[8];
  $title = "Estimate test model m\@+0t\@+";

  print CSH "psbasemap $Bg $R_poly $J_poly -K -O -V $shift >> $psfile\n";
  print CSH "awk '{print \$1 / $normgx,\$2 / $normgy}' $poly_curve | psxy $c_info_rd $J_poly $R_poly -K -O -V >> $psfile\n";
  print CSH "awk '{print \$1 / $normgx,\$4 / $normgy}' $poly_curve | psxy $c_info_bd $J_poly $R_poly -K -O -V >> $psfile\n";
  print CSH "psxy $J_poly $R_poly $c_info_ks -K -O -V >>$psfile<<EOF\n$x1 $a1y1 \n$x1 $y1 \n$a1x1 $y1 \nEOF\n";
  print CSH "psxy $J_poly $R_poly $p_info_k -K -O -V >>$psfile<<EOF\n$x1 $y1 \nEOF\n";
  print CSH "psxy $J_poly $R_poly $p_info_w -K -O -V >>$psfile<<EOF\n$x2 $a1y1 \nEOF\n";
  print CSH "pstext -N $J_poly $R_title -K -O -V >>$psfile<<EOF\n $x_title $z_title $fsize_title 0 $fontno CM $title \nEOF\n";

  #-----------------------------
  print CSH "pstext -N $J_title $R_title -O -V >>$psfile<<EOF\n 10 10 $fsize_title 0 $fontno CM $title \nEOF\n"; # FINISH
  print CSH "convert $psfile -rotate 90 $jpgfile\n";
  if($ieps==1){print CSH "convert $psfile $epsfile\n"}

  #=======================================================================
  # file names for figures
  $name = "talk_cg01_5"; $psfile  = "$name.ps"; $jpgfile = "$name.jpg"; $epsfile = "$name.eps";

  #=============================================
  # model for the data

  $B = "$B0".$Bopts[15];
  $title = "Phase speed for data";

  # phase speed map
  print CSH "psbasemap $B $R $J -K -V $origin > $psfile\n"; # START
  if($icolor==1) {print CSH "awk '{print \$1,\$2,\$4*100}' $file1dat | pscontour $R $J -A- -C$cpt_vel -I -K -O -V  >> $psfile\n"}
  print CSH "pscoast $J $R $coast_info -K -O -V >> $psfile\n";
  if($ishelf) {print CSH "awk '{print \$2,\$1}' $shelf_file |psxy $J $R $Wshelf -K -O -V >> $psfile\n"}
  print CSH "awk '\$1 == \"R\" {print \$2,\$3}' $recfile |psxy -N $J $R -K -O -V $rec >> $psfile\n";
  print CSH "awk '{print \$1,\$2}' $evefile |psxy -N $J $R -K -O -V $src >> $psfile\n";
  print CSH "psscale -C$cpt_vel $Dscale $Bscale1 -K -O -V >> $psfile \n";
  print CSH "pstext -N $J_title $R_title -K -O -V >>$psfile<<EOF\n $x_title $z_title $fsize_title 0 $fontno CM $title \nEOF\n";

  #=============================================
  # chi as a function of receiver

  # total misfit for the model
  open(IN,"$chifile");
  $chi = <IN>;
  $schia = "\@~\143\@~ ( $smod )";
  $schi = sprintf("${schia}  =  %3.3e",$chi);
  $title = "Summed misfit : $schi";
  $B = "$B0".$Bopts[15];
  $shift = $shift1;

  print CSH "psbasemap $B $R $J -K -O -V $shift >> $psfile\n";
  #print CSH "awk '{print \$1,\$2,\$3 / $norm3 }' $file > temp\n";
  print CSH "awk '{print \$1,\$2,\$3 / $norm3 }' $file3 | pscontour $J $R -A- -C$cpt_chi -I -K -O -V >> $psfile\n";
  print CSH "pscoast $J $R $coast_info -K -O -V >> $psfile\n";
  if($ishelf) {print CSH "awk '{print \$2,\$1}' $shelf_file |psxy $J $R $Wshelf -K -O -V >> $psfile\n"}
  print CSH "psscale -C$cpt_chi $Dscale $Bscale3 -K -O -V >> $psfile \n";
  print CSH "awk '{print \$1,\$2}' $evefile |psxy -N $J $R -K -O -V $src >> $psfile\n";
  print CSH "awk '\$1 == \"R\" {print \$2,\$3}' $recfile |psxy -N $J $R -K -O -V $rec >> $psfile\n";
  print CSH "pstext -N $J_title $R_title -K -O -V >>$psfile<<EOF\n $x_title $z_title $fsize_title 0 $fontno CM $title \nEOF\n";

  #=============================================
  # graph

  $Bg = "$B1".$Bopts[8];
  $title = "Plot  $schia";
  $shift = $shift2;

  print CSH "psbasemap $Bg $R_poly $J_poly -K -O -V $shift >> $psfile\n";
  print CSH "awk '{print \$1 / $normgx,\$2 / $normgy}' $poly_curve | psxy $c_info_rd $J_poly $R_poly -K -O -V >> $psfile\n";
  print CSH "awk '{print \$1 / $normgx,\$4 / $normgy}' $poly_curve | psxy $c_info_bd $J_poly $R_poly -K -O -V >> $psfile\n";
  print CSH "psxy $J_poly $R_poly $c_info_ks -K -O -V >>$psfile<<EOF\n$x1 $a1y1 \n$x1 $y1 \n$a1x1 $y1 \nEOF\n";
  print CSH "psxy $J_poly $R_poly $c_info_ks -K -O -V >>$psfile<<EOF\n$x3 $a1y1 \n$x3 $y3 \n$a1x1 $y3 \nEOF\n";
  print CSH "psxy $J_poly $R_poly $p_info_k -K -O -V >>$psfile<<EOF\n$x1 $y1 \n$x3 $y3\nEOF\n";
  print CSH "psxy $J_poly $R_poly $p_info_w -K -O -V >>$psfile<<EOF\n$x2 $a1y1 \nEOF\n";
  print CSH "pstext -N $J_poly $R_title -K -O -V >>$psfile<<EOF\n $x_title $z_title $fsize_title 0 $fontno CM $title \nEOF\n";

  #-----------------------------
  print CSH "pstext -N $J_poly $R_title -O -V >>$psfile<<EOF\n 10 10 $fsize_title 0 $fontno CM $title \nEOF\n"; # FINISH
  print CSH "convert $psfile -rotate 90 $jpgfile\n";
  if($ieps==1){print CSH "convert $psfile $epsfile\n"}

  #=======================================================================
  # file names for figures
  $name = "talk_cg01_6"; $psfile  = "$name.ps"; $jpgfile = "$name.jpg"; $epsfile = "$name.eps";

  #=============================================
  # model for the data

  $B = "$B0".$Bopts[15];
  $title = "Phase speed for data";

  # phase speed map
  print CSH "psbasemap $B $R $J -K -V $origin > $psfile\n"; # START
  if($icolor==1) {print CSH "awk '{print \$1,\$2,\$4*100}' $file1dat | pscontour $R $J -A- -C$cpt_vel -I -K -O -V  >> $psfile\n"}
  print CSH "pscoast $J $R $coast_info -K -O -V >> $psfile\n";
  if($ishelf) {print CSH "awk '{print \$2,\$1}' $shelf_file |psxy $J $R $Wshelf -K -O -V >> $psfile\n"}
  print CSH "awk '\$1 == \"R\" {print \$2,\$3}' $recfile |psxy -N $J $R -K -O -V $rec >> $psfile\n";
  print CSH "awk '{print \$1,\$2}' $evefile |psxy -N $J $R -K -O -V $src >> $psfile\n";
  print CSH "psscale -C$cpt_vel $Dscale $Bscale1 -K -O -V >> $psfile \n";
  print CSH "pstext -N $J_title $R_title -K -O -V >>$psfile<<EOF\n $x_title $z_title $fsize_title 0 $fontno CM $title \nEOF\n";

  #=============================================
  # kernel

  $shift = $shift1;
  $title = "Gradient for $smod";
  $B = "$B0".$Bopts[15];

  print CSH "psbasemap $B $R $J -K -O -V $shift >> $psfile\n";
  if($icolor==1) {print CSH "awk '{print \$1,\$2,\$7 / $norm2 }' $file2 | pscontour $J $R -A- -C$cpt_ker -I -K -O -V >> $psfile\n"}
  #if($icolor==1) {print CSH "awk '{print \$1,\$2,\$3 / $norm2 }' $file2 | pscontour $J $R -A- -C$cpt_ker -I -K -O -V >> $psfile\n"}
  print CSH "pscoast $J $R $coast_info -K -O -V >> $psfile\n";
  if($ishelf) {print CSH "awk '{print \$2,\$1}' $shelf_file |psxy $J $R $Wshelf -K -O -V >> $psfile\n"}
  print CSH "psscale -C$cpt_ker $Dscale $Bscale2 -K -O -V >> $psfile \n";
  print CSH "awk '{print \$1,\$2}' $evefile |psxy -N $J $R -K -O -V $src >> $psfile\n";
  print CSH "pstext -N $J_title $R_title -K -O -V >>$psfile<<EOF\n $x_title $z_title $fsize_title 0 $fontno CM $title \nEOF\n";
  print CSH "awk '\$1 == \"R\" {print \$2,\$3}' $recfile |psxy -N $J $R -K -O -V $rec >> $psfile\n";

  #=============================================
  # graph

  $shift = $shift2;
  $Bg = "$B1".$Bopts[8];
  $title = "Plot  g\@+0t\@+ = d\@~\143\@~ / dm($smod)";

  print CSH "psbasemap $Bg $R_poly $J_poly -K -O -V $shift >> $psfile\n";
  print CSH "awk '{print \$1 / $normgx,\$2 / $normgy}' $poly_curve | psxy $c_info_rd $J_poly $R_poly -K -O -V >> $psfile\n";
  print CSH "awk '{print \$1 / $normgx,\$3 / $normgy}' $poly_curve | psxy $c_info_rd $J_poly $R_poly -K -O -V >> $psfile\n";
  #print CSH "awk '{print \$1 / $normgx,\$4 / $normgy}' $poly_curve | psxy $c_info_bd $J_poly $R_poly -K -O -V >> $psfile\n";
  print CSH "psxy $J_poly $R_poly $c_info_ks -K -O -V >>$psfile<<EOF\n$x1 $a1y1 \n$x1 $y1 \n$a1x1 $y1 \nEOF\n";
  print CSH "psxy $J_poly $R_poly $c_info_ks -K -O -V >>$psfile<<EOF\n$x3 $a1y1 \n$x3 $y3 \n$a1x1 $y3 \nEOF\n";
  print CSH "psxy $J_poly $R_poly $p_info_k -K -O -V >>$psfile<<EOF\n$x1 $y1 \n$x3 $y3\nEOF\n";
  print CSH "psxy $J_poly $R_poly $p_info_w -K -O -V >>$psfile<<EOF\n$x2 $a1y1 \nEOF\n";
  print CSH "pstext -N $J_poly $R_title -K -O -V >>$psfile<<EOF\n $x_title $z_title $fsize_title 0 $fontno CM $title \nEOF\n";

#-----------------------------
  print CSH "pstext -N $J_title $R_title -O -V >>$psfile<<EOF\n 10 10 $fsize_title 0 $fontno CM $title \nEOF\n"; # FINISH
  print CSH "convert $psfile -rotate 90 $jpgfile\n";
  if($ieps==1){print CSH "convert $psfile $epsfile\n"}

  #=======================================================================
  # file names for figures
  $name = "talk_cg01_7"; $psfile  = "$name.ps"; $jpgfile = "$name.jpg"; $epsfile = "$name.eps";

  #=============================================
  # model for the data

  $B = "$B0".$Bopts[15];
  $title = "Phase speed for data";

  # phase speed map
  print CSH "psbasemap $B $R $J -K -V $origin > $psfile\n"; # START
  if($icolor==1) {print CSH "awk '{print \$1,\$2,\$4*100}' $file1dat | pscontour $R $J -A- -C$cpt_vel -I -K -O -V  >> $psfile\n"}
  print CSH "pscoast $J $R $coast_info -K -O -V >> $psfile\n";
  if($ishelf) {print CSH "awk '{print \$2,\$1}' $shelf_file |psxy $J $R $Wshelf -K -O -V >> $psfile\n"}
  print CSH "awk '\$1 == \"R\" {print \$2,\$3}' $recfile |psxy -N $J $R -K -O -V $rec >> $psfile\n";
  print CSH "awk '{print \$1,\$2}' $evefile |psxy -N $J $R -K -O -V $src >> $psfile\n";
  print CSH "psscale -C$cpt_vel $Dscale $Bscale1 -K -O -V >> $psfile \n";
  print CSH "pstext -N $J_title $R_title -K -O -V >>$psfile<<EOF\n $x_title $z_title $fsize_title 0 $fontno CM $title \nEOF\n";

#=================================================================
# NEW MODEL

  # wave2d run number
  $strun = sprintf("%4.4i",$irun + 2);
  $dir = "${odir}$strun";

  $mod = $mods[2];
  $smod = "m\@+$mod\@+";

  # get the summed kernel and chi files
  $file1syn = "${dir}/$mfile_syn";
  $file2 = "${dir}/fun_smooth.dat";       # smoothed kernel
  #$file2 = "${dir}/summed_ker.dat";
  $file3 = "${dir}/summed_chi_r.dat";
  $recfile = "${dir}/${edir}001/sr.txt";  # src-rec for first event
  $chifile = "${dir}/summed_chi_all.dat";

  if (not -f $file1syn){ die("Check if $file1syn exist or not\n") }
  if (not -f $file2)   { die("Check if $file2 exist or not\n") }
  if (not -f $file3)   { die("Check if $file3 exist or not\n") }
  if (not -f $recfile) { die("Check if $recfile exist or not\n") }
  if (not -f $chifile) { die("Check if $chifile exist or not\n") }
#=================================================================

  #=============================================
  # model for synthetics

  $title = "Phase speed model $smod";
  $shift = $shift1;

  # phase speed map
  print CSH "psbasemap $B $R $J -K -O -V $shift >> $psfile\n";
  if($icolor==1) {print CSH "awk '{print \$1,\$2,\$4*100}' $file1syn | pscontour $R $J -A- -C$cpt_vel -I -O -K -V >> $psfile\n"}
  print CSH "pscoast $J $R $coast_info -K -O -V >> $psfile\n";
  if($ishelf) {print CSH "awk '{print \$2,\$1}' $shelf_file |psxy $J $R $Wshelf -K -O -V >> $psfile\n"}
  print CSH "awk '{print \$1,\$2}' $evefile |psxy -N $J $R -K -O -V $src >> $psfile\n";
  print CSH "awk '\$1 == \"R\" {print \$2,\$3}' $recfile |psxy -N $J $R -K -O -V $rec >> $psfile\n";
  print CSH "psscale -C$cpt_vel $Dscale $Bscale1 -K -O -V >> $psfile \n";
  print CSH "pstext -N $J_title $R_title -K -O -V >>$psfile<<EOF\n $x_title $z_title $fsize_title 0 $fontno CM $title \nEOF\n";

  #=============================================
  # graph

  $shift = $shift2;
  $Bg = "$B1".$Bopts[8];
  $title = "Cubic interpolation to get m\@+1\@+";

  print CSH "psbasemap $Bg $R_poly $J_poly -K -O -V $shift >> $psfile\n";
  print CSH "awk '{print \$1 / $normgx,\$2 / $normgy}' $poly_curve | psxy $c_info_rd $J_poly $R_poly -K -O -V >> $psfile\n";
  print CSH "awk '{print \$1 / $normgx,\$3 / $normgy}' $poly_curve | psxy $c_info_rd $J_poly $R_poly -K -O -V >> $psfile\n";
  print CSH "awk '{print \$1 / $normgx,\$5 / $normgy}' $poly_curve | psxy $c_info_ks2 $J_poly $R_poly -K -O -V >> $psfile\n";
  print CSH "psxy $J_poly $R_poly $c_info_ks -K -O -V >>$psfile<<EOF\n$x1 $a1y1 \n$x1 $y1 \n$a1x1 $y1 \nEOF\n";
  print CSH "psxy $J_poly $R_poly $c_info_ks -K -O -V >>$psfile<<EOF\n$x3 $a1y1 \n$x3 $y3 \n$a1x1 $y3 \nEOF\n";
  print CSH "psxy $J_poly $R_poly $c_info_kd -K -O -V >>$psfile<<EOF\n$x5 $a1y1 \n$x5 $y5 \n$a1x1 $y5 \nEOF\n";
  print CSH "psxy $J_poly $R_poly $p_info_k -K -O -V >>$psfile<<EOF\n$x1 $y1 \n$x3 $y3\nEOF\n";
  print CSH "psxy $J_poly $R_poly $p_info_w -K -O -V >>$psfile<<EOF\n$x2 $a1y1 \n$x4 $y4 \n$x5 $y5 \nEOF\n";
  print CSH "pstext -N $J_poly $R_title -K -O -V >>$psfile<<EOF\n $x_title $z_title $fsize_title 0 $fontno CM $title \nEOF\n";
  #-----------------------------
  print CSH "pstext -N $J_title $R_title -O -V >>$psfile<<EOF\n 10 10 $fsize_title 0 $fontno CM $title \nEOF\n"; # FINISH
  print CSH "convert $psfile -rotate 90 $jpgfile\n";
  if($ieps==1){print CSH "convert $psfile $epsfile\n"}

  #=======================================================================
  # file names for figures
  $name = "talk_cg01_8"; $psfile  = "$name.ps"; $jpgfile = "$name.jpg"; $epsfile = "$name.eps";

  #=============================================
  # model for the data

  $B = "$B0".$Bopts[15];
  $title = "Phase speed for data";

  # phase speed map
  print CSH "psbasemap $B $R $J -K -V $origin > $psfile\n"; # START
  if($icolor==1) {print CSH "awk '{print \$1,\$2,\$4*100}' $file1dat | pscontour $R $J -A- -C$cpt_vel -I -K -O -V  >> $psfile\n"}
  print CSH "pscoast $J $R $coast_info -K -O -V >> $psfile\n";
  if($ishelf) {print CSH "awk '{print \$2,\$1}' $shelf_file |psxy $J $R $Wshelf -K -O -V >> $psfile\n"}
  print CSH "awk '\$1 == \"R\" {print \$2,\$3}' $recfile |psxy -N $J $R -K -O -V $rec >> $psfile\n";
  print CSH "awk '{print \$1,\$2}' $evefile |psxy -N $J $R -K -O -V $src >> $psfile\n";
  print CSH "psscale -C$cpt_vel $Dscale $Bscale1 -K -O -V >> $psfile \n";
  print CSH "pstext -N $J_title $R_title -K -O -V >>$psfile<<EOF\n $x_title $z_title $fsize_title 0 $fontno CM $title \nEOF\n";

  #=============================================
  # chi as a function of receiver

  # total misfit for the model
  open(IN,"$chifile");
  $chi = <IN>;
  $schia = "\@~\143\@~ ( $smod )";
  $schi = sprintf("${schia}  =  %3.3e",$chi);
  $title = "Summed misfit : $schi";
  $B = "$B0".$Bopts[15];
  $shift = $shift1;

  print CSH "psbasemap $B $R $J -K -O -V $shift >> $psfile\n";
  #print CSH "awk '{print \$1,\$2,\$3 / $norm3 }' $file > temp\n";
  print CSH "awk '{print \$1,\$2,\$3 / $norm3 }' $file3 | pscontour $J $R -A- -C$cpt_chi -I -K -O -V >> $psfile\n";
  print CSH "pscoast $J $R $coast_info -K -O -V >> $psfile\n";
  if($ishelf) {print CSH "awk '{print \$2,\$1}' $shelf_file |psxy $J $R $Wshelf -K -O -V >> $psfile\n"}
  print CSH "psscale -C$cpt_chi $Dscale $Bscale3 -K -O -V >> $psfile \n";
  print CSH "awk '{print \$1,\$2}' $evefile |psxy -N $J $R -K -O -V $src >> $psfile\n";
  print CSH "awk '\$1 == \"R\" {print \$2,\$3}' $recfile |psxy -N $J $R -K -O -V $rec >> $psfile\n";
  print CSH "pstext -N $J_title $R_title -K -O -V >>$psfile<<EOF\n $x_title $z_title $fsize_title 0 $fontno CM $title \nEOF\n";

  #=============================================
  # graph

  $shift = $shift2;
  $Bg = "$B1".$Bopts[8];
  $title = "Plot  $schia";

  print CSH "psbasemap $Bg $R_poly $J_poly -K -O -V $shift >> $psfile\n";
  print CSH "awk '{print \$1 / $normgx,\$2 / $normgy}' $poly_curve | psxy $c_info_rd $J_poly $R_poly -K -O -V >> $psfile\n";
  print CSH "awk '{print \$1 / $normgx,\$3 / $normgy}' $poly_curve | psxy $c_info_rd $J_poly $R_poly -K -O -V >> $psfile\n";
  print CSH "awk '{print \$1 / $normgx,\$5 / $normgy}' $poly_curve | psxy $c_info_ks2 $J_poly $R_poly -K -O -V >> $psfile\n";

  print CSH "psxy $J_poly $R_poly $c_info_ks -K -O -V >>$psfile<<EOF\n$x1 $a1y1 \n$x1 $y1 \n$a1x1 $y1 \nEOF\n";
  print CSH "psxy $J_poly $R_poly $c_info_ks -K -O -V >>$psfile<<EOF\n$x3 $a1y1 \n$x3 $y3 \n$a1x1 $y3 \nEOF\n";
  print CSH "psxy $J_poly $R_poly $c_info_kd -K -O -V >>$psfile<<EOF\n$x5 $a1y1 \n$x5 $y5 \n$a1x1 $y5 \nEOF\n";
  print CSH "psxy $J_poly $R_poly $c_info_ks -K -O -V >>$psfile<<EOF\n$x5 $a1y1 \n$x6 $y6 \n$a1x1 $y6 \nEOF\n";

  print CSH "psxy $J_poly $R_poly $p_info_w -K -O -V >>$psfile<<EOF\n$x2 $a1y1 \n$x4 $y4 \n$x5 $y5 \nEOF\n";
  print CSH "psxy $J_poly $R_poly $p_info_k -K -O -V >>$psfile<<EOF\n$x1 $y1 \n$x3 $y3 \n $x6 $y6 \nEOF\n";
  print CSH "pstext -N $J_poly $R_title -K -O -V >>$psfile<<EOF\n $x_title $z_title $fsize_title 0 $fontno CM $title \nEOF\n";
  #-----------------------------
  print CSH "pstext -N $J_title $R_title -O -V >>$psfile<<EOF\n 10 10 $fsize_title 0 $fontno CM $title \nEOF\n"; # FINISH
  print CSH "convert $psfile -rotate 90 $jpgfile\n";
  if($ieps==1){print CSH "convert $psfile $epsfile\n"}

  #=======================================================================
  # file names for figures
  $name = "talk_cg01_9"; $psfile  = "$name.ps"; $jpgfile = "$name.jpg"; $epsfile = "$name.eps";

  #=============================================
  # model for the data

  $B = "$B0".$Bopts[15];
  $title = "Phase speed for data";

  # phase speed map
  print CSH "psbasemap $B $R $J -K -V $origin > $psfile\n"; # START
  if($icolor==1) {print CSH "awk '{print \$1,\$2,\$4*100}' $file1dat | pscontour $R $J -A- -C$cpt_vel -I -K -O -V  >> $psfile\n"}
  print CSH "pscoast $J $R $coast_info -K -O -V >> $psfile\n";
  if($ishelf) {print CSH "awk '{print \$2,\$1}' $shelf_file |psxy $J $R $Wshelf -K -O -V >> $psfile\n"}
  print CSH "awk '\$1 == \"R\" {print \$2,\$3}' $recfile |psxy -N $J $R -K -O -V $rec >> $psfile\n";
  print CSH "awk '{print \$1,\$2}' $evefile |psxy -N $J $R -K -O -V $src >> $psfile\n";
  print CSH "psscale -C$cpt_vel $Dscale $Bscale1 -K -O -V >> $psfile \n";
  print CSH "pstext -N $J_title $R_title -K -O -V >>$psfile<<EOF\n $x_title $z_title $fsize_title 0 $fontno CM $title \nEOF\n";

  #=============================================
  # model for synthetics

  $title = "Phase speed model $smod";
  $shift = $shift1;

  # phase speed map
  print CSH "psbasemap $B $R $J -K -O -V $shift >> $psfile\n";
  if($icolor==1) {print CSH "awk '{print \$1,\$2,\$4*100}' $file1syn | pscontour $R $J -A- -C$cpt_vel -I -O -K -V >> $psfile\n"}
  print CSH "pscoast $J $R $coast_info -K -O -V >> $psfile\n";
  if($ishelf) {print CSH "awk '{print \$2,\$1}' $shelf_file |psxy $J $R $Wshelf -K -O -V >> $psfile\n"}
  print CSH "awk '{print \$1,\$2}' $evefile |psxy -N $J $R -K -O -V $src >> $psfile\n";
  print CSH "awk '\$1 == \"R\" {print \$2,\$3}' $recfile |psxy -N $J $R -K -O -V $rec >> $psfile\n";
  print CSH "psscale -C$cpt_vel $Dscale $Bscale1 -K -O -V >> $psfile \n";
  print CSH "pstext -N $J_title $R_title -K -O -V >>$psfile<<EOF\n $x_title $z_title $fsize_title 0 $fontno CM $title \nEOF\n";

  #=============================================
  # graph

  $shift = $shift2;
  $Bg = "$B2".$Bopts[8];
  $title = "New line search";

  print CSH "psbasemap $Bg $R_poly $J_poly -K -O -V $shift >> $psfile\n";
  print CSH "psxy $J_poly $R_poly $c_info_ks -K -O -V >>$psfile<<EOF\n$x1 $a1y1 \n$x1 $y6 \n$a1x1 $y6 \nEOF\n";
  print CSH "psxy $J_poly $R_poly $p_info_k -K -O -V >>$psfile<<EOF\n$x1 $y6 \nEOF\n";
  print CSH "pstext -N $J_poly $R_title -K -O -V >>$psfile<<EOF\n $x_title $z_title $fsize_title 0 $fontno CM $title \nEOF\n";
  #-----------------------------
  print CSH "pstext -N $J_title $R_title -O -V >>$psfile<<EOF\n 10 10 $fsize_title 0 $fontno CM $title \nEOF\n"; # FINISH
  print CSH "convert $psfile -rotate 90 $jpgfile\n";
  if($ieps==1){print CSH "convert $psfile $epsfile\n"}

  close (CSH); system("csh -f $cshfile");
  if($ixv==1) {system("xv $jpgfile &")}
}



#===========================================================================
#===========================================================================

if ($ifig02 == 1) {

  # shifting the subplots
  $xfac = 1.20;
  $yfac = 1.20;
  $dX = $xfac*$Jwid; $dY = 0; $shift1 = "-X$dX -Y$dY";
  $dX = -2*$dX; $dY = -$yfac*$Jwid; $shift2 = "-X$dX -Y$dY";

  #===============================================
  print "\nWriting CSH file...\n";

  open(CSH,">$cshfile");
  print CSH "gmtset BASEMAP_TYPE plain TICK_LENGTH $tick LABEL_FONT_SIZE $fsize2 ANOT_FONT_SIZE $fsize2 PLOT_DEGREE_FORMAT D HEADER_FONT $fontno ANOT_FONT $fontno LABEL_FONT $fontno HEADER_FONT_SIZE $fsize1 FRAME_PEN $fpen TICK_PEN $tpen\n";
  #===============================================

  $niter = 8;

  for ($k = 0; $k <= $niter; $k = $k+1) {
    $irun = $irun0 + 2*$k;              # wave2d run number
    $strun = sprintf("%4.4i",$irun);
    $dir = "${odir}$strun";
    $file1dat = "${dir}/$mfile_dat";
    $file1syn = "${dir}/$mfile_syn";
    $file1ker = "${dir}/summed_ker.dat";
    $ker_files[$k] = $file1ker;
    $mod_files[$k] = $file1syn;
    if (not -f $file1syn)   { die("Check if $file1syn exist or not\n") }

    # load chi values
    $chifile = "${dir}/summed_chi_all.dat";
    open(IN,"$chifile");
    $chi = <IN>;
    $it_vals[$k] = $k;
    $chi_vals[$k] = $chi;
    #$schi = sprintf("\@~\143\@~ ( $smod )  =  %3.3e",$chi);
  }

  #die("testing");

  # get the receivers and sources
  $recfile = "${dir0}/${edir}001/sr.txt";  # src-rec for first event
  $evefile = "${dir0}/events_lonlat.dat";
  if (not -f $recfile) { die("Check if $recfile exist or not\n") }
  if (not -f $evefile) { die("Check if $evefile exist or not\n") }

  #=============================================
  # model for the data

  $B = "$B0".$Bopts[15];

for ($k = 0; $k <= $niter; $k = $k+1) {

  # file names for figures
  $name = "talk_cg02_${k}"; $psfile  = "$name.ps"; $jpgfile = "$name.jpg"; $epsfile = "$name.eps";
  $schi = "\@~\143\@~ ( m\@+${k}\@+ )";

  chomp($it_vals[$k]); chomp($chi_vals[$k]);

  # phase speed map for data
  $title = "Phase speed for data";
  print CSH "psbasemap $B $R $J -K -V $origin > $psfile\n"; # START
  if($icolor==1) {print CSH "awk '{print \$1,\$2,\$4*100}' $file1dat | pscontour $R $J -A- -C$cpt_vel -I -K -O -V  >> $psfile\n"}
  print CSH "pscoast $J $R $coast_info -K -O -V >> $psfile\n";
  if($ishelf) {print CSH "awk '{print \$2,\$1}' $shelf_file |psxy $J $R $Wshelf -K -O -V >> $psfile\n"}
  print CSH "awk '\$1 == \"R\" {print \$2,\$3}' $recfile |psxy -N $J $R -K -O -V $rec >> $psfile\n";
  print CSH "awk '{print \$1,\$2}' $evefile |psxy -N $J $R -K -O -V $src >> $psfile\n";
  #print CSH "psscale -C$cpt_vel $Dscale $Bscale1 -K -O -V >> $psfile \n";
  print CSH "pstext -N $J_title $R_title -K -O -V >>$psfile<<EOF\n $x_title $z_title $fsize_title 0 $fontno CM $title \nEOF\n";

  # phase speed map for synthetics
  $mod = $mods[$k*2]; $smod = "m\@+$mod\@+"; $title = "Phase speed model $smod";
  print CSH "psbasemap $B $R $J -K -O -V $shift1 >> $psfile\n";
  if($icolor==1) {print CSH "awk '{print \$1,\$2,\$4*100}' $mod_files[$k] | pscontour $R $J -A- -C$cpt_vel -I -K -O -V  >> $psfile\n"}
  print CSH "pscoast $J $R $coast_info -K -O -V >> $psfile\n";
  if($ishelf) {print CSH "awk '{print \$2,\$1}' $shelf_file |psxy $J $R $Wshelf -K -O -V >> $psfile\n"}
  print CSH "awk '\$1 == \"R\" {print \$2,\$3}' $recfile |psxy -N $J $R -K -O -V $rec >> $psfile\n";
  print CSH "awk '{print \$1,\$2}' $evefile |psxy -N $J $R -K -O -V $src >> $psfile\n";
  #print CSH "psscale -C$cpt_vel $Dscale $Bscale1 -K -O -V >> $psfile \n";
  print CSH "pstext -N $J_title $R_title -K -O -V >>$psfile<<EOF\n $x_title $z_title $fsize_title 0 $fontno CM $title \nEOF\n";

  # misfit vs iteration
  $shift = $shift1;
  $Bg = "$B3b".$Bopts[8];
  $title = "Misfit,  \@~\143\@~ ( m )  ( $utype )";

  print CSH "psbasemap $Bg $R_chi $J_chi -K -O -V $shift >> $psfile\n";
  print CSH "awk '{print \$1,\$2}' $chi_curve | psxy $c_info_rd $J_chi $R_chi -K -O -V >> $psfile\n";  # red dot
  for ($j = 0; $j <= $niter; $j = $j+1) {
    chomp($it_vals[$j]); chomp($chi_vals[$j]);
    print CSH "psxy $J_chi $R_chi $p_info_k -K -O -V >>$psfile<<EOF\n $it_vals[$j] $chi_vals[$j]\nEOF\n";
  }
  $xpt = $it_vals[$k] + 0.5;
  print CSH "psxy $J_chi $R_chi $p_info_r -K -O -V >>$psfile<<EOF\n $it_vals[$k] $chi_vals[$k]\nEOF\n";
  print CSH "pstext -N $J_chi $R_chi -K -O -V >>$psfile<<EOF\n $xpt $chi_vals[$k] $fsize1 0 $fontno LM $schi \nEOF\n";
  print CSH "pstext -N $J_title $R_title -K -O -V >>$psfile<<EOF\n $x_title $z_title $fsize_title 0 $fontno CM $title \nEOF\n";

#-----------------------------
  print CSH "pstext -N $J_title $R_title -O -V >>$psfile<<EOF\n 10 10 $fsize_title 0 $fontno CM $title \nEOF\n"; # FINISH
  print CSH "convert $psfile -rotate 90 $jpgfile\n";
  if($ieps==1){print CSH "convert $psfile $epsfile\n"}

}

  close (CSH);
  system("csh -f $cshfile");
  system("xv $jpgfile &");

}

#===========================================================================
#===========================================================================

if ($ifig03 == 1) {

  # shifting the subplots
  $xfac = 1.20;
  $yfac = 1.20;
  $dX = $xfac*$Jwid; $dY = 0; $shift1 = "-X$dX -Y$dY";
  $dX = -2*$dX; $dY = -$yfac*$Jwid; $shift2 = "-X$dX -Y$dY";

  #===============================================
  print "\nWriting CSH file...\n";

  open(CSH,">$cshfile");
  print CSH "gmtset BASEMAP_TYPE plain TICK_LENGTH $tick LABEL_FONT_SIZE $fsize2 ANOT_FONT_SIZE $fsize2 PLOT_DEGREE_FORMAT D HEADER_FONT $fontno ANOT_FONT $fontno LABEL_FONT $fontno HEADER_FONT_SIZE $fsize1 FRAME_PEN $fpen TICK_PEN $tpen\n";
  #===============================================

  $dir = "/home/carltape/sem2d/2d_adjoint_banana/mat_SAVE/Gik_examples";
  #$dir  = "/home/carltape/sem2d/2d_adjoint_banana";
  $B = "$B0".$Bopts[15];

  # load index file into sample kernels and splines
  $index_file = "${dir}/sample_index.dat";
  if (not -f $index_file) { die("Check if $index_file exist or not\n") }
  open(IN,"$index_file"); @index = <IN>; $nfig = @index;

  # get the splines
  $spl_file = "${dir}/con_lonlat_q08.dat";
  if (not -f $spl_file) { die("Check if $spl_file exist or not\n") }

  # get the receivers and sources
  $recfile = "${dir}/recs_lonlat.dat";
  $evefile = "${dir}/events_lonlat.dat";
  if (not -f $recfile) { die("Check if $recfile exist or not\n") }
  if (not -f $evefile) { die("Check if $evefile exist or not\n") }
  open(IN,"$recfile"); @recs = <IN>; $nrec = @recs;
  open(IN,"$evefile"); @eves = <IN>; $neve = @eves;
  $nmeas = $nrec*$neve;

  print "\n $nfig figures showing Gik elements \n";
  print "\n $nrec receivers and $neve events = $nmeas paths \n";

  #$nfig = 1;                # use for testing one only
  $shift = $shift1;

  # color scaling factor
  @facs  = (0.5,0.5,0.5,0.15,0.5,0.5);
  $spl_fk = "-Sc3p -G0";

  #for ($k = $nfig-1; $k <= $nfig-1; $k = $k+1) {
  for ($k = 0; $k <= $nfig-1; $k = $k+1) {

    ($ieve,$irec,$iker,$ispl,$Gik) = split(" ",$index[$k]);
    $seve = sprintf("%4.4i",$ieve);
    $srec = sprintf("%4.4i",$irec);
    $sker = sprintf("%4.4i",$iker);
    $sspl = sprintf("%4.4i",$ispl);
    $sgik = sprintf("%.2f",$Gik);

    $fileG = "${dir}/Aij_e${seve}_r${srec}_k${sker}_s${sspl}";
    if (not -f $fileG) { die("Check if $fileG exist or not\n") }
    print "\n $k, $iker, $ispl, $Gik";

    # lat-lon for target source and receiver
    #$ieve = $ieves[$k];
    #$irec = $irecs[$k];
    ($slon,$slat,$junk) = split(" ",$eves[$ieve-1]);
    ($rlon,$rlat,$junk) = split(" ",$recs[$irec-1]);

    # lat-lon for target spline
    open(IN,"$spl_file");
    @splines = <IN>;
    ($blon,$blat) = split(" ",$splines[$ispl-1]);

    #----------------------------
    # make an enhanced color scale for the K*B plots

    $fac = $facs[$k];
    $ss = $kmax*$fac;
    $ds = 2*$ss/$scale_color;
    $T2b = sprintf("-T%3.3e/%3.3e/%3.3e",-$ss,$ss,$ds);
    print "\nT2b = $T2b\n";
    $bs2b = sprintf("%2.2e",0.9*$ss);
    $Bscale2c  = sprintf("-B${bs2b}:\" K * B  ( 10\@+%2.2i\@+ ) \": -E10p",$kpwr);

    $cpt_kerb = "color1b.cpt";
    print CSH "makecpt -C$colorbar $T2b > temp1\n";
    print CSH "sed 's/^B.*/B       170     0      0  /' temp1 >  temp2\n";
    print CSH "sed 's/^F.*/F         0     0    255  /' temp2 > $cpt_kerb\n";
    #----------------------------

    # file names for figures
    $name = "talk_Gik03_${k}"; $psfile  = "$name.ps"; $jpgfile = "$name.jpg"; $epsfile = "$name.eps";

  # kernel
  $title = "Kernel ${iker} / ${nmeas} (event ${ieve}, rec ${irec})";
  print CSH "psbasemap $B $R $J -K -V $origin > $psfile\n"; # START
  if($icolor==1) {print CSH "awk '{print \$1,\$2,\$3 / $norm2 }' $fileG | pscontour $R $J -A- -C$cpt_ker -I -K -O -V  >> $psfile\n"}
  print CSH "pscoast $J $R $coast_info -K -O -V >> $psfile\n";
  if($ishelf) {print CSH "awk '{print \$2,\$1}' $shelf_file |psxy $J $R $Wshelf -K -O -V >> $psfile\n"}
  print CSH "psxy $J $R -K -O -V $src1 >>$psfile<<EOF\n $slon $slat \nEOF\n";
  print CSH "psxy $J $R -K -O -V $rec1 >>$psfile<<EOF\n $rlon $rlat \nEOF\n";
  print CSH "psxy $J $R -K -O -V $rec_fw >>$psfile<<EOF\n $blon $blat \nEOF\n";
  #print CSH "awk '{print \$1,\$2}' $recfile |psxy -N $J $R -K -O -V $rec >> $psfile\n";
  #print CSH "awk '{print \$1,\$2}' $evefile |psxy -N $J $R -K -O -V $src >> $psfile\n";
  print CSH "psscale -C$cpt_ker $Dscale $Bscale2 -K -O -V >> $psfile \n";
  print CSH "pstext -N $J_title $R_title -K -O -V >>$psfile<<EOF\n $x_title $z_title $fsize_title 0 $fontno CM $title \nEOF\n";

  # spline
  $title = "Spline ${ispl}";
  print CSH "psbasemap $B $R $J -K -O -V $shift >> $psfile\n";
  if($icolor==1) {print CSH "awk '{print \$1,\$2,\$4 }' $fileG | pscontour $R $J -A- -C$cpt_spl -I -K -O -V >> $psfile\n"}
  print CSH "pscoast $J $R $coast_info -K -O -V >> $psfile\n";
  if($ishelf) {print CSH "awk '{print \$2,\$1}' $shelf_file |psxy $J $R $Wshelf -K -O -V >> $psfile\n"}
  print CSH "psxy $J $R -K -O -V $src1 >>$psfile<<EOF\n $slon $slat \nEOF\n";
  print CSH "psxy $J $R -K -O -V $rec1 >>$psfile<<EOF\n $rlon $rlat \nEOF\n";
  print CSH "awk '{print \$1,\$2}' $spl_file |psxy -N $J $R -K -O -V $spl_fk >> $psfile\n";
  print CSH "psxy $J $R -K -O -V $rec_fw >>$psfile<<EOF\n $blon $blat \nEOF\n";
  #print CSH "awk '{print \$1,\$2}' $recfile |psxy -N $J $R -K -O -V $rec >> $psfile\n";
  #print CSH "awk '{print \$1,\$2}' $evefile |psxy -N $J $R -K -O -V $src >> $psfile\n";
  print CSH "psscale -C$cpt_spl $Dscale $Bscale_spl -K -O -V >> $psfile \n";
  print CSH "pstext -N $J_title $R_title -K -O -V >>$psfile<<EOF\n $x_title $z_title $fsize_title 0 $fontno CM $title \nEOF\n";

  # kernel * spline
  $title = "Kernel * Spline ( ${sgik} )";
  print CSH "psbasemap $B $R $J -K -O -V $shift >> $psfile\n";
  if($icolor==1) {print CSH "awk '{print \$1,\$2,\$5 / $norm2 }' $fileG | pscontour $R $J -A- -C$cpt_kerb -I -K -O -V  >> $psfile\n"}
  print CSH "pscoast $J $R $coast_info -K -O -V >> $psfile\n";
  if($ishelf) {print CSH "awk '{print \$2,\$1}' $shelf_file |psxy $J $R $Wshelf -K -O -V >> $psfile\n"}
  print CSH "psxy $J $R -K -O -V $src1 >>$psfile<<EOF\n $slon $slat \nEOF\n";
  print CSH "psxy $J $R -K -O -V $rec1 >>$psfile<<EOF\n $rlon $rlat \nEOF\n";
  print CSH "psxy $J $R -K -O -V $rec_fw >>$psfile<<EOF\n $blon $blat \nEOF\n";
  #print CSH "awk '{print \$1,\$2}' $recfile |psxy -N $J $R -K -O -V $rec >> $psfile\n";
  #print CSH "awk '{print \$1,\$2}' $evefile |psxy -N $J $R -K -O -V $src >> $psfile\n";
  print CSH "psscale -C$cpt_kerb $Dscale $Bscale2c -K -O -V >> $psfile \n";
  print CSH "pstext -N $J_title $R_title -K -O -V >>$psfile<<EOF\n $x_title $z_title $fsize_title 0 $fontno CM $title \nEOF\n";

#-----------------------------
  print CSH "pstext -N $J_title $R_title -O -V >>$psfile<<EOF\n 10 10 $fsize_title 0 $fontno CM $title \nEOF\n"; # FINISH
  print CSH "convert $psfile -rotate 90 $jpgfile\n";
  if($ieps==1){print CSH "convert $psfile $epsfile\n"}

  }

  close (CSH);
  system("csh -f $cshfile");
  system("xv $jpgfile &");

}

#===========================================================================
if ($ifig03p == 1) {

  $origin = "-X0.5 -Y4.75";
  $Jwid = 2.25;
  $J = "-JM${Jwid}i";      # in lat-lon
  $J_title = "-JX${Jwid}";

  # shifting the subplots
  $xfac = 1.20;
  $yfac = 1.20;
  $dX = $xfac*$Jwid; $dY = 0; $shift1 = "-X$dX -Y$dY";
  $dX = -2*$dX; $dY = 1.2*$yfac*$Jwid; $shift2 = "-X$dX -Y$dY";

  # colorbar
  $Dlen = 0.75*$Jwid; $Dx = $Jwid/2; $Dy = -0.08*$Jwid;
  $Dscale = "-D$Dx/$Dy/$Dlen/0.12h";

  $fontno = 4;
  $fsize_title = 9;

  #===============================================
  print "\nWriting CSH file...\n";

  open(CSH,">$cshfile");
  print CSH "gmtset BASEMAP_TYPE plain TICK_LENGTH $tick LABEL_FONT_SIZE 9 ANOT_FONT_SIZE 9 PLOT_DEGREE_FORMAT D HEADER_FONT $fontno ANOT_FONT $fontno LABEL_FONT $fontno HEADER_FONT_SIZE 12 FRAME_PEN $fpen TICK_PEN $tpen\n";
  #===============================================

  $dir = "/home/carltape/sem2d/2d_adjoint_banana/mat_SAVE/Gik_examples";
  #$dir  = "/home/carltape/sem2d/2d_adjoint_banana";

  # load index file into sample kernels and splines
  $index_file = "${dir}/sample_index.dat";
  if (not -f $index_file) { die("Check if $index_file exist or not\n") }
  open(IN,"$index_file"); @index = <IN>; $nfig = @index;
  print "\n index of Gik examples to plot for kernels \n @index \n";

  # get the splines
  $spl_file = "${dir}/con_lonlat_q08.dat";
  if (not -f $spl_file) { die("Check if $spl_file exist or not\n") }
  $ngrid = 286;

  # get the receivers and sources
  $recfile = "${dir}/recs_lonlat.dat";
  $evefile = "${dir}/events_lonlat.dat";
  if (not -f $recfile) { die("Check if $recfile exist or not\n") }
  if (not -f $evefile) { die("Check if $evefile exist or not\n") }
  open(IN,"$recfile"); @recs = <IN>; $nrec = @recs;
  open(IN,"$evefile"); @eves = <IN>; $neve = @eves;
  $nmeas = $nrec*$neve;

  print "\n $nfig figures showing Gik elements \n";
  print "\n $nrec receivers and $neve events = $nmeas paths \n";

  #$nfig = 1;                # use for testing one only
  $shift = $shift1;

  # color scaling factor
  @facs  = (0.5,0.5,0.5,0.15,0.5,0.5);
  $spl_fk = "-Sc3p -G0";

  # KEY COMMAND: pick the spline from the index file
  $k = $nfig-1;

    ($ieve,$irec,$iker,$ispl,$Gik) = split(" ",$index[$k]);
    $seve = sprintf("%4.4i",$ieve);
    $srec = sprintf("%4.4i",$irec);
    $sker = sprintf("%4.4i",$iker);
    $sspl = sprintf("%4.4i",$ispl);
    $sgik = sprintf("%.2f",$Gik);

    $fileG = "${dir}/Aij_ker_e${seve}_r${srec}_k${sker}_s${sspl}";
    if (not -f $fileG) { die("Check if $fileG exist or not\n") }
    print "\n $k, $iker, $ispl, $Gik";

    # lat-lon for target source and receiver
    #$ieve = $ieves[$k];
    #$irec = $irecs[$k];
    ($slon,$slat,$junk) = split(" ",$eves[$ieve-1]);
    ($rlon,$rlat,$junk) = split(" ",$recs[$irec-1]);

    # lat-lon for target spline
    open(IN,"$spl_file");
    @splines = <IN>;
    ($blon,$blat) = split(" ",$splines[$ispl-1]);

    #----------------------------
    # make an enhanced color scale for the K*B plots

    $fac = $facs[$k];
    $ss = $kmax*$fac;
    $ds = 2*$ss/$scale_color;
    $T2b = sprintf("-T%3.3e/%3.3e/%3.3e",-$ss,$ss,$ds);
    print "\nT2b = $T2b\n";
    $bs2b = sprintf("%2.2e",0.9*$ss);
    $Bscale2c  = sprintf("-B${bs2b}:\" K * B  ( 10\@+%2.2i\@+  m\@+-2\@+ s ) \": -E10p",$kpwr);

    $cpt_kerb = "color1b.cpt";
    print CSH "makecpt -C$colorbar $T2b > temp1\n";
    print CSH "sed 's/^B.*/B       170     0      0  /' temp1 >  temp2\n";
    print CSH "sed 's/^F.*/F         0     0    255  /' temp2 > $cpt_kerb\n";
    #----------------------------

    # file names for figures
    $name = "talk_Gik03_${k}_ker_paper"; $psfile  = "$name.ps"; $jpgfile = "$name.jpg"; $epsfile = "$name.eps";

  @labs = ("d","e","f","a","b","c");
  @Baxes = (15,15,15,7,4,99);

  # kernel
  $title = "($labs[0])  Kernel : i = ${iker} / ${nmeas} (event ${ieve}, rec ${irec})";
  $B = "$B0".$Bopts[$Baxes[0]];
  print CSH "psbasemap $B $R $J -K -V -P $origin > $psfile\n"; # START
  if($icolor==1) {print CSH "awk '{print \$1,\$2,\$3 / $norm2 }' $fileG | pscontour $R $J -A- -C$cpt_ker -I -K -O -V -P  >> $psfile\n"}
  print CSH "pscoast $J $R $coast_info -K -O -V -P >> $psfile\n";
  if($ishelf) {print CSH "awk '{print \$2,\$1}' $shelf_file |psxy $J $R $Wshelf -K -O -V -P >> $psfile\n"}
  print CSH "psxy $J $R -K -O -V -P $src1 >>$psfile<<EOF\n $slon $slat \nEOF\n";
  print CSH "psxy $J $R -K -O -V -P $rec1 >>$psfile<<EOF\n $rlon $rlat \nEOF\n";
  print CSH "psxy $J $R -K -O -V -P $rec_fw >>$psfile<<EOF\n $blon $blat \nEOF\n";
  #print CSH "awk '{print \$1,\$2}' $recfile |psxy -N $J $R -K -O -V -P $rec >> $psfile\n";
  #print CSH "awk '{print \$1,\$2}' $evefile |psxy -N $J $R -K -O -V -P $src >> $psfile\n";
  print CSH "psscale -C$cpt_ker $Dscale $Bscale2b -K -O -V -P >> $psfile \n";
  print CSH "pstext -N $J_title $R_title -K -O -V -P >>$psfile<<EOF\n $x_title $z_title $fsize_title 0 $fontno CM $title \nEOF\n";

  # spline
  $title = "($labs[1])  Spline : k = ${ispl} / $ngrid";
  $B = "$B0".$Bopts[$Baxes[1]];
  print CSH "psbasemap $B $R $J -K -O -V -P $shift >> $psfile\n";
  if($icolor==1) {print CSH "awk '{print \$1,\$2,\$4 }' $fileG | pscontour $R $J -A- -C$cpt_spl -I -K -O -V -P >> $psfile\n"}
  print CSH "pscoast $J $R $coast_info -K -O -V -P >> $psfile\n";
  if($ishelf) {print CSH "awk '{print \$2,\$1}' $shelf_file |psxy $J $R $Wshelf -K -O -V -P >> $psfile\n"}
  print CSH "psxy $J $R -K -O -V -P $src1 >>$psfile<<EOF\n $slon $slat \nEOF\n";
  print CSH "psxy $J $R -K -O -V -P $rec1 >>$psfile<<EOF\n $rlon $rlat \nEOF\n";
  print CSH "awk '{print \$1,\$2}' $spl_file |psxy -N $J $R -K -O -V -P $spl_fk >> $psfile\n";
  print CSH "psxy $J $R -K -O -V -P $rec_fw >>$psfile<<EOF\n $blon $blat \nEOF\n";
  print CSH "psscale -C$cpt_spl $Dscale $Bscale_spl -K -O -V -P >> $psfile \n";
  print CSH "pstext -N $J_title $R_title -K -O -V -P >>$psfile<<EOF\n $x_title $z_title $fsize_title 0 $fontno CM $title \nEOF\n";

  # kernel * spline
  $title = "($labs[2])  Kernel * Spline ( G\@-ik\@- = ${sgik} )";
  $B = "$B0".$Bopts[$Baxes[2]];
  print CSH "psbasemap $B $R $J -K -O -V -P $shift >> $psfile\n";
  if($icolor==1) {print CSH "awk '{print \$1,\$2,\$5 / $norm2 }' $fileG | pscontour $R $J -A- -C$cpt_kerb -I -K -O -V -P  >> $psfile\n"}
  print CSH "pscoast $J $R $coast_info -K -O -V -P >> $psfile\n";
  if($ishelf) {print CSH "awk '{print \$2,\$1}' $shelf_file |psxy $J $R $Wshelf -K -O -V -P >> $psfile\n"}
  print CSH "psxy $J $R -K -O -V -P $src1 >>$psfile<<EOF\n $slon $slat \nEOF\n";
  print CSH "psxy $J $R -K -O -V -P $rec1 >>$psfile<<EOF\n $rlon $rlat \nEOF\n";
  print CSH "psxy $J $R -K -O -V -P $rec_fw >>$psfile<<EOF\n $blon $blat \nEOF\n";
  print CSH "psscale -C$cpt_kerb $Dscale $Bscale2c -K -O -V -P >> $psfile \n";
  print CSH "pstext -N $J_title $R_title -K -O -V -P >>$psfile<<EOF\n $x_title $z_title $fsize_title 0 $fontno CM $title \nEOF\n";

  #-----------------------------
  #print CSH "pstext -N $J_title $R_title -O -V -P >>$psfile<<EOF\n 10 10 $fsize_title 0 $fontno CM $title \nEOF\n"; # FINISH
  #print CSH "convert $psfile $jpgfile\n";
  #if($ieps==1){print CSH "convert $psfile $epsfile\n"}

  # KERNEL
  #================================
  # RAY PATH

  # file names for figures
    #$name = "talk_Gik03_${k}_ray_paper"; $psfile  = "$name.ps"; $jpgfile = "$name.jpg"; $epsfile = "$name.eps";

   $iextra = 0;
  ($ieve,$irec,$iray,$ispl,$ispl2) = (1,126,126,203,111);

   $seve = sprintf("%4.4i",$ieve);
   $srec = sprintf("%4.4i",$irec);
   $sker = sprintf("%4.4i",$iray);
   $sspl = sprintf("%4.4i",$ispl);
   $sspl2 = sprintf("%4.4i",$ispl2);

  $ray_dir = "${dir}/rays";
  $ray_file = "${ray_dir}/ray_path.dat";
  $profile_file = "${ray_dir}/gaussian.dat";
  $profile_file2 = "${ray_dir}/gaussian_extra.dat";
  $int_file = "${ray_dir}/gaussian_int.dat";
  $basis_file = "${ray_dir}/Aij_ray_e${seve}_r${srec}_k${sker}_s${sspl}.dat";
  if (not -f $ray_file) { die("Check if $ray_file exist or not\n") }
  if (not -f $profile_file) { die("Check if $profile_file exist or not\n") }
  if (not -f $profile_file2) { die("Check if $profile_file2 exist or not\n") }
  if (not -f $int_file) { die("Check if $int_file exist or not\n") }
  if (not -f $basis_file) { die("Check if $basis_file exist or not\n") }

  # Gik value and integral of Gaussian
  open(IN,"$int_file"); @lims = <IN>;
  ($int,$Gik) = split(" ",$lims[0]);
  #print "\n $int, $Gik \n\n"; die("testing");
  $sgik = sprintf("%.2f",$Gik);

  $src_col = "-G255/0/0";
  $rec_col = "-G0/255/255";

  # lat-lon for target spline
  ($blon,$blat) = split(" ",$splines[$ispl-1]);
  ($blon2,$blat2) = split(" ",$splines[$ispl2-1]);

  # ray path
  $title = "($labs[3])  Ray : i = ${iray} / ${nmeas} (event ${ieve}, rec ${irec})";
  $B = "$B0".$Bopts[$Baxes[3]];
  #print CSH "psbasemap $B $R $J -K -V -P $origin > $psfile\n"; # START
  print CSH "psbasemap $B $R $J -K -O -V -P $shift2 >> $psfile\n";
  print CSH "psxy $ray_file $J $R -K -O -V -P -W1.0p >>$psfile\n";
  print CSH "pscoast $J $R $coast_info -K -O -V -P >> $psfile\n";
  if($ishelf) {print CSH "awk '{print \$2,\$1}' $shelf_file |psxy $J $R $Wshelf -K -O -V -P >> $psfile\n"}
  print CSH "psxy $J $R -K -O -V -P -G255/0/0 $src1 $src_col >>$psfile<<EOF\n $slon $slat \nEOF\n";
  print CSH "psxy $J $R -K -O -V -P $rec1 $rec_col >>$psfile<<EOF\n $rlon $rlat \nEOF\n";
  print CSH "psxy $J $R -K -O -V -P $rec_fw >>$psfile<<EOF\n $blon $blat \nEOF\n";
  if($iextra==1) {print CSH "psxy $J $R -K -O -V -P $rec_fw >>$psfile<<EOF\n $blon2 $blat2 \nEOF\n"}
  #print CSH "psscale -C$cpt_ker $Dscale $Bscale2 -K -O -V -P >> $psfile \n";
  print CSH "pstext -N $J_title $R_title -K -O -V -P >>$psfile<<EOF\n $x_title $z_title $fsize_title 0 $fontno CM $title \nEOF\n";

  # spline
  $title = "($labs[4])  Spline : k = ${ispl} / $ngrid";
  $B = "$B0".$Bopts[$Baxes[4]];
  print CSH "psbasemap $B $R $J -K -O -V -P $shift >> $psfile\n";
  if($icolor==1) {print CSH "pscontour $basis_file $R $J -A- -C$cpt_spl -I -K -O -V -P >> $psfile\n"}
  print CSH "psxy $ray_file $J $R -K -O -V -P -W1.0p >>$psfile\n";
  print CSH "pscoast $J $R $coast_info -K -O -V -P >> $psfile\n";
  if($ishelf) {print CSH "awk '{print \$2,\$1}' $shelf_file |psxy $J $R $Wshelf -K -O -V -P >> $psfile\n"}
  print CSH "psxy $J $R -K -O -V -P $src1 $src_col >>$psfile<<EOF\n $slon $slat \nEOF\n";
  print CSH "psxy $J $R -K -O -V -P $rec1 $rec_col >>$psfile<<EOF\n $rlon $rlat \nEOF\n";
  print CSH "awk '{print \$1,\$2}' $spl_file |psxy -N $J $R -K -O -V -P $spl_fk >> $psfile\n";
  print CSH "psxy $J $R -K -O -V -P $rec_fw >>$psfile<<EOF\n $blon $blat \nEOF\n";
  if($iextra==1) {print CSH "psxy $J $R -K -O -V -P $rec_fw >>$psfile<<EOF\n $blon2 $blat2 \nEOF\n"}
  #print CSH "psscale -C$cpt_spl $Dscale $Bscale_spl -K -O -V -P >> $psfile \n";
  print CSH "pstext -N $J_title $R_title -K -O -V -P >>$psfile<<EOF\n $x_title $z_title $fsize_title 0 $fontno CM $title \nEOF\n";

   # axes limits for gaussian plot
  $axes_file = "${ray_dir}/axes.dat";
  if (not -f $axes_file)   { die("Check if $axes_file exist or not\n") }
  open(IN,"$axes_file");
  @ax_lims = <IN>;
  ($a1x1,$a1x2,$a1y1,$a1y2) = split(" ",$ax_lims[0]);  # inner box
  ($a2x1,$a2x2,$a2y1,$a2y2) = split(" ",$ax_lims[1]);  # outer box

  # find the power of the distance-along-ray
  $lpwr_x = int( log($a2x2) / log(10) );
  $normx = "1e$lpwr_x";

  $xmin1 = $a1x1 / $normx;
  $xmax1 = $a1x2 / $normx;
  $xmin2 = $a2x1 / $normx;
  $xmax2 = $a2x2 / $normx;

  $spl_tick = 0.2;
  $Jgaus = "-JX${Jwid}";
  $Rgaus = "-R$xmin2/$xmax2/$a2y1/$a2y2";
  $Bgaus = sprintf("-B1:\" Distance along ray path  (10\@+%2.2i\@+ m) \":/${spl_tick}".$Bopts[4],$lpwr_x);

  $Bscale_spl  = "-B${spl_tick} -A";
  $Dlen = 0.84 * $Jwid; $Dx = 0; $Dy = $Jwid/2;
  $Dscale = "-D$Dx/$Dy/$Dlen/0.12";

  print "\n $a2x2 $lpwr_x $normx \n $Jgaus $Rgaus \n\n";
  #die("testing");

  # ray * spline
  $title = "($labs[5])  Spline along ray ( G\@-ik\@- = ${sgik} )";
  print CSH "psbasemap $Bgaus $Rgaus $Jgaus -K -O -V -P $shift >> $psfile\n";

  print CSH "awk '{print \$1 / $normx,\$2}' $profile_file | psxy $Jgaus $Rgaus -K -O -V -P -W1.0p >>$psfile\n";
  if($iextra==1) {print CSH "awk '{print \$1 / $normx,\$2}' $profile_file2 | psxy $Jgaus $Rgaus -K -O -V -P -W1.0tap >>$psfile\n"}
  print CSH "psscale -C$cpt_spl2 $Dscale $Bscale_spl -K -O -V -P >> $psfile \n";
  print CSH "psxy $Jgaus $Rgaus -K -O -V -P $src1 $src_col >>$psfile<<EOF\n $xmin1 0 \nEOF\n";
  print CSH "psxy $Jgaus $Rgaus -K -O -V -P $rec1 $rec_col >>$psfile<<EOF\n $xmax1 0 \nEOF\n";
  print CSH "pstext -N $J_title $R_title -K -O -V -P >>$psfile<<EOF\n $x_title $z_title $fsize_title 0 $fontno CM $title \nEOF\n";

#-----------------------------
  print CSH "pstext -N $J_title $R_title -O -V -P >>$psfile<<EOF\n 10 10 $fsize_title 0 $fontno CM $title \nEOF\n"; # FINISH
  print CSH "convert $psfile $jpgfile\n";
  if($ieps==1){print CSH "convert $psfile $epsfile\n"}

  close (CSH);
  system("csh -f $cshfile");
  system("xv $jpgfile &");
}

#===========================================================================
if ($ifig04 == 1) {

  $Jwid = 3;
  $J = "-JM${Jwid}i";      # in lat-lon
  $origin = "-X0.75 -Y3.0";

  # colorbar
  $Dlen = 0.75*${Jwid}; $Dx = $Jwid/2; $Dy = -0.10*$Jwid;
  $Dscale = "-D$Dx/$Dy/$Dlen/0.20h";

  # shifting the subplots
  $xfac = 1.15;
  $yfac = 1.35;
  $dX = $xfac*$Jwid; $dY = 0; $shift1 = "-X$dX -Y$dY";
  $dX = -2*$dX; $dY = -$yfac*$Jwid; $shift2 = "-X$dX -Y$dY";

  #===============================================
  print "\nWriting CSH file...\n";

  open(CSH,">$cshfile");
  print CSH "gmtset BASEMAP_TYPE plain TICK_LENGTH $tick LABEL_FONT_SIZE $fsize2 ANOT_FONT_SIZE $fsize2 PLOT_DEGREE_FORMAT D HEADER_FONT $fontno ANOT_FONT $fontno LABEL_FONT $fontno HEADER_FONT_SIZE $fsize1 FRAME_PEN $fpen TICK_PEN $tpen\n";
  #===============================================

  if($irun0 == 0) { $iev = 5 } else { $iev = 1 }
  $stev = sprintf("%3.3i",$iev);

  # get the summed kernel and chi files
  $file2 = "${dir}/summed_ker.dat";
  #$file2 = "${dir}/fun_smooth.dat";       # smoothed kernel
  $recfile = "${dir}/${edir}${stev}/sr.txt";  # src-rec for event 001

  if (not -f $file2)   { die("Check if $file2 exist or not\n") }
  if (not -f $recfile) { die("Check if $recfile exist or not\n") }

  # number of receivers
  open(IN,$recfile); @temp = <IN>; $nrec = @temp - 1;

  # KEY: make TWO figures, one having the event kernel
  for ($k = 0; $k <= 1; $k = $k+1) {

  if($k == 0) { $slab = "mod" } else { $slab = "mod_ker" }

  # file names for figures
  $name = "talk_${slab}_${strun0}";
  $psfile = "$name.ps"; $jpgfile = "$name.jpg"; $epsfile = "$name.eps";

  #=============================================
  # model for the data

  $B = "$B0".$Bopts[15];
  $title = "Phase speed for data";

  # phase speed map
  print CSH "psbasemap $B $R $J -K -V $origin > $psfile\n"; # START
  #if($icolor==1) {print CSH "awk '{print \$1,\$2,\$3/1000}' $file1dat | pscontour $R $J -A- -C$cpt_vel -I -K -O -V >> $psfile\n"}
  if($icolor==1) {print CSH "awk '{print \$1,\$2,\$4*100}' $file1dat | pscontour $R $J -A- -C$cpt_vel -I -K -O -V  >> $psfile\n"}
  print CSH "pscoast $J $R $coast_info -K -O -V >> $psfile\n";
  if($ishelf) {print CSH "awk '{print \$2,\$1}' $shelf_file |psxy $J $R $Wshelf -K -O -V >> $psfile\n"}
  print CSH "awk '\$1 == \"R\" {print \$2,\$3}' $recfile |psxy -N $J $R -K -O -V $rec_fw >> $psfile\n";
  print CSH "awk '\$1 == \"S\" {print \$2,\$3}' $recfile |psxy -N $J $R -K -O -V $src_fr >> $psfile\n";
  #print CSH "awk '{print \$1,\$2}' $evefile |psxy -N $J $R -K -O -V $src >> $psfile\n";
  print CSH "psscale -C$cpt_vel $Dscale $Bscale1 -K -O -V >> $psfile \n";
  print CSH "pstext -N $J_title $R_title -K -O -V >>$psfile<<EOF\n $x_title $z_title $fsize_title 0 $fontno CM $title \nEOF\n";

  #=============================================
  # model for synthetics

  $title = "Phase speed for synthetics";
  $shift = $shift1;

  # phase speed map
  print CSH "psbasemap $B $R $J -K -O -V $shift >> $psfile\n";
  #if($icolor==1) {print CSH "awk '{print \$1,\$2,\$3/1000}' $file1syn | pscontour $R $J -A- -C$cpt_vel -I -O -K -V >> $psfile\n"}
  if($icolor==1) {print CSH "awk '{print \$1,\$2,\$4*100}' $file1syn | pscontour $R $J -A- -C$cpt_vel -I -O -K -V >> $psfile\n"}
  print CSH "pscoast $J $R $coast_info -K -O -V >> $psfile\n";
  if($ishelf) {print CSH "awk '{print \$2,\$1}' $shelf_file |psxy $J $R $Wshelf -K -O -V >> $psfile\n"}
  #print CSH "awk '{print \$1,\$2}' $evefile |psxy -N $J $R -K -O -V $src >> $psfile\n";
  print CSH "awk '\$1 == \"R\" {print \$2,\$3}' $recfile |psxy -N $J $R -K -O -V $rec_fw >> $psfile\n";
  print CSH "awk '\$1 == \"S\" {print \$2,\$3}' $recfile |psxy -N $J $R -K -O -V $src_fr >> $psfile\n";
  print CSH "psscale -C$cpt_vel $Dscale $Bscale1 -K -O -V >> $psfile \n";
  print CSH "pstext -N $J_title $R_title -K -O -V >>$psfile<<EOF\n $x_title $z_title $fsize_title 0 $fontno CM $title \nEOF\n";

  #=============================================
  # kernel -- UNSMOOTHED VERSION

  if ($k == 1) {

  $shift = $shift1;
  $title = "Event kernel ($nrec receivers)";
  $B = "$B0".$Bopts[15];

  print CSH "psbasemap $B $R $J -K -O -V $shift >> $psfile\n";

  if($icolor==1) {print CSH "awk '{print \$1,\$2,\$3 / $norm2 }' $file2 | pscontour $J $R -A- -C$cpt_ker -I -K -O -V >> $psfile\n"}
  #if($icolor==1) {print CSH "awk '{print \$1,\$2,\$7 / $norm2 }' $file2 | pscontour $J $R -A- -C$cpt_ker -I -K -O -V >> $psfile\n"}

  print CSH "pscoast $J $R $coast_info -K -O -V >> $psfile\n";
  if($ishelf) {print CSH "awk '{print \$2,\$1}' $shelf_file |psxy $J $R $Wshelf -K -O -V >> $psfile\n"}
  print CSH "psscale -C$cpt_ker $Dscale $Bscale2 -K -O -V >> $psfile \n";
  #print CSH "awk '{print \$1,\$2}' $evefile |psxy -N $J $R -K -O -V $src >> $psfile\n";
  print CSH "awk '\$1 == \"R\" {print \$2,\$3}' $recfile |psxy -N $J $R -K -O -V $rec >> $psfile\n";
  print CSH "awk '\$1 == \"S\" {print \$2,\$3}' $recfile |psxy -N $J $R -K -O -V $src_fr >> $psfile\n";
  print CSH "pstext -N $J_title $R_title -K -O -V >>$psfile<<EOF\n $x_title $z_title $fsize_title 0 $fontno CM $title \nEOF\n";

}

#-----------------------------
  print CSH "pstext -N $J_title $R_title -O -V >>$psfile<<EOF\n 10 10 $fsize_title 0 $fontno CM $title \nEOF\n"; # FINISH
  print CSH "convert $psfile -rotate 90 $jpgfile\n";
  if($ieps==1){print CSH "convert $psfile $epsfile\n"}

  }  # $k

  close (CSH); system("csh -f $cshfile");
  system("xv $jpgfile &");

}


#===========================================================================
if ($ifig04p == 1) {

  $Jwid = 3;
  $J = "-JM${Jwid}i";      # in lat-lon
  $origin = "-X0.75 -Y3.0";
  $fontno = 4;

  # colorbar
  $Dlen = 0.75*${Jwid}; $Dx = $Jwid/2; $Dy = -0.08*$Jwid;
  $Dscale = "-D$Dx/$Dy/$Dlen/0.12h";

  # shifting the subplots
  $xfac = 1.15;
  $yfac = 1.35;
  $dX = $xfac*$Jwid; $dY = 0; $shift1 = "-X$dX -Y$dY";
  $dX = -2*$dX; $dY = -$yfac*$Jwid; $shift2 = "-X$dX -Y$dY";

  #===============================================
  print "\nWriting CSH file...\n";

  open(CSH,">$cshfile");
  print CSH "gmtset BASEMAP_TYPE plain TICK_LENGTH $tick LABEL_FONT_SIZE $fsize2 ANOT_FONT_SIZE $fsize2 PLOT_DEGREE_FORMAT D HEADER_FONT $fontno ANOT_FONT $fontno LABEL_FONT $fontno HEADER_FONT_SIZE $fsize1 FRAME_PEN $fpen TICK_PEN $tpen\n";
  #===============================================

  if($irun0 == 0) { $iev = 5 } else { $iev = 1 }
  $stev = sprintf("%3.3i",$iev);

  # get the summed kernel and chi files
  $file2 = "${dir}/summed_ker.dat";
  #$file2 = "${dir}/fun_smooth.dat";       # smoothed kernel
  $recfile = "${dir}/${edir}${stev}/sr.txt";  # src-rec for event 001

  if (not -f $file2)   { die("Check if $file2 exist or not\n") }
  if (not -f $recfile) { die("Check if $recfile exist or not\n") }

  # number of receivers
  open(IN,$recfile); @temp = <IN>; $nrec = @temp - 1;

  $k = 0;
  if($k == 0) { $slab = "mod" } else { $slab = "mod_ker" }

  # file names for figures
  $name = "talk_${slab}_${strun0}_paper";
  $psfile = "$name.ps"; $jpgfile = "$name.jpg"; $epsfile = "$name.eps";

  #=============================================
  # model for the data

  $B = "$B0".$Bopts[15];
  $title = "Target phase speed model (data)";

  # phase speed map
  print CSH "psbasemap $B $R $J -K -P -V $origin > $psfile\n"; # START
  #if($icolor==1) {print CSH "awk '{print \$1,\$2,\$3/1000}' $file1dat | pscontour $R $J -A- -C$cpt_vel -I -K -O -P -V >> $psfile\n"}
  if($icolor==1) {print CSH "awk '{print \$1,\$2,\$4*100}' $file1dat | pscontour $R $J -A- -C$cpt_vel -I -K -O -P -V  >> $psfile\n"}
  print CSH "pscoast $J $R $coast_info -K -O -P -V >> $psfile\n";
  if($ishelf) {print CSH "awk '{print \$2,\$1}' $shelf_file |psxy $J $R $Wshelf -K -O -P -V >> $psfile\n"}
  print CSH "awk '\$1 == \"R\" {print \$2,\$3}' $recfile |psxy -N $J $R -K -O -P -V $rec_fw >> $psfile\n";
  print CSH "awk '\$1 == \"S\" {print \$2,\$3}' $recfile |psxy -N $J $R -K -O -P -V $src_fr >> $psfile\n";
  #print CSH "awk '{print \$1,\$2}' $evefile |psxy -N $J $R -K -O -P -V $src >> $psfile\n";
  print CSH "psscale -C$cpt_vel $Dscale $Bscale1 -K -O -P -V >> $psfile \n";
  print CSH "pstext -N $J_title $R_title -K -O -P -V >>$psfile<<EOF\n $x_title $z_title $fsize2 0 $fontno CM $title \nEOF\n";

  #=============================================
  # model for synthetics

  $title = "Initial phase speed model (synthetics)";
  $shift = $shift1;

  # phase speed map
  print CSH "psbasemap $B $R $J -K -O -P -V $shift >> $psfile\n";
  #if($icolor==1) {print CSH "awk '{print \$1,\$2,\$3/1000}' $file1syn | pscontour $R $J -A- -C$cpt_vel -I -O -K -P -V >> $psfile\n"}
  if($icolor==1) {print CSH "awk '{print \$1,\$2,\$4*100}' $file1syn | pscontour $R $J -A- -C$cpt_vel -I -O -K -P -V >> $psfile\n"}
  print CSH "pscoast $J $R $coast_info -K -O -P -V >> $psfile\n";
  if($ishelf) {print CSH "awk '{print \$2,\$1}' $shelf_file |psxy $J $R $Wshelf -K -O -P -V >> $psfile\n"}
  #print CSH "awk '{print \$1,\$2}' $evefile |psxy -N $J $R -K -O -P -V $src >> $psfile\n";
  print CSH "awk '\$1 == \"R\" {print \$2,\$3}' $recfile |psxy -N $J $R -K -O -P -V $rec_fw >> $psfile\n";
  print CSH "awk '\$1 == \"S\" {print \$2,\$3}' $recfile |psxy -N $J $R -K -O -P -V $src_fr >> $psfile\n";
  print CSH "psscale -C$cpt_vel $Dscale $Bscale1 -K -O -P -V >> $psfile \n";
  print CSH "pstext -N $J_title $R_title -K -O -P -V >>$psfile<<EOF\n $x_title $z_title $fsize2 0 $fontno CM $title \nEOF\n";

#-----------------------------
  print CSH "pstext -N $J_title $R_title -O -P -V >>$psfile<<EOF\n 10 10 $fsize_title 0 $fontno CM $title \nEOF\n"; # FINISH
  #print CSH "convert $psfile -rotate 90 $jpgfile\n";
  print CSH "convert $psfile $jpgfile\n";
  if($ieps==1){print CSH "convert $psfile $epsfile\n"}

  close (CSH); system("csh -f $cshfile");
  system("xv $jpgfile &");
}

#===========================================================================
if ($ifig05 == 1) {

  $Jwid = 1.5;
  $J = "-JM${Jwid}i";      # in lat-lon
  $origin = "-X1.5 -Y1.0";
  $J_title = "-JX${Jwid}";

  # colorbar
  $Dlen = 0.75*${Jwid}; $Dx = $Jwid/2; $Dy = -0.15*$Jwid;
  $Dscale = "-D$Dx/$Dy/$Dlen/0.10h";

  # shifting the subplots
  $xfac = 1.15;
  $yfac = 1.15;
  $dX = 0; $dY = $yfac*$Jwid; $shift1 = "-X$dX -Y$dY";
  $dX = -2*$dX; $dY = -$yfac*$Jwid; $shift2 = "-X$dX -Y$dY";

  #===============================================
  print "\nWriting CSH file...\n";

  open(CSH,">$cshfile");
  print CSH "gmtset BASEMAP_TYPE plain TICK_LENGTH 0.1c LABEL_FONT_SIZE $fsize3 ANOT_FONT_SIZE $fsize3 PLOT_DEGREE_FORMAT D HEADER_FONT $fontno ANOT_FONT $fontno LABEL_FONT $fontno HEADER_FONT_SIZE $fsize2 FRAME_PEN $fpen TICK_PEN $tpen\n";
  #===============================================

  if($irun0 == 0) { $iev = 5 } else { $iev = 1 }
  $stev = sprintf("%3.3i",$iev);

  $dir1 = "${dir}/${edir}${stev}";

  # get src-rec file
  $recfile = "${dir1}/sr.txt";  # src-rec for event 001
  if (not -f $recfile) { die("Check if $recfile exist or not\n") }
  open(IN,$recfile); @temp = <IN>; $nrec = @temp - 1;

  # pick frames
  $dt = 0.06;
  @frames  = (400,1200,2000,2800);
  $numf = @frames;

  #@files = glob("${dir1}/forward*");
  #$nfile = @files;

  # file names for figures
  $name = "talk_forward_${strun0}";
  $psfile = "$name.ps"; $jpgfile = "$name.jpg"; $epsfile = "$name.eps";
  $title = " Regular Wavefield";

  #------------------------------------
  # color for forward wavefield
  $pwr = 3;
  $cmax = 1.0;
  $norm = "1e-$pwr";

  $ss  = $cmax;
  $ds  = 2*$ss/$scale_color;
  $bs  = sprintf("%3.3e",0.9*$ss); # colorbar
  $T   = sprintf("-T%3.3e/%3.3e/%3.3e",-$ss,$ss,$ds);

  $cpt_wave = "color4.cpt";
  print CSH "makecpt -C$colorbar $T > temp1\n";
  print CSH "sed 's/^B.*/B       170     0      0  /' temp1 >  temp2\n";
  print CSH "sed 's/^F.*/F         0     0    255  /' temp2 > $cpt_wave\n";

  $BscaleS1 = sprintf("-B%2.2e:\" s ( x, y, t )  ( 10\@+-%2.2i\@+  m )\": -E10p",$bs,$pwr);  # $k = 0
  #------------------------------------

for ($i = 0; $i < $numf; $i++) {

   $j1 = $frames[$i];           # forward frame
   $time = sprintf("%04d",$j1*$dt);
   $snapshot_f = sprintf("${dir1}/forward_%05d",$j1);
   if (not -f $snapshot_f) {die("Check if $snapshot_f exist or not\n");}

   print CSH "echo $psfile\n";
   print CSH "echo $snapshot_f\n";

   $B = "-B1/1:\"t = $time s\"::.\"  \":wsne";
   $B_row1 = "-B1/1:\"t = $time s\"::.\"  \":wsne";
   if ($i == 0) { $B = $B_row1}
   #$B = "-B2/2:.\" \":wesn";

   if ($i == 0) {print CSH "psbasemap $B $R $J -K -V $origin > $psfile\n"}
   else         {print CSH "psbasemap $B $R $J -K -O -V $shift1 >> $psfile\n"}

   # PLOT THE FORWARD WAVEFIELD
   if($icolor==1) {print CSH "awk '{print \$1,\$2,\$3 / $norm}' $snapshot_f | pscontour $J $R $B -A- -C${cpt_wave} -I -K -O -V >> $psfile\n"}
   print CSH "pscoast $J $R $B $coast_info -K -O -V >> $psfile\n";
   if ($i == 0) {print CSH "psscale -C${cpt_wave} $Dscale $BscaleS1 -K -O -V >> $psfile \n";}
   print CSH "awk '\$1 == \"R\" {print \$2,\$3}' $recfile |psxy -N $J $R -K -O -V $rec >> $psfile\n";
   print CSH "awk '\$1 == \"S\" {print \$2,\$3}' $recfile |psxy -N $J $R -K -O -V $src_fw >> $psfile\n";

   # plot the time of the snapshot
   $tstr = "t = $time s";
   print CSH "pstext -N $J_title $R_title -K -O -V >>$psfile<<EOF\n -0.2 0.5 $fsize2 90 $fontno CM $tstr \nEOF\n";

}
  print CSH "pstext -N $J_title $R_title -K -O -V >>$psfile<<EOF\n 0.5 1.2 $fsize2 0 $fontno CM $title \nEOF\n";

#-----------------------------

  if($irun == 1 || $irun == 2) {      # data file generated in gji_figs.m (11-28-05)

    # get axes limits
    $ax_file     = "${dir}/time_series_axes.dat";
    if (not -f $ax_file) { die("Check if $ax_file exist or not\n") }
    open(IN,$ax_file); @axlims = <IN>;
    @pwr_vec = split(" ",$axlims[0]);
    @cmx_vec = split(" ",$axlims[1]);
    ($tmin,$tmax,$junk,$junk) = split(" ",$axlims[2]);
    print "\n $tmin $tmax \n @pwr_vec \n @cmx_vec \n\n";
    @norm = ("1e$pwr_vec[0]","1e$pwr_vec[1]","1e$pwr_vec[2]","1e$pwr_vec[3]");

    # make plot bounds
    $R1 = "-R${tmin}/${tmax}/-${cmx_vec[0]}/${cmx_vec[0]}";
    $R2 = "-R${tmin}/${tmax}/-${cmx_vec[1]}/${cmx_vec[1]}";
    $R3 = "-R${tmin}/${tmax}/-${cmx_vec[2]}/${cmx_vec[2]}";
    $R4 = "-R${tmin}/${tmax}/-${cmx_vec[3]}/${cmx_vec[3]}";
    #print "\n $R1 \n $R2 \n $R3 \n $R4 \n";

    # get time series
    $series_file = "${dir}/time_series.dat";
    if (not -f $series_file) { die("Check if $series_file exist or not\n") }

    # get measurement
    $meas_file = "${dir}/measure_vec.dat";
    if (not -f $meas_file) { die("Check if $meas_file exist or not\n") }
    open(IN,$meas_file); $meas = <IN>; chomp($meas);
    $stmeas = sprintf(" \@~\104\@~T = %.2f s",$meas);

    print "\n $stmeas \n";

    # plot dimensions
    $Jseis = "-JX4/${Jwid}";
    $dY = -2.7*$yfac*$Jwid;
    $shift = "-X3 -Y${dY}";

    if($irun == 2) {$units = "m\@+-1\@+"} else {$units = "m\@+-1\@+ s"}

    # titles (note the units)
    $title1 = sprintf("Regular Source Function  (10\@+%2.2i\@+ kg s\@+-2\@+)",$pwr_vec[0]);
    $title2 = sprintf("Displacement  (10\@+%2.2i\@+ m)",$pwr_vec[1]);
    $title3 = sprintf("Velocity  (10\@+%2.2i\@+ m s\@+-1\@+)",$pwr_vec[2]);
    $title4 = sprintf("Adjoint Source Function  (10\@+%2.2i\@+ $units)",$pwr_vec[3]);

    $B1 = sprintf("-B%3.3f:\"Time  (s)\":/%3.3f:\" \":WeSn",20,2);
    $B2 = sprintf("-B%3.3f:\"Time  (s)\":/%3.3f:\" \":Wesn",20,2);
    $B3 = sprintf("-B%3.3f:\"Time  (s)\":/%3.3f:\" \":Wesn",20,1);
    $B4 = sprintf("-B%3.3f:\"Time  (s)\":/%3.3f:\" \":Wesn",20,1);

    $dY = 1.2 * $yfac * $Jwid;
    $shifty = "-Y${dY}";
    $x_tit = 1.5;
    #$y_tit = 0.9;
    $y_tit = 1.15;
    print CSH "psbasemap $B1 $R1 $Jseis -K -O -V $shift >> $psfile\n";
    print CSH "awk '{print \$1,\$2/$norm[0]}' $series_file | psxy $Jseis $R1 -K -O -V -W1p >> $psfile\n";
    print CSH "pstext -N $J_title $R_title -K -O -V >>$psfile<<EOF\n $x_tit $y_tit $fsize2 0 $fontno CM $title1 \nEOF\n";

    #print CSH "psbasemap $B2 $R2 $Jseis -K -O -V $shifty >> $psfile\n";
    #print CSH "awk '{print \$1,\$3/$norm[1]}' $series_file | psxy $Jseis $R2 -K -O -V -W1p >> $psfile\n";
    #print CSH "awk '{print \$1,\$4/$norm[1]}' $series_file | psxy $Jseis $R2 -K -O -V -W1/255/0/0tap >> $psfile\n";
    #print CSH "pstext -N $J_title $R_title -K -O -V >>$psfile<<EOF\n $x_tit $y_tit $fsize2 0 $fontno CM $title2 \nEOF\n";
    #print CSH "pstext -N $Jseis $R2 -K -O -V >>$psfile<<EOF\n 150 -3 $fsize2 0 $fontno CM $stmeas \nEOF\n";

    print CSH "psbasemap $B3 $R3 $Jseis -K -O -V $shifty >> $psfile\n";
    print CSH "pstext -N $J_title $R_title -K -O -V >>$psfile<<EOF\n $x_tit $y_tit $fsize2 0 $fontno CM $title3 \nEOF\n";
    if($irun == 1) {
      print CSH "awk '{print \$1,\$5/$norm[2]}' $series_file | psxy $Jseis $R3 -K -O -V -W1p >> $psfile\n";
      print CSH "awk '{print \$1,\$6/$norm[2]}' $series_file | psxy $Jseis $R3 -K -O -V -W1/255/0/0tap >> $psfile\n";
      print CSH "pstext -N $Jseis $R3 -K -O -V >>$psfile<<EOF\n 150 1.5 $fsize2 0 $fontno CM $stmeas \nEOF\n";
      print CSH "pstext -N $Jseis $R3 -K -O -V           >>$psfile<<EOF\n 143 -1.5 $fsize3 0 $fontno RM DATA \nEOF\n";
      print CSH "pstext -N $Jseis $R3 -G255/0/0 -K -O -V >>$psfile<<EOF\n 161 -1.5 $fsize3 0 $fontno LM SYN \nEOF\n";
    } else {
       print CSH "awk '{print \$1,\$6/$norm[2]}' $series_file | psxy $Jseis $R3 -K -O -V -W1p >> $psfile\n";
    }

    print CSH "psbasemap $B4 $R4 $Jseis -K -O -V $shifty >> $psfile\n";
    print CSH "awk '{print \$1,\$7/$norm[3]}' $series_file | psxy $Jseis $R4 -K -O -V -W1p >> $psfile\n";
    print CSH "pstext -N $J_title $R_title -K -O -V >>$psfile<<EOF\n $x_tit $y_tit $fsize2 0 $fontno CM $title4 \nEOF\n";

  }  # $irun == 1

#-----------------------------
  print CSH "pstext -N $J_title $R_title -O -V >>$psfile<<EOF\n 10 10 $fsize_title 0 $fontno CM $title \nEOF\n"; # FINISH
  print CSH "convert $psfile -rotate 90 $jpgfile\n";
  if($ieps==1){print CSH "convert $psfile $epsfile\n"}

  close (CSH); system("csh -f $cshfile");
  system("xv $jpgfile &");

}

#===========================================================================
if ($ifig06 == 1) {

  $Jwid = 5;
  $J = "-JX${Jwid}/-${Jwid}";  # flip the y-axis
  $origin = "-X3.5 -Y2.0";
  $origin = "-X1.5 -Y2.0";
  $J_title = "-JX${Jwid}";

  # colorbar
  $Dlen = 0.75*${Jwid}; $Dx = $Jwid/2; $Dy = -0.1*$Jwid;
  $Dscale = "-D$Dx/$Dy/$Dlen/0.15h";

  # shifting the subplots
  $xfac = 1.15;
  $yfac = 1.15;
  $dX = 0; $dY = $yfac*$Jwid; $shift1 = "-X$dX -Y$dY";
  $dX = -2*$dX; $dY = -$yfac*$Jwid; $shift2 = "-X$dX -Y$dY";

  #===============================================
  print "\nWriting CSH file...\n";

  open(CSH,">$cshfile");
  print CSH "gmtset BASEMAP_TYPE plain TICK_LENGTH 0.2c LABEL_FONT_SIZE $fsize1 ANOT_FONT_SIZE $fsize1 PLOT_DEGREE_FORMAT D HEADER_FONT $fontno ANOT_FONT $fontno LABEL_FONT $fontno HEADER_FONT_SIZE $fsize0 FRAME_PEN $fpen TICK_PEN $tpen\n";
  #===============================================

  $slab = "ker";
  $slab = "ray";

  $dir = "/home/carltape/sem2d/2d_adjoint_banana/mat_SAVE/design_matrix/run_0000";
  $file_hess = "${dir}/AtA_${slab}_gmt.dat";
  if (not -f $file_hess) {die("Check if $file_hess exist or not\n")}

  $ngrid = 286;

  # file names for figures
  $name = "hessian_${slab}_000";
  $psfile = "$name.ps"; $jpgfile = "$name.jpg"; $epsfile = "$name.eps";
  $title = "The Hessian  (G\@+T\@+G)";

  # find the mins and maxs of Gik
  ($imin,$imax,$jmin,$jmax,$Gmin,$Gmax) = split(" ",`minmax -C $file_hess`);

   # find the power of the max value of AtA
  ($Gmax_abs) = sort{$b <=> $a} (abs($Gmin),abs($Gmax));
  $pwr = int( log($Gmax_abs) / log(10) );
  if($pwr < 0) {$pwr = -1 * ( abs($pwr)+1 )}

  #------------------------------------
  # color for hessian
  $norm = "1e$pwr";
  $fac = 0.5;
  $cmax = $fac*$Gmax_abs/$norm;
  $cmax = 1;

  print "\n $imin, $imax, $jmin, $jmax \n $Gmin, $Gmax, $Gmax_abs, $pwr \n $norm $fac $cmax \n";
  #die("testing");

  $ss  = $cmax;
  $ds  = 2*$ss/$scale_color;
  $T   = sprintf("-T%3.3e/%3.3e/%3.3e",-$ss,$ss,$ds);

  $cpt_hess = "color.cpt";
  print CSH "makecpt -C$colorbar $T > temp1\n";
  print CSH "sed 's/^B.*/B       170     0      0  /' temp1 >  temp2\n";
  print CSH "sed 's/^F.*/F         0     0    255  /' temp2 > $cpt_hess\n";

  $bs  = 0.5*$ss;
  $BscaleH = sprintf("-B%.2f:\" G\@-ik\@- for $slab  ( 10\@+%2.2i\@+ )\": -Ef10p",$bs,$pwr);  # $k = 0
  #------------------------------------

  print "\n $BscaleH \n $T  \n";

  $R = "-R1/${ngrid}/1/${ngrid}";
  $B = "-B40:\" column index\":/40:\" row index\"::.\"  \":WsNe";

  print CSH "psbasemap $B $R $J -K -P -V $origin > $psfile\n";

  #if($icolor==1) {print CSH "awk '{print \$1,\$2,\$3 / $norm}' $file_hess | pscontour $J $R $B -A- -C${cpt_hess} -I -K -O -P -V >> $psfile\n"}

  # make grd file, then plot
  $grdfile = "temp.grd";
  $grd_info = "-I1/1";
  print CSH "awk '{print \$1,\$2,\$3 / $norm}' $file_hess | xyz2grd -G$grdfile $grd_info $R \n";
  print CSH "grdimage $grdfile $R $J $B -C${cpt_hess} -T -K -O -P -V >> $psfile\n";

  print CSH "psscale -C${cpt_hess} $Dscale $BscaleH -K -O -P -V >> $psfile \n";
  print CSH "psbasemap $B $R $J -K -O -P -V >> $psfile\n";
  #print CSH "pstext -N $J_title $R_title -K -O -P -V >>$psfile<<EOF\n 0.5 1.2 $fsize_title 0 $fontno CM $title \nEOF\n";

#-----------------------------
  print CSH "pstext -N $J_title $R_title -O -P -V >>$psfile<<EOF\n 10 10 $fsize_title 0 $fontno CM $title \nEOF\n"; # FINISH
  #print CSH "convert $psfile -rotate 90 $jpgfile\n";
  print CSH "convert $psfile $jpgfile\n";
  if($ieps==1){print CSH "convert $psfile $epsfile\n"}

  close (CSH); system("csh -f $cshfile");
  system("xv $jpgfile &");

}

#===========================================================================
# Hessian plot copied from ifig06 on 07-March-2006

if ($ifig07 == 1) {

  $Jwid = 2.75;
  $J = "-JM${Jwid}i";           # in lat-lon
  $Jb = "-JX${Jwid}/-${Jwid}";  # flip the y-axis
  $origin = "-X1.0 -Y3.0";
  $J_title = "-JX${Jwid}";
  $x_title = 0.05;
  $z_title = 1.15;

  # colorbar
  $Dlen = 0.6*${Jwid}; $Dx = $Jwid/2; $Dy = -0.08*$Jwid;
  $Dscale = "-D$Dx/$Dy/$Dlen/0.15h";

  # shifting the subplots
  $xfac = 1.2;
  $yfac = 1.15;
  $dX = $xfac*$Jwid; $dY = 0; $shift1 = "-X$dX -Y$dY";

  #===============================================
  print "\nWriting CSH file...\n";

  open(CSH,">$cshfile");
  print CSH "gmtset BASEMAP_TYPE plain TICK_LENGTH 0.15c LABEL_FONT_SIZE $fsize2 ANOT_FONT_SIZE 10 PLOT_DEGREE_FORMAT D HEADER_FONT $fontno ANOT_FONT $fontno LABEL_FONT $fontno HEADER_FONT_SIZE $fsize1 FRAME_PEN $fpen TICK_PEN $tpen\n";
  #===============================================

  $slab = "ker";

  # locate Hessian data file
  $dirb = "/home/carltape/sem2d/2d_adjoint_banana/mat_SAVE/design_matrix/run_0000";
  $file_hess = "${dirb}/AtA_${slab}_gmt.dat";
  if (not -f $file_hess) {die("Check if $file_hess exist or not\n")}

  # locate Hessian diagonal data file
  $dirb = "/home/carltape/sem2d/2d_adjoint_banana/mat_SAVE/design_vector/run_0020";
  $file_diag = "${dirb}/m_IMODEL_3_rayleigh_id01_plot";
  if (not -f $file_diag) {die("Check if $file_diag exist or not\n")}

  # locate gradient data file
  $dirb = "/home/carltape/sem2d/2d_adjoint_banana/mat_SAVE/gradient_vector/run_0020";
  $file_grad = "${dirb}/m_IMODEL_3_rayleigh_plot";
  if (not -f $file_grad) {die("Check if $file_grad exist or not\n")}

  $ngrid = 286;

  # file names for figures
  $name = "hessian_and_gradient_ker";
  $psfile = "$name.ps"; $jpgfile = "$name.jpg"; $epsfile = "$name.eps";
  $title = "The Hessian  (G\@+T\@+G)";

  # find the mins and maxs of Hij
  ($imin,$imax,$jmin,$jmax,$Gmin,$Gmax) = split(" ",`minmax -C $file_hess`);

   # find the power of the max value of AtA
  ($Gmax_abs) = sort{$b <=> $a} (abs($Gmin),abs($Gmax));
  $pwr = int( log($Gmax_abs) / log(10) );
  if($pwr < 0) {$pwr = -1 * ( abs($pwr)+1 )}

  #------------------------------------
  # color for hessian
  $normH = "1e$pwr";
  $fac = 0.5;
  $cmax = $fac*$Gmax_abs/$normH;
  $cmax = 1.01;

  print "\n $imin, $imax, $jmin, $jmax \n $Gmin, $Gmax, $Gmax_abs, $pwr \n $normH $fac $cmax \n";
  #die("testing");

  $ss  = $cmax;
  $ds  = 2*$ss/$scale_color;
  $T   = sprintf("-T%3.3e/%3.3e/%3.3e",-$ss,$ss,$ds);

  $cpt_hess = "color.cpt";
  print CSH "makecpt -C$colorbar $T > temp1\n";
  print CSH "sed 's/^B.*/B       170     0      0  /' temp1 >  temp2\n";
  print CSH "sed 's/^F.*/F         0     0    255  /' temp2 > $cpt_hess\n";

  $bs  = 0.5;
  #$BscaleH = sprintf("-B%.2f:\" Hessian matrix H\@-ij\@-  (10\@+%2.2i\@+ s\@+2\@+)\": -Ef10p",$bs,$pwr);  # $k = 0
  $BscaleH = sprintf("-B%.2f:\" Hessian matrix H  (10\@+%2.2i\@+ s\@+2\@+)\": -Ef10p",$bs,$pwr);
  print "\n $BscaleH \n $T  \n";

  #------------------------------------
  # color for gradient

  $pwr = 4;
  $normG = "1e$pwr";
  $cmax = 1.01;

  print "\n $imin, $imax, $jmin, $jmax \n $Gmin, $Gmax, $Gmax_abs, $pwr \n $normG $fac $cmax \n";
  #die("testing");

  $ss  = $cmax;
  $ds  = 2*$ss/$scale_color;
  $T   = sprintf("-T%3.3e/%3.3e/%3.3e",-$ss,$ss,$ds);

  $cpt_grad = "color_grad.cpt";
  print CSH "makecpt -C$colorbar $T > temp1\n";
  print CSH "sed 's/^B.*/B       170     0      0  /' temp1 >  temp2\n";
  print CSH "sed 's/^F.*/F         0     0    255  /' temp2 > $cpt_grad\n";

  $bs  = 0.5;
  #$BscaleG = sprintf("-B%.2f:\" gradient vector g\@-k\@-  (10\@+%2.2i\@+ s\@+2\@+)\": -Ef10p",$bs,$pwr);  # $k = 0
  $BscaleG = sprintf("-B%.2f:\" gradient vector g  (10\@+%2.2i\@+ s\@+2\@+)\": -Ef10p",$bs,$pwr);  # $k = 0
  print "\n $BscaleG \n $T  \n";

  #------------------------------------
  # color for Hessian diagonal

  $pwr = 4;
  $normD = "1e$pwr";
  $cmax = 4.01;

  $ss  = $cmax;
  $ds  = $ss/$scale_color;
  $T   = sprintf("-T%3.3e/%3.3e/%3.3e",0,$ss,$ds);

  $cpt_diag = "color_diag.cpt";
  #print CSH "makecpt -Chot $T > $cpt_diag\n";
  print CSH "makecpt -Crainbow $T > temp1\n";
  #print CSH "sed 's/^B.*/B       255   255    255  /' temp1 >  temp2\n";
  #print CSH "sed 's/^F.*/F         0     0      0  /' temp2 > $cpt_diag\n";
  print CSH "sed 's/^B.*/B       255     0    255  /' temp1 >  temp2\n";
  print CSH "sed 's/^F.*/F       255     0      0  /' temp2 > $cpt_diag\n";

  $bs  = 1.0;
  $BscaleD = sprintf("-B%.2f:\" Hessian diagonal H\@-ii\@-   (10\@+%2.2i\@+ s\@+2\@+)\": -Ef10p",$bs,$pwr);  # $k = 0
  print "\n $BscaleD \n $T  \n";

  #------------------------------------
  # plot Hessian

  $Rhess = "-R1/${ngrid}/1/${ngrid}";
  $Bhess = "-B40:\" column index\":/40:\" row index\"::.\"  \":WsNe";
  $title = "(a)";

  print CSH "psbasemap $Bhess $Rhess $Jb -K -V $origin > $psfile\n";

  #if($icolor==1) {print CSH "awk '{print \$1,\$2,\$3 / $normH}' $file_hess | pscontour $Jb $Rhess $Bhess -A- -C${cpt_hess} -I -K -O -V >> $psfile\n"}

  # make grd file, then plot
  $grdfile = "temp.grd";
  $grd_info = "-I1/1";
  print CSH "awk '{print \$1,\$2,\$3 / $normH}' $file_hess | xyz2grd -G$grdfile $grd_info $Rhess \n";
  print CSH "grdimage $grdfile $Rhess $Jb $Bhess -C${cpt_hess} -T -K -O -V >> $psfile\n";
  print CSH "psscale -C${cpt_hess} $Dscale $BscaleH -K -O -V >> $psfile \n";
  print CSH "psbasemap $Bhess $Rhess $Jb -K -O -V >> $psfile\n";
  print CSH "pstext -N $J_title $R_title -K -O -V >>$psfile<<EOF\n $x_title $z_title $fsize_title 0 $fontno CM $title \nEOF\n";

  #------------------------------------
  # plot Hessian diagonal

  $recfile = "${dir}/${edir}001/sr.txt";  # src-rec for first event
  $evefile = "${dir0}/events_lonlat.dat";
  if (not -f $recfile) { die("Check if $recfile exist or not\n") }
  if (not -f $evefile) { die("Check if $evefile exist or not\n") }

  $B = "$B0".$Bopts[5];
  $title = "(b)";

  print CSH "psbasemap $B $R $J -K -O -V $shift1 >> $psfile\n";
  if($icolor==1) {print CSH "awk '{print \$1,\$2,\$3 / $normD}' ${file_diag} | pscontour $R $J -A- -C${cpt_diag} -I -K -O -V >> $psfile\n"}
  print CSH "pscoast $J $R $coast_info -K -O -V >> $psfile\n";
  print CSH "awk '\$1 == \"R\" {print \$2,\$3}' $recfile |psxy -N $J $R -K -O -V $rec >> $psfile\n";
  print CSH "awk '{print \$1,\$2}' $evefile |psxy -N $J $R -K -O -V $src -W0.5p >> $psfile\n";
  print CSH "psscale -C${cpt_diag} $Dscale $BscaleD -K -O -V >> $psfile \n";
  print CSH "pstext -N $J_title $R_title -K -O -V >>$psfile<<EOF\n $x_title $z_title $fsize_title 0 $fontno CM $title \nEOF\n";

  #------------------------------------
  # plot gradient

  $B = "$B0".$Bopts[5];
  $title = "(c)";

  print CSH "psbasemap $B $R $J -K -O -V $shift1 >> $psfile\n";
  if($icolor==1) {print CSH "awk '{print \$1,\$2,\$3 / $normG}' ${file_grad} | pscontour $R $J -A- -C${cpt_grad} -I -K -O -V >> $psfile\n"}
  print CSH "pscoast $J $R $coast_info -K -O -V >> $psfile\n";
  print CSH "awk '\$1 == \"R\" {print \$2,\$3}' $recfile |psxy -N $J $R -K -O -V $rec >> $psfile\n";
  print CSH "awk '{print \$1,\$2}' $evefile |psxy -N $J $R -K -O -V $src -W0.5p >> $psfile\n";
  print CSH "psscale -C${cpt_grad} $Dscale $BscaleG -K -O -V >> $psfile \n";
  print CSH "pstext -N $J_title $R_title -K -O -V >>$psfile<<EOF\n $x_title $z_title $fsize_title 0 $fontno CM $title \nEOF\n";

#-----------------------------
  print CSH "pstext -N $J_title $R_title -O -V >>$psfile<<EOF\n 10 10 $fsize_title 0 $fontno CM $title \nEOF\n"; # FINISH
  print CSH "convert $psfile -rotate 90 $jpgfile\n";
  #print CSH "convert $psfile $jpgfile\n";
  if($ieps==1){print CSH "convert $psfile $epsfile\n"}

  close (CSH); system("csh -f $cshfile");
  system("xv $jpgfile &");

}

#===========================================================================
if ($ifig05p == 1) {

  $Jwid = 1.5;
  $J = "-JM${Jwid}i";      # in lat-lon
  $origin = "-X1.5 -Y6.0";
  $J_title = "-JX${Jwid}";

  # colorbar
  $Dlen = 0.75*${Jwid}; $Dx = $Jwid/2; $Dy = -0.15*$Jwid;
  $Dscale = "-D$Dx/$Dy/$Dlen/0.10h";

  # shifting the subplots
  $xfac = 1.25;
  $yfac = 1.45;
  $dX = 0; $dY = -$yfac*$Jwid; $shift1 = "-X$dX -Y$dY";

  #===============================================
  print "\nWriting CSH file...\n";

  open(CSH,">$cshfile");
  print CSH "gmtset BASEMAP_TYPE plain TICK_LENGTH 0.1c LABEL_FONT_SIZE $fsize2 ANOT_FONT_SIZE $fsize3 PLOT_DEGREE_FORMAT D HEADER_FONT $fontno ANOT_FONT $fontno LABEL_FONT $fontno HEADER_FONT_SIZE $fsize1 FRAME_PEN $fpen TICK_PEN $tpen\n";
  #===============================================

  $dir = "/home/carltape/2d_adjoint/OUTPUT/run_0001";
  $dir1 = "${dir}/event_001";

  # get src-rec file
  $recfile = "${dir1}/sr.txt";  # src-rec for event 001
  if (not -f $recfile) { die("Check if $recfile exist or not\n") }
  open(IN,$recfile); @temp = <IN>; $nrec = @temp - 1;

  # file names for figures
  $name = "talk_seis_01";
  $psfile = "$name.ps"; $jpgfile = "$name.jpg"; $epsfile = "$name.eps";
  $title = " Regular Wavefield";

    # get axes limits
    $ax_file = "${dir}/time_series_axes.dat";
    if (not -f $ax_file) { die("Check if $ax_file exist or not\n") }
    open(IN,$ax_file); @axlims = <IN>;
    @pwr_vec = split(" ",$axlims[0]);
    @cmx_vec = split(" ",$axlims[1]);
    ($tmin,$tmax,$junk,$junk) = split(" ",$axlims[2]);
    print "\n $tmin $tmax \n @pwr_vec \n @cmx_vec \n\n";
    @norm = ("1e$pwr_vec[0]","1e$pwr_vec[1]","1e$pwr_vec[2]","1e$pwr_vec[3]");

    # make plot bounds
    $R1 = "-R${tmin}/${tmax}/-${cmx_vec[0]}/${cmx_vec[0]}";
    $R2 = "-R${tmin}/${tmax}/-${cmx_vec[1]}/${cmx_vec[1]}";
    $R3 = "-R${tmin}/${tmax}/-${cmx_vec[2]}/${cmx_vec[2]}";
    $R4 = "-R${tmin}/${tmax}/-${cmx_vec[3]}/${cmx_vec[3]}";
    #print "\n $R1 \n $R2 \n $R3 \n $R4 \n";

    # get time series
    $series_file = "${dir}/time_series.dat";
    if (not -f $series_file) { die("Check if $series_file exist or not\n") }

    # get measurement
    $meas_file = "${dir}/measure_vec.dat";
    if (not -f $meas_file) { die("Check if $meas_file exist or not\n") }
    open(IN,$meas_file); $meas = <IN>; chomp($meas);
    $stmeas = sprintf(" \@~\104\@~T = %.2f s",$meas);
    print "\n $stmeas \n";

    # plot dimensions
    $Jseis = "-JX4/${Jwid}";
    $units = "m\@+3\@+ s";

    $title1 = "(a)  Regular Source";
    #$title2 = "(b)  Displacement at Receiver";
    $title3 = "(b)  Velocity at Receiver";
    $title4 = "(c)  Adjoint Source";

    $title1y = sprintf("Amplitude  (10\@+%2.2i\@+ kg s\@+-2\@+)",$pwr_vec[0]);
    #$title2y = sprintf("Displacement  (10\@+%2.2i\@+ m)",$pwr_vec[1]);
    $title3y = sprintf("Speed  (10\@+%2.2i\@+ m s\@+-1\@+)",$pwr_vec[2]);
    $title4y = sprintf("Amplitude  (10\@+%2.2i\@+ m\@+-1\@+ s)",$pwr_vec[3]);

    $B1 = sprintf("-B%3.3f:\" \":/%3.3f:\"$title1y\":WeSn",20,2);
    #$B2 = sprintf("-B%3.3f:\" \":/%3.3f:\"$title2y\":WeSn",20,2);
    $B3 = sprintf("-B%3.3f:\" \":/%3.3f:\"$title3y\":WeSn",20,1);
    $B4 = sprintf("-B%3.3f:\"Time  (s)\":/%3.3f:\"$title4y\":WeSn",20,1);

    $x_tit = 0.1;
    $y_tit = 1.15;
    $shift = $shift1;

    print CSH "psbasemap $B1 $R1 $Jseis -K -P -V $origin > $psfile\n";  # START
    print CSH "awk '{print \$1,\$2/$norm[0]}' $series_file | psxy $Jseis $R1 -K -O -P -V -W1p >> $psfile\n";
    print CSH "pstext -N $J_title $R_title -K -O -P -V >>$psfile<<EOF\n $x_tit $y_tit $fsize2 0 $fontno LM $title1 \nEOF\n";

    print CSH "psbasemap $B3 $R3 $Jseis -K -O -P -V $shift >> $psfile\n";
    print CSH "pstext -N $J_title $R_title -K -O -P -V >>$psfile<<EOF\n $x_tit $y_tit $fsize2 0 $fontno LM $title3 \nEOF\n";
    print CSH "awk '{print \$1,\$5/$norm[2]}' $series_file | psxy $Jseis $R3 -K -O -P -V -W1p >> $psfile\n";
    print CSH "awk '{print \$1,\$6/$norm[2]}' $series_file | psxy $Jseis $R3 -K -O -P -V -W1/255/0/0tap >> $psfile\n";
    print CSH "pstext -N $Jseis $R3 -K -O -P -V >>$psfile<<EOF\n 150 1.5 $fsize2 0 $fontno CM $stmeas \nEOF\n";
    print CSH "pstext -N $Jseis $R3 -K -O -P -V           >>$psfile<<EOF\n 143 -1.5 $fsize3 0 $fontno RM DATA \nEOF\n";
    print CSH "pstext -N $Jseis $R3 -G255/0/0 -K -O -P -V >>$psfile<<EOF\n 161 -1.5 $fsize3 0 $fontno LM SYN \nEOF\n";

    print CSH "psbasemap $B4 $R4 $Jseis -K -O -P -V $shift >> $psfile\n";
    print CSH "awk '{print \$1,\$7/$norm[3]}' $series_file | psxy $Jseis $R4 -K -O -P -V -W1p >> $psfile\n";
    print CSH "pstext -N $J_title $R_title -K -O -P -V >>$psfile<<EOF\n $x_tit $y_tit $fsize2 0 $fontno LM $title4 \nEOF\n";

#-----------------------------
  print CSH "pstext -N $J_title $R_title -O -P -V >>$psfile<<EOF\n 10 10 $fsize_title 0 $fontno CM $title \nEOF\n"; # FINISH
  #print CSH "convert $psfile -rotate 90 $jpgfile\n";
  print CSH "convert $psfile $jpgfile\n";
  if($ieps==1){print CSH "convert $psfile $epsfile\n"}

  close (CSH); system("csh -f $cshfile");
  system("xv $jpgfile &");

}


#===========================================================================
#===========================================================================
if ($ifig10 == 1) {

  $Jwid = 2;
  $J = "-JM${Jwid}i";      # in lat-lon
  $origin = "-X1 -Y6.0";

  $J_poly_wid = $Jwid;
  $J_poly = "-JX${J_poly_wid}/${Jwid}";
  $J_title = "-JX${Jwid}";
  $J_chi = "-JX${Jwid}i/${Jwid}il";        # note log scale on y-axis

  $fsize_title = 10;

  # shifting the subplots
  $xfac = 1.2;
  $yfac = 1.35;
  $dX = $xfac*$Jwid; $dY = 0; $shift1 = "-X$dX -Y$dY";
  $dX = -2*$dX; $dY = -$yfac*$Jwid; $shift2 = "-X$dX -Y$dY";

  # colorbar
  $Dlen = 1.2; $Dx = $Jwid/2; $Dy = -0.05*$Jwid;
  $Dscale = "-D$Dx/$Dy/$Dlen/0.10h";

  #===============================================
  print "\nWriting CSH file...\n";

  open(CSH,">$cshfile");
  print CSH "gmtset BASEMAP_TYPE plain TICK_LENGTH $tick LABEL_FONT_SIZE $fsize3 ANOT_FONT_SIZE $fsize3 PLOT_DEGREE_FORMAT D HEADER_FONT $fontno ANOT_FONT $fontno LABEL_FONT $fontno HEADER_FONT_SIZE $fsize2 FRAME_PEN $fpen TICK_PEN $tpen\n";
  #===============================================

  # get the summed kernel and chi files
  #$file2 = "${dir}/summed_ker.dat";
  $file2 = "${dir}/fun_smooth.dat";       # smoothed kernel
  $file3 = "${dir}/summed_chi_r.dat";
  $recfile = "${dir}/${edir}001/sr.txt";  # src-rec for first event
  $evefile = "${dir0}/events_lonlat.dat";
  $chifile = "${dir}/summed_chi_all.dat";

  if (not -f $file2)   { die("Check if $file2 exist or not\n") }
  if (not -f $file3)   { die("Check if $file3 exist or not\n") }
  if (not -f $recfile) { die("Check if $recfile exist or not\n") }
  if (not -f $chifile) { die("Check if $chifile exist or not\n") }

  # file names for figures
  $name = "talk_cg10_${strun0}"; $psfile = "$name.ps"; $jpgfile = "$name.jpg"; $epsfile = "$name.eps";

  # number of receivers
  open(IN,$file3); @temp = <IN>; $nrec = @temp;

  #=============================================
  # model for the data

  $B = "$B0".$Bopts[15];
  $title = "Phase speed for data";

  # phase speed map
  print CSH "psbasemap $B $R $J -K -V $origin > $psfile\n"; # START
  #if($icolor==1) {print CSH "awk '{print \$1,\$2,\$3/1000}' $file1dat | pscontour $R $J -A- -C$cpt_vel -I -K -O -V >> $psfile\n"}
  if($icolor==1) {print CSH "awk '{print \$1,\$2,\$4*100}' $file1dat | pscontour $R $J -A- -C$cpt_vel -I -K -O -V  >> $psfile\n"}
  print CSH "pscoast $J $R $coast_info -K -O -V >> $psfile\n";
  if($ishelf) {print CSH "awk '{print \$2,\$1}' $shelf_file |psxy $J $R $Wshelf -K -O -V >> $psfile\n"}
  print CSH "awk '\$1 == \"R\" {print \$2,\$3}' $recfile |psxy -N $J $R -K -O -V $rec >> $psfile\n";
  print CSH "awk '{print \$1,\$2}' $evefile |psxy -N $J $R -K -O -V $src >> $psfile\n";
  print CSH "psscale -C$cpt_vel $Dscale $Bscale1b -K -O -V >> $psfile \n";
  print CSH "pstext -N $J_title $R_title -K -O -V >>$psfile<<EOF\n $x_title $z_title $fsize_title 0 $fontno CM $title \nEOF\n";

  #=============================================
  # model for synthetics

  $title = "Phase speed model $smod";
  $shift = $shift1;

  # phase speed map
  print CSH "psbasemap $B $R $J -K -O -V $shift >> $psfile\n";
  #if($icolor==1) {print CSH "awk '{print \$1,\$2,\$3/1000}' $file1syn | pscontour $R $J -A- -C$cpt_vel -I -O -K -V >> $psfile\n"}
  if($icolor==1) {print CSH "awk '{print \$1,\$2,\$4*100}' $file1syn | pscontour $R $J -A- -C$cpt_vel -I -O -K -V >> $psfile\n"}
  print CSH "pscoast $J $R $coast_info -K -O -V >> $psfile\n";
  if($ishelf) {print CSH "awk '{print \$2,\$1}' $shelf_file |psxy $J $R $Wshelf -K -O -V >> $psfile\n"}
  print CSH "awk '{print \$1,\$2}' $evefile |psxy -N $J $R -K -O -V $src >> $psfile\n";
  print CSH "awk '\$1 == \"R\" {print \$2,\$3}' $recfile |psxy -N $J $R -K -O -V $rec >> $psfile\n";
  #print CSH "psscale -C$cpt_vel $Dscale $Bscale1 -K -O -V >> $psfile \n";
  print CSH "pstext -N $J_title $R_title -K -O -V >>$psfile<<EOF\n $x_title $z_title $fsize_title 0 $fontno CM $title \nEOF\n";

  #=============================================
  # kernel -- SMOOTHED VERSION

  $shift = $shift1;
  $title = "Gradient (smoothed) for $smod";
  $B = "$B0".$Bopts[15];

  print CSH "psbasemap $B $R $J -K -O -V $shift >> $psfile\n";

  #if($icolor==1) {print CSH "awk '{print \$1,\$2,\$3 / $norm2 }' $file2 | pscontour $J $R -A- -C$cpt_ker -I -K -O -V >> $psfile\n"}
  if($icolor==1) {print CSH "awk '{print \$1,\$2,\$7 / $norm2 }' $file2 | pscontour $J $R -A- -C$cpt_ker -I -K -O -V >> $psfile\n"}

  print CSH "pscoast $J $R $coast_info -K -O -V >> $psfile\n";
  if($ishelf) {print CSH "awk '{print \$2,\$1}' $shelf_file |psxy $J $R $Wshelf -K -O -V >> $psfile\n"}
  print CSH "psscale -C$cpt_ker $Dscale $Bscale2c -K -O -V >> $psfile \n";
  print CSH "awk '{print \$1,\$2}' $evefile |psxy -N $J $R -K -O -V $src >> $psfile\n";
  print CSH "pstext -N $J_title $R_title -K -O -V >>$psfile<<EOF\n $x_title $z_title $fsize_title 0 $fontno CM $title \nEOF\n";
  print CSH "awk '\$1 == \"R\" {print \$2,\$3}' $recfile |psxy -N $J $R -K -O -V $rec >> $psfile\n";

  #=============================================
  # graph

  $shift = $shift2;
  $Bg = "$B1".$Bopts[8];
  $title = "Estimate test model m\@+0t\@+";

  print CSH "psbasemap $Bg $R_poly $J_poly -K -O -V $shift >> $psfile\n";
  print CSH "awk '{print \$1 / $normgx,\$2 / $normgy}' $poly_curve | psxy $c_info_rd $J_poly $R_poly -K -O -V >> $psfile\n";
  print CSH "awk '{print \$1 / $normgx,\$4 / $normgy}' $poly_curve | psxy $c_info_bd $J_poly $R_poly -K -O -V >> $psfile\n";
  print CSH "psxy $J_poly $R_poly $c_info_ks -K -O -V >>$psfile<<EOF\n$x1 $a1y1 \n$x1 $y1 \n$a1x1 $y1 \nEOF\n";
  print CSH "psxy $J_poly $R_poly $p_info_k -K -O -V >>$psfile<<EOF\n$x1 $y1 \nEOF\n";
  print CSH "psxy $J_poly $R_poly $p_info_w -K -O -V >>$psfile<<EOF\n$x2 $a1y1 \nEOF\n";
  print CSH "pstext -N $J_title $R_title -K -O -V >>$psfile<<EOF\n $x_title $z_title $fsize_title 0 $fontno CM $title \nEOF\n";

#=================================================================
# NEW MODEL

  # wave2d run number
  $strun = sprintf("%4.4i",$irun + 1);
  $dir = "${odir}$strun";

  $mod = $mods[1];
  $smod = "m\@+$mod\@+";

  # get the summed kernel and chi files
  $file1syn = "${dir}/$mfile_syn";
  #$file2 = "${dir}/summed_ker.dat";
  $file2 = "${dir}/fun_smooth.dat";
  $file3 = "${dir}/summed_chi_r.dat";
  $recfile = "${dir}/${edir}001/sr.txt";  # src-rec for first event
  $chifile = "${dir}/summed_chi_all.dat";

  if (not -f $file1syn){ die("Check if $file1syn exist or not\n") }
  if (not -f $file2)   { die("Check if $file2 exist or not\n") }
  if (not -f $file3)   { die("Check if $file3 exist or not\n") }
  if (not -f $recfile) { die("Check if $recfile exist or not\n") }
  if (not -f $chifile) { die("Check if $chifile exist or not\n") }
#=================================================================

  #=============================================
  # model for synthetics

  $title = "Phase speed model $smod";
  $shift = $shift1;

  # phase speed map
  print CSH "psbasemap $B $R $J -K -O -V $shift >> $psfile\n";
  if($icolor==1) {print CSH "awk '{print \$1,\$2,\$4*100}' $file1syn | pscontour $R $J -A- -C$cpt_vel -I -O -K -V >> $psfile\n"}
  print CSH "pscoast $J $R $coast_info -K -O -V >> $psfile\n";
  if($ishelf) {print CSH "awk '{print \$2,\$1}' $shelf_file |psxy $J $R $Wshelf -K -O -V >> $psfile\n"}
  print CSH "awk '{print \$1,\$2}' $evefile |psxy -N $J $R -K -O -V $src >> $psfile\n";
  print CSH "awk '\$1 == \"R\" {print \$2,\$3}' $recfile |psxy -N $J $R -K -O -V $rec >> $psfile\n";
  #print CSH "psscale -C$cpt_vel $Dscale $Bscale1 -K -O -V >> $psfile \n";
  print CSH "pstext -N $J_title $R_title -K -O -V >>$psfile<<EOF\n $x_title $z_title $fsize_title 0 $fontno CM $title \nEOF\n";

  #=============================================
  # kernel -- SMOOTHED VERSION

  $shift = $shift1;
  $title = "Gradient (smoothed) for $smod";
  $B = "$B0".$Bopts[15];

  print CSH "psbasemap $B $R $J -K -O -V $shift >> $psfile\n";
  #if($icolor==1) {print CSH "awk '{print \$1,\$2,\$3 / $norm2 }' $file2 | pscontour $J $R -A- -C$cpt_ker -I -K -O -V >> $psfile\n"}
  if($icolor==1) {print CSH "awk '{print \$1,\$2,\$7 / $norm2 }' $file2 | pscontour $J $R -A- -C$cpt_ker -I -K -O -V >> $psfile\n"}
  print CSH "pscoast $J $R $coast_info -K -O -V >> $psfile\n";
  if($ishelf) {print CSH "awk '{print \$2,\$1}' $shelf_file |psxy $J $R $Wshelf -K -O -V >> $psfile\n"}
  #print CSH "psscale -C$cpt_ker $Dscale $Bscale2 -K -O -V >> $psfile \n";
  print CSH "awk '{print \$1,\$2}' $evefile |psxy -N $J $R -K -O -V $src >> $psfile\n";
  print CSH "pstext -N $J_title $R_title -K -O -V >>$psfile<<EOF\n $x_title $z_title $fsize_title 0 $fontno CM $title \nEOF\n";
  print CSH "awk '\$1 == \"R\" {print \$2,\$3}' $recfile |psxy -N $J $R -K -O -V $rec >> $psfile\n";

  #=============================================
  # graph

  $shift = $shift2;
  $Bg = "$B1".$Bopts[8];
  $title = "Cubic interpolation to get m\@+1\@+";

  print CSH "psbasemap $Bg $R_poly $J_poly -K -O -V $shift >> $psfile\n";
  print CSH "awk '{print \$1 / $normgx,\$2 / $normgy}' $poly_curve | psxy $c_info_rd $J_poly $R_poly -K -O -V >> $psfile\n";
  print CSH "awk '{print \$1 / $normgx,\$3 / $normgy}' $poly_curve | psxy $c_info_rd $J_poly $R_poly -K -O -V >> $psfile\n";
  print CSH "awk '{print \$1 / $normgx,\$5 / $normgy}' $poly_curve | psxy $c_info_ks2 $J_poly $R_poly -K -O -V >> $psfile\n";
  print CSH "psxy $J_poly $R_poly $c_info_ks -K -O -V >>$psfile<<EOF\n$x1 $a1y1 \n$x1 $y1 \n$a1x1 $y1 \nEOF\n";
  print CSH "psxy $J_poly $R_poly $c_info_ks -K -O -V >>$psfile<<EOF\n$x3 $a1y1 \n$x3 $y3 \n$a1x1 $y3 \nEOF\n";
  print CSH "psxy $J_poly $R_poly $c_info_kd -K -O -V >>$psfile<<EOF\n$x5 $a1y1 \n$x5 $y5 \n$a1x1 $y5 \nEOF\n";
  print CSH "psxy $J_poly $R_poly $p_info_k -K -O -V >>$psfile<<EOF\n$x1 $y1 \n$x3 $y3\nEOF\n";
  print CSH "psxy $J_poly $R_poly $p_info_w -K -O -V >>$psfile<<EOF\n$x2 $a1y1 \n$x4 $y4 \n$x5 $y5 \nEOF\n";
  print CSH "pstext -N $J_title $R_title -K -O -V >>$psfile<<EOF\n $x_title $z_title $fsize_title 0 $fontno CM $title \nEOF\n";

#=================================================================
# NEW MODEL

  # wave2d run number
  $strun = sprintf("%4.4i",$irun + 2);
  $dir = "${odir}$strun";

  $mod = $mods[2];
  $smod = "m\@+$mod\@+";

  # get the summed kernel and chi files
  $file1syn = "${dir}/$mfile_syn";
  $file2 = "${dir}/summed_ker.dat";
  $file3 = "${dir}/summed_chi_r.dat";
  $recfile = "${dir}/${edir}001/sr.txt";  # src-rec for first event
  $chifile = "${dir}/summed_chi_all.dat";

  if (not -f $file1syn){ die("Check if $file1syn exist or not\n") }
  if (not -f $file2)   { die("Check if $file2 exist or not\n") }
  if (not -f $file3)   { die("Check if $file3 exist or not\n") }
  if (not -f $recfile) { die("Check if $recfile exist or not\n") }
  if (not -f $chifile) { die("Check if $chifile exist or not\n") }
#=================================================================

  #=============================================
  # model for synthetics

  $title = "Phase speed model $smod";
  $shift = $shift1;

  # phase speed map
  print CSH "psbasemap $B $R $J -K -O -V $shift >> $psfile\n";
  if($icolor==1) {print CSH "awk '{print \$1,\$2,\$4*100}' $file1syn | pscontour $R $J -A- -C$cpt_vel -I -O -K -V >> $psfile\n"}
  print CSH "pscoast $J $R $coast_info -K -O -V >> $psfile\n";
  if($ishelf) {print CSH "awk '{print \$2,\$1}' $shelf_file |psxy $J $R $Wshelf -K -O -V >> $psfile\n"}
  print CSH "awk '{print \$1,\$2}' $evefile |psxy -N $J $R -K -O -V $src >> $psfile\n";
  print CSH "awk '\$1 == \"R\" {print \$2,\$3}' $recfile |psxy -N $J $R -K -O -V $rec >> $psfile\n";
  #print CSH "psscale -C$cpt_vel $Dscale $Bscale1 -K -O -V >> $psfile \n";
  print CSH "pstext -N $J_title $R_title -K -O -V >>$psfile<<EOF\n $x_title $z_title $fsize_title 0 $fontno CM $title \nEOF\n";

  #=============================================
  # misfit vs iteration

  $shift = $shift1;
  $Bg = "$B3b".$Bopts[8];
  $title = "Misfit for first two models";

  print CSH "psbasemap $Bg $R_chi $J_chi -K -O -V $shift >> $psfile\n";
  print CSH "awk '{print \$1,\$2}' $chi_curve | psxy $c_info_rd $J_chi $R_chi -K -O -V >> $psfile\n";
  print CSH "psxy $J_chi $R_chi $p_info_k -K -O -V >>$psfile<<EOF\n0 $y1t \n1 $y7t \nEOF\n";
  print CSH "pstext -N $J_chi $R_chi -K -O -V >>$psfile<<EOF\n 0.75 $y1t $fsize2 0 $fontno LM $schi0 \nEOF\n";
  print CSH "pstext -N $J_chi $R_chi -K -O -V >>$psfile<<EOF\n 1.75 $y7t $fsize2 0 $fontno LM $schi1 \nEOF\n";
  print CSH "pstext -N $J_title $R_title -K -O -V >>$psfile<<EOF\n $x_title $z_title $fsize_title 0 $fontno CM $title \nEOF\n";
#-----------------------------
  print CSH "pstext -N $J_title $R_title -O -V >>$psfile<<EOF\n 10 10 $fsize_title 0 $fontno CM $title \nEOF\n"; # FINISH
  print CSH "convert $psfile -rotate 90 $jpgfile\n";
  if($ieps==1){print CSH "convert $psfile $epsfile\n"}

  close (CSH); system("csh -f $cshfile");
  system("xv $jpgfile &");
}

#===========================================================================


#=================================================================
