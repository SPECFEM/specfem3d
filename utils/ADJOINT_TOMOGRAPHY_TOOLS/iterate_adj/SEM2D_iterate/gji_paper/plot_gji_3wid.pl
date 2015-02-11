#!/usr/bin/perl -w

#==========================================================
#
#  plot_gji_3wid.pl
#  Carl Tape
#  04-Nov-2005
#
#  First run wave2d.f90, then model_optimize_poly.m to get axes limits, etc.
#
#  ISTRUCTURE:
#    1. checker, Nfac=3
#    2. checker, Nfac=2
#    3. checker, Nfac=1
#    4. rayleigh wave (T=20s, smoothed with sig=10km)
#
#    plot_gji_3wid.pl 25 60000 -6/3.0/0/80/08/03 1 0460 0  8 1 1 1 0   # big checker: ifig01,2,3,4
#    plot_gji_3wid.pl 25 30000 -6/3.0/0/80/10/05 1 0780 0  8 1 1 3 0   # sm  checker: ifig02,4
#    plot_gji_3wid.pl 25 15000 -6/3.0/0/80/20/10 1 4000 0 16 1 1 3 0   # sm  checker, 20 percent: ifig04
#    plot_gji_3wid.pl 25 30000 -6/3.0/0/20/08/03 1 0280 0  8 0 1 4 0   # rayleigh   : ifig02,4
#
#    plot_gji_3wid.pl 25 30000 -6/3.0/0/20/08/03 1 0100 0  0 1 4 0 0   # rayleigh: ifig04
#    plot_gji_3wid.pl 25 30000 -6/3.0/0/20/08/03 1 0140 0  0 1 4 0 0   # rayleigh: ifig04
#
#    plot_gji_3wid.pl 25 60000 -6/3.0/0/80/08/03 1 2500 0  1 1 2 0 0
#
#    plot_gji_3wid.pl 25 90000 -7/6.0/0/80/08/03 1 0500 1  1 0 1 0 0   # ifig04b (kernels only)
#    plot_gji_3wid.pl 25 75000 -7/6.0/0/80/08/03 1 0480 1  1 0 1 0 0
#    plot_gji_3wid.pl 25 60000 -7/6.0/0/80/08/03 1 0460 1  1 0 1 0 0
#    plot_gji_3wid.pl 25 45000 -7/6.0/0/80/08/03 1 0440 1  1 0 1 0 0
#    plot_gji_3wid.pl 25 30000 -7/6.0/0/80/08/03 1 0420 1  1 0 1 0 0
#    plot_gji_3wid.pl 25 15000 -7/6.0/0/80/08/03 1 0400 1  1 0 1 0 0
#
#    plot_gji_3wid.pl 25 60000 -7/3.0/0/20/08/03 1 0460 0  8 1 1 1 0    # big checker, 04b
#    plot_gji_3wid.pl 25 30000 -7/3.0/0/80/08/03 1 0780 0  8 1 1 3 0    # sm checker, 04b
#    plot_gji_3wid.pl 25 30000 -7/1.0/0/20/08/03 1 0280 0  8 0 1 4 0    # rayleigh, 04b
#
#    plot_gji_3wid.pl 25 60000 -6/2.0/0/20/08/03 1 0460 0  8 1 0 1 0   # ifig05
#    plot_gji_3wid.pl 25 60000 -6/1.5/0/20/08/03 1 0460 0  8 1 0 1 0   # ifig07 -- big checker
#    plot_gji_3wid.pl 25 60000 -6/3.0/0/20/08/03 1 0280 0  8 0 0 4 0   # ifig07 -- rayleigh
#    plot_gji_3wid.pl 25 60000 -7/7.0/0/20/08/03 1 0460 0  8 1 0 1 0   # 6,8,8b,9
#
#    plot_gji_3wid.pl 25 30000 -7/1.0/0/20/08/03 1 0280 0  8 0 0 4 0   # rayleigh, ifig10
#    plot_gji_3wid.pl 25 30000 -7/1.0/0/20/08/03 1 0520 0  8 0 0 4 0   # rayleigh, ifig10 (pert source)
#    plot_gji_3wid.pl 25 30000 -7/1.0/0/20/08/03 1 0540 0  8 0 0 4 0   # rayleigh, ifig10 (pert dT 50%)
#    plot_gji_3wid.pl 25 30000 -7/1.0/0/20/08/03 1 0560 0  8 0 0 4 0   # rayleigh, ifig10 (pert dT 20%)
#
#    plot_gji_3wid.pl 25 30000 -7/1.0/0/20/08/03 1 0280 0  8 0 1 4 0   # rayleigh, ifig10
#    plot_gji_3wid.pl 25 30000 -7/1.0/0/20/08/03 1 0100 0 16 0 1 4 0   # rayleigh, ifig10 -- with CG cubic-quad
#    plot_gji_3wid.pl 25 30000 -7/1.0/0/20/08/03 1 0100 0 16 0 1 4 0   # rayleigh, ifig10b
#
#    plot_gji_3wid.pl 25 30000 -7/1.0/0/20/09/03 1 0005 0  0 1 0 1 1   # ifig20
#    plot_gji_3wid.pl 25 30000 -7/1.0/0/20/10/05 1 0010 0  0 1 0 3 1   # ifig20, 07/24/06
#    plot_gji_3wid.pl 25 30000 -7/1.0/0/20/20/10 1 0011 0  0 1 0 3 1   # ifig20, 07/24/06
#    plot_gji_3wid.pl 25 30000 -7/1.0/0/20/09/03 1 0020 0  0 0 0 4 1   # ifig20, 07/19/06
#    plot_gji_3wid.pl 25 30000 -7/1.0/0/20/09/03 1 0021 0  0 0 0 4 1   # ifig20
#
#    plot_gji_3wid.pl 25 30000 -7/1.0/0/20/0.8 1  20 0  0 0 0 4 1   # ifig20b
#
#==========================================================

if (@ARGV < 8) {die("Usage: plot_kernels.pl nevent Gamma colors iker irun0 iter niter ichecker ipoly istructure ibanana \n");}
($nevent,$Gamma,$colors,$iker,$irun0,$iter,$niter,$ichecker,$ipoly,$istructure,$ibanana) = @ARGV;

if($ibanana==1) {$odir      = "../../../2d_adjoint_banana/OUTPUT_banana/run_"}
else            {$odir      = "../../OUTPUT/run_"}

$edir      = "event_001";  # event for first event (sr.txt)

$mfile_dat = "socal_vel_dat.dat";
$mfile_syn = "socal_vel_syn.dat";
$kerfile   = "kernel_basis";

$cshfile = "plot_gji_3wid.csh";

if($istructure==1){$Nfac=3}
if($istructure==2){$Nfac=2}
if($istructure==3){$Nfac=1}
$sNfac = sprintf("%1i",$Nfac);

# boolean commands for plotting
$icolor = 1;   # ccc

$ifig01  = 0;   # original CG 2-figure, 3x3 subplots
$ifig02  = 0;   # modified CG 1-figure, 3x3 subplots -- THIS MUST NOT BE RUN AFTER ifig01
$ifig03  = 0;   # phase vel models (6), 3x3 subplots
$ifig04  = 0;   # phase vel models (8), 3x3 subplots
$ifig04b = 0;   # kernels (9), 3x3 subplots
$ifig05  = 0;   # smoothing process for kernels
$ifig06  = 0;   # smoothing vs scalelength, 4x3 subplots
$ifig07  = 1;   # event kernels --> summed event kernel, 3x3 subplots
$ifig08  = 0;   # maps    : map vs Nevent, 4x3 subplots
$ifig08b = 0;   # maps    : map vs Nevent, 3x1 subplots
$ifig09  = 0;   # kernels : map vs Nevent, 4x3 subplots
$ifig10  = 0;
$ifig10b = 0;

# L-curve plots with recovered models
$ifig20  = 0;   # recovered images for classical tomography (hessian)
$ifig20b = 0;   # paper

$irun = $irun0 + $iter;
@mods = ("0","0t","1","1t","2","2t","3","3t","4","4t","5","5t","6","6t","7","7t","8","8t","9","9t","10","10t","11","11t","12","12t","13","13t","14","14t","15","15t","16","16t");
$mod = $mods[$iter];
$smod = "m\@+$mod\@+";

# label for the type of map used to generate the data
$smap = sprintf("m%2.2i",$istructure);

# wave2d run number
$strun0 = sprintf("%4.4i",$irun0);
$strun = sprintf("%4.4i",$irun);
$dir0 = "$odir$strun0";
$dir = "$odir$strun";

$stGam  = sprintf("\@~\107\@~ = %.1f km",$Gamma/1000);
$stGam2 = sprintf("%3.3i",$Gamma/1000);

print "\n $dir, $edir, $mfile_dat, $mfile_syn, $kerfile";
print "\n $colors, $iker, $mod, $stGam2 \n";
#die("testing");

# colors for the kernel and misfit function
@cols = split("/",$colors);
$kpwr  = $cols[0];  # power for kernels
$kmax  = $cols[1];  # value for kernels
$opwr  = $cols[2];  # power for misfit maps
$omax  = $cols[3];  # value for misfit maps
$cmax  = $cols[4];  # value for phase speed maps (percent pert)
$ctick = $cols[5];  # value for phase speed maps (percent pert tick)

#@files = glob("$dir/$kerfile");
#$numk = @files;

@titles = ("Waveform","Traveltime (xcorr), misfit","Amplitude (xcorr), misfit","Traveltime (MT), misfit","Amplitude (MT), misfit","Traveltime (xcorr), sampling","Amplitude (xcorr), sampling");
@units = ("m\@+2\@+ s","s\@+2\@+","xxx","xxx","xxx","xxx","xxx");
$ktype = $titles[$iker];
$utype = $units[$iker];

$plabel = "/home/carltape/2d_adjoint/scripts/plot_ker_mod.pl";

# data files
$dir_gmt = "/home/carltape/gmt";
$ishelf = 0;
$shelf_file = "../../INPUT/oms_shelf";
$plate_file = "${dir_gmt}/plate_boundaries";
$fault_file = "${dir_gmt}/faults/jennings.xy";

# plotting specifications
$fsize0 = "14";
$fsize1 = "11";
$fsize2 = "8";
$fsize3 = "5";
$fontno = "4";    # 1 or 4
$tick   = "0.1c";
$fpen   = "1.5p";
$tpen   = "1.0p";

# plot symbols for sources, receivers, and shelf
$src    = "-W0.5p -Sa8p";
$rec    = "-Sc2p -G0";
$rec0   = "-Sc10p -W0.5p";
$src_ev = "-W0.5p -Sa12p -G255";
$Wshelf = "-W0.5/0/0/0tap";

$coast_info = "-W1p -Na/1p -Df";
$coast_info2 = "-W0.5p -Na/0.5p -Df";

$plate_info_w = "-M -W2.5p/255/255/255";
$plate_info_r = "-M -W2.5p/255/0/0";
$plate_info_b = "-M -W2.5p/0/0/255";
$plate_info_k = "-M -W2.5p";

$fault_info_r = "-M -W1.5p/255/0/0";
$fault_info_k = "-M -W1.5p/0/0/0";
$fault_info_w = "-M -W1.5p/255/255/255";

#-------------------------
# color for kernels

# KEY: scaling for color
$scale_color = 21.0;
$colorbar = "seis";
#open(IN,"$dir/$cfile");
#@colors = <IN>;

$norm2 = "1e$kpwr";
$ss = $kmax;
$ds = 2*$ss/$scale_color;
$bs2 = sprintf("%3.3e",0.9*$ss); # colorbar
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

open(IN,"$dir0/socal_vel_c0.dat");
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
$Jwid = 2;
$J = "-JM${Jwid}i";      # in lat-lon
$origin = "-X0.8 -Y8.25";

# which borders to plot the lat-lon
# 1 four sides, 4 single sides, 6 two sides, 4 three sides, 1 zero sides
@Bopts = ("WESN","Wesn","wesN","wEsn","weSn","WesN","wEsN","wESn","WeSn","WEsn","weSN","wESN","WESn","WeSN","WEsN","wesn");
$B0 = "-B1:.\" \":";

# axes scale for phase speed maps: c(th, ph)
#$bs1 = 0.5;
#$Bscale1d  = sprintf("-B%2.2e:\" Phase Speed for data ( km s\@+-1\@+ )\":",$bs1);
#$Bscale1s  = sprintf("-B%2.2e:\" Phase Speed for model $smod ( km s\@+-1\@+ )\":",$bs1);
$bs1 = $ctick;
$Bscale1  = sprintf("-B%2.2ef1:\" \@~\045\@~ pert. from %2.2f km/s\": -E7p",$bs1,$c0/1000);
$Bscale1b = sprintf("-B%2.2ef1:\" \": -E7p",$bs1);

# axes scale for kernels: K(th, ph)
# [\@~\143\@~] --> s
$tp = "\@~\146\@~, \@~\161\@~";
$Bscale2  = sprintf("-B%2.2e:\" K ( $tp )  ( 10\@+%2.2i\@+  m\@+-2\@+ s )\": -E7p",$bs2,$kpwr);
$Bscale2b = sprintf("-B%2.2e:\" \": -E7p",$bs2);

# axes scale for chi_plots: chi(th_r, ph_r)
$Bscale3  = sprintf("-B%2.2e:\" \@~\143\@~ ( \@~\161\@~\@-r\@- , \@~\146\@~\@-r\@- )  ( 10\@+%2.2i\@+ )\": -Ef7p",$bs3,$opwr);

#-------------------------
# phase speed model
$file1dat = "$dir/$mfile_dat";
$file1syn = "$dir/$mfile_syn";

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
$fsize_title = 9;

#=================================================================
# files and parameters related to polynomial plots and misfit-vs-iteration plots

  $p_info_k   = "-N -Sc6p -W0.5p -G0";
  $p_info_w   = "-N -Sc6p -W0.5p -G255";
  $c_info_ks  = "-W0.5p";
  $c_info_ks2 = "-W1.5p";
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
  $J_poly = "-JX$Jwid";

  # polynomial lines and curves
  $poly_curve = "poly_curve_${strun0}.dat";
  $poly_pts   = "poly_points_${strun0}.dat";
  if (not -f ${poly_curve}) { die("Check if ${poly_curve} exist or not\n") }
  if (not -f ${poly_pts})   { die("Check if ${poly_pts} exist or not\n") }
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
  $chi_ran = $xmaxg - $xming;
  $R_chi = "-R$xming/$xmaxg/$yming/$ymaxg";
  $J_chi = "-JX${Jwid}i/${Jwid}il";        # note log scale on y-axis

  # text labels for chi-vs-m plots
  $schi0 = "\@~\143\@~(m\@+0\@+)";
  $schi1 = "\@~\143\@~(m\@+1\@+)";
  $x0chit = 0 + 0.05*$chi_ran;
  $x1chit = 1 + 0.05*$chi_ran;

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
  $B1 = sprintf("-B${xtick}:\"\@~\156\@~   ( 10\@+%2.2i\@+ ) \":/${ytick}:\" \@~\143\@~\@+0\@+ [ m (\@~\156\@~) ]   ( 10\@+%2.2i\@+  $utype ) \":",$lpwr_x,$lpwr_y);
  $B2 = sprintf("-B${xtick}:\"\@~\156\@~   ( 10\@+%2.2i\@+ ) \":/${ytick}:\" \@~\143\@~\@+1\@+ [ m (\@~\156\@~) ]   ( 10\@+%2.2i\@+  $utype ) \":",$lpwr_x,$lpwr_y);

  # scale for chi-vs-m plots
  $iter_tick = 2;
  $B3a    = "-Ba${iter_tick}f1:\" k, model number \":/a1f2g1p:\" \@~\143\@~ (m)   ( $utype ) \":";
  $B3b    = "-Ba${iter_tick}f1:\" k, model number \":/a1f2g1p:\" \":";
  $B4     = "-Ba${iter_tick}f1:\" k, model number \":/20:\" \":";

}

#===========================================================================
# create colorpoint files

  open(CSH,">$cshfile");

  $cpt_vel_map = "../../model_files/socal_color.cpt";

  # make uniform colorpoint file
  #$T1 = "-T3/4/0.1";
  $dc = $cmax/10;
  $T1 = "-T-$cmax/$cmax/$dc";
  $cpt_vel_uniform = "color0.cpt";
  print CSH "makecpt -C$colorbar $T1 > temp1\n";
  print CSH "sed 's/^B.*/B       170     0      0  /' temp1 >  temp2\n";
  print CSH "sed 's/^F.*/F         0     0    255  /' temp2 > $cpt_vel_uniform\n";

  # phase speed model
  if($ichecker==0) {$cpt_vel = $cpt_vel_map}
  else             {$cpt_vel = $cpt_vel_uniform}

  # kernel
  $cpt_ker = "color1.cpt";
  print CSH "makecpt -C$colorbar $T2 > temp1\n";
  print CSH "sed 's/^B.*/B       170     0      0  /' temp1 >  temp2\n";
  print CSH "sed 's/^F.*/F         0     0    255  /' temp2 > $cpt_ker\n";

  # misfit (as a function of receiver)
  $cpt_chi = "color2.cpt";
  print CSH "makecpt -Chot $T3 -I > $cpt_chi\n";

  close (CSH);
  system("csh -f $cshfile");
  #die("testing");

#===========================================================================
if ($ifig01 == 1) {

  # shifting the subplots
  $xfac = 1.25;
  $yfac = 1.6;
  $dX = $xfac*$Jwid; $dY = 0; $shift1 = "-X$dX -Y$dY";
  $dX = -2*$dX; $dY = -$yfac*$Jwid; $shift2 = "-X$dX -Y$dY";

  # colorbar
  $Dlen = 1.2; $Dx = $Jwid/2; $Dy = -0.10*$Jwid;
  $Dscale = "-D$Dx/$Dy/$Dlen/0.10h";

  #===============================================
  print "\nWriting CSH file...\n";

  open(CSH,">$cshfile");
  print CSH "gmtset BASEMAP_TYPE plain PLOT_DEGREE_FORMAT D TICK_LENGTH $tick LABEL_FONT_SIZE $fsize2 ANOT_FONT_SIZE $fsize2  HEADER_FONT $fontno ANOT_FONT $fontno LABEL_FONT $fontno HEADER_FONT_SIZE $fsize1 FRAME_PEN $fpen TICK_PEN $tpen \n";
  #===============================================

  # get the summed kernel and chi files
  $file2 = "$dir/summed_ker.dat";
  $file3 = "$dir/summed_chi_r.dat";
  $recfile = "$dir/${edir}/sr.txt";            # src-rec for first event
  $evefile = "$dir0/events_lonlat.dat";
  $chifile = "$dir/summed_chi_all.dat";

  if (not -f $file2)   { die("Check if $file2 exist or not\n") }
  if (not -f $file3)   { die("Check if $file3 exist or not\n") }
  if (not -f $recfile) { die("Check if $recfile exist or not\n") }
  if (not -f $chifile) { die("Check if $chifile exist or not\n") }

  # file names for figures
  $name    = "${smap}_cg01a_${strun0}";
  $psfile  = "$name.ps";
  $jpgfile = "$name.jpg";

  # number of receivers
  open(IN,$file3); @temp = <IN>; $nrec = @temp;
  # $npath = $nrec * $nevent;

  #=============================================
  # model for the data

  $B = "$B0".$Bopts[15];
  $title = "Target phase speed model (data)";

  # phase speed map
  print CSH "psbasemap $B $R $J -P -K -V $origin > $psfile\n"; # START
  #if($icolor==1) {print CSH "awk '{print \$1,\$2,\$3/1000}' $file1dat | pscontour $R $J -A- -C$cpt_vel -I -P -K -O -V >> $psfile\n"}
  if($icolor==1) {print CSH "awk '{print \$1,\$2,\$4*100}' $file1dat | pscontour $R $J -A- -C$cpt_vel -I -P -K -O -V  >> $psfile\n"}
  print CSH "pscoast $J $R $coast_info -P -K -O -V >> $psfile\n";
  if($ishelf) {print CSH "awk '{print \$2,\$1}' $shelf_file |psxy $J $R $Wshelf -K -O -P -V >> $psfile\n"}
  print CSH "awk '\$1 == \"R\" {print \$2,\$3}' $recfile |psxy -N $J $R -K -O -P -V $rec >> $psfile\n";
  print CSH "awk '{print \$1,\$2}' $evefile |psxy -N $J $R -K -O -P -V $src >> $psfile\n";
  print CSH "psscale -C$cpt_vel $Dscale $Bscale1 -P -K -O -V >> $psfile \n";
  print CSH "pstext -N $J_title $R_title -K -O -V -P >>$psfile<<EOF\n $x_title $z_title $fsize_title 0 $fontno CM $title \nEOF\n";

  #=============================================
  # model for synthetics

  $title = "Phase speed model $smod";
  $shift = $shift1;

  # phase speed map
  print CSH "psbasemap $B $R $J -P -K -O -V $shift >> $psfile\n";
  #if($icolor==1) {print CSH "awk '{print \$1,\$2,\$3/1000}' $file1syn | pscontour $R $J -A- -C$cpt_vel -I -P -O -K -V >> $psfile\n"}
  if($icolor==1) {print CSH "awk '{print \$1,\$2,\$4*100}' $file1syn | pscontour $R $J -A- -C$cpt_vel -I -P -O -K -V >> $psfile\n"}
  print CSH "pscoast $J $R $coast_info -P -K -O -V >> $psfile\n";
  if($ishelf) {print CSH "awk '{print \$2,\$1}' $shelf_file |psxy $J $R $Wshelf -K -O -P -V >> $psfile\n"}
  print CSH "awk '{print \$1,\$2}' $evefile |psxy -N $J $R -K -O -P -V $src >> $psfile\n";
  print CSH "awk '\$1 == \"R\" {print \$2,\$3}' $recfile |psxy -N $J $R -K -O -P -V $rec >> $psfile\n";
  print CSH "psscale -C$cpt_vel $Dscale $Bscale1 -P -K -O -V >> $psfile \n";
  print CSH "pstext -N $J_title $R_title -K -O -V -P >>$psfile<<EOF\n $x_title $z_title $fsize_title 0 $fontno CM $title \nEOF\n";

  #=============================================
  # chi as a function of receiver

  # total misfit for the model
  open(IN,"$chifile");
  $chi = <IN>;
  $schi = sprintf("\@~\143\@~ ( $smod )  =  %3.3e",$chi);
  $title = "Summed misfit : $schi";
  $B = "$B0".$Bopts[15];
  $shift = $shift1;

  print CSH "psbasemap $B $R $J -P -K -O -V $shift >> $psfile\n";
  #print CSH "awk '{print \$1,\$2,\$3 / $norm3 }' $file > temp\n";
  print CSH "awk '{print \$1,\$2,\$3 / $norm3 }' $file3 | pscontour $J $R -A- -C$cpt_chi -I -K -O -P -V >> $psfile\n";
  print CSH "pscoast $J $R $coast_info -K -O -P -V >> $psfile\n";
  if($ishelf) {print CSH "awk '{print \$2,\$1}' $shelf_file |psxy $J $R $Wshelf -K -O -P -V >> $psfile\n"}
  print CSH "psscale -C$cpt_chi $Dscale $Bscale3 -K -O -P -V >> $psfile \n";
  print CSH "awk '{print \$1,\$2}' $evefile |psxy -N $J $R -K -O -P -V $src >> $psfile\n";
  print CSH "awk '\$1 == \"R\" {print \$2,\$3}' $recfile |psxy -N $J $R -K -O -P -V $rec >> $psfile\n";
  print CSH "pstext -N $J_title $R_title -K -O -V -P >>$psfile<<EOF\n $x_title $z_title $fsize_title 0 $fontno CM $title \nEOF\n";

  #=============================================
  # graph

  $shift = $shift2;
  $title = "Plot  $schi";
  $Bg = "$B1".$Bopts[8];

  print CSH "psbasemap $Bg $R_poly $J_poly -P -K -O -V $shift >> $psfile\n";
  print CSH "psxy $J_poly $R_poly $p_info_k -K -O -V -P >>$psfile<<EOF\n$x1 $y1 \nEOF\n";
  print CSH "psxy $J_poly $R_poly $c_info_ks -K -O -V -P >>$psfile<<EOF\n$x1 $a1y1 \n$x1 $y1 \n$a1x1 $y1 \nEOF\n";
  print CSH "pstext -N $J_title $R_title -K -O -V -P >>$psfile<<EOF\n $x_title $z_title $fsize_title 0 $fontno CM $title \nEOF\n";

  #=============================================
  # kernel

  $shift = $shift1;
  $title = "Gradient for $smod";
  $B = "$B0".$Bopts[15];

  print CSH "psbasemap $B $R $J -P -K -O -V $shift >> $psfile\n";
  if($icolor==1) {print CSH "awk '{print \$1,\$2,\$3 / $norm2 }' $file2 | pscontour $J $R -A- -C$cpt_ker -I -K -O -P -V >> $psfile\n"}
  print CSH "pscoast $J $R $coast_info -K -O -P -V >> $psfile\n";
  if($ishelf) {print CSH "awk '{print \$2,\$1}' $shelf_file |psxy $J $R $Wshelf -K -O -P -V >> $psfile\n"}
  print CSH "psscale -C$cpt_ker $Dscale $Bscale2 -K -O -P -V >> $psfile \n";
  print CSH "awk '{print \$1,\$2}' $evefile |psxy -N $J $R -K -O -P -V $src >> $psfile\n";
  print CSH "pstext -N $J_title $R_title -K -O -V -P >>$psfile<<EOF\n $x_title $z_title $fsize_title 0 $fontno CM $title \nEOF\n";
  print CSH "awk '\$1 == \"R\" {print \$2,\$3}' $recfile |psxy -N $J $R -K -O -P -V $rec >> $psfile\n";

  #=============================================
  # graph

  $shift = $shift1;
  $Bg = "$B1".$Bopts[8];
  $title = "Plot  g\@+0\@+ = d\@~\143\@~ / dm($smod)";

  print CSH "psbasemap $Bg $R_poly $J_poly -P -K -O -V $shift >> $psfile\n";
  print CSH "awk '{print \$1 / $normgx,\$2 / $normgy}' $poly_curve | psxy $c_info_rd $J_poly $R_poly -K -O -P -V >> $psfile\n";
  print CSH "psxy $J_poly $R_poly $c_info_ks -K -O -V -P >>$psfile<<EOF\n$x1 $a1y1 \n$x1 $y1 \n$a1x1 $y1 \nEOF\n";
  print CSH "psxy $J_poly $R_poly $p_info_k -K -O -V -P >>$psfile<<EOF\n$x1 $y1 \nEOF\n";
  print CSH "pstext -N $J_title $R_title -K -O -V -P >>$psfile<<EOF\n $x_title $z_title $fsize_title 0 $fontno CM $title \nEOF\n";

  #=============================================
  # graph

  $shift = $shift2;
  $Bg = "$B1".$Bopts[8];
  $title = "Estimate test model m\@+0t\@+";

  print CSH "psbasemap $Bg $R_poly $J_poly -P -K -O -V $shift >> $psfile\n";
  print CSH "awk '{print \$1 / $normgx,\$2 / $normgy}' $poly_curve | psxy $c_info_rd $J_poly $R_poly -K -O -P -V >> $psfile\n";
  print CSH "awk '{print \$1 / $normgx,\$4 / $normgy}' $poly_curve | psxy $c_info_bd $J_poly $R_poly -K -O -P -V >> $psfile\n";
  print CSH "psxy $J_poly $R_poly $c_info_ks -K -O -V -P >>$psfile<<EOF\n$x1 $a1y1 \n$x1 $y1 \n$a1x1 $y1 \nEOF\n";
  print CSH "psxy $J_poly $R_poly $p_info_k -K -O -V -P >>$psfile<<EOF\n$x1 $y1 \nEOF\n";
  print CSH "psxy $J_poly $R_poly $p_info_w -K -O -V -P >>$psfile<<EOF\n$x2 $a1y1 \nEOF\n";
  print CSH "pstext -N $J_title $R_title -K -O -V -P >>$psfile<<EOF\n $x_title $z_title $fsize_title 0 $fontno CM $title \nEOF\n";

#=================================================================
# NEW MODEL

  # wave2d run number
  $strun = sprintf("%4.4i",$irun + 1);
  $dir = "$odir$strun";

  $mod = $mods[1];
  $smod = "m\@+$mod\@+";

  # get the summed kernel and chi files
  $file1syn = "$dir/$mfile_syn";
  $file2 = "$dir/summed_ker.dat";
  $file3 = "$dir/summed_chi_r.dat";
  $recfile = "$dir/${edir}/sr.txt";  # src-rec for first event
  $chifile = "$dir/summed_chi_all.dat";

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
  print CSH "psbasemap $B $R $J -P -K -O -V $shift >> $psfile\n";
  if($icolor==1) {print CSH "awk '{print \$1,\$2,\$4*100}' $file1syn | pscontour $R $J -A- -C$cpt_vel -I -P -O -K -V >> $psfile\n"}
  print CSH "pscoast $J $R $coast_info -P -K -O -V >> $psfile\n";
  if($ishelf) {print CSH "awk '{print \$2,\$1}' $shelf_file |psxy $J $R $Wshelf -K -O -P -V >> $psfile\n"}
  print CSH "awk '{print \$1,\$2}' $evefile |psxy -N $J $R -K -O -P -V $src >> $psfile\n";
  print CSH "awk '\$1 == \"R\" {print \$2,\$3}' $recfile |psxy -N $J $R -K -O -P -V $rec >> $psfile\n";
  print CSH "psscale -C$cpt_vel $Dscale $Bscale1 -P -K -O -V >> $psfile \n";
  print CSH "pstext -N $J_title $R_title -K -O -V -P >>$psfile<<EOF\n $x_title $z_title $fsize_title 0 $fontno CM $title \nEOF\n";

  #=============================================
  # chi as a function of receiver

  # total misfit for the model
  open(IN,"$chifile");
  $chi = <IN>;
  $schi = sprintf("\@~\143\@~ ( $smod )  =  %3.3e",$chi);
  $title = "Summed misfit : $schi";
  $B = "$B0".$Bopts[15];
  $shift = $shift1;

  print CSH "psbasemap $B $R $J -P -K -O -V $shift >> $psfile\n";
  #print CSH "awk '{print \$1,\$2,\$3 / $norm3 }' $file > temp\n";
  print CSH "awk '{print \$1,\$2,\$3 / $norm3 }' $file3 | pscontour $J $R -A- -C$cpt_chi -I -K -O -P -V >> $psfile\n";
  print CSH "pscoast $J $R $coast_info -K -O -P -V >> $psfile\n";
  if($ishelf) {print CSH "awk '{print \$2,\$1}' $shelf_file |psxy $J $R $Wshelf -K -O -P -V >> $psfile\n"}
  print CSH "psscale -C$cpt_chi $Dscale $Bscale3 -K -O -P -V >> $psfile \n";
  print CSH "awk '{print \$1,\$2}' $evefile |psxy -N $J $R -K -O -P -V $src >> $psfile\n";
  print CSH "awk '\$1 == \"R\" {print \$2,\$3}' $recfile |psxy -N $J $R -K -O -P -V $rec >> $psfile\n";
  print CSH "pstext -N $J_title $R_title -K -O -V -P >>$psfile<<EOF\n $x_title $z_title $fsize_title 0 $fontno CM $title \nEOF\n";

#-----------------------------
  print CSH "pstext -N $J_title $R_title -O -V -P >>$psfile<<EOF\n 10 10 $fsize_title 0 $fontno CM $title \nEOF\n"; # FINISH
  print CSH "convert $psfile $jpgfile\n";
  print CSH "echo done with $psfile\n";

  close (CSH);
  system("csh -f $cshfile");
  system("xv $jpgfile &");
}

#=================================================================
# SECOND FIGURE
#=================================================================

#===========================================================================
if ($ifig01 == 1) {

  #===============================================
  print "\nWriting CSH file...\n";
  open(CSH,">$cshfile");
  print CSH "gmtset BASEMAP_TYPE plain TICK_LENGTH $tick LABEL_FONT_SIZE $fsize2 ANOT_FONT_SIZE $fsize2 PLOT_DEGREE_FORMAT D HEADER_FONT $fontno ANOT_FONT $fontno LABEL_FONT $fontno HEADER_FONT_SIZE $fsize1 FRAME_PEN $fpen TICK_PEN $tpen\n";
  #===============================================

  # file names for figures
  $name    = "${smap}_cg01b_${strun0}";
  $psfile  = "$name.ps";
  $jpgfile = "$name.jpg";

  #=============================================
  # graph

  $Bg = "$B1".$Bopts[8];
  $title = "Plot  $schi";

  print CSH "psbasemap $Bg $R_poly $J_poly -P -K -V $origin > $psfile\n";  # START
  print CSH "awk '{print \$1 / $normgx,\$2 / $normgy}' $poly_curve | psxy $c_info_rd $J_poly $R_poly -K -O -P -V >> $psfile\n";
  print CSH "awk '{print \$1 / $normgx,\$4 / $normgy}' $poly_curve | psxy $c_info_bd $J_poly $R_poly -K -O -P -V >> $psfile\n";
  print CSH "psxy $J_poly $R_poly $c_info_ks -K -O -V -P >>$psfile<<EOF\n$x1 $a1y1 \n$x1 $y1 \n$a1x1 $y1 \nEOF\n";
  print CSH "psxy $J_poly $R_poly $c_info_ks -K -O -V -P >>$psfile<<EOF\n$x3 $a1y1 \n$x3 $y3 \n$a1x1 $y3 \nEOF\n";
  print CSH "psxy $J_poly $R_poly $p_info_k -K -O -V -P >>$psfile<<EOF\n$x1 $y1 \n$x3 $y3\nEOF\n";
  print CSH "psxy $J_poly $R_poly $p_info_w -K -O -V -P >>$psfile<<EOF\n$x2 $a1y1 \nEOF\n";
  print CSH "pstext -N $J_title $R_title -K -O -V -P >>$psfile<<EOF\n $x_title $z_title $fsize_title 0 $fontno CM $title \nEOF\n";

  #=============================================
  # kernel

  $shift = $shift1;
  $title = "Gradient for $smod";
  $B = "$B0".$Bopts[15];

  print CSH "psbasemap $B $R $J -P -K -O -V $shift >> $psfile\n";
  if($icolor==1) {print CSH "awk '{print \$1,\$2,\$3 / $norm2 }' $file2 | pscontour $J $R -A- -C$cpt_ker -I -K -O -P -V >> $psfile\n"}
  print CSH "pscoast $J $R $coast_info -K -O -P -V >> $psfile\n";
  if($ishelf) {print CSH "awk '{print \$2,\$1}' $shelf_file |psxy $J $R $Wshelf -K -O -P -V >> $psfile\n"}
  print CSH "psscale -C$cpt_ker $Dscale $Bscale2 -K -O -P -V >> $psfile \n";
  print CSH "awk '{print \$1,\$2}' $evefile |psxy -N $J $R -K -O -P -V $src >> $psfile\n";
  print CSH "pstext -N $J_title $R_title -K -O -V -P >>$psfile<<EOF\n $x_title $z_title $fsize_title 0 $fontno CM $title \nEOF\n";
  print CSH "awk '\$1 == \"R\" {print \$2,\$3}' $recfile |psxy -N $J $R -K -O -P -V $rec >> $psfile\n";

  #=============================================
  # graph

  $shift = $shift1;
  $Bg = "$B1".$Bopts[8];
  $title = "Plot  g\@+0t\@+ = d\@~\143\@~ / dm($smod)";

  print CSH "psbasemap $Bg $R_poly $J_poly -P -K -O -V $shift >> $psfile\n";
  print CSH "awk '{print \$1 / $normgx,\$2 / $normgy}' $poly_curve | psxy $c_info_rd $J_poly $R_poly -K -O -P -V >> $psfile\n";
  print CSH "awk '{print \$1 / $normgx,\$3 / $normgy}' $poly_curve | psxy $c_info_rd $J_poly $R_poly -K -O -P -V >> $psfile\n";
  #print CSH "awk '{print \$1 / $normgx,\$4 / $normgy}' $poly_curve | psxy $c_info_bd $J_poly $R_poly -K -O -P -V >> $psfile\n";
  print CSH "psxy $J_poly $R_poly $c_info_ks -K -O -V -P >>$psfile<<EOF\n$x1 $a1y1 \n$x1 $y1 \n$a1x1 $y1 \nEOF\n";
  print CSH "psxy $J_poly $R_poly $c_info_ks -K -O -V -P >>$psfile<<EOF\n$x3 $a1y1 \n$x3 $y3 \n$a1x1 $y3 \nEOF\n";
  print CSH "psxy $J_poly $R_poly $p_info_k -K -O -V -P >>$psfile<<EOF\n$x1 $y1 \n$x3 $y3\nEOF\n";
  print CSH "psxy $J_poly $R_poly $p_info_w -K -O -V -P >>$psfile<<EOF\n$x2 $a1y1 \nEOF\n";
  print CSH "pstext -N $J_title $R_title -K -O -V -P >>$psfile<<EOF\n $x_title $z_title $fsize_title 0 $fontno CM $title \nEOF\n";

  #=============================================
  # graph

  $shift = $shift2;
  $Bg = "$B1".$Bopts[8];
  $title = "Cubic interpolation to get m\@+1\@+";

  print CSH "psbasemap $Bg $R_poly $J_poly -P -K -O -V $shift >> $psfile\n";
  print CSH "awk '{print \$1 / $normgx,\$2 / $normgy}' $poly_curve | psxy $c_info_rd $J_poly $R_poly -K -O -P -V >> $psfile\n";
  print CSH "awk '{print \$1 / $normgx,\$3 / $normgy}' $poly_curve | psxy $c_info_rd $J_poly $R_poly -K -O -P -V >> $psfile\n";
  print CSH "awk '{print \$1 / $normgx,\$5 / $normgy}' $poly_curve | psxy $c_info_ks2 $J_poly $R_poly -K -O -P -V >> $psfile\n";
  print CSH "psxy $J_poly $R_poly $c_info_ks -K -O -V -P >>$psfile<<EOF\n$x1 $a1y1 \n$x1 $y1 \n$a1x1 $y1 \nEOF\n";
  print CSH "psxy $J_poly $R_poly $c_info_ks -K -O -V -P >>$psfile<<EOF\n$x3 $a1y1 \n$x3 $y3 \n$a1x1 $y3 \nEOF\n";
  print CSH "psxy $J_poly $R_poly $c_info_kd -K -O -V -P >>$psfile<<EOF\n$x5 $a1y1 \n$x5 $y5 \n$a1x1 $y5 \nEOF\n";
  print CSH "psxy $J_poly $R_poly $p_info_k -K -O -V -P >>$psfile<<EOF\n$x1 $y1 \n$x3 $y3\nEOF\n";
  print CSH "psxy $J_poly $R_poly $p_info_w -K -O -V -P >>$psfile<<EOF\n$x2 $a1y1 \n$x4 $y4 \n$x5 $y5 \nEOF\n";
  print CSH "pstext -N $J_title $R_title -K -O -V -P >>$psfile<<EOF\n $x_title $z_title $fsize_title 0 $fontno CM $title \nEOF\n";

#=================================================================
# NEW MODEL

  # wave2d run number
  $strun = sprintf("%4.4i",$irun + 2);
  $dir = "$odir$strun";

  $mod = $mods[2];
  $smod = "m\@+$mod\@+";

  # get the summed kernel and chi files
  $file1syn = "$dir/$mfile_syn";
  $file2 = "$dir/summed_ker.dat";
  $file3 = "$dir/summed_chi_r.dat";
  $recfile = "$dir/${edir}/sr.txt";  # src-rec for first event
  $chifile = "$dir/summed_chi_all.dat";

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
  print CSH "psbasemap $B $R $J -P -K -O -V $shift >> $psfile\n";
  if($icolor==1) {print CSH "awk '{print \$1,\$2,\$4*100}' $file1syn | pscontour $R $J -A- -C$cpt_vel -I -P -O -K -V >> $psfile\n"}
  print CSH "pscoast $J $R $coast_info -P -K -O -V >> $psfile\n";
  if($ishelf) {print CSH "awk '{print \$2,\$1}' $shelf_file |psxy $J $R $Wshelf -K -O -P -V >> $psfile\n"}
  print CSH "awk '{print \$1,\$2}' $evefile |psxy -N $J $R -K -O -P -V $src >> $psfile\n";
  print CSH "awk '\$1 == \"R\" {print \$2,\$3}' $recfile |psxy -N $J $R -K -O -P -V $rec >> $psfile\n";
  print CSH "psscale -C$cpt_vel $Dscale $Bscale1 -P -K -O -V >> $psfile \n";
  print CSH "pstext -N $J_title $R_title -K -O -V -P >>$psfile<<EOF\n $x_title $z_title $fsize_title 0 $fontno CM $title \nEOF\n";

  #=============================================
  # chi as a function of receiver

  # total misfit for the model
  open(IN,"$chifile");
  $chi = <IN>;
  $schi = sprintf("\@~\143\@~ ( $smod )  =  %3.3e",$chi);
  $title = "Summed misfit : $schi";
  $B = "$B0".$Bopts[15];
  $shift = $shift1;

  print CSH "psbasemap $B $R $J -P -K -O -V $shift >> $psfile\n";
  #print CSH "awk '{print \$1,\$2,\$3 / $norm3 }' $file > temp\n";
  print CSH "awk '{print \$1,\$2,\$3 / $norm3 }' $file3 | pscontour $J $R -A- -C$cpt_chi -I -K -O -P -V >> $psfile\n";
  print CSH "pscoast $J $R $coast_info -K -O -P -V >> $psfile\n";
  if($ishelf) {print CSH "awk '{print \$2,\$1}' $shelf_file |psxy $J $R $Wshelf -K -O -P -V >> $psfile\n"}
  print CSH "psscale -C$cpt_chi $Dscale $Bscale3 -K -O -P -V >> $psfile \n";
  print CSH "awk '{print \$1,\$2}' $evefile |psxy -N $J $R -K -O -P -V $src >> $psfile\n";
  print CSH "awk '\$1 == \"R\" {print \$2,\$3}' $recfile |psxy -N $J $R -K -O -P -V $rec >> $psfile\n";
  print CSH "pstext -N $J_title $R_title -K -O -V -P >>$psfile<<EOF\n $x_title $z_title $fsize_title 0 $fontno CM $title \nEOF\n";

  #=============================================
  # graph

  $shift = $shift2;
  $Bg = "$B1".$Bopts[8];
  $title = "Plot  $schi";

  print CSH "psbasemap $Bg $R_poly $J_poly -P -K -O -V $shift >> $psfile\n";
  print CSH "awk '{print \$1 / $normgx,\$2 / $normgy}' $poly_curve | psxy $c_info_rd $J_poly $R_poly -K -O -P -V >> $psfile\n";
  print CSH "awk '{print \$1 / $normgx,\$3 / $normgy}' $poly_curve | psxy $c_info_rd $J_poly $R_poly -K -O -P -V >> $psfile\n";
  print CSH "awk '{print \$1 / $normgx,\$5 / $normgy}' $poly_curve | psxy $c_info_ks2 $J_poly $R_poly -K -O -P -V >> $psfile\n";

  print CSH "psxy $J_poly $R_poly $c_info_ks -K -O -V -P >>$psfile<<EOF\n$x1 $a1y1 \n$x1 $y1 \n$a1x1 $y1 \nEOF\n";
  print CSH "psxy $J_poly $R_poly $c_info_ks -K -O -V -P >>$psfile<<EOF\n$x3 $a1y1 \n$x3 $y3 \n$a1x1 $y3 \nEOF\n";
  print CSH "psxy $J_poly $R_poly $c_info_kd -K -O -V -P >>$psfile<<EOF\n$x5 $a1y1 \n$x5 $y5 \n$a1x1 $y5 \nEOF\n";
  print CSH "psxy $J_poly $R_poly $c_info_ks -K -O -V -P >>$psfile<<EOF\n$x5 $a1y1 \n$x6 $y6 \n$a1x1 $y6 \nEOF\n";

  print CSH "psxy $J_poly $R_poly $p_info_w -K -O -V -P >>$psfile<<EOF\n$x2 $a1y1 \n$x4 $y4 \n$x5 $y5 \nEOF\n";
  print CSH "psxy $J_poly $R_poly $p_info_k -K -O -V -P >>$psfile<<EOF\n$x1 $y1 \n$x3 $y3 \n $x6 $y6 \nEOF\n";
  print CSH "pstext -N $J_title $R_title -K -O -V -P >>$psfile<<EOF\n $x_title $z_title $fsize_title 0 $fontno CM $title \nEOF\n";

  #=============================================
  # graph

  $shift = $shift1;
  $Bg = "$B2".$Bopts[8];
  $title = "New line search";

  print CSH "psbasemap $Bg $R_poly $J_poly -P -K -O -V $shift >> $psfile\n";

  print CSH "psxy $J_poly $R_poly $c_info_ks -K -O -V -P >>$psfile<<EOF\n$x1 $a1y1 \n$x1 $y6 \n$a1x1 $y6 \nEOF\n";
  print CSH "psxy $J_poly $R_poly $p_info_k -K -O -V -P >>$psfile<<EOF\n$x1 $y6 \nEOF\n";
  print CSH "pstext -N $J_title $R_title -K -O -V -P >>$psfile<<EOF\n $x_title $z_title $fsize_title 0 $fontno CM $title \nEOF\n";

  #=============================================
  # misfit vs iteration

  $shift = $shift1;
  $Bg = "$B3a".$Bopts[8];
  $title = "Misfit for first two models";

  print CSH "psbasemap $Bg $R_chi $J_chi -P -K -O -V $shift >> $psfile\n";
  print CSH "awk '{print \$1,\$2}' $chi_curve | psxy $c_info_rd $J_chi $R_chi -K -O -P -V >> $psfile\n";
  print CSH "psxy $J_chi $R_chi $p_info_k -K -O -V -P >>$psfile<<EOF\n0 $y1t \n1 $y7t \nEOF\n";
  print CSH "pstext -N $J_chi $R_chi -K -O -V -P >>$psfile<<EOF\n 0.75 $y1t $fsize1 0 $fontno LM $schi0 \nEOF\n";
  print CSH "pstext -N $J_chi $R_chi -K -O -V -P >>$psfile<<EOF\n 1.75 $y7t $fsize1 0 $fontno LM $schi1 \nEOF\n";
  print CSH "pstext -N $J_title $R_title -K -O -V -P >>$psfile<<EOF\n $x_title $z_title $fsize_title 0 $fontno CM $title \nEOF\n";

#-----------------------------
  print CSH "pstext -N $J_title $R_title -O -V -P >>$psfile<<EOF\n 10 10 $fsize_title 0 $fontno CM $title \nEOF\n"; # FINISH
  print CSH "convert $psfile $jpgfile\n";
  print CSH "echo done with $psfile\n";

  close (CSH);
  system("csh -f $cshfile");
  system("xv $jpgfile &");
}

#===========================================================================
#===========================================================================
if ($ifig02 == 1) {

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
  print CSH "gmtset BASEMAP_TYPE plain TICK_LENGTH $tick LABEL_FONT_SIZE $fsize2 ANOT_FONT_SIZE $fsize2 PLOT_DEGREE_FORMAT D HEADER_FONT $fontno ANOT_FONT $fontno LABEL_FONT $fontno HEADER_FONT_SIZE $fsize1 FRAME_PEN $fpen TICK_PEN $tpen\n";
  #===============================================

  # get the summed kernel and chi files
  #$file2 = "$dir/summed_ker.dat";
  $file2 = "$dir/fun_smooth.dat";       # smoothed kernel
  $file3 = "$dir/summed_chi_r.dat";
  $recfile = "$dir/${edir}/sr.txt";  # src-rec for first event
  $evefile = "$dir0/events_lonlat.dat";
  $chifile = "$dir/summed_chi_all.dat";

  if (not -f $file2)   { die("Check if $file2 exist or not\n") }
  if (not -f $file3)   { die("Check if $file3 exist or not\n") }
  if (not -f $recfile) { die("Check if $recfile exist or not\n") }
  if (not -f $chifile) { die("Check if $chifile exist or not\n") }

  # file names for figures
  $name    = "${smap}_cg02_${strun0}";
  $psfile  = "$name.ps";
  $jpgfile = "$name.jpg";

  # number of receivers
  open(IN,$file3); @temp = <IN>; $nrec = @temp;

  #=============================================
  # model for the data

  $B = "$B0".$Bopts[15];
  $title = "(a)  Target model (data)";

  # phase speed map
  print CSH "psbasemap $B $R $J -P -K -V $origin > $psfile\n"; # START
  #if($icolor==1) {print CSH "awk '{print \$1,\$2,\$3/1000}' $file1dat | pscontour $R $J -A- -C$cpt_vel -I -P -K -O -V >> $psfile\n"}
  if($icolor==1) {print CSH "awk '{print \$1,\$2,\$4*100}' $file1dat | pscontour $R $J -A- -C$cpt_vel -I -P -K -O -V  >> $psfile\n"}
  print CSH "pscoast $J $R $coast_info -P -K -O -V >> $psfile\n";
  if($ishelf) {print CSH "awk '{print \$2,\$1}' $shelf_file |psxy $J $R $Wshelf -K -O -P -V >> $psfile\n"}
  print CSH "awk '\$1 == \"R\" {print \$2,\$3}' $recfile |psxy -N $J $R -K -O -P -V $rec >> $psfile\n";
  print CSH "awk '{print \$1,\$2}' $evefile |psxy -N $J $R -K -O -P -V $src >> $psfile\n";
  print CSH "psscale -C$cpt_vel $Dscale $Bscale1b -P -K -O -V >> $psfile \n";
  print CSH "pstext -N $J_title $R_title -K -O -V -P >>$psfile<<EOF\n $x_title $z_title $fsize_title 0 $fontno CM $title \nEOF\n";

  #=============================================
  # model for synthetics

  $title = "(b)  Initial model $smod (synthetics)";
  $shift = $shift1;

  # phase speed map
  print CSH "psbasemap $B $R $J -P -K -O -V $shift >> $psfile\n";
  #if($icolor==1) {print CSH "awk '{print \$1,\$2,\$3/1000}' $file1syn | pscontour $R $J -A- -C$cpt_vel -I -P -O -K -V >> $psfile\n"}
  if($icolor==1) {print CSH "awk '{print \$1,\$2,\$4*100}' $file1syn | pscontour $R $J -A- -C$cpt_vel -I -P -O -K -V >> $psfile\n"}
  print CSH "pscoast $J $R $coast_info -P -K -O -V >> $psfile\n";
  if($ishelf) {print CSH "awk '{print \$2,\$1}' $shelf_file |psxy $J $R $Wshelf -K -O -P -V >> $psfile\n"}
  print CSH "awk '{print \$1,\$2}' $evefile |psxy -N $J $R -K -O -P -V $src >> $psfile\n";
  print CSH "awk '\$1 == \"R\" {print \$2,\$3}' $recfile |psxy -N $J $R -K -O -P -V $rec >> $psfile\n";
  #print CSH "psscale -C$cpt_vel $Dscale $Bscale1 -P -K -O -V >> $psfile \n";
  print CSH "pstext -N $J_title $R_title -K -O -V -P >>$psfile<<EOF\n $x_title $z_title $fsize_title 0 $fontno CM $title \nEOF\n";

  #=============================================
  # kernel -- SMOOTHED VERSION

  $shift = $shift1;
  $title = "(c)  Gradient (smoothed) for $smod";
  $B = "$B0".$Bopts[15];

  print CSH "psbasemap $B $R $J -P -K -O -V $shift >> $psfile\n";

  #if($icolor==1) {print CSH "awk '{print \$1,\$2,\$3 / $norm2 }' $file2 | pscontour $J $R -A- -C$cpt_ker -I -K -O -P -V >> $psfile\n"}
  if($icolor==1) {print CSH "awk '{print \$1,\$2,\$7 / $norm2 }' $file2 | pscontour $J $R -A- -C$cpt_ker -I -K -O -P -V >> $psfile\n"}

  print CSH "pscoast $J $R $coast_info -K -O -P -V >> $psfile\n";
  if($ishelf) {print CSH "awk '{print \$2,\$1}' $shelf_file |psxy $J $R $Wshelf -K -O -P -V >> $psfile\n"}
  print CSH "psscale -C$cpt_ker $Dscale $Bscale2b -K -O -P -V >> $psfile \n";
  print CSH "awk '{print \$1,\$2}' $evefile |psxy -N $J $R -K -O -P -V $src >> $psfile\n";
  print CSH "pstext -N $J_title $R_title -K -O -V -P >>$psfile<<EOF\n $x_title $z_title $fsize_title 0 $fontno CM $title \nEOF\n";
  print CSH "awk '\$1 == \"R\" {print \$2,\$3}' $recfile |psxy -N $J $R -K -O -P -V $rec >> $psfile\n";

  #=============================================
  # graph

  $shift = $shift2;
  $Bg = "$B1".$Bopts[8];
  $title = "(d)  Estimate test model m\@+0t\@+";

  print CSH "psbasemap $Bg $R_poly $J_poly -P -K -O -V $shift >> $psfile\n";
  print CSH "awk '{print \$1 / $normgx,\$2 / $normgy}' $poly_curve | psxy $c_info_rd $J_poly $R_poly -K -O -P -V >> $psfile\n";
  print CSH "awk '{print \$1 / $normgx,\$4 / $normgy}' $poly_curve | psxy $c_info_bd $J_poly $R_poly -K -O -P -V >> $psfile\n";
  print CSH "psxy $J_poly $R_poly $c_info_ks -K -O -V -P >>$psfile<<EOF\n$x1 $a1y1 \n$x1 $y1 \n$a1x1 $y1 \nEOF\n";
  print CSH "psxy $J_poly $R_poly $p_info_k -K -O -V -P >>$psfile<<EOF\n$x1 $y1 \nEOF\n";
  print CSH "psxy $J_poly $R_poly $p_info_w -K -O -V -P >>$psfile<<EOF\n$x2 $a1y1 \nEOF\n";
  print CSH "pstext -N $J_title $R_title -K -O -V -P >>$psfile<<EOF\n $x_title $z_title $fsize_title 0 $fontno CM $title \nEOF\n";

#=================================================================
# NEW MODEL

  # wave2d run number
  $strun = sprintf("%4.4i",$irun + 1);
  $dir = "$odir$strun";

  $mod = $mods[1];
  $smod = "m\@+$mod\@+";

  # get the summed kernel and chi files
  $file1syn = "$dir/$mfile_syn";
  #$file2 = "$dir/summed_ker.dat";
  $file2 = "$dir/fun_smooth.dat";
  $file3 = "$dir/summed_chi_r.dat";
  $recfile = "$dir/${edir}/sr.txt";  # src-rec for first event
  $chifile = "$dir/summed_chi_all.dat";

  if (not -f $file1syn){ die("Check if $file1syn exist or not\n") }
  if (not -f $file2)   { die("Check if $file2 exist or not\n") }
  if (not -f $file3)   { die("Check if $file3 exist or not\n") }
  if (not -f $recfile) { die("Check if $recfile exist or not\n") }
  if (not -f $chifile) { die("Check if $chifile exist or not\n") }
#=================================================================

  #=============================================
  # model for synthetics

  $title = "(e)  Test model $smod";
  $shift = $shift1;

  # phase speed map
  print CSH "psbasemap $B $R $J -P -K -O -V $shift >> $psfile\n";
  if($icolor==1) {print CSH "awk '{print \$1,\$2,\$4*100}' $file1syn | pscontour $R $J -A- -C$cpt_vel -I -P -O -K -V >> $psfile\n"}
  print CSH "pscoast $J $R $coast_info -P -K -O -V >> $psfile\n";
  if($ishelf) {print CSH "awk '{print \$2,\$1}' $shelf_file |psxy $J $R $Wshelf -K -O -P -V >> $psfile\n"}
  print CSH "awk '{print \$1,\$2}' $evefile |psxy -N $J $R -K -O -P -V $src >> $psfile\n";
  print CSH "awk '\$1 == \"R\" {print \$2,\$3}' $recfile |psxy -N $J $R -K -O -P -V $rec >> $psfile\n";
  #print CSH "psscale -C$cpt_vel $Dscale $Bscale1 -P -K -O -V >> $psfile \n";
  print CSH "pstext -N $J_title $R_title -K -O -V -P >>$psfile<<EOF\n $x_title $z_title $fsize_title 0 $fontno CM $title \nEOF\n";

  #=============================================
  # kernel -- SMOOTHED VERSION

  $shift = $shift1;
  $title = "(f)  Gradient (smoothed) for $smod";
  $B = "$B0".$Bopts[15];

  print CSH "psbasemap $B $R $J -P -K -O -V $shift >> $psfile\n";
  #if($icolor==1) {print CSH "awk '{print \$1,\$2,\$3 / $norm2 }' $file2 | pscontour $J $R -A- -C$cpt_ker -I -K -O -P -V >> $psfile\n"}
  if($icolor==1) {print CSH "awk '{print \$1,\$2,\$7 / $norm2 }' $file2 | pscontour $J $R -A- -C$cpt_ker -I -K -O -P -V >> $psfile\n"}
  print CSH "pscoast $J $R $coast_info -K -O -P -V >> $psfile\n";
  if($ishelf) {print CSH "awk '{print \$2,\$1}' $shelf_file |psxy $J $R $Wshelf -K -O -P -V >> $psfile\n"}
  #print CSH "psscale -C$cpt_ker $Dscale $Bscale2 -K -O -P -V >> $psfile \n";
  print CSH "awk '{print \$1,\$2}' $evefile |psxy -N $J $R -K -O -P -V $src >> $psfile\n";
  print CSH "pstext -N $J_title $R_title -K -O -V -P >>$psfile<<EOF\n $x_title $z_title $fsize_title 0 $fontno CM $title \nEOF\n";
  print CSH "awk '\$1 == \"R\" {print \$2,\$3}' $recfile |psxy -N $J $R -K -O -P -V $rec >> $psfile\n";

  #=============================================
  # graph

  $shift = $shift2;
  $Bg = "$B1".$Bopts[8];
  $title = "(g)  Cubic interpolation to get m\@+1\@+";

  print CSH "psbasemap $Bg $R_poly $J_poly -P -K -O -V $shift >> $psfile\n";
  print CSH "awk '{print \$1 / $normgx,\$2 / $normgy}' $poly_curve | psxy $c_info_rd $J_poly $R_poly -K -O -P -V >> $psfile\n";
  print CSH "awk '{print \$1 / $normgx,\$3 / $normgy}' $poly_curve | psxy $c_info_rd $J_poly $R_poly -K -O -P -V >> $psfile\n";
  print CSH "awk '{print \$1 / $normgx,\$5 / $normgy}' $poly_curve | psxy $c_info_ks2 $J_poly $R_poly -K -O -P -V >> $psfile\n";
  print CSH "psxy $J_poly $R_poly $c_info_ks -K -O -V -P >>$psfile<<EOF\n$x1 $a1y1 \n$x1 $y1 \n$a1x1 $y1 \nEOF\n";
  print CSH "psxy $J_poly $R_poly $c_info_ks -K -O -V -P >>$psfile<<EOF\n$x3 $a1y1 \n$x3 $y3 \n$a1x1 $y3 \nEOF\n";
  print CSH "psxy $J_poly $R_poly $c_info_kd -K -O -V -P >>$psfile<<EOF\n$x5 $a1y1 \n$x5 $y5 \n$a1x1 $y5 \nEOF\n";
  print CSH "psxy $J_poly $R_poly $p_info_k -K -O -V -P >>$psfile<<EOF\n$x1 $y1 \n$x3 $y3\nEOF\n";
  print CSH "psxy $J_poly $R_poly $p_info_w -K -O -V -P >>$psfile<<EOF\n$x2 $a1y1 \n$x4 $y4 \n$x5 $y5 \nEOF\n";
  print CSH "pstext -N $J_title $R_title -K -O -V -P >>$psfile<<EOF\n $x_title $z_title $fsize_title 0 $fontno CM $title \nEOF\n";

#=================================================================
# NEW MODEL

  # wave2d run number
  $strun = sprintf("%4.4i",$irun + 2);
  $dir = "$odir$strun";

  $mod = $mods[2];
  $smod = "m\@+$mod\@+";

  # get the summed kernel and chi files
  $file1syn = "$dir/$mfile_syn";
  $file2 = "$dir/summed_ker.dat";
  $file3 = "$dir/summed_chi_r.dat";
  $recfile = "$dir/${edir}/sr.txt";  # src-rec for first event
  $chifile = "$dir/summed_chi_all.dat";

  if (not -f $file1syn){ die("Check if $file1syn exist or not\n") }
  if (not -f $file2)   { die("Check if $file2 exist or not\n") }
  if (not -f $file3)   { die("Check if $file3 exist or not\n") }
  if (not -f $recfile) { die("Check if $recfile exist or not\n") }
  if (not -f $chifile) { die("Check if $chifile exist or not\n") }
#=================================================================

  #=============================================
  # model for synthetics

  $title = "(h)  Updated model $smod";
  $shift = $shift1;

  # phase speed map
  print CSH "psbasemap $B $R $J -P -K -O -V $shift >> $psfile\n";
  if($icolor==1) {print CSH "awk '{print \$1,\$2,\$4*100}' $file1syn | pscontour $R $J -A- -C$cpt_vel -I -P -O -K -V >> $psfile\n"}
  print CSH "pscoast $J $R $coast_info -P -K -O -V >> $psfile\n";
  if($ishelf) {print CSH "awk '{print \$2,\$1}' $shelf_file |psxy $J $R $Wshelf -K -O -P -V >> $psfile\n"}
  print CSH "awk '{print \$1,\$2}' $evefile |psxy -N $J $R -K -O -P -V $src >> $psfile\n";
  print CSH "awk '\$1 == \"R\" {print \$2,\$3}' $recfile |psxy -N $J $R -K -O -P -V $rec >> $psfile\n";
  #print CSH "psscale -C$cpt_vel $Dscale $Bscale1 -P -K -O -V >> $psfile \n";
  print CSH "pstext -N $J_title $R_title -K -O -V -P >>$psfile<<EOF\n $x_title $z_title $fsize_title 0 $fontno CM $title \nEOF\n";

  #=============================================
  # misfit vs iteration

  $shift = $shift1;
  $Bg = "$B3b".$Bopts[8];
  $title = "(i)  Misfit for first two models";

  print CSH "psbasemap $Bg $R_chi $J_chi -P -K -O -V $shift >> $psfile\n";
  print CSH "awk '{print \$1,\$2}' $chi_curve | psxy $c_info_rd $J_chi $R_chi -K -O -P -V >> $psfile\n";
  print CSH "psxy $J_chi $R_chi $p_info_k -K -O -V -P >>$psfile<<EOF\n0 $y1t \n1 $y7t \nEOF\n";
  print CSH "pstext -N $J_chi $R_chi -K -O -V -P >>$psfile<<EOF\n 0.75 $y1t $fsize1 0 $fontno LM $schi0 \nEOF\n";
  print CSH "pstext -N $J_chi $R_chi -K -O -V -P >>$psfile<<EOF\n 1.75 $y7t $fsize1 0 $fontno LM $schi1 \nEOF\n";
  print CSH "pstext -N $J_title $R_title -K -O -V -P >>$psfile<<EOF\n $x_title $z_title $fsize_title 0 $fontno CM $title \nEOF\n";


#-----------------------------

  print CSH "pstext -N $J_title $R_title -O -V -P >>$psfile<<EOF\n 10 10 $fsize_title 0 $fontno CM $title \nEOF\n"; # FINISH
  print CSH "convert $psfile $jpgfile\n";
  print CSH "echo done with $psfile\n";

  close (CSH);
  system("csh -f $cshfile");
  system("xv $jpgfile &");
}


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

  $niter = 8;

  for ($k = 0; $k <= $niter; $k = $k+1) {
    $irun = $irun0 + 2*$k;              # wave2d run number
    $strun = sprintf("%4.4i",$irun);
    $dir = "$odir$strun";
    $file1syn = "$dir/$mfile_syn";
    $file1ker = "$dir/summed_ker.dat";
    $ker_files[$k] = $file1ker;
    $mod_files[$k] = $file1syn;
    if (not -f $file1syn)   { die("Check if $file1syn exist or not\n") }

    # load chi values
    $chifile = "$dir/summed_chi_all.dat";
    open(IN,"$chifile");
    $chi = <IN>;
    $it_vals[$k] = $k;
    $chi_vals[$k] = $chi;
    #$schi = sprintf("\@~\143\@~ ( $smod )  =  %3.3e",$chi);
  }

  #die("testing");

  # get the receivers and sources
  $recfile = "$dir0/${edir}/sr.txt";  # src-rec for first event
  $evefile = "$dir0/events_lonlat.dat";
  if (not -f $recfile) { die("Check if $recfile exist or not\n") }
  if (not -f $evefile) { die("Check if $evefile exist or not\n") }

  # file names for figures
  $name    = "${smap}_cg03_${strun0}";
  $psfile  = "$name.ps";
  $jpgfile = "$name.jpg";

  #=============================================
  # model for the data

  $B = "$B0".$Bopts[15];
  $title = "Target phase speed model (data)";

  # phase speed map
  print CSH "psbasemap $B $R $J -P -K -V $origin > $psfile\n"; # START
  if($icolor==1) {print CSH "awk '{print \$1,\$2,\$4*100}' $file1dat | pscontour $R $J -A- -C$cpt_vel -I -P -K -O -V  >> $psfile\n"}
  print CSH "pscoast $J $R $coast_info -P -K -O -V >> $psfile\n";
  if($ishelf) {print CSH "awk '{print \$2,\$1}' $shelf_file |psxy $J $R $Wshelf -K -O -P -V >> $psfile\n"}
  print CSH "awk '\$1 == \"R\" {print \$2,\$3}' $recfile |psxy -N $J $R -K -O -P -V $rec >> $psfile\n";
  print CSH "awk '{print \$1,\$2}' $evefile |psxy -N $J $R -K -O -P -V $src >> $psfile\n";
  #print CSH "psscale -C$cpt_vel $Dscale $Bscale1 -P -K -O -V >> $psfile \n";
  print CSH "pstext -N $J_title $R_title -K -O -V -P >>$psfile<<EOF\n $x_title $z_title $fsize_title 0 $fontno CM $title \nEOF\n";

  #=============================================
  # model for synthetics

  $niter_max = 4;
  for ($k = 0; $k <= $niter_max; $k = $k+1) {

    $mod = $mods[$k*2]; $smod = "m\@+$mod\@+"; $title = "Phase speed model $smod";
    if ($k % 3 == 2) {$shift = $shift2} else {$shift = $shift1}

    # phase speed map
    print CSH "psbasemap $B $R $J -P -K -O -V $shift >> $psfile\n";
    if ($icolor==1) {
      print CSH "awk '{print \$1,\$2,\$4*100}' $mod_files[$k] | pscontour $R $J -A- -C$cpt_vel -I -P -O -K -V >> $psfile\n";
    }
    print CSH "pscoast $J $R $coast_info -P -K -O -V >> $psfile\n";
    if($ishelf) {print CSH "awk '{print \$2,\$1}' $shelf_file |psxy $J $R $Wshelf -K -O -P -V >> $psfile\n"}
    print CSH "awk '{print \$1,\$2}' $evefile |psxy -N $J $R -K -O -P -V $src >> $psfile\n";
    print CSH "awk '\$1 == \"R\" {print \$2,\$3}' $recfile |psxy -N $J $R -K -O -P -V $rec >> $psfile\n";
    #print CSH "psscale -C$cpt_vel $Dscale $Bscale1 -P -K -O -V >> $psfile \n";
    print CSH "pstext -N $J_title $R_title -K -O -V -P >>$psfile<<EOF\n $x_title $z_title $fsize_title 0 $fontno CM $title \nEOF\n";
  }

  # plot final model
  $shift = $shift2;
  $mod = $mods[$niter*2]; $smod = "m\@+$mod\@+"; $title = "Phase speed model $smod";
  print CSH "psbasemap $B $R $J -P -K -O -V $shift >> $psfile\n";
  if ($icolor==1) {
    print CSH "awk '{print \$1,\$2,\$4*100}' $mod_files[$niter] | pscontour $R $J -A- -C$cpt_vel -I -P -O -K -V >> $psfile\n";
  }
  print CSH "pscoast $J $R $coast_info -P -K -O -V >> $psfile\n";
  if($ishelf) {print CSH "awk '{print \$2,\$1}' $shelf_file |psxy $J $R $Wshelf -K -O -P -V >> $psfile\n"}
  print CSH "awk '{print \$1,\$2}' $evefile |psxy -N $J $R -K -O -P -V $src >> $psfile\n";
  print CSH "awk '\$1 == \"R\" {print \$2,\$3}' $recfile |psxy -N $J $R -K -O -P -V $rec >> $psfile\n";
  #print CSH "psscale -C$cpt_vel $Dscale $Bscale1 -P -K -O -V >> $psfile \n";
  print CSH "pstext -N $J_title $R_title -K -O -V -P >>$psfile<<EOF\n $x_title $z_title $fsize_title 0 $fontno CM $title \nEOF\n";

  #=============================================
  # misfit vs iteration

  $shift = $shift1;
  $Bg = "$B3b".$Bopts[8];
  $title = "Misfit,  \@~\143\@~ (m)  ( $utype )";

  print CSH "psbasemap $Bg $R_chi $J_chi -P -K -O -V $shift >> $psfile\n";
  print CSH "awk '{print \$1,\$2}' $chi_curve | psxy $c_info_rd $J_chi $R_chi -K -O -P -V >> $psfile\n";
  for ($k = 0; $k <= $niter; $k = $k+1) {
    print CSH "psxy $J_chi $R_chi $p_info_k -K -O -V -P >>$psfile<<EOF\n $it_vals[$k] $chi_vals[$k]\nEOF\n";
  }
  print CSH "pstext -N $J_chi $R_chi -K -O -V -P >>$psfile<<EOF\n 1 $y1t $fsize1 0 $fontno LM $schi0 \nEOF\n";
  print CSH "pstext -N $J_chi $R_chi -K -O -V -P >>$psfile<<EOF\n 2 $y7t $fsize1 0 $fontno LM $schi1 \nEOF\n";
  print CSH "pstext -N $J_title $R_title -K -O -V -P >>$psfile<<EOF\n $x_title $z_title $fsize_title 0 $fontno CM $title \nEOF\n";

  #=============================================
  # var-red vs iteration

  $shift = $shift1;
  $Bg = "$B4".$Bopts[8];
  $title = "Variance reduction";

  print CSH "psbasemap $Bg $R_var $J_var -P -K -O -V $shift >> $psfile\n";

  print CSH "awk '{print \$1,\$2}' $var_fit | psxy $c_info_rd $J_var $R_var -K -O -P -V >> $psfile\n";
  print CSH "awk '{print \$1,\$2}' $var_pts | psxy $p_info_k $J_var $R_var -K -O -P -V >> $psfile\n";

  print CSH "pstext -N $J_title $R_title -K -O -V -P >>$psfile<<EOF\n $x_title $z_title $fsize_title 0 $fontno CM $title \nEOF\n";

#-----------------------------
  print CSH "pstext -N $J_title $R_title -O -V -P >>$psfile<<EOF\n 10 10 $fsize_title 0 $fontno CM $title \nEOF\n"; # FINISH
  print CSH "convert $psfile $jpgfile\n";
  print CSH "echo done with $psfile\n";

  close (CSH);
  system("csh -f $cshfile");
  system("xv $jpgfile &");
}


#===========================================================================
if ($ifig04 == 1) {

  # shifting the subplots
  $xfac = 1.20;
  $yfac = 1.20;
  $dX = $xfac*$Jwid; $dY = 0; $shift1 = "-X$dX -Y$dY";
  $dX2 = -2*$dX; $dY2 = -$yfac*$Jwid; $shift2 = "-X$dX2 -Y$dY2";

  # colorbar
  $Dlen = 1.2; $Dx = $Jwid/2; $Dy = -0.10*$Jwid;
  $Dscale = "-D$Dx/$Dy/$Dlen/0.10h";

  #===============================================
  print "\nWriting CSH file...\n";

  open(CSH,">$cshfile");
  print CSH "gmtset BASEMAP_TYPE plain TICK_LENGTH $tick LABEL_FONT_SIZE $fsize2 ANOT_FONT_SIZE $fsize2 PLOT_DEGREE_FORMAT D HEADER_FONT $fontno ANOT_FONT $fontno LABEL_FONT $fontno HEADER_FONT_SIZE $fsize1 FRAME_PEN $fpen TICK_PEN $tpen\n";
  #===============================================

  # iteration index (subset of all possible)
  #@kvec = (0,1,2,3);
  #@kvec = (0,1,2,3,8,12,16);
  #@kvec = (0,1,2,3,8);
  @kvec = (1,2,3,4,8);
  $numk = @kvec;
  $niter_max = $kvec[$numk-1];  # max iteration is the last in the list

  if($numk % 3 == 1) {$shift3 = "-X$dX"}
  if($numk % 3 == 2) {$shift3 = "-Y$dY2"}
  if($numk % 3 == 0) {$shift3 = "-X-$dX -Y$dY2"}

  # load all possible files
  for ($k = 0; $k <= $niter_max; $k = $k+1) {

    $irun = $irun0 + 2*$k;              # wave2d run number

    $strun = sprintf("%4.4i",$irun);
    $dir = "$odir$strun";
    $file1syn = "$dir/$mfile_syn";
    $file1ker = "$dir/summed_ker.dat";
    if (not -f $file1syn)   { die("Check if $file1syn exist or not\n") }
    #if (not -f $file1ker)   { die("Check if $file1ker exist or not\n") }
    $ker_files[$k] = $file1ker;
    $mod_files[$k] = $file1syn;

    # load chi values
    $chifile = "$dir/summed_chi_all.dat";
    open(IN,"$chifile");
    $chi = <IN>;
    chomp($chi);
    $it_vals[$k] = $k;
    $chi_vals[$k] = $chi;
    #$schi = sprintf("\@~\143\@~ ( $smod )  =  %3.3e",$chi);
  }
  #print "\n @it_vals \n @chi_vals \n\n"; die("testing");

  # get the receivers and sources
  $recfile = "$dir0/${edir}/sr.txt";  # src-rec for first event
  $evefile = "$dir0/events_lonlat.dat";
  if (not -f $recfile) { die("Check if $recfile exist or not\n") }
  if (not -f $evefile) { die("Check if $evefile exist or not\n") }

  # file names for figures
  $name    = "${smap}_cg04_${strun0}";
  $psfile  = "$name.ps";
  $jpgfile = "$name.jpg";

  $B = "$B0".$Bopts[15];

  #=============================================
  # phase speed map for m1

  #$k = 1;
  #$mod = $mods[$k*2]; $smod = "m\@+$mod\@+"; $title = "Phase speed model $smod";

  #print CSH "psbasemap $B $R $J -P -K -V $origin > $psfile\n"; # START
  #if($icolor==1) {print CSH "awk '{print \$1,\$2,\$4*100}' $mod_files[$k] | pscontour $R $J -A- -C$cpt_vel -I -P -K -O -V  >> $psfile\n"}
  #print CSH "pscoast $J $R $coast_info -P -K -O -V >> $psfile\n";
  #if($ishelf) {print CSH "awk '{print \$2,\$1}' $shelf_file |psxy $J $R $Wshelf -K -O -P -V >> $psfile\n"}
  #print CSH "awk '\$1 == \"R\" {print \$2,\$3}' $recfile |psxy -N $J $R -K -O -P -V $rec >> $psfile\n";
  #print CSH "awk '{print \$1,\$2}' $evefile |psxy -N $J $R -K -O -P -V $src >> $psfile\n";
  ##print CSH "psscale -C$cpt_vel $Dscale $Bscale1 -P -K -O -V >> $psfile \n";
  #print CSH "pstext -N $J_title $R_title -K -O -V -P >>$psfile<<EOF\n $x_title $z_title $fsize_title 0 $fontno CM $title \nEOF\n";

  #=============================================
  # recovered models for synthetics

  @labs  = ("a","b","c","d","e","f","g","h","i","j");

  #$numk = @kvec;
  for ($i = 0; $i < $numk; $i = $i+1) {
    $k = $kvec[$i];

    $mod = $mods[2*$k]; $smod = "m\@+${mod}\@+"; $title = "($labs[$i])  Phase speed model ${smod}";
    if ($i % 3 == 0) {$shift = $shift2} else {$shift = $shift1}

    # phase speed map
    if($i == 0) { print CSH "psbasemap $B $R $J -P -K -V $origin > $psfile\n"} # START
    else        { print CSH "psbasemap $B $R $J -P -K -O -V $shift >> $psfile\n"}
    if ($icolor==1) {
      print CSH "awk '{print \$1,\$2,\$4*100}' $mod_files[$k] | pscontour $R $J -A- -C$cpt_vel -I -P -O -K -V >> $psfile\n";
    }
    print CSH "pscoast $J $R $coast_info -P -K -O -V >> $psfile\n";
    if($ishelf) {print CSH "awk '{print \$2,\$1}' $shelf_file |psxy $J $R $Wshelf -K -O -P -V >> $psfile\n"}
    print CSH "awk '{print \$1,\$2}' $evefile |psxy -N $J $R -K -O -P -V $src >> $psfile\n";
    print CSH "awk '\$1 == \"R\" {print \$2,\$3}' $recfile |psxy -N $J $R -K -O -P -V $rec >> $psfile\n";
    #print CSH "psscale -C$cpt_vel $Dscale $Bscale1 -P -K -O -V >> $psfile \n";
    print CSH "pstext -N $J_title $R_title -K -O -V -P >>$psfile<<EOF\n $x_title $z_title $fsize_title 0 $fontno CM $title \nEOF\n";
  }

  #=============================================
  # misfit vs iteration

  # whether to plot the phase speed map for data
  $idatamap = 0;    # =0 for GJI figure
  if($idatamap==1) {$shift = $shift3} else {$shift = $shift1}

  $Bg = "$B3b".$Bopts[8];
  $title = "($labs[$i])  Misfit,  \@~\143\@~ (m\@+k\@+)  ($utype)";

  print CSH "psbasemap $Bg $R_chi $J_chi -P -K -O -V $shift >> $psfile\n";
  print CSH "awk '{print \$1,\$2}' $chi_curve | psxy $c_info_rd $J_chi $R_chi -K -O -P -V >> $psfile\n";
  for ($k = 0; $k <= $niter_max; $k = $k+1) {
    print CSH "psxy $J_chi $R_chi $p_info_k -K -O -V -P >>$psfile<<EOF\n $it_vals[$k] $chi_vals[$k]\nEOF\n";  # plot point
  }
  # labels for chi0 and chi1
  print CSH "pstext -N $J_chi $R_chi -K -O -V -P >>$psfile<<EOF\n $x0chit $y1t $fsize1 0 $fontno LM $schi0 \nEOF\n";
  print CSH "pstext -N $J_chi $R_chi -K -O -V -P >>$psfile<<EOF\n $x1chit $y7t $fsize1 0 $fontno LM $schi1 \nEOF\n";
  print CSH "pstext -N $J_title $R_title -K -O -V -P >>$psfile<<EOF\n $x_title $z_title $fsize_title 0 $fontno CM $title \nEOF\n";

  #=============================================
  # phase speed map for data

  if($idatamap==1) {

  $shift = $shift1;
  $title = "Target phase speed model (data)";
  print CSH "psbasemap $B $R $J -P -K -O -V $shift >> $psfile\n";
  if($icolor==1) {print CSH "awk '{print \$1,\$2,\$4*100}' $file1dat | pscontour $R $J -A- -C$cpt_vel -I -P -K -O -V  >> $psfile\n"}
  print CSH "pscoast $J $R $coast_info -P -K -O -V >> $psfile\n";
  if($ishelf) {print CSH "awk '{print \$2,\$1}' $shelf_file |psxy $J $R $Wshelf -K -O -P -V >> $psfile\n"}
  print CSH "awk '\$1 == \"R\" {print \$2,\$3}' $recfile |psxy -N $J $R -K -O -P -V $rec >> $psfile\n";
  print CSH "awk '{print \$1,\$2}' $evefile |psxy -N $J $R -K -O -P -V $src >> $psfile\n";
  print CSH "psscale -C$cpt_vel $Dscale $Bscale1 -P -K -O -V >> $psfile \n";
  print CSH "pstext -N $J_title $R_title -K -O -V -P >>$psfile<<EOF\n $x_title $z_title $fsize_title 0 $fontno CM $title \nEOF\n";

  }

#-----------------------------
  print CSH "pstext -N $J_title $R_title -O -V -P >>$psfile<<EOF\n 10 10 $fsize_title 0 $fontno CM $title \nEOF\n"; # FINISH
  print CSH "convert $psfile $jpgfile\n";
  print CSH "echo done with $psfile\n";

  close (CSH);
  system("csh -f $cshfile");
  system("xv $jpgfile &");
}


#===========================================================================
if ($ifig04b == 1) {

  # shifting the subplots
  $xfac = 1.15;
  $yfac = 1.15;
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
    $dir = "$odir$strun";
    $file1syn = "$dir/$mfile_syn";
    $file1ker = "$dir/summed_ker.dat";
    $ker_files[$k] = $file1ker;
    $mod_files[$k] = $file1syn;
    if (not -f $file1syn)   { die("Check if $file1syn exist or not\n") }
  }

  #die("testing");

  # get the receivers and sources
  $recfile = "$dir0/${edir}/sr.txt";  # src-rec for first event
  $evefile = "$dir0/events_lonlat.dat";
  if (not -f $recfile) { die("Check if $recfile exist or not\n") }
  if (not -f $evefile) { die("Check if $evefile exist or not\n") }

  # file names for figures
  $name    = "${smap}_ker9_${stGam2}_${strun0}";
  $psfile  = "$name.ps";
  $jpgfile = "$name.jpg";

  #=============================================
  # kernels

  $B = "$B0".$Bopts[15];
  $niter_max = $niter;

  for ($k = 0; $k <= $niter_max; $k = $k+1) {

    $mod = $mods[$k*2]; $smod = "m\@+$mod\@+"; $title = "Kernel for $smod";
    if ($k % 3 == 0) {$shift = $shift2} else {$shift = $shift1}

    # phase speed map
    if ($k == 0) {print CSH "psbasemap $B $R $J -P -K -V $origin > $psfile\n"}
    else         {print CSH "psbasemap $B $R $J -P -K -O -V $shift >> $psfile\n"}

    if($icolor==1) {print CSH "awk '{print \$1,\$2,\$3 / $norm2 }' $ker_files[$k] | pscontour $J $R -A- -C$cpt_ker -I -K -O -P -V >> $psfile\n"}
    print CSH "pscoast $J $R $coast_info -P -K -O -V >> $psfile\n";
    if($ishelf) {print CSH "awk '{print \$2,\$1}' $shelf_file |psxy $J $R $Wshelf -K -O -P -V >> $psfile\n"}
    print CSH "awk '{print \$1,\$2}' $evefile |psxy -N $J $R -K -O -P -V $src >> $psfile\n";
    print CSH "awk '\$1 == \"R\" {print \$2,\$3}' $recfile |psxy -N $J $R -K -O -P -V $rec >> $psfile\n";
    #print CSH "psscale -C$cpt_vel $Dscale $Bscale1 -P -K -O -V >> $psfile \n";
    print CSH "pstext -N $J_title $R_title -K -O -V -P >>$psfile<<EOF\n $x_title $z_title $fsize_title 0 $fontno CM $title \nEOF\n";
  }

#-----------------------------
  print CSH "pstext -N $J_title $R_title -O -V -P >>$psfile<<EOF\n 10 10 $fsize_title 0 $fontno CM $title \nEOF\n"; # FINISH
  print CSH "convert $psfile $jpgfile\n";
  print CSH "echo done with $psfile\n";

  close (CSH);
  system("csh -f $cshfile");
  system("xv $jpgfile &");
}


#===========================================================================
if ($ifig05 == 1) {

  # shifting the subplots
  $xfac = 1.15;
  $yfac = 1.20;
  $dX = $xfac*$Jwid; $dY = 0; $shift1 = "-X$dX -Y$dY";
  $dX = -$dX; $dY = -$yfac*$Jwid; $shift2 = "-X$dX -Y$dY";

  # colorbar
  $Dlen = 1.2; $Dx = $Jwid/2; $Dy = -0.45;
  $Dscale = "-D$Dx/$Dy/$Dlen/0.10h";

  #===============================================
  print "\nWriting CSH file...\n";

  open(CSH,">$cshfile");
  print CSH "gmtset BASEMAP_TYPE plain TICK_LENGTH $tick LABEL_FONT_SIZE $fsize2 ANOT_FONT_SIZE $fsize2 PLOT_DEGREE_FORMAT D HEADER_FONT $fontno ANOT_FONT $fontno LABEL_FONT $fontno HEADER_FONT_SIZE $fsize1 FRAME_PEN $fpen TICK_PEN $tpen\n";
  #===============================================

  # file with desired functions
  $file_smooth = "$dir0/fun_smooth.dat";
  if (not -f $file_smooth) { die("Check if $file_smooth exist or not\n") }

  # get the receivers and sources
  $recfile = "$dir0/${edir}/sr.txt";  # src-rec for first event
  $evefile = "$dir0/events_lonlat.dat";
  if (not -f $recfile) { die("Check if $recfile exist or not\n") }
  if (not -f $evefile) { die("Check if $evefile exist or not\n") }

  # file names for figures
  $name    = "${smap}_smooth_${strun0}";
  $psfile  = "$name.ps";
  $jpgfile = "$name.jpg";

  #=============================================
  # model for the data

  $B = "$B0".$Bopts[8]; $title = "(a)  Unsmoothed kernel for m\@+0\@+"; $shift = $shift1;

  # phase speed map
  print CSH "psbasemap $B $R $J -P -K -V $origin > $psfile\n"; # START
  if($icolor==1) {print CSH "awk '{print \$1,\$2,\$5 / $norm2 }' $file_smooth | pscontour $R $J -A- -C$cpt_ker -I -P -K -O -V  >> $psfile\n"}
  print CSH "pscoast $J $R $coast_info2 -P -K -O -V >> $psfile\n";
  #if($ishelf) {print CSH "awk '{print \$2,\$1}' $shelf_file |psxy $J $R $Wshelf -K -O -P -V >> $psfile\n"}
  print CSH "awk '\$1 == \"R\" {print \$2,\$3}' $recfile |psxy -N $J $R -K -O -P -V $rec >> $psfile\n";
  print CSH "awk '{print \$1,\$2}' $evefile |psxy -N $J $R -K -O -P -V $src >> $psfile\n";
  print CSH "psscale -C$cpt_ker $Dscale $Bscale2 -P -K -O -V >> $psfile \n";
  print CSH "pstext -N $J_title $R_title -K -O -V -P >>$psfile<<EOF\n $x_title $z_title $fsize_title 0 $fontno CM $title \nEOF\n";

  #$stsig = sprintf("\@~\163\@~ = %.1f km",$sigma/1000);
  #$stGam = sprintf("\@~\107\@~ = %.1f km",$Gamma/1000);
  $B = "$B0".$Bopts[15]; $title = "(b)  Smoothed kernel,  $stGam"; $shift = $shift1;

  # smooth function
  print CSH "psbasemap $B $J $R -P -K -O -V $shift >> $psfile\n";
  if($icolor==1) {print CSH "awk '{print \$1,\$2,\$7 / $norm2}' $file_smooth | pscontour $R $J -A- -C$cpt_ker -I -P -K -O -V  >> $psfile\n"}
  print CSH "pscoast $J $R $coast_info2 -P -K -O -V >> $psfile\n";
  #if($ishelf) {print CSH "awk '{print \$2,\$1}' $shelf_file |psxy $J $R $Wshelf -K -O -P -V >> $psfile\n"}
  print CSH "awk '\$1 == \"R\" {print \$2,\$3}' $recfile |psxy -N $J $R -K -O -P -V $rec >> $psfile\n";
  print CSH "awk '{print \$1,\$2}' $evefile |psxy -N $J $R -K -O -P -V $src >> $psfile\n";
  print CSH "pstext -N $J_title $R_title -K -O -V -P >>$psfile<<EOF\n $x_title $z_title $fsize_title 0 $fontno CM $title \nEOF\n";

  #------------------------------------------
  # make colorpoint file for Gaussian
  #$cmax2 = 1/(2*3.14159*$sigma*$sigma);
  $cmax2 = 4/(3.14159*$Gamma*$Gamma);
  $dc2 = $cmax2/10;
  $T2  = "-T-${cmax2}/${cmax2}/$dc2";
  $cpt_gaus = "color_g1.cpt";
  print CSH "makecpt -C$colorbar $T2 > temp1\n";
  print CSH "sed 's/^B.*/B       170     0      0  /' temp1  >  temp2\n";
  print CSH "sed 's/^F.*/F         0     0    255  /' temp2 > $cpt_gaus\n";

  # gaussian function bounds
  # BE CAREFUL THAT THE GAUSSIAN IS SCALED PROPERLY
  $xmin0 = 0; $xmax0 = 480; $zmin0 = 0; $zmax0 = 480;
  $xceng = ($xmax0-$xmin0)/2;
  $zceng = ($zmax0-$zmin0)/2;
  $wfac = 130;   # width of inset box
  $xmin = $xceng - $wfac/2; $xmax = $xceng + $wfac/2;
  $zmin = $zceng - $wfac/2; $zmax = $zceng + $wfac/2;
  $Rgaus = "-R$xmin/$xmax/$zmin/$zmax";
  $Jwid_gaus = $wfac/($xmax0-$xmin0) * $Jwid;
  $Jgaus = "-JX${Jwid_gaus}i";
  #print "\n $Jgaus $Rgaus \n $J $R \n"; die("testing");

  print CSH "gmtset TICK_LENGTH 0\n";
  print CSH "awk '{print \$3/1000,\$4/1000,\$6}' $file_smooth | pscontour $Rgaus $Jgaus -A- -C$cpt_gaus -I -P -K -O -V  >> $psfile\n";
  #print CSH "awk '{print \$1,\$2,\$6}' $file_smooth | pscontour $Rgaus $Jgaus -A- -C$cpt_gaus -I -P -K -O -V  >> $psfile\n";

  # plot a scale bar that is sfac times sigma in length
  $zsig1 = $zceng - $Gamma/1000 / 2;
  $zsig2 = $zceng + $Gamma/1000 / 2;
  $xsig  = $xceng - $wfac/2 * 0.7;
  print CSH "psxy $Jgaus $Rgaus -W1.0p -K -O -V -P >>$psfile<<EOF\n$xsig $zsig1 \n$xsig $zsig2 \nEOF\n";

  print CSH "psbasemap -B $Jgaus $Rgaus -P -K -O -V >> $psfile\n";
  print CSH "gmtset TICK_LENGTH $tick\n";
  #------------------------------------------

  $B = "$B0".$Bopts[15]; $title = "(c)  Residual  =  (a) - (b)"; $shift = $shift1;

  # residual function
  print CSH "psbasemap $B $J $R -P -K -O -V $shift >> $psfile\n";
  if($icolor==1) {print CSH "awk '{print \$1,\$2,\$8 / $norm2}' $file_smooth | pscontour $R $J -A- -C$cpt_ker -I -P -K -O -V  >> $psfile\n"}
  print CSH "pscoast $J $R $coast_info2 -P -K -O -V >> $psfile\n";
  #if($ishelf) {print CSH "awk '{print \$2,\$1}' $shelf_file |psxy $J $R $Wshelf -K -O -P -V >> $psfile\n"}
  #print CSH "awk '\$1 == \"R\" {print \$2,\$3}' $recfile |psxy -N $J $R -K -O -P -V $rec >> $psfile\n";
  #print CSH "awk '{print \$1,\$2}' $evefile |psxy -N $J $R -K -O -P -V $src >> $psfile\n";
  print CSH "pstext -N $J_title $R_title -K -O -V -P >>$psfile<<EOF\n $x_title $z_title $fsize_title 0 $fontno CM $title \nEOF\n";

  #-----------------------------

  # wave2d run number
  $irun = 400; $strun = sprintf("%4.4i",$irun); $dir = "$odir$strun";
  $file_smooth = "$dir/fun_smooth.dat";
  if (not -f $file_smooth) { die("Check if $file_smooth exist or not\n") }

  $Gamma = 15000;
  #$stsig = sprintf("\@~\163\@~ = %.1f km",$sigma/1000);
  $stGam = sprintf("\@~\107\@~ = %.1f km",$Gamma/1000);
  $B = "$B0".$Bopts[15]; $title = "(d)  Smoothed kernel,  $stGam"; $shift = $shift2;

  # smooth function
  print CSH "psbasemap $B $J $R -P -K -O -V $shift >> $psfile\n";
  if($icolor==1) {print CSH "awk '{print \$1,\$2,\$7 / $norm2}' $file_smooth | pscontour $R $J -A- -C$cpt_ker -I -P -K -O -V  >> $psfile\n"}
  print CSH "pscoast $J $R $coast_info2 -P -K -O -V >> $psfile\n";
  #if($ishelf) {print CSH "awk '{print \$2,\$1}' $shelf_file |psxy $J $R $Wshelf -K -O -P -V >> $psfile\n"}
  print CSH "awk '\$1 == \"R\" {print \$2,\$3}' $recfile |psxy -N $J $R -K -O -P -V $rec >> $psfile\n";
  print CSH "awk '{print \$1,\$2}' $evefile |psxy -N $J $R -K -O -P -V $src >> $psfile\n";
  print CSH "pstext -N $J_title $R_title -K -O -V -P >>$psfile<<EOF\n $x_title $z_title $fsize_title 0 $fontno CM $title \nEOF\n";

  #------------------------------------------
  # make colorpoint file for Gaussian
  #$cmax2 = 1/(2*3.14159*$sigma*$sigma);
  $cmax2 = 4/(3.14159*$Gamma*$Gamma);
  $dc2 = $cmax2/10;
  $T2  = "-T-${cmax2}/${cmax2}/$dc2";
  $cpt_gaus = "color_g2.cpt";
  print CSH "makecpt -C$colorbar $T2 > temp1\n";
  print CSH "sed 's/^B.*/B       170     0      0  /' temp1  >  temp2\n";
  print CSH "sed 's/^F.*/F         0     0    255  /' temp2 > $cpt_gaus\n";

  print CSH "gmtset TICK_LENGTH 0\n";
  print CSH "awk '{print \$3/1000,\$4/1000,\$6}' $file_smooth | pscontour $Rgaus $Jgaus -A- -C$cpt_gaus -I -P -K -O -V  >> $psfile\n";
  #print CSH "awk '{print \$1,\$2,\$6}' $file_smooth | pscontour $Rgaus $Jgaus -A- -C$cpt_gaus -I -P -K -O -V  >> $psfile\n";

  # plot a scale bar that is sfac times sigma in length
  $zsig1 = $zceng - $Gamma/1000 / 2;
  $zsig2 = $zceng + $Gamma/1000 / 2;
  print CSH "psxy $Jgaus $Rgaus -W1.0p -K -O -V -P >>$psfile<<EOF\n$xsig $zsig1 \n$xsig $zsig2 \nEOF\n";

  print CSH "psbasemap -B $Jgaus $Rgaus -P -K -O -V >> $psfile\n";
  print CSH "gmtset TICK_LENGTH $tick\n";
  #------------------------------------------

  $B = "$B0".$Bopts[15]; $title = "(e)  Residual  =  (a) - (d)"; $shift = $shift1;

  # residual function
  print CSH "psbasemap $B $J $R -P -K -O -V $shift >> $psfile\n";
  if($icolor==1) {print CSH "awk '{print \$1,\$2,\$8 / $norm2}' $file_smooth | pscontour $R $J -A- -C$cpt_ker -I -P -K -O -V  >> $psfile\n"}
  print CSH "pscoast $J $R $coast_info2 -P -K -O -V >> $psfile\n";
  #if($ishelf) {print CSH "awk '{print \$2,\$1}' $shelf_file |psxy $J $R $Wshelf -K -O -P -V >> $psfile\n"}
  #print CSH "awk '\$1 == \"R\" {print \$2,\$3}' $recfile |psxy -N $J $R -K -O -P -V $rec >> $psfile\n";
  #print CSH "awk '{print \$1,\$2}' $evefile |psxy -N $J $R -K -O -P -V $src >> $psfile\n";
  print CSH "pstext -N $J_title $R_title -K -O -V -P >>$psfile<<EOF\n $x_title $z_title $fsize_title 0 $fontno CM $title \nEOF\n";

#-----------------------------
  print CSH "pstext -N $J_title $R_title -O -V -P >>$psfile<<EOF\n 10 10 $fsize_title 0 $fontno CM $title \nEOF\n"; # FINISH
  print CSH "convert $psfile $jpgfile\n";
  print CSH "echo done with $psfile\n";

  close (CSH);
  system("csh -f $cshfile");
  system("xv $jpgfile &");
}


#===========================================================================
if ($ifig06 == 1) {

  # shifting the subplots
  $xfac = 1.15;
  $yfac = 1.15;
  $dX = $xfac*$Jwid; $dY = 0; $shift1 = "-X$dX -Y$dY";
  $dX = -2*$dX; $dY = -$yfac*$Jwid; $shift2 = "-X$dX -Y$dY";

  # colorbar
  $Dlen = 1.2; $Dx = $Jwid/2; $Dy = -0.45;
  $Dscale = "-D$Dx/$Dy/$Dlen/0.10h";

  # all models are m8
  $mod = $mods[16]; $smod = "m\@+$mod\@+";

  #===============================================
  print "\nWriting CSH file...\n";

  open(CSH,">$cshfile");
  print CSH "gmtset BASEMAP_TYPE plain TICK_LENGTH $tick LABEL_FONT_SIZE $fsize2 ANOT_FONT_SIZE $fsize2 PLOT_DEGREE_FORMAT D HEADER_FONT $fontno ANOT_FONT $fontno LABEL_FONT $fontno HEADER_FONT_SIZE $fsize1 FRAME_PEN $fpen TICK_PEN $tpen\n";
  #===============================================

  # get the receivers and sources
  $recfile = "$dir0/${edir}/sr.txt";  # src-rec for first event
  $evefile = "$dir0/events_lonlat.dat";
  if (not -f $recfile) { die("Check if $recfile exist or not\n") }
  if (not -f $evefile) { die("Check if $evefile exist or not\n") }

  # file names for figures
  $name    = "smooth_vs_scale";
  $psfile  = "$name.ps";
  $jpgfile = "$name.jpg";

  #=============================================
  # model for the data

  @runs  = (420,800,780,460,820,900,500,840,920);                       # irun0
  @ifile = (1,1,1,1,1,1,1,1,1);                                         # boolean: which runs are done
  @imods = (8,8,8,8,8,8,8,8,7);                                         # which model to plot
  @gams  = (30000,30000,30000,60000,60000,60000,90000,90000,90000);     # Gamma smoothing
  @nfacs = (3,2,1,3,2,1,3,2,1);                                         # scalelength of cheker map
  @labs  = ("a","e","i","b","f","j","c","g","k","d","h","l");

  # phase speed model
  $k = 0;
  $irun0 = $runs[$k]; $strun = sprintf("%4.4i",$irun0); $dir = "$odir$strun"; $file1dat = "$dir/$mfile_dat";
  if (not -f $file1dat) { die("Check if $file1dat exist or not\n") }
  $B = "$B0".$Bopts[15]; $title = "($labs[$k]) Target phase speed model  (n = $nfacs[$k])";

  print CSH "psbasemap $B $R $J -P -K -V $origin > $psfile\n"; # START
  if($icolor==1) {print CSH "awk '{print \$1,\$2,\$4*100}' $file1dat | pscontour $R $J -A- -C$cpt_vel -I -P -K -O -V  >> $psfile\n"}
  print CSH "pscoast $J $R $coast_info -P -K -O -V >> $psfile\n";
  if($ishelf) {print CSH "awk '{print \$2,\$1}' $shelf_file |psxy $J $R $Wshelf -K -O -P -V >> $psfile\n"}
  print CSH "awk '\$1 == \"R\" {print \$2,\$3}' $recfile |psxy -N $J $R -K -O -P -V $rec >> $psfile\n";
  print CSH "awk '{print \$1,\$2}' $evefile |psxy -N $J $R -K -O -P -V $src >> $psfile\n";
  #print CSH "psscale -C$cpt_vel $Dscale $Bscale1 -P -K -O -V >> $psfile \n";
  print CSH "pstext -N $J_title $R_title -K -O -V -P >>$psfile<<EOF\n $x_title $z_title $fsize_title 0 $fontno CM $title \nEOF\n";

  # phase speed model
  $k = $k + 1;
  $irun0 = $runs[$k]; $strun = sprintf("%4.4i",$irun0); $dir = "$odir$strun"; $file1dat = "$dir/$mfile_dat";
  if (not -f $file1dat) { die("Check if $file1dat exist or not\n") }
  $B = "$B0".$Bopts[15]; $title = "($labs[$k]) Target phase speed model  (n = $nfacs[$k])";  $shift = $shift1;

  print CSH "psbasemap $B $R $J -P -K -O -V $shift >> $psfile\n";
  if($icolor==1) {print CSH "awk '{print \$1,\$2,\$4*100}' $file1dat | pscontour $R $J -A- -C$cpt_vel -I -P -K -O -V  >> $psfile\n"}
  print CSH "pscoast $J $R $coast_info -P -K -O -V >> $psfile\n";
  if($ishelf) {print CSH "awk '{print \$2,\$1}' $shelf_file |psxy $J $R $Wshelf -K -O -P -V >> $psfile\n"}
  print CSH "awk '\$1 == \"R\" {print \$2,\$3}' $recfile |psxy -N $J $R -K -O -P -V $rec >> $psfile\n";
  print CSH "awk '{print \$1,\$2}' $evefile |psxy -N $J $R -K -O -P -V $src >> $psfile\n";
  print CSH "pstext -N $J_title $R_title -K -O -V -P >>$psfile<<EOF\n $x_title $z_title $fsize_title 0 $fontno CM $title \nEOF\n";

  # phase speed model
  $k = $k + 1;
  $irun0 = $runs[$k]; $strun = sprintf("%4.4i",$irun0); $dir = "$odir$strun"; $file1dat = "$dir/$mfile_dat";
  if (not -f $file1dat) { die("Check if $file1dat exist or not\n") }
  $B = "$B0".$Bopts[15]; $title = "($labs[$k]) Target phase speed model  (n = $nfacs[$k])";  $shift = $shift1;

  print CSH "psbasemap $B $R $J -P -K -O -V $shift >> $psfile\n";
  if($icolor==1) {print CSH "awk '{print \$1,\$2,\$4*100}' $file1dat | pscontour $R $J -A- -C$cpt_vel -I -P -K -O -V  >> $psfile\n"}
  print CSH "pscoast $J $R $coast_info -P -K -O -V >> $psfile\n";
  if($ishelf) {print CSH "awk '{print \$2,\$1}' $shelf_file |psxy $J $R $Wshelf -K -O -P -V >> $psfile\n"}
  print CSH "awk '\$1 == \"R\" {print \$2,\$3}' $recfile |psxy -N $J $R -K -O -P -V $rec >> $psfile\n";
  print CSH "awk '{print \$1,\$2}' $evefile |psxy -N $J $R -K -O -P -V $src >> $psfile\n";
  print CSH "pstext -N $J_title $R_title -K -O -V -P >>$psfile<<EOF\n $x_title $z_title $fsize_title 0 $fontno CM $title \nEOF\n";

  #=========================

  $kmax = 8;
  for ($k = 0; $k <= $kmax; $k = $k+1) {

    if ($k % 3 == 0) {$shift = $shift2} else {$shift = $shift1}

    # recovered phase speed model
    $irun0 = $runs[$k]; $Gam = $gams[$k];
    #$stsig = sprintf("\@~\163\@~ = %.1f km",$sigma/1000);
    $stGam = sprintf("\@~\107\@~ = %.1f km",$Gam/1000);
    $strun = sprintf("%4.4i",$irun0+2*$imods[$k]); $dir = "$odir$strun"; $file1syn = "$dir/$mfile_syn";
    $B = "$B0".$Bopts[15]; $title = "($labs[$k+3]) Recovered  (m\@+$imods[$k]\@+,  $stGam)";

    print CSH "psbasemap $B $R $J -P -K -O -V $shift >> $psfile\n";
    if($ifile[$k]==1) {
       if (not -f $file1syn) { die("Check if $file1syn exist or not\n") }
       if($icolor==1) {print CSH "awk '{print \$1,\$2,\$4*100}' $file1syn | pscontour $R $J -A- -C$cpt_vel -I -P -K -O -V  >> $psfile\n"}
     }
    print CSH "pscoast $J $R $coast_info -P -K -O -V >> $psfile\n";
    if($ishelf) {print CSH "awk '{print \$2,\$1}' $shelf_file |psxy $J $R $Wshelf -K -O -P -V >> $psfile\n"}
    print CSH "awk '\$1 == \"R\" {print \$2,\$3}' $recfile |psxy -N $J $R -K -O -P -V $rec >> $psfile\n";
    print CSH "awk '{print \$1,\$2}' $evefile |psxy -N $J $R -K -O -P -V $src >> $psfile\n";
    print CSH "pstext -N $J_title $R_title -K -O -V -P >>$psfile<<EOF\n $x_title $z_title $fsize_title 0 $fontno CM $title \nEOF\n";
  }

#-----------------------------
  print CSH "pstext -N $J_title $R_title -O -V -P >>$psfile<<EOF\n 10 10 $fsize_title 0 $fontno CM $title \nEOF\n"; # FINISH
  print CSH "convert $psfile $jpgfile\n";
  print CSH "echo done with $psfile\n";

  close (CSH);
  system("csh -f $cshfile");
  system("xv $jpgfile &");
}

#===========================================================================
if ($ifig07 == 1) {

  # shifting the subplots
  $xfac = 1.15;
  $yfac = 1.15;
  $dX = $xfac*$Jwid; $dY = 0; $shift1 = "-X$dX -Y$dY";
  $dX = -2*$dX; $dY = -$yfac*$Jwid; $shift2 = "-X$dX -Y$dY";

  # colorbar
  $Dlen = 1.2; $Dx = $Jwid/2; $Dy = -0.2;
  $Dscale = "-D$Dx/$Dy/$Dlen/0.10h";

  # all models are m8
  $mod = $mods[16]; $smod = "m\@+$mod\@+";

  # get the receivers and sources
  $recfile = "$dir0/${edir}/sr.txt";  # src-rec for first event
  $evefile = "$dir0/events_lonlat.dat";
  if (not -f $recfile) { die("Check if $recfile exist or not\n") }
  if (not -f $evefile) { die("Check if $evefile exist or not\n") }

  # file names for figures
  $name    = "${smap}_event_kernel_${strun0}";
  $psfile  = "$name.ps";
  $jpgfile = "$name.jpg";

  #===============================================
  print "\nWriting CSH file...\n";

  open(CSH,">$cshfile");
  print CSH "gmtset BASEMAP_TYPE plain TICK_LENGTH $tick LABEL_FONT_SIZE $fsize2 ANOT_FONT_SIZE $fsize2 PLOT_DEGREE_FORMAT D HEADER_FONT $fontno ANOT_FONT $fontno LABEL_FONT $fontno HEADER_FONT_SIZE $fsize1 FRAME_PEN $fpen TICK_PEN $tpen\n";
  #===============================================
  # model for the data

  # base directory
  #$irun0 = 40;
  $strun = sprintf("%4.4i",$irun0); $dir = "$odir$strun";
  $B = "$B0".$Bopts[15];

  #-------------------------------------
  # color scale for EVENT kernels
  $kpwr = -7;
  $kmax = 2;
  $norm_ev = "1e$kpwr";
  $ss = $kmax;
  $ds = 2*$ss/$scale_color;
  $bs2 = sprintf("%3.3e",0.9*$ss); # colorbar
  $Tev = sprintf("-T%3.3e/%3.3e/%3.3e",-$ss,$ss,$ds);
  print "Tev = $Tev\n";

  $Bscale_ev  = sprintf("-B%2.2e:\" K ( \@~\161\@~, \@~\146\@~ )  ( 10\@+%2.2i\@+  m\@+-2\@+ s)\": -E7p",$bs2,$kpwr);

  $cpt_ev = "color5.cpt";
  print CSH "makecpt -C$colorbar $Tev > temp1\n";
  print CSH "sed 's/^B.*/B       170     0      0  /' temp1 >  temp2\n";
  print CSH "sed 's/^F.*/F         0     0    255  /' temp2 > $cpt_ev\n";
  #-------------------------------------

# these events follow the ORIGINAL ordering of events
#@ev_picks = (1,2,3,4,6,11,5);
@ev_picks = (1,2,3,4,5,6,7);
$nevent_plot = @ev_picks;
@labs = ("a","b","c","d","e","f","g","h","i");

# loop over events
for ($k = 0; $k < $nevent_plot; $k = $k+1) {

  # get the event kernel
  $ev = $k+1;
  $stev = sprintf("%2.2i",$ev);

  $ev_dir = sprintf("$dir/event_%3.3i",$ev_picks[$k]);
  $file2 = "$ev_dir/$kerfile";
  $srcfile = "$ev_dir/sr.txt";
  if (not -f $file2) { die("Check if $file2 exist or not\n") }
  if (not -f $srcfile) { die("Check if $srcfile exist or not\n") }

  # number of receivers
  $file3 = "$dir0/summed_chi_r.dat";
  if (not -f $file3) { die("Check if $file3 exist or not\n") }
  open(IN,$file3); @temp = <IN>; $nrec = @temp;
  $npath = $nrec * 1;
  print "\n nrec = $nrec, npath = $npath";

  $title = "($labs[$k])  Event kernel  $stev / $nevent";
  $mod = $mods[$k*2];
  if ($k % 3 == 0) {$shift = $shift2} else {$shift = $shift1}

  if($k==0) {print CSH "psbasemap $B $R $J -P -K -V $origin > $psfile\n"}
  else      {print CSH "psbasemap $B $R $J -P -K -O -V $shift >> $psfile\n"}
  if($icolor==1) {print CSH "awk '{print \$1,\$2,\$3 / $norm_ev }' $file2 | pscontour $J $R -A- -C$cpt_ev -I -K -O -P -V >> $psfile\n";}
  print CSH "pscoast $J $R $coast_info -P -K -O -V >> $psfile\n";
  if($ishelf) {print CSH "awk '{print \$2,\$1}' $shelf_file |psxy $J $R $Wshelf -K -O -P -V >> $psfile\n"}
  print CSH "awk '\$1 == \"R\" {print \$2,\$3}' $recfile |psxy -N $J $R -K -O -P -V $rec >> $psfile\n";
  print CSH "awk '\$1 == \"S\" {print \$2,\$3}' $srcfile |psxy -N $J $R -K -O -P -V $src_ev >> $psfile\n";
  print CSH "pstext -N $J_title $R_title -K -O -V -P >>$psfile<<EOF\n $x_title $z_title $fsize_title 0 $fontno CM $title \nEOF\n";
  if($k==$nevent_plot-1) {print CSH "psscale -C$cpt_ev $Dscale $Bscale_ev -K -O -P -V >> $psfile \n"}

}  # for k < nevent

  # summed event kernel

  # get the summed kernel and chi files
  $file2 = "$dir/summed_ker.dat";
  $recfile = "$dir/${edir}/sr.txt";     # src-rec for first event
  $evefile = "$dir0/events_lonlat.dat";
  if (not -f $file2)   { die("Check if $file2 exist or not\n") }
  if (not -f $recfile) { die("Check if $recfile exist or not\n") }

  # number of receivers
  #open(IN,$file3); @temp = <IN>; $nrec = @temp;
  $npath = $nrec * $nevent;
  $title = "($labs[$k])  Misfit kernel ($nevent events)";
  $shift = $shift1;

  print CSH "psbasemap $B $R $J -P -K -O -V $shift >> $psfile\n";
  if($icolor==1) {print CSH "awk '{print \$1,\$2,\$3 / $norm2 }' $file2 | pscontour $J $R -A- -C$cpt_ker -I -K -O -P -V >> $psfile\n"}
  print CSH "pscoast $J $R $coast_info -K -O -P -V >> $psfile\n";
  if($ishelf) {print CSH "awk '{print \$2,\$1}' $shelf_file |psxy $J $R $Wshelf -K -O -P -V >> $psfile\n"}
  print CSH "psscale -C$cpt_ker $Dscale $Bscale2 -K -O -P -V >> $psfile \n";
  print CSH "awk '{print \$1,\$2}' $evefile |psxy -N $J $R -K -O -P -V $src >> $psfile\n";
  print CSH "pstext -N $J_title $R_title -K -O -V -P >>$psfile<<EOF\n $x_title $z_title $fsize_title 0 $fontno CM $title \nEOF\n";
  print CSH "awk '\$1 == \"R\" {print \$2,\$3}' $recfile |psxy -N $J $R -K -O -P -V $rec >> $psfile\n";

  # phase speed for data
  $file1dat = "$dir/$mfile_dat";
  if (not -f $file1dat) { die("Check if $file1dat exist or not\n") }
  $title = "($labs[$k+1])  Target phase speed model (data)"; $shift = $shift1;

  print CSH "psbasemap $B $R $J -P -K -O -V $shift >> $psfile\n";
  if($icolor==1) {print CSH "awk '{print \$1,\$2,\$4*100}' $file1dat | pscontour $R $J -A- -C$cpt_vel -I -P -K -O -V  >> $psfile\n"}
  print CSH "pscoast $J $R $coast_info -P -K -O -V >> $psfile\n";
  if($ishelf) {print CSH "awk '{print \$2,\$1}' $shelf_file |psxy $J $R $Wshelf -K -O -P -V >> $psfile\n"}
  print CSH "awk '\$1 == \"R\" {print \$2,\$3}' $recfile |psxy -N $J $R -K -O -P -V $rec >> $psfile\n";
  print CSH "awk '{print \$1,\$2}' $evefile |psxy -N $J $R -K -O -P -V $src >> $psfile\n";
  print CSH "psscale -C$cpt_vel $Dscale $Bscale1 -P -K -O -V >> $psfile \n";
  print CSH "pstext -N $J_title $R_title -K -O -V -P >>$psfile<<EOF\n $x_title $z_title $fsize_title 0 $fontno CM $title \nEOF\n";

#-----------------------------
  print CSH "pstext -N $J_title $R_title -O -V -P >>$psfile<<EOF\n 10 10 $fsize_title 0 $fontno CM $title \nEOF\n"; # FINISH

  print CSH "convert $psfile $jpgfile\n";
  print CSH "echo done with $psfile\n";
  close (CSH);
  system("csh -f $cshfile");
  system("xv $jpgfile &");
}


#===========================================================================
if ($ifig08 == 1) {

  # shifting the subplots
  $xfac = 1.15;
  $yfac = 1.15;
  $dX = $xfac*$Jwid; $dY = 0; $shift1 = "-X$dX -Y$dY";
  $dX = -2*$dX; $dY = -$yfac*$Jwid; $shift2 = "-X$dX -Y$dY";

  # colorbar
  $Dlen = 1.2; $Dx = $Jwid/2; $Dy = -0.45;
  $Dscale = "-D$Dx/$Dy/$Dlen/0.10h";

  # all models are m8
  $mod = $mods[16]; $smod = "m\@+$mod\@+";

  #===============================================
  print "\nWriting CSH file...\n";

  open(CSH,">$cshfile");
  print CSH "gmtset BASEMAP_TYPE plain TICK_LENGTH $tick LABEL_FONT_SIZE $fsize2 ANOT_FONT_SIZE $fsize2 PLOT_DEGREE_FORMAT D HEADER_FONT $fontno ANOT_FONT $fontno LABEL_FONT $fontno HEADER_FONT_SIZE $fsize1 FRAME_PEN $fpen TICK_PEN $tpen\n";
  #===============================================

  # get the receivers and sources
  $recfile = "$dir0/${edir}/sr.txt";  # src-rec for first event
  $evefile = "$dir0/events_lonlat.dat";
  if (not -f $recfile) { die("Check if $recfile exist or not\n") }
  if (not -f $evefile) { die("Check if $evefile exist or not\n") }

  # file names for figures
  $name    = "nevent_maps";
  $psfile  = "$name.ps";
  $jpgfile = "$name.jpg";

  #=============================================
  # model for the data

  @runs  = (300,700,200,340,740,240,380,780,280);                       # irun0
  @ifile = (1,1,1,1,1,1,1,1,1);                                         # boolean: which runs are done
  @imods = (8,8,8,8,8,8,8,8,8);                                         # which model to plot
  @gams  = (30000,30000,30000,30000,30000,30000,30000,30000,30000);     # Gamma smoothing
  @nev   = (5,5,5,15,15,15,25,25,25);                                   # Nevent

  # phase speed model
  $k = 0;
  $irun0 = $runs[$k]; $strun = sprintf("%4.4i",$irun0); $dir = "$odir$strun"; $file1dat = "$dir/$mfile_dat";
  if (not -f $file1dat) { die("Check if $file1dat exist or not\n") }
  $B = "$B0".$Bopts[15]; $title = "Phase speed for data  (n = 3)";

  print CSH "psbasemap $B $R $J -P -K -V $origin > $psfile\n"; # START
  if($icolor==1) {print CSH "awk '{print \$1,\$2,\$4*100}' $file1dat | pscontour $R $J -A- -C$cpt_vel -I -P -K -O -V  >> $psfile\n"}
  print CSH "pscoast $J $R $coast_info -P -K -O -V >> $psfile\n";
  if($ishelf) {print CSH "awk '{print \$2,\$1}' $shelf_file |psxy $J $R $Wshelf -K -O -P -V >> $psfile\n"}
  print CSH "awk '\$1 == \"R\" {print \$2,\$3}' $recfile |psxy -N $J $R -K -O -P -V $rec >> $psfile\n";
  print CSH "awk '{print \$1,\$2}' $evefile |psxy -N $J $R -K -O -P -V $src >> $psfile\n";
  #print CSH "psscale -C$cpt_vel $Dscale $Bscale1 -P -K -O -V >> $psfile \n";
  print CSH "pstext -N $J_title $R_title -K -O -V -P >>$psfile<<EOF\n $x_title $z_title $fsize_title 0 $fontno CM $title \nEOF\n";

  # phase speed model
  $k = $k + 1;
  $irun0 = $runs[$k]; $strun = sprintf("%4.4i",$irun0); $dir = "$odir$strun"; $file1dat = "$dir/$mfile_dat";
  if (not -f $file1dat) { die("Check if $file1dat exist or not\n") }
  $B = "$B0".$Bopts[15]; $title = "Phase speed for data  (n = 1)";  $shift = $shift1;

  print CSH "psbasemap $B $R $J -P -K -O -V $shift >> $psfile\n";
  if($icolor==1) {print CSH "awk '{print \$1,\$2,\$4*100}' $file1dat | pscontour $R $J -A- -C$cpt_vel -I -P -K -O -V  >> $psfile\n"}
  print CSH "pscoast $J $R $coast_info -P -K -O -V >> $psfile\n";
  if($ishelf) {print CSH "awk '{print \$2,\$1}' $shelf_file |psxy $J $R $Wshelf -K -O -P -V >> $psfile\n"}
  print CSH "awk '\$1 == \"R\" {print \$2,\$3}' $recfile |psxy -N $J $R -K -O -P -V $rec >> $psfile\n";
  print CSH "awk '{print \$1,\$2}' $evefile |psxy -N $J $R -K -O -P -V $src >> $psfile\n";
  print CSH "pstext -N $J_title $R_title -K -O -V -P >>$psfile<<EOF\n $x_title $z_title $fsize_title 0 $fontno CM $title \nEOF\n";

  # phase speed model
  $k = $k + 1;
  $irun0 = $runs[$k]; $strun = sprintf("%4.4i",$irun0); $dir = "$odir$strun"; $file1dat = "$dir/$mfile_dat";
  if (not -f $file1dat) { die("Check if $file1dat exist or not\n") }
  $B = "$B0".$Bopts[15]; $title = "Phase speed for data  (Rayleigh, T = 20s)";  $shift = $shift1;

  print CSH "psbasemap $B $R $J -P -K -O -V $shift >> $psfile\n";
  if($icolor==1) {print CSH "awk '{print \$1,\$2,\$4*100}' $file1dat | pscontour $R $J -A- -C$cpt_vel_map -I -P -K -O -V  >> $psfile\n"}
  print CSH "pscoast $J $R $coast_info -P -K -O -V >> $psfile\n";
  if($ishelf) {print CSH "awk '{print \$2,\$1}' $shelf_file |psxy $J $R $Wshelf -K -O -P -V >> $psfile\n"}
  print CSH "awk '\$1 == \"R\" {print \$2,\$3}' $recfile |psxy -N $J $R -K -O -P -V $rec >> $psfile\n";
  print CSH "awk '{print \$1,\$2}' $evefile |psxy -N $J $R -K -O -P -V $src >> $psfile\n";
  print CSH "pstext -N $J_title $R_title -K -O -V -P >>$psfile<<EOF\n $x_title $z_title $fsize_title 0 $fontno CM $title \nEOF\n";

  #=========================

  $kmax = 8;
  for ($k = 0; $k <= $kmax; $k = $k+1) {

    if ($k % 3 == 0) {$shift = $shift2} else {$shift = $shift1}
    if ($k % 3 == 2) {$cpt = $cpt_vel_map} else {$cpt = $cpt_vel}

    # recovered phase speed model
    $irun0 = $runs[$k]; $Gam = $gams[$k];
    #$stsig = sprintf("\@~\163\@~ = %.1f km",$sigma/1000);
    $stGam = sprintf("\@~\107\@~ = %.1f km",$Gam/1000);
    $strun = sprintf("%4.4i",$irun0+2*$imods[$k]); $dir  = "$odir$strun";
    $strun0 = sprintf("%4.4i",$irun0);             $dir0 = "$odir$strun0";
    $file1syn = "$dir/$mfile_syn";
    $B = "$B0".$Bopts[15];
    $title = "Recovered  (m\@+$imods[$k]\@+,  Nevents = $nev[$k])";

    print CSH "psbasemap $B $R $J -P -K -O -V $shift >> $psfile\n";
    if($ifile[$k]==1) {
       if (not -f $file1syn) { die("Check if $file1syn exist or not\n") }
       if($icolor==1) {print CSH "awk '{print \$1,\$2,\$4*100}' $file1syn | pscontour $R $J -A- -C$cpt -I -P -K -O -V  >> $psfile\n"}

       $evefile = "$dir0/events_lonlat.dat";
       if (not -f $evefile) {die("Check if $evefile exist or not\n")}
     }
    print CSH "pscoast $J $R $coast_info -P -K -O -V >> $psfile\n";
    if($ishelf) {print CSH "awk '{print \$2,\$1}' $shelf_file |psxy $J $R $Wshelf -K -O -P -V >> $psfile\n"}
    print CSH "awk '\$1 == \"R\" {print \$2,\$3}' $recfile |psxy -N $J $R -K -O -P -V $rec >> $psfile\n";
    print CSH "awk '{print \$1,\$2}' $evefile |psxy -N $J $R -K -O -P -V $src >> $psfile\n";
    print CSH "pstext -N $J_title $R_title -K -O -V -P >>$psfile<<EOF\n $x_title $z_title $fsize_title 0 $fontno CM $title \nEOF\n";
  }

#-----------------------------
  print CSH "pstext -N $J_title $R_title -O -V -P >>$psfile<<EOF\n 10 10 $fsize_title 0 $fontno CM $title \nEOF\n"; # FINISH
  print CSH "convert $psfile $jpgfile\n";
  print CSH "echo done with $psfile\n";

  close (CSH);
  system("csh -f $cshfile");
  system("xv $jpgfile &");
}

#===========================================================================
if ($ifig08b == 1) {

  # shifting the subplots
  $xfac = 1.15;
  $yfac = 1.15;
  $dX = $xfac*$Jwid; $dY = 0; $shift1 = "-X$dX -Y$dY";
  $dX = -2*$dX; $dY = -$yfac*$Jwid; $shift2 = "-X$dX -Y$dY";

  # colorbar
  $Dlen = 1.2; $Dx = $Jwid/2; $Dy = -0.45;
  $Dscale = "-D$Dx/$Dy/$Dlen/0.10h";

  # all models are m8
  $mod = $mods[16]; $smod = "m\@+$mod\@+";

  #===============================================
  print "\nWriting CSH file...\n";

  open(CSH,">$cshfile");
  print CSH "gmtset BASEMAP_TYPE plain TICK_LENGTH $tick LABEL_FONT_SIZE $fsize2 ANOT_FONT_SIZE $fsize2 PLOT_DEGREE_FORMAT D HEADER_FONT $fontno ANOT_FONT $fontno LABEL_FONT $fontno HEADER_FONT_SIZE $fsize1 FRAME_PEN $fpen TICK_PEN $tpen\n";
  #===============================================

  # get the receivers and sources
  $recfile = "$dir0/${edir}/sr.txt";  # src-rec for first event
  $evefile = "$dir0/events_lonlat.dat";
  if (not -f $recfile) { die("Check if $recfile exist or not\n") }
  if (not -f $evefile) { die("Check if $evefile exist or not\n") }

  # file names for figures
  $name    = "nevent_maps_red";
  $psfile  = "$name.ps";
  $jpgfile = "$name.jpg";

  #=============================================
  # model for the data

  @runs  = (700,740,780);         # irun0
  @ifile = (1,1,1);               # boolean: which runs are done
  @imods = (8,8,8);               # which model to plot
  @gams  = (30000,30000,30000);   # Gamma smoothing
  @nev   = (5,15,25);             # Nevent
  $numk = @runs;

  #=========================

  $kmax = $numk;
  for ($k = 0; $k < $kmax; $k = $k+1) {

    # recovered phase speed model
    $irun0 = $runs[$k]; $Gam = $gams[$k];
    #$stsig = sprintf("\@~\163\@~ = %.1f km",$sigma/1000);
    $stGam = sprintf("\@~\107\@~ = %.1f km",$Gam/1000);
    $strun = sprintf("%4.4i",$irun0+2*$imods[$k]); $dir  = "$odir$strun";
    $strun0 = sprintf("%4.4i",$irun0);             $dir0 = "$odir$strun0";
    $file1syn = "$dir/$mfile_syn";

    # get number of sources and receivers
    $file3 = "$dir/summed_chi_r.dat";
    if (not -f $file3) {die("Check if $file3 exist or not\n")}
    open(IN,$file3); @temp = <IN>; $nrec = @temp;
    $nevent = $nev[$k];
    $N = $nevent * $nrec;

    print "\m nevent = $nevent, nrec = $nrec, N = $N \n";

    # load misfit associated with model
    $chifile = "$dir/summed_chi_all.dat";
    if (not -f $chifile) {die("Check if $chifile exist or not\n")}
    open(IN,"$chifile");
    $chi = <IN>;
    $schi = sprintf("\@~\143\@~ = %2.2e",$chi);
    #print "\n chi is $chi, misfit function is $schi\n";

    # compute average traveltime anomaly
    $del = sqrt( 2*$chi / $N );
    $sdel = sprintf("\@~\104\@~T = %.3f s",$del);
    #print "\n del is $del, average traveltime anomaly is $sdel seconds\n";

    $B = "$B0".$Bopts[15];
    $title = "m\@+$imods[$k]\@+,  Nevents = $nevent,  $sdel";

    if($k == 0){ print CSH "psbasemap $B $R $J -P -K -V $origin > $psfile\n"}
    else       { print CSH "psbasemap $B $R $J -P -K -O -V $shift1 >> $psfile\n"}
    if($ifile[$k]==1) {
       if (not -f $file1syn) { die("Check if $file1syn exist or not\n") }
       if($icolor==1) {print CSH "awk '{print \$1,\$2,\$4*100}' $file1syn | pscontour $R $J -A- -C$cpt_vel -I -P -K -O -V  >> $psfile\n"}

       $evefile = "$dir0/events_lonlat.dat";
       if (not -f $evefile) {die("Check if $evefile exist or not\n")}
     }
    print CSH "pscoast $J $R $coast_info -P -K -O -V >> $psfile\n";
    if($ishelf) {print CSH "awk '{print \$2,\$1}' $shelf_file |psxy $J $R $Wshelf -K -O -P -V >> $psfile\n"}
    print CSH "awk '\$1 == \"R\" {print \$2,\$3}' $recfile |psxy -N $J $R -K -O -P -V $rec >> $psfile\n";
    print CSH "awk '{print \$1,\$2}' $evefile |psxy -N $J $R -K -O -P -V $src >> $psfile\n";
    print CSH "pstext -N $J_title $R_title -K -O -V -P >>$psfile<<EOF\n $x_title $z_title $fsize_title 0 $fontno CM $title \nEOF\n";
  }

#-----------------------------
  print CSH "pstext -N $J_title $R_title -O -V -P >>$psfile<<EOF\n 10 10 $fsize_title 0 $fontno CM $title \nEOF\n"; # FINISH
  print CSH "convert $psfile $jpgfile\n";
  print CSH "echo done with $psfile\n";

  close (CSH);
  system("csh -f $cshfile");
  system("xv $jpgfile &");
}


#===========================================================================
if ($ifig09 == 1) {

  # shifting the subplots
  $xfac = 1.15;
  $yfac = 1.15;
  $dX = $xfac*$Jwid; $dY = 0; $shift1 = "-X$dX -Y$dY";
  $dX = -2*$dX; $dY = -$yfac*$Jwid; $shift2 = "-X$dX -Y$dY";

  # colorbar
  $Dlen = 1.2; $Dx = $Jwid/2; $Dy = -0.45;
  $Dscale = "-D$Dx/$Dy/$Dlen/0.10h";

  # all models are m8
  $mod = $mods[16]; $smod = "m\@+$mod\@+";

  #===============================================
  print "\nWriting CSH file...\n";

  open(CSH,">$cshfile");
  print CSH "gmtset BASEMAP_TYPE plain TICK_LENGTH $tick LABEL_FONT_SIZE $fsize2 ANOT_FONT_SIZE $fsize2 PLOT_DEGREE_FORMAT D HEADER_FONT $fontno ANOT_FONT $fontno LABEL_FONT $fontno HEADER_FONT_SIZE $fsize1 FRAME_PEN $fpen TICK_PEN $tpen\n";
  #===============================================

  # get the receivers and sources
  $recfile = "$dir0/${edir}/sr.txt";  # src-rec for first event
  $evefile = "$dir0/events_lonlat.dat";
  if (not -f $recfile) { die("Check if $recfile exist or not\n") }
  if (not -f $evefile) { die("Check if $evefile exist or not\n") }

  # file names for figures
  $name    = "nevent_kernels";
  $psfile  = "$name.ps";
  $jpgfile = "$name.jpg";

  #=============================================
  # model for the data

  @runs  = (300,700,200,340,740,240,380,780,280);                       # irun0
  @ifile = (1,1,1,1,1,1,1,1,1);                                         # boolean: which runs are done
  @imods = (0,0,0,0,0,0,0,0,0);                                         # which model to plot
  @gams  = (30000,30000,30000,30000,30000,30000,30000,30000,30000);     # Gamma smoothing
  @nev   = (5,5,5,15,15,15,25,25,25);                                   # Nevent

  # phase speed model
  $k = 0;
  $irun0 = $runs[$k]; $strun = sprintf("%4.4i",$irun0); $dir = "$odir$strun"; $file1dat = "$dir/$mfile_dat";
  if (not -f $file1dat) { die("Check if $file1dat exist or not\n") }
  $B = "$B0".$Bopts[15]; $title = "Phase speed for data  (N = 3)";

  print CSH "psbasemap $B $R $J -P -K -V $origin > $psfile\n"; # START
  if($icolor==1) {print CSH "awk '{print \$1,\$2,\$4*100}' $file1dat | pscontour $R $J -A- -C$cpt_vel -I -P -K -O -V  >> $psfile\n"}
  print CSH "pscoast $J $R $coast_info -P -K -O -V >> $psfile\n";
  if($ishelf) {print CSH "awk '{print \$2,\$1}' $shelf_file |psxy $J $R $Wshelf -K -O -P -V >> $psfile\n"}
  print CSH "awk '\$1 == \"R\" {print \$2,\$3}' $recfile |psxy -N $J $R -K -O -P -V $rec >> $psfile\n";
  print CSH "awk '{print \$1,\$2}' $evefile |psxy -N $J $R -K -O -P -V $src >> $psfile\n";
  #print CSH "psscale -C$cpt_vel $Dscale $Bscale1 -P -K -O -V >> $psfile \n";
  print CSH "pstext -N $J_title $R_title -K -O -V -P >>$psfile<<EOF\n $x_title $z_title $fsize_title 0 $fontno CM $title \nEOF\n";

  # phase speed model
  $k = $k + 1;
  $irun0 = $runs[$k]; $strun = sprintf("%4.4i",$irun0); $dir = "$odir$strun"; $file1dat = "$dir/$mfile_dat";
  if (not -f $file1dat) { die("Check if $file1dat exist or not\n") }
  $B = "$B0".$Bopts[15]; $title = "Phase speed for data  (N = 1)";  $shift = $shift1;

  print CSH "psbasemap $B $R $J -P -K -O -V $shift >> $psfile\n";
  if($icolor==1) {print CSH "awk '{print \$1,\$2,\$4*100}' $file1dat | pscontour $R $J -A- -C$cpt_vel -I -P -K -O -V  >> $psfile\n"}
  print CSH "pscoast $J $R $coast_info -P -K -O -V >> $psfile\n";
  if($ishelf) {print CSH "awk '{print \$2,\$1}' $shelf_file |psxy $J $R $Wshelf -K -O -P -V >> $psfile\n"}
  print CSH "awk '\$1 == \"R\" {print \$2,\$3}' $recfile |psxy -N $J $R -K -O -P -V $rec >> $psfile\n";
  print CSH "awk '{print \$1,\$2}' $evefile |psxy -N $J $R -K -O -P -V $src >> $psfile\n";
  print CSH "pstext -N $J_title $R_title -K -O -V -P >>$psfile<<EOF\n $x_title $z_title $fsize_title 0 $fontno CM $title \nEOF\n";

  # phase speed model
  $k = $k + 1;
  $irun0 = $runs[$k]; $strun = sprintf("%4.4i",$irun0); $dir = "$odir$strun"; $file1dat = "$dir/$mfile_dat";
  if (not -f $file1dat) { die("Check if $file1dat exist or not\n") }
  $B = "$B0".$Bopts[15]; $title = "Phase speed for data  (Rayleigh, T = 20s)";  $shift = $shift1;

  print CSH "psbasemap $B $R $J -P -K -O -V $shift >> $psfile\n";
  if($icolor==1) {print CSH "awk '{print \$1,\$2,\$4*100}' $file1dat | pscontour $R $J -A- -C$cpt_vel_map -I -P -K -O -V  >> $psfile\n"}
  print CSH "pscoast $J $R $coast_info -P -K -O -V >> $psfile\n";
  if($ishelf) {print CSH "awk '{print \$2,\$1}' $shelf_file |psxy $J $R $Wshelf -K -O -P -V >> $psfile\n"}
  print CSH "awk '\$1 == \"R\" {print \$2,\$3}' $recfile |psxy -N $J $R -K -O -P -V $rec >> $psfile\n";
  print CSH "awk '{print \$1,\$2}' $evefile |psxy -N $J $R -K -O -P -V $src >> $psfile\n";
  print CSH "pstext -N $J_title $R_title -K -O -V -P >>$psfile<<EOF\n $x_title $z_title $fsize_title 0 $fontno CM $title \nEOF\n";

  #=========================

  $kmax = 8;
  for ($k = 0; $k <= $kmax; $k = $k+1) {

    if ($k % 3 == 0) {$shift = $shift2} else {$shift = $shift1}

    # recovered phase speed model
    $irun0 = $runs[$k]; $Gam = $gams[$k];
    #$stsig = sprintf("\@~\163\@~ = %.1f km",$sigma/1000);
    $stGam = sprintf("\@~\107\@~ = %.1f km",$Gam/1000);
    $strun = sprintf("%4.4i",$irun0+2*$imods[$k]); $dir = "$odir$strun"; $file1ker = "$dir/summed_ker.dat";
    $B = "$B0".$Bopts[15];
    $title = "$Ktp for m\@+$imods[$k]\@+,  Nevents = $nev[$k]";

    print CSH "psbasemap $B $R $J -P -K -O -V $shift >> $psfile\n";
    if($ifile[$k]==1) {
       if (not -f $file1ker) { die("Check if $file1ker exist or not\n") }
       if($icolor==1) {print CSH "awk '{print \$1,\$2,\$3 / $norm2 }' $file1ker | pscontour $J $R -A- -C$cpt_ker -I -K -O -P -V >> $psfile\n"}

       $evefile = "$dir/events_lonlat.dat";
       if (not -f $evefile) {die("Check if $evefile exist or not\n")}
     }
    print CSH "pscoast $J $R $coast_info -P -K -O -V >> $psfile\n";
    if($ishelf) {print CSH "awk '{print \$2,\$1}' $shelf_file |psxy $J $R $Wshelf -K -O -P -V >> $psfile\n"}
    print CSH "awk '\$1 == \"R\" {print \$2,\$3}' $recfile |psxy -N $J $R -K -O -P -V $rec >> $psfile\n";
    print CSH "awk '{print \$1,\$2}' $evefile |psxy -N $J $R -K -O -P -V $src >> $psfile\n";
    print CSH "pstext -N $J_title $R_title -K -O -V -P >>$psfile<<EOF\n $x_title $z_title $fsize_title 0 $fontno CM $title \nEOF\n";
  }

#-----------------------------
  print CSH "pstext -N $J_title $R_title -O -V -P >>$psfile<<EOF\n 10 10 $fsize_title 0 $fontno CM $title \nEOF\n"; # FINISH
  print CSH "convert $psfile $jpgfile\n";
  print CSH "echo done with $psfile\n";

  close (CSH);
  system("csh -f $cshfile");
  system("xv $jpgfile &");
}


#===========================================================================
if ($ifig10 == 1) {

  # shifting the subplots
  $xfac = 1.20;
  $yfac = 1.50;
  $dX = $xfac*$Jwid; $dY = 0; $shift1 = "-X$dX -Y$dY";
  $dX = -2*$dX; $dY = -$yfac*$Jwid; $shift2 = "-X$dX -Y$dY";

  $origin2 = "-X0.8 -Y7.25";
  $x_title2 = 0.5;
  $z_title2 = 1.3;

  # colorbar
  $Dlen = 1.2; $Dx = $Jwid/2; $Dy = -0.10*$Jwid;
  $Dscale = "-D$Dx/$Dy/$Dlen/0.10h";

  # KEY: zoomed-in region
  $R2 = "-R-119/-116/33/36";

  #===============================================
  print "\nWriting CSH file...\n";

  open(CSH,">$cshfile");
  print CSH "gmtset BASEMAP_TYPE plain TICK_LENGTH $tick LABEL_FONT_SIZE $fsize2 ANOT_FONT_SIZE $fsize2 PLOT_DEGREE_FORMAT D HEADER_FONT $fontno ANOT_FONT $fontno LABEL_FONT $fontno HEADER_FONT_SIZE $fsize1 FRAME_PEN $fpen TICK_PEN $tpen\n";
  #===============================================

  for ($k= 0; $k <= $niter; $k = $k+1) {
    $irun = $irun0 + 2*$k;              # wave2d run number
    $strun = sprintf("%4.4i",$irun);
    $dir = "$odir$strun";
    $file1syn = "$dir/$mfile_syn";
    $mod_files[$k] = $file1syn;
    if (not -f $file1syn)   { die("Check if $file1syn exist or not\n") }
  }

  # get the receivers and sources
  $recfile = "$dir0/${edir}/sr.txt";  # src-rec for first event
  $evefile = "$dir0/events_lonlat.dat";
  if (not -f $recfile) { die("Check if $recfile exist or not\n") }
  if (not -f $evefile) { die("Check if $evefile exist or not\n") }

  # file names for figures
  $name    = "${smap}_cg04zoom_${strun0}";
  $psfile  = "$name.ps";
  $jpgfile = "$name.jpg";

  $B = "$B0".$Bopts[15];

  #=============================================
  # recovered phase speed map

  $niter_max = 3;
  @kvec = (0,1,2,3,$niter);
  @labs = ("a","b","c","d","e");
  $numk = @kvec;

  for ($i = 0; $i < $numk; $i = $i+1) {

    $k = $kvec[$i];
    $mod = $mods[$k*2]; $smod = "m\@+$mod\@+";
    $title = "($labs[$i])  Phase speed model $smod";
    if ($k== 0) {$title = "($labs[$i])  Initial phase speed model $smod"}
    if ($k % 3 == 0) {$shift = $shift2} else {$shift = $shift1}

    if ($k==0) {print CSH "psbasemap $B $R2 $J -P -K -V $origin2 > $psfile\n"}
    else       {print CSH "psbasemap $B $R2 $J -P -K -O -V $shift >> $psfile\n"}
    if ($icolor==1) {print CSH "awk '{print \$1,\$2,\$4*100}' $mod_files[$k] | pscontour $R2 $J -A- -C$cpt_vel -I -P -O -K -V >> $psfile\n"}
    print CSH "pscoast $J $R2 $coast_info -P -K -O -V >> $psfile\n";

    if ($k==0) {
      print CSH "psxy $plate_file $J $R2 $plate_info_k -K -V -O >> $psfile \n";
      print CSH "psxy $fault_file $J $R2 $fault_info_k -K -V -O >> $psfile \n";
      print CSH "awk '{print \$1,\$2}' $evefile |psxy $J $R2 -K -O -P -V $src -G255 >> $psfile\n";
      print CSH "pstext -N $J $R2 -K -V -O >>$psfile<<EOF
      -118.3 34.98 7  30 1 LM Garlock
      -118.5 34.82 7 -25 1 LM San Andreas
      -118.4 33.12 7  0 1 LM Pacific Ocean
EOF\n";
    } else {
      print CSH "awk '{print \$1,\$2}' $evefile |psxy $J $R2 -K -O -P -V $src >> $psfile\n";
    }

    print CSH "awk '{print \$2,\$1}' $shelf_file |psxy $J $R2 $Wshelf -K -O -P -V >> $psfile\n";
    print CSH "awk '\$1 == \"R\" {print \$2,\$3}' $recfile |psxy $J $R2 -K -O -P -V $rec >> $psfile\n";
    if($i == $numk-1) {print CSH "psscale -C$cpt_vel $Dscale $Bscale1 -P -K -O -V >> $psfile \n";}
    print CSH "pstext -N $J_title $R_title -K -O -V -P >>$psfile<<EOF\n $x_title2 $z_title2 $fsize_title 0 $fontno CM $title \nEOF\n";

    print CSH "psbasemap $B $R2 $J -P -K -O -V >> $psfile\n";
  }

  #=============================================
  # chi-vs-m curve

  $shift = $shift1;
  $Bg = "$B3b".$Bopts[8];
  $Jywid = $Jwid * 1.2;
  $J_chi = "-JX${Jwid}i/${Jywid}il";
  $J_title = "-JX${Jwid}i/${Jywid}i";
  $title = "(f)  Misfit,  \@~\143\@~ (m\@+k\@+)  ($utype)";

  for ($k = 0; $k <= $niter; $k = $k+1) {
    $irun = $irun0 + 2*$k;
    $strun = sprintf("%4.4i",$irun);
    $dir = "$odir$strun";
    $chifile = "$dir/summed_chi_all.dat";
    open(IN,"$chifile"); $chi = <IN>; chomp($chi);
    $it_vals[$k] = $k; $chi_vals[$k] = $chi;
  }

  print CSH "psbasemap $Bg $R_chi $J_chi -P -K -O -V $shift >> $psfile\n";
  if($irun0 != 100) {print CSH "awk '{print \$1,\$2}' $chi_curve | psxy $c_info_rd $J_chi $R_chi -K -O -P -V >> $psfile\n";}
  for ($k = 0; $k <= $niter; $k = $k+1) {
    print CSH "psxy $J_chi $R_chi $p_info_k -K -O -V -P >>$psfile<<EOF\n $it_vals[$k] $chi_vals[$k]\nEOF\n";  # plot point
  }

  if($irun0 == 100) {

    $irun0 = 20;
    for ($k = 0; $k <= $niter; $k = $k+1) {
      $irun = $irun0 + 2*$k;
      $strun = sprintf("%4.4i",$irun);
      $dir = "$odir$strun";
      $chifile = "$dir/summed_chi_all.dat";
      open(IN,"$chifile"); $chi = <IN>; chomp($chi);
      print CSH "psxy $J_chi $R_chi $p_info_w -K -O -V -P >>$psfile<<EOF\n $k $chi\nEOF\n";  # plot point
    }

 # legend
  $origin_scale = "-Xa0.6 -Ya1.95";            # position of origin (inches)
  $bxmin = 0; $bymin = 0;                     # location of lower left box corner
  $bwid = 0.70;                               # width of box
  $nrow = 2;                                  # KEY: number of rows in legend
  $dy = 0.07;                                 # spacing for each row
  $dx1 = 0.05;                                # spacing between left box edge and label
  $dx2 = 0.05;                                # spacing between label and text
  $dy1 = $dx1;                                # spacing between bottom box edge and label
  $bxmax = $bxmin + $bwid;                    # position of right edge
  $bymax = $bymin + 2*$dy1 + ($nrow-1)*$dy;   # position of top edge
  $lx = $bxmin + $dx1;                        # x position of label
  $tx = $lx + $dx2;                           # x position of text
  $ry1 = $bymin + $dy1;                       # y position of bottom row
  $ry2 = $bymin + $dy1 + $dy;
  $ry3 = $bymin + $dy1 + 2*$dy;
  $ry4 = $bymin + $dy1 + 3*$dy;

  # legend labels (bottom row to top)
  @legendlabels = ("cubic interpolation","quadratic interpolation");

  print CSH "psxy -W1.0p $J_title $R_title -K -O -V -P $origin_scale >>$psfile<<EOF\n$bxmin $bymin\n$bxmax $bymin\n$bxmax $bymax\n$bxmin $bymax\n$bxmin $bymin\nEOF\n";
  print CSH "psxy $p_info_w $J_title $R_title -K -O -V -P $origin_scale >>$psfile<<EOF\n $lx $ry2 \nEOF\n";
  print CSH "psxy $p_info_k $J_title $R_title -K -O -V -P $origin_scale >>$psfile<<EOF\n $lx $ry1 \nEOF\n";
  print CSH "pstext -N $J_title $R_title -K -O -V -P $origin_scale >>$psfile<<EOF\n $tx $ry2 9 0 $fontno LM $legendlabels[1] \nEOF\n";
  print CSH "pstext -N $J_title $R_title -K -O -V -P $origin_scale >>$psfile<<EOF\n $tx $ry1 9 0 $fontno LM $legendlabels[0] \nEOF\n";

  }

  print CSH "pstext -N $J_chi $R_chi -K -O -V -P >>$psfile<<EOF\n $x0chit $y1t $fsize1 0 $fontno BL $schi0 \nEOF\n";
  print CSH "pstext -N $J_chi $R_chi -K -O -V -P >>$psfile<<EOF\n $x1chit $y7t $fsize1 0 $fontno BL $schi1 \nEOF\n";
  print CSH "pstext -N $J_title $R_title -K -O -V -P >>$psfile<<EOF\n $x_title $z_title $fsize_title 0 $fontno CM $title \nEOF\n";

#  #=============================================
#  # phase speed map for data

#  #$B = "$B0".$Bopts[15];
#  $title = "(f)  Phase speed for data";
#  $shift = $shift1;

#  print CSH "psbasemap $B $R2 $J -P -K -O -V $shift >> $psfile\n";
#  if($icolor==1) {print CSH "awk '{print \$1,\$2,\$4*100}' $file1dat | pscontour $R2 $J -A- -C$cpt_vel -I -P -K -O -V  >> $psfile\n"}
#  print CSH "pscoast $J $R2 $coast_info -P -K -O -V >> $psfile\n";
#  print CSH "awk '{print \$2,\$1}' $shelf_file |psxy $J $R2 $Wshelf -K -O -P -V >> $psfile\n";
#  print CSH "awk '\$1 == \"R\" {print \$2,\$3}' $recfile |psxy $J $R2 -K -O -P -V $rec >> $psfile\n";
#  print CSH "awk '{print \$1,\$2}' $evefile |psxy $J $R2 -K -O -P -V $src >> $psfile\n";
#  print CSH "pstext -N $J_title $R_title -K -O -V -P >>$psfile<<EOF\n $x_title2 $z_title2 $fsize_title 0 $fontno CM $title \nEOF\n";
#  print CSH "psbasemap $B $R2 $J -P -K -O -V >> $psfile\n";

#-----------------------------
  print CSH "pstext -N $J_title $R_title -O -V -P >>$psfile<<EOF\n 10 10 $fsize_title 0 $fontno CM $title \nEOF\n"; # FINISH
  print CSH "convert $psfile $jpgfile\n";
  print CSH "echo done with $psfile\n";

  close (CSH);
  system("csh -f $cshfile");
  system("xv $jpgfile &");
}



#===========================================================================
if ($ifig10b == 1) {

  # shifting the subplots
  $xfac = 1.20;
  $yfac = 1.50;
  $dX1 = $xfac*$Jwid; $dY1 = 0; $shift1 = "-X$dX1 -Y$dY1";
  $dX1b = (1+($xfac-1)*2)*$Jwid; $dY1b = 0; $shift1b = "-X$dX1b -Y$dY1b";
  $dX2 = -2*$dX1; $dY2 = -$yfac*$Jwid; $shift2 = "-X$dX2 -Y$dY2";

  $origin2 = "-X0.8 -Y7.25";
  $x_title2 = 0.5;
  $z_title2 = 1.3;

  # colorbar
  $Dlen = 1.2; $Dx = $Jwid/2; $Dy = -0.10*$Jwid;
  $Dscale = "-D$Dx/$Dy/$Dlen/0.10h";

  # KEY: zoomed-in region
  $R2 = "-R-119/-116/33/36";

  #===============================================
  print "\nWriting CSH file...\n";

  open(CSH,">$cshfile");
  print CSH "gmtset BASEMAP_TYPE plain TICK_LENGTH $tick LABEL_FONT_SIZE $fsize1 ANOT_FONT_SIZE $fsize2 PLOT_DEGREE_FORMAT D HEADER_FONT $fontno ANOT_FONT $fontno LABEL_FONT $fontno HEADER_FONT_SIZE $fsize1 FRAME_PEN $fpen TICK_PEN $tpen\n";
  #===============================================

  $k = $niter;
  $irun = $irun0 + 2*$k;
  $strun = sprintf("%4.4i",$irun);
  $dir = "$odir$strun";
  $file1syn = "$dir/$mfile_syn";
  $mod_file = $file1syn;
  if (not -f $file1syn)   { die("Check if $file1syn exist or not\n") }

  # get the receivers and sources
  $recfile = "$dir0/${edir}/sr.txt";  # src-rec for first event
  $evefile = "$dir0/events_lonlat.dat";
  if (not -f $recfile) { die("Check if $recfile exist or not\n") }
  if (not -f $evefile) { die("Check if $evefile exist or not\n") }

  # file names for figures
  $name    = "${smap}_zoom_gji";
  $psfile  = "$name.ps";
  $jpgfile = "$name.jpg";

  $B = "$B0".$Bopts[15];

  #=============================================
  # phase speed map for data

  $B = "$B0".$Bopts[15];
  $title = "(a)  Target model (data)";

  print CSH "psbasemap $B $R2 $J -P -K -V $origin2 > $psfile\n";  # START
  if($icolor==1) {print CSH "awk '{print \$1,\$2,\$4*100}' $file1dat | pscontour $R2 $J -A- -C$cpt_vel -I -P -K -O -V  >> $psfile\n"}
  print CSH "pscoast $J $R2 $coast_info -P -K -O -V >> $psfile\n";
  print CSH "awk '\$1 == \"R\" {print \$2,\$3}' $recfile |psxy $J $R2 -K -O -P -V $rec >> $psfile\n";
  if(0 == 1) {
     print CSH "psxy $plate_file $J $R2 $plate_info_k -K -V -O >> $psfile \n";
     print CSH "psxy $fault_file $J $R2 $fault_info_k -K -V -O >> $psfile \n";
     print CSH "awk '{print \$1,\$2}' $evefile |psxy $J $R2 -K -O -P -V $src -G255 >> $psfile\n";
     print CSH "pstext -N $J $R2 -K -V -O >>$psfile<<EOF
     -118.3 34.98 7  30 1 LM Garlock
     -118.5 34.82 7 -25 1 LM San Andreas
     -118.4 33.12 7  0 1 LM Pacific Ocean
EOF\n";
  } else {
     print CSH "awk '{print \$1,\$2}' $evefile |psxy $J $R2 -K -O -P -V $src >> $psfile\n";
  }
  print CSH "pstext -N $J_title $R_title -K -O -V -P >>$psfile<<EOF\n $x_title2 $z_title2 $fsize_title 0 $fontno CM $title \nEOF\n";
  print CSH "psbasemap $B $R2 $J -P -K -O -V >> $psfile\n";

  #=============================================
  # phase speed map from rays

  # both have the same damping value
  $dir_hess = "/home/carltape/wave2d/2d_adjoint_banana/mat_SAVE/recovered_model_vector";
  $stgam = sprintf("\@~\147\@~ = %.2f", 10.0);

  $ray_file = "${dir_hess}/run_0020/lcurve50/m_IMODEL_3_rayleigh_id03_plot";
  if (not -f $ray_file)   { die("Check if $ray_file exist or not\n") }
  $title = "(b)  Model from rays, $stgam";
  $shift = $shift1;

  print CSH "psbasemap $B $R2 $J -P -K -O -V $shift >> $psfile\n";
  if ($icolor==1) {print CSH "awk '{print \$1,\$2,\$3*100}' $ray_file | pscontour $R2 $J -A- -C$cpt_vel -I -P -O -K -V >> $psfile\n"}
  print CSH "pscoast $J $R2 $coast_info -P -K -O -V >> $psfile\n";
  print CSH "awk '{print \$1,\$2}' $evefile |psxy $J $R2 -K -O -P -V $src >> $psfile\n";
  print CSH "awk '\$1 == \"R\" {print \$2,\$3}' $recfile |psxy $J $R2 -K -O -P -V $rec >> $psfile\n";
  print CSH "pstext -N $J_title $R_title -K -O -V -P >>$psfile<<EOF\n $x_title2 $z_title2 $fsize_title 0 $fontno CM $title \nEOF\n";
  print CSH "psbasemap $B $R2 $J -P -K -O -V >> $psfile\n";

  #=============================================
  # phase speed map from kernels

  $ker_file = "${dir_hess}/run_0020/lcurve05/m_IMODEL_3_rayleigh_id03_plot";
  if (not -f $ker_file)   { die("Check if $ker_file exist or not\n") }

  $title = "(c)  Model from kernels, $stgam";
  $shift = $shift1;

  print CSH "psbasemap $B $R2 $J -P -K -O -V $shift >> $psfile\n";
  if ($icolor==1) {print CSH "awk '{print \$1,\$2,\$3*100}' $ker_file | pscontour $R2 $J -A- -C$cpt_vel -I -P -O -K -V >> $psfile\n"}
  print CSH "pscoast $J $R2 $coast_info -P -K -O -V >> $psfile\n";
  print CSH "awk '{print \$1,\$2}' $evefile |psxy $J $R2 -K -O -P -V $src >> $psfile\n";
  print CSH "awk '\$1 == \"R\" {print \$2,\$3}' $recfile |psxy $J $R2 -K -O -P -V $rec >> $psfile\n";
  print CSH "pstext -N $J_title $R_title -K -O -V -P >>$psfile<<EOF\n $x_title2 $z_title2 $fsize_title 0 $fontno CM $title \nEOF\n";
  print CSH "psbasemap $B $R2 $J -P -K -O -V >> $psfile\n";

  #=============================================
  # phase speed map from adjoint-CG

  $mod = $mods[$k*2]; $smod = "m\@+$mod\@+";
  $title = "(d)  Model from adjoint, $smod";
  $shift = $shift2;

  print CSH "psbasemap $B $R2 $J -P -K -O -V $shift >> $psfile\n";
  if ($icolor==1) {print CSH "awk '{print \$1,\$2,\$4*100}' $mod_file | pscontour $R2 $J -A- -C$cpt_vel -I -P -O -K -V >> $psfile\n"}
  print CSH "pscoast $J $R2 $coast_info -P -K -O -V >> $psfile\n";
  print CSH "awk '{print \$1,\$2}' $evefile |psxy $J $R2 -K -O -P -V $src >> $psfile\n";
  print CSH "awk '\$1 == \"R\" {print \$2,\$3}' $recfile |psxy $J $R2 -K -O -P -V $rec >> $psfile\n";
  #print CSH "psscale -C$cpt_vel $Dscale $Bscale1 -P -K -O -V >> $psfile \n";
  print CSH "psscale -C$cpt_vel $Dscale $Bscale1b -P -K -O -V >> $psfile \n";
  print CSH "pstext -N $J $R -K -O -V -P >>$psfile<<EOF\n-115.2 30.6 11 0 $fontno LM \@~\045\@~\nEOF\n";
  print CSH "pstext -N $J_title $R_title -K -O -V -P >>$psfile<<EOF\n $x_title2 $z_title2 $fsize_title 0 $fontno CM $title \nEOF\n";
  print CSH "psbasemap $B $R2 $J -P -K -O -V >> $psfile\n";

  #=============================================
  # chi-vs-m curve

  $shift = $shift1b;

  $Bg = "-Ba2f1:\" k, model number \":/a1f2p:\"  Misfit,  \@~\143\@~ (m\@+k\@+)  (s\@+2\@+)  \":".$Bopts[8];

  $Jxwid = $Jwid * 2.0;
  $Jywid = $Jwid * 1.2;
  $J_chi = "-JX${Jxwid}i/${Jywid}il";
  $J_title = "-JX${Jxwid}i/${Jywid}i";
  $title = "(e)  Misfit comparison: classical vs adjoint";

  # load ray chi
  $dir_ray = "/home/store2/carltape/OUTPUT_banana/run_0023";
  $chifile = "${dir_ray}/summed_chi_all.dat";
  open(IN,"$chifile"); $chi_ray = <IN>; chomp($chi_ray);
  print "\n chi_ray = ${chi_ray} \n";

  # load kernel chi
  $dir_ker = "/home/store2/carltape/OUTPUT_banana/run_0022";
  $chifile = "${dir_ker}/summed_chi_all.dat";
  open(IN,"$chifile"); $chi_ker = <IN>; chomp($chi_ker);
  print "\n chi_ker = ${chi_ker} \n";

  for ($k = 0; $k <= $niter; $k = $k+1) {
    $irun = $irun0 + 2*$k;
    $strun = sprintf("%4.4i",$irun);
    $dir = "$odir$strun";
    $chifile = "$dir/summed_chi_all.dat";
    open(IN,"$chifile"); $chi = <IN>; chomp($chi);
    $it_vals[$k] = $k; $chi_vals[$k] = $chi;
  }

  print CSH "psbasemap $Bg $R_chi $J_chi -P -K -O -V $shift >> $psfile\n";

  # rays and kernels
  print CSH "psxy  -W1tap   $J_chi $R_chi -K -O -V -P >>$psfile<<EOF\n$a2x1 $chi_ker\n$a2x2 $chi_ker\nEOF\n";
  print CSH "psxy  -W1p     $J_chi $R_chi -K -O -V -P >>$psfile<<EOF\n$a2x1 $chi_ray\n$a2x2 $chi_ray\nEOF\n";

  # text labels
  print CSH "pstext -N $J_title $R_title -K -O -V -P >>$psfile<<EOF\n 0.02 0.40 10 0 $fontno LM (b) Hessian - Ray \nEOF\n";
  print CSH "pstext -N $J_title $R_title -K -O -V -P >>$psfile<<EOF\n 0.02 0.29 10 0 $fontno LM (c) Hessian - Kernel \nEOF\n";
  print CSH "pstext -N $J_title $R_title -K -O -V -P >>$psfile<<EOF\n 0.86 0.08 10 0 $fontno LM (d)  $smod  \nEOF\n";

  print CSH "awk '{print \$1,\$2}' $chi_curve | psxy $c_info_rd $J_chi $R_chi -K -O -P -V >> $psfile\n";
  for ($k = 0; $k <= $niter; $k = $k+1) {
    print CSH "psxy $J_chi $R_chi $p_info_k -K -O -V -P >>$psfile<<EOF\n $it_vals[$k] $chi_vals[$k]\nEOF\n";  # plot point
  }
  #print CSH "pstext -N $J_chi $R_chi -K -O -V -P >>$psfile<<EOF\n $x0chit $y1t $fsize1 0 $fontno BL $schi0 \nEOF\n";
  #print CSH "pstext -N $J_chi $R_chi -K -O -V -P >>$psfile<<EOF\n $x1chit $y7t $fsize1 0 $fontno BL $schi1 \nEOF\n";
  print CSH "pstext -N $J_title $R_title -K -O -V -P >>$psfile<<EOF\n $x_title $z_title $fsize_title 0 $fontno CM $title \nEOF\n";

#-----------------------------
  print CSH "pstext -N $J_title $R_title -O -V -P >>$psfile<<EOF\n 10 10 $fsize_title 0 $fontno CM $title \nEOF\n"; # FINISH
  print CSH "convert $psfile $jpgfile\n";
  print CSH "echo done with $psfile\n";

  close (CSH);
  system("csh -f $cshfile");
  system("xv $jpgfile &");
}




#===========================================================================
if ($ifig20 == 1) {

  # shifting the subplots
  $xfac = 1.15;
  $yfac = 1.15;
  $dX0 = $xfac*$Jwid; $dY0 = -$yfac*$Jwid;
  $dX1 = $dX0;         $dY1 = 0;     $shift1  = "-X$dX1  -Y$dY1";
  $dX2 = -2*$dX1;      $dY2 = $dY0;  $shift2  = "-X$dX2  -Y$dY2";

  # colorbar
  $Dlen = 1.2; $Dx = $Jwid/2; $Dy = -0.08*$Jwid;
  $Dscale = "-D$Dx/$Dy/$Dlen/0.10h";

  #===============================================
  print "\nWriting CSH file...\n";

  open(CSH,">$cshfile");
  print CSH "gmtset BASEMAP_TYPE plain TICK_LENGTH $tick LABEL_FONT_SIZE $fsize2 ANOT_FONT_SIZE $fsize2 PLOT_DEGREE_FORMAT D HEADER_FONT $fontno ANOT_FONT $fontno LABEL_FONT $fontno HEADER_FONT_SIZE $fsize1 FRAME_PEN $fpen TICK_PEN $tpen\n";
  #===============================================

  # KEY COMMAND
  $ilcurve = 50;   # 50-51 (rays), 5 (kernels)
  $stil = sprintf("%2.2i",$ilcurve);

  $recover_dir = "/home/carltape/wave2d/2d_adjoint_banana/mat_SAVE/recovered_model_vector/run_${strun0}/lcurve${stil}";
  $model_lab   = "IMODEL_2_Nfac${sNfac}";
  if($istructure==4) {$model_lab = "IMODEL_3_rayleigh"}

  # find in ALL plotting files
  @mod_files = glob("${recover_dir}/*_plot");
  $nplot = @mod_files;
  print "\n $nplot plotting files";
  print "\n first file : ${mod_files[0]} \n";

  # load in ALL the damping values
  open(IN,"${recover_dir}/lamvec"); @lams = <IN>;
  $ndamp = @lams;
  print "\n $ndamp damping values for the recovered models \n";

  # for the ray-based inversions, the no-damping case does not exist
  $nbad = $ndamp - $nplot;
  print "\n $nbad files that we do not use \n";
  #die("testing");

  #------------------
  # load specs for L-curve plots

  # axes limits
  $axes_file = "${recover_dir}/lcurve_axis";
  if (not -f ${axes_file}) { die("Check if ${axes_file} exist or not\n") }
  open(IN,"$axes_file"); @ax0 = <IN>;
  ($a1x1,$a1x2,$a1y1,$a1y2) = split(" ",$ax0[0]);

  $normgx = 1;
  $normgy = 1;
  $xming = $a1x1/$normgx;
  $xmaxg = $a1x2/$normgx;
  $yming = $a1y1/$normgy;
  $ymaxg = $a1y2/$normgy;
  $R_Lcurve  = "-R$xming/$xmaxg/$yming/$ymaxg";
  $J_Lcurve  = "-JX$Jwid";

  $dx = $xmaxg - $xming;
  $dy = $ymaxg - $yming;
  $xtick = 0.25*$dx;
  $ytick = 0.25*$dy;
  #$B_Lcurve  = sprintf("-B%.2f:\" log10 ( Misfit Norm ) \":/%.2f:\" log10 ( Model Norm )\":",$xtick,$ytick);
  $B_Lcurve  = sprintf("-B%.2f:\" log10 ( Model Norm ) \":/%.2f:\" log10 ( Misfit Norm )\":",$xtick,$ytick);

  # polynomial lines and curves
  $Lcurve_file = "${recover_dir}/lcurve_pts";
  if (not -f ${Lcurve_file}) { die("Check if ${Lcurve_file} exist or not\n") }
  open(IN,"${Lcurve_file}"); @Lcurve_pts = <IN>;

  #------------------

  #for ($k = 1; $k <= $ndamp; $k = $k+1) {
  #  $stid = sprintf("%2.2i",$k);
  #  $file1syn = "${recover_dir}m_${model_lab}_id${stid}_plot";
  #  $mod_files[$k-1] = $file1syn;
  #  if (not -f $file1syn)   { die("Check if $file1syn exist or not\n") }
  #}

  # get the receivers and sources
  $recfile = "$dir0/${edir}/sr.txt";  # src-rec for first event
  $evefile = "$dir0/events_lonlat.dat";
  if (not -f $recfile) { die("Check if $recfile exist or not\n") }
  if (not -f $evefile) { die("Check if $evefile exist or not\n") }

  # file names for figures
  # hess refers to the fact that we have used the Hessian (GtG) in the computations,
  # following the classical tomography approach
  $name    = "${smap}_hess01_${strun0}_ilcurve${stil}";
  $psfile  = "$name.ps";
  $jpgfile = "$name.jpg";

  $B = "$B0".$Bopts[15];

  #=============================================
  # recovered phase speed map

  # pick a subset of the possible figures
  print "\n pick a subset of the models to plot \n";
  @kvec = (0,2,4,6,8,10);
  #@kvec = (0,4,8,12,16,20);
  #@kvec = (0,1,2,4,6,$nplot-1);
  #@kvec = (0,1,2,3,4,5);
  $numk = @kvec;

  for ($i = 0; $i < $numk; $i = $i+1) {

    $k = $kvec[$i];
    $stgam = sprintf("\@~\147\@~ = %.2f",$lams[$k+$nbad]);
    $title = "Model for $stgam";
    if ($i % 3 == 0) {$shift = $shift2} else {$shift = $shift1}
    if (not -f $mod_files[$k]) { die("Check if $mod_files[$k] exist or not\n") }

    if ($k==0) {
       print CSH "psbasemap $B $R $J -P -K -V $origin > $psfile\n"; # START
    } else {
       print CSH "psbasemap $B $R $J -P -K -O -V $shift >> $psfile\n";
    }
    # note 3rd column, not 4th
    if ($icolor==1) {print CSH "awk '{print \$1,\$2,\$3*100}' $mod_files[$k] | pscontour $R $J -A- -C$cpt_vel -I -P -O -K -V >> $psfile\n"}
    print CSH "pscoast $J $R $coast_info -P -K -O -V >> $psfile\n";
    if($ishelf) {print CSH "awk '{print \$2,\$1}' $shelf_file |psxy $J $R $Wshelf -K -O -P -V >> $psfile\n"}
    print CSH "awk '{print \$1,\$2}' $evefile |psxy $J $R -K -O -P -V $src >> $psfile\n";
    print CSH "awk '\$1 == \"R\" {print \$2,\$3}' $recfile |psxy $J $R -K -O -P -V $rec >> $psfile\n";
    #print CSH "psscale -C$cpt_vel $Dscale $Bscale1 -P -K -O -V >> $psfile \n";
    print CSH "pstext -N $J_title $R_title -K -O -V -P >>$psfile<<EOF\n $x_title $z_title $fsize_title 0 $fontno CM $title \nEOF\n";

    print CSH "psbasemap $B $R $J -P -K -O -V >> $psfile\n";
  }

  #=============================================
  # Lcurve

  $B = "$B_Lcurve".$Bopts[8];  # 15 for no axes
  $title = "L-curve for the models";

  $col = $numk % 3;
  if($col == 1) {$shift3 = "-X$dX1 -Y$dY2"}
  if($col == 2) {$shift3 = "-Y$dY2"}
  if($col == 0) {$shift3 = "-X-$dX1 -Y$dY2"}
  $shift = $shift3;

  print CSH "psbasemap $B $R_Lcurve $J_Lcurve -P -K -O -V $shift >> $psfile\n";
  print CSH "psxy $Lcurve_file $R_Lcurve $J_Lcurve $p_info_k -P -K -O -V >> $psfile\n";
  print CSH "psxy $Lcurve_file $R_Lcurve $J_Lcurve -W0.75p/0/0/255 -P -K -O -V >> $psfile\n";

  # print a couple labels
  #@kvec2 = (0,1,4);
  @kvec2 = @kvec;
  $numk2 = @kvec2;
  for ($i = 0; $i < $numk2; $i = $i+1) {
    $k = $kvec2[$i];
    #$stgam = sprintf("\@~\147\@~ = %.2f",$lams[$k]);
    $stgam = sprintf("%.2f",$lams[$k]);
    ($xpt,$ypt) = split(" ",$Lcurve_pts[$k]);
    $xpt = $xpt + 0.03*$dx;
    $ypt = $ypt + 0.03*$dy;
    print CSH "pstext -N $J_Lcurve $R_Lcurve -K -O -V -P >>$psfile<<EOF\n $xpt $ypt $fsize_title 0 $fontno LM $stgam \nEOF\n";
  }

  print CSH "pstext -N $J_title $R_title -K -O -V -P >>$psfile<<EOF\n $x_title $z_title $fsize_title 0 $fontno CM $title \nEOF\n";

  #=============================================
  # phase speed map for data

  $B = "$B0".$Bopts[15];
  $title = "Target phase speed model (data)";
  $shift = $shift1;

  print CSH "psbasemap $B $R $J -P -K -O -V $shift >> $psfile\n";
  if($icolor==1) {print CSH "awk '{print \$1,\$2,\$4*100}' $file1dat | pscontour $R $J -A- -C$cpt_vel -I -P -K -O -V  >> $psfile\n"}
  print CSH "pscoast $J $R $coast_info -P -K -O -V >> $psfile\n";
  if($ishelf) {print CSH "awk '{print \$2,\$1}' $shelf_file |psxy $J $R $Wshelf -K -O -P -V >> $psfile\n"}
  print CSH "awk '\$1 == \"R\" {print \$2,\$3}' $recfile |psxy $J $R -K -O -P -V $rec >> $psfile\n";
  print CSH "awk '{print \$1,\$2}' $evefile |psxy $J $R -K -O -P -V $src >> $psfile\n";
  print CSH "psscale -C$cpt_vel $Dscale $Bscale1b -P -K -O -V >> $psfile \n";
  print CSH "pstext -N $J_title $R_title -K -O -V -P >>$psfile<<EOF\n $x_title $z_title $fsize_title 0 $fontno CM $title \nEOF\n";
  print CSH "psbasemap $B $R $J -P -K -O -V >> $psfile\n";

#-----------------------------
  print CSH "pstext -N $J_title $R_title -O -V -P >>$psfile<<EOF\n 10 10 $fsize_title 0 $fontno CM $title \nEOF\n"; # FINISH
  print CSH "convert $psfile $jpgfile\n";
  print CSH "echo done with $psfile\n";

  close (CSH);
  system("csh -f $cshfile");
  system("xv $jpgfile &");
}

#=================================================================
if ($ifig20b == 1) {

  # shifting the subplots
  $xfac = 1.35;
  $yfac = 1.25;
  $dX0 = $xfac*$Jwid; $dY0 = -$yfac*$Jwid;
  $dX = $dX0;         $dY = 0;     $shift1  = "-X$dX  -Y$dY";
  $dX = -2*$dX;       $dY = $dY0;  $shift2  = "-X$dX  -Y$dY";
  $dXb = -1.5*$dX0;   $dY = $dY0;  $shift2b = "-X$dXb -Y$dY";
  $origin = "-X0.5 -Y8";

  # colorbar
  $Dlen = 1.2; $Dx = $Jwid/2; $Dy = -0.05*$Jwid;
  $Dscale = "-D$Dx/$Dy/$Dlen/0.10h";

  #===============================================
  print "\nWriting CSH file...\n";

  open(CSH,">$cshfile");
  print CSH "gmtset BASEMAP_TYPE plain TICK_LENGTH $tick LABEL_FONT_SIZE $fsize2 ANOT_FONT_SIZE $fsize2 PLOT_DEGREE_FORMAT D HEADER_FONT $fontno ANOT_FONT $fontno LABEL_FONT $fontno HEADER_FONT_SIZE $fsize1 FRAME_PEN $fpen TICK_PEN $tpen\n";
  #===============================================

  $ilcurve = 11;
  $stil = sprintf("%2.2i",$ilcurve);

  $recover_dir = "/home/carltape/wave2d/2d_adjoint_banana/mat_SAVE/recovered_model_vector/run_${strun0}/lcurve${stil}";
  $model_lab   = "IMODEL_2_Nfac${sNfac}";
  if($istructure==4) {$model_lab = "IMODEL_3_rayleigh"}

  # load in ALL the damping values
  open(IN,"${recover_dir}/lamvec");
  @lams = <IN>;
  $ndamp = @lams;
  print "\n $ndamp damping values for the recovered models \n";

  #------------------
  # L-curve plots

  # axes limits
  $axes_file = "${recover_dir}/lcurve_axis";
  if (not -f $axes_file)   { die("Check if $axes_file exist or not\n") }
  open(IN,"$axes_file");
  @ax0 = <IN>;
  ($a1x1,$a1x2,$a1y1,$a1y2) = split(" ",$ax0[0]);

  $normgx = 1;
  $normgy = 1;
  $xming = $a1x1/$normgx;
  $xmaxg = $a1x2/$normgx;
  $yming = $a1y1/$normgy;
  $ymaxg = $a1y2/$normgy;
  $R_Lcurve  = "-R$xming/$xmaxg/$yming/$ymaxg";
  $J_Lcurve  = "-JX$Jwid";

  $dx = $xmaxg - $xming;
  $dy = $ymaxg - $yming;
  $xtick = 0.25*$dx;
  $ytick = 0.25*$dy;

  $B_Lcurve  = sprintf("-B%.2f:\" log10 ( Model Norm ) \":/%.2f:\" log10 ( Misfit Norm )\":",$xtick,$ytick);
  #$B_Lcurve  = sprintf("-B%.2f:\" log10 ( Misfit Norm ) \":/%.2f:\" log10 ( Model Norm )\":",$xtick,$ytick);

  # polynomial lines and curves
  $Lcurve_file = "${recover_dir}/lcurve_pts";
  if (not -f $Lcurve_file)   { die("Check if $Lcurve_file exist or not\n") }
  open(IN,"$Lcurve_file");
     @Lcurve_pts = <IN>;

  #------------------

  for ($k = 1; $k <= $ndamp; $k = $k+1) {
    $stid = sprintf("%2.2i",$k);
    $file1syn = "${recover_dir}/m_${model_lab}_id${stid}_plot";
    $mod_files[$k-1] = $file1syn;
    if (not -f $file1syn)   { die("Check if $file1syn exist or not\n") }
  }

  # get the receivers and sources
  $recfile = "$dir0/${edir}/sr.txt";  # src-rec for first event
  $evefile = "$dir0/events_lonlat.dat";
  if (not -f $recfile) { die("Check if $recfile exist or not\n") }
  if (not -f $evefile) { die("Check if $evefile exist or not\n") }

  # file names for figures
  $name    = "${smap}_hess01_${strun0}_ilcurve${stil}_paper";
  $psfile  = "$name.ps";
  $jpgfile = "$name.jpg";

  $B = "$B0".$Bopts[15];

  #=============================================
  # recovered phase speed map

  #@kvec = (0,1,2,3,4,5,5);  # pick a subset of the possible figures
  @kvec = (0,8,12,14,16,20,20);
  $numk = @kvec;
  @labs = ("a","b","c","d","e","f","g","h","i","j");

  for ($i = 0; $i < $numk; $i = $i+1) {

    $k = $kvec[$i];
    $stgam = sprintf("\@~\147\@~ = %.2f",$lams[$k]);
    $title = "($labs[$i])  Model for $stgam";
    if ($i % 3 == 0) {$shift = $shift2} else {$shift = $shift1}

    if ($k==0) {
       print CSH "psbasemap $B $R $J -P -K -V $origin > $psfile\n"; # START
    } else {
       print CSH "psbasemap $B $R $J -P -K -O -V $shift >> $psfile\n";
    }
    # note 3rd column, not 4th
    if ($icolor==1) {print CSH "awk '{print \$1,\$2,\$3*100}' $mod_files[$k] | pscontour $R $J -A- -C$cpt_vel -I -P -O -K -V >> $psfile\n"}

    if ($icolor==1 && $i==$numk-1)
    {print CSH "awk '{print \$1,\$2,\$3*100}' $mod_files[$k] | pscontour $R $J -A- -C$cpt_vel_uniform -I -P -O -K -V >> $psfile\n"}

    print CSH "pscoast $J $R $coast_info -P -K -O -V >> $psfile\n";
    if($ishelf) {print CSH "awk '{print \$2,\$1}' $shelf_file |psxy $J $R $Wshelf -K -O -P -V >> $psfile\n"}
    print CSH "awk '{print \$1,\$2}' $evefile |psxy $J $R -K -O -P -V $src >> $psfile\n";
    print CSH "awk '\$1 == \"R\" {print \$2,\$3}' $recfile |psxy $J $R -K -O -P -V $rec >> $psfile\n";

    print CSH "pstext -N $J_title $R_title -K -O -V -P >>$psfile<<EOF\n $x_title $z_title $fsize_title 0 $fontno CM $title \nEOF\n";

    print CSH "psbasemap $B $R $J -P -K -O -V >> $psfile\n";

    # colorbar for exaggerted 'gradient' map
    if($i == $numk-1) {
      $bs1 = 0.4;
      $Bscaleg = sprintf("-B%2.2ef1:\" \": -E7p",$bs1);
      print CSH "psscale -C$cpt_vel_uniform $Dscale $Bscaleg -P -K -O -V >> $psfile \n";
      print CSH "pstext -N $J $R -K -O -V -P >>$psfile<<EOF\n-115.2 30.95 11 0 $fontno LM \@~\045\@~\nEOF\n";
    }
  }

  #=============================================
  # Lcurve

  $B = "$B_Lcurve".$Bopts[8];  # 15 for no axes
  $title = "($labs[$i])  L-curve for the models";
  $shift = $shift1;

  print CSH "psbasemap $B $R_Lcurve $J_Lcurve -P -K -O -V $shift >> $psfile\n";
  print CSH "psxy $Lcurve_file $R_Lcurve $J_Lcurve -W0.75p/0/0/255 -P -K -O -V >> $psfile\n";
  print CSH "psxy $Lcurve_file $R_Lcurve $J_Lcurve $p_info_w -P -K -O -V >> $psfile\n";

  # print a couple labels
  #@kvec = (0,1,2,3,4,5);   # from lowest to highest damping values
  #$numk = @kvec;

  @xshifts = (-0.035,-0.035,-0.06,-0.11,-0.11,0.04);
  @yshifts = (-0.06,-0.06,-0.05,0.02,0.02,0.04);
  @stgams  = ("(a)","(b)","(c)","(d)","(e)",sprintf("\@~\147\@~ = %.2f  (f, g)",$lams[$kvec[$numk-1]]));

  for ($j = 0; $j < $numk; $j = $j+1) {
    $k = $kvec[$j];

    # plot gamma = xx next to the point
    #if($j==5) {$stgam = sprintf("\@~\147\@~ = %.2f",$lams[$k])}
    #else {$stgam = sprintf("%.2f",$lams[$k])}
    $stgam = $stgams[$j];

    ($xpt0,$ypt0) = split(" ",$Lcurve_pts[$k]);
    $xpt = $xpt0 + $xshifts[$j]*$dx;
    $ypt = $ypt0 + $yshifts[$j]*$dy;
    print CSH "pstext -N $J_Lcurve $R_Lcurve -K -O -V -P >>$psfile<<EOF\n $xpt $ypt $fsize_title 0 $fontno LM $stgam \nEOF\n";
    #print CSH "psxy -N $J_Lcurve $R_Lcurve $p_info_w -K -O -V -P >>$psfile<<EOF\n $xpt0 $ypt0\nEOF\n";
  }
  print CSH "pstext -N $J_title $R_title -K -O -V -P >>$psfile<<EOF\n $x_title $z_title $fsize_title 0 $fontno CM $title \nEOF\n";

  #=============================================
  # phase speed for data

  $B = "$B0".$Bopts[15];
  $title = "($labs[$i+1])  Target phase speed model (data)";
  $shift = $shift1;

  print CSH "psbasemap $B $R $J -P -K -O -V $shift >> $psfile\n";
  if($icolor==1) {print CSH "awk '{print \$1,\$2,\$4*100}' $file1dat | pscontour $R $J -A- -C$cpt_vel -I -P -K -O -V  >> $psfile\n"}
  print CSH "pscoast $J $R $coast_info -P -K -O -V >> $psfile\n";
  if($ishelf) {print CSH "awk '{print \$2,\$1}' $shelf_file |psxy $J $R $Wshelf -K -O -P -V >> $psfile\n"}
  print CSH "awk '\$1 == \"R\" {print \$2,\$3}' $recfile |psxy $J $R -K -O -P -V $rec >> $psfile\n";
  print CSH "awk '{print \$1,\$2}' $evefile |psxy $J $R -K -O -P -V $src >> $psfile\n";
  print CSH "psscale -C$cpt_vel $Dscale $Bscale1b -P -K -O -V >> $psfile \n";
  print CSH "pstext -N $J $R -K -O -V -P >>$psfile<<EOF\n-115.2 30.95 11 0 $fontno LM \@~\045\@~\nEOF\n";
  print CSH "pstext -N $J_title $R_title -K -O -V -P >>$psfile<<EOF\n $x_title $z_title $fsize_title 0 $fontno CM $title \nEOF\n";
  print CSH "psbasemap $B $R $J -P -K -O -V >> $psfile\n";

#-----------------------------
  print CSH "pstext -N $J_title $R_title -O -V -P >>$psfile<<EOF\n 10 10 $fsize_title 0 $fontno CM $title \nEOF\n"; # FINISH
  print CSH "convert $psfile $jpgfile\n";
  print CSH "echo done with $psfile\n";

  close (CSH);
  system("csh -f $cshfile");
  system("xv $jpgfile &");
}

#=================================================================
