#!/usr/bin/perl -w

#==========================================================
#
#  joint_subspace.pl
#  Carl Tape
#  22-Nov-2008
#
#  This script plots figures for the joint source-structure inversions
#  produced by the 2D SEM wave propagation code.
#
#  EXAMPLE:
#    ../joint_subspace.pl -6/3.0/0/80/0.08 1 700 600 1  2 0  # structure (ifig2)
#
#  PREVIOUS EXAMPLES:
#    ../joint_subspace.pl -6/3.0/0/80/0.08 1 8350 8350 1  2 0  # source (ifig1)
#    ../joint_subspace.pl -6/3.0/0/80/0.08 1 8450 8450 1  4 0  # structure (ifig2)
#    ../joint_subspace.pl -6/3.0/0/80/0.08 1 8550 8550 1  1 0  # joint (ifig3, ifig0)
#
#==========================================================

if (@ARGV < 7) {die("Usage: plot_kernels.pl colors iker irun0 irun0_cg ipoly qmax HESSIAN\n");}
($colors,$iker,$irun0,$irun0_cg,$ipoly,$qmax,$ihessian) = @ARGV;
$iter = 0;

$sth = sprintf("%1i",$ihessian);

# USER change
$dirq = "/home/carltape/wave2d/SEM2D_iterate";

# directories
$pdir = "${dirq}/PLOTTING";
$ddir = "${pdir}/DATA_FILES";
$idir = "${dirq}/INPUT";
$odir = "${dirq}/OUTPUT/run_";
if (not -e $ddir) {die("Check if ddir $ddir exist or not\n");}
$edir = "event_005";  # event for first event (sr.txt)

# misfit function variable
$misfitvar = "\@~\143\@~";
$misfitvar = "S";

# resolution of color plots
$interp = "-I0.5m/0.5m -S4m";   # key information
$grdfile = "temp.grd";

$mfile_dat_str = "structure_dat.dat";
$mfile_syn_str = "structure_syn.dat";
$mfile_dat_src = "src_dat.dat";
$mfile_syn_src = "src_syn.dat";

$cshfile = "joint_subspace.csh";

#if($istructure==1){$Nfac=3}
#if($istructure==2){$Nfac=2}
#if($istructure==3){$Nfac=1}
#$sNfac = sprintf("%1i",$Nfac);

# boolean commands for plotting
$icolor = 1;   # ccc

$irun = $irun0_cg + $iter;
@mods = ("0","0t","1","1t","2","2t","3","3t","4","4t","5","5t","6","6t","7","7t","8","8t","9","9t","10","10t","11","11t","12","12t","13","13t","14","14t","15","15t","16","16t");
$mod = $mods[$iter];
$smod = "m\@+$mod\@+";

# label for the type of map used to generate the data
#$smap = sprintf("m%2.2i",$istructure);

# wave2d run number
$strun0_cg = sprintf("%4.4i",$irun0_cg);
$strun0 = sprintf("%4.4i",$irun0);
$strun = sprintf("%4.4i",$irun);
$dir0 = "${odir}${strun0}";
$dir = "${odir}${strun}";

#$stgam  = sprintf("\@~\147\@~ = %.1f km",$gamma/1000);
#$stgam2 = sprintf("%3.3i",$gamma/1000);

print "\n $dir, $edir, $mfile_dat_str, $mfile_syn_str";
print "\n $colors, $iker, $mod\n";
#die("testing");

# colors for the kernel and misfit function
@cols = split("/",$colors);
$kpwr = $cols[0];  # power for kernels
$kmax = $cols[1];  # value for kernels
$opwr = $cols[2];  # power for misfit maps
$omax = $cols[3];  # value for misfit maps
$cmax = $cols[4];  # value for phase velocity maps (percent pert)

#@files = glob("$dir/$kerfile");
#$numk = @files;

@titles = ("Waveform","Traveltime (xcorr), misfit","Amplitude (xcorr), misfit","Traveltime (MT), misfit","Amplitude (MT), misfit","Traveltime (xcorr), sampling","Amplitude (xcorr), sampling");
@units = ("m\@+2\@+ s","s\@+2\@+","xxx","xxx","xxx","xxx","xxx");
$ktype = $titles[$iker];
$utype = $units[$iker];

$plabel = "/home/carltape/wave2d/2d_adjoint/scripts/joint_subspace.pl";

# data files
#$dir_gmt = "/home/carltape/gmt";
$ishelf = 0;
$shelf_file = "$idir/BOUNDARIES/oms_shelf";
$plate_file = "$idir/BOUNDARIES/plate_boundaries";
$fault_file = "$idir/BOUNDARIES/faults/jennings.xy";

# plotting specifications
$fsize0 = "14";
$fsize1 = "11";
$fsize2 = "9";
$fsize3 = "5";
$fontno = "1";    # 1 or 4
$tick   = "0.15c";
$fpen   = "1.5p";
$tpen   = "1.0p";

# KEY: parameters for the source error STARS
$ref_rad = 6;
$fac = 3;

# plot symbols for sources, receivers, and shelf
$src    = "-W0.5p -Sa${ref_rad}p";
$rec    = "-Sc2p -G0";
$rec0   = "-Sc10p -W0.5p";
$src_ev = "-W0.5p -Sa12p -G255";
$Wshelf = "-W0.5p,0/0/0,--";

$dot_info = "-W0.5p -Sap";

$coast_info = "-W1p -Na/1p -Df";
$coast_info2 = "-W0.5p -Na/0.5p -Df";

$plate_info_w = "-M -W2.5p,255/255/255";
$plate_info_r = "-M -W2.5p,255/0/0";
$plate_info_b = "-M -W2.5p,0/0/255";
$plate_info_k = "-M -W2.5p";

$fault_info_r = "-M -W1.5p,255/0/0";
$fault_info_k = "-M -W1.5p,0/0/0";
$fault_info_w = "-M -W1.5p,255/255/255";

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

open(IN,"$dir0/reference_values.dat"); @lines = <IN>; ($alpha0,$beta0) = split(" ",$lines[0]);
open(IN,"$dir0/reference_period.dat"); @lines = <IN>; ($per) = split(" ",$lines[0]);
$lam = $beta0*$per;
$per_lab = sprintf("T = %3.3f s",$per);
$beta0_lab  = sprintf("beta0 = %3.3f km/s",$beta0/1000);
$lam_lab  = sprintf("\@~\154\@~ = %.0f km",$lam/1000);
print "\n$per_lab, $beta0_lab, $lam_lab\n";
#-------------------------

# write plotting scripts
$Jwid = 2;
$J = "-JM${Jwid}i";      # in lat-lon
$origin = "-X0.8 -Y8.25";

# which borders to plot the lat-lon
# 1 four sides, 4 single sides, 6 two sides, 4 three sides, 1 zero sides
@Bopts = ("WESN","Wesn","wesN","wEsn","weSn","WesN","wEsN","wESn","WeSn","WEsn","weSN","wESN","WESn","WeSN","WEsN","wesn");
$B0 = "-B1:.\" \":";

# axes scale for phase velocity maps: c(th, ph)
#$bs1 = 0.5;
#$Bscale1d  = sprintf("-B%2.2e:\" Phase Velocity for data ( km s\@+-1\@+ )\":",$bs1);
#$Bscale1s  = sprintf("-B%2.2e:\" Phase Velocity for model $smod ( km s\@+-1\@+ )\":",$bs1);
$bs1 = $cmax;
$bs2 = $bs1/4;
#$Bscale1  = sprintf("-B%2.2ef1:\" \@~\045\@~ pert. from %2.2f km/s\": -E5p",$bs1,$beta0/1000);
$Bscale1  = sprintf("-Ba%2.2ef%2.2e:\" ln(c / c0), c0 = %2.2f km/s\": -E5p",$bs1,$bs2,$beta0/1000);
$Bscale1b = sprintf("-Ba%2.2ef%2.2e:\" \": -E5p",$bs1,$bs2);

# axes scale for kernels: K(th, ph)
# [$misfitvar] --> s
$tp = "\@~\146\@~, \@~\161\@~";
$Bscale2  = sprintf("-B%2.2e:\" K ( $tp )  ( 10\@+%2.2i\@+  m\@+-2\@+ s )\": -E5p",$bs2,$kpwr);
$Bscale2b = sprintf("-B%2.2e:\" \": -E5p",$bs2);

# axes scale for chi_plots: chi(th_r, ph_r)
#$Bscale3  = sprintf("-B%2.2e:\" $misfitvar ( \@~\161\@~\@-r\@- , \@~\146\@~\@-r\@- )  ( 10\@+%2.2i\@+ )\": -Ef5p",$bs3,$opwr);

#-------------------------

# model files from the final directory
$file1dat_str = "${dir}/${mfile_dat_str}";
$file1syn_str = "${dir}/${mfile_syn_str}";
$file1dat_src = "${dir}/${mfile_dat_src}";
$file1syn_src = "${dir}/${mfile_syn_src}";

# set bounds for the plotting
#$xmin = -121; $xmax = -115.0; $zmin = 31.75; $zmax = 36.5;
($xmin,$xmax,$zmin,$zmax) = split(" ",`minmax -C $file1dat_str`);
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
$c_info_rd  = "-W1.0p,255/0/0,--";
$c_info_bd  = "-W1.0p,0/0/255,--";
$c_info_kd  = "-W0.5p,0/0/0,--";

if ($ipoly==1) {

  # axes limits for polynomial plots
  $axes_file = "${ddir}/axes_${strun0_cg}.dat";
  if (not -f $axes_file) {die("Check if axes_file $axes_file exist or not\n");}
  open(IN,"$axes_file");
  @ax_lims = <IN>;
  ($a1x1,$a1x2,$a1y1,$a1y2) = split(" ",$ax_lims[0]); # polynomial
  ($a2x1,$a2x2,$a2y1,$a2y2) = split(" ",$ax_lims[1]); # chi vs iteration
  ($a3x1,$a3x2,$a3y1,$a3y2) = split(" ",$ax_lims[2]); # var-red vs iteration

  # chi-vs-iteration
  $xming = $a2x1; $xmaxg = $a2x2; $yming = $a2y1; $ymaxg = $a2y2;
  $chi_ran = $xmaxg - $xming;
  $R_chi = "-R$xming/$xmaxg/$yming/$ymaxg";
  $J_chi = "-JX${Jwid}i/${Jwid}il"; # note log scale on y-axis

  # text labels for chi-vs-m plots
  $schi0 = "$misfitvar ( m\@+0\@+ )";
  $schi1 = "$misfitvar ( m\@+1\@+ )";
  $x0chit = 0 + 0.05*$chi_ran;
  $x1chit = 1 + 0.05*$chi_ran;

  $chi_curve = "${ddir}/chi_curve_${strun0_cg}.dat";
  if (not -f $chi_curve) {die("Check if $chi_curve exist or not\n");}

  # scale for chi-vs-m plots (a1f2g1p)
  $iter_tick = 2;
  $B3a    = "-B${iter_tick}f1g2:\" k, model number \":/a1f3g3p:\" $misfitvar ( m )   ( $utype ) \":";
  $B3b    = "-B${iter_tick}f1g2:\" k, model number \":/a1f3g3p:\" \":";
}

#===========================================================================
# create colorpoint files

open(CSH,">$cshfile");

# phase velocity maps
$dc = $cmax/10;
$T1 = "-T-$cmax/$cmax/$dc";
$cpt_vel = "color0.cpt";
print CSH "makecpt -C$colorbar $T1 -D > $cpt_vel\n";

# origin time
$cmax = 1;
$dc = $cmax/40;
$T2 = "-T-$cmax/$cmax/$dc";
$cptfile = "color_src.cpt";
print CSH "makecpt -Cpolar $T2 -I -D > $cptfile\n";

close (CSH);
system("csh -f $cshfile");
#die("testing");

#===========================================================================

# shifting the subplots
$xfac = 1.20;
$yfac = 1.20;
$dX1p = $xfac*$Jwid; $dX1m = -$xfac*$Jwid;
$dY1p = $yfac*$Jwid; $dY1m = -$yfac*$Jwid;
$dX2p = 2*$xfac*$Jwid; $dX2m = -2*$xfac*$Jwid;
$dY2p = 2*$yfac*$Jwid; $dY2m = -2*$yfac*$Jwid;

# shifting the subplots
$yfac = 1.50;
$dY1p = $yfac*$Jwid; $dY1m = -$yfac*$Jwid;
$dY2p = 2*$yfac*$Jwid; $dY2m = -2*$yfac*$Jwid;
$yfac = 1.20; $dY1mB = -$yfac*$Jwid;

# colorbar
$Dlen = 1.2; $Dx = $Jwid/2; $Dy = -0.10*$Jwid;
$Dscale = "-D$Dx/$Dy/$Dlen/0.10h";

#===============================================
print "\nWriting CSH file...\n";

open(CSH,">$cshfile");
print CSH "gmtset BASEMAP_TYPE plain PAPER_MEDIA letter TICK_LENGTH $tick LABEL_FONT_SIZE $fsize2 ANOT_FONT_SIZE $fsize2 PLOT_DEGREE_FORMAT D HEADER_FONT $fontno ANOT_FONT $fontno LABEL_FONT $fontno HEADER_FONT_SIZE $fsize1 FRAME_PEN $fpen TICK_PEN $tpen\n";
#===============================================

# iteration index (subset of all possible)
@kvec = (0,$qmax);
$numk = @kvec;
$niter_max = $kvec[$numk-1];  # max iteration is the last in the list

# load all possible files
for ($k = 0; $k <= $niter_max; $k = $k+1) {

  $irun = $irun0_cg + 2*$k; # wave2d run number
  $strun = sprintf("%4.4i",$irun);
  $dir = "${odir}${strun}";

  # if the simulation stopped on a test model, use it instead
  if ( ($k == $niter_max) && (not -e $dir) ) {
     $irun = $irun - 1;
     $strun = sprintf("%4.4i",$irun);
     $dir = "${odir}${strun}";
   }

  $file1syn_str = "${dir}/$mfile_syn_str";
  $file1syn_src = "${dir}/$mfile_syn_src";
  if (not -f $file1syn_str) {die("Check if file1syn_str $file1syn_str exist or not\n");}
  if (not -f $file1syn_src) {die("Check if file1syn_src $file1syn_src exist or not\n");}
  $str_files[$k] = $file1syn_str;
  $src_files[$k] = $file1syn_src;

  # load chi values
  $chifile = "${dir}/chi.dat";
  if (not -f $chifile) {die("Check if chifile $chifile exist or not\n");}
  open(IN,"$chifile");
  $chi = <IN>; chomp($chi);
  $it_vals[$k] = $k;
  $chi_vals[$k] = $chi;
  #$schi = sprintf("$misfitvar ( $smod )  =  %3.3e",$chi);
}

# load files from Hessian run
$it_vals_hess[0] = $it_vals[0];
$chi_vals_hess[0] = $chi_vals[0];
$str_files_hess[0] = $str_files[0];
$src_files_hess[0] = $src_files[0];

for ($h = 1; $h <= $qmax; $h = $h+1) {
  $dir = sprintf("$dir0/HESSIAN/model_m%2.2i",$h);
  $file1syn_str = "${dir}/$mfile_syn_str";
  $file1syn_src = "${dir}/$mfile_syn_src";
  if (not -f $file1syn_str) {die("Check if file1syn_str $file1syn_str exist or not\n")}
  if (not -f $file1syn_src) {die("Check if file1syn_src $file1syn_src exist or not\n")}
  $str_files_hess[$h] = $file1syn_str;
  $src_files_hess[$h] = $file1syn_src;

  # load chi values
  $chifile = "${dir}/chi.dat";
  if (not -f $chifile) {die("Check if chifile $chifile exist or not\n");}
  open(IN,"$chifile");
  $chi = <IN>; chomp($chi);
  $it_vals_hess[$h] = $h;
  $chi_vals_hess[$h] = $chi;
}

print "\n @str_files \n @src_files \n @str_files_hess \n @src_files_hess \n";

# get the receivers and sources
$recfile = "$dir0/${edir}/sr.txt"; # src-rec for first event
$evefile = "$dir0/events_dat_lonlat.dat"; # sources for DATA
if (not -f $recfile) {die("Check if $recfile exist or not\n");}
if (not -f $evefile) {die("Check if $evefile exist or not\n");}

@labs  = ("a","b","c","d","e","f","g","h","i","j");

$ifig0 = 0;    # initial and target models
$ifig1 = 0;    # source inversion
$ifig2 = 1;    # structure inversion

#=============================================
# INITIAL AND TARGET MODELS
#=============================================

if($ifig0 == 1) {

# file names for figures
$name    = "joint_subspace_00"; $psfile  = "$name.ps"; $jpgfile = "$name.jpg";

@bvec = (15,9,15);

# target structure
$B = "$B0".$Bopts[9];
$title = "(a)  Target structure";
$str_file = $file1dat_str;
print CSH "psbasemap $B $R $J -K -V -P $origin > $psfile\n";  # START
if ($icolor==1) {
  print CSH "awk '{print \$1,\$2,\$7}' $str_file | nearneighbor -G$grdfile $R $interp\n";
  print CSH "grdimage $grdfile -C${cpt_vel} $J -K -O -V -Q >> $psfile\n";
}
print CSH "pscoast $J $R $coast_info -K -O -V >> $psfile\n";
print CSH "awk '{print \$1,\$2}' $evefile |psxy -N $J $R -K -O -V $src >> $psfile\n";
print CSH "awk '\$1 == \"R\" {print \$2,\$3}' $recfile |psxy -N $J $R -K -O -V $rec >> $psfile\n";
print CSH "pstext -N $J_title $R_title -K -O -V >>$psfile<<EOF\n $x_title $z_title $fsize_title 0 $fontno CM $title \nEOF\n";

# target source
$B = "$B0".$Bopts[15];
$title = "(b)  Target source";
$shift = "-X$dX1p";
$src_file = $file1dat_src;
$dots_file = "temp";
print CSH "awk '{print \$1,\$2,\$8,${ref_rad} + sqrt(\$6*\$6 + \$7*\$7)*$fac/1000}' $src_file > $dots_file\n";
print CSH "psbasemap $B $R $J -K -O -V $shift >> $psfile\n";
print CSH "pscoast $J $R $coast_info2 -K -O -V >> $psfile\n";
print CSH "awk '\$1 == \"R\" {print \$2,\$3}' $recfile |psxy -N $J $R -K -O -V $rec >> $psfile\n";
print CSH "psxy $dots_file $J $R $dot_info -C$cptfile -K -O -V >> $psfile \n";
print CSH "pstext -N $J_title $R_title -K -O -V >>$psfile<<EOF\n $x_title $z_title $fsize_title 0 $fontno CM $title \nEOF\n";

# initial structure
$B = "$B0".$Bopts[9];
$title = "(c)  Initial structure";
$shift = "-X$dX1m -Y$dY1mB";
$str_file = $str_files[0];
print CSH "psbasemap $B $R $J -K -V -O $shift >> $psfile\n";
if ($icolor==1) {
  print CSH "awk '{print \$1,\$2,\$7}' $str_file | nearneighbor -G$grdfile $R $interp\n";
  print CSH "grdimage $grdfile -C${cpt_vel} $J -K -O -V -Q >> $psfile\n";
}
print CSH "pscoast $J $R $coast_info -K -O -V >> $psfile\n";
print CSH "awk '{print \$1,\$2}' $evefile |psxy -N $J $R -K -O -V $src >> $psfile\n";
print CSH "awk '\$1 == \"R\" {print \$2,\$3}' $recfile |psxy -N $J $R -K -O -V $rec >> $psfile\n";
print CSH "pstext -N $J_title $R_title -K -O -V >>$psfile<<EOF\n $x_title $z_title $fsize_title 0 $fontno CM $title \nEOF\n";
print CSH "psscale -C$cpt_vel $Dscale $Bscale1 -K -O -V >> $psfile \n";

# initial source
$B = "$B0".$Bopts[15];
$shift = "-X$dX1p";
$title = "(d)  Initial source";
$src_file = $src_files[0];
$dots_file = "temp";
print CSH "awk '{print \$1,\$2,\$8,${ref_rad} + sqrt(\$6*\$6 + \$7*\$7)*$fac/1000}' $src_file > $dots_file\n";
print CSH "psbasemap $B $R $J -K -O -V $shift >> $psfile\n";
print CSH "pscoast $J $R $coast_info2 -K -O -V >> $psfile\n";
print CSH "awk '\$1 == \"R\" {print \$2,\$3}' $recfile |psxy -N $J $R -K -O -V $rec >> $psfile\n";
print CSH "psxy $dots_file $J $R $dot_info -C$cptfile -K -O -V >> $psfile \n";
print CSH "pstext -N $J_title $R_title -K -O -V >>$psfile<<EOF\n $x_title $z_title $fsize_title 0 $fontno CM $title \nEOF\n";

# legend for sources
$Dx = 0.2; $Dy = -0.2;
$Dscale_otime = "-D$Dx/$Dy/0.8/0.10h";;
$Bscale_otime = "-B1f0.25:\"Origin time error (s)\": -E5p";
print CSH "psscale -C$cptfile $Dscale_otime $Bscale_otime -K -O -V >> $psfile\n";

# convert source errors to dot size
@misloc_dots = (0,2,4);
$ndot = @misloc_dots;
$xlon0 = -117; $dlon = 1.0;
for ($j = 0; $j < $ndot; $j = $j+1) {
  $misloc = $misloc_dots[$j];
  $misloc_size[$j] = $ref_rad + $misloc*$fac; # KEY: use same formula as above
  $xlon[$j] = $xlon0 + $j*$dlon;
}

# source error scale -- symbols
$origin_box = "-Xa0 -Ya0"; # location for mislocation legend
$yp0 = $zmin - 2.0;
$yp1 = $yp0 + 0.5;
$yp2 = $yp1 + 0.5;
$yp3 = $yp2 + 0.5;
print CSH "psxy $J $R $dot_info -N -K -V -O $origin_box >> $psfile <<EOF
  $xlon[2] $yp3 $misloc_size[2]
  $xlon[1] $yp3 $misloc_size[1]
  $xlon[0] $yp3 $misloc_size[0]
EOF\n";

# source error scale -- text
print CSH "pstext $J $R -N -K -V -O $origin_box >> $psfile <<EOF
  $xlon[2] $yp2 $fsize2 0 $fontno CM 4
  $xlon[1] $yp2 $fsize2 0 $fontno CM 2
  $xlon[0] $yp2 $fsize2 0 $fontno CM 0
  $xlon[1] $yp1 $fsize2 0 $fontno CM Source
  $xlon[1] $yp0 $fsize2 0 $fontno CM mislocation (km)
EOF\n";

#-----------------------------
print CSH "pstext -N $J_title $R_title -O -V >>$psfile<<EOF\n 10 10 $fsize_title 0 $fontno CM $title \nEOF\n"; # FINISH
print CSH "convert $psfile $jpgfile\n";
print CSH "echo done with $psfile\n";
print CSH "ghostview $psfile &\n";

}

#=============================================
# SOURCE INVERSION
#=============================================

if($ifig1 == 1) {

# file names for figures
$name    = "joint_subspace_01"; $psfile  = "$name.ps"; $jpgfile = "$name.jpg";

@tlabs = ("CG: Error in initial source","CG: Error in source","CG: Error in final source");
@kvec = (0,1,2); $numk = @kvec;
@bvec = (15,9,15);
$t = 0;
for ($i = 0; $i < $numk; $i = $i+1) {

  $k = $kvec[$i];
  $t = $t+1;
  $b = $bvec[$i];
  $B = "$B0".$Bopts[$b];

  $mod = $mods[2*$k];
  $smod = "m${mod}";
  #$smod = "m\@+${mod}\@+";
  $title = "($labs[$t-1])  $tlabs[$i] ${smod}";

  $source_error_file = $src_files[$k];
  $dots_file = "temp";
  print CSH "awk '{print \$1,\$2,\$8,${ref_rad} + sqrt(\$6*\$6 + \$7*\$7)*$fac/1000}' $source_error_file > $dots_file\n";
  print CSH "awk '{print \$1,\$2,\$8,sqrt(\$6*\$6 + \$7*\$7)/1000}' $source_error_file > source_error\n";
  #print "\n $source_error_file \n $dots_file \n"; die("testing");
  #-------------------

  if ($i==0) {
     $shift = "-X$dX2m -Y$dY1mB";
  } elsif ($i==1 || $i==2) {
     $shift = "-X$dX1p";
  }

  if($i==0) {
     print CSH "psbasemap $B $R $J -K -V -P $origin > $psfile\n";  # START
  } else {
     print CSH "psbasemap $B $R $J -K -O -V $shift >> $psfile\n";
  }
  print CSH "pscoast $J $R $coast_info2 -K -O -V >> $psfile\n";

  print CSH "awk '\$1 == \"R\" {print \$2,\$3}' $recfile |psxy -N $J $R -K -O -V $rec >> $psfile\n";

  # plot source errors as colored and sized stars
  $dot_info = "-W0.5p -Sap";
  print CSH "psxy $dots_file $J $R $dot_info -C$cptfile -K -O -V >> $psfile \n";

  print CSH "pstext -N $J_title $R_title -K -O -V >>$psfile<<EOF\n $x_title $z_title $fsize_title 0 $fontno CM $title \nEOF\n";
}

#----------------------------

@tlabs = ("Misfit,  S (mk)","Subspace: Error in source","Subspace: Error in final source");
@kvec = (0,1,2); $numk = @kvec;
@bvec = (8,1,8);
for ($i = 0; $i < $numk; $i = $i+1) {

  $k = $kvec[$i];
  $t = $t+1;
  $b = $bvec[$i];
  $B = "$B0".$Bopts[$b];

  $mod = $mods[2*$k];
  $smod = "m${mod}";
  #$smod = "m\@+${mod}\@+";
  $title = "($labs[$t-1])  $tlabs[$i] ${smod}";
  $source_error_file = $src_files_hess[$k];

  $dots_file = "temp";
  print CSH "awk '{print \$1,\$2,\$8,${ref_rad} + sqrt(\$6*\$6 + \$7*\$7)*$fac/1000}' $source_error_file > $dots_file\n";
  print CSH "awk '{print \$1,\$2,\$8,sqrt(\$6*\$6 + \$7*\$7)/1000}' $source_error_file > source_error\n";
  #print "\n $source_error_file \n $dots_file \n"; die("testing");
  #-------------------

  if ($i==0) {
    $shift = "-X$dX2m -Y$dY1mB";
  } elsif ($i==1 || $i==2) {
    $shift = "-X$dX1p";
  }

  if ($i == 0) {

    # plot the misfit function
    $Bg = "$B3b".$Bopts[$b];
    $title = "($labs[$t-1])  $tlabs[$i]";
    print CSH "psbasemap $Bg $R_chi $J_chi -K -O -V $shift >> $psfile\n";
    #print CSH "awk '{print \$1,\$2}' $chi_curve | psxy $c_info_rd $J_chi $R_chi -K -O -V >> $psfile\n";
    for ($p = 0; $p <= $niter_max; $p = $p+1) {
      print CSH "psxy $J_chi $R_chi $p_info_w -K -O -V >>$psfile<<EOF\n $it_vals[$p] $chi_vals[$p]\nEOF\n";
      print CSH "psxy $J_chi $R_chi $p_info_k -K -O -V >>$psfile<<EOF\n $it_vals_hess[$p] $chi_vals_hess[$p]\nEOF\n";
    }
    print CSH "pstext -N $J_title $R_title -K -O -V >>$psfile<<EOF\n $x_title $z_title $fsize_title 0 $fontno CM $title \nEOF\n";

  } else {

    # plot the source errors
    print CSH "psbasemap $B $R $J -K -O -V $shift >> $psfile\n";
    print CSH "pscoast $J $R $coast_info2 -K -O -V >> $psfile\n";
    print CSH "awk '\$1 == \"R\" {print \$2,\$3}' $recfile |psxy -N $J $R -K -O -V $rec >> $psfile\n";

    # plot source errors as colored and sized stars
    $dot_info = "-W0.5p -Sap";
    print CSH "psxy $dots_file $J $R $dot_info -C$cptfile -K -O -V >> $psfile \n";

    if ($i==1) {
      # depth scale
      $Dx = 0.3; $Dy = -0.2;

      $Dscale_otime = "-D$Dx/$Dy/0.8/0.10h";;
      $Bscale_otime = "-B1f0.25:\"Origin time error (s)\": -E5p";
      print CSH "psscale -C$cptfile $Dscale_otime $Bscale_otime -K -O -V >> $psfile\n";

      # convert source errors to dot size
      @misloc_dots = (0,2,4);
      $ndot = @misloc_dots;
      $xlon0 = -117; $dlon = 1.0;
      for ($j = 0; $j < $ndot; $j = $j+1) {
  $misloc = $misloc_dots[$j];
  $misloc_size[$j] = $ref_rad + $misloc*$fac; # KEY: use same formula as above
  $xlon[$j] = $xlon0 + $j*$dlon;
      }

      # source error scale -- symbols
      $origin_box = "-Xa0.3 -Ya-0.13"; # location for mislocation legend
      $yp0 = $zmin - 2.0;
      $yp1 = $yp0 + 0.5;
      $yp2 = $yp1 + 0.5;
      $yp3 = $yp2 + 0.5;
      print CSH "psxy $J $R $dot_info -N -K -V -O $origin_box >> $psfile <<EOF
  $xlon[2] $yp3 $misloc_size[2]
  $xlon[1] $yp3 $misloc_size[1]
  $xlon[0] $yp3 $misloc_size[0]
EOF\n";

      # source error scale -- text
      print CSH "pstext $J $R -N -K -V -O $origin_box >> $psfile <<EOF
  $xlon[2] $yp2 $fsize2 0 $fontno CM 4
  $xlon[1] $yp2 $fsize2 0 $fontno CM 2
  $xlon[0] $yp2 $fsize2 0 $fontno CM 0
  $xlon[1] $yp1 $fsize2 0 $fontno CM Source mislocation (km)
EOF\n";
    }
    print CSH "pstext -N $J_title $R_title -K -O -V >>$psfile<<EOF\n $x_title $z_title $fsize_title 0 $fontno CM $title \nEOF\n";
  }
}


#-----------------------------
print CSH "pstext -N $J_title $R_title -O -V >>$psfile<<EOF\n 10 10 $fsize_title 0 $fontno CM $title \nEOF\n"; # FINISH
print CSH "convert $psfile $jpgfile\n";
print CSH "echo done with $psfile\n";
print CSH "ghostview $psfile &\n";

}

#=============================================
#  STRUCTURE INVERSION
#=============================================

if($ifig2 == 1) {

# file names for figures
$name    = "joint_subspace_02"; $psfile  = "$name.ps"; $jpgfile = "$name.jpg";

@tlabs = ("CG: Structure","CG: Structure","CG: Final structure");
@kvec = (1,2,$qmax); $numk = @kvec;
@bvec = (15,9,15);
$t = 0;
for ($i = 0; $i < $numk; $i = $i+1) {

  $k = $kvec[$i];
  $t = $t+1;
  $b = $bvec[$i];
  $B = "$B0".$Bopts[$b];

  $mod = $mods[2*$k];
  $smod = "m${mod}";
  #$smod = "m\@+${mod}\@+";
  $title = "($labs[$t-1])  $tlabs[$i] ${smod}";
  $str_file = $str_files[$k];

  if ($i==0) {
     $shift = "-X$dX2m -Y$dY1mB";
  } elsif ($i==1 || $i==2) {
     $shift = "-X$dX1p";
  }

  if($i==0) {
     print CSH "psbasemap $B $R $J -K -V -P $origin > $psfile\n";  # START
  } else {
     print CSH "psbasemap $B $R $J -K -O -V $shift >> $psfile\n";
  }

 # plot the structure
  if ($icolor==1) {
    #print CSH "awk '{print \$1,\$2,\$6}' $str_file | pscontour $R $J -A- -C$cpt_vel -I -O -K -V >> $psfile\n";
    print CSH "awk '{print \$1,\$2,\$7}' $str_file | nearneighbor -G$grdfile $R $interp\n";
    print CSH "grdimage $grdfile -C${cpt_vel} $J -K -O -V -Q >> $psfile\n";
  }
  print CSH "pscoast $J $R $coast_info -K -O -V >> $psfile\n";
  print CSH "awk '{print \$1,\$2}' $evefile |psxy -N $J $R -K -O -V $src >> $psfile\n";
  print CSH "awk '\$1 == \"R\" {print \$2,\$3}' $recfile |psxy -N $J $R -K -O -V $rec >> $psfile\n";

  print CSH "pstext -N $J_title $R_title -K -O -V >>$psfile<<EOF\n $x_title $z_title $fsize_title 0 $fontno CM $title \nEOF\n";
}

#----------------------------

@tlabs = ("Subspace: Structure","Subspace: Structure","Subspace: Final structure");
@kvec = (1,2,$qmax); $numk = @kvec;
@bvec = (15,9,15);
for ($i = 0; $i < $numk; $i = $i+1) {

  $k = $kvec[$i];
  $t = $t+1;
  $b = $bvec[$i];
  $B = "$B0".$Bopts[$b];

  $mod = $mods[2*$k];
  $smod = "m${mod}";
  #$smod = "m\@+${mod}\@+";
  $title = "($labs[$t-1])  $tlabs[$i] ${smod}";
  $str_file = $str_files_hess[$k];

  if ($i==0) {
    $shift = "-X$dX2m -Y$dY1mB";
  } elsif ($i==1 || $i==2) {
    $shift = "-X$dX1p";
  }

  # plot the structure
  print CSH "psbasemap $B $R $J -K -O -V $shift >> $psfile\n";
  if ($icolor==1) {
    #print CSH "awk '{print \$1,\$2,\$6}' $str_file | pscontour $R $J -A- -C$cpt_vel -I -O -K -V >> $psfile\n";
    print CSH "awk '{print \$1,\$2,\$7}' $str_file | nearneighbor -G$grdfile $R $interp\n";
    print CSH "grdimage $grdfile -C${cpt_vel} $J -K -O -V -Q >> $psfile\n";
  }
  print CSH "pscoast $J $R $coast_info -K -O -V >> $psfile\n";
  print CSH "awk '{print \$1,\$2}' $evefile |psxy -N $J $R -K -O -V $src >> $psfile\n";
  print CSH "awk '\$1 == \"R\" {print \$2,\$3}' $recfile |psxy -N $J $R -K -O -V $rec >> $psfile\n";
  if ($i == 0) {
    print CSH "psscale -C$cpt_vel $Dscale $Bscale1 -K -O -V >> $psfile \n";
  }

  print CSH "pstext -N $J_title $R_title -K -O -V >>$psfile<<EOF\n $x_title $z_title $fsize_title 0 $fontno CM $title \nEOF\n";

}

# plot the target map
$B = "$B0".$Bopts[8];
$t = $t+1;
$shift = "-X$dX1m -Y$dY1mB";
$title = "($labs[$t-1])  Target structure";
print CSH "psbasemap $B $R $J -K -O -V $shift >> $psfile\n";
if ($icolor==1) {
  print CSH "awk '{print \$1,\$2,\$7}' $file1dat_str | nearneighbor -G$grdfile $R $interp\n";
  print CSH "grdimage $grdfile -C${cpt_vel} $J -K -O -V -Q >> $psfile\n";
}
print CSH "pscoast $J $R $coast_info -K -O -V >> $psfile\n";
print CSH "awk '{print \$1,\$2}' $evefile |psxy -N $J $R -K -O -V $src >> $psfile\n";
print CSH "awk '\$1 == \"R\" {print \$2,\$3}' $recfile |psxy -N $J $R -K -O -V $rec >> $psfile\n";
print CSH "pstext -N $J_title $R_title -K -O -V >>$psfile<<EOF\n $x_title $z_title $fsize_title 0 $fontno CM $title \nEOF\n";

# plot the misfit function
$Bg = "$B3b".$Bopts[8];
$t = $t+1;
$title = "($labs[$t-1])  Misfit,  S (mk)";
$shift = "-X$dX1p";
print CSH "psbasemap $Bg $R_chi $J_chi -K -O -V $shift >> $psfile\n";
#print CSH "awk '{print \$1,\$2}' $chi_curve | psxy $c_info_rd $J_chi $R_chi -K -O -V >> $psfile\n";
for ($p = 0; $p <= $niter_max; $p = $p+1) {
  print CSH "psxy $J_chi $R_chi $p_info_w -K -O -V >>$psfile<<EOF\n $it_vals[$p] $chi_vals[$p]\nEOF\n";
  print CSH "psxy $J_chi $R_chi $p_info_k -K -O -V >>$psfile<<EOF\n $it_vals_hess[$p] $chi_vals_hess[$p]\nEOF\n";
}
print CSH "pstext -N $J_title $R_title -K -O -V >>$psfile<<EOF\n $x_title $z_title $fsize_title 0 $fontno CM $title \nEOF\n";


#-----------------------------
print CSH "pstext -N $J_title $R_title -O -V >>$psfile<<EOF\n 10 10 $fsize_title 0 $fontno CM $title \nEOF\n"; # FINISH
print CSH "convert $psfile $jpgfile\n";
print CSH "echo done with $psfile\n";
print CSH "ghostview $psfile &\n";

}

#=============================================

close (CSH);
system("csh -f $cshfile");

#=================================================================
