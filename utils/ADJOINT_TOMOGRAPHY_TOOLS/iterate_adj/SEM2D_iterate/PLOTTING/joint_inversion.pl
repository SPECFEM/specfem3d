#!/usr/bin/perl -w

#==========================================================
#
#  joint_inversion.pl
#  Carl Tape
#  18-Nov-2009
#
#  This script plots figures for the joint source-structure inversions
#  produced by the 2D SEM wave propagation code.
#  It also accomodates the structure-only and source-only inversions.
#  The user must first have run the matlab scripts wave2d_cg_poly.m and wave2d_cg_figs.m.
#
#  EXAMPLES:
#    ../joint_inversion.pl -6/3.0/0/80/0.08 1 100 1 15 0    # 1 source, structure only
#    ../joint_inversion.pl -6/3.0/0/80/0.08 1 200 1  3 0    # 1 source, source only
#    ../joint_inversion.pl -6/3.0/0/80/0.08 1 300 1 10 0    # 1 source, joint
#    ../joint_inversion.pl -6/3.0/0/80/0.08 1 400 1 16 0    # 5 sources, joint
#    ../joint_inversion.pl -6/3.0/0/80/0.08 1 500 1 16 0    # 25 sources, joint
#    ../joint_inversion.pl -6/3.0/0/80/0.08 1 600 1 16 0    # 25 sources, structure only
#    ../joint_inversion.pl -6/3.0/0/80/0.08 1 700 1  2 1    # 25 sources, structure only, SUBSPACE
#
#    ../joint_inversion.pl -6/3.0/0/80/0.08 1 800 1 16 0    # 1 sources, structure
#    ../joint_inversion.pl -6/3.0/0/80/0.08 1 XXX 1 16 0    # 5 sources, joint
#    ../joint_inversion.pl -6/3.0/0/80/0.08 1 XXX 1 16 0    # 25 sources, joint
#
#==========================================================

if (@ARGV < 6) {die("Usage: plot_kernels.pl colors iker irun0 ipoly qmax READ_IN\n");}
($colors,$iker,$irun0,$ipoly,$qmax,$ihessian) = @ARGV;
$iter = 0;

# base directory
$pwd = $ENV{PWD};
$dirplot = `dirname $pwd`; chomp($dirplot);
$basedir = `dirname $dirplot`; chomp($basedir);
#$basedir = "/data1/cig/seismo/3D/ADJOINT_TOMO/iterate_adj/SEM2D_iterate_work";

#-----------------------------------
# USER INPUT

# OPTIONAL: directory with the "reference" conjugate gradient simulation is
$irun0_cg = $irun0;
#$irun0_cg = 600;

# resolution of color plots
#$interp = "-I0.5m/0.5m -S4m";
$interp = "-I2m/2m -S4m";

$edir = "event_005";  # event for first event (sr.txt)

#-----------------------------------

$sth = sprintf("%1i",$ihessian);

# directories
$pdir = "${basedir}/PLOTTING";
$ddir = "${pdir}/DATA_FILES";
$idir = "${basedir}/INPUT";
$odir = "${basedir}/OUTPUT/run_";
#$odir = "/data1/cig/seismo/3D/ADJOINT_TOMO/iterate_adj/SEM2D_iterate_OUTPUT_OLD_2/run_";  # CHT
if (not -e $ddir) {die("Check if ddir $ddir exist or not\n");}

# misfit function variable
$misfitvar = "\@~\143\@~";
$misfitvar = "S";

# grd file for interpolation
$grdfile = "temp.grd";

$mfile_dat_str = "structure_dat.dat";
$mfile_syn_str = "structure_syn.dat";
$mfile_dat_src = "src_dat.dat";
$mfile_syn_src = "src_syn.dat";

$cshfile = "joint_inversion.csh";

#if($istructure==1){$Nfac=3}
#if($istructure==2){$Nfac=2}
#if($istructure==3){$Nfac=1}
#$sNfac = sprintf("%1i",$Nfac);

# boolean commands for plotting
$icolor = 1;   # ccc

$irun = $irun0_cg + $iter;
@mods = ("0","0t","1","1t","2","2t","3","3t","4","4t","5","5t","6","6t","7","7t","8","8t","9","9t","10","10t","11","11t","12","12t","13","13t","14","14t","15","15t","16","16t","17","17t","18","18t","19","19t","20","20t","21","21t","22","22t","23","23t","24","24t","25","25t","26","26t","27","27t","28","28t","29","29t","30","30t","31","31t","32");
$mod = $mods[$iter];
$smod = "m\@+$mod\@+";

# label for the type of map used to generate the data
#$smap = sprintf("m%2.2i",$istructure);

# wave2d run number
$strun0_cg = sprintf("%4.4i",$irun0_cg);
$strun0 = sprintf("%4.4i",$irun0);
$strun = sprintf("%4.4i",$irun);
$dir0 = "$odir$strun0";
$dir = "$odir$strun";

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
#@units = ("m\@+2\@+ s","s\@+2\@+","xxx","xxx","xxx","xxx","xxx");
@units = ("  ","  ","xxx","xxx","xxx","xxx","xxx");
$ktype = $titles[$iker];
$utype = $units[$iker];

$plabel = "${pdir}/joint_inversion.pl";

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

##-------------------------
## color for chi (objective function)

#$norm3 = "1e$opwr";
#$ss = $omax;
#$ds = $ss/$scale_color;
#$bs3 = sprintf("%3.3e",0.9*$ss/3); # colorbar
#$T3 = sprintf("-T%3.3e/%3.3e/%3.3e",0,$ss,$ds);
#print "T3 = $T3\n";

#-------------------------
# get reference phase velocity and period

#$c0 = 3500; $per = 20;

#open(IN,"$dir0/socal_vel_c0.dat");
#@vals = <IN>;
#$c0       = $vals[0];
#$per      = $vals[1];
#$lam      = $c0*$per;
#$per_lab  = sprintf("T = %.1f s",$per);
#$c0_lab   =  sprintf("c0 = %.2f km/s",$c0/1000);
#$lam_lab  = sprintf("\@~\154\@~ = %.0f km",$lam/1000);
#print "\n$per_lab, $c0_lab, $lam_lab \n";

$fileR1 = "$dir0/reference_values.dat";
$fileR2 = "$dir0/reference_period.dat";
$fileR3 = "$dir0/sigma_values.dat";
if (not -f $fileR1) {die("Check if fileR1 $fileR1 exist or not\n");}
if (not -f $fileR2) {die("Check if fileR2 $fileR2 exist or not\n");}
if (not -f $fileR3) {die("Check if fileR3 $fileR3 exist or not\n");}

open(IN,$fileR1); @lines = <IN>; ($alpha0,$beta0) = split(" ",$lines[0]);
open(IN,$fileR2); $per = <IN>; chomp($per);
$lam = $beta0*$per;
$per_lab = sprintf("T = %3.3f s",$per);
$beta0_lab  = sprintf("beta0 = %3.3f km/s",$beta0/1000);
$lam_lab  = sprintf("\@~\154\@~ = %.0f km",$lam/1000);
print "\n$per_lab, $beta0_lab, $lam_lab\n";

open(IN,$fileR3); @lines = <IN>;
($sigma_B,undef) = split(" ",$lines[0]);
($sigma_ts,undef) = split(" ",$lines[1]);
($sigma_xs,undef) = split(" ",$lines[2]);
($sigma_zs,undef) = split(" ",$lines[3]);
print "\n SIGMAS: $sigma_B, $sigma_ts, $sigma_xs, $sigma_zs\n";
#die("TESTING");

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
$file1dat_str = "$dir/${mfile_dat_str}";
$file1syn_str = "$dir/${mfile_syn_str}";
$file1dat_src = "$dir/${mfile_dat_src}";
$file1syn_src = "$dir/${mfile_syn_src}";

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
  if ($a2x2 > 16) {$iter_tick = 4}
  #$B3a    = "-B${iter_tick}f1g2:\" k, model number \":/a1f3g3p:\" $misfitvar ( m )   ( $utype ) \":";
  $B3a    = "-B${iter_tick}f1g2:\" k, model number \":/a1f3g3p:\" $misfitvar ( m )\":";
  $B3b    = "-B${iter_tick}f1g2:\" k, model number \":/a1f3g3p:\" \":";
}

#===========================================================================
# create colorpoint files

open(CSH,">$cshfile");

$cpt_vel_map = "../../model_files/socal_color.cpt";

# phase velocity model
$ichecker = 1;
if ($ichecker==0) {
  $cpt_vel = $cpt_vel_map;
} else {
  # make colorpoint file
  #$T1 = "-T3/4/0.1";
  $dc = $cmax/10;
  $T1 = "-T-$cmax/$cmax/$dc";
  $cpt_vel = "color0.cpt";
  print CSH "makecpt -C$colorbar $T1 -D > $cpt_vel\n";
  #print CSH "makecpt -C$colorbar $T1 > temp1\n";
  #print CSH "sed 's/^B.*/B       170     0      0  /' temp1 >  temp2\n";
  #print CSH "sed 's/^F.*/F         0     0    255  /' temp2 > $cpt_vel\n";
}

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
  $dir = "$odir$strun";

  # if the simulation stopped on a test model, use it instead
  if ( ($k == $niter_max) && (not -e $dir) ) {
     $irun = $irun - 1;
     $strun = sprintf("%4.4i",$irun);
     $dir = "$odir$strun";
   }

  $file1syn_str = "$dir/$mfile_syn_str";
  $file1syn_src = "$dir/$mfile_syn_src";
  if (not -f $file1syn_str) {die("Check if $file1syn_str exist or not\n");}
  if (not -f $file1syn_src) {die("Check if $file1syn_src exist or not\n");}
  $str_files[$k] = $file1syn_str;
  $src_files[$k] = $file1syn_src;

  # load chi values
  $chifile = "$dir/chi.dat";
  open(IN,"$chifile");
  $chi = <IN>; chomp($chi);
  $it_vals[$k] = $k;
  $chi_vals[$k] = $chi;
  #$schi = sprintf("$misfitvar ( $smod )  =  %3.3e",$chi);
}
#print "\n @it_vals \n @chi_vals \n\n"; die("testing");

# replace source and structure files with those from the Hessian run
if($ihessian == 1) {
  for ($h = 1; $h <= $qmax; $h = $h+1) {
    $dir = sprintf("$dir0/READ_IN/model_m%4.4i",$h);
    $file1syn_str = "$dir/$mfile_syn_str";
    $file1syn_src = "$dir/$mfile_syn_src";
    if (not -f $file1syn_str) {die("Check if $file1syn_str exist or not\n");}
    if (not -f $file1syn_src) {die("Check if $file1syn_src exist or not\n");}
    $str_files[$h] = $file1syn_str;
    $src_files[$h] = $file1syn_src;

    # load chi values
    $chifile = "$dir/chi.dat";
    open(IN,"$chifile");
    $chi = <IN>; chomp($chi);
    $it_vals[$h] = $h;
    $chi_vals[$h] = $chi;
  }
}

print "\n @str_files \n @src_files \n";

# get the receivers and sources
$recfile = "$dir0/${edir}/sr.txt"; # src-rec for first event
$evefile = "$dir0/events_dat_lonlat.dat"; # sources for DATA
if (not -f $recfile) {die("Check if $recfile exist or not\n");}
if (not -f $evefile) {die("Check if $evefile exist or not\n");}

# file names for figures
$name    = "joint_inversion_${strun0}_h${sth}";
$psfile  = "$name.ps";
$jpgfile = "$name.jpg";

@labs  = ("a","b","c","d","e","f","g","h","i","j");
@tlabs = ("Initial structure","Structure","Final structure","Misfit","Target structure");
if ($qmax==1) {
  @tlabs = ("Initial structure","Final structure","Misfit","Target structure");
}

#=============================================
# TOP ROW : STRUCTURE

$shift = "-X$dX1p";
$B = "$B0".$Bopts[1];
if($qmax > 1) { @kvec = (0,1,$qmax); $numk = @kvec; }

for ($i = 0; $i < $numk; $i = $i+1) {

  $t = $i;
  $k = $kvec[$i];
  $mod = $mods[2*$k];
  $smod = "m${mod}";
  #$smod = "m\@+${mod}\@+";
  $title = "($labs[$t])  $tlabs[$t] ${smod}";

  # phase velocity map
  if ($i == 0) {
    print CSH "psbasemap $B $R $J -P -K -V $origin > $psfile\n"; # START
  } else {
    print CSH "psbasemap $B $R $J -K -O -V $shift >> $psfile\n";
  }
  if ($icolor==1) {
    #print CSH "awk '{print \$1,\$2,\$6}' $str_files[$k] | pscontour $R $J -A- -C$cpt_vel -I -O -K -V >> $psfile\n";
    print CSH "awk '{print \$1,\$2,\$6}' $str_files[$k] | nearneighbor -G$grdfile $R $interp\n";
    print CSH "grdimage $grdfile -C${cpt_vel} $J -K -O -V -Q >> $psfile\n";
  }

  print "\n  $str_files[$k]\n";

  print CSH "pscoast $J $R $coast_info -K -O -V >> $psfile\n";
  print CSH "awk '{print \$1,\$2}' $evefile | psxy -N $J $R -K -O -V $src >> $psfile\n";
  print CSH "awk '\$1 == \"R\" {print \$2,\$3}' $recfile | psxy -N $J $R -K -O -V $rec >> $psfile\n";
  if ($i == 0) {
    print CSH "psscale -C$cpt_vel $Dscale $Bscale1 -K -O -V >> $psfile \n";
  }
  print CSH "pstext -N $J_title $R_title -K -O -V >>$psfile<<EOF\n $x_title $z_title $fsize_title 0 $fontno CM $title \nEOF\n";
}

#---------------------------------------------
# (2,2) subplot : misfit vs iteration

if($qmax > 1) {$shift = "-X$dX1m -Y$dY1m";} else {$shift = "-Y$dY1m";}
$Bg = "$B3b".$Bopts[8];
$t = $t+1;
$title = "($labs[$t])  $tlabs[$t],  $misfitvar (mk)";
#$title = "($labs[$t])  $tlabs[$t],  $misfitvar (mk)  (s\@+2\@+)";

print CSH "psbasemap $Bg $R_chi $J_chi -K -O -V $shift >> $psfile\n"; # START

#if ($irun0_cg != 2500) {
#  # plot curve for basic structure inversion
#  print "\n plotting reference chi curve...\n";
#  $chi_curve_ref = "chi_curve_2500.dat";
#  if (not -f $chi_curve_ref) {
#    die("Check if $chi_curve_ref exist or not\n");
#  }
#  print CSH "awk '{print \$1,\$2}' $chi_curve_ref | psxy $c_info_rd $J_chi $R_chi -K -O -V >> $psfile\n";
#}

print CSH "awk '{print \$1,\$2}' $chi_curve | psxy $c_info_rd $J_chi $R_chi -K -O -V >> $psfile\n";
for ($k = 0; $k <= $niter_max; $k = $k+1) {
  print CSH "psxy $J_chi $R_chi $p_info_k -K -O -V >>$psfile<<EOF\n $it_vals[$k] $chi_vals[$k]\nEOF\n"; # plot point
}
print CSH "pstext -N $J_title $R_title -K -O -V >>$psfile<<EOF\n $x_title $z_title $fsize_title 0 $fontno CM $title \nEOF\n";

print CSH "pstext -N $J_title $R_title -Xa-1.75 -Ya1.25 -K -O -V >>$psfile<<EOF\n 0 0 16 0 $fontno CM run_$irun0 (H-${sth}) \nEOF\n";

#---------------------------------------------
# (2,3) subplot :  target data
$B = "$B0".$Bopts[8];
$t = $t+1;
$shift = "-X$dX1p";
$title = "($labs[$t])  $tlabs[$t]";

print CSH "psbasemap $B $R $J -K -O -V $shift >> $psfile\n";
if ($icolor==1) {
  print CSH "awk '{print \$1,\$2,\$6}' $file1dat_str | nearneighbor -G$grdfile $R $interp\n";
  print CSH "grdimage $grdfile -C${cpt_vel} $J -K -O -V -Q >> $psfile\n";
}
print CSH "pscoast $J $R $coast_info -K -O -V >> $psfile\n";
print CSH "awk '{print \$1,\$2}' $evefile |psxy -N $J $R -K -O -V $src >> $psfile\n";
print CSH "awk '\$1 == \"R\" {print \$2,\$3}' $recfile |psxy -N $J $R -K -O -V $rec >> $psfile\n";
#print CSH "psscale -C$cpt_vel $Dscale $Bscale1 -K -O -V >> $psfile \n";
print CSH "pstext -N $J_title $R_title -K -O -V >>$psfile<<EOF\n $x_title $z_title $fsize_title 0 $fontno CM $title \nEOF\n";

#=============================================
# error in source parameters

# shifting the subplots
$yfac = 1.50;
$dY1p = $yfac*$Jwid; $dY1m = -$yfac*$Jwid;
$dY2p = 2*$yfac*$Jwid; $dY2m = -2*$yfac*$Jwid;
$yfac = 1.20; $dY1mB = -$yfac*$Jwid;   # same as above

$B = "$B0".$Bopts[1];

# make colorpoint file (-Z for continuous; -I to invert)
$cmax = 2*$sigma_ts;
$dc = $cmax/40;
$T2 = "-T-$cmax/$cmax/$dc";
$cptfile = "color_src.cpt";
print CSH "makecpt -Cpolar $T2 -I -D > $cptfile\n";

#---------------

@tlabs = ("Error in initial source","Error in source","Error in final source","Error in target source");
if ($qmax==1) {
  @tlabs = ("Error in initial source","Error in final source","Error in target source");
}

#@kvec = (0,$qmax); $numk = @kvec;
#for ($i = 0; $i < $numk; $i = $i+1) {
for ($i = 0; $i < $numk; $i = $i+1) {

  $k = $kvec[$i];
  #$k = $i;

  $mod = $mods[2*$k];
  $smod = "m${mod}";
  #$smod = "m\@+${mod}\@+";
  $t = $t+1;
  $title = "($labs[$t])  $tlabs[$i] ${smod}";

  #-------------------
  # KEY COMMAND -- formula to scale from mislocation to dot size
  # THIS MUST BE REPEATED FOR THE KEY BELOW!

#  $irun = $irun0_cg + 2*$k;  # wave2d run number (CG)
#  $strun = sprintf("%4.4i",$irun);
#  $dir = "$odir$strun";

#  if ( ($k == $niter_max) && (not -e $dir) ) {
#     $irun = $irun - 1;
#     $strun = sprintf("%4.4i",$irun);
#     $dir = "$odir$strun";
#   }

#  if($ihessian==0) {
#     $source_error_file = "$dir/$mfile_syn_src";
#  } else {
#     $source_error_file = $src_files[$i];
#  }

  $source_error_file = $src_files[$k];
  $dots_file = "temp";
  print CSH "awk '{print \$1,\$2,\$6,${ref_rad} + sqrt(\$7*\$7 + \$8*\$8)*$fac/1000}' $source_error_file > $dots_file\n";
  print CSH "awk '{print \$1,\$2,\$6,sqrt(\$7*\$7 + \$8*\$8)/1000}' $source_error_file > source_error\n";
  #print "\n $source_error_file \n $dots_file \n"; die("testing");
  #-------------------

  if ($i==0) {
     $shift = "-X$dX2m -Y$dY1m";
  } elsif ($i==1 || $i==2) {
     $shift = "-X$dX1p";
  }

  print CSH "psbasemap $B $R $J -K -O -V $shift >> $psfile\n";
  print CSH "pscoast $J $R $coast_info2 -K -O -V >> $psfile\n";

  print CSH "awk '\$1 == \"R\" {print \$2,\$3}' $recfile |psxy -N $J $R -K -O -V $rec >> $psfile\n";

  # plot source errors as colored and sized stars
  $dot_info = "-W0.5p -Sap";
  print CSH "psxy $dots_file $J $R $dot_info -C$cptfile -K -O -V >> $psfile \n";

  if ($i==0) {
    # origin time scale
    #$Dx = 0.25; $Dy = -0.2;
    $Dx = 0.25; $Dy = 3.0;
    $Dscale_otime = "-D$Dx/$Dy/0.8/0.10h";
    $Bscale_otime = sprintf("-B%.2ff0.1:\"Origin time error (s)\": -E8p",2*$sigma_ts);
    print CSH "psscale -C$cptfile $Dscale_otime $Bscale_otime -K -O -V >> $psfile\n";

    # convert source errors to dot size
    #$sigma_xs = 2.0;
    @misloc_dots = (0,$sigma_xs/1000,2*$sigma_xs/1000);
    $ndot = @misloc_dots;
    $xlon0 = -117; $dlon = 1.0;
    for ($j = 0; $j < $ndot; $j = $j+1) {
      $misloc = $misloc_dots[$j];
      $misloc_size[$j] = $ref_rad + $misloc*$fac; # KEY: use same formula as above
      $xlon[$j] = $xlon0 + $j*$dlon;
    }

    # source error scale -- symbols
    #$origin_box = "-Xa0 -Ya0";   # location for mislocation legend
    $origin_box = "-Xa0.0 -Ya3.25"; # location for mislocation legend
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
  $xlon[2] $yp2 $fsize2 0 $fontno CM $misloc_dots[2]
  $xlon[1] $yp2 $fsize2 0 $fontno CM $misloc_dots[1]
  $xlon[0] $yp2 $fsize2 0 $fontno CM $misloc_dots[0]
  $xlon[1] $yp1 $fsize2 0 $fontno CM Source
  $xlon[1] $yp0 $fsize2 0 $fontno CM mislocation (km)
EOF\n";
  }

  print CSH "pstext -N $J_title $R_title -K -O -V >>$psfile<<EOF\n $x_title $z_title $fsize_title 0 $fontno CM $title \nEOF\n";
}

# third column is the target source
$shift = "-Y$dY1mB";
$t = $t+1;
$title = "($labs[$t])  $tlabs[$i]";
$source_error_file = $file1dat_src;
print CSH "awk '{print \$1,\$2,\$6,${ref_rad} + sqrt(\$7*\$7 + \$8*\$8)*$fac/1000}' $source_error_file > $dots_file\n";

print CSH "psbasemap $B $R $J -K -O -V $shift >> $psfile\n";
print CSH "pscoast $J $R $coast_info2 -K -O -V >> $psfile\n";
print CSH "awk '\$1 == \"R\" {print \$2,\$3}' $recfile |psxy -N $J $R -K -O -V $rec >> $psfile\n";
print CSH "psxy $dots_file $J $R $dot_info -C$cptfile -K -O -V >> $psfile \n";
print CSH "pstext -N $J_title $R_title -K -O -V >>$psfile<<EOF\n $x_title $z_title $fsize_title 0 $fontno CM $title \nEOF\n";

#-----------------------------
print CSH "pstext -N $J_title $R_title -O -V >>$psfile<<EOF\n 10 10 $fsize_title 0 $fontno CM $title \nEOF\n"; # FINISH
print CSH "convert $psfile $jpgfile\n";
print CSH "echo done with $psfile\n";
print CSH "gv $psfile &\n";

close (CSH);
system("csh -f $cshfile");
#system("xv $jpgfile &");

#=================================================================
