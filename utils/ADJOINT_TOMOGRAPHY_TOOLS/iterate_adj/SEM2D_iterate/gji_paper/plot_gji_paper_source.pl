#!/usr/bin/perl -w

#==========================================================
#
#  plot_gji_paper_source.pl
#  Carl Tape
#  27-Feb-2006
#
#  copied from plot_gji_3wid.pl on 27-Feb-2006
#
#  ISTRUCTURE:
#    1. checker, Nfac=3
#    2. checker, Nfac=2
#    3. checker, Nfac=1
#    4. rayleigh wave (T=20s, smoothed with sig=10km)
#
#    plot_gji_paper_source.pl 25 60000 -6/3.0/0/80/8 1 2500 0 1 1 2 0 16
#    plot_gji_paper_source.pl 25 60000 -6/3.0/0/80/8 1 2550 0 1 1 2 0 16
#    plot_gji_paper_source.pl 25 60000 -6/3.0/0/80/8 1 2600 0 1 1 2 0 14
#    plot_gji_paper_source.pl 25 60000 -6/3.0/0/80/8 1 2650 0 1 1 2 0 3
#
#==========================================================

if (@ARGV < 8) {die("Usage: plot_kernels.pl nevent gamma colors iker irun0 iter ichecker ipoly istructure ibanana qmax \n");}
($nevent,$gamma,$colors,$iker,$irun0,$iter,$ichecker,$ipoly,$istructure,$ibanana,$qmax) = @ARGV;

if($ibanana==1) {$odir      = "../../../2d_adjoint_banana/OUTPUT_banana/run_"}
else            {$odir      = "../../OUTPUT/run_"}

$edir      = "event_001";  # event for first event (sr.txt)

$mfile_dat = "socal_vel_dat.dat";
$mfile_syn = "socal_vel_syn.dat";
$kerfile   = "kernel_basis";

$cshfile = "plot_gji_paper_source.csh";

if($istructure==1){$Nfac=3}
if($istructure==2){$Nfac=2}
if($istructure==3){$Nfac=1}
$sNfac = sprintf("%1i",$Nfac);

# boolean commands for plotting
$icolor = 1;   # ccc

$ifig01 = 0;   # chi(mk), recovered model, source error
$ifig02 = 1;

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
$cmax = $cols[4];  # value for phase velocity maps (percent pert)

#@files = glob("$dir/$kerfile");
#$numk = @files;

@titles = ("Waveform","Traveltime (xcorr), misfit","Amplitude (xcorr), misfit","Traveltime (MT), misfit","Amplitude (MT), misfit","Traveltime (xcorr), sampling","Amplitude (xcorr), sampling");
@units = ("m\@+2\@+ s","s\@+2\@+","xxx","xxx","xxx","xxx","xxx");
$ktype = $titles[$iker];
$utype = $units[$iker];

$plabel = "/home/carltape/sem2d/2d_adjoint/scripts/plot_ker_mod.pl";

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
$ref_rad = 6;
$src    = "-W0.5p -Sa${ref_rad}p";
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
# get reference phase velocity and period

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

# axes scale for phase velocity maps: c(th, ph)
#$bs1 = 0.5;
#$Bscale1d  = sprintf("-B%2.2e:\" Phase Velocity for data ( km s\@+-1\@+ )\":",$bs1);
#$Bscale1s  = sprintf("-B%2.2e:\" Phase Velocity for model $smod ( km s\@+-1\@+ )\":",$bs1);
$bs1 = 3;
$Bscale1  = sprintf("-B%2.2ef1:\" \@~\045\@~ pert. from %2.2f km/s\": -E5p",$bs1,$c0/1000);
$Bscale1b = sprintf("-B%2.2ef1:\" \": -E5p",$bs1);

# axes scale for kernels: K(th, ph)
# [\@~\143\@~] --> s
$tp = "\@~\146\@~, \@~\161\@~";
$Bscale2  = sprintf("-B%2.2e:\" K ( $tp )  ( 10\@+%2.2i\@+  m\@+-2\@+ s )\": -E5p",$bs2,$kpwr);
$Bscale2b = sprintf("-B%2.2e:\" \": -E5p",$bs2);

# axes scale for chi_plots: chi(th_r, ph_r)
$Bscale3  = sprintf("-B%2.2e:\" \@~\143\@~ ( \@~\161\@~\@-r\@- , \@~\146\@~\@-r\@- )  ( 10\@+%2.2i\@+ )\": -Ef5p",$bs3,$opwr);

#-------------------------
# phase velocity model
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

  # chi-vs-iteration
  $xming = $a2x1; $xmaxg = $a2x2; $yming = $a2y1; $ymaxg = $a2y2;
  $chi_ran = $xmaxg - $xming;
  $R_chi = "-R$xming/$xmaxg/$yming/$ymaxg";
  $J_chi = "-JX${Jwid}i/${Jwid}il";        # note log scale on y-axis

  # text labels for chi-vs-m plots
  $schi0 = "\@~\143\@~ ( m\@+0\@+ )";
  $schi1 = "\@~\143\@~ ( m\@+1\@+ )";
  $x0chit = 0 + 0.05*$chi_ran;
  $x1chit = 1 + 0.05*$chi_ran;

  $chi_curve = "chi_curve_${strun0}.dat";
  if (not -f $chi_curve)   { die("Check if $chi_curve exist or not\n") }

  # scale for chi-vs-m plots
  $iter_tick = 2;
  $B3a    = "-B${iter_tick}:\" k, model number \":/a1f2g1p:\" \@~\143\@~ ( m )   ( $utype ) \":";
  $B3b    = "-B${iter_tick}:\" k, model number \":/a1f2g1p:\" \":";
}

#===========================================================================
# create colorpoint files

  open(CSH,">$cshfile");

  $cpt_vel_map = "../../model_files/socal_color.cpt";

  # phase velocity model
  if($ichecker==0) {$cpt_vel = $cpt_vel_map}
  else {
     # make colorpoint file
     #$T1 = "-T3/4/0.1";
     $dc = $cmax/10;
     $T1 = "-T-$cmax/$cmax/$dc";
     $cpt_vel = "color0.cpt";
     print CSH "makecpt -C$colorbar $T1 > temp1\n";
     print CSH "sed 's/^B.*/B       170     0      0  /' temp1 >  temp2\n";
     print CSH "sed 's/^F.*/F         0     0    255  /' temp2 > $cpt_vel\n";
  }

#  # kernel
#  $cpt_ker = "color1.cpt";
#  print CSH "makecpt -C$colorbar $T2 > temp1\n";
#  print CSH "sed 's/^B.*/B       170     0      0  /' temp1 >  temp2\n";
#  print CSH "sed 's/^F.*/F         0     0    255  /' temp2 > $cpt_ker\n";

#  # misfit (as a function of receiver)
#  $cpt_chi = "color2.cpt";
#  print CSH "makecpt -Chot $T3 -I > $cpt_chi\n";

  close (CSH);
  system("csh -f $cshfile");
  #die("testing");

#===========================================================================
if ($ifig01 == 1) {

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
  @kvec = (0,$qmax);
  $numk = @kvec;
  $niter_max = $kvec[$numk-1];  # max iteration is the last in the list

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
    $chi = <IN>; chomp($chi);
    $it_vals[$k] = $k;
    $chi_vals[$k] = $chi;
    #$schi = sprintf("\@~\143\@~ ( $smod )  =  %3.3e",$chi);
  }
  #print "\n @it_vals \n @chi_vals \n\n"; die("testing");

  # get the receivers and sources
  $recfile = "$dir0/${edir}/sr.txt";          # src-rec for first event
  $evefile = "$dir0/events_lonlat_dat.dat";   # sources for DATA
  if (not -f $recfile) { die("Check if $recfile exist or not\n") }
  if (not -f $evefile) { die("Check if $evefile exist or not\n") }

  # file names for figures
  $name    = "plot_gji_paper_source_${strun0}";
  $psfile  = "$name.ps";
  $jpgfile = "$name.jpg";

  #=============================================
  # misfit vs iteration

  $shift = $shift1;
  $Bg = "$B3b".$Bopts[8];
  $title = "(a)  Misfit,  \@~\143\@~ (m\@+k\@+)  (s\@+2\@+)";

  print CSH "psbasemap $Bg $R_chi $J_chi -P -K -V $origin > $psfile\n"; # START

  if($irun0 != 2500) {
    # plot curve for basic structure inversion
    print "\n plotting reference chi curve...\n";
    $chi_curve_ref = "chi_curve_2500.dat";
    if (not -f $chi_curve_ref)   { die("Check if $chi_curve_ref exist or not\n") }
    print CSH "awk '{print \$1,\$2}' $chi_curve_ref | psxy $c_info_rd $J_chi $R_chi -K -O -P -V >> $psfile\n";
  }

  print CSH "awk '{print \$1,\$2}' $chi_curve | psxy $c_info_rd $J_chi $R_chi -K -O -P -V >> $psfile\n";
  for ($k = 0; $k <= $niter_max; $k = $k+1) {
    print CSH "psxy $J_chi $R_chi $p_info_k -K -O -V -P >>$psfile<<EOF\n $it_vals[$k] $chi_vals[$k]\nEOF\n";  # plot point
  }
  print CSH "pstext -N $J_title $R_title -K -O -V -P >>$psfile<<EOF\n $x_title $z_title $fsize_title 0 $fontno CM $title \nEOF\n";

  #=============================================
  # recovered models for synthetics

  $k = $qmax;
  $mod = $mods[2*$k]; $smod = "m\@+${mod}\@+"; $title = "(b)  Structure model ${smod}";
  $shift = $shift1;
  $B = "$B0".$Bopts[1];

  # phase velocity map
  print CSH "psbasemap $B $R $J -P -K -O -V $shift >> $psfile\n";
  if ($icolor==1) {
    print CSH "awk '{print \$1,\$2,\$4*100}' $mod_files[$k] | pscontour $R $J -A- -C$cpt_vel -I -P -O -K -V >> $psfile\n";
  }
  print CSH "pscoast $J $R $coast_info -P -K -O -V >> $psfile\n";
  if($ishelf) {print CSH "awk '{print \$2,\$1}' $shelf_file |psxy $J $R $Wshelf -K -O -P -V >> $psfile\n"}
  print CSH "awk '{print \$1,\$2}' $evefile |psxy -N $J $R -K -O -P -V $src >> $psfile\n";
  print CSH "awk '\$1 == \"R\" {print \$2,\$3}' $recfile |psxy -N $J $R -K -O -P -V $rec >> $psfile\n";
  print CSH "psscale -C$cpt_vel $Dscale $Bscale1 -P -K -O -V >> $psfile \n";
  print CSH "pstext -N $J_title $R_title -K -O -V -P >>$psfile<<EOF\n $x_title $z_title $fsize_title 0 $fontno CM $title \nEOF\n";

  #=============================================
  # error in source parameters

  $shift = $shift1;
  $B = "$B0".$Bopts[1];
  $title = "(c)  Error in source parameters";

  # KEY COMMAND -- formula to scale from magnitude to dot size
  # THIS MUST BE REPEATED FOR THE KEY BELOW!
  $fac = 3;
  $source_error_file = "gji_source_error_run_${strun0}.dat";
  $dots_file = "temp";
  if (not -f $source_error_file) { die("Check if $source_error_file exist or not\n") }
  print CSH "awk '{print \$1,\$2,\$6,${ref_rad} + \$5*$fac}' $source_error_file > $dots_file\n";

  # make colorpoint file (-Z for continuous; -I to invert)
  $cmax = 1;
  $dc = $cmax/40;
  $T2 = "-T-$cmax/$cmax/$dc";
  $cptfile = "color_src.cpt";
  print CSH "makecpt -Cpolar $T2 -I > temp1\n";
  print CSH "sed 's/^B.*/B       255     0      0  /' temp1 >  temp2\n";
  print CSH "sed 's/^F.*/F         0     0    255  /' temp2 > $cptfile\n";

  #---------------

  print CSH "psbasemap $B $R $J -P -K -O -V $shift >> $psfile\n";

  #print CSH "pscoast $J $R $coast_info -P -K -O -V >> $psfile\n";
  if($ishelf) {print CSH "awk '{print \$2,\$1}' $shelf_file |psxy $J $R $Wshelf -K -O -P -V >> $psfile\n"}
  print CSH "awk '\$1 == \"R\" {print \$2,\$3}' $recfile |psxy -N $J $R -K -O -P -V $rec >> $psfile\n";
  #print CSH "awk '{print \$1,\$2}' $evefile |psxy -N $J $R -K -O -P -V $src -G255 >> $psfile\n";

  # plot source errors as colored and sized stars
  $dot_info = "-W0.5p -Sap";
  print CSH "psxy $dots_file $J $R $dot_info -C$cptfile -K -O -V -P >> $psfile \n";

  # depth scale
  $Dscale_dep = "-D0.25/-0.2/0.8/0.10h";;
  $Bscale_dep = "-B1f0.25:\"Origin time error (s)\": -E5p";
  print CSH "psscale -C$cptfile $Dscale_dep $Bscale_dep -K -O -V -P >> $psfile\n";

  # convert source errors to dot size
  @mag_dots = (0,2,4);
  $ndot = @mag_dots;
  $xlon0 = -117; $dlon = 1.0;
  for ($k = 0; $k < $ndot; $k = $k+1) {
    $mag = $mag_dots[$k];
    $mag_size[$k] = $ref_rad + $mag*$fac;   # KEY: use same formula as above
    $xlon[$k] = $xlon0 + $k*$dlon;
  }

  # source error scale -- symbols
  $origin_box = "-Xa0 -Ya0";
  $yp1 = $zmin - 1.6;
  $yp2 = $yp1 + 0.3;
  $yp3 = $yp2 + 0.5;
  print CSH "psxy $J $R $dot_info -N -K -V -O -P $origin_box >> $psfile <<EOF
  $xlon[2] $yp3 $mag_size[2]
  $xlon[1] $yp3 $mag_size[1]
  $xlon[0] $yp3 $mag_size[0]
EOF\n";

  # source error scale -- text
  print CSH "pstext $J $R -N -K -V -O -P $origin_box >> $psfile <<EOF
  $xlon[2] $yp2 $fsize2 0 $fontno CM 4
  $xlon[1] $yp2 $fsize2 0 $fontno CM 2
  $xlon[0] $yp2 $fsize2 0 $fontno CM 0
  $xlon[1] $yp1 $fsize2 0 $fontno CM Source mislocation (km)
EOF\n";

  print CSH "pstext -N $J_title $R_title -K -O -V -P >>$psfile<<EOF\n $x_title $z_title $fsize_title 0 $fontno CM $title \nEOF\n";

#-----------------------------
  print CSH "pstext -N $J_title $R_title -O -V -P >>$psfile<<EOF\n 10 10 $fsize_title 0 $fontno CM $title \nEOF\n"; # FINISH
  print CSH "convert $psfile $jpgfile\n";

  close (CSH);
  system("csh -f $cshfile");
  system("xv $jpgfile &");
}

#===========================================================================
#===========================================================================

if ($ifig02 == 1) {

  $ipick = 2;

  # file name for figure
  $name    = sprintf("plot_gji_paper_source_f%2.2i",$ipick);
  $psfile  = "$name.ps";
  $jpgfile = "$name.jpg";

  # shifting the subplots
  $xfac = 1.20;
  $yfac = 1.40;
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

  for ($ii = 0; $ii <= 1; $ii = $ii+1) {

    if($ipick == 1) {

      if($ii == 0) {$irun0 = 2500; $qmax = 16; @labs = ("a","b","c"); @tags = (" ","(zero)");}
      else         {$irun0 = 2550; $qmax = 16; @labs = ("d","e","f"); @tags = (" "," "); }
    } else {
      @tags = ("zero"," ");
      if($ii == 0) {$irun0 = 2600; $qmax = 14; @labs = ("a","b","c"); @tags = (" ","(fixed)");}
      else         {$irun0 = 2650; $qmax = 3; @labs = ("d","e","f");  @tags = ("(fixed)"," "); }
    }

    $irun = $irun0 + $iter;
    $strun0 = sprintf("%4.4i",$irun0);
    $strun = sprintf("%4.4i",$irun);
    $dir0 = "$odir$strun0";
    $dir = "$odir$strun";

    print "\n irun0 = $irun0, irun = $irun \n";
    #die("testing");

    $chi_curve = "chi_curve_${strun0}.dat";
    if (not -f $chi_curve)   { die("Check if $chi_curve exist or not\n") }

  # iteration index (subset of all possible)
  @kvec = (0,$qmax);
  $numk = @kvec;
  $niter_max = $kvec[$numk-1];  # max iteration is the last in the list

  # load all possible files
  for ($k = 0; $k <= $niter_max; $k = $k+1) {

    $irun = $irun0 + 2*$k;              # wave2d run number

    $strun = sprintf("%4.4i",$irun);
    $dir = "$odir$strun";
    $file1syn = "$dir/$mfile_syn";
    if (not -f $file1syn)   { die("Check if $file1syn exist or not\n") }
    $mod_files[$k] = $file1syn;

    # load chi values
    $chifile = "$dir/summed_chi_all.dat";
    open(IN,"$chifile");
    $chi = <IN>; chomp($chi);
    $it_vals[$k] = $k;
    $chi_vals[$k] = $chi;
    #$schi = sprintf("\@~\143\@~ ( $smod )  =  %3.3e",$chi);
  }

  # get the receivers and sources
  $recfile = "$dir0/${edir}/sr.txt";          # src-rec for first event
  $evefile = "$dir0/events_lonlat_dat.dat";   # sources for DATA
  if (not -f $recfile) { die("Check if $recfile exist or not\n") }
  if (not -f $evefile) { die("Check if $evefile exist or not\n") }

  #=============================================
  # misfit vs iteration

  $shift = $shift1;
  $Bg = "$B3b".$Bopts[8];
  $title = "($labs[0])  Misfit,  \@~\143\@~ (m\@+k\@+)  (s\@+2\@+)";

  if($ii == 0) {
    print CSH "psbasemap $Bg $R_chi $J_chi -P -K -V $origin > $psfile\n"; # START
  } else {
    print CSH "psbasemap $Bg $R_chi $J_chi -P -K -O -V $shift2 >> $psfile\n";
  }

  if($irun0 != 2500) {
    # plot curve for basic structure inversion
    print "\n plotting reference chi curve...\n";
    $chi_curve_ref = "chi_curve_2500.dat";
    if (not -f $chi_curve_ref)   { die("Check if $chi_curve_ref exist or not\n") }
    print CSH "awk '{print \$1,\$2}' $chi_curve_ref | psxy $c_info_rd $J_chi $R_chi -K -O -P -V >> $psfile\n";
  }

  print CSH "awk '{print \$1,\$2}' $chi_curve | psxy $c_info_rd $J_chi $R_chi -K -O -P -V >> $psfile\n";
  for ($k = 0; $k <= $niter_max; $k = $k+1) {
    print CSH "psxy $J_chi $R_chi $p_info_k -K -O -V -P >>$psfile<<EOF\n $it_vals[$k] $chi_vals[$k]\nEOF\n";  # plot point
  }
  print CSH "pstext -N $J_title $R_title -K -O -V -P >>$psfile<<EOF\n $x_title $z_title $fsize_title 0 $fontno CM $title \nEOF\n";

  #=============================================
  # recovered models for synthetics

  $k = $qmax;
  $mod = $mods[2*$k]; $smod = "m\@+${mod}\@+"; $title = "($labs[1])  Structure model ${smod}  $tags[0]";
  $shift = $shift1;
  if($ii == 0) {$B = "$B0".$Bopts[8]} else {$B = "$B0".$Bopts[1]}

  # phase velocity map
  print CSH "psbasemap $B $R $J -P -K -O -V $shift >> $psfile\n";
  if ($icolor==1) {
    print CSH "awk '{print \$1,\$2,\$4*100}' $mod_files[$k] | pscontour $R $J -A- -C$cpt_vel -I -P -O -K -V >> $psfile\n";
  }
  print CSH "pscoast $J $R $coast_info -P -K -O -V >> $psfile\n";
  if($ishelf) {print CSH "awk '{print \$2,\$1}' $shelf_file |psxy $J $R $Wshelf -K -O -P -V >> $psfile\n"}
  print CSH "awk '{print \$1,\$2}' $evefile |psxy -N $J $R -K -O -P -V $src >> $psfile\n";
  print CSH "awk '\$1 == \"R\" {print \$2,\$3}' $recfile |psxy -N $J $R -K -O -P -V $rec >> $psfile\n";
  if($ii == 1) {print CSH "psscale -C$cpt_vel $Dscale $Bscale1 -P -K -O -V >> $psfile \n";}
  print CSH "pstext -N $J_title $R_title -K -O -V -P >>$psfile<<EOF\n $x_title $z_title $fsize_title 0 $fontno CM $title \nEOF\n";

  #=============================================
  # error in source parameters

  $shift = $shift1;
  if($ii == 0) {$B = "$B0".$Bopts[8]} else {$B = "$B0".$Bopts[1]}
  $title = "($labs[2])  Error in source parameters  $tags[1]";

  # KEY COMMAND -- formula to scale from magnitude to dot size
  # THIS MUST BE REPEATED FOR THE KEY BELOW!
  $fac = 3;
  $source_error_file = "gji_source_error_run_${strun0}.dat";
  $dots_file = "temp";
  if (not -f $source_error_file) { die("Check if $source_error_file exist or not\n") }
  print CSH "awk '{print \$1,\$2,\$6,${ref_rad} + \$5*$fac}' $source_error_file > $dots_file\n";

  # make colorpoint file (-Z for continuous; -I to invert)
  $cmax = 1;
  $dc = $cmax/40;
  $T2 = "-T-$cmax/$cmax/$dc";
  $cptfile = "color_src.cpt";
  print CSH "makecpt -Cpolar $T2 -I > temp1\n";
  print CSH "sed 's/^B.*/B       255     0      0  /' temp1 >  temp2\n";
  print CSH "sed 's/^F.*/F         0     0    255  /' temp2 > $cptfile\n";

  #---------------

  print CSH "psbasemap $B $R $J -P -K -O -V $shift >> $psfile\n";

  #print CSH "pscoast $J $R $coast_info -P -K -O -V >> $psfile\n";
  if($ishelf) {print CSH "awk '{print \$2,\$1}' $shelf_file |psxy $J $R $Wshelf -K -O -P -V >> $psfile\n"}
  print CSH "awk '\$1 == \"R\" {print \$2,\$3}' $recfile |psxy -N $J $R -K -O -P -V $rec >> $psfile\n";
  #print CSH "awk '{print \$1,\$2}' $evefile |psxy -N $J $R -K -O -P -V $src -G255 >> $psfile\n";

  # plot source errors as colored and sized stars
  $dot_info = "-W0.5p -Sap";
  print CSH "psxy $dots_file $J $R $dot_info -C$cptfile -K -O -V -P >> $psfile \n";

  if($ii == 1) {

  # depth scale
  $Dscale_dep = "-D0.25/-0.2/0.8/0.10h";;
  $Bscale_dep = "-B1f0.25:\"Origin time error (s)\": -E5p";
  print CSH "psscale -C$cptfile $Dscale_dep $Bscale_dep -K -O -V -P >> $psfile\n";

  # convert source errors to dot size
  @mag_dots = (0,2,4);
  $ndot = @mag_dots;
  $xlon0 = -117; $dlon = 1.0;
  for ($k = 0; $k < $ndot; $k = $k+1) {
    $mag = $mag_dots[$k];
    $mag_size[$k] = $ref_rad + $mag*$fac;   # KEY: use same formula as above
    $xlon[$k] = $xlon0 + $k*$dlon;
  }

  # source error scale -- symbols
  $origin_box = "-Xa0 -Ya0";
  $yp1 = $zmin - 1.6;
  $yp2 = $yp1 + 0.3;
  $yp3 = $yp2 + 0.5;
  print CSH "psxy $J $R $dot_info -N -K -V -O -P $origin_box >> $psfile <<EOF
  $xlon[2] $yp3 $mag_size[2]
  $xlon[1] $yp3 $mag_size[1]
  $xlon[0] $yp3 $mag_size[0]
EOF\n";

  # source error scale -- text
  print CSH "pstext $J $R -N -K -V -O -P $origin_box >> $psfile <<EOF
  $xlon[2] $yp2 $fsize2 0 $fontno CM 4
  $xlon[1] $yp2 $fsize2 0 $fontno CM 2
  $xlon[0] $yp2 $fsize2 0 $fontno CM 0
  $xlon[1] $yp1 $fsize2 0 $fontno CM Source mislocation (km)
EOF\n";

  }

  if($ii == 1 && $ipick == 1) {
    # plot a symbol next to a particular event (GJI figure)
    print "\n dir is $dir0 \n";
    $evefile = "$dir0/events_lonlat_dat.dat";
    if (not -f $evefile) { die("Check if $evefile exist or not\n") }
    open(IN,"$evefile"); @events = <IN>;
    $ievent = 25;
    ($elon,$elat,$num) = split(" ",$events[$ievent-1]);
    print "\n $elon $elat $num \n";

    $xtxt = $elon-0.05;
    $ytxt = $elat+0.3;
    print CSH "pstext -N $J $R -K -O -V -P >>$psfile<<EOF\n $xtxt $ytxt 10 0 1 CM S \nEOF\n";
    #print CSH "psxy -N $J $R -K -O -V -P $src -G0/255/0 >>$psfile<<EOF\n $xtxt $ytxt\nEOF\n";
  }

  print CSH "pstext -N $J_title $R_title -K -O -V -P >>$psfile<<EOF\n $x_title $z_title $fsize_title 0 $fontno CM $title \nEOF\n";

} # for ii

#-----------------------------
  print CSH "pstext -N $J_title $R_title -O -V -P >>$psfile<<EOF\n 10 10 $fsize_title 0 $fontno CM $title \nEOF\n"; # FINISH
  print CSH "convert $psfile $jpgfile\n";

  close (CSH);
  system("csh -f $cshfile");
  system("xv $jpgfile &");
}

#=================================================================
