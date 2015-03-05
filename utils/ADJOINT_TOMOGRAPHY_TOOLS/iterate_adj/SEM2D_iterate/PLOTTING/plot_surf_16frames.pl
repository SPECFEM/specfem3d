#!/usr/bin/perl -w

#==========================================================
#
#  plot_surf_16frames.pl
#  Carl Tape
#  26-Aug-2006
#
#  Perl script to plot snapshots of the forward wavefield,
#  adjoint wavefield, interaction field, and kernel.
#
#  copied from plot_oms_kernel.pl and plot_for_adj_oms.pl on 09-June-2005
#
#-------------------------
#  EXAMPLES (GJI figures):
#    (27-Nov-2005, traveltime) -- socal kernels
#    [000] 132 recs (checker) -- GJI Fig. 8
#    [001] single src-rec (homogeneous) -- GJI Fig. 5
#    [002] single src-rec (banana-doughnut)
#    [003] 132 recs (rayleigh)
#
#  EXECUTE FROM /PLOTTING/FIGURES:
#    ../plot_surf_16frames.pl /net/denali/scratch1/carltape/OUTPUT 0 5 0/4000/400 400/1200/2000/2400 3/8/7 1.0/2.5/0.5 0.06 0 1 1
#    ../plot_surf_16frames.pl /net/denali/scratch1/carltape/OUTPUT 1 1 0/4000/400 400/1200/2000/2800 3/8/7 1.0/2.5/0.5 0.06 0 1 1
#    ../plot_surf_16frames.pl /net/denali/scratch1/carltape/OUTPUT 2 1 0/4000/400 400/1200/2000/2800 3/9/9 1.0/3.0/4.0 0.06 0 5 1
#    ../plot_surf_16frames.pl /net/denali/scratch1/carltape/OUTPUT 3 5 0/4000/400 400/1200/2000/2800 3/8/7 1.0/2.5/0.5 0.06 0 1 1
#-------------------------
#
#  01-Feb-2006
#    ../plot_surf_16frames.pl OUTPUT/run_0012/event_001 0/4000/400 800/1600/2400/3200 3/8/7 1.0/10/1    0.06 0 1 2
#    ../plot_surf_16frames.pl OUTPUT/run_1000/event_005 0/4000/400 800/1600/2400/3200 3/9/8 1.0/05/1    0.06 0 1 2
#
#    ../plot_surf_16frames.pl OUTPUT/run_3100/event_005 0/4000/400 800/1200/1600/2000 3/9/8 1.0/05/1    0.06 0 1 2
#    ../plot_surf_16frames.pl OUTPUT/run_3101/event_005 0/4000/400 800/1200/1600/2000 3/5/5 1.0/01/1    0.06 0 3 2
#
#   AGU 2006
#    ../plot_surf_16frames.pl OUTPUT 4700 1 0/4000/400 400/1200/2000/2400 3/8/8 1.0/2/5 0.06 0 1 0
#    ../plot_surf_16frames.pl OUTPUT 4750 1 0/4000/400 400/1200/2000/2400 3/8/8 1.0/2/5 0.06 0 1 0
#
#-------------------------
#  04-Jan-2008
#    ../plot_surf_16frames.pl OUTPUT_2 6000 5 0/4000/400 400/1200/2000/2400 3/8/8 1.0/2/5 0.06 0 1 2
#
#==========================================================

if (@ARGV < 9) {die("Usage: plot_surf_16frames.pl basedir f-first/f-last/f-inc f1/f2/f3 pwr1/pwr2/pwr3 c1/c2/c3 DT bool_finite_source IKER \n");}
($idir,$irun,$ievent,$range,$ftemp,$ptemp,$ctemp,$dt,$ifinite,$iker,$ilabs) = @ARGV;
($first,$end,$nint) = split("/",$range);
@frames  = split("/",$ftemp);
@pwr     = split("/",$ptemp);     # PWR  : increase (larger negative power) for more contrast
@cmax    = split("/",$ctemp);     # CMAX : decrease for more contrast
$numf = @frames;

# base directory
$strun = sprintf("%4.4i",$irun);
$stev  = sprintf("%3.3i",$ievent);
$basedir = "$idir/run_$strun/event_$stev";

# ilabs: how many labels you want on the plots
#  [0] for talks or posters
#  [1] for publication
#  [2] for notes

# header for the kernel column
$col4 = "Event  Kernel";
#$col4 = "Sensitivity  Kernel";
#$col4 = "Kernel  K ( x, y, t )";

@fig_title = ("Waveform","Traveltime (xcorr), misfit","Amplitude (xcorr), misfit","Traveltime (MT), misfit","Amplitude (MT), misfit","Traveltime (xcorr), sampling","Amplitude (xcorr), sampling");

# plot the color frames or not
$icolor = 1;    # ccc

$comp = 1;  # y-component only (surface waves)

# KEY: resolution of color plots
#$interp = "-I0.5m/0.5m -S4m";
$interp = "-I1m/1m -S4m";   # key information
$interp = "-I2m/2m -S4m";   # key information
$grdfile = "temp.grd";

# plotting specifications
$fsize0 = "18";
$fsize1 = "10";
$fsize2 = "8";
$fsize3 = "4";
$fontno = "1";
#if($ilabs == 0) {$fontno = "1"} else {$fontno = "4"}
$tick   = "0.1c";

# plot symbols for sources, receivers, and shelf
if($ifinite==0) {$src = "-W0.5p -Sa0.10";} else {$src = "-Sc0.05";}
$rec = "-W0.5p,0/0/0 -St0.08";
$Wshelf = "-W1.0p,0/0/0,--";

$src = "-W0.5p -Sa0.2 -G255";   # nice white source (GJI figures)

$cshfile = "plot_surf_16frames.csh";
open(CSH,">$cshfile");
print CSH "gmtset BASEMAP_TYPE plain PAPER_MEDIA letter PLOT_DEGREE_FORMAT D TICK_LENGTH $tick LABEL_FONT_SIZE $fsize2 ANOT_FONT_SIZE $fsize2  HEADER_FONT $fontno ANOT_FONT $fontno LABEL_FONT $fontno HEADER_FONT_SIZE $fsize1 CHAR_ENCODING Standard+\n";


# KEY: scaling for color
$scale_color = 21.0;
$colorbar = "seis";
@norm = ("1e-$pwr[0]","1e-$pwr[1]","1e-$pwr[2]");
$fac = 5;      # KEY: enhance the interaction field (increase for more contrast)

print "@norm \n";

#die("testing");

$pt = "\@~f\@~, \@~q\@~";  # phi, theta
#$pt = "\@~\146\@~, \@~\161\@~";  # phi, theta

# NOTE the syn (or data) and abr (alpha-beta-rho; or cbr or kmr) suffixes.
@wavefield = ("forward_syn","adjoint","kernel_abr");
@wavefield = ("forward","adjoint","kernel");      # GJI figures
@titles    = ("Regular  Wavefield","Adjoint  Wavefield","Interaction  Field",$col4);
$numw = @wavefield;

#$numw = 2;

for ($k = 0; $k < $numw; $k ++ ){
  $ss = 0;
  for ($j = $first; $j <= $end; $j = $j + $nint) {
    $isnap = ($j - $first) / $nint + 1;
    $snap = sprintf("%05d",$j);
    $snapshot = "${basedir}/$wavefield[$k]_$snap";
    if (not -f $snapshot) {die("check if snapshot $snapshot exist or not\n");}

    # determine abs(max) for three fields
    if($k < $numw-1) {
      ($xmin,$xmax,$zmin,$zmax,$smin,$smax) = split(" ",`minmax -C $snapshot`);
      ($ss) = sort{$b <=> $a} ($ss,abs($smin),abs($smax));
    } else {
      ($xmin,$xmax,$zmin,$zmax,$imin,$imax,$smin,$smax) = split(" ",`minmax -C $snapshot`);
      ($ss) = sort{$b <=> $a} ($ss,abs($smin),abs($smax));
    }
  }
  print "\n $type Wavefield : $wavefield[$k]\n";
  print "smax of all = $ss \n";
  print "xmin = $xmin; xmax = $xmax\n";
  print "zmin = $zmin; zmax = $zmax\n";

  $dinc = 0.25;  # buffer around min/max
  $xmin = $xmin-$dinc;  $xmax = $xmax+$dinc;
  $zmin = $zmin-$dinc;  $zmax = $zmax+$dinc;
  $region = "-R$xmin/$xmax/$zmin/$zmax";
  print "\nregion is $region\n";

  $Rc[$k] = "-R$xmin/$xmax/$zmin/$zmax";
  print "R = $Rc[$k]\n";

  $ss[$k] = $cmax[$k];
  $ds[$k] = 2*$ss[$k]/$scale_color;
  $bs[$k] = sprintf("%3.3e",0.9*$ss[$k]);  # colorbar
  $Ts[$k] = sprintf("-T%3.3e/%3.3e/%3.3e",-$ss[$k],$ss[$k],$ds[$k]);
  print "Ts = $Ts[$k]\n";

  print CSH "makecpt -C$colorbar $Ts[$k] -D > color_${k}.cpt\n";
  #print CSH "makecpt -C$colorbar $Ts[$k] > color.cpt\n";
  #print CSH "sed 's/^B.*/B       170     0      0  /' color.cpt  >  color1.cpt\n";
  #print CSH "sed 's/^F.*/F         0    34    226  /' color1.cpt >  color_${k}.cpt\n";
}

#die("testing");

#-------------------------
# color for kernels

# color for the kernel
$ss = $cmax[2];
$ds = 2*$ss/$scale_color;
$bs = sprintf("%3.3e",0.9*$ss);  # colorbar
$TsK = sprintf("-T%3.3e/%3.3e/%3.3e",-$ss,$ss,$ds);
print "TsK = $TsK\n";

# color for the interaction
$ss2 = $cmax[2] / $fac;
$ds2 = 2*$ss2/$scale_color;
$bs2 = sprintf("%3.3e",0.9*$ss2);  # colorbar
$TsI = sprintf("-T%3.3e/%3.3e/%3.3e",-$ss2,$ss2,$ds2);
print "TsI = $TsI\n";

print CSH "makecpt -C$colorbar $TsK -D > color_K.cpt\n";
#print CSH "makecpt -C$colorbar $TsK > color.cpt\n";
#print CSH "sed 's/^B.*/B       170     0      0  /' color.cpt  >  color1.cpt\n";
#print CSH "sed 's/^F.*/F         0    34    226  /' color1.cpt >  color_K.cpt\n";

print CSH "makecpt -C$colorbar $TsI -D > color_I.cpt\n";
#print CSH "makecpt -C$colorbar $TsI > color.cpt\n";
#print CSH "sed 's/^B.*/B       170     0      0  /' color.cpt  >  color1.cpt\n";
#print CSH "sed 's/^F.*/F         0    34    226  /' color1.cpt >  color_I.cpt\n";

# color bars
# NOTE: there seems to be a max LENGTH of the label string for the color bars
$BscaleS1 = sprintf("-B%2.2e:\"s ($pt, t), 10\@+-%2.2i\@+  m\":",$bs[0],$pwr[0]);  # $k = 0
$BscaleS2 = sprintf("-B%2.2e:\"s\@+\262\@+ ($pt, t), 10\@+-%2.2i\@+ kg\@+-1\@+ s  [F]\":",$bs[1],$pwr[1]);   # $k = 1
$BscaleI  = sprintf("-B%2.2e:\"K\302 ($pt, t), 10\@+-%2.2i\@+ m\@+-2\@+ s\@+-1\@+ [F]\":",$bs2,$pwr[2]);
$BscaleK  = sprintf("-B%2.2e:\"K ($pt, t), 10\@+-%2.2i\@+  m\@+-2\@+ [F]\":",$bs,$pwr[2]);
#$BscaleS1 = sprintf("-B%2.2e:\" s ( $pt, t )  ( 10\@+-%i\@+  m )\":",$bs[0],$pwr[0]);  # $k = 0
#$BscaleS2 = sprintf("-B%2.2e:\" s\@+\262\@+ ( $pt, t )  ( 10\@+-%i\@+ )\":",$bs[1],$pwr[1]);   # $k = 1
#$BscaleI  = sprintf("-B%2.2e:\" K\302 ( $pt, t )  ( 10\@+-%i\@+ )\":",$bs2,$pwr[2]);
#$BscaleK  = sprintf("-B%2.2e:\" K ( $pt, t )  ( 10\@+-%i\@+ )\":",$bs,$pwr[2]);

#$Bscale = sprintf("-B%2.2e",$bs[$k]);

#print "\n $BscaleS1 \n $BscaleS2 \n $BscaleI \n $BscaleK\n"; die("testing\n");

#-------------------------

# write plotting scripts
$J = "-JM1.6i";      # in lat-lon

$origin = "-X0.75 -Y2.00";
$dX     = "-X1.9";
$mdX    = "-X-5.7";
$dY     = "-Y1.9";

$Dscale = "-D0.85/-0.40/1.25/0.1h -E8p";            # colorbar

$name = "plot_16_${strun}_${stev}_${iker}_${ilabs}";
$psfile  = "$name.ps";
$jpgfile = "$name.jpg";

print "\nWriting CSH file...\n";

#$numf = 1;

#for ($j = $first; $j <= $end; $j = $j + $nint) {
for ($i = 0; $i < $numf; $i++) {

   $j1 = $frames[$i];           # forward frame
   $j2 = $end + $first - $j1;   # corresponding adjoint frame
   $snap1 = sprintf("%05d",$j1);
   $snap2 = sprintf("%05d",$j2);
   $time = sprintf("%04d",$j1*$dt);

   $snapshot_f = "${basedir}/$wavefield[0]_${snap1}";
   $snapshot_a = "${basedir}/$wavefield[1]_${snap2}";
   $snapshot_k = "${basedir}/$wavefield[2]_${snap2}";

   if (not -f $snapshot_f) {die("check if snapshot_f $snapshot_f exist or not\n");}
   if (not -f $snapshot_a) {die("check if snapshot_a $snapshot_a exist or not\n");}
   if (not -f $snapshot_k) {die("check if snapshot_k $snapshot_k exist or not\n");}

   $k = 0;
   $R = $Rc[$k];

   print CSH "echo $psfile\n";
   print CSH "echo $snapshot_f\n";

   $B = "-B1/1:\"t = $time s\"::.\"  \":Wsne";
   $B_row1 = "-B1/1:\"t = $time s\"::.\"  \":WSne";
   if ($i == 0) { $B = $B_row1;}
   #$B = "-B2/2:.\" \":wesn";

   if ($i == 0) {print CSH "psbasemap $J $R $B -K -P -V $origin > $psfile\n"} # START
   else {print CSH "psbasemap $J $R $B -K -O -P -V $dY $mdX >> $psfile\n"}

   # PLOT THE FORWARD WAVEFIELD
   if($icolor==1) {
      print CSH "awk '{print \$1,\$2,\$($comp+2) / $norm[$k]}' $snapshot_f | nearneighbor -G$grdfile $R $interp\n";
      print CSH "grdimage $grdfile -Ccolor_${k}.cpt $J -K -O -P -V -Q >> $psfile\n";
   }
   print CSH "pscoast $J $R $B -W1p -Na/1p -Dh -K -O -P -V >> $psfile\n";
   #print CSH "awk '{print \$2,\$1}' INPUT/oms_shelf |psxy $J $R $Wshelf -K -O -P -V >> $psfile\n";
   if ($i == 0 && $ilabs > 0) {print CSH "psscale -Ccolor_${k}.cpt $Dscale $BscaleS1 -K -O -P -V >> $psfile \n";}
   print CSH "awk '\$1 == \"R\" {print \$2,\$3}' ${basedir}/sr.txt |psxy -N $J $R -K -O -P -V $rec >> $psfile\n";
   print CSH "awk '\$1 == \"S\" {print \$2,\$3}' ${basedir}/sr.txt |psxy -N $J $R -K -O -P -V $src >> $psfile\n";

   # plot the time of the snapshot (for some reason, it won't work inside the B command)
   $xtext = $xmin-0.3*($xmax-$xmin);
   $ztext = $zmin+0.5*($zmax-$zmin);
   $tstr = "t = $time s";
   print CSH "pstext -N $J $R -K -O -P -V >>$psfile<<EOF\n $xtext $ztext $fsize1 90 $fontno CM $tstr\nEOF\n";

   # plot title
   $xtx = $xmin+0.5*($xmax-$xmin); $ztx = $zmin+1.15*($zmax-$zmin);
   if ($i == $numf-1) {print CSH "pstext -N $J $R -K -O -P -V >>$psfile<<EOF\n $xtx $ztx $fsize1 0 $fontno CM $titles[0]\nEOF\n";}

   #-------------------------

   $k = 1;
   $R = $Rc[$k];

   print CSH "echo $psfile\n";
   print CSH "echo $snapshot_a\n";

   $B = "-B1/1:.\" \":wesn";
   $B_row1 = "-B1/1:.\" \":weSn";
   if ($i == 0) { $B = $B_row1;}

   # PLOT THE ADJOINT WAVEFIELD
   print CSH "psbasemap $J $R $B -K -O -P -V $dX >> $psfile\n";
   if($icolor==1) {
      print CSH "awk '{print \$1,\$2,\$($comp+2) / $norm[$k]}' $snapshot_a | nearneighbor -G$grdfile $R $interp\n";
      print CSH "grdimage $grdfile -Ccolor_${k}.cpt $J -K -O -P -V -Q >> $psfile\n";
   }
   print CSH "pscoast $J $R $B -W1p -Na/1p -Dh -K -O -P -V >> $psfile\n";
   if ($i == 0 && $ilabs > 0) {print CSH "psscale -Ccolor_${k}.cpt $Dscale $BscaleS2 -K -O -P -V >> $psfile \n";}
   print CSH "awk '\$1 == \"R\" {print \$2,\$3}' ${basedir}/sr.txt |psxy -N $J $R -K -O -P -V $rec >> $psfile\n";
   print CSH "awk '\$1 == \"S\" {print \$2,\$3}' ${basedir}/sr.txt |psxy -N $J $R -K -O -P -V $src >> $psfile\n";
   if ($i == $numf-1) {print CSH "pstext -N $J $R -K -O -P -V >>$psfile<<EOF\n $xtx $ztx $fsize1 0 $fontno CM $titles[1]\nEOF\n";}

   #-------------------------

   $k = 2;

   print CSH "echo $psfile\n";
   print CSH "echo $snapshot_k\n";

   # PLOT THE INTERACTION FIELD
   print CSH "psbasemap $J $R $B -K -O -P -V $dX >> $psfile\n";
   if($icolor==1) {
      print CSH "awk '{print \$1,\$2,\$3 / $norm[$k]}' $snapshot_k | nearneighbor -G$grdfile $R $interp\n";
      #print CSH "awk '{print \$1,\$2,\$7 / $norm[$k]}' $snapshot_k | nearneighbor -G$grdfile $R $interp\n";   # GJI figures
      print CSH "grdimage $grdfile -Ccolor_I.cpt $J -K -O -P -V -Q >> $psfile\n";
   }
   print CSH "pscoast $J $R $B -W1p -Na/1p -Dh -K -O -P -V >> $psfile\n";
   if ($i == 0 && $ilabs > 0) {print CSH "psscale -Ccolor_I.cpt $Dscale $BscaleI -K -O -P -V >> $psfile \n";}
   print CSH "awk '\$1 == \"R\" {print \$2,\$3}' ${basedir}/sr.txt |psxy -N $J $R -K -O -P -V $rec >> $psfile\n";
   print CSH "awk '\$1 == \"S\" {print \$2,\$3}' ${basedir}/sr.txt |psxy -N $J $R -K -O -P -V $src >> $psfile\n";
   if ($i == $numf-1) {print CSH "pstext -N $J $R -K -O -P -V >>$psfile<<EOF\n $xtx $ztx $fsize1 0 $fontno CM $titles[2]\nEOF\n";}

   # PLOT THE KERNEL
   print CSH "psbasemap $J $R $B -K -O -P -V $dX >> $psfile\n";
   if($icolor==1) {
      print CSH "awk '{print \$1,\$2,\$4 / $norm[$k]}' $snapshot_k | nearneighbor -G$grdfile $R $interp\n";
      print CSH "grdimage $grdfile -Ccolor_K.cpt $J -K -O -P -V -Q >> $psfile\n";
   }
   print CSH "pscoast $J $R $B -W1p -Na/1p -Dh -K -O -P -V >> $psfile\n";
   if ($i == 0 && $ilabs > 0) {print CSH "psscale -Ccolor_K.cpt $Dscale $BscaleK -K -O -P -V >> $psfile \n";}
   print CSH "awk '\$1 == \"R\" {print \$2,\$3}' ${basedir}/sr.txt |psxy -N $J $R -K -O -P -V $rec >> $psfile\n";
   print CSH "awk '\$1 == \"S\" {print \$2,\$3}' ${basedir}/sr.txt |psxy -N $J $R -K -O -P -V $src >> $psfile\n";
   if ($i == $numf-1) {print CSH "pstext -N $J $R -K -O -P -V >>$psfile<<EOF\n $xtx $ztx $fsize1 0 $fontno CM $titles[3]\nEOF\n";}
}

# plot title and GMT header
if ($ilabs == 2) {
  $plabel = "/home/carltape/wave2d/2d_adjoint/plot_surf_16frames.pl";
  $Utag = "-U/-3.75/1.95/$plabel";  # GMT header
  $shift = "-X-2.0i -Y0.7i";
  $title = "Kernel construction -- $fig_title[$iker]";
  print CSH "pstext -N $J $R $Utag -K -O -P -V $shift >>$psfile<<EOF\n $xmin $zmax $fsize0 0 $fontno CM $title\nEOF\n";
}

#-------------------------
print CSH "pstext $J -R0/1/0/1 -O -P -V >>$psfile<<EOF\n 10 10 $fsize0 0 $fontno CM junk \nEOF\n";  # FINISH
#print CSH "convert $psfile $jpgfile\n";

close (CSH);
system("csh -f $cshfile");
system("ghostview $psfile &");
