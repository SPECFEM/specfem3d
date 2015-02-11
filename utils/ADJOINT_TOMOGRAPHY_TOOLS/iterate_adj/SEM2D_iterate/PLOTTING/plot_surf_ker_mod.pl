#!/usr/bin/perl -w

#==========================================================
#
#  plot_ker_mod.pl
#  Carl Tape
#  04-Aug-2006
#
#  This script generates a figure with four subplots:
#     1. target model (data)
#     2. current model (synthetics)
#     3. misfit map
#     4. event kernel
#  Option to plot event kernels (isum=0) or summed event kernels (isum=1)
#
#  EXAMPLES:
#------------------
# traveltime (irun = 600--616), Rayleigh wave, sigma = 10km
#    scripts/plot_ker_mod.pl OUTPUT/run_ event_ socal_vel_dat.dat socal_vel_syn.dat kernel_basis -6/3.0/0/20 1 26 1 600 0
#    scripts/plot_ker_mod.pl OUTPUT/run_ event_ socal_vel_dat.dat socal_vel_syn.dat kernel_basis -6/1.0/0/20 1 26 1 600 1
#    scripts/plot_ker_mod.pl OUTPUT/run_ event_ socal_vel_dat.dat socal_vel_syn.dat kernel_basis -6/1.0/0/20 1 26 1 600 2
#    scripts/plot_ker_mod.pl OUTPUT/run_ event_ socal_vel_dat.dat socal_vel_syn.dat kernel_basis -6/1.0/0/10 1 26 1 600 3
#    scripts/plot_ker_mod.pl OUTPUT/run_ event_ socal_vel_dat.dat socal_vel_syn.dat kernel_basis -7/4.0/0/10 1 26 1 600 4
#    scripts/plot_ker_mod.pl OUTPUT/run_ event_ socal_vel_dat.dat socal_vel_syn.dat kernel_basis -7/4.0/0/10 1 26 1 600 16
#
# traveltime: INDIVIDUAL kernels and summed kernel
# scripts/plot_ker_mod.pl OUTPUT/run_ event_ socal_vel_dat.dat socal_vel_syn.dat kernel_basis -7/3.0/0/80 1 5 0 700 0
#
# VARY THE SCALELENGTH AND SMOOTHING
# scripts/plot_ker_mod.pl OUTPUT/run_ event_ socal_vel_dat.dat socal_vel_syn.dat kernel_basis -7/3.0/0/10 1 25 1 420 16  # 1
# scripts/plot_ker_mod.pl OUTPUT/run_ event_ socal_vel_dat.dat socal_vel_syn.dat kernel_basis -7/3.0/0/10 1 25 1 800 16  # 2
# scripts/plot_ker_mod.pl OUTPUT/run_ event_ socal_vel_dat.dat socal_vel_syn.dat kernel_basis -7/3.0/0/10 1 25 1 780 16  # 3
# scripts/plot_ker_mod.pl OUTPUT/run_ event_ socal_vel_dat.dat socal_vel_syn.dat kernel_basis -7/3.0/0/10 1 25 1 460 16  # 4 -- main
# scripts/plot_ker_mod.pl OUTPUT/run_ event_ socal_vel_dat.dat socal_vel_syn.dat kernel_basis -7/3.0/0/10 1 25 1 820 16  # 5
# scripts/plot_ker_mod.pl OUTPUT/run_ event_ socal_vel_dat.dat socal_vel_syn.dat kernel_basis -7/3.0/0/10 1 25 1 900 16  # 6
# scripts/plot_ker_mod.pl OUTPUT/run_ event_ socal_vel_dat.dat socal_vel_syn.dat kernel_basis -7/3.0/0/10 1 25 1 500 16  # 7
# scripts/plot_ker_mod.pl OUTPUT/run_ event_ socal_vel_dat.dat socal_vel_syn.dat kernel_basis -7/3.0/0/10 1 25 1 840 16  # 8
# scripts/plot_ker_mod.pl OUTPUT/run_ event_ socal_vel_dat.dat socal_vel_syn.dat kernel_basis -7/3.0/0/10 1 25 1 920 14  # 9
#
# VARY THE NUMBER OF EVENTS (Nfac=1, gamma=30)
# scripts/plot_ker_mod.pl OUTPUT/run_ event_ socal_vel_dat.dat socal_vel_syn.dat kernel_basis -7/3.0/0/10 1  5 1 700 16
# scripts/plot_ker_mod.pl OUTPUT/run_ event_ socal_vel_dat.dat socal_vel_syn.dat kernel_basis -7/3.0/0/10 1 10 1 720 16
# scripts/plot_ker_mod.pl OUTPUT/run_ event_ socal_vel_dat.dat socal_vel_syn.dat kernel_basis -7/3.0/0/10 1 15 1 740 16
# scripts/plot_ker_mod.pl OUTPUT/run_ event_ socal_vel_dat.dat socal_vel_syn.dat kernel_basis -7/3.0/0/10 1 20 1 760 16
# scripts/plot_ker_mod.pl OUTPUT/run_ event_ socal_vel_dat.dat socal_vel_syn.dat kernel_basis -7/3.0/0/10 1 25 1 780 16
#
# VARY THE NUMBER OF EVENTS (Nfac=3, gamma=30)
# scripts/plot_ker_mod.pl OUTPUT/run_ event_ socal_vel_dat.dat socal_vel_syn.dat kernel_basis -7/3.0/0/10 1  5 1 300 16
# scripts/plot_ker_mod.pl OUTPUT/run_ event_ socal_vel_dat.dat socal_vel_syn.dat kernel_basis -7/3.0/0/10 1 10 1 320 16
# scripts/plot_ker_mod.pl OUTPUT/run_ event_ socal_vel_dat.dat socal_vel_syn.dat kernel_basis -7/3.0/0/10 1 15 1 340 16
# scripts/plot_ker_mod.pl OUTPUT/run_ event_ socal_vel_dat.dat socal_vel_syn.dat kernel_basis -7/3.0/0/10 1 20 1 360 16
# scripts/plot_ker_mod.pl OUTPUT/run_ event_ socal_vel_dat.dat socal_vel_syn.dat kernel_basis -7/3.0/0/10 1 25 1 380 16
#
# VARY THE SMOOTHING (Nfac=3, gamma = 15,30,45,60,75,90)
# scripts/plot_ker_mod.pl OUTPUT/run_ event_ socal_vel_dat.dat socal_vel_syn.dat kernel_basis -7/3.0/0/10 1 1/25 1 400 16
# scripts/plot_ker_mod.pl OUTPUT/run_ event_ socal_vel_dat.dat socal_vel_syn.dat kernel_basis -7/3.0/0/10 1 1/25 1 420 16
# scripts/plot_ker_mod.pl OUTPUT/run_ event_ socal_vel_dat.dat socal_vel_syn.dat kernel_basis -7/3.0/0/10 1 1/25 1 440 16
# scripts/plot_ker_mod.pl OUTPUT/run_ event_ socal_vel_dat.dat socal_vel_syn.dat kernel_basis -7/3.0/0/10 1 1/25 1 460 16
# scripts/plot_ker_mod.pl OUTPUT/run_ event_ socal_vel_dat.dat socal_vel_syn.dat kernel_basis -7/3.0/0/10 1 1/25 1 480 16
# scripts/plot_ker_mod.pl OUTPUT/run_ event_ socal_vel_dat.dat socal_vel_syn.dat kernel_basis -7/3.0/0/10 1 1/25 1 500 16
#
# PERTURBED CONDITIONS (Rayleigh, gamma = 30)
# scripts/plot_ker_mod.pl OUTPUT/run_ event_ socal_vel_dat.dat socal_vel_syn.dat kernel_basis -7/3.0/0/10 1 1/25 1 280 16  # reference
# scripts/plot_ker_mod.pl OUTPUT/run_ event_ socal_vel_dat.dat socal_vel_syn.dat kernel_basis -6/1.0/1/2  1 1/25 1 520 16  # pert source
#
# scripts/plot_ker_mod.pl OUTPUT/run_ event_ socal_vel_dat.dat socal_vel_syn.dat kernel_basis -7/3.0/0/10 1 1/25 1 540 16
# scripts/plot_ker_mod.pl OUTPUT/run_ event_ socal_vel_dat.dat socal_vel_syn.dat kernel_basis -7/3.0/0/10 1 1/25 1 560 16
#
# (wave2d_3.f90, 19-Nov-2005)
# scripts/plot_ker_mod.pl OUTPUT/run_ event_ socal_vel_dat.dat socal_vel_syn.dat kernel_basis -8/1.0/0/1 1 1 0 580 0
# scripts/plot_ker_mod.pl OUTPUT/run_ event_ socal_vel_dat.dat socal_vel_syn.dat kernel_basis -8/1.0/0/1 1 1 0 581 0
# scripts/plot_ker_mod.pl OUTPUT/run_ event_ socal_vel_dat.dat socal_vel_syn.dat kernel_basis -8/1.0/0/1 1 1 0 582 0
#
# GJI figure
# scripts/plot_ker_mod.pl OUTPUT/run_ event_ socal_vel_dat.dat socal_vel_syn.dat kernel_basis -8/5.0/0/8 1 5/5 0 000 0
# scripts/plot_ker_mod.pl OUTPUT/run_ event_ socal_vel_dat.dat socal_vel_syn.dat kernel_basis -8/5.0/0/8 1 1/1 0 001 0
#
# (31-Jan-2006, single kernel)
# scripts/plot_ker_mod.pl OUTPUT/run_ event_ socal_vel_dat.dat socal_vel_syn.dat kernel_basis -8/5.0/0/1 1 1/1 0 012 0
#
#==========================================================

if (@ARGV < 11) {die("Usage: plot_kernels.pl out_dir_pre event_dir_pre modelfile_dat modelfile_syn kernelfile kpwr/kmax/oprw/omax IKER iev1/iev2 isum irun0 iter model \n");}
($odir,$edir,$mfile_dat,$mfile_syn,$kerfile,$colors,$iker,$ievents,$isum,$irun0,$iter) = @ARGV;

$irun = $irun0 + $iter;

$ipvel = 0;  # use a fixed color map if useing Rayleigh wave phase velocity maps

@mods = ("0","0t","1","1t","2","2t","3","3t","4","4t","5","5t","6","6t","7","7t","8","8t");
$mod = $mods[$iter];

# event indices
@ievs = split("/",$ievents);
$ieventmin = $ievs[0];
$ieventmax = $ievs[1];
$nevent = $ieventmax - $ieventmin + 1;
$stev1 = sprintf("%3.3i",$ieventmin);     # first event dir
print "\n $ievents, $ieventmin, $ieventmax, $nevent\n";

# model iteration
$smod = "m\@+$mod\@+";

# wave2d run number
$strun0 = sprintf("%4.4i",$irun0);
$strun = sprintf("%4.4i",$irun);
$dir0 = "$odir$strun0";
$dir = "$odir$strun";

print "\n $dir, $edir, $mfile_dat, $mfile_syn, $kerfile";
print "\n $colors, $iker, $nevent, $isum, $mod \n";
#die("testing");

@cols = split("/",$colors);
$kpwr = $cols[0];
$kmax = $cols[1];
$opwr = $cols[2];
$omax = $cols[3];

#@files = glob("$dir/$kerfile");
#$numk = @files;

@titles = ("Waveform","Traveltime (xcorr), misfit","Amplitude (xcorr), misfit","Traveltime (MT), misfit","Amplitude (MT), misfit","Traveltime (xcorr), sampling","Amplitude (xcorr), sampling");
$ktype = $titles[$iker];

$plabel = "/home/carltape/sem2d/2d_adjoint/scripts/plot_ker_mod.pl";

# plotting specifications
$fsize0 = "14";
$fsize1 = "11";
$fsize2 = "9";
$fsize3 = "5";
$fontno = "1";    # 1 or 4
$tick   = "0.15c";

# plot symbols for sources, receivers, and shelf
$src = "-W1.0p -Sa0.20 -G255/0/0";
$rec = "-W0.5p/0/0/0 -St0.10";
$rec0 = "-Sc10p -W0.5p";
$Wshelf = "-W1.0/0/0/0tap";

$shelf_file = "INPUT/oms_shelf";

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
$bs3 = sprintf("%3.3e",0.9*$ss/4); # colorbar
$T3 = sprintf("-T%3.3e/%3.3e/%3.3e",0,$ss,$ds);
print "T3 = $T3\n";

#-------------------------
# get reference phase velocity and period

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
$Jwid = 3.25;
$J = "-JM${Jwid}i";      # in lat-lon
$origin = "-X0.6 -Y6.25";

$B1dat = "-B1:.\" \":WeSN";
$B1syn = "-B1:.\" \":wESN";
$B2 = "-B1:.\" \":wESN";
$B3 = "-B1:.\" \":WeSN";

$Dlen = 2.0;
$Dx = $Jwid/2;
$Dy = -0.35;

$Dscale1 = "-D$Dx/$Dy/$Dlen/0.10h";
$Dscale2 = "-D$Dx/$Dy/$Dlen/0.10h";
$Dscale3 = "-D$Dx/$Dy/$Dlen/0.10h";

#$bs1 = 0.5;
#$Bscale1d  = sprintf("-B%2.2e:\" Phase Velocity for data ( km s\@+-1\@+ )\":",$bs1);
#$Bscale1s  = sprintf("-B%2.2e:\" Phase Velocity for model $smod ( km s\@+-1\@+ )\":",$bs1);
$cmax = 9; $bs1 = 3;
$Bscale1d  = sprintf("-B%2.2ef1:\" Structure for data  (\@~\045\@~ pert. from %2.2f km/s)\": -E10p",$bs1,$c0/1000);
$Bscale1s  = sprintf("-B%2.2ef1:\" Structure for model $smod  (\@~\045\@~ pert from %2.2f km/s)\": -E10p",$bs1,$c0/1000);
$Bscale2  = sprintf("-B%2.2e:\" K ( x, y )  ( 10\@+%2.2i\@+  m\@+-2\@+ [\@~d\143\@~] )\": -E10p",$bs2,$kpwr);
$Bscale3  = sprintf("-B%2.2e:\" \@~\143\@~ ( x\@-r\@- , y\@-r\@- )  ( 10\@+%2.2i\@+ )\": -Ef10p",$bs3,$opwr);

#-------------------------
# phase velocity model
$file1dat = "$dir/$mfile_dat";
$file1syn = "$dir/$mfile_syn";

$title = "Rayleigh Wave Phase Velocity";

# set bounds for the plotting
#$xmin = -121; $xmax = -115.0; $zmin = 31.75; $zmax = 36.5;
($xmin,$xmax,$zmin,$zmax) = split(" ",`minmax -C $file1dat`);
$dinc = 0.25;  # buffer around min/max
$xmin = $xmin-$dinc;  $xmax = $xmax+$dinc;
$zmin = $zmin-$dinc;  $zmax = $zmax+$dinc;
$R = "-R$xmin/$xmax/$zmin/$zmax";

# plot title
$xtx = $xmin+0.5*($xmax-$xmin);
$ztx = $zmin+1.10*($zmax-$zmin);

#-------------------------
# 3-plot figure:
#   (1) upper left : model
#   (2) lower left : kernel
#   (3) middle right : data misfit (chi)
#-------------------------

#===========================================================================
if ($isum == 0) {
print "\n Plotting event kernels...\n";
#===========================================================================

for ($k = $ieventmin-1; $k <= $ieventmax-1; $k = $k+1) {
#for ($k = 0; $k <= $nevent-1; $k = $k+1) {

  #===============================================
  print "\nWriting CSH file...\n";
  $cshfile = "plot_ker_mod.csh";
  open(CSH,">$cshfile");
  print CSH "gmtset BASEMAP_TYPE plain TICK_LENGTH $tick LABEL_FONT_SIZE $fsize2 ANOT_FONT_SIZE $fsize2 PLOT_DEGREE_FORMAT D HEADER_FONT $fontno ANOT_FONT $fontno LABEL_FONT $fontno HEADER_FONT_SIZE $fsize1\n";
  #===============================================

  # get the kernel and chi files
  $ev = $k+1;
  $stev = sprintf("%3.3i",$ev);
  $ev_dir = "$dir/$edir$stev";
  $file2 = "$ev_dir/$kerfile";
  $file3 = "$ev_dir/chi_r.dat";
  $srcfile = "$ev_dir/sr.txt";
  $chifile = "$ev_dir/single_chi_r.dat";
  if (not -f $file2) { die("Check if $file2 exist or not\n") }
  if (not -f $file3) { die("Check if $file3 exist or not\n") }
  if (not -f $srcfile) { die("Check if $srcfile exist or not\n") }
  if (not -f $chifile) { die("Check if $chifile exist or not\n") }

  # file names for figures
  $name    = "$dir/event_$stev";
  $psfile  = "$name.ps";
  $jpgfile = "$name.jpg";

  # number of receivers
  open(IN,$file3); @temp = <IN>; $nrec = @temp;

  $npath = $nrec * 1;
  $plot_title = "\"Actual\" structure, model, misfit, and summed kernel ($ktype, EVENT $stev) -- irun = $strun";
  $plot_title2 = "1 event, $nrec receivers, $npath paths (seismograms)";

#}
#die("testing");
#for ($k = 0; $k <= $nevent-1; $k = $k+1) {

  #=============================================
  # model for data

  if($ipvel==1) {$cpt_vel = "model_files/socal_color.cpt"}
  else {
     # make colorpoint file
     #$T1 = "-T3/4/0.1";
     $dc = $cmax/10;
     $T1 = "-T-$cmax/$cmax/$dc";
     $cpt_vel = "color0.cpt";
     print CSH "makecpt -C$colorbar $T1 > temp1\n";
     print CSH "sed 's/^B.*/B       170     0      0  /' temp1  >  temp2\n";
     print CSH "sed 's/^F.*/F         0     0    255  /' temp2 > $cpt_vel\n";
  }

  # phase velocity map
  #print CSH "awk '{print \$1,\$2,\$3/1000}' $file1dat | pscontour $R $J -A- -C$cpt_vel -I -P -K -V $origin > $psfile\n"; # START
  print CSH "awk '{print \$1,\$2,\$4*100}' $file1dat | pscontour $R $J -A- -C$cpt_vel -I -P -K -V $origin > $psfile\n"; # START
  print CSH "pscoast $J $R $B1dat -W1p -Na/1p -Dh -P -K -O -V >> $psfile\n";
  #print CSH "awk '{print \$2,\$1}' ../../INPUT/oms_shelf |psxy $J $R $Wshelf -K -O -P -V >> $psfile\n";

  # plot receivers with numbered label
  print CSH "awk '\$1 == \"R\" {print \$2,\$3}' $srcfile |psxy -N $J $R -K -O -P -V $rec0 >> $psfile\n";
  $rec_file = text_rec; $angle = 0; $just = "CM";
  print CSH "awk '\$1 == \"R\" {print \$2,\$3,$fsize3,$angle,$fontno,\"$just\",\$4}' $srcfile > $rec_file\n";
  print CSH "pstext $rec_file -N $J $R -K -O -P -V >> $psfile\n";

  print CSH "awk '\$1 == \"S\" {print \$2,\$3}' $srcfile |psxy -N $J $R -K -O -P -V $src >> $psfile\n";
  print CSH "psscale -C$cpt_vel $Dscale1 $Bscale1d -P -K -O -V >> $psfile \n"; # -I gives fancy shading
  #print CSH "pstext -N $J $R -K -O -P -V >>$psfile<<EOF\n $xtx $ztx $fsize1 0 $fontno CM $title\nEOF\n";

  # plot title and GMT header
  $Utag = "-U/0/0.25/$plabel";
  $shift = "-Xa0 -Ya4";
  $shift2 = "-Xa0 -Ya3.75";
  print CSH "pstext -N $J $R $Utag -K -O -P -V $shift >>$psfile<<EOF\n $xmin $zmin $fsize1 0 $fontno LM $plot_title\nEOF\n";
  print CSH "pstext -N $J $R -K -O -P -V $shift2 >>$psfile<<EOF\n $xmin $zmin $fsize1 0 $fontno LM $plot_title2\nEOF\n";

  #--------------------------------------------
  # model for synthetics

  $dX = 1.1*$Jwid;
  $shift = "-X$dX";

  # phase velocity map
  #print CSH "awk '{print \$1,\$2,\$3/1000}' $file1syn | pscontour $R $J -A- -C$cpt_vel $shift -I -P -O -K -V >> $psfile\n";
  print CSH "awk '{print \$1,\$2,\$4*100}' $file1syn | pscontour $R $J -A- -C$cpt_vel $shift -I -P -O -K -V >> $psfile\n";
  print CSH "pscoast $J $R $B1syn -W1p -Na/1p -Dh -P -K -O -V >> $psfile\n";
  #print CSH "awk '{print \$2,\$1}' ../../INPUT/oms_shelf |psxy $J $R $Wshelf -K -O -P -V >> $psfile\n";

  # plot source, receivers with numbered label
  print CSH "awk '\$1 == \"R\" {print \$2,\$3}' $srcfile |psxy -N $J $R -K -O -P -V $rec0 >> $psfile\n";
  $rec_file = text_rec; $angle = 0; $just = "CM";
  print CSH "awk '\$1 == \"R\" {print \$2,\$3,$fsize3,$angle,$fontno,\"$just\",\$4}' $srcfile > $rec_file\n";
  print CSH "pstext $rec_file -N $J $R -K -O -P -V >> $psfile\n";
  print CSH "awk '\$1 == \"S\" {print \$2,\$3}' $srcfile |psxy -N $J $R -K -O -P -V $src >> $psfile\n";
  print CSH "psscale -C$cpt_vel $Dscale1 $Bscale1s -P -K -O -V >> $psfile \n";

  #=============================================
  # chi as a function of receiver

  $fac = 0.40; $dY = -(1+$fac)*$Jwid;
  $shift = "-X-$dX -Y$dY";

  #$dX = 1.1*$Jwid;
  #$dY = $Jwid/2 + $fac/2*$Jwid;
  #$shift = "-X$dX -Y$dY";

  # total misfit for the event
  open(IN,"$chifile");
  $chi = <IN>;
  $schi = sprintf("\@~\143\@~\@-$stev\@- ( $smod )  =  %3.3e",$chi);

  $cptfile3 = "color_3.cpt";
  print CSH "makecpt -Chot $T3 -I > $cptfile3\n";
  #print CSH "sed 's/^B.*/B       170     0      0  /' color.cpt  >  color1.cpt\n";
  #print CSH "sed 's/^F.*/F         0    34    226  /' color1.cpt > $cptfile3\n";

  print CSH "psbasemap $J $R $B3 $shift -K -O -P -V >> $psfile\n";
  if($nrec >= 10) {
    #print CSH "awk '{print \$1,\$2,\$3 / $norm3 }' $file3 > temp\n";
    print CSH "awk '{print \$1,\$2,\$3 / $norm3 }' $file3 | pscontour $J $R $B3 -A- -C$cptfile3 -I -K -O -P -V >> $psfile\n";
    print CSH "psscale -C$cptfile3 $Dscale3 $Bscale3 -K -O -P -V >> $psfile \n";
  }
  print CSH "pscoast $J $R $B3 -W1p -Na/1p -Dh -K -O -P -V >> $psfile\n";
  #print CSH "pstext -N $J $R -K -O -P -V >>$psfile<<EOF\n $xtx $ztx $fsize1 0 $fontno CM $titles[$k]\nEOF\n";
  print CSH "awk '\$1 == \"R\" {print \$2,\$3}' $srcfile |psxy -N $J $R -K -O -P -V $rec >> $psfile\n";
  print CSH "awk '\$1 == \"S\" {print \$2,\$3}' $srcfile |psxy -N $J $R -K -O -P -V $src >> $psfile\n";

  $chiloc = "-Xa$Dx -Ya-1.1";
  print CSH "pstext -N $J $R -K -O -P -V $chiloc >>$psfile<<EOF\n $xmin $zmin $fsize1 0 $fontno CM $schi \nEOF\n";

  #=============================================
  # kernel

  #$fac = 0.35;
  #$dY = -(1+$fac)*$Jwid;
  #$shift = "-Y$dY";
  $shift = "-X$dX";

  print CSH "\n echo here we are $shift \n";

  $cptfile2 = "color_2.cpt";
  print CSH "makecpt -C$colorbar $T2 > color.cpt\n";
  print CSH "sed 's/^B.*/B       170     0      0  /' color.cpt  >  color1.cpt\n";
  print CSH "sed 's/^F.*/F         0     0    255  /' color1.cpt > $cptfile2\n";

  #($pwr,$cmax) = split(" ",$colors[$k]);
  #$norm = "1e-$kpwr";
  #print CSH "echo color parameters for $file: norm = $norm, cmax = $cmax\n";

  #print CSH "awk '{print \$1,\$2,\$3 / $norm2 }' $file > temp\n";
  print CSH "awk '{print \$1,\$2,\$3 / $norm2 }' $file2 | pscontour $J $R $B2 -A- -C$cptfile2 $shift -I -K -O -P -V >> $psfile\n";
  print CSH "pscoast $J $R $B2 -W1p -Na/1p -Dh -K -O -P -V >> $psfile\n";
  print CSH "psscale -C$cptfile2 $Dscale2 $Bscale2 -K -O -P -V >> $psfile \n";
  #print CSH "pstext -N $J $R -K -O -P -V >>$psfile<<EOF\n $xtx $ztx $fsize1 0 $fontno CM $titles[$k]\nEOF\n";
  print CSH "awk '\$1 == \"R\" {print \$2,\$3}' $srcfile |psxy -N $J $R -K -O -P -V $rec >> $psfile\n";
  print CSH "awk '\$1 == \"S\" {print \$2,\$3}' $srcfile |psxy -N $J $R -O -P -V $src >> $psfile\n";  # FINISH

  #----------------------------
  print CSH "convert $psfile $jpgfile\n";
  print CSH "echo done with $psfile\n";

  close (CSH);
  system("csh -f $cshfile");
  system("xv $jpgfile &");

}  # for k < nevent

#=================================================================
} else {      # summed event kernel (misfit kernel)
#=================================================================

  #===============================================
  print "\nWriting CSH file...\n";
  $cshfile = "plot_ker_mod.csh";
  open(CSH,">$cshfile");
  print CSH "gmtset BASEMAP_TYPE plain TICK_LENGTH $tick LABEL_FONT_SIZE $fsize2 ANOT_FONT_SIZE $fsize2 PLOT_DEGREE_FORMAT D HEADER_FONT $fontno ANOT_FONT $fontno LABEL_FONT $fontno HEADER_FONT_SIZE $fsize1\n";
  #===============================================

  $src = "-W0.5p -Sa0.15";

  # get the summed kernel and chi files
  $file2 = "$dir/summed_ker.dat";
  $file3 = "$dir/summed_chi_r.dat";
  $recfile = "$dir/${edir}$stev1/sr.txt";  # src-rec for first event
  $evefile = "$dir0/events_lonlat.dat";
  $chifile = "$dir/summed_chi_all.dat";

  if (not -f $file2)   { die("Check if $file2 exist or not\n") }
  if (not -f $file3)   { die("Check if $file3 exist or not\n") }
  if (not -f $recfile) { die("Check if $recfile exist or not\n") }
  if (not -f $chifile) { die("Check if $chifile exist or not\n") }

  # file names for figures
  $name    = "$dir/event_sum_$strun";
  $psfile  = "$name.ps";
  $jpgfile = "$name.jpg";

  # number of receivers
  open(IN,$file3); @temp = <IN>; $nrec = @temp;

  $npath = $nrec * $nevent;
  $plot_title = "\"Actual\" structure, model, misfit, and summed kernel ($ktype) -- irun = $strun";
  $plot_title2 = "$nevent events, $nrec receivers, $npath paths (seismograms) -- $c0_lab,  $per_lab,  $lam_lab";

  #=============================================
  # model for the data

  if($ipvel==1) {$cpt_vel = "model_files/socal_color.cpt"}
  else {
     # make colorpoint file
     #$T1 = "-T3/4/0.1";
     $dc = $cmax/10;
     $T1 = "-T-$cmax/$cmax/$dc";
     $cpt_vel = "color0.cpt";
     print CSH "makecpt -C$colorbar $T1 > temp1\n";
     print CSH "sed 's/^B.*/B       170     0      0  /' temp1  >  temp2\n";
     print CSH "sed 's/^F.*/F         0     0    255  /' temp2 > $cpt_vel\n";
  }

  # phase velocity map
  #print CSH "awk '{print \$1,\$2,\$3/1000}' $file1dat | pscontour $R $J -A- -C$cpt_vel -I -P -K -V $origin > $psfile\n"; # START
  print CSH "awk '{print \$1,\$2,\$4*100}' $file1dat | pscontour $R $J -A- -C$cpt_vel -I -P -K -V $origin > $psfile\n"; # START
  print CSH "pscoast $J $R $B1dat -W1p -Na/1p -Dh -P -K -O -V >> $psfile\n";
  print CSH "awk '{print \$2,\$1}' $shelf_file |psxy $J $R $Wshelf -K -O -P -V >> $psfile\n";

  # plot receivers with numbered label
  print CSH "awk '\$1 == \"R\" {print \$2,\$3}' $recfile |psxy -N $J $R -K -O -P -V $rec0 >> $psfile\n";
  $rec_file = text_rec; $angle = 0; $just = "CM";
  print CSH "awk '\$1 == \"R\" {print \$2,\$3,$fsize3,$angle,$fontno,\"$just\",\$4}' $recfile > $rec_file\n";
  print CSH "pstext $rec_file -N $J $R -K -O -P -V >> $psfile\n";

  #print CSH "awk '\$1 == \"S\" {print \$2,\$3}' $recfile |psxy -N $J $R -K -O -P -V $src >> $psfile\n";
  print CSH "awk '{print \$1,\$2}' $evefile |psxy -N $J $R -K -O -P -V $src >> $psfile\n";
  print CSH "psscale -C$cpt_vel $Dscale1 $Bscale1d -P -K -O -V >> $psfile \n";

  # plot title and GMT header
  $Utag = "-U/0/0.25/$plabel";
  $shift = "-Xa0 -Ya4";
  $shift2 = "-Xa0 -Ya3.75";
  print CSH "pstext -N $J $R $Utag -K -O -P -V $shift >>$psfile<<EOF\n $xmin $zmin $fsize1 0 $fontno LM $plot_title\nEOF\n";
  print CSH "pstext -N $J $R -K -O -P -V $shift2 >>$psfile<<EOF\n $xmin $zmin $fsize1 0 $fontno LM $plot_title2\nEOF\n";

  #--------------------------------------------
  # model for synthetics

  $dX = 1.1*$Jwid;
  $shift = "-X$dX";

  # phase velocity map
  #print CSH "awk '{print \$1,\$2,\$3/1000}' $file1syn | pscontour $R $J -A- -C$cpt_vel $shift -I -P -O -K -V >> $psfile\n";
  print CSH "awk '{print \$1,\$2,\$4*100}' $file1syn | pscontour $R $J -A- -C$cpt_vel $shift -I -P -O -K -V >> $psfile\n";
  print CSH "pscoast $J $R $B1syn -W1p -Na/1p -Dh -P -K -O -V >> $psfile\n";
  print CSH "awk '{print \$2,\$1}' $shelf_file |psxy $J $R $Wshelf -K -O -P -V >> $psfile\n";

  # plot receivers with numbered label
  print CSH "awk '\$1 == \"R\" {print \$2,\$3}' $recfile |psxy -N $J $R -K -O -P -V $rec0 >> $psfile\n";
  $rec_file = text_rec; $angle = 0; $just = "CM";
  print CSH "awk '\$1 == \"R\" {print \$2,\$3,$fsize3,$angle,$fontno,\"$just\",\$4}' $recfile > $rec_file\n";
  print CSH "pstext $rec_file -N $J $R -K -O -P -V >> $psfile\n";
  print CSH "awk '{print \$1,\$2}' $evefile |psxy -N $J $R -K -O -P -V $src >> $psfile\n";
  print CSH "psscale -C$cpt_vel $Dscale1 $Bscale1s -P -K -O -V >> $psfile \n";

  #=============================================
  # chi as a function of receiver

  $fac = 0.40; $dY = -(1+$fac)*$Jwid;
  #$shift = "-Y$dY";
  $shift = "-X-$dX -Y$dY";

  # total misfit for the model
  open(IN,"$chifile");
  $chi = <IN>;
  $schi = sprintf("\@~\143\@~ ( $smod )  =  %3.3e",$chi);

  $cptfile3 = "color_3.cpt";
  print CSH "makecpt -Chot $T3 -I > $cptfile3\n";

  #print CSH "awk '{print \$1,\$2,\$3 / $norm3 }' $file > temp\n";
  print CSH "awk '{print \$1,\$2,\$3 / $norm3 }' $file3 | pscontour $J $R $B3 -A- -C$cptfile3 $shift -I -K -O -P -V >> $psfile\n";
  print CSH "pscoast $J $R $B3 -W1p -Na/1p -Dh -K -O -P -V >> $psfile\n";
  print CSH "awk '{print \$2,\$1}' $shelf_file |psxy $J $R $Wshelf -K -O -P -V >> $psfile\n";

  print CSH "psscale -C$cptfile3 $Dscale3 $Bscale3 -K -O -P -V >> $psfile \n";
  print CSH "awk '{print \$1,\$2}' $evefile |psxy -N $J $R -K -O -P -V $src >> $psfile\n";
  print CSH "awk '\$1 == \"R\" {print \$2,\$3}' $recfile |psxy -N $J $R -K -O -P -V $rec >> $psfile\n";
  $chiloc = "-Xa$Dx -Ya-1.1";
  print CSH "pstext -N $J $R -K -O -P -V $chiloc >>$psfile<<EOF\n $xmin $zmin $fsize1 0 $fontno CM $schi \nEOF\n";

  #=============================================
  # kernel

  $shift = "-X$dX";

  $cptfile2 = "color_2.cpt";
  print CSH "makecpt -C$colorbar $T2 > color.cpt\n";
  print CSH "sed 's/^B.*/B       170     0      0  /' color.cpt  >  color1.cpt\n";
  print CSH "sed 's/^F.*/F         0     0    255  /' color1.cpt > $cptfile2\n";

  print CSH "awk '{print \$1,\$2,\$3 / $norm2 }' $file2 | pscontour $J $R $B2 -A- -C$cptfile2 $shift -I -K -O -P -V >> $psfile\n";
  print CSH "pscoast $J $R $B2 -W1p -Na/1p -Dh -K -O -P -V >> $psfile\n";
  print CSH "awk '{print \$2,\$1}' $shelf_file |psxy $J $R $Wshelf -K -O -P -V >> $psfile\n";

  print CSH "psscale -C$cptfile2 $Dscale2 $Bscale2 -K -O -P -V >> $psfile \n";
  print CSH "awk '{print \$1,\$2}' $evefile |psxy -N $J $R -K -O -P -V $src >> $psfile\n";
  print CSH "awk '\$1 == \"R\" {print \$2,\$3}' $recfile |psxy -N $J $R -O -P -V $rec >> $psfile\n"; # FINISH

#-----------------------------
  print CSH "convert $psfile $jpgfile\n";
  print CSH "echo done with $psfile\n";

  close (CSH);
  system("csh -f $cshfile");
  system("xv $jpgfile &");
  #system("rm $psfile &");
}

#=================================================================
