#!/usr/bin/perl -w

#==========================================================
#
#  plot_gji_2wid.pl
#  Carl Tape
#  15-March-2006
#
#    plot_gji_2wid.pl ../../OUTPUT/run_ event_ socal_vel_dat.dat socal_vel_syn.dat kernel_basis -6/3.0/0/80 1 280 0
#
#==========================================================

if (@ARGV < 9) {die("Usage: plot_gji_2wid.pl out_dir_pre event_dir_pre modelfile_dat modelfile_syn kernelfile kpwr/kmax/oprw/omax IKER irun0 iter model \n");}
($odir,$edir,$mfile_dat,$mfile_syn,$kerfile,$colors,$iker,$irun0,$iter) = @ARGV;

$ishelf = 0;   # plot the shelf as a dashed line

# boolean commands for plotting
$icolor = 1;

$ifig00 = 0;   # title page figure: SoCal topo/bath, phase speed map, faults
$ifig01 = 0;   # SoCal topo/bath, phase speed map, faults
$ifig02 = 1;   # GJI paper: SoCal topo/bath, phase speed map, src-rec geometry
$ifig03 = 0;   # src-rec geometry, ray coverage

$ifig04 = 0;   # topo/bath, including Baja CA (Ge161)

$cshfile = "plot_gji_2wid.csh";

$irun = $irun0 + $iter;

@mods = ("0","0t","1","1t","2","2t","3","3t","4","4t","5","5t","6","6t","7","7t","8","8t");
$mod = $mods[$iter];

# model iteration
$smod = "m\@+$mod\@+";

# topo files
$dir_topo = "/home/datalib/Topography/Globe_Becker/data";
$grdfile  = "$dir_topo/west_data/w140n40.Bathmetry.srtm.swap.grd";
$gradfile = "$dir_topo/west_data/w140n40.Bathmetry.srtm.swap.grad";

# boundary files
$dir0 = "/home/carltape/gmt";
$plate_file_nuvel   = "${dir0}/nuvel_1_plates";
$plate_file   = "${dir0}/plate_boundaries";
$fault_file   = "${dir0}/faults/jennings.xy";
$kcf_file     = "${dir0}/faults/kcf.xy";
$shelf_file = "../../INPUT/oms_shelf";

# wave2d run number
$strun0 = sprintf("%4.4i",$irun0);
$strun = sprintf("%4.4i",$irun);
$dir0 = "$odir$strun0";
$dir = "$odir$strun";

print "\n $dir, $edir, $mfile_dat, $mfile_syn, $kerfile";
print "\n $colors, $iker, $mod \n";
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
$fsize2 = "10";
$fsize3 = "5";
$fontno = "4";    # 1 or 4
$tick   = "0.15c";
$fpen   = "1.5p";
$tpen   = "1.0p";

# plot symbols for sources, receivers, and shelf
$src          = "-W1p -Sa12p";
$src_r        = "-W1p/255/0/0 -Sa12p";
$src_b        = "-W1p/0/0/255 -Sa12p";
$rec          = "-Sc5p -G0";
$rec_r        = "-W0.5p/0/0/0 -Sc5p -G255/0/0";
$rec_w        = "-W0.5p/0/0/0 -Sc5p -G255";
#$rec0 = "-Sc10p -W0.5p";
$Wshelf       = "-W1.0/0/0/0tap";

$plate_info_w = "-M -W2.5p/255/255/255";
$plate_info_r = "-M -W2.5p/255/0/0";
$plate_info_b = "-M -W2.5p/0/0/255";
$plate_info_k = "-M -W2.5p";

$coast_info_k  = "-W1p -Na/1p -Df -A10";
$coast_info_w  = "-W1p/255/255/255 -Na/1p -Df -A10";
$coast_info2  = "$coast_info_k -I1/1.0p";
$coast_info3  = "$coast_info_k -S220/255/255 -C220/255/255 -G200";
$coast_info_inset  = "-W0.5p -Dl -Na/0.5p -A1000 -S220/255/255 -C220/255/255 -G200";

$fault_info_r     = "-M -W1.5p/255/0/0";
$fault_info_k    = "-M -W1.5p/0/0/0";
$fault_info_w    = "-M -W1.5p/255/255/255";

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

$dir = "/home/carltape/sem2d/2d_adjoint/OUTPUT/run_0280";

open(IN,"$dir/socal_vel_c0.dat");
@vals = <IN>;
$c0       = $vals[0];
$per      = $vals[1];
$lam      = $c0*$per;
$per_lab  = sprintf("T = %.1f s",$per);
$c0_lab   =  sprintf("c0 = %.2f km/s",$c0/1000);
$lam_lab  = sprintf("\@~\154\@~ = %.0f km",$lam/1000);
print "\n$per_lab, $c0_lab, $lam_lab \n";

#die("testing");

#-------------------------

# write plotting scripts
$Jwid = 3.25;
$J = "-JM${Jwid}i";      # in lat-lon
$origin = "-X0.8 -Y6.25";

# which borders to plot the lat-lon
# 1 four sides, 4 single sides, 6 two sides, 4 three sides, 1 zero sides
@Bopts = ("WESN","Wesn","wesN","wEsn","weSn","WesN","wEsN","wESn","WeSn","WEsn","weSN","wESN","WESn","WeSN","WEsN","wesn");
$B0 = "-B1:.\" \":";

## polynomial plots
#$lpwr_x = 4;
#$lpwr_y = 3;
#$normgx = "1e$lpwr_x";
#$normgy = "1e$lpwr_y";
#$B1 = sprintf("-B0.5:\"\@~\154\@~   ( 10\@+%2.2i\@+ ) \":/1:\" \@~\143\@~\@+0\@+ [ m (\@~\154\@~) ]   ( 10\@+%2.2i\@+  s\@+2\@+ ) \":",$lpwr_x,$lpwr_y);
#$B2 = sprintf("-B0.5:\"\@~\154\@~   ( 10\@+%2.2i\@+ ) \":/1:\" \@~\143\@~\@+1\@+ [ m (\@~\154\@~) ]   ( 10\@+%2.2i\@+  s\@+2\@+ ) \":",$lpwr_x,$lpwr_y);
#$B3a = "-B1:\" model number (iteration) \":/a1f2p:\" \@~\143\@~ ( m )   ( s\@+2\@+ ) \":";
#$B3b = "-B1:\" model number (iteration) \":/a1f2p:\" \":";
#$B4  = "-B1:\" model number (iteration) \":/20:\" \":";

# colorbar
$Dlen = $Jwid*0.7;
$Dx = $Jwid/2;
$Dy = -0.12*$Jwid;
$Dscale = "-D$Dx/$Dy/$Dlen/0.15h";

#$bs1 = 0.5;
#$Bscale1d  = sprintf("-B%2.2e:\" Phase Speed for data ( km s\@+-1\@+ )\":",$bs1);
#$Bscale1s  = sprintf("-B%2.2e:\" Phase Speed for model $smod ( km s\@+-1\@+ )\":",$bs1);
$bs1 = 1; $cmax = 8;
$Bscale1a  = sprintf("-Ba%2.2e:\" \@~\045\@~ pert. from %2.2f km/s\": -E10p -A",$bs1,$c0/1000);
$Bscale1b  = sprintf("-Ba%2.2e:\" \": -E10p -A",$bs1,$c0/1000);
$bs2 = 0.05;
$Bscale1c = sprintf("-B%2.2ef1:\" Phase Speed\": -E10p",$bs2);
$Bscale2   = sprintf("-B%2.2e:\" K ( x, y )  ( 10\@+%2.2i\@+  m\@+-2\@+ [\@~d\143\@~] )\": -E10p",$bs2,$kpwr);
$Bscale3   = sprintf("-B%2.2e:\" \@~\143\@~ ( x\@-r\@- , y\@-r\@- )  ( 10\@+%2.2i\@+ )\": -Ef10p",$bs3,$opwr);

$Bscale_topo = "-B2000f500:\"Elevation (km)\":";

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
$z_title = 1.15;
$fsize_title = 12;

# shifting the subplots
$xfac = 1.1;
$yfac = 1.6;
$dX = $xfac*$Jwid; $dY = 0; $shift1 = "-X$dX -Y$dY";
$dX = -2*$dX; $dY = -$yfac*$Jwid; $shift2 = "-X$dX -Y$dY";

#===========================================================================
# create colorpoint files

  open(CSH,">$cshfile");

  # phase speed model
  if(1==1) {
     $cpt_vel  = "../../model_files/socal_color.cpt";      # ranges from -5% to 3%
     $cpt_vel2 = "../../model_files/socal_color_vel.cpt";

  } else {
     # make colorpoint file
     #$T1 = "-T3/4/0.1";
     $dc = 2*$cmax/20;
     $T1 = "-T-$cmax/$cmax/$dc";
     $cpt_vel = "color0.cpt";
     print CSH "makecpt -C$colorbar $T1 > temp1\n";
     print CSH "sed 's/^B.*/B       170     0      0  /' temp1  >  temp2\n";
     print CSH "sed 's/^F.*/F         0     0    255  /' temp2 > $cpt_vel\n";

     # make colorpoint file
     $pmin = -$cmax;
     $pmax = $cmax;
     $cmin = ($pmin/100 + 1)*$c0 / 1000;
     $cmax = ($pmax/100 + 1)*$c0 / 1000;
     $dc = ($cmax-$cmin)/20;
     $T1 = "-T$cmin/$cmax/$dc";
     $cpt_vel2 = "color0.cpt";
     print CSH "makecpt -C$colorbar $T1 > temp1\n";
     print CSH "sed 's/^B.*/B       170     0      0  /' temp1  >  temp2\n";
     print CSH "sed 's/^F.*/F         0     0    255  /' temp2 > $cpt_vel2\n";
  }

  # topography and bathymetry
  $cpt_topo = "color.cpt";
  $cmax = 5000; $dc = 2*$cmax/20;
  $T = "-T-$cmax/$cmax/$dc";
  print CSH "makecpt -Cglobe $T -Z > $cpt_topo \n";

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
if ($ifig00 == 1 || $ifig01 == 1) {

  #===============================================
  print "\nWriting CSH file...\n";

  open(CSH,">$cshfile");
  print CSH "gmtset BASEMAP_TYPE plain PLOT_DEGREE_FORMAT D TICK_LENGTH $tick LABEL_FONT_SIZE $fsize2 ANOT_FONT_SIZE $fsize2  HEADER_FONT $fontno ANOT_FONT $fontno LABEL_FONT $fontno HEADER_FONT_SIZE $fsize1 FRAME_PEN $fpen TICK_PEN $tpen \n";
  #===============================================

  # file names for figures
  if($ifig01==1) {$name = "socal_01"} else {$name = "socal_00"}
  $psfile  = "$name.ps";
  $jpgfile = "$name.jpg";

  #=============================================
  # topography and bathymetry
  if($ifig01==1) {$B = "$B0".$Bopts[13]} else {$B = "$B0".$Bopts[15]}
  $title = "Southern California";

  print CSH "psbasemap $J $R $B -K -V -P $origin > $psfile\n";  # START
  if ($icolor==1) {print CSH "grdimage $grdfile -C$cpt_topo $J $R -I$gradfile -K -O -V -P >> $psfile\n"}
  print CSH "pscoast $J $R $coast_info_k -K -O -V -P  >> $psfile \n";
  print CSH "psbasemap $J $R $B -K -V -P -O >> $psfile\n";
  print CSH "psxy $plate_file $J $R $plate_info_k -K -V -O -P >> $psfile \n";
  print CSH "psxy $fault_file $J $R $fault_info_k -K -V -O >> $psfile \n";
  print CSH "psxy $kcf_file $J $R $fault_info_k -K -V -O >> $psfile \n";
  if($ishelf == 1) {print CSH "awk '{print \$2,\$1}' $shelf_file |psxy $J $R $Wshelf -K -O -P -V >> $psfile\n"}
  if($ifig01 == 1) {
    print CSH "psscale -C$cpt_topo $Dscale $Bscale_topo -V -P -K -O >> $psfile \n";
    print CSH "pstext -N $J_title $R_title -K -O -V -P >>$psfile<<EOF\n $x_title $z_title $fsize_title 0 $fontno CM $title \nEOF\n";
  }

  #=============================================
  # phase speed map

  $title = "Rayleigh wave phase speed map";
  if($ifig01==1) {$B = "$B0".$Bopts[6]} else {$B = "$B0".$Bopts[15]}
  $shift = $shift1;

  $file01 = "$dir/socal_vel_dat.dat";

  print CSH "psbasemap $B $R $J -P -K -O -V $shift >> $psfile\n";
  if ($icolor==1) {
    print CSH "awk '{print \$1,\$2,\$4*100}' $file01 | pscontour $R $J -A- -C$cpt_vel -I -P -O -K -V >> $psfile\n";
  }
  print CSH "pscoast $J $R $coast_info_k -P -K -O -V >> $psfile\n";
  if($ishelf == 1) {print CSH "awk '{print \$2,\$1}' $shelf_file |psxy $J $R $Wshelf -K -O -P -V >> $psfile\n"}

  print CSH "psxy $plate_file $J $R $plate_info_k -K -V -O -P >> $psfile \n";
  print CSH "psxy $fault_file $J $R $fault_info_k -K -V -O >> $psfile \n";
  print CSH "psxy $kcf_file $J $R $fault_info_k -K -V -O >> $psfile \n";

  if($ifig01 == 1) {
  # scale bar
  print CSH "psscale -C$cpt_vel $Dscale $Bscale1b -P -K -O -V >> $psfile \n";
  print CSH "psscale -C$cpt_vel2 $Dscale $Bscale1c -P -K -O -V >> $psfile \n";
  print CSH "pstext -N $J $R -K -O -V -P >>$psfile<<EOF
-115.2 31.42 11 0 $fontno LM \@~\045\@~
-115.2 30.65 11 0 $fontno LM km/s
EOF\n";

  print CSH "pstext -N $J_title $R_title -K -O -V -P >>$psfile<<EOF\n $x_title $z_title $fsize_title 0 $fontno CM $title \nEOF\n";
}

  #=============================================
  # inset map
  $iwid = 0.8;
  $Jinset = "-JM$iwid";
  $x1q = -130; $x2q = -110;
  $y1q = 25; $y2q = 50;
  $Rinset = "-R$x1q/$x2q/$y1q/$y2q";
  $Binset = "-B1::wesn";
  $origin_inset = "-Xa-0.9 -Ya1.85";
  $coast_info_inset  = "-W0.5p -Dl -Na/0.5p -A1000 -S220/255/255 -C220/255/255 -G200";

  print CSH "pscoast $Jinset $Rinset $coast_info_inset $origin_inset -K -O -V -P >> $psfile \n";
  #print CSH "psxy $plate_file_nuvel $Jinset $Rinset -M -W1.5p/0/0/255 $origin_inset -K -V -O -P >> $psfile \n";
  print CSH "psxy $Jinset $Rinset -W2.0p/255/0/0 -N -K -O -V -P $origin_inset >>$psfile<<EOF\n $xmin $zmin\n$xmax $zmin\n$xmax $zmax\n$xmin $zmax\n$xmin $zmin\nEOF\n";
  #print CSH "psxy $Jinset $Rinset -W2.0p -N -K -O -V $origin_inset >>$psfile<<EOF\n $x1q $y1q\n$x2q $y1q\n$x2q $y2q\n$x1q $y2q\n$x1q $y1q\nEOF\n";

  print CSH "gmtset TICK_LENGTH 0\n";
  print CSH "psbasemap $Jinset $Rinset $Binset $origin_inset -K -O -V -P >> $psfile \n";
  print CSH "gmtset TICK_LENGTH $tick\n";

#-----------------------------
  print CSH "pstext -N $J_title $R_title -O -V -P >>$psfile<<EOF\n 10 10 $fsize_title 0 $fontno CM $title \nEOF\n"; # FINISH
  print CSH "convert $psfile $jpgfile\n";

  close (CSH);
  system("csh -f $cshfile");
  system("xv $jpgfile &");
}

#===========================================================================
if ($ifig02 == 1) {

  #===============================================
  print "\nWriting CSH file...\n";

  open(CSH,">$cshfile");
  print CSH "gmtset BASEMAP_TYPE plain PLOT_DEGREE_FORMAT D TICK_LENGTH $tick LABEL_FONT_SIZE $fsize2 ANOT_FONT_SIZE $fsize2  HEADER_FONT $fontno ANOT_FONT $fontno LABEL_FONT $fontno HEADER_FONT_SIZE $fsize1 FRAME_PEN $fpen TICK_PEN $tpen \n";
  #===============================================

  # get the receivers and sources
  $recfile = "$dir0/${edir}001/sr.txt";  # src-rec for first event
  $evefile = "$dir0/events_lonlat.dat";
  if (not -f $recfile) { die("Check if $recfile exist or not\n") }
  if (not -f $evefile) { die("Check if $evefile exist or not\n") }

  # file names for figures
  $name    = "socal_02";
  $psfile  = "$name.ps";
  $jpgfile = "$name.jpg";

  #=============================================
  # topography and bathymetry

  $B = "$B0".$Bopts[13];
  $title = "(a)  Source-receiver geometry";

  print CSH "psbasemap $J $R $B -K -V -P $origin > $psfile\n";  # START
  if ($icolor==1) {print CSH "grdimage $grdfile -C$cpt_topo $J $R -I$gradfile -K -O -V -P >> $psfile\n"}
  print CSH "pscoast $J $R $coast_info2 -K -O -V -P  >> $psfile \n";
  print CSH "psscale -C$cpt_topo $Dscale $Bscale_topo -V -P -K -O >> $psfile \n";
  print CSH "psbasemap $J $R $B -K -V -P -O >> $psfile\n";
  #print CSH "psxy $plate_file $J $R $plate_info_r -K -V -O -P >> $psfile \n";
  if($ishelf == 1) {print CSH "awk '{print \$2,\$1}' $shelf_file |psxy $J $R $Wshelf -K -O -P -V >> $psfile\n"}
  print CSH "awk '\$1 == \"R\" {print \$2,\$3}' $recfile |psxy -N $J $R -K -O -P -V $rec_w >> $psfile\n";
  print CSH "awk '{print \$1,\$2}' $evefile |psxy -N $J $R -K -O -P -V $src >> $psfile\n";
  print CSH "pstext -N $J_title $R_title -K -O -V -P >>$psfile<<EOF\n $x_title $z_title $fsize_title 0 $fontno CM $title \nEOF\n";

  #=============================================
  # phase speed map

  $title = "(b)  Rayleigh wave phase speed map";
  $shift = $shift1;
  $B = "$B0".$Bopts[6];

  $file01 = "$dir/socal_vel_dat.dat";

  print CSH "psbasemap $B $R $J -P -K -O -V $shift >> $psfile\n";
  if ($icolor==1) {
    print CSH "awk '{print \$1,\$2,\$4*100}' $file01 | pscontour $R $J -A- -C$cpt_vel -I -P -O -K -V >> $psfile\n";
  }
  print CSH "pscoast $J $R $coast_info_k -P -K -O -V >> $psfile\n";
  #print CSH "psxy $plate_file $J $R $plate_info_k -K -V -O -P >> $psfile \n";
  if($ishelf == 1) {print CSH "awk '{print \$2,\$1}' $shelf_file |psxy $J $R $Wshelf -K -O -P -V >> $psfile\n"}
  print CSH "awk '\$1 == \"R\" {print \$2,\$3}' $recfile |psxy -N $J $R -K -O -P -V $rec_w >> $psfile\n";
  print CSH "awk '{print \$1,\$2}' $evefile |psxy -N $J $R -K -O -P -V $src >> $psfile\n";

  # scale bar
  print CSH "psscale -C$cpt_vel $Dscale $Bscale1b -P -K -O -V >> $psfile \n";
  print CSH "psscale -C$cpt_vel2 $Dscale $Bscale1c -P -K -O -V >> $psfile \n";
  print CSH "pstext -N $J $R -K -O -V -P >>$psfile<<EOF
-115.2 31.42 11 0 $fontno LM \@~\045\@~
-115.2 30.65 11 0 $fontno LM km/s
EOF\n";

  print CSH "pstext -N $J_title $R_title -K -O -V -P >>$psfile<<EOF\n $x_title $z_title $fsize_title 0 $fontno CM $title \nEOF\n";

  #=============================================
  # inset map
  $iwid = 0.8;
  $Jinset = "-JM$iwid";
  $x1q = -130; $x2q = -110;
  $y1q = 25; $y2q = 50;
  $Rinset = "-R$x1q/$x2q/$y1q/$y2q";
  $Binset = "-B1::wesn";
  $origin_inset = "-Xa-0.9 -Ya1.85";
  $coast_info_inset  = "-W0.5p -Dl -Na/0.5p -A1000 -S220/255/255 -C220/255/255 -G200";

  print CSH "pscoast $Jinset $Rinset $coast_info_inset $origin_inset -K -O -V -P >> $psfile \n";
  #print CSH "psxy $plate_file_nuvel $Jinset $Rinset -M -W1.5p/0/0/255 $origin_inset -K -V -O -P >> $psfile \n";
  print CSH "psxy $Jinset $Rinset -W2.0p/255/0/0 -N -K -O -V -P $origin_inset >>$psfile<<EOF\n $xmin $zmin\n$xmax $zmin\n$xmax $zmax\n$xmin $zmax\n$xmin $zmin\nEOF\n";
  #print CSH "psxy $Jinset $Rinset -W2.0p -N -K -O -V $origin_inset >>$psfile<<EOF\n $x1q $y1q\n$x2q $y1q\n$x2q $y2q\n$x1q $y2q\n$x1q $y1q\nEOF\n";

  print CSH "gmtset TICK_LENGTH 0\n";
  print CSH "psbasemap $Jinset $Rinset $Binset $origin_inset -K -O -V -P >> $psfile \n";
  print CSH "gmtset TICK_LENGTH $tick\n";

#-----------------------------
  print CSH "pstext -N $J_title $R_title -O -V -P >>$psfile<<EOF\n 10 10 $fsize_title 0 $fontno CM $title \nEOF\n"; # FINISH
  print CSH "convert $psfile $jpgfile\n";

  close (CSH);
  system("csh -f $cshfile");
  system("xv $jpgfile &");
}

#===========================================================================
if ($ifig03 == 1) {

  #===============================================
  print "\nWriting CSH file...\n";

  open(CSH,">$cshfile");
  print CSH "gmtset BASEMAP_TYPE plain TICK_LENGTH $tick LABEL_FONT_SIZE $fsize2 ANOT_FONT_SIZE $fsize2 PLOT_DEGREE_FORMAT D HEADER_FONT $fontno ANOT_FONT $fontno LABEL_FONT $fontno HEADER_FONT_SIZE $fsize1 FRAME_PEN $fpen TICK_PEN $tpen\n";
  #===============================================

  # get the receivers and sources
  $recfile = "$dir0/${edir}001/sr.txt";  # src-rec for first event
  $evefile = "$dir0/events_lonlat.dat";
  if (not -f $recfile) { die("Check if $recfile exist or not\n") }
  if (not -f $evefile) { die("Check if $evefile exist or not\n") }

  # file names for figures
  $name    = "socal_03";
  $psfile  = "$name.ps";
  $jpgfile = "$name.jpg";

  #=============================================
  # source-receiver geometry

  $B = "$B0".$Bopts[13];
  $title = "Source-receiver geometry";

  $Bscale = "-B2000f500:\"Elevation (km)\": -A";

  print CSH "psbasemap $J $R $B -K -V -P $origin > $psfile\n";  # START
  print CSH "pscoast $J $R $coast_info3 -K -O -V -P  >> $psfile \n";
  print CSH "psbasemap $J $R $B -K -V -P -O >> $psfile\n";
  #print CSH "psxy $plate_file $J $R $plate_info_r -K -V -O -P >> $psfile \n";
  if($ishelf == 1) {print CSH "awk '{print \$2,\$1}' $shelf_file |psxy $J $R $Wshelf -K -O -P -V >> $psfile\n"}
  print CSH "awk '\$1 == \"R\" {print \$2,\$3}' $recfile |psxy -N $J $R -K -O -P -V $rec_r >> $psfile\n";
  print CSH "awk '{print \$1,\$2}' $evefile |psxy -N $J $R -K -O -P -V $src_b >> $psfile\n";

  # plot receivers with numbered label
  #print CSH "awk '\$1 == \"R\" {print \$2,\$3}' $recfile |psxy -N $J $R -K -O -P -V  $rec0 >> $psfile\n";
  #$rec_file = text_rec; $angle = 0; $just = "CM";
  #print CSH "awk '\$1 == \"R\" {print \$2,\$3,$fsize3,$angle,$fontno,\"$just\",\$4}' $recfile > $rec_file\n";
  #print CSH "pstext $rec_file -N $J $R -K -O -P -V >> $psfile\n";

  print CSH "pstext -N $J_title $R_title -K -O -V -P >>$psfile<<EOF\n $x_title $z_title $fsize_title 0 $fontno CM $title \nEOF\n";

  #=============================================
  # phase speed map

  $shift = $shift1;
  $B = "$B0".$Bopts[6];

  $file01 = "$dir/socal_vel_dat.dat";

  print CSH "psbasemap $B $R $J -P -K -O -V $shift >> $psfile\n";
  if ($icolor==1) {
    print CSH "awk '{print \$1,\$2,\$4*100}' $file01 | pscontour $R $J -A- -C$cpt_vel -I -P -O -K -V >> $psfile\n";
  }
  print CSH "pscoast $J $R $coast)info_k -P -K -O -V >> $psfile\n";
  #print CSH "psxy $plate_file $J $R $plate_info_k -K -V -O -P >> $psfile \n";
  if($ishelf == 1) {print CSH "awk '{print \$2,\$1}' $shelf_file |psxy $J $R $Wshelf -K -O -P -V >> $psfile\n"}
  #print CSH "awk '\$1 == \"R\" {print \$2,\$3}' $recfile |psxy -N $J $R -K -O -P -V $rec_w >> $psfile\n";
  #print CSH "awk '{print \$1,\$2}' $evefile |psxy -N $J $R -K -O -P -V $src >> $psfile\n";

  # plot all the possible ray paths
  print "\n $evefile $recfile \n";
  print CSH "awk '{print \$1,\$2}' $evefile > src_temp \n";
  print CSH "awk '\$1 == \"R\" {print \$2,\$3}' $recfile > rec_temp \n";
  open(IN,"rec_temp"); @recs = <IN>; $nrec = @recs;
  open(IN,"src_temp"); @srcs = <IN>; $nsrc = @srcs;
  print "\n nrec = $nrec, nsrc = $nsrc \n";
  if($icolor==1) {
    for ($isrc = 0; $isrc <= $nsrc-1; $isrc = $isrc+1) {
       for ($irec = 0; $irec <= $nrec-1; $irec = $irec+1) {
          chomp($srcs[$isrc]); chomp($recs[$irec]);
          #print "\n $isrc $irec $srcs[$isrc] $recs[$irec] ";
          print CSH "psxy -W0.25p $J $R -K -O -V -P >>$psfile<<EOF\n $srcs[$isrc] \n $recs[$irec] \nEOF\n";
       }
    }
  }  # icolor

  $npath = $nsrc * $nrec;
  $title = "$nsrc events  \@~\264\@~  $nrec receivers  =  $npath ray-paths";

  #die("testing");

  # scale bar
  print CSH "psscale -C$cpt_vel $Dscale $Bscale1b -P -K -O -V >> $psfile \n";
  print CSH "psscale -C$cpt_vel2 $Dscale $Bscale1c -P -K -O -V >> $psfile \n";
  print CSH "pstext -N $J $R -K -O -V -P >>$psfile<<EOF
-115.2 31.42 11 0 $fontno LM \@~\045\@~
-115.2 30.65 11 0 $fontno LM km/s
EOF\n";

  print CSH "pstext -N $J_title $R_title -K -O -V -P >>$psfile<<EOF\n $x_title $z_title $fsize_title 0 $fontno CM $title \nEOF\n";

#  #=============================================
#  # inset map
#  $iwid = 0.8;
#  $Jinset = "-JM$iwid";
#  $x1q = -130; $x2q = -110;
#  $y1q = 25; $y2q = 50;
#  $Rinset = "-R$x1q/$x2q/$y1q/$y2q";
#  $Binset = "-B1::wesn";
#  $origin_inset = "-Xa-0.9 -Ya1.85";

#  print CSH "pscoast $Jinset $Rinset $coast_info_inset $origin_inset -K -O -V >> $psfile \n";
#  print CSH "psxy $Jinset $Rinset -W2.0p/255/0/0 -N -K -O -V $origin_inset >>$psfile<<EOF\n $xmin $zmin\n$xmax $zmin\n$xmax $zmax\n$xmin $zmax\n$xmin $zmin\nEOF\n";
#  #print CSH "psxy $Jinset $Rinset -W2.0p -N -K -O -V $origin_inset >>$psfile<<EOF\n $x1q $y1q\n$x2q $y1q\n$x2q $y2q\n$x1q $y2q\n$x1q $y1q\nEOF\n";

#  print CSH "gmtset TICK_LENGTH 0\n";
#  print CSH "psbasemap $Jinset $Rinset $Binset $origin_inset -K -O -V >> $psfile \n";
#  print CSH "gmtset TICK_LENGTH $tick\n";

#-----------------------------
  print CSH "pstext -N $J_title $R_title -O -V -P >>$psfile<<EOF\n 10 10 $fsize_title 0 $fontno CM $title \nEOF\n"; # FINISH
  print CSH "convert $psfile $jpgfile\n";

  close (CSH);
  system("csh -f $cshfile");
  system("xv $jpgfile &");
}


#===========================================================================
if ($ifig04 == 1) {

  #===============================================
  print "\nWriting CSH file...\n";

  open(CSH,">$cshfile");
  print CSH "gmtset BASEMAP_TYPE plain TICK_LENGTH $tick LABEL_FONT_SIZE $fsize2 ANOT_FONT_SIZE $fsize2 PLOT_DEGREE_FORMAT D HEADER_FONT $fontno ANOT_FONT $fontno LABEL_FONT $fontno HEADER_FONT_SIZE $fsize1 FRAME_PEN $fpen TICK_PEN $tpen\n";
  #===============================================

  # file names for figures
  $name    = "socal_04";
  $psfile  = "$name.ps";
  $jpgfile = "$name.jpg";

  #=============================================
  # source-receiver geometry

  $R = "-R-123/-112/29/36.5";
  $Jwid = 8;
  $J = "-JM${Jwid}i";      # in lat-lon
  $origin = "-X1 -Y1";

  $B = "$B0".$Bopts[0];
  $title = "Source-receiver geometry";

  $Bscale = "-B2000f500:\"Elevation (km)\": -A";

  print CSH "psbasemap $J $R $B -K -V $origin > $psfile\n";  # START
  if ($icolor==1) {print CSH "grdimage $grdfile -C$cpt_topo $J $R -I$gradfile -K -O -V >> $psfile\n"}

  #print CSH "pscoast $J $R $coast_info3 -K -O -V  >> $psfile \n";
  print CSH "pscoast $J $R $coast_info_k -K -O -V  >> $psfile \n";
  print CSH "psbasemap $J $R $B -K -V -O >> $psfile\n";
  print CSH "psxy $plate_file $J $R $plate_info_k -K -V -O >> $psfile \n";
  print CSH "psxy $fault_file $J $R $fault_info_k -K -V -O >> $psfile \n";
  print CSH "psxy $kcf_file $J $R $fault_info_k -K -V -O >> $psfile \n";
  if($ishelf == 1) {print CSH "awk '{print \$2,\$1}' $shelf_file |psxy $J $R -W1.5tap -K -O -V >> $psfile\n"}

#-----------------------------
  print CSH "pstext -N $J_title $R_title -O -V >>$psfile<<EOF\n 10 10 $fsize_title 0 $fontno CM $title \nEOF\n"; # FINISH
  print CSH "convert $psfile -rotate 90 $jpgfile\n";

  close (CSH);
  system("csh -f $cshfile");
  system("xv $jpgfile &");
}

#=================================================================
