#!/usr/bin/perl -w

#==========================================================
#
#  plot_gauss_mean.pl
#  Carl Tape
#  15-June-2010
#
#  Figure to present at Albert Tarantola symposium
#
#    ../plot_gauss_mean.pl
#
#==========================================================

$imodel = 0;
$im = 1;
$irun0 = $im;

# base directory
$pwd = $ENV{PWD};
$dirplot = `dirname $pwd`; chomp($dirplot);
$basedir = `dirname $dirplot`; chomp($basedir);
#$basedir = "/data1/cig/seismo/3D/ADJOINT_TOMO/iterate_adj/SEM2D_iterate_work";

$mfile_syn = "structure_syn.dat";

# directories
$plotdir = "${basedir}/PLOTTING";
$odir = "${basedir}/OUTPUT";
$idir = "${basedir}/INPUT";
$figdir = "${plotdir}/FIGURES";
if (not -e $figdir) {die("Check if figdir $figdir exist or not\n");}

$stimodel = sprintf("m%.1f",$imodel/2);
$stirun0 = sprintf("%4.4i",$irun0);
$stirun = sprintf("%4.4i",$irun0+$imodel);

$dir0 = "$odir/run_${stirun0}";

if(1==1) {
   $dir1 = "$odir/run_${stirun}";
} else {
   $dir1 = sprintf("$dir0/READ_IN_CG/model_m%4.4i",$imodel);
}

$file1syn = "$dir1/${mfile_syn}";
if (not -f $file1syn) { die("Check if $file1syn exist or not\n") }
#if (not -f $evefile) { die("Check if $evefile exist or not\n") }

# read in reference values: alpha0, beta0, per
if (1==1) {
  $fileref1 = "$dir0/reference_values.dat";
  $fileref2 = "$dir0/reference_period.dat";
  open(IN,$fileref1); @lines = <IN>; ($alpha0,$beta0) = split(" ",$lines[0]);
  open(IN,$fileref2); $per = <IN>; chomp($per);
} else {
  $per = 20;
  $beta0 = 3500;
}
$per_lab = sprintf("T = %3.3f s",$per);
$beta0_lab  = sprintf("beta0 = %3.3f km/s",$beta0/1000);
#print "\n-- $per -- $per_lab -- $beta0_lab --\n"; die("testing");

#$edir  = sprintf("$dir/event_%3.3i",$event);
#print "\n $dir,$edir,$mfile_dat,$mfile_syn,$beta0,$per,$ifinite,$iopt \n";

$plabel = "${plotdir}/plot_gauss_mean.pl";

#die("\ntesting\n");

# plotting specifications
$fsize0 = "14";
$fsize1 = "11";
$fsize2 = "9";
$fsize3 = "5";
$fontno = "1";    # 1 or 4
$tick   = "0.1c";

$tinc = 50;   # tick increment (in seconds) for source time function

# plot symbols for sources, receivers, and shelf
$src0 = "-W0.75p -Sa0.20 -G255/0/0";
$src = "-W0.75p -Sa0.20";
$rec = "-W0.5p/0/0/0 -St0.10";
$rec0 = "-Sc10p -W0.5p";
$Wshelf = "-W1.0/0/0/0tap";

# resolution of color plots
$interp = "-I0.5m/0.5m -S4m";   # key information
$grdfile = "temp.grd";

#-------------------------
# color

# KEY: scaling for color
$scale_color = 21.0;
$colorbar = "seis";

#-------------------------

# write plotting scripts
$wid = 2.5;
$J1 = "-JM${wid}i";      # in lat-lon
$origin = "-X1.5 -Y2.0";

$B1dat = "-B1:.\" \":WeSN";
$B1syn = "-B1:.\" \":wESN";

$Dlen = 2.0;
$Dx = $wid/2;
$Dy = -0.35;

$Dscale1 = "-D$Dx/$Dy/$Dlen/0.10h";

$ptick = 0.05;
#$Bscale1  = sprintf("-B%2.2e:\" Phase Velocity ( km s\@+-1\@+ )\": -E10p",$ptick);
$Bscale1  = sprintf("-B%2.2e:\" Perturbation from %3.3f  km s\@+-1\@+ \": -E10p",$ptick,$beta0/1000);

#-------------------------

$title = "Membrane Wave Speed";
#$title = "Rayleigh Wave Phase Velocity";

# set bounds for the plotting
#$xmin = -121; $xmax = -115.0; $zmin = 31.75; $zmax = 36.5;
($xmin,$xmax,$zmin,$zmax,$smin,$smax,$tmin,$tmax) = split(" ",`minmax -C $file1syn`);
#($tt) = sort{$b <=> $a} ($tt,abs($tmin),abs($tmax));
$dinc = 0.25;  # buffer around min/max
$xmin = $xmin-$dinc;  $xmax = $xmax+$dinc;
$zmin = $zmin-$dinc;  $zmax = $zmax+$dinc;
$R = "-R$xmin/$xmax/$zmin/$zmax";

print "\n $R \n $smin $smax \n";

# plot title
$xtx = $xmin+0.5*($xmax-$xmin);
$ztx = $zmin+1.10*($zmax-$zmin);

  #===============================================
  $cshfile = "plot_gauss_mean.csh";
  print "\nWriting to CSH file ${cshfile}...\n";
  open(CSH,">$cshfile");
  print CSH "gmtset BASEMAP_TYPE plain PAPER_MEDIA letter TICK_LENGTH $tick LABEL_FONT_SIZE $fsize2 ANOT_FONT_SIZE $fsize2 PLOT_DEGREE_FORMAT D HEADER_FONT $fontno ANOT_FONT $fontno LABEL_FONT $fontno HEADER_FONT_SIZE $fsize1\n";
  #===============================================

  # make colorpoint file
  $pmax = 0.1;
  $dc = $pmax/10;
  $T1 = "-T-$pmax/$pmax/$dc";
  $cptfile1 = "color1.cpt";
  print CSH "makecpt -C$colorbar -D $T1 > $cptfile1\n";

  $pmax = 0.1/10;
  $dc = $pmax/10;
  $T1 = "-T-$pmax/$pmax/$dc";
  $cptfile2 = "color2.cpt";
  print CSH "makecpt -C$colorbar -D $T1 > $cptfile2\n";

  $pmax = 0.0025;
  $dc = $pmax/20;
  $T1 = "-T0/$pmax/$dc";
  $cptfile3 = "color3.cpt";
  print CSH "makecpt -Chot -D -I $T1 > $cptfile3\n";

  $Cfac = 10;
  $pmax = 0.0025/$Cfac;
  $dc = $pmax/20;
  $T1 = "-T0/$pmax/$dc";
  $cptfile4 = "color4.cpt";
  print CSH "makecpt -Chot -D -I $T1 > $cptfile4\n";

#=============================================

  $origin = "-X1 -Y7.5";
  $B1syn = "-B1:.\" \":WESN";

#=============================================

$idir = "/data1/cig/seismo/3D/ADJOINT_TOMO/iterate_adj/SEM2D_iterate_INPUT/random_fields/model_set_0001";

$name    = "gaussian_mean";
$psfile  = "$name.ps";
$jpgfile = "$name.jpg";

  $B1syn = "-B1:.\" \":WesN";
  $file1syn = "$idir/mean.dat";
  if (not -f $file1syn) {die("Check if file1syn $file1syn exist or not\n");}
   $olab = "-Xa1.3 -Ya-0.15";

  # mean model
  print CSH "psbasemap $B1syn $R $J1 -P -K -V $origin > $psfile\n";  # START
  print CSH "awk '{print \$1,\$2,\$6}' $file1syn | nearneighbor -G$grdfile $R $interp\n";
  print CSH "grdimage $grdfile -C$cptfile1 $J1 -K -O -V -Q >> $psfile\n";
  print CSH "pscoast $J1 $R $B1syn -W1p -Na/1p -Dh -K -O -V >> $psfile\n";
  #print CSH "psscale -C$cptfile1 $Dscale1 $Bscale1 -K -O -V >> $psfile \n";
  $faultfile = "/home/carltape/gmt/faults/jennings_more.xy";
  if (not -f $faultfile) {die("Check if faultfile $faultfile exist or not\n");}
  print CSH "psxy $faultfile $J1 $R -W0.75p -M -K -V -O >> $psfile\n";
  print CSH "pstext -N $J1 $R -K -O -V $olab >>$psfile<<EOF\n $xmin $zmin $fsize0 0 $fontno CM Mean Model\nEOF\n";

  $B1syn = "-B1:.\" \":Wesn";
  $file1syn = "$idir/mean_est.dat";
  if (not -f $file1syn) {die("Check if file1syn $file1syn exist or not\n");}
  $dY0 = 3.2;
  $shift = "-Y-$dY0";

  # estimated mean model
  print CSH "psbasemap $B1syn $R $J1 -K -O -V $shift >> $psfile\n";
  print CSH "awk '{print \$1,\$2,\$6}' $file1syn | nearneighbor -G$grdfile $R $interp\n";
  print CSH "grdimage $grdfile -C$cptfile1 $J1 -K -O -V -Q >> $psfile\n";
  print CSH "pscoast $J1 $R $B1syn -W1p -Na/1p -Dh -K -O -V >> $psfile\n";
  #print CSH "psscale -C$cptfile1 $Dscale1 $Bscale1 -K -O -V >> $psfile \n";
  print CSH "psxy $faultfile $J1 $R -W0.75p -M -K -V -O >> $psfile\n";
  print CSH "pstext -N $J1 $R -K -O -V $olab >>$psfile<<EOF\n $xmin $zmin $fsize0 0 $fontno CM Estimated Mean Model\nEOF\n";

  $B1syn = "-B1:.\" \":Wesn";

  # estimated mean model -- zoomed
  print CSH "psbasemap $B1syn $R $J1 -K -O -V $shift >> $psfile\n";
  print CSH "awk '{print \$1,\$2,\$6}' $file1syn | nearneighbor -G$grdfile $R $interp\n";
  print CSH "grdimage $grdfile -C$cptfile2 $J1 -K -O -V -Q >> $psfile\n";
  print CSH "pscoast $J1 $R $B1syn -W1p -Na/1p -Dh -K -O -V >> $psfile\n";
  #print CSH "psscale -C$cptfile1 $Dscale1 $Bscale1 -K -O -V >> $psfile \n";
  print CSH "psxy $faultfile $J1 $R -W0.75p -M -K -V -O >> $psfile\n";
  print CSH "pstext -N $J1 $R -K -O -V $olab >>$psfile<<EOF\n $xmin $zmin $fsize0 0 $fontno CM Estimated Mean Model (x10)\nEOF\n";

#------------------------------------------------------------------------
# covariance matrices -- NOTE plotting j,i,Cij

  $file1syn = "$idir/C_Cest.dat";
  if (not -f $file1syn) {die("Check if file1syn $file1syn exist or not\n");}
  $B1syn = "-Ba256f32:.\" \":WesN";
  $J1 = "-JX${wid}i/-${wid}i";
  $dY2 = 2*$dY0;
  $shift = "-X3.5i -Y${dY2}i";
  $xmin = 0; $xmax = 1025; $zmin = 0; $zmax = 1025;
  $R = "-R$xmin/$xmax/$zmin/$zmax";
  $interp = "-I2/2 -S2";   # key information

  print CSH "psbasemap $B1syn $R $J1 -K -O -V $shift >> $psfile\n";
  #print CSH "awk '{print \$1,\$2,\$3}' $file1syn | psxy -C$cptfile3 -Sc1p $R $J1 -K -O -V >> $psfile\n";
  print CSH "awk '{print \$2,\$1,\$3}' $file1syn | nearneighbor -G$grdfile $R $interp\n";
  print CSH "grdimage $grdfile -C$cptfile3 $J1 -K -O -V -Q >> $psfile\n";
  print CSH "pstext -N $J1 $R -K -O -V $olab >>$psfile<<EOF\n $xmin $zmax $fsize0 0 $fontno CM Covariance Matrix\nEOF\n";
  print CSH "psbasemap $B1syn $R $J1 -K -O -V >> $psfile\n";

  $shift = "-Y-$dY0";

  print CSH "psbasemap $B1syn $R $J1 -K -O -V $shift >> $psfile\n";
  #print CSH "awk '{print \$1,\$2,\$3}' $file1syn | psxy -C$cptfile3 -Sc1p $R $J1 -K -O -V >> $psfile\n";
  print CSH "awk '{print \$2,\$1,\$4}' $file1syn | nearneighbor -G$grdfile $R $interp\n";
  print CSH "grdimage $grdfile -C$cptfile3 $J1 -K -O -V -Q >> $psfile\n";
  print CSH "pstext -N $J1 $R -K -O -V $olab >>$psfile<<EOF\n $xmin $zmax $fsize0 0 $fontno CM Maximum Likelihood Matrix\nEOF\n";
  print CSH "psbasemap $B1syn $R $J1 -K -O -V >> $psfile\n";

  print CSH "psbasemap $B1syn $R $J1 -K -O -V $shift >> $psfile\n";
  #print CSH "awk '{print \$1,\$2,\$3}' $file1syn | psxy -C$cptfile4 -Sc1p $R $J1 -K -O -V >> $psfile\n";
  print CSH "awk '{print \$2,\$1,\$4-\$3}' $file1syn | nearneighbor -G$grdfile $R $interp\n";
  print CSH "grdimage $grdfile -C$cptfile4 $J1 -K -O -V -Q >> $psfile\n";
  print CSH "pstext -N $J1 $R -K -O -V $olab >>$psfile<<EOF\n $xmin $zmax $fsize0 0 $fontno CM Residual Matrix (x$Cfac)\nEOF\n";
  print CSH "psbasemap $B1syn $R $J1 -K -O -V >> $psfile\n";


#-----------------------------
#  print CSH "convert $psfile $jpgfile\n";

  print CSH "pstext -N $J1 $R -O -V -Xa-0 -Ya-0.65 >>$psfile<<EOF\n $xmin $zmin $fsize0 0 $fontno CM \nEOF\n";  # FINISH

  close (CSH);
  system("csh -f $cshfile");
  system("gv $psfile &")
  #if($iopt <= 2) {system("xv $jpgfile &")}

#=================================================================
