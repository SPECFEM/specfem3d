#!/usr/bin/perl -w

#==========================================================
#
#  plot_surf_model.pl
#  Carl Tape
#  15-Jan-2010
#
#  Plot the source time function, the velocity model for the synthetics,
#  and the velocity model for the data (optional).
#  This program can executed from wave2d.f90.
#
#  iopt = 1 : reference model
#  iopt = 2 : data model, reference model
#  iopt = 3 : reference model with source + receivers, source time function
#  iopt = 4 : data model and reference model with source + receivers, source time function
#
#  Examples:
#    ../plot_surf_model.pl 4 100 0 0/1/5 0.1/0.05 stf_00001_1 0
#    ../plot_surf_model.pl 4 200 0 0/1/5 0.1/0.05 stf_00001_1 0
#    ../plot_surf_model.pl 4 300 0 0/1/5 0.1/0.05 stf_00001_1 0
#    ../plot_surf_model.pl 4 400 0 1/0/5 0.1/0.05 stf_00001_1 0
#    ../plot_surf_model.pl 4 500 0 1/0/5 0.1/0.05 stf_00001_1 0
#    ../plot_surf_model.pl 4 600 0 1/0/5 0.1/0.05 stf_00001_1 0
#    ../plot_surf_model.pl 4 700 0 1/0/5 0.1/0.05 stf_00001_1 0
#    ../plot_surf_model.pl 4 800 0 0/1/5 0.1/0.05 stf_00001_1 0
#
#    ../plot_surf_model.pl 4 5000 0 1/0/5 0.1/0.05 stf_00001_1 0
#    ../plot_surf_model.pl 4 5000 2 1/0/5 0.1/0.05 stf_00001_1 0
#
#  Set of 1000 models
#    ../plot_surf_model.pl 1 0001 0 1/0/1 0.1/0.05 stf_00001_1 0
#
#==========================================================

if (@ARGV < 7) {die("Usage: plot_surf_model.pl xxx\n");}
($iopt,$irun0,$imodel,$ievents,$pbar,$stf_file1,$ifinite) = @ARGV;
$comp = 1;

# base directory
$pwd = $ENV{PWD};
$dirplot = `dirname $pwd`; chomp($dirplot);
$basedir = `dirname $dirplot`; chomp($basedir);
#$basedir = "/data1/cig/seismo/3D/ADJOINT_TOMO/iterate_adj/SEM2D_iterate_work";

($pmax,$ptick) = split("/",$pbar);
($ievent_all,$ievent_one,$ievent) = split("/",$ievents);
$rec_file1 = "sr.txt";
$efile_syn = "events_syn_lonlat.dat";
$mfile_syn = "structure_syn.dat";
$mfile_dat = "structure_dat.dat";

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

$dir2 = sprintf("$dir1/event_%3.3i",$ievent);
$stf_file = "$dir2/${stf_file1}";
$rec_file = "$dir2/${rec_file1}";
$file1dat = "$dir1/${mfile_dat}";
$file1syn = "$dir1/${mfile_syn}";
$evefile = "$dir0/${efile_syn}";     # event file from the base directory
if (not -f $file1dat) { die("Check if $file1dat exist or not\n") }
if (not -f $file1syn) { die("Check if $file1syn exist or not\n") }
if (not -f $evefile) { die("Check if $evefile exist or not\n") }

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

$plabel = "${plotdir}/plot_surf_model.pl";

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
if($ifinite==0) {$src0 = "-W0.75p -Sa0.20 -G255/0/0"} else {$src0 = "-Sc0.05"}
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
$wid = 3.5;
$J1 = "-JM${wid}i";      # in lat-lon
$origin = "-X1.5 -Y2.0";

$B1dat = "-B1:.\" \":WeSN";
$B1syn = "-B1:.\" \":wESN";

$Dlen = 2.0;
$Dx = $wid/2;
$Dy = -0.35;

$Dscale1 = "-D$Dx/$Dy/$Dlen/0.10h";

#$Bscale1  = sprintf("-B%2.2e:\" Phase Velocity ( km s\@+-1\@+ )\": -E10p",$ptick);
$Bscale1  = sprintf("-B%2.2e:\" Perturbation from %3.3f  km s\@+-1\@+ \": -E10p",$ptick,$beta0/1000);

#-------------------------

$title = "Membrane Wave Speed";
#$title = "Rayleigh Wave Phase Velocity";

# set bounds for the plotting
#$xmin = -121; $xmax = -115.0; $zmin = 31.75; $zmax = 36.5;
($xmin,$xmax,$zmin,$zmax,$smin,$smax,$tmin,$tmax) = split(" ",`minmax -C $file1dat`);
($tt) = sort{$b <=> $a} ($tt,abs($tmin),abs($tmax));
$dinc = 0.25;  # buffer around min/max
$xmin = $xmin-$dinc;  $xmax = $xmax+$dinc;
$zmin = $zmin-$dinc;  $zmax = $zmax+$dinc;
$R = "-R$xmin/$xmax/$zmin/$zmax";

print "\n $R \n $smin $smax \n -$tt $tt \n";

# plot title
$xtx = $xmin+0.5*($xmax-$xmin);
$ztx = $zmin+1.10*($zmax-$zmin);

$name    = "model_${stirun}";
$psfile  = "$name.ps";
$jpgfile = "$name.jpg";

  #===============================================
  $cshfile = "plot_surf_model.csh";
  print "\nWriting to CSH file ${cshfile}...\n";
  open(CSH,">$cshfile");
  print CSH "gmtset BASEMAP_TYPE plain PAPER_MEDIA letter TICK_LENGTH $tick LABEL_FONT_SIZE $fsize2 ANOT_FONT_SIZE $fsize2 PLOT_DEGREE_FORMAT D HEADER_FONT $fontno ANOT_FONT $fontno LABEL_FONT $fontno HEADER_FONT_SIZE $fsize1\n";
  #===============================================

  # make colorpoint file
#  $T1 = "-T3/4/0.1";
#  $pmax = 10;
  $dc = $pmax/10;
  $T1 = "-T-$pmax/$pmax/$dc";
  $cptfile1 = "color0.cpt";
  print CSH "makecpt -C$colorbar -D $T1 > $cptfile1\n";
  #print CSH "makecpt -C$colorbar -D $T1 > temp1\n";
  #print CSH "sed 's/^B.*/B       170     0      0  /' temp1  >  temp2\n";
  #print CSH "sed 's/^F.*/F         0     0    255  /' temp2 > $cptfile1\n";

#=============================================
if($iopt == 1) {
#=============================================

  $plot_title = "Model for synthetics";
  $origin = "-X1.5 -Y5.0";
  $B1syn = "-B1:.\" \":WESN";

  # phase velocity map
  print CSH "psbasemap $B1syn $R $J1 -P -K -V $origin > $psfile\n";  # START
  #print CSH "awk '{print \$1,\$2,\$6}' $file1syn | pscontour $R $J1 -A- -C$cptfile1 -I -K -O -V >> $psfile\n";
  print CSH "awk '{print \$1,\$2,\$6}' $file1syn | nearneighbor -G$grdfile $R $interp\n";
  print CSH "grdimage $grdfile -C$cptfile1 $J1 -K -O -V -Q >> $psfile\n";
  print CSH "pscoast $J1 $R $B1syn -W1p -Na/1p -Dh -K -O -V >> $psfile\n";
  #print CSH "awk '{print \$2,\$1}' $idir/BOUNDARIES/oms_shelf |psxy $J1 $R $Wshelf -K -O -V >> $psfile\n";
  print CSH "psscale -C$cptfile1 $Dscale1 $Bscale1 -K -O -V >> $psfile \n";

  if (1==1) {
    $faultfile = "/home/carltape/gmt/faults/jennings_more.xy";
    if (not -f $faultfile) {die("Check if faultfile $faultfile exist or not\n");}
    print CSH "psxy $faultfile $J1 $R -W0.75p -M -K -V -O >> $psfile\n";
  }

  # plot title and GMT header
  $Utag = "-U/0/0.25/$plabel";
  $shift = "-Xa0 -Ya4";
  print CSH "pstext -N $J1 $R $Utag -O -V $shift >>$psfile<<EOF\n $xmin $zmin $fsize0 0 $fontno LM $plot_title\nEOF\n";  # FINISH

#=============================================
} if($iopt == 2) {
#=============================================

  $plot_title = "Model for data, model for synthetics";
  $origin = "-X0.6 -Y5.0";

  # phase velocity map -- data
  print CSH "psbasemap $B1dat $R $J1 -P -K -V $origin > $psfile\n";  # START
  #print CSH "awk '{print \$1,\$2,\$6}' $file1dat | pscontour $R $J1 -A- -C$cptfile1 -I -V -K -O >> $psfile\n";
  print CSH "awk '{print \$1,\$2,\$6}' $file1dat | nearneighbor -G$grdfile $R $interp\n";
  print CSH "grdimage $grdfile -C$cptfile1 $J1 -K -O -V -Q >> $psfile\n";
  print CSH "pscoast $J1 $R $B1dat -W1p -Na/1p -Dh -V -K -O >> $psfile\n";
  #print CSH "awk '{print \$2,\$1}' $idir/BOUNDARIES/oms_shelf |psxy $J1 $R $Wshelf -V -K -O >> $psfile\n";
  print CSH "psscale -C$cptfile1 $Dscale1 $Bscale1 -V -K -O >> $psfile \n";

  $dX = 1.1*$wid; $dY = 0;
  #$dX = 0; $dY = 4.0;
  $shift = "-X$dX -Y$dY";

  # phase velocity map -- synthetics
  print CSH "psbasemap $B1syn $R $J1 -K -O -V $shift >> $psfile\n";
  #print CSH "awk '{print \$1,\$2,\$6}' $file1syn | pscontour $R $J1 -A- -C$cptfile1 -I -V -K -O >> $psfile\n";
  print CSH "awk '{print \$1,\$2,\$6}' $file1syn | nearneighbor -G$grdfile $R $interp\n";
  print CSH "grdimage $grdfile -C$cptfile1 $J1 -K -O -V -Q >> $psfile\n";
  print CSH "pscoast $J1 $R $B1syn -W1p -Na/1p -Dh -V -K -O >> $psfile\n";
  print CSH "psscale -C$cptfile1 $Dscale1 $Bscale1 -V -K -O >> $psfile \n";

  # plot title and GMT header
  $Utag = "-U/0/0.25/$plabel";
  $shift = "-Xa-$dX -Ya4";
  print CSH "pstext -N $J1 $R $Utag -V -O $shift >>$psfile<<EOF\n $xmin $zmin $fsize0 0 $fontno LM $plot_title\nEOF\n";  # FINISH

#=============================================
} if($iopt > 2) {
#=============================================

  #$rec_file = "$edir/sr.txt";
  if (not -f $rec_file) { die("Check if rec_file $rec_file exist or not\n") }

  # source time function
  #$stf_file = "$edir/stf_00001_$comp";
  print "\n Source time function is $stf_file\n";
  ($tmin,$tmax,$smin,$smax) = split(" ",`minmax -C $stf_file`);
  ($ss) = sort{$b <=> $a} (abs($smin),abs($smax));

  if ($ss==0) {$ss = 1}  # do not normalize by zero
  print "smax = $ss \n tmin = $tmin tmax = $tmax\n";

  #$stf_data = "temp";
  $Bs = sprintf("-B%3.3f:\"Time  (s)\":/%3.3f:\" h (t) \"::.\"Source time function (hdur = %3.1f s)\":WeSn",$tinc,0.5,$per/2);
  $Rs = "-R$tmin/$tmax/-1/1";
  $Js = "-JX5/1.5";
  #print CSH "awk '{print \$1,\$2/$ss}' $stf_file > $stf_data\n";
  #print CSH "psxy $stf_data $Rs $Js $Bs -W1p -K -O -V -Y-3.5 >> $psfile\n";

  print CSH "psbasemap $Rs $Js $Bs -P -K -V $origin > $psfile\n"; # START
  print CSH "awk '{print \$1,\$2/$ss}' $stf_file | psxy $Rs $Js $Bs -W1p -K -O -V >> $psfile\n";

if($iopt == 4) {

  # model for data (optional)
  $plot_title = "Target model, model $stimodel, and source time function (run_${stirun})";
  $shift = "-Y4.0 -X-0.9";

  # phase velocity map
  print CSH "psbasemap $B1dat $R $J1 -K -O -V $shift >> $psfile\n";
  #print CSH "awk '{print \$1,\$2,\$6}' $file1dat | pscontour $R $J1 -A- -C$cptfile1 -I -K -O -V >> $psfile\n";
  print CSH "awk '{print \$1,\$2,\$6}' $file1dat | nearneighbor -G$grdfile $R $interp\n";
  print CSH "grdimage $grdfile -C$cptfile1 $J1 -K -O -V -Q >> $psfile\n";
  print CSH "pscoast $J1 $R $B1dat -W1p -Na/1p -Dh -K -O -V >> $psfile\n";

  # plot receivers with numbered label
  print CSH "awk '\$1 == \"R\" {print \$2,\$3}' $rec_file | psxy -N $J1 $R -K -O -V $rec0 >> $psfile\n";
  $rec_file2 = text_rec; $angle = 0; $just = "CM";
  print CSH "awk '\$1 == \"R\" {print \$2,\$3,$fsize3,$angle,$fontno,\"$just\",\$4}' $rec_file > $rec_file2\n";
  print CSH "pstext $rec_file2 -N $J1 $R -K -O -V >> $psfile\n";

  if($ievent_one) {print CSH "awk '\$1 == \"S\" {print \$2,\$3}' $rec_file | psxy $src0 -N $J1 $R -K -O -V >> $psfile\n";}
  if($ievent_all) {print CSH "awk '{print \$1,\$2}' $evefile | psxy $src -N $J1 $R -K -O -V >> $psfile\n";}
  print CSH "psscale -C$cptfile1 $Dscale1 $Bscale1 -K -O -V >> $psfile \n";

  $dX = 1.1*$wid; $dY = 0;

} else {

  $plot_title = "Velocity model and source time function";
  $dX = 0; $dY = 4.0;
}

  # model for synthetics
  $shift = "-X$dX -Y$dY";

  # phase velocity map
  print CSH "psbasemap $B1syn  $R $J1 -K -O -V $shift >> $psfile\n";
  #print CSH "awk '{print \$1,\$2,\$6}' $file1syn | pscontour $R $J1 -A- -C$cptfile1 $shift -I -O -K -V >> $psfile\n";
  print CSH "awk '{print \$1,\$2,\$6}' $file1syn | nearneighbor -G$grdfile $R $interp\n";
  print CSH "grdimage $grdfile -C$cptfile1 $J1 -K -O -V -Q >> $psfile\n";
  print CSH "pscoast $J1 $R $B1syn -W1p -Na/1p -Dh -K -O -V >> $psfile\n";

  # plot receivers with numbered label
  print CSH "awk '\$1 == \"R\" {print \$2,\$3}' $rec_file | psxy -N $J1 $R -K -O -V $rec0 >> $psfile\n";
  $rec_file2 = text_rec; $angle = 0; $just = "CM";
  print CSH "awk '\$1 == \"R\" {print \$2,\$3,$fsize3,$angle,$fontno,\"$just\",\$4}' $rec_file > $rec_file2\n";
  print CSH "pstext $rec_file2 -N $J1 $R -K -O -V >> $psfile\n";
  if($ievent_one) {print CSH "awk '\$1 == \"S\" {print \$2,\$3}' $rec_file | psxy $src0 -N $J1 $R -K -O -V >> $psfile\n";}
  if($ievent_all) {print CSH "awk '{print \$1,\$2}' $evefile | psxy $src -N $J1 $R -K -O -V >> $psfile\n";}
  print CSH "psscale -C$cptfile1 $Dscale1 $Bscale1 -K -O -V >> $psfile \n";

  # plot title and GMT header
  $Utag = "-U/0/0.25/$plabel";
  $shift = "-Xa-$dX -Ya4";
  print CSH "pstext -N $J1 $R $Utag -O -V $shift >>$psfile<<EOF\n $xmin $zmin $fsize0 0 $fontno LM $plot_title\nEOF\n";  # FINISH

}  # $iopt > 2

#-----------------------------
  print CSH "convert $psfile $jpgfile\n";

  close (CSH);
  system("csh -f $cshfile");
  system("gv $psfile &")
  #if($iopt <= 2) {system("xv $jpgfile &")}

#=================================================================
