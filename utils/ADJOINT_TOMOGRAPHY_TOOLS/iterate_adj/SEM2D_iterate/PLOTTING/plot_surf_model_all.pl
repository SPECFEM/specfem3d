#!/usr/bin/perl -w

#==========================================================
#
#  plot_surf_model_all.pl
#  Carl Tape
#  15-Jan-2010
#
#  Full display of models and model perturbations.
#  FIGURE 1: 6-panel figure
#  FIGURE 2: 9-panel figure (emphasizing CG test model)
#  Input allows for a range of models as well.
#
#  Example:
#    ../plot_surf_model_all.pl 800 0/32 0/1/5 0.1/0.05 0
#    ../plot_surf_model_all.pl 7200 64/64 1/0/5 0.1/0.05 0
#
#==========================================================

if (@ARGV < 5) {die("Usage: plot_surf_model_all.pl xxx\n");}
($irun0,$imodels,$ievents,$pbar,$ifinite) = @ARGV;
$comp = 1;

# USER INPUT
$ifig0 = 0;
$ifig1 = 0;
$ifig2 = 1;
$itest = 0;      # plot CG test models or not

$icolor = 1;
$ijpg = 0;
$ipdf = 1;

# base directory
$pwd = $ENV{PWD};
$dirplot = `dirname $pwd`; chomp($dirplot);
$basedir = `dirname $dirplot`; chomp($basedir);
#$basedir = "/data1/cig/seismo/3D/ADJOINT_TOMO/iterate_adj/SEM2D_iterate_work";

($imodelmin,$imodelmax) = split("/",$imodels);
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

$stirun0 = sprintf("%4.4i",$irun0);

$dir0 = "$odir/run_${stirun0}";
$dir2 = sprintf("$dir0/event_%3.3i",$ievent);
$rec_file = "$dir2/${rec_file1}";
$evefile = "$dir0/${efile_syn}";     # event file from the base directory
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

$plabel = "${plotdir}/plot_surf_model_all.pl";

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
$interp = "-I2m/2m -S4m";   # key information
$grdfile = "temp.grd";

#-------------------------
# color

# KEY: scaling for color
$scale_color = 21.0;
$colorbar = "seis";

#-------------------------

$wid1 = 2.5;
$wid2 = 2.0;

$B1dat = "-B1:.\" \":WeSN";
$B1syn = "-B1:.\" \":wESN";

#===============================================
$cshfile = "plot_surf_model_all.csh";
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

#$rec_file = "$edir/sr.txt";
if (not -f $rec_file) {die("Check if rec_file $rec_file exist or not\n");}

# sides to show the lat-lon tick marks
$B0 = "-B1:.\" \":";
@Bopts1 = ("WesN","wEsN","Wesn","wEsn","Wesn","wEsn");
@Bopts2 = ("WesN","wesN","wEsN","Wesn","wesn","wEsn","Wesn","wesn","wEsn");

# figure 1
$dX = 1.1*$wid1; $dY = 0; $shift1 = "-X$dX -Y$dY";
$dX = -1.1*$wid1; $dY = -1.1*$wid1; $shift2 = "-X$dX -Y$dY";
@shifts1 = ($origin,$shift1,$shift2,$shift1,$shift2,$shift1);

# figure 2
$dX = 1.1*$wid2; $dY = 0; $shift1 = "-X$dX -Y$dY";
$dX = -2.2*$wid2; $dY = -1.1*$wid2; $shift2 = "-X$dX -Y$dY";
@shifts2 = ($origin,$shift1,$shift1,$shift2,$shift1,$shift1,$shift2,$shift1,$shift1);

$fcol = 6;
#@finds = (7,7,7,7,7,7,7,7,7);

$J_title = "-JM1";
$R_title = "-R0/1/0/1";

#================

# loop over models
for ($m = $imodelmin; $m <= $imodelmax; $m = $m+1) {

  $imodel = $m;     # current model
  $mtest = $m % 2;    # =1 for a test model

  $stirun = sprintf("%4.4i",$irun0+$imodel);
  if ($mtest==1) {
    $stmod = sprintf("m%2.2it",$imodel/2);
  } else {
    $stmod = sprintf("m%2.2i",$imodel/2);
  }

  # target model
  $filedat = "$dir0/${mfile_dat}";

  # directories (current model and previos model)
  $stirun  = sprintf("%4.4i",$irun0+$imodel);
  $dir1    = "$odir/run_${stirun}";

  #-----------------------------
  # FIGURE 0

  if ($ifig0 == 1) {

    $name = "model_all_${stirun}_f0";
    $psfile0 = "$name.ps";
    $jpgfile = "$name.jpg";

    # write plotting scripts
    $wid = $wid1;
    $J1 = "-JM${wid}i";   # in lat-lon
    $origin = "-X1.5 -Y7.25";
    $Dlen = 2.0; $Dx = $wid/2; $Dy = -0.35;
    $Dscale1 = "-D$Dx/$Dy/$Dlen/0.15h";
    $Bscale1  = sprintf("-B%2.2e:\" Perturbation from %3.3f  km s\@+-1\@+ \": -E10p",$ptick,$beta0/1000);
    $Bscale2  = sprintf("-B%2.2e:\" Perturbation from mtar\": -E10p",$ptick);
    $xlab = $wid/2;
    $ylab = -0.02*$wid;
    $olab = "-Xa$xlab -Ya$ylab";

    print "\n----FIGURE 1--model $stmod----\n";

    @files = ("$dir0/${mfile_syn}","$dir0/${mfile_syn}",
        "$dir1/${mfile_syn}","$dir1/${mfile_syn}",
        "$dir0/${mfile_dat}","$dir0/${mfile_dat}");
    @stlabs = ("m00","mtar - m00",$stmod,"mtar - $stmod","mtar","mtar - mtar");

    # only plot even models (m04 = iteration 8)
    $kmaxt = 6;
    $kmaxm = 6*(1 - $mtest);  # =0 for test model
    if ($itest == 1) {
      $kmax = $kmaxt;
    } else {
      $kmax = $kmaxm;
    }

    # loop over subplots
    for ($k = 0; $k < $kmax; $k = $k+1) {

      $iB = $Bopts1[$k];
      $shift = $shifts1[$k];
      $B = "${B0}${iB}";
      $find = $fcol;
      $find2 = 2*$fcol;

      $file1 = $files[$k];
      if (not -f $file1) {
  die("Check if $file1 exist or not\n");
      }

      # column-concatenate target model and current model
      if ($k % 2 == 1) {
  $file2 = "ftemp";
  print CSH "paste $filedat $file1 > $file2\n";
      }

      # phase velocity map
      if ($k==0 && $m==$imodelmin) {
  # set bounds for the plotting
  #$xmin = -121; $xmax = -115.0; $zmin = 31.75; $zmax = 36.5;
  ($xmin,$xmax,$zmin,$zmax,$smin,$smax,$tmin,$tmax) = split(" ",`minmax -C $file1`);
  $dinc = 0.25;   # buffer around min/max
  $xmin = $xmin-$dinc;  $xmax = $xmax+$dinc;
  $zmin = $zmin-$dinc;  $zmax = $zmax+$dinc;
  $R = "-R$xmin/$xmax/$zmin/$zmax";
  print "\n$R\n";
      }

      if ($k==0) {
  print CSH "psbasemap $B $R $J1 -K -V -P $origin > $psfile0\n"; # START
      } else {
  print CSH "psbasemap $B $R $J1 -K -O -V $shift >> $psfile0\n";
      }

      if ($icolor == 1) {
  #print CSH "awk '{print \$1,\$2,\$7}' $file1 | pscontour $R $J1 -A- -C$cptfile1 -I -K -O -V >> $psfile0\n";
  if ($k % 2 == 0) {
    print CSH "awk '{print \$1,\$2,\$${find}}' $file1 | nearneighbor -G$grdfile $R $interp\n";
  } else {
    print CSH "awk '{print \$1,\$2,\$${find}-\$${find2}}' $file2 | nearneighbor -G$grdfile $R $interp\n";
  }
  print CSH "grdimage $grdfile -C$cptfile1 $J1 -K -O -V -Q >> $psfile0\n";
      }
      print CSH "pscoast $J1 $R $B -W1p -Na/1p -Dh -K -O -V >> $psfile0\n";

      # plot receivers
      if (0==1) {   # numbered large circles
  print CSH "awk '\$1 == \"R\" {print \$2,\$3}' $rec_file | psxy -N $J1 $R -K -O -V $rec0 >> $psfile0\n";
  $rec_file2 = text_rec; $angle = 0; $just = "CM";
  print CSH "awk '\$1 == \"R\" {print \$2,\$3,$fsize3,$angle,$fontno,\"$just\",\$4}' $rec_file > $rec_file2\n";
  print CSH "pstext $rec_file2 -N $J1 $R -K -O -V >> $psfile0\n";

      } else {      # small circles
  print CSH "awk '\$1 == \"R\" {print \$2,\$3}' $rec_file | psxy -N -Sc5p -W0.5p $J1 $R -K -O -V >> $psfile0\n";

      }

      # plot sources
      if ($ievent_one) {
  print CSH "awk '\$1 == \"S\" {print \$2,\$3}' $rec_file | psxy $src0 -N $J1 $R -K -O -V >> $psfile0\n";
      }
      if ($ievent_all) {
  print CSH "awk '{print \$1,\$2}' $evefile | psxy $src -N $J1 $R -K -O -V >> $psfile0\n";
      }

      # plot label
      $stlab = $stlabs[$k];
      $textinfo = "-N -C4p -W255/255/255o,1.0p,0/0/0,solid -G0/0/0";
      print CSH "pstext $R_title $J_title $textinfo -K -O -V $olab >>$psfile0<<EOF\n0 0 12 0 $fontno CM $stlab\nEOF\n";

      # colorscale bars
      if ($k==4) {
  print CSH "psscale -C$cptfile1 $Dscale1 $Bscale1 -K -O -V >> $psfile0\n";
      }
      if ($k==5) {
  print CSH "psscale -C$cptfile1 $Dscale1 $Bscale2 -K -O -V >> $psfile0\n";
      }

      # plot title and GMT header
      if ($k==5) {
  $plot_title = "Model m00, model $stmod, and target model (run_${stirun})";
  $Utag = "-U/0/0.25/$plabel";
  $shift = "-Xa-3.5 -Ya8.5";
  print CSH "pstext -N $J1 $R $Utag -O -V $shift >>$psfile0<<EOF\n $xmin $zmin $fsize0 0 $fontno LM $plot_title\nEOF\n"; # FINISH

  if ($ijpg == 1) {
    print CSH "convert $psfile0 $jpgfile\n";
  }
  if ($ipdf == 1) {
    print CSH "ps2pdf $psfile0\n";
  }

      }
    }       # loop over subplots
  }

  #-----------------------------
  # FIGURE 1

  if ($ifig1 == 1) {

    $name = "model_all_${stirun}_f1";
    $psfile1 = "$name.ps";
    $jpgfile = "$name.jpg";

    # write plotting scripts
    $wid = $wid1;
    $J1 = "-JM${wid}i";   # in lat-lon
    $origin = "-X1.5 -Y7.25";
    $Dlen = 2.0; $Dx = $wid/2; $Dy = -0.35;
    $Dscale1 = "-D$Dx/$Dy/$Dlen/0.15h";
    $Bscale1  = sprintf("-B%2.2e:\" Perturbation from %3.3f  km s\@+-1\@+ \": -E10p",$ptick,$beta0/1000);
    $Bscale2  = sprintf("-B%2.2e:\" Perturbation from m00 \": -E10p",$ptick);
    $xlab = $wid/2;
    $ylab = -0.02*$wid;
    $olab = "-Xa$xlab -Ya$ylab";

    print "\n----FIGURE 1--model $stmod----\n";

    @files = ("$dir0/${mfile_syn}","$dir0/${mfile_syn}",
        "$dir1/${mfile_syn}","$dir1/${mfile_syn}",
        "$dir0/${mfile_dat}","$dir0/${mfile_dat}");
    @stlabs = ("m00","m00 - m00",$stmod,"$stmod - m00","mtar","mtar - m00");

    # only plot even models (m04 = iteration 8)
    $kmaxt = 6;
    $kmaxm = 6*(1 - $mtest);  # =0 for test model
    if ($itest == 1) {
      $kmax = $kmaxt;
    } else {
      $kmax = $kmaxm;
    }

    # loop over subplots
    for ($k = 0; $k < $kmax; $k = $k+1) {

      $iB = $Bopts1[$k];
      $shift = $shifts1[$k];
      $B = "${B0}${iB}";
      $find = $fcol;
      $find2 = 2*$fcol;

      $file0 = $files[0]; # initial model
      $file1 = $files[$k];
      if (not -f $file1) {
  die("Check if $file1 exist or not\n");
      }

      # column-concatenate current model and initial model
      if ($k % 2 == 1) {
  $file2 = "ftemp";
  print CSH "paste $file1 $file0 > $file2\n";
      }

      # phase velocity map
      if ($k==0 && $m==$imodelmin) {
  # set bounds for the plotting
  #$xmin = -121; $xmax = -115.0; $zmin = 31.75; $zmax = 36.5;
  ($xmin,$xmax,$zmin,$zmax,$smin,$smax,$tmin,$tmax) = split(" ",`minmax -C $file1`);
  $dinc = 0.25;   # buffer around min/max
  $xmin = $xmin-$dinc;  $xmax = $xmax+$dinc;
  $zmin = $zmin-$dinc;  $zmax = $zmax+$dinc;
  $R = "-R$xmin/$xmax/$zmin/$zmax";
  print "\n$R\n";
      }

      if ($k==0) {
  print CSH "psbasemap $B $R $J1 -K -V -P $origin > $psfile1\n"; # START
      } else {
  print CSH "psbasemap $B $R $J1 -K -O -V $shift >> $psfile1\n";
      }

      if ($icolor == 1) {
  #print CSH "awk '{print \$1,\$2,\$7}' $file1 | pscontour $R $J1 -A- -C$cptfile1 -I -K -O -V >> $psfile1\n";
  if ($k % 2 == 0) {
    print CSH "awk '{print \$1,\$2,\$${find}}' $file1 | nearneighbor -G$grdfile $R $interp\n";
  } else {
    print CSH "awk '{print \$1,\$2,\$${find}-\$${find2}}' $file2 | nearneighbor -G$grdfile $R $interp\n";
  }
  print CSH "grdimage $grdfile -C$cptfile1 $J1 -K -O -V -Q >> $psfile1\n";
      }
      print CSH "pscoast $J1 $R $B -W1p -Na/1p -Dh -K -O -V >> $psfile1\n";

      # plot receivers
      if (0==1) {   # numbered large circles
  print CSH "awk '\$1 == \"R\" {print \$2,\$3}' $rec_file | psxy -N $J1 $R -K -O -V $rec0 >> $psfile1\n";
  $rec_file2 = text_rec; $angle = 0; $just = "CM";
  print CSH "awk '\$1 == \"R\" {print \$2,\$3,$fsize3,$angle,$fontno,\"$just\",\$4}' $rec_file > $rec_file2\n";
  print CSH "pstext $rec_file2 -N $J1 $R -K -O -V >> $psfile1\n";

      } else {      # small circles
  print CSH "awk '\$1 == \"R\" {print \$2,\$3}' $rec_file | psxy -N -Sc5p -W0.5p $J1 $R -K -O -V >> $psfile1\n";

      }

      # plot sources
      if ($ievent_one) {
  print CSH "awk '\$1 == \"S\" {print \$2,\$3}' $rec_file | psxy $src0 -N $J1 $R -K -O -V >> $psfile1\n";
      }
      if ($ievent_all) {
  print CSH "awk '{print \$1,\$2}' $evefile | psxy $src -N $J1 $R -K -O -V >> $psfile1\n";
      }

      # plot label
      $stlab = $stlabs[$k];
      $textinfo = "-N -C4p -W255/255/255o,1.0p,0/0/0,solid -G0/0/0";
      print CSH "pstext $R_title $J_title $textinfo -K -O -V $olab >>$psfile1<<EOF\n0 0 12 0 $fontno CM $stlab\nEOF\n";

      # colorscale bars
      if ($k==4) {
  print CSH "psscale -C$cptfile1 $Dscale1 $Bscale1 -K -O -V >> $psfile1\n";
      }
      if ($k==5) {
  print CSH "psscale -C$cptfile1 $Dscale1 $Bscale2 -K -O -V >> $psfile1\n";
      }

      # plot title and GMT header
      if ($k==5) {
  $plot_title = "Model m00, model $stmod, and target model (run_${stirun})";
  $Utag = "-U/0/0.25/$plabel";
  $shift = "-Xa-3.5 -Ya8.5";
  print CSH "pstext -N $J1 $R $Utag -O -V $shift >>$psfile1<<EOF\n $xmin $zmin $fsize0 0 $fontno LM $plot_title\nEOF\n"; # FINISH

  if ($ijpg == 1) {
    print CSH "convert $psfile1 $jpgfile\n";
  }
  if ($ipdf == 1) {
    print CSH "ps2pdf $psfile1\n";
  }

      }
    }       # loop over subplots
  }

  #-----------------------------
  # FIGURE 2

  # if you want fig 2 and the model is m01,m02,m03,etc (iterations 2,3,4,etc)
  if ( ($ifig2 == 1) && (($imodel >= 1) && ($mtest==0)) ) {

    # previous test model
    $imodelt = $imodel - 1;
    $stmodt  = sprintf("m%2.2it",$imodelt/2);
    $stirunt = sprintf("%4.4i",$irun0+$imodelt);
    $dir1t   = "$odir/run_${stirunt}";

    # previous model
    $imodelp = $imodel - 2;
    $stmodp  = sprintf("m%2.2i",$imodelp/2);
    $stirunp = sprintf("%4.4i",$irun0+$imodelp);
    $dir1p   = "$odir/run_${stirunp}";

    print "\n--------FIGURE 2---models $stmodp, $stmodt, $stmod---\n";

    $name2 = "model_all_${stirun}_f2";
    $psfile2 = "$name2.ps";
    $jpgfile2 = "$name2.jpg";

    # write plotting scripts
    $wid = $wid2;
    $J1 = "-JM${wid}i";   # in lat-lon
    $origin = "-X1 -Y7.25";
    $Dlen = 1.5; $Dx = $wid/2; $Dy = -0.35;
    $Dscale1 = "-D$Dx/$Dy/$Dlen/0.15h";
    $Bscale1  = sprintf("-B%2.2e:\" Perturbation from %3.3f  km s\@+-1\@+ \": -E10p",$ptick,$beta0/1000);
    $Bscale2  = sprintf("-B%2.2e:\" Perturbation from $stmodp \": -E10p",$ptick);
    $Bscale3  = sprintf("-B%2.2e:\" Perturbation from mtar\": -E10p",$ptick);
    $xlab = $wid/2;
    $ylab = -0.02*$wid;
    $olab = "-Xa$xlab -Ya$ylab";

    # nine subplots
    @files0 = ("$dir1p/${mfile_syn}","$dir1p/${mfile_syn}","$dir1p/${mfile_syn}",
               "$dir1t/${mfile_syn}","$dir1p/${mfile_syn}","$dir1t/${mfile_syn}",
               "$dir1/${mfile_syn}","$dir1p/${mfile_syn}","$dir1/${mfile_syn}");
    @files1 = ("$dir1p/${mfile_syn}","$dir1p/${mfile_syn}","$dir0/${mfile_dat}",
               "$dir1t/${mfile_syn}","$dir1t/${mfile_syn}","$dir0/${mfile_dat}",
               "$dir1/${mfile_syn}","$dir1/${mfile_syn}","$dir0/${mfile_dat}");
    @stlabs = ($stmodp,"$stmodp - $stmodp","mtar - $stmodp",
               $stmodt,"$stmodt - $stmodp","mtar - $stmodt",
               $stmod,"$stmod - $stmodp","mtar - $stmod");

    $kmax = 9;

    # loop over subplots
    for ($k = 0; $k < $kmax; $k = $k+1) {

      $iB = $Bopts2[$k];
      $shift = $shifts2[$k];
      $B = "${B0}${iB}";
      $find = $fcol;
      $find2 = 2*$fcol;

      $kmod = $k % 3;
      #print "\n -- k = $k -- kmod = $kmod --\n";

      # two files to concatenate
      $file0 = $files0[$k];
      $file1 = $files1[$k];
      if (not -f $file0) {
  die("Check if file0 $file0 exist or not\n");
      }

      # column-concatenate two files (columns 2 and 3)
      if ($kmod > 0) {
  if (not -f $file1) {
    die("Check if file1 $file1 exist or not\n");
  }
  $file2 = "ftemp";
        print "concatenating two files:\n  $file0\n  $file1\n";
  print CSH "paste $file1 $file0 > $file2\n";
      }

      # phase velocity map
      if ($k==0) {
  # set bounds for the plotting
  #$xmin = -121; $xmax = -115.0; $zmin = 31.75; $zmax = 36.5;
  ($xmin,$xmax,$zmin,$zmax,$smin,$smax,$tmin,$tmax) = split(" ",`minmax -C $file1`);
  $dinc = 0.25;   # buffer around min/max
  $xmin = $xmin-$dinc;  $xmax = $xmax+$dinc;
  $zmin = $zmin-$dinc;  $zmax = $zmax+$dinc;
  $R = "-R$xmin/$xmax/$zmin/$zmax";
  print "\n$R\n";
      }

      if ($k==0) {
  print CSH "psbasemap $B $R $J1 -K -V -P $origin > $psfile2\n"; # START
      } else {
  print CSH "psbasemap $B $R $J1 -K -O -V $shift >> $psfile2\n";
      }

      if ($icolor == 1) {
  #print CSH "awk '{print \$1,\$2,\$7}' $file1 | pscontour $R $J1 -A- -C$cptfile1 -I -K -O -V >> $psfile2\n";
  if ($kmod == 0) {
    print CSH "awk '{print \$1,\$2,\$${find}}' $file1 | nearneighbor -G$grdfile $R $interp\n";
  } else {
    print CSH "awk '{print \$1,\$2,\$${find}-\$${find2}}' $file2 | nearneighbor -G$grdfile $R $interp\n";
  }
  print CSH "grdimage $grdfile -C$cptfile1 $J1 -K -O -V -Q >> $psfile2\n";
      }
      print CSH "pscoast $J1 $R $B -W1p -Na/1p -Dh -K -O -V >> $psfile2\n";

      # plot receivers
      if (0==1) {   # numbered large circles
  print CSH "awk '\$1 == \"R\" {print \$2,\$3}' $rec_file | psxy -N $J1 $R -K -O -V $rec0 >> $psfile2\n";
  $rec_file2 = text_rec; $angle = 0; $just = "CM";
  print CSH "awk '\$1 == \"R\" {print \$2,\$3,$fsize3,$angle,$fontno,\"$just\",\$4}' $rec_file > $rec_file2\n";
  print CSH "pstext $rec_file2 -N $J1 $R -K -O -V >> $psfile2\n";

      } else {      # small circles
  print CSH "awk '\$1 == \"R\" {print \$2,\$3}' $rec_file | psxy -N -Sc5p -W0.5p $J1 $R -K -O -V >> $psfile2\n";

      }

      # plot sources
      if ($ievent_one) {
  print CSH "awk '\$1 == \"S\" {print \$2,\$3}' $rec_file | psxy $src0 -N $J1 $R -K -O -V >> $psfile2\n";
      }
      if ($ievent_all) {
  print CSH "awk '{print \$1,\$2}' $evefile | psxy $src -N $J1 $R -K -O -V >> $psfile2\n";
      }

      # plot label
      $stlab = $stlabs[$k];
      $textinfo = "-N -C4p -W255/255/255o,1.0p,0/0/0,solid -G0/0/0";
      print CSH "pstext $R_title $J_title $textinfo -K -O -V $olab >>$psfile2<<EOF\n0 0 10 0 $fontno CM $stlab\nEOF\n";

      # colorscale bars
      if ($k==6) {
  print CSH "psscale -C$cptfile1 $Dscale1 $Bscale1 -K -O -V >> $psfile2\n";
      }
      if ($k==7) {
  print CSH "psscale -C$cptfile1 $Dscale1 $Bscale2 -K -O -V >> $psfile2\n";
      }
      if ($k==8) {
  print CSH "psscale -C$cptfile1 $Dscale1 $Bscale3 -K -O -V >> $psfile2\n";
      }

      # plot title and GMT header
      if ($k==8) {
  $plot_title = "Models $stmodp, $stmodt, $stmod (run_${stirun})";
  $Utag = "-U/0/0.25/$plabel";
  $shift = "-Xa-4.5 -Ya7";
  print CSH "pstext -N $J1 $R $Utag -O -V $shift >>$psfile2<<EOF\n $xmin $zmin $fsize0 0 $fontno LM $plot_title\nEOF\n"; # FINISH

  if ($ijpg == 1) {
    print CSH "convert $psfile2 $jpgfile\n";
  }
  if ($ipdf == 1) {
    print CSH "ps2pdf $psfile2\n";
  }

      }
    }       # loop over subplots

  }

  #-----------------------------

}       # loop over models

close (CSH);
system("csh -f $cshfile");
if($ifig0==1) {system("gv $psfile0 &")}
if($ifig1==1) {system("gv $psfile1 &")}
if($ifig2==1) {system("gv $psfile2 &")}

#=================================================================
