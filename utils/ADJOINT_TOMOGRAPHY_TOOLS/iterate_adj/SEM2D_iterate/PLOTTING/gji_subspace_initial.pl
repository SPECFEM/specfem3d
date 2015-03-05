#!/usr/bin/perl -w

#==========================================================
#
#  gji_subspace_initial.pl
#  Carl Tape
#  15-Jan-2010
#
#  Copied from plot_surf_model_all.pl
#  This plots mprior, m00, mtar, and ln(mtar/m00)
#
#  Example:
#    ../gji_subspace_initial.pl 7200 0/0 1/0/5 0.1/0.05 0
#
#==========================================================

if (@ARGV < 5) {die("Usage: gji_subspace_initial.pl xxx\n");}
($irun0,$imodels,$ievents,$pbar,$ifinite) = @ARGV;
$comp = 1;

# USER INPUT
$ifig1 = 1;
$ifig2 = 0;
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

$plabel = "${plotdir}/gji_subspace_initial.pl";

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
$cshfile = "gji_subspace_initial.csh";
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

  # directories (current model and previos model)
  $stirun  = sprintf("%4.4i",$irun0+$imodel);
  $dir1    = "$odir/run_${stirun}";

  $name = "gji_subspace_initial";
  $psfile = "$name.ps";
  $jpgfile = "$name.jpg";

# write plotting scripts
$wid = $wid1;
$J1 = "-JM${wid}i";      # in lat-lon
$origin = "-X1.5 -Y7.25";
$Dlen = 2.0; $Dx = $wid/2; $Dy = -0.35;
$Dscale1 = "-D$Dx/$Dy/$Dlen/0.15h";
#$Bscale1  = sprintf("-B%2.2e:\" Perturbation from %3.3f  km s\@+-1\@+ \": -E10p",$ptick,$beta0/1000);
#$Bscale2  = sprintf("-B%2.2e:\" Perturbation from m00 \": -E10p",$ptick);
$Bscale  = sprintf("-B%2.2e:\" \": -E10p",$ptick);
$xlab = $wid/2;
$ylab = -0.02*$wid;
$olab = "-Xa$xlab -Ya$ylab";

    @files = ("$dir0/${mfile_syn}","$dir0/${mfile_syn}",
             "$dir1/${mfile_dat}","$dir0/${mfile_dat}");
    @stlabs = ("(a)  mprior","(b)  $stmod","(c)  mtar","(d)  ln(mtar/$stmod)");

    # only plot even models (m04 = iteration 8)
    $kmaxt = 4;
    $kmaxm = 4*(1 - $mtest);  # =0 for test model
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
  print CSH "psbasemap $B $R $J1 -K -V -P $origin > $psfile\n"; # START
      } else {
  print CSH "psbasemap $B $R $J1 -K -O -V $shift >> $psfile\n";
      }

      if ($icolor == 1) {
  #print CSH "awk '{print \$1,\$2,\$7}' $file1 | pscontour $R $J1 -A- -C$cptfile1 -I -K -O -V >> $psfile\n";
  if ($k == 0) {
          print CSH "awk '{print \$1,\$2,0}' $file0 | nearneighbor -G$grdfile $R $interp\n";
  } elsif ($k == 1) {
    print CSH "awk '{print \$1,\$2,\$${find}}' $file0 | nearneighbor -G$grdfile $R $interp\n";
  } elsif ($k == 2) {
          print CSH "awk '{print \$1,\$2,\$${find}}' $file1 | nearneighbor -G$grdfile $R $interp\n";
        } elsif ($k == 3) {
          $file2 = "ftemp";
    print CSH "paste $file1 $file0 > $file2\n";
    print CSH "awk '{print \$1,\$2,\$${find}-\$${find2}}' $file2 | nearneighbor -G$grdfile $R $interp\n";
        } else {
          die("k ($k) must be 0,1,2,3");
  }
  print CSH "grdimage $grdfile -C$cptfile1 $J1 -K -O -V -Q >> $psfile\n";
      }
      print CSH "pscoast $J1 $R $B -W1p -Na/1p -Dh -K -O -V >> $psfile\n";

      # plot receivers
      if (0==1) {   # numbered large circles
  print CSH "awk '\$1 == \"R\" {print \$2,\$3}' $rec_file | psxy -N $J1 $R -K -O -V $rec0 >> $psfile\n";
  $rec_file2 = text_rec; $angle = 0; $just = "CM";
  print CSH "awk '\$1 == \"R\" {print \$2,\$3,$fsize3,$angle,$fontno,\"$just\",\$4}' $rec_file > $rec_file2\n";
  print CSH "pstext $rec_file2 -N $J1 $R -K -O -V >> $psfile\n";

      } else {      # small circles
  print CSH "awk '\$1 == \"R\" {print \$2,\$3}' $rec_file | psxy -N -Sc5p -W0.5p $J1 $R -K -O -V >> $psfile\n";

      }

      # plot sources
      if ($ievent_one) {
  print CSH "awk '\$1 == \"S\" {print \$2,\$3}' $rec_file | psxy $src0 -N $J1 $R -K -O -V >> $psfile\n";
      }
      if ($ievent_all) {
  print CSH "awk '{print \$1,\$2}' $evefile | psxy $src -N $J1 $R -K -O -V >> $psfile\n";
      }

      # plot label
      $stlab = $stlabs[$k];
      $textinfo = "-N -C4p -W255/255/255o,1.0p,0/0/0,solid -G0/0/0";
      print CSH "pstext $R_title $J_title $textinfo -K -O -V $olab >>$psfile<<EOF\n0 0 12 0 $fontno CM $stlab\nEOF\n";

      # colorscale bars
      #if ($k==2) {print CSH "psscale -C$cptfile1 $Dscale1 $Bscale1 -K -O -V >> $psfile\n";}
      #if ($k==3) {print CSH "psscale -C$cptfile1 $Dscale1 $Bscale2 -K -O -V >> $psfile\n";}
      if ($k==2 || $k==3) {print CSH "psscale -C$cptfile1 $Dscale1 $Bscale -K -O -V >> $psfile\n";}

      # plot title and GMT header
      if ($k==3) {
  $plot_title = "Model m00, model $stmod, and target model (run_${stirun})";
  $Utag = "-U/0/0.25/$plabel";
  $shift = "-Xa-3.5 -Ya8.5";
  print CSH "pstext -N $J1 $R $Utag -O -V $shift >>$psfile<<EOF\n $xmin $zmin $fsize0 0 $fontno LM $plot_title\nEOF\n"; # FINISH

  if ($ijpg == 1) {
    print CSH "convert $psfile $jpgfile\n";
  }
  if ($ipdf == 1) {
    print CSH "ps2pdf $psfile\n";
  }

      }
    }       # loop over subplots

}       # loop over models

close (CSH);
system("csh -f $cshfile");
if($ifig1==1) {system("gv $psfile &")}
if($ifig2==1) {system("gv $psfile2 &")}

#=================================================================
