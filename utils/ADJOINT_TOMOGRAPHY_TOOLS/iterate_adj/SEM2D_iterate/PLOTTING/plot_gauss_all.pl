#!/usr/bin/perl -w

#==========================================================
#
#  plot_gauss_all.pl
#  Carl Tape
#  15-Jan-2010
#
#
#  Set of 1000 models
#    ../plot_gauss_all.pl
#
#==========================================================

#if (@ARGV < 5) {die("Usage: plot_gauss_all.pl xxx\n");}
#($im,$imodel,$pbar) = @ARGV;

$imodel = 0;
$im = 1;
$irun0 = $im;

# base directory
$pwd = $ENV{PWD};
$dirplot = `dirname $pwd`; chomp($dirplot);
$basedir = `dirname $dirplot`; chomp($basedir);
#$basedir = "/data1/cig/seismo/3D/ADJOINT_TOMO/iterate_adj/SEM2D_iterate_work";

$pbar = "0.1/0.05";
($pmax,$ptick) = split("/",$pbar);
$mfile_syn = "structure_syn.dat";

# directories
$plotdir = "${basedir}/PLOTTING";
$odir = "${basedir}/OUTPUT";    # might need to be modified
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

# source and receivers
$rfile = "$dir0/recs_lonlat.dat";
$efile = "$dir0/events_syn_lonlat.dat";
if (not -f $efile) { die("Check if efile $efile exist or not\n") }
if (not -f $rfile) { die("Check if rfile $rfile exist or not\n") }

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

$plabel = "${plotdir}/plot_gauss_all.pl";

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
$src = "-W0.75p -Sa0.20";
$rec = "-W0.5p,0/0/0 -Sc3p";
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
  $cshfile = "plot_gauss_all.csh";
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

  $origin = "-X1.5 -Y5.0";
  $B1syn = "-B1:.\" \":WESN";

#=============================================

$pmin = 1; $pmax = 1000;
#$pmin = 95; $pmax = 97;

$ibest = 1; $ik = 0;

$idir = "/data1/cig/seismo/3D/ADJOINT_TOMO/iterate_adj/SEM2D_iterate_INPUT/random_fields/model_set_0001";

for ($p = $pmin; $p <= $pmax; $p ++ ) {

  $im = $p;
  $stim = sprintf("%4.4i",$im);
  $file1syn = "$idir/structure_model_${stim}.dat";
  if (not -f $file1syn) {
    die("Check if file1syn $file1syn exist or not\n");
  }

  $stirun = sprintf("%4.4i",$irun0+$im-1);
  $chifile = "$odir/run_${stirun}/chi.dat";
  if (not -f $chifile) {die("Check if chifile $chifile exist or not\n");}
  open(IN,"$chifile");
  $chi = <IN>;
  if ($im==933) {$chi = 0}   # target model
  #$schi = sprintf("F(m) = %.1f",$chi);
  $schi = sprintf("%.1f",$chi);
  close(IN);

  if ( $ibest==0 || ($ibest==1 && $chi <= 100) ) {
    #print "$p, $schi\n";
    #print CSH "echo $p, $schi\n";

    if($ibest==1) {
      $ik = $ik + 1;
      $name = "gaussian_cov_model_ibest_${stim}";
      #$name = "gaussian_cov_model_ibest_$ik";
    } else {
      $name = "gaussian_cov_model_${stim}";
    }
    $psfile  = "$name.ps";
    $jpgfile = "$name.jpg";

    # phase velocity map
    print CSH "psbasemap $B1syn $R $J1 -P -K -V $origin > $psfile\n"; # START
    #print CSH "awk '{print \$1,\$2,\$6}' $file1syn | pscontour $R $J1 -A- -C$cptfile1 -I -K -O -V >> $psfile\n";
    print CSH "awk '{print \$1,\$2,\$6}' $file1syn | nearneighbor -G$grdfile $R $interp\n";
    print CSH "grdimage $grdfile -C$cptfile1 $J1 -K -O -V -Q >> $psfile\n";
    print CSH "pscoast $J1 $R $B1syn -W1p -Na/1p -Dh -K -O -V >> $psfile\n";
    #print CSH "awk '{print \$2,\$1}' $idir/BOUNDARIES/oms_shelf |psxy $J1 $R $Wshelf -K -O -V >> $psfile\n";
    print CSH "psscale -C$cptfile1 $Dscale1 $Bscale1 -K -O -V >> $psfile \n";

    if (1==1) {
      # faults
      $faultfile = "/home/carltape/gmt/faults/jennings_more.xy";
      if (not -f $faultfile) {
  die("Check if faultfile $faultfile exist or not\n");
      }
      print CSH "psxy $faultfile $J1 $R -W0.75p -M -K -V -O >> $psfile\n";

      # sources and receivers, plus chi label
      if ($ibest==1) {
  print CSH "pstext -N $J1 $R -K -O -V -Xa-0 -Ya-0.9 >>$psfile<<EOF\n $xmin $zmin $fsize0 0 $fontno CM $schi\nEOF\n";
        print CSH "awk '{print \$1,\$2}' $efile | psxy $R $J1 $src -K -O -V >> $psfile\n";
        print CSH "awk '{print \$1,\$2}' $rfile | psxy $R $J1 $rec -K -O -V >> $psfile\n";
      }
    }

    # plot title and GMT header
    $Utag = "-U/0/0.25/$plabel";
    $shift = "-Xa-0.5 -Ya-0.5";

    $stim = sprintf("%4.4i",$im);
    print CSH "pstext -N $J1 $R -K -O -V -Xa-0 -Ya-0.4 >>$psfile<<EOF\n $xmin $zmin $fsize0 0 $fontno CM Model\nEOF\n";
    print CSH "pstext -N $J1 $R -O -V -Xa-0 -Ya-0.65 >>$psfile<<EOF\n $xmin $zmin $fsize0 0 $fontno CM $stim\nEOF\n"; # FINISH

  }  # ibest = 1
}

#-----------------------------
#  print CSH "convert $psfile $jpgfile\n";

  close (CSH);
  system("csh -f $cshfile");
  system("gv $psfile &")
  #if($iopt <= 2) {system("xv $jpgfile &")}

#=================================================================
