#!/usr/bin/perl -w

#==========================================================
#
#  plot_body_wavefield.pl
#  Carl Tape
#  17-Jan-2007
#
#-------------------------
#  EXAMPLES:
#
#    scripts/plot_body_wavefield.pl OUTPUT_body 1 1 5100 1 parameters1.log 1/2/3 1200/1600/2000/2400/2800/3200 3/3/3 3.0/7.0/2.0 2 events_dat_xz.dat recs_xz.dat
#    scripts/plot_body_wavefield.pl OUTPUT_body 1 0 5100 1 parameters1.log 1/2/3 1200/1600/2400/3200/4000/6000/8000 3/3/3 3.0/7.0/2.0 2 events_dat_xz.dat recs_xz.dat
#    scripts/plot_body_wavefield.pl OUTPUT_body 0 0 5100 1 parameters1.log 1/2/3 1200/1600/2400/3200/4000/6000/8000 3/3/3 3.0/7.0/2.0 2 events_dat_xz.dat recs_xz.dat
#
#    scripts/plot_body_wavefield.pl OUTPUT_body 5011 1 parameters1.log 2 1200/1600/2000/2400/3200/5200/6400 3 10.0 2 events_xy_dat.dat recs_xz.dat
#    scripts/plot_body_wavefield.pl OUTPUT_body 5011 1 parameters1.log 1/2/3 1200/1600/2000/2800/3600/4800 3/3/3 3.0/7.0/2.0 2 events_xy_dat.dat recs_xz.dat
#    scripts/plot_body_wavefield.pl OUTPUT_body 5011 2 parameters1.log 1/2/3 1200/1600/2000/2800/3600/4800 3/3/3 3.0/10.0/3.0 2 events_xy_dat.dat recs_xz.dat
#    scripts/plot_body_wavefield.pl OUTPUT_body 5011 3 parameters1.log 1/2/3 1200/1600/2000/2800/3600/4800 3/3/3 3.0/7.0/2.0 2 events_xy_dat.dat recs_xz.dat
#    scripts/plot_body_wavefield.pl OUTPUT_body 5011 4 parameters1.log 1/2/3 1200/1600/2000/2800/3600/4800 3/3/3 2.5/7.0/1.5 2 events_xy_dat.dat recs_xz.dat
#
#    scripts/plot_body_wavefield.pl OUTPUT_body 5010 1 parameters1.log 2 1200/1600/2000/2400/3200/5200/6400 3 10.0 2
#    scripts/plot_body_wavefield.pl OUTPUT_body 5010 1 parameters1.log 1/2/3 1200/1600/2000/2800/3600/4800 3/3/3 3.0/7.0/2.0 2
#    scripts/plot_body_wavefield.pl OUTPUT_body 5010 2 parameters1.log 1/2/3 1200/1600/2000/2800/3600/4800 3/3/3 3.0/10.0/3.0 2
#    scripts/plot_body_wavefield.pl OUTPUT_body 5010 3 parameters1.log 1/2/3 1200/1600/2000/2800/3600/4800 3/3/3 3.0/7.0/2.0 2
#    scripts/plot_body_wavefield.pl OUTPUT_body 5010 4 parameters1.log 1/2/3 1200/1600/2000/2800/3600/4800 3/3/3 2.5/7.0/1.5 2
#
#    scripts/plot_body_wavefield.pl OUTPUT_body 5000 2 parameters1.log 2 1200/1600/2000/2400/3200/5200/6400 3 10.0 2
#    scripts/plot_body_wavefield.pl OUTPUT_body 5000 1 parameters1.log 1/2/3 800/1200/1600/2000/2400/2800 3/3/3 4.0/10.0/3.0 2
#    scripts/plot_body_wavefield.pl OUTPUT_body 5000 2 parameters1.log 1/2/3 800/1200/1600/2000/2400/2800 3/3/3 4.0/10.0/3.0 2
#    scripts/plot_body_wavefield.pl OUTPUT_body 5000 3 parameters1.log 1/2/3 800/1200/1600/2000/2400/2800 3/3/3 4.0/10.0/3.0 2
#    scripts/plot_body_wavefield.pl OUTPUT_body 5001 1 parameters1.log 1/2/3 800/1200/1600/2000/2400/2800 3/3/3 3.5/10.0/7.0 2
#    scripts/plot_body_wavefield.pl OUTPUT_body 5001 2 parameters1.log 1/2/3 800/1200/1600/2000/2400/2800 3/3/3 3.5/10.0/7.0 2
#    scripts/plot_body_wavefield.pl OUTPUT_body 5001 3 parameters1.log 1/2/3 800/1200/1600/2000/2400/2800 3/3/3 3.5/10.0/7.0 2
#
#==========================================================

if (@ARGV < 3) {die("Usage: plot_body_wavefield.pl basedir ievent parm_file f1/f2/f3 pwr1/pwr2/pwr3 c1/c2/c3 \n");}
($odir,$idata,$iportrait,$irun,$ievent,$pfile,$wtemp,$ftemp,$ptemp,$ctemp,$ilabs,$efile,$rfile) = @ARGV;

$strun = sprintf("%4.4i",$irun);
$stev  = sprintf("%3.3i",$ievent);

$basedir = "$odir/run_$strun";
$dir = "$basedir/event_$stev";

$evefile = "$basedir/$efile";
$recfile = "$basedir/$rfile";

@comps   = split("/",$wtemp);     # indices of 1 or 2 components to plot
#($first,$end,$nint) = split("/",$range);
@frames  = split("/",$ftemp);
@pwr     = split("/",$ptemp);     # PWR  : increase (larger negative power) for more contrast
@cmax    = split("/",$ctemp);     # CMAX : decrease for more contrast
$numf = @frames;
$ncomp = @comps;

# suffixes for the wavefield files
@labs = ("_syn","_dat");
$suffix = $labs[$idata];
#$suffix = "";

# get parameters
$parmfile = "$basedir/$pfile";
if (not -f $parmfile) {die("Check if $parmfile exist or not\n");}
open(IN,"$parmfile");
@lines = <IN>;
($first,$end,$nint) = split($lines[0]);
($XFOR,$YFOR,$ZFOR) = split(" ",$lines[1]);
$stsrc = "source XYZ = ${XFOR}${YFOR}${ZFOR}";
$dt = $lines[2];
$hdur = $lines[3];
print "\n frame $first to $end in increment of $nint \n DT = $dt, HDUR =  $hdur \n $stsrc \n";

# get event
if (not -f $evefile) {die("Check if $evefile exist or not\n")}
open(IN,"$evefile"); @lines = <IN>;
($xsrc,$zsrc,$junk1,$namesrc) = split(" ",$lines[$ievent-1]);
print "\n source $namesrc is at (x, z) = ($xsrc, $zsrc) \n";

# ilabs: how many labels you want on the plots
#  [0] for talks
#  [1] for publication
#  [2] for notes

# plot the color frames or not
$icolor = 1;    # ccc

if($iportrait == 1){
  $orient = "-P";
  $wid = 2.25;
  $origin = "-X0.7 -Y8.70";
} else {
  $orient = " ";
  $wid = 3.25;
  $origin = "-X0.7 -Y6.70";
}

# plotting specifications
$fsize0 = "18";
$fsize1 = "10";
$fsize2 = "8";
$fsize3 = "4";
if($ilabs == 0) {$fontno = "1"} else {$fontno = "4"}
$tick   = "0.1c";

# plot symbols for sources, receivers, and shelf
if($ifinite==0) {$src = "-W0.5p -Sa0.10";} else {$src = "-Sc0.05";}
$rec = "-Si6p -W0.5p,0/0/0";
$src = "-Sa10p -W0.5p -G255/255/255";   # nice white source (GJI figures)

@Bopts = ("WESN","Wesn","wesN","wEsn","weSn","WesN","wEsN","wESn","WeSn","WEsn","weSN","wESN","WESn","WeSN","WEsN","wesn");

#=================
$cshfile = "plot_body_wavefield.csh";
open(CSH,">$cshfile");
print CSH "gmtset BASEMAP_TYPE plain MEASURE_UNIT inch PLOT_DEGREE_FORMAT D TICK_LENGTH $tick LABEL_FONT_SIZE $fsize2 ANOT_FONT_SIZE $fsize2  HEADER_FONT $fontno ANOT_FONT $fontno LABEL_FONT $fontno HEADER_FONT_SIZE $fsize1\n";

# KEY: scaling for color
$scale_color = 21.0;
$colorbar = "seis";
@norm = ("1e-$pwr[0]","1e-$pwr[1]","1e-$pwr[2]");

print "@norm \n";

@clabs = ("x","y","z");
@wavefield=("forward${suffix}","adjoint","kernel");
$numw = @wavefield;

if($ncomp==1) { $stcomp = $clabs[$comps[0]-1] }
if($ncomp==2) { $stcomp = $clabs[$comps[0]-1].$clabs[$comps[1]-1] }
if($ncomp==3) { $stcomp = $clabs[$comps[0]-1].$clabs[$comps[1]-1].$clabs[$comps[2]-1] }

# bounds for plotting (we could also use the file global_mesh.dat)
$firstframe = "$dir/forward${suffix}_00000";
if (not -f $firstframe) {die("Check if $firstframe exist or not\n");}
#($xmin,$xmax,$zmin,$zmax) = split(" ",`minmax -C -I1/1 $firstframe`);
($xmin,$xmax,$zmin,$zmax) = split(" ",`minmax -C $firstframe`);
$xmin = $xmin/1000; $xmax = $xmax/1000;
$zmin = $zmin/1000; $zmax = $zmax/1000;
$pinc = 0.0;  # buffer around min/max
$xran = $xmax - $xmin; $zran = $zmax - $zmin;
$xmin = $xmin-$pinc*$xran;  $xmax = $xmax+$pinc*$xran;
$zmin = $zmin-$pinc*$zran;  $zmax = $zmax+$pinc*$zran;
$R = "-R$xmin/$xmax/$zmin/$zmax";

$k = 0;

# loop over two selected components
for ($j = 0; $j < $ncomp; $j++) {

  $ss[$j] = $cmax[$j];
  $ds[$j] = 2*$ss[$j]/$scale_color;
  $bs[$j] = sprintf("%3.3e",0.9*$ss[$j]);  # colorbar
  $Ts[$j] = sprintf("-T%3.3e/%3.3e/%3.3e",-$ss[$j],$ss[$j],$ds[$j]);
  print "Ts = $Ts[$j]\n";
  print CSH "makecpt -C$colorbar $Ts[$j] -D > color_${j}.cpt\n";
}

#-------------------------

# write plotting scripts
$hwid = $wid/2;
$dx0 = $wid;
$dy0 = $dx0*($zran/$xran);
$J = "-JX${dx0}i/${dy0}i";

$dx1 = $dx0*1.05;
$dy1 = $dy0*1.5;
$dy2 = ($numf-1) * $dy1;

$shift1 = "-X0 -Y-$dy1";
$shift2 = "-X$dx1 -Y$dy2";

$Dscale = "-D$hwid/-0.40/$hwid/0.1h";            # colorbar

$name = "wavefield${suffix}_${strun}_${stev}_${stcomp}";
$psfile  = "$name.ps";
$jpgfile = "$name.jpg";

print "\nWriting CSH file...\n";

#$numf = 1;

$nplots = $numf*$ncomp;

# loop over two selected components
for ($j = 0; $j < $ncomp; $j++) {

  $comp = $clabs[$comps[$j]-1];
  print "\n component to plot is $comp ($comps[$j]) \n";

  if($j==0){$locate = $origin}
  if($j==1){$locate = $shift2}
  if($j==2){$locate = $shift2}

  # color bar
  $BscaleS1 = sprintf("-B%2.2e:\" s\@-$comp\@- ( x, z, t )  ( 10\@+-%2.2i\@+  m )\": -E10p",$bs[$j]/2,$pwr[$j]);

# loop over time frames
for ($i = 0; $i < $numf; $i++) {

   $j1 = $frames[$i];           # forward frame
   if($i > 0){$locate = $shift1}

   $time = sprintf("%04d",$j1*$dt);

   $snapshot_f = sprintf("$dir/$wavefield[0]_%05d",$j1);
   if (not -f $snapshot_f) {die("Check if $snapshot_f exist or not\n");}

   $iplot = $j*$numf + ($i+1);
   print CSH "echo plot $iplot out of $nplots \n";

   $Bfull = "-B50/50:\"t = $time s\"::.\"  \":";
   $B = $Bfull.$Bopts[15];
   if ($j == 0) {$B = $Bfull.$Bopts[1]}
   if ($i==$numf-1) {$B = $Bfull.$Bopts[4]}
   if ($i==$numf-1 && $j==0) {$B = $Bfull.$Bopts[8]}

   if ($i == 0 && $j == 0) {print CSH "psbasemap $J $R $B -K $orient -V $locate > $psfile\n"} # START
   else {print CSH "psbasemap $J $R $B -K -O $orient -V $locate >> $psfile\n"}

   # PLOT THE FORWARD WAVEFIELD
   #if($icolor==1) {print CSH "awk '{print \$1/1000,\$2/1000,\$($comps[$j]+2) / $norm[$j]}' $snapshot_f | pscontour $J $R $B -A- -Ccolor_${j}.cpt -I -K -O $orient -V >> $psfile\n"}

   if($icolor==1) {
      $grdfile = "temp.grd";
      $grdinfo = "-I1/1 -S2";   # key information
      #print CSH "awk '{print \$1/1000,\$2/1000,\$($comps[$j]+2) / $norm[$j]}' $snapshot_f > test\n";
      print CSH "awk '{print \$1/1000,\$2/1000,\$($comps[$j]+2) / $norm[$j]}' $snapshot_f | nearneighbor -G$grdfile $R $grdinfo\n";
      print CSH "grdimage $grdfile -Ccolor_${j}.cpt $J -K -O $orient -V >> $psfile\n";
   }

  # plot discontinutities
  print CSH "psxy -N $J $R -K -O -P -V -W0.5p >> $psfile<<EOF\n 0 45 \n 400 45 \nEOF\n";
  print CSH "psxy -N $J $R -K -O -P -V -W0.5p >> $psfile<<EOF\n 0 65 \n 400 65 \nEOF\n";
  print CSH "psxy -N $J $R -K -O -P -V -W0.5p >> $psfile<<EOF\n 0 75 \n 400 75 \nEOF\n";

   if ($i == $numf-1 && $ilabs > 0) {print CSH "psscale -Ccolor_${j}.cpt $Dscale $BscaleS1 -K -O $orient -V >> $psfile \n"}

   # plot source and receivers
   $xs = $xsrc/1000; $zs = $zsrc/1000;
   print CSH "psxy -N $J $R -K -O -V $orient $src >> $psfile<<EOF\n $xs $zs \nEOF\n";
   print CSH "awk '{print \$1/1000,\$2/1000}' $recfile |psxy -N $J $R -K -O $orient -V $rec >> $psfile\n";
   #print CSH "awk '\$1 == \"R\" {print \$2/1000,\$3/1000}' $dir/sr.txt |psxy -N $J $R -K -O $orient -V $rec >> $psfile\n";
   #print CSH "awk '\$1 == \"S\" {print \$2/1000,\$3/1000}' $dir/sr.txt |psxy -N $J $R -K -O $orient -V $src >> $psfile\n";

   # plot title
   $xtx = $xmin+0.5*($xmax-$xmin); $ztx = $zmin+1.15*($zmax-$zmin);
   #if ($i == $numf-1) {print CSH "pstext -N $J $R -K -O $orient -V >>$psfile<<EOF\n $xtx $ztx $fsize1 0 $fontno CM $titles[0]\nEOF\n";}

# plot title and GMT header
if ($ilabs == 2 && $i == 0) {
  $plabel = "/home/denali2/carltape/wave2d/2d_adjoint/plot_body_wavefield.pl";
  $Utag = "-U/0/0.2i/$plabel";  # GMT header
  if($j > 0) {$Utag = " "}
  $shift0 = "-Xa0i -Ya${dy1}i";
  $title = "$comp comp of disp, $stsrc";
  print CSH "pstext -N $J $R $Utag -K -O $orient -V $shift0 >>$psfile<<EOF\n 0 0 $fsize1 0 $fontno LM $title\nEOF\n";
}

 }  # for $i

}  # for $j

#-------------------------
print CSH "pstext $J -R0/1/0/1 -O $orient -V >>$psfile<<EOF\n 10 10 $fsize0 0 $fontno CM junk \nEOF\n";  # FINISH
print CSH "echo done with $psfile $snapshot_f\n";
if($iportrait == 0){print CSH "convert $psfile -rotate 90 $jpgfile\n"}
if($iportrait == 1){print CSH "convert $psfile $jpgfile\n"}

close (CSH);
system("csh -f $cshfile");
system("xv $jpgfile &");
