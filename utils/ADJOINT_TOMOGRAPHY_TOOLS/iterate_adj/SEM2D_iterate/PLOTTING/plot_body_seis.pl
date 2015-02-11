#!/usr/bin/perl -w

#==========================================================
#
#  plot_body_seis.pl
#  Carl Tape
#  07-Aug-2006
#
#-------------------------
#  EXAMPLES:
#    scripts/plot_body_seis.pl OUTPUT_body 0 5000 2 parameters1.log events_xz.dat recs_xz.dat 2
#    scripts/plot_body_seis.pl OUTPUT_body 0 5000 4 parameters1.log events_xz.dat recs_xz.dat 1/2/3
#    scripts/plot_body_seis.pl OUTPUT_body 0 5001 2 parameters1.log events_xy.dat recs_xz.dat 1/2/3
#    scripts/plot_body_seis.pl OUTPUT_body 0 5011 2 parameters1.log events_xy.dat recs_xz.dat 1/2/3
#    scripts/plot_body_seis.pl OUTPUT_body 0 5016 1 parameters1.log events_xz.dat recs_xz.dat 1/2/3
#
#    scripts/plot_body_seis.pl OUTPUT_body 1 1 5100 1  1 14 parameters1.log events_dat_xz.dat recs_xz.dat 1/2/3
#    scripts/plot_body_seis.pl OUTPUT_body 1 1 5100 1 16 30 parameters1.log events_dat_xz.dat recs_xz.dat 1/2/3
#
#    scripts/plot_body_seis.pl OUTPUT_body 1 1 5300 1  1  8 parameters1.log events_dat_xz.dat recs_xz.dat 1/2/3
#    scripts/plot_body_seis.pl OUTPUT_body 1 1 5300 1  9 16 parameters1.log events_dat_xz.dat recs_xz.dat 1/2/3
#
#==========================================================

if (@ARGV < 3) {die("Usage: plot_body_seis.pl basedir ievent parm_file xxx \n");}
($odir,$idata,$isyn,$irun,$ievent,$rmin,$rmax,$pfile,$efile,$rfile,$wtemp) = @ARGV;

$strun = sprintf("%4.4i",$irun);
$stev  = sprintf("%3.3i",$ievent);
$basedir = "$odir/run_$strun";
$dir = "$basedir/event_$stev";

$evefile = "$basedir/$efile";
$recfile = "$basedir/$rfile";

@comps   = split("/",$wtemp);     # indices of 1 or 2 components to plot
$ncomp = @comps;

# get parameters
$parmfile = "$basedir/$pfile";
if (not -f $parmfile) {die("Check if $parmfile exist or not\n")}
open(IN,"$parmfile"); @lines = <IN>;
($first,$end,$nint) = split(" ",$lines[0]);
($XFOR,$YFOR,$ZFOR) = split(" ",$lines[1]);
$stsrc = "source XYZ = ${XFOR}${YFOR}${ZFOR}";
$dt = split(" ",$lines[2]);
$hdur = split(" ",$lines[3]);
print "\n frame $first to $end in increment of $nint \n DT = $dt, HDUR =  $hdur \n $stsrc \n";

# get event
if (not -f $evefile) {die("Check if $evefile exist or not\n")}
open(IN,"$evefile"); @lines = <IN>;
($xsrc,$zsrc,$junk1,$namesrc) = split(" ",$lines[$ievent-1]);
print "\n source $namesrc is at (x, z) = ($xsrc, $zsrc) \n";

# get number of receivers
#open(IN,$recfile); @temp = <IN>; $nrec = @temp;
$nrec = $rmax-$rmin+1;
$strmin = sprintf("%2.2i",$rmin);
$strmax = sprintf("%2.2i",$rmax);

$iportrait = 1;
if($iportrait == 1){$orient = "-P"} else {$orient = " "}

# plotting specifications
$fsize0 = "18";
$fsize1 = "10";
$fsize2 = "8";
$fsize3 = "4";
if($ilabs == 0) {$fontno = "1"} else {$fontno = "4"}
$tick   = "0.1c";

@Bopts = ("WESN","Wesn","wesN","wEsn","weSn","WesN","wEsN","wESn","WeSn","WEsn","weSN","wESN","WESn","WeSN","WEsN","wesn");

#=================
$cshfile = "plot_body_seis.csh";
open(CSH,">$cshfile");
print CSH "gmtset BASEMAP_TYPE plain MEASURE_UNIT inch PLOT_DEGREE_FORMAT D TICK_LENGTH $tick LABEL_FONT_SIZE $fsize2 ANOT_FONT_SIZE $fsize2  HEADER_FONT $fontno ANOT_FONT $fontno LABEL_FONT $fontno HEADER_FONT_SIZE $fsize1\n";

@clabs = ("x","y","z");
if($ncomp==1) { $stcomp = $clabs[$comps[0]-1] }
if($ncomp==2) { $stcomp = $clabs[$comps[0]-1].$clabs[$comps[1]-1] }
if($ncomp==3) { $stcomp = $clabs[$comps[0]-1].$clabs[$comps[1]-1].$clabs[$comps[2]-1] }
print "\n $ncomp $stcomp \n";

#-------------------------

# write plotting scripts
$xwid = 2.5;
$ywid = 9/$nrec;
$J = "-JX$xwid/$ywid";

$dx0 = $xwid;
$dy0 = $ywid;
$dx1 = $dx0 + 0.1;
$dy1 = $dy0 + 0.;
$dy2 = ($nrec-1) * $dy1;

$shift1 = "-X0 -Y-$dy1";
$shift2 = "-X$dx1 -Y$dy2";

$origin = "-X0.5 -Y9.25";
$origin = "-X0.5 -Y8.75";

$Dscale = "-D$hwid/-0.40/$hwid/0.1h";            # colorbar

$name = "seis_${strun}_${stev}_${stcomp}_r${strmin}_r${strmax}";
$psfile  = "$name.ps";
$jpgfile = "$name.jpg";

print "\nWriting CSH file...\n";

#$nrec = 1;

$nplots = $nrec*$ncomp;

# loop over two selected components
for ($j = 0; $j < $ncomp; $j++) {

  $icomp = $comps[$j];
  $comp = $clabs[$icomp-1];
  print "\n component to plot is $comp ($icomp) \n";

  if($j==0){$locate = $origin}
  if($j==1){$locate = $shift2}
  if($j==2){$locate = $shift2}

# loop over seismograms
#for ($i = 0; $i < $nrec; $i++) {
for ($i = $rmin-1; $i <= $rmax-1; $i++) {

   $irec = $i+1;
   if($i > $rmin-1){$locate = $shift1}
   $iplot = $j*$nrec + $i+1 - ($rmin-1);

   if($idata==1){
     $seisfile = sprintf("$dir/dat_%05d_%1d",$irec,$icomp);
     if (not -f $seisfile) {die("Check if $seisfile exist or not\n")}
   }
   if($isyn == 1) {
     $seisfile2 = sprintf("$dir/syn_%05d_%1d",$irec,$icomp);
     if (not -f $seisfile2) {die("Check if $seisfile exist or not\n")}
   }

   print CSH "echo seismogram $iplot out of $nplots, $seisfile \n";

   # get the bounds for the seismogram
   $ss = 0;
   ($tmin,$tmax,$smin,$smax) = split(" ",`minmax -C $seisfile`);
   ($ss) = sort{$b <=> $a} ($ss,abs($smin),abs($smax));
   $ymax = $ss*1.1;
   $tmin = 0;
   $R = "-R$tmin/$tmax/-$ymax/$ymax";

   $tinc1 = 20;
   $tinc2 = 10;

   $B0 = sprintf("-Ba{%3.3f}f{%3.3f}:\" \":/%3.3f:\" \"::.\" \":wesn",$tinc1,$tinc2,$ss/2);
   $Bfull = sprintf("-Ba{%3.3f}f{%3.3f}:\"Time  (s)\":/%3.3f:\" z, h (t) \"::.\" \":weSn",$tinc1,$tinc2,$ss/2);
   if ($i == $rmax-1) {$B = $Bfull.$Bopts[4]} else {$B = $B0.$Bopts[15]}

   if ($i == $rmin-1 && $j == 0) {print CSH "psbasemap $J $R $B -K -V $orient $locate > $psfile\n"} # START
   else {print CSH "psbasemap $J $R $B -K -O -V $orient $locate >> $psfile\n"}

   # PLOT THE SEISMOGRAMS
   if( $isyn == 1) {print CSH "awk '{print \$1,\$2}' $seisfile2 | psxy $R $J $B -W0.5p,255/0/0,-- $orient -K -O -V >> $psfile\n"}
   if($idata == 1) {print CSH "awk '{print \$1,\$2}' $seisfile | psxy $R $J $B -W0.5p $orient -K -O -V >> $psfile\n"}

   $label = "$irec($comp)";
   $txtinfo = "-N -C3p -W255o";
   #print CSH "pstext -N -JX1i/1i -R0/1/0/1 -K -O -V $orient >>$psfile<<EOF\n 0.1 0.5 $fsize2 0 $fontno LM $label \nEOF\n";
   print CSH "pstext $txtinfo $J -R0/1/0/1 -K -O -V $orient >>$psfile<<EOF\n -0.05 0.5 $fsize2 0 $fontno LM $label \nEOF\n";

# plot title and GMT header
if ($i == $rmin-1 && $j == 0) {
  $plabel = "/home/denali2/carltape/wave2d/2d_adjoint/plot_body_seis.pl";
  $Utag = "-U/0/0.5i/$plabel";  # GMT header
  $shift0 = "-Xa0i -Ya${dy1}i";
  $title = "Seismograms for $dir/, $stsrc";
  print CSH "pstext -N -JX1i/1i -R0/1/0/1 $Utag -K -O $orient -V $shift0 >>$psfile<<EOF\n 0 0.25 $fsize1 0 $fontno LM $title\nEOF\n";
}

 }  # for $i

}  # for $j

#-------------------------
print CSH "pstext $J -R0/1/0/1 -O $orient -V >>$psfile<<EOF\n 10 10 $fsize0 0 $fontno CM junk \nEOF\n";  # FINISH
print CSH "echo done with $psfile $seisfile\n";
if($iportrait == 0){print CSH "convert $psfile -rotate 90 $jpgfile\n"}
if($iportrait == 1){print CSH "convert $psfile $jpgfile\n"}

close (CSH);
system("csh -f $cshfile");
system("xv $jpgfile &");
