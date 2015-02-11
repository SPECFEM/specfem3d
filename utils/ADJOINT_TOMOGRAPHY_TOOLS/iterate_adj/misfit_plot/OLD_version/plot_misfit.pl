#!/usr/bin/perl -w

#-----------------------------------
# Carl Tape, 02-Sept-2008
# plot_misfit.pl
#
#
# EXAMPLE:
#    plot_misfit.pl 9818433 USC/CI 6/40
#
#-----------------------------------

#use Getopt::Std;
#use POSIX;

if (@ARGV < 3) {die("Usage: plot_misfit.pl eid sta/net Tmin/Tmax\n");}
($eid,$stanet,$Ts) = @ARGV;

#($sta,$net) = split(/\//,$stanet);
($sta,$net) = split("/",$stanet);

# bandpass range
($Tmin,$Tmax) = split("/",$Ts);
$sTmin = sprintf("T%3.3i",$Tmin);
$sTmax = sprintf("T%3.3i",$Tmax);
$Ttag = "${sTmin}_${sTmax}";

#---------------------------
# USER INPUT

# models to compare
$minm = 0;
$maxm = 7;
$smodel_min = sprintf("m%2.2i",$minm);
$smodel_max = sprintf("m%2.2i",$maxm);

# base directory (contains STATIONS and eid_list)
#$dir0 = "/net/denali/scratch1/carltape/svn/cig/seismo/3D/ADJOINT_TOMO/iterate_adj_work/matlab/OUTPUT_MISFIT";

# CMTSOLUTION file
$cmtfile = "../CMTSOLUTION";
@cmtvals = &get_cmt($cmtfile);   # see subroutine below

# STATIONS file
$stafile = "../STATIONS_ADJOINT";

# plot faults (southern California)
$ifault = 1;

# file to use to extract the SAC headers
$Zsyn = "${sta}.${net}.BHZ.semd.sac.d.${Ttag}.${smodel_max}";

$saclst = "/opt/seismo-util/bin/saclst";

# axes positions
$x0_title = 0.5; $y0_title = 10.5;
$dY_under_title = 0.4;
$dY_under_subtitle = 0.4;
$dX_between_col = 0.2;
$dY_under_seis = 0.3;

# dimensions of seismogram subplots
$xwid = 2.4;
$ywid = 1.0;
$J = "-JX${xwid}/${ywid}";

# subtitles
@subs = ("(a)  Vertical component, Z","(b)  Radial component, R","(c)  Transverse component, T","(d)  Variance reduction for each component","(e) Source-receiver geometry");

@slabs = ("d and s($smodel_min)","d and s($smodel_max)","s($smodel_min) - d","s($smodel_max) - d","s($smodel_max) - s($smodel_min)");
$nrows = 5;    # <= 5

$x0_map = 5.0; $y0_map = 2.3 + (5-$nrows)*$ywid;
$x0_misfit = 1.5; $y0_misfit = 2.9 + (5-$nrows)*$ywid;

#---------------------------

# create psmeca file from CMTSOLUTION file
@M[0..5] = @cmtvals[12..17];
$year = $cmtvals[0]; $month = $cmtvals[1]; $day = $cmtvals[2];
$elat = $cmtvals[9]; $elon = $cmtvals[10]; $edep = $cmtvals[11];
$eid = $cmtvals[18];
$cmtpsmeca = "temp.txt";
open(MIJ,">$cmtpsmeca");
#printf MIJ "$elon $elat $edep $M[0] $M[1] $M[2] $M[3] $M[4] $M[5] 1 X Y $eid";
printf MIJ "$elon $elat $edep $M[0] $M[1] $M[2] $M[3] $M[4] $M[5] 1";
close(MIJ);

# compute moment using Dahlen and Tromp, Eq. (5.91)
# compute moment magnitude using Kanamori (1977) -- see Carl's Latex notes cmt.pdf
$M0 = 1/sqrt(2) * sqrt($M[0]**2 + $M[1]**2 + $M[2]**2 + 2*($M[3]**2 + $M[4]**2 + $M[5]**2 ) );
$k  = -11.2 + (2/3)*log(5)/log(10);
$Mw = (2/3)*log($M0)/log(10) + $k;
$stM0 = sprintf("M0 = %.3e dyne-cm",$M0);
$stMw = sprintf("Mw = %.2f",$Mw);
print "\n           moment : $stM0 \n moment magnitude : $stMw\n";

# get the receiver location, azimuth, distance (from the $Zsyn file)
# We want the azimuth for the file name, so that all the PDFs will be sorted by azimuth.
(undef,$tmin0,$tmax0,$slon,$slat,$dist,$az,$Tchan) = split(" ",`$saclst b e stlo stla dist az kcmpnm f $Zsyn`);
$edate  = sprintf("%4.4i . %2.2i . %2.2i",$year,$month,$day);
$stedep = sprintf("depth = %.1f km",$edep);
$strec  = "$sta . $net";
$stdist = sprintf("dist = %.1f km",$dist);
$staz   = sprintf("az = %.1f deg",$az);

#---------------------------

$fontno = 1;  # font number to use for text

# set GMT defaults
`gmtset BASEMAP_TYPE plain PAPER_MEDIA letter ANNOT_FONT_SIZE = 10p LABEL_FONT_SIZE = 10p PAGE_ORIENTATION = portrait PLOT_DEGREE_FORMAT D TICK_LENGTH 4p`;

# name of the output files
$azlabel = sprintf("%3.3i",$az);
$name = "${azlabel}_${sta}_${net}_${Ttag}";
$psfile = "$name.ps";
$jpg_file = "$name.jpg";

# overall title
$title = "Event $eid,  station ${sta}.${net},  bandpass ${Tmin}-${Tmax} s,  models ${smodel_min} and ${smodel_max}";
`pstext -R0/1/0/1 -JX1 -N -Xa${x0_title} -Ya${y0_title} -K -P > $psfile <<EOF\n 0 0 14 0 1 LM $title \nEOF\n`;  # START (no -O)

# subtitles
$x0_sub = $x0_title;
$y0_sub = $y0_title - $dY_under_title;

for ($k = 1; $k <= 3; $k = $k+1) {
  $x_sub = $x0_sub + ($k-1)*($xwid+$dX_between_col);
  $origin_sub = "-Xa${x_sub} -Ya${y0_sub}";
  `pstext -R0/1/0/1 -JX1 -N $origin_sub -K -O >> $psfile <<EOF\n 0 0 10 0 1 LM $subs[$k-1] \nEOF\n`;
}

# subtitle for misfit plot
$x0_sub2 = $x0_misfit - 0.8;
$y0_sub2 = $y0_title - $dY_under_title - $dY_under_subtitle - $nrows*$ywid - $dY_under_seis;
$origin_sub = "-Xa${x0_sub2} -Ya${y0_sub2}";
`pstext -R0/1/0/1 -JX1 -N $origin_sub -K -O >> $psfile <<EOF\n 0 0 10 0 1 LM $subs[3] \nEOF\n`;

# subtitle for map
$x0_sub3 = $x0_map - 0.8;
$y0_sub3 = $y0_sub2;
$origin_sub = "-Xa${x0_sub3} -Ya${y0_sub3}";
`pstext -R0/1/0/1 -JX1 -N $origin_sub -K -O >> $psfile <<EOF\n 0 0 10 0 1 LM $subs[4] \nEOF\n`;

#---------------------------------------------------------------------------
# PLOT SEISMOGRAMS AND RESIDUALS

$x0_seis = $x0_title;
$y0_seis = $y0_title - $dY_under_title - $dY_under_subtitle - $ywid;

$B0 = "-Ba50f10g50:\" \":/a1f1::wesN";
$B = "-Ba50f10g50:\" \":/a1f1::wesn";

#$slabinfo = "-C2p -W0/255/255o,0/255/255 -G0/0/0 -N";
$slabinfo = "-C3p -W255/255/255o,0/0/0 -G0/0/0 -N";

@comps = ("Z","R","T");
for ($k = 1; $k <= 3; $k = $k+1) {
  $comp = $comps[$k-1];

  $Pfile = "GMT_time_series_${Ttag}_${comp}.dat";
  $bounds = "GMT_axes_${Ttag}_${comp}.dat";
  $Nfile = "GMT_norms_${Ttag}_${comp}.dat";

  # if the seismograms exist, then plot them
  if ( (-f $Pfile) && (-f $bounds) ) {

    open(IN,"$bounds"); @lines = <IN>;
    ($tmin,$tmax,$ymin,$ymax) = split(" ",$lines[0]);
    $tmin = -10;
    if($Tmin==2) {$tmax = $tmax/2;}    # USER: reduce the records for 2s
    $R = "-R$tmin/$tmax/$ymin/$ymax";

    # plot seismograms
    $x_seis = $x0_seis + ($k-1)*($xwid+$dX_between_col);
    $x_seis_lab = $x_seis + $xwid - 0.05;

    for ($j = 1; $j <= $nrows; $j = $j+1) {
      $y_seis = $y0_seis - ($j-1)*$ywid;
      $originP = "-Xa${x_seis} -Ya${y_seis}";

      $y_seis_lab = $y_seis + $ywid - 0.05;
      $originPlab = "-Xa${x_seis_lab} -Ya${y_seis_lab}";

      if($j==1) {`psbasemap $J $R $B0 -O -K $originP >> $psfile`}
      else {`psbasemap $J $R $B -O -K $originP >> $psfile`}
      if ($j==1) {
  `awk \'{print \$1,\$2}\' $Pfile | psxy $J $R -W0.5p,0/0/0   -O -K $originP >> $psfile`;
  `awk \'{print \$1,\$3}\' $Pfile | psxy $J $R -W0.5p,255/0/0 -O -K $originP >> $psfile`;
        `pstext -R0/1/0/1 -JX1 $slabinfo -K -O $originPlab >> $psfile <<EOF\n0 0 8 0 $fontno RT $slabs[0]\nEOF\n`;
      } elsif ($j==2) {
  `awk \'{print \$1,\$2}\' $Pfile | psxy $J $R -W0.5p,0/0/0   -O -K $originP >> $psfile`;
  `awk \'{print \$1,\$4}\' $Pfile | psxy $J $R -W0.5p,255/0/0 -O -K $originP >> $psfile`;
        `pstext -R0/1/0/1 -JX1 $slabinfo -K -O $originPlab >> $psfile <<EOF\n0 0 8 0 $fontno RT $slabs[1]\nEOF\n`;
      } elsif ($j==3) {
  `awk \'{print \$1,\$5}\' $Pfile | psxy $J $R -W0.5p,0/0/0 -O -K $originP >> $psfile`;
        `pstext -R0/1/0/1 -JX1 $slabinfo -K -O $originPlab >> $psfile <<EOF\n0 0 8 0 $fontno RT $slabs[2]\nEOF\n`;
      } elsif ($j==4) {
  `awk \'{print \$1,\$6}\' $Pfile | psxy $J $R -W0.5p,0/0/0 -O -K $originP >> $psfile`;
        `pstext -R0/1/0/1 -JX1 $slabinfo -K -O $originPlab >> $psfile <<EOF\n0 0 8 0 $fontno RT $slabs[3]\nEOF\n`;
      } elsif ($j==5) {
  `awk \'{print \$1,\$7}\' $Pfile | psxy $J $R -W0.5p,0/0/0 -O -K $originP >> $psfile`;
         `pstext -R0/1/0/1 -JX1 $slabinfo -K -O $originPlab >> $psfile <<EOF\n0 0 8 0 $fontno RT $slabs[4]\nEOF\n`;
      }
    }

  } else {
    print "seismogram files $Pfile and $bounds do not both exist\n";
  }

}

#---------------------------------------------------------------------------
# PLOT THE MISFIT

$originN = "-Xa${x0_misfit} -Ya${y0_misfit}";

$xwid = 1.5;
$ywid = 1.3;
$J = "-JX${xwid}/${ywid}";
$B = "-Ba1f1:\"Model index\":/a50f10:\"Variance reduction\":WeSn";
$kmax = $maxm + 3;
$R = "-R-0.5/$kmax/-50/100";

$Zinfo = "-Sc12p -W0.5p,0/0/0 -G255 -N";

print "J option --> $J\n";
print "R option --> $R\n";
print "B option --> $B\n";

@comps = ("Z","R","T");
for ($k = 1; $k <= 3; $k = $k+1) {
  $comp = $comps[$k-1];

  $Nfile = "GMT_norms_${Ttag}_${comp}.dat";

  # if the norms file exists
  if ( -f $Nfile ) {
    # reduction in the misfit
    open(IN,"$Nfile"); @Nlines = <IN>;
    #$yred = 100*(1.0 - $Nlines[4]/$Nlines[3]);
    $yred = $Nlines[8]; chomp($yred);
    print "yred --> $yred\n";

    # plot the misfit
    if($k==1) {`psbasemap $J $R $B -O -K $originN >> $psfile`}
    `psxy $R $J -W1p -K -O $originN >> $psfile <<EOF\n$minm 0\n$maxm $yred\nEOF\n`;
    `psxy $R $J $Zinfo -K -O $originN >> $psfile <<EOF\n$minm 0\nEOF\n`;
    `psxy $R $J $Zinfo -K -O $originN >> $psfile <<EOF\n$maxm $yred\nEOF\n`;
    `pstext $R $J -N -K -O $originN >> $psfile <<EOF\n$maxm $yred 8 0 $fontno CM $comp\nEOF\n`;

  } else {
    print "norm file $Nfile does not exist\n";
  }

}

#---------------------------------------------------------------------------
# PLOT BASE MAP with stations, faults, CMT

$originM = "-Xa${x0_map} -Ya${y0_map}";
$proj = "-JM2.5";
$textinfo = "-G0 -S1p,255/255/255";

$xmin = -122; $xmax = -114;
$ymin = 32; $ymax = 37;
$bounds="-R$xmin/$xmax/$ymin/$ymax ";
$xtick = 1; $ytick = 1;
$tick="-B${xtick}/${ytick}::wesn";
$cinfo="-W0.5 -Dh -A100";

# plot the base map
`psbasemap $bounds $proj $tick $originM -K -O >> $psfile`;
`pscoast $bounds $proj $cinfo -S150/200/255 -W2 -G200/255/150 -N1/1p -N2/0.5p -O -K $originM >> $psfile`;

# plot southern California faults
if($ifault==1) {
   $dir1 = "/home/carltape/gmt";
   $fault_file = "$dir1/faults/jennings.xy";
   $kcf_file   = "$dir1/faults/kcf.xy";
   $breck_file = "$dir1/faults/breck.xy";
   $fault_infoK = "-M -W0.5p,0/0/0";
   `psxy $fault_file $bounds $proj $fault_infoK -O -K $originM >> $psfile`;
   `psxy $kcf_file $bounds $proj $fault_infoK -O -K $originM >> $psfile`;
   `psxy $breck_file $bounds $proj $fault_infoK -O -K $originM >> $psfile`;
}

# plot the stations
`awk \'{print \$4,\$3,0,4,B,C}\' $stafile | psxy $proj $bounds -M -St0.10 -W0.25p -G200/200/200 -N -O -K $originM >> $psfile`;

#------------

# plot labels next to the map
@mapstrings = ($edate,$eid,$stMw,$stedep,$strec,$stdist,$staz);
$nstrings = @mapstrings;
$x1 = $xmin - 4*$xtick;
for ($i=0; $i < $nstrings; $i++) {
  $y1 = $ymax - 0.5 - $i*$ytick*0.7;
  `pstext $bounds $proj $tick -N -K -O $originM >> $psfile <<EOF\n $x1 $y1 10 0 $fontno LM $mapstrings[$i] \nEOF\n`;
}

# plot ray path, station, and CMT source
`psxy $bounds $proj -W1p -K -O $originM >> $psfile <<EOF\n$elon $elat\n$slon $slat\nEOF\n`;  # ray path
#$cmtfont = 10; $cmt_scale = 0.3; $cmtinfo = "-Sm${cmt_scale}/$cmtfont -L0.5p/0/0/0 -G255/0/0";
#`psmeca $cmtpsmeca $bounds $proj $cmtinfo -O -K $originM >> $psfile`;  # CMT source
$cmt_scale = 0.3; $cmtinfo = "-Sm${cmt_scale} -L0.5p/0/0/0 -G255/0/0";
`psmeca $cmtpsmeca $bounds $proj $cmtinfo -O -K $originM >> $psfile`;  # CMT source
$elat_shift = $elat + 0.6;
`pstext $bounds $proj $textinfo -K -O $originM >> $psfile <<EOF\n $elon $elat_shift 9 0 $fontno CM $eid\nEOF\n`;  # station label
`rm $cmtpsmeca`;

`psxy $bounds $proj -St0.15 -W2 -G255/0/0 -K -O $originM >> $psfile <<EOF\n$slon $slat\nEOF\n`;  # red station
$slat_shift = $slat - 0.3;
`pstext $bounds $proj $textinfo -K -O $originM >> $psfile <<EOF\n $slon $slat_shift 8 0 $fontno CM $sta \nEOF\n`;  # station label

#---------------------------------------------------------------------------

`pstext -R0/1/0/1 -JX1 -N -Xa${x0_title} -Ya${y0_title} -O >> $psfile <<EOF\n 10 10 14 0 1 LM $title\nEOF\n`;  # FINISH (no -K)

#`convert $psfile $jpg_file`; `xv $jpg_file &`;
`ps2pdf $psfile`;
#`rm $psfile`;       # remove PS file
#system("xv $jpg_file &");

#================================================
sub get_cmt {
  my ($cmtfile)=@_;
  open(CMT,"$cmtfile") or die("Error opening cmtfile $cmtfile\n");
  @cmt = <CMT>;
  my($pde,$oyear,$omonth,$oday,$ohr,$omin,$osec1)=split(" ",$cmt[0]);
  my($osec,$omsec)=split(/\./,$osec1); $omsec=$omsec*10;
  my(undef,undef,$eid)=split(" ",$cmt[1]);
  my(undef,undef,$tshift)=split(" ",$cmt[2]);
  my(undef,undef,$hdur)=split(" ",$cmt[3]);
  my(undef,$elat)=split(" ",$cmt[4]);
  my(undef,$elon)=split(" ",$cmt[5]);
  my(undef,$edep)=split(" ",$cmt[6]);
  my(undef,$Mrr)=split(" ",$cmt[7]);
  my(undef,$Mtt)=split(" ",$cmt[8]);
  my(undef,$Mpp)=split(" ",$cmt[9]);
  my(undef,$Mrt)=split(" ",$cmt[10]);
  my(undef,$Mrp)=split(" ",$cmt[11]);
  my(undef,$Mtp)=split(" ",$cmt[12]);
  close(CMT);
  return ($oyear,$omonth,$oday,$ohr,$omin,$osec,$omsec,$tshift,$hdur,$elat,$elon,$edep,$Mrr,$Mtt,$Mpp,$Mrt,$Mrp,$Mtp,$eid);
}

#================================================
