#!/usr/bin/perl 
#
#-----------------------------------
# plot_win_adj.pl
#
# This script plots a single figure for a single source-receiver pair
# showing data, synthetics, windows, and adjoint sources.
#
# NOTE: USER must change the file suffixes below, depending on the dataset.
#       (Search for "USER".)
# 
# EXAMPLE:
#   plot_win_adj.pl -M ../CMTSOLUTION_9818433 -n MPM/CI/BH -b 0 -l -10/200 -k 7/1 -a STATIONS_ADJOINT -d DATA -s SYN -c RECON -w MEASUREMENT.WINDOWS -i m16 -j 6/30
#
#-----------------------------------

use Getopt::Std;
use POSIX;
use List::Util qw[min max];

sub Usage{
  print STDERR <<END;

  plot_win_adj.pl -M CMTFILE -n station/network -b iboth -r -a STATION_FILE -d data_dir -s syn_dir -c recon_dir -w winfile
  with

       -M -- CMTSOLUTION file to plot event focal mechanism (beach ball) on the map
       -l -- Start/End of trace to be cut
       -b -- plot MT and CC adjoint sources (=1)
       -k -- type of measurement (1-8) / plot adjoint source(0 or 1)
       -n -- station/network/chan
       -r -- plot 3-component seismograms and adjoint sources in relatively scaled amplitude
       -a -- station file to plot all the stations on the map
       -d -- data directory
       -s -- synthetics directory
       -c -- reconstructed directory
       -w -- MEASUREMENT.WINDOWS file
       -i -- label for the model iteration (for title)
       -j -- Tmin/Tmax (for title)
END
  exit(1);
}

if (@ARGV == 0) { Usage(); }
if (!getopts('k:m:l:n:b:ra:d:s:c:w:i:j:')) {die('Check input arguments\n');}
print "STATIONS file is $opt_a\n";
if ($opt_m) {@cmtvals = &get_cmt($opt_m);}
if ($opt_k) {($imeas,$iplot_adj) = split(/\//,$opt_k);}
if ($opt_l) {($tmin,$tmax) = split(/\//,$opt_l);}
if ($opt_w) {$Winfile=$opt_w;}
if ($opt_n) {($sta,$net,$chan)=split(/\//,$opt_n);}
if(!$opt_b) {$iboth = 0} else {$iboth = $opt_b}
if (!$opt_d) {$opt_d = "DATA"}
if (!$opt_s) {$opt_s = "SYN"}
if (!$opt_c) {$opt_c = "RECON"}
if (!$opt_i) {$opt_i = "m00"}
if ($opt_j) {($Tmin,$Tmax) = split(/\//,$opt_j);}
$smodel = $opt_i;

# suffix for adjoint sources
$klab = sprintf('iker%2.2i',$imeas);
$klab2 = sprintf('iker%2.2i',$imeas-2);
@ktitles = ("waveform, -d","waveform, s-d","banana-doughnut CC-TT","banana-doughnut dlnA",
            "cross-correlation TT","amplitude dlnA","multitaper TT","multitaper dlnA");
$ktitle = $ktitles[$imeas-1];

#if($imeas == 0)    {$klab = "ik"; $ktitle = "waveform";}
#elsif($imeas == 1) {$klab = "mtm"; $ktitle = "multitaper TT";}
#elsif($imeas == 2) {$klab = "cc"; $ktitle = "cross-correlation TT";}
#elsif($imeas == 3) {$klab = "bdcc"; $ktitle = "banana-doughnut CC-TT";}
#else {die("ERROR: imeas ($imeas) must be 0,1, 2, or 3\n");}

$sTmin = sprintf("T%3.3i",$Tmin);
$sTmax = sprintf("T%3.3i",$Tmax);
$Ttag = "${sTmin}_${sTmax}";
$Ttag_title = sprintf("bandpass %i s to %i s",$Tmin,$Tmax);
$Ttag_plot = sprintf("bp [%i s, %i s]",$Tmin,$Tmax);

#---------------------------
# ADDITIONAL USER INPUT

# name of file containing measurements (see mt_sub.f90)
$meas_file = "window_chi";
if (not -f $meas_file) {die("\n check if $meas_file exists\n")}

# plot faults (southern California)
$ifault = 1;

# with iboth = 1, you need to have made both the CC and MT adjoint sources
if($iboth==1) {
   if( ($imeas != 7) && ($imeas != 8) ) {die("\nFor plotting CC and MT adjoint sources, use imeas = 7 or 8.\n")}
   $ktitle = "TT -- MT = black ; CC = red";
}

# plotting preferences
#$iplot_adj = 1;      # plot adjoint sources
$iplot_win = 1;      # plot windows
$iplot_CClabel = 1;  # plot CC measurement labels for each window
$iplot_recon = 1;    # plot reconstructed synthetics
if($iplot_CClabel==1) {$iplot_win = 1;}

# no reconstructed records or labels for waveform difference misfit
if ($imeas < 5) {$iplot_recon = 0; $iplot_CClabel = 0;}

# adjust plot size based on whether you have adjoint sources
if($iplot_adj==1) {$swid = 3.5; $fsizetitle = 10;}
else {$swid = 7.4; $fsizetitle = 12;};

#---------------------------

# extract values from the array cmtvals
@M[0..5] = @cmtvals[12..17];
$year = $cmtvals[0]; $month = $cmtvals[1]; $day = $cmtvals[2];
$elat = $cmtvals[9]; $elon = $cmtvals[10]; $edep = $cmtvals[11];
$eid = $cmtvals[18];
open(MIJ,">temp.txt");
printf MIJ "$elon $elat $edep $M[0] $M[1] $M[2] $M[3] $M[4] $M[5] 1 X Y $eid";
close(MIJ);

# compute moment using Dahlen and Tromp, Eq. (5.91)
# compute moment magnitude using Kanamori (1977) -- see Carl's Latex notes cmt.pdf
$M0 = 1/sqrt(2) * sqrt($M[0]**2 + $M[1]**2 + $M[2]**2 + 2*($M[3]**2 + $M[4]**2 + $M[5]**2 ) );
$k  = -11.2 + (2/3)*log(5)/log(10);
$Mw = (2/3)*log($M0)/log(10) + $k;
$stM0 = sprintf("M0 = %.3e dyne-cm",$M0);
$stMw = sprintf("Mw = %.2f",$Mw);
print "\n           moment : $stM0 \n moment magnitude : $stMw\n";

$fontno = 1;  # font number to use for text

# set GMT defaults
`gmtset BASEMAP_TYPE plain PAPER_MEDIA letter MEASURE_UNIT inch ANNOT_FONT_SIZE = 12p LABEL_FONT_SIZE = 12p PAGE_ORIENTATION = portrait PLOT_DEGREE_FORMAT D`;

#$name = "${sta}_${net}_win_adj_${klab}";
$name = "${eid}_${Ttag}_${sta}_${net}_${smodel}_${klab}_win_adj";
$ps_file = "$name.ps";
$jpg_file = "$name.jpg";

#$saclst="/opt/seismo-util/bin/saclst";
$saclst="saclst";

#------------
$xmin = -122; $xmax = -114;
$ymin = 32; $ymax = 37;
$bounds="-R$xmin/$xmax/$ymin/$ymax ";
#$bounds="-R-120/-115/32/36.5";
$xtick = 1; $ytick = 1;
$tick="-B$xtick/$ytick";
$cinfo="-W0.5 -Dh -A100";
$proj="-JM4";

# plot the base map
$origin = "-X2.75";
`psbasemap $bounds $proj $tick $origin -K -P > $ps_file`;   # START
`pscoast $bounds $proj $cinfo -S150/200/255 -W2 -G200/255/150 -N1/1p -N2/0.5p -O -K >> $ps_file`;

# plot southern California faults
if($ifault==1) {
   $dir0 = "/opt/seismo/datalib/SC/";
   $fault_file = "$dir0/jennings.xy";
#   $kcf_file   = "$dir0/faults/kcf.xy";
#   $breck_file = "$dir0/faults/breck.xy";
   $fault_infoK = "-m -W0.5p,0/0/0";
   `psxy $fault_file $bounds $proj $fault_infoK -O -K >> $ps_file`;
#   `psxy $kcf_file $bounds $proj $fault_infoK -O -K >> $ps_file`;
#   `psxy $breck_file $bounds $proj $fault_infoK -O -K >> $ps_file`;
}

# plot the stations
`awk \'{print \$4,\$3,0,4,B,C}\' $opt_a | psxy $proj $bounds -m -St0.15 -W2 -G200/200/200 -N -O -K >>$ps_file`;

#------------
# data, synthetics, and (optional) reconstructed synthetics

$suffix = ".d";     # USER MUST CHANGE

# DATA
$d_tag = "${opt_d}/$eid.$net.$sta";
$Tdat = `ls ${d_tag}.*T.sac*`; chomp($Tdat);
$Rdat = `ls ${d_tag}.*R.sac*`; chomp($Rdat);
$Zdat = `ls ${d_tag}.*Z.sac*`; chomp($Zdat);
#$Tdat = "$opt_d/$eid.$net.$sta.BHT.sac${suffix}";
#$Rdat = "$opt_d/$eid.$net.$sta.BHR.sac${suffix}";
#$Zdat = "$opt_d/$eid.$net.$sta.BHZ.sac${suffix}";
if (0==1) {
  ($Tdat,$Tchan) = split(" ",`$saclst kcmpnm f ${d_tag}*T.sac${suffix}`);
  ($Rdat,$Rchan) = split(" ",`$saclst kcmpnm f ${d_tag}*R.sac${suffix}`);
  ($Zdat,$Zchan) = split(" ",`$saclst kcmpnm f ${d_tag}*Z.sac${suffix}`);
  print "\n $Tdat $Tchan \n $Rdat $Rchan \n $Zdat $Zchan \n";
}
if(-f $Tdat) {$Tflag = 1} else {$Tflag = 0}
if(-f $Rdat) {$Rflag = 1} else {$Rflag = 0}
if(-f $Zdat) {$Zflag = 1} else {$Zflag = 0}
#print "\n Tflag : $Tflag, Rflag : $Rflag, Zflag : $Zflag \n";

# SYNTHETICS -- component label always has the form BH_ or LH_,
# even if the data file is something else (HH_, LH_).
# We assume that ALL components exist.
$Tlab = "$sta.$net.${chan}T";
$Rlab = "$sta.$net.${chan}R";
$Zlab = "$sta.$net.${chan}Z";
$Tsyn = `ls ${opt_s}/${Tlab}*`; chomp($Tsyn);
$Rsyn = `ls ${opt_s}/${Rlab}*`; chomp($Rsyn);
$Zsyn = `ls ${opt_s}/${Zlab}*`; chomp($Zsyn);
#$Tsyn = "${opt_s}/${Tlab}.semd.sac${suffix}";
#$Rsyn = "${opt_s}/${Rlab}.semd.sac${suffix}";
#$Zsyn = "${opt_s}/${Zlab}.semd.sac${suffix}";

# RECONSTRUCTED SYNTHETICS
$suffix = ".recon.sac";     # USER MUST CHANGE
$Trecon = "${opt_c}/${Tlab}${suffix}";
$Rrecon = "${opt_c}/${Rlab}${suffix}";
$Zrecon = "${opt_c}/${Zlab}${suffix}";

#------------

#print "\n checking files: \n -- $Tdat -- $Rdat -- $Zdat -- \n -- $Tsyn -- $Rsyn -- $Zsyn -- \n -- $Trecon -- $Rrecon -- $Zrecon -- \n"; die("testing");

# get the receiver location, azimuth, distance (from the $Zsyn file)
(undef,$tmin0,$tmax0,$slon,$slat,$dist,$az,$Tchan) = split(" ",`$saclst b e stlo stla dist az kcmpnm f $Zsyn`);
$edate  = sprintf("%4.4i . %2.2i . %2.2i",$year,$month,$day);
$stedep = sprintf("depth = %.2f km",$edep);
$strec  = "$sta . $net";
$stdist = sprintf("dist = %.1f km",$dist);
$staz   = sprintf("az = %.1f deg",$az);
$stmodel  = "model $smodel";

# if the time cut is not specified, then use the whole record
if (!$opt_l) {$tmin = $tmin0; $tmax = $tmax0;}
$trange = $tmax - $tmin;

# plot labels next to the map
@mapstrings = ($edate,$eid,$stMw,$stedep,$strec,$stdist,$staz,"--",${Ttag_plot},$stmodel);
$nstrings = @mapstrings;
$x1 = $xmin - 4*$xtick;
for ($i=0; $i < $nstrings; $i++) {
  $y1 = $ymax - $i*$ytick/2;
  `pstext $bounds $proj $tick -N -K -O >> $ps_file <<EOF\n $x1 $y1 12 0 $fontno LM $mapstrings[$i] \nEOF\n`;
}

#`pstext $bounds $proj $tick -N -K -O >> $ps_file <<EOF\n $x1 $y0 12 0 $fontno LM $eid \nEOF\n`;
#`pstext $bounds $proj $tick -N -K -O >> $ps_file <<EOF\n $x1 $y1 12 0 $fontno LM $stMw \nEOF\n`;
#`pstext $bounds $proj $tick -N -K -O >> $ps_file <<EOF\n $x1 $y2 12 0 $fontno LM $stedep \nEOF\n`;
#`pstext $bounds $proj $tick -N -K -O >> $ps_file <<EOF\n $x1 $y3 12 0 $fontno LM $strec \nEOF\n`;
#`pstext $bounds $proj $tick -N -K -O >> $ps_file <<EOF\n $x1 $y4 12 0 $fontno LM $stdist \nEOF\n`;
#`pstext $bounds $proj $tick -N -K -O >> $ps_file <<EOF\n $x1 $y5 12 0 $fontno LM $staz \nEOF\n`;

# plot ray path, station, and CMT source
`psxy $bounds $proj -W1p -K -O >> $ps_file <<EOF\n$elon $elat\n$slon $slat\nEOF\n`;  # ray path
$cmtfont = 10; $cmt_scale = 0.4; $cmtinfo = "-Sm${cmt_scale}/$cmtfont -L0.5p/0/0/0 -G255/0/0";
`psmeca temp.txt $bounds $proj $cmtinfo -O -K >> $ps_file`;  # CMT source
`rm temp.txt`;
`psxy $bounds $proj -St0.15 -W2 -G255/0/0 -K -O>> $ps_file <<EOF\n$slon $slat\nEOF\n`;  # red station
$slat_shift=$slat+0.2;
`pstext $bounds $proj -K -O>> $ps_file <<EOF\n $slon $slat_shift 8 0 $fontno CM $sta \nEOF\n`;  # station label

#-----------------------------------
# limits for seismogram axes
$Smin = 0;
$Smax = 2;
$Srange = $Smax - $Smin;

$xtextpos1 = $tmin - 0.03*$trange;      # xmax label
$xtextpos2 = $tmin + 0.05*$trange;      # first window measurement label
$xtextpos3 = $tmin + 0.50*$trange;      # title position
$xtextpos4 = $tmin + 0.95*$trange;      # final window measurement label

$ytextpos5 = 1.35*$Srange;    # titles
$ytextpos4 = 1.03*$Srange;    # dT label
$ytextpos3 = 0.90*$Srange;    # dlnA label
$ytextpos2 = 0.50*$Srange;    # Z, R, T label
$ytextpos1 = 0.15*$Srange;    # ymax label

# pstext info for measurement labels
$textinfo = "-C2p -W0/255/255o,0/255/255 -G0/0/0";

# dimension of seismograms
$proj = "-JX${swid}/1";
$dX = $swid + 0.5;
$dY = 1.5;

$tick = "-Ba50f10:\"Time (s)\":/a1f1S";
$tick2 = "-Ba50f10/a1f1S";

$pen = "0.7p";      # line thickness for seismograms
$red_pen = "-W${pen},250/0/0"; # syn
$black_pen = "-W${pen},0/0/0"; # data
$recon_pen = "-W${pen},0/255/255"; # recon_syn

my(@Twinb); my(@Twine);
my(@Rwinb); my(@Rwine);
my(@Zwinb); my(@Zwine);
my($undef);

# windows and measurements: transverse component
$ntline=`grep -n \"$Tsyn\" $Winfile | awk -F: '{ print \$1 }'`;
chomp($ntline);
if ($ntline) {
  $nline=$ntline;
  $nline=$nline+1;
  $nT=`sed -n ${nline}p $Winfile`;
  @mlines = `grep $Tlab ${meas_file}`; # measurements
  for ($i=0; $i < $nT; $i++) {
    $nline=$nline+1;
    $win=`sed -n ${nline}p $Winfile`; chomp($win);
    ($winb,$wine)=split(" ",$win);
    push @Twinb, "$winb";
    push @Twine, "$wine";

    (undef,undef,undef,undef,undef,$kplotT,undef,undef,$chi2T,$chi3T,$chi4T,$chi5T,$meas2T,$meas3T,$meas4T,$meas5T,$sigma2T,$sigma3T,$sigma4T,$sigma5T) = split(" ",$mlines[$i]);
    $stlabT_dT[$i] = sprintf("dT = %.2f \\261 %.2f",$meas4T,$sigma4T);
    $stlabT_dA[$i] = sprintf("dA = %.2f \\261 %.2f",$meas5T,$sigma5T);
  }
}

# windows and measurements: radial component
$nrline=`grep -n \"$Rsyn\" $Winfile | awk -F: '{ print \$1 }'`;
chomp($nrline);
#print "$nrline\n";
if ($nrline) {
  $nline=$nrline;
  $nline=$nline+1;
  $nR=`sed -n ${nline}p $Winfile`;
  @mlines = `grep $Rlab ${meas_file}`; # measurements
  for ($i=0; $i < $nR; $i++) {
    $nline=$nline+1;
    $win=`sed -n ${nline}p $Winfile`;
    chomp($win);
    #  print "$win\n";
    ($winb,$wine)=split(" ",$win);
    #  print "$winb\n";
    #  print "$wine\n";
    push @Rwinb, "$winb";
    push @Rwine, "$wine";

    (undef,undef,undef,undef,undef,$kplotR,undef,undef,$chi2R,$chi3R,$chi4R,$chi5R,$meas2R,$meas3R,$meas4R,$meas5R,$sigma2R,$sigma3R,$sigma4R,$sigma5R) = split(" ",$mlines[$i]);
    $stlabR_dT[$i] = sprintf("dT = %.2f \\261 %.2f",$meas4R,$sigma4R);
    $stlabR_dA[$i] = sprintf("dA = %.2f \\261 %.2f",$meas5R,$sigma5R);
  }
}

# windows and measurements: vertical component
$nzline=`grep -n \"$Zsyn\" $Winfile | awk -F: '{ print \$1 }'`;
chomp($nzline);
if ($nzline) {
  $nline=$nzline;
  $nline=$nline+1;
  $nZ=`sed -n ${nline}p $Winfile`;
  @mlines = `grep $Zlab ${meas_file}`; # measurements

  for ($i=0; $i < $nZ; $i++) {
    $nline=$nline+1;
    $win=`sed -n ${nline}p $Winfile`;
    chomp($win);
     # print "$win\n";
    ($winb,$wine)=split(" ",$win);
     # print "$winb\n";
     # print "$wine\n";
    push @Zwinb, "$winb";
    push @Zwine, "$wine";

    (undef,undef,undef,undef,undef,$kplotZ,undef,undef,$chi2Z,$chi3Z,$chi4Z,$chi5Z,$meas2Z,$meas3Z,$meas4Z,$meas5Z,$sigma2Z,$sigma3Z,$sigma4Z,$sigma5Z) = split(" ",$mlines[$i]);
    $stlabZ_dT[$i] = sprintf("dT = %.2f \\261 %.2f",$meas4Z,$sigma4Z);
    $stlabZ_dA[$i] = sprintf("dA = %.2f \\261 %.2f",$meas5Z,$sigma5Z);
  }
}

#==========================================

# pssac2 scaling
#      -M vertical scaling in sacfile_unit/MEASURE_UNIT = size<required>
#          size: each trace will normalized to size (in MEASURE_UNIT)
#              scale =  (1/size) * [data(max) - data(min)]

# limits for data records
$maxT = 0; $maxR = 0; $maxZ = 0;
if($Tflag==1) {(undef,$minT,$maxT)=split(" ",`$saclst depmin depmax f $Tdat`); $maxT=max(abs($minT),abs($maxT));}
if($Rflag==1) {(undef,$minR,$maxR)=split(" ",`$saclst depmin depmax f $Rdat`); $maxR=max(abs($minR),abs($maxR));}
if($Zflag==1) {(undef,$minZ,$maxZ)=split(" ",`$saclst depmin depmax f $Zdat`); $maxZ=max(abs($minZ),abs($maxZ));}

# limits for synthetic records
(undef,$minT2,$maxT2) = split(" ",`$saclst depmin depmax f $Tsyn`);
$maxT2 = max(abs($minT2),abs($maxT2));
(undef,$minR2,$maxR2) = split(" ",`$saclst depmin depmax f $Rsyn`);
$maxR2 = max(abs($minR2),abs($maxR2));
(undef,$minZ2,$maxZ2) = split(" ",`$saclst depmin depmax f $Zsyn`);
$maxZ2 = max(abs($minZ2),abs($maxZ2));

# overall max values
$maxTDS = max($maxT,$maxT2);
$maxRDS = max($maxR,$maxR2);
$maxZDS = max($maxZ,$maxZ2);
$maxD = max($maxT,max($maxR,$maxZ));
$maxS = max($maxT2,max($maxR2,$maxZ2));
$max = max($maxD,$maxS);
if ($maxTDS == 0) {$maxTDS=$maxD;}
if ($maxRDS == 0) {$maxRDS=$maxD;}
if ($maxZDS == 0) {$maxZDS=$maxD;}

#print "\n-- $Tsyn -- $Tdat\n -- $Rsyn -- $Rdat\n -- $Zsyn -- $Zdat\n";
#print "\n Overall max values: $maxD (data), $maxS (syn) $max (overall)\n";

$bounds = "-R$tmin/$tmax/$Smin/$Smax";
$scale = 1.5;       # scaling for plotting seismograms
$size = 1/$scale;

if ($opt_r) {
  # normalize by the same value (max)
  $sizeT  = $size*$maxT/$max;
  $sizeT2 = $size*$maxT2/$max;
  $sizeR  = $size*$maxR/$max;
  $sizeR2 = $size*$maxR2/$max;
  $sizeZ  = $size*$maxZ/$max;
  $sizeZ2 = $size*$maxZ2/$max;

} else {
  # normalize by different values
  $sizeT  = $size*$maxT/$maxTDS;
  $sizeT2 = $size*$maxT2/$maxTDS;
  $sizeR  = $size*$maxR/$maxRDS;
  $sizeR2 = $size*$maxR2/$maxRDS;
  $sizeZ  = $size*$maxZ/$maxZDS;
  $sizeZ2 = $size*$maxZ2/$maxZDS;
}
$maxsize=max($sizeT,$sizeR,$sizeZ);
$maxsize2=max($sizeT2,$sizeR2,$sizeZ2);

if ($sizeT == 0) {$sizeT = $maxsize;}
if ($sizeR == 0) {$sizeR = $maxsize;}
if ($sizeZ == 0) {$sizeZ = $maxsize;}
if ($sizeT2 == 0) {$sizeT2 = $maxsize2;}
if ($sizeR2 == 0) {$sizeR2 = $maxsize2;}
if ($sizeZ2 == 0) {$sizeZ2 = $maxsize2;}

# limits for reconstructed records
if ($Tflag && $ntline) {
  (undef,$minT3,$maxT3) = split(" ",`$saclst depmin depmax f $Trecon`);
  $maxT3 = max(abs($minT3),abs($maxT3));
  if ($opt_r) {
    $sizeT3 = $size*$maxT3/$max;
  } else {
    $sizeT3 = $size*$maxT3/$maxTDS;
  }
}
if ($Rflag && $nrline) {
  (undef,$minR3,$maxR3) = split(" ",`$saclst depmin depmax f $Rrecon`);
  $maxR3 = max(abs($minR3),abs($maxR3));
   if ($opt_r) {
    $sizeR3 = $size*$maxR3/$max;
  } else {
    $sizeR3 = $size*$maxR3/$maxRDS;
  }
}
if ($Zflag && $nzline) {
  (undef,$minZ3,$maxZ3) = split(" ",`$saclst depmin depmax f $Zrecon`);
  $maxZ3 = max(abs($minZ3),abs($maxZ3));
  if ($opt_r) {
    $sizeZ3 = $size*$maxZ3/$max;
  } else {
    $sizeZ3 = $size*$maxZ3/$maxZDS;
  }
}

# strings for plotting -- show the yaxis limits for each trace (based on synthetics)
$strT = sprintf('%.2e',$maxT2/$sizeT2);
$strR = sprintf('%.2e',$maxR2/$sizeR2);
$strZ = sprintf('%.2e',$maxZ2/$sizeZ2);

#==========================================

# TRANSVERSE component: data, synthetics, and windows
#print "pssac2 $Tsyn -X-2 -Y4.5 $proj $bounds -Ent-3 -M${sizeT2} $red_pen -N -K -O $tick \n";
if (-f $Tsyn){
`pssac2 $Tsyn -X-2 -Y4.5 $proj $bounds -Ent-3 -M${sizeT2} $red_pen -N -K -O $tick >> $ps_file`;	# synthetics
if ($ntline) {
  if ($iplot_win==1) {
    for ($i=0; $i<$nT; $i++) {
      open(GMT,"|psxy $proj $bounds -W2 -G220/220/220 -O -K >> $ps_file");
      print GMT "$Twinb[$i] 0\n $Twine[$i] 0\n $Twine[$i] 2\n $Twinb[$i] 2\n";
      #print GMT "$Twinb[$i] 0\n";
      #print GMT "$Twine[$i] 0\n";
      #print GMT "$Twine[$i] 2\n";
      #print GMT "$Twinb[$i] 2\n";
      close(GMT);
    }
  }
  if ($iplot_CClabel) {
    for ($i=0; $i<$nT; $i++) {
      $xtext = $Twinb[$i]; $just = "LB";
      if ($i == 0) {
	$xtext = $xtextpos2; $just = "LB";
      }
      if ($i == $nT-1) {
	$xtext = $xtextpos4; $just = "RB";
      }
      `echo \"$xtext $ytextpos4 7 0 1 $just $stlabT_dT[$i]\" | pstext $textinfo $proj $bounds -N -O -K >> $ps_file`;
      `echo \"$xtext $ytextpos3 7 0 1 $just $stlabT_dA[$i]\" | pstext $textinfo $proj $bounds -N -O -K >> $ps_file`;
    }
  }
}}
if ($Tflag==1) {`pssac2 $Tdat $proj $bounds -Ent-3 -M${sizeT} $black_pen -N -K -O >> $ps_file`;}    # data
if (-f $Tsyn){
`pssac2 $Tsyn $proj $bounds -Ent-3 -M${sizeT2} $red_pen -N -K -O >> $ps_file`;}                     # synthetics
if ($Tflag && $ntline && ($iplot_recon==1)) {`pssac2 $Trecon $proj $bounds -Ent-3 -M${sizeT3} $recon_pen -N -K -O >> $ps_file`;}   # reconstructed
`echo \"$xtextpos1 $ytextpos2 14 0 1 BC T\" | pstext $proj $bounds -N -O -K >> $ps_file`;
`echo \"$xtextpos1 $ytextpos1 9 0 1 CM $strT\" | pstext $proj $bounds -N -O -K >> $ps_file`;

# RADIAL component: data, synthetics, and windows
if (-f $Rsyn) {
`pssac2 $Rsyn -Y$dY $proj $bounds -Ent-3 -M${sizeR2} $red_pen -N -K -O $tick2 >> $ps_file`;}
if ($nrline) {
  if ($iplot_win==1) {
    for ($i=0; $i<$nR; $i++) {
      open(GMT,"|psxy $proj $bounds -W2 -G220/220/220 -O -K >> $ps_file");
      print GMT "$Rwinb[$i] 0\n $Rwine[$i] 0\n $Rwine[$i] 2\n $Rwinb[$i] 2\n";
      close(GMT);
    }
  }
  if ($iplot_CClabel==1) {
    for ($i=0; $i<$nR; $i++) {
      $xtext = $Rwinb[$i]; $just = "LB";
      if ($i == 0) {
	$xtext = $xtextpos2; $just = "LB";
      }
      if ($i == $nR-1) {
	$xtext = $xtextpos4; $just = "RB";
      }
      `echo \"$xtext $ytextpos4 7 0 1 $just $stlabR_dT[$i]\" | pstext $textinfo $proj $bounds -N -O -K >> $ps_file`;
      `echo \"$xtext $ytextpos3 7 0 1 $just $stlabR_dA[$i]\" | pstext $textinfo $proj $bounds -N -O -K >> $ps_file`;
    }
  }
}
if ($Rflag==1) {`pssac2 $Rdat $proj $bounds -Ent-3 -M${sizeR} $black_pen -N -K -O >> $ps_file`;}
if (-f $Rsyn){
`pssac2 $Rsyn $proj $bounds -Ent-3 -M${sizeR2} $red_pen -N -K -O >> $ps_file`;}
if ($Rflag && $nrline && ($iplot_recon==1)) {`pssac2 $Rrecon $proj $bounds -Ent-3 -M${sizeR3} $recon_pen -N -K -O >> $ps_file`;}
`echo \"$xtextpos1 $ytextpos2 14 0 1 BC R\" | pstext $proj $bounds -N -O -K >> $ps_file`;
`echo \"$xtextpos1 $ytextpos1 9 0 1 CM $strR\" | pstext $proj $bounds -N -O -K >> $ps_file`;

# VERTICAL component: data, synthetics, and windows
if (-f $Zsyn) {
`pssac2 $Zsyn -Y$dY $proj $bounds -Ent-3 -M${sizeZ2} $red_pen -N -K -O  $tick2 >> $ps_file`;
if ($nzline) {
  if ($iplot_win==1) {
    for ($i=0; $i<$nZ; $i++) {
      open(GMT,"|psxy $proj $bounds -W2 -G220/220/220 -O -K >> $ps_file");
      print GMT "$Zwinb[$i] 0\n $Zwine[$i] 0\n $Zwine[$i] 2\n $Zwinb[$i] 2\n";
      close(GMT);
      #print "*** $Zwinb[$i] $Zwine[$i]\n";
    }
  }
  if ($iplot_CClabel==1) {
    for ($i=0; $i<$nZ; $i++) {
      $xtext = $Zwinb[$i]; $just = "LB";
      if ($i == 0) {
	$xtext = $xtextpos2; $just = "LB";
      }
      if ($i == $nZ-1) {
	$xtext = $xtextpos4; $just = "RB";
      }
      `echo \"$xtext $ytextpos4 7 0 1 $just $stlabZ_dT[$i]\" | pstext $textinfo $proj $bounds -N -O -K >> $ps_file`;
      `echo \"$xtext $ytextpos3 7 0 1 $just $stlabZ_dA[$i]\" | pstext $textinfo $proj $bounds -N -O -K >> $ps_file`;
    }
  }
}}
if ($Zflag==1) {`pssac2 $Zdat $proj $bounds -Ent-3 -M${sizeZ}  $black_pen -N -K -O $tick2 >> $ps_file`;}
if (-f $Zsyn){
`pssac2 $Zsyn $proj $bounds -Ent-3 -M${sizeZ2} $red_pen -N -O -K >> $ps_file`;}
if ($Zflag && $nzline && ($iplot_recon==1)) {`pssac2 $Zrecon $proj $bounds -Ent-3 -M${sizeZ3} $recon_pen -N -K -O >> $ps_file`;}
`echo \"$xtextpos1 $ytextpos2 14 0 1 BC Z\" | pstext $proj $bounds -N -K -O >> $ps_file`;
`echo \"$xtextpos1 $ytextpos1 9 0 1 CM $strZ\" | pstext $proj $bounds -N -O -K >> $ps_file`;

# plot titles
`echo \"$xtextpos3 $ytextpos5 $fsizetitle 0 $fontno BC Data (Black) and Synthetics (Red) -- ${Ttag_title}\" | pstext $proj $bounds -N -K -O >> $ps_file`;
if ($iplot_adj==1) {
   `echo \"$xtextpos3 $ytextpos5 $fsizetitle 0 $fontno BC Adjoint Source ($ktitle)\" | pstext $proj $bounds -N -X$dX -K -O >> $ps_file`;
}

#----------------------------------------
# PLOTTING ADJOINT SOURCES

if($iplot_adj == 1) {

  $tick="-Ba50f10:\"Time (s)\":/a1f1S";
  $tick2="-Ba50f10/a1f1S";
  $proj="-JX${swid}/1";
  #$red_pen = "-W1p,250/0/0";
  #$black_pen = "-W1p,0/0/0";

  # adjoint sources (always have the BH_ or LH_ channel, even if the data do not)
  $Zadj="ADJOINT_SOURCES/${sta}.${net}.${chan}Z.${klab}.adj";
  $Radj="ADJOINT_SOURCES/${sta}.${net}.${chan}R.${klab}.adj";
  $Tadj="ADJOINT_SOURCES/${sta}.${net}.${chan}T.${klab}.adj";

  #$maxZ=1.0; $maxR=1.0; $maxT=1.0;

  $scale=1.2;			# adjust for plotting (this times max gives the axis limits)
  if (-f $Zadj) {
    ($minZ,$maxZ)=split(/\//,`minmax $Zadj | awk '{print \$6}'`);
    (undef,$minZ)=split(/</,$minZ);($maxZ,undef)=split(/>/,$maxZ);
    $maxZ=$scale*max(abs($maxZ),abs($minZ));
  }
  if (-f $Radj) {
    ($minR,$maxR)=split(/\//,`minmax $Radj | awk '{print \$6}'`);
    (undef,$minR)=split(/</,$minR);($maxR,undef)=split(/>/,$maxR);
    $maxR=$scale*max(abs($maxR),abs($minR));
  }
  if (-f $Tadj) {
    ($minT,$maxT)=split(/\//,`minmax $Tadj | awk '{print \$6}'`);
    (undef,$minT)=split(/</,$minT);($maxT,undef)=split(/>/,$maxT);
    $maxT=$scale*max(abs($maxT),abs($minT));
  }

  # overall max value
  $maxA = max($maxT,max($maxR,$maxZ));
  #print "\n Overall adjoint source max value: $maxA\n";
  if ($opt_r) {
    $maxZ=$maxA;$maxR=$maxA;$maxT=$maxA;
  }

  # CHT: 10e-3 tolerance too large for waveform difference
  $tol = 10e-8;
  if ($maxZ<$tol) {
    $maxZ=1.0;
  }
  if ($maxR<$tol) {
    $maxR=1.0;
  }
  if ($maxT<$tol) {
    $maxT=1.0;
  }

  #$maxZ=1.0; $maxR=1.0; $maxT=1.0;    # testing

  # strings for plotting
  $strT = sprintf('%.2e',$maxT);
  $strR = sprintf('%.2e',$maxR);
  $strZ = sprintf('%.2e',$maxZ);

  # adjoint sources for cross-correlation plot in red
  if ($iboth==1) {
    $Zadj2="ADJOINT_SOURCES/$sta.$net.${chan}Z.${klab2}.adj";
    $Radj2="ADJOINT_SOURCES/$sta.$net.${chan}R.${klab2}.adj";
    $Tadj2="ADJOINT_SOURCES/$sta.$net.${chan}T.${klab2}.adj";
    
  }

  # VERTICAL component: adjoint source and windows
  $bounds="-R$tmin/$tmax/-$maxZ/$maxZ";
  `psbasemap $proj $bounds $tick2 -K -O >> $ps_file`;
  if ($iplot_win==1) {
    if ($nzline) {
      for ($i =0; $i<$nZ; $i++) {
	open(GMT,"|psxy $proj $bounds -W2 -G220/220/220 -O -K >> $ps_file");
	print GMT "$Zwinb[$i] -$maxZ\n";
	print GMT "$Zwine[$i] -$maxZ\n";
	print GMT "$Zwine[$i] $maxZ\n";
	print GMT "$Zwinb[$i] $maxZ\n";
	close(GMT);
      }
    }
  }
  if (-f $Zadj) {
    $ytextpos2 = 0; $ytextpos1 = -$maxZ + 0.15*(2*$maxZ);
    `psxy $proj $Zadj $bounds $black_pen -N -K -O >> $ps_file`;
    if ($iboth==1) {
      `psxy $proj $Zadj2 $bounds $red_pen -N -K -O >> $ps_file`;
    }
    `echo \"$xtextpos1 $ytextpos2 14 0 1 BC Z\" | pstext $proj $bounds -N -K -O >> $ps_file`;
    `echo \"$xtextpos1 $ytextpos1 9 0 $fontno CM $strZ\" | pstext $proj $bounds -N -K -O >> $ps_file`;
  }

  # RADIAL component: adjoint source and windows
  $bounds="-R$tmin/$tmax/-$maxR/$maxR";
  `psbasemap -Y-$dY $proj $bounds $tick2 -K -O >> $ps_file`;
  if ($iplot_win==1) {
    if ($nrline) {
      for ($i =0; $i<$nR; $i++) {
	open(GMT,"|psxy $proj $bounds -W2 -G220/220/220 -O -K >> $ps_file");
	print GMT "$Rwinb[$i] -$maxR\n";
	print GMT "$Rwine[$i] -$maxR\n";
	print GMT "$Rwine[$i] $maxR\n";
	print GMT "$Rwinb[$i] $maxR\n";
	close(GMT);
      }
    }
  }
  if (-f $Radj) {
    $ytextpos2 = 0; $ytextpos1 = -$maxR + 0.15*(2*$maxR);
    `psxy $proj $Radj $bounds $black_pen -N -K -O >> $ps_file`;
    `echo \"$xtextpos1 $ytextpos2 14 0 1 BC R\" | pstext $proj $bounds -N -K -O >> $ps_file`;
    if ($iboth==1) {
      `psxy $proj $Radj2 $bounds $red_pen -N -K -O >> $ps_file`;
    }
    `echo \"$xtextpos1 $ytextpos1 9 0 $fontno CM $strR\" | pstext $proj $bounds -N -K -O >> $ps_file`;
  }

  # TRANSVESE component: adjoint source and windows
  $bounds="-R$tmin/$tmax/-$maxT/$maxT";
  `psbasemap -Y-$dY $proj $bounds $tick -K -O >> $ps_file`;
  if ($iplot_win==1) {
    if ($ntline) {
      for ($i =0; $i<$nT; $i++) {
	open(GMT,"|psxy $proj $bounds -W2 -G220/220/220 -O -K >> $ps_file");
	print GMT "$Twinb[$i] -$maxT\n";
	print GMT "$Twine[$i] -$maxT\n";
	print GMT "$Twine[$i] $maxT\n";
	print GMT "$Twinb[$i] $maxT\n";
	close(GMT);
      }
    }
  }
  if (-f $Tadj) {
    $ytextpos2 = 0; $ytextpos1 = -$maxT + 0.15*(2*$maxT);
    `psxy $proj $Tadj $bounds $black_pen -N -K -O >> $ps_file`;
    if ($iboth==1) {
      `psxy $proj $Tadj2 $bounds $red_pen -N -K -O >> $ps_file`;
    }
    `echo \"$xtextpos1 $ytextpos2 14 0 1 BC T\" | pstext $proj $bounds -N -K -O >> $ps_file`;
    `echo \"$xtextpos1 $ytextpos1 9 0 $fontno CM $strT\" | pstext $proj $bounds -N -K -O >> $ps_file`;
  }

}

#-----------------------

# plot titles
#`echo \"$xtextpos3 $maxT 12 0 $fontno BC Adjoint Source ($ktitle)\" | pstext $proj $bounds -N -Y3.2 -K -O >> $ps_file`;
#`echo \"$xtextpos3 $maxT 12 0 $fontno BC Data and Synthetics\" | pstext $proj $bounds -N -O -X-4 >> $ps_file`;
`echo \"0 0\" | psxy $proj $bounds -N -O -X15 >> $ps_file`;   # FINISH

#`convert $ps_file $jpg_file ; xv $jpg_file &`;
`ps2pdf $ps_file`;
`rm $ps_file`;       # remove PS file
#system("xv $jpg_file &");

#================================================

# subroutine name: max
# Input: number1, number2
# returns greater of 2 numbers
#sub max {
#	if ($_[0]<$_[1]) {return $_[1]} else {return $_[0]};
#}

sub get_cmt {
  my ($cmt_file)=@_;
  open(CMT, "$cmt_file") or die("Error opening $cmt_file\n");
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
