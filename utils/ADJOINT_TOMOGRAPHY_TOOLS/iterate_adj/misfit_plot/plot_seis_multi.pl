#!/usr/bin/perl -w

#-----------------------------------
# Carl Tape, 17-Nov-2008
# plot_seis_multi.pl
#
# Text here.
#
# 
# EXAMPLE:
#
# plot_seis_multi.pl -m CMTSOLUTION_3321590 -n PHL/CI -l 10/190 -a 3321590_T006_T030_m16_STATIONS_ADJOINT -d DATA -s SYN -w MEASUREMENT_WINDOWS_3321590_T006_T030_m16/MEASUREMENT_WINDOWS_3321590_T006_T030_m00 -x 3321590_T006_T030_m16_window_chi/3321590_T006_T030_m00_window_chi -i 3321590/m16/m00/6/30
#
# plot_seis_multi.pl -m CMTSOLUTION_3321590 -n PHL/CI -l 10/190 -a 3321590_T006_T030_m16_STATIONS_ADJOINT -d DATA -s SYN -i 3321590/m16/m00/6/30 MEASUREMENT_WINDOWS_3321590_T006_T030_m16 MEASUREMENT_WINDOWS_3321590_T006_T030_m00 3321590_T006_T030_m16_window_chi 3321590_T006_T030_m00_window_chi
#
#-----------------------------------

use Getopt::Std;
use POSIX;

sub Usage{
  print STDERR <<END;

  plot_seis_multi.pl -m CMTFILE -n station/network -r -d data_dir -s syn_dir -w winfile
  with

       -m -- CMTSOLUTION file to plot event focal mechanism (beach ball) on the map
       -l -- Start/End of trace to be cut
       -n -- station/network
       -r -- plot 3-compoment seismograms and adjoint sources in relatively scaled amplitude
       -d -- data directory
       -s -- synthetics directory
       -w -- window file for model 1 / window file for model 2
       -x -- measurement file for model 1 / measurement file for model 2
       -i -- event ID / Tmin / Tmax
END
  exit(1);
}

if (@ARGV == 0) { Usage(); }
if (!getopts('m:l:n:r:d:s:i:')) {die('Check input arguments\n');}
if ($opt_m) {@cmtvals = &get_cmt($opt_m);}
if ($opt_l) {($tmin,$tmax) = split(/\//,$opt_l);$trange = $tmax-$tmin;}   # search for opt_l
#if ($opt_w) {($win_file1,$win_file2) = split(/\//,$opt_w);}
#if ($opt_x) {($meas_file1,$meas_file2) = split(/\//,$opt_x);}
if ($opt_n) {($sta,$net)=split(/\//,$opt_n);}
#if (!$opt_a) {$opt_a = "STATIONS_ADJOINT"}
if (!$opt_d) {$opt_d = "DATA"}
if (!$opt_s) {$opt_s = "SYN"}
if ($opt_i) {($eid,$Tmin,$Tmax) = split(/\//,$opt_i);}

# final four files are: window file 1, window file 2, measurement file 1, measurement file 2
if(@ARGV < 3) {die("EXIT: must have at least three files -- smodel, window file, measurement file");}
@allfiles = @ARGV;
$nfiles = @allfiles;
$nseis = $nfiles/3;     # number of models for which we want to compare seismograms
#print "-- nfiles -- $nfiles -- nseis -- $nseis --\n";
if($nfiles % 3 != 0) {die("files must be a factor of three\n")}
for ($kk = 1; $kk <= $nseis; $kk++) {
   $smodels[$kk-1] = $allfiles[$kk - 1];
   $winfiles[$kk-1] = $allfiles[$kk+$nseis-1];
   $measfiles[$kk-1] = $allfiles[$kk+2*$nseis-1];
}
#print "-- @smodels --\n-- @winfiles --\n-- @measfiles --\n";

# check that at least one set of input files is there
$smodel1 = $smodels[0];
$win_file1 = $winfiles[0];
$meas_file1 = $measfiles[0];
if (not -f $meas_file1) {die("\n check if measurment file $meas_file1 exists\n")}
if (not -f $win_file1) {die("\n check if window file $win_file1 exists\n")}

$sTmin = sprintf("T%3.3i",$Tmin);
$sTmax = sprintf("T%3.3i",$Tmax);
$Ttag = "${sTmin}_${sTmax}";
$Ttag_title = sprintf("bandpass %i s to %i s",$Tmin,$Tmax);
$Ttag_plot = sprintf("bp [%i s, %i s]",$Tmin,$Tmax);

print "\nCHECKING USER INPUT:\n";
print "Event ID: ${eid}\n";
print "Station $sta network $net\n";
#print "Models $smodel1 and $smodel2\n";
print "Period range $Ttag\n";
print "Plotting interval is $trange seconds: $tmin to $tmax\n";
print "Window file 1 is $win_file1\n";
#print "Window file 2 is $win_file2\n";
print "Measurement file 1 is $meas_file1\n";
#print "Measurement file 2 is $meas_file2\n";
print "CMTSOLUTION: @cmtvals\n";
#die("CHECKING USER INPUT");

#---------------------------
# ADDITIONAL USER INPUT

$iplot_win = 1;      # plot windows
$iplot_dTlabel = 1;  # plot CC measurement for dT
$iplot_dAlabel = 0;  # plot CC measurement for dlnA
$iplot_sigma = 0;
$dTfsize = 9;

$ilims = 0;          # plot limits next to axis label
$iletter = 0;        # plot a letter for a publication figure
$letter = "B";

$sdX = 3.0; $sdY = 1.2; $origin = "-X7.5 -Y0.5";
$labfsize = 10;          # font size for model tag (e.g., "Z-m16")
$tpen = "1.0p";
$fpen = "1.0p";
$pen = "0.5p";           # line thickness for seismograms
$winfo = "-W1p -G230/230/230"; # specs for gray windows

# time ticks for seismograms
if ($trange >= 100) {$ttick1 = 20; $ttick2 = 10;}
elsif ($trange < 100 && $trange > 60) {$ttick1 = 10; $ttick2 = 10;}
else {$ttick1 = 10; $ttick2 = 2;}

# science paper figure
$sdY = 1.0; $iplot_win = 1; $iplot_dTlabel = 1; $iplot_dAlabel = 0; $iplot_sigma = 0; $iletter = 1; $letter = "B"; $pen = "1p";

# GJI paper figure
$sdY = 1.2; $iplot_win = 0; $iplot_dTlabel = 0; $iplot_dAlabel = 0; $iplot_sigma = 0; $iletter = 1; $letter = "b"; $pen = "1p";

# simple version
#$sdY = 1.0; $iplot_win = 0; $iplot_dTlabel = 1; $iplot_dAlabel = 0; $iplot_sigma = 0; $iletter = 0; $pen = "1p";

#---------------------------

# # extract values from the array cmtvals
# @M[0..5] = @cmtvals[12..17];
# $year = $cmtvals[0]; $month = $cmtvals[1]; $day = $cmtvals[2];
# $elat = $cmtvals[9]; $elon = $cmtvals[10]; $edep = $cmtvals[11];
# $eid = $cmtvals[18];
# $cmtpsmeca = "temp.txt";
# open(MIJ,">$cmtpsmeca");
# #printf MIJ "$elon $elat $edep $M[0] $M[1] $M[2] $M[3] $M[4] $M[5] 1 X Y $eid";
# printf MIJ "$elon $elat $edep $M[0] $M[1] $M[2] $M[3] $M[4] $M[5] 1";
# close(MIJ);

# # compute moment using Dahlen and Tromp, Eq. (5.91)
# # compute moment magnitude using Kanamori (1977) -- see Carl's Latex notes cmt.pdf
# $M0 = 1/sqrt(2) * sqrt($M[0]**2 + $M[1]**2 + $M[2]**2 + 2*($M[3]**2 + $M[4]**2 + $M[5]**2 ) );
# $k  = -11.2 + (2/3)*log(5)/log(10);
# $Mw = (2/3)*log($M0)/log(10) + $k;
# $stM0 = sprintf("M0 = %.3e dyne-cm",$M0);
# $stMw = sprintf("Mw = %.2f",$Mw);
# print "\n           moment : $stM0 \n moment magnitude : $stMw\n";

$fontno = 1;  # font number to use for text

# pstext info for measurement labels
$color = "0/255/255";
$color = "255/255/255";
$textinfo1 = "-C2p -W${color}o,${color} -G0/0/0";
$textinfo2 = "-C4p -W255/255/255o,1.0p,0/0/0,solid -G0/0/0";
$textinfo3 = "-C4p -W255/255/255o,1.0p,255/255/255,solid -G0/0/0";
$textinfo4 = "-C2p -W255/255/255o,1.0p,255/255/255,solid -G0/0/0";

#---------------------------

# set GMT defaults
`gmtset BASEMAP_TYPE plain PAPER_MEDIA letter ANNOT_FONT_SIZE = 9p LABEL_FONT_SIZE = 9p PAGE_ORIENTATION = landscape PLOT_DEGREE_FORMAT D FRAME_PEN $fpen TICK_PEN $tpen`;

#$name = "${sta}_${net}_win_adj_${klab}";
$name = "${eid}_${Ttag}_${sta}_${net}_${smodel1}_seis_${nseis}_multi";
$ps_file = "$name.ps";
$jpg_file = "$name.jpg";

$saclst="/opt/seismo-util/bin/saclst";

#------------
# data, synthetics

# DATA
$d_tag = "$eid.$net.$sta";
$Tdat = `ls ${opt_d}/${d_tag}.*T.sac*`; chomp($Tdat);
$Rdat = `ls ${opt_d}/${d_tag}.*R.sac*`; chomp($Rdat);
$Zdat = `ls ${opt_d}/${d_tag}.*Z.sac*`; chomp($Zdat);
if(-f $Tdat) {$Tflag = 1} else {$Tflag = 0}
if(-f $Rdat) {$Rflag = 1} else {$Rflag = 0}
if(-f $Zdat) {$Zflag = 1} else {$Zflag = 0}
#print "\n Tflag : $Tflag, Rflag : $Rflag, Zflag : $Zflag \n";

$smodel = $smodel1;

# SYNTHETICS -- component label always has the form BH_,
# even if the data file is something else (HH_, LH_).
# We assume that ALL componenets exist.
$s_tag = "semd.sac.${smodel}.${Ttag}";
$Tlab = "$sta.$net.BHT";
$Rlab = "$sta.$net.BHR";
$Zlab = "$sta.$net.BHZ";
$Tsyn = `ls ${opt_s}/${Tlab}*${s_tag}`; chomp($Tsyn);
$Rsyn = `ls ${opt_s}/${Rlab}*${s_tag}`; chomp($Rsyn);
$Zsyn = `ls ${opt_s}/${Zlab}*${s_tag}`; chomp($Zsyn);

#------------

#-----------------------------------
# limits for seismogram axes
$Smin = 0;
$Smax = 2;
$Srange = $Smax - $Smin;

$xtextpos1 = $tmin - 0.00*$trange;     # xmax label
$xtextpos2 = $tmin + 0.02*$trange;      # first window measurement label
$xtextpos3 = $tmin + 0.50*$trange;      # title position
$xtextpos4 = $tmin + 0.98*$trange;      # final window measurement label

$ytextpos5 = 1.35*$Srange;	# titles
#$ytextpos4 = 1.03*$Srange;    # dT label
#$ytextpos3 = 0.90*$Srange;    # dlnA label
$ytextpos4 = 0.95*$Srange;	# dT label
$ytextpos3 = 0.80*$Srange;	# dlnA label
$ytextpos2 = 0.35*$Srange;	# Z, R, T label
$ytextpos1 = 0.10*$Srange;	# ymax label

$ytextpos4 = 0.85*$Srange;	# lower dT label

# dimension of seismograms
$proj = "-JX${sdX}/${sdY}";
$dYseis = $sdY + 0.00;

$dXseis = $sdX + 0.5;
$dXseis2 = 2*$dXseis;

$Bseis0 = "-Ba${ttick1}f${ttick2}:\"  \":/a1f1";
#$tick = "-Ba50f10:\"Time (s)\":/a1f1S";
#$tick2 = "-Ba50f10/a1f1S";

$red_pen = "-W${pen},250/0/0";
$black_pen = "-W${pen},0/0/0";
$recon_pen = "-W${pen},0/255/255";

#============================================================================

$kmin = 1;
$kmax = $nseis;

$shift1 = "-X-${dXseis}";
$shift2 = "-X${dXseis2} -Y${dYseis}";

# loop over sets of seismograms (
for ($kk = $kmin; $kk <= $kmax; $kk++) {

  my(@Twinb); my(@Twine);
  my(@Rwinb); my(@Rwine);
  my(@Zwinb); my(@Zwine);
  my($undef);

  $smodel = $smodels[$kk-1];
  $winfile = $winfiles[$kk-1];
  $measfile = $measfiles[$kk-1];
  $smodeltag = $smodel;
  if($smodel eq "m99") {$smodeltag = "1D";}

  if($kk == 1) {$Bseis = "${Bseis0}::weSn";}
  elsif($kk == $kmax) {$Bseis = "${Bseis0}::wesN";}
  else {$Bseis = "${Bseis0}::wesn";}


# SYNTHETICS -- component label always has the form BH_,
# even if the data file is something else (HH_, LH_).
# We assume that ALL componenets exist.
$s_tag = "semd.sac.${smodel}.${Ttag}";
$Tlab = "$sta.$net.BHT";
$Rlab = "$sta.$net.BHR";
$Zlab = "$sta.$net.BHZ";
$Tsyn = `ls ${opt_s}/${Tlab}*${s_tag}`; chomp($Tsyn);
$Rsyn = `ls ${opt_s}/${Rlab}*${s_tag}`; chomp($Rsyn);
$Zsyn = `ls ${opt_s}/${Zlab}*${s_tag}`; chomp($Zsyn);

# windows and measurements: transverse component
$ntline=`grep -n \"$Tsyn\" $winfile | awk -F: '{ print \$1 }'`;
chomp($ntline);
if ($ntline) {
  $nline=$ntline;
  $nline=$nline+1;
  $nT=`sed -n ${nline}p $winfile`;
  @mlines = `grep $Tlab $measfile`; # measurements
  for ($i=0; $i < $nT; $i++) {
    $nline=$nline+1;
    $win=`sed -n ${nline}p $winfile`; chomp($win);
    (undef,$winb,$wine)=split(/ +/,$win);
    push @Twinb, "$winb";
    push @Twine, "$wine";

    (undef,undef,undef,undef,undef,$kplotT,undef,undef,$chi2T,$chi3T,$chi4T,$chi5T,$meas2T,$meas3T,$meas4T,$meas5T,$sigma2T,$sigma3T,$sigma4T,$sigma5T) = split(" ",$mlines[$i]);
    if ($iplot_sigma==1) {
      $stlabT_dT[$i] = sprintf("\@~D\@~T = %.2f \\261 %.2f s",$meas4T,$sigma4T);
      $stlabT_dA[$i] = sprintf("\@~D\@~A = %.2f \\261 %.2f",$meas5T,$sigma5T);
    } else {
      $stlabT_dT[$i] = sprintf("\@~D\@~T = %.2f s",$meas4T);
      $stlabT_dA[$i] = sprintf("\@~D\@~A = %.2f",$meas5T);
    }
  }
}
#print "Transverse component: ntline = $ntline, nT = $nT\n";

# windows and measurements: radial component
  $nrline=`grep -n \"$Rsyn\" $winfile | awk -F: '{ print \$1 }'`;
  chomp($nrline);
  #print "$nrline\n";
  if ($nrline) {
    $nline=$nrline;
    $nline=$nline+1;
    $nR=`sed -n ${nline}p $winfile`;
    @mlines = `grep $Rlab $measfile`; # measurements
    for ($i=0; $i < $nR; $i++) {
      $nline=$nline+1;
      $win=`sed -n ${nline}p $winfile`;
      chomp($win);
      #  print "$win\n";
      (undef,$winb,$wine)=split(/ +/,$win);
      #  print "$winb\n";
      #  print "$wine\n";
      push @Rwinb, "$winb";
      push @Rwine, "$wine";

      (undef,undef,undef,undef,undef,$kplotR,undef,undef,$chi2R,$chi3R,$chi4R,$chi5R,$meas2R,$meas3R,$meas4R,$meas5R,$sigma2R,$sigma3R,$sigma4R,$sigma5R) = split(" ",$mlines[$i]);
      if ($iplot_sigma==1) {
	$stlabR_dT[$i] = sprintf("\@~D\@~T = %.2f \\261 %.2f s",$meas4R,$sigma4R);
	$stlabR_dA[$i] = sprintf("\@~D\@~A = %.2f \\261 %.2f",$meas5R,$sigma5R);
      } else {
	$stlabR_dT[$i] = sprintf("\@~D\@~T = %.2f s",$meas4R);
	$stlabR_dA[$i] = sprintf("\@~D\@~A = %.2f",$meas5R);
      }
    }
  }
#print "Radial component: nrline = $nrline, nR = $nR\n";

# windows and measurements: vertical component
  $nzline=`grep -n \"$Zsyn\" $winfile | awk -F: '{ print \$1 }'`;
  chomp($nzline);
  #print "$nzline\n";
  if ($nzline) {
    $nline=$nzline;
    $nline=$nline+1;
    $nZ=`sed -n ${nline}p $winfile`;
    @mlines = `grep $Zlab $measfile`; # measurements
    for ($i=0; $i < $nZ; $i++) {
      $nline=$nline+1;
      $win=`sed -n ${nline}p $winfile`;
      chomp($win);
      #  print "$win\n";
      (undef,$winb,$wine)=split(/ +/,$win);
      #  print "$winb\n";
      #  print "$wine\n";
      push @Zwinb, "$winb";
      push @Zwine, "$wine";

      (undef,undef,undef,undef,undef,$kplotZ,undef,undef,$chi2Z,$chi3Z,$chi4Z,$chi5Z,$meas2Z,$meas3Z,$meas4Z,$meas5Z,$sigma2Z,$sigma3Z,$sigma4Z,$sigma5Z) = split(" ",$mlines[$i]);
      if ($iplot_sigma==1) {
	$stlabZ_dT[$i] = sprintf("\@~D\@~T = %.2f \\261 %.2f s",$meas4Z,$sigma4Z);
	$stlabZ_dA[$i] = sprintf("\@~D\@~A = %.2f \\261 %.2f",$meas5Z,$sigma5Z);
      } else {
	$stlabZ_dT[$i] = sprintf("\@~D\@~T = %.2f s",$meas4Z);
	$stlabZ_dA[$i] = sprintf("\@~D\@~A = %.2f",$meas5Z);
      }
    }
  }
#print "Vertical component: nrline = $nrline, nZ = $nZ\n";

#==========================================

# pssac2 scaling
#      -M vertical scaling in sacfile_unit/MEASURE_UNIT = size<required>
#          size: each trace will normalized to size (in MEASURE_UNIT)
#              scale =  (1/size) * [data(max) - data(min)]

# limits for synthetic records
(undef,$minT2,$maxT2) = split(" ",`$saclst depmin depmax f $Tsyn`);
$maxT2 = max(abs($minT2),abs($maxT2));
(undef,$minR2,$maxR2) = split(" ",`$saclst depmin depmax f $Rsyn`);
$maxR2 = max(abs($minR2),abs($maxR2));
(undef,$minZ2,$maxZ2) = split(" ",`$saclst depmin depmax f $Zsyn`);
$maxZ2 = max(abs($minZ2),abs($maxZ2));

# use the same vertical scale for m00 and m16 Z, the same vertical scale for m00 and m16 R, etc
if ($kk == 1) {

  # limits for data records
  $maxT = 0; $maxR = 0; $maxZ = 0;
  if ($Tflag==1) {(undef,$minT,$maxT)=split(" ",`$saclst depmin depmax f $Tdat`); $maxT=max(abs($minT),abs($maxT));}
  if ($Rflag==1) {(undef,$minR,$maxR)=split(" ",`$saclst depmin depmax f $Rdat`); $maxR=max(abs($minR),abs($maxR));}
  if ($Zflag==1) {(undef,$minZ,$maxZ)=split(" ",`$saclst depmin depmax f $Zdat`); $maxZ=max(abs($minZ),abs($maxZ));}

  # overall max values
  $maxTDS = max($maxT,$maxT2);
  $maxRDS = max($maxR,$maxR2);
  $maxZDS = max($maxZ,$maxZ2);
  $maxD = max($maxT,max($maxR,$maxZ));
  $maxS = max($maxT2,max($maxR2,$maxZ2));
  $max = max($maxD,$maxS);
  #print "\n Overall max values: $maxD (data), $maxS (syn) $max (overall)\n";

  $bounds = "-R$tmin/$tmax/$Smin/$Smax";
  $scale = 1.2;		       # KEY: scaling for plotting seismograms
  $size = 1/$scale;

}

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

# strings for plotting -- show the yaxis limits for each trace (based on synthetics)
$strT = sprintf('%.2e',$maxT2/$sizeT2);
$strR = sprintf('%.2e',$maxR2/$sizeR2);
$strZ = sprintf('%.2e',$maxZ2/$sizeZ2);

#==========================================

# TRANSVERSE component: data, synthetics, and windows
# synthetics
if ($kk==1) {
  `psbasemap $origin $proj $bounds -K $Bseis > $ps_file`;   # START
} else {
  `psbasemap $shift2 $proj $bounds -K -O $Bseis >> $ps_file`;
}
`pssac2 $proj $Tsyn $bounds -Ent-3 -M${sizeT2} $red_pen -N -K -O $Bseis >> $ps_file`;
if ($ntline) {
  if ($iplot_win==1) {
    for ($i=0; $i<$nT; $i++) {
      open(GMT,"|psxy $proj $bounds $winfo -O -K >> $ps_file");
      print GMT "$Twinb[$i] 0\n $Twine[$i] 0\n $Twine[$i] 2\n $Twinb[$i] 2\n";
      #print GMT "$Twinb[$i] 0\n";
      #print GMT "$Twine[$i] 0\n";
      #print GMT "$Twine[$i] 2\n";
      #print GMT "$Twinb[$i] 2\n";
      close(GMT);
    }
  }
  if ($iplot_dAlabel || $iplot_dTlabel) {
    for ($i=0; $i<$nT; $i++) {
      $xtext = $Twinb[$i]; $just = "LB";
      if ($i == 0) {
	$xtext = $xtextpos2; $just = "LB";
      }
      if ($i == $nT-1) {
	$xtext = $xtextpos4; $just = "RB";
      }
      if($iplot_dTlabel) {`echo \"$xtext $ytextpos4 $dTfsize 0 1 $just $stlabT_dT[$i]\" | pstext $textinfo1 $proj $bounds -N -O -K >> $ps_file`;}
      if($iplot_dAlabel) {`echo \"$xtext $ytextpos3 $dTfsize 0 1 $just $stlabT_dA[$i]\" | pstext $textinfo1 $proj $bounds -N -O -K >> $ps_file`;}
    }
  }
}
if ($Tflag==1) {`pssac2 $proj $Tdat $bounds -Ent-3 -M${sizeT} $black_pen -N -K -O >> $ps_file`;}    # data
`pssac2 $proj $Tsyn $bounds -Ent-3 -M${sizeT2} $red_pen -N -K -O >> $ps_file`;                      # synthetics
#if ($ntline && ($iplot_recon==1)) {`pssac2 $proj $Trecon $bounds -Ent-3 -M${sizeT3} $recon_pen -N -K -O >> $ps_file`;}   # reconstructed
#`psbasemap $proj $bounds -K -O $Bseis >> $ps_file`;  # replot basemap
$Tlabel = "T-${smodeltag}";
if($ilims==1) {`echo \"$xtextpos1 $ytextpos2 9 0 1 CM $strT\" | pstext $proj $bounds $textinfo3 -N -O -K >> $ps_file`;}
`echo \"$xtextpos1 $ytextpos1 $labfsize 0 1 BC $Tlabel\" | pstext $proj $bounds $textinfo2 -N -O -K >> $ps_file`;

#-----------------------

# RADIAL component: data, synthetics, and windows
#`pssac2 -Y$dYseis $proj $Rsyn $bounds -Ent-3 -M${sizeR2} $red_pen -N -K -O $Bseis >> $ps_file`;
`psbasemap $shift1 $proj $bounds -K -O $Bseis >> $ps_file`;
`pssac2 $proj $Rsyn $bounds -Ent-3 -M${sizeR2} $red_pen -N -K -O $Bseis >> $ps_file`;
if ($nrline) {
  if ($iplot_win==1) {
    for ($i=0; $i<$nR; $i++) {
      open(GMT,"|psxy $proj $bounds $winfo -O -K >> $ps_file");
      print GMT "$Rwinb[$i] 0\n $Rwine[$i] 0\n $Rwine[$i] 2\n $Rwinb[$i] 2\n";
      close(GMT);
    }
  }
  if ($iplot_dTlabel==1 || $iplot_dAlabel==1) {
    for ($i=0; $i<$nR; $i++) {
      $xtext = $Rwinb[$i]; $just = "LB";
      if ($i == 0) {
	$xtext = $xtextpos2; $just = "LB";
      }
      if ($i == $nR-1) {
	$xtext = $xtextpos4; $just = "RB";
      }
      if($iplot_dTlabel==1) {`echo \"$xtext $ytextpos4 $dTfsize 0 1 $just $stlabR_dT[$i]\" | pstext $textinfo1 $proj $bounds -N -O -K >> $ps_file`;}
      if($iplot_dAlabel==1) {`echo \"$xtext $ytextpos3 $dTfsize 0 1 $just $stlabR_dA[$i]\" | pstext $textinfo1 $proj $bounds -N -O -K >> $ps_file`;}
    }
  }
}
if ($Rflag==1) {`pssac2 $proj $Rdat $bounds -Ent-3 -M${sizeR} $black_pen -N -K -O >> $ps_file`;}
`pssac2 $proj $Rsyn $bounds -Ent-3 -M${sizeR2} $red_pen -N -K -O >> $ps_file`;
#if ($nrline && ($iplot_recon==1)) {`pssac2 $proj $Rrecon $bounds -Ent-3 -M${sizeR3} $recon_pen -N -K -O >> $ps_file`;}
#`psbasemap $proj $bounds -K -O $Bseis >> $ps_file`;  # replot basemap
$Rlabel = "R-${smodeltag}";
if($ilims==1) {`echo \"$xtextpos1 $ytextpos2 9 0 1 CM $strR\" | pstext $proj $bounds $textinfo3 -N -O -K >> $ps_file`;}
`echo \"$xtextpos1 $ytextpos1 $labfsize 0 1 BC $Rlabel\" | pstext $proj $bounds $textinfo2 -N -O -K >> $ps_file`;

#-----------------------

# VERTICAL component: data, synthetics, and windows
#`pssac2 -Y$dYseis $proj $Zsyn $bounds -Ent-3 -M${sizeZ2} $red_pen -N -K -O $Bseis >> $ps_file`;
`psbasemap $shift1 $proj $bounds -K -O $Bseis >> $ps_file`;
`pssac2 $proj $Zsyn $bounds -Ent-3 -M${sizeZ2} $red_pen -N -K -O $Bseis >> $ps_file`;
if ($nzline) {
  if ($iplot_win==1) {
    for ($i=0; $i<$nZ; $i++) {
      open(GMT,"|psxy $proj $bounds $winfo -O -K >> $ps_file");
      print GMT "$Zwinb[$i] 0\n $Zwine[$i] 0\n $Zwine[$i] 2\n $Zwinb[$i] 2\n";
      close(GMT);
    }
  }
  if ($iplot_dTlabel==1 || $iplot_dAlabel==1) {
    for ($i=0; $i<$nZ; $i++) {
      $xtext = $Zwinb[$i]; $just = "LB";
      if ($i == 0) {
	$xtext = $xtextpos2; $just = "LB";
      }
      if ($i == $nZ-1) {
	$xtext = $xtextpos4; $just = "RB";
      }
      if($iplot_dTlabel==1) {`echo \"$xtext $ytextpos4 $dTfsize 0 1 $just $stlabZ_dT[$i]\" | pstext $textinfo1 $proj $bounds -N -O -K >> $ps_file`;}
      if($iplot_dAlabel==1) {`echo \"$xtext $ytextpos3 $dTfsize 0 1 $just $stlabZ_dA[$i]\" | pstext $textinfo1 $proj $bounds -N -O -K >> $ps_file`;}
    }
  }
}
if ($Zflag==1) {`pssac2 $proj $Zdat $bounds -Ent-3 -M${sizeZ}  $black_pen -N -K -O $Bseis >> $ps_file`;}
`pssac2 $proj $Zsyn $bounds -Ent-3 -M${sizeZ2} $red_pen -N -O -K >> $ps_file`;
#if ($nzline && ($iplot_recon==1)) {`pssac2 $proj $Zrecon $bounds -Ent-3 -M${sizeZ3} $recon_pen -N -K -O >> $ps_file`;}
#`psbasemap $proj $bounds -K -O $Bseis >> $ps_file`;  # replot basemap
$Zlabel = "Z-${smodeltag}";
if($ilims==1) {`echo \"$xtextpos1 $ytextpos2 9 0 1 CM $strZ\" | pstext $proj $bounds $textinfo3 -N -O -K >> $ps_file`;}
`echo \"$xtextpos1 $ytextpos1 $labfsize 0 1 BC $Zlabel\" | pstext $proj $bounds $textinfo2 -N -K -O >> $ps_file`;

# plot overall label (for publication)
if ($iletter==1 && $kk==$kmax) {
  $textinfo = "-N -C4p -W255/255/255o,1.0p,0/0/0,solid -G0/0/0";
  `echo \"-0.15 1.2 20 0 $fontno TL $letter\" | pstext -R0/1/0/1 -JM1 $textinfo -K -O -V >>$ps_file`;
}


# plot titles
#`echo \"$xtextpos3 $ytextpos5 12 0 $fontno BC Data (Black) and Synthetics (Red) -- ${Ttag_title}\" | pstext $proj $bounds -N -K -O >> $ps_file`;


}   # for loop

#----------------------------------------

# plot titles
#`echo \"$xtextpos3 $maxT 12 0 $fontno BC Adjoint Source ($ktitle)\" | pstext $proj $bounds -N -Y3.2 -K -O >> $ps_file`;
#`echo \"$xtextpos3 $maxT 12 0 $fontno BC Data and Synthetics\" | pstext $proj $bounds -N -O -X-4 >> $ps_file`;
`echo \"0 0\" | psxy $proj $bounds -N -O -X15 >> $ps_file`;   # FINISH

#`convert $ps_file $jpg_file ; xv $jpg_file &`;
`ps2pdf $ps_file`;
#`rm $ps_file`;       # remove PS file
system("gv ${ps_file} &");

print "done with ${name}.pdf\n";

#================================================

# subroutine name: max
# Input: number1, number2
# returns greater of 2 numbers
sub max {
	if ($_[0]<$_[1]) {return $_[1]} else {return $_[0]};
}

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
