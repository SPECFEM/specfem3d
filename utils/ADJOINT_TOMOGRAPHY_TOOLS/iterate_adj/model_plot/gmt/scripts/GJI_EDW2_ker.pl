#!/usr/bin/perl -w

#==========================================================
#
#  GJI_EDW2_ker.pl
#  Carl Tape
#  26-May-2009
#
#  This script reads in a horizontal cross-section data file with five columns:
#     lon, lat, model-1, model-2, ln(model-2 / model-1)
#  and plots the three cross-sections.
#
#    GJI_EDW2_ker.pl -40/5 1 0.10 1.5 1 0.75
#
#==========================================================

if (@ARGV < 6) {
  die("Usage: GJI_EDW2.pl xxx\n");
}
($zran,$ivel,$cpert2,$vexag,$irun,$heightxc0) = @ARGV;

$icolor = 1;
$ivertlab = 0;      # vertical exaggeration label
$isimbox = 1;     # box around simulation region
$iletter = 1;     # label for publication figures
$letter = "a";

$pmin = 11; $pmax = 11;
$smodel1 = "m00"; $smodel2 = "m16";

#($pmin,$pmax) = split(/\//,$pinds);
#($smodel1,$smodel2) = split(/\//,$smodels);
$smodeltag = "${smodel2}_${smodel1}";
if ($ivel==1) {
  $modlab = "vs"; $mtit = "Vs"; $vmin = 2.0; $vmax = 4.21; $icolv = 3; $vtick = 1;
} elsif ($ivel==2) {
  $modlab = "vb"; $mtit = "Vb"; $vmin = 2.0; $vmax = 5.5; $icolv = 6; $vtick = 1;
} elsif ($ivel==3) {
  $modlab = "vp"; $mtit = "Vp";  $vmin = 2.5; $vmax = 7.5; $icolv = 9; $vtick = 1;
}
($zmin0,$zmax0) = split(/\//,$zran);

# directory containing the data files
$stirun = sprintf("%2.2i",$irun);
$dirdat = "INPUT/vert_${stirun}";

#---------------------------------------------------------
# absolute origins

$widmap = 2.5;      # with of horizontal cross sections
$height = 2.1;      # approximate (depends on projection)
$dxgap = 0.5;
$dygap = 0.5;
$x1 = 1.5; $y1 = 7.5;   # Vs m16 horizontal cross section
$or1 = "-Xa${x1} -Ya${y1}";

$x2 = $x1 + $widmap + $dxgap; $y2 = $y1; $or2 = "-Xa${x2} -Ya${y2}";
$x3 = $x1; $y3 = $y2 - $ygap - $height; $or3 = "-Xa${x3} -Ya${y3}";
$x4 = $x2; $y4 = $y3; $or4 = "-Xa${x4} -Ya${y4}";
$x5 = $x1; $y5 = $y3 - $ygap - $height; $or5 = "-Xa${x5} -Ya${y5}";
$x6 = $x2; $y6 = $y5; $or6 = "-Xa${x6} -Ya${y6}";
@omaps = ($or1,$or2,$or3,$or4,$or5,$or6);

$Dthick = 0.1;      # thickness of colorbar
$xgap = 0.2;
$ygap = 0.5;
$xo1 = $x1 - $xgap - $Dthick; $yo1 = $y1 + $ygap; $ocb1 = "-Xa${xo1} -Ya${yo1}";
$xo2 = $x2 + $widmap + $xgap; $yo2 = $yo1; $ocb2 = "-Xa${xo2} -Ya${yo2}";
$xo3 = $xo1; $yo3 = $y3 + $ygap; $ocb3 = "-Xa${xo3} -Ya${yo3}";
$xo4 = $xo2; $yo4 = $yo3; $ocb4 = "-Xa${xo4} -Ya${yo4}";
$xo5 = $xo1; $yo5 = $y5 + $ygap; $ocb5 = "-Xa${xo5} -Ya${yo5}";
$xo6 = $xo2; $yo6 = $yo5; $ocb6 = "-Xa${xo6} -Ya${yo6}";
@ocbs = ($ocb1,$ocb2,$ocb3,$ocb4,$ocb5,$ocb6);

# vertical cross section
$xvert = $x1 + 0.0;
$yvert = $y1 + 1.9;
$overt = "-Xa${xvert} -Ya${yvert}";

#---------------------------------------------------------

# subtitles
@titles = ("$mtit  $smodel1","$mtit  $smodel2","ln(${smodel2} / ${smodel1})");

# parameters for color scales
#$cpert1 = 0.2;
#$cpert2 = 0.2;
$scale_color = 65.0;    # 65.0
$colorbar = "seis";

# KEY: resolution of color plots
$interp = "-I0.5/0.5 -S1.5";  # nearneighbor
#$interp = "-I1/1 -S2";   # nearneighbor
#$interp = "-I1.0";            # xyz2grd
$grdfile = "temp.grd";

@Bopts = ("WESN","Wesn","wesN","wEsn","weSn","WesN","wEsN","wESn","WeSn","WEsn","weSN","wESN","WESn","WeSN","WEsN","wesn");

$dir0 = "/home/carltape/gmt";
$dir = "socal_2005";

$fault_file   = "$dir0/faults/jennings_more.xy";
$fault_file2  = "$dir0/faults/scitex.lin";
#$kcf_file     = "$dir0/faults/kcf.xy";
#$breck_file   = "$dir0/faults/breck.xy";

$plate_file  = "$dir0/plates/plate_boundaries/bird_boundaries";
$plate_labels = "$dir0/socal_2005/socal_plate_labs.xyz";

# continental shelf
$ishelf = 0;
$shelf_file = "/home/carltape/gmt/coastlines/socal_shelf";
$Wshelf = "-W0.75p,0/0/0,--";

# CMT sources
$dir_sources = "/home/carltape/results/SOURCES/socal_16";
$outer_boundary = "${dir_sources}/SPECFEM_outer_points.dat";
$inner_boundary = "${dir_sources}/SPECFEM_inner_points.dat";

# cities
$cities = "$dir0/cities/socal_tomo_cities";

# socal map
#$tick1 = 100; $tick2 = 100;  # no ticks for map
$tick1 = 1; $tick2 = 0.5;
$tlen_map = "3p";
#$Rmap = "-R-122/-114/32/37";
$Rmap = "-R-119/-115/33/35.2";

$topo_labels  = "$dir0/socal_2005/socal_topo_labs.xyz";
$fault_labels = "$dir0/socal_2005/socal_fault_labs.xyz";
$coast_res = "-Df -A0/0/4";

$B = "$B0".$Bopts[0];

$iportrait = 1;
if ($iportrait == 1) {
  $orient = "-P"; $rotangle = 0;
} else {
  $orient = " "; $rotangle = 90;
}

# plotting specifications
$fsize0 = "14";
$fsize1 = "12";
$fsize2 = "9";
$fsize3 = "6";
$fontno = "1";      # 1 or 4
$tlen   = "6p";
$fpen   = "1p";
$tpen   = "1p";

# plotting specificiations
# NOTE: pen attribute (-W) does not work for psvelo
#$coast_info     = "-A5000 -Dl -W1.0p -S220/255/255 -C220/255/255 -G200";  # A : smallest feature plotted, in km^2; D : resolution
$coast_res      = "-Df -A0/0/4";
$coast_info     = "${coast_res} -W0.75p -Na/0.75p";
#$coast_info     = "${coast_res} -W1p -Na/1.0p -S220/255/255 -C220/255/255 -G200";  # -I1/1.0p for rivers
$coast_info2     = "${coast_res} -S220/255/255 -C220/255/255 -G200"; # -I1/1.0p for rivers
$fault_info     = "-M -W0.75p,0/0/0";

$sky_color = "100/100/100";
$city_info = "-Sc10p -W1.0p/0/0/0 -G255/255/0";

#$colinfo = "-W1p,0/0/0 -G255/255/255";
$colinfo_map = "-W1.0p,0/0/0 -G0/255/255";
$srcinfo_map = "-Sa15p $colinfo_map -N";
$recinfo_map = "-Si15p $colinfo_map -N";
$colinfo_xc = "-W1.0p,0/0/0 -G0/255/255";
$srcinfo_xc  = "-Sa15p $colinfo_xc -N";
$recinfo_xc  = "-Si15p $colinfo_map -N";

$textinfo       = "-G255 -S1.5p";
$textinfo2      = "-G0/0/0 -S2p,255/255/255";
$textinfo3      = "-G0/0/0 -S2p,255/255/0";

# BOOLEAN: WHICH FIGURES TO PLOT
#$ixv = 1; $ipdf = 0;
$ixv = 0; $ipdf = 1;
$ijpg = 0;

# plot title
$J_title = "-JM1";
$R_title = "-R0/1/0/1";
$fsize_title = $fsize0;

#---------------------------------------------------------

# load values and check that the files exist
($nump,undef,undef) = split(" ",`ls -1 ${dirdat}/*path | wc`);
if ($nump == 0) {
  die("No ray paths available");
}
if ($pmax > $nump) {
  $pmax = $nump;
}
if ($pmin < 1) {
  $pmin = 1;
}
print "\nNumber of cuts is $nump\n";

# load fault positions along each profile
#$ffile = "$dirdat/ALL_${stirun}_rays_fault_positions";
$ffile = "$dirdat/ALL_${stirun}_rays_fault_positions_mod";
if (not -f $ffile) {
  die("Check if ffile $ffile exist or not\n");
}
open(IN,$ffile); @faults = <IN>; close(IN);

$cshfile = "GJI_EDW2.csh";
open(CSH,">$cshfile");
print CSH "gmtset COLOR_NAN $sky_color PAPER_MEDIA letter MEASURE_UNIT inch BASEMAP_TYPE plain PLOT_DEGREE_FORMAT D TICK_LENGTH $tlen LABEL_FONT_SIZE $fsize2 ANOT_FONT_SIZE $fsize2  HEADER_FONT $fontno ANOT_FONT $fontno LABEL_FONT $fontno HEADER_FONT_SIZE $fsize1 FRAME_PEN $fpen TICK_PEN $tpen \n";

#@iBs = (9,3,14,6,14,6);   # what sides of each cross section to show tick marks
@iBs = (7,4,7,4,7,4);

# modify the height of each cross-section based on the length
$p = 11;
$stip = sprintf("%3.3i",$p);

# datafile -- COLUMNS: lon, lat, m00, m16
$dfile = "${dirdat}/vert_${stirun}_xc_vsvbvp_${smodeltag}_${stip}.dat";
if (not -f $dfile) {
  die("Check if dfile $dfile exist or not\n");
}
print "Data file is $dfile\n";

# source and receiver points -- NOT necessarily the endpoints of the ray path
# COLUMNS: slons(ii),slats(ii),0,rlons(ii),rlats(ii),deg2km(dq),azq,iflip
$rfile = "${dirdat}/vert_${stirun}_xc_${stip}_A_B";
if (not -f $rfile) {
  die("Check if rfile $rfile exist or not\n");
}
open(IN,$rfile); @line = <IN>; close(IN);
($slon,$slat,$sdist,$sdep,$rlon,$rlat,$rdist,$rdep,$raz,$iflip) = split(" ",$line[0]);
#print "-- $slon -- $slat -- $sdist -- $sdep --\n";
#print "$rlon -- $rlat -- $rdist -- $rdep -- $raz -- $iflip--\n";

# positions of six faults along each profile
# (1) SAF, (2) GF, (3) SGF, (4) MCF, (5) SYF, (6) CRF, (7) EF, (8) SN, (9) KC
(undef,$eid,$stanet,$azi1,$azi2,$saf,$gf,$sgf,$mcf,$syf,$crf,$ef,$snf,$kcf) = split(" ",$faults[$p-1]);
($sta,$net) = split("\\.",$stanet); # only for irun = 1
$sazi1 = sprintf("%3.3i",$azi1);
$sazi2 = sprintf("%3.3i",$azi2);
print "-- $stanet -- $sta -- $net -- $eid -- $sazi1 -- $sazi2 --\n";

# axes limits for cross sections
($dmin,$dmax,$zmin,$zmax) = split(" ",`minmax -C -I2 $dfile`);
$dran = $dmax - $dmin; $zran = $zmax - $zmin;
$heightxc = $heightxc0;

# compute the height, including vertical exaggeration
$wid = $heightxc*$dran/$zran / $vexag;
$dfile = $dfile;
$R = "-R$dmin/$dmax/$zmin0/$zmax0";
$J = "-JX$wid/$heightxc";
if ($iflip==1) {
  $J = "-JX-$wid/$heightxc";
}

#$widmap = $heightxc*(1.9/1.5); # width of map -- tied to height of xc
$Jmap = "-JM${widmap}";
$B0map = "-Ba${tick1}f${tick2}d:.\" \":";
#$Bmap = "$B0map".$Bopts[7];
$Bmap = "$B0map".$Bopts[15];

#---------------------------------------------------------
# USER INPUT

# parameters controlling the placement of subplots and labels
$x0_map = 1;
$y0_map = 8;
$x0_xc = $x0_map + 3.5; $y0_xc = $y0_map;
$ygap_lab = 0.1; $xgap_lab = 0.6*$ygap_lab; # labels on the xc
$xgap_cbar = 0.2;   # gap between plot and cbar

# factors (times $heightxc) controlling y position of colorbar and labels
$fygap_bot = 0.0;
$f_cbar = 0.9;
$fygap_top = 0.0;
$f_lab = $fygap_bot + $f_cbar + $fygap_top;
$Dlen = $heightxc*$f_cbar;  # length of colorbar
$Dscale0 = "-D0/0/${Dlen}/${Dthick}";

#---------------------------------------------------------

# scale for plots
$tick1x = 50; $tick2x = 10;
$tick1z = 20; $tick2z = 10;
@Bopts = ("WESN","Wesn","wesN","wEsn","weSn","WesN","wEsN","wESn","WeSn","WEsn","weSN","wESN","WESn","WeSN","WEsN","wesn");
$B0 = sprintf("-Ba%3.3ff%3.3fd:\"   \":/a%3.3ff%3.3fd:\"   \"::.\"  \":WesN",$tick1x,$tick2x,$tick1z,$tick2z);
#$B0 = sprintf("-Ba%3.3ff%3.3fg%3.3f:\"   \":/a%3.3ff%3.3fg%3.3f:\"   \"::.\"  \":WesN",$tick1x,$tick2x,$tick1x,$tick1z,$tick2z,$tick1z);   # with gridlines

#==================================================

#-------------------
# absolute origins -- these depend on the WIDTH of the cross section

# overall label
$xgap = 0.1; $ygap = -0.1;
$x0_label = $xvert + $xgap;
$y0_label = $yvert + $heightxc0 + $ygap;
$olab = "-Xa${x0_label} -Ya${y0_label}";

# colorbar
$fac = $heightxc*$fygap_bot + $Dlen/2;
$xcb1 = $xvert + $wid + $xgap_cbar; $ycb1 = $yvert + $fac;
$ocb = "-Xa${xcb1} -Ya${ycb1}";

# label for colorbar
$xcblab1 = $xcb1 + 0.6; $ycblab1 = $ycb1 + 0.5;
$ocblab1 = "-Xa${xcblab1} -Ya${ycblab1}";
$xcblab2 = $xcblab1; $ycblab2 = $ycblab1 - 0.2;
$ocblab2 = "-Xa${xcblab2} -Ya${ycblab2}";

#-------------------

#$dtitle = sprintf("Depth  =  %.1f km",$zdep/1000);

$cpert1 = $vpert;

# colorpoint file for m00 and m16 -- VARIES with depth -- in KILOMETERS
$cptfile1 = "color1.cpt";
$dc = ($vmax-$vmin)/${scale_color};
$T = sprintf("-T%3.3e/%3.3e/%3.3e",$vmin,$vmax,$dc);
print CSH "makecpt -C$colorbar $T -D -Z > color.cpt\n"; # -Z for continuous (not for pscontour)
print CSH "sed 's/^F.*/F       200     200     200/' color.cpt > $cptfile1\n";

# colorbars for cross sections
$Bscale1 = sprintf("-B%2.2ef0.5:\"  \": -E8p",$vtick);

#==================================================

$stip = sprintf("%3.3i",$p);
$fname = "GJI_EDW2_ker";
$psfile = "$fname.ps"; $jpgfile = "$fname.jpg";

#---------------------
# 1: VERTICAL CROSS SECTION OF MODEL

# datafile -- COLUMNS: lon, lat, m00, m16
$dfile = "${dirdat}/vert_${stirun}_xc_vsvbvp_${smodeltag}_${stip}.dat";
if (not -f $dfile) {
  die("Check if dfile $dfile exist or not\n");
}
print "Data file is $dfile\n";

# source and receiver points -- NOT necessarily the endpoints of the ray path
# COLUMNS: slons(ii),slats(ii),0,rlons(ii),rlats(ii),deg2km(dq),azq,iflip
$rfile = "${dirdat}/vert_${stirun}_xc_${stip}_A_B";
if (not -f $rfile) {
  die("Check if rfile $rfile exist or not\n");
}
open(IN,$rfile); @line = <IN>; close(IN);
($slon,$slat,$sdist,$sdep,$rlon,$rlat,$rdist,$rdep,$raz,$iflip) = split(" ",$line[0]);

# axes limits for cross sections
($dmin,$dmax,$zmin,$zmax) = split(" ",`minmax -C -I2 $dfile`);
$dran = $dmax - $dmin; $zran = $zmax - $zmin;

# compute the height, including vertical exaggeration
#$heightxc = $zran*$wid/$dran * $vexag;
$wid = $heightxc*$dran/$zran / $vexag;
$R = "-R$dmin/$dmax/$zmin0/$zmax0";
$J = "-JX$wid/$heightxc";
if ($iflip==1) {
  $J = "-JX-$wid/$heightxc";
}

#-----------------------------------------------

$k = 2;

#$shift = $shifts[$k-1];
$B = "$B0".$Bopts[5];
$title = $titles[$k-1];

#------------------
# get absolute origins for various pieces

$Dscale = $Dscale0;

#------------------

print CSH "psbasemap $J $R $B -K -V $overt $orient > $psfile\n"; # START
print CSH "psbasemap $J $R $B -G${sky_color} -K -O -V $overt >> $psfile\n";

if ($icolor==1) {
  print CSH "awk '{print \$1,\$2,\$($icolv+$k-1)}' $dfile | nearneighbor -G$grdfile $R $interp\n";
  #print CSH "awk '{print \$1,\$2,\$4}' $dfile | xyz2grd -G$grdfile $R $interp\n";
  print CSH "grdimage $grdfile -C$cptfile1 $J -Q -K -O -V $overt >> $psfile\n";
}

print CSH "psscale -C${cptfile1} $Dscale $Bscale1 -K -O -V $ocb >> $psfile\n";
$slab1 = "Vs  m16";
print CSH "pstext -N $R_title $J_title -K -O -V $ocblab1 >>$psfile<<EOF\n0 0 $fsize2 0 $fontno CM $slab1\nEOF\n";
$slab2 = "km/s";
print CSH "pstext -N $R_title $J_title -K -O -V $ocblab2 >>$psfile<<EOF\n0 0 $fsize2 0 $fontno CM $slab2\nEOF\n";

# horizontal line
$finfo = "-W1.0p,0/0/0,--";
$zdep = -4; $x1 = -500; $x2 = 500;
print CSH "psxy $J $R $finfo -K -O -V $overt >>$psfile<<EOF\n$x1 $zdep\n$x2 $zdep\nEOF\n";

# vertical lines
$zminf0 = -5;
$finfo = "-W1.5p,0/0/0";
print CSH "psxy $J $R $finfo -K -O -V $overt >>$psfile<<EOF\n$saf $zmax0\n$saf $zminf0\nEOF\n";
print CSH "psxy $J $R $finfo -K -O -V $overt >>$psfile<<EOF\n$gf $zmax0\n$gf $zminf0\nEOF\n";
print CSH "psxy $J $R $finfo -K -O -V $overt >>$psfile<<EOF\n$sgf $zmax0\n$sgf $zminf0\nEOF\n";
#print CSH "psxy $J $R $finfo -K -O -V $overt >>$psfile<<EOF\n$mcf $zmax0\n$mcf $zminf0\nEOF\n";
print CSH "psxy $J $R $finfo -K -O -V $overt >>$psfile<<EOF\n$syf $zmax0\n$syf $zminf0\nEOF\n";
print CSH "psxy $J $R $finfo -K -O -V $overt >>$psfile<<EOF\n$crf $zmax0\n$crf $zminf0\nEOF\n";
print CSH "psxy $J $R $finfo -K -O -V $overt >>$psfile<<EOF\n$ef $zmax0\n$ef $zminf0\nEOF\n";
print CSH "psxy $J $R $finfo -K -O -V $overt >>$psfile<<EOF\n$snf $zmax0\n$snf $zminf0\nEOF\n";
print CSH "psxy $J $R $finfo -K -O -V $overt >>$psfile<<EOF\n$kcf $zmax0\n$kcf $zminf0\nEOF\n";

# plot source and receiver
if ($irun == 1) {
  print CSH "awk '{print \$3,-\$4}' $rfile | psxy $J $R $srcinfo_xc -K -O -V $overt >> $psfile\n";
  print CSH "awk '{print \$7,$zmax0}' $rfile | psxy $J $R $recinfo_xc -K -O -V $overt >> $psfile\n";
}

# details for specific cross-sections
$talign = "CB";
$flsize = 9;
$textinfo = "-N -C3p -W255/255/255o,1.0p,0/0/0,solid -G0/0/0";
if ($saf != 9999) {
  print CSH "pstext $J $R $textinfo -K -O -V $overt >>$psfile<<EOF\n$saf $zmax0 $flsize 0 $fontno $talign SA\nEOF\n";
}
if ($gf != 9999) {
  print CSH "pstext $J $R $textinfo -K -O -V $overt >>$psfile<<EOF\n$gf $zmax0 $flsize 0 $fontno $talign G\nEOF\n";
}
if ($sgf != 9999) {
  print CSH "pstext $J $R $textinfo -K -O -V $overt >>$psfile<<EOF\n$sgf $zmax0 $flsize 0 $fontno $talign SG\nEOF\n";
}
#if($mcf != 9999) {print CSH "pstext $J $R $textinfo -K -O -V $overt >>$psfile<<EOF\n$mcf $zmax0 $flsize 0 $fontno $talign MC\nEOF\n";}
if ($syf != 9999) {
  print CSH "pstext $J $R $textinfo -K -O -V $overt >>$psfile<<EOF\n$syf $zmax0 $flsize 0 $fontno $talign SY\nEOF\n";
}
if ($crf != 9999) {
  print CSH "pstext $J $R $textinfo -K -O -V $overt >>$psfile<<EOF\n$crf $zmax0 $flsize 0 $fontno $talign CR\nEOF\n";
}
if ($ef != 9999) {
  print CSH "pstext $J $R $textinfo -K -O -V $overt >>$psfile<<EOF\n$ef $zmax0 $flsize 0 $fontno $talign E\nEOF\n";
}
if ($snf != 9999) {
  print CSH "pstext $J $R $textinfo -K -O -V $overt >>$psfile<<EOF\n$snf $zmax0 $flsize 0 $fontno $talign SN\nEOF\n";
}
if ($kcf != 9999) {
  print CSH "pstext $J $R $textinfo -K -O -V $overt >>$psfile<<EOF\n$kcf $zmax0 $flsize 0 $fontno $talign KC\nEOF\n";
}

# plot overall label (for publication)
if ($iletter==1) {
  $textinfo = "-N -C3p -W255/255/255o,1.0p,0/0/0,solid -G0/0/0";
  print CSH "pstext $R_title $J_title $textinfo -K -O -V $olab >>$psfile<<EOF\n0 0 14 0 $fontno TL $letter\nEOF\n";
}

#    # plot label on each cross section
#    for ($k = 2; $k <= 2; $k ++ ) {
#      $title = $titles[$k-1];
#      $talign = "LB";
#      #$otitle = "-Xa${x3} -Ya${y3}";
#      $otitle = $overtlabs[$k-1];
#      $textinfo = "-N -C4p -W255/255/255o,1.0p,0/0/0,solid -G0/0/0";
#      print CSH "pstext $R_title $J_title $textinfo -K -O -V $otitle >>$psfile<<EOF\n0 0 12 0 $fontno $talign $title\nEOF\n";
#    }

#---------------------
# 2: HORIZONTAL CROSS SECTIONS OF MODEL

$icut = 13;     # depth to plot
$sticut = sprintf("%3.3i",$icut);
$irunh = 2; $stirunh = sprintf("%2.2i",$irunh);
$dirdath = "INPUT/horz_${stirunh}";

# KEY: file showing the cuts and the plotting range for the values
$fcuts = "${dirdath}/horz_${stirunh}_xc_cuts_mod";
if (not -f $fcuts) {
  die("Check if fcuts $fcuts exist or not\n");
}
open(IN,$fcuts); @lines = <IN>; $nump = @lines;

# factors (times $heightxc) controlling y position of colorbar and labels
$Dthick = 0.1;      # thickness of colorbar
$Dlen = 1.0;  # length of colorbar
$Dscale = "-D0/0/${Dlen}/${Dthick}";

#---------------------
# 3: HORIZONTAL CROSS SECTIONS OF KERNELS

# parameters controlling which kernel slices to plot
@icols = (4,4, 3,3,4,5);  # column in the datafile
@icuts = (1,13, 1,13,13,13);    # index into depth slice
@letters = ("b","c","d","e","f","g");

# KEY: data files
@dtags = ("${modlab}_${smodeltag}","${modlab}_${smodeltag}","ker_h000km_v000km","ker_h000km_v000km","ker_h003km_v001km","ker_h003km_v001km");

# color scales
@norms = (1e-3,1e-3, 1e12,1e12,1e12,1e12);
@kvals = (99,99, 1.21,0.81,1.01,0.81); # max value for color scale
@labels = ("Vs m16, depth 0 km","Vs m16, depth 4 km","K1, depth 0 km","K1, depth 4 km","K2, depth 4 km","K3, depth 4 km");
@clabels1 = (" "," ","10\@+-12\@+","10\@+-12\@+","10\@+-12\@+","10\@+-12\@+",);
$uker = "m\@+-3\@+ s";      # BDK units
@clabels2 = ("km/s","km/s",$uker,$uker,$uker,$uker);
@tick1s = (0.4,0.2,0.8,0.8,0.8,0.8);
@tick2s = (0.2,0.2,0.2,0.2,0.2,0.2);
@Aflip = ("-A"," ","-A"," ","-A"," ",);
@xtxs = (-0.1,0.2,-0.1,0.2,-0.1,0.2);

#for ($ik = 1; $ik <= 5; $ik ++ ) {
for ($ik = 1; $ik <= 6; $ik ++ ) {

  $origin = $omaps[$ik-1];
  $iB = $iBs[$ik-1];
  $Bmap = "$B0map".$Bopts[$iB];
  $Bscale1 = sprintf("-B%2.2ef%2.2e:\"  \": -E8p %s",$tick1s[$ik-1],$tick2s[$ik-1],$Aflip[$ik-1]);

  print CSH "psbasemap $Jmap $Rmap $Bmap -K -O -V $origin >> $psfile\n";

  $dtag = $dtags[$ik-1];
  $norm = $norms[$ik-1];
  $icol = $icols[$ik-1];
  $cmax = $kvals[$ik-1];
  $icut = $icuts[$ik-1]; $sticut = sprintf("%3.3i",$icut);

  # minmax for color scale
  if ($ik <= 2) {
    # load values and check that the files exist
    (undef,$zcut,$vs_norm,$vb_norm,$vpert) = split(" ",$lines[$icut-1]);
    $cnorm = $vs_norm;
    $cpert = $vpert;
    if ($ik==2) {
      $cpert = 0.10;
    }
    $cmin = $cnorm*exp(-$cpert)*$norm;
    $cmax = $cnorm*exp($cpert)*$norm;

  } else {
    $cmin = -$cmax;
  }

  # data file
  if($ik <= 2) {
     $dfileh = "${dirdath}/horz_${stirunh}_xc_${dtag}_${sticut}_mask1.dat";
  } else {
     $dfileh = "${dirdath}/horz_${stirunh}_xc_${dtag}_${sticut}.dat";
  }
  if (not -f ${dfileh}) {
    die("Check if dfileh ${dfileh} exist or not\n");
  }

  # colorpoint file
  $cptfileh = "colorK.cpt";
  $dc = ($cmax-$cmin)/${scale_color};
  $T = sprintf("-T%3.3e/%3.3e/%3.3e",$cmin,$cmax,$dc);
  print CSH "makecpt -C$colorbar $T -D -Z > $cptfileh\n";

  if ($icolor==1) {
    $interph = "-I1m -S2m";
    $grdfileh = "temph.grd";
    #print CSH "awk '{print \$1,\$2,\$${icol}*1e12}' $dfileh | xyz2grd -G$grdfile $Rmap $interph\n";
    print CSH "awk '{print \$1,\$2,\$${icol}*$norm}' $dfileh | nearneighbor -G$grdfileh $Rmap $interph\n";
    print CSH "grdimage $grdfileh -C$cptfileh $Jmap -Q -K -O -V $origin >> $psfile\n";
  }

  print CSH "pscoast $Jmap $Rmap $coast_info -K -O -V $origin >> $psfile \n";
  if ($ishelf==1) {
    print CSH "awk '{print \$2,\$1}' ${shelf_file} | psxy $Jmap $Rmap $Wshelf -K -O -V $origin >> $psfile\n";
  }


# colorbar
$ocb = $ocbs[$ik-1];

# color scale
print CSH "psscale -C${cptfileh} $Dscale $Bscale1 -K -O -V $ocb >> $psfile\n";
$slab1 = $clabels1[$ik-1];
print CSH "pstext -N $R_title $J_title -K -O -V $ocb >>$psfile<<EOF\n$xtxs[$ik-1] 1.0 $fsize2 0 $fontno CM $slab1\nEOF\n";
$slab2 = $clabels2[$ik-1];
print CSH "pstext -N $R_title $J_title -K -O -V $ocb >>$psfile<<EOF\n$xtxs[$ik-1] 0.8 $fsize2 0 $fontno CM $slab2\nEOF\n";

  print CSH "psxy ${fault_file} $Jmap $Rmap $fault_info -K -V -O $origin >> $psfile \n";
  #print CSH "psxy ${kcf_file} $Jmap $Rmap $fault_info -K -V -O $origin >> $psfile \n";

  # plot ray path
  $pfile = "${dirdat}/vert_${stirun}_xc_${stip}_ray_path";
  if (not -f $pfile) {
    die("Check if pfile $pfile exist or not\n");
  }
  print CSH "awk '{print \$1,\$2}' ${pfile} | psxy $Jmap $Rmap -W1.0p,0/0/0,-- -K -O -V $origin >> $psfile\n";
  #print CSH "awk '{print \$1,\$2}' ${pfile} | psxy $Jmap $Rmap -W6p,0/0/0 -K -O -V $origin >> $psfile\n";
  #print CSH "awk '{print \$1,\$2}' ${pfile} | psxy $Jmap $Rmap -W3p,0/255/255 -K -O -V $origin >> $psfile\n";

  # plot source and receiver lon-lat
  if ($irun==1) {
    $rfile = "${dirdat}/vert_${stirun}_xc_${stip}_A_B";
    if (not -f $rfile) {
      die("Check if rfile $rfile exist or not\n");
    }
    print CSH "awk '{print \$1,\$2}' ${rfile} | psxy $Jmap $Rmap $srcinfo_map -K -O -V $origin >> $psfile\n";
    print CSH "awk '{print \$5,\$6}' ${rfile} | psxy $Jmap $Rmap $recinfo_map -K -O -V $origin >> $psfile\n";
  }
  print CSH "psbasemap $Jmap $Rmap $Bmap -K -O -V $origin >> $psfile\n";

  # plot overall label (for publication)
  if ($iletter==1) {
    $letter = @letters[$ik-1];
    #$textinfo = "-N -C3p -W255/255/255o,1.0p,0/0/0,solid -G0/0/0";
    #print CSH "pstext $R_title $J_title $textinfo -K -O -V $olab >>$psfile<<EOF\n0 0 14 0 $fontno TL $letter\nEOF\n";
    print CSH "pstext $R_title $J_title $textinfo -K -O -V $origin >>$psfile<<EOF\n0.1 1.6 14 0 $fontno TL $letter\nEOF\n";
    $label = @labels[$ik-1];
    print CSH "pstext $R_title $J_title $textinfo -K -O -V $origin >>$psfile<<EOF\n2.4 1.6 10 0 $fontno TR $label\nEOF\n";
  }

}

#------------------------------------

print CSH "pstext -N $R_title $J_title -O -V >>$psfile<<EOF\n 10 10 16 0 $fontno CM \nEOF\n"; # FINISH
#if($ixv==1) {print CSH "convert $psfile -rotate $rotangle $jpgfile\n";}
if ($ixv==1) {print CSH "gv $psfile &\n";}
if ($ipdf==1) {print CSH "ps2pdf $psfile\n";}
if ($ijpg==1) {print CSH "convert $psfile $jpgfile\n";}

print CSH "gv $psfile &\n";

#==================================================
close (CSH);
system("csh -f $cshfile");
#==================================================
