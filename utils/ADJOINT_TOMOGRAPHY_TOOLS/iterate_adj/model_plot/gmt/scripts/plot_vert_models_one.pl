#!/usr/bin/perl -w

#==========================================================
#
#  plot_vert_models_one.pl
#  Carl Tape
#  07-Dec-2008
#
#  This script reads in a horizontal cross-section data file with five columns:
#     lon, lat, model-1, model-2, ln(model-2 / model-1)
#  and plots the three cross-sections.
#
#  EXAMPLE:
#    plot_vert_models_one.pl 1/86 m00/m16 -40/5 1 0.10 3.0 1 2.2 # Vs -- src-rec paths
#    plot_vert_models_one.pl 4/32 m00/m16 -40/5 1 0.10 3.0 2 2.2 # Vs -- SAF-perp profiles and three more (cut 14,25)
#    plot_vert_models_one.pl 1/10 m00/m16 -40/5 1 0.10 3.0 3 2.2 # Vs -- Garlock-perp profiles
#    plot_vert_models_one.pl 1/13 m00/m16 -40/5 1 0.10 3.0 4 2.2 # Vs -- N-S (constant longitude)
#    plot_vert_models_one.pl 1/9  m00/m16 -40/5 1 0.10 3.0 5 2.2 # Vs -- E-W (constant latitude)
#
#
#==========================================================

if (@ARGV < 8) {die("Usage: plot_vert_models_one.pl xxx\n");}
($pinds,$smodels,$zran,$ivel,$cpert2,$vexag,$irun,$heightxc0) = @ARGV;

$icolor = 1;
$ivertlab = 0;    # vertical exaggeration label
$isimbox = 1;     # box around simulation region
$iletter = 0;     # label for publication figures
$letter = "a";
$imfinal = 1;     # plot an additional figure with only m16 and the map

($pmin,$pmax) = split(/\//,$pinds);
($smodel1,$smodel2) = split(/\//,$smodels);
$smodeltag = "${smodel2}_${smodel1}";
if($ivel==1) {$modlab = "vs"; $mtit = "Vs"; $vmin = 2.0; $vmax = 4.21; $icol = 3; $vtick = 1;}
elsif($ivel==2) {$modlab = "vb"; $mtit = "Vb"; $vmin = 2.0; $vmax = 5.5; $icol = 6; $vtick = 1;}
elsif($ivel==3) {$modlab = "vp"; $mtit = "Vp";  $vmin = 2.5; $vmax = 7.5; $icol = 9; $vtick = 1;}
($zmin0,$zmax0) = split(/\//,$zran);

# directory containing the data files
$stirun = sprintf("%2.2i",$irun);
$dirdat = "INPUT/vert_${stirun}";

#---------------------------------------------------------

# subtitles
@titles = ("$mtit  $smodel1","$mtit  $smodel2","ln(${smodel2} / ${smodel1})");

# parameters for color scales
#$cpert1 = 0.2;
#$cpert2 = 0.2;
$scale_color = 65.0;   # 65.0
$colorbar = "seis";

# KEY: resolution of color plots
$interp = "-I0.5/0.5 -S1.5";   # nearneighbor
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
$tick1 = 1; $tick2 = 1;
$tlen_map = "3p";
#$Rmap = "-R-121.2/-114.8/32.3/36.7";
$Rmap = "-R-122/-114/32/37";

$topo_labels  = "$dir0/socal_2005/socal_topo_labs.xyz";
$fault_labels = "$dir0/socal_2005/socal_fault_labs.xyz";
$coast_res = "-Df -A0/0/4";

$B = "$B0".$Bopts[0];

$iportrait = 0;
if($iportrait == 1){$orient = "-P"; $rotangle = 0}
else {$orient = " "; $rotangle = 90}

# plotting specifications
$fsize0 = "18";
$fsize1 = "16";
$fsize2 = "14";
$fsize3 = "12";
$fontno = "1";    # 1 or 4
$tlen   = "6p";
$fpen   = "1p";
$tpen   = "1p";

# plotting specificiations
# NOTE: pen attribute (-W) does not work for psvelo
#$coast_info     = "-A5000 -Dl -W1.0p -S220/255/255 -C220/255/255 -G200";  # A : smallest feature plotted, in km^2; D : resolution
$coast_res      = "-Df -A0/0/4";
#$coast_info     = "${coast_res} -W0.75p -Na/0.75p";
$coast_info     = "${coast_res} -W1p -Na/1.0p -S220/255/255 -C220/255/255 -G200";  # -I1/1.0p for rivers
$coast_info2     = "${coast_res} -S220/255/255 -C220/255/255 -G200";  # -I1/1.0p for rivers
$fault_info     = "-M -W0.5p,255/0/0";
$fault_infoK    = "-M -W0.5p,0/0/0";

$sky_color = "100/100/100";
$city_info = "-Sc10p -W1.0p/0/0/0 -G255/255/0";

#$colinfo = "-W1p,0/0/0 -G255/255/255";
$colinfo_map = "-W1.5p,0/0/0 -G0/255/255";
$colinfo_xc = "-W2p,0/0/0";
$srcinfo_xc  = "-Sa30p $colinfo_xc -N";
$recinfo_xc  = "-Si30p $colinfo_map -N";

$textinfo       = "-G255 -S1.5p";
$textinfo2      = "-G0/0/0 -S2p,255/255/255";
$textinfo3      = "-G0/0/0 -S2p,255/255/0";

# BOOLEAN: WHICH FIGURES TO PLOT
$ixv = 0;
$ipdf = 0;
$ieps = 1;
$ijpg = 1;

# plot title
$J_title = "-JM1";
$R_title = "-R0/1/0/1";
$fsize_title = $fsize0;

#---------------------------------------------------------

# load values and check that the files exist
($nump,undef,undef) = split(" ",`ls -1 ${dirdat}/*path | wc`);
if($nump == 0) {die("No ray paths available")}
if($pmax > $nump) {$pmax = $nump;}
if($pmin < 1) {$pmin = 1;}
print "\nNumber of cuts is $nump\n";

# load fault positions along each profile
  $ffile = "$dirdat/ALL_${stirun}_rays_fault_positions_mod";
  if (not -f $ffile) {die("Check if ffile $ffile exist or not\n");}
  open(IN,$ffile); @faults = <IN>; close(IN);

$cshfile = "plot_vert_models_one.csh";
open(CSH,">$cshfile");
print CSH "gmtset COLOR_NAN $sky_color PAPER_MEDIA letter MEASURE_UNIT inch BASEMAP_TYPE plain PLOT_DEGREE_FORMAT D TICK_LENGTH $tlen LABEL_FONT_SIZE $fsize2 ANOT_FONT_SIZE $fsize2  HEADER_FONT $fontno ANOT_FONT $fontno LABEL_FONT $fontno HEADER_FONT_SIZE $fsize1 FRAME_PEN $fpen TICK_PEN $tpen \n";

@iBs = (5,8,4);  # what sides of each cross section to show tick marks

# modify the height of each cross-section based on the length
for ($p = $pmin; $p <= $pmax; $p ++ ) {
  $stip = sprintf("%3.3i",$p);

  # datafile -- COLUMNS: lon, lat, m00, m16
  $dfile = "${dirdat}/vert_${stirun}_xc_vsvbvp_${smodeltag}_${stip}.dat";
  if (not -f $dfile) {die("Check if dfile $dfile exist or not\n");}
  print "Data file is $dfile\n";

  # source and receiver points -- NOT necessarily the endpoints of the ray path
  # COLUMNS: slons(ii),slats(ii),0,rlons(ii),rlats(ii),deg2km(dq),azq,iflip
  $rfile = "${dirdat}/vert_${stirun}_xc_${stip}_A_B";
  if (not -f $rfile) {die("Check if rfile $rfile exist or not\n");}
  open(IN,$rfile); @line = <IN>; close(IN);
  ($slon,$slat,$sdist,$sdep,$rlon,$rlat,$rdist,$rdep,$raz,$iflip) = split(" ",$line[0]);
   #print "-- $slon -- $slat -- $sdist -- $sdep --\n";
   #print "$rlon -- $rlat -- $rdist -- $rdep -- $raz -- $iflip--\n";

   # positions of six faults along each profile
   # (1) SAF, (2) GF, (3) SGF, (4) MCF, (5) SYF, (6) CRF, (7) EF, (8) SN, (9) KC
   (undef,$eid,$stanet,$azi1,$azi2,$saf,$gf,$sgf,$mcf,$syf,$crf,$ef,$snf,$kcf) = split(" ",$faults[$p-1]);
   ($sta,$net) = split("\\.",$stanet);    # only for irun = 1
   $sazi1 = sprintf("%3.3i",$azi1);
   $sazi2 = sprintf("%3.3i",$azi2);
   print "-- $stanet -- $sta -- $net -- $eid -- $sazi1 -- $sazi2 --\n";

  # axes limits for cross sections
  ($dmin,$dmax,$zmin,$zmax) = split(" ",`minmax -C -I2 $dfile`);
  $dran = $dmax - $dmin; $zran = $zmax - $zmin;
  $heightxc = $heightxc0;
  #if ($dran > 450) {
  #  $heightxc = 1/($dran/700 + 0.03); # empirical relationship
  #}

  # compute the height, including vertical exaggeration
  $wid = $heightxc*$dran/$zran / $vexag;
  $dfile = $dfile;
  $R = "-R$dmin/$dmax/$zmin0/$zmax0";
  $J = "-JX$wid/$heightxc";
  if ($iflip==1) {$J = "-JX-$wid/$heightxc";}

  $widmap = $heightxc*(1.9/1.5); # width of map -- tied to height of xc
  $Jmap = "-JM${widmap}";
  $B0map = "-Ba${tick1}f${tick2}d:.\" \":";
  #$Bmap = "$B0map".$Bopts[7];
  $Bmap = "$B0map".$Bopts[15];

  #---------------------------------------------------------
  # USER INPUT

  # parameters controlling the placement of subplots and labels
  $x0_xc = 0.75; $y0_xc = 3;  # lower left of m00 cross section
  #$x0_map = 0.6;
  #$y0_map = 5;
  #$xgap_map_xc = 0.8;
  #$x0_xc = $x0_map + $widmap + $xgap_map_xc;
  #$y0_xc = $y0_map;
  $ygap_lab = 0.1; $xgap_lab = 0.6*$ygap_lab; # labels on the xc
  $xgap_cbar = 0.2;   # gap between plot and cbar
  #$omap = "-Xa${x0_map} -Ya${y0_map}";
  #$oxc  = "-Xa${x0_xc} -Ya${y0_xc}";

  # factors (times $heightxc) controlling y position of colorbar and labels
  $fygap_bot = 0.1;
  $f_cbar = 0.6;
  $fygap_top = 0.20;
  $f_lab = $fygap_bot + $f_cbar + $fygap_top;
  $Dthick = 0.20;   # thickness of colorbar
  $Dlen = $heightxc*$f_cbar;  # length of colorbar
  $Dscale0 = "-D0/0/${Dlen}/${Dthick}";
  #@Dscales = ($Dscale0,$Dscale0,"${Dscale0}h");
  @Dscales = ($Dscale0,$Dscale0,$Dscale0);

  # factors controlling position of the vertical exaggeration label
  $xv = $x0_xc;
  $yv = $y0_xc + $heightxc + 0.5;
  $ov = "-Xa${xv} -Ya${yv}";

  #---------------------------------------------------------

  # scale for plots
  $tick1x = 50; $tick2x = 10;
  $tick1z = 10; $tick2z = 5;
  if ($dran > 550) {
    $tick1x = 100; $tick1z = 20;
  }
  @Bopts = ("WESN","Wesn","wesN","wEsn","weSn","WesN","wEsN","wESn","WeSn","WEsn","weSN","wESN","WESn","WeSN","WEsN","wesn");
  $B0 = sprintf("-Ba%3.3ff%3.3fd:\"   \":/a%3.3ff%3.3fd:\"   \"::.\"  \":WesN",$tick1x,$tick2x,$tick1z,$tick2z);
  #$B0 = sprintf("-Ba%3.3ff%3.3fg%3.3f:\"   \":/a%3.3ff%3.3fg%3.3f:\"   \"::.\"  \":WesN",$tick1x,$tick2x,$tick1x,$tick1z,$tick2z,$tick1z);   # with gridlines


  #==================================================

  #-------------------
  # absolute origins -- these depend on the WIDTH of the cross section

  # cross sections
  $xgap_xc = 0.5;   # hspace between cross sections
  #if($wid < 3.5) {$xgap_xc = 1.0;}
  $x1 = $x0_xc; $y1 = $y0_xc; $oxc1 = "-Xa${x1} -Ya${y1}";
  $x2 = $x0_xc; $y2 = $y0_xc; $oxc2 = "-Xa${x2} -Ya${y2}";
  $x3 = $x0_xc; $y3 = $y0_xc; $oxc3 = "-Xa${x3} -Ya${y3}";
  @oxcs = ($oxc1,$oxc2,$oxc3);

  # cross-section label
  $x1_label = $x1 + $xgap_lab; $y1_label = $y1 + $ygap_lab; $oxclab1 = "-Xa${x1_label} -Ya${y1_label}";
  $x2_label = $x2 + $xgap_lab; $y2_label = $y2 + $ygap_lab; $oxclab2 = "-Xa${x2_label} -Ya${y2_label}";
  $x3_label = $x3 + $xgap_lab; $y3_label = $y3 + $ygap_lab; $oxclab3 = "-Xa${x3_label} -Ya${y3_label}";
  @oxclabs = ($oxclab1,$oxclab2,$oxclab3);

  # colorbar
  $fac = $heightxc*$fygap_bot + $Dlen/2;
  $xcb1 = $x1 + $wid + $xgap_cbar; $ycb1 = $y1 + $fac; $ocb1 = "-Xa${xcb1} -Ya${ycb1}";
  $xcb2 = $x2 + $wid + $xgap_cbar; $ycb2 = $y2 + $fac; $ocb2 = "-Xa${xcb2} -Ya${ycb2}";
  $xcb3 = $x3 - $xgap_cbar - $Dthick; $ycb3 = $y3 + $fac; $ocb3 = "-Xa${xcb3} -Ya${ycb3}";
  #$xcb3 = $x0_map + $widmap + $Dlen/2 + 0.4; $ycb3 = $y3 + $heightxc + 1.0; $ocb3 = "-Xa${xcb3} -Ya${ycb3}";
  #$xcb3 = $x3 + $wid - $Dlen/2; $ycb3 = $y3 + $heightxc + 1.0; $ocb3 = "-Xa${xcb3} -Ya${ycb3}";
  @ocbs = ($ocb1,$ocb2,$ocb3);

  # label for colorbar
  $xcblab1 = $xcb1 + $xgap_cbar; $ycblab1 = $y1 + $heightxc*$f_lab; $ocblab1 = "-Xa${xcblab1} -Ya${ycblab1}";
  $xcblab2 = $xcb2 + $xgap_cbar; $ycblab2 = $y2 + $heightxc*$f_lab; $ocblab2 = "-Xa${xcblab2} -Ya${ycblab2}";
  $ocblab3 = "-Xa0 -Ya0";
  @ocbarlabs = ($ocblab1,$ocblab2,$ocblab3);

  #-------------------

  #$dtitle = sprintf("Depth  =  %.1f km",$zdep/1000);

  $cpert1 = $vpert;

  # colorpoint file for m00 and m16 -- VARIES with depth -- in KILOMETERS
  $cptfile1 = "color1.cpt";
  $dc = ($vmax-$vmin)/${scale_color};
  $T = sprintf("-T%3.3e/%3.3e/%3.3e",$vmin,$vmax,$dc);
  print CSH "makecpt -C$colorbar $T -D -Z > color.cpt\n";   # -Z for continuous (not for pscontour)
  print CSH "sed 's/^F.*/F       200     200     200/' color.cpt > $cptfile1\n";

  # colorpoint file for perturbation ln(m16/m00) -- FIXED with depth
  $cptfile2 = "color2.cpt";
  $cmin = -$cpert2; $cmax = $cpert2;
  $dc = ($cmax-$cmin)/${scale_color};
  $T = sprintf("-T%3.3e/%3.3e/%3.3e",$cmin*1.01,$cmax*1.01,$dc);
  print CSH "makecpt -C$colorbar $T -D > $cptfile2\n";

  # colorbars for cross sections
  $Bscale1 = sprintf("-B%2.2ef0.5:\"  \": -E10p",$vtick);
  #$Bscale2  = sprintf("-Ba%2.2ef0.05:\"$titles[2]\": -E10p",$cpert2);
  $Bscale2  = sprintf("-Ba%2.2ef0.05:\" \": -E10p -A",$cpert2); # -A for ticks on other side

  #------------------------------------

  # plot ray path
  $pfile = "${dirdat}/vert_${stirun}_xc_${stip}_ray_path";
  if (not -f $pfile) {
    die("Check if pfile $pfile exist or not\n");
  }

  # plot source and receiver lon-lat
  $rfile = "${dirdat}/vert_${stirun}_xc_${stip}_A_B";
  if (not -f $rfile) {die("Check if rfile $rfile exist or not\n");}

  #-----------------------------------------------

  $k = 2;

  #$shift = $shifts[$k-1];
  $iB = $iBs[$k-1];
  $B = "$B0".$Bopts[$iB];
  $title = $titles[$k-1];

  #------------------

  # vertical lines
  $zminf0 = -5;
  $finfo = "-W1.5p,0/0/0";

  # details for specific cross-sections
  #if ($irun==1 && $k == 1) {
  (undef,undef,undef,undef,undef,$saf,$gf,$sgf,$mcf,$syf) = split(" ",$faults[$p-1]);
  $talign = "CB";
  $flsize = 9;
  $textinfo = "-N -C3p -W255/255/255o,1.0p,0/0/0,solid -G0/0/0";

#==================================================
# figure with map with m16 beneath

if ($imfinal == 1) {

  @iBs = (5,8,4);

  #for ($p = $pmin; $p <= $pmax; $p ++ ) {

    $stip = sprintf("%3.3i",$p);

    $fname = "vert_${stirun}_xc_${modlab}_${smodeltag}_${stip}_single";
    if($irun==1) {$fname = "vert_${stirun}_xc_${sta}_${net}_${sazi1}_${eid}_${modlab}_${smodeltag}_${stip}";}
    $psfile = "$fname.ps"; $jpgfile = "$fname.jpg"; $epsfile = "$fname.eps";

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

    #-------------------
    # absolute origins -- these depend on the WIDTH of the cross section

    # overall label
    $x0_label = $x0_xc;
    $y0_label = $y0_xc + $heightxc + 1;
    $olab = "-Xa${x0_label} -Ya${y0_label}";
    $x0_label2 = $x0_label;
    $y0_label2 = $y0_xc + 0.8*$heightxc + 1;
    $olab2 = "-Xa${x0_label2} -Ya${y0_label2}";

    #-----------------------------------------------

      #$shift = $shifts[$k-1];
      $iB = $iBs[$k-1];
      $B = "$B0".$Bopts[$iB];
      $title = $titles[$k-1];

      #------------------
      # get absolute origins for various pieces

      $oxc = $oxcs[$k-1];
      $ocb = $ocbs[$k-1];
      $ocbarlab = $ocbarlabs[$k-1];
      $Dscale = $Dscales[$k-1];

      #------------------

      print CSH "psbasemap $J $R $B -K -V $oxc > $psfile\n";
      print CSH "psbasemap $J $R $B -G${sky_color} -K -O -V $oxc >> $psfile\n";

      if ($icolor==1) {
  print CSH "awk '{print \$1,\$2,\$($icol+$k-1)}' $dfile | nearneighbor -G$grdfile $R $interp\n";
        #print CSH "awk '{print \$1,\$2,\$4}' $dfile | xyz2grd -G$grdfile $R $interp\n";
  print CSH "grdimage $grdfile -C$cptfile1 $J -Q -K -O -V $oxc >> $psfile\n";

  print CSH "grdcontour $grdfile -Cconfile1 -A- -W1.5p -Q100 $J -K -O -V $oxc >> $psfile\n";
      }

      print CSH "psscale -C${cptfile1} $Dscale $Bscale1 -K -O -V $ocb >> $psfile\n";
      $slab1 = "km/s";
      print CSH "pstext -N $R_title $J_title -K -O -V $ocbarlab >>$psfile<<EOF\n0 0 $fsize2 0 $fontno CM $slab1\nEOF\n";

      # details for specific cross-sections
      #if ($irun==1) {
  # positions of five faults along each profile
  # (1) SAF, (2) GF, (3) SGF, (4) MCF, (5) SYF, (6) CRF
  #(undef,undef,undef,$saf,$gf,$sgf,$mcf,$syf,$crf) = split(" ",$faults[$p-1]);
  #print "$p -- SAF $saf -- GF -- $gf -- SGF $sgf -- MCF $mcf -- SYF $syf\n"; die("TESTING");

  # vertical lines
      print CSH "psxy $J $R $finfo -K -O -V $oxc >>$psfile<<EOF\n$saf $zmax0\n$saf $zminf0\nEOF\n";
      print CSH "psxy $J $R $finfo -K -O -V $oxc >>$psfile<<EOF\n$gf $zmax0\n$gf $zminf0\nEOF\n";
      print CSH "psxy $J $R $finfo -K -O -V $oxc >>$psfile<<EOF\n$sgf $zmax0\n$sgf $zminf0\nEOF\n";
      print CSH "psxy $J $R $finfo -K -O -V $oxc >>$psfile<<EOF\n$mcf $zmax0\n$mcf $zminf0\nEOF\n";
      print CSH "psxy $J $R $finfo -K -O -V $oxc >>$psfile<<EOF\n$syf $zmax0\n$syf $zminf0\nEOF\n";
      print CSH "psxy $J $R $finfo -K -O -V $oxc >>$psfile<<EOF\n$crf $zmax0\n$crf $zminf0\nEOF\n";
      print CSH "psxy $J $R $finfo -K -O -V $oxc >>$psfile<<EOF\n$ef $zmax0\n$ef $zminf0\nEOF\n";
      print CSH "psxy $J $R $finfo -K -O -V $oxc >>$psfile<<EOF\n$snf $zmax0\n$snf $zminf0\nEOF\n";
      print CSH "psxy $J $R $finfo -K -O -V $oxc >>$psfile<<EOF\n$kcf $zmax0\n$kcf $zminf0\nEOF\n";

      #print CSH "psxy $J $R -W1p,0/0/0,-- -K -O -V $oxc >>$psfile<<EOF\n$dmin -10\n$dmax -10\nEOF\n";

      # plot source and receiver
      if ($irun == 1) {
  print CSH "awk '{print \$3,-\$4}' $rfile | psxy $J $R $srcinfo_xc -K -O -V $oxc >> $psfile\n";
  print CSH "awk '{print \$7,$zmax0}' $rfile | psxy $J $R $recinfo_xc -K -O -V $oxc >> $psfile\n";
      }

      # details for specific cross-sections
      #if ($irun==1) {
  #(undef,undef,undef,$saf,$gf,$sgf,$mcf,$syf) = split(" ",$faults[$p-1]);
      $talign = "CB";
      $flsize = $fsize1;
      $textinfo = "-N -C4p -W255/255/255o,1.0p,0/0/0,solid -G0/0/0";
      if($saf != 9999) {print CSH "pstext $J $R $textinfo -K -O -V $oxc >>$psfile<<EOF\n$saf $zmax0 $flsize 0 $fontno $talign SA\nEOF\n";}
      if($gf != 9999) {print CSH "pstext $J $R $textinfo -K -O -V $oxc >>$psfile<<EOF\n$gf $zmax0 $flsize 0 $fontno $talign G\nEOF\n";}
      if($sgf != 9999) {print CSH "pstext $J $R $textinfo -K -O -V $oxc >>$psfile<<EOF\n$sgf $zmax0 $flsize 0 $fontno $talign SG\nEOF\n";}
      if($mcf != 9999) {print CSH "pstext $J $R $textinfo -K -O -V $oxc >>$psfile<<EOF\n$mcf $zmax0 $flsize 0 $fontno $talign MC\nEOF\n";}
      if($syf != 9999) {print CSH "pstext $J $R $textinfo -K -O -V $oxc >>$psfile<<EOF\n$syf $zmax0 $flsize 0 $fontno $talign SY\nEOF\n";}
      if($crf != 9999) {print CSH "pstext $J $R $textinfo -K -O -V $oxc >>$psfile<<EOF\n$crf $zmax0 $flsize 0 $fontno $talign CR\nEOF\n";}
      if($ef != 9999) {print CSH "pstext $J $R $textinfo -K -O -V $oxc >>$psfile<<EOF\n$ef $zmax0 $flsize 0 $fontno $talign E\nEOF\n";}
      if($snf != 9999) {print CSH "pstext $J $R $textinfo -K -O -V $oxc >>$psfile<<EOF\n$snf $zmax0 $flsize 0 $fontno $talign SN\nEOF\n";}
      if($kcf != 9999) {print CSH "pstext $J $R $textinfo -K -O -V $oxc >>$psfile<<EOF\n$kcf $zmax0 $flsize 0 $fontno $talign KC\nEOF\n";}
      #}

    # plot overall label (for publication)
    if ($iletter==1) {
      $textinfo = "-N -C4p -W255/255/255o,1.0p,0/0/0,solid -G0/0/0";
      print CSH "pstext $R_title $J_title $textinfo -K -O -V $olab >>$psfile<<EOF\n0 0 18 0 $fontno TL $letter\nEOF\n";
    } else {
      if ($irun==1) {
  $textinfo = "-N";
  print CSH "pstext $R_title $J_title $textinfo -K -O -V $olab >>$psfile<<EOF\n0 0 $fsize1 0 $fontno TC $eid\nEOF\n";
  print CSH "pstext $R_title $J_title $textinfo -K -O -V $olab2 >>$psfile<<EOF\n0 0 $fsize1 0 $fontno TC $stanet\nEOF\n";
      }
    }

    # plot label on each cross section
    for ($k = 2; $k <= 2; $k ++ ) {
      $title = $titles[$k-1];
      $talign = "LB";
      #$otitle = "-Xa${x3} -Ya${y3}";
      $otitle = $oxclabs[$k-1];
      $textinfo = "-N -C4p -W255/255/255o,1.0p,0/0/0,solid -G0/0/0";
      print CSH "pstext $R_title $J_title $textinfo -K -O -V $otitle >>$psfile<<EOF\n0 0 $fsize0 0 $fontno $talign $title\nEOF\n";
    }

    #------------------------------------

    print CSH "pstext -N $R_title $J_title -O -V >>$psfile<<EOF\n 10 10 16 0 $fontno CM \nEOF\n"; # FINISH
    #if($ixv==1) {print CSH "convert $psfile -rotate $rotangle $jpgfile\n";}
    if ($ixv==1) {print CSH "gv $psfile &\n";}
    #if ($ipdf==1) {print CSH "ps2pdf $psfile\n";}
    if ($ieps==1) {print CSH "ps2eps -f -R + -l $psfile\n";}
    #if ($ijpg==1) {print CSH "convert -rotate $rotangle $psfile $jpgfile\n";}
    if ($ijpg==1) {print CSH "convert $epsfile $jpgfile\n";}

  }       # imfinal == 1
}       # loop over p

#==================================================

#------------------------------------

close (CSH);
system("csh -f $cshfile");

#   print "convert $psfile $jpgfile \n";
#   system("convert $psfile -rotate $rotangle $jpgfile");
#   if ($ipdf==1) {system("ps2pdf $psfile")}
#   if ($ixv==1) {system("gv $psfile &");}
system("gv $psfile &");

#==================================================
