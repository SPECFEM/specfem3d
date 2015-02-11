#!/usr/bin/perl -w

#==========================================================
#
#  plot_vert_models_basemap.pl
#  Carl Tape
#  02-April-2009
#
#  This script reads in a horizontal cross-section data file with five columns:
#     lon, lat, model-1, model-2, ln(model-2 / model-1)
#  and plots the three cross-sections.
#
#  EXAMPLE:
#    plot_vert_models_basemap.pl 1 2 3
#    plot_vert_models_basemap.pl 4 5 6 7 8 9 10 12 13 15 17 18 19 20 21 23 24 26 27 28 29 30 31 32   # SAF
#    plot_vert_models_basemap.pl 1 2 3 4 5 6 7 8 9 10  # Garlock
#    plot_vert_models_basemap.pl 1 2 3 4 5 6 7 8 9 10 11 12 13  # lon lines
#    plot_vert_models_basemap.pl 1 2 3 4 5 6 7 8 9  # lat lines
#
#    plot_vert_models_basemap.pl 4 5 10 25
#
#==========================================================

if (@ARGV < 1) {die("Usage: plot_vert_models_basemap.pl xxx\n");}
@pinds = @ARGV;
$nump = @pinds;
$stp = sprintf("%3.3i",$nump);

$irun = 5;        # KEY: 1, 2 (SAF), 3 (Garlock)
$icolor = 1;
$ivertlab = 0;    # vertical exaggeration label
$isimbox = 1;     # box around simulation region
$iletter = 0;        # label for publication figures
$letter = "A";
$imfinal = 0;     # plot a figure with only m16 and the map

#($pmin,$pmax) = split(/\//,$pinds);
#($smodel1,$smodel2) = split(/\//,$smodels);
#$smodeltag = "${smodel2}_${smodel1}";
#if($ivel==1) {$modlab = "vs"; $mtit = "Vs"; $vmin = 2.0; $vmax = 4.21; $icol = 3; $vtick = 1;}
#elsif($ivel==2) {$modlab = "vb"; $mtit = "Vb"; $vmin = 2.0; $vmax = 5.5; $icol = 6; $vtick = 1;}
#elsif($ivel==3) {$modlab = "vp"; $mtit = "Vp";  $vmin = 2.5; $vmax = 7.5; $icol = 9; $vtick = 1;}
#($zmin0,$zmax0) = split(/\//,$zran);

# directory containing the data files
$stirun = sprintf("%2.2i",$irun);
$dirdat = "INPUT/vert_${stirun}";


#---------------------------------------------------------

# subtitles
#@titles = ("$mtit  $smodel1","$mtit  $smodel2","ln(${smodel2} / ${smodel1})");

# parameters for color scales
#$cpert1 = 0.2;
#$cpert2 = 0.2;
$scale_color = 65.0;
$colorbar = "seis";

# KEY: resolution of color plots
$interp = "-I0.5/0.5 -S1.5";   # nearneighbor
#$interp = "-I1.0";            # xyz2grd
$grdfile = "temp.grd";

@Bopts = ("WESN","Wesn","wesN","wEsn","weSn","WesN","wEsN","wESn","WeSn","WEsn","weSN","wESN","WESn","WeSN","WEsN","wesn");

$dir0 = "/home/carltape/gmt";
$dir = "socal_2005";

$fault_file   = "$dir0/faults/jennings_more.xy";
$fault_file2  = "$dir0/faults/socal_faults_xc.lonlat";
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
$tick1 = 1; $tick2 = 0.2;
#$Rmap = "-R-121.2/-114.8/32.3/36.7";
$Rmap = "-R-122/-114/32/37";

$topo_labels  = "$dir0/socal_2005/socal_topo_labs.xyz";
$fault_labels = "$dir0/socal_2005/socal_fault_labs.xyz";
$coast_res = "-Df -A0/0/4";

$iportrait = 1;
if($iportrait == 1){$orient = "-P"; $rotangle = 0}
else {$orient = " "; $rotangle = 90}

# plotting specifications
$fsize0 = "14";
$fsize1 = "12";
$fsize2 = "11";
$fsize3 = "10";
$fontno = "1";    # 1 or 4
$tlen   = "6p";
$fpen   = "1.5p";
$tpen   = "1.5p";

# plotting specificiations
# NOTE: pen attribute (-W) does not work for psvelo
#$coast_info     = "-A5000 -Dl -W1.0p -S220/255/255 -C220/255/255 -G200";  # A : smallest feature plotted, in km^2; D : resolution
$coast_res      = "-Df -A0/0/4";
#$coast_info     = "${coast_res} -W0.75p -Na/0.75p";
$coast_info     = "${coast_res} -W1p -Na/1.0p -S220/255/255 -C220/255/255 -G200";  # -I1/1.0p for rivers
$coast_info2     = "${coast_res} -S220/255/255 -C220/255/255 -G200";  # -I1/1.0p for rivers
$fault_info     = "-M -W0.5p,255/0/0";
$fault_infoK    = "-M -W0.5p,0/0/0";
$fault_info2     = "-M -W2p,255/0/0";

$sky_color = "100/100/100";
$city_info = "-Sc10p -W1.0p/0/0/0 -G255/255/0";

#$colinfo = "-W1p,0/0/0 -G255/255/255";
$colinfo_map = "-W1.5p,0/0/0 -G0/255/255";
$srcinfo_map = "-Sa20p $colinfo_map -N";
$recinfo_map = "-Si20p $colinfo_map -N";
$srcinfo_mapinset = "-Sa10p -W0.5p,0/0/0 -G0/255/255 -N";
$recinfo_mapinset = "-Si10p -W0.5p,0/0/0 -G0/255/255 -N";
$colinfo_xc = "-W2p,0/0/0";
$srcinfo_xc  = "-Sa20p $colinfo_xc -N";
$recinfo_xc  = "-Si20p $colinfo_map -N";

$textinfo       = "-G255 -S1.5p";
$textinfo2      = "-G0/0/0 -S2p,255/255/255";
$textinfo3      = "-G0/0/0 -S2p,255/255/0";

# BOOLEAN: WHICH FIGURES TO PLOT
#$ixv = 1; $ipdf = 0;
$ixv = 0; $ipdf = 1;

# plot title
$J_title = "-JM1";
$R_title = "-R0/1/0/1";
$fsize_title = $fsize0;

#---------------------------------------------------------

# load values and check that the files exist
($nump0,undef,undef) = split(" ",`ls -1 ${dirdat}/*path | wc`);
if($nump0 == 0) {die("No ray paths available")}
#if($pmax > $nump) {$pmax = $nump;}
#if($pmin < 1) {$pmin = 1;}
print "\nNumber of total cuts is $nump0\n";

# load ALL fault positions along each profile
  $ffile = "$dirdat/ALL_${stirun}_rays_fault_positions_mod";
  if (not -f $ffile) {die("Check if ffile $ffile exist or not\n");}
  open(IN,$ffile); @faults = <IN>; close(IN);

  $widmap = 6.5; # width of map
  $Jmap = "-JM${widmap}";
  $B0map = "-Ba${tick1}f${tick2}d:.\" \":";
  #$Bmap = "$B0map".$Bopts[7];
  $Bmap = "$B0map".$Bopts[0];

  # parameters controlling the placement of subplots and labels
  $omap = "-Xa1 -Ya3";

  $fname = "vert_${stirun}_basemap_${stp}";
  $psfile = "$fname.ps"; $jpgfile = "$fname.jpg";

#==============================================

$cshfile = "plot_vert_models_basemap.csh";
open(CSH,">$cshfile");
print CSH "gmtset COLOR_NAN $sky_color PAPER_MEDIA letter BASEMAP_TYPE plain PLOT_DEGREE_FORMAT D TICK_LENGTH $tlen LABEL_FONT_SIZE $fsize2 ANOT_FONT_SIZE $fsize2  HEADER_FONT $fontno ANOT_FONT $fontno LABEL_FONT $fontno HEADER_FONT_SIZE $fsize1 FRAME_PEN $fpen TICK_PEN $tpen \n";

  print CSH "psbasemap $Jmap $Rmap $Bmap -K -V $orient $omap > $psfile\n"; # START

  print CSH "pscoast $Jmap $Rmap $coast_info -K -O -V $omap >> $psfile \n";
  if ($ishelf==1) {
    print CSH "awk '{print \$2,\$1}' ${shelf_file} | psxy $Jmap $Rmap $Wshelf -K -O -V $omap >> $psfile\n";
  }

  print CSH "psxy ${fault_file} $Jmap $Rmap $fault_info -K -V -O $omap >> $psfile \n";
  print CSH "psxy ${fault_file2} $Jmap $Rmap $fault_info2 -K -V -O $omap >> $psfile \n";

  if ($isimbox==1) {
    $outer_boundary = "${dir_sources}/SPECFEM_outer_points.dat";
    if (not -f ${outer_boundary}) {
      die("Check if outer_boundary ${outer_boundary} exist or not\n");
    }
    print CSH "psxy ${outer_boundary} $Jmap $Rmap -W2p,0/0/0 -K -O -V $omap >> $psfile\n";
  }

#==============================================

# modify the height of each cross-section based on the length
#for ($p = $pmin; $p <= $pmax; $p ++ ) {
for ($ip = 1; $ip <= $nump; $ip ++ ) {
  $p = $pinds[$ip-1];
  $stip = sprintf("%3.3i",$p);
  print "-- $ip -- $nump -- $p --\n";

  # source and receiver points -- NOT necessarily the endpoints of the ray path
  # COLUMNS: slons(ii),slats(ii),0,rlons(ii),rlats(ii),deg2km(dq),azq,iflip
  $rfile = "${dirdat}/vert_${stirun}_xc_${stip}_A_B";
  if (not -f $rfile) {die("Check if rfile $rfile exist or not\n");}
  open(IN,$rfile); @line = <IN>; close(IN);
  ($slon,$slat,$sdist,$sdep,$rlon,$rlat,$rdist,$rdep,$raz,$iflip) = split(" ",$line[0]);

   # positions of six faults along each profile
   # (1) SAF, (2) GF, (3) SGF, (4) MCF, (5) SYF, (6) CRF
   (undef,$eid,$stanet,$azi1,$azi2,$saf,$gf,$sgf,$mcf,$syf,$crf) = split(" ",$faults[$p-1]);
   ($sta,$net) = split("\\.",$stanet);
   $sazi1 = sprintf("%3.3i",$azi1);
   $sazi2 = sprintf("%3.3i",$azi2);
   print "-- $stanet -- $sta -- $net -- $eid -- $sazi1 -- $sazi2 --\n";

  #---------------------------------------------------------

  # plot ray path
  $pfile = "${dirdat}/vert_${stirun}_xc_${stip}_ray_path";
  if (not -f $pfile) {
    die("Check if pfile $pfile exist or not\n");
  }
  print CSH "awk '{print \$1,\$2}' ${pfile} | psxy $Jmap $Rmap -W6p,0/0/0 -K -O -V $omap >> $psfile\n";
  print CSH "awk '{print \$1,\$2}' ${pfile} | psxy $Jmap $Rmap -W3p,0/255/255 -K -O -V $omap >> $psfile\n";

  # plot source and receiver lon-lat
  $rfile = "${dirdat}/vert_${stirun}_xc_${stip}_A_B";
  if (not -f $rfile) {
    die("Check if rfile $rfile exist or not\n");
  }

  if ($irun==1) {
  # plot source
  print CSH "awk '{print \$1,\$2}' ${rfile} | psxy $Jmap $Rmap $srcinfo_map -K -O -V $omap >> $psfile\n";

  # plot receiver
  print CSH "awk '{print \$5,\$6}' ${rfile} | psxy $Jmap $Rmap $recinfo_map -K -O -V $omap >> $psfile\n";
}

  # plot cities
  ##print CSH "psxy $Jmap $Rmap $city_info -K -O -V  >>$psfile<<EOF\n -118.1717 34.1483 \nEOF\n";   # Pasadena
  #print CSH "awk '{print \$2,\$1}' $cities | psxy ${city_info} $Jmap $Rmap -K -V -O >> $psfile \n";
  #print CSH "awk '{print \$2,\$1,12,0,1,\"LB\",\$3}' $cities | pstext $textinfo3 -D0.05/0.05 $Jmap $Rmap -K -V -O >> $psfile \n";

  #print CSH "pstext ${topo_labels} $Jmap $Rmap $textinfo -N -K -V -O >> $psfile\n";
  #print CSH "pstext ${fault_labels} $textinfo2 $Jmap $Rmap -K -V -O >> $psfile \n";

#     # plot overall label (for publication)
#     if ($iletter==1) {
#       $textinfo = "-N -C4p -W255/255/255o,1.0p,0/0/0,solid -G0/0/0";
#       print CSH "pstext $R_title $J_title $textinfo -K -O -V $olab >>$psfile<<EOF\n0 0 18 0 $fontno TL $letter\nEOF\n";
#     } else {
#       $textinfo = "-N";
#       print CSH "pstext $R_title $J_title $textinfo -K -O -V $olab >>$psfile<<EOF\n0 0 10 0 $fontno TC $eid\nEOF\n";
#       print CSH "pstext $R_title $J_title $textinfo -K -O -V $olab2 >>$psfile<<EOF\n0 0 10 0 $fontno TC $stanet\nEOF\n";
#     }

}       # loop over p

  # fault labels
  $fault_labs = "$dir0/faults/socal_fault_labs.lonlat";
  $textinfo = "-N -C4p -W255/255/255o,1.0p,255/0/0,solid -G0/0/0";
  print CSH "awk '{print \$1,\$2,10,0,1,\"CM\",\$3}' ${fault_labs} | pstext $textinfo $Jmap $Rmap -K -V -O $omap >> $psfile\n";

  print CSH "pstext -N $R_title $J_title -K -O -V >>$psfile<<EOF\n 10 10 16 0 $fontno CM \nEOF\n";
  print CSH "psbasemap $Jmap $Rmap $Bmap -O -V $omap >> $psfile\n"; # FINISH (no K)

  #if($ixv==1) {print CSH "convert $psfile -rotate $rotangle $jpgfile\n";}
  if ($ixv==1) {print CSH "gv $psfile &\n";}
  if ($ipdf==1) {print CSH "ps2pdf $psfile\n";}

#------------------------------------

close (CSH);
system("csh -f $cshfile");

system("gv $psfile &")

#==================================================
