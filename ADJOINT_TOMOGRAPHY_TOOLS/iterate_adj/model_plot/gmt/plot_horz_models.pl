#!/usr/bin/perl -w

#==========================================================
#
#  plot_horz_models.pl
#  Carl Tape
#  07-Dec-2008 (last checked 10-Jan-2013 with GMT 4.5.3)
#
#  This script reads in a horizontal cross-section data file with five columns:
#     lon, lat, model-1, model-2, ln(model-2 / model-1)
#  and plots the three cross-sections.
#
#  EXAMPLE:
#    plot_horz_models.pl 1/81 m00/m16 1 1/1/1/0.2      # Vs
#    plot_horz_models.pl 1/81 m00/m16 0 1/1/1/0.3      # Vb
#
#  FIGURES:
#    plot_horz_models.pl 5/5   m00/m16 1 1/0/1/0.2      # Vs
#    plot_horz_models.pl 21/21 m00/m16 1 0/0/1/0.2      # Vs TEST
#    plot_horz_models.pl 41/41 m00/m16 1 0/0/1/0.2      # Vs
#    plot_horz_models.pl 5/5   m00/m16 0 1/0/1/0.3      # Vb
#    plot_horz_models.pl 21/21 m00/m16 0 0/0/1/0.3      # Vb
#    plot_horz_models.pl 41/41 m00/m16 0 0/0/1/0.3      # Vb
#
#    plot_horz_models.pl 3/3   m00/m01 1 1/0/1/0.2      # Vs
#
#==========================================================

if (@ARGV < 4) {die("Usage: plot_horz_models.pl xxx\n");}
($pinds,$smodels,$ivs,$parms) = @ARGV;

$icolor = 1;
($pmin,$pmax) = split(/\//,$pinds);
($smodel1,$smodel2) = split(/\//,$smodels);
$smodeltag = "${smodel2}_${smodel1}";
if($ivs==1) {$modlab = "vs"; $mtit = "Vs";} else {$modlab = "vb"; $mtit = "Vb";}
($ilonlab,$itoplab,$isidelab,$ctick1b_km) = split(/\//,$parms);

# detailed flags for publication figures
#$ilonlab = 0;      # display longitude tick labels at the top
#$itoplab = 0;      # label the model above the plot
$iinsetlab = 1;     # label the model inan inset
$iletter = 0;       # label for publication figures
$letter = "C";
$iextra = 0;        # extra information on figures
$irun = 2;          # default 2; 1 for m01/m00
$imask = 1;
$isingle = 0;       # single figure showing only final model

$cgray = "-G220/220/220";

$stirun = sprintf("%2.2i",$irun);

#---------------------------------------------------------
# USER INPUT

# directory containing the data files
$dirdat = "INPUT/horz_${stirun}";

# KEY: file showing the cuts and the plotting range for the values
$fcuts = "${dirdat}/horz_xc_cuts_mod";
if (not -f $fcuts) {die("Check if fcuts $fcuts exist or not\n")}

# load values and check that the files exist
open(IN,$fcuts); @lines = <IN>; $nump = @lines;
for ($p = 1; $p <= $nump; $p ++ ) {
   $stip = sprintf("%3.3i",$p);
   (undef,$zcut,$vs_norm,$vb_norm,$vpert) = split(" ",$lines[$p-1]);
   $vperts[$p-1] = $vpert;
   $vsnorms[$p-1] = $vs_norm;
   $vbnorms[$p-1] = $vb_norm;
   $zcuts[$p-1]  = -$zcut;
   $dfiles[$p-1] = "${dirdat}/horz_xc_${modlab}_${smodeltag}_${stip}_mask${imask}.dat";
   #if (not -f $dfile) {die("Check if dfile $dfile exist or not\n")}
}
print "\nNumber of cuts is $nump\n";
if($pmax > $nump) {$pmax = $nump;}
if($pmin < 1) {$pmin = 1;}

# subtitles
@titles = ("$mtit  model  $smodel1","$mtit  model  $smodel2","ln ( ${smodel2} / ${smodel1} )");
@labs = ("$smodel1","$smodel2","ln(${smodel2}/${smodel1})");

# parameters for color scales
$cpert1 = 0.2;
#$cpert2 = 0.2;
$scale_color = 65.0;
$colorbar = "seis";

# resolution of color plots
#$interp = "-I2m/2m -S4m";   # key information
$interp = "-I2m"; 
$grdfile = "temp.grd";

# position of titles and labels
$xtx1 = 0.5;   $ytx1 = 1.00;
$xtx2 = -0.15; $ytx2 = 0.4;
$xtx3 = 0.7;   $ytx3 = -0.19;    # cbar 1
$xtx4 = $xtx3; $ytx4 = -0.1;    # cbar 2
$xtx5 = -0.15; $ytx5 = 0.85;    # A, B, C, etc
$xtx6 = 0.0; $ytx6 = 0.05;    # m00, m16, ln(m16/m00)

#print "$dtitle\n"; die("TESTING");

#---------------------------------------------------------
# DATA FILES (MODIFY FOR EACH USER)

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

# check that files exist
#if (not -f $xxx)        { die("Check if $xxx exist or not\n") }

# bounds and dimensions
# 0 = LA
# 1 = southern California
$iregion = 1;

if ($iregion==1) {
  # southern California
  $wid = 2.8;			# width of figure (inches)
  $xmin = -122; $xmax = -114; $ymin = 31.5; $ymax = 37;
  $xmin = -121.2; $xmax = -114.8; $ymin = 32.3; $ymax = 36.7;   # SPECFEM
  $tick1 = 1; $tick2 = 0.5;
  $cmax = 5000; $cmin = -$cmax;
  $y_title = 0.98; $y_title2 = 0.92; $x_title = 0.5; $x_title2 = 0.5;
  $iportrait = 0;
  $stag = "Southern California";
  $topo_labels  = "$dir0/socal_2005/socal_topo_labs.xyz";
  $fault_labels = "$dir0/socal_2005/socal_fault_labs.xyz";
  $name = "plot_horz_models";
  $origin = "-X1i -Y5i"; 

  # inset map
  $iinset = 1;
  $origin_inset = "-Xa7.5 -Ya4.5";
  $Jinset = "-JM1.5";
  $Rinset = "-R-132/-110/25/50";
  $Binset = "-B100wesn"; 
  $coast_res = "-Df -A0/0/4";
}

$dX = $wid + 0.5;
$hwid = $wid/2;
#$R = "-Rg";
$R = "-R$xmin/$xmax/$ymin/$ymax";
#$J = "-JK180/${wid}i";
$J = "-JM${wid}";
#$J = "-Jm0.70";

# scale for plots
@Bopts = ("WESN","Wesn","wesN","wEsn","weSn","WesN","wEsN","wESn","WeSn","WEsn","weSN","wESN","WESn","WeSN","WEsN","wesn");
$B0 = "-Ba${tick1}f${tick2}d:.\" \":";
#$B0 = "-Ba${tick1}f${tick1}g${tick1}:.\" \":";   # with gridlines

$B = "$B0".$Bopts[0];

if($iportrait == 1){$orient = "-P"; $rotangle = 0}
else {$orient = " "; $rotangle = 90}

# color bar
$Dwid = 0.15;
$Dlen = $wid*0.6;
$Dx = $Dlen*0.5;
$Dy = -0.2;
$Dy2 = $Dy-$Dwid;
$Dy2 = $Dy;
$Dscale = "-D$Dx/$Dy/$Dlen/${Dwid}h";
$Dscale2 = "-D$Dx/$Dy2/$Dlen/${Dwid}h";

# plotting specifications
$fsize0 = "14";
$fsize1 = "12";
$fsize2 = "10";
$fsize3 = "8";
$fontno = "1";    # 1 or 4
$tick   = "6p";
$fpen   = "1p";
$tpen   = "1p";

# plotting specificiations
# NOTE: pen attribute (-W) does not work for psvelo
#$coast_info     = "-A5000 -Dl -W1.0p -S220/255/255 -C220/255/255 -G200";  # A : smallest feature plotted, in km^2; D : resolution

$coast_info     = "${coast_res} -W0.75p -Na/0.75p";
$fault_info     = "-m -W0.75p,255/0/0";
$fault_infoK    = "-m -W0.75p,0/0/0";

$city_info = "-Sc10p -W1.0p/0/0/0 -G255/255/0";

$textinfo       = "-G255 -S1.5p";
$textinfo2      = "-G0/0/0 -S2p,255/255/255";
$textinfo3      = "-G0/0/0 -S2p,255/255/0";

# BOOLEAN: WHICH FIGURES TO PLOT
$ixv = 1; $ipdf = 1;
#$ixv = 1; $ipdf = 0;

# plot title
$J_title = "-JM${wid}";
$R_title = "-R0/1/0/1";
$fsize_title = $fsize0;

#==================================================

$cshfile = "plot_horz_models.csh";
open(CSH,">$cshfile");
print CSH "gmtset PAPER_MEDIA letter MEASURE_UNIT inch BASEMAP_TYPE plain PLOT_DEGREE_FORMAT D TICK_LENGTH $tick LABEL_FONT_SIZE $fsize2 ANOT_FONT_SIZE $fsize2  HEADER_FONT $fontno ANOT_FONT $fontno LABEL_FONT $fontno HEADER_FONT_SIZE $fsize1 FRAME_PEN $fpen TICK_PEN $tpen \n";

@shifts = ($origin,"-X$dX","-X$dX");
if($ilonlab==1) {@iBs = (2,14,2);} else {@iBs = (15,9,15);}

for ($p = $pmin; $p <= $pmax; $p ++ ) {

  $stip = sprintf("%3.3i",$p);
  $fname0 = "horz_xc_${modlab}_${smodeltag}_${stip}_extra${iextra}";
  $fname = $fname0;
  $psfile = "$fname.ps"; $jpgfile = "$fname.jpg";

  $vsnorm = $vsnorms[$p-1];	# reference velocity for m00 and m16
  $vbnorm = $vbnorms[$p-1];	# reference velocity for m00 and m16
  $vpert = $vperts[$p-1];       # perturbation from reference velocity to plot
  $zdep = $zcuts[$p-1];		# depth of the cross-section (m)

  # make less of a range for Vb compared with Vs -- factor
  if($ivs == 0) {$vpert = 1.0*$vpert;}

  $dtitle = sprintf("Depth  =  %.1f km",$zdep/1000);

  # datafile -- COLUMNS: lon, lat, m00, m16
  $dfile = $dfiles[$p-1];
  if (not -f $dfile) {die("Check if dfile $dfile exist or not\n");}
  print "DATA FILE: $dfile\n";

  if($ivs == 1) {$cnorm = $vsnorm;} else {$cnorm = $vbnorm;}
  $cpert1 = $vpert;
  $cpert2 = $vpert;    # NOTE: over-rule the input value

  # colorpoint file for m00 and m16 -- VARIES with depth -- in ln(m00/cnorm)
  $cptfile1a = "color1a.cpt";
  $cmin = -$cpert1*1.01; $cmax = $cpert1*1.01; 
  $dc = ($cmax-$cmin)/${scale_color};
  $T = sprintf("-T%3.3e/%3.3e/%3.3e",$cmin,$cmax,$dc);
  print CSH "makecpt -C$colorbar $T -D > $cptfile1a\n";

  # colorpoint file for m00 and m16 -- VARIES with depth -- in KILOMETERS
  # THIS IS FOR THE COLORBAR ONLY
  $cptfile1b = "color1b.cpt";
  $cmin = $cnorm*exp(-$cpert1)*1e-3;
  $cmax = $cnorm*exp($cpert1)*1e-3;
  $dc = ($cmax-$cmin)/${scale_color};
  $T = sprintf("-T%3.3e/%3.3e/%3.3e",$cmin,$cmax,$dc);
  print CSH "makecpt -C$colorbar $T -D > $cptfile1b\n";
  $Bscale1b = sprintf("-B%2.2ef0.1:\"  \": -E10p",$ctick1b_km);

  # colorpoint file for perturbation ln(m16/m00) -- FIXED with depth
  $cptfile2 = "color2.cpt";
  $cmin = -$cpert2*1.01; $cmax = $cpert2*1.01; 
  $dc = ($cmax-$cmin)/${scale_color};
  $T = sprintf("-T%3.3e/%3.3e/%3.3e",$cmin,$cmax,$dc);
  print CSH "makecpt -C$colorbar $T -D > $cptfile2\n";

  # ticks for colorbars
  $Bscale1a  = sprintf("-Ba%2.2ef0.05:\" \": -E10p -A",$cpert1);
  $Bscale2  = sprintf("-Ba%2.2ef0.05:\" \": -E10p",$cpert2);

# attempting to get shaded relief superimposed
#$gradfile0 = "/home/datalib/Topography/GLOBAL/ETOPO1/ETOPO1_Ice_g.grad";
#if (not -f $gradfile0) {die("Check if gradfile $gradfile0 exist or not\n");}

  # LOOP over initial model, final model, and logarithmic difference
  for ($k = 1; $k <= 3; $k ++ ) {
  #for ($k = 1; $k <= 1; $k ++ ) {

    $shift = $shifts[$k-1];
    $iB = $iBs[$k-1];
    $B = "$B0".$Bopts[$iB];
    $title = $titles[$k-1];

    if ($k==1) {
      print CSH "psbasemap $J $R $B -K -V $orient $origin > $psfile\n";	# START
    } else {
      print CSH "psbasemap $J $R $B -K -O -V $shift >> $psfile\n";
    }
    print CSH "psbasemap $J $R $B $cgray -K -O -V >> $psfile\n";  # mask color
    if ($icolor==1) {
      if ($k == 1) {
	##print CSH "awk '{print \$1,\$2,log(\$3/$cnorm)}' $dfile | nearneighbor -G$grdfile $R $interp\n";
        #print CSH "awk '{print \$1,\$2,log(\$3/$cnorm)}' $dfile | xyz2grd -G$grdfile $R $interp\n";
	#print CSH "grdimage $grdfile -C$cptfile1a $J -Q -K -O -V >> $psfile\n";
        #print CSH "awk '{print \$1,\$2,\$3/1000}' $dfile | xyz2grd -G$grdfile $R $interp\n";
        print CSH "grep -v NaN $dfile | awk '{print \$1,\$2,\$3/1000}' | xyz2grd -G$grdfile $R $interp\n";
	print CSH "grdimage $grdfile -C$cptfile1b $J -Q -K -O -V >> $psfile\n";

        #print CSH "psscale -C${cptfile1a} $Dscale $Bscale1a -K -O -V >> $psfile\n";
        print CSH "psscale -C${cptfile1b} $Dscale2 $Bscale1b -K -O -V >> $psfile\n";
        #$slab1 = sprintf("ln(%s / %.2f)",$mtit,$cnorm/1000);
        $slab1 = sprintf("(%.2f \\261 %.0f \\045)",$cnorm/1000,$cpert1*100);

        $slab2 = "$mtit  km/s";
        print CSH "pstext -N $R_title $J_title -K -O -V >>$psfile<<EOF\n$xtx3 $ytx3 $fsize2 0 $fontno LM $slab1\nEOF\n";
        print CSH "pstext -N $R_title $J_title -K -O -V >>$psfile<<EOF\n$xtx4 $ytx4 $fsize2 0 $fontno LM $slab2\nEOF\n";

      } elsif ($k == 2) {
	##print CSH "awk '{print \$1,\$2,log(\$4/$cnorm)}' $dfile | nearneighbor -G$grdfile $R $interp\n";
        #print CSH "awk '{print \$1,\$2,log(\$4/$cnorm)}' $dfile | xyz2grd -G$grdfile $R $interp\n";
	#print CSH "grdimage $grdfile -C$cptfile1a $J -Q -K -O -V >> $psfile\n";
        print CSH "grep -v NaN $dfile | awk '{print \$1,\$2,\$4/1000}' | xyz2grd -G$grdfile $R $interp\n";
	print CSH "grdimage $grdfile -C$cptfile1b $J -Q -K -O -V >> $psfile\n";

        #print CSH "psscale -C${cptfile1a} $Dscale $Bscale1a -K -O -V >> $psfile\n";
        print CSH "psscale -C${cptfile1b} $Dscale2 $Bscale1b -K -O -V >> $psfile\n";
        #$slab1 = sprintf("ln(%s / %.2f)",$mtit,$cnorm/1000);
        $slab1 = sprintf("(%.2f \\261 %.0f \\045)",$cnorm/1000,$cpert1*100);
        $slab2 = "$mtit  km/s";
        print CSH "pstext -N $R_title $J_title -K -O -V >>$psfile<<EOF\n$xtx3 $ytx3 $fsize2 0 $fontno LM $slab1\nEOF\n";
        print CSH "pstext -N $R_title $J_title -K -O -V >>$psfile<<EOF\n$xtx4 $ytx4 $fsize2 0 $fontno LM $slab2\nEOF\n";

      } elsif ($k == 3) {
	#print CSH "awk '{print \$1,\$2,log(\$4/\$3)}' $dfile | nearneighbor -G$grdfile $R $interp\n";
        #print CSH "awk '{print \$1,\$2,\$5}' $dfile | nearneighbor -G$grdfile $R $interp\n";
        print CSH "awk '{print \$1,\$2,\$5}' $dfile | xyz2grd -G$grdfile $R $interp\n";
	print CSH "grdimage $grdfile -C$cptfile2 $J -K -O -V -Q >> $psfile\n";

        print CSH "psscale -C${cptfile2} $Dscale $Bscale2 -K -O -V >> $psfile\n";
        $slab1 = "ln($smodel2 / $smodel1)";
        print CSH "pstext -N $R_title $J_title -K -O -V >>$psfile<<EOF\n$xtx4 $ytx4 $fsize2 0 $fontno LM $slab1\nEOF\n";
      }
    }
    print CSH "pscoast $J $R $coast_info -K -O -V >> $psfile \n";
    if($ishelf==1) {print CSH "awk '{print \$2,\$1}' ${shelf_file} |psxy $J $R $Wshelf -K -O -P -V >> $psfile\n";}

    #--------------------------------

    #print CSH "psxy ${plate_file} $J $R $plate_info -K -V -O >> $psfile \n";
    print CSH "psxy ${fault_file} $J $R $fault_infoK -K -V -O >> $psfile \n";
    #print CSH "psxy ${kcf_file} $J $R $fault_infoK -K -V -O >> $psfile \n";
    #print CSH "psxy ${breck_file} $J $R $fault_infoK -K -V -O >> $psfile \n";
    #print CSH "psxy ${fault_file2} $J $R $fault_infoK -K -V -O >> $psfile \n";

    if (0==1) {
      # boundaries of simulation
      print CSH "psxy ${outer_boundary} $J $R -W2p,0/0/0 -K -O -V >>$psfile\n";  
      print CSH "psxy ${inner_boundary} $J $R -W1.5p,0/0/0,-- -K -O -V >>$psfile\n";  

      $ibox = 0;
      if ($ibox==1) {
	$boxinfo = "-W1.5p,0/0/255"; # -A : suppress drawing line segments as great circle arcs

	# Lin model (2007) box
	$lin_boundary = "/net/denali/home2/carltape/gmt/tomography/lin_2007/lin_boundary_points.dat";
	print CSH "psxy ${lin_boundary} $J $R $boxinfo -K -O -V >>$psfile\n";  
      }
    }

    #--------------------------------

    if ($iextra==1) {
      $finfo1 = "-W3p,0/0/0";
      $finfo2 = "-W1.5p,0/255/255";

      # plot a line at -119 longitude
      if ($p==11 && $ivs==1) {
	$xmark = -119;
	print CSH "psxy $J $R $finfo1 -K -O -V >>$psfile<<EOF\n$xmark $ymin\n$xmark $ymax\nEOF\n";
	print CSH "psxy $J $R $finfo2 -K -O -V >>$psfile<<EOF\n$xmark $ymin\n$xmark $ymax\nEOF\n";
      }

      # plot two cross section ray paths
      if ($p==21 && $ivs==1) {
	@irs = (2,5); $nray = @irs ;
	for ($ik = 1; $ik <= $nray; $ik ++ ) {
	  $ir = $irs[$ik-1];
	  $stir = sprintf("%3.3i",$ir);
	  $pfile = "./INPUT/vert_01/vert_xc_${stir}_ray_path";
	  if (not -f $pfile) {
	    die("Check if pfile $pfile exist or not\n");
	  }
	  print CSH "awk '{print \$1,\$2}' ${pfile} | psxy $J $R $finfo1 -K -O -V >> $psfile\n";
	  print CSH "awk '{print \$1,\$2}' ${pfile} | psxy $J $R $finfo2 -K -O -V >> $psfile\n";
	}
      }
    }

    #--------------------------------

    print CSH "psbasemap $J $R $B -K -V -O >> $psfile\n";

    # plot cities
    ##print CSH "psxy $J $R $city_info -K -O -V  >>$psfile<<EOF\n -118.1717 34.1483 \nEOF\n";   # Pasadena
    #print CSH "awk '{print \$2,\$1}' $cities | psxy ${city_info} $J $R -K -V -O >> $psfile \n";
    #print CSH "awk '{print \$2,\$1,12,0,1,\"LB\",\$3}' $cities | pstext $textinfo3 -D0.05/0.05 $J $R -K -V -O >> $psfile \n";

    #print CSH "pstext ${topo_labels} $J $R $textinfo -N -K -V -O >> $psfile\n";
    #print CSH "pstext ${fault_labels} $textinfo2 $J $R -K -V -O >> $psfile \n";

    #----------------------------

    # plot horizontal title above each column
    if ($itoplab==1) {
      print CSH "pstext -N $R_title $J_title -K -O -V >>$psfile<<EOF\n$xtx1 $ytx1 $fsize0 0 $fontno CM $title\nEOF\n";
    }
 
    # inset label for each plot
    if ($iinsetlab==1) {
      $lab = $labs[$k-1];
      $textinfo = "-N -C4p -W255/255/255o,1.0p,0/0/0,solid -G0/0/0";
      print CSH "pstext $R_title $J_title $textinfo -K -O -V >>$psfile<<EOF\n$xtx6 $ytx6 11 0 $fontno LB $lab\nEOF\n";
    }

    # plot vertical title left of the row
    if ($isidelab==1 && $k==1) {
      print CSH "pstext -N $R_title $J_title -K -O -V >>$psfile<<EOF\n$xtx2 $ytx2 $fsize0 90 $fontno CM $dtitle\nEOF\n";
    } 

    # plot overall label (for publication)
    if ($iletter==1 && $k==1) {
      $textinfo = "-N -C4p -W255/255/255o,1.0p,0/0/0,solid -G0/0/0";
      print CSH "pstext $R_title $J_title $textinfo -K -O -V >>$psfile<<EOF\n$xtx5 $ytx5 18 0 $fontno TL $letter\nEOF\n";
    }

  }				# loop over k

  print CSH "pstext -N $R_title $J_title -O -V >>$psfile<<EOF\n $x_title $y_title 16 0 $fontno CM \nEOF\n"; # FINISH 
  #if($ixv==1) {print CSH "convert $psfile -rotate $rotangle $jpgfile\n";}
  #if($ixv==1) {print CSH "gv $psfile &\n";}
  if($ipdf==1) {print CSH "ps2pdf $psfile\n";}

  #-------------------------------------------

  if ($isingle == 1) {

  $fname = "${fname0}_single";
  $psfile = "$fname.ps"; $jpgfile = "$fname.jpg";

  $k = 2;

    $shift = $shifts[$k-1];
    $B = "$B0"."wENs";
    $title = $titles[$k-1];

    print CSH "psbasemap $J $R $B -K -V $orient $origin > $psfile\n";	# START
   
    print CSH "psbasemap $J $R $B $cgray -K -O -V >> $psfile\n";
    if ($icolor==1) {
        print CSH "awk '{print \$1,\$2,\$4/1000}' $dfile | xyz2grd -G$grdfile $R $interp\n";
	print CSH "grdimage $grdfile -C$cptfile1b $J -Q -K -O -V >> $psfile\n";
        print CSH "psscale -C${cptfile1b} $Dscale2 $Bscale1b -K -O -V >> $psfile\n";
        $slab1 = sprintf("(%.2f \\261 %.0f \\045)",$cnorm/1000,$cpert1*100);
        $slab2 = "$mtit  km/s";
        print CSH "pstext -N $R_title $J_title -K -O -V >>$psfile<<EOF\n$xtx3 $ytx3 $fsize2 0 $fontno LM $slab1\nEOF\n";
        print CSH "pstext -N $R_title $J_title -K -O -V >>$psfile<<EOF\n$xtx4 $ytx4 $fsize2 0 $fontno LM $slab2\nEOF\n";
    }
    print CSH "pscoast $J $R $coast_info -K -O -V >> $psfile \n";
    if($ishelf==1) {print CSH "awk '{print \$2,\$1}' ${shelf_file} |psxy $J $R $Wshelf -K -O -P -V >> $psfile\n";}
    print CSH "psxy ${fault_file} $J $R $fault_infoK -K -V -O >> $psfile \n";

    print CSH "psbasemap $J $R $B -K -V -O >> $psfile\n";

    # plot cities
    ##print CSH "psxy $J $R $city_info -K -O -V  >>$psfile<<EOF\n -118.1717 34.1483 \nEOF\n";   # Pasadena
    #print CSH "awk '{print \$2,\$1}' $cities | psxy ${city_info} $J $R -K -V -O >> $psfile \n";
    #print CSH "awk '{print \$2,\$1,12,0,1,\"LB\",\$3}' $cities | pstext $textinfo3 -D0.05/0.05 $J $R -K -V -O >> $psfile \n";

    #print CSH "pstext ${topo_labels} $J $R $textinfo -N -K -V -O >> $psfile\n";
    #print CSH "pstext ${fault_labels} $textinfo2 $J $R -K -V -O >> $psfile \n";

     # fault labels
     $fault_labs = "$dir0/faults/socal_fault_labs_red.lonlat";
     $textinfo = "-N -C2p -W255/255/255o,1.0p,0/0/0,solid -G0/0/0";
     print CSH "awk '{print \$1,\$2,10,0,1,\"CM\",\$3}' ${fault_labs} | pstext $textinfo $J $R -K -V -O >> $psfile\n";

    #----------------------------

    # plot horizontal title above each column
    if ($itoplab==1) {
      print CSH "pstext -N $R_title $J_title -K -O -V >>$psfile<<EOF\n$xtx1 $ytx1 $fsize0 0 $fontno CM $title\nEOF\n";
    }
 
    # inset label for each plot
    if ($iinsetlab==1) {
      $lab = $labs[$k-1];
      $textinfo = "-N -C4p -W255/255/255o,1.0p,0/0/0,solid -G0/0/0";
      print CSH "pstext $R_title $J_title $textinfo -K -O -V >>$psfile<<EOF\n$xtx6 $ytx6 11 0 $fontno LB $lab\nEOF\n";
    }

    # plot vertical title left of the row
    if ($isidelab==1) {
      print CSH "pstext -N $R_title $J_title -K -O -V >>$psfile<<EOF\n$xtx2 $ytx2 $fsize0 90 $fontno CM $dtitle\nEOF\n";
    } 

    # plot overall label (for publication)
    if ($iletter==1) {
      $textinfo = "-N -C4p -W255/255/255o,1.0p,0/0/0,solid -G0/0/0";
      print CSH "pstext $R_title $J_title $textinfo -K -O -V >>$psfile<<EOF\n$xtx5 $ytx5 18 0 $fontno TL $letter\nEOF\n";
    }

  print CSH "pstext -N $R_title $J_title -O -V >>$psfile<<EOF\n $x_title $y_title 16 0 $fontno CM \nEOF\n"; # FINISH 
  if($ixv==1) {print CSH "convert $psfile -rotate $rotangle $jpgfile\n";}
  if($ixv==1) {print CSH "gv $psfile &\n";}
  if($ipdf==1) {print CSH "ps2pdf $psfile\n";}

  }

}				# loop over p
  
#------------------------------------
close (CSH);
system("csh -f $cshfile");

#   print "convert $psfile $jpgfile \n";
#   system("convert $psfile -rotate $rotangle $jpgfile");
#   if ($ipdf==1) {system("ps2pdf $psfile")}
#   if ($ixv==1) {system("gv $psfile &");}

system("gv $psfile &");

#==================================================
