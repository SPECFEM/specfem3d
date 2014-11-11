#!/usr/bin/perl -w

#==========================================================
#
#  plot_horz_coverage.pl
#  Carl Tape
#  16-Nov-2009
#
#  horz_xc_vs_coverage_003.dat
#  This script reads in a horizontal cross-section data file with five columns:
#     lon, lat, log10(Ksum), log10(Ksum)+mask
#  and plots the three cross-sections.
#
#  EXAMPLE (SCIENCE PAPER):
#    plot_horz_coverage.pl 0 5/5 0/4 0/1/6    # no mask
#    plot_horz_coverage.pl 1 5/5 0/4 0/0/6    # with mask
#    plot_horz_coverage.pl 0 21/21 0/4 0/1/3
#    plot_horz_coverage.pl 1 21/21 0/4 0/0/3
#    plot_horz_coverage.pl 0 41/41 0/4 0/1/3
#    plot_horz_coverage.pl 1 41/41 0/4 0/0/3
#
#==========================================================

if (@ARGV < 4) {die("Usage: plot_horz_coverage.pl xxx\n");}
($imask,$pinds,$cvals,$parms) = @ARGV;

$icolor = 1;
($pmin,$pmax) = split(/\//,$pinds);
($cmean,$cpert) = split(/\//,$cvals);
($itoplab,$isidelab,$iB) = split(/\//,$parms);

# KEY COMMAND: HOW DID YOU WEIGHT THE SUMMED KERNELS
$smodeltag = "coverage_sum_abs_nwin";
#$smodeltag = "coverage_sum_abs";

# detailed flags for publication figures
#$ilonlab = 0;      # display longitude tick labels at the top
#$itoplab = 1;      # label the model above the plot
#$isidelab = 1;      # label of the depth on the left of the plot
$ivs = 1;           # KEY COMMAND
$iinsetlab = 0;    # label the model inan inset
$iletter = 0;         # label for publication figures
$letter = "C";

$cgray = "-G220/220/220";

if($ivs==1) {$modlab = "vs"; $mtit = "K_Vs";} else {$modlab = "vb"; $mtit = "K_Vb";}

#---------------------------------------------------------
# USER INPUT

# directory containing the data files
$stirun = "02";
$dirdat = "INPUT/horz_${stirun}";

# KEY: file showing the cuts and the plotting range for the values
$fcuts = "${dirdat}/horz_${stirun}_xc_cuts_mod";
if (not -f $fcuts) {die("Check if fcuts $fcuts exist or not\n")}

# load values and check that the files exist
open(IN,$fcuts); @lines = <IN>; $nump = @lines; close(IN);
for ($p = 1; $p <= $nump; $p ++ ) {
   $stip = sprintf("%3.3i",$p);
   (undef,$zcut,$vs_norm,$vb_norm,$vpert) = split(" ",$lines[$p-1]);
   $vperts[$p-1] = $vpert;
   $vsnorms[$p-1] = $vs_norm;
   $vbnorms[$p-1] = $vb_norm;
   $zcuts[$p-1]  = -$zcut;
   $dfiles[$p-1] = "${dirdat}/horz_${stirun}_xc_${modlab}_m16_${smodeltag}_${stip}.dat";
   #if (not -f $dfile) {die("Check if dfile $dfile exist or not\n")}
}
print "\nNumber of cuts is $nump\n";
if($pmax > $nump) {$pmax = $nump;}
if($pmin < 1) {$pmin = 1;}

# subtitles
@titles = ("$mtit  coverage");
@labs = ("$mtit");

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
$xtx3 = 0.5;   $ytx3 = -0.1;
$xtx4 = $xtx3; $ytx4 = -0.3;
$xtx5 = -0.15; $ytx5 = 0.85;    # A, B, C, etc
$xtx6 = 0.0; $ytx6 = 0.05;    # m00, m16, ln(m16/m00)

#print "$dtitle\n"; die("TESTING");

#---------------------------------------------------------
# DATA FILES

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

 # $xmin = -121.2; $xmax = -114.8; $ymin = 32.3; $ymax = 36.7; $tick1 = 1; $tick2 = 0.5;
  $xmin = -122; $xmax = -114; $ymin = 32; $ymax = 37; 
  $xtick1 = 2; $xtick2 = 0.5; $ytick1 = 1; $ytick2 = 0.5;

  $cmax = 5000; $cmin = -$cmax;
  $y_title = 0.98; $y_title2 = 0.92; $x_title = 0.5; $x_title2 = 0.5;
  $iportrait = 0;
  $stag = "Southern California";
  $topo_labels  = "$dir0/socal_2005/socal_topo_labs.xyz";
  $fault_labels = "$dir0/socal_2005/socal_fault_labs.xyz";
  $name = "plot_horz_coverage";
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

$B0 = "-Ba${xtick1}f${xtick2}d:\" \":/a${ytick1}f${ytick2}d:\" \":";
#$B0 = "-Ba${xtick1}f${xtick2}d:.\" \":";
#$B0 = "-Ba${tick1}f${tick1}g${tick1}:.\" \":";   # with gridlines

$B = "$B0".$Bopts[0];

if($iportrait == 1){$orient = "-P"; $rotangle = 0}
else {$orient = " "; $rotangle = 90}

# color bar
$Dlen = $wid*0.4;
$Dx = $Dlen*0.5;
$Dy = -0.2;
$Dscale = "-D$Dx/$Dy/$Dlen/0.15h";

# plotting specifications
$fsize0 = "14";
$fsize1 = "12";
$fsize2 = "10";
$fsize3 = "6";
$fontno = "1";    # 1 or 4
$tick   = "6p";
$fpen   = "1p";
$tpen   = "1p";

# plotting specificiations
# NOTE: pen attribute (-W) does not work for psvelo
#$coast_info     = "-A5000 -Dl -W1.0p -S220/255/255 -C220/255/255 -G200";  # A : smallest feature plotted, in km^2; D : resolution

$coast_info     = "${coast_res} -W0.75p -Na/0.75p";
$fault_info     = "-M -W0.75p,255/0/0";
$fault_infoK    = "-M -W0.75p,0/0/0";

$city_info = "-Sc10p -W1.0p/0/0/0 -G255/255/0";

$textinfo       = "-G255 -S1.5p";
$textinfo2      = "-G0/0/0 -S2p,255/255/255";
$textinfo3      = "-G0/0/0 -S2p,255/255/0";

# BOOLEAN: WHICH FIGURES TO PLOT
$ixv = 0; $ipdf = 1;
#$ixv = 1; $ipdf = 0;

# plot title
$J_title = "-JM${wid}";
$R_title = "-R0/1/0/1";
$fsize_title = $fsize0;

#==================================================

$cshfile = "plot_horz_coverage.csh";
open(CSH,">$cshfile");
print CSH "gmtset PAPER_MEDIA letter BASEMAP_TYPE plain PLOT_DEGREE_FORMAT D TICK_LENGTH $tick LABEL_FONT_SIZE $fsize2 ANOT_FONT_SIZE $fsize2  HEADER_FONT $fontno ANOT_FONT $fontno LABEL_FONT $fontno HEADER_FONT_SIZE $fsize1 FRAME_PEN $fpen TICK_PEN $tpen \n";

@shifts = ($origin,"-X$dX","-X$dX");
#if($ilonlab==1) {@iBs = (2,14,2);} else {@iBs = (15,9,15);}

for ($p = $pmin; $p <= $pmax; $p ++ ) {

  $stip = sprintf("%3.3i",$p);
  $fname = "horz_${stirun}_xc_${modlab}_m16_${smodeltag}_${stip}_${imask}";
  $psfile = "$fname.ps"; $jpgfile = "$fname.jpg";

  $vsnorm = $vsnorms[$p-1];	# reference velocity for m00 and m16
  $vbnorm = $vbnorms[$p-1];	# reference velocity for m00 and m16
  $vpert = $vperts[$p-1];       # perturbation from reference velocity to plot
  $zdep = $zcuts[$p-1];		# depth of the cross-section (m)

   # zero-level for mask normalization
   $cfile = "${dirdat}/horz_${stirun}_xc_${modlab}_m16_${smodeltag}_${stip}_mask_value.dat";
   if (not -f $cfile) {die("Check if cfile $cfile exist or not\n")}
   open(IN,$cfile); $cnorm = <IN>; close(IN); chomp($cnorm);
   #$cnorms[$p-1] = $cnorm;
   #print "-- $cnorm --\n"; die("TESTING:");

  $dtitle = sprintf("Depth  =  %.1f km",$zdep/1000);

  # datafile -- COLUMNS: lon, lat, m00, m16
  $dfile = $dfiles[$p-1];
  if (not -f $dfile) {die("Check if dfile $dfile exist or not\n");}
  print "Data file is $dfile\n";

  # colorpoint file for m00 and m16 -- VARIES with depth -- in ln(m00/cnorm)
  $cptfile = "color.cpt";
#   if($cmean >= 0) {
#     $cmin = $cmean * exp(-$cpert);
#     $cmax = $cmean * exp($cpert);
#   } else {
#     $cmin = $cmean * exp($cpert);
#     $cmax = $cmean * exp(-$cpert);
#   }
  $cmin = -$cpert;
  $cmax = $cpert;
  $dc = ($cmax-$cmin)/${scale_color};
  $T = sprintf("-T%3.3e/%3.3e/%3.3e",$cmin,$cmax,$dc);
  #print CSH "makecpt -C$colorbar $T -D -Z > $cptfile\n";
  print CSH "makecpt -Chot $T -D > $cptfile\n";
  #print "\n colorbar is -- $T -- $cmin, $cmean, $cmax \n"; die("TESTING");

  $ctick = $cmax;
  $Bscale  = sprintf("-Ba%2.2ef%2.2e:\" \": -E10p",$ctick,$ctick/2);

  for ($k = 1; $k <= 1; $k ++ ) {

    $shift = $shifts[$k-1];
    #$iB = $iBs[$k-1];
    $B = "$B0".$Bopts[$iB];
    $title = $titles[$k-1];
    $kcol = 3 + $imask;

    if ($k==1) {
      print CSH "psbasemap $J $R $B -K -V $orient $origin > $psfile\n";	# START
    } else {
      print CSH "psbasemap $J $R $B -K -O -V $shift >> $psfile\n";
    }
    print CSH "psbasemap $J $R $B $cgray -K -O -V >> $psfile\n";
    if ($icolor==1) {
	#print CSH "awk '{print \$1,\$2,log(\$3/$cnorm)}' $dfile | nearneighbor -G$grdfile $R $interp\n";
        print CSH "awk '{print \$1,\$2,log(\$${kcol}/$cnorm)}' $dfile | xyz2grd -G$grdfile $R $interp\n";
	print CSH "grdimage $grdfile -C${cptfile} $J -Q -K -O -V >> $psfile\n";
        print CSH "psscale -C${cptfile} $Dscale $Bscale -K -O -V >> $psfile\n";

        #$slab1 = sprintf("ln ( %s / %.1e )",$mtit,$cnorm);
        $slab1 = sprintf("ln ( %s / 10\@+%2.2i\@+ )",$mtit,log($cnorm)/log(10));

        print CSH "pstext -N $R_title $J_title -K -O -V >>$psfile<<EOF\n$xtx3 $ytx3 $fsize2 0 $fontno LM $slab1\nEOF\n";
    }
    print CSH "pscoast $J $R $coast_info -K -O -V >> $psfile \n";
    if($ishelf==1) {print CSH "awk '{print \$2,\$1}' ${shelf_file} |psxy $J $R $Wshelf -K -O -P -V >> $psfile\n";}

    #--------------------------------

    #print CSH "psxy ${plate_file} $J $R $plate_info -K -V -O >> $psfile \n";
    print CSH "psxy ${fault_file} $J $R $fault_infoK -K -V -O >> $psfile \n";
    #print CSH "psxy ${kcf_file} $J $R $fault_infoK -K -V -O >> $psfile \n";
    #print CSH "psxy ${breck_file} $J $R $fault_infoK -K -V -O >> $psfile \n";
    #print CSH "psxy ${fault_file2} $J $R $fault_infoK -K -V -O >> $psfile \n";

    if ($imask==1) {
      # boundaries of simulation
      print CSH "psxy ${outer_boundary} $J $R -W2p,0/0/0 -K -O -V >>$psfile\n";  
      #print CSH "psxy ${inner_boundary} $J $R -W1.5p,0/0/0,-- -K -O -V >>$psfile\n";  

      $ibox = 0;
      if ($ibox==1) {
	$boxinfo = "-W1.5p,0/0/255"; # -A : suppress drawing line segments as great circle arcs

	# Lin model (2007) box
	$lin_boundary = "/net/denali/home2/carltape/gmt/tomography/lin_2007/lin_boundary_points.dat";
	print CSH "psxy ${lin_boundary} $J $R $boxinfo -K -O -V >>$psfile\n";  
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
  if($ixv==1) {print CSH "ghostview $psfile &\n";}
  if($ipdf==1) {print CSH "ps2pdf $psfile\n";}

}				# loop over p
  
#------------------------------------
# print CSH "pstext -N $R_title $J_title -O -V >>$psfile<<EOF\n $x_title $y_title 16 0 $fontno CM \nEOF\n"; # FINISH 

close (CSH);
system("csh -f $cshfile");

#   print "convert $psfile $jpgfile \n";
#   system("convert $psfile -rotate $rotangle $jpgfile");
#   if ($ipdf==1) {system("ps2pdf $psfile")}
#   if ($ixv==1) {system("ghostview $psfile &");}

system("gv $psfile &")

#==================================================
