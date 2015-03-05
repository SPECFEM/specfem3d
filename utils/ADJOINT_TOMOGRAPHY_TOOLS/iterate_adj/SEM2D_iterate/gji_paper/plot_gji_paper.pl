#!/usr/bin/perl -w
#==========================================================
#
#  plot_gji_paper.pl
#  Carl Tape
#  10-Jan-2006
#
#  Figure 1: plot of misfit vs iteration (irun0 = 100, 20)
#
#==========================================================

$cshfile = "plot_gji_paper.csh";

# boolean commands for plotting
$ixv    = 1;
$ieps   = 0;

# plotting specifications
$fsize0 = "24";
$fsize1 = "18";
$fsize2 = "16";
$fsize3 = "12";
$fontno = "4";    # 1 or 4
$tick   = "0.2c";
$fpen   = "1.5p";
$tpen   = "1.0p";

# which borders to plot the lat-lon
# 1 four sides, 4 single sides, 6 two sides, 4 three sides, 1 zero sides
@Bopts = ("WESN","Wesn","wesN","wEsn","weSn","WesN","wEsN","wESn","WeSn","WEsn","weSN","wESN","WESn","WeSN","WEsN","wesn");

#=================================================================
# files and parameters related to polynomial plots and misfit-vs-iteration plots

  $pmark      = "-Sc10p";
  $p_info_k   = "-N $pmark -W1.0p -G0";
  $p_info_w   = "-N $pmark -W1.0p -G255";
  $p_info_r   = "-N $pmark -W1.0p -G255/0/0";
  $p_info_b   = "-N $pmark -W1.0p -G0/0/255";
  $p_info_g   = "-N $pmark -W1.0p -G0/0/255";

  $c_info_ks  = "-W1.0p";
  $c_info_ks2 = "-W2.0p";
  $c_info_rd  = "-W1.0/255/0/0tap";
  $c_info_bd  = "-W1.0/0/0/255tap";
  $c_info_kd  = "-W0.5/0/0/0tap";

  # input files to read
  $chi_curve = "chi_cubic_quad.dat";
  $axes_file = "chi_cubic_quad_axes.dat";
  if (not -f $axes_file)   { die("Check if $axes_file exist or not\n") }
  if (not -f $chi_curve)   { die("Check if $chi_curve exist or not\n") }

  # get bounds
  open(IN,"$axes_file"); @ax_lims = <IN>;
  ($a1x1,$a1x2,$a1y1,$a1y2) = split(" ",$ax_lims[0]);

  # chi-vs-iteration
  $xming = $a1x1;
  $xmaxg = $a1x2;
  $yming = $a1y1;
  $ymaxg = $a1y2;
  $R_chi = "-R$xming/$xmaxg/$yming/$ymaxg";

  # scale for chi-vs-m plots
  $B3a    = "-B2:\" model number (iteration) \":/a1f2p:\" \@~\143\@~ ( m )   ( $utype ) \":";
  $B3b    = "-B2:\" model number (iteration) \":/a1f2p:\" \":";
  $B4     = "-B2:\" model number (iteration) \":/20:\" \":";

#===========================================================================
#===========================================================================

  $origin = "-X1.25 -Y3";                   # origin
  $Jwid = 6;                               # width of plot
  $J_chi = "-JX${Jwid}i/${Jwid}il";        # note log scale on y-axis

  # plot title
  $J_title = "-JX${Jwid}";
  $R_title = "-R0/1/0/1";
  $x_title = 0.5;
  $z_title = 1.08;
  $z_title = 10;
  $fsize_title = $fsize1;

  # file names
  $name = "plot_gji_paper_01";
  $psfile = "${name}.ps"; $jpgfile = "${name}.jpg"; $epsfile = "${name}.eps";
  $title = " test";

  #===============================================
  print "\nWriting CSH file...\n";

  open(CSH,">$cshfile");
  print CSH "gmtset BASEMAP_TYPE plain TICK_LENGTH $tick LABEL_FONT_SIZE $fsize1 ANOT_FONT_SIZE $fsize2 PLOT_DEGREE_FORMAT D HEADER_FONT $fontno ANOT_FONT $fontno LABEL_FONT $fontno HEADER_FONT_SIZE $fsize1 FRAME_PEN $fpen TICK_PEN $tpen\n";
  #===============================================

  # misfit vs iteration
  $shift = $shift1;
  $Bg = "-Ba2f1:\" k, model number \":/a1f2p:\"  Misfit,  \@~\143\@~ (m\@+k\@+)  (s\@+2\@+)  \":".$Bopts[8];
  $title = " Recovered models for classical tomography and CG methods";

  $lab1 = "\@~\143\@~ (m\@+0\@+)";

  print CSH "psbasemap $Bg $R_chi $J_chi -K -V -P $origin > $psfile\n";  # START
  print CSH "awk '{print \$1,\$4}' $chi_curve | psxy -W1p $J_chi $R_chi -K -O -V -P >> $psfile\n";      # ray line
  print CSH "awk '{print \$1,\$5}' $chi_curve | psxy -W1tap $J_chi $R_chi -K -O -V -P >> $psfile\n";    # kernel line
  print CSH "awk '{print \$1,\$3}' $chi_curve | psxy $p_info_w $J_chi $R_chi -K -O -V -P >> $psfile\n"; # cubic points
  print CSH "awk '{print \$1,\$2}' $chi_curve | psxy $p_info_r $J_chi $R_chi -K -O -V -P >> $psfile\n"; # quad points

  # text labels
  print CSH "pstext -N $J_title $R_title -K -O -V -P >>$psfile<<EOF\n 0.08 0.38 14 0 $fontno LM Hessian - Ray \nEOF\n";
  print CSH "pstext -N $J_title $R_title -K -O -V -P >>$psfile<<EOF\n 0.08 0.31 14 0 $fontno LM Hessian - Kernel \nEOF\n";
  print CSH "pstext -N $J_title $R_title -K -O -V -P >>$psfile<<EOF\n 0.08 0.82 14 0 $fontno LM $lab1 \nEOF\n";

  #----------------------------------------------------------------------------
  # legend
  $origin_scale = "-Xa2.8 -Ya4.8";            # position of origin (inches)
  $bxmin = 0; $bymin = 0;                     # location of lower left box corner
  $bwid = 0.47;                               # width of box
  $nrow = 2;                                  # KEY: number of rows in legend
  $dy = 0.05;                                 # spacing for each row
  $dx1 = 0.04;                                # spacing between left box edge and label
  $dx2 = 0.03;                                # spacing between label and text
  $dy1 = $dx1;                                # spacing between bottom box edge and label
  $bxmax = $bxmin + $bwid;                    # position of right edge
  $bymax = $bymin + 2*$dy1 + ($nrow-1)*$dy;   # position of top edge
  $lx = $bxmin + $dx1;                        # x position of label
  $tx = $lx + $dx2;                           # x position of text
  $ry1 = $bymin + $dy1;                       # y position of bottom row
  $ry2 = $bymin + $dy1 + $dy;
  $ry3 = $bymin + $dy1 + 2*$dy;
  $ry4 = $bymin + $dy1 + 3*$dy;

  # legend labels (bottom row to top)
  @legendlabels = ("using a cubic polynomial","using a quadratic polynomial");

  print CSH "psxy -W1.0p $J_title $R_title -K -O -V -P $origin_scale >>$psfile<<EOF\n$bxmin $bymin\n$bxmax $bymin\n$bxmax $bymax\n$bxmin $bymax\n$bxmin $bymin\nEOF\n";
  print CSH "psxy $p_info_r $J_title $R_title -K -O -V -P $origin_scale >>$psfile<<EOF\n $lx $ry2 \nEOF\n";
  print CSH "psxy $p_info_w $J_title $R_title -K -O -V -P $origin_scale >>$psfile<<EOF\n $lx $ry1 \nEOF\n";
  print CSH "pstext -N $J_title $R_title -K -O -V -P $origin_scale >>$psfile<<EOF\n $tx $ry2 14 0 $fontno LM $legendlabels[1] \nEOF\n";
  print CSH "pstext -N $J_title $R_title -K -O -V -P $origin_scale >>$psfile<<EOF\n $tx $ry1 14 0 $fontno LM $legendlabels[0] \nEOF\n";
  #----------------------------------------------------------------------------

  print CSH "pstext -N $J_title $R_title -O -V -P >>$psfile<<EOF\n $x_title $z_title $fsize_title 0 $fontno CM $title \nEOF\n";  # FINISH

  print CSH "convert $psfile $jpgfile\n";
  if($ieps==1){print CSH "convert $psfile $epsfile\n"}
  close (CSH);
  system("csh -f $cshfile");
  system("xv $jpgfile &");

#=================================================================
