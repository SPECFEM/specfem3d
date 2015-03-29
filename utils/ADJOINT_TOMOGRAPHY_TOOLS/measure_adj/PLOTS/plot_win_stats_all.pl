#!/usr/bin/perl -w
#
#-----------------------------------
# plot_win_stats_all.pl
#
# This script calls plot_win_stats.pl to plot histograms of measurements for each event.
# It also compiles all the measurements into one histogram, with the option of plotting here.
#
# The set of final histograms are copied to an output directory.
#
# EXAMPLE:
#    plot_win_stats_all.pl 6/30 m16 0
#    plot_win_stats_all.pl 6/30 m16 1 (plot histogram)
#
#-----------------------------------

if (@ARGV < 3) {die("Usage: plot_win_stats_all.pl Tmin/Tmax model iplot\n")}
($Ts,$smodel,$iplot) = @ARGV;

$pwd = $ENV{PWD};

($Tmin,$Tmax) = split("/",$Ts);
$sTmin = sprintf("T%3.3i",$Tmin);
$sTmax = sprintf("T%3.3i",$Tmax);
$Ttag = "${sTmin}_${sTmax}";
$sTper = sprintf("T = %i - %i s",$Tmin,$Tmax);
$ftag = "${Ttag}_${smodel}";

$mid = 16;         # index for directory with earthquake sources
$dir_source = "/home/carltape/results/SOURCES/socal_${mid}";
$dir_source_text = "${dir_source}/v${mid}_files_text";

#$file_eids = "${dir_source}/SOCAL_FINAL_CMT_v${mid}_eid";
#$file_eids = "/home/carltape/results/EID_LISTS/eids_simulation";
$file_eids = "/home/carltape/results/EID_LISTS/eids_tomo";
#$file_eids = "/home/carltape/results/EID_LISTS/eids_extra";
if (not -f $file_eids) {die("\n check if $file_eids exists\n")}
open(IN,$file_eids); @eids = <IN>; $nevent = @eids;

# directories
$odir    = "/home/carltape/results/MEASUREMENTS";
$dir_all = "$odir/$smodel/PDF_EID_${Ttag}";
$dir_run = "/home/carltape/socal/socal_3D/RUNS";
if (not -e $dir_all) {die("check if $dir_all exist or not\n")}
if (not -e $dir_run) {die("check if $dir_run exist or not\n")}

# name of file containing ALL source focal mechanisms in PSMECA format (see socal_tomo_sources.m)
$cmtall_psmeca = "${dir_source}/SOCAL_FINAL_CMT_v${mid}_psmeca_eid";
if (not -f $cmtall_psmeca) {die("\n check if $cmtall_psmeca exists\n")}

#$nevent = 5;

# make a directory to dump all the individual measurement files into
$dir_meas_collect = "meas_files_${ftag}";
`rm -rf ${dir_meas_collect} ; mkdir ${dir_meas_collect}`;

# remove plot-related files from local directory
`rm *dAsigma *jpg *ps *pdf`;

print "\n\n Running plot_win_stats_all.pl for $nevent possible events...\n";

if (0==1) {
  for ($i = 1; $i <= $nevent; $i = $i+1) {
    $eid = $eids[$i-1]; chomp($eid);
    print "-- $i, $eid --\n";
  }
  die("LISTING EIDS ONLY");
}

# loop over all events
$k = 0;
$imin = 1; $imax = $nevent;   # default
#$imin = 1; $imax = 20;
#$imin = 81; $imax = $imin;

for ($i = $imin; $i <= $imax; $i = $i+1) {

  $eid = $eids[$i-1]; chomp($eid);
  $ftag = "${eid}_${Ttag}_${smodel}";

  print "\n------------------------------";
  print "\n $i, $eid \n";

  # name of file containing measurements (see ma_sub.f90)
  $dir_meas = "${dir_run}/${eid}/${smodel}/MEASURE_${Ttag}";
  $meas_file = "${dir_meas}/${ftag}_window_chi";
  if (-f $meas_file) {

    # name of stations file
    $sta_file = "${dir_meas}/ADJOINT_${Ttag}/STATIONS_ADJOINT";

    # file containing text for plotting here
    $eid_text = "${dir_source_text}/${eid}_text";
    if (not -f $eid_text) {
      die("\n check if eid_text $eid_text exists\n");
    }

    # plot figure
    $name = sprintf("%3.3i_${ftag}_win_stats",$i);
    print "\n plot_win_stats.pl $Ts $eid $name $meas_file $sta_file $eid_text ${cmtall_psmeca}\n\n";
    #die("TESTING");
    `plot_win_stats.pl $Ts $eid $name $meas_file $sta_file $eid_text ${cmtall_psmeca}`;

    # check if the figure was made
    if (-f "${name}.pdf") {
      $k = $k+1;
      print "$i, $eid, ${name}.pdf was made\n";
    } else {
      print "$i, $eid, ${name}.pdf was NOT made\n";
    }

  } else {
    print "meas_file $meas_file does not exist\n";
  }
}

#----------------------------------------------

#$file1 = "${dir_meas_collect}/ALL_dT_cc";
#$file2 = "${dir_meas_collect}/ALL_dA_cc";
#$file3 = "${dir_meas_collect}/ALL_dT_cc_sigma";
#$file4 = "${dir_meas_collect}/ALL_dA_cc_sigma";
#`cat ${dir_meas_collect}/*dT_cc > $file1`;
#`cat ${dir_meas_collect}/*dA_cc > $file2`;
#`cat ${dir_meas_collect}/*dT_cc_sigma > $file3`;
#`cat ${dir_meas_collect}/*dA_cc_sigma > $file4`;

`mv *dAsigma ${dir_meas_collect}`;
$file0 = "${dir_meas_collect}/ALL_cc_dT_dA_dTsigma_dAsigma";
`cat ${dir_meas_collect}/*dAsigma > $file0`;

# copy measurement files to output directory
`cp -r ${dir_meas_collect} $dir_all`;

#================================================
# PLOT SUMMARY HISTOGRAMS OF ALL WINDOWS

if ($iplot==1) {

$fontno = 1;

$name = "ALL_${Ttag}_${smodel}_win_stats";
$psfile = "$name.ps";
$jpgfile = "$name.jpg";
$pdffile = "$name.pdf";

#=============================================
# write the C-shell script to file
$cshfile = "plot_win_stats_all.csh";
print "\nWriting to $cshfile ...\n";
open(CSH,">$cshfile");

# set GMT defaults
print CSH "gmtset BASEMAP_TYPE plain PAPER_MEDIA letter MEASURE_UNIT inch TICK_LENGTH 5p ANNOT_FONT_SIZE = 10p LABEL_FONT_SIZE = 12p PLOT_DEGREE_FORMAT D\n";

#-----------------------------------------
# plot histograms

# get the number of events in the list
($nwin,undef,undef,undef) = split(" ",`wc $file0`);
print "\n -- $nwin total windows -- \n";

if($nwin < 1) {die("\n No windows -- stop here.\n\n")}

#($min1,$max1) = split(" ",`minmax -C $file1`); $st1 = sprintf("min/max = %.2f / %.2f",$min1,$max1);
#($min2,$max2) = split(" ",`minmax -C $file2`); $st2 = sprintf("min/max = %.2f / %.2f",$min2,$max2);
#($min3,$max3) = split(" ",`minmax -C $file3`); $st3 = sprintf("min/max = %.2f / %.2f",$min3,$max3);
#($min4,$max4) = split(" ",`minmax -C $file4`); $st4 = sprintf("min/max = %.2f / %.2f",$min4,$max4);

($min1,$max1,$min2,$max2,$min3,$max3,$min4,$max4,$min5,$max5,$min6,$max6) = split(" ",`minmax -C $file0`);
$st1 = sprintf("min/max = %.2f / %.2f",$min1,$max1);
$st2 = sprintf("min/max = %.2f / %.2f",$min2,$max2);
$st3 = sprintf("min/max = %.2f / %.2f",$min3,$max3);
$st4 = sprintf("min/max = %.2f / %.2f",$min4,$max4);
$st5 = sprintf("min/max = %.2f / %.2f",$min5,$max5);
$st6 = sprintf("min/max = %.2f / %.2f",$min6,$max6);

print "\n $st1 \n $st2 \n $st3 \n $st4 \n $st5 \n $st6 \n";

# dimensions of histogram subplots
$widhist = 2.5; $heighthist = 2.0;

# KEY: commands for histograms axes
$Nbin = 20;
$max_DT = 8.0; $min_DT = -$max_DT;
$max_DA = 3.0; $min_DA = -$max_DA;
$max_sigma_DT = 2.5; $min_sigma_DT = 0;
$max_sigma_DA = $max_sigma_DT; $min_sigma_DA = $min_sigma_DT ;
$max_DTerr = $max_DT; $min_DTerr = -$max_DTerr;
$max_DAerr = $max_DA; $min_DAerr = -$max_DAerr;
if($Tmin == 2 || $Tmin == 3) {
  $max_DT = 4; $min_DT = -$max_DT;
  $max_DA = 1.5; $min_DA = -$max_DA;
  $max_sigma_DT = 1.5; $min_sigma_DT = 0;
  $max_sigma_DA = $max_sigma_DT; $min_sigma_DA = $min_sigma_DT ;
  $max_DTerr = 10; $min_DTerr = -$max_DTerr;
  $max_DAerr = 4; $min_DAerr = -$max_DAerr;
}
$DT_bin = ($max_DT - $min_DT)/$Nbin;
$DA_bin = ($max_DA - $min_DA)/$Nbin;
$DT_sigma_bin = ($max_sigma_DT - $min_sigma_DT)/$Nbin;
$DA_sigma_bin = ($max_sigma_DA - $min_sigma_DA)/$Nbin;
$DTerr_bin = ($max_DTerr - $min_DTerr)/$Nbin;
$DAerr_bin = ($max_DAerr - $min_DAerr)/$Nbin;
#print "\n $DT_bin $DA_bin $DT_sigma_bin $DA_sigma_bin $DTerr_bin $DAerr_bin \n"; die("testing");

# locations for plots
$X0 = 1.25; $Y0 = 8;
$dX = $widhist+1; $dY = $heighthist+1;
@ohistX = ($X0, $X0+$dX, $X0, $X0+$dX, $X0, $X0+$dX);
@ohistY = ($Y0, $Y0, $Y0-$dY, $Y0-$dY, $Y0-2*$dY, $Y0-2*$dY);
$origin_hist1 = "-Xa$ohistX[0] -Ya$ohistY[0]";
$origin_hist2 = "-Xa$ohistX[1] -Ya$ohistY[1]";
$origin_hist3 = "-Xa$ohistX[2] -Ya$ohistY[2]";
$origin_hist4 = "-Xa$ohistX[3] -Ya$ohistY[3]";
$origin_hist5 = "-Xa$ohistX[4] -Ya$ohistY[4]";
$origin_hist6 = "-Xa$ohistX[5] -Ya$ohistY[5]";

  # plot histograms of DA and DT
  $max_count = 40; $count_tick = 5;

  $Rhist1 = "-R$min_DT/$max_DT/0/$max_count";
  $Rhist2 = "-R$min_DA/$max_DA/0/$max_count";
  $Rhist3 = "-R$min_sigma_DA/$max_sigma_DA/0/$max_count";
  $Rhist4 = "-R$min_sigma_DA/$max_sigma_DA/0/$max_count";
  $Rhist5 = "-R$min_DTerr/$max_DTerr/0/$max_count";
  $Rhist6 = "-R$min_DAerr/$max_DAerr/0/$max_count";

  $Jhist = "-JX${widhist}i/${heighthist}";

  $sty = "Percent  (N = $nwin)";
  $Bhist1 = "-Ba4f1:\"dT (s)  ($st1)\":/a${count_tick}:\"$sty\":WeSn";
  $Bhist2 = "-Ba1f0.25:\"dlnA  ($st2)\":/a${count_tick}:\"$sty\":WeSn";
  $Bhist3 = "-Ba0.5f0.1:\"sigma-dT (s)  ($st3)\":/a${count_tick}:\"$sty\":WeSn";
  $Bhist4 = "-Ba0.5f0.1:\"sigma-dlnA  ($st4)\":/a${count_tick}:\"$sty\":WeSn";
  $Bhist5 = "-Ba4f1:\"dT / sigma-dT  ($st5)\":/a${count_tick}:\"$sty\":WeSn";
  $Bhist6 = "-Ba1f0.25:\"dlnA / sigma-dlnA  ($st6)\":/a${count_tick}:\"$sty\":WeSn";

  $hist_info = "-W5 -G0/255/255 -L0.5p,0/0/0 -Z1";

  print CSH "pshistogram $file0 $hist_info $Jhist $Bhist1 $Rhist1 -T0 -W${DT_bin} -K -V $origin_hist1 > $psfile\n";   # START
  print CSH "pshistogram $file0 $hist_info $Jhist $Bhist2 $Rhist2 -T1 -W${DA_bin} -K -O -V $origin_hist2 >> $psfile\n";
  print CSH "pshistogram $file0 $hist_info $Jhist $Bhist3 $Rhist3 -T2 -W${DT_sigma_bin} -K -O -V $origin_hist3 >> $psfile\n";
  print CSH "pshistogram $file0 $hist_info $Jhist $Bhist4 $Rhist4 -T3 -W${DA_sigma_bin} -K -O -V $origin_hist4 >> $psfile\n";
  print CSH "pshistogram $file0 $hist_info $Jhist $Bhist5 $Rhist5 -T4 -W${DTerr_bin} -K -O -V $origin_hist5 >> $psfile\n";
  print CSH "pshistogram $file0 $hist_info $Jhist $Bhist6 $Rhist6 -T5 -W${DAerr_bin} -K -O -V $origin_hist6 >> $psfile\n";
#  print CSH "awk '{print \$1/\$3}' $file0 | pshistogram $hist_info $Jhist $Bhist5 $Rhist5 -W${DTerr_bin} -K -O -V $origin_hist5 >> $psfile\n";
#  print CSH "awk '{print \$2/\$4}' $file0 | pshistogram $hist_info $Jhist $Bhist6 $Rhist6 -W${DAerr_bin} -K -O -V $origin_hist6 >> $psfile\n";

  # plot lines
  print CSH "psxy -W1.5p,0/0/0 $Jhist $Rhist1 -K -O -V $origin_hist1 >>$psfile<<EOF\n0 0\n0 $max_count\nEOF\n";
  print CSH "psxy -W1.5p,0/0/0 $Jhist $Rhist2 -K -O -V $origin_hist2 >>$psfile<<EOF\n0 0\n0 $max_count\nEOF\n";
  print CSH "psxy -W1.5p,0/0/0 $Jhist $Rhist5 -K -O -V $origin_hist5 >>$psfile<<EOF\n0 0\n0 $max_count\nEOF\n";
  print CSH "psxy -W1.5p,0/0/0 $Jhist $Rhist6 -K -O -V $origin_hist6 >>$psfile<<EOF\n0 0\n0 $max_count\nEOF\n";

# show the number of windows
$x1 = -0.5; $y1 = $heighthist + 0.5;
print CSH "pstext -R0/1/0/1 -JX1i/1i -N -K -O -V $origin_hist1 >> $psfile <<EOF\n $x1 $y1 14 0 $fontno LM $nwin windows picked for period range $sTper for $k out of $nevent events\nEOF\n";

#-----------------------------------------
# FINISH
print CSH "psxy -R0/1/0/1 -JX1i/1i -X15 -O -V >>$psfile<<EOF\n 0 0\nEOF\n";

print CSH "convert $psfile $jpgfile\n";
print CSH "ps2pdf $psfile $pdffile\n";
#print CSH "rm $psfile\n";
print CSH "ghostview $psfile &\n";

#-----------------------------------------
close(CSH);
system("csh -f $cshfile");

}
#================================================

# concatenate pdf files into one
$ofile = "ALL_${smodel}_${Ttag}_win_stats_full.pdf";
`/home/carltape/bin/pdcat -r $pdffile \[0-9\]*${Ttag}*.pdf $ofile`;

# copy composite file to output directory
`cp $ofile $odir/$ofile`;

print "\n\n";

#================================================
