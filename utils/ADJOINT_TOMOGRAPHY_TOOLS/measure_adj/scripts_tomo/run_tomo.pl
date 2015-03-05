#!/usr/bin/perl -w

#==========================================================
#
#  run_tomo.pl
#  Carl Tape
#  20-Oct-2009
#
#  This is the master program for running the measurement code on all the
#  events for which the windowing code has been run.
#
#  NOTE: The command run_measure_adj.pl below contains additional input parameters.
#
#  EXAMPLES (19-Oct-2009):
#    run_tomo.pl m16 0.011 1.0 5 BH 6/30 0/1/1/0 -5.0/5.0/-1.5/1.5/0.7 1/1.0/0.5 3/0.02/2.5/2.0/2.5/3.5/1.5  # CC: no adjoint source
#    run_tomo.pl m16 0.011 1.0 1 BH 6/30 0/1/1/1 -5.0/5.0/-1.5/1.5/0.7 1/1.0/0.5 3/0.02/2.5/2.0/2.5/3.5/1.5
#    run_tomo.pl m16 0.011 1.0 2 BH 6/30 0/1/1/1 -5.0/5.0/-1.5/1.5/0.7 1/1.0/0.5 3/0.02/2.5/2.0/2.5/3.5/1.5
#    run_tomo.pl m16 0.011 1.0 3 BH 6/30 0/1/1/1 -5.0/5.0/-1.5/1.5/0.7 1/1.0/0.5 3/0.02/2.5/2.0/2.5/3.5/1.5
#    run_tomo.pl m16 0.011 1.0 4 BH 6/30 0/1/1/1 -5.0/5.0/-1.5/1.5/0.7 1/1.0/0.5 3/0.02/2.5/2.0/2.5/3.5/1.5
#    run_tomo.pl m16 0.011 1.0 5 BH 6/30 0/1/1/1 -5.0/5.0/-1.5/1.5/0.7 1/1.0/0.5 3/0.02/2.5/2.0/2.5/3.5/1.5  # CC: adjoint source
#    run_tomo.pl m16 0.011 1.0 6 BH 6/30 0/1/1/1 -5.0/5.0/-1.5/1.5/0.7 1/1.0/0.5 3/0.02/2.5/2.0/2.5/3.5/1.5
#    run_tomo.pl m16 0.011 1.0 7 BH 6/30 0/1/1/1 -5.0/5.0/-1.5/1.5/0.7 1/1.0/0.5 1/0.02/2.5/2.0/2.5/3.5/1.5
#    run_tomo.pl m16 0.011 1.0 8 BH 6/30 0/1/1/1 -5.0/5.0/-1.5/1.5/0.7 1/1.0/0.5 1/0.02/2.5/2.0/2.5/3.5/1.5
#    run_tomo.pl m16 0.011 1.0 7 BH 3/30 0/1/1/1 -3.0/3.0/-1.5/1.5/0.7 1/1.0/0.5 1/0.02/2.5/0.25/0.5/3.5/1.5
#    run_tomo.pl m16 0.011 1.0 7 BH 2/30 0/1/1/1 -2.0/2.0/-1.5/1.5/0.7 1/1.0/0.5 1/0.02/2.5/0.25/0.5/3.5/1.5
#
#  socal m00 - m04, 0.2/0.2
#  socal m05 - m08, 0.4/0.4
#  socal m09, MAX: 4.0 pm 0.6 (2-30s), 4.0 pm 0.8 (3-30s), 5.0 pm 0.8 (6-30s)
#  socal m10, MAX: 2.0 pm 0.6 (2-30s), 3.0 pm 0.8 (3-30s), 5.0 pm 0.8 (6-30s)
#  socal m12, MAX: 2.0 pm 0.6 (2-30s), 3.0 pm 0.8 (3-30s), 5.0 pm 1.0 (6-30s)
#  socal m13 --> m16, MAX: 2.0 pm 1.0 (2-30s), 3.0 pm 1.0 (3-30s), 5.0 pm 1.0 (6-30s)
#
#==========================================================

if (@ARGV < 10) {die("Usage: run_tomo.pl xxx\n")}
($smodel,$dt,$lcut_frac,$imeas,$chan,$Ts,$iparbools,$par1,$par2,$par3) = @ARGV;
#($smodel,$ibp,$Ts,$dt,$iparbools,$par1,$par2,$par3,$lcut_frac) = @ARGV;

$pwd = $ENV{PWD};

($Tmin,$Tmax) = split("/",$Ts);
$sTmin = sprintf("T%3.3i",$Tmin);
$sTmax = sprintf("T%3.3i",$Tmax);
$Ttag = "${sTmin}_${sTmax}";

# cut times for plotting purposes only
#($lcut_min,$lcut_max) = split("/",$lcut);

#================================================================
# USER INPUT

$dir0 = "/home/carltape";
#$dir0 = "/media/raid/carltape/SOCAL_ADJOINT";
$dir_out = "$dir0/RUNS";
$dir_pdf_eid = "$dir0/results/MEASUREMENTS/${smodel}/PDF_EID_${Ttag}";
$dir_pdf_all = "$dir0/results/MEASUREMENTS/${smodel}/PDF_ALL_${Ttag}";
if (not -e $dir_out) {die("check if dir_out $dir_out exist or not\n")}
if (not -e $dir_pdf_eid) {die("check if dir_pdf_eid $dir_pdf_eid exist or not\n")}
if (not -e $dir_pdf_all) {die("check if dir_pdf_all $dir_pdf_all exist or not\n")}

# file containing event IDs, half-durations for source, and number of time-steps for simulations
$file_eid = "$dir0/results/SOURCES/socal_16/EIDs_simulation_duration";
if (not -f ${file_eid}) {die("check if ${file_eid} exist or not\n")}
open(IN,"${file_eid}"); @elines = <IN>; close(IN);
$nevent = @elines;
print "\n $nevent events in the list\n";

# full stations file (copy into PLOTS directory)
$file_stations = "/home/carltape/gmt/stations/seismic/Matlab_output/STATIONS_CALIFORNIA_TOMO_INNER_specfem";
if (not -f ${file_stations}) {die("check if file_stations ${file_stations} exist or not\n")}
`cp $file_stations PLOTS/STATIONS_TOMO`;

# sorted stations by distance (1) or azimuth (2)
$isort = 2;
if($isort == 1) {$sortlab = "dist"} else {$sortlab = "az"}

#================================================================

for ($ievent = 1; $ievent <= $nevent; $ievent++) {
   ($eindex,$eid,$Ldur_round,$nstep,$hdur,$Ldur,$mw) = split(" ",$elines[$ievent-1]);
   print "\n $ievent : $eid";
}
#die("testing");

$imin = 1; $imax = $nevent;
#$imin = 1; $imax = 265;
#$imin = 261; $imax = $imin;    # 183 Yorba Linda

# EXAMPLE lines from EIDs_simulation_duration
#      2   9967025     200.000     18200   0.29     4.084
#      3   9967137     200.000     18200   0.33     4.211
#      4   9967249     200.000     18200   0.32     4.191
#      5   9967541     200.000     18200   0.20     3.790

#=============================================
# write the C-shell script to file
$cshfile = "run_tomo.csh";
print "\nWriting to $cshfile ...\n";
open(CSH,">$cshfile");

for ($ievent = $imin; $ievent <= $imax; $ievent++) {

   # read in a row of numbers
   ($eindex,$eid,$Ldur_round,$nstep,$hdur,$mw) = split(" ",$elines[$ievent-1]);

   # get sorted stations file for ordering the plots
   $file_stations_sort = "$dir0/results/SOURCES/EID_STATION_LISTS/STATIONS_by_${sortlab}_from_${eid}";
   if (not -f ${file_stations_sort}) {die("check if file_stations_sort ${file_stations_sort} exist or not\n")}

   print "$ievent -- $imin to $imax : $eid\n";
   print CSH "echo $ievent -- $imin to $imax : $eid\n";

   # duration of records for plotting
   $lcut_min = -10;
   $lcut_max = $lcut_frac*$Ldur_round;
   $lcut = "${lcut_min}/${lcut_max}";

   # time vector for adjoint sources
   # KEY: The factor -1.5 is because this factor is used in SPECFEM3D.
   #      When the kernel simulation is run, it will pick the start time
   #      to be -1.5 hdur, where hdur comes from the CMTSOLUTION.
   $tstart = -1.5*$hdur;
   $tvec = "$tstart/$dt/$nstep";

   $ftag = "${eid}_${Ttag}_${smodel}";

   # check if window file exists
   $dir_eid = "${dir_out}/$eid/$smodel";
   $dir1 = "${dir_eid}/WINDOW_${Ttag}";
   $dir_win = "${dir_eid}/WINDOW_${Ttag}";
   $dir_meas = "${dir_eid}/MEASURE_${Ttag}";
   $dir_adj  = "${dir_meas}/ADJOINT_${Ttag}";

   $window_file = "${dir_win}/MEASUREMENT_WINDOWS_${ftag}";

   #if (-e $dir_adj) {
   if (-e $dir_meas) {
     #print "--> DO NOT RUN: adjoint source directory already exists ($dir_adj)\n";
     print "--> DO NOT RUN: measurement directory already exists ($dir_meas)\n";

   } else {

     if (not -f ${window_file}) {
     print "--> DO NOT RUN: no window file exists\n${window_file}\n";

     } else {

       # copy sorted STATIONS file
       print CSH "cp ${file_stations_sort} PLOTS/STATIONS_sort\n";
       open(IN,${file_stations_sort}); @stalines = <IN>; $nrec = @stalines; close(IN);
       print "Total number of stations sorted by $sortlab is $nrec\n";

       print "--> delete individual PDF files in ${dir_pdf_all}\n";
       print CSH "rm -rf ${dir_pdf_all}/*${eid}*\n";
       print CSH "sleep 3s\n";

       # this also compiles the measurement code, which is not needed for each new event
       print "--> run the measurement code and make adjoint sources\n";
       print CSH "prepare_measure_adj.pl $smodel 0 1 0 $Ts $eid\n";

       #------------------------------------

       $itest = 0;
       $iplot = 1;
       $iboth = 0;

       print CSH "run_measure_adj.pl $smodel $lcut $itest $iplot $iboth $tvec $imeas $chan $Ts $iparbools $par1 $par2 $par3\n";

       #print CSH "run_mt_cc_plot.pl $smodel $lcut $itest $tvec $imeas $chan $Ts $iparbools $par1 $par2 $par3\n";

       #------------------------------------

       # copy run files into RUN directory
       $infile1 = "window_chi";           $outfile1 = "${dir_meas}/${ftag}_${infile1}";
       $infile2 = "run_file";             $outfile2 = "${dir_meas}/${ftag}_${infile2}";
       $infile3 = "MEASUREMENT.PAR";      $outfile3 = "${dir_meas}/${ftag}_${infile3}";
       $infile4 = "MEASUREMENT.WINDOWS";  $outfile4 = "${dir_meas}/${ftag}_${infile4}";
                                          $outfile5 = "${dir_meas}/${ftag}_${infile4}_orig";
       print CSH "rm -rf ${dir_meas}; mkdir ${dir_meas}\n";
       print CSH "cp $infile1 $outfile1\n";
       print CSH "cp $infile2 $outfile2\n";
       print CSH "cp $infile3 $outfile3\n";
       print CSH "cp $infile4 $outfile4\n";

       # make a single PDF, with pages sorted by azimuth or distance
       print CSH "cd PLOTS\n";
       print CSH "make_pdf_by_event.pl\n";
       $infile5 = "${ftag}_ALL.pdf";
       $outfile5 = "${dir_meas}/${infile5}";
       print CSH "cp $infile5 $outfile5\n";        # composite PDF into RUN
       if (0==1) {
   print CSH "mv $infile5 ${dir_pdf_eid}\n"; # composite PDF into single directory
   print CSH "mv *.pdf ${dir_pdf_all}\n";    # individual PDFs into single directory
       }
       print CSH "cd $pwd\n";

       # move adjoint sources into RUN directory
       print CSH "rm -rf ${dir_adj}\nsleep 5s\n";
       print CSH "mv ADJOINT_SOURCES ${dir_adj}; mkdir ADJOINT_SOURCES\n";

       print CSH "echo done with Event $eid -- $ievent out of $imax/$nevent\n";
     }
   }

}

#-----------------------------------------
close(CSH);
#system("csh -f $cshfile");

#=================================================================
