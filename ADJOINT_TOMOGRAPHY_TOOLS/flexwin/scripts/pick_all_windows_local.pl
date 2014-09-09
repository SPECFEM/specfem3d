#!/usr/bin/perl -w

#==========================================================
#
#  pick_all_windows_local.pl
#  Carl Tape
#  20-Ocit-2009
#
#  This is the "local version" of pick_all_windows_local.pl
#  This is run from the directory containing MEASURE, DATA, SYN,
#  and assumes you have ALREADY run prepare_seis and prepare_input.
#
#  INPUT:
#    input         input file for windowing code
#    dir           directory containing windowing code output files
#    Tmin/Tmax     bandpass periods for filtering data and synthetics
#    idebug        boolean: debugging information
#    iplot         boolean: make composite plots for all windows
#    imeas         boolean: make *mt_input files for measurement code
#    ibody         boolean: turn body-wave mode ON
#
#  CALLS:
#     plot_windows_all.pl
#     prepare_meas_all.pl
#
#  EXAMPLES:
#     /data1/cig/seismo/3D/ADJOINT_TOMO/flexwin_work/scripts/pick_all_windows_local.pl
#
#     pick_all_windows_local.pl 2/30 1/1/1/1 input MEASURE
#     pick_all_windows_local.pl 6/30 1/1/1/0 input MEASURE
#
#     pick_all_windows_local.pl 6/30 1/1/1/0 input_test MEASURE_TEST    # plots AND window file
#     pick_all_windows_local.pl 6/30 0/0/1/0 input_test MEASURE_TEST    # window file only
#     pick_all_windows_local.pl 6/30 1/1/0/0 input_test MEASURE_TEST    # plots only
#     pick_all_windows_local.pl 2/30 1/1/0/1 input_test MEASURE_TEST    # plots only
#
#     pick_all_windows_local.pl 3/10 1/1/1/0 input_test MEASURE_TEST    # plots AND window file
#
#==========================================================

if (@ARGV < 4) {die("Usage: pick_all_windows_local.pl  Tmin/Tmax idebug/iplot/imeas/ibody input_file MEASURE_dir\n");}
($Ts,$ibools,$input,$dir) = @ARGV;

$pwd = $ENV{PWD};

#-------------------------------------
# USER INPUT

# source code directory
$dir_win_code = "/home/carltape/ADJOINT_TOMO/flexwin_work";

#-------------------------------------

$dir_scripts = "${dir_win_code}/scripts";
if (not -e ${dir_win_code}) {die("check if ${dir_win_code} exist or not\n")}
if (not -e ${dir_scripts}) {die("check if ${dir_scripts} exist or not\n")}

# split input entries
($Tmin,$Tmax)  = split("/",$Ts);
($imin,$imax)  = split("/",$inds);

# split boolean entries
($idebug,$iplot,$imeas,$ibody)  = split("/",$ibools);
if($idebug==1) {$bdebug = ".true."} else {$bdebug = ".false."}
if($iplot==1) {$bplot = ".true."} else {$bplot = ".false."}
if($imeas==1) {$bmeas = ".true."} else {$bmeas = ".false."}
if($ibody==1) {$bbody = ".true."} else {$bbody = ".false."}

# string for PDF
$sTmin = sprintf("%3.3i",$Tmin);
$sTmax = sprintf("%3.3i",$Tmax);
$sTminD = sprintf("%.2f",$Tmin);
$sTmaxD = sprintf("%.2f",$Tmax);

# user parameters for windowing code
$par_file = "${dir_win_code}/PAR_FILE";

# name of executable windowing code (in $dir_win_code)
$win_execute = "flexwin";

# make command for windowing code
$make = "make -f make_gfortran";

# location of plot_windows_all.pl
$plot_windows_perl = "${dir_scripts}/plot_windows_all.pl";
if (not -f ${plot_windows_perl}) {die("check if ${plot_windows_perl} exist or not\n")}

# location of prepare_meas_all.pl
$prepare_meas_perl = "${dir_scripts}/prepare_meas_all.pl";
if (not -f ${prepare_meas_perl}) {die("check if ${prepare_meas_perl} exist or not\n")}

# script to summarize window picks
$idataset = 1;
if ($idataset == 1) {
   $extract_script = "${dir_scripts}/extract_event_windowing_stats_carl.sh";
} elsif ($idataset == 2) {
   $extract_script = "${dir_scripts}/extract_event_windowing_stats_min.sh";
}
if (not -f ${extract_script}) {die("check if ${extract_script} exist or not\n")}

# sorted receiver files
$reclist0      = "${dir}/rectext";
$reclist_dist  = "${reclist0}_dist";
$reclist_picks = "${reclist0}_picks";

#===========================================================================
# write the C-shell script to file
$cshfile = "pick_all_windows_local.csh";
print "\nWriting to $cshfile ...\n";
open(CSH,">$cshfile");

# for the first event ($imin), change the user parameters file and compile
print CSH "cd ${dir_win_code}\n";
print CSH "sed '/^WIN_MIN_PERIOD                  =/s/^.*\$/WIN_MIN_PERIOD                  = $sTminD/' ${par_file} > temp1\n";
print CSH "sed '/^WIN_MAX_PERIOD                  =/s/^.*\$/WIN_MAX_PERIOD                  = $sTmaxD/' temp1 > temp2\n";
print CSH "sed '/^DEBUG                           =/s/^.*\$/DEBUG                           = $bdebug/' temp2 > temp3\n";
print CSH "sed '/^MAKE_SEISMO_PLOTS               =/s/^.*\$/MAKE_SEISMO_PLOTS               = $bplot/' temp3 > temp4\n";
print CSH "sed '/^MAKE_WINDOW_FILES               =/s/^.*\$/MAKE_WINDOW_FILES               = $bmeas/' temp4 > temp5\n";
print CSH "sed '/^BODY_WAVE_ONLY                  =/s/^.*\$/BODY_WAVE_ONLY                  = $bbody/' temp5 > temp6\n";
print CSH "\\mv temp6 ${par_file}\n \\rm temp1 temp2 temp3 temp4 temp5\n";

# compile the windowing code
print CSH "$make clean\n $make ${win_execute}\n";

# remove contents in target directory
# NOTE: This copies in the PAR_FILE from the main directory -- it is not the local version!
print CSH "cd $pwd\n";
print CSH "\\cp ${par_file} .\n";
print CSH "\\rm -rf ${dir}\n mkdir ${dir}\n";

# output file tag
$eout = "test_T${sTmin}";
 
# run the windowing code
print CSH "${dir_win_code}/${win_execute} < $input\n";

if ($iplot == 1) {
  # run extract_event_windowing_stats_carl.sh to get rectext and windows summary figure
  print CSH "${extract_script} ${dir}\n";

  # sort rectext into rectext_dist and rectext_picks
  print CSH "sort -g -k 2 $reclist0 > ${reclist_dist}\n";

  # generate composite PDF file using plot_windows_all.pl
  # make sure Carl's settings are on in plot_seismos_gmt.sh
  $pdffile = "$dir/${eout}_all.pdf";
  print CSH "${plot_windows_perl} ${dir_win_code} ${dir} ${reclist_dist} $pdffile\n";
}

if ($imeas == 1) {
  # generate MEASUREMENT.WINDOWS file from *mt_input files
  $ofile = "MEASUREMENT_WINDOWS";
  print CSH "${prepare_meas_perl} ${dir} $ofile\n";
}

#-----------------------------------------
close(CSH);
#system("csh -f $cshfile");

print "\n";

#=================================================================
