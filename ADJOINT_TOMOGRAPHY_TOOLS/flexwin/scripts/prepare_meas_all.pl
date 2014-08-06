#!/usr/bin/perl -w

#==========================================================
#
#  prepare_meas_all.pl
#  Carl Tape
#  01-Nov-2007
#
#  This file creates a MEASUREMENT.WINDOWS file that can be used in the
#  measurement code mt_measure_adj.f90.
#
#  NOTE: these lines should simply be incorporated into pick_all_windows.pl.
#
#  INPUT:
#     dir_win_run_meas    directory with MEASURE
#     ofile               MEASUREMENT.WINDOWS output file name
#
#  EXAMPLES:
#    prepare_meas_all.pl MEASURE MEASUREMENT_WINDOWS_9818433_T06
#
#  CALLED BY:
#    pick_all_windows.pl
#
#==========================================================

if (@ARGV < 2) {die("Usage: prepare_meas_all.pl dir_win_run_meas ofile\n")}
($dir_win_run_meas,$ofile) = @ARGV;

#===========================================================================

$ofile2 = "${dir_win_run_meas}/$ofile";

# create the composite windowing file
# STOP if there are no window files in the MEASURE directory
$files = "${dir_win_run_meas}/*mt*";
$file_windows0 = "${dir_win_run_meas}/${ofile}_0";
$file_windows1 = "${dir_win_run_meas}/${ofile}_1";
$file_windows2 = "${dir_win_run_meas}/${ofile}_2";
@windowfiles = glob($files); $numwin = @windowfiles;
if ($numwin == 0) {die("no windows so check $files (prepare_meas_all.pl)")}

# write the C-shell script to file
$cshfile = "prepare_meas_all.csh";
print "\nWriting to $cshfile ...\n";
open(CSH,">$cshfile");

print CSH "echo $numwin seismograms -- BHZ,BHR,BHT\n";
print CSH "cat $files > ${file_windows1}\n";
print CSH "seq $numwin $numwin > ${file_windows0}\n";
print CSH "cat ${file_windows0} ${file_windows1} > ${file_windows2}\n";
print CSH "cp ${file_windows2} $ofile2\n";
print CSH "\\rm ${file_windows0} ${file_windows1} ${file_windows2}\n";

#-----------------------------------------
close(CSH);
system("csh -f $cshfile");

#=================================================================
