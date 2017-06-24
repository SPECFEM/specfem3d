#!/usr/bin/perl -w

# adapted from Carl Tape's version
#

if (@ARGV < 2) {die("Usage: prepare_meas_all.pl dir_win_run_meas ofile\n")}
($dir_win_run_meas,$ofile) = @ARGV;

@windowfiles = glob("${dir_win_run_meas}/*mt*");
$numfile = @windowfiles;
if ($numfile== 0) {die("no windows so check $dir_win_run_meas/*mt*\n")}

# write the C-shell script to file
$cshfile="prepare_meas_all.bash";
open(CSH,">$cshfile");
print CSH "echo $numfile > $ofile \n";

foreach $file (@windowfiles) {print CSH "cat $file >> $ofile \n";}

close(CSH);
system("csh -f $cshfile");
