#!/usr/bin/perl -w

#==========================================================
#
#  run_mt_cc_plot.pl
#  Carl Tape
#  20-Oct-2009
#
#  This runs mt_measure_adj.f90 first in cross-correlation mode,
#  and then in multitaper mode.  The purpose is to have
#  a direct comparison between the two types of adjoint sources.
#  We do not need both for the adjoint tomography inversion.
#
#  EXAMPLE:
#     run_mt_cc_plot.pl m16 -10/200 0 -0.585/0.011/18200 7 BH 6/30 0/1/1/1 -5.0/5.0/-1.5/1.5/0.7 1/1.0/0.5 1/0.02/2.5/2.0/2.5/3.5/1.5
#
#==========================================================

if (@ARGV < 10) {die("Usage: run_mt_cc_plot.pl xxx\n")}
($smodel,$lcut,$itest,$tvec,$imeas,$chan,$Ts,$iparbools,$par1,$par2,$par3) = @ARGV;
#($ibp,$Ts,$lcut,$tvec,$itest,$iparbools,$par1,$par2,$par3) = @ARGV;

$pwd = $ENV{PWD};

($Tmin,$Tmax) = split("/",$Ts);

# time interval to cut for plotting
#$Lstart = -10;
#$Lend = 180;
#if($Tmin==2) {$Lend = 120;}
#$lcut = "$Lstart/$Lend";

if ($itest == 0) {
  $dir_syn  = "SYN"; $dir_data = "DATA"; $dir_meas = "MEASURE"; $dir_recon = "RECON";
} else {
  $dir_syn  = "SYN_TEST"; $dir_data = "DATA_TEST"; $dir_meas = "MEASURE_TEST"; $dir_recon = "RECON_TEST";
}
$plot_recon_dir = "PLOTS/${dir_recon}";

# delete all figures in PLOTS
`rm PLOTS/*ps PLOTS/*jpg PLOTS/*pdf`;

# first run in cross-correlation TT mode
$imeas = 5;
$itaper = 3;

# replace itaper in par3 with this value
(undef,$v2,$v3,$v4,$v5,$v6,$v7) = split("/",$par3);
$par3 = "$itaper/$v2/$v3/$v4/$v5/$v6/$v7";

$plot_adj_dir = "PLOTS/ADJOINT_SOURCES";

#$cmtfile = "CMTSOLUTION_9818433";
#$stafile1 = "EXAMPLE_2/STATIONS_CHT";
$stafile1 = "STATIONS_TOMO";
$stafile2 = "STATIONS_ADJOINT";

# a CMTSOLUTION file should have been copied here from prepare_mt_measure_adj.pl
$cmtfile = glob("CMTSOLUTION_*");

print "\n-- $cmtfile --\n";

#=============================================
# write the C-shell script to file
$cshfile = "run_mt_cc_plot.csh";
print "\nWriting to $cshfile ...\n";
open(CSH,">$cshfile");

# empty the output files directory and compile the measurement code
# (this is done in prepare_measure_adj.pl)
#print CSH "\\rm -rf OUTPUT_FILES ; mkdir OUTPUT_FILES ; make clean\n make\n";

# create MEASUREMENT.PAR file
#$wtr = 0.02; $npi = 2.5; # "default" values
#print CSH "write_par_file.pl OUTPUT_FILES $ibp $Ts $tvec $itaper $imeas $wtr/$npi $iparbools $par1 $par2 $par3\n";
print CSH "write_par_file.pl $tvec $imeas $chan $Ts $iparbools $par1 $par2 $par3\n";

#---------------------------------------------

# run the measurement code
print CSH "measure_adj\n";

# copy adjoint sources into PLOTS
print CSH "\\rm -rf ${plot_adj_dir}\n mkdir ${plot_adj_dir}\n";
print CSH "cp OUTPUT_FILES/*adj ${plot_adj_dir}\n";

# copy reconstructed records into PLOTS
  print CSH "cp OUTPUT_FILES/*recon.cc* ${plot_recon_dir}\n";

#////////////////////////////////////////////////////////

# second run in multitaper TT mode
$imeas = 7;
$itaper = 1;

# replace itaper in par3 with this value
(undef,$v2,$v3,$v4,$v5,$v6,$v7) = split("/",$par3);
$par3 = "$itaper/$v2/$v3/$v4/$v5/$v6/$v7";

# empty the output files directory
print CSH "\\rm -rf OUTPUT_FILES ; mkdir OUTPUT_FILES\n";

# copy chi values into PLOTS
print CSH "cp window_chi PLOTS\n"; 

# create MEASUREMENT.PAR file
print CSH "write_par_file.pl $tvec $imeas $chan $Ts $iparbools $par1 $par2 $par3\n";

# run the measurement code and save output file
# (re-compilation is not necessary)
print CSH "measure_adj > run_file\n";
print CSH "cp OUTPUT_FILES/*adj ${plot_adj_dir}\n";

# create adjoint sources and STATIONS_ADJOINT file for SPECFEM3D
# prepare_adj_src.pl dumps the ZEN adjoint sources into $adj_dir
$adj_dir = "ADJOINT_SOURCES";
print CSH "\\rm -rf ${adj_dir}\n mkdir ${adj_dir}\n";
print CSH "prepare_adj_src.pl -m $cmtfile -z $chan -s PLOTS/$stafile1 -o ${adj_dir} OUTPUT_FILES/*adj\n";
print CSH "cp STATIONS_ADJOINT ${adj_dir}\n";
print CSH "\\mv STATIONS_ADJOINT PLOTS\n";

print CSH "cd PLOTS\n";
print CSH "plot_win_adj_all.pl -l $lcut -m ../$cmtfile -n $chan -b 1 -k $imeas -a $stafile2 -d $dir_data -s $dir_syn -c $dir_recon -w MEASUREMENT.WINDOWS -i $smodel -j $Ts\n";

print CSH "\\rm test*.pdf\n";

# make a single PDF file
#$isort = 3;
#$otag = "mt_cc_all";
#print CSH "make_pdf_by_event.pl $isort $otag\n";

#print CSH "/home/carltape/bin/pdcat -r *.pdf mt_cc_all.pdf\n";     # replace with make_pdf_by_event.pl
print CSH "cd $pwd\n";

#-----------------------------------------
close(CSH);
system("csh -f $cshfile");

#=================================================================
