#!/usr/bin/perl -w

#==========================================================
#
#  prepare_measure_adj.pl
#  Carl Tape
#  19-Oct-2009
#
#  This script prepares the data and synthetics for running in the measurement code.
#  At the moment, it is somewhat tailored to the SoCal dataset.
#
#  INPUT:
#    itest          MEASURE_TEST or MEASURE
#    iplot          process data and synthetics for plots
#    iwindow        grab records from windowing directory
#    Tmin/Tmax      bandpass periods for filtering data and synthetics
#
#  EXAMPLE:
#    prepare_measure_adj.pl m16 0 1 0 6/30 9818433
#
#==========================================================

if (@ARGV < 6) {die("Usage: prepare_measure_adj.pl smodel itest iplot iwindow Tmin/Tmax eid\n");}
($smodel,$itest,$iplot,$iwindow,$Ts,$eid) = @ARGV;

# iwindow: controls where the data, synthetics, and windows files are
# 1: grab everything from the windowing directory
# 0: grab everything from a general output directory
#$iwindow = 0;

# label for directory
($Tmin,$Tmax) = split("/",$Ts);
$sTmin = sprintf("T%3.3i",$Tmin);
$sTmax = sprintf("T%3.3i",$Tmax);
$Ttag = "${sTmin}_${sTmax}";

if ($itest == 0) {
  $dir_syn  = "SYN"; $dir_data = "DATA"; $dir_meas = "MEASURE";
} else {
  $dir_syn  = "SYN_TEST"; $dir_data = "DATA_TEST"; $dir_meas = "MEASURE_TEST";
  if($iwindow==0) {die("Exit: must have iwindow = 1 for TEST option.");}
}

#============================================
# USER INPUT

# CMTSOLUTION file
$dir_cmt = "/home/carltape/results/SOURCES/socal_16/v16_files";
$cmtfile = "${dir_cmt}/CMTSOLUTION_${eid}";
if (not -f $cmtfile) {die("check if cmtfile $cmtfile exist or not\n")}

# directories for data, synthetics, and window files : MUST BE CHANGED FOR EACH USER
if ($iwindow==1) {
  # windowing code run directory
  $dir_win_run = "/data2/SVN/seismo/3D/flexwin_run";
  $dir_win_run_syn  = "${dir_win_run}/${dir_syn}";
  $dir_win_run_data = "${dir_win_run}/${dir_data}";
  $win_in = "${dir_win_run}/${dir_meas}/MEASUREMENT_WINDOWS";
  $par_in = "${dir_win_run}/${dir_meas}/PAR_FILE";
  #$dir_win_run_meas = "${dir_win_run}/${dir_meas}";

} else {
  $dir0 = "/home/carltape/RUNS/$eid/$smodel";

  $dir1 = "$dir0/WINDOW_${Ttag}";
  $dir_win_run_syn  = "$dir1/SYN";
  $dir_win_run_data = "$dir1/DATA";
  $win_in = "$dir1/MEASUREMENT_WINDOWS_${eid}_${Ttag}_${smodel}";
  $par_in = "$dir1/PAR_FILE";

  #$dir0 = "/net/sierra/raid1/carltape/socal/socal_3D";
  #$dir_win_run_syn  = "$dir0/SYN/model_m0/${eid}/PROCESSED";
  #$dir_win_run_data = "$dir0/DATA/FINAL/${eid}/PROCESSED";
  #$dir_win_out = "/net/sierra/raid1/carltape/results/WINDOWS/model_m0/";
}

#============================================

if (not -e ${dir_win_run_syn}) {die("check if ${dir_win_run_syn} exist or not\n");}
if (not -e ${dir_win_run_data}) {die("check if ${dir_win_run_data} exist or not\n");}
if (not -e ${win_in}) {die("check if ${win_in} exist or not\n");}
if (not -e ${par_in}) {die("check if ${par_in} exist or not\n");}

# directories for measurement code
$dir_meas_syn  = ${dir_syn};
$dir_meas_data = ${dir_data};
#$dir_meas_meas = ${dir_meas};

#=============================================
# write the C-shell script to file
$cshfile = "prepare_measure_adj.csh";
print "\nWriting to $cshfile ...\n";
open(CSH,">$cshfile");

# remove folders in the measurement directory
print CSH "\\rm -rf ${dir_meas_syn} ${dir_meas_data}\n";
#print CSH "mkdir ${dir_meas_syn} ${dir_meas_data}\n";

# copy ALL data and synthetics into measurement directory (BHZ,BHR,BHT)
# NOTE: copy the whole directory, because the list might be too long for using cp
#       Another option is to LINK the data and syn directories to the local run directory.
print CSH "echo copying files into the measurement directory...\n";
print CSH "cp -r ${dir_win_run_syn} ${dir_meas_syn}\n sleep 2s\n";
print CSH "cp -r ${dir_win_run_data} ${dir_meas_data}\n sleep 2s\n";
#print CSH "cp ${dir_win_run_syn}/* ${dir_meas_syn}\n";
#print CSH "cp ${dir_win_run_data}/* ${dir_meas_data}\n";

# CMTSOLUTION file
print CSH "\\rm CMTSOLUTION*\n";
print CSH "echo $cmtfile\n";
print CSH "cp $cmtfile .\n";

# process data and synthetics for the PLOTTING figures (NOT required for adjoint sources)
# We checked that the processing used in the Perl scripts matches what is used in
# both the windowing code and in the measurement code.
# NOTE: NO EXTENSION IS PUT ON THE BAND-PASSED VERSIONS
if ($iplot == 1) { 
  # synthetics
  $syn_dir_plot = "PLOTS/${dir_syn}";
  print CSH "\\rm -rf ${syn_dir_plot}\n";
  print CSH "cp -r ${dir_meas_syn} ${syn_dir_plot}\n";
  #print CSH "cd ${syn_dir_plot} ; mv_files.pl -x \"*sac*\" \"*sac\" ; cd -\n";  # remove extensions
  #print CSH "mkdir ${syn_dir_plot}\n";
  #print CSH "cd ${dir_meas_syn}\n";      
  #print CSH "process_trinet_syn_new.pl -S -t $Ts -d ../${syn_dir_plot} *sac.d\n"; 
  #print CSH "cd -\n";

  # data
  $data_dir_plot = "PLOTS/${dir_data}";
  print CSH "\\rm -rf ${data_dir_plot}\n";
  print CSH "cp -r ${dir_meas_data} ${data_dir_plot}\n";
  #print CSH "cd ${data_dir_plot} ; mv_files.pl -x \"*sac*\" \"*sac\" ; cd -\n";  # remove extensions
  #print CSH "mkdir ${data_dir_plot}\n";
  #print CSH "cd ${dir_meas_data}\n";
  #print CSH "process_cal_data.pl -t $Ts -d ../${data_dir_plot} *sac.d\n";
  #print CSH "cd -\n";

  # remove any old plots in PLOTS directory
  print CSH "\\rm PLOTS/*pdf PLOTS/*jpg PLOTS/*ps\n";

#  # make a set of SYNTHETIC stations lists sorted by: station, distance, azimuth
#  @tags = ("sta","dist","az");
#  $fall = "${syn_dir_plot}/*";
#  for ($ij = 1; $ij <= 3; $ij = $ij+1) {
#    $tag = $tags[$ij-1];
#    $ofile = "PLOTS/STATIONS_${tag}";
#    print CSH "\\rm -rf $ofile\n";
#    print CSH "saclst kstnm knetwk dist az f $fall | awk '{print \$2\".\"\$3,\$4,\$5}' | uniq | sort -g -k $ij > $ofile\n";
#  }

}

## this version allows for flexibility for the name of the window file
#if ($iwindow==1) {
#  $win_in = "${dir_win_run_meas}/MEASUREMENT*";
#  @files = glob("${win_in}"); $nfiles = @files;
#  if ($nfiles != 1) {die("check if ${win_in} exist or not\n");}
#  print CSH "\\cp @files ${win_out}\n";
#}

# check that the number of records matches the first line of the file
# NOTE: this is important when manually removing windows
($nseis,undef,undef) = split(" ",`grep DATA ${win_in} | wc`);                    # number of listed data records
open(IN,${win_in}); @lines = <IN>; close(IN); $nlist = $lines[0]; chomp($nlist); # header number
print "$eid -- nlist = $nlist, nseis = $nseis\n";
if($nlist != $nseis) {
   print "${win_in}\n";
   die("check window files\n");
}

# copy windowing file and PAR file
$win_out = "./MEASUREMENT.WINDOWS";
$par_out = "./PAR_FILE";
print CSH "\\rm ${win_out}\n";
print CSH "\\cp ${win_in} ${win_out}\n";
print CSH "\\rm ${par_out}\n";
print CSH "\\cp ${par_in} ${par_out}\n";

# copy the windows file to PLOTS
if ($iplot == 1) {print CSH "cp $win_out PLOTS\n"}

# empty the output files directory
print CSH "\\rm -rf OUTPUT_FILES\n";
print CSH "mkdir OUTPUT_FILES\n";

# compile the measurement code
print CSH "make clean\n make\n";

#-----------------------------------------
close(CSH);
system("csh -f $cshfile");
#=================================================================
