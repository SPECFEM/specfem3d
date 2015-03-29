#!/usr/bin/perl -w

#==========================================================
#
#  run_measure_adj.pl
#  Carl Tape
#  19-Oct-2009
#
#  This script runs following prepare_measure.pl.
#  It executes measure_adj and then ouputs a complete set of
#  un-rotated adjoint sources for SPECFEM3D as well as STATIONS_ADJOINT.
#  It also provides the option of making a set of plots by
#  calling plot_win_adj_all.pl.
#
#  INPUT:
#     smodel           index label for model
#     lcut             cut times for plotting records
#     itest            use DATA_TEST/SYN_TEST or DATA/SYN
#     iplot            make plots of data, syn, windows, and adjoint sources
#     iboth            plot both MT and CC adjoint sources (=1) or not (=0)
#     ...              see write_par_file.pl
#
#  EXAMPLE (T = [6s, 30s]):
#     run_measure_adj.pl m16 -10/200 0 1 0 -0.585/0.011/18200 7 BH 6/30 0/1/1/1 -5.0/5.0/-1.5/1.5/0.7 1/1.0/0.5 1/0.02/2.5/2.0/2.5/3.5/1.5  # data + syn + adj
#     run_measure_adj.pl m16 -10/200 0 1 0 -0.585/0.011/18200 7 BH 6/30 0/1/1/0 -5.0/5.0/-1.5/1.5/0.7 1/1.0/0.5 1/0.02/2.5/2.0/2.5/3.5/1.5  # data + syn
#
#==========================================================

if (@ARGV < 13) {die("Usage: run_measure_adj.pl xxx\n")}
($smodel,$lcut,$itest,$iplot,$iboth,$tvec,$imeas,$chan,$Ts,$iparbools,$par1,$par2,$par3) = @ARGV;

$pwd = $ENV{PWD};

$plot_adj_dir = "PLOTS/ADJOINT_SOURCES";
#$iboth = 0;   # option to plot both cross-correlation and multitaper adjoint sources

($Tmin,$Tmax) = split("/",$Ts);

# whether to compute (and plot) adjoint sources
# note: If you want to compute adjoint sources but only plot the data+syn,
#       then you can make manual changes to $iplot_adj in plot_win_adj.pl.
(undef,undef,undef,$iadj) = split("/",$iparbools);

#$cmtfile = "CMTSOLUTION_9818433";
#$stafile1 = "EXAMPLE_2/STATIONS_CHT";
$stafile1 = "STATIONS_TOMO";
$stafile2 = "STATIONS_ADJOINT";

# a CMTSOLUTION file should have been copied here from prepare_measure_adj.pl
$cmtfile = glob("CMTSOLUTION_*");

print "\n-- $cmtfile --\n";

if ($itest == 0) {
   $dir_syn  = "SYN"; $dir_data = "DATA"; $dir_meas = "MEASURE"; $dir_recon = "RECON";
} else {
   $dir_syn  = "SYN_TEST"; $dir_data = "DATA_TEST"; $dir_meas = "MEASURE_TEST"; $dir_recon = "RECON_TEST";
}
$plot_recon_dir = "PLOTS/${dir_recon}";

#=============================================
# write the C-shell script to file
$cshfile = "run_measure_adj.csh";
print "\nWriting to $cshfile ...\n";
open(CSH,">$cshfile");

# empty the output files directory and compile the measurement code
# (this is done in prepare_measure_adj.pl)
#print CSH "\\rm -rf OUTPUT_FILES ; mkdir OUTPUT_FILES ; make clean\n make\n";

#----------------------------
# choose ipar1 based on the values used in the windowing code PAR_FILE (Feb 2009)
# This is a more sensible option, since you guarantee consistency with FLEXWIN.
if(1==1) {
   # Allow for some "loosening" of the parameters, since the exact nature
   # of the measurement in the measurement code may exclude some input windows.
   $dt_fac = 0.5;      # in seconds
   $cc_fac = 0.02;     # from 0.0 to 1.0

   (undef,undef,$tshift_base) = split(" ",`grep TSHIFT_BASE PAR_FILE`);
   (undef,undef,$tshift_ref) = split(" ",`grep TSHIFT_REFERENCE PAR_FILE`);
   (undef,undef,$dlna_base) = split(" ",`grep DLNA_BASE PAR_FILE`);
   (undef,undef,$dlna_ref) = split(" ",`grep DLNA_REFERENCE PAR_FILE`);
   (undef,undef,$cc_base) = split(" ",`grep CC_BASE PAR_FILE`);

   $TSHIFT_MIN  = $tshift_ref - $tshift_base - $dt_fac;
   $TSHIFT_MAX  = $tshift_ref + $tshift_base + $dt_fac;
   $DLNA_MIN    = $dlna_ref - $dlna_base;
   $DLNA_MAX    = $dlna_ref + $dlna_base;
   $CC_MIN = $cc_base - $cc_fac;
   if($CC_MIN < 0) {$CC_MIN = 0}

   $par1_old = $par1;
   $par1 = "${TSHIFT_MIN}/${TSHIFT_MAX}/${DLNA_MIN}/${DLNA_MAX}/${CC_MIN}";
   print "Updating values from PAR_FILE:\n";
   print " old -- $par1_old --\n new -- $par1 --\n";
}
#----------------------------

# create MEASUREMENT.PAR file
#$wtr = 0.02; $npi = 2.5;   # "default" values
#$itaper = 3;
#print CSH "write_par_file.pl OUTPUT_FILES $ibp $Ts $tvec $itaper $imeas $wtr/$npi $iparbools $par1 $par2 $par3\n";
print CSH "write_par_file.pl $tvec $imeas $chan $Ts $iparbools $par1 $par2 $par3\n";

#---------------------------------------------

# run the measurement code
print CSH "measure_adj > run_file\n";    # output to file
#print CSH "measure_adj\n";

if ($iplot == 1) {
  # copy adjoint sources into PLOTS
  if($iboth==0) {print CSH "\\rm -rf ${plot_adj_dir}\n mkdir ${plot_adj_dir}\n";}
  print CSH "cp OUTPUT_FILES/*adj ${plot_adj_dir}\n";

  # copy reconstructed records into PLOTS
  print CSH "cp OUTPUT_FILES/*recon.cc* ${plot_recon_dir}\n";

  # copy chi values into PLOTS
  print CSH "cp window_chi PLOTS\n";

  # convert to SAC files (creates *.sac files)
  #print CSH "ascii2sac.csh ${plot_adj_dir}*.adj\n";
}

# create adjoint sources and STATIONS_ADJOINT file for SPECFEM3D
# prepare_adj_src.pl dumps the ZEN adjoint sources into $adj_dir
# NOTE: You may want plots of data+syn+windows even if you do not compute adjoint sources;
#       hence you can use the *recon.cc.sac to write a STATIONS_ADJOINT file for plot_win_adj.pl.
if ($iadj == 1) {
  $tfiles = "-i OUTPUT_FILES/*adj";
} else {
  $tfiles = "OUTPUT_FILES/*recon.cc.sac";
}
$adj_dir = "ADJOINT_SOURCES";
print CSH "\\rm -rf ${adj_dir}\n mkdir ${adj_dir}\n";
print CSH "prepare_adj_src.pl -m $cmtfile -s PLOTS/$stafile1 -o ${adj_dir} $tfiles\n";
print CSH "cp STATIONS_ADJOINT ${adj_dir}\n";  # only for iadj=1
print CSH "\\mv STATIONS_ADJOINT PLOTS\n";

# make plots of (filtered) data, synthetics, windows, and adjoint sources
if ($iplot == 1) {
  print CSH "cd PLOTS\n";
  print CSH "plot_win_adj_all.pl -l $lcut -m ../$cmtfile -n $chan -b $iboth -k $imeas/$iadj -a $stafile2 -d $dir_data -s $dir_syn -c $dir_recon -w MEASUREMENT.WINDOWS -i $smodel -j $Ts\n";
  print CSH "cd $pwd\n";
}

#-----------------------------------------
close(CSH);
system("csh -f $cshfile");

#=================================================================
