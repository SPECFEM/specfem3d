#!/usr/bin/perl -w

#==========================================================
#
#  multi_plot_rs.pl
#  Carl Tape, 22-Sept-2005
#
#  This script repeatedly calls plot_rs.pl to make sub record sections.
#  It then concatenates the record sections into a single pdf document.
#
#  1.  Run plot_model.pl
#  2.  Run make_sac_files.pl (generates STATIONS).
#  3.  Run multi_plot_rs.pl
#
#  EXAMPLES (execute from main directory):
#    (02-Aug-2005, waveform)
#    scripts/multi_plot_rs.pl out01 OUTPUT 1 sacdat sacsyn 0/150 10/280 -0.001 -0.001
#
#    (02-Aug-2005, traveltime, xcorr, misfit)
#    scripts/multi_plot_rs.pl out01 OUTPUT 1 sacdat sacsyn 0/150 10/280 -0.001 -5000
#
#    (02-Aug-2005, amplitude, xcorr, misfit)
#    scripts/multi_plot_rs.pl out01 OUTPUT 1 sacdat sacsyn 0/150 10/280 -0.001 -100
#
#    (02-Aug-2005, traveltime, mtm, misfit)
#    scripts/multi_plot_rs.pl out01 OUTPUT 1 sacdat sacsyn 0/150 10/280 -0.001 -3e6
#
#    (02-Aug-2005, amplitude, mtm, misfit)
#    scripts/multi_plot_rs.pl out01 OUTPUT 1 sacdat sacsyn 0/150 10/280 -0.001 -1e5
#
#    (02-Aug-2005, traveltime, xcorr, sampling)
#    scripts/multi_plot_rs.pl out01 OUTPUT 0 junk sacsyn 0/150 10/280 -0.001 -2000
#
#    (02-Aug-2005, amplitude, xcorr, sampling)
#    scripts/multi_plot_rs.pl out01 OUTPUT 0 junk sacsyn 0/150 10/280 -0.001 -700
#
#-----------------------
#
#    (22-Sept-2005, waveform)
#    scripts/multi_plot_rs.pl out01 OUTPUT/run_02/event_01 1 sacdat sacsyn 0/240 10/500 20/50 -0.008 -0.008
#
#    (22-Sept-2005, traveltime, xcorr, misfit)
#    scripts/multi_plot_rs.pl out01 OUTPUT/run_01/event_01 1 sacdat sacsyn 0/240 10/500 20/50 -0.008 -20000
#
#    (06-Oct-2005, traveltime, xcorr, misfit)
#    scripts/multi_plot_rs.pl out01 OUTPUT/run_99/event_11 1 sacdat sacsyn 0/240 10/320 20/50 -0.006 -15000
#
#    (27-Nov-2005, traveltime, xcorr, misfit)
#    scripts/multi_plot_rs.pl out01 OUTPUT/run_000/event_005 1 sacdat sacsyn 0/160 0/300 20/50 -0.006 -15000
#    scripts/multi_plot_rs.pl out01 OUTPUT/run_001/event_001 1 sacdat sacsyn 0/180 300/400 20/10 -0.001 -15000
#
#    (23-May-2006)
#    scripts/multi_plot_rs.pl record_section OUTPUT/run_3100/event_005 1 sacdat sacsyn 0/160 0/300 20/50 -0.001 -5000
#    scripts/multi_plot_rs.pl record_section OUTPUT/run_3101/event_005 1 sacdat sacsyn 0/160 0/300 20/50 -0.001 -3e6
#
#    (19-July-2006, traveltime, xcorr, misfit)
#    scripts/multi_plot_rs.pl record_section_e01 OUTPUT/run_4200/event_001 1 sacdat sacsyn 0/240 0/500/100 20/10 -0.002 -15000
#    scripts/multi_plot_rs.pl record_section_e05 OUTPUT/run_4200/event_005 1 sacdat sacsyn 0/240 0/500/100 20/10 -0.002 -15000
#
#    scripts/multi_plot_rs.pl record_section_e01 OUTPUT/run_4400/event_001 1 sacdat sacsyn 0/240 0/500/100 20/10 -0.002 -15000
#    scripts/multi_plot_rs.pl record_section_e05 OUTPUT/run_4400/event_005 1 sacdat sacsyn 0/240 0/500/100 20/10 -0.002 -15000
#    scripts/multi_plot_rs.pl record_section_e05 OUTPUT/run_4480/event_005 1 sacdat sacsyn 0/240 0/500/100 20/10 -0.002 -1e7
#
#  Type pssac2 for info on scaling the records:
#      -M vertical scaling in sacfile_unit/MEASURE_UNIT = size<required>
#          size: each trace will normalized to size (in MEASURE_UNIT)
#              scale =  (1/size) * [data(max) - data(min)]
#          size/alpha: plot absolute amplitude multiplied by (1/size)*r^alpha
#              where r is the distance range in km across surface
#              specifying alpha = 0.0 will give absolute amplitudes
#              scale = (1/size) * r^alpha
#          size/s: plot absolute amplitude multiplied by (1/size)*sqrt(sin(gcarc))
#              where gcarc is the distance in degrees.
#              scale = (1/size) * sqrt(sin(gcarc))
#
#==========================================================

if (@ARGV < 10) {die("Usage: multi_plot_rs.pl outfile datadir idat dat_suffix syn_suffix tmin/tmax dmin/dmax/dinc ttick/dtick Msyn Madj \n");}
($outfile,$datadir,$idat,$dat_suffix,$syn_suffix,$tran,$dran,$ticks,$Msyn,$Madj) = @ARGV;

$comp = 1;
$alpha = 0;

($ttick,$dtick) = split("/",$ticks);
($tmin,$tmax) = split("/",$tran);
($dmin,$dmax,$dinc) = split("/",$dran);
if(not $dinc) {$dinc = $dmax-$dmin;}

$B = "$ttick/$dtick";

$nump = ($dmax-$dmin)/$dinc;

print "\n Range for time axis     : $tmin $tmax ($ttick)";
print "\n Range for distance axis : $dmin $dmax ($dtick)";
print "\n $nump separate record sections (dinc = $dinc) \n";

# write the C-shell script to file
$cfile = "multi_plot_rs.csh";
print "\nWriting CSH file ($cfile)...\n";
open(CSH,">$cfile");

$t1 = $tmin;
$t2 = $tmax;
$d1 = $dmin;
$d2 = $dmin + $dinc;
for ($i = 1; $i <= $nump; $i = $i+1) {

  $R        = "-R$t1/$t2/$d1/$d2";

  $pslab    = "_${t1}_${t2}_${d1}_${d2}.ps";
  $title1   = "socal_seis";
  $title2   = "socal_stf_adjoint";
  $psfile1 = $title1.$pslab;
  $psfile2 = $title2.$pslab;

  print "\n idat = $idat \n";

  if ($idat == 0) {

     # record section of synthetics only
     print CSH "scripts/plot_rs.pl $R -B$B -Ekt0 -M$Msyn/$alpha -W0.5p $datadir/*$syn_suffix \n";

  } else {

     # record section of data + synthetics
     print CSH "scripts/plot_rs.pl $R -B$B -m $title1 -Ekt0 -M$Msyn/$alpha -W0.5p -w0.5/255/0/0tap -S $datadir,$syn_suffix -D $datadir,$dat_suffix,1,STATIONS \n";
  }

  # record section of adjoint source time functions (NOTE: sac suffix is NOT user input)
  print CSH "scripts/plot_rs.pl $R -B$B -m $title2 -Ekt0 -n \"Socal adjoint source time functions\" -M$Madj/$alpha -W0.5p $datadir/*sacstfadj \n";

  print CSH "ps2pdf $psfile1\n";
  print CSH "ps2pdf $psfile2\n";

  $d1 = $d2;
  $d2 = $d2 + $dinc;
}

print CSH "ps2pdf model.ps\n";
close (CSH);

#die("testing");

system("csh -f -v $cfile");
# system("csh -f $cfile");

#---------------

# bring it all into one pdf file
@files = glob("*.pdf");
if (@files == 0){die("\n STOP: No pdf record sections available to combine\n\n"); }

#print CSH "pdcat @files $outfile.pdf \n\n";
print "\n pdcat -r @files $outfile.pdf \n\n";
system("/home/wagholi2/lqy/local/PDF-CLE/bin/pdcat -r @files $outfile.pdf");

#------------------------------------------------------
