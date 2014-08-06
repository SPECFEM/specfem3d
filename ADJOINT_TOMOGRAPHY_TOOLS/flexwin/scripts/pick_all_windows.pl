#!/usr/bin/perl -w

#==========================================================
#
#  pick_all_windows.pl
#  Carl Tape
#  01-Feb-2009
#
#  The script reads in a set of data directories, each labeled as an event ID,
#  and it copies the data and synthetics into the windowing code, runs the code,
#  and outputs a PDF file of the window picks.
#
#  For working on a single event, it is easier to first run this scripts,
#  but for subsequent events, run pick_all_windows_local.pl.
#
#  INPUT:
#    Tmin/Tmax     bandpass periods for filtering data and synthetics
#                  --> as of 07-June-2008, this data and syn can be PRE-BANDPASS FILTERED
#    imin/imax     index of first event and last event to run
#
#    idebug        boolean: debugging information
#    iplot         boolean: make composite plots for all windows
#    imeas         boolean: make *mt_input files for measurement code
#    ibody         boolean: turn body-wave mode ON
#
#    idataset      integer: which dataset to use (1 for socal; 2 for japan)
#    iexecute      boolean: whether to execute the Shell script
#
#  DIRECTORIES:
#    dir_data      directory containing data directories
#    dir_syn       directory containing synthetics directories
#    dir_win_code  directory for compiling the windowing code
#    dir_win_run   directory for running the windowing code
#    odir          output directory for figures and other files
#    
#  CALLS:
#     plot_windows_all.pl
#     prepare_meas_all.pl
#
#  EXAMPLES:
#    pick_all_windows.pl m16 0 6/30   1/300 1/1/1/0 1 0     # make plots and WINDOWS file, T = 6-30s
#    pick_all_windows.pl m16 0 3/30   1/300 1/1/1/0 1 0     # make plots and WINDOWS file, T = 3-30s
#    pick_all_windows.pl m16 0 2/30   1/300 1/1/1/1 1 0     # make plots and WINDOWS file, T = 2-30s
#
#    pick_all_windows.pl m00 0 6/30 179/179 1/1/0/0 1 0     # make plots only, T = 6-30s
#    pick_all_windows.pl m00 0 6/30 179/179 0/0/1/0 1 0     # make WINDOWS file, T = 6-30s
#
#    pick_all_windows.pl m16 0 3/10 198/198 1/1/1/0 1 0     # socal 14236768
#
#==========================================================

if (@ARGV < 7) {die("Usage: pick_all_windows.pl model Tmin/Tmax imin/imax idebug/iplot/imeas/ibody idataset iexecute\n");}
($smodel,$ibp,$Ts,$inds,$ibools,$idataset,$iexecute) = @ARGV;

# plot a figure for every record -- this takes a long time (iplot = 1 for this option)
$iplotall = 0;

# split input entries
($Tmin,$Tmax)  = split("/",$Ts);
($imin,$imax)  = split("/",$inds);

# string for PDF
$sTmin = sprintf("T%3.3i",$Tmin);
$sTmax = sprintf("T%3.3i",$Tmax);
$Ttag = "${sTmin}_${sTmax}";
$sTminD = sprintf("%.2f",$Tmin);   # double precision for PAR_FILE
$sTmaxD = sprintf("%.2f",$Tmax);

# split boolean entries
($idebug,$iplot,$imeas,$ibody)  = split("/",$ibools);
if($idebug==1) {$bdebug = ".true."} else {$bdebug = ".false."}
if($iplot==1) {$bplot = ".true."} else {$bplot = ".false."}
if($imeas==1) {$bmeas = ".true."} else {$bmeas = ".false."}
if($ibody==1) {$bbody = ".true."} else {$bbody = ".false."}
if($ibp==1) {$bbp = ".true."} else {$bbp = ".false."}

#-------------------------------------
# USER INPUT -- CHANGE FOR EACH USER

# directories for data and synthetics
# 1=socal; 2=Japan
if($idataset == 1) {
  #$dir_data  = "/home/carltape/SOCAL_ADJOINT/DATA/FINAL";
  #$dir_syn  = "/home/carltape/SOCAL_ADJOINT/SYN/model_${smodel}";
  #$dir_data  = "/data2/data/TEST";
  $dir_data  = "/data2/data/calif/FINAL";  
  $dir_syn  = "/data2/syn/socal/model_${smodel}";  
  $sdataset = "socal";
} elsif ($idataset == 2) {
  $dir_data  = "/home/mchen/DATA/TEST";
  $dir_syn  = "/home/mchen/SEM/TEST";
  $sdataset = "japan";
} else {
  die("\n idataset must be 1 (socal) or 2 (japan)\n");
}

# specify various directories (MUST BE MODIFIED FOR EACH USER)
# note: may want a second copy to run simultaneously
$dir0 = "/data2/SVN/seismo/3D";    # path to CIG SVN repository
$dir_win_code = "$dir0/ADJOINT_TOMO/flexwin_work";
$dir_win_run  = "$dir0/flexwin_run";
#$dir_win_code = "$dir0/ADJOINT_TOMO/flexwin_work_copy2";
#$dir_win_run  = "$dir0/flexwin_run_copy2";
$dir_scripts  = "${dir_win_code}/scripts";

# directory to collect COPIES of various output files
$odir = "/home/carltape/results/WINDOWS/model_${smodel}";

# run directory for windows, measurements, adjoint sources, and kernels
#$rundir = "/home/carltape/SOCAL_ADJOINT/RUNS";
$rundir = "/data2/RUNS";

# channels to use for data and synthetics
# NOTE: update prepare_input also ($pinput)
#@chans = ("BHZ","BHT","BHR","HHZ","HHR","HHT");
@dchans = ("BHZ","BHT","BHR","HHZ","HHR","HHT");
#@schans = ("BHZ","BHT","BHR","BHZ","BHR","BHT");   # original
@schans = ("HXZ","HXT","HXR","HXZ","HXR","HXT");
$nchan = @dchans;

#-------------------------------------

# scripts for preparing data and synthetics
$dir_prepare = "${dir_scripts}/prepare_scripts/${sdataset}";

# scripts for user files
$dir_user = "${dir_win_code}/user_files";

# check that directories exist
if (not -e ${dir_data}) {die("check if ${dir_data} exist or not\n")}
if (not -e ${dir_syn}) {die("check if ${dir_syn} exist or not\n")}
if (not -e ${dir_win_code}) {die("check if ${dir_win_code} exist or not\n")}
if (not -e ${dir_win_run}) {die("check if ${dir_win_run} exist or not\n")}
if (not -e ${dir_scripts}) {die("check if ${dir_scripts} exist or not\n")}
if (not -e ${odir}) {die("check if ${odir} exist or not\n")}
if (not -e ${dir_prepare}) {die("check if ${dir_prepare} exist or not\n")}
if (not -e ${rundir}) {die("check if ${rundir} exist or not\n")}

# data and synthetic files (MUST BE MODIFIED FOR EACH USER)
if ($idataset == 1) {
   #$suffix_syn = "semd.sac.d.${Ttag}";
   $suffix_syn = "semd.sac.${smodel}.${Ttag}";
   $suffix_dat = "sac.d.${Ttag}";

} elsif ($idataset == 2) {
  $suffix_syn = "bp12_150s";
  $suffix_dat = "bp12_150s";
}

# directories for windowing code
$dir_win_run_syn  = "${dir_win_run}/SYN";
$dir_win_run_data = "${dir_win_run}/DATA";
$dir_win_run_meas = "${dir_win_run}/MEASURE";

# EVENT LIST
$eid_list = "/home/carltape/results/EID_LISTS/syn_run_${smodel}";
#$eid_list = "/home/carltape/results/EID_LISTS/syn_run_iterate";
if (not -f $eid_list) {die("check if eid_list ${eid_list} exist or not\n")}
open(IN,$eid_list); @eids = <IN>; close(IN);
$nevent0 = @eids;

## obtain list of events
##@datadirs = `ls -1 -d ${dir_data}/[1-9]* | sort -g`;
#@datadirs = glob("${dir_data}/[1-9]*");
#$nevent0 = @datadirs;

# adjust indices for events
print "\n imin = $imin ; imax = $imax ; nevent = $nevent0";
if($imin < 1) {$imin = 1}
if($imin > $nevent0) {$imin = $nevent0}
if($imax > $nevent0) {$imax = $nevent0}
$nevent = $imax - $imin + 1;
print "\n imin = $imin ; imax = $imax ; nevent = $nevent\n";
print "\n $nevent events to use in the windowing file\n";

for ($ievent = $imin; $ievent <= $imax; $ievent++) {
  #$datadir1 = $datadirs[$ievent-1]; chomp($datadir1);
  #$eid = `basename $datadir1`; chomp($eid);   # event ID
  $eid = $eids[$ievent-1]; chomp($eid);
  print "\n Event $ievent : $eid";
}
#die("testing");

#--------------------------------------------------------------

# user parameters for windowing code
# NOTE: user must first copy their PAR_FILE from the user_files directory
$par_file = "${dir_win_code}/PAR_FILE";
$userfun_file = "${dir_win_code}/user_functions.f90";

# copy PAR_FILE into local directory
# directories for data and synthetics
if($idataset == 1) {
   $par_file_in = "${dir_user}/socal_3D/PAR_FILE_${Ttag}_${smodel}";
   $userfun_file_in = "${dir_user}/socal_3D/user_functions_${smodel}.f90";

} elsif ($idataset == 2) {
   $par_file_in = "${dir_user}/japan_3D/PAR_FILE_${Ttag}";
   $userfun_file_in = "${dir_user}/japan_3D/user_functions.f90";
}
if (not -f ${par_file_in}) {die("check if PAR_FILE ${par_file_in} exist or not\n")}
if (not -f ${userfun_file_in}) {die("check if PAR_FILE ${userfun_file_in} exist or not\n")}

# copy user files to code directory
`cp ${par_file_in} ${par_file}`;
`cp ${userfun_file_in} ${userfun_file}`;

# copy IASP files into run directory
$iasp_in1 = "${dir_win_code}/ttimes_mod/iasp91.tbl";
$iasp_in2 = "${dir_win_code}/ttimes_mod/iasp91.hed";
$iasp_out1 = "${dir_win_run}/iasp91.tbl";
$iasp_out2 = "${dir_win_run}/iasp91.hed";
`cp $iasp_in1 $iasp_out1`;
`cp $iasp_in2 $iasp_out2`;

# name of executable windowing code (in $dir_win_code)
$win_execute = "flexwin";

# make command for windowing code
$make = "make -f make_gfortran";   # change as of 9-30-08
#$make = "make -f make_intel_caltech";

# location of plot_windows_all.pl
$plot_windows_perl = "${dir_scripts}/plot_windows_all.pl";
if (not -f ${plot_windows_perl}) {die("check if ${plot_windows_perl} exist or not\n")}

# location of prepare_meas_all.pl
$prepare_meas_perl = "${dir_scripts}/prepare_meas_all.pl";
if (not -f ${prepare_meas_perl}) {die("check if ${prepare_meas_perl} exist or not\n")}

# script to summarize window picks
if ($idataset == 1) {
   $extract_script = "${dir_scripts}/extract_event_windowing_stats_carl.sh";
} elsif ($idataset == 2) {
   $extract_script = "${dir_scripts}/extract_event_windowing_stats_min.sh";
} 
if (not -f ${extract_script}) {die("check if ${extract_script} exist or not\n")}

# sorted receiver files
$reclist0      = "${dir_win_run_meas}/rectext";
$reclist_dist  = "${reclist0}_dist";
$reclist_picks = "${reclist0}_picks";

#$nevent = 2;     # testing

#===========================================================================
# write the C-shell script to file
$cshfile = "pick_all_windows.csh";
print "\nWriting to $cshfile ...\n";
open(CSH,">$cshfile");

# for the first event ($imin), change the user parameters file and compile
print CSH "cd ${dir_win_code}\n";
print CSH "sed '/^RUN_BANDPASS                    =/s/^.*\$/RUN_BANDPASS                    = $bbp/' ${par_file} > temp0\n";
print CSH "sed '/^WIN_MIN_PERIOD                  =/s/^.*\$/WIN_MIN_PERIOD                  = $sTminD/' temp0 > temp1\n";
print CSH "sed '/^WIN_MAX_PERIOD                  =/s/^.*\$/WIN_MAX_PERIOD                  = $sTmaxD/' temp1 > temp2\n";
print CSH "sed '/^DEBUG                           =/s/^.*\$/DEBUG                           = $bdebug/' temp2 > temp3\n";
print CSH "sed '/^MAKE_SEISMO_PLOTS               =/s/^.*\$/MAKE_SEISMO_PLOTS               = $bplot/' temp3 > temp4\n";
print CSH "sed '/^MAKE_WINDOW_FILES               =/s/^.*\$/MAKE_WINDOW_FILES               = $bmeas/' temp4 > temp5\n";
print CSH "sed '/^BODY_WAVE_ONLY                  =/s/^.*\$/BODY_WAVE_ONLY                  = $bbody/' temp5 > temp6\n";
print CSH "\\mv temp6 ${par_file}\n \\rm temp0 temp1 temp2 temp3 temp4 temp5\n";

# compile the windowing code
print CSH "cd ${dir_win_code}\n $make clean\n $make ${win_execute}\n";

# copy prepare_seis, prepare_input, prepare_input_test, PAR_FILE into window directory
$pseis       = "prepare_seis.pl";
$pinput      = "prepare_input";
$pinput_test = "prepare_input_test";

print CSH "\\cp ${dir_prepare}/${pseis} ${dir_win_run}/${pseis}\n";
print CSH "\\cp ${dir_prepare}/${pinput} ${dir_win_run}/${pinput}\n";
print CSH "\\cp ${dir_prepare}/${pinput_test} ${dir_win_run}/${pinput_test}\n";
print CSH "\\cp ${par_file} ${dir_win_run}\n";

# copy PAR_FILE and user_functions.f90 to output directory
print CSH "\\cp ${par_file} ${userfun_file} ${dir_prepare}/${pseis} $odir\n";

#die("TESTING");

#===========================================================================

$kk = 0;

#foreach $file (@datadirs) {
for ($ievent = $imin; $ievent <= $imax; $ievent++) {

  # directories
  $eid = $eids[$ievent-1]; chomp($eid);
  $datadir1 = "${dir_data}/$eid";
  #$datadir1 = $datadirs[$ievent-1]; chomp($datadir1);
  #$eid = `basename $datadir1`; chomp($eid); # event ID

  $kk = $kk + 1;
  print "Event ID is $eid -- $kk out of $nevent\n";

  if ($idataset == 1) {
    $datadir1 = "${dir_data}/${eid}/PROCESSED/PROCESSED_${Ttag}";
    $syndir1 = "${dir_syn}/${eid}/PROCESSED/PROCESSED_${Ttag}";
  } elsif ($idataset == 2) {
    $datadir1 = "$datadir1/bp12_150s";
    $syndir1 = "${dir_syn}/${eid}/bp12_150s";
  }

  # First, check if the data and synthetics directories exist.
  # If so, then check if the windowing code has already been run.
  if ( (not -e ${datadir1}) || (not -e ${syndir1}) ) {
    if (not -e ${datadir1}) {
      print "--> check if ${datadir1} exist\n";
    }
    if (not -e ${syndir1}) {
      print "--> check if ${syndir1} exist\n";
    }

  } else {	# both data and synthetic directories exist

    print CSH "echo Event ID is $eid -- $kk out of $nevent\n";

    $eout = "${eid}_${Ttag}_${smodel}";	                # output file tag
    $savedir = "$rundir/$eid/$smodel/WINDOW_${Ttag}";   # main run directory

    # If the windowing code has already been run, then move to the next event.
    # NOTE: to overwrite, just comment out this option
    if (-e $savedir) {
    #if (-f "$savedir/${eout}_runfile") {
      print "--> $savedir/${eout}_runfile exists -- on to next event\n";

    } else {
      print "--> do this event\n";
      print CSH "echo Data dir is $datadir1\n";
      print CSH "echo  Syn dir is $syndir1\n";

      # remove folders in the windowing directory
      print CSH "\\rm -rf ${dir_win_run_syn} ${dir_win_run_data} ${dir_win_run_meas} \n";
      print CSH "mkdir ${dir_win_run_syn} ${dir_win_run_data} ${dir_win_run_meas} \n";
      

      # copy data and synthetics into windowing directory (BHZ,BHR,BHT)
      print CSH "echo copying data and syn files to ${dir_win_run}\n";
      for ($j = 1; $j <= $nchan; $j++) {
	$dchan = $dchans[$j-1];
	$schan = $schans[$j-1];
	print CSH "cp $syndir1/*${schan}.${suffix_syn} ${dir_win_run_syn}\n"; # synthetics
	print CSH "cp $datadir1/*${dchan}.${suffix_dat} ${dir_win_run_data}\n"; # data
      }

      #  print CSH "cp $syndir1/$sfile1  ${dir_win_run_syn}\n";
      #  print CSH "cp $syndir1/$sfile2 ${dir_win_run_syn}\n";
      #  print CSH "cp $syndir1/$sfile3 ${dir_win_run_syn}\n";
      #  print CSH "cp $datadir1/$dfile1 ${dir_win_run_data}\n";
      #  print CSH "cp $datadir1/$dfile2 ${dir_win_run_data}\n";
      #  print CSH "cp $datadir1/$dfile3 ${dir_win_run_data}\n";

      # prepare sac files for windowing code (BHZ,BHR,BHT)
      # NOTE: the file extensions must be checked within prepare_seis.pl and prepare_input
      print CSH "cd ${dir_win_run}\n";
      #print CSH "${pseis}\n";    # now this is done in PRE-PROCESSING
      print CSH "${pinput}\n";

      # copy prepared data and synthetic files into a directory for measurement code
      if ($imeas == 1) {
	#$savedir = "$rundir/$eid/$smodel/WINDOW_${Ttag}";
	print CSH "echo copying prepared DATA and SYN to ${savedir}...\n";
	print CSH "mkdir -p $rundir/$eid\n";
	print CSH "mkdir -p $rundir/$eid/$smodel/\n";
	print CSH "mkdir -p $savedir\n";
	print CSH "\\rm -rf $savedir/SYN ; cp -r ${dir_win_run_syn} $savedir\n";
	print CSH "\\rm -rf $savedir/DATA ; cp -r ${dir_win_run_data} $savedir\n";
      }
 
      # copy input file into output directory
      $input_in = "${dir_win_run}/input";
      $input_out = "${odir}/${eout}_input";
      print CSH "\\cp ${input_in} ${input_out}\n"; # list of matching records
      print CSH "\\cp ${input_in}_bad ${input_out}_bad\n"; # list of probably bad datafiles
      print CSH "\\cp ${input_in}_no_syn ${input_out}_no_syn\n"; # list of missing synthetics

      # run the windowing code
      # NOTE: The >& command will record the premature exit messages in the run_file.
      $runfile = "${dir_win_run_meas}/${eout}_runfile";
      print CSH "echo running windowing code ${win_execute} -- output file is $runfile\n";
      print CSH "${dir_win_code}/${win_execute} < input >& $runfile\n";

      # copy run file to output directory
      print CSH "\\cp $runfile $odir\n";

      if ($iplot == 1) {
	# run extract_event_windowing_stats_carl.sh to get rectext and windows summary figure
	print CSH "${extract_script} MEASURE\n";

	# sort rectext into rectext_dist and rectext_picks
	print CSH "sort -g -k 2 $reclist0 > ${reclist_dist}\n";

	# copy code output and receiver list output to output directory
	$reclist_out = "$odir/${eout}_rectext_dist";
	print CSH "\\cp ${reclist_dist} ${reclist_out}\n";

	# generate composite PDF file using plot_windows_all.pl
        # NOTE: This takes a long time, and is not really necessary
        #       if have already determined the windowing code parameters
        #       and you are planning to plot them in mt_measure_adj.
        if($iplotall == 1) {
	   $pdffile = "$odir/${eout}_all.pdf";
	   print CSH "${plot_windows_perl} ${dir_win_code} ${dir_win_run_meas} ${reclist_dist} $pdffile\n";
	 }

	# copy record section to output directory
	$rs_in = "${dir_win_run_meas}/event_recordsection.pdf";
	$rs_out = "${odir}/${eout}_rs.pdf";
	print CSH "\\cp ${rs_in} ${rs_out}\n";

	# copy window statistics figure to output directory
	$stats_in = "${dir_win_run_meas}/event_winstats.pdf";
	$stats_out = "${odir}/${eout}_winstats.pdf";
	print CSH "\\cp ${stats_in} ${stats_out}\n";
      }

      if ($imeas == 1) {
	# generate MEASUREMENT.WINDOWS file from *mt_input files
	$ofile = "MEASUREMENT_WINDOWS_${eout}";
	print CSH "${prepare_meas_perl} ${dir_win_run_meas} $ofile\n";

	# copy window file to another location
	print CSH "\\cp ${dir_win_run_meas}/${ofile} ${odir}/${ofile}\n";

	# leave copy in MEASURE directory (use the generic name called by mt_measure_adj.f90)
	print CSH "\\cp ${dir_win_run_meas}/${ofile} ${dir_win_run_meas}/MEASUREMENT_WINDOWS\n";

	# copy all window output files, including figures, into RUN directory
	print CSH "\\cp ${odir}/*${eout}* $savedir\n";
	print CSH "\\cp ${savedir}/${ofile} ${savedir}/${ofile}_orig\n";    # extra copy
	print CSH "\\cp ${par_file} ${userfun_file} ${dir_prepare}/${pseis} $savedir\n";
      }

    }

  }
}

#-----------------------------------------
close(CSH);
if($iexecute==1) {system("csh -f $cshfile");}
print "csh -f pick_all_windows.csh\n\n";

#=================================================================
