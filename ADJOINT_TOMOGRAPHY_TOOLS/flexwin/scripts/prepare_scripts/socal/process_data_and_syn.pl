#!/usr/bin/perl -w

#==========================================================
#
#  Carl Tape
#  19-Oct-2009
#  process_data_and_syn.pl
#
#  This script processes data and 3D-Socal-SEM synthetics for Southern California.
#  It calls two Caltech in-house SAC-based Perl processing scripts
#       process_cal_data.pl 
#       process_trinet_syn_new.pl
#  The general procedure is advised for pre-processing data and synthetics for
#  seismic tomography applications.
#
#  INPUT:
#    CMT_list  list of event ID labels in a file
#    dir_syn   directory for synthetics
#    dir_dat   directory for data
#    stafile   STATIONS file ("none" for default)
#    hdur      shortest period for filtering data and synthetics
#
#    syn_ext   extension for processed synthetics
#    dat_ext   extension for processed data
#    Trange    band-pass filter
#
#  ORDER OF OPERATIONS:
#    ~/UTILS/process_data_and_syn.pl 1 m16 1 0 d 6/30 # data, initial pre-processing
#    ~/UTILS/process_data_and_syn.pl 1 m16 0 1 d 6/30 # syn, initial pre-processing (REPEAT)
#    ~/UTILS/process_data_and_syn.pl 2 m16 1 1 d 6/30 # both, create cut files
#    ~/UTILS/process_data_and_syn.pl 3 m16 1 0 d 6/30 # data, execute cut file
#    ~/UTILS/process_data_and_syn.pl 3 m16 0 1 d 6/30 # syn, execute cut file
#
#    ~/UTILS/process_data_and_syn.pl 4 m16 1 0 d 6/30 # data, bandpass T=6-30
#    ~/UTILS/process_data_and_syn.pl 4 m16 0 1 d 6/30 # syn, bandpass T=6-30
#
#    ~/UTILS/process_data_and_syn.pl 4 m16 1 0 d 3/30 # data, bandpass T=3-30
#    ~/UTILS/process_data_and_syn.pl 4 m16 0 1 d 3/30 # syn, bandpass T=3-30
#
#    ~/UTILS/process_data_and_syn.pl 4 m16 1 0 d 2/30 # data, bandpass T=2-30
#    ~/UTILS/process_data_and_syn.pl 4 m16 0 1 d 2/30 # syn, bandpass T=2-30
#
#  FILE-NAMING CONVENTIONS (example for event 9818433):
#    DATA                                 -- main data directory
#       9818433                           -- event directory
#          9818433.CI.PAS.BHZ.sac         -- raw record
#          9818433_dat_cut_done           -- list of data files to cut
#          9818433_syn_cut                -- list of syn files to cut
#          PROCESSED                      -- initial processing directory
#             9818433.CI.PAS.BHZ.sac.d    -- sac headers, interpolated, cut
#             PROCESSED_T006_T030         -- bandpass directory
#                 9818433.CI.PAS.BHZ.sac.d.T006_T030  -- bandpass (and rotated)
#
#    SYN                                  -- main synthetics directory
#       model_m00                         -- model iteration
#          9818433                        -- event directory
#             PAS.CI.BHZ.semd.sac         -- "raw" synthetic record
#             9818433_syn_cut_done        -- list of synthetics that were cut
#             PROCESSED                   -- initial processing directory
#                PAS.CI.BHZ.semd.sac.m00  -- sac headers, hdur convolution, interpolated, cut
#                PROCESSED_T006_T030      -- bandpass directory
#                    PAS.CI.BHZ.semd.sac.m00.T006_T030 -- bandpass (and rotated)
#
#==========================================================

if (@ARGV < 6) {die("Usage: process_data_and_syn.pl iprocess smodel idata isyn dat_ext Tmin/Tmax\n")}
($iprocess,$smodel,$idata,$isyn,$dat_ext,$Trange) = @ARGV;

#-------------------------------------
# USER INPUT

# USER PARAMETERS
$sps = 20;               # samples per second for interpolation/sub-sampling
$tfac = 1.0;             # factor to extend lengths of records (should be > 1.0)
$bmin = -40;             # minimum allowable time before origin time for cut records
$chdur = 0;              # convolve synthetics with Gaussian half duration (default = 1)

$pdir = "PROCESSED";     # tag for processed directories
$syn_ext = $smodel;      # KEY: include model index for comparison among synthetics from different models
$syn_suffix0 = "semd.sac";                  # suffix for synthetic files
$syn_suffix = "${syn_suffix0}.${syn_ext}";
$dat_suffix0 = "sac";                       # suffix for data files
$dat_suffix = "${dat_suffix0}.${dat_ext}"; 

$iexecute = 0;           # execute CSH file immediately
$ilist = 1;              # read EIDs from a list (=1) or grab all CMTSOLUTION files (=0)

# EID list
$CMT_list = "/home/carltape/SOCAL_ADJOINT/SYN/model_${smodel}";
$dir_source = "/home/carltape/results/SOURCES/socal_16";
$dirCMT    = "${dir_source}/v16_files";
$CMT_list  = "${dir_source}/EIDs_only_loc";
#$CMT_list  = "/home/carltape/results/EID_LISTS/syn_run_${smodel}";
#$CMT_list  = "/home/carltape/results/EID_LISTS/syn_run_m12";

# directories
$dirdat0    = "/home/carltape/SOCAL_ADJOINT/DATA/FINAL";
$dirsyn0    = "/home/carltape/SOCAL_ADJOINT/SYN/model_${smodel}";
#$dirsyn0 = "/home/carltape/SOCAL_ADJOINT/SYN/model_pre_${smodel}";
# STATIONS file
#$stafile = "/home/carltape/gmt/stations/seismic/Matlab_output/STATIONS";
$stafile = "/home/carltape/gmt/stations/seismic/Matlab_output/STATIONS_CALIFORNIA_TOMO_INNER_specfem";
#$stafile = "/home/carltape/gmt/stations/seismic/Matlab_output/STATIONS_CALIFORNIA_TOMO_OUTER_specfem";

#-------------------------------------

# check if the files exist
if (not -f $stafile) {print "\n check if $stafile exists";}
if (not -f ${CMT_list}) {print "\n check if CMT_list ${CMT_list} exists";}

# period range
($Tmin,$Tmax) = split("/",$Trange);
$sTmin = sprintf("T%3.3i",$Tmin);
$sTmax = sprintf("T%3.3i",$Tmax);
$Ttag = "${sTmin}_${sTmax}";
$pdirbpass = "${pdir}_${Ttag}";
$ftag = "${Ttag}_${smodel}";

# grab all the CMT solution files or read in the list of event IDs
if($ilist == 1) {
  open(IN,"${CMT_list}"); @cmtfiles = <IN>; close(IN);
} else {
  @cmtfiles = glob("$dirCMT/CMTSOLUTION*");
}
$ncmt = @cmtfiles;

print "\n $dirCMT \n $dirsyn0 \n $dirdat0\n";
print "\n $ncmt CMTSOLUTION files\n";
if($ncmt == 0) {die("\n No CMTSOLUTION files")}

#$ncmt = 1;   # testing
if($ilist == 1) {
  for ($ievent = 1; $ievent <= $ncmt; $ievent++) {
    $cmtfile = $cmtfiles[$ievent-1]; chomp($cmtfile);
    $eid = $cmtfile;
    #$file0 = `basename $cmtfile`; ($junk,$eid) = split("_",$file0); chomp($eid);
    print "\n $ievent, Event $eid";
  }
}
#die("testing");

# write the C-shell script to file
if($idata==1 && $isyn==1) {
   $cshfile = "process_data_and_syn_${ftag}.csh";
} elsif($idata==1 && $isyn==0) {
   $cshfile = "process_data_${ftag}.csh";
} elsif($idata==0 && $isyn==1) {
   $cshfile = "process_syn_${ftag}.csh";
} else {
   die("check idata and isyn\n");
}
print "\nWriting to $cshfile ...\n";
open(CSH,">$cshfile");

if($iprocess==0) {
  $sfile = "nprocess_syn";
  open(SYN,">$sfile");
}

$imin = 1; $imax = $ncmt;  # default
#$imin = 1; $imax = 10;
$imin = 76; $imax = $imin;

#----------------------------------------------------------------------

#foreach $file (@cmtfiles) {
#for ($ievent = 1; $ievent <= $ncmt; $ievent++) {
for ($ievent = $imin; $ievent <= $imax; $ievent++) {

  # get the event ID
  if ($ilist == 1) {
    $eid = $cmtfiles[$ievent-1]; chomp($eid);
    $cmtfile = "$dirCMT/CMTSOLUTION_$eid";

  } else {
    $cmtfile = $cmtfiles[$ievent-1];
    $file0 = `basename $cmtfile`; ($junk,$eid) = split("_",$file0); chomp($eid);
  }
  print "$ievent, $imax, Event $eid\n";
  print CSH "echo $ievent, $imax, Event $eid\n";

  # data and synthetics directories
  $dirsyn = "${dirsyn0}/${eid}";
  $dirdat = "${dirdat0}/${eid}";
  $dirdat_pro_1 = "${dirdat0}/${eid}/$pdir";
  $dirsyn_pro_1 = "${dirsyn0}/${eid}/$pdir";
  $dirdat_pro_2 = "${dirdat_pro_1}/$pdirbpass";
  $dirsyn_pro_2 = "${dirsyn_pro_1}/$pdirbpass";

  # cut times files
  $cutfile_dat      = "$dirdat/${eid}_dat_cut";
  $cutfile_dat_done = "${cutfile_dat}_done";
  $cutfile_syn      = "$dirdat/${eid}_syn_cut";       # note: data directory
  $cutfile_syn_done = "$dirsyn/${eid}_syn_cut_done";  # note: syn directory

  # optional -- delete pre-processed directories
  #print CSH "rm -rf $dirsyn/PRO*\n";
  #print CSH "rm -rf $dirdat/PRO*\n";

  #----------------------------------------------------------------------
  # PROCESSING PART 0: check the number of processed synthetic files for each event

  if ($iprocess == 0) {
     if(-e ${dirsyn_pro_1}) {
       ($nfile,undef,undef) = split(" ",`ls -1 ${dirsyn_pro_1}/* | wc`);
       print SYN "$eid $nfile\n";
     }
  }

  #----------------------------------------------------------------------
  # PROCESSING PART 1: assign SAC headers, interpolate, and pick P and S arrivals (based on a 1D socal model)

  if ($iprocess == 1) {

    # synthetics -- this will convolve with the source half-duration (prior to interpolating)
    if ($isyn == 1) {

      if (-e $dirsyn) {

        # this block is required if you are also converting ASCII files to SAC
 	($nsemd,undef,undef) = split(" ",`ls -1 $dirsyn/*semd | wc`); # number of ascii files
 	($nsacd,undef,undef) = split(" ",`ls -1 $dirsyn/*sac | wc`); # number of sac files
 	print "nsemd = $nsemd, nsacd = $nsacd\n";
 	if ($nsemd*$nsacd > 0) {
 	  if ($nsemd == $nsacd) {
 	    print CSH "rm $dirsyn/*semd\n";
 	    print CSH "sleep 5s\n";
 	  } else {
 	    die("error in converting to SAC files\n");
 	  }
 	}

	if (not -e ${dirsyn_pro_1}) {
	  print CSH "cd $dirsyn\n";
          # converting to SAC is so slow on the cluster -- even in parallel -- that we must do it here
          if($nsemd > 0 && $nsacd==0) {
	     print CSH "process_trinet_syn_new.pl -m $cmtfile -a $stafile *semd\n";
	     print CSH "sleep 5s\n";
	  } else {
	    if ($chdur == 1) {
	      print CSH "process_trinet_syn_new.pl -S -m $cmtfile -h -a $stafile -s $sps -p -d $pdir -x ${syn_ext} *.${syn_suffix0}\n";
	    } else {
	      print CSH "process_trinet_syn_new.pl -S -m $cmtfile -a $stafile -s $sps -p -d $pdir -x ${syn_ext} *.${syn_suffix0}\n";
	    }
          }
	} else {
	  print "dir ${dirsyn_pro_1} already exists\n";
	}
      } else {
	print "$dirsyn does not exist\n";
      }
    }	  # isyn

    # data
    if ($idata == 1) {
      if (-e $dirdat) {
	if (not -e ${dirdat_pro_1}) {
	print CSH "cd $dirdat\n";
        #print CSH "mv $pdir ${pdir}_OLD\n";
	#print CSH "\\rm -rf $pdir; mkdir $pdir\n";
	print CSH "process_cal_data.pl -m $cmtfile -p -s $sps -d $pdir -x ${dat_ext} *.${dat_suffix0}\n";
	} else {
	  print "dir ${dirdat_pro_1} already exists\n";
	}
      } else {
	print "$dirdat does not exist\n";
      }
    }	  # idata

  }  # iprocess = 1

  #----------------------------------------------------------------------
  # PROCESSING PART 2: getting the cut times for the records

  if ($iprocess == 2) {

    # BOTH the INITIALLY processed synthetics and data directories must exist,
    # even if you only want to process synthetics.

    if ( (not -e ${dirdat_pro_1}) || (not -e ${dirsyn_pro_1}) ) {
      print "--> dirdat ${dirdat_pro_1} and dirsyn ${dirsyn_pro_1} do not both exist\n";

    } elsif ( ((-f $cutfile_syn) || (-f $cutfile_syn_done)) || ((-f $cutfile_dat) || (-f $cutfile_dat_done)) ) {

      if (-f $cutfile_syn) {
	print "cutfile_syn ${cutfile_syn} already exists\n";
      }
      if (-f $cutfile_syn_done) {
	print "cutfile_syn_done ${cutfile_syn_done} already exists\n";
      }
      if (-f $cutfile_dat) {
	print "cutfile_dat ${cutfile_dat} already exists\n";
      }
      if (-f $cutfile_dat_done) {
	print "cutfile_dat_done ${cutfile_dat_done} already exists\n";
      }
      print "--> you are probably ready for cutting or bandpassing...\n";

    } else {

       # check that the number of unprocessed files matches the number of processed files
       $ns1 = `ls -1 ${dirsyn}/*.sac | wc | awk '{print \$1}'`; chomp($ns1);
       $nd1 = `ls -1 ${dirdat}/*.sac | wc | awk '{print \$1}'`; chomp($nd1);
       $ns2 = `ls -1 ${dirsyn_pro_1}/*.sac.${smodel} | wc | awk '{print \$1}'`; chomp($ns2);
       $nd2 = `ls -1 ${dirdat_pro_1}/*.sac.d | wc | awk '{print \$1}'`; chomp($nd2);
       print "-- ndata1 $nd1 -- ndata2 $nd2 -- nsyn1 $ns1 -- nsyn2 $ns2 --\n";
       if( ($nd1 != $nd2) || ($ns1 != $ns2)) {
          print "mismatch of expected records\n";
          die("RESOLVE MISMATCH\n");
       }

      print "\nWriting to cutfiles ...\n";
      open(CUTDAT,">${cutfile_dat}");
      open(CUTSYN,">${cutfile_syn}");

      # grab all the initially processed DATA files
      @files = glob("${dirdat_pro_1}/*");
      $nfile = @files;
      print "\n $nfile data files to line up with synthetics\n";

      foreach $datfile (@files) { 
	# read the sac headers -- network, station, component
	(undef,$net,$sta,$chan) = split(" ",`saclst knetwk kstnm kcmpnm f $datfile`);
	$comp = `echo $chan | awk '{print substr(\$1,3,1)}'`;
	chomp($comp);

	# synthetics are always BH_ component
	$synfile = "${dirsyn_pro_1}/${sta}.${net}.BH${comp}.${syn_suffix}";

	# if the synthetic file exists, then go on
	if (-f $synfile) {
          # only list files that are pairs
          # only base name is used for syn file
          $synfile_base = `basename $synfile`; chomp($synfile_base);
	  print "$datfile $synfile_base\n";

	  # get info on data and synthetics
	  (undef,$bd,$ed,$deltad,$nptd) = split(" ",`saclst b e delta npts f $datfile`);
	  (undef,$bs,$es,$deltas,$npts) = split(" ",`saclst b e delta npts f $synfile`);
	  $tlend = $ed - $bd;
	  $tlens = $es - $bs;
    
          # if the end time of the data is less than the origin time, STOP and move the record to REJECTED
          if ($ed < 0) {
             close(CUTDAT); close(CUTSYN); `rm ${cutfile_dat} ${cutfile_syn}`;
             print "move unprocessed data file to REJECTED\n";
             print "end time is $ed\n";
             die("data record ends before 0 -- REJECTED\n")
          }

	  # dt should be the same for both records ALREADY
	  if (log($deltad/$deltas) > 0.01) {
            close(CUTDAT); close(CUTSYN); `rm ${cutfile_dat} ${cutfile_syn}`;
	    print "$datfile $synfile\n";
	    print "DT values are not close enough: $deltad, $deltas\n";
	    die("fix the DT values\n");
	  } else {
	    $dt = $deltad;
	  }

	  # determine the cut times for the records
          $b0 = $bd;
	  #if ($bd < $bs) {$b0 = $bd} else {$b0 = $bs}
	  if ($b0 < $bmin) {$b0 = $bmin}
	  if ($ed < $es) {$e0 = $ed} else {$e0 = $es}

	  $b = $b0;
	  $tlen0 = $e0 - $b;
	  $tlen = $tfac * $tlen0; # extend record length (if desired)
	  $e = $b0 + $tlen;
	  $npt = int( ($e-$b)/$dt );

          # print the cut times if the begin time is before the end time
          if ($b < $e) {
	    print CUTDAT "$datfile $b $e $npt $dt\n";
	    print CUTSYN "$synfile_base $b $e $npt $dt\n";
	  }

	  if (0==1) {
	    print "\n Data : $bd $ed $deltad $nptd -- $tlend";
	    print "\n Syn  : $bs $es $deltas $npts -- $tlens";
	    print "\n b0, e0, tlen0 : $b0, $e0, $tlen0 ";
	    print "\n   b : $bd, $bs, $b ";
	    print "\n   e : $ed, $es, $e ";
	    print "\n npt : $nptd, $npts, $npt ";
	    print "\n $tlen = $tfac * ($e0 - $b)";
	    print "\n $e = $b0 + $tlen \n";
	  }
	}			# if synfile exists
      }				# for all data files
      print "\n Done making cutfile $cutfile_dat\n";
      print "\n Done making cutfile $cutfile_syn\n";
      close(CUTDAT);
      close(CUTSYN);
 
    }				# dirdat_pro_1 and dirsyn_pro_1 exist
  }				# iprocess=2

  #----------------------------------------------------------------------
  # PROCESSING PART 3: cutting records and padding zeros

  if ($iprocess == 3) {

    if ($idata == 1) {
      if (-e ${dirdat_pro_2}) {
	print "--> dirdat ${dirdat_pro_2} already exists\n";

      } else {
	if (-f $cutfile_dat_done) {
	  print "cutfile_dat_done ${cutfile_dat_done} already exists\n";

	} else {

	  if (not -f $cutfile_dat) {
	    print "cutfile_dat $cutfile_dat does not exist -- try iprocess = 2\n";

	  } else {

	    # read cut file
	    open(IN,"${cutfile_dat}"); @lines = <IN>; close(IN); $nlines = @lines;

	    $sacfile = "sacdat_${ftag}.mac";
	    `echo echo on > $sacfile`;
	    `echo readerr badfile fatal >> $sacfile`;

	    for ($j = 1; $j <= $nlines; $j++) {

	      $line = $lines[$j-1]; chomp($line);
	      ($datfile,$b,$e,$npt,$dt) = split(" ",$line);
	      print "$j out of $nlines -- $datfile\n";
	      #print "-- $datfile -- $synfile -- $b -- $e -- $npt -- $dt -- \n";

	      # cut records and fill zeros
	      `echo r $datfile >> $sacfile`;
	      `echo cuterr fillz >> $sacfile`;
	      `echo "cut $b n $npt" >> $sacfile`;
	      `echo r $datfile >> $sacfile`;
	      `echo cut off >> $sacfile`;
	      `echo w over >> $sacfile`;
	    } 
	    `echo quit >> $sacfile`;

	    # KEY: execute SAC command
	    `sac $sacfile`;
	    `sleep 5s`;
	    `rm $sacfile`;
	    print "\n Done cutting pre-processed data files\n";

	    # rename cut file in data directory
	    `mv ${cutfile_dat} ${cutfile_dat_done}`;

	  }			# cutfile exist
	}			# cutfile_done exist
      }				# bandpass dir exist
    }				# idata

    #------------------

    if ($isyn == 1) {
      if ( (-e ${dirsyn_pro_2}) || (not -e ${dirsyn_pro_1}) ) {
        if ( -e ${dirsyn_pro_2}) {print "--> dirsyn ${dirsyn_pro_2} already exists\n";}
        if ( not -e ${dirsyn_pro_1} ) {print "--> dirsyn ${dirsyn_pro_1} does not exist\n";}

      } else {
	if (-f $cutfile_syn_done) {
	  print "cutfile_syn_done ${cutfile_syn_done} already exists\n";

	} else {

	  if (not -f $cutfile_syn) {
	    print "cutfile_syn $cutfile_syn does not exist -- try iprocess = 2\n";

	  } else {

	    # read cut file
	    open(IN,"${cutfile_syn}"); @lines = <IN>; close(IN); $nlines = @lines;

	    $sacfile = "sacsyn_${ftag}.mac";
	    `echo echo on > $sacfile`;
	    `echo readerr badfile fatal >> $sacfile`;

	    for ($j = 1; $j <= $nlines; $j++) {

	      $line = $lines[$j-1]; chomp($line);
	      ($synfile_base,$b,$e,$npt,$dt) = split(" ",$line);

              # KEY: use the suffix of the present model -- the cut file may have been
              #      generated using a different set of synthetics
              #      THIS FILE FORMAT MIGHT DIFFER FOR DIFFERENT USERS
              #      HERE: PHOB.NC.BHE.semd.sac.m12
              ($tag1,$tag2,$tag3,$tag4,$tag5,$tag6) = split("\\.",$synfile_base);
              $synfile_base_new = "${tag1}.${tag2}.${tag3}.${tag4}.${tag5}.${smodel}";

              # KEY: indicate the base directory
              $synfile = "${dirsyn_pro_1}/${synfile_base_new}";
              if (not -f $synfile) {
                 print "synfile $synfile does not exist\n";
                 #die("synfile $synfile does not exist");

              } else {
		print "$j out of $nlines -- $synfile\n";
		#print "-- $datfile -- $synfile -- $b -- $e -- $npt -- $dt -- \n";

		# cut records and fill zeros
		`echo r $synfile >> $sacfile`;
		`echo cuterr fillz >> $sacfile`;
		`echo "cut $b n $npt" >> $sacfile`;
		`echo r $synfile >> $sacfile`;
		`echo cut off >> $sacfile`;
		`echo w over >> $sacfile`;
              }
	    } 
	    `echo quit >> $sacfile`;

	    # KEY: execute SAC command
	    `sac $sacfile`;
	    `sleep 5s`;
	    `rm $sacfile`;
	    print "\n Done cutting pre-processed syn files\n";

	    # copy syn cut file into syn directory
	    `cp ${cutfile_syn} ${cutfile_syn_done}`;

	  }			# cutfile exist
	}			# cutfile_done exist
      }				# bandpass dir exist
    }				# isyn
 
  }				# iprocess = 3


  #----------------------------------------------------------------------
  # PROCESSING PART 4: bandpass

  if ($iprocess == 4) {

    if ($isyn == 1) {
      if (-e ${dirsyn_pro_2}) {
	print "--> dirsyn ${dirsyn_pro_2} already exists\n";

      } else {
	if (not -f ${cutfile_syn_done}) {
	  print "cutfile_syn_done ${cutfile_syn_done} does not exist\n";

	} else {
	  if (not -e ${dirsyn_pro_1}) {
	    print "${dirsyn_pro_1} does not exist\n";

	  } else {
	    print CSH "cd ${dirsyn_pro_1}\n";
	    #print CSH "\\rm -rf $pdirbpass\n";
	    print CSH "process_trinet_syn_new.pl -S -t $Trange -d $pdirbpass -x $Ttag *.${syn_suffix} \n";
	    print CSH "cd $pdirbpass\n";
	    print CSH "rotate.pl *E.${syn_suffix}.${Ttag}\n";
	  }
	}
      }
    }				# isyn

    #-----------

    if ($idata == 1) {
      if (-e ${dirdat_pro_2}) {
	print "--> dirdat ${dirdat_pro_2} already exists\n";

      } else {
	if (not -f ${cutfile_dat_done}) {
	  print "cutfile_dat_done ${cutfile_dat_done} does not exist\n";

	} else {
	  if (not -e ${dirdat_pro_1}) {
	    print "${dirdat_pro_1} does not exist\n";

	  } else {
	    $ofile1 = "${eid}_ofile";
	    $ofile2 = "${eid}_no_pz_files";
	    $odir = "${dirdat0}/CHECK_DIR";
	    print CSH "mkdir -p $odir\n";

	    print CSH "cd ${dirdat_pro_1}\n";
	    #print CSH "\\rm -rf $pdirbpass\n";
	    print CSH "process_cal_data.pl -i none -t $Trange -d $pdirbpass -x $Ttag *.${dat_suffix} > $ofile1\n";
	    print CSH "grep skip $ofile1 | awk '{print \$7}' > $ofile2\n";
	    print CSH "\\cp $ofile1 $ofile2 $odir\n";

	    print CSH "cd $pdirbpass\n";
	    print CSH "rotate.pl *E.${dat_suffix}.${Ttag}\n";
	  }
	}
      }
    }				# idata

  }				# iprocess==4

}  # END OF LOOP OVER EVENTS

if($iprocess==0) {close(SYN);}

#======================
print CSH "echo done with $cshfile\n";
close(CSH);
print "closing $cshfile\n";
if(($iprocess==1) || ($iprocess==4)) {print "csh -f $cshfile\n";}
if($iexecute==1) {system("csh -f $cshfile");}

print "\n ";
#=================================================================
