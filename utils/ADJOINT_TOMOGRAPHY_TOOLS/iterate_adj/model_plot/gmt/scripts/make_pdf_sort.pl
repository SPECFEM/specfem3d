#!/usr/bin/perl -w

#==========================================================
#
#  make_pdf_sort.pl
#  Carl Tape
#  28-March-2009
#
#  This sorts all PDF files by
#    1  azimuth or distance
#    2  of stations or of events
#
#  EXAMPLES:
#    make_pdf_sort.pl 1/1/0 T006_T030   # for each station, sort events by azimuth
#    make_pdf_sort.pl 1/0/1 T006_T030   # for each event, sort station by azimuth
#    make_pdf_sort.pl 1/1/1 T006_T030   # both at once
#
#==========================================================

if (@ARGV < 2) {die("Usage: make_pdf_sort.pl iaz/istation/ievent Ttag\n")}
($ibool,$Ttag) = @ARGV;

($iaz,$istation,$ievent) = split("\\/",$ibool);
if($iaz == 1) {$sortlab = "az"} else {$sortlab = "dist"}

#================================================================
# USER INPUT

#$emin = 1;         # number of events recorded by a station required to make a composite PDF
#$smodel = "m16";
#$Ttag = "T003_T030";
#$isort = 2;         # sorted stations by distance (1) or azimuth (2)

# full stations file
$file_stations = "/home/carltape/gmt/stations/seismic/Matlab_output/STATIONS_CALIFORNIA_TOMO_INNER_gmt";
if (not -f ${file_stations}) {die("check if file_stations ${file_stations} exist or not\n")}
open(IN,"${file_stations}"); @stalines0 = <IN>; close(IN); $nrec = @stalines0;

# directory containing sorted stations lists
# --> to make these files, see /ADJOINT_TOMO/iterate_adj/UTILS/station_lists/
$dir0 = "/home/carltape/results/SOURCES";
$edir = "${dir0}/EID_STATION_LISTS";

# list of all possible event IDs
$file_eid = "${dir0}/socal_16/EIDs_only_loc";
if (not -f ${file_eid}) {die("check if ${file_eid} exist or not\n")}
open(IN,"${file_eid}"); @elines0 = <IN>; close(IN);
$nevent = @elines0;

#================================================================

print "\n nrec -- $nrec -- nevent -- $nevent --\n";

if (0==1) {
  for ($ik = 1; $ik <= $nrec; $ik = $ik+1) {
    ($stalon,$stalat,$station,$network,undef,undef) = split(" ",$stalines0[$ik-1]);
    print "$ik out of $nrec : station $station.$network\n";
  }
  die("TESTING: listing only the stations\n");
}

#--------------------------------------

if ($istation==1) {

  $imin = 1; $imax = $nrec;
  #$imin = 120; $imax = 200;
  #$imin = 120; $imax = $imin;

  # NOTE: make sure that the pdcat executable is there
  @pdcat = "/home/carltape/bin/pdcat -r";
  $k = 1;
  $ofile = "ALL_xc_and_seis_${Ttag}_station.pdf";

  $ffile = "sort_stations_by_${sortlab}";
  open(OUT,">$ffile");

  # loop over all the stations in the sorted list
  for ($ik = $imin; $ik <= $imax; $ik = $ik+1) {

    # get station name from list
    ($stalon,$stalat,$station,$network,undef,undef) = split(" ",$stalines0[$ik-1]);
    print "$ik out of $imin to $imax : station $station.$network\n";

    # get list of all possible events sorted by azimuth or distance
    # EXAMPLE: EIDS_by_az_from_LT2.NP
    $event_sort = "${edir}/EIDS_by_${sortlab}_from_${station}.${network}";
    if (not -f "${event_sort}") {die("check if event_sort ${event_sort} exist or not\n");}
    open(IN,"${event_sort}"); @elines = <IN>; close(IN); $nevent = @elines;
    print "$nevent events in the list \n";

    # loop over all possible events
    for ($j = 1; $j <= $nevent; $j = $j+1) {

      # get the event ID in the sorted list
      ($eid,$elon,$elat,undef,undef,undef,undef,undef) = split(" ",$elines[$j-1]);

      # EXAMPLE: xc_and_seis_14418600_HEC_CI_T006_T030.pdf
      $ftag = "${eid}_${station}_${network}_${Ttag}";
      @files = glob("*${ftag}.pdf");
      $numf = @files;
      if ($numf == 0) {
  #print "$j -- $eid --> no pdf file exists\n";
      } elsif ($numf == 1) {
  #print "$j -- $eid --> pdf file exists\n";
  $pdffile = $files[0]; chomp($pdffile);
  $pdcat[$k] = $pdffile;
        $k = $k+1;
        print OUT "$ftag\n";
      } else {
  print "$eid\n";
  die("more than one pdf file exists\n");
      }
    }       # for
  }

  # make composite PDF file
  print "output file is $ofile\n";
  $pdcat[$k] = "$ofile";
  print "@pdcat\n";
  `@pdcat`;
  `sleep 3s`;
}

#=================================================================

if ($ievent==1) {

  $jmin = 1; $jmax = $nevent;
  #$jmin = 120; $jmax = 200;
  #$jmin = 120; $jmax = $jmin;

  # NOTE: make sure that the pdcat executable is there
  @pdcat = "/home/carltape/bin/pdcat -r";
  $k = 1;
  $ofile = "ALL_xc_and_seis_${Ttag}_event.pdf";

  $ffile = "sort_events_by_${sortlab}";
  open(OUT,">$ffile");

  # loop over all possible events
  for ($j = $jmin; $j <= $jmax; $j = $j+1) {

    $eid = $elines0[$j-1]; chomp($eid);
    print "$j out of $jmin to $jmax : eid $eid\n";

    # get list of all possible stations sorted by azimuth or distance
    # EXAMPLE: STATIONS_by_az_from_9818433
    #            U10A.TA -116.3297   36.4193     9818433 -117.7840   33.9133  308.4102   24.9641
    #             AMD.SN -116.2809   36.4526     9818433 -117.7840   33.9133  313.6531   25.3876
    $station_sort = "${edir}/STATIONS_by_${sortlab}_from_${eid}";
    if (not -f "${station_sort}") {die("check if station_sort ${station_sort} exist or not\n");}
    open(IN,"${station_sort}"); @stalines = <IN>; close(IN); $nrec = @stalines;
    print "$nrec stations in the list \n";

    # loop over all the stations in the sorted list
    for ($ik = 1; $ik <= $nrec; $ik = $ik+1) {

      # get station name from list
      ($stanet,undef,undef,undef,undef,undef,undef,undef,undef,undef) = split(" ",$stalines[$ik-1]);
      ($station,$network) = split("\\.",$stanet);
      #print "$ik out of $nrec : station $station.$network\n";

      # EXAMPLE: xc_and_seis_14418600_HEC_CI_T006_T030.pdf
      $ftag = "${eid}_${station}_${network}_${Ttag}";
      @files = glob("*${ftag}.pdf");
      $numf = @files;
      if ($numf == 0) {
  #print "$j -- $eid --> no pdf file exists\n";
      } elsif ($numf == 1) {
  #print "$j -- $eid --> pdf file exists\n";
  $pdffile = $files[0]; chomp($pdffile);
  $pdcat[$k] = $pdffile;
        $k = $k+1;
        print OUT "$ftag\n";
      } else {
  print "$eid\n";
  die("more than one pdf file exists\n");
      }
    }       # for
  }

  # make composite PDF file
  print "output file is $ofile\n";
  $pdcat[$k] = "$ofile";
  print "@pdcat\n";
  `@pdcat`;
}

#=================================================================
