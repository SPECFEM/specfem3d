#!/usr/bin/perl -w
#
#==========================================================
#
#  make_pdf_by_station.pl
#
#  This sorts the output PDF files in PLOTS into order by
#    1  distance
#    2  azimuth
#
#  EXAMPLE (from /home/carltape/results/MEASUREMENTS/m16/PDF_ALL_T006_T030):
#    make_pdf_by_station.pl
#
#==========================================================

#if (@ARGV < 3) {die("Usage: make_pdf_by_station.pl xxx\n")}
#($isort,$otag) = @ARGV;

#================================================================
# USER INPUT

$emin = 1;         # number of events recorded by a station required to make a composite PDF
$smodel = "m16";
$Ttag = "T003_T030";
$isort = 2;         # sorted stations by distance (1) or azimuth (2)

# full stations file
$file_stations = "/home/carltape/gmt/stations/seismic/Matlab_output/STATIONS_CALIFORNIA_TOMO_INNER_gmt";
if (not -f ${file_stations}) {die("check if file_stations ${file_stations} exist or not\n")}
open(IN,"${file_stations}"); @stalines = <IN>; close(IN); $nrec = @stalines;

# directory containing sorted stations lists
# (See /ADJOINT_TOMO/iterate_adj/UTILS/station_lists/ )
$edir = "/home/carltape/results/SOURCES/EID_STATION_LISTS";

# subset list of events that you want to use
#$file_eid_sub = "/home/carltape/results/SOURCES/socal_16/EIDs_only_loc";
$file_eid_sub = "/home/carltape/results/EID_LISTS/eids_simulation";
if (not -f ${file_eid_sub}) {die("check if ${file_eid_sub} exist or not\n")}

# list of all possible event IDs
#$file_eid = "/home/carltape/results/EID_LISTS/syn_run_${smodel}";
#if (not -f ${file_eid}) {die("check if ${file_eid} exist or not\n")}
#open(IN,"${file_eid}"); @elines = <IN>; close(IN);
#$nevent = @elines;
#print "\n $nevent events in the list \n";

if($isort == 1) {$sortlab = "dist"} else {$sortlab = "az"}

#================================================================

print "\n $nrec \n";

if (0==1) {
  for ($ik = 1; $ik <= $nrec; $ik = $ik+1) {
    ($stalon,$stalat,$station,$network,undef,undef) = split(" ",$stalines[$ik-1]);
    print "$ik out of $nrec : station $station.$network\n";
  }
  die("TESTING: listing only the stations\n");
}

#--------------------------------------

$imin = 1; $imax = $nrec;
#$imin = 120; $imax = $imin;

# loop over all the stations in the sorted list
for ($ik = $imin; $ik <= $imax; $ik = $ik+1) {

  # NOTE: make sure that the executable is there
  @pdcat = "/home/carltape/bin/pdcat -r"; $k = 1;

  # get station name from list
  ($stalon,$stalat,$station,$network,undef,undef) = split(" ",$stalines[$ik-1]);
  print "$ik out of $nrec : station $station.$network\n";

  # get list of all possible events sorted by azimuth or distance
  # EXAMPLE: EIDS_by_az_from_LT2.NP
  $event_sort = "${edir}/EIDS_by_${sortlab}_from_${station}.${network}";
  if (not -f "${event_sort}") {die("check if event_sort ${event_sort} exist or not\n")}
  open(IN,"${event_sort}"); @elines = <IN>; close(IN); $nevent = @elines;
  print "$nevent events in the list \n";

  # output file
  $ofile = "${station}_${network}_${Ttag}_${smodel}_ALL.pdf";

  if (-f $ofile) {
    print "--> $ofile already exists\n";

  } else {
    # loop over all possible events
    for ($j = 1; $j <= $nevent; $j = $j+1) {

      # get the event ID in the sorted list
      ($eid,$elon,$elat,undef,undef,undef,undef,undef) = split(" ",$elines[$j-1]);

      # check if event is in the subset list
      ($nmatch,undef,undef) = split(" ",`grep $eid ${file_eid_sub} | wc`);

      if ($nmatch == 1) {
  # EXAMPLE: 9828889_T006_T030_GSC_CI_m16_cc_win_adj.pdf
  @files = glob("${eid}_${Ttag}_${station}_${network}_${smodel}*pdf");
  $numf = @files;
  if ($numf == 0) {
    #print "$j -- $eid --> no pdf file exists\n";

  } elsif ($numf == 1) {
    #print "$j -- $eid --> pdf file exists\n";
    $pdffile = $files[0]; chomp($pdffile);
    $pdcat[$k] = $pdffile; $k = $k+1;

  } else {
    print "$eid\n";
    die("more than one pdf file exists\n");
  }
      }
    }       # for

    # if there is at least one file, then make the composite PDF
    if ($k > $emin+1) {
      print "output file is $ofile\n";
      $pdcat[$k] = "./$ofile";
      print "@pdcat\n";
      `@pdcat`;     # execute
      `sleep 5s`;

    } else {
      print "--> Fewer than $emin events out of $nevent for this station\n";
    }

  }

}

#=================================================================
