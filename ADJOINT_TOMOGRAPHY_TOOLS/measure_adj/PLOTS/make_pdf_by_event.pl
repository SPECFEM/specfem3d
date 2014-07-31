#!/usr/bin/perl -w
#
#==========================================================
#  make_pdf_by_event.pl
#
#  This sorts the output PDF files in PLOTS into order by
#    1  distance
#    2  azimuth
#
#  EXAMPLE (from PLOTS): make_pdf_by_event.pl 3 mt_cc_all
#
#==========================================================

#if (@ARGV < 3) {die("Usage: make_pdf_by_event.pl xxx\n")}
#($isort,$otag) = @ARGV;

$station_sort = "STATIONS_sort";
if (not -f "${station_sort}") {die("check if station_sort ${station_sort} exist or not\n")}
open(IN,"${station_sort}"); @stalines = <IN>; close(IN); $nrec = @stalines;

print "\n $nrec \n";

# NOTE: make sure that the executable is there
@pdcat = "/home/carltape/bin/pdcat -r"; $k = 1;

# loop over all the stations in the sorted list
for ($ik = 1; $ik <= $nrec; $ik = $ik+1) {

  ($stanet,$stalon,$stlat,undef,undef,undef,undef,undef) = split(" ",$stalines[$ik-1]);
  ($station,$network) = split("\\.",$stanet);
  #print "$ik out of $nrec : $station $network\n";
  
  @files = glob("*${station}_${network}*pdf");
  $numf = @files;
  if ($numf == 0) {
     #print "$ik out of $nrec : $station $network --> no pdf file exists\n";

  } elsif ($numf == 1) {
     print "$ik out of $nrec : $station $network --> one pdf file exists\n";
     $pdffile = $files[0]; chomp($pdffile);
     $pdcat[$k] = $pdffile; $k = $k+1;

  } else {
     die("more than one pdf file exists\n");
  }

}

print "@pdcat\n";

# construct the output file name from the first file in the list
($flabel,undef) = split("\\.",$pdcat[1]);
($eid,$sTmin,$sTmax,$sta,$net,$smodel,$tag1,$tag2,$tag3) = split("_",$flabel);
$ofile = "${eid}_${sTmin}_${sTmax}_${smodel}_ALL.pdf";
print "output file is $ofile\n";

$pdcat[$k] = "./$ofile";

# execute pdf command
`@pdcat`;
`sleep 3s`;

#-----------------
# OBSOLETE VERSION AS OF 02-OCT-2008

#if (@ARGV < 2) {die("Usage: make_pdf_by_event.pl xxx\n")}
#($isort,$otag) = @ARGV;

#@tags = ("sta","dist","az");
#$tag = $tags[$isort-1];

#$ofile = "${otag}_${tag}.pdf";

#$station_sort = "STATIONS_${tag}";
#if (not -f "${station_sort}") {die("check if station_sort ${station_sort} exist or not\n");}
#@pdffiles = glob("*pdf");
#$npdf = @pdffiles;

#print "\n $npdf files\n";

## loop over all the PDF files
#for ($ik = 1; $ik <= $npdf; $ik = $ik+1) {

#  $pdffile0 = $pdffiles[$ik-1];
#  ($pdffile) = split(" ",`basename $pdffile0`);
#  ($sta,$net,undef,undef,undef) = split("_",$pdffile);
#  $stanet = "$sta.$net";
#  $iline=`grep -n \"$stanet\" ${station_sort} | awk -F: '{print \$1}'`; chomp($iline);
#  print "$ik $stanet $iline\n";

#  $sti = sprintf("%3.3i",$iline);
#  $pdffile_out = "${sti}_${pdffile}";
#  `\\mv $pdffile0 ${pdffile_out}`;
#  #print "move $iline $sti $pdffile0 ${pdffile_out}\n";
#}

## concatenate into one file
#`/home/carltape/bin/pdcat -r *.pdf $ofile`;

#=================================================================
