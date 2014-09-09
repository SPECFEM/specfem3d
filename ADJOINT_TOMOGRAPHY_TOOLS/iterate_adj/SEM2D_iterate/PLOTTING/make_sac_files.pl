#!/usr/bin/perl -w

#==========================================================
#
#  make_sac_files.pl
#  Carl Tape
#  23-May-2006
#
#  This script converts ascii seismograms to sac files,
#  and also includes the src and rec lat-lon.
#
#  This generates the STATIONS file, which is needed for multi_plot_rs.pl
#
#  EXAMPLE (execute from main directory):
#    scripts/make_sac_files.pl OUTPUT/ dat 1 48.0
#    scripts/make_sac_files.pl OUTPUT/ syn 1 48.0
#    scripts/make_sac_files.pl OUTPUT/ stfadj 1 48.0
#
#    scripts/make_sac_files.pl OUTPUT/run_3100/event_005 dat 1 48.0
#    scripts/make_sac_files.pl OUTPUT/run_3100/event_005 syn 1 48.0
#    scripts/make_sac_files.pl OUTPUT/run_3100/event_005 stfadj 1 48.0
#
#    scripts/make_sac_files.pl OUTPUT/run_4600/event_005 stfadj 1 48.0; scripts/make_sac_files.pl OUTPUT/run_4600/event_005 syn 1 48.0; scripts/make_sac_files.pl OUTPUT/run_4600/event_005 dat 1 48.0
#
#    scripts/make_sac_files.pl OUTPUT_body/run_5000/event_001 dat 1 16.0
#    scripts/make_sac_files.pl OUTPUT_body/run_5000/event_001 dat 2 16.0
#    scripts/make_sac_files.pl OUTPUT_body/run_5000/event_001 dat 3 16.0
#
#    scripts/make_sac_files.pl OUTPUT_body/run_5000/event_001 dat 3 16.0
#
#==========================================================

if (@ARGV < 3) {die("Usage: make_sac_files.pl datadir file_prefix comp tshift\n");}
($datadir,$file_prefix,$ncomp,$tshift) = @ARGV;

# origin time for records
$mtshift = -$tshift;

for($i = 1; $i <= $ncomp; $i++) {

# get ascii seismograms
@files = glob("$datadir/${file_prefix}*$i");
$numf = 1 + $#files;
print "\n $numf files in $datadir with label $datafile \n";
if (@files == 0) {die("No data files available to plot\n");}

# create source and receiver lists from wave2d.f90 output file [sr.txt]
$src_file = "src_list";
$rec_file = "rec_list";
system("awk '\$1 == \"S\" {print \$2,\$3}' ${datadir}/sr.txt > $src_file");
system("awk '\$1 == \"R\" {print \$2,\$3}' ${datadir}/sr.txt > $rec_file");
open(IN,"$src_file"); @src = <IN>;
open(IN,"$rec_file"); @rec = <IN>;

# get the (first) source
($evlo,$evla)=split(" ",$src[0]);

# make station file for plot_rs.pl
$station_file = "STATIONS";
open(OUT, ">$station_file");

$net = "CI";
$ref_time = 0;

#foreach $file (@files) {
for($k = 0; $k < $numf; $k++){

   # labels for data files
   $dir_file_ascii = $files[$k];
   ($dir,$file_suffix) = split($file_prefix,$dir_file_ascii);
   $file_ascii = $file_prefix.$file_suffix;
   #print "\n $dir_file_ascii $file_ascii";

   # extract station name (assume <1000 receivers)
   ($lab,$sta_temp,$comp) = split("_",$file_ascii);
   $sta = substr($sta_temp, 2, 3);
   #@chars = split(//,$sta_temp); $sta = join("",@chars[2..4]);
   #print "\n $file_ascii to $sta_temp to $sta. \n";

   # convert ascii file to sac file
   #$file_sac = $lab."."."$sta".".".$comp.".sac";
   system("ascii2sac.csh $dir_file_ascii");
   $dir_file_sac = "$datadir/$sta.$net.$comp.sac$lab";
   system("mv $dir_file_ascii.sac $dir_file_sac");

   # extract station longitude-latitude
   # THIS IS NOT A FOOL-PROOF LINK, SINCE THE GLOB COMMAND ABOVE
   # MUST SELECT THE FILES IN THE SAME ORDER IN WHICH THE STATION
   # LIST WAS MADE.
   ($stlo,$stla)=split(" ",$rec[$k]);

   # update labels on the sac file
   #system("sacch b $mtshift knetwk $net kcmpnm $comp kstnm $sta evlo $evlo evla $evla stlo $stlo stla $stla f $dir_file_sac");
   system("sacch t0 $ref_time knetwk $net kcmpnm $comp kstnm $sta evlo $evlo evla $evla stlo $stlo stla $stla f $dir_file_sac");

   print OUT "$sta\n";
}
close(OUT)

}
#==================================
