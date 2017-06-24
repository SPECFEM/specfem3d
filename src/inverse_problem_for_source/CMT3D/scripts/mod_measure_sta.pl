#!/usr/bin/perl

# this program removes the associated measurement files from MEASURE_FILE
# according to a list of stations

if (@ARGV != 2) {die("mod_measure_sta.pl MEASURE syn-file-list\n");}

%synfile=();
open(STA,"$ARGV[1]");
@msyn=<STA>;
for $syn (@msyn) {
  chomp($syn);
  ($sta,$net,$comp)=split(/\./,$syn);
  if ($comp=~/BHZ/) {
    $synfile{"$sta.$net.$comp"} = 1; }
  else {
    $synfile{"$sta.$net.BHR"} = 1; $synfile{"$sta.$net.BHT"} = 1; }
}

@syns=keys %synfile;
print "Missing synthetics: @syns\n";

open(FIN,"$ARGV[0]") || die("Check if $ARGV[0] exists or not\n");

$tmp_file="measure.tmp";
open(FOUT,">$tmp_file");

($nfile) = split(" ",<FIN>);

$nn=0;
for ($i=0;$i<$nfile;$i++) {
  $data[$i] = <FIN>; chomp($data[$i]);
  $syn[$i] = <FIN>; chomp($syn[$i]);
  $nwin[$i] = <FIN>; chomp($nwin[$i]);
   for ($j=0;$j<$nwin[$i];$j++) {
     $tt[$j] = <FIN>;}
   for ($k=0;$k<@syns;$k++) {
     if ($syn[$i] =~/$syns[$k]/) {last;}}
#  print "processing $data[$i] -- $k file\n";
   if ($k==@syns) { # no matching
     $nn++;
     print FOUT "$data[$i]\n$syn[$i]\n$nwin[$i]\n";
     for ($j=0;$j<$nwin[$i];$j++) {print FOUT "$tt[$j]";}
   } else{
     print "Matching: $data[$i] and $syn[$i]\n";}
 }
close(FIN);
close(FOUT);
print "Orignal number of files $nfile; New number of files $nn\n";

system("echo $nn > $ARGV[0]; cat $tmp_file >> $ARGV[0]; rm -f $tmp_file");




